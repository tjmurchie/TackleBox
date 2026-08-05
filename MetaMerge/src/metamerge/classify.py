"""Core merge and DNA-support classification logic for MetaMerge.

This module contains the main ``build_merge`` function, which drives the
per-taxon classification pipeline. The underlying ``classify_status``/
``classify_status_v2`` classifiers live in ``evidence.py``, re-exported here for
backward compatibility:

  1. For each taxon in the MEGAN count matrix:
     a. Compute real-library and blank-library count statistics.
     b. Look up exact Holi/metaDMG damage support for each real library and
        each blank library.
     c. Assess lineage support from other exact-damage-supported taxa in the
        same libraries (optional, conservative).
     d. Determine the QC label for the best exact-match row.
     e. Classify into one of the six DNA-support categories.
     f. Resolve a common name (optional).
     g. Emit a merged record.

DNA-support categories
-----------------------
Very high confidence
    Exact Holi/metaDMG damage support (damage > threshold, significance >
    threshold, N_reads ≥ high_confidence_n_reads_min) + strong conservative
    MEGAN count support + no strong QC caution.

High confidence
    Exact Holi/metaDMG damage support + strong conservative MEGAN count support,
    but without the extra ≥100-read criterion and/or with mild QC caution.

Supported
    Exact damage support and/or strong MEGAN count support and/or low-rank
    lineage-consistent support, but below the confidence thresholds.

Tentative
    Weak or incomplete DNA support that still exceeds the weak-support floor
    and is not blank-associated.

Weak support
    Minimal non-control evidence, below the Tentative floor.

Blank-associated
    Substantial blank overlap and weak real-library support.

Design principles
------------------
- No fuzzy taxon matching.
- Classification is based on DNA evidence only (damage signals, read counts,
  lineage consistency).  Additional lines of evidence such as ecological
  plausibility, biogeographic range, or macrofossil data can be applied as a
  separate interpretive layer on top of these DNA-based categories.
- Each step is documented so the classification can be audited row-by-row.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from .common_names import (
    load_common_name_cache,
    load_common_name_overrides,
    resolve_common_name,
    save_common_name_cache,
)
from .evidence import SourceSignals, TaxonEvidence, classify_status, classify_status_v2
from .fillet_matching import (
    build_fillet_taxonomy_lookup,
    choose_best_fillet_row,
    make_fillet_exact_index,
    row_has_fillet_strong_count_support,
    row_is_fillet_authenticated,
)
from .holi import (
    build_holi_taxonomy_lookup,
    build_lineage_indexes,
    choose_best_exact_row,
    compute_qc_label,
    is_meaningful_low_rank_lineage_support,
    make_holi_exact_index,
    row_has_exact_damage_support,
    summarize_lineage_support_for_taxon,
)
from .utils import STATUS_PRIORITY, normalize_name, normalize_rank, select_best_status


def _broad_group(tax_path: str) -> str:
    """Infer a broad taxonomic group label from a lineage path.

    Used for grouping taxa in heatmap reports.

    Args:
        tax_path: Semicolon-delimited lineage path from metaDMG.

    Returns:
        One of ``"fungi"``, ``"plants"``, ``"animals"``, ``"microbes"``, or
        ``"other"``.
    """
    path_lower = tax_path.lower()
    if "fungi" in path_lower:
        return "fungi"
    if any(
        x in path_lower
        for x in ["viridiplantae", "streptophyta", "embryophyta", "chlorophyta", "tracheophyta"]
    ):
        return "plants"
    if any(
        x in path_lower
        for x in ["metazoa", "animalia", "chordata", "arthropoda", "mollusca",
                   "annelida", "nematoda", "cnidaria"]
    ):
        return "animals"
    if "bacteria" in path_lower or "archaea" in path_lower:
        return "microbes"
    return "other"


def build_merge(
    metadata: pd.DataFrame,
    megan_df: pd.DataFrame,
    holi_df: pd.DataFrame,
    config: dict,
    common_name_overrides: dict | None = None,
    cache_path=None,
    fillet_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Merge MEGAN + Holi/metaDMG (+ optionally Fillet) + metadata into one row
    per taxon.

    This is the project-neutral classification engine.  It classifies taxa
    based on DNA evidence (damage signals, read counts, lineage consistency),
    and -- when ``fillet_df`` is supplied -- Fillet's own independent
    multi-proxy authentication (composite_authenticity/authenticity_tier plus
    eco/pal/fos ecological support lines).

    Backward compatibility: when ``fillet_df`` is ``None`` (the historical
    call signature), every taxon's DNA-support status is produced by calling
    the original, untouched ``classify_status()`` directly (via
    ``evidence.classify_status_v2``'s MEGAN+Holi path) -- this function's
    output is unchanged from before Fillet support existed. This is currently
    the MEGAN-anchored path only: a real Holi+Fillet-only run with no MEGAN
    count matrix at all is not yet supported by this function (tracked in
    FILLET_METAMERGE_INTEGRATION_PLAN.md as the next increment) -- pass an
    empty/dummy ``megan_df`` is not a substitute, this simply isn't wired up
    yet for that case.

    Args:
        metadata: Validated library-linker DataFrame (from ``load_metadata``).
        megan_df: Loaded MEGAN count matrix (from ``load_megan_counts``).
        holi_df: Loaded metaDMG/Holi table (from ``load_holi``).
        config: Full MetaMerge config dict.
        common_name_overrides: Optional dict from ``load_common_name_overrides``.
        cache_path: Optional Path for persisting online common-name lookups.
        fillet_df: Optional loaded Fillet MetaMerge evidence table (from
            ``io.load_fillet``). When supplied, ``metadata`` must have a
            ``fillet_library_name`` column (raises ``ValueError`` otherwise).

    Returns:
        Tuple of ``(merged_df, summary_dict)`` where:
          - ``merged_df`` has one row per MEGAN taxon, sorted by DNA support
            status, MEGAN max count, and scientific name.
          - ``summary_dict`` contains ``n_taxa``, ``status_counts``, and
            ``unmatched_taxa_examples``.
    """
    thresholds  = config["thresholds"]
    exact_by_id, exact_by_name_rank = make_holi_exact_index(holi_df)
    tax_lookup  = build_holi_taxonomy_lookup(holi_df)
    common_name_overrides = common_name_overrides or {}
    common_name_cache = load_common_name_cache(cache_path) if cache_path else {}

    have_fillet = fillet_df is not None
    if have_fillet and "fillet_library_name" not in metadata.columns:
        raise ValueError(
            "build_merge() was given fillet_df but metadata has no "
            "fillet_library_name column -- add one to the linker file."
        )
    if have_fillet:
        fillet_exact_by_id, fillet_exact_by_name_rank = make_fillet_exact_index(fillet_df)
        fillet_tax_lookup = build_fillet_taxonomy_lookup(fillet_df)
        # Per-sample index of every Fillet-authenticated row, regardless of
        # which taxon it's for -- used by the discordance check below to ask
        # "did Fillet authenticate something ELSE, unrelated, for this same
        # library" without re-scanning the whole table per taxon.
        fillet_authenticated_by_sample: dict = defaultdict(list)
        for _row in fillet_df.to_dict(orient="records"):
            if row_is_fillet_authenticated(_row, thresholds):
                fillet_authenticated_by_sample[_row["sample"]].append(_row)
        fillet_authenticated_by_sample = dict(fillet_authenticated_by_sample)

    is_pos_ctrl = metadata.get("is_positive_control",      pd.Series(False, index=metadata.index))
    is_env_ctrl = metadata.get("is_environmental_control", pd.Series(False, index=metadata.index))
    real_meta  = metadata.loc[~metadata["is_negative_control"] & ~is_env_ctrl].copy()
    blank_meta = metadata.loc[ metadata["is_negative_control"]].copy()

    megan_library_cols = metadata["megan_library_name"].tolist()
    real_cols          = real_meta["megan_library_name"].tolist()
    blank_cols         = blank_meta["megan_library_name"].tolist()

    # Pre-build lineage indexes once; reused for every MEGAN taxon.
    _, supported_name_index_by_sample, descendant_index_by_sample = build_lineage_indexes(
        holi_df, thresholds
    )

    records        = []
    unmatched_taxa = []

    for row in tqdm(megan_df.to_dict(orient="records"), desc="Classifying taxa", unit="taxon"):
        tax_id     = row.get("tax_id")
        tax_id_str = row.get("tax_id_str") or ""
        tax_name   = normalize_name(row.get("tax_name"))
        tax_rank   = normalize_rank(row.get("tax_rank"))

        counts      = {lib: float(row.get(lib, 0) or 0) for lib in megan_library_cols}
        real_counts = {lib: counts[lib] for lib in real_cols}
        blank_counts= {lib: counts[lib] for lib in blank_cols}

        real_positive_libs  = [lib for lib, v in real_counts.items()  if v > 0]
        blank_positive_libs = [lib for lib, v in blank_counts.items() if v > 0]

        n_real_positive = len(real_positive_libs)
        max_real_count  = max(real_counts.values())  if real_counts  else 0.0
        max_blank_count = max(blank_counts.values()) if blank_counts else 0.0
        blank_ratio     = (max_blank_count / max_real_count) if max_real_count > 0 else np.nan

        strong_count_support = (
            n_real_positive >= thresholds["strong_count_min_libraries"]
            and max_real_count >= thresholds["strong_count_min_reads"]
        )
        some_count_support = n_real_positive >= 1

        # ── Exact Holi matching for real libraries ──────────────────────────
        exact_rows_real             = []
        exact_damage_supported_real = []

        for _, meta_row in real_meta.iterrows():
            megan_lib = meta_row["megan_library_name"]
            holi_lib  = meta_row["holi_library_name"]
            if counts.get(megan_lib, 0) <= 0:
                continue
            rows = []
            if tax_id_str:
                rows.extend(exact_by_id.get((holi_lib, tax_id_str), []))
            if not rows and tax_name:
                rows.extend(exact_by_name_rank.get((holi_lib, tax_name, tax_rank), []))
                if not rows:
                    rows.extend(exact_by_name_rank.get((holi_lib, tax_name, ""), []))
            best = choose_best_exact_row(rows, thresholds)
            if best:
                augmented = dict(best)
                augmented["megan_library_name"] = megan_lib
                exact_rows_real.append(augmented)
                if row_has_exact_damage_support(best, thresholds):
                    exact_damage_supported_real.append(augmented)

        # ── Exact Holi matching for blank libraries ──────────────────────────
        exact_rows_blank             = []
        exact_damage_supported_blank = []

        for _, meta_row in blank_meta.iterrows():
            megan_lib = meta_row["megan_library_name"]
            holi_lib  = meta_row["holi_library_name"]
            if counts.get(megan_lib, 0) <= 0:
                continue
            rows = []
            if tax_id_str:
                rows.extend(exact_by_id.get((holi_lib, tax_id_str), []))
            if not rows and tax_name:
                rows.extend(exact_by_name_rank.get((holi_lib, tax_name, tax_rank), []))
                if not rows:
                    rows.extend(exact_by_name_rank.get((holi_lib, tax_name, ""), []))
            best = choose_best_exact_row(rows, thresholds)
            if best:
                augmented = dict(best)
                augmented["megan_library_name"] = megan_lib
                exact_rows_blank.append(augmented)
                if row_has_exact_damage_support(best, thresholds):
                    exact_damage_supported_blank.append(augmented)

        exact_damage_support     = len(exact_damage_supported_real) > 0
        exact_damage_support_ge100 = any(
            float(x.get("N_reads") or 0) >= thresholds["high_confidence_n_reads_min"]
            for x in exact_damage_supported_real
        )
        blank_damage_support = len(exact_damage_supported_blank) > 0
        blank_caution = (
            max_blank_count >= thresholds["blank_absolute_min"]
            or (max_real_count > 0 and max_blank_count / max_real_count >= thresholds["blank_relative_min"])
        )
        blank_associated = blank_damage_support or (
            blank_caution and not strong_count_support and not exact_damage_support
        )

        # Best real-library exact row (highest scoring across all real libs).
        best_exact = None
        if exact_rows_real:
            best_exact = max(
                exact_rows_real,
                key=lambda x: (
                    int(row_has_exact_damage_support(x, thresholds)),
                    float(x.get("significance") or -np.inf),
                    float(x.get("N_reads") or -np.inf),
                    float(x.get("damage") or -np.inf),
                ),
            )

        # ── Exact Fillet matching for real libraries (only if fillet_df given) ──
        fillet_rows_real         = []
        fillet_authenticated_real = []
        best_fillet              = None
        if have_fillet:
            for _, meta_row in real_meta.iterrows():
                megan_lib  = meta_row["megan_library_name"]
                fillet_lib = meta_row.get("fillet_library_name")
                if pd.isna(fillet_lib) or not fillet_lib or counts.get(megan_lib, 0) <= 0:
                    continue
                rows = []
                if tax_id_str:
                    rows.extend(fillet_exact_by_id.get((fillet_lib, tax_id_str), []))
                if not rows and tax_name:
                    rows.extend(fillet_exact_by_name_rank.get((fillet_lib, tax_name, tax_rank), []))
                    if not rows:
                        rows.extend(fillet_exact_by_name_rank.get((fillet_lib, tax_name, ""), []))
                best = choose_best_fillet_row(rows, thresholds)
                if best:
                    augmented = dict(best)
                    augmented["megan_library_name"] = megan_lib
                    fillet_rows_real.append(augmented)
                    if row_is_fillet_authenticated(best, thresholds):
                        fillet_authenticated_real.append(augmented)
            if fillet_rows_real:
                best_fillet = max(
                    fillet_rows_real,
                    key=lambda x: (
                        int(row_is_fillet_authenticated(x, thresholds)),
                        float(x.get("composite_authenticity") or -np.inf),
                        float(x.get("direct_hard_reads") or -np.inf),
                    ),
                )

        fillet_authenticated = len(fillet_authenticated_real) > 0
        fillet_strong_count_support = any(
            row_has_fillet_strong_count_support(x, thresholds) for x in fillet_rows_real
        )
        fillet_eco_support = any(bool(x.get("eco_support")) for x in fillet_authenticated_real)
        fillet_pal_support = any(bool(x.get("pal_support")) for x in fillet_authenticated_real)
        fillet_fos_support = any(bool(x.get("fos_support")) for x in fillet_authenticated_real)

        # ── Canonical taxon info from the Holi taxonomy lookup ───────────────
        focal_info = None
        if tax_id_str and tax_id_str in tax_lookup:
            focal_info = tax_lookup[tax_id_str]
        elif tax_name and tax_name in tax_lookup:
            focal_info = tax_lookup[tax_name]
        else:
            unmatched_taxa.append(tax_name or tax_id_str or "unknown")
            focal_info = {"tax_name": tax_name, "tax_rank": tax_rank, "tax_path": ""}

        # ── Lineage support ──────────────────────────────────────────────────
        lineage_support_libraries = []
        lineage_support_examples  = []
        if config["lineage"]["enabled"] and focal_info:
            for _, meta_row in real_meta.iterrows():
                megan_lib = meta_row["megan_library_name"]
                holi_lib  = meta_row["holi_library_name"]
                if counts.get(megan_lib, 0) <= 0:
                    continue
                libs, examples = summarize_lineage_support_for_taxon(
                    focal_tax_name=focal_info.get("tax_name", tax_name),
                    focal_tax_rank=focal_info.get("tax_rank", tax_rank),
                    focal_tax_path=focal_info.get("tax_path", ""),
                    holi_sample=holi_lib,
                    supported_name_index_by_sample=supported_name_index_by_sample,
                    descendant_index_by_sample=descendant_index_by_sample,
                    package_cfg=config,
                )
                lineage_support_libraries.extend(libs)
                lineage_support_examples.extend(examples)

        lineage_support = len(lineage_support_libraries) > 0

        # ── Fillet vs. Holi/MEGAN taxonomic discordance ──────────────────────
        # True if, for some real library where MEGAN/Holi consider this taxon
        # present, Fillet's own authenticated call for that same library is a
        # DIFFERENT, taxonomically-unrelated taxon -- i.e. Fillet didn't just
        # fail to corroborate this taxon, it positively authenticated
        # something else entirely for the same sample. Deliberately
        # conservative, mirroring the existing lineage-support philosophy:
        # only flags discordance when both taxa's lineage paths are actually
        # resolvable (via Holi's own taxonomy lookup, since Fillet's export
        # carries no tax_path of its own) -- if we can't confidently tell
        # whether two taxa are related, this does NOT count as discordant.
        discordant = False
        if have_fillet and focal_info and focal_info.get("tax_path"):
            fillet_authenticated_taxa_here = {
                (x.get("tax_id_str"), x.get("tax_name")) for x in fillet_authenticated_real
            }
            for _, meta_row in real_meta.iterrows():
                megan_lib  = meta_row["megan_library_name"]
                fillet_lib = meta_row.get("fillet_library_name")
                if pd.isna(fillet_lib) or not fillet_lib or counts.get(megan_lib, 0) <= 0:
                    continue
                for other in fillet_authenticated_by_sample.get(fillet_lib, []):
                    other_key = (other.get("tax_id_str"), other.get("tax_name"))
                    if other_key in fillet_authenticated_taxa_here:
                        continue
                    other_info = (
                        tax_lookup.get(other.get("tax_id_str"))
                        or tax_lookup.get(other.get("tax_name"))
                    )
                    if other_info is None or not other_info.get("tax_path"):
                        continue
                    if is_meaningful_low_rank_lineage_support(
                        focal_tax_name=focal_info.get("tax_name", tax_name),
                        focal_tax_rank=focal_info.get("tax_rank", tax_rank),
                        focal_path=focal_info.get("tax_path", ""),
                        candidate_tax_name=other_info.get("tax_name", ""),
                        candidate_tax_rank=other_info.get("tax_rank", ""),
                        candidate_path=other_info.get("tax_path", ""),
                        cfg=config,
                    ):
                        continue
                    discordant = True
                    break
                if discordant:
                    break

        # ── QC label for the best exact row ──────────────────────────────────
        if best_exact:
            qc_label, align_ratio = compute_qc_label(
                n_reads=best_exact.get("N_reads"),
                n_alignments=best_exact.get("N_alignments"),
                map_valid=bool(best_exact.get("MAP_valid")),
                rho_ac=best_exact.get("rho_Ac"),
                thresholds=thresholds,
            )
        else:
            qc_label    = "not-applicable"
            align_ratio = np.nan

        # ── Final classification ─────────────────────────────────────────────
        # Routed through the source-agnostic classify_status_v2() rather than
        # calling classify_status() directly. When fillet_df is None (the
        # historical call signature), ev.fillet.present is False, and
        # classify_status_v2's MEGAN+Holi path calls the original
        # classify_status() directly and unmodified -- so this produces
        # byte-for-byte identical (status, basis) output to before Fillet
        # support existed (see evidence.py's TestBackwardCompatParity).
        evidence = TaxonEvidence(
            megan=SourceSignals(
                present=True,
                strong_count_support=strong_count_support,
                some_count_support=some_count_support,
                blank_associated=blank_associated,
            ),
            holi=SourceSignals(
                present=True,
                exact_damage_support=exact_damage_support,
                exact_damage_support_ge100=exact_damage_support_ge100,
                qc_label=qc_label,
                blank_associated=blank_associated,
            ),
            fillet=SourceSignals(
                present=have_fillet,
                exact_damage_support=fillet_authenticated if have_fillet else None,
                strong_count_support=fillet_strong_count_support if have_fillet else None,
                authenticated=fillet_authenticated if have_fillet else None,
                eco_support=fillet_eco_support if have_fillet else None,
                pal_support=fillet_pal_support if have_fillet else None,
                fos_support=fillet_fos_support if have_fillet else None,
            ),
            lineage_support=lineage_support,
            max_real_count=max_real_count,
            discordant=discordant,
        )
        status, basis, ensemble_support_score = classify_status_v2(evidence, config)

        # ── Common name ───────────────────────────────────────────────────────
        common_name, common_name_source = resolve_common_name(
            scientific_name=focal_info.get("tax_name", tax_name),
            tax_id=tax_id,
            builtin_map=config["taxonomy"]["builtin_common_name_map"],
            override_map=common_name_overrides,
            online=bool(config["taxonomy"].get("online_common_names")),
            cache=common_name_cache,
            language=config["taxonomy"].get("common_name_language", "eng"),
        )

        tax_path   = focal_info.get("tax_path", "")
        broad_grp  = _broad_group(tax_path)

        # ── Assemble output record ────────────────────────────────────────────
        record = {
            "tax_id":                            tax_id,
            "tax_id_str":                        tax_id_str,
            "scientific_name":                   focal_info.get("tax_name", tax_name),
            "tax_rank":                          focal_info.get("tax_rank", tax_rank),
            "tax_path":                          tax_path,
            "common_name":                       common_name,
            "common_name_source":                common_name_source,
            "broad_group":                       broad_grp,
            "aDNA_support_status":               status,
            "support_basis_summary":             basis,
            "megan_positive_libraries_n":        n_real_positive,
            "megan_positive_libraries":          "; ".join(real_positive_libs),
            "megan_max_count":                   max_real_count,
            "megan_blank_positive_libraries_n":  len(blank_positive_libs),
            "megan_blank_positive_libraries":    "; ".join(blank_positive_libs),
            "megan_blank_max_count":             max_blank_count,
            "blank_ratio":                       blank_ratio,
            "strong_count_support":              strong_count_support,
            "some_count_support":                some_count_support,
            "Holi_exact_damage_sig_libraries_n":    len(exact_damage_supported_real),
            "Holi_exact_damage_sig_libraries":      "; ".join(
                sorted({x["sample"] for x in exact_damage_supported_real})
            ),
            "Holi_exact_damage_sig_ge100_libraries_n": sum(
                1 for x in exact_damage_supported_real
                if float(x.get("N_reads") or 0) >= thresholds["high_confidence_n_reads_min"]
            ),
            "Low_rank_lineage_support_libraries_n":   len(sorted(set(lineage_support_libraries))),
            "Low_rank_lineage_support_libraries":     "; ".join(sorted(set(lineage_support_libraries))),
            "Low_rank_lineage_support_examples":      "; ".join(sorted(set(lineage_support_examples))[:10]),
            "Holi_blank_exact_damage_sig":            blank_damage_support,
            "Holi_blank_best_N_reads":    max(
                (float(x.get("N_reads") or 0) for x in exact_rows_blank), default=np.nan
            ),
            "Holi_blank_best_damage":     max(
                (float(x.get("damage") or np.nan) for x in exact_rows_blank), default=np.nan
            ),
            "Holi_blank_best_significance": max(
                (float(x.get("significance") or np.nan) for x in exact_rows_blank),
                default=np.nan,
            ),
            # ── Ensemble/Fillet fields (present regardless of fillet_df, so the
            # output schema is stable whether or not a given run supplies
            # Fillet evidence -- Fillet-specific fields are simply NaN/False
            # when fillet_df is None) ──────────────────────────────────────
            "ensemble_support_score":       round(ensemble_support_score, 4),
            "fillet_authenticated":         fillet_authenticated if have_fillet else pd.NA,
            "fillet_composite_authenticity": (
                float(best_fillet.get("composite_authenticity"))
                if have_fillet and best_fillet is not None and pd.notna(best_fillet.get("composite_authenticity"))
                else np.nan
            ),
            "fillet_authenticity_tier": (
                int(best_fillet.get("authenticity_tier"))
                if have_fillet and best_fillet is not None and pd.notna(best_fillet.get("authenticity_tier"))
                else pd.NA
            ),
            "fillet_direct_reads": (
                float(best_fillet.get("direct_hard_reads"))
                if have_fillet and best_fillet is not None and pd.notna(best_fillet.get("direct_hard_reads"))
                else np.nan
            ),
            "fillet_eco_support":  fillet_eco_support if have_fillet else pd.NA,
            "fillet_pal_support":  fillet_pal_support if have_fillet else pd.NA,
            "fillet_fos_support":  fillet_fos_support if have_fillet else pd.NA,
            "fillet_holi_megan_discordant": discordant if have_fillet else pd.NA,
        }

        # Per-library MEGAN counts as separate columns (prefixed with count__).
        for _, meta_row in metadata.iterrows():
            record["count__" + meta_row["merged_library_name"]] = counts.get(
                meta_row["megan_library_name"], 0.0
            )

        # Per-library aDNA support status (prefixed with aDNA_support_lib__).
        lineage_support_holi_libs = set(lineage_support_libraries)
        for _, meta_row in metadata.iterrows():
            megan_lib  = meta_row["megan_library_name"]
            holi_lib   = meta_row["holi_library_name"]
            merged_lib = meta_row["merged_library_name"]
            count_val  = counts.get(megan_lib, 0)
            is_pos = bool(is_pos_ctrl.get(meta_row.name, False))

            if meta_row["is_negative_control"]:
                lib_adna = "Blank-library" if count_val > 0 else "Not detected"
            elif bool(is_env_ctrl.get(meta_row.name, False)):
                lib_adna = "Environmental-control" if count_val > 0 else "Not detected"
            elif count_val <= 0:
                lib_adna = "Not detected"
            else:
                lib_dmg = [
                    x for x in exact_damage_supported_real
                    if x.get("megan_library_name") == megan_lib
                ]
                if lib_dmg:
                    has_ge100 = any(
                        float(x.get("N_reads") or 0) >= thresholds["high_confidence_n_reads_min"]
                        for x in lib_dmg
                    )
                    if has_ge100:
                        # Only award the top label when the best-supported row also has
                        # acceptable QC — a "strong caution" flag (e.g. extreme
                        # multi-mapping) makes the damage estimate unreliable regardless
                        # of read count.
                        best_lib = max(
                            lib_dmg,
                            key=lambda x: float(x.get("N_reads") or 0),
                        )
                        lib_qc, _ = compute_qc_label(
                            n_reads=best_lib.get("N_reads"),
                            n_alignments=best_lib.get("N_alignments"),
                            map_valid=bool(best_lib.get("MAP_valid")),
                            rho_ac=best_lib.get("rho_Ac"),
                            thresholds=thresholds,
                        )
                        lib_adna = (
                            "Damage-supported (>=100 reads)"
                            if lib_qc != "strong caution"
                            else "Damage-supported"
                        )
                    else:
                        lib_adna = "Damage-supported"
                elif holi_lib in lineage_support_holi_libs:
                    lib_adna = "Lineage-supported"
                else:
                    lib_adna = "Count-only"

            record[f"aDNA_support_lib__{merged_lib}"] = lib_adna

        # Per-library Holi damage values (NaN for controls or libraries with no Holi match).
        for _, meta_row in metadata.iterrows():
            megan_lib  = meta_row["megan_library_name"]
            merged_lib = meta_row["merged_library_name"]
            is_pos = bool(is_pos_ctrl.get(meta_row.name, False))
            if meta_row["is_negative_control"] or is_pos or bool(is_env_ctrl.get(meta_row.name, False)):
                record[f"Holi_damage_lib__{merged_lib}"]       = np.nan
                record[f"Holi_significance_lib__{merged_lib}"] = np.nan
            else:
                dmg_row = next(
                    (x for x in exact_rows_real if x.get("megan_library_name") == megan_lib),
                    None,
                )
                if dmg_row is not None:
                    record[f"Holi_damage_lib__{merged_lib}"] = (
                        float(dmg_row.get("damage") or np.nan)
                    )
                    record[f"Holi_significance_lib__{merged_lib}"] = (
                        float(dmg_row.get("significance") or np.nan)
                    )
                else:
                    record[f"Holi_damage_lib__{merged_lib}"]       = np.nan
                    record[f"Holi_significance_lib__{merged_lib}"] = np.nan

        # Best real-library Holi fields.
        if best_exact:
            record.update({
                "Holi_best_library":             best_exact.get("sample"),
                "Holi_best_damage_source_taxon": best_exact.get("tax_name"),
                "Holi_best_tax_rank":            best_exact.get("tax_rank"),
                "Holi_best_tax_path_short":      " > ".join(
                    [x for x in best_exact.get("tax_path", "").split(";")[-6:] if x]
                ),
                "Holi_best_N_reads":             float(best_exact.get("N_reads") or np.nan),
                "Holi_best_N_alignments":        float(best_exact.get("N_alignments") or np.nan),
                "Holi_best_damage":              float(best_exact.get("damage") or np.nan),
                "Holi_best_significance":        float(best_exact.get("significance") or np.nan),
                "Holi_best_abs_rho_Ac": (
                    abs(float(best_exact.get("rho_Ac")))
                    if best_exact.get("rho_Ac") is not None and pd.notna(best_exact.get("rho_Ac"))
                    else np.nan
                ),
                "Holi_best_MAP_valid":               bool(best_exact.get("MAP_valid")),
                "Holi_best_multimapping_fit_qc":     qc_label,
                "Holi_best_alignments_per_read":     align_ratio,
            })
        else:
            record.update({
                "Holi_best_library":             pd.NA,
                "Holi_best_damage_source_taxon": pd.NA,
                "Holi_best_tax_rank":            pd.NA,
                "Holi_best_tax_path_short":      pd.NA,
                "Holi_best_N_reads":             np.nan,
                "Holi_best_N_alignments":        np.nan,
                "Holi_best_damage":              np.nan,
                "Holi_best_significance":        np.nan,
                "Holi_best_abs_rho_Ac":          np.nan,
                "Holi_best_MAP_valid":           pd.NA,
                "Holi_best_multimapping_fit_qc": "not-applicable",
                "Holi_best_alignments_per_read": np.nan,
            })

        # Best blank-library Holi fields.
        if exact_rows_blank:
            best_blank = max(
                exact_rows_blank,
                key=lambda x: (
                    int(row_has_exact_damage_support(x, thresholds)),
                    float(x.get("significance") or -np.inf),
                    float(x.get("N_reads") or -np.inf),
                ),
            )
            n_reads_blank = best_blank.get("N_reads")
            record["Holi_blank_best_library"]         = best_blank.get("sample")
            record["Holi_blank_best_N_alignments"]    = float(best_blank.get("N_alignments") or np.nan)
            record["Holi_blank_best_alignments_per_read"] = (
                float(best_blank.get("N_alignments") or np.nan) / float(n_reads_blank)
                if n_reads_blank and float(n_reads_blank) > 0
                else np.nan
            )
        else:
            record["Holi_blank_best_library"]              = pd.NA
            record["Holi_blank_best_N_alignments"]         = np.nan
            record["Holi_blank_best_alignments_per_read"]  = np.nan

        records.append(record)

    if cache_path:
        save_common_name_cache(cache_path, common_name_cache)

    merged = pd.DataFrame.from_records(records).sort_values(
        by=["aDNA_support_status", "megan_max_count", "scientific_name"],
        key=lambda s: (
            s.map(STATUS_PRIORITY).fillna(99)
            if s.name == "aDNA_support_status"
            else s
        ),
        ascending=[True, False, True],
    )

    summary = {
        "n_taxa":                int(len(merged)),
        "status_counts":         merged["aDNA_support_status"].value_counts(dropna=False).to_dict(),
        "unmatched_taxa_examples": sorted(set(x for x in unmatched_taxa if x))[:20],
    }
    return merged, summary

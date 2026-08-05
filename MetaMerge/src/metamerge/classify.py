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


def _distinct_taxa(df: pd.DataFrame | None, id_col: str = "tax_id_str",
                    name_col: str = "tax_name", rank_col: str = "tax_rank") -> dict:
    """Return every distinct taxon appearing anywhere in ``df``, keyed by a
    canonical identity: ``tax_id_str`` when present, else ``(tax_name, tax_rank)``.

    Used to build the taxon universe for the no-MEGAN merge path, where there
    is no single count-matrix row-per-taxon to anchor on -- Holi and Fillet
    each carry one row per (taxon, sample), so the same taxon appears many
    times and must be deduplicated.
    """
    out: dict = {}
    if df is None or df.empty:
        return out
    cols = [c for c in (id_col, name_col, rank_col) if c in df.columns]
    for _, row in df[cols].drop_duplicates().iterrows():
        tid  = str(row.get(id_col)) if id_col in cols and pd.notna(row.get(id_col)) and row.get(id_col) != "" else ""
        name = normalize_name(row.get(name_col)) if name_col in cols else ""
        rank = normalize_rank(row.get(rank_col)) if rank_col in cols else ""
        if not tid and not name:
            continue
        key = tid if tid else (name, rank)
        if key not in out:
            out[key] = {"tax_id_str": tid, "tax_name": name, "tax_rank": rank}
    return out


def _real_count_proxy(holi_rows: list[dict], fillet_rows: list[dict]) -> float:
    """Best available "how much real read support does this library have"
    signal when there is no MEGAN count matrix to read it from directly.

    Prefers Holi's own ``N_reads`` (the product's existing "real read count"
    signal, already compared against ``high_confidence_n_reads_min`` etc.),
    falling back to Fillet's ``direct_hard_reads`` when Holi has no exact
    match for this taxon/library.
    """
    if holi_rows:
        n = float(max(holi_rows, key=lambda x: float(x.get("N_reads") or 0)).get("N_reads") or 0)
        if n > 0:
            return n
    if fillet_rows:
        return float(max(fillet_rows, key=lambda x: float(x.get("direct_hard_reads") or 0)).get("direct_hard_reads") or 0)
    return 0.0


def _per_library_adna_label(
    lib_holi_damage_rows: list[dict],
    lib_fillet_auth_rows: list[dict],
    thresholds: dict,
) -> tuple[str | None, int]:
    """Per-library aDNA_support_lib__ label, Fillet+Holi combined when both
    have a signal for this SPECIFIC library (Tyler, 2026-08-05: "the
    authenticity badge is the authenticity score expanded to be Fillet+Holi,
    given both can give their own damage/significance/likelihood metrics").

    Dual corroboration in the same library -- Holi's own damage signal AND
    Fillet's own authentication BOTH present for this library -- reaches the
    top label even when Holi's own read count alone doesn't independently
    clear the >=100-read bar, mirroring the taxon-level cascade's own
    "2 of 3 corroborating signals reaches the top tier" logic. Deliberately
    reuses the SAME label strings regardless of which source(s) contributed
    (rather than inventing new category strings per combination), since the
    per-source breakdown is already visible elsewhere in the output
    (fillet_authenticated, Holi_damage_lib__, etc.) -- this keeps the R
    heatmap's existing shape legend unchanged.

    Returns:
        Tuple of ``(label_or_None, per_library_methods_agreement_count)``.
        ``label`` is ``None`` when neither source has a signal for this
        library, so the caller falls back to its own Lineage-supported/
        Count-only logic.
    """
    holi_ok   = len(lib_holi_damage_rows) > 0
    fillet_ok = len(lib_fillet_auth_rows) > 0
    agreement = int(holi_ok) + int(fillet_ok)

    if not holi_ok and not fillet_ok:
        return None, agreement

    has_ge100 = holi_ok and any(
        float(x.get("N_reads") or 0) >= thresholds["high_confidence_n_reads_min"]
        for x in lib_holi_damage_rows
    )
    qc_ok = True
    if holi_ok:
        best_lib = max(lib_holi_damage_rows, key=lambda x: float(x.get("N_reads") or 0))
        lib_qc, _ = compute_qc_label(
            n_reads=best_lib.get("N_reads"),
            n_alignments=best_lib.get("N_alignments"),
            map_valid=bool(best_lib.get("MAP_valid")),
            rho_ac=best_lib.get("rho_Ac"),
            thresholds=thresholds,
        )
        qc_ok = lib_qc != "strong caution"

    if (has_ge100 or (holi_ok and fillet_ok)) and qc_ok:
        return "Damage-supported (>=100 reads)", agreement
    return "Damage-supported", agreement


def _build_merge_no_megan(
    metadata: pd.DataFrame,
    holi_df: pd.DataFrame | None,
    fillet_df: pd.DataFrame | None,
    config: dict,
    common_name_overrides: dict | None,
    cache_path,
) -> tuple[pd.DataFrame, dict]:
    """Merge path for when no MEGAN count matrix is supplied at all (Holi-only,
    Fillet-only, or Holi+Fillet -- Tyler's expected primary real-world use
    case going forward, per FILLET_METAMERGE_INTEGRATION_PLAN.md decision 1:
    "any subset of {MEGAN, Holi, Fillet} must work, not just 2 fixed + 1
    optional").

    There is no MEGAN count matrix to anchor the taxon universe or the
    per-library "how much real support" signal on, so both are rebuilt from
    whichever of Holi/Fillet are present: the taxon universe is the union of
    every distinct taxon in Holi and/or Fillet (see ``_distinct_taxa``), and
    the per-library count proxy is Holi's N_reads, falling back to Fillet's
    direct_hard_reads (see ``_real_count_proxy``).

    Deliberately kept as a separate function from the MEGAN-anchored merge
    loop above, rather than unified into one taxon-union loop, so that the
    existing, tested MEGAN-anchored path (build_merge's main body) is not put
    at risk of a subtle behavior change -- e.g. this function's taxon universe
    legitimately includes taxa MEGAN never called at all, which the
    MEGAN-anchored path must NOT start doing (that would break the "output is
    byte-for-byte unchanged when fillet_df is None" backward-compatibility
    guarantee tested in TestBackwardCompatParity/test_smoke.py).

    Produces the SAME output record schema as the MEGAN-anchored path (see
    build_merge's docstring) so workbook.py/report.py need no changes. Fields
    literally named ``megan_*``/``count__`` are populated from the fallback
    proxy described above rather than true MEGAN counts -- this is
    deliberately documented here rather than solved by renaming the schema, to
    avoid a second breaking change to every downstream consumer. Heatmap/plot
    rendering (``write_plot_inputs``, the R script) still assumes MEGAN-style
    counts are the primary visualization signal; that remains a known,
    separate follow-up, not attempted here.
    """
    thresholds = config["thresholds"]
    have_holi   = holi_df is not None
    have_fillet = fillet_df is not None

    if have_holi:
        exact_by_id, exact_by_name_rank = make_holi_exact_index(holi_df)
        tax_lookup = build_holi_taxonomy_lookup(holi_df)
        _, supported_name_index_by_sample, descendant_index_by_sample = build_lineage_indexes(
            holi_df, thresholds
        )
    else:
        exact_by_id, exact_by_name_rank, tax_lookup = {}, {}, {}
        supported_name_index_by_sample, descendant_index_by_sample = {}, {}

    if have_fillet:
        fillet_exact_by_id, fillet_exact_by_name_rank = make_fillet_exact_index(fillet_df)
        fillet_authenticated_by_sample: dict = defaultdict(list)
        for _row in fillet_df.to_dict(orient="records"):
            if row_is_fillet_authenticated(_row, thresholds):
                fillet_authenticated_by_sample[_row["sample"]].append(_row)
        fillet_authenticated_by_sample = dict(fillet_authenticated_by_sample)
    else:
        fillet_exact_by_id, fillet_exact_by_name_rank, fillet_authenticated_by_sample = {}, {}, {}

    common_name_overrides = common_name_overrides or {}
    common_name_cache = load_common_name_cache(cache_path) if cache_path else {}

    is_pos_ctrl = metadata.get("is_positive_control",      pd.Series(False, index=metadata.index))
    is_env_ctrl = metadata.get("is_environmental_control", pd.Series(False, index=metadata.index))
    real_meta  = metadata.loc[~metadata["is_negative_control"] & ~is_env_ctrl].copy()

    taxon_union: dict = {}
    if have_holi:
        for key, info in _distinct_taxa(holi_df).items():
            taxon_union.setdefault(key, info)
    if have_fillet:
        for key, info in _distinct_taxa(fillet_df).items():
            taxon_union.setdefault(key, info)

    records = []
    unmatched_taxa = []

    for taxon_info in tqdm(taxon_union.values(), desc="Classifying taxa (no MEGAN)", unit="taxon"):
        tax_id_str = taxon_info["tax_id_str"]
        tax_name   = taxon_info["tax_name"]
        tax_rank   = taxon_info["tax_rank"]
        tax_id     = int(tax_id_str) if tax_id_str.isdigit() else pd.NA

        # ── One pass over every library: exact Holi/Fillet matches + the
        # real-count proxy, split into real/blank/env by control status ──────
        exact_rows_real, exact_damage_supported_real = [], []
        exact_rows_blank, exact_damage_supported_blank = [], []
        fillet_rows_real, fillet_authenticated_real = [], []
        all_lib_counts: dict = {}
        real_positive_libs, blank_positive_libs = [], []
        real_counts, blank_counts = {}, {}

        for _, meta_row in metadata.iterrows():
            merged_lib = meta_row["merged_library_name"]
            is_blank   = bool(meta_row["is_negative_control"])
            is_env     = bool(is_env_ctrl.get(meta_row.name, False))
            holi_lib   = meta_row.get("holi_library_name")
            fillet_lib = meta_row.get("fillet_library_name")

            lib_holi_rows = []
            if have_holi and pd.notna(holi_lib) and holi_lib:
                if tax_id_str:
                    lib_holi_rows.extend(exact_by_id.get((holi_lib, tax_id_str), []))
                if not lib_holi_rows and tax_name:
                    lib_holi_rows.extend(exact_by_name_rank.get((holi_lib, tax_name, tax_rank), []))
                    if not lib_holi_rows:
                        lib_holi_rows.extend(exact_by_name_rank.get((holi_lib, tax_name, ""), []))
            best_holi = choose_best_exact_row(lib_holi_rows, thresholds)

            lib_fillet_rows = []
            if have_fillet and pd.notna(fillet_lib) and fillet_lib:
                if tax_id_str:
                    lib_fillet_rows.extend(fillet_exact_by_id.get((fillet_lib, tax_id_str), []))
                if not lib_fillet_rows and tax_name:
                    lib_fillet_rows.extend(fillet_exact_by_name_rank.get((fillet_lib, tax_name, tax_rank), []))
                    if not lib_fillet_rows:
                        lib_fillet_rows.extend(fillet_exact_by_name_rank.get((fillet_lib, tax_name, ""), []))
            best_fillet_row = choose_best_fillet_row(lib_fillet_rows, thresholds)

            proxy = _real_count_proxy([best_holi] if best_holi else [], [best_fillet_row] if best_fillet_row else [])
            all_lib_counts[merged_lib] = proxy

            if not is_blank and not is_env:
                real_counts[merged_lib] = proxy
                if proxy > 0:
                    real_positive_libs.append(merged_lib)
                if best_holi:
                    augmented = dict(best_holi)
                    augmented["merged_library_name"] = merged_lib
                    exact_rows_real.append(augmented)
                    if row_has_exact_damage_support(best_holi, thresholds):
                        exact_damage_supported_real.append(augmented)
                if best_fillet_row:
                    augmented = dict(best_fillet_row)
                    augmented["merged_library_name"] = merged_lib
                    fillet_rows_real.append(augmented)
                    if row_is_fillet_authenticated(best_fillet_row, thresholds):
                        fillet_authenticated_real.append(augmented)
            elif is_blank:
                blank_counts[merged_lib] = proxy
                if proxy > 0:
                    blank_positive_libs.append(merged_lib)
                if best_holi:
                    augmented = dict(best_holi)
                    augmented["merged_library_name"] = merged_lib
                    exact_rows_blank.append(augmented)
                    if row_has_exact_damage_support(best_holi, thresholds):
                        exact_damage_supported_blank.append(augmented)

        n_real_positive = len(real_positive_libs)
        max_real_count  = max(real_counts.values())  if real_counts  else 0.0
        max_blank_count = max(blank_counts.values()) if blank_counts else 0.0
        blank_ratio     = (max_blank_count / max_real_count) if max_real_count > 0 else np.nan

        strong_count_support = (
            n_real_positive >= thresholds["strong_count_min_libraries"]
            and max_real_count >= thresholds["strong_count_min_reads"]
        )
        some_count_support = n_real_positive >= 1

        exact_damage_support = len(exact_damage_supported_real) > 0
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

        fillet_authenticated = len(fillet_authenticated_real) > 0
        fillet_strong_count_support = any(
            row_has_fillet_strong_count_support(x, thresholds) for x in fillet_rows_real
        )
        fillet_eco_support = any(bool(x.get("eco_support")) for x in fillet_authenticated_real)
        fillet_pal_support = any(bool(x.get("pal_support")) for x in fillet_authenticated_real)
        fillet_fos_support = any(bool(x.get("fos_support")) for x in fillet_authenticated_real)
        best_fillet = None
        if fillet_rows_real:
            best_fillet = max(
                fillet_rows_real,
                key=lambda x: (
                    int(row_is_fillet_authenticated(x, thresholds)),
                    float(x.get("composite_authenticity") or -np.inf),
                    float(x.get("direct_hard_reads") or -np.inf),
                ),
            )

        # ── Canonical taxon info (only available via Holi's tax_path) ────────
        focal_info = None
        if tax_id_str and tax_id_str in tax_lookup:
            focal_info = tax_lookup[tax_id_str]
        elif tax_name and tax_name in tax_lookup:
            focal_info = tax_lookup[tax_name]
        else:
            unmatched_taxa.append(tax_name or tax_id_str or "unknown")
            focal_info = {"tax_name": tax_name, "tax_rank": tax_rank, "tax_path": ""}

        # ── Lineage support (only possible when Holi is present) ─────────────
        lineage_support_libraries, lineage_support_examples = [], []
        if have_holi and config["lineage"]["enabled"] and focal_info:
            for _, meta_row in real_meta.iterrows():
                holi_lib = meta_row.get("holi_library_name")
                if pd.isna(holi_lib) or not holi_lib:
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

        # ── Fillet vs. Holi taxonomic discordance (needs both present) ───────
        discordant = False
        if have_fillet and have_holi and focal_info and focal_info.get("tax_path"):
            fillet_authenticated_taxa_here = {
                (x.get("tax_id_str"), x.get("tax_name")) for x in fillet_authenticated_real
            }
            for _, meta_row in real_meta.iterrows():
                fillet_lib = meta_row.get("fillet_library_name")
                if pd.isna(fillet_lib) or not fillet_lib:
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

        # ── QC label for the best exact Holi row ──────────────────────────────
        if best_exact:
            qc_label, align_ratio = compute_qc_label(
                n_reads=best_exact.get("N_reads"),
                n_alignments=best_exact.get("N_alignments"),
                map_valid=bool(best_exact.get("MAP_valid")),
                rho_ac=best_exact.get("rho_Ac"),
                thresholds=thresholds,
            )
        else:
            qc_label, align_ratio = "not-applicable", np.nan

        # ── Final classification ──────────────────────────────────────────────
        evidence = TaxonEvidence(
            megan=SourceSignals(present=False),
            holi=SourceSignals(
                present=have_holi,
                exact_damage_support=exact_damage_support if have_holi else None,
                exact_damage_support_ge100=exact_damage_support_ge100 if have_holi else None,
                qc_label=qc_label if have_holi else None,
                # blank_associated is a merge-level determination (from
                # whichever real-count proxy is available), not something only
                # Holi personally detects -- always set here regardless of
                # have_holi, matching the existing MEGAN-anchored path's own
                # convention of setting the same computed value on both its
                # megan and holi slots unconditionally.
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
        status, basis, ensemble_support_score, methods_agreement_count, methods_agreement_fraction = (
            classify_status_v2(evidence, config)
        )

        # ── Common name ────────────────────────────────────────────────────────
        common_name, common_name_source = resolve_common_name(
            scientific_name=focal_info.get("tax_name", tax_name),
            tax_id=int(tax_id) if pd.notna(tax_id) else None,
            builtin_map=config["taxonomy"]["builtin_common_name_map"],
            override_map=common_name_overrides,
            online=bool(config["taxonomy"].get("online_common_names")),
            cache=common_name_cache,
            language=config["taxonomy"].get("common_name_language", "eng"),
        )

        tax_path  = focal_info.get("tax_path", "")
        broad_grp = _broad_group(tax_path)

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
            "ensemble_support_score":       round(ensemble_support_score, 4),
            "methods_agreement_count":      methods_agreement_count,
            "methods_agreement_fraction":   round(methods_agreement_fraction, 4),
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

        # Per-library counts (count__), keyed by merged_library_name -- the
        # only universally-present identifier when MEGAN is absent.
        for _, meta_row in metadata.iterrows():
            merged_lib = meta_row["merged_library_name"]
            record["count__" + merged_lib] = all_lib_counts.get(merged_lib, 0.0)

        # Per-library aDNA support status (aDNA_support_lib__), Fillet+Holi
        # combined when both have a signal for this specific library -- plus
        # a per-library methods_agreement_lib__ count driving the R
        # heatmap's separate ensemble-tier badge (see _per_library_adna_label).
        lineage_support_holi_libs = set(lineage_support_libraries)
        for _, meta_row in metadata.iterrows():
            merged_lib = meta_row["merged_library_name"]
            holi_lib   = meta_row.get("holi_library_name")
            count_val  = all_lib_counts.get(merged_lib, 0)
            lib_agreement = 0

            if meta_row["is_negative_control"]:
                lib_adna = "Blank-library" if count_val > 0 else "Not detected"
            elif bool(is_env_ctrl.get(meta_row.name, False)):
                lib_adna = "Environmental-control" if count_val > 0 else "Not detected"
            elif count_val <= 0:
                lib_adna = "Not detected"
            else:
                lib_dmg = [
                    x for x in exact_damage_supported_real
                    if x.get("merged_library_name") == merged_lib
                ]
                lib_fillet_auth = [
                    x for x in fillet_authenticated_real
                    if x.get("merged_library_name") == merged_lib
                ]
                lib_adna, lib_agreement = _per_library_adna_label(lib_dmg, lib_fillet_auth, thresholds)
                if lib_adna is None:
                    lib_adna = "Lineage-supported" if holi_lib in lineage_support_holi_libs else "Count-only"

            record[f"aDNA_support_lib__{merged_lib}"]        = lib_adna
            record[f"methods_agreement_lib__{merged_lib}"]   = lib_agreement

        # Per-library Holi damage values.
        for _, meta_row in metadata.iterrows():
            merged_lib = meta_row["merged_library_name"]
            is_pos = bool(is_pos_ctrl.get(meta_row.name, False))
            if meta_row["is_negative_control"] or is_pos or bool(is_env_ctrl.get(meta_row.name, False)):
                record[f"Holi_damage_lib__{merged_lib}"]       = np.nan
                record[f"Holi_significance_lib__{merged_lib}"] = np.nan
            else:
                dmg_row = next(
                    (x for x in exact_rows_real if x.get("merged_library_name") == merged_lib),
                    None,
                )
                if dmg_row is not None:
                    record[f"Holi_damage_lib__{merged_lib}"]       = float(dmg_row.get("damage") or np.nan)
                    record[f"Holi_significance_lib__{merged_lib}"] = float(dmg_row.get("significance") or np.nan)
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

    merged = pd.DataFrame.from_records(records)
    if not merged.empty:
        merged = merged.sort_values(
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
        "status_counts":         merged["aDNA_support_status"].value_counts(dropna=False).to_dict() if not merged.empty else {},
        "unmatched_taxa_examples": sorted(set(x for x in unmatched_taxa if x))[:20],
    }
    return merged, summary


def _build_merge_megan_fillet_no_holi(
    metadata: pd.DataFrame,
    megan_df: pd.DataFrame,
    fillet_df: pd.DataFrame,
    config: dict,
    common_name_overrides: dict | None,
    cache_path,
) -> tuple[pd.DataFrame, dict]:
    """MEGAN+Fillet merge path with no Holi at all.

    Per Tyler's design (2026-08-05): by default Fillet is the primary
    counts+support source (its own composite_authenticity/authenticity_tier +
    eco/pal/fos lines already combine multiple evidence dimensions on their
    own), MEGAN's own counts are additional corroborating support -- see
    ``evidence._classify_megan_fillet`` and
    ``config["source_roles"]["megan_fillet_primary"]`` to swap the roles.

    MEGAN still anchors the taxon universe here (unlike
    ``_build_merge_no_megan``, which needs a real taxon-union since there is
    no count matrix to anchor on at all) -- this keeps the well-tested
    "iterate megan_df's own taxa" structure from the MEGAN-anchored body
    above, just without any Holi lookups, plus real Fillet exact-matching.

    No lineage support or discordance detection is possible in this path:
    both require a taxonomy/lineage source with real ``tax_path`` data, which
    only Holi supplies today (Fillet's own export carries no ``tax_path``).
    Conservatively left as False/empty rather than guessed at, matching the
    existing "if we can't confidently tell, don't flag" philosophy.
    """
    thresholds = config["thresholds"]
    common_name_overrides = common_name_overrides or {}
    common_name_cache = load_common_name_cache(cache_path) if cache_path else {}

    if "fillet_library_name" not in metadata.columns:
        raise ValueError(
            "build_merge() was given fillet_df but metadata has no "
            "fillet_library_name column -- add one to the linker file."
        )
    fillet_exact_by_id, fillet_exact_by_name_rank = make_fillet_exact_index(fillet_df)
    fillet_tax_lookup = build_fillet_taxonomy_lookup(fillet_df)

    is_pos_ctrl = metadata.get("is_positive_control",      pd.Series(False, index=metadata.index))
    is_env_ctrl = metadata.get("is_environmental_control", pd.Series(False, index=metadata.index))
    real_meta  = metadata.loc[~metadata["is_negative_control"] & ~is_env_ctrl].copy()
    blank_meta = metadata.loc[ metadata["is_negative_control"]].copy()

    megan_library_cols = metadata["megan_library_name"].tolist()
    real_cols          = real_meta["megan_library_name"].tolist()
    blank_cols         = blank_meta["megan_library_name"].tolist()

    records        = []
    unmatched_taxa = []

    for row in tqdm(megan_df.to_dict(orient="records"), desc="Classifying taxa (no Holi)", unit="taxon"):
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

        # Blank detection is purely MEGAN-count-driven here -- no Holi damage
        # signal exists to corroborate blanks with, unlike the MEGAN-anchored
        # body's blank_damage_support check.
        blank_caution = (
            max_blank_count >= thresholds["blank_absolute_min"]
            or (max_real_count > 0 and max_blank_count / max_real_count >= thresholds["blank_relative_min"])
        )
        blank_associated = blank_caution and not strong_count_support

        # ── Exact Fillet matching for real libraries ──────────────────────────
        fillet_rows_real, fillet_authenticated_real = [], []
        best_fillet = None
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

        # ── Canonical taxon info: no Holi taxonomy lookup available; fall
        # back to Fillet's own name/rank lookup (still no tax_path, matching
        # Fillet's export schema) ──────────────────────────────────────────
        focal_info = None
        if tax_id_str and tax_id_str in fillet_tax_lookup:
            focal_info = fillet_tax_lookup[tax_id_str]
        elif tax_name and tax_name in fillet_tax_lookup:
            focal_info = fillet_tax_lookup[tax_name]
        else:
            unmatched_taxa.append(tax_name or tax_id_str or "unknown")
            focal_info = {"tax_name": tax_name, "tax_rank": tax_rank, "tax_path": ""}

        # No lineage support or discordance detection possible without a
        # tax_path source (see docstring).
        lineage_support = False
        discordant       = False

        evidence = TaxonEvidence(
            megan=SourceSignals(
                present=True,
                strong_count_support=strong_count_support,
                some_count_support=some_count_support,
                blank_associated=blank_associated,
            ),
            holi=SourceSignals(present=False, blank_associated=blank_associated),
            fillet=SourceSignals(
                present=True,
                exact_damage_support=fillet_authenticated,
                strong_count_support=fillet_strong_count_support,
                authenticated=fillet_authenticated,
                eco_support=fillet_eco_support,
                pal_support=fillet_pal_support,
                fos_support=fillet_fos_support,
            ),
            lineage_support=lineage_support,
            max_real_count=max_real_count,
            discordant=discordant,
        )
        status, basis, ensemble_support_score, methods_agreement_count, methods_agreement_fraction = (
            classify_status_v2(evidence, config)
        )

        common_name, common_name_source = resolve_common_name(
            scientific_name=focal_info.get("tax_name", tax_name),
            tax_id=tax_id,
            builtin_map=config["taxonomy"]["builtin_common_name_map"],
            override_map=common_name_overrides,
            online=bool(config["taxonomy"].get("online_common_names")),
            cache=common_name_cache,
            language=config["taxonomy"].get("common_name_language", "eng"),
        )
        tax_path  = focal_info.get("tax_path", "")
        broad_grp = _broad_group(tax_path)

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
            "Holi_exact_damage_sig_libraries_n":       0,
            "Holi_exact_damage_sig_libraries":          "",
            "Holi_exact_damage_sig_ge100_libraries_n": 0,
            "Low_rank_lineage_support_libraries_n":    0,
            "Low_rank_lineage_support_libraries":       "",
            "Low_rank_lineage_support_examples":        "",
            "Holi_blank_exact_damage_sig":              False,
            "Holi_blank_best_N_reads":                  np.nan,
            "Holi_blank_best_damage":                   np.nan,
            "Holi_blank_best_significance":             np.nan,
            "ensemble_support_score":       round(ensemble_support_score, 4),
            "methods_agreement_count":      methods_agreement_count,
            "methods_agreement_fraction":   round(methods_agreement_fraction, 4),
            "fillet_authenticated":         fillet_authenticated,
            "fillet_composite_authenticity": (
                float(best_fillet.get("composite_authenticity"))
                if best_fillet is not None and pd.notna(best_fillet.get("composite_authenticity"))
                else np.nan
            ),
            "fillet_authenticity_tier": (
                int(best_fillet.get("authenticity_tier"))
                if best_fillet is not None and pd.notna(best_fillet.get("authenticity_tier"))
                else pd.NA
            ),
            "fillet_direct_reads": (
                float(best_fillet.get("direct_hard_reads"))
                if best_fillet is not None and pd.notna(best_fillet.get("direct_hard_reads"))
                else np.nan
            ),
            "fillet_eco_support":  fillet_eco_support,
            "fillet_pal_support":  fillet_pal_support,
            "fillet_fos_support":  fillet_fos_support,
            "fillet_holi_megan_discordant": discordant,
        }

        for _, meta_row in metadata.iterrows():
            record["count__" + meta_row["merged_library_name"]] = counts.get(
                meta_row["megan_library_name"], 0.0
            )

        # Per-library aDNA support status -- no Holi damage tiers possible
        # here (lib_dmg is always empty), so _per_library_adna_label
        # collapses to Fillet-authenticated/Count-only/blank/env-control/
        # not-detected. Still reuses the shared helper for consistency with
        # the other two build paths' methods_agreement_lib__ column.
        for _, meta_row in metadata.iterrows():
            megan_lib  = meta_row["megan_library_name"]
            merged_lib = meta_row["merged_library_name"]
            count_val  = counts.get(megan_lib, 0)
            lib_agreement = 0
            if meta_row["is_negative_control"]:
                lib_adna = "Blank-library" if count_val > 0 else "Not detected"
            elif bool(is_env_ctrl.get(meta_row.name, False)):
                lib_adna = "Environmental-control" if count_val > 0 else "Not detected"
            elif count_val <= 0:
                lib_adna = "Not detected"
            else:
                lib_fillet_auth = [x for x in fillet_authenticated_real if x.get("megan_library_name") == megan_lib]
                lib_adna, lib_agreement = _per_library_adna_label([], lib_fillet_auth, thresholds)
                if lib_adna is None:
                    lib_adna = "Count-only"
            record[f"aDNA_support_lib__{merged_lib}"]      = lib_adna
            record[f"methods_agreement_lib__{merged_lib}"] = lib_agreement

        for _, meta_row in metadata.iterrows():
            merged_lib = meta_row["merged_library_name"]
            record[f"Holi_damage_lib__{merged_lib}"]       = np.nan
            record[f"Holi_significance_lib__{merged_lib}"] = np.nan

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
        record["Holi_blank_best_library"]             = pd.NA
        record["Holi_blank_best_N_alignments"]        = np.nan
        record["Holi_blank_best_alignments_per_read"] = np.nan

        records.append(record)

    if cache_path:
        save_common_name_cache(cache_path, common_name_cache)

    merged = pd.DataFrame.from_records(records)
    if not merged.empty:
        merged = merged.sort_values(
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
        "status_counts":         merged["aDNA_support_status"].value_counts(dropna=False).to_dict() if not merged.empty else {},
        "unmatched_taxa_examples": sorted(set(x for x in unmatched_taxa if x))[:20],
    }
    return merged, summary


def build_merge(
    metadata: pd.DataFrame,
    megan_df: pd.DataFrame | None,
    holi_df: pd.DataFrame | None,
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

    Any non-empty subset of {MEGAN, Holi, Fillet} is supported: pass ``None``
    for whichever source(s) are unavailable. When ``megan_df`` is ``None``,
    this delegates to ``_build_merge_no_megan``, which builds the taxon
    universe and per-library "real count" signal from whichever of
    Holi/Fillet are present instead of a MEGAN count matrix -- this is
    Tyler's expected primary real-world case going forward (Holi+Fillet, no
    MEGAN). When ``holi_df`` is ``None`` but ``megan_df``/``fillet_df`` are
    both supplied, this delegates to ``_build_merge_megan_fillet_no_holi``,
    which keeps MEGAN's own taxon-anchored structure but adds real Fillet
    exact-matching/authentication with no Holi lookups at all (no lineage
    support or discordance detection is possible in that path, since both
    require a real ``tax_path`` source, which only Holi supplies today). At
    least one of the three must be supplied.

    Backward compatibility: when ``megan_df`` is supplied and ``fillet_df`` is
    ``None`` (the historical call signature), every taxon's DNA-support status
    is produced by calling the original, untouched ``classify_status()``
    directly (via ``evidence.classify_status_v2``'s MEGAN+Holi path) -- this
    function's output is unchanged from before Fillet support existed.

    Args:
        metadata: Validated library-linker DataFrame (from ``load_metadata``).
        megan_df: Loaded MEGAN count matrix (from ``load_megan_counts``), or
            ``None`` when no MEGAN count matrix is available for this run.
        holi_df: Loaded metaDMG/Holi table (from ``load_holi``), or ``None``.
        config: Full MetaMerge config dict.
        common_name_overrides: Optional dict from ``load_common_name_overrides``.
        cache_path: Optional Path for persisting online common-name lookups.
        fillet_df: Optional loaded Fillet MetaMerge evidence table (from
            ``io.load_fillet``). When supplied, ``metadata`` must have a
            ``fillet_library_name`` column (raises ``ValueError`` otherwise).

    Returns:
        Tuple of ``(merged_df, summary_dict)`` where:
          - ``merged_df`` has one row per taxon (from the MEGAN count matrix
            when supplied, otherwise the union of Holi/Fillet taxa), sorted by
            DNA support status, real-count signal, and scientific name.
          - ``summary_dict`` contains ``n_taxa``, ``status_counts``, and
            ``unmatched_taxa_examples``.
    """
    if megan_df is None and holi_df is None and fillet_df is None:
        raise ValueError(
            "build_merge() needs at least one of megan_df, holi_df, fillet_df -- "
            "all three were None."
        )

    if megan_df is None:
        return _build_merge_no_megan(
            metadata=metadata,
            holi_df=holi_df,
            fillet_df=fillet_df,
            config=config,
            common_name_overrides=common_name_overrides,
            cache_path=cache_path,
        )

    if holi_df is None:
        if fillet_df is None:
            raise ValueError(
                "build_merge() was given megan_df but neither holi_df nor "
                "fillet_df -- the MEGAN-anchored merge path needs at least "
                "one of Holi or Fillet alongside MEGAN. Pass megan_df=None "
                "instead for a Holi/Fillet-only merge."
            )
        return _build_merge_megan_fillet_no_holi(
            metadata=metadata,
            megan_df=megan_df,
            fillet_df=fillet_df,
            config=config,
            common_name_overrides=common_name_overrides,
            cache_path=cache_path,
        )

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
        status, basis, ensemble_support_score, methods_agreement_count, methods_agreement_fraction = (
            classify_status_v2(evidence, config)
        )

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
            "methods_agreement_count":      methods_agreement_count,
            "methods_agreement_fraction":   round(methods_agreement_fraction, 4),
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

        # Per-library aDNA support status (prefixed with aDNA_support_lib__),
        # Fillet+Holi combined when both have a signal for this specific
        # library -- plus a per-library methods_agreement_lib__ count driving
        # the R heatmap's separate ensemble-tier badge.
        lineage_support_holi_libs = set(lineage_support_libraries)
        for _, meta_row in metadata.iterrows():
            megan_lib  = meta_row["megan_library_name"]
            holi_lib   = meta_row["holi_library_name"]
            merged_lib = meta_row["merged_library_name"]
            count_val  = counts.get(megan_lib, 0)
            is_pos = bool(is_pos_ctrl.get(meta_row.name, False))
            lib_agreement = 0

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
                lib_fillet_auth = [
                    x for x in fillet_authenticated_real
                    if x.get("megan_library_name") == megan_lib
                ]
                lib_adna, lib_agreement = _per_library_adna_label(lib_dmg, lib_fillet_auth, thresholds)
                if lib_adna is None:
                    lib_adna = "Lineage-supported" if holi_lib in lineage_support_holi_libs else "Count-only"

            record[f"aDNA_support_lib__{merged_lib}"]      = lib_adna
            record[f"methods_agreement_lib__{merged_lib}"] = lib_agreement

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

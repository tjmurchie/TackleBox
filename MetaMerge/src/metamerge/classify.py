"""Core merge and DNA-support classification logic for MetaMerge.

This module contains the main ``build_merge`` function and the ``classify_status``
helper.  Together they implement the per-taxon classification pipeline:

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
from .holi import (
    build_holi_taxonomy_lookup,
    build_lineage_indexes,
    choose_best_exact_row,
    compute_qc_label,
    make_holi_exact_index,
    row_has_exact_damage_support,
    summarize_lineage_support_for_taxon,
)
from .utils import normalize_name, normalize_rank, select_best_status


def classify_status(
    exact_damage_support: bool,
    exact_damage_support_ge100: bool,
    strong_count_support: bool,
    some_count_support: bool,
    lineage_support: bool,
    blank_associated: bool,
    qc_label: str,
    max_real_count: float,
    thresholds: dict,
) -> tuple[str, str]:
    """Assign a DNA-support category and short basis summary.

    Applies the classification rules in priority order.  Blank-associated always
    wins; Very high confidence requires the strictest combination of signals.

    Args:
        exact_damage_support: True if any real library has exact Holi damage+sig.
        exact_damage_support_ge100: True if any supporting row has N_reads ≥ threshold.
        strong_count_support: True if count evidence meets the strong-count thresholds.
        some_count_support: True if at least one real library is positive.
        lineage_support: True if any real library has lineage-consistent support.
        blank_associated: True if blank evidence dominates.
        qc_label: QC label from compute_qc_label ("clean", "caution",
            "strong caution").
        max_real_count: Maximum MEGAN count across all real libraries.
        thresholds: The "thresholds" sub-dict from the MetaMerge config.

    Returns:
        Tuple of ``(status_string, basis_summary_string)``.
    """
    if blank_associated:
        return "Blank-associated", "blank-associated"

    if exact_damage_support_ge100 and strong_count_support and qc_label != "strong caution":
        return (
            "Very high confidence",
            "exact damage-supported; strong counts; >=100 Holi reads",
        )

    if exact_damage_support and strong_count_support:
        if qc_label == "clean":
            return "High confidence", "exact damage-supported; strong counts"
        return "High confidence", f"exact damage-supported; strong counts; {qc_label}"

    if (
        exact_damage_support
        or strong_count_support
        or (lineage_support and some_count_support)
    ):
        basis = []
        if exact_damage_support:
            basis.append("damage-supported (exact)")
        if lineage_support:
            basis.append("lineage-supported")
        if strong_count_support:
            basis.append("strong counts")
        elif some_count_support:
            basis.append("some counts")
        if qc_label != "clean" and exact_damage_support:
            basis.append(qc_label)
        return "Supported", "; ".join(basis) if basis else "supported"

    if max_real_count >= thresholds["tentative_min_reads"]:
        return "Tentative", "weak/incomplete DNA support"

    if max_real_count >= thresholds["weak_support_min_reads"]:
        return "Weak support", "very weak DNA support"

    return "Weak support", "no non-control support"


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
) -> tuple[pd.DataFrame, dict]:
    """Merge MEGAN + Holi/metaDMG + metadata into one row per taxon.

    This is the project-neutral classification engine.  It classifies taxa
    based on DNA evidence (damage signals, read counts, lineage consistency).

    Args:
        metadata: Validated library-linker DataFrame (from ``load_metadata``).
        megan_df: Loaded MEGAN count matrix (from ``load_megan_counts``).
        holi_df: Loaded metaDMG/Holi table (from ``load_holi``).
        config: Full MetaMerge config dict.
        common_name_overrides: Optional dict from ``load_common_name_overrides``.
        cache_path: Optional Path for persisting online common-name lookups.

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
        status, basis = classify_status(
            exact_damage_support=exact_damage_support,
            exact_damage_support_ge100=exact_damage_support_ge100,
            strong_count_support=strong_count_support,
            some_count_support=some_count_support,
            lineage_support=lineage_support,
            blank_associated=blank_associated,
            qc_label=qc_label,
            max_real_count=max_real_count,
            thresholds=thresholds,
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

    STATUS_SORT_MAP = {
        "Very high confidence": 0,
        "High confidence": 1,
        "Supported": 2,
        "Tentative": 3,
        "Weak support": 4,
        "Blank-associated": 5,
    }

    merged = pd.DataFrame.from_records(records).sort_values(
        by=["aDNA_support_status", "megan_max_count", "scientific_name"],
        key=lambda s: (
            s.map(STATUS_SORT_MAP).fillna(99)
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

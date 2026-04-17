"""Holi/metaDMG matching and lineage-support logic for MetaMerge.

This module contains all the logic for:
  - Building fast lookup indexes from the metaDMG output.
  - Exact matching of MEGAN taxa to metaDMG rows (by tax_id or tax_name+rank).
  - QC label assignment for the metaDMG fit quality.
  - Lineage-support detection (conservative ancestor/descendant logic).

Exact vs. lineage matching
---------------------------
*Exact matching* requires a taxon to appear in the Holi/metaDMG output with
the same tax_id (preferred) or the same tax_name+rank (fallback), and to pass
the configured damage and significance thresholds.

*Lineage support* is a supplementary signal: if a lower-rank relative of the
focal taxon is exact-damage-supported in the same library, the focal taxon
receives indirect support.  This is deliberately conservative:
  - Only direct ancestor/descendant relationships count.
  - Broad clades (e.g. Mammalia, Boreoeutheria) are explicitly forbidden.
  - Only family rank or lower (configurable via ``lineage_max_rank_level``).
  - Within a limited number of rank steps (configurable via ``lineage_max_steps``).

QC labels
----------
MetaMerge assigns a QC label to the best Holi row for each taxon/library:
  - ``"clean"``         — mapping and fit quality within acceptable bounds.
  - ``"caution"``       — one QC metric is slightly elevated; use with care.
  - ``"strong caution"``— ``MAP_valid`` is False, or QC metrics exceed the
                          caution thresholds; treat damage estimate unreliable.
  - ``"not-applicable"``— no Holi/metaDMG data found for this taxon/library.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

from .utils import normalize_name, normalize_rank, rank_to_level


def build_holi_taxonomy_lookup(holi_df: pd.DataFrame) -> dict:
    """Build a canonical taxon-info lookup dict from the Holi/metaDMG table.

    Used to retrieve the authoritative tax_path and tax_rank for taxa that are
    found in the MEGAN count matrix but whose path is absent from MEGAN output.

    Args:
        holi_df: The loaded metaDMG DataFrame.

    Returns:
        Dict mapping either a ``tax_id_str`` or a normalized ``tax_name`` to a
        sub-dict ``{"tax_name": ..., "tax_rank": ..., "tax_path": ...}``.
    """
    info = {}
    group_cols = ["tax_id_str", "tax_name", "tax_rank", "tax_path"]
    for _, row in holi_df[group_cols].drop_duplicates().iterrows():
        key_id = row["tax_id_str"] if pd.notna(row["tax_id_str"]) else ""
        if key_id and key_id not in info:
            info[key_id] = {
                "tax_name": row["tax_name"],
                "tax_rank": row["tax_rank"],
                "tax_path": row["tax_path"],
            }
        name = row["tax_name"]
        if name and name not in info:
            info[name] = {
                "tax_name": row["tax_name"],
                "tax_rank": row["tax_rank"],
                "tax_path": row["tax_path"],
            }
    return info


def make_holi_exact_index(holi_df: pd.DataFrame) -> tuple[dict, dict]:
    """Create exact-match indexes for O(1) lookup by sample and taxon.

    Builds two indexes:
      - ``by_id[(sample, tax_id_str)]``         → list of matching rows
      - ``by_name_rank[(sample, name, rank)]``  → list of matching rows
      - ``by_name_rank[(sample, name, "")]``    → list of matching rows (rank-agnostic)

    Args:
        holi_df: The loaded metaDMG DataFrame.

    Returns:
        Tuple of (by_id, by_name_rank) dicts.
    """
    by_id        = defaultdict(list)
    by_name_rank = defaultdict(list)
    for row in holi_df.to_dict(orient="records"):
        sample = row["sample"]
        tid    = row.get("tax_id_str", "")
        name   = row.get("tax_name",   "")
        rank   = row.get("tax_rank",   "")
        if tid:
            by_id[(sample, tid)].append(row)
        if name:
            by_name_rank[(sample, name, rank)].append(row)
            by_name_rank[(sample, name, "")].append(row)
    return dict(by_id), dict(by_name_rank)


def compute_qc_label(
    n_reads: float,
    n_alignments: float,
    map_valid: bool,
    rho_ac: float,
    thresholds: dict,
) -> tuple[str, float]:
    """Assign a QC label to a metaDMG row based on mapping and fit quality.

    Checks three independent QC signals:
      1. ``MAP_valid`` — if False, the MAP optimisation did not converge.
      2. ``N_alignments / N_reads`` — high ratios indicate excessive multi-mapping.
      3. ``|rho_Ac|`` — high values indicate poor correlation between the fit and
         the data.

    Args:
        n_reads: metaDMG N_reads value.
        n_alignments: metaDMG N_alignments value.
        map_valid: Whether the MAP fit is valid (metaDMG MAP_valid field).
        rho_ac: Fit-correlation metric rho_Ac from metaDMG.
        thresholds: The "thresholds" sub-dict of the MetaMerge config.

    Returns:
        Tuple of ``(label_string, alignment_ratio)`` where ``label_string`` is
        one of ``"clean"``, ``"caution"``, or ``"strong caution"``.
    """
    align_ratio = (
        float(n_alignments) / float(n_reads)
        if n_reads and n_reads > 0
        else np.nan
    )
    abs_rho = (
        abs(float(rho_ac))
        if rho_ac is not None and not pd.isna(rho_ac)
        else np.nan
    )

    if not map_valid:
        return "strong caution", align_ratio

    if (
        (not np.isnan(align_ratio) and align_ratio > thresholds["qc_alignments_per_read_caution_max"])
        or (not np.isnan(abs_rho) and abs_rho > thresholds["qc_abs_rho_caution_max"])
    ):
        return "strong caution", align_ratio

    if (
        (not np.isnan(align_ratio) and align_ratio > thresholds["qc_alignments_per_read_clean_max"])
        or (not np.isnan(abs_rho) and abs_rho > thresholds["qc_abs_rho_clean_max"])
    ):
        return "caution", align_ratio

    return "clean", align_ratio


def choose_best_exact_row(rows: list[dict], thresholds: dict) -> dict | None:
    """Select the single best Holi exact-match row for a taxon/library pair.

    Scoring priority (lexicographic):
      1. Whether the row passes the damage+significance thresholds.
      2. Significance value (higher is better).
      3. N_reads (higher is better).
      4. Damage value (higher is better).

    Args:
        rows: Candidate metaDMG rows (all from the same sample for one taxon).
        thresholds: The "thresholds" sub-dict of the MetaMerge config.

    Returns:
        The best-scoring row dict, or ``None`` if ``rows`` is empty.
    """
    if not rows:
        return None

    def _score(row: dict) -> tuple:
        supported = int(
            pd.notna(row.get("damage"))
            and pd.notna(row.get("significance"))
            and row["damage"] > thresholds["damage_min"]
            and row["significance"] > thresholds["significance_min"]
        )
        return (
            supported,
            float(row.get("significance") or -np.inf),
            float(row.get("N_reads") or -np.inf),
            float(row.get("damage") or -np.inf),
        )

    return max(rows, key=_score)


def row_has_exact_damage_support(row: dict, thresholds: dict) -> bool:
    """Return ``True`` if a metaDMG row passes the damage/significance thresholds.

    Args:
        row: A single metaDMG row dict.
        thresholds: The "thresholds" sub-dict of the MetaMerge config.

    Returns:
        ``True`` if ``damage > damage_min`` AND ``significance > significance_min``.
    """
    if row is None:
        return False
    damage       = row.get("damage")
    significance = row.get("significance")
    return (
        pd.notna(damage)
        and pd.notna(significance)
        and damage > thresholds["damage_min"]
        and significance > thresholds["significance_min"]
    )


def lineage_path_names(path: str) -> list[str]:
    """Split a semicolon-delimited taxonomic path into normalized taxon names.

    Args:
        path: Semicolon-delimited lineage string from metaDMG tax_path.

    Returns:
        List of normalized, non-empty taxon name strings.
    """
    if not path:
        return []
    return [normalize_name(x) for x in str(path).split(";") if normalize_name(x)]


def is_meaningful_low_rank_lineage_support(
    focal_tax_name: str,
    focal_tax_rank: str,
    focal_path: str,
    candidate_tax_name: str,
    candidate_tax_rank: str,
    candidate_path: str,
    cfg: dict,
    forbid_set: set[str] | None = None,
    ordered_ranks: list[str] | None = None,
    max_level: int | None = None,
) -> bool:
    """Decide whether a candidate taxon provides meaningful lineage support.

    Applies conservative rules to prevent over-counting lineage support:
      - Both taxa must be below the ``lineage_max_rank_level``.
      - One taxon must be a direct ancestor/descendant of the other.
      - The rank gap must not exceed ``lineage_max_steps``.
      - Neither taxon may be in the ``forbid_broad_names`` list.
      - The two taxa must not be the same name.

    Args:
        focal_tax_name: Name of the taxon being evaluated (from MEGAN).
        focal_tax_rank: Rank of the focal taxon.
        focal_path: Semicolon-delimited lineage path of the focal taxon.
        candidate_tax_name: Name of the potential lineage-support taxon (from Holi).
        candidate_tax_rank: Rank of the candidate taxon.
        candidate_path: Semicolon-delimited lineage path of the candidate taxon.
        cfg: Full MetaMerge config dict.
        forbid_set: Pre-computed set of forbidden broad names (optional; built
            from cfg if not provided).
        ordered_ranks: Pre-computed ordered rank list (optional).
        max_level: Pre-computed maximum allowed rank level (optional).

    Returns:
        ``True`` if the candidate provides valid lineage support.
    """
    focal_tax_name     = normalize_name(focal_tax_name)
    candidate_tax_name = normalize_name(candidate_tax_name)
    focal_tax_rank     = normalize_rank(focal_tax_rank)
    candidate_tax_rank = normalize_rank(candidate_tax_rank)

    if not focal_tax_name or not candidate_tax_name:
        return False
    if focal_tax_name == candidate_tax_name:
        return False

    forbid = forbid_set or {normalize_name(x) for x in cfg["lineage"]["forbid_broad_names"]}
    if focal_tax_name in forbid or candidate_tax_name in forbid:
        return False

    ordered_ranks = ordered_ranks or [normalize_rank(x) for x in cfg["lineage"]["rank_levels"]]
    max_rank  = normalize_rank(cfg["thresholds"]["lineage_max_rank_level"])
    max_level = rank_to_level(max_rank, ordered_ranks) if max_level is None else max_level

    focal_level = rank_to_level(focal_tax_rank, ordered_ranks)
    cand_level  = rank_to_level(candidate_tax_rank, ordered_ranks)

    if focal_level is None or cand_level is None or max_level is None:
        return False
    if focal_level > max_level or cand_level > max_level:
        return False

    focal_lineage = lineage_path_names(focal_path)
    cand_lineage  = lineage_path_names(candidate_path)
    if not focal_lineage or not cand_lineage:
        return False

    focal_in_candidate  = focal_tax_name in cand_lineage
    candidate_in_focal  = candidate_tax_name in focal_lineage
    if not (focal_in_candidate or candidate_in_focal):
        return False

    if abs(focal_level - cand_level) > int(cfg["thresholds"]["lineage_max_steps"]):
        return False

    return True


def build_lineage_indexes(
    holi_df: pd.DataFrame,
    thresholds_cfg: dict,
) -> tuple[dict[str, list[dict]], dict[str, dict[str, list[dict]]], dict[str, dict[str, list[dict]]]]:
    """Pre-compute per-sample lineage-support indexes for fast lookup.

    Building these indexes once at the start of the merge loop (rather than
    re-scanning all Holi rows for every MEGAN taxon) reduces time complexity
    from O(taxa × holi_rows) to O(taxa + holi_rows).

    Args:
        holi_df: The loaded metaDMG DataFrame.
        thresholds_cfg: The "thresholds" sub-dict of the MetaMerge config.

    Returns:
        Tuple of three dicts:
          - ``exact_supported_rows_by_sample``:
              sample → list of all exact-damage-supported rows.
          - ``supported_name_index_by_sample``:
              sample → tax_name → list of rows for that name.
          - ``descendant_index_by_sample``:
              sample → lineage_ancestor_name → list of rows whose tax_path
              contains that ancestor name.
    """
    exact_by_sample          = defaultdict(list)
    name_index_by_sample     = defaultdict(lambda: defaultdict(list))
    descendant_index_by_sample = defaultdict(lambda: defaultdict(list))

    for row in holi_df.to_dict(orient="records"):
        if not row_has_exact_damage_support(row, thresholds_cfg):
            continue
        sample = row["sample"]
        exact_by_sample[sample].append(row)
        name = normalize_name(row.get("tax_name", ""))
        if name:
            name_index_by_sample[sample][name].append(row)
        for lineage_name in lineage_path_names(row.get("tax_path", "")):
            descendant_index_by_sample[sample][lineage_name].append(row)

    return (
        dict(exact_by_sample),
        {k: dict(v) for k, v in name_index_by_sample.items()},
        {k: dict(v) for k, v in descendant_index_by_sample.items()},
    )


def summarize_lineage_support_for_taxon(
    focal_tax_name: str,
    focal_tax_rank: str,
    focal_tax_path: str,
    holi_sample: str,
    supported_name_index_by_sample: dict[str, dict[str, list[dict]]],
    descendant_index_by_sample: dict[str, dict[str, list[dict]]],
    package_cfg: dict,
) -> tuple[list[str], list[str]]:
    """Find lineage-consistent exact-damage-supported taxa for a focal taxon.

    Searches for both:
      - *Ancestor support*: a taxon named in the focal taxon's lineage is
        exact-damage-supported in the same library.
      - *Descendant support*: a taxon whose lineage includes the focal taxon
        name is exact-damage-supported in the same library.

    Args:
        focal_tax_name: Name of the focal taxon (from MEGAN).
        focal_tax_rank: Rank of the focal taxon.
        focal_tax_path: Semicolon-delimited lineage path of the focal taxon.
        holi_sample: Holi/metaDMG sample name for this library.
        supported_name_index_by_sample: Pre-built index from build_lineage_indexes.
        descendant_index_by_sample: Pre-built index from build_lineage_indexes.
        package_cfg: Full MetaMerge config dict.

    Returns:
        Tuple of ``(libraries, examples)`` where:
          - ``libraries`` is a list of Holi sample names that provide lineage support.
          - ``examples`` is a list of human-readable description strings.
    """
    supported_libraries = []
    examples = []

    if not focal_tax_name:
        return supported_libraries, examples

    forbid_set    = {normalize_name(x) for x in package_cfg["lineage"]["forbid_broad_names"]}
    ordered_ranks = [normalize_rank(x) for x in package_cfg["lineage"]["rank_levels"]]
    max_rank      = normalize_rank(package_cfg["thresholds"]["lineage_max_rank_level"])
    max_level     = rank_to_level(max_rank, ordered_ranks)

    focal_lineage   = lineage_path_names(focal_tax_path)
    candidate_rows  = []

    # Ancestor support: candidate is named directly in the focal taxon's lineage.
    for lineage_name in focal_lineage:
        candidate_rows.extend(
            supported_name_index_by_sample.get(holi_sample, {}).get(lineage_name, [])
        )

    # Descendant support: the focal taxon name appears in the candidate's lineage.
    candidate_rows.extend(
        descendant_index_by_sample.get(holi_sample, {}).get(
            normalize_name(focal_tax_name), []
        )
    )

    # Deduplicate.
    seen    = set()
    deduped = []
    for row in candidate_rows:
        key = (row.get("sample"), row.get("tax_name"), row.get("tax_rank"))
        if key not in seen:
            seen.add(key)
            deduped.append(row)

    for row in deduped:
        if is_meaningful_low_rank_lineage_support(
            focal_tax_name=focal_tax_name,
            focal_tax_rank=focal_tax_rank,
            focal_path=focal_tax_path,
            candidate_tax_name=row.get("tax_name", ""),
            candidate_tax_rank=row.get("tax_rank", ""),
            candidate_path=row.get("tax_path", ""),
            cfg=package_cfg,
            forbid_set=forbid_set,
            ordered_ranks=ordered_ranks,
            max_level=max_level,
        ):
            lib = row["sample"]
            supported_libraries.append(lib)
            examples.append(f"{lib}: {row.get('tax_name')} ({row.get('tax_rank')})")

    return supported_libraries, examples

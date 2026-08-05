"""Fillet matching logic for MetaMerge -- the Fillet-specific analogue of
holi.py's Holi/metaDMG matching and index-building functions.

Fillet's own authentication signals (composite_authenticity, authenticity_tier,
confidence_tier, eco/pal/fos_support) are structurally different from Holi's
damage/significance/rho_Ac triple, so this is a genuinely separate module
rather than a variant of holi.py -- mirrors that module's own index-building /
exact-matching / per-taxon-signal-extraction shape so the two stay easy to
compare side by side, without conflating two different tools' evidence models
into one function.
"""

from __future__ import annotations

from collections import defaultdict

import pandas as pd

from .utils import normalize_name, normalize_rank


def make_fillet_exact_index(fillet_df: pd.DataFrame) -> tuple[dict, dict]:
    """Create exact-match indexes for O(1) lookup by sample and taxon.

    Mirrors holi.make_holi_exact_index()'s shape exactly.

    Args:
        fillet_df: The loaded Fillet MetaMerge evidence DataFrame (from
            io.load_fillet()).

    Returns:
        Tuple of (by_id, by_name_rank) dicts, same structure as
        holi.make_holi_exact_index().
    """
    by_id        = defaultdict(list)
    by_name_rank = defaultdict(list)
    for row in fillet_df.to_dict(orient="records"):
        sample = row["sample"]
        tid    = row.get("tax_id_str", "")
        name   = row.get("tax_name", "")
        rank   = row.get("tax_rank", "")
        if tid:
            by_id[(sample, tid)].append(row)
        if name:
            by_name_rank[(sample, name, rank)].append(row)
            by_name_rank[(sample, name, "")].append(row)
    return dict(by_id), dict(by_name_rank)


def build_fillet_taxonomy_lookup(fillet_df: pd.DataFrame) -> dict:
    """Build a canonical taxon-info lookup dict from Fillet's evidence table.

    Fillet's export has no tax_path column (unlike Holi's), so this lookup
    only carries tax_name/tax_rank -- callers needing a lineage path for a
    Fillet-only taxon should fall back to the Holi taxonomy lookup or leave
    tax_path blank, same as MEGAN-only taxa do today when Holi has no
    corresponding row. Mirrors holi.build_holi_taxonomy_lookup()'s shape.
    """
    info = {}
    group_cols = ["tax_id_str", "tax_name", "tax_rank"]
    for _, row in fillet_df[group_cols].drop_duplicates().iterrows():
        key_id = row["tax_id_str"] if pd.notna(row["tax_id_str"]) else ""
        if key_id and key_id not in info:
            info[key_id] = {"tax_name": row["tax_name"], "tax_rank": row["tax_rank"], "tax_path": ""}
        name = row["tax_name"]
        if name and name not in info:
            info[name] = {"tax_name": row["tax_name"], "tax_rank": row["tax_rank"], "tax_path": ""}
    return info


def row_is_fillet_authenticated(row: dict, thresholds: dict) -> bool:
    """Return True if a Fillet row clears Fillet's own authentication bar.

    Requires BOTH composite_authenticity >= fillet_composite_min AND
    authenticity_tier <= fillet_authenticity_tier_max (lower tier number =
    higher confidence, 1 = highest, 5 = lowest, 0 = rejected) -- matching
    Fillet's own "never collapse continuous and discrete evidence into one
    number" design: a taxon must clear both independently, not just a
    weighted average of the two.
    """
    if row is None:
        return False
    composite = row.get("composite_authenticity")
    tier      = row.get("authenticity_tier")
    if pd.isna(composite) or pd.isna(tier):
        return False
    return (
        composite >= thresholds["fillet_composite_min"]
        and 1 <= tier <= thresholds["fillet_authenticity_tier_max"]
    )


def choose_best_fillet_row(rows: list[dict], thresholds: dict) -> dict | None:
    """Select the single best Fillet exact-match row for a taxon/library pair.

    Scoring priority (lexicographic), mirroring holi.choose_best_exact_row():
      1. Whether the row passes Fillet's own authentication bar.
      2. composite_authenticity (higher is better).
      3. direct_hard_reads (higher is better).
    """
    if not rows:
        return None

    def _score(row: dict) -> tuple:
        authenticated = int(row_is_fillet_authenticated(row, thresholds))
        composite = row.get("composite_authenticity")
        reads = row.get("direct_hard_reads")
        return (
            authenticated,
            float(composite) if pd.notna(composite) else -1.0,
            float(reads) if pd.notna(reads) else -1.0,
        )

    return max(rows, key=_score)


def row_has_fillet_strong_count_support(row: dict, thresholds: dict) -> bool:
    """Return True if a Fillet row's own read count clears
    fillet_strong_count_min_reads -- Fillet's analogue of MEGAN's
    strong_count_support, used when no MEGAN count matrix is available
    (the Holi+Fillet-only case)."""
    if row is None:
        return False
    reads = row.get("direct_hard_reads")
    if pd.isna(reads):
        return False
    return float(reads) >= thresholds["fillet_strong_count_min_reads"]

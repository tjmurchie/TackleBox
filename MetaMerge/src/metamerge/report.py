"""Reporting helpers: long-format plotting tables for MetaMerge.

This module prepares plot-friendly TSVs for the R heatmap / stacked-bar layer.
Compared with earlier iterations, the logic here now does more of the heavy
lifting before the data ever reach R:

- technical libraries are collapsed to biological-sample plotting units
- per-library support is aggregated to the strongest support observed within a
  plotting unit
- linker metadata columns are carried forward automatically
- spiderweb / blanks / undetermined groups are separated more cleanly
- a family-like display group is derived directly from tax_path where possible
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any

import pandas as pd

from .utils import STATUS_PRIORITY

# All ranks that can ever appear in report tables.
REPORT_ALLOWED_RANKS = {
    "subspecies", "species", "species group", "genus", "subgenus",
    "subtribe", "tribe", "subfamily", "family",
}

# Depth map: 0 = most specific (subspecies), higher = more inclusive.
# min_rank_for_reports selects all ranks with depth <= this taxon's depth.
_RANK_DEPTH: dict[str, int] = {
    "subspecies": 0, "species": 1, "species group": 2,
    "subgenus": 3, "genus": 4, "subtribe": 5, "tribe": 6,
    "subfamily": 7, "family": 8, "superfamily": 9,
    "infraorder": 10, "suborder": 11, "order": 12,
    "superorder": 13, "class": 14, "subphylum": 15,
    "phylum": 16, "kingdom": 17, "domain": 18,
}

TAX_RANK_ORDER = {
    "family": 0, "subfamily": 1, "tribe": 2, "subtribe": 3,
    "genus": 4, "subgenus": 5, "species group": 6, "species": 7,
    "subspecies": 8,
}

LIB_SUPPORT_ORDER = {
    "Damage-supported (>=100 reads)": 0,
    "Damage-supported": 1,
    "Lineage-supported": 2,
    "Count-only": 3,
    "Blank-library": 4,
    "Not detected": 5,
}


def _short_plot_label(value: str) -> str:
    """Create a compact x-axis label for plotting units while preserving uniqueness as much as possible."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    s = str(value).strip().replace("_", "-")
    s = re.sub(r"-(?:E-PNmp|E-PNm|Shg)$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"-(?:Alpha|Beta|Gamma)$", "", s, flags=re.IGNORECASE)
    return s


def _parse_tax_path(path: str | float | None) -> list[tuple[str, str]]:
    """Parse metaDMG/Holi tax_path into (name, rank) tuples."""
    if path is None or (isinstance(path, float) and pd.isna(path)):
        return []
    out: list[tuple[str, str]] = []
    for part in str(path).split("	"):
        m = re.match(r'^[^:]+:([^:]+):"?([^":]+)"?$', part)
        if not m:
            continue
        out.append((m.group(1).strip(), m.group(2).strip().lower()))
    return out


def _display_group_from_path(path: str | float | None, broad_group: str | None) -> tuple[str, str]:
    """Pick a human-readable grouping label from tax_path."""
    lineage = _parse_tax_path(path)
    if not lineage:
        return (broad_group or "Other", "")

    prefs = {
        "animals": ["family", "subfamily", "tribe", "order", "class", "phylum"],
        "plants": ["family", "subfamily", "tribe", "order", "class", "phylum"],
        "fungi": ["family", "order", "class", "phylum"],
        "microbes": ["phylum", "class", "order", "family"],
        "other": ["family", "order", "class", "phylum"],
    }
    ranks = prefs.get(str(broad_group or "").lower(), prefs["other"])
    for pref in ranks:
        for name, rank in lineage:
            if rank == pref:
                return (name, pref)
    return (lineage[0][0], lineage[0][1])


def _normalize_plot_sample_id(sample_id: Any, biological_sample_id: Any, holi_library_name: Any, merged_library_name: Any, sample_type: Any, is_negative_control: Any) -> str:
    """Collapse technical-library names to a biological plotting unit."""
    def _s(x: Any) -> str:
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        return str(x).strip()

    neg = str(is_negative_control).lower() in {"true", "1", "yes"}
    stype = _s(sample_type).lower()
    sid = _s(sample_id)
    base_sid = _s(biological_sample_id)
    holi_name = _s(holi_library_name)
    merged = _s(merged_library_name)

    if base_sid:
        base = base_sid
    elif not neg and holi_name and holi_name not in {"blanks", "NO_HOLI_DATA"} and not holi_name.startswith("REVIEW_NEEDED"):
        base = holi_name
    else:
        base = sid or merged

    base = base.replace("_", "-")
    base = re.sub(r'-(?:Alpha|Beta|Gamma)$', '', base, flags=re.IGNORECASE)
    base = re.sub(r'([A-Za-z0-9])(?:[cp])$', r'\\1', base)
    base = base.replace('-Cob-', '-Cob').replace('-Web-', '-Web')
    base = re.sub(r'\s+', ' ', base).strip('- ')
    if "undetermined" in stype or base.lower().startswith("undetermined"):
        return "Undetermined"
    return base or merged or sid or holi_name or 'unknown'


def _derive_plot_group(row: pd.Series) -> str:
    """Choose a plotting subset/group label."""
    site = str(row.get("site") or "").strip()
    sample_type = str(row.get("sample_type") or "").strip().lower()
    group = str(row.get("group") or "").strip()
    sample_id = str(row.get("plot_sample_id") or row.get("sample_id") or "")
    neg = str(row.get("is_negative_control")).lower() in {"true", "1", "yes"}
    pos = str(row.get("is_positive_control", "false")).lower() in {"true", "1", "yes"}
    env = str(row.get("is_environmental_control", "false")).lower() in {"true", "1", "yes"}

    if env:
        return "environmental_controls"
    if neg or "blank" in sample_type or ("negative" in sample_type and "positive" not in sample_type):
        return "negative_controls"
    if pos or ("positive" in sample_type and "control" in sample_type):
        return "positive_controls"
    if "undetermined" in sample_type or sample_id.lower().startswith("undetermined"):
        return "undetermined"
    if "web" in sample_type:
        return f"{site}_webs" if site else "webs"
    if group and group.lower() not in {"nan", "none", "", "blanks", "positive_controls"}:
        return group
    if site and sample_type:
        st = sample_type.replace(" ", "_")
        if not st.endswith("s"):
            st = st + "s"
        return f"{site}_{st}"
    if site:
        return site
    return "subset"


def _first_non_null(series: pd.Series):
    for value in series:
        if pd.notna(value) and str(value).strip() != "":
            return value
    return pd.NA


def _collapse_metadata(series: pd.Series):
    vals = [str(v).strip() for v in series if pd.notna(v) and str(v).strip() != ""]
    if not vals:
        return pd.NA
    uniq = []
    for v in vals:
        if v not in uniq:
            uniq.append(v)
    if len(uniq) == 1:
        return uniq[0]
    return "; ".join(uniq)


def _best_support(values: pd.Series) -> str:
    vals = [str(v) for v in values.dropna().tolist() if str(v).strip() != ""]
    if not vals:
        return "Not detected"
    return sorted(set(vals), key=lambda x: LIB_SUPPORT_ORDER.get(x, 999))[0]


def filter_for_report(df: pd.DataFrame, broad_group: str, config: dict) -> pd.DataFrame:
    """Return all eligible taxa for a broad group.

    Rank filtering: ``config["report"]["max_rank_for_reports"]`` controls the
    *most inclusive* (broadest) rank allowed.  All ranks at that level or more
    specific (lower depth) are included.  Accepted values:

    * ``"species"`` — species, species group, subspecies only
    * ``"genus"`` — genus and below
    * ``"family"`` — family and below  *(default)*
    * ``"order"`` — order and below (includes family-level and more specific)
    * ``"class"`` — class and below
    * any rank string from the _RANK_DEPTH map

    Examples:  ``"order"`` includes Proboscidea alongside its constituent
    families (Elephantidae etc.).  ``"family"`` is the most useful default
    for publication figures.
    """
    min_rank = str(config["report"]["max_rank_for_reports"]).lower().strip()
    min_depth = _RANK_DEPTH.get(min_rank, _RANK_DEPTH["family"])
    allowed_ranks = {r for r, d in _RANK_DEPTH.items() if d <= min_depth}

    out = df.loc[df["broad_group"] == broad_group].copy()
    if out.empty:
        return out

    out["_status_priority"] = out["aDNA_support_status"].map(STATUS_PRIORITY).fillna(99)
    out["_rank_allowed"] = out["tax_rank"].fillna("").str.lower().isin(allowed_ranks)
    out = out.loc[out["_rank_allowed"]].copy()

    display = out.apply(lambda row: _display_group_from_path(row.get("tax_path"), row.get("broad_group")), axis=1)
    out["display_group"] = [x[0] for x in display]
    out["display_group_rank"] = [x[1] for x in display]

    out = out.sort_values(
        ["_status_priority", "megan_max_count", "display_group", "scientific_name"],
        ascending=[True, False, True, True],
    )
    return out.drop(columns=["_status_priority", "_rank_allowed"], errors="ignore")


def make_plot_input(df: pd.DataFrame, metadata: pd.DataFrame, broad_group: str, config: dict) -> pd.DataFrame:
    """Convert merged output into long-format plot input for one broad group."""
    subset = filter_for_report(df, broad_group=broad_group, config=config)
    if subset.empty:
        return subset

    count_cols   = [c for c in subset.columns if c.startswith("count__")]
    support_cols = [c for c in subset.columns if c.startswith("aDNA_support_lib__")]
    damage_cols  = [c for c in subset.columns if c.startswith("Holi_damage_lib__")]
    id_cols = [c for c in subset.columns if c not in count_cols and c not in support_cols and c not in damage_cols]

    long = subset.melt(id_vars=id_cols, value_vars=count_cols, var_name="plot_library", value_name="count")
    long["merged_library_name"] = long["plot_library"].str.replace(r"^count__", "", regex=True)

    if support_cols:
        sup_long = subset[["scientific_name"] + support_cols].melt(
            id_vars=["scientific_name"],
            value_vars=support_cols,
            var_name="_sup_col",
            value_name="library_adna_support",
        )
        sup_long["merged_library_name"] = sup_long["_sup_col"].str.replace(r"^aDNA_support_lib__", "", regex=True)
        long = long.merge(
            sup_long[["scientific_name", "merged_library_name", "library_adna_support"]],
            on=["scientific_name", "merged_library_name"],
            how="left",
        )
    else:
        long["library_adna_support"] = pd.NA

    if damage_cols:
        dmg_long = subset[["scientific_name"] + damage_cols].melt(
            id_vars=["scientific_name"],
            value_vars=damage_cols,
            var_name="_dmg_col",
            value_name="library_damage",
        )
        dmg_long["merged_library_name"] = dmg_long["_dmg_col"].str.replace(r"^Holi_damage_lib__", "", regex=True)
        long = long.merge(
            dmg_long[["scientific_name", "merged_library_name", "library_damage"]],
            on=["scientific_name", "merged_library_name"],
            how="left",
        )
    else:
        long["library_damage"] = pd.NA

    # carry all linker metadata columns forward
    meta_cols = [c for c in metadata.columns if c != "megan_library_name"]
    meta_for_merge = metadata[["merged_library_name"] + [c for c in meta_cols if c != "merged_library_name"]].copy()
    if not meta_for_merge.empty:
        long = long.merge(meta_for_merge, on="merged_library_name", how="left")

    long["count"] = pd.to_numeric(long["count"], errors="coerce").fillna(0)
    long = long.loc[long["count"] > 0].copy()
    if long.empty:
        return long
    long["plot_sample_id"] = long.apply(
        lambda row: _normalize_plot_sample_id(
            row.get("sample_id"),
            row.get("biological_sample_id"),
            row.get("holi_library_name"),
            row.get("merged_library_name"),
            row.get("sample_type"),
            row.get("is_negative_control"),
        ),
        axis=1,
    )
    long["plot_group"] = long.apply(_derive_plot_group, axis=1)

    # collapse technical libraries / chemistries to plotting units
    fixed_cols = [
        c for c in long.columns
        if c not in {"plot_library", "count", "library_adna_support", "library_damage",
                     "merged_library_name", "plot_library_label", "support_priority"}
    ]

    rows = []
    for (sci, plot_sample_id, plot_group), group in long.groupby(["scientific_name", "plot_sample_id", "plot_group"], dropna=False):
        dmg_vals = pd.to_numeric(group["library_damage"], errors="coerce")
        out = {
            "scientific_name": sci,
            "plot_sample_id": plot_sample_id,
            "plot_group": plot_group,
            "count": pd.to_numeric(group["count"], errors="coerce").fillna(0).sum(),
            "library_adna_support": _best_support(group["library_adna_support"]),
            "plot_damage": float(dmg_vals.max()) if dmg_vals.notna().any() else float("nan"),
            "merged_library_name": plot_sample_id,
            "plot_library_label": _short_plot_label(plot_sample_id),
            "n_technical_libraries": group["merged_library_name"].nunique(),
            "technical_libraries": "; ".join(sorted(group["merged_library_name"].dropna().astype(str).unique().tolist())),
        }
        for col in fixed_cols:
            if col in {"scientific_name", "plot_sample_id", "plot_group", "count", "library_adna_support", "merged_library_name", "plot_library_label"}:
                continue
            if col not in group.columns:
                out[col] = pd.NA
            elif col in {"technical_libraries", "notes", "metadata_notes", "context", "sediment_zone", "site", "group", "sample_type", "label", "label_base", "putative_taxon_common", "putative_taxon_scientific", "element", "portion", "museum", "extraction_method"}:
                out[col] = _collapse_metadata(group[col])
            else:
                out[col] = _first_non_null(group[col])
        rows.append(out)

    collapsed = pd.DataFrame(rows)
    if collapsed.empty:
        return collapsed

    collapsed["support_priority"] = collapsed["aDNA_support_status"].map(STATUS_PRIORITY).fillna(99)
    collapsed["library_support_priority"] = collapsed["library_adna_support"].map(LIB_SUPPORT_ORDER).fillna(99)

    # helpful ranking metrics for the R layer
    collapsed["overall_taxon_sum"] = collapsed.groupby("scientific_name")["count"].transform("sum")
    collapsed["overall_taxon_max"] = collapsed.groupby("scientific_name")["count"].transform("max")
    collapsed["subset_taxon_sum"] = collapsed.groupby(["plot_group", "scientific_name"])["count"].transform("sum")
    collapsed["subset_taxon_max"] = collapsed.groupby(["plot_group", "scientific_name"])["count"].transform("max")
    collapsed["display_group_sum"] = collapsed.groupby(["plot_group", "display_group"])["count"].transform("sum")
    collapsed["display_group_max"] = collapsed.groupby(["plot_group", "display_group"])["count"].transform("max")
    collapsed["display_group_status_priority"] = collapsed.groupby(["plot_group", "display_group"])["support_priority"].transform("min")
    collapsed["tax_rank_order"] = collapsed["tax_rank"].fillna("").str.lower().map(TAX_RANK_ORDER).fillna(999).astype(int)

    collapsed = collapsed.sort_values(
        [
            "plot_group",
            "display_group_status_priority",
            "display_group_sum",
            "display_group_max",
            "support_priority",
            "library_support_priority",
            "subset_taxon_sum",
            "subset_taxon_max",
            "tax_rank_order",
            "display_group",
            "scientific_name",
            "plot_sample_id",
        ],
        ascending=[True, True, False, False, True, True, False, False, True, True, True, True],
    )
    return collapsed


def write_plot_inputs(merged: pd.DataFrame, metadata: pd.DataFrame, outdir: Path, config: dict) -> list[Path]:
    """Write one long-format TSV per broad taxonomic group."""
    paths: list[Path] = []
    report_dir = outdir / "report_inputs"
    report_dir.mkdir(parents=True, exist_ok=True)

    for broad_group in ["animals", "plants", "fungi", "microbes"]:
        table = make_plot_input(merged, metadata, broad_group=broad_group, config=config)
        if table.empty:
            continue
        path = report_dir / f"{broad_group}_heatmap_input.tsv"
        table.to_csv(path, sep="	", index=False)
        paths.append(path)

    return paths

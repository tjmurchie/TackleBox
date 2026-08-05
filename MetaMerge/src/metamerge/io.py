"""Input reading for MEGAN count matrices and Holi/metaDMG tables.

MEGAN count matrix format
--------------------------
MetaMerge expects a **wide** count matrix with:
  - One row per taxon.
  - One column per sequencing library (labeled by the library name used in the
    MEGAN run).
  - An optional ``tax_id`` column (NCBI taxonomy ID).
  - An optional ``tax_name`` column (taxon name string).
  - An optional ``tax_rank`` column (e.g. species, genus, family).

MEGAN's TSV export places the taxon ID in a column named ``#Datasets``.
The package recognises this automatically via the column-alias list in
``defaults.py``.

Multiple MEGAN files
---------------------
``load_megan_counts`` accepts either a single file path or a list of paths.
When multiple files are provided (common when different sample groups were run
through MEGAN separately), they are merged on tax_id using a full outer join:
taxa present in some files but not others get a count of 0 for the missing
libraries.  All files must use the same taxon-ID column convention.

Holi / metaDMG CSV format
---------------------------
MetaMerge reads only the columns it needs from the (potentially very large)
metaDMG CSV using ``usecols``, so even 1 GB+ files load without memory issues.

Required metaDMG columns: sample, tax_id, tax_name, tax_rank, N_reads,
N_alignments, damage, significance, rho_Ac, MAP_valid, tax_path.

Fillet evidence table format
------------------------------
``load_fillet`` reads Fillet's own native MetaMerge export
(``fillet_metamerge_evidence.tsv``, written by
``fillet.report.write_fillet_metamerge_input``) -- a long table, one row per
taxon per sample, carrying Fillet's aggregated authenticity assessment
(composite_authenticity, authenticity_tier, confidence_tier) alongside the
individual evidence dimensions worth reasoning about independently
(damage, reference breadth/stacking, blank_fraction, eco/pal/fos support).
Fillet is an independent classifier, not a second copy of Holi or MEGAN, so
this is a genuinely separate input path, not a variant of load_holi/
load_megan_counts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import pandas as pd

from .utils import find_first_matching_column, normalize_name, normalize_rank, read_table


def _load_single_megan(path: Path, aliases: dict, sheet_name=None) -> pd.DataFrame:
    """Load one MEGAN TSV/CSV and standardise its taxon-ID columns.

    Returns a DataFrame with columns tax_id, tax_name (optional), tax_rank
    (optional), plus all library count columns from that file.
    """
    df = read_table(path, sheet_name=sheet_name)

    tax_id_col   = find_first_matching_column(df, aliases["tax_id"])
    tax_name_col = find_first_matching_column(df, aliases["tax_name"])
    tax_rank_col = find_first_matching_column(df, aliases["tax_rank"])

    if not tax_id_col and not tax_name_col:
        raise ValueError(
            f"Could not identify a taxon-ID or taxon-name column in {path.name}.  "
            "Searched for: " + str(aliases["tax_id"] + aliases["tax_name"]) + ".  "
            "Either rename the column or add an alias override to your config YAML."
        )

    rename = {}
    if tax_id_col:
        rename[tax_id_col] = "tax_id"
    if tax_name_col:
        rename[tax_name_col] = "tax_name"
    if tax_rank_col:
        rename[tax_rank_col] = "tax_rank"

    return df.rename(columns=rename).copy()


def load_megan_counts(
    path: Union[str, list[str]],
    metadata: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    """Load one or more MEGAN wide count matrices and standardise taxon columns.

    When multiple files are supplied they are merged on tax_id (full outer join)
    so that taxa appearing in only some files receive a count of 0 for the
    libraries in other files.  This handles the common case where different
    sample groups were run through MEGAN separately.

    The function:
      1. Reads each file (xlsx/csv/tsv auto-detected).
      2. Detects taxon-ID, taxon-name, and taxon-rank columns via alias lists.
      3. Merges multiple files on tax_id (outer join, counts filled with 0).
      4. Checks that every library referenced in the linker metadata has a
         matching column in the merged matrix.
      5. Normalises taxon names and ranks.
      6. Returns a tidy DataFrame with only the columns needed for the merge.

    Args:
        path: Path to the MEGAN count matrix file, or a list of paths when
            libraries are spread across multiple MEGAN TSV outputs.
        metadata: Validated linker DataFrame (from ``load_metadata``).
        config: MetaMerge config dict.

    Returns:
        DataFrame with columns: tax_id, tax_name, tax_rank, tax_id_str, and
        one column per library listed in ``metadata["megan_library_name"]``.

    Raises:
        ValueError: If no taxon-ID or taxon-name column can be found, or if a
            library referenced in the linker is missing from the merged matrix.
    """
    aliases    = config["column_aliases"]
    sheet_name = config["io"].get("counts_sheet")

    paths = [path] if isinstance(path, str) else list(path)

    # Load each file independently and collect the DataFrames.
    frames = [_load_single_megan(Path(p), aliases, sheet_name) for p in paths]

    if len(frames) == 1:
        df = frames[0]
    else:
        # Merge on tax_id with a full outer join.  Fill missing library counts
        # with 0 (taxon not detected in that MEGAN run → 0 reads).
        df = frames[0]
        for other in frames[1:]:
            # Identify taxon key columns present in both frames.
            merge_on = [c for c in ["tax_id", "tax_name", "tax_rank"]
                        if c in df.columns and c in other.columns]
            if not merge_on:
                merge_on = ["tax_id"] if "tax_id" in df.columns else ["tax_name"]
            df = pd.merge(df, other, on=merge_on, how="outer")

    # Keep only relevant count columns (fill NaN from outer join with 0).
    for col in metadata["megan_library_name"]:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    # Validate that every linker library has a matching column.
    missing_libs = [col for col in metadata["megan_library_name"] if col not in df.columns]
    if missing_libs:
        raise ValueError(
            f"The following {len(missing_libs)} linker library name(s) could not be "
            f"found in any of the supplied MEGAN count file(s):\n"
            + "\n".join(f"  {lib}" for lib in missing_libs[:10])
            + ("\n  …" if len(missing_libs) > 10 else "")
            + "\n\nCheck that megan_library_name values in your linker file "
            "exactly match the column headers in the MEGAN TSV(s)."
        )

    # Keep only the taxon columns plus the library count columns.
    keep = [c for c in ["tax_id", "tax_name", "tax_rank"] if c in df.columns]
    keep += metadata["megan_library_name"].tolist()
    df = df[keep].copy()

    # Ensure all three taxon columns exist (fill missing ones with pd.NA).
    for col in ["tax_id", "tax_name", "tax_rank"]:
        if col not in df.columns:
            df[col] = pd.NA

    df["tax_name"]   = df["tax_name"].map(normalize_name)
    df["tax_rank"]   = df["tax_rank"].map(normalize_rank)
    df["tax_id_str"] = df["tax_id"].astype("string").str.strip()

    # Ensure library count columns are numeric (coerce errors to 0.0).
    for col in metadata["megan_library_name"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(float)

    return df


def load_fillet(path: str, config: dict) -> pd.DataFrame:
    """Load Fillet's native MetaMerge evidence table (fillet_metamerge_evidence.tsv,
    written by fillet.report.write_fillet_metamerge_input).

    Unlike load_holi() (which reads a large third-party metaDMG CSV and therefore
    uses usecols for memory efficiency), Fillet's export is already a small,
    purpose-built table -- read via the generic read_table() helper (tsv/csv/xlsx
    auto-detected) rather than a hardcoded pd.read_csv, so it tolerates the same
    format flexibility as the MEGAN/linker loaders.

    Args:
        path: Path to Fillet's exported evidence table.
        config: MetaMerge config dict (provides fillet_required_columns).

    Returns:
        DataFrame with the required columns present, normalised and typed the
        same way load_holi() normalises Holi's sample/tax_name/tax_rank/tax_id_str.

    Raises:
        ValueError: If any required column is absent from the file.
    """
    required = config["fillet_required_columns"]
    df = read_table(Path(path))

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "Fillet MetaMerge export is missing required columns: "
            + ", ".join(missing)
            + ".  Required columns: " + ", ".join(required)
        )

    df["tax_name"]   = df["tax_name"].map(normalize_name)
    df["tax_rank"]   = df["tax_rank"].map(normalize_rank)
    df["tax_id_str"] = df["tax_id"].astype("string").str.strip()
    df["sample"]     = df["sample"].astype(str).str.strip()

    for col in ("eco_support", "pal_support", "fos_support"):
        df[col] = df[col].map(
            lambda v: str(v).strip().lower() in ("true", "1", "yes") if pd.notna(v) else False
        )

    numeric_cols = [
        "direct_hard_reads", "direct_weighted_reads",
        "cumulative_hard_reads", "cumulative_weighted_reads",
        "composite_authenticity", "authenticity_tier",
        "mean_damage_score", "max_damage_score", "terminus_ct_5p", "terminus_ga_3p",
        "best_reference_breadth", "mean_reference_breadth",
        "stack_concentration", "n_covered_windows", "coverage_gini",
        "blank_fraction", "n_support_lines", "barcode_fraction", "mean_hit_uniqueness",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "confidence_tier" in df.columns:
        df["confidence_tier"] = df["confidence_tier"].fillna("").astype(str)
    if "flags" in df.columns:
        df["flags"] = df["flags"].fillna("").astype(str)

    return df


def load_holi(path: str, config: dict) -> pd.DataFrame:
    """Load a Holi/metaDMG table, reading only the columns MetaMerge needs.

    Using ``usecols`` makes this fast even for multi-gigabyte metaDMG CSVs.

    Args:
        path: Path to the metaDMG/Holi CSV.
        config: MetaMerge config dict (provides the required column list).

    Returns:
        DataFrame with exactly the required columns, normalised and typed.

    Raises:
        ValueError: If any required column is absent from the file.
    """
    usecols = config["holi_required_columns"]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)

    missing = [c for c in usecols if c not in df.columns]
    if missing:
        raise ValueError(
            "Holi/metaDMG file is missing required columns: "
            + ", ".join(missing)
            + ".  Required columns: " + ", ".join(usecols)
        )

    df["tax_name"]   = df["tax_name"].map(normalize_name)
    df["tax_rank"]   = df["tax_rank"].map(normalize_rank)
    df["tax_id_str"] = df["tax_id"].astype("string").str.strip()
    df["sample"]     = df["sample"].astype(str).str.strip()
    df["tax_path"]   = df["tax_path"].fillna("").astype(str)
    df["MAP_valid"]  = df["MAP_valid"].fillna(False).astype(bool)

    for col in ["N_reads", "N_alignments", "damage", "significance", "rho_Ac"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df

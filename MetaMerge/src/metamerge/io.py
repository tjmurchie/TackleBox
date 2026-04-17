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

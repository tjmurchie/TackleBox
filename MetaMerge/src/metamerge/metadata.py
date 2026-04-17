"""Metadata / library-linker parsing and validation for MetaMerge.

The *library linker* (also called the metadata file) is what makes MetaMerge
project-agnostic.  It tells the merge engine:

  - Which MEGAN count-matrix column corresponds to which Holi/metaDMG sample.
  - Which libraries are negative controls.
  - Which libraries are positive controls.
  - What project metadata (site, context, group, etc.) to carry through to the
    merged workbook.

Required columns in the linker file
-------------------------------------
megan_library_name
    The *exact* column header used in the MEGAN count matrix.  This must match
    character-for-character, including the long file-path suffix that MEGAN
    sometimes writes (e.g. ``PFX-SITE-S01-c-Shg_S60_PE.mapped.min24...``).

holi_library_name
    The sample name as it appears in the ``sample`` column of the Holi/metaDMG
    CSV.  One Holi sample may cover multiple MEGAN libraries (e.g., all
    extraction replicates of one biological sample are pooled on the Holi side).

merged_library_name
    A human-readable, unique name used to label this library in workbook output
    and heatmaps.  Must be unique within the linker file.

Recommended optional columns
------------------------------
sample_id, sample_type, is_negative_control, is_positive_control,
site, depth, age, group, notes

Any of these not present in the linker will simply be absent from the merged
workbook.  The package recognises them using the alias lists in defaults.py so
minor column-name differences are handled automatically.

How environmental controls are detected
-----------------------------------------
A library is treated as an environmental control (separate from both
biological samples and extraction blanks) if *any* of these are true:

  1. The ``is_environmental_control`` column (if present) is truthy.
  2. The ``sample_type`` value contains "environmental".
  3. The ``sample_type`` value contains "spiderweb".
  4. The ``sample_type`` value contains "env_control" or "env_ctrl".

Environmental controls are detected before negative controls so that
samples with "environmental_control" in sample_type (which contains
"control") are not misclassified as extraction blanks.

How negative controls are detected
-------------------------------------
A library is treated as a negative control if *any* of these are true:

  1. The ``is_negative_control`` column (if present) is truthy (True, 1, yes, t).
  2. The ``sample_type`` value (if present) contains "blank" (catches air blanks,
     extraction blanks, library blanks, etc.).
  3. The ``sample_type`` value (if present) contains "negative".
  4. The ``sample_type`` value (if present) contains "control" but NOT "positive".

Note: positive controls and environmental controls are never flagged as negative.

Positive controls are NOT treated as negative controls and will never drive
Blank-associated classification.  They are tracked separately via
``is_positive_control`` so they can be distinguished from true biological
samples in outputs and reports.

How positive controls are detected
--------------------------------------
A library is treated as a positive control if *any* of these are true:

  1. The ``is_positive_control`` column (if present) is truthy.
  2. The ``sample_type`` value contains "positive" (catches "positive control",
     "positive ctrl", etc., case-insensitive).
  3. The ``control_type`` column (written by ``metamerge linker``) equals
     "positive_control".
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .utils import find_first_matching_column, read_table, safe_bool


# These three columns are the absolute minimum required in every linker file.
REQUIRED_METADATA_FIELDS = ["megan_library_name", "holi_library_name", "merged_library_name"]


def load_metadata(path: str, config: dict) -> pd.DataFrame:
    """Load and validate the MetaMerge library linker file.

    Reads the linker CSV/Excel, applies column-name aliases, normalises the
    ``is_negative_control`` flag, and checks that all required fields are
    present and that ``merged_library_name`` values are unique.

    Rows where all fields are NaN (e.g. trailing blank rows in an Excel file)
    are silently dropped before validation.

    Args:
        path: Path to the linker file (CSV, TSV, or Excel).
        config: MetaMerge config dict (used for column aliases and sheet name).

    Returns:
        Validated DataFrame with at minimum the three required columns plus
        ``is_negative_control`` and ``sample_type``.

    Raises:
        ValueError: If required columns are missing or ``merged_library_name``
            contains duplicates.
    """
    df = read_table(Path(path), sheet_name=config["io"].get("metadata_sheet"))

    # Drop completely-empty rows (common in Excel files with trailing blank rows).
    df = df.dropna(how="all").reset_index(drop=True)

    aliases = config["column_aliases"]

    # Build rename map from alias lists.
    rename = {}
    for target in REQUIRED_METADATA_FIELDS + [
        "sample_id", "sample_type", "is_negative_control", "is_positive_control",
        "site", "age", "depth", "group", "notes",
    ]:
        col = find_first_matching_column(df, aliases.get(target, [target]))
        if col:
            rename[col] = target

    df = df.rename(columns=rename).copy()

    # Check required fields.
    missing = [field for field in REQUIRED_METADATA_FIELDS if field not in df.columns]
    if missing:
        raise ValueError(
            "Metadata/linker file is missing required columns: "
            + ", ".join(missing)
            + ".  Required columns are: "
            + ", ".join(REQUIRED_METADATA_FIELDS)
            + ".  See docs/metadata_linker.md for the full linker specification."
        )

    # Ensure optional columns have defaults.
    if "is_negative_control" not in df.columns:
        df["is_negative_control"] = False

    if "is_positive_control" not in df.columns:
        df["is_positive_control"] = False

    if "sample_type" not in df.columns:
        df["sample_type"] = ""

    # Normalise is_positive_control using the explicit flag, sample_type, and
    # the control_type column written by `metamerge linker` when present.
    def _is_positive(row_flag, row_type: str, ctrl_type: str) -> bool:
        row_type_lc = str(row_type).lower()
        ctrl_lc = str(ctrl_type).lower()
        return (
            safe_bool(row_flag)
            or "positive" in row_type_lc
            or ctrl_lc == "positive_control"
        )

    ctrl_type_col = df["control_type"] if "control_type" in df.columns else [""] * len(df)
    df["is_positive_control"] = [
        _is_positive(flag, stype, ctype)
        for flag, stype, ctype in zip(df["is_positive_control"], df["sample_type"], ctrl_type_col)
    ]

    # Detect environmental controls BEFORE negative controls, so that samples
    # with "environmental_control" in their sample_type (which contains "control")
    # are not misclassified as extraction/library blanks.
    if "is_environmental_control" not in df.columns:
        df["is_environmental_control"] = False

    def _is_environmental(row_flag, row_type: str, is_pos: bool) -> bool:
        """Environmental controls (field samples, spiderwebs, water controls, etc.)
        are neither biological samples nor extraction blanks."""
        if is_pos:
            return False
        row_type_lc = str(row_type).lower()
        return (
            safe_bool(row_flag)
            or "environmental" in row_type_lc
            or "spiderweb" in row_type_lc
            or "env_control" in row_type_lc
            or "env_ctrl" in row_type_lc
        )

    df["is_environmental_control"] = [
        _is_environmental(flag, stype, is_pos)
        for flag, stype, is_pos in zip(
            df["is_environmental_control"], df["sample_type"], df["is_positive_control"]
        )
    ]

    # Normalise is_negative_control using both the explicit flag and sample_type.
    # A positive control or environmental control must never be flagged as a negative control.
    def _is_negative(row_flag, row_type: str, is_pos: bool, is_env: bool) -> bool:
        if is_pos or is_env:
            return False
        row_type_lc = str(row_type).lower()
        return (
            safe_bool(row_flag)
            or "blank" in row_type_lc
            or "negative" in row_type_lc
            or ("control" in row_type_lc and "positive" not in row_type_lc)
        )

    df["is_negative_control"] = [
        _is_negative(flag, stype, is_pos, is_env)
        for flag, stype, is_pos, is_env in zip(
            df["is_negative_control"], df["sample_type"],
            df["is_positive_control"], df["is_environmental_control"]
        )
    ]

    # Validate uniqueness of merged_library_name.
    if df["merged_library_name"].duplicated().any():
        dups = df.loc[df["merged_library_name"].duplicated(), "merged_library_name"].tolist()
        raise ValueError(
            "merged_library_name values must be unique across the linker file.  "
            "Duplicates found: " + ", ".join(map(str, dups[:10]))
        )

    return df

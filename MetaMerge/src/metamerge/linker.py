#!/usr/bin/env python3
"""
make_metamerge_linker.py — MetaMerge library linker generator
==============================================================

Generates the ``library_linker.csv`` required by MetaMerge's ``check`` and
``run`` commands from the project metadata spreadsheet and the MEGAN count-
matrix TSV file(s).

The *linker* is the cross-walk that tells MetaMerge:
  1. Which MEGAN column (long file-path name) belongs to which biological sample.
  2. Which Holi/metaDMG sample name corresponds to each MEGAN library.
  3. Which libraries are negative controls.
  4. What project metadata (site, sample type, context …) to carry through to
     the merged workbook.

MINIMUM REQUIRED INPUTS (for any new project)
----------------------------------------------
1. A metadata table containing at least:
     * seq_library_id_Shg   Stem name of the first-chemistry sequenced library
       (or equivalent column; see --shg-col).
     * seq_library_id_Enr   Stem name of the second-chemistry library (optional;
       see --enr-col).  Not all libraries have both.
     * sample_id            Biological sample identifier.  This must match or
       closely correspond to the ``sample`` column in the metaDMG/Holi CSV
       (the script will attempt to strip trailing -c / -p extraction suffixes
       automatically).
     * sample_type          Sample category string.  Any value containing
       "blank" is automatically treated as a negative control.

2. The MEGAN count-matrix TSV/CSV file(s) whose column headers are the full
   library file names.

OPTIONAL INPUTS
---------------
* --holi        Holi/metaDMG CSV (for validation — confirms that derived holi_library_name values appear in the Holi sample column).
* --seq-summary Sequencing summary TSV(s) with per-library QC metrics.
* --holi-map    A two-column CSV (megan_stem, holi_sample_name) for manual
                overrides when automatic derivation fails.
* --out         Output linker CSV path (default: library_linker.csv).

HOW holi_library_name IS DERIVED
---------------------------------
1. If the sample_type contains "blank" (case-insensitive) → "blanks"
   (all blanks are pooled in the metaDMG step)
2. If a manual override is provided via --holi-map → use it
3. If sample_id matches a known Holi sample name exactly → use it
4. Strip trailing extraction suffix (-c / -p) and any resulting trailing
   hyphen → check (e.g. PROJ-S01c → PROJ-S01)
5. Normalise underscores to hyphens → check (e.g. PROJ-La_S1 → PROJ-La-S1)
6. Strip leading zeros from numeric suffix → check (e.g. SITE-04 → SITE-4)
7. Combinations of steps 4+5, 4+6, 5+6
8. Otherwise → "REVIEW_NEEDED:<sample_id>" (flags for manual review)

Note: Undetermined libraries and Beta-replicate libraries are handled before
this derivation step (see above).

HOW MEGAN STEMS ARE MATCHED TO METADATA
-----------------------------------------
MEGAN file-path column headers are stripped to a library stem (e.g.
PFX-SITE-S01c-Shg after lane tag removal) and looked up in the metadata by
seq_library_id_Shg / seq_library_id_Enr.  If no match is found, the script
tries the following fallbacks in order:

  1. Convert underscores to hyphens (La_S1c → La-S1c).
  2. Strip a "-Beta" extraction-replicate tag from the stem.  Some sequencing
     runs produce duplicate extractions labelled "-Beta-" that share a Holi
     sample with the primary extraction (e.g. PFX-SITE-S01-p-Beta-Shg is
     treated identically to PFX-SITE-S01-p-Shg because metaDMG merged them).

UNDETERMINED INDEXES
---------------------
Libraries whose column header contains "Undetermined" are automatically assigned
sample_type="undetermined_indexes" with no Holi data (holi_library_name set to
"NO_HOLI_DATA").  These contribute to MEGAN count columns in the merge output
but will never receive Holi damage support.  Remove them from the linker
entirely before running MetaMerge if you prefer to exclude undetermined reads.

The REVIEW_NEEDED entries must be resolved before running MetaMerge.

HOW merged_library_name IS DERIVED
------------------------------------
The library stem (MEGAN column name after stripping the long file-path suffix)
with the batch prefix (AiL-, BiL-, iL-) removed.  This is unique per library
and human-readable.

USAGE EXAMPLES
--------------
  # Typical project with MEGAN directory and Holi/metaDMG CSV:
  python make_metamerge_linker.py \\
      --metadata project_metadata.xlsx \\
      --megan    megan_counts/ \\
      --holi     metadmg_output.csv \\
      --out      library_linker.csv

  # Single MEGAN file and no metaDMG for validation:
  python make_metamerge_linker.py \\
      --metadata project_metadata.xlsx \\
      --megan    megan_counts.tsv \\
      --out      library_linker.csv

  # With manual holi-name overrides for edge cases:
  python make_metamerge_linker.py \\
      --metadata project_metadata.xlsx \\
      --megan    megan_counts.tsv \\
      --holi-map manual_holi_overrides.csv \\
      --out      library_linker.csv


CONTROL / BLANK IDENTIFICATION
------------------------------
MetaMerge linker distinguishes four control types:

  negative_control      — wet-lab blanks: air blanks, extraction blanks,
                          library blanks.  Any sample_type containing
                          "blank", "negative", "air blank",
                          "extraction blank", or "library blank"
                          (case-insensitive).  These drive the
                          Blank-associated classification in the merge.

  positive_control      — known-positive reference samples (e.g. a bear
                          specimen of known species, a mock community).
                          Any sample_type containing "positive" or
                          "positive control".  These are NOT treated as
                          blanks and will NEVER drive Blank-associated
                          calls.  They are tracked separately so they
                          can be distinguished from true biological
                          samples in outputs and reports.

  environmental_control — webs, washes, swabs, or surface samples used
                          as environmental references.  Not blank-driving.

  control_unspecified   — any other sample_type containing "control"
                          without a positive/negative qualifier.
                          Treated conservatively as a negative control.

The generated ``control_type`` column records the exact classification.
``is_negative_control`` is True only for negative_control and
control_unspecified rows.

NAME INTERPRETATION / HOLI MATCHING
-----------------------------------
The linker does **not** assume all projects name libraries the same way.
It derives ``holi_library_name`` using a small sequence of conservative
normalizations so the audit trail stays transparent:

1. Exact match between metadata ``sample_id`` and a Holi/metaDMG
   ``sample`` name.
2. Strip common extraction/library suffixes (for example trailing
   ``-c`` / ``-p`` or similar chemistry suffixes).
3. Convert underscores to hyphens.
4. Normalize some leading-zero patterns (for example ``SITE-04`` → ``SITE-4``).
5. Fall back to ``REVIEW_NEEDED`` if no safe automatic match is found.

The generated ``holi_derivation`` / ``how_matched`` column records which
rule was used. If your project uses different naming conventions, the
recommended solution is to:
- standardize ``sample_id`` upstream, or
- edit the linker after generation, or
- supply an override/holi-map file when supported by your workflow.

OUTPUT COLUMNS
--------------
megan_library_name    Exact MEGAN column header (must match character for character).
holi_library_name     Corresponding sample name in the metaDMG/Holi ``sample`` column.
merged_library_name   Human-readable unique name for the merged workbook and plots.
sample_id             Biological sample ID from the metadata.
sample_type           sediment / bone / blank / control / etc.
is_negative_control   True / False.
site                  Site/location code (e.g. SITE1, SITE2).
context               Archaeological context string (e.g. cave, marsh).
group                 Broad sample group for report ordering.
chemistry             Sequencing chemistry (Shg, E-PNmp, E-PNm, etc.).
notes                 Any notes from the metadata.
holi_derivation       How the holi_library_name was derived (for audit).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


# ─── Column-name configuration ───────────────────────────────────────────────
# Adjust these if your metadata uses different column names.

DEFAULT_SHG_COL   = "seq_library_id_Shg"    # First-chemistry library ID column
DEFAULT_ENR_COL   = "seq_library_id_Enr"    # Second-chemistry library ID column
DEFAULT_SAMPLE_ID = "sample_id"             # Biological sample ID column
DEFAULT_SITE_COL  = "site"                  # Site/location column
DEFAULT_TYPE_COL  = "sample_type"           # Sample type column
DEFAULT_CONTEXT   = "context"               # Archaeological context column
DEFAULT_GREEK     = "greek_replicate"       # Extraction replicate (alpha/beta/gamma)

# Pattern that strips the long MEGAN file-path suffix from a column name.
# Matches everything from _PE.mapped (or .mapped for PE/SE) onwards.
MEGAN_SUFFIX_RE = re.compile(
    r"[_.]PE\.mapped.*|[_.]SE\.mapped.*|[_.]mapped\..*", re.IGNORECASE
)

# Pattern to strip the lane/index tag (_S{digits}) from the end of a library stem.
LANE_TAG_RE = re.compile(r"_S\d+$")

# Batch prefixes to strip from library stems when building merged_library_name.
# Adjust this pattern for your project's batch-prefix convention, or set to
# r"^$" to disable stripping entirely.
BATCH_PREFIX_RE = re.compile(r"^(?:AiL|BiL|iL)-")

# Extraction-method suffixes appended to sample_id.
# Typical values: -c (cold-spin), -p (PowerLyzer), or similar single-letter codes.
# The suffix may be preceded by a hyphen/underscore separator OR follow directly
# after a digit (e.g. PROJ-S01c, PROJ-La_S1c → both end with 'c' after a digit).
EXTRACTION_SUFFIX_RE = re.compile(r"[cp]$", re.IGNORECASE)


# ─── Helper functions ─────────────────────────────────────────────────────────

def strip_megan_suffix(col: str) -> str:
    """Strip the MEGAN file-path suffix from a column header.

    Turns ``PFX-SITE-S01-c-Shg_S60_PE.mapped.min24.adapt-filtered...`` into
    ``PFX-SITE-S01-c-Shg_S60``.

    Args:
        col: Full MEGAN column header string.

    Returns:
        Library stem with suffix removed.
    """
    return MEGAN_SUFFIX_RE.sub("", col).rstrip(".")


def strip_lane_tag(stem: str) -> str:
    """Strip the lane/index number tag from a library stem.

    Turns ``PFX-SITE-S01-c-Shg_S60`` into ``PFX-SITE-S01-c-Shg``.

    Args:
        stem: Library stem (after stripping the MEGAN suffix).

    Returns:
        Library ID without the ``_S{N}`` lane tag.
    """
    return LANE_TAG_RE.sub("", stem)


def clean_merged_name(stem: str) -> str:
    """Derive a clean human-readable merged_library_name from a library stem.

    Strips batch prefixes (configured via BATCH_PREFIX_RE) to reduce
    verbosity while keeping the name unique and informative.

    Args:
        stem: Library stem (``PFX-SITE-S01-c-Shg_S60`` or similar).

    Returns:
        Clean merged library name (``SITE-S01-c-Shg_S60``).
    """
    return BATCH_PREFIX_RE.sub("", stem)


def detect_chemistry(lib_id: str) -> str:
    """Detect the sequencing chemistry from a library ID string.

    Args:
        lib_id: Library ID as it appears in the metadata (e.g. ``PFX-SITE-S01-c-Shg``).

    Returns:
        Chemistry short string: ``"Shg"``, ``"E-PNmp"``, ``"E-PNm"``, or ``"unknown"``.
    """
    lib_lower = lib_id.lower()
    if "shg" in lib_lower:
        return "Shg"
    if "e-pnmp" in lib_lower:
        return "E-PNmp"
    if "e-pnm" in lib_lower:
        return "E-PNm"
    return "unknown"


def _strip_leading_zeros(s: str) -> str:
    """Strip leading zeros from hyphen-separated numeric suffixes.

    Converts ``SITE-04`` → ``SITE-4``, ``SITE-09`` → ``SITE-9``, leaving
    non-numeric or already-minimal strings unchanged.

    Args:
        s: Sample ID string.

    Returns:
        String with leading zeros removed from trailing numeric component.
    """
    return re.sub(r"-0+(\d)", r"-\1", s)




def normalize_biological_sample_id(sample_id: str) -> str:
    """Normalise a metadata sample_id to a biological plotting/aggregation ID."""
    if not sample_id:
        return ""
    s = str(sample_id).strip().replace("_", "-")
    s = re.sub(r"-(?:Alpha|Beta|Gamma)$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"([A-Za-z0-9])(?:[cp])$", r"\1", s)
    return s.strip("- ")


def classify_control_type(sample_id: str, sample_type: str, site: str = "", context: str = "") -> str:
    """Classify whether a record is a negative, positive, environmental, or non-control sample.

    The previous logic treated any plain ``sample_type == control`` as a negative control.
    That is too aggressive because some projects use generic labels like ``control`` for
    positive controls (for example a known-positive reference sample).

    Returns one of: ``negative_control``, ``positive_control``, ``environmental_control``,
    ``control_unspecified``, or ``not_control``.
    """
    sid = str(sample_id or "").strip().lower()
    stype = str(sample_type or "").strip().lower()
    site = str(site or "").strip().lower()
    context = str(context or "").strip().lower()
    hay = " | ".join([sid, stype, site, context])

    # Positive controls: explicit wording or common positive-control naming patterns.
    if (
        "positive" in hay
        or re.search(r"(?:pos(?:itive)?[-_ ]?ctrl|positive[-_ ]?control)", hay)
    ):
        return "positive_control"

    # Environmental/context controls such as webs, washes, swabs, etc.  These are not
    # negative wet-lab blanks and should not drive Blank-associated calls.
    if any(tok in hay for tok in ["web", "wash", "swab", "surface", "environmental control"]):
        return "environmental_control"

    # Negative controls / blanks.
    if any(tok in hay for tok in ["blank", "negative", "air blank", "extraction blank", "library blank"]):
        return "negative_control"
    if "control" in stype or "control" in context:
        return "control_unspecified"

    return "not_control"


def infer_plot_group(site: str, sample_type: str, context: str, biological_sample_id: str, control_type: str | None = None) -> str:
    """Derive a plotting subset/group label."""
    site = (site or "").strip()
    stype = (sample_type or "").strip().lower()
    sid = (biological_sample_id or "").strip()
    control_type = control_type or classify_control_type(sid, sample_type, site, context)

    if control_type == "negative_control" or control_type == "control_unspecified":
        return "negative_controls"
    if control_type == "positive_control":
        return "positive_controls"
    if control_type == "environmental_control" and ("web" in stype or "web" in sid.lower()):
        return f"{site}_webs" if site else "webs"
    if "undetermined" in stype or sid.lower().startswith("undetermined"):
        return "undetermined"
    if "web" in stype:
        return f"{site}_webs" if site else "webs"
    return _infer_group(site, sample_type, context)

def load_metadata_table(
    path: Path,
    shg_col: str,
    enr_col: str,
    sample_id_col: str,
    site_col: str,
    type_col: str,
    context_col: str,
    greek_col: str,
) -> tuple[pd.DataFrame, dict[str, dict]]:
    """Load the project metadata spreadsheet and build a library-stem lookup.

    Args:
        path: Path to the metadata xlsx/csv.
        shg_col: Column name for first-chemistry library IDs.
        enr_col: Column name for second-chemistry library IDs.
        sample_id_col: Column name for biological sample IDs.
        site_col: Column name for site/location.
        type_col: Column name for sample type.
        context_col: Column name for archaeological context.
        greek_col: Column name for extraction replicate (alpha/beta/gamma).

    Returns:
        Tuple of:
          - The raw metadata DataFrame (after dropping all-NaN rows).
          - A dict mapping library-ID-stem → row info dict.
    """
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        df = pd.read_excel(path, sheet_name=0)
    else:
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(path, sep=sep)

    # Drop trailing blank rows (common in Excel files).
    df = df.dropna(how="all").reset_index(drop=True)

    lookup: dict[str, dict] = {}

    def _to_str(val) -> str:
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return ""
        return str(val).strip()

    for _, row in df.iterrows():
        sample_id = _to_str(row.get(sample_id_col, ""))
        sample_type = _to_str(row.get(type_col, ""))
        site = _to_str(row.get(site_col, ""))
        context = _to_str(row.get(context_col, ""))
        greek = _to_str(row.get(greek_col, ""))

        info = {
            "sample_id":   sample_id,
            "sample_type": sample_type,
            "site":        site,
            "context":     context,
            "greek":       greek,
        }

        # Carry all remaining metadata fields through so optional project
        # information (age, depth, zone, museum ID, etc.) can flow into the
        # linker, merged workbook, and plotting layer automatically.
        for col, val in row.items():
            if col in {shg_col, enr_col}:
                continue
            if col not in info:
                info[col] = val

        # Index by both chemistry library stems.
        for col in [shg_col, enr_col]:
            lib_id = _to_str(row.get(col, ""))
            if lib_id:
                lookup[lib_id] = info

    return df, lookup


def load_megan_columns(paths: list[Path]) -> list[tuple[str, str, Path]]:
    """Read the header row of each MEGAN TSV to collect all library columns.

    Args:
        paths: List of MEGAN count-matrix TSV/CSV file paths.

    Returns:
        List of tuples: ``(full_column_name, library_stem, source_file_path)``.
        The first column (tax_id / ``#Datasets``) is excluded.
    """
    entries = []
    seen_stems = set()
    for fpath in paths:
        with fpath.open("r", encoding="utf-8-sig") as handle:
            header = handle.readline().rstrip("\n")

        sep = "\t" if "\t" in header else ","
        cols = header.split(sep)
        # Skip the first column (#Datasets / tax_id).
        for col in cols[1:]:
            col = col.strip()
            if not col:
                continue
            stem = strip_megan_suffix(col)
            if stem not in seen_stems:
                seen_stems.add(stem)
                entries.append((col, stem, fpath))

    return entries


def derive_holi_name(
    sample_id: str,
    sample_type: str,
    known_holi_samples: set[str],
    manual_overrides: dict[str, str],
) -> tuple[str, str]:
    """Derive the metaDMG/Holi sample name for a library.

    Rules applied in priority order (see module docstring for details):
      1. Blank → "blanks"
      2. Manual override → use it
      3. Exact match against known Holi samples → use sample_id
      4. Strip -c/-p extraction suffix (and trailing hyphen) and check
      5. Normalise underscores to hyphens (PROJ-La_S1 → PROJ-La-S1)
      6. Strip leading zeros from numeric suffix (SITE-04 → SITE-4)
      7. Combinations of 4+5, 4+6, 5+6
      8. Positive control with no Holi match → "NO_HOLI_DATA"
         (positive controls may not be in metaDMG output; MEGAN counts
         are still included in the merge but no damage support is
         applied — this is expected behaviour, not an error)
      9. Fallback → "REVIEW_NEEDED:<sample_id>"

    Args:
        sample_id: Biological sample ID from the metadata.
        sample_type: Sample type string (e.g. "sediment", "blank", "bone").
        known_holi_samples: Set of sample names present in the metaDMG CSV.
        manual_overrides: Dict mapping sample_id → holi_sample_name for cases
            that need manual specification.

    Returns:
        Tuple of ``(holi_library_name, derivation_method_string)``.
    """
    type_lower = str(sample_type).lower()
    is_positive_ctrl = "positive" in type_lower

    # Rule 1 — blanks are always pooled.
    if "blank" in type_lower:
        return "blanks", "blank-sample-type"

    # Rule 2 — manual override.
    if sample_id in manual_overrides:
        return manual_overrides[sample_id], "manual-override"

    # Rule 3 — exact match.
    if sample_id in known_holi_samples:
        return sample_id, "exact-match"

    # Build a list of candidate Holi names by applying normalizations in
    # priority order and checking each against known_holi_samples.
    #
    # Transformations applied (alone or combined):
    #   suffix  — strip trailing extraction letter (-c / -p), then strip any
    #             trailing hyphen left behind (e.g. "PROJ-S01c" → "PROJ-S01").
    #   hyphen  — convert underscores to hyphens (e.g. "PROJ-La_S1" → "PROJ-La-S1").
    #   nozero  — strip leading zeros from numeric suffix (SITE-04 → SITE-4).

    # --- Suffix-stripped form ---
    stripped = EXTRACTION_SUFFIX_RE.sub("", sample_id).rstrip("-")
    stripped = stripped if stripped != sample_id else None  # None if unchanged

    # --- Underscore-normalised form ---
    hyphened = sample_id.replace("_", "-")
    hyphened = hyphened if hyphened != sample_id else None

    # --- Leading-zero-stripped form ---
    nozero = _strip_leading_zeros(sample_id)
    nozero = nozero if nozero != sample_id else None

    # Candidate list: (candidate_string, method_label)
    candidates: list[tuple[str, str]] = []

    # Rule 4: extraction suffix stripped
    if stripped is not None:
        candidates.append((stripped, "stripped-extraction-suffix"))

    # Rule 5: underscore → hyphen
    if hyphened is not None:
        candidates.append((hyphened, "normalised-underscores"))

    # Rule 6: leading zeros stripped
    if nozero is not None:
        candidates.append((nozero, "stripped-leading-zeros"))

    # Combined: suffix + hyphen normalisation (covers PROJ-La_S1c → PROJ-La-S1)
    if stripped is not None:
        sh = stripped.replace("_", "-")
        if sh not in {sample_id, stripped}:
            candidates.append((sh, "stripped-suffix-normalised-underscores"))

    # Combined: suffix + leading zeros
    if stripped is not None:
        sn = _strip_leading_zeros(stripped)
        if sn not in {sample_id, stripped}:
            candidates.append((sn, "stripped-suffix-and-leading-zeros"))

    # Combined: hyphen normalisation + leading zeros
    if hyphened is not None:
        hn = _strip_leading_zeros(hyphened)
        if hn not in {sample_id, hyphened}:
            candidates.append((hn, "normalised-underscores-and-leading-zeros"))

    for candidate, method in candidates:
        if candidate and candidate in known_holi_samples:
            return candidate, method

    # Rule 8 — positive controls with no Holi match get NO_HOLI_DATA, not
    # REVIEW_NEEDED.  Positive controls are often not run through metaDMG
    # (no expected ancient damage signal), so a missing match is expected
    # behaviour.  MEGAN counts are still included in the merge.
    if is_positive_ctrl:
        return "NO_HOLI_DATA", "positive-control-no-holi-match"

    # Rule 9 — flag for manual review.
    return f"REVIEW_NEEDED:{sample_id}", "no-match-found"


def load_holi_samples(path: Path) -> set[str]:
    """Read unique sample names from a metaDMG/Holi CSV.

    Only reads the ``sample`` column to avoid loading the full (potentially
    very large) CSV into memory.

    Args:
        path: Path to the metaDMG CSV.

    Returns:
        Set of unique sample name strings.
    """
    try:
        df = pd.read_csv(path, usecols=["sample"], low_memory=False)
        return set(df["sample"].dropna().astype(str).str.strip().unique())
    except Exception as exc:
        print(f"Warning: could not read metaDMG file for validation: {exc}", file=sys.stderr)
        return set()


def load_holi_map(path: Optional[Path]) -> dict[str, str]:
    """Load a manual holi-name override CSV (two columns: megan_stem, holi_sample_name).

    Args:
        path: Path to the override CSV, or None to skip.

    Returns:
        Dict mapping library stem or sample_id → holi_sample_name.
    """
    if path is None:
        return {}
    df = pd.read_csv(path)
    if "megan_stem" not in df.columns or "holi_sample_name" not in df.columns:
        print(
            "Warning: --holi-map CSV must have columns 'megan_stem' and 'holi_sample_name'.  "
            "File will be ignored.",
            file=sys.stderr,
        )
        return {}
    return dict(zip(df["megan_stem"].astype(str), df["holi_sample_name"].astype(str)))


# ─── Main ─────────────────────────────────────────────────────────────────────

def build_linker(
    meta_path: Path,
    megan_paths: list[Path],
    metadmg_path: Optional[Path],
    holi_map_path: Optional[Path],
    shg_col: str,
    enr_col: str,
    sample_id_col: str,
    site_col: str,
    type_col: str,
    context_col: str,
    greek_col: str,
    verbose: bool = True,
) -> pd.DataFrame:
    """Build the MetaMerge library linker DataFrame.

    Args:
        meta_path: Path to the project metadata spreadsheet.
        megan_paths: List of MEGAN TSV file paths.
        metadmg_path: Optional path to the metaDMG CSV for validation.
        holi_map_path: Optional path to a manual holi-name override CSV.
        shg_col: Metadata column for first-chemistry library IDs.
        enr_col: Metadata column for second-chemistry library IDs.
        sample_id_col: Metadata column for biological sample IDs.
        site_col: Metadata column for site/location.
        type_col: Metadata column for sample type.
        context_col: Metadata column for archaeological context.
        greek_col: Metadata column for extraction replicate.
        verbose: Whether to print progress messages.

    Returns:
        Complete linker DataFrame ready to write as CSV.
    """
    if verbose:
        print(f"Reading metadata from: {meta_path}")
    _, stem_lookup = load_metadata_table(
        meta_path, shg_col, enr_col, sample_id_col,
        site_col, type_col, context_col, greek_col,
    )

    if verbose:
        print(f"Library stems in metadata: {len(stem_lookup)}")

    # Build a secondary lookup keyed by the extraction-suffix-stripped sample_id.
    # Used to resolve Beta-replicate libraries that share a Holi sample with the
    # primary extraction but may use a different extraction letter.
    # e.g. PFX-SITE-S01-c-Shg has sample_id=SITE-S01-c → stripped key = "SITE-S01"
    _stripped_sample_lookup: dict[str, dict] = {}
    for _info in stem_lookup.values():
        _sid = _info.get("sample_id", "")
        _stripped_sid = EXTRACTION_SUFFIX_RE.sub("", _sid).rstrip("-") if _sid else ""
        if _stripped_sid and _stripped_sid not in _stripped_sample_lookup:
            _stripped_sample_lookup[_stripped_sid] = _info

    # Regex for stripping chemistry suffixes so we can recover the bare sample_id.
    _CHEM_SUFFIX_RE = re.compile(r"[-_](?:Shg|E-PNmp?|E-PNm)$", re.IGNORECASE)

    known_holi: set[str] = set()
    if metadmg_path:
        if verbose:
            print(f"Reading Holi/metaDMG sample names from: {metadmg_path}")
        known_holi = load_holi_samples(metadmg_path)
        if verbose:
            print(f"  Found {len(known_holi)} unique metaDMG samples.")

    manual_overrides = load_holi_map(holi_map_path)
    if verbose and manual_overrides:
        print(f"  Loaded {len(manual_overrides)} manual holi-name overrides.")

    if verbose:
        print(f"Reading MEGAN column headers from {len(megan_paths)} file(s)…")
    megan_entries = load_megan_columns(megan_paths)
    if verbose:
        print(f"  Found {len(megan_entries)} unique MEGAN library columns.")

    rows = []
    unmatched_libs = []

    for full_col, stem, source_file in megan_entries:
        lib_id = strip_lane_tag(stem)  # e.g. PFX-SITE-S01-c-Shg

        # ── Undetermined indexes ────────────────────────────────────────────
        # Libraries whose column header includes "Undetermined" are sequencing
        # artefacts (unassigned index reads).  They have no Holi data, so we
        # assign them a sentinel holi_library_name and skip the metadata lookup.
        if "undetermined" in lib_id.lower() or "undetermined" in full_col.lower():
            merged_name = clean_merged_name(stem)
            rows.append({
                "megan_library_name":  full_col,
                "holi_library_name":   "NO_HOLI_DATA",
                "merged_library_name": merged_name,
                "sample_id":           "Undetermined",
                "sample_type":         "undetermined_indexes",
                "is_negative_control": "false",
                "site":                "NA",
                "context":             "",
                "group":               "undetermined",
                "chemistry":           detect_chemistry(lib_id),
                "notes":               "Unassigned index reads — no Holi/metaDMG data",
                "holi_derivation":     "undetermined-indexes",
            })
            continue

        # ── Normal metadata lookup ─────────────────────────────────────────
        info = stem_lookup.get(lib_id)

        if info is None:
            # Fallback 1: MEGAN file paths sometimes use underscores where the
            # metadata spreadsheet uses hyphens (e.g. La_S11c vs La-S11c).
            alt_id = lib_id.replace("_", "-")
            if alt_id != lib_id:
                info = stem_lookup.get(alt_id)

        if info is None and "-Beta-" in lib_id:
            # Fallback 2: Strip a "-Beta" extraction-replicate tag.  Some
            # projects produce Beta-replicate libraries that are merged with
            # the primary extraction on the Holi/metaDMG side.
            # e.g. PFX-SITE-S01-p-Beta-Shg → PFX-SITE-S01-p-Shg
            beta_stripped = lib_id.replace("-Beta-", "-")
            info = stem_lookup.get(beta_stripped)
            if info is None:
                info = stem_lookup.get(beta_stripped.replace("_", "-"))
            if info is None:
                # The Beta library may use a different extraction letter than the
                # metadata entry.  Derive the core sample_id by stripping the
                # batch prefix, chemistry suffix, and extraction letter, then
                # look up against the stripped-sample-id secondary table.
                no_prefix  = BATCH_PREFIX_RE.sub("", beta_stripped)
                no_chem    = _CHEM_SUFFIX_RE.sub("", no_prefix)
                core_sid   = EXTRACTION_SUFFIX_RE.sub("", no_chem).rstrip("-")
                if core_sid:
                    info = _stripped_sample_lookup.get(core_sid)

        if info is None:
            unmatched_libs.append((full_col, lib_id, source_file.name))
            # Still include the library with placeholder values so the user can
            # see it and fill in the metadata manually.
            info = {
                "sample_id":   "UNKNOWN",
                "sample_type": "UNKNOWN",
                "site":        "UNKNOWN",
                "context":     "",
                "greek":       "",
            }

        sample_id   = info["sample_id"]
        sample_type = info["sample_type"]
        site        = info["site"]
        context     = info["context"]
        chemistry   = detect_chemistry(lib_id)

        # Derive holi_library_name using the four-rule hierarchy.
        holi_name, derivation = derive_holi_name(
            sample_id=sample_id,
            sample_type=sample_type,
            known_holi_samples=known_holi,
            manual_overrides=manual_overrides,
        )

        control_type = classify_control_type(sample_id, sample_type, site, context)
        is_negative_control = control_type == "negative_control"

        biological_sample_id = normalize_biological_sample_id(sample_id)
        group = infer_plot_group(site, sample_type, context, biological_sample_id, control_type=control_type)

        merged_name = clean_merged_name(stem)

        row_out = {
            "megan_library_name":  full_col,
            "holi_library_name":   holi_name,
            "merged_library_name": merged_name,
            "sample_id":           sample_id,
            "biological_sample_id": biological_sample_id,
            "sample_type":         sample_type,
            "control_type":        control_type,
            "is_negative_control": str(is_negative_control).lower(),
            "site":                site,
            "context":             context,
            "group":               group,
            "chemistry":           chemistry,
            "notes":               "",
            "holi_derivation":     derivation,
        }

        for key, val in info.items():
            if key in {"sample_id", "sample_type", "site", "context", "greek"}:
                continue
            if key not in row_out:
                row_out[key] = val

        rows.append(row_out)

    linker = pd.DataFrame(rows)

    # Report issues.
    if unmatched_libs:
        print(
            f"\nWARNING: {len(unmatched_libs)} MEGAN library column(s) could not be "
            f"matched to the metadata.  These rows have sample_id=UNKNOWN and will need "
            f"manual correction before running MetaMerge:",
            file=sys.stderr,
        )
        for col, lib_id, fname in unmatched_libs[:20]:
            print(f"  [{fname}]  stem: {lib_id}", file=sys.stderr)

    review_needed = linker[linker["holi_library_name"].str.startswith("REVIEW_NEEDED")]
    if not review_needed.empty:
        print(
            f"\nWARNING: {len(review_needed)} row(s) could not have a Holi sample name "
            f"derived automatically.  Edit the 'holi_library_name' column in the linker "
            f"file before running MetaMerge.  Affected sample_ids:",
            file=sys.stderr,
        )
        for sid in review_needed["sample_id"].unique()[:20]:
            print(f"  {sid}", file=sys.stderr)

    if verbose:
        n_neg  = linker["is_negative_control"].eq("true").sum()
        n_pos  = (linker.get("control_type", pd.Series(dtype=str)) == "positive_control").sum()
        n_env  = (linker.get("control_type", pd.Series(dtype=str)) == "environmental_control").sum()
        n_real = len(linker) - n_neg - n_pos - n_env
        print(f"\nLinker rows: {len(linker)}")
        print(f"  Real samples:       {n_real}")
        print(f"  Negative controls:  {n_neg}  (air blanks / extraction blanks / library blanks)")
        if n_pos:
            print(f"  Positive controls:  {n_pos}")
        if n_env:
            print(f"  Environmental ctrl: {n_env}")
        if not review_needed.empty:
            print(f"  REVIEW_NEEDED rows: {len(review_needed)}")

    return linker


def _infer_group(site: str, sample_type: str, context: str) -> str:
    """Infer a broad group string for report ordering from site/type/context.

    Derives a ``{site}_{type}`` label from the sample_type string.  The
    result is used as the ``group`` column in the linker and drives plot
    subsetting in the R report scripts.

    Args:
        site: Site/location code from the metadata (e.g. ``"SITE1"``).
        sample_type: Sample type string (e.g. ``"sediment"``, ``"bone"``).
        context: Archaeological context string.

    Returns:
        A group string such as ``"SITE1_sediments"``, ``"SITE2_bones"``,
        ``"blanks"``, ``"positive_controls"``, or ``"other"``.
    """
    type_lower = str(sample_type).lower()
    site_str   = str(site).strip() if site and str(site).strip() not in {"", "NA", "nan"} else ""

    if "blank" in type_lower:
        return "blanks"
    if "control" in type_lower and "positive" in type_lower:
        return "positive_controls"
    if "control" in type_lower:
        return "blanks"

    if "bone" in type_lower:
        return f"{site_str}_bones" if site_str else "bones"
    if "sediment" in type_lower or "web" in type_lower:
        return f"{site_str}_sediments" if site_str else "sediments"

    return f"{site_str}" if site_str else "other"


def resolve_megan_paths(inputs: list[str]) -> list[Path]:
    """Resolve a list of file/directory paths to MEGAN TSV file paths.

    Directories are expanded to all ``.tsv`` files they contain.

    Args:
        inputs: List of file or directory path strings from the CLI.

    Returns:
        Sorted list of resolved MEGAN TSV file Paths.

    Raises:
        SystemExit: If any path does not exist.
    """
    paths = []
    for inp in inputs:
        p = Path(inp)
        if not p.exists():
            print(f"Error: path does not exist: {p}", file=sys.stderr)
            sys.exit(1)
        if p.is_dir():
            found = sorted(p.glob("*.tsv")) + sorted(p.glob("*.csv"))
            if not found:
                print(f"Warning: no .tsv/.csv files found in directory: {p}", file=sys.stderr)
            paths.extend(found)
        else:
            paths.append(p)
    return paths


def _trunc(val: object, width: int) -> str:
    """Return a string representation of val, truncated with … if over width."""
    s = str(val) if val is not None else ""
    return s if len(s) <= width else s[: width - 1] + "…"


def print_linker_report(linker: pd.DataFrame, out_path: Path) -> None:
    """Print a human-readable summary of the generated linker.

    Covers:
      - Overall row counts and resolution status
      - Details on any REVIEW_NEEDED / UNKNOWN rows
      - A short preview table of representative rows
      - Guidance on what to fix and how, plus next steps
    """
    SEP  = "=" * 66
    SEP2 = "-" * 66

    n_total   = len(linker)
    n_blanks  = linker["is_negative_control"].eq("true").sum()
    n_pos_ctrl = (linker.get("control_type", pd.Series(dtype=str)) == "positive_control").sum()
    n_env_ctrl = (linker.get("control_type", pd.Series(dtype=str)) == "environmental_control").sum()
    n_real    = n_total - n_blanks - n_pos_ctrl - n_env_ctrl

    review_rows  = linker[linker["holi_library_name"].str.startswith("REVIEW_NEEDED", na=False)]
    unknown_rows = linker[linker["sample_id"].eq("UNKNOWN")]
    no_holi_rows = linker[linker["holi_library_name"].eq("NO_HOLI_DATA")]
    n_review     = len(review_rows)
    n_unknown    = len(unknown_rows)

    print(f"\n{SEP}")
    print(f"  Linker written to : {out_path}")
    print(f"  Total rows        : {n_total}")
    print(f"    Real samples    : {n_real}")
    print(f"    Negative ctrl   : {n_blanks}  (air / extraction / library blanks)")
    if n_pos_ctrl:
        print(f"    Positive ctrl   : {n_pos_ctrl}")
    if n_env_ctrl:
        print(f"    Environ ctrl    : {n_env_ctrl}")
    if no_holi_rows is not None and len(no_holi_rows) > 0:
        undetermined = no_holi_rows[no_holi_rows["holi_derivation"].eq("undetermined-indexes")]
        other_no_holi = no_holi_rows[~no_holi_rows["holi_derivation"].eq("undetermined-indexes")]
        if len(undetermined) > 0:
            print(f"  No Holi data      : {len(undetermined)} undetermined-index librar{'y' if len(undetermined)==1 else 'ies'} (expected)")
        if len(other_no_holi) > 0:
            print(f"  No Holi data      : {len(other_no_holi)} other librar{'y' if len(other_no_holi)==1 else 'ies'} (check these)")
    print(SEP)

    # ── Resolution status ─────────────────────────────────────────────────────
    print()
    if n_review == 0 and n_unknown == 0:
        print("  STATUS: All libraries resolved automatically.")
        print("          No manual edits needed before running MetaMerge.")
    else:
        issues = []
        if n_review:
            issues.append(f"{n_review} holi_library_name value(s) flagged REVIEW_NEEDED")
        if n_unknown:
            issues.append(f"{n_unknown} row(s) with sample_id=UNKNOWN (metadata not found)")
        print(f"  STATUS: {len(issues)} issue type(s) need attention before running MetaMerge.")
        for issue in issues:
            print(f"          - {issue}")

    # ── REVIEW_NEEDED details ─────────────────────────────────────────────────
    if n_review > 0:
        print(f"\n{SEP2}")
        print(f"  REVIEW NEEDED — {n_review} row(s) could not be auto-resolved")
        print(SEP2)
        print(
            "  These libraries could not be matched to a Holi/metaDMG sample name.\n"
            "  Each is stored as REVIEW_NEEDED:<sample_id> in the holi_library_name\n"
            "  column and must be corrected before running MetaMerge.\n"
        )
        for i, (_, r) in enumerate(review_rows.head(15).iterrows(), start=1):
            sid  = r["sample_id"]
            holi = r["holi_library_name"]
            mlib = r["merged_library_name"]
            stype = r["sample_type"]
            print(f"  {i}. merged_library_name : {mlib}")
            print(f"     sample_id           : {sid}")
            print(f"     sample_type         : {stype}")
            print(f"     holi_library_name   : {holi}  ← needs replacing")
            print()
        if n_review > 15:
            print(f"  … and {n_review - 15} more (see the CSV for the full list).")

    # ── UNKNOWN sample_id details ─────────────────────────────────────────────
    if n_unknown > 0:
        print(f"\n{SEP2}")
        print(f"  UNKNOWN METADATA — {n_unknown} MEGAN column(s) not in metadata spreadsheet")
        print(SEP2)
        print(
            "  The MEGAN column stem could not be matched to any row in your metadata.\n"
            "  Check that the library name in MEGAN exactly matches the seq_library_id\n"
            "  column in your spreadsheet (hyphens vs underscores, batch prefixes, etc.).\n"
        )
        for _, r in unknown_rows.head(10).iterrows():
            print(f"  - {r['merged_library_name']}")
        if n_unknown > 10:
            print(f"  … and {n_unknown - 10} more.")

    # ── Sample preview table ──────────────────────────────────────────────────
    print(f"\n{SEP2}")
    print("  Sample preview — representative rows (megan_library_name omitted,")
    print("  it is stored correctly in the CSV but is too long to display here)")
    print(SEP2)

    # Select a cross-section of rows: one per site/group, one blank, REVIEW rows.
    selected = []
    seen_merged = set()

    def _add(row):
        if row["merged_library_name"] not in seen_merged:
            seen_merged.add(row["merged_library_name"])
            selected.append(row)

    # One real-sample row per group
    real = linker[linker["is_negative_control"].eq("false") &
                  ~linker["holi_library_name"].str.startswith("REVIEW_NEEDED", na=False) &
                  ~linker["sample_id"].eq("UNKNOWN") &
                  ~linker["holi_library_name"].eq("NO_HOLI_DATA")]
    for group in real["group"].dropna().unique():
        g_rows = real[real["group"] == group]
        if not g_rows.empty:
            _add(g_rows.iloc[0])

    # One blank
    blank_sample = linker[linker["is_negative_control"].eq("true")]
    if not blank_sample.empty:
        _add(blank_sample.iloc[0])

    # One positive control (if any)
    if "control_type" in linker.columns:
        pos_ctrl = linker[linker["control_type"].eq("positive_control")]
        if not pos_ctrl.empty:
            _add(pos_ctrl.iloc[0])

    # One undetermined (if any)
    undet = linker[linker["holi_derivation"].eq("undetermined-indexes")]
    if not undet.empty:
        _add(undet.iloc[0])

    # Any REVIEW_NEEDED rows (up to 3)
    for _, r in review_rows.head(3).iterrows():
        _add(r)

    # Column widths for display
    W = {"merged": 26, "sample_id": 16, "type": 14, "holi": 22, "ctrl": 5, "how": 30}
    header = (
        f"  {'merged_library_name':<{W['merged']}}  "
        f"{'sample_id':<{W['sample_id']}}  "
        f"{'sample_type':<{W['type']}}  "
        f"{'holi_library_name':<{W['holi']}}  "
        f"{'ctrl':<{W['ctrl']}}  "
        f"{'how_matched'}"
    )
    print(header)
    print("  " + "-" * (W["merged"] + W["sample_id"] + W["type"] + W["holi"] + W["ctrl"] + W["how"] + 12))

    for row in selected:
        if str(row["is_negative_control"]).lower() == "true":
            ctrl = "NEG"
        elif str(row.get("control_type", "")).lower() == "positive_control":
            ctrl = "POS"
        else:
            ctrl = "-"
        print(
            f"  {_trunc(row['merged_library_name'], W['merged']):<{W['merged']}}  "
            f"{_trunc(row['sample_id'], W['sample_id']):<{W['sample_id']}}  "
            f"{_trunc(row['sample_type'], W['type']):<{W['type']}}  "
            f"{_trunc(row['holi_library_name'], W['holi']):<{W['holi']}}  "
            f"{ctrl:<{W['ctrl']}}  "
            f"{_trunc(row['holi_derivation'], W['how'])}"
        )

    # ── Next steps ────────────────────────────────────────────────────────────
    print(f"\n{SEP2}")
    print("  Next steps")
    print(SEP2)

    if n_review > 0 or n_unknown > 0:
        print(
            "  !! Fix the issues above before running MetaMerge.\n"
        )
        print(
            "  How to edit the linker\n"
            "  ----------------------\n"
            f"  Open {out_path} in Excel, LibreOffice Calc, or a text editor.\n"
            "\n"
            "  For each REVIEW_NEEDED row:\n"
            "    - Find the row by its sample_id (column D).\n"
            "    - Replace the holi_library_name value (column B) with the exact\n"
            "      sample name as it appears in the 'sample' column of your\n"
            "      metaDMG/Holi CSV.\n"
            "    - If the sample has no metaDMG data (e.g. a positive control that\n"
            "      was not processed by metaDMG), set holi_library_name to:\n"
            "        NO_HOLI_DATA\n"
            "      The library will still contribute MEGAN counts to the merge but\n"
            "      will receive no damage support.  Note: positive controls whose\n"
            "      sample_type contains 'positive' are assigned NO_HOLI_DATA\n"
            "      automatically rather than REVIEW_NEEDED.\n"
            "\n"
            "  For each UNKNOWN sample_id row:\n"
            "    - The MEGAN column stem could not be found in your metadata.\n"
            "    - Check that the library ID in the metadata matches the MEGAN\n"
            "      column header (after stripping the file-path suffix).\n"
            "    - Fill in sample_id, sample_type, site, and is_negative_control\n"
            "      manually, then set holi_library_name as above.\n"
            "\n"
            "  You can also supply a --holi-map CSV with two columns\n"
            "  (megan_stem, holi_sample_name) to provide overrides without\n"
            "  editing the linker by hand — re-run this script to rebuild.\n"
        )
        print(
            "  Once all issues are resolved, validate with:\n"
            f"    metamerge check --megan-counts <megan_dir_or_file> --holi <metadmg.csv> --linker {out_path}\n"
        )
    else:
        print(
            "  Everything looks good.  You can proceed directly to the merge.\n"
            "\n"
            "  (Optional) Spot-check a few rows in the preview above:\n"
            "    - holi_library_name should match an actual sample name from your\n"
            "      metaDMG CSV's 'sample' column.\n"
            "    - ctrl=NEG for negative controls (blanks), ctrl=POS for positive\n"
            "      controls, ctrl=- for real samples.  Negative controls drive the\n"
            "      Blank-associated classification so incorrect flags here will\n"
            "      affect results.  Positive controls are treated as real samples\n"
            "      for count evidence but are never blank-associated.\n"
            "    - how_matched shows how the Holi name was derived; 'exact-match' is\n"
            "      best, 'stripped-extraction-suffix' and similar are also reliable.\n"
            "    - Positive controls with no Holi match get 'NO_HOLI_DATA'\n"
            "      automatically (not REVIEW_NEEDED) — this is expected.\n"
            "\n"
            "  Step 1 — validate inputs (fast, nothing written):\n"
            f"    metamerge check --megan-counts <megan_dir_or_file> --holi <metadmg.csv> --linker {out_path}\n"
            "\n"
            "  Step 2 — run the full merge:\n"
            f"    metamerge run --megan-counts <megan_dir_or_file> --holi <metadmg.csv> --linker {out_path} --outdir <output_directory/> --online-common-names --render-graphs\n"
            "\n"
            "  `--yes` skips the interactive/manual validation stop and is mainly useful\n"
            "  for automated runs or when you intentionally want to continue despite warnings.\n"
        )


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the MetaMerge linker subcommand parser."""
    parser = argparse.ArgumentParser(
        prog="metamerge linker",
        description=__doc__,
        epilog=(
            "Control detection rules\n"
            "-----------------------\n"
            "Rows are classified into four control types using sample_id + sample_type:\n\n"
            "  negative_control      — sample_type contains 'blank', 'negative',\n"
            "                          'air blank', 'extraction blank', or\n"
            "                          'library blank' (case-insensitive).\n"
            "                          These drive Blank-associated calls.\n\n"
            "  positive_control      — sample_type contains 'positive' or\n"
            "                          'positive control'.  These are NOT blanks\n"
            "                          and will NEVER drive Blank-associated calls.\n"
            "                          Positive controls with no Holi match are\n"
            "                          assigned NO_HOLI_DATA automatically.\n\n"
            "  environmental_control — sample_type contains 'web', 'wash', 'swab',\n"
            "                          'surface', or 'environmental control'.\n\n"
            "  control_unspecified   — sample_type contains 'control' without a\n"
            "                          positive/negative qualifier; treated as\n"
            "                          a negative control conservatively.\n\n"
            "How Holi names are interpreted\n"
            "------------------------------\n"
            "The linker first tries an exact sample_id <-> Holi sample match, then\n"
            "applies conservative normalizations (strip extraction suffix, convert\n"
            "underscores to hyphens, remove leading zeros).  Positive controls with\n"
            "no Holi match get NO_HOLI_DATA (not REVIEW_NEEDED).  All other\n"
            "unresolved rows are flagged REVIEW_NEEDED for manual correction."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version="metamerge linker 0.5.6")

    parser.add_argument(
        "--metadata", required=True,
        help=(
            "Project metadata spreadsheet (.xlsx, .csv, or .tsv).  "
            "Must contain columns for library IDs (--shg-col, --enr-col), "
            "sample_id, sample_type, and optionally site, context."
        ),
    )
    parser.add_argument(
        "--megan-counts", "--megan", dest="megan_counts", nargs="+", required=True,
        help=(
            "One or more MEGAN count-matrix TSV/CSV files, or a directory "
            "containing them.  Column headers must be the full MEGAN library "
            "file names."
        ),
    )
    parser.add_argument(
        "--holi", "--metadmg", dest="holi",
        help=(
            "Holi/metaDMG CSV file (optional but recommended). Used for validation only — "
            "confirms that derived holi_library_name values exist as sample "
            "names in the metaDMG output."
        ),
    )
    parser.add_argument(
        "--holi-map",
        help=(
            "Manual holi-name override CSV with columns 'megan_stem' and "
            "'holi_sample_name'.  Use this for edge cases where the automatic "
            "derivation fails (e.g. positive controls with non-standard naming)."
        ),
    )
    parser.add_argument(
        "--out", default="library_linker.csv",
        help="Output CSV path (default: library_linker.csv).",
    )
    parser.add_argument(
        "--shg-col", default=DEFAULT_SHG_COL,
        help=f"Metadata column for first-chemistry library IDs (default: {DEFAULT_SHG_COL}).",
    )
    parser.add_argument(
        "--enr-col", default=DEFAULT_ENR_COL,
        help=f"Metadata column for second-chemistry library IDs (default: {DEFAULT_ENR_COL}).",
    )
    parser.add_argument(
        "--sample-id-col", default=DEFAULT_SAMPLE_ID,
        help=f"Metadata column for biological sample IDs (default: {DEFAULT_SAMPLE_ID}).",
    )
    parser.add_argument(
        "--site-col", default=DEFAULT_SITE_COL,
        help=f"Metadata column for site/location (default: {DEFAULT_SITE_COL}).",
    )
    parser.add_argument(
        "--type-col", default=DEFAULT_TYPE_COL,
        help=f"Metadata column for sample type (default: {DEFAULT_TYPE_COL}).",
    )
    parser.add_argument(
        "--context-col", default=DEFAULT_CONTEXT,
        help=f"Metadata column for archaeological context (default: {DEFAULT_CONTEXT}).",
    )
    parser.add_argument(
        "--greek-col", default=DEFAULT_GREEK,
        help=f"Metadata column for extraction replicate (default: {DEFAULT_GREEK}).",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress messages (warnings are still shown).",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    """CLI entry point for the MetaMerge linker generator."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    megan_paths = resolve_megan_paths(args.megan_counts)
    if not megan_paths:
        print("Error: no MEGAN TSV files found.", file=sys.stderr)
        sys.exit(1)

    linker = build_linker(
        meta_path=Path(args.metadata),
        megan_paths=megan_paths,
        metadmg_path=Path(args.holi) if args.holi else None,
        holi_map_path=Path(args.holi_map) if args.holi_map else None,
        shg_col=args.shg_col,
        enr_col=args.enr_col,
        sample_id_col=args.sample_id_col,
        site_col=args.site_col,
        type_col=args.type_col,
        context_col=args.context_col,
        greek_col=args.greek_col,
        verbose=not args.quiet,
    )

    out_path = Path(args.out)
    linker.to_csv(out_path, index=False)

    if not args.quiet:
        print_linker_report(linker, out_path)


if __name__ == "__main__":
    main()

"""General utility helpers shared across MetaMerge modules.

This module provides lightweight, dependency-free helpers for:
  - taxon name and rank normalization
  - table reading (xlsx / csv / tsv)
  - column-name detection via alias lists
  - status priority ordering

All normalization is deliberately conservative: we only strip whitespace and
compress internal spaces. Fuzzy or phonetic matching is intentionally excluded
because false positive taxon merges are worse than unmatched taxa.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd


# DNA-support status hierarchy, most supported → least supported.
STATUS_ORDER = [
    "Very high confidence",
    "High confidence",
    "Supported",
    "Tentative",
    "Weak support",
    "Blank-associated",
]

STATUS_PRIORITY = {name: i for i, name in enumerate(STATUS_ORDER)}


def normalize_name(value: object) -> str:
    """Normalize taxon names conservatively for exact matching.

    Operations performed:
      - Convert to string and strip surrounding whitespace
      - Compress internal whitespace runs to a single space
      - Replace curly/backtick apostrophes with a straight apostrophe

    We intentionally avoid stemming, lower-casing, or diacritic removal.
    If taxon names differ more substantially between workflows, users should
    provide a manual common-name override or fix the names upstream.

    Args:
        value: Any value — None, float NaN, int, or string.

    Returns:
        Normalized string, or empty string if the input is None/NaN.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    text = text.replace("\u2019", "'").replace("`", "'")
    text = re.sub(r"\s+", " ", text)
    return text


def normalize_rank(value: object) -> str:
    """Normalize taxonomic rank strings to lowercase with collapsed whitespace.

    Args:
        value: Any value — None, float NaN, or string.

    Returns:
        Lowercase normalized rank string, or empty string if None/NaN.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return re.sub(r"\s+", " ", str(value).strip().lower())


def slugify(text: str) -> str:
    """Create a filesystem-safe ASCII slug from a taxon or library name.

    Args:
        text: Arbitrary string to slugify.

    Returns:
        Lowercase alphanumeric slug with non-alphanumeric characters replaced
        by underscores.  Returns ``"unnamed"`` if the result would be empty.
    """
    text = normalize_name(text).lower()
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text or "unnamed"


def detect_delimiter(path: Path) -> str:
    """Heuristically detect whether a plain-text table uses tabs or commas.

    Reads only the first line of the file, so this is fast even for very large
    tables.

    Args:
        path: Path to the file.

    Returns:
        ``"\\t"`` if tabs outnumber commas on the first line, else ``","``.
    """
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        first = handle.readline()
    return "\t" if first.count("\t") > first.count(",") else ","


def safe_bool(value: object) -> bool:
    """Convert a wide range of text or numeric values to a Python bool.

    Handles the typical spreadsheet representations of True/False that pandas
    may not parse automatically when reading metadata tables.

    Args:
        value: Any value to coerce.

    Returns:
        ``True`` for ``1``, ``"true"``, ``"t"``, ``"yes"``, ``"y"`` (all case-
        insensitive); ``False`` for everything else, including ``None``.
    """
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def select_best_status(series: Iterable[str]) -> Optional[str]:
    """Return the highest-priority DNA-support status from an iterable.

    Args:
        series: Iterable of status strings.

    Returns:
        The string with the lowest STATUS_PRIORITY value (i.e., most supported),
        or ``None`` if no recognized status is present.
    """
    items = [x for x in series if x in STATUS_PRIORITY]
    if not items:
        return None
    return min(items, key=lambda x: STATUS_PRIORITY[x])


def rank_to_level(rank: str, ordered_ranks: list[str]) -> Optional[int]:
    """Map a taxonomic rank name to a numeric level in a lineage hierarchy.

    Args:
        rank: Rank string (will be normalized before lookup).
        ordered_ranks: Ordered list of normalized rank strings, from most
            specific (index 0) to broadest.

    Returns:
        Integer index of the rank, or ``None`` if not found.
    """
    rank = normalize_rank(rank)
    if rank not in ordered_ranks:
        return None
    return ordered_ranks.index(rank)


def read_table(path: Path, sheet_name: Optional[str] = None) -> pd.DataFrame:
    """Read a tabular file (xlsx, csv, or tsv) into a pandas DataFrame.

    For Excel files the first sheet is used unless ``sheet_name`` is provided.
    For CSV/TSV files the delimiter is detected automatically.

    Args:
        path: Path to the file.
        sheet_name: Sheet name for Excel files; ignored for CSV/TSV.

    Returns:
        DataFrame with all rows and columns from the file.

    Raises:
        ValueError: If the file extension is not supported.
    """
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(path, sheet_name=0 if sheet_name is None else sheet_name)
    if suffix in {".csv", ".tsv", ".txt"}:
        sep = detect_delimiter(path)
        return pd.read_csv(path, sep=sep)
    raise ValueError(f"Unsupported input format: {path}")


def find_first_matching_column(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    """Return the first column name from ``candidates`` that exists in ``df``.

    Used by the IO layer to map project-specific column names to the
    standardized internal names via alias lists.

    Args:
        df: DataFrame to search.
        candidates: Ordered list of candidate column names.

    Returns:
        The first matching column name, or ``None`` if none match.
    """
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    return None

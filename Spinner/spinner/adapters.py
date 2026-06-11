"""Adapter / primer contamination scanning.

Adapters are loaded from a TSV with columns: name, sequence, max_mismatch, action.
Scanning uses a sliding-window Hamming-distance search on both the forward and
(optionally) the reverse-complement of each adapter.  Terminal vs. internal
classification uses a configurable endpoint window in bp.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .seq_utils import revcomp


@dataclass
class Adapter:
    name: str
    sequence: str
    max_mismatch: int
    action: str = "reject"


def load_adapters(path: str) -> List[Adapter]:
    out: List[Adapter] = []
    if not path:
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                name = (row.get("name") or "").strip()
                seq = re.sub(r"\s+", "", (row.get("sequence") or "").upper())
                if name and seq:
                    out.append(
                        Adapter(
                            name=name,
                            sequence=seq,
                            max_mismatch=int(row.get("max_mismatch") or 0),
                            action=(row.get("action") or "reject").lower(),
                        )
                    )
    except FileNotFoundError:
        pass
    return out


def hamming_leq(a: str, b: str, max_mm: int) -> bool:
    """Return True if strings *a* and *b* differ in at most *max_mm* positions."""
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > max_mm:
                return False
    return True


def find_adapter(
    seq: str, ad: Adapter, cfg: dict
) -> Optional[Tuple[int, str]]:
    """Search *seq* for adapter *ad*, allowing mismatches.

    Returns ``(position, strand)`` of the first hit (0-based start on *seq*),
    or ``None`` if not found.  Strand is ``"+"`` for forward or ``"-"`` for
    reverse complement.
    """
    acfg = cfg.get("adapter_screen", {})
    min_len = int(acfg.get("min_adapter_match_length", 12))
    scan_rc = acfg.get("scan_reverse_complement", True)

    patterns: List[Tuple[str, str]] = [(ad.sequence, "+")]
    if scan_rc:
        patterns.append((revcomp(ad.sequence), "-"))

    for pat, strand in patterns:
        if len(pat) < min_len or len(seq) < len(pat):
            continue
        # Exact search first (fast path).
        idx = seq.find(pat)
        if idx >= 0:
            return idx, strand
        # Hamming sliding window for mismatches.
        if ad.max_mismatch > 0:
            for i in range(len(seq) - len(pat) + 1):
                if hamming_leq(seq[i : i + len(pat)], pat, ad.max_mismatch):
                    return i, strand
    return None


def is_terminal(pos: int, match_len: int, seq_len: int, win: int) -> bool:
    """Return True if a hit at *pos* with length *match_len* is within *win*
    bases of either end of a sequence of length *seq_len*."""
    return pos <= win or (pos + match_len) >= (seq_len - win)

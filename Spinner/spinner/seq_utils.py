"""Low-level sequence utilities: IUPAC, revcomp, GC, entropy, homopolymer."""
from __future__ import annotations

import hashlib
import math
from collections import Counter

IUPAC_DNA: set = set("ACGTRYSWKMBDHVN")

_COMP_TABLE = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def revcomp(s: str) -> str:
    return s.translate(_COMP_TABLE)[::-1]


def seq_hash(s: str) -> str:
    return hashlib.sha256(s.upper().encode()).hexdigest()


def gc_frac(s: str) -> float:
    bases = [x for x in s.upper() if x in "ACGT"]
    return sum(1 for x in bases if x in "GC") / len(bases) if bases else 0.0


def entropy(s: str) -> float:
    counts = Counter(x for x in s.upper() if x in "ACGT")
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counts.values())


def max_hpoly(s: str) -> int:
    if not s:
        return 0
    best = cur = 1
    last = s[0].upper()
    for b in s[1:].upper():
        if b == last:
            cur += 1
            if cur > best:
                best = cur
        else:
            last = b
            cur = 1
    return best

"""Bad-keyword header scanning.

Keywords are loaded from a TSV with columns: keyword, action, reason.
Action should be ``reject`` or ``review``.  Case-insensitive matching is
the default (configurable via ``bad_keywords.case_sensitive`` in config).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import List


@dataclass
class Keyword:
    keyword: str
    action: str   # "reject" | "review"
    reason: str = ""


def load_keywords(path: str) -> List[Keyword]:
    out: List[Keyword] = []
    if not path:
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                kw = (row.get("keyword") or "").strip()
                if kw:
                    out.append(
                        Keyword(
                            keyword=kw,
                            action=(row.get("action") or "review").lower(),
                            reason=row.get("reason") or "",
                        )
                    )
    except FileNotFoundError:
        pass
    return out

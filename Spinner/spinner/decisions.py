"""Scoring, decision assignment, and decisions TSV writing.

score_decide() is called once before capping (to compute initial decisions)
and once after capping (to incorporate cap_exceeded reasons).
"""
from __future__ import annotations

import csv
import dataclasses
from typing import Dict

from .annotation import Annotation


def score_decide(ann: Dict[str, Annotation], cfg: dict) -> None:
    """Assign decision_score and decision for every record in *ann*.

    Hard-reject reasons override score thresholds.  Review reasons prevent a
    record from being promoted to KEEP even if the score is high enough.
    """
    rules = cfg.get("decision_rules", {})
    hard_reject = set(rules.get("hard_reject_reasons", []))
    review_reasons = set(rules.get("review_reasons", []))
    scoring = cfg.get("scoring", {})
    keep_min = int(rules.get("score_thresholds", {}).get("keep_min", 80))
    rev_min = int(rules.get("score_thresholds", {}).get("review_min", 50))

    for a in ann.values():
        score = int(scoring.get("start", 100))
        for r in a.reasons:
            score += int(scoring.get(r, 0))
        a.decision_score = score
        rs = set(a.reasons)

        if rs & hard_reject:
            a.decision = "REJECT"
        elif score >= keep_min and not (rs & review_reasons):
            a.decision = "KEEP"
        elif score >= rev_min:
            a.decision = "REVIEW"
        else:
            a.decision = "REJECT"


def write_decisions(ann: Dict[str, Annotation], path: str) -> None:
    """Write decisions TSV with one row per input record (including duplicates)."""
    rows = [a.as_dict() for a in ann.values()]
    if not rows:
        return
    fields = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

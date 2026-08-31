"""Scoring, decision assignment, and decisions TSV writing.

score_decide() is called once before capping (to compute initial decisions)
and once after capping (to incorporate cap_exceeded reasons).
"""
from __future__ import annotations

import csv
from typing import Dict

from .annotation import Annotation


def score_decide(ann: Dict[str, Annotation], cfg: dict) -> None:
    """Assign decision_score and decision for every record in *ann*.

    Hard-reject reasons always override score thresholds -- unaffected by
    ``three_state_mode``.

    Default (``decision_rules.three_state_mode: false``, real design change 2026-08-30):
    accept/reject only. Every reason in ``review_reasons`` still applies its own score
    penalty, but no longer blocks KEEP on its own -- a single ``score_thresholds.keep_min``
    cutoff decides. This reflects a real shift in Spinner's role: Fillet's ``bait-eval``
    (a separate, comprehensive, final identity/capture-competition check run against the
    actual finished baits) is now the authoritative last word on contamination, so Spinner
    itself only needs to do cheap, fast triage -- it no longer needs a middle "REVIEW" tier
    forcing manual review before trusting its own calls. See CHANGELOG.md for the full
    real-data motivation.

    Set ``three_state_mode: true`` to restore the old KEEP/REVIEW/REJECT behavior
    (``score_thresholds.review_min`` only has an effect in that mode): hard-reject reasons
    override thresholds; ``review_reasons`` additionally block KEEP even when the score is
    high enough; a mid-range score is REVIEW rather than REJECT.
    """
    rules = cfg.get("decision_rules", {})
    hard_reject = set(rules.get("hard_reject_reasons", []))
    review_reasons = set(rules.get("review_reasons", []))
    three_state_mode = bool(rules.get("three_state_mode", False))
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
        elif not three_state_mode:
            a.decision = "KEEP" if score >= keep_min else "REJECT"
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

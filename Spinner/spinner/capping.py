"""Species × marker capping to limit over-represented taxa.

cap_refs() runs after score_decide() so that it can use existing decisions
and decision_scores for ranking.  Records that exceed the cap get the
cap_exceeded reason added; score_decide() is called again afterward to
incorporate that reason into the final decision.

rescue_sole_representatives() runs after the final score_decide() and promotes
the best available REVIEW record to KEEP for any species that would otherwise
have zero KEEP records — critical for ancient DNA work where the only reference
for a rare taxon may be imperfect.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict

from .annotation import Annotation


def cap_refs(ann: Dict[str, Annotation], cfg: dict) -> None:
    """Cap the number of references per species × marker class combination."""
    if not cfg.get("steps", {}).get("cap_references", True):
        return

    c = cfg.get("capping", {})
    mode = c.get("mode", "species_marker")
    by_class = c.get("max_per_species_marker", {})
    max_total = int(c.get("max_per_species_total", 10**9))
    cap_action = c.get("cap_action", "review").upper()
    uncapped = set(c.get("uncapped_classes", []))

    groups: Dict[tuple, list] = defaultdict(list)
    for a in ann.values():
        # Already hard-rejected records are excluded from counting.
        if a.decision == "REJECT":
            continue
        if a.marker_class in uncapped:
            continue
        species = a.species_guess or a.genus_guess or "Unknown"
        key = (species, a.marker_class) if mode == "species_marker" else (species, "ALL")
        groups[key].append(a)

    for (sp, klass), items in groups.items():
        cap = int(by_class.get(klass, max_total)) if mode == "species_marker" else max_total
        # Sort: prefer KEEP > REVIEW, then higher score, longer, lower N.
        items.sort(
            key=lambda a: (
                a.decision != "KEEP",
                -a.decision_score,
                -a.length,
                a.n_fraction,
                a.accession,
            )
        )
        for rank, a in enumerate(items, 1):
            a.cap_rank = rank
            if rank > cap:
                a.add_reason("cap_exceeded")
                if cap_action == "REJECT":
                    a.decision = "REJECT"
                elif a.decision == "KEEP":
                    a.decision = "REVIEW"


def rescue_sole_representatives(ann: Dict[str, Annotation], cfg: dict) -> int:
    """Promote the best REVIEW record to KEEP for species with zero KEEP records.

    For ancient DNA metagenomics, the only NCBI reference for a rare or extinct
    species may be imperfect (high N, no BLAST match, etc.).  Losing it entirely
    leaves a hole in the mapping database.  This function rescues one record per
    species × marker_class group that would otherwise have zero KEEP records.

    Returns the count of records rescued.

    Activated by setting ``capping.rescue_sole_representatives: true`` in config.
    The rescued record receives reason ``sole_representative`` and its decision
    is set to KEEP regardless of score.
    """
    if not cfg.get("capping", {}).get("rescue_sole_representatives", False):
        return 0

    # Group by (species, marker_class) — same grouping as cap_refs().
    groups: Dict[tuple, list] = defaultdict(list)
    for a in ann.values():
        species = a.species_guess or a.genus_guess or ""
        if not species or species.lower() in ("unknown", ""):
            continue  # no species info — can't meaningfully rescue
        groups[(species, a.marker_class)].append(a)

    rescued = 0
    for (species, klass), items in groups.items():
        keep_count = sum(1 for a in items if a.decision == "KEEP")
        if keep_count > 0:
            continue  # already has at least one KEEP — no rescue needed

        candidates = [a for a in items if a.decision == "REVIEW"]
        if not candidates:
            continue  # all REJECTED — nothing safe to rescue

        # Pick the best candidate: highest score, then longest, then lowest N.
        best = max(
            candidates,
            key=lambda a: (a.decision_score, a.length, -a.n_fraction),
        )
        best.decision = "KEEP"
        best.add_reason("sole_representative")
        rescued += 1

    return rescued

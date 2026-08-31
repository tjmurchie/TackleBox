"""Species × marker capping to limit over-represented taxa.

cap_refs() runs after score_decide() so that it can use existing decisions
and decision_scores for ranking.  Records that exceed the cap get the
cap_exceeded reason added; score_decide() is called again afterward to
incorporate that reason into the final decision.

Note on cap_action: cap_refs() sets a.decision directly, but that is purely
advisory/for-immediate-inspection — the *authoritative* decision is always
whatever the subsequent score_decide() call computes from cfg["decision_rules"]
purely off of a.reasons (score_decide() never reads or preserves a pre-existing
a.decision).  This means cap_action == "review" naturally has real, lasting
effect (it only ever adds the soft "cap_exceeded" review reason).  But
cap_action == "reject" would otherwise be silently thrown away by the next
score_decide() call, because plain "cap_exceeded" is a *review* reason, not a
hard-reject reason, in cfg["decision_rules"].  To make cap_action == "reject"
actually reject, cap_refs() additionally adds the distinct "cap_exceeded_reject"
reason (kept in cfg["decision_rules"]["hard_reject_reasons"]) so the following
score_decide() call forces REJECT with no special-casing needed here.

rescue_sole_representatives() runs after the final score_decide() and promotes
the best available REVIEW record to KEEP for any species that would otherwise
have zero KEEP records — critical for ancient DNA work where the only reference
for a rare taxon may be imperfect.  Because a hard-capped-to-zero record now
ends up with decision == "REJECT" (not "REVIEW") once cap_action == "reject",
it is automatically excluded from rescue_sole_representatives()'s REVIEW-only
candidate list — no separate cap-awareness is needed there.
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
                    # "cap_exceeded" alone is only a *review* reason in the
                    # default decision_rules, so the score_decide() call that
                    # runs right after cap_refs() would otherwise silently
                    # undo this REJECT.  Add a distinct hard-reject reason so
                    # cap_action: "reject" survives that recomputation.
                    a.add_reason("cap_exceeded_reject")
                    a.decision = "REJECT"
                # cap_action == "review" (default): soft-penalty-only, 2026-08-30 --
                # adding the plain "cap_exceeded" reason above is sufficient; the
                # following score_decide() call is the sole authority on the real
                # decision either way (three_state_mode blocks KEEP via review_reasons
                # membership; accept/reject-only mode just applies cap_exceeded's score
                # penalty), so no advisory decision override belongs here at all.


def rescue_sole_representatives(ann: Dict[str, Annotation], cfg: dict) -> int:
    """Promote the best available non-hard-rejected record to KEEP for species with
    zero KEEP records.

    For ancient DNA metagenomics, the only NCBI reference for a rare or extinct
    species may be imperfect (high N, no BLAST match, etc.).  Losing it entirely
    leaves a hole in the mapping database.  This function rescues one record per
    species × marker_class group that would otherwise have zero KEEP records.

    Returns the count of records rescued.

    Activated by setting ``capping.rescue_sole_representatives: true`` in config.
    The rescued record receives reason ``sole_representative`` and its decision
    is set to KEEP regardless of score.

    Candidate pool depends on ``decision_rules.three_state_mode``:

    - **Three-state mode** (legacy): candidates are records already sitting at
      ``REVIEW`` -- unchanged from this function's original behavior. A record
      hard-capped to REJECT (``cap_action: "reject"``, via the distinct
      ``cap_exceeded_reject`` hard-reject reason) is already excluded, since it
      never reaches ``REVIEW`` in the first place.
    - **Accept/reject-only mode** (default since 2026-08-30): there is no REVIEW
      state to rescue FROM, so a species whose sole candidate merely scored below
      ``keep_min`` would otherwise be lost entirely. Candidates here are records at
      ``REJECT`` that were NOT hard-rejected (``decision_rules.hard_reject_reasons``)
      -- i.e. rejected purely on score, not for a structural reason (contamination,
      chimera, duplicate, adapter/vector hit) that a human wouldn't want rescued
      regardless of scarcity.
    """
    if not cfg.get("capping", {}).get("rescue_sole_representatives", False):
        return 0

    three_state_mode = bool(cfg.get("decision_rules", {}).get("three_state_mode", False))
    hard_reject = set(cfg.get("decision_rules", {}).get("hard_reject_reasons", []))

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

        if three_state_mode:
            candidates = [a for a in items if a.decision == "REVIEW"]
        else:
            candidates = [
                a for a in items
                if a.decision == "REJECT" and not (set(a.reasons) & hard_reject)
            ]
        if not candidates:
            continue  # nothing safe to rescue (all hard-rejected, or none at all)

        # Pick the best candidate: highest score, then longest, then lowest N.
        best = max(
            candidates,
            key=lambda a: (a.decision_score, a.length, -a.n_fraction),
        )
        best.decision = "KEEP"
        best.add_reason("sole_representative")
        rescued += 1

    return rescued

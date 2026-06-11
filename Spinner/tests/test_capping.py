"""Tests for species × marker capping (spinner.capping)."""
from __future__ import annotations

import copy

import pytest

from spinner.annotation import Annotation
from spinner.capping import cap_refs, rescue_sole_representatives
from spinner.config import DEFAULT_CONFIG
from spinner.decisions import score_decide


def _make_mito(key: str, species: str = "Rangifer tarandus", score: int = 100) -> Annotation:
    a = Annotation(
        accession=key,
        record_key=key,
        source_file="",
        header=f">{key} {species} mitochondrion",
        length=200,
        seq_sha256=key,
        species_guess=species,
        genus_guess=species.split()[0],
        kingdom="Animal",
        marker_class="Mito",
    )
    a.decision = "KEEP"
    a.decision_score = score
    return a


def _cfg(**overrides):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["steps"]["cap_references"] = True  # capping tests need capping enabled
    cfg["capping"].update(overrides)
    return cfg


def make_ann_dict(*anns):
    return {a.record_key: a for a in anns}


# ---------------------------------------------------------------------------
# Basic capping
# ---------------------------------------------------------------------------

def test_cap_exceeded_above_limit():
    """Records beyond the cap should get cap_exceeded reason."""
    # Cap Mito at 2 for tests.
    cfg = _cfg(max_per_species_marker={"Mito": 2, "Plastid": 20, "NucMark": 10, "Other": 5})
    recs = [_make_mito(f"R{i}") for i in range(5)]
    ann = make_ann_dict(*recs)
    cap_refs(ann, cfg)
    exceeded = [a for a in ann.values() if "cap_exceeded" in a.reasons]
    assert len(exceeded) == 3


def test_cap_not_exceeded_within_limit():
    """Records within the cap should not be flagged."""
    cfg = _cfg(max_per_species_marker={"Mito": 5, "Plastid": 20, "NucMark": 10, "Other": 5})
    recs = [_make_mito(f"R{i}") for i in range(3)]
    ann = make_ann_dict(*recs)
    cap_refs(ann, cfg)
    exceeded = [a for a in ann.values() if "cap_exceeded" in a.reasons]
    assert len(exceeded) == 0


def test_cap_respects_score_order():
    """Higher-scoring records should be kept; lower-scoring ones capped first."""
    cfg = _cfg(max_per_species_marker={"Mito": 1, "Plastid": 20, "NucMark": 10, "Other": 5})
    high = _make_mito("HIGH", score=120)
    low = _make_mito("LOW", score=60)
    ann = make_ann_dict(high, low)
    cap_refs(ann, cfg)
    assert "cap_exceeded" not in ann["HIGH"].reasons
    assert "cap_exceeded" in ann["LOW"].reasons


def test_cap_skips_already_rejected():
    """Records that are already REJECT should not be counted toward the cap."""
    cfg = _cfg(max_per_species_marker={"Mito": 1, "Plastid": 20, "NucMark": 10, "Other": 5})
    rejected = _make_mito("REJECT_REC")
    rejected.decision = "REJECT"
    kept = _make_mito("KEPT_REC")
    ann = make_ann_dict(rejected, kept)
    cap_refs(ann, cfg)
    assert "cap_exceeded" not in ann["KEPT_REC"].reasons


def test_cap_different_species_independent():
    """Records from different species are capped independently."""
    cfg = _cfg(max_per_species_marker={"Mito": 1, "Plastid": 20, "NucMark": 10, "Other": 5})
    sp1 = _make_mito("SP1A", species="Rangifer tarandus")
    sp2 = _make_mito("SP2A", species="Betula papyrifera")
    ann = make_ann_dict(sp1, sp2)
    cap_refs(ann, cfg)
    # Each species has only 1 record, which is within their cap of 1.
    assert "cap_exceeded" not in ann["SP1A"].reasons
    assert "cap_exceeded" not in ann["SP2A"].reasons


def test_cap_action_changes_decision_to_review():
    """cap_action=review should change KEEP -> REVIEW for cap-exceeded records."""
    cfg = _cfg(
        max_per_species_marker={"Mito": 1, "Plastid": 20, "NucMark": 10, "Other": 5},
        cap_action="review",
    )
    recs = [_make_mito(f"R{i}") for i in range(3)]
    ann = make_ann_dict(*recs)
    cap_refs(ann, cfg)
    exceeded = [a for a in ann.values() if "cap_exceeded" in a.reasons]
    for a in exceeded:
        assert a.decision == "REVIEW"


def test_cap_uncapped_classes_not_capped():
    """Classes listed in uncapped_classes should not be capped."""
    cfg = _cfg(
        max_per_species_marker={"Mito": 1, "Plastid": 20, "NucMark": 10, "Other": 5},
        uncapped_classes=["Mito"],
    )
    recs = [_make_mito(f"R{i}") for i in range(5)]
    ann = make_ann_dict(*recs)
    cap_refs(ann, cfg)
    exceeded = [a for a in ann.values() if "cap_exceeded" in a.reasons]
    assert len(exceeded) == 0


# ---------------------------------------------------------------------------
# Sole representative rescue
# ---------------------------------------------------------------------------

def _make_review_mito(key: str, species: str = "Rangifer tarandus", score: int = 65) -> Annotation:
    a = _make_mito(key, species, score)
    a.decision = "REVIEW"
    return a


def _rescue_cfg(**overrides):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["capping"]["rescue_sole_representatives"] = True
    cfg["capping"].update(overrides)
    return cfg


def test_rescue_promotes_best_review_when_no_keep():
    """Species with 0 KEEP and ≥1 REVIEW should have its best record rescued."""
    review1 = _make_review_mito("R1", score=65)
    review2 = _make_review_mito("R2", score=55)
    ann = make_ann_dict(review1, review2)
    cfg = _rescue_cfg()
    rescued = rescue_sole_representatives(ann, cfg)
    assert rescued == 1
    assert ann["R1"].decision == "KEEP"
    assert "sole_representative" in ann["R1"].reasons
    assert ann["R2"].decision == "REVIEW"


def test_rescue_does_not_promote_when_keep_exists():
    """Species already with a KEEP record should not trigger rescue."""
    kept = _make_mito("K1")  # decision=KEEP
    review = _make_review_mito("R1")
    ann = make_ann_dict(kept, review)
    cfg = _rescue_cfg()
    rescued = rescue_sole_representatives(ann, cfg)
    assert rescued == 0
    assert ann["K1"].decision == "KEEP"
    assert ann["R1"].decision == "REVIEW"


def test_rescue_does_not_promote_from_reject():
    """All-REJECT species cannot be rescued (no safe candidate)."""
    rejected = _make_mito("REJ1")
    rejected.decision = "REJECT"
    ann = make_ann_dict(rejected)
    cfg = _rescue_cfg()
    rescued = rescue_sole_representatives(ann, cfg)
    assert rescued == 0
    assert ann["REJ1"].decision == "REJECT"


def test_rescue_disabled_when_explicitly_false():
    """rescue_sole_representatives should be a no-op when config flag is explicitly False."""
    review = _make_review_mito("R1")
    ann = make_ann_dict(review)
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["capping"]["rescue_sole_representatives"] = False  # explicitly disable
    rescued = rescue_sole_representatives(ann, cfg)
    assert rescued == 0
    assert ann["R1"].decision == "REVIEW"


def test_rescue_independent_per_species():
    """Each species is rescued independently; one with KEEP should not prevent rescue of another."""
    sp1_keep = _make_mito("SP1_K", species="Rangifer tarandus")
    sp2_review = _make_review_mito("SP2_R", species="Betula papyrifera")
    ann = make_ann_dict(sp1_keep, sp2_review)
    cfg = _rescue_cfg()
    rescued = rescue_sole_representatives(ann, cfg)
    assert rescued == 1
    assert ann["SP2_R"].decision == "KEEP"
    assert ann["SP1_K"].decision == "KEEP"

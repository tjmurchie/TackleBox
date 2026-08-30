"""Tests for scoring and decision assignment (spinner.decisions)."""
from __future__ import annotations

import copy


from spinner.annotation import Annotation, annotate
from spinner.config import DEFAULT_CONFIG
from spinner.decisions import score_decide
from spinner.fasta import parse_fasta


def make_ann(**kwargs) -> Annotation:
    defaults = dict(
        accession="TEST.1",
        record_key="TEST.1",
        source_file="",
        header=">TEST.1 test",
        length=200,
        seq_sha256="abc",
    )
    defaults.update(kwargs)
    return Annotation(**defaults)


def make_ann_dict(*anns):
    return {a.record_key: a for a in anns}


def _cfg():
    return copy.deepcopy(DEFAULT_CONFIG)


# ---------------------------------------------------------------------------
# Hard-reject reasons -> REJECT
# ---------------------------------------------------------------------------

def test_hard_reject_adapter_internal():
    a = make_ann()
    a.add_reason("adapter_internal")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "REJECT"


def test_hard_reject_duplicate_accession():
    a = make_ann()
    a.add_reason("duplicate_accession")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "REJECT"


def test_n_fraction_high_is_soft_penalty():
    """n_fraction_high is NOT a hard reject — it is a score penalty only.
    A high-N record from a rare species should be kept, not hard-rejected.
    Score: 100 + (-20) = 80 >= keep_min(65) and no review_reason -> KEEP.
    """
    a = make_ann()
    a.add_reason("n_fraction_high")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision_score == 80
    assert ann["TEST.1"].decision == "KEEP"


def test_hard_reject_bad_keyword_reject():
    a = make_ann()
    a.add_reason("bad_keyword_reject")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "REJECT"


def test_hard_reject_length_below_min():
    a = make_ann()
    a.add_reason("length_below_min")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "REJECT"


# ---------------------------------------------------------------------------
# Review reasons -> REVIEW (not KEEP, even if score is high)
# ---------------------------------------------------------------------------

def test_review_from_adapter_terminal():
    a = make_ann()
    a.add_reason("adapter_terminal")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    # Score starts at 100, adapter_terminal = -40 -> 60; review reason present.
    assert ann["TEST.1"].decision == "REVIEW"


def test_review_from_bad_keyword_review():
    a = make_ann()
    a.add_reason("bad_keyword_review")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "REVIEW"


def test_review_from_cap_exceeded():
    a = make_ann()
    a.add_reason("cap_exceeded")
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "REVIEW"


# ---------------------------------------------------------------------------
# Clean record -> KEEP
# ---------------------------------------------------------------------------

def test_clean_record_keep():
    a = make_ann()  # no reasons
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision == "KEEP"
    assert ann["TEST.1"].decision_score == 100


# ---------------------------------------------------------------------------
# Score arithmetic
# ---------------------------------------------------------------------------

def test_score_decreases_on_penalty():
    a = make_ann()
    a.add_reason("adapter_terminal")  # -40
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision_score == 60


def test_score_increases_on_bonus():
    a = make_ann()
    a.add_reason("taxonomy_same_species")  # +20
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["TEST.1"].decision_score == 120


def test_score_below_review_min_becomes_reject():
    """If score drops below review_min (30) without any hard-reject or review reason,
    the record should be REJECT.  Clear review_reasons so the decision falls through
    to the score path: 100 - 40 - 30 - 20 = 10 < review_min(30) -> REJECT.
    """
    cfg = _cfg()
    cfg["decision_rules"]["review_reasons"] = []  # test pure score logic
    a = make_ann()
    a.add_reason("adapter_terminal")    # -40
    a.add_reason("bad_keyword_review")  # -30
    a.add_reason("n_fraction_high")     # -20 (soft penalty, not hard reject)
    # Score = 100 - 40 - 30 - 20 = 10, which is below review_min=30.
    ann = make_ann_dict(a)
    score_decide(ann, cfg)
    assert ann["TEST.1"].decision_score == 10
    assert ann["TEST.1"].decision == "REJECT"


# ---------------------------------------------------------------------------
# End-to-end annotation + score_decide integration
# ---------------------------------------------------------------------------

def test_annotation_to_decision_clean(tmp_path):
    p = tmp_path / "t.fasta"
    p.write_text(">CLEAN.1 Rangifer tarandus mitochondrion complete genome\n" + "ACGT" * 30 + "\n")
    recs = parse_fasta([str(p)])
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["steps"]["adapter_screen"] = False
    cfg["steps"]["bad_keyword_screen"] = False
    ann = annotate(recs, cfg)
    score_decide(ann, cfg)
    assert ann["CLEAN.1"].decision == "KEEP"


def test_annotation_to_decision_all_n(tmp_path):
    """An all-N record is not hard-rejected: n_fraction_high is a score penalty.
    The 80-N sequence also triggers homopolymer_long (run of 80 > max 60), which
    is a review_reason -> REVIEW.  The record is preserved for manual inspection,
    not silently discarded.
    """
    p = tmp_path / "t.fasta"
    p.write_text(">ALLN.1 Salix arctica ITS1\n" + "N" * 80 + "\n")
    recs = parse_fasta([str(p)])
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["steps"]["adapter_screen"] = False
    cfg["steps"]["bad_keyword_screen"] = False
    ann = annotate(recs, cfg)
    score_decide(ann, cfg)
    assert "n_fraction_high" in ann["ALLN.1"].reasons
    assert ann["ALLN.1"].decision == "REVIEW"


# ---------------------------------------------------------------------------
# taxonomy_exempt_length fix (2026-08-29): a clean, length-exempted complete
# genome must reach KEEP on its own merit, unlike the "taxonomy_not_checked"
# reason it replaces for this case.
# ---------------------------------------------------------------------------

def test_length_exempt_complete_genome_reaches_keep():
    """Pinned real scenario: EU153454.1 (16,480 bp complete mitogenome, real cluster
    centroid in a live PalaeoSCOPE run) had reasons
    `complete_organelle;taxonomy_not_checked;cluster_representative` and
    decision_score=110 (well above keep_min=65) but stayed stuck at REVIEW purely
    because `taxonomy_not_checked` is in decision_rules.review_reasons. With the fix,
    the exact same score-contributing reasons, but tagged `taxonomy_exempt_length`
    instead of `taxonomy_not_checked` (since this record was excluded from the search
    by max_query_length, not genuinely searched-and-not-found), must reach KEEP."""
    a = make_ann(accession="EU153454.1", record_key="EU153454.1", length=16480)
    a.add_reason("complete_organelle")       # +5
    a.add_reason("taxonomy_exempt_length")   # 0
    a.add_reason("cluster_representative")   # +5
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["EU153454.1"].decision_score == 110
    assert ann["EU153454.1"].decision == "KEEP"


def test_same_score_with_old_taxonomy_not_checked_reason_still_forces_review():
    """Contrast case, same file: a record that was genuinely SUBMITTED to the taxonomy
    search and got no hit (real "taxonomy_not_checked", not a length exemption) must
    still be correctly blocked from auto-KEEP at the identical score -- the fix must not
    weaken review-forcing for genuinely-unverified records, only for confidently-exempt
    ones."""
    a = make_ann(accession="SHORT.1", record_key="SHORT.1", length=500)
    a.add_reason("complete_organelle")     # +5
    a.add_reason("taxonomy_not_checked")   # 0
    a.add_reason("cluster_representative")  # +5
    ann = make_ann_dict(a)
    score_decide(ann, _cfg())
    assert ann["SHORT.1"].decision_score == 110
    assert ann["SHORT.1"].decision == "REVIEW"

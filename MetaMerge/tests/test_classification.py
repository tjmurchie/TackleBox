"""Tests for the MetaMerge classification logic (classify_status).

These tests verify that the DNA-support category assignment in
``metamerge.classify.classify_status`` works correctly for key combinations of
inputs, including the priority ordering (blank-associated always wins) and the
boundary between Very high confidence and High confidence.
"""

from metamerge.classify import classify_status


THRESHOLDS = {"tentative_min_reads": 5, "weak_support_min_reads": 1}


def test_very_high_confidence():
    """Exact damage + ge100 reads + strong counts + clean QC → Very high confidence."""
    status, basis = classify_status(
        exact_damage_support=True,
        exact_damage_support_ge100=True,
        strong_count_support=True,
        some_count_support=True,
        lineage_support=False,
        blank_associated=False,
        qc_label="clean",
        max_real_count=200,
        thresholds=THRESHOLDS,
    )
    assert status == "Very high confidence"
    assert "damage-supported" in basis.lower() or "exact" in basis.lower()


def test_very_high_confidence_blocked_by_strong_caution():
    """Strong caution QC prevents Very high confidence even with ge100 reads."""
    status, _ = classify_status(
        exact_damage_support=True,
        exact_damage_support_ge100=True,
        strong_count_support=True,
        some_count_support=True,
        lineage_support=False,
        blank_associated=False,
        qc_label="strong caution",
        max_real_count=200,
        thresholds=THRESHOLDS,
    )
    # Should fall back to High confidence rather than Very high confidence.
    assert status == "High confidence"


def test_high_confidence():
    """Exact damage + strong counts, but without ge100 reads → High confidence."""
    status, basis = classify_status(
        exact_damage_support=True,
        exact_damage_support_ge100=False,
        strong_count_support=True,
        some_count_support=True,
        lineage_support=False,
        blank_associated=False,
        qc_label="clean",
        max_real_count=100,
        thresholds=THRESHOLDS,
    )
    assert status == "High confidence"


def test_blank_associated_wins():
    """Blank-associated always takes priority regardless of damage/count support."""
    status, _ = classify_status(
        exact_damage_support=True,
        exact_damage_support_ge100=True,
        strong_count_support=True,
        some_count_support=True,
        lineage_support=True,
        blank_associated=True,
        qc_label="clean",
        max_real_count=200,
        thresholds=THRESHOLDS,
    )
    assert status == "Blank-associated"


def test_tentative():
    """No damage or count support but above tentative floor → Tentative."""
    status, _ = classify_status(
        exact_damage_support=False,
        exact_damage_support_ge100=False,
        strong_count_support=False,
        some_count_support=False,
        lineage_support=False,
        blank_associated=False,
        qc_label="not-applicable",
        max_real_count=10,
        thresholds=THRESHOLDS,
    )
    assert status == "Tentative"


def test_weak_support():
    """Below tentative floor but above 0 → Weak support."""
    status, _ = classify_status(
        exact_damage_support=False,
        exact_damage_support_ge100=False,
        strong_count_support=False,
        some_count_support=False,
        lineage_support=False,
        blank_associated=False,
        qc_label="not-applicable",
        max_real_count=2,
        thresholds=THRESHOLDS,
    )
    assert status == "Weak support"

"""Tests for the generalized, source-agnostic classify_status_v2()/
compute_ensemble_support_score() in metamerge.evidence.

The most important test class here is TestBackwardCompatParity: it proves
classify_status_v2() reduces to byte-for-byte identical (status, basis) output
to the original classify_status() for every MEGAN+Holi input combination
test_classification.py already covers -- the load-bearing guarantee that
existing 2-source MEGAN+Holi users see zero behavior change from this
integration.
"""

from metamerge.classify import classify_status
from metamerge.evidence import (
    SourceSignals,
    TaxonEvidence,
    classify_status_v2,
    compute_ensemble_support_score,
)

THRESHOLDS = {"tentative_min_reads": 5, "weak_support_min_reads": 1}
CONFIG = {"thresholds": THRESHOLDS, "ensemble_score": {}}


def _megan_holi_evidence(
    exact_damage_support=False,
    exact_damage_support_ge100=False,
    strong_count_support=False,
    some_count_support=False,
    lineage_support=False,
    blank_associated=False,
    qc_label="clean",
    max_real_count=0,
) -> TaxonEvidence:
    return TaxonEvidence(
        megan=SourceSignals(
            present=True,
            strong_count_support=strong_count_support,
            some_count_support=some_count_support,
            blank_associated=blank_associated,
            max_count=max_real_count,
        ),
        holi=SourceSignals(
            present=True,
            exact_damage_support=exact_damage_support,
            exact_damage_support_ge100=exact_damage_support_ge100,
            qc_label=qc_label,
        ),
        lineage_support=lineage_support,
        max_real_count=max_real_count,
    )


class TestBackwardCompatParity:
    """For every (exact_damage_support, exact_damage_support_ge100,
    strong_count_support, some_count_support, lineage_support, blank_associated,
    qc_label, max_real_count) combination exercised by test_classification.py,
    classify_status_v2() must produce the identical status and basis as calling
    classify_status() directly."""

    CASES = [
        dict(exact_damage_support=True, exact_damage_support_ge100=True,
             strong_count_support=True, some_count_support=True,
             qc_label="clean", max_real_count=200),
        dict(exact_damage_support=True, exact_damage_support_ge100=True,
             strong_count_support=True, some_count_support=True,
             qc_label="strong caution", max_real_count=200),
        dict(exact_damage_support=True, exact_damage_support_ge100=False,
             strong_count_support=True, some_count_support=True,
             qc_label="clean", max_real_count=100),
        dict(exact_damage_support=True, exact_damage_support_ge100=False,
             strong_count_support=False, some_count_support=True,
             qc_label="clean", max_real_count=10),
        dict(exact_damage_support=False, exact_damage_support_ge100=False,
             strong_count_support=True, some_count_support=True,
             qc_label="not-applicable", max_real_count=60),
        dict(exact_damage_support=False, exact_damage_support_ge100=False,
             strong_count_support=False, some_count_support=True,
             lineage_support=True, qc_label="not-applicable", max_real_count=3),
        dict(exact_damage_support=False, exact_damage_support_ge100=False,
             strong_count_support=False, some_count_support=False,
             qc_label="not-applicable", max_real_count=6),
        dict(exact_damage_support=False, exact_damage_support_ge100=False,
             strong_count_support=False, some_count_support=False,
             qc_label="not-applicable", max_real_count=1),
        dict(exact_damage_support=False, exact_damage_support_ge100=False,
             strong_count_support=False, some_count_support=False,
             qc_label="not-applicable", max_real_count=0),
        dict(blank_associated=True, exact_damage_support=True,
             exact_damage_support_ge100=True,
             strong_count_support=True, some_count_support=True,
             qc_label="clean", max_real_count=500),
    ]

    def test_parity_across_all_cases(self):
        for raw_case in self.CASES:
            case = dict(raw_case)  # don't mutate the shared class-level CASES list
            lineage_support = case.pop("lineage_support", False)
            blank_associated = case.pop("blank_associated", False)

            old_status, old_basis = classify_status(
                exact_damage_support=case["exact_damage_support"],
                exact_damage_support_ge100=case["exact_damage_support_ge100"],
                strong_count_support=case["strong_count_support"],
                some_count_support=case["some_count_support"],
                lineage_support=lineage_support,
                blank_associated=blank_associated,
                qc_label=case["qc_label"],
                max_real_count=case["max_real_count"],
                thresholds=THRESHOLDS,
            )

            ev = _megan_holi_evidence(
                lineage_support=lineage_support, blank_associated=blank_associated, **case,
            )
            new_status, new_basis, score, agreement_count, agreement_fraction = classify_status_v2(ev, CONFIG)

            assert new_status == old_status, f"status mismatch for {case}"
            assert new_basis == old_basis, f"basis mismatch for {case}"
            assert 0.0 <= score <= 1.0
            assert 0 <= agreement_count <= 2
            assert 0.0 <= agreement_fraction <= 1.0


class TestHoliFilletOnly:
    """Holi+Fillet, no MEGAN -- Tyler's expected primary real-world case."""

    def test_ge100_damage_plus_authenticated_eco_support_is_very_high(self):
        ev = TaxonEvidence(
            holi=SourceSignals(present=True, exact_damage_support=True, exact_damage_support_ge100=True),
            fillet=SourceSignals(present=True, authenticated=True, eco_support=True),
            max_real_count=0,
        )
        status, basis, score, agreement_count, agreement_fraction = classify_status_v2(ev, CONFIG)
        assert status == "Very high confidence"
        assert "eco" in basis.lower()
        assert score > 0.5
        assert agreement_count == 2  # holi + fillet, matches sources_present here

    def test_damage_supported_no_fillet_authentication_is_supported(self):
        ev = TaxonEvidence(
            holi=SourceSignals(present=True, exact_damage_support=True, exact_damage_support_ge100=False),
            fillet=SourceSignals(present=True, authenticated=False),
        )
        status, basis, score, _, _ = classify_status_v2(ev, CONFIG)
        assert status == "Supported"

    def test_no_support_falls_to_weak(self):
        ev = TaxonEvidence(
            holi=SourceSignals(present=True, exact_damage_support=False),
            fillet=SourceSignals(present=True, authenticated=False),
            max_real_count=0,
        )
        status, basis, score, agreement_count, agreement_fraction = classify_status_v2(ev, CONFIG)
        assert status == "Weak support"
        assert score == 0.0
        assert agreement_count == 0

    def test_tentative_from_read_count_alone(self):
        ev = TaxonEvidence(
            holi=SourceSignals(present=True, exact_damage_support=False),
            fillet=SourceSignals(present=True, authenticated=False),
            max_real_count=6,
        )
        status, _, _, _, _ = classify_status_v2(ev, CONFIG)
        assert status == "Tentative"


class TestMeganFilletOnly:
    def test_fillet_damage_plus_strong_counts_is_high_confidence(self):
        ev = TaxonEvidence(
            megan=SourceSignals(present=True, strong_count_support=True, some_count_support=True),
            fillet=SourceSignals(present=True, exact_damage_support=True, authenticated=True),
            max_real_count=80,
        )
        status, _, _, _, _ = classify_status_v2(ev, CONFIG)
        assert status == "High confidence"


class TestSingleSourceOnly:
    def test_fillet_only_authenticated_caps_at_supported(self):
        """Single-source cases never reach High/Very high confidence -- those
        tiers require convergent evidence from more than one source."""
        ev = TaxonEvidence(
            fillet=SourceSignals(present=True, authenticated=True, eco_support=True, fos_support=True),
        )
        status, basis, score, agreement_count, agreement_fraction = classify_status_v2(ev, CONFIG)
        assert status == "Supported"
        assert "fillet only" in basis.lower()
        assert agreement_count == 1
        assert agreement_fraction == 1.0  # 1 of 1 present sources agrees

    def test_holi_only_exact_damage_caps_at_supported(self):
        ev = TaxonEvidence(holi=SourceSignals(present=True, exact_damage_support=True))
        status, _, _, _, _ = classify_status_v2(ev, CONFIG)
        assert status == "Supported"

    def test_megan_only_strong_counts_caps_at_supported(self):
        ev = TaxonEvidence(megan=SourceSignals(present=True, strong_count_support=True))
        status, _, _, _, _ = classify_status_v2(ev, CONFIG)
        assert status == "Supported"

    def test_no_sources_present_is_weak_support(self):
        ev = TaxonEvidence()
        status, basis, score, agreement_count, agreement_fraction = classify_status_v2(ev, CONFIG)
        assert status == "Weak support"
        assert score == 0.0
        assert agreement_count == 0
        assert agreement_fraction == 0.0


class TestThreeSourceCorroboration:
    def test_all_three_agreeing_upgrades_to_3source_corroborated(self):
        ev = TaxonEvidence(
            megan=SourceSignals(present=True, strong_count_support=True, some_count_support=True),
            holi=SourceSignals(present=True, exact_damage_support=True, exact_damage_support_ge100=True, qc_label="clean"),
            fillet=SourceSignals(present=True, authenticated=True, eco_support=True),
            max_real_count=200,
            discordant=False,
        )
        status, basis, score, agreement_count, agreement_fraction = classify_status_v2(ev, CONFIG)
        assert status == "Very high confidence (3-source corroborated)"
        assert agreement_count == 3
        assert agreement_fraction == 1.0
        assert score > compute_ensemble_support_score(
            TaxonEvidence(
                megan=SourceSignals(present=True, strong_count_support=True),
                holi=SourceSignals(present=True, exact_damage_support_ge100=True),
            )
        )

    def test_discordant_does_not_upgrade_even_if_otherwise_qualifying(self):
        ev = TaxonEvidence(
            megan=SourceSignals(present=True, strong_count_support=True, some_count_support=True),
            holi=SourceSignals(present=True, exact_damage_support=True, exact_damage_support_ge100=True, qc_label="clean"),
            fillet=SourceSignals(present=True, authenticated=True, eco_support=True),
            max_real_count=200,
            discordant=True,
        )
        status, basis, score, _, _ = classify_status_v2(ev, CONFIG)
        assert status != "Very high confidence (3-source corroborated)"
        assert "discordant" in basis.lower() or "disagree" in basis.lower()


class TestEnsembleSupportScore:
    def test_all_three_sources_supporting_scores_highest(self):
        ev = TaxonEvidence(
            megan=SourceSignals(present=True, strong_count_support=True),
            holi=SourceSignals(present=True, exact_damage_support_ge100=True),
            fillet=SourceSignals(present=True, authenticated=True, eco_support=True, pal_support=True, fos_support=True),
        )
        score = compute_ensemble_support_score(ev)
        assert score == 1.0  # clamped; raw sum exceeds 1.0 with all bonuses

    def test_single_source_scores_lower_than_two_sources(self):
        single = compute_ensemble_support_score(
            TaxonEvidence(fillet=SourceSignals(present=True, authenticated=True))
        )
        double = compute_ensemble_support_score(
            TaxonEvidence(
                holi=SourceSignals(present=True, exact_damage_support_ge100=True),
                fillet=SourceSignals(present=True, authenticated=True),
            )
        )
        assert single < double

    def test_discordance_reduces_score(self):
        agreeing = compute_ensemble_support_score(
            TaxonEvidence(
                holi=SourceSignals(present=True, exact_damage_support_ge100=True),
                fillet=SourceSignals(present=True, authenticated=True),
                discordant=False,
            )
        )
        discordant = compute_ensemble_support_score(
            TaxonEvidence(
                holi=SourceSignals(present=True, exact_damage_support_ge100=True),
                fillet=SourceSignals(present=True, authenticated=True),
                discordant=True,
            )
        )
        assert discordant < agreeing

    def test_no_evidence_scores_zero(self):
        assert compute_ensemble_support_score(TaxonEvidence()) == 0.0

    def test_score_never_negative(self):
        ev = TaxonEvidence(discordant=True)
        assert compute_ensemble_support_score(ev) == 0.0

    def test_custom_weights_respected(self):
        ev = TaxonEvidence(megan=SourceSignals(present=True, strong_count_support=True))
        default_score = compute_ensemble_support_score(ev)
        custom_score = compute_ensemble_support_score(ev, {"megan_weight": 0.9})
        assert custom_score != default_score
        assert custom_score == 0.9

"""Tests for fillet_matching.py -- the Fillet-specific analogue of holi.py's
matching/index-building functions."""

import pandas as pd
import pytest

from metamerge.config import load_config
from metamerge.fillet_matching import (
    build_fillet_taxonomy_lookup,
    choose_best_fillet_row,
    make_fillet_exact_index,
    row_has_fillet_strong_count_support,
    row_is_fillet_authenticated,
)


@pytest.fixture
def thresholds():
    return load_config(None)["thresholds"]


def _row(sample="lib1", tax_id_str="9696", tax_name="Puma", tax_rank="genus",
         composite=0.5, tier=2, direct_hard_reads=100):
    return {
        "sample": sample, "tax_id_str": tax_id_str, "tax_name": tax_name, "tax_rank": tax_rank,
        "composite_authenticity": composite, "authenticity_tier": tier,
        "direct_hard_reads": direct_hard_reads,
    }


class TestRowIsFilletAuthenticated:
    def test_high_composite_and_low_tier_authenticated(self, thresholds):
        row = _row(composite=0.8, tier=1)
        assert row_is_fillet_authenticated(row, thresholds) is True

    def test_low_composite_not_authenticated(self, thresholds):
        row = _row(composite=0.05, tier=1)
        assert row_is_fillet_authenticated(row, thresholds) is False

    def test_high_tier_number_not_authenticated(self, thresholds):
        """tier=5 is the LOWEST confidence tier -- must fail even with high composite."""
        row = _row(composite=0.9, tier=5)
        assert row_is_fillet_authenticated(row, thresholds) is False

    def test_tier_zero_rejected_not_authenticated(self, thresholds):
        row = _row(composite=0.9, tier=0)
        assert row_is_fillet_authenticated(row, thresholds) is False

    def test_none_row_not_authenticated(self, thresholds):
        assert row_is_fillet_authenticated(None, thresholds) is False

    def test_nan_composite_not_authenticated(self, thresholds):
        row = _row(composite=float("nan"), tier=1)
        assert row_is_fillet_authenticated(row, thresholds) is False


class TestRowHasStrongCountSupport:
    def test_above_threshold_true(self, thresholds):
        row = _row(direct_hard_reads=100)
        assert row_has_fillet_strong_count_support(row, thresholds) is True

    def test_below_threshold_false(self, thresholds):
        row = _row(direct_hard_reads=2)
        assert row_has_fillet_strong_count_support(row, thresholds) is False

    def test_none_row_false(self, thresholds):
        assert row_has_fillet_strong_count_support(None, thresholds) is False


class TestChooseBestFilletRow:
    def test_empty_list_returns_none(self, thresholds):
        assert choose_best_fillet_row([], thresholds) is None

    def test_authenticated_row_beats_unauthenticated(self, thresholds):
        weak = _row(composite=0.05, tier=5, direct_hard_reads=1000)
        strong = _row(composite=0.8, tier=1, direct_hard_reads=10)
        assert choose_best_fillet_row([weak, strong], thresholds) is strong

    def test_higher_composite_wins_among_authenticated(self, thresholds):
        a = _row(composite=0.5, tier=2)
        b = _row(composite=0.9, tier=1)
        assert choose_best_fillet_row([a, b], thresholds) is b


class TestMakeFilletExactIndex:
    def test_index_by_id_and_name_rank(self):
        df = pd.DataFrame([_row(), _row(sample="lib2", tax_id_str="9903", tax_name="Bison bison", tax_rank="species")])
        by_id, by_name_rank = make_fillet_exact_index(df)
        assert ("lib1", "9696") in by_id
        assert ("lib2", "Bison bison", "species") in by_name_rank
        assert ("lib2", "Bison bison", "") in by_name_rank


class TestBuildFilletTaxonomyLookup:
    def test_lookup_by_id_and_name(self):
        df = pd.DataFrame([_row()])
        lookup = build_fillet_taxonomy_lookup(df)
        assert lookup["9696"]["tax_name"] == "Puma"
        assert lookup["Puma"]["tax_rank"] == "genus"
        assert lookup["9696"]["tax_path"] == ""

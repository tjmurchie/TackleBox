"""Tests for build_merge()'s no-MEGAN path (_build_merge_no_megan in
classify.py) -- Holi+Fillet (Tyler's expected primary real-world case), plus
Holi-only and Fillet-only single-source runs. Mirrors test_smoke.py's shape
but exercises megan_df=None instead of the MEGAN-anchored path."""

import pandas as pd
import pytest

from metamerge.classify import build_merge
from metamerge.config import load_config


def _metadata(rows):
    return pd.DataFrame(rows)


def _base_metadata():
    """2 real libraries (enough for strong_count_min_libraries=2) + 1 blank,
    with holi_library_name AND fillet_library_name crosswalk columns but no
    megan_library_name column at all."""
    return _metadata({
        "merged_library_name": ["Real1", "Real2", "Blank1"],
        "holi_library_name":   ["Real1_holi", "Real2_holi", "Blank1_holi"],
        "fillet_library_name": ["Real1_fillet", "Real2_fillet", ""],
        "is_negative_control": [False, False, True],
        "sample_type":         ["sediment", "sediment", "extraction blank"],
    })


def _holi_df(n_reads=(150.0, 110.0, 0.0), damage=(0.12, 0.10, 0.0), sig=(3.5, 3.0, 0.0)):
    return pd.DataFrame({
        "sample":       ["Real1_holi", "Real2_holi", "Blank1_holi"],
        "tax_id":       [9696, 9696, 9696],
        "tax_id_str":   ["9696", "9696", "9696"],
        "tax_name":     ["Puma", "Puma", "Puma"],
        "tax_rank":     ["genus", "genus", "genus"],
        "N_reads":      list(n_reads),
        "N_alignments": [x * 2 for x in n_reads],
        "damage":       list(damage),
        "significance": list(sig),
        "rho_Ac":       [0.05, 0.04, 0.0],
        "MAP_valid":    [True, True, True],
        "tax_path":     ["Animalia;Chordata;Mammalia;Carnivora;Felidae;Puma"] * 3,
    })


def _fillet_df(composite=(0.8, 0.75), tier=(1, 1), reads=(200.0, 180.0)):
    return pd.DataFrame({
        "sample":                 ["Real1_fillet", "Real2_fillet"],
        "tax_id_str":             ["9696", "9696"],
        "tax_name":               ["Puma", "Puma"],
        "tax_rank":               ["genus", "genus"],
        "composite_authenticity": list(composite),
        "authenticity_tier":      list(tier),
        "direct_hard_reads":      list(reads),
        "eco_support":            [True, True],
        "pal_support":            [False, False],
        "fos_support":            [True, True],
    })


@pytest.fixture
def config():
    return load_config(None)


class TestHoliFilletNoMegan:
    def test_reaches_high_confidence(self, config):
        merged, summary = build_merge(_base_metadata(), None, _holi_df(), config, fillet_df=_fillet_df())
        assert len(merged) == 1
        row = merged.iloc[0]
        assert row["scientific_name"] == "Puma"
        assert row["aDNA_support_status"] in {"Very high confidence", "High confidence"}
        assert row["ensemble_support_score"] > 0.5
        assert bool(row["fillet_authenticated"]) is True
        assert summary["n_taxa"] == 1

    def test_count_columns_populated_from_fallback_proxy(self, config):
        merged, _ = build_merge(_base_metadata(), None, _holi_df(), config, fillet_df=_fillet_df())
        row = merged.iloc[0]
        assert row["count__Real1"] == pytest.approx(150.0)
        assert row["count__Real2"] == pytest.approx(110.0)
        assert row["megan_max_count"] == pytest.approx(150.0)

    def test_blank_associated_without_megan(self, config):
        holi = _holi_df(n_reads=(0.0, 0.0, 120.0), damage=(0.0, 0.0, 0.15), sig=(0.0, 0.0, 3.2))
        merged, _ = build_merge(_base_metadata(), None, holi, config, fillet_df=None)
        assert merged.iloc[0]["aDNA_support_status"] == "Blank-associated"

    def test_taxon_present_only_in_fillet_still_appears(self, config):
        """A taxon Holi never saw at all (e.g. Holi's alignment missed it) but
        Fillet independently authenticated should still surface -- this is a
        real capability the old MEGAN-anchored loop never had (it could only
        ever iterate MEGAN's own taxon list)."""
        holi = _holi_df()  # only has Puma (9696)
        fillet = pd.concat([
            _fillet_df(),
            pd.DataFrame({
                # Deliberately weak on every OTHER signal (low read count,
                # no eco/pal/fos lines) so this exercises the single-source
                # ceiling cleanly -- Fillet's ability to reach higher tiers
                # on its own strength (composite+tier+eco/pal/fos) when
                # Holi has nothing for this taxon is covered separately in
                # test_evidence.py's TestHoliFilletOnly.
                "sample": ["Real1_fillet"], "tax_id_str": ["9999"], "tax_name": ["Urocitellus parryii"],
                "tax_rank": ["species"], "composite_authenticity": [0.6], "authenticity_tier": [2],
                "direct_hard_reads": [20.0], "eco_support": [False], "pal_support": [False], "fos_support": [False],
            }),
        ], ignore_index=True)
        merged, summary = build_merge(_base_metadata(), None, holi, config, fillet_df=fillet)
        assert summary["n_taxa"] == 2
        assert "Urocitellus parryii" in set(merged["scientific_name"])
        squirrel = merged.loc[merged["scientific_name"] == "Urocitellus parryii"].iloc[0]
        # Single-source (Fillet only, for this taxon -- Holi has no rows for it,
        # and Fillet's own corroboration is weak) is capped at "Supported".
        assert squirrel["aDNA_support_status"] == "Supported"


class TestHoliOnlyNoMeganNoFillet:
    def test_single_source_capped_at_supported(self, config):
        merged, summary = build_merge(_base_metadata(), None, _holi_df(), config, fillet_df=None)
        assert summary["n_taxa"] == 1
        row = merged.iloc[0]
        assert row["aDNA_support_status"] == "Supported"
        assert pd.isna(row["fillet_authenticated"])


class TestFilletOnlyNoMeganNoHoli:
    def test_single_source_capped_at_supported(self, config):
        merged, summary = build_merge(_base_metadata(), None, None, config, fillet_df=_fillet_df())
        assert summary["n_taxa"] == 1
        row = merged.iloc[0]
        assert row["aDNA_support_status"] == "Supported"
        assert bool(row["fillet_authenticated"]) is True
        # No Holi at all -- no lineage/taxonomy source, so tax_path stays empty
        # and the taxon is recorded as unmatched (not an error, just no lineage info).
        assert row["tax_path"] == ""
        assert "Puma" in summary["unmatched_taxa_examples"]


def _megan_fillet_metadata():
    """2 real libraries + 1 blank, with megan_library_name AND
    fillet_library_name crosswalk columns but no holi_library_name at all."""
    return _metadata({
        "merged_library_name": ["Real1", "Real2", "Blank1"],
        "megan_library_name":  ["Real1", "Real2", "Blank1"],
        "fillet_library_name": ["Real1_fillet", "Real2_fillet", ""],
        "is_negative_control": [False, False, True],
        "sample_type":         ["sediment", "sediment", "extraction blank"],
    })


def _megan_df(counts=(100.0, 90.0, 0.0)):
    return pd.DataFrame({
        "tax_id":     [9696],
        "tax_id_str": ["9696"],
        "tax_name":   ["Puma"],
        "tax_rank":   ["genus"],
        "Real1":      [counts[0]],
        "Real2":      [counts[1]],
        "Blank1":     [counts[2]],
    })


class TestMeganFilletNoHoli:
    """MEGAN+Fillet, no Holi at all -- previously unsupported, now real."""

    def test_reaches_high_confidence(self, config):
        merged, summary = build_merge(
            _megan_fillet_metadata(), _megan_df(), None, config, fillet_df=_fillet_df()
        )
        assert summary["n_taxa"] == 1
        row = merged.iloc[0]
        assert row["scientific_name"] == "Puma"
        assert row["aDNA_support_status"] in {"Very high confidence", "High confidence"}
        assert bool(row["fillet_authenticated"]) is True
        assert row["ensemble_support_score"] > 0.0
        # No Holi at all -- no lineage/taxonomy source, so no lineage/discordance signal.
        assert row["tax_path"] == ""
        assert pd.isna(row["Holi_best_damage"])

    def test_count_columns_come_from_megan(self, config):
        merged, _ = build_merge(
            _megan_fillet_metadata(), _megan_df(), None, config, fillet_df=_fillet_df()
        )
        row = merged.iloc[0]
        assert row["count__Real1"] == pytest.approx(100.0)
        assert row["count__Real2"] == pytest.approx(90.0)
        assert row["megan_max_count"] == pytest.approx(100.0)

    def test_blank_associated_from_megan_counts_alone(self, config):
        megan = _megan_df(counts=(2.0, 1.0, 100.0))
        merged, _ = build_merge(
            _megan_fillet_metadata(), megan, None, config, fillet_df=_fillet_df()
        )
        assert merged.iloc[0]["aDNA_support_status"] == "Blank-associated"

    def test_fillet_alone_no_megan_corroboration_still_supported(self, config):
        """Fillet authenticated but MEGAN counts too weak to corroborate --
        still reaches at least 'Supported' via Fillet's own strength."""
        megan = _megan_df(counts=(1.0, 0.0, 0.0))  # below strong_count thresholds
        merged, _ = build_merge(
            _megan_fillet_metadata(), megan, None, config, fillet_df=_fillet_df()
        )
        assert merged.iloc[0]["aDNA_support_status"] in {
            "Supported", "High confidence", "Very high confidence",
        }


class TestValidation:
    def test_all_three_none_raises(self, config):
        with pytest.raises(ValueError, match="at least one"):
            build_merge(_base_metadata(), None, None, config, fillet_df=None)

    def test_megan_alone_with_no_holi_and_no_fillet_raises(self, config):
        """MEGAN alone (no Holi, no Fillet at all) is still not supported --
        MEGAN+Fillet-without-Holi now IS supported (see TestMeganFilletNoHoli
        below), so this specifically covers megan_df being the ONLY source."""
        megan = pd.DataFrame({
            "tax_id": [9696], "tax_name": ["Puma"], "tax_rank": ["genus"], "tax_id_str": ["9696"],
        })
        with pytest.raises(ValueError, match="neither holi_df nor fillet_df"):
            build_merge(_base_metadata(), megan, None, config, fillet_df=None)

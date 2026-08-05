"""Tests for report.py's make_plot_input() -- specifically the
plot_methods_agreement column that drives the R heatmap's new ensemble-tier
badge (Tyler, 2026-08-05)."""

import pandas as pd
import pytest

from metamerge.config import load_config
from metamerge.report import make_plot_input


@pytest.fixture
def config():
    return load_config(None)


def _merged_row(**overrides):
    row = {
        "tax_id": 9696, "scientific_name": "Puma", "tax_rank": "genus",
        "tax_path": "Animalia;Chordata;Mammalia;Carnivora;Felidae;Puma",
        "common_name": "Puma", "broad_group": "animals",
        "aDNA_support_status": "Very high confidence",
        "megan_max_count": 100.0,
        "count__Real1": 100.0, "count__Blank1": 0.0,
        "aDNA_support_lib__Real1": "Damage-supported (>=100 reads)",
        "aDNA_support_lib__Blank1": "Not detected",
        "methods_agreement_lib__Real1": 2,
        "methods_agreement_lib__Blank1": 0,
        "Holi_damage_lib__Real1": 0.12, "Holi_damage_lib__Blank1": float("nan"),
    }
    row.update(overrides)
    return row


def _metadata():
    return pd.DataFrame({
        "merged_library_name": ["Real1", "Blank1"],
        "megan_library_name":  ["Real1", "Blank1"],
        "holi_library_name":   ["Real1_holi", "Blank1_holi"],
        "is_negative_control": [False, True],
        "sample_type":         ["sediment", "extraction blank"],
    })


class TestPlotMethodsAgreement:
    def test_agreement_column_present_and_correct(self, config):
        merged = pd.DataFrame([_merged_row()])
        table = make_plot_input(merged, _metadata(), broad_group="animals", config=config)
        assert not table.empty
        assert "plot_methods_agreement" in table.columns
        # Only Real1 (the real library) has count > 0, so exactly one row survives.
        assert len(table) == 1
        assert table.iloc[0]["plot_methods_agreement"] == 2

    def test_missing_agreement_columns_default_to_zero(self, config):
        """A merged_df from an older MetaMerge version with no
        methods_agreement_lib__ columns at all should still produce valid
        plot input, just with agreement=0 (no badge drawn)."""
        row = _merged_row()
        del row["methods_agreement_lib__Real1"]
        del row["methods_agreement_lib__Blank1"]
        merged = pd.DataFrame([row])
        table = make_plot_input(merged, _metadata(), broad_group="animals", config=config)
        assert not table.empty
        assert (table["plot_methods_agreement"] == 0).all()

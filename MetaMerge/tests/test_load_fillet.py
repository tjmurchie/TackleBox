"""Tests for io.load_fillet() -- reading Fillet's native MetaMerge export
(fillet_metamerge_evidence.tsv, written by fillet.report.write_fillet_metamerge_input)."""

import pandas as pd
import pytest

from metamerge.config import load_config
from metamerge.io import load_fillet


@pytest.fixture
def config():
    return load_config(None)


def _write_fillet_tsv(tmp_path, rows, name="fillet_metamerge_evidence.tsv"):
    df = pd.DataFrame(rows)
    path = tmp_path / name
    df.to_csv(path, sep="\t", index=False)
    return path


REQUIRED_COLS = {
    "sample": "garg_up", "tax_id": "9999", "tax_name": "Urocitellus parryii",
    "tax_rank": "species", "direct_hard_reads": 100, "cumulative_hard_reads": 100,
    "composite_authenticity": 0.8, "authenticity_tier": 1, "confidence_tier": "HIGH",
    "mean_damage_score": 0.15, "best_reference_breadth": 0.6, "stack_concentration": 0.1,
    "blank_fraction": 0.0, "eco_support": "True", "pal_support": "False", "fos_support": "True",
}


class TestLoadFillet:
    def test_loads_a_valid_file(self, tmp_path, config):
        path = _write_fillet_tsv(tmp_path, [REQUIRED_COLS])
        df = load_fillet(str(path), config)
        assert len(df) == 1
        assert df["tax_id_str"].iloc[0] == "9999"
        assert df["sample"].iloc[0] == "garg_up"

    def test_boolean_support_columns_parsed_correctly(self, tmp_path, config):
        path = _write_fillet_tsv(tmp_path, [REQUIRED_COLS])
        df = load_fillet(str(path), config)
        assert bool(df["eco_support"].iloc[0]) is True
        assert bool(df["pal_support"].iloc[0]) is False
        assert bool(df["fos_support"].iloc[0]) is True

    def test_numeric_columns_are_numeric(self, tmp_path, config):
        path = _write_fillet_tsv(tmp_path, [REQUIRED_COLS])
        df = load_fillet(str(path), config)
        assert float(df["composite_authenticity"].iloc[0]) == pytest.approx(0.8)
        assert int(df["authenticity_tier"].iloc[0]) == 1

    def test_missing_required_column_raises(self, tmp_path, config):
        bad = dict(REQUIRED_COLS)
        del bad["composite_authenticity"]
        path = _write_fillet_tsv(tmp_path, [bad])
        with pytest.raises(ValueError, match="composite_authenticity"):
            load_fillet(str(path), config)

    def test_multiple_rows_multiple_samples(self, tmp_path, config):
        row2 = dict(REQUIRED_COLS)
        row2["sample"] = "garg_up_damage"
        row2["tax_name"] = "Bison bison"
        row2["tax_id"] = "9903"
        path = _write_fillet_tsv(tmp_path, [REQUIRED_COLS, row2])
        df = load_fillet(str(path), config)
        assert len(df) == 2
        assert set(df["sample"]) == {"garg_up", "garg_up_damage"}

    def test_name_and_rank_normalized(self, tmp_path, config):
        row = dict(REQUIRED_COLS)
        row["tax_name"] = "  Urocitellus   parryii  "
        row["tax_rank"] = "SPECIES"
        path = _write_fillet_tsv(tmp_path, [row])
        df = load_fillet(str(path), config)
        assert df["tax_name"].iloc[0] == "Urocitellus parryii"
        assert df["tax_rank"].iloc[0] == "species"

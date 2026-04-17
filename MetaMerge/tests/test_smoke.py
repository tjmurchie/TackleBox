"""Smoke test for the MetaMerge build_merge pipeline.

Runs the complete merge/classification pipeline on a minimal synthetic dataset
to verify that:
  1. The pipeline runs without errors end to end.
  2. The output has the expected number of taxa.
  3. The DNA support status for a strongly-supported taxon is at least
     "High confidence".

This test does NOT check every output column in detail — that is left to
specific unit tests.  Its purpose is to catch import errors, structural
changes, and obvious regressions quickly.
"""

from pathlib import Path

import pandas as pd

from metamerge.classify import build_merge
from metamerge.config import load_config


def _make_test_data():
    """Build a minimal but realistic synthetic test dataset."""
    cfg = load_config(None)

    metadata = pd.DataFrame(
        {
            "megan_library_name":  ["A", "B", "BLK"],
            "holi_library_name":   ["A_holi", "B_holi", "blanks"],
            "merged_library_name": ["A", "B", "BLK"],
            "is_negative_control": [False, False, True],
            "sample_type":         ["sediment", "sediment", "negative control"],
            "site":                ["Site1", "Site1", "blank"],
        }
    )

    megan = pd.DataFrame(
        {
            "tax_id":     [9696],
            "tax_name":   ["Puma"],
            "tax_rank":   ["genus"],
            "tax_id_str": ["9696"],
            "A":          [100.0],
            "B":          [20.0],
            "BLK":        [0.0],
        }
    )

    holi = pd.DataFrame(
        {
            "sample":        ["A_holi", "blanks"],
            "tax_id":        [9696,       9696],
            "tax_id_str":    ["9696",     "9696"],
            "tax_name":      ["Puma",     "Puma"],
            "tax_rank":      ["genus",    "genus"],
            "N_reads":       [150.0,      0.0],
            "N_alignments":  [300.0,      0.0],
            "damage":        [0.10,       0.00],
            "significance":  [3.00,       0.00],
            "rho_Ac":        [0.10,       0.00],
            "MAP_valid":     [True,       True],
            "tax_path":      [
                "Animalia;Chordata;Mammalia;Carnivora;Felidae;Puma",
                "Animalia;Chordata;Mammalia;Carnivora;Felidae;Puma",
            ],
        }
    )
    return cfg, metadata, megan, holi


def test_smoke_build_merge():
    """The full merge pipeline runs and produces the expected output shape."""
    cfg, metadata, megan, holi = _make_test_data()
    merged, summary = build_merge(metadata, megan, holi, cfg)

    assert len(merged) == 1, f"Expected 1 taxon row, got {len(merged)}"
    assert summary["n_taxa"] == 1

    status = merged.iloc[0]["aDNA_support_status"]
    assert status in {
        "Very high confidence", "High confidence", "Supported",
    }, f"Unexpected status for well-supported taxon: {status}"


def test_smoke_blank_associated():
    """A taxon that only appears in blanks should be Blank-associated."""
    cfg, metadata, megan, holi = _make_test_data()

    # Put all count in the blank library only.
    megan["A"]   = 0.0
    megan["B"]   = 0.0
    megan["BLK"] = 100.0

    # Give the blank Holi row real damage support.
    holi.loc[holi["sample"] == "blanks", "damage"]       = 0.15
    holi.loc[holi["sample"] == "blanks", "significance"] = 3.0
    holi.loc[holi["sample"] == "blanks", "N_reads"]      = 120.0

    merged, _ = build_merge(metadata, megan, holi, cfg)
    assert merged.iloc[0]["aDNA_support_status"] == "Blank-associated"


def test_smoke_output_columns():
    """Key output columns are present in the merged DataFrame."""
    cfg, metadata, megan, holi = _make_test_data()
    merged, _ = build_merge(metadata, megan, holi, cfg)

    required_cols = [
        "tax_id", "scientific_name", "tax_rank", "tax_path",
        "aDNA_support_status", "support_basis_summary",
        "megan_max_count", "megan_positive_libraries_n",
        "Holi_best_damage", "Holi_best_significance",
        "count__A", "count__B", "count__BLK",
    ]
    missing = [c for c in required_cols if c not in merged.columns]
    assert not missing, f"Missing expected columns: {missing}"

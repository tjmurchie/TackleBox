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


def test_smoke_blank_associated_not_overridden_by_strong_caution_multimapping():
    """Regression test for the A51-LBL-1 Homo sapiens misclassification:
    an exact-damage-support call built on a "strong caution" (heavily
    multi-mapping) Holi fit must not be allowed to unconditionally override
    a substantial blank overlap. Modeled directly on the real numbers:
    N_reads=28012/N_alignments=1140245 (40.7 alignments/read, over double
    qc_alignments_per_read_caution_max=20.0) and blank MEGAN count 12.8x the
    real-sample count, with the blank's own Holi row deliberately kept below
    the damage/significance thresholds so only the QC-label branch is
    exercised (not the pre-existing blank_damage_support branch)."""
    cfg, metadata, megan, holi = _make_test_data()

    # Weak, single-library real MEGAN support (mirrors LBL-1: only 1 of 3
    # libraries positive, 6 reads) and a blank count that dominates it.
    megan["A"]   = 6.0
    megan["B"]   = 0.0
    megan["BLK"] = 77.0

    # Real library's Holi row: clears damage/significance (exact_damage_support
    # True) but is heavily multi-mapping -> qc_label "strong caution".
    holi.loc[holi["sample"] == "A_holi", "damage"]       = 0.0201
    holi.loc[holi["sample"] == "A_holi", "significance"] = 2.46
    holi.loc[holi["sample"] == "A_holi", "N_reads"]      = 28012.0
    holi.loc[holi["sample"] == "A_holi", "N_alignments"] = 1140245.0

    # Blank's own Holi row stays below threshold -- isolates the new branch.
    holi.loc[holi["sample"] == "blanks", "damage"]       = 0.011
    holi.loc[holi["sample"] == "blanks", "significance"] = 1.21

    merged, _ = build_merge(metadata, megan, holi, cfg)
    row = merged.iloc[0]
    assert row["Holi_best_multimapping_fit_qc"] == "strong caution"
    assert row["aDNA_support_status"] == "Blank-associated"


def test_smoke_blank_associated_depth_normalized():
    """A taxon with a small RAW blank count -- below BOTH blank_absolute_min
    and blank_relative_min, so the pre-existing raw-count checks alone would
    NOT flag it -- should still be Blank-associated if it is proportionally
    far more concentrated in a shallowly-sequenced blank than in a deeply-
    sequenced real library (reads-per-million of each library's own total
    classified-read count). Motivating real case: A51 project blanks were
    sequenced far more shallowly than real samples, and several background
    bacteria (e.g. Cutibacterium acnes) had small raw blank counts that
    looked negligible by raw-count comparison but were orders of magnitude
    more concentrated in the blanks once normalized for sequencing depth."""
    cfg, metadata, megan, holi = _make_test_data()

    # Focal taxon (row 0): real library "A" has a raw count comfortably above
    # blank_absolute_min and >10x the blank's raw count (5), so the OLD
    # raw-count-only checks would NOT flag blank_caution at all:
    #   blank_absolute_min=10  -> 5 < 10, no trigger
    #   blank_relative_min=0.10 -> 5/100 = 0.05 < 0.10, no trigger
    # "B" stays at 0 (non-positive) so strong_count_support requires >=2
    # positive real libraries and is correctly False here regardless of A's
    # count, letting blank_associated actually fire.
    megan["A"]   = 100.0
    megan["B"]   = 0.0
    megan["BLK"] = 5.0

    # Filler taxon (a second megan row) inflates each library's TOTAL
    # classified-read depth independently of the focal taxon's own count --
    # library "A" is deeply sequenced (2,000,000 total), "BLK" is shallow
    # (50 total), mirroring the real A51 blank-vs-sample depth gap.
    filler = pd.DataFrame({
        "tax_id": [1], "tax_name": ["Filler"], "tax_rank": ["genus"],
        "tax_id_str": ["1"], "A": [1_999_900.0], "B": [0.0], "BLK": [45.0],
    })
    megan = pd.concat([megan, filler], ignore_index=True)

    # Real library's Holi row stays below damage/significance thresholds so
    # exact_damage_support is False -- isolating the depth-normalized branch
    # from the pre-existing blank_damage_support / strong-caution branches.
    holi.loc[holi["sample"] == "A_holi", "damage"]       = 0.0
    holi.loc[holi["sample"] == "A_holi", "significance"] = 0.0
    holi.loc[holi["sample"] == "blanks", "damage"]       = 0.0
    holi.loc[holi["sample"] == "blanks", "significance"] = 0.0

    merged, _ = build_merge(metadata, megan, holi, cfg)
    row = merged[merged["tax_id_str"] == "9696"].iloc[0]

    # real RPM = 100/2,000,000 * 1e6 = 50; blank RPM = 5/50 * 1e6 = 100,000.
    # real_rpm (50) is nowhere near 10x blank_rpm (100,000) -> depth-normalized
    # blank concern should fire even though raw counts alone would not.
    assert row["blank_caution_depth_normalized"], (
        f"expected depth-normalized blank concern to fire: "
        f"real_rpm={row['real_rpm']}, blank_rpm={row['blank_rpm']}"
    )
    assert row["aDNA_support_status"] == "Blank-associated"


def test_smoke_blank_depth_normalize_can_be_disabled():
    """--disable-blank-depth-normalize (blank_depth_normalize=False in config)
    restores the old raw-count-only behavior for the same scenario as
    test_smoke_blank_associated_depth_normalized -- with it off, a small raw
    blank count that fails BOTH raw-count checks should NOT be Blank-
    associated on its own (no other blank/damage signal present)."""
    cfg, metadata, megan, holi = _make_test_data()

    megan["A"]   = 100.0
    megan["B"]   = 0.0
    megan["BLK"] = 5.0
    filler = pd.DataFrame({
        "tax_id": [1], "tax_name": ["Filler"], "tax_rank": ["genus"],
        "tax_id_str": ["1"], "A": [1_999_900.0], "B": [0.0], "BLK": [45.0],
    })
    megan = pd.concat([megan, filler], ignore_index=True)
    holi.loc[holi["sample"] == "A_holi", "damage"]       = 0.0
    holi.loc[holi["sample"] == "A_holi", "significance"] = 0.0
    holi.loc[holi["sample"] == "blanks", "damage"]       = 0.0
    holi.loc[holi["sample"] == "blanks", "significance"] = 0.0

    cfg["thresholds"]["blank_depth_normalize"] = False
    merged, _ = build_merge(metadata, megan, holi, cfg)
    row = merged[merged["tax_id_str"] == "9696"].iloc[0]
    assert row["aDNA_support_status"] != "Blank-associated"


def test_smoke_blank_depth_normalized_overrides_strong_count_support():
    """PI decision (2026-10-21, real A51 case: Cutibacterium acnes): the
    depth-normalized blank check is an UNCONDITIONAL override -- same
    priority tier as an exact blank-damage match -- that beats even
    strong_count_support (>=2 real libraries each with >=50 RAW reads,
    which has no depth-awareness of its own). Real libraries "A" and "B"
    both clear strong_count_min_reads/min_libraries here, which would have
    exempted this taxon from Blank-associated under the OLD (conditional)
    design; with the depth-normalized check promoted to unconditional, it
    must still be forced to Blank-associated when real_rpm fails to clear
    blank_depth_fold_threshold x blank_rpm."""
    cfg, metadata, megan, holi = _make_test_data()

    # Both real libraries positive and >=50 reads -> strong_count_support
    # would be True under the old design.
    megan["A"]   = 100.0
    megan["B"]   = 60.0
    megan["BLK"] = 5.0
    filler = pd.DataFrame({
        "tax_id": [1], "tax_name": ["Filler"], "tax_rank": ["genus"],
        "tax_id_str": ["1"], "A": [499_900.0], "B": [599_940.0], "BLK": [45.0],
    })
    megan = pd.concat([megan, filler], ignore_index=True)

    # Keep exact_damage_support False on the real side too, so this test
    # isolates the strong_count_support interaction specifically (a separate
    # test below covers the exact-damage-support interaction).
    holi.loc[holi["sample"] == "A_holi", "damage"]       = 0.0
    holi.loc[holi["sample"] == "A_holi", "significance"] = 0.0
    holi.loc[holi["sample"] == "blanks", "damage"]       = 0.0
    holi.loc[holi["sample"] == "blanks", "significance"] = 0.0

    merged, _ = build_merge(metadata, megan, holi, cfg)
    row = merged[merged["tax_id_str"] == "9696"].iloc[0]

    assert row["strong_count_support"], "test setup should produce strong_count_support=True"
    # real RPM = max(100/500,000, 60/600,000) * 1e6 = 200; blank RPM = 5/50 * 1e6 = 100,000.
    assert row["blank_caution_depth_normalized"]
    assert row["aDNA_support_status"] == "Blank-associated"


def test_smoke_blank_depth_normalized_overrides_exact_damage_support():
    """Same unconditional-override principle, but isolating the OTHER
    pre-existing escape hatch: genuine exact Holi/metaDMG damage support
    with a clean QC fit. Even a real library with authentic-looking ancient
    damage must still be forced to Blank-associated if its depth-normalized
    rate fails to clear the real-vs-blank fold threshold."""
    cfg, metadata, megan, holi = _make_test_data()

    megan["A"]   = 100.0
    megan["B"]   = 0.0
    megan["BLK"] = 5.0
    filler = pd.DataFrame({
        "tax_id": [1], "tax_name": ["Filler"], "tax_rank": ["genus"],
        "tax_id_str": ["1"], "A": [1_999_900.0], "B": [0.0], "BLK": [45.0],
    })
    megan = pd.concat([megan, filler], ignore_index=True)

    # Real library "A" gets genuine, clean-fit exact damage support.
    holi.loc[holi["sample"] == "A_holi", "damage"]       = 0.10
    holi.loc[holi["sample"] == "A_holi", "significance"] = 3.0
    holi.loc[holi["sample"] == "A_holi", "N_reads"]      = 150.0
    holi.loc[holi["sample"] == "blanks", "damage"]       = 0.0
    holi.loc[holi["sample"] == "blanks", "significance"] = 0.0

    merged, _ = build_merge(metadata, megan, holi, cfg)
    row = merged[merged["tax_id_str"] == "9696"].iloc[0]

    assert row["blank_caution_depth_normalized"]
    assert row["aDNA_support_status"] == "Blank-associated"


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
        "real_rpm", "blank_rpm", "blank_caution_depth_normalized",
    ]
    missing = [c for c in required_cols if c not in merged.columns]
    assert not missing, f"Missing expected columns: {missing}"

"""Tests for metadata.load_metadata()'s N-of-3 classifier crosswalk relaxation.

Historically the linker file required exactly megan_library_name +
holi_library_name + merged_library_name. Since a real run may now use any
non-empty subset of {MEGAN, Holi, Fillet} -- e.g. Holi+Fillet with no MEGAN at
all -- only merged_library_name is unconditionally required, and at least one
of the three classifier crosswalk columns must be present (but not all three).
"""

import pandas as pd
import pytest

from metamerge.config import load_config
from metamerge.metadata import load_metadata


@pytest.fixture
def config():
    return load_config(None)


def _write_linker(tmp_path, df: pd.DataFrame, name="linker.csv"):
    path = tmp_path / name
    df.to_csv(path, index=False)
    return path


class TestLegacyTwoSourceLinkerStillWorks:
    """The original megan_library_name + holi_library_name + merged_library_name
    linker shape must keep working unchanged."""

    def test_classic_megan_holi_linker_loads(self, tmp_path, config):
        df = pd.DataFrame({
            "megan_library_name": ["lib1_S1"],
            "holi_library_name": ["lib1"],
            "merged_library_name": ["Sample1"],
        })
        path = _write_linker(tmp_path, df)
        loaded = load_metadata(str(path), config)
        assert "megan_library_name" in loaded.columns
        assert "holi_library_name" in loaded.columns
        assert loaded["merged_library_name"].iloc[0] == "Sample1"


class TestHoliFilletOnlyLinker:
    """A linker with holi_library_name + fillet_library_name but no
    megan_library_name column at all must load successfully."""

    def test_holi_fillet_linker_loads_without_megan_column(self, tmp_path, config):
        df = pd.DataFrame({
            "holi_library_name": ["lib1"],
            "fillet_library_name": ["garg_up"],
            "merged_library_name": ["Sample1"],
        })
        path = _write_linker(tmp_path, df)
        loaded = load_metadata(str(path), config)
        assert "megan_library_name" not in loaded.columns
        assert "holi_library_name" in loaded.columns
        assert "fillet_library_name" in loaded.columns
        assert loaded["fillet_library_name"].iloc[0] == "garg_up"


class TestFilletOnlyLinker:
    def test_fillet_only_linker_loads(self, tmp_path, config):
        df = pd.DataFrame({
            "fillet_library_name": ["garg_up"],
            "merged_library_name": ["Sample1"],
        })
        path = _write_linker(tmp_path, df)
        loaded = load_metadata(str(path), config)
        assert "fillet_library_name" in loaded.columns


class TestMissingRequiredColumns:
    def test_missing_merged_library_name_raises(self, tmp_path, config):
        df = pd.DataFrame({
            "holi_library_name": ["lib1"],
            "fillet_library_name": ["garg_up"],
        })
        path = _write_linker(tmp_path, df)
        with pytest.raises(ValueError, match="merged_library_name"):
            load_metadata(str(path), config)

    def test_no_classifier_crosswalk_columns_at_all_raises(self, tmp_path, config):
        df = pd.DataFrame({
            "merged_library_name": ["Sample1"],
            "site": ["Somewhere"],
        })
        path = _write_linker(tmp_path, df)
        with pytest.raises(ValueError, match="crosswalk"):
            load_metadata(str(path), config)


class TestDuplicateMergedLibraryNameStillValidated:
    def test_duplicate_merged_library_name_still_raises(self, tmp_path, config):
        df = pd.DataFrame({
            "fillet_library_name": ["garg_up", "garg_up_damage"],
            "merged_library_name": ["Sample1", "Sample1"],
        })
        path = _write_linker(tmp_path, df)
        with pytest.raises(ValueError, match="unique"):
            load_metadata(str(path), config)

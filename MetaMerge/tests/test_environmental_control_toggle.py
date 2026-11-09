"""Tests for the ``controls.enable_environmental_control`` config toggle.

Verifies that disabling the environmental-control category (e.g. via the
``--disable-environmental-control`` CLI flag) makes ``load_metadata`` treat
spiderweb/wash/swab/surface libraries as ordinary samples, and makes
``_derive_plot_group`` fall through to the normal linker ``group`` value
instead of a separate "environmental_controls" / "<site>_webs" bucket.
"""

import copy

import pandas as pd

from metamerge.defaults import DEFAULT_CONFIG
from metamerge.metadata import load_metadata
from metamerge.report import _derive_plot_group


def _write_linker(tmp_path):
    path = tmp_path / "linker.csv"
    pd.DataFrame([
        {
            "megan_library_name": "Web-c-Shg",
            "holi_library_name": "Site-Web",
            "merged_library_name": "Web-c-Shg",
            "sample_id": "Site-Web-c",
            "sample_type": "spiderweb",
            "group": "Site_sediments",
            "site": "Site",
        },
        {
            "megan_library_name": "Sed-c-Shg",
            "holi_library_name": "Site-Sed",
            "merged_library_name": "Sed-c-Shg",
            "sample_id": "Site-Sed-c",
            "sample_type": "sediment",
            "group": "Site_sediments",
            "site": "Site",
        },
    ]).to_csv(path, index=False)
    return path


def test_environmental_control_enabled_by_default(tmp_path):
    config = copy.deepcopy(DEFAULT_CONFIG)
    metadata = load_metadata(str(_write_linker(tmp_path)), config)
    web_row = metadata.loc[metadata["megan_library_name"] == "Web-c-Shg"].iloc[0]
    assert bool(web_row["is_environmental_control"]) is True
    assert bool(web_row["is_negative_control"]) is False


def test_environmental_control_disabled_treats_web_as_ordinary_sample(tmp_path):
    config = copy.deepcopy(DEFAULT_CONFIG)
    config["controls"]["enable_environmental_control"] = False
    metadata = load_metadata(str(_write_linker(tmp_path)), config)
    web_row = metadata.loc[metadata["megan_library_name"] == "Web-c-Shg"].iloc[0]
    # Disabled means: not flagged as environmental, and NOT accidentally
    # reclassified as a negative control either (sample_type "spiderweb"
    # contains no "blank"/"negative"/"control" substring).
    assert bool(web_row["is_environmental_control"]) is False
    assert bool(web_row["is_negative_control"]) is False


def test_plot_group_environmental_controls_bucket_when_enabled():
    row = pd.Series({
        "site": "Site", "sample_type": "spiderweb", "group": "Site_sediments",
        "plot_sample_id": "Site-Web-c", "is_negative_control": "False",
        "is_positive_control": "False", "is_environmental_control": "True",
    })
    assert _derive_plot_group(row, enable_environmental_control=True) == "environmental_controls"


def test_plot_group_falls_back_to_linker_group_when_disabled():
    row = pd.Series({
        "site": "Site", "sample_type": "spiderweb", "group": "Site_sediments",
        "plot_sample_id": "Site-Web-c", "is_negative_control": "False",
        "is_positive_control": "False", "is_environmental_control": "True",
    })
    # Even though is_environmental_control is (still) True on the row, the
    # disabled flag must override both the environmental_controls bucket AND
    # the "web"-substring fallback, landing on the ordinary linker group.
    assert _derive_plot_group(row, enable_environmental_control=False) == "Site_sediments"


def test_negative_control_gets_own_bucket_by_default():
    row = pd.Series({
        "site": "Site", "sample_type": "blank", "group": "Site_sediments",
        "plot_sample_id": "Blank-1", "is_negative_control": "True",
        "is_positive_control": "False", "is_environmental_control": "False",
    })
    assert _derive_plot_group(row) == "negative_controls"


def test_negative_control_falls_back_to_linker_group_when_kept_with_group():
    row = pd.Series({
        "site": "Site", "sample_type": "blank", "group": "Site_sediments",
        "plot_sample_id": "Blank-1", "is_negative_control": "True",
        "is_positive_control": "False", "is_environmental_control": "False",
    })
    assert _derive_plot_group(row, keep_negative_controls_with_group=True) == "Site_sediments"

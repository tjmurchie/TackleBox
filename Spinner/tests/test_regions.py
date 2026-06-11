"""Tests for marker/region classification and taxonomy guessing (spinner.regions)."""
from __future__ import annotations

import copy
from pathlib import Path

import pytest

from spinner.config import DEFAULT_CONFIG
from spinner.regions import classify, guess_species, load_regions, load_species_kingdom, norm_kingdom


# ---------------------------------------------------------------------------
# guess_species
# ---------------------------------------------------------------------------

def test_guess_species_ncbi_format():
    header = ">NC_012345.1 Rangifer tarandus mitochondrion complete genome"
    sp, gen = guess_species(header)
    assert sp == "Rangifer tarandus"
    assert gen == "Rangifer"


def test_guess_species_bracket_format():
    # Use "16S" as the second token (starts with digit) so the first-format
    # heuristic doesn't match and we fall through to the bracket regex.
    header = ">AX123456.1 16S rRNA gene partial [Betula papyrifera]"
    sp, gen = guess_species(header)
    assert sp == "Betula papyrifera"
    assert gen == "Betula"


def test_guess_species_no_match():
    header = ">UNKN.1 some sequence"
    sp, gen = guess_species(header)
    assert sp == ""
    assert gen == ""


# ---------------------------------------------------------------------------
# norm_kingdom
# ---------------------------------------------------------------------------

def test_norm_kingdom_animal():
    assert norm_kingdom("Animalia") == "Animal"
    assert norm_kingdom("Metazoa") == "Animal"


def test_norm_kingdom_plant():
    assert norm_kingdom("Plantae") == "Plant"
    assert norm_kingdom("Viridiplantae") == "Plant"


def test_norm_kingdom_fungi():
    assert norm_kingdom("Fungi") == "Fungi"


def test_norm_kingdom_unknown():
    assert norm_kingdom("") == "Unknown"
    assert norm_kingdom(None) == "Unknown"


# ---------------------------------------------------------------------------
# classify — built-in heuristics
# ---------------------------------------------------------------------------

def _cfg():
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    return cfg


def test_classify_mito_builtin():
    header = ">NC_001234.1 Rangifer tarandus mitochondrion complete genome"
    klass, rid = classify(header, [], _cfg())
    assert klass == "Mito"
    assert rid == "MITO_BUILTIN"


def test_classify_plastid_builtin():
    header = ">MK123456.1 Betula papyrifera chloroplast rbcL gene"
    klass, rid = classify(header, [], _cfg())
    assert klass == "Plastid"
    assert rid == "PLASTID_BUILTIN"


def test_classify_nucmark_18s_builtin():
    header = ">AY123456.1 Salix arctica 18S ribosomal RNA gene"
    klass, rid = classify(header, [], _cfg())
    assert klass == "NucMark"
    assert rid == "NUCMARK_BUILTIN"


def test_classify_nucmark_its_builtin():
    header = ">AY999999.1 Salix arctica internal transcribed spacer 1"
    klass, rid = classify(header, [], _cfg())
    assert klass == "NucMark"


def test_classify_other_builtin():
    header = ">ZZ000001.1 Unknown organism partial sequence"
    klass, rid = classify(header, [], _cfg())
    assert klass == "Other"


# ---------------------------------------------------------------------------
# classify — from regions_config
# ---------------------------------------------------------------------------

def test_classify_from_regions_config(regions_file):
    rules = load_regions(str(regions_file))
    header = ">NC_001234.1 Rangifer tarandus mitochondrion complete genome"
    klass, rid = classify(header, rules, _cfg())
    assert klass == "Mito"
    assert rid == "MITO"


def test_classify_plastid_from_regions_config(regions_file):
    rules = load_regions(str(regions_file))
    header = ">MK123456.1 Betula papyrifera chloroplast genome"
    klass, rid = classify(header, rules, _cfg())
    assert klass == "Plastid"
    assert rid == "PLASTID"


def test_classify_nucmark_from_regions_config(regions_file):
    rules = load_regions(str(regions_file))
    header = ">AY123456.1 Salix arctica 18S ribosomal RNA gene"
    klass, rid = classify(header, rules, _cfg())
    assert klass == "NucMark"
    assert rid == "NUCRDNA_18S"


# ---------------------------------------------------------------------------
# load_species_kingdom
# ---------------------------------------------------------------------------

def test_load_species_kingdom_header(species_kingdom_file):
    sp2k, g2k = load_species_kingdom(str(species_kingdom_file))
    assert sp2k.get("rangifer tarandus") == "Animal"
    assert sp2k.get("betula papyrifera") == "Plant"
    assert sp2k.get("salix arctica") == "Plant"
    assert g2k.get("rangifer") == "Animal"


def test_load_species_kingdom_missing_file():
    sp2k, g2k = load_species_kingdom("/no/such/file.tsv")
    assert sp2k == {}
    assert g2k == {}


def test_load_species_kingdom_empty_path():
    sp2k, g2k = load_species_kingdom("")
    assert sp2k == {}
    assert g2k == {}

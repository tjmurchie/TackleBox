"""Tests for marker/region classification and taxonomy guessing (spinner.regions)."""
from __future__ import annotations

import copy


from spinner.config import DEFAULT_CONFIG
from spinner.regions import (
    build_genus_abbrev_index,
    classify,
    guess_species,
    load_regions,
    load_species_kingdom,
    norm_kingdom,
)


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
# guess_species — abbreviated-genus headers (regression: Bug A)
#
# Real accession Z54097.1: "M.primigenius" is a legacy NCBI abbreviated
# Genus.species token smashed into ONE field.  The naive parser used to treat
# fields[1]="M.primigenius" as a genus and fields[2]="mitochondrial" (an
# ordinary description word) as the species epithet, fabricating the bogus
# pseudo-species "M.primigenius mitochondrial" -- a singleton that bypassed
# real clustering/capping entirely against the correctly-parsed 544+ other
# "Mammuthus primigenius" records for the same real species.
# ---------------------------------------------------------------------------

_Z54097_HEADER = (
    ">Z54097.1 M.primigenius mitochondrial partial 16S rRNA gene "
    "(Khatanga mammoth; 93 bp)"
)


def test_guess_species_abbreviated_genus_no_fabrication_without_index():
    """Without a known-species index, must not fabricate a pseudo-species."""
    sp, gen = guess_species(_Z54097_HEADER)
    assert sp == ""
    assert gen == ""
    # In particular, must never fabricate this bogus two-word combination.
    assert sp != "M.primigenius mitochondrial"
    assert gen != "M.primigenius"


def test_guess_species_abbreviated_genus_resolved_with_index():
    """With a genus_abbrev_index, the abbreviation resolves to the real species."""
    index = {"m.primigenius": "Mammuthus"}
    sp, gen = guess_species(_Z54097_HEADER, index)
    assert sp == "Mammuthus primigenius"
    assert gen == "Mammuthus"


def test_guess_species_abbreviated_genus_unresolvable_key_still_safe():
    """An index that doesn't contain the abbreviation still must not fabricate."""
    index = {"r.tarandus": "Rangifer"}
    sp, gen = guess_species(_Z54097_HEADER, index)
    assert sp == ""
    assert gen == ""


def test_guess_species_full_genus_unaffected_by_abbrev_check():
    """A normal fully-spelled-genus header must be completely unaffected."""
    header = ">NC_012345.1 Rangifer tarandus mitochondrion complete genome"
    index = {"m.primigenius": "Mammuthus"}
    sp, gen = guess_species(header, index)
    assert sp == "Rangifer tarandus"
    assert gen == "Rangifer"


# ---------------------------------------------------------------------------
# build_genus_abbrev_index
# ---------------------------------------------------------------------------

def test_build_genus_abbrev_index_basic():
    pairs = [
        ("Mammuthus primigenius", "Mammuthus"),
        ("Rangifer tarandus", "Rangifer"),
    ]
    index = build_genus_abbrev_index(pairs)
    assert index == {"m.primigenius": "Mammuthus", "r.tarandus": "Rangifer"}


def test_build_genus_abbrev_index_skips_empty_and_abbreviated_entries():
    pairs = [
        ("", ""),                                  # no guess at all
        ("M.primigenius", "M.primigenius"),        # already-abbreviated genus — must be skipped
        ("Mammuthus primigenius", "Mammuthus"),    # genuine full pair
    ]
    index = build_genus_abbrev_index(pairs)
    assert index == {"m.primigenius": "Mammuthus"}


def test_build_genus_abbrev_index_first_seen_wins():
    pairs = [
        ("Mammuthus primigenius", "Mammuthus"),
        ("Mustela primigenius", "Mustela"),  # hypothetical clash on the same abbreviation
    ]
    index = build_genus_abbrev_index(pairs)
    assert index["m.primigenius"] == "Mammuthus"


def test_guess_species_abbreviated_genus_end_to_end_via_index_from_real_dataset():
    """End-to-end: build the index the way annotate() does, from a mix of
    headers including the correctly-parsed full-genus siblings of Z54097.1,
    and confirm the abbreviated record resolves to the identical species
    string used by its correctly-parsed siblings (required for correct
    clustering/capping grouping).
    """
    sibling_header = ">OQ000001.1 Mammuthus primigenius mitochondrion, complete genome"
    all_headers = [sibling_header, _Z54097_HEADER]
    index = build_genus_abbrev_index(guess_species(h) for h in all_headers)
    sp, gen = guess_species(_Z54097_HEADER, index)
    sibling_sp, sibling_gen = guess_species(sibling_header)
    assert sp == sibling_sp == "Mammuthus primigenius"
    assert gen == sibling_gen == "Mammuthus"


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

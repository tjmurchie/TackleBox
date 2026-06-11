"""Shared pytest fixtures for Spinner tests."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Inline FASTA content used across multiple test modules.
# ---------------------------------------------------------------------------

CLEAN_MITO = ">NC_CLEAN.1 Rangifer tarandus mitochondrion complete genome\n" + "ACGT" * 30 + "\n"
ADAPTER_INTERNAL = (
    ">ADAPTER_INT.1 Rangifer tarandus COI partial sequence\n"
    # Internal TruSeq adapter: positioned well away from both ends (120 padding + adapter + 120)
    + "A" * 120
    + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    + "A" * 120
    + "\n"
)
ADAPTER_TERMINAL = (
    ">ADAPTER_TERM.1 Rangifer tarandus COI gene\n"
    # Adapter starts at position 2 — within the 25 bp terminal window
    + "ACAGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    + "A" * 80
    + "\n"
)
DUP_ACCESSION_A = ">DUPSEQ.1 Betula papyrifera chloroplast rbcL gene\n" + "ACGT" * 18 + "\n"
DUP_ACCESSION_B = ">DUPSEQ.1 Betula papyrifera chloroplast rbcL gene (duplicate accession)\n" + "GGGG" * 18 + "\n"
DUP_SEQUENCE = (
    ">DUPSEQ2.1 Betula papyrifera chloroplast rbcL gene (duplicate sequence)\n"
    + "ACGT" * 18  # same sequence as DUP_ACCESSION_A
    + "\n"
)
HIGH_N = ">HIGHN.1 Salix arctica internal transcribed spacer 1\n" + "N" * 80 + "\n"
BAD_KEYWORD = ">BADVEC.1 Synthetic construct vector pUC19 sequence\n" + "ACGT" * 20 + "\n"
CLEAN_NUCLEAR = ">NUC18S.1 Salix arctica 18S ribosomal RNA gene\n" + "ACGTACGT" * 10 + "\n"
CLEAN_PLASTID = ">PLAST.1 Betula papyrifera chloroplast rbcL gene\n" + "GCTAGCTA" * 10 + "\n"


MINIMAL_FASTA = (
    CLEAN_MITO
    + ADAPTER_INTERNAL
    + ADAPTER_TERMINAL
    + DUP_ACCESSION_A
    + DUP_ACCESSION_B
    + DUP_SEQUENCE
    + HIGH_N
    + BAD_KEYWORD
    + CLEAN_NUCLEAR
    + CLEAN_PLASTID
)

SPECIES_KINGDOM_TSV = textwrap.dedent("""\
    species\tgenus\tkingdom
    Rangifer tarandus\tRangifer\tAnimalia
    Betula papyrifera\tBetula\tPlantae
    Salix arctica\tSalix\tPlantae
""")

ADAPTERS_TSV = textwrap.dedent("""\
    name\tsequence\tmax_mismatch\taction
    Illumina_TruSeq_Universal\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\t2\treject
    Illumina_short\tAGATCGGAAGAGC\t1\treject
""")

BAD_KEYWORDS_TSV = textwrap.dedent("""\
    keyword\taction\treason
    UNVERIFIED\treject\tNCBI record flagged as unverified
    synthetic construct\treject\tSynthetic/artificial construct
    vector\treject\tHeader suggests vector sequence
    clone\treview\tClone-derived record
    uncultured\treview\tUncultured source
""")

REGIONS_TSV = textwrap.dedent("""\
    region_id\tclass\tenabled_default\tregex\tncbi_title_clause
    MITO\tMito\t1\tmitochondri[on]|mitochondrial\t
    PLASTID\tPlastid\t1\tchloroplast|plastid|plastome\t
    NUCRDNA_18S\tNucMark\t1\t\\b18S\\b|small subunit ribosomal\t
    ITS\tNucMark\t1\tinternal transcribed spacer|\\bITS\\b\t
""")


@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
    p = tmp_path / "test.fasta"
    p.write_text(MINIMAL_FASTA, encoding="utf-8")
    return p


@pytest.fixture
def species_kingdom_file(tmp_path: Path) -> Path:
    p = tmp_path / "species_kingdom.tsv"
    p.write_text(SPECIES_KINGDOM_TSV, encoding="utf-8")
    return p


@pytest.fixture
def adapters_file(tmp_path: Path) -> Path:
    p = tmp_path / "adapters.tsv"
    p.write_text(ADAPTERS_TSV, encoding="utf-8")
    return p


@pytest.fixture
def bad_keywords_file(tmp_path: Path) -> Path:
    p = tmp_path / "bad_keywords.tsv"
    p.write_text(BAD_KEYWORDS_TSV, encoding="utf-8")
    return p


@pytest.fixture
def regions_file(tmp_path: Path) -> Path:
    p = tmp_path / "regions.tsv"
    p.write_text(REGIONS_TSV, encoding="utf-8")
    return p


@pytest.fixture
def default_cfg():
    """Return a clean DEFAULT_CONFIG copy."""
    from spinner.config import DEFAULT_CONFIG
    import copy
    return copy.deepcopy(DEFAULT_CONFIG)

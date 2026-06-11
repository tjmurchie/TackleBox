"""Tests for basic QC and duplicate detection via annotate() (spinner.annotation)."""
from __future__ import annotations

import copy
from pathlib import Path

import pytest

from spinner.config import DEFAULT_CONFIG
from spinner.fasta import parse_fasta
from spinner.annotation import annotate


def make_cfg(**qc_overrides):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["basic_qc"].update(qc_overrides)
    # Disable adapter/keyword scan to isolate QC logic.
    cfg["steps"]["adapter_screen"] = False
    cfg["steps"]["bad_keyword_screen"] = False
    return cfg


def parse_inline(fasta_text: str, tmp_path: Path):
    p = tmp_path / "t.fasta"
    p.write_text(fasta_text, encoding="utf-8")
    return parse_fasta([str(p)])


# ---------------------------------------------------------------------------
# Length filters
# ---------------------------------------------------------------------------

def test_length_below_min(tmp_path):
    recs = parse_inline(">A.1 test\nACGT\n", tmp_path)
    cfg = make_cfg(min_length=50)
    ann = annotate(recs, cfg)
    assert "length_below_min" in ann["A.1"].reasons


def test_length_above_max(tmp_path):
    recs = parse_inline(">A.1 test\n" + "A" * 100 + "\n", tmp_path)
    cfg = make_cfg(max_length=50)
    ann = annotate(recs, cfg)
    assert "length_above_max" in ann["A.1"].reasons


def test_length_within_bounds_no_reason(tmp_path):
    recs = parse_inline(">A.1 test\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg(min_length=50, max_length=200)
    ann = annotate(recs, cfg)
    assert "length_below_min" not in ann["A.1"].reasons
    assert "length_above_max" not in ann["A.1"].reasons


# ---------------------------------------------------------------------------
# N fraction
# ---------------------------------------------------------------------------

def test_n_fraction_high(tmp_path):
    recs = parse_inline(">A.1 test\n" + "N" * 80 + "\n", tmp_path)
    cfg = make_cfg(max_n_fraction=0.05)
    ann = annotate(recs, cfg)
    assert "n_fraction_high" in ann["A.1"].reasons
    assert ann["A.1"].n_fraction == pytest.approx(1.0)


def test_n_fraction_ok(tmp_path):
    recs = parse_inline(">A.1 test\n" + "N" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg(max_n_fraction=0.05)
    ann = annotate(recs, cfg)
    assert "n_fraction_high" not in ann["A.1"].reasons


# ---------------------------------------------------------------------------
# Non-IUPAC fraction
# ---------------------------------------------------------------------------

def test_non_iupac_fraction_high(tmp_path):
    # 'Z' and 'X' are not IUPAC DNA.
    recs = parse_inline(">A.1 test\n" + "Z" * 50 + "ACGT" * 10 + "\n", tmp_path)
    cfg = make_cfg(max_non_iupac_fraction=0.01)
    ann = annotate(recs, cfg)
    assert "non_iupac_fraction_high" in ann["A.1"].reasons


# ---------------------------------------------------------------------------
# Low complexity / Shannon entropy
# ---------------------------------------------------------------------------

def test_low_complexity_all_same_base(tmp_path):
    recs = parse_inline(">A.1 test\n" + "A" * 100 + "\n", tmp_path)
    cfg = make_cfg(min_shannon_entropy=1.2)
    ann = annotate(recs, cfg)
    assert "low_complexity" in ann["A.1"].reasons
    assert ann["A.1"].shannon_entropy == pytest.approx(0.0)


def test_normal_complexity_no_flag(tmp_path):
    recs = parse_inline(">A.1 test\n" + "ACGT" * 25 + "\n", tmp_path)
    cfg = make_cfg(min_shannon_entropy=1.2)
    ann = annotate(recs, cfg)
    assert "low_complexity" not in ann["A.1"].reasons


# ---------------------------------------------------------------------------
# Homopolymer
# ---------------------------------------------------------------------------

def test_homopolymer_long(tmp_path):
    # Use "T" * 70 so the run doesn't merge with surrounding ACGT sequences.
    recs = parse_inline(">A.1 test\n" + "ACGT" * 5 + "T" * 70 + "ACGT" * 5 + "\n", tmp_path)
    cfg = make_cfg(max_homopolymer_run=60)
    ann = annotate(recs, cfg)
    assert "homopolymer_long" in ann["A.1"].reasons
    assert ann["A.1"].max_homopolymer >= 70


def test_homopolymer_ok(tmp_path):
    recs = parse_inline(">A.1 test\n" + "ACGT" * 25 + "\n", tmp_path)
    cfg = make_cfg(max_homopolymer_run=60)
    ann = annotate(recs, cfg)
    assert "homopolymer_long" not in ann["A.1"].reasons


# ---------------------------------------------------------------------------
# Duplicate accession detection
# ---------------------------------------------------------------------------

def test_duplicate_accession_second_flagged(tmp_path):
    content = ">DUP.1 first\nACGT\n>DUP.1 second\nGGGG\n"
    recs = parse_inline(content, tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert len(ann) == 2
    first = ann["DUP.1"]
    second = ann["DUP.1__dup2"]
    assert not first.duplicate_accession
    assert second.duplicate_accession
    assert "duplicate_accession" in second.reasons


def test_duplicate_accession_first_not_flagged(tmp_path):
    content = ">DUP.1 first\nACGT\n>DUP.1 second\nGGGG\n"
    recs = parse_inline(content, tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert not ann["DUP.1"].duplicate_accession


# ---------------------------------------------------------------------------
# Duplicate sequence detection
# ---------------------------------------------------------------------------

def test_duplicate_sequence_different_accession(tmp_path):
    seq = "ACGT" * 20
    content = f">A.1 first\n{seq}\n>B.1 second\n{seq}\n"
    recs = parse_inline(content, tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    # First is not flagged, second is.
    assert not ann["A.1"].duplicate_sequence
    assert ann["B.1"].duplicate_sequence
    assert "duplicate_sequence" in ann["B.1"].reasons


def test_duplicate_sequence_same_accession_both_flagged(tmp_path):
    """Duplicate accession with identical sequences should flag both dup_acc and dup_seq on second."""
    seq = "ACGT" * 10
    content = f">X.1 a\n{seq}\n>X.1 b\n{seq}\n"
    recs = parse_inline(content, tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    second = ann["X.1__dup2"]
    assert second.duplicate_accession
    assert second.duplicate_sequence


# ---------------------------------------------------------------------------
# Per-class minimum length
# ---------------------------------------------------------------------------

def test_length_below_class_min(tmp_path):
    """A NucMark record shorter than min_length_by_class should get length_below_class_min."""
    recs = parse_inline(">A.1 Homo sapiens 18S ribosomal RNA\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg(min_length_by_class={"NucMark": 500})
    ann = annotate(recs, cfg)
    assert "length_below_class_min" in ann["A.1"].reasons


def test_length_below_class_min_other_class_unaffected(tmp_path):
    """A Mito record should not be affected by a NucMark-specific minimum."""
    recs = parse_inline(">A.1 Homo sapiens mitochondrion complete genome\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg(min_length=50, min_length_by_class={"NucMark": 500})
    ann = annotate(recs, cfg)
    assert "length_below_class_min" not in ann["A.1"].reasons


def test_length_above_class_min_no_flag(tmp_path):
    """A NucMark record meeting the per-class minimum should not be flagged."""
    recs = parse_inline(">A.1 Homo sapiens 18S ribosomal RNA\n" + "ACGT" * 200 + "\n", tmp_path)
    cfg = make_cfg(min_length_by_class={"NucMark": 50})
    ann = annotate(recs, cfg)
    assert "length_below_class_min" not in ann["A.1"].reasons


# ---------------------------------------------------------------------------
# Score bonus reasons: refseq_preferred, voucher_keyword, complete_organelle
# ---------------------------------------------------------------------------

def test_refseq_preferred_nc_prefix(tmp_path):
    recs = parse_inline(">NC_001234.1 Homo sapiens mitochondrion\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "refseq_preferred" in ann["NC_001234.1"].reasons


def test_refseq_preferred_nr_prefix(tmp_path):
    recs = parse_inline(">NR_024540.1 Homo sapiens 18S ribosomal RNA\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "refseq_preferred" in ann["NR_024540.1"].reasons


def test_refseq_preferred_not_triggered_for_genbank(tmp_path):
    recs = parse_inline(">MK123456.1 Rangifer tarandus COI\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "refseq_preferred" not in ann["MK123456.1"].reasons


def test_voucher_keyword_in_header(tmp_path):
    recs = parse_inline(">MK123456.1 Rangifer tarandus voucher TMu-001 COI gene\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "voucher_keyword" in ann["MK123456.1"].reasons


def test_voucher_keyword_holotype(tmp_path):
    recs = parse_inline(">MK123456.1 Betula papyrifera holotype rbcL gene\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "voucher_keyword" in ann["MK123456.1"].reasons


def test_voucher_keyword_not_triggered(tmp_path):
    recs = parse_inline(">MK123456.1 Rangifer tarandus cytochrome b gene\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "voucher_keyword" not in ann["MK123456.1"].reasons


def test_complete_organelle_complete_genome(tmp_path):
    recs = parse_inline(">NC_001234.1 Homo sapiens mitochondrion complete genome\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "complete_organelle" in ann["NC_001234.1"].reasons


def test_complete_organelle_not_triggered(tmp_path):
    recs = parse_inline(">MK123456.1 Betula papyrifera rbcL partial sequence\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert "complete_organelle" not in ann["MK123456.1"].reasons


# ---------------------------------------------------------------------------
# All records present in decisions even with duplicates
# ---------------------------------------------------------------------------

def test_all_records_in_annotation(tmp_path):
    content = ">A.1 a\nAAAA\n>A.1 b\nCCCC\n>B.1 c\nGGGG\n"
    recs = parse_inline(content, tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg)
    assert len(ann) == 3
    assert "A.1" in ann
    assert "A.1__dup2" in ann
    assert "B.1" in ann

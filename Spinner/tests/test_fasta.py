"""Tests for FASTA parsing and writing (spinner.fasta)."""
from __future__ import annotations

import gzip
import textwrap
from pathlib import Path

import pytest

from spinner.fasta import FastaRecord, parse_fasta, write_fasta


def write_tmp(tmp_path: Path, content: str) -> Path:
    p = tmp_path / "test.fasta"
    p.write_text(content, encoding="utf-8")
    return p


def test_parse_single_record(tmp_path):
    p = write_tmp(tmp_path, ">ACC1.1 Homo sapiens mitochondrion\nACGTACGT\n")
    recs = parse_fasta([str(p)])
    assert len(recs) == 1
    assert recs[0].accession == "ACC1.1"
    assert recs[0].seq_upper == "ACGTACGT"
    assert recs[0].id == "ACC1.1"


def test_parse_multiline_sequence(tmp_path):
    p = write_tmp(tmp_path, ">ACC2.1 desc\nACGT\nACGT\nACGT\n")
    recs = parse_fasta([str(p)])
    assert recs[0].seq_upper == "ACGTACGTACGT"


def test_parse_lowercase_sequence(tmp_path):
    p = write_tmp(tmp_path, ">ACC3.1 desc\nacgtacgt\n")
    recs = parse_fasta([str(p)])
    assert recs[0].seq_upper == "ACGTACGT"


def test_parse_empty_file(tmp_path):
    p = write_tmp(tmp_path, "")
    recs = parse_fasta([str(p)])
    assert recs == []


def test_parse_multiple_records(tmp_path):
    content = ">A.1 desc\nAAAA\n>B.1 desc\nCCCC\n>C.1 desc\nGGGG\n"
    p = write_tmp(tmp_path, content)
    recs = parse_fasta([str(p)])
    assert len(recs) == 3
    assert [r.accession for r in recs] == ["A.1", "B.1", "C.1"]


def test_duplicate_accession_preserved(tmp_path):
    """Both records with the same accession must appear in the output."""
    content = ">DUP.1 seq A\nAAAA\n>DUP.1 seq B\nCCCC\n"
    p = write_tmp(tmp_path, content)
    recs = parse_fasta([str(p)])
    assert len(recs) == 2
    keys = [r.id for r in recs]
    assert "DUP.1" in keys
    assert "DUP.1__dup2" in keys


def test_duplicate_accession_third_copy(tmp_path):
    content = ">X.1 a\nAAAA\n>X.1 b\nCCCC\n>X.1 c\nGGGG\n"
    p = write_tmp(tmp_path, content)
    recs = parse_fasta([str(p)])
    assert len(recs) == 3
    assert [r.id for r in recs] == ["X.1", "X.1__dup2", "X.1__dup3"]


def test_write_fasta_roundtrip(tmp_path):
    """Records written by write_fasta should be re-parseable with identical sequences."""
    src = write_tmp(tmp_path, ">ACC.1 test\nACGTACGT\n")
    recs = parse_fasta([str(src)])
    dst = tmp_path / "out.fasta"
    write_fasta(recs, str(dst))
    recs2 = parse_fasta([str(dst)])
    assert len(recs2) == 1
    assert recs2[0].seq_upper == "ACGTACGT"


def test_header_preserved_in_write(tmp_path):
    content = ">NC_001234.1 Rangifer tarandus mitochondrion complete genome\nACGT\n"
    p = write_tmp(tmp_path, content)
    recs = parse_fasta([str(p)])
    out = tmp_path / "out.fasta"
    write_fasta(recs, str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == ">NC_001234.1 Rangifer tarandus mitochondrion complete genome"


def test_parse_multiple_files(tmp_path):
    p1 = tmp_path / "a.fasta"
    p2 = tmp_path / "b.fasta"
    p1.write_text(">A.1 desc\nAAAA\n")
    p2.write_text(">B.1 desc\nCCCC\n")
    recs = parse_fasta([str(p1), str(p2)])
    assert len(recs) == 2
    assert {r.accession for r in recs} == {"A.1", "B.1"}


def test_parse_gzipped_fasta(tmp_path):
    """parse_fasta should transparently handle .fasta.gz input."""
    content = ">GZ.1 Homo sapiens mitochondrion\nACGTACGT\n"
    gz_path = tmp_path / "test.fasta.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as f:
        f.write(content)
    recs = parse_fasta([str(gz_path)])
    assert len(recs) == 1
    assert recs[0].accession == "GZ.1"
    assert recs[0].seq_upper == "ACGTACGT"


def test_parse_gzipped_multiple_records(tmp_path):
    """Gzipped FASTA with multiple records should parse all records."""
    content = ">A.1 desc\nACGT\n>B.1 desc\nGGGG\n"
    gz_path = tmp_path / "multi.fasta.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as f:
        f.write(content)
    recs = parse_fasta([str(gz_path)])
    assert len(recs) == 2
    assert {r.accession for r in recs} == {"A.1", "B.1"}

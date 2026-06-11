"""Tests for bad-keyword header scanning (spinner.keywords + annotation)."""
from __future__ import annotations

import copy
from pathlib import Path

import pytest

from spinner.config import DEFAULT_CONFIG
from spinner.fasta import parse_fasta
from spinner.annotation import annotate


def make_cfg():
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["steps"]["adapter_screen"] = False
    cfg["steps"]["basic_qc"] = False
    return cfg


def parse_inline(fasta_text: str, tmp_path: Path):
    p = tmp_path / "t.fasta"
    p.write_text(fasta_text, encoding="utf-8")
    return parse_fasta([str(p)])


def test_keyword_reject(tmp_path, bad_keywords_file):
    recs = parse_inline(">BADVEC.1 Synthetic construct vector sequence\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg, bad_keywords_tsv=str(bad_keywords_file))
    a = ann["BADVEC.1"]
    assert a.bad_keyword_hit
    assert a.bad_keyword_action == "reject"
    assert "bad_keyword_reject" in a.reasons


def test_keyword_review(tmp_path, bad_keywords_file):
    recs = parse_inline(">CLON.1 Some clone from library\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg, bad_keywords_tsv=str(bad_keywords_file))
    a = ann["CLON.1"]
    assert a.bad_keyword_hit
    assert a.bad_keyword_action == "review"
    assert "bad_keyword_review" in a.reasons
    assert "bad_keyword_reject" not in a.reasons


def test_keyword_case_insensitive(tmp_path, bad_keywords_file):
    # Keyword is "UNVERIFIED" (uppercase in TSV), header has lowercase version.
    recs = parse_inline(">BADX.1 unverified sequence from NCBI\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    cfg["bad_keywords"]["case_sensitive"] = False
    ann = annotate(recs, cfg, bad_keywords_tsv=str(bad_keywords_file))
    a = ann["BADX.1"]
    assert a.bad_keyword_hit


def test_keyword_case_sensitive_no_match(tmp_path, bad_keywords_file):
    recs = parse_inline(">BADX.1 unverified sequence\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    cfg["bad_keywords"]["case_sensitive"] = True  # "UNVERIFIED" vs "unverified" -> no match
    ann = annotate(recs, cfg, bad_keywords_tsv=str(bad_keywords_file))
    a = ann["BADX.1"]
    assert not a.bad_keyword_hit


def test_clean_record_no_keyword(tmp_path, bad_keywords_file):
    recs = parse_inline(">CLEAN.1 Rangifer tarandus mitochondrion complete genome\n" + "ACGT" * 30 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg, bad_keywords_tsv=str(bad_keywords_file))
    assert not ann["CLEAN.1"].bad_keyword_hit


def test_keyword_first_match_only(tmp_path, bad_keywords_file):
    """Only the first matching keyword should be reported."""
    # Header contains both "synthetic construct" and "vector".
    recs = parse_inline(">X.1 synthetic construct vector pUC19\n" + "ACGT" * 20 + "\n", tmp_path)
    cfg = make_cfg()
    ann = annotate(recs, cfg, bad_keywords_tsv=str(bad_keywords_file))
    a = ann["X.1"]
    # Only one keyword should be recorded.
    kw_reasons = [r for r in a.reasons if r.startswith("bad_keyword")]
    assert len(kw_reasons) == 1

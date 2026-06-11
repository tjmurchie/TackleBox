"""Tests for adapter/primer contamination scanning (spinner.adapters)."""
from __future__ import annotations

import copy

import pytest

from spinner.adapters import Adapter, find_adapter, hamming_leq, is_terminal, load_adapters
from spinner.config import DEFAULT_CONFIG


def make_cfg(**adapter_overrides):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["adapter_screen"].update(adapter_overrides)
    return cfg


# ---------------------------------------------------------------------------
# hamming_leq
# ---------------------------------------------------------------------------

def test_hamming_exact_match():
    assert hamming_leq("AAAA", "AAAA", 0) is True


def test_hamming_one_mismatch_allowed():
    assert hamming_leq("AAAA", "AATA", 1) is True


def test_hamming_one_mismatch_not_allowed():
    assert hamming_leq("AAAA", "AATA", 0) is False


def test_hamming_exceeds_max():
    assert hamming_leq("AAAA", "TTTT", 2) is False


def test_hamming_all_match_zero_allowed():
    assert hamming_leq("GCGC", "GCGC", 0) is True


# ---------------------------------------------------------------------------
# is_terminal
# ---------------------------------------------------------------------------

def test_is_terminal_at_left_edge():
    # Position 0, match length 10, seq length 100, window 25 -> terminal
    assert is_terminal(0, 10, 100, 25) is True


def test_is_terminal_within_left_window():
    # Position 20, match 10, seq 100, window 25 -> position 20 <= 25 -> terminal
    assert is_terminal(20, 10, 100, 25) is True


def test_is_terminal_at_right_edge():
    # Position 80, match 20, seq 100, window 25 -> 80+20=100 >= 75 -> terminal
    assert is_terminal(80, 20, 100, 25) is True


def test_is_internal():
    # Position 40, match 10, seq 200, window 25 -> 40 > 25 and 50 < 175 -> internal
    assert is_terminal(40, 10, 200, 25) is False


# ---------------------------------------------------------------------------
# find_adapter
# ---------------------------------------------------------------------------

TRUSEQ = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
AD_TRUSEQ = Adapter("TruSeq", TRUSEQ, max_mismatch=0, action="reject")
AD_TRUSEQ_MM2 = Adapter("TruSeq_mm2", TRUSEQ, max_mismatch=2, action="reject")


def test_find_adapter_exact_internal():
    # Adapter placed in the middle of a long sequence -> internal hit.
    seq = "A" * 80 + TRUSEQ + "A" * 80
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=False)
    result = find_adapter(seq, AD_TRUSEQ, cfg)
    assert result is not None
    pos, strand = result
    assert pos == 80
    assert strand == "+"


def test_find_adapter_exact_terminal():
    # Adapter at position 2 (within 25 bp window) -> terminal.
    seq = "AC" + TRUSEQ + "A" * 80
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=False)
    result = find_adapter(seq, AD_TRUSEQ, cfg)
    assert result is not None
    pos, strand = result
    assert pos == 2
    assert is_terminal(pos, len(TRUSEQ), len(seq), 25) is True


def test_find_adapter_with_mismatch():
    # One mismatch in the adapter sequence.
    mutated = TRUSEQ[:5] + "X" + TRUSEQ[6:]
    seq = "A" * 50 + mutated + "A" * 50
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=False)
    result = find_adapter(seq, AD_TRUSEQ_MM2, cfg)
    assert result is not None


def test_find_adapter_no_match():
    seq = "A" * 200
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=False)
    result = find_adapter(seq, AD_TRUSEQ, cfg)
    assert result is None


def test_find_adapter_revcomp():
    """Reverse complement of adapter should be detected on minus strand."""
    from spinner.seq_utils import revcomp
    rc = revcomp(TRUSEQ)
    seq = "A" * 50 + rc + "A" * 50
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=True)
    result = find_adapter(seq, AD_TRUSEQ, cfg)
    assert result is not None
    pos, strand = result
    assert strand == "-"


def test_find_adapter_too_short_adapter():
    """Adapter shorter than min_adapter_match_length should not match."""
    short_ad = Adapter("short", "AGATCG", max_mismatch=0)
    seq = "ACGT" * 20 + "AGATCG" + "ACGT" * 20
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=False)
    result = find_adapter(seq, short_ad, cfg)
    assert result is None


def test_find_adapter_seq_shorter_than_adapter():
    seq = "AGATCGG"  # shorter than TRUSEQ
    cfg = make_cfg(min_adapter_match_length=12, scan_reverse_complement=False)
    result = find_adapter(seq, AD_TRUSEQ, cfg)
    assert result is None


# ---------------------------------------------------------------------------
# load_adapters
# ---------------------------------------------------------------------------

def test_load_adapters(adapters_file):
    ads = load_adapters(str(adapters_file))
    assert len(ads) >= 1
    names = [a.name for a in ads]
    assert "Illumina_TruSeq_Universal" in names


def test_load_adapters_missing_file():
    ads = load_adapters("/no/such/file.tsv")
    assert ads == []


def test_load_adapters_empty_path():
    ads = load_adapters("")
    assert ads == []

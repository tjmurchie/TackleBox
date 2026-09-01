"""Tests for self_repeat_softmask() and the --skip-self-mask CLI default.

Real finding, 2026-09-01, from a live ~18,500-reference multi-species PalaeoSCOPE
Phase B panel: self_repeat_softmask() counts k-mers ACROSS THE ENTIRE input FASTA, not
per-reference -- correct for the small/few-reference case this was originally built for
(detecting genuine within-genome structural repeats), but at real multi-species panel
scale it stops detecting genuine repeats and instead flags ordinary, biologically
CONSERVED single-copy coding sequence -- shared across species precisely because it's a
good marker -- as "repetitive". Confirmed on a real 18,541-reference panel: 93.7% of all
bases masked overall, and a real 651bp COI barcode gene masked 100%. `--skip-self-mask`
now defaults to True (masking skipped) instead of False; `--no-skip-self-mask` restores
the original masked-by-default behavior for anyone who still wants it (a genuinely
small/few-genome reference set is exactly where within-genome repeat detection is still
the right call). This module previously had zero test coverage of any kind.
"""
from __future__ import annotations

import FlyForge as ff


def _make_parser():
    return ff.build_arg_parser()


def _parse(parser, extra_args):
    return parser.parse_args(["-i", "x.fa", "--prefix", "x", *extra_args])


class TestSkipSelfMaskDefault:
    def test_bare_invocation_defaults_to_skipped(self):
        args = _parse(_make_parser(), [])
        assert args.skip_self_mask is True

    def test_explicit_skip_flag_still_works(self):
        """Backward compatible: anyone already passing --skip-self-mask explicitly gets
        the exact same True value as before this default changed."""
        args = _parse(_make_parser(), ["--skip-self-mask"])
        assert args.skip_self_mask is True

    def test_explicit_no_skip_restores_old_masked_default(self):
        args = _parse(_make_parser(), ["--no-skip-self-mask"])
        assert args.skip_self_mask is False

    def test_repeat_k_and_threshold_defaults_unchanged(self):
        """Only the skip/apply default changed -- the masking algorithm's own parameters,
        used whenever masking IS enabled via --no-skip-self-mask, are untouched."""
        args = _parse(_make_parser(), ["--no-skip-self-mask"])
        assert args.repeat_k == 15
        assert args.repeat_threshold == 3


class TestSelfRepeatSoftmask:
    """Pure-function correctness of self_repeat_softmask() itself -- unchanged by the
    CLI default flip, but had zero prior test coverage."""

    def test_kmer_below_threshold_is_not_masked(self):
        # A 15-mer appearing only twice (threshold defaults to 3) must survive unmasked.
        seq = "ACGTACGTACGTACG" * 2  # 30bp, repeats one 15-mer exactly twice
        masked = ff.self_repeat_softmask({"a": seq, "b": "TTTTTTTTTTTTTTT" * 2}, k=15, threshold=3)
        assert masked["a"] == seq.upper()  # unchanged, no lowercase introduced

    def test_kmer_at_or_above_threshold_is_masked(self):
        shared_kmer = "ACGTACGTACGTACG"  # 15bp
        # Non-repetitive flanks (no internal homopolymer/low-complexity runs, and no
        # accidental cross-sequence collisions with each other or the shared prefix) --
        # isolates the shared-prefix signal from any other source of "repeat" count.
        seqs = {
            "a": shared_kmer + "CAGATTTTCATATTATGCAG",
            "b": shared_kmer + "AAAGCGGCACTTGTGAAGTG",
            "c": shared_kmer + "CCGTAATGCCTTTCCCTAAC",
        }
        masked = ff.self_repeat_softmask(seqs, k=15, threshold=3)
        # The shared 15bp prefix (appearing in all 3 sequences -- at/above threshold)
        # must be lowercased in every sequence that contains it.
        for sid in seqs:
            assert masked[sid][:15] == shared_kmer.lower()
            assert masked[sid][15:] == seqs[sid][15:].upper()  # rest untouched

    def test_kmers_containing_n_are_never_counted_or_masked(self):
        seqs = {
            "a": "ACGTNCGTACGTACG" + "CAGATTTTCATATTATGCAG",
            "b": "ACGTNCGTACGTACG" + "AAAGCGGCACTTGTGAAGTG",
            "c": "ACGTNCGTACGTACG" + "CCGTAATGCCTTTCCCTAAC",
        }
        masked = ff.self_repeat_softmask(seqs, k=15, threshold=3)
        # Despite appearing 3x, every window overlapping the N stays uncounted/unmasked
        # -- the N-containing prefix must survive entirely unmasked (still uppercase).
        for sid in seqs:
            assert masked[sid][:15] == seqs[sid][:15].upper()

    def test_case_insensitive_counting_real_world_mixed_case_input(self):
        """Real GenBank-derived FASTA can contain mixed-case sequence -- counting must
        treat 'acgt' and 'ACGT' as the same k-mer (the function upper()s internally)."""
        shared_kmer_upper = "ACGTACGTACGTACG"
        seqs = {
            "a": shared_kmer_upper + "CAGATTTTCATATTATGCAG",
            "b": shared_kmer_upper.lower() + "AAAGCGGCACTTGTGAAGTG",  # same k-mer, lowercase input
            "c": shared_kmer_upper + "CCGTAATGCCTTTCCCTAAC",
        }
        masked = ff.self_repeat_softmask(seqs, k=15, threshold=3)
        for sid in seqs:
            assert masked[sid][:15] == shared_kmer_upper.lower()
            assert masked[sid][15:] == seqs[sid][15:].upper()  # flank untouched, isolates the claim

    def test_empty_input_returns_empty_dict(self):
        assert ff.self_repeat_softmask({}, k=15, threshold=3) == {}

    def test_real_conserved_marker_gene_gets_masked_across_many_species(self):
        """Directly demonstrates the real finding this whole change is about: a single
        conserved region shared across many otherwise-distinct sequences (standing in
        for a real conserved coding gene shared across many species) gets masked in
        EVERY one of them, not retained as a single representative -- this is exactly
        why --skip-self-mask now defaults to True for multi-species panel use."""
        import random
        conserved = "ATGGCCATTGTAATG"  # 15bp, a stand-in for a conserved codon stretch
        # 20 "species", each with its own genuinely distinct, non-repetitive divergent
        # flank (not a homopolymer -- a run of one repeated base would itself trigger
        # real repeat detection independent of the conserved locus, confounding the
        # claim this test is isolating) but the SAME conserved 15-mer.
        seqs = {
            f"species_{i}": conserved + "".join(random.Random(i).choice("ACGT") for _ in range(20))
            for i in range(20)
        }
        masked = ff.self_repeat_softmask(seqs, k=15, threshold=3)
        n_masked_at_locus = sum(1 for s in masked.values() if s[:15] == conserved.lower())
        assert n_masked_at_locus == 20  # every species loses this locus, not just the "extra" copies


class TestFlyForgeAuditSkipSelfMaskDefault:
    """FlyForgeAudit.py is a separate tool but reuses FlyForge's own
    self_repeat_softmask() and had its own independent --skip-self-mask flag/default --
    fixed for consistency alongside FlyForge's own default change, via the same shared
    add_shared_filtering_args() helper both the `audit` and `augment` subcommands use."""

    def _parser(self):
        import argparse
        import FlyForgeAudit as ffa
        p = argparse.ArgumentParser()
        ffa.add_shared_filtering_args(p)
        return p

    def test_bare_invocation_defaults_to_skipped(self):
        assert self._parser().parse_args([]).skip_self_mask is True

    def test_explicit_skip_flag_still_works(self):
        assert self._parser().parse_args(["--skip-self-mask"]).skip_self_mask is True

    def test_explicit_no_skip_restores_old_masked_default(self):
        assert self._parser().parse_args(["--no-skip-self-mask"]).skip_self_mask is False

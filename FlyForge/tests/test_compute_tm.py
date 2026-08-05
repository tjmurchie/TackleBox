"""Tests for FlyForge.compute_tm() -- the melting-temperature calculation
that feeds filter_melting_temp()'s pass/fail decision for every bait.

compute_tm() wraps Bio.SeqUtils.MeltingTemp.Tm_NN in a bare
``except Exception: return 0.0``. That is a real, silent failure mode: a
computation error (empty/degenerate sequence, unexpected characters) produces
the exact same output (0.0) as a bait that was measured and genuinely has a
very low Tm. Downstream, filter_melting_temp's default min_tm=50.0 then
discards both cases identically -- a bait can vanish from the design because
its Tm calculation failed, not because its Tm is actually low, with no error
raised anywhere. These tests lock in and document that real behavior (not
propose changing it) so a future refactor can't silently make it worse
(e.g. by letting a raised exception propagate and crash a multi-hour run
over one bad reference sequence), and so anyone reading the test suite sees
the risk explicitly rather than having to rediscover it by reading the
source.
"""

import pytest

import FlyForge as ff


class TestComputeTmNormalInputs:
    def test_real_sequence_returns_plausible_positive_tm(self):
        # A real, non-degenerate 80 nt sequence -- physically, nearest-
        # neighbor Tm for a fragment this long and GC-balanced should land
        # comfortably in a normal PCR/hybridization range. Checking a sanity
        # RANGE rather than an exact literal, since the exact value is
        # Bio.SeqUtils.MeltingTemp's own internal computation, not something
        # this test should re-derive by hand.
        seq = ("ACGT" * 20)[:80]
        assert len(seq) == 80
        tm = ff.compute_tm(seq)
        assert 40.0 < tm < 100.0

    def test_gc_rich_sequence_has_higher_tm_than_at_rich(self):
        # A basic physical-plausibility check on the real thermodynamic
        # model, independent of the exact numbers: more G/C pairing (3
        # hydrogen bonds) should raise Tm relative to an A/T-heavy sequence
        # of the same length (2 hydrogen bonds).
        at_rich = ("AT" * 40)[:80]
        gc_rich = ("GC" * 40)[:80]
        assert len(at_rich) == len(gc_rich) == 80
        assert ff.compute_tm(gc_rich) > ff.compute_tm(at_rich)

    def test_lowercase_input_handled_same_as_uppercase(self):
        seq_upper = "ACGTACGATCGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGATCGATCGTAGCTA"
        seq_lower = seq_upper.lower()
        assert ff.compute_tm(seq_upper) == ff.compute_tm(seq_lower)


class TestComputeTmFailureIsSilent:
    """Confirms real inputs that make Bio.SeqUtils.MeltingTemp.Tm_NN raise
    (IndexError, verified directly against biopython) are swallowed and
    turned into 0.0 rather than propagating -- the documented, real
    consequence flagged to Tyler: this looks identical to "measured, and
    it's genuinely very low" to every downstream consumer."""

    @pytest.mark.parametrize("bad_seq", [
        "",                     # empty string
        "N" * 80,               # all-ambiguous, no real bases for the NN model
        "X" * 80,               # not a valid base at all
    ])
    def test_pathological_sequences_return_zero_not_raise(self, bad_seq):
        assert ff.compute_tm(bad_seq) == 0.0

    def test_zero_tm_is_indistinguishable_from_a_real_low_tm_bait(self):
        """Documents the concrete downstream consequence: a Tm-computation
        failure and a real (hypothetical) 0-degree bait both fail
        filter_melting_temp's default cutoff identically -- there is no
        signal anywhere in the output that distinguishes "this bait's Tm
        could not be computed" from "this bait was measured and is bad"."""
        failed_bait = _make_bait("ref|fail", "N" * 80, tm=ff.compute_tm("N" * 80))
        assert failed_bait.tm == 0.0

        kept, dropped_count = ff.filter_melting_temp([failed_bait], min_tm=50.0)
        assert kept == []
        assert dropped_count == 1


def _make_bait(bait_id: str, seq: str, tm: float) -> "ff.Bait":
    return ff.Bait(
        bait_id=bait_id, seq=seq, ref_id="ref", ref_start=1, ref_end=len(seq),
        gc_frac=0.5, masked_frac=0.0, ambiguous_count=0, tm=tm,
    )

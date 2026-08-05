"""Tests for FlyForge.compute_tm() -- the melting-temperature calculation
that feeds filter_melting_temp()'s pass/fail decision for every bait.

compute_tm() wraps Bio.SeqUtils.MeltingTemp.Tm_NN in a bare
``except Exception: return 0.0``. That is still the same value a genuinely
very-low-Tm bait would have -- a Tm-calculation failure and a real low-Tm
bait are excluded identically by filter_melting_temp's default min_tm=50.0,
and that filtering behavior is deliberately UNCHANGED (a failed bait should
still not be trusted for synthesis). What changed is visibility: compute_tm()
now accepts an optional ``on_failure`` callback, and run_pipeline() uses it
to count failures across a whole run and print an explicit warning if any
occurred -- so the difference between "this bait's Tm couldn't be computed"
and "this bait was measured and is bad" is no longer invisible, without
raising and aborting a run over one bad reference sequence.
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


class TestComputeTmFailureIsSilentToTheFilter:
    """Confirms real inputs that make Bio.SeqUtils.MeltingTemp.Tm_NN raise
    (IndexError, verified directly against biopython) are swallowed and
    turned into 0.0 rather than propagating -- filter_melting_temp still
    can't tell a failure from a genuinely bad bait, by design (see module
    docstring). What CAN tell them apart now is on_failure (tested below)."""

    @pytest.mark.parametrize("bad_seq", [
        "",                     # empty string
        "N" * 80,               # all-ambiguous, no real bases for the NN model
        "X" * 80,               # not a valid base at all
    ])
    def test_pathological_sequences_return_zero_not_raise(self, bad_seq):
        assert ff.compute_tm(bad_seq) == 0.0

    def test_zero_tm_is_indistinguishable_from_a_real_low_tm_bait_to_the_filter(self):
        """Documents the concrete downstream consequence: a Tm-computation
        failure and a real (hypothetical) 0-degree bait both fail
        filter_melting_temp's default cutoff identically -- filter_melting_
        temp itself has no way to tell "couldn't be computed" apart from
        "measured and is bad" (that's what on_failure is for instead)."""
        failed_bait = _make_bait("ref|fail", "N" * 80, tm=ff.compute_tm("N" * 80))
        assert failed_bait.tm == 0.0

        kept, dropped_count = ff.filter_melting_temp([failed_bait], min_tm=50.0)
        assert kept == []
        assert dropped_count == 1


class TestComputeTmFailureVisibility:
    """compute_tm()'s on_failure callback -- how run_pipeline() actually
    surfaces these failures as a count + warning instead of leaving them
    invisible in the output."""

    def test_on_failure_called_once_per_failed_sequence(self):
        calls = []
        assert ff.compute_tm("N" * 80, on_failure=lambda: calls.append(1)) == 0.0
        assert len(calls) == 1

    def test_on_failure_not_called_for_a_successful_calculation(self):
        calls = []
        seq = ("ACGT" * 20)[:80]
        tm = ff.compute_tm(seq, on_failure=lambda: calls.append(1))
        assert tm != 0.0
        assert calls == []

    def test_on_failure_none_is_a_safe_default(self):
        """Backward compatible: existing callers that don't pass on_failure
        at all must keep working exactly as before."""
        assert ff.compute_tm("N" * 80) == 0.0

    def test_tile_sequence_forwards_on_tm_failure_for_every_bait(self):
        """tile_sequence() calls compute_tm() at two different call sites
        (the short-sequence padding branch and the main tiling loop) --
        confirm both are wired to the same counter, using a sequence long
        enough to produce several tiled windows."""
        calls = []
        seq = "N" * 300  # every window will fail Tm calculation
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0,
                                  on_tm_failure=lambda: calls.append(1))
        assert len(baits) > 1
        assert len(calls) == len(baits)

    def test_tile_sequence_short_sequence_padding_branch_forwards_too(self, monkeypatch):
        """The pad_min branch always appends real "T" bases to reach
        bait_len, which in practice is enough for Tm_NN to compute SOMETHING
        (even a physically implausible negative Tm) rather than raise --
        confirmed directly, an all-N short sequence does NOT reach compute_tm's
        except branch once padded. So this specifically tests the WIRING
        (does tile_sequence forward on_tm_failure through this call site at
        all) via a forced failure, rather than relying on finding real
        degenerate content that happens to make biopython raise here."""
        def fake_compute_tm(seq, on_failure=None):
            if on_failure is not None:
                on_failure()
            return 0.0

        monkeypatch.setattr(ff, "compute_tm", fake_compute_tm)
        calls = []
        seq = "N" * 71  # hits the pad_min branch (single padded bait)
        baits = ff.tile_sequence("ref", seq, bait_len=80,
                                  on_tm_failure=lambda: calls.append(1))
        assert len(baits) == 1
        assert len(calls) == 1

    def test_tile_sequence_no_failures_when_all_windows_are_valid(self):
        calls = []
        seq = ("ACGT" * 100)[:300]
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0,
                                  on_tm_failure=lambda: calls.append(1))
        assert len(baits) > 1
        assert calls == []


def _make_bait(bait_id: str, seq: str, tm: float) -> "ff.Bait":
    return ff.Bait(
        bait_id=bait_id, seq=seq, ref_id="ref", ref_start=1, ref_end=len(seq),
        gc_frac=0.5, masked_frac=0.0, ambiguous_count=0, tm=tm,
    )

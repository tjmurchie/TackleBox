"""Tests for FlyForge.tile_sequence() -- the tiling window/step arithmetic,
including the circular-genome wraparound case used for mitochondrial and
other circular targets.

This is classic off-by-one-prone interval math: computing a step size from a
density, generating start positions, guaranteeing full 3' coverage even when
the step doesn't evenly divide the sequence, wrapping windows across the
origin for circular references, and handling several short-sequence edge
cases (padding vs. dropping) at specific length thresholds. None of it had
any regression protection before this suite.
"""

import random

import pytest

import FlyForge as ff


def _rand_seq(n: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestShortSequenceThresholds:
    """FlyForge's default thresholds are omit_short_leq=70, pad_min=71 --
    i.e. a sequence must be AT LEAST 71 nt to produce any bait at all, even
    though the "too short" cutoff is 70. There's no accidental silent-drop
    gap at these specific defaults (they're adjacent), but the boundary
    itself is exactly the kind of off-by-one detail worth pinning down."""

    def test_length_at_or_below_omit_short_leq_returns_nothing(self):
        seq = _rand_seq(70)
        assert ff.tile_sequence("ref", seq, bait_len=80) == []

    def test_length_at_pad_min_is_padded_to_bait_len(self):
        seq = _rand_seq(71)
        baits = ff.tile_sequence("ref", seq, bait_len=80)
        assert len(baits) == 1
        bait = baits[0]
        assert len(bait.seq) == 80
        # Padding is T, appended after the real sequence.
        assert bait.seq[:71] == seq
        assert bait.seq[71:] == "T" * 9
        # ref_end reflects the REAL sequence length, not the padded length --
        # a subtle distinction: reporting ref_end=80 here would silently
        # claim reference support for 9 nt that don't exist in the source.
        assert bait.ref_start == 1
        assert bait.ref_end == 71

    def test_length_between_pad_min_and_bait_len_still_single_padded_bait(self):
        seq = _rand_seq(75)
        baits = ff.tile_sequence("ref", seq, bait_len=80)
        assert len(baits) == 1
        assert len(baits[0].seq) == 80
        assert baits[0].ref_end == 75

    def test_length_just_below_pad_min_with_custom_thresholds_is_dropped(self):
        """With NON-default thresholds (a real gap between omit_short_leq and
        pad_min), a sequence longer than the "too short" cutoff but shorter
        than pad_min is silently dropped rather than padded -- documenting a
        real, if narrow, coverage-loss edge case."""
        seq = _rand_seq(65)
        baits = ff.tile_sequence(
            "ref", seq, bait_len=80, omit_short_leq=50, pad_min=71,
        )
        assert baits == []


class TestLinearTilingCoverage:
    def test_step_matches_bait_len_over_density(self):
        seq = _rand_seq(500)
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=4.0)
        expected_step = round(80 / 4.0)
        starts = [b.ref_start - 1 for b in baits]
        # Every consecutive gap should be the computed step, EXCEPT possibly
        # the final one (which is forced to align with the sequence end).
        gaps = [starts[i + 1] - starts[i] for i in range(len(starts) - 1)]
        assert all(g == expected_step for g in gaps[:-1])

    def test_final_window_always_reaches_the_last_base(self):
        """The tiling step (bait_len / density) will rarely divide evenly
        into (L - bait_len) -- the code explicitly appends a final window
        aligned to the true end so 3' sequence is never left uncovered.
        This is exactly the kind of guarantee that silently breaks under an
        unrelated refactor without a test."""
        for length in (137, 250, 401, 999):
            seq = _rand_seq(length, seed=length)
            baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0)
            last_bait = max(baits, key=lambda b: b.ref_end)
            assert last_bait.ref_end == length, f"length={length}"

    def test_no_window_exceeds_sequence_bounds(self):
        seq = _rand_seq(333)
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0)
        for b in baits:
            assert 1 <= b.ref_start
            assert b.ref_end <= 333
            assert b.ref_end - b.ref_start + 1 == 80
            assert len(b.seq) == 80

    def test_bait_sequence_matches_source_at_reported_coordinates(self):
        seq = _rand_seq(400)
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0)
        for b in baits:
            assert b.seq == seq[b.ref_start - 1: b.ref_end]


class TestCircularWraparound:
    def test_a_window_actually_wraps_across_the_origin(self):
        """With a short-ish circular sequence and default tiling density,
        at least one generated window should start near the end and wrap
        to the beginning -- otherwise the circular-handling code path isn't
        being exercised at all."""
        seq = _rand_seq(200, seed=7)
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0, circular=True)
        wrapping = [b for b in baits if b.ref_start + 80 - 1 > 200]
        assert wrapping, "expected at least one wrapped window for a 200 nt circular sequence with bait_len=80"
        wrapped = wrapping[0]
        wrap_len = (wrapped.ref_start + 80 - 1) - 200
        # The wrapped fragment must be the tail of the sequence followed by
        # its own head -- exactly what a circular molecule looks like when
        # linearized at a different origin.
        assert wrapped.seq == seq[wrapped.ref_start - 1:] + seq[:wrap_len]
        assert len(wrapped.seq) == 80

    def test_circular_mode_covers_every_position_including_the_origin(self):
        """The whole point of circular mode (per the function's own
        docstring) is that origin-adjacent bases aren't under-enriched the
        way they would be if the molecule were treated as linear. Confirm
        every position 0..L-1 is covered by at least one 80 nt window."""
        length = 200
        seq = _rand_seq(length, seed=11)
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0, circular=True)

        covered = [False] * length
        for b in baits:
            start0 = b.ref_start - 1
            for offset in range(80):
                covered[(start0 + offset) % length] = True
        assert all(covered), "circular tiling left at least one position (possibly at the origin) uncovered"

    def test_non_circular_tiling_of_the_same_sequence_does_not_wrap(self):
        seq = _rand_seq(200, seed=7)
        baits = ff.tile_sequence("ref", seq, bait_len=80, tiling_density=3.0, circular=False)
        for b in baits:
            assert b.ref_end <= 200

"""Tests for FlyForge.parse_circular_id_set() -- resolving --circular-ids hints against
the actual set of reference IDs present after tiling input is finalized.

Real bug, found 2026-08-31 via a live ~10,000-species PalaeoSCOPE Phase B run:
`--remove-complements` can chop a single input record into 2+ linear sub-fragments (each
renamed with a `_[start:end]` suffix) when it hits its own self-BLAST minus-strand
region -- if that same accession was also requested as `--circular-ids`, resolving the
bare accession against the now-fragmented ID set is genuinely ambiguous (it matches
multiple fragments), and this used to raise RuntimeError, aborting the entire FlyForge
run and losing every other reference's already-completed tiling work along with it.
"""

import FlyForge as ff


def test_exact_id_match_resolves_normally():
    resolved = ff.parse_circular_id_set("NC_011818.1", ["NC_011818.1", "OTHER.1"])
    assert resolved == {"NC_011818.1"}


def test_single_chopped_fragment_still_resolves_via_fuzzy_contains_match():
    """Exactly one fragment survived (or the accession was merely renamed, not split) --
    still unambiguous, so this should keep working exactly as before."""
    resolved = ff.parse_circular_id_set(
        "NC_011818.1", ["NC_011818.1_[0:16775]", "OTHER.1"],
    )
    assert resolved == {"NC_011818.1_[0:16775]"}


def test_ambiguous_multi_fragment_hint_is_skipped_not_raised():
    """The real bug case: --remove-complements split NC_011818.1 into two fragments, so
    the bare accession hint now matches BOTH -- ambiguous. Must degrade to 'skip this
    hint' rather than crash the whole run."""
    available = ["NC_011818.1_[0:5000]", "NC_011818.1_[5000:16775]", "OTHER.1"]
    resolved = ff.parse_circular_id_set("NC_011818.1", available)
    assert resolved == set()


def test_unresolvable_hint_is_skipped_not_raised():
    resolved = ff.parse_circular_id_set("NOT_PRESENT_AT_ALL.1", ["OTHER.1", "OTHER.2"])
    assert resolved == set()


def test_ambiguous_and_unresolvable_hints_are_logged_as_warnings():
    warnings = []
    available = ["NC_011818.1_[0:5000]", "NC_011818.1_[5000:16775]", "OTHER.1"]
    resolved = ff.parse_circular_id_set("NC_011818.1,NOT_PRESENT.1", available, log_fn=warnings.append)
    assert resolved == set()
    assert len(warnings) == 2
    assert all("could not resolve circular reference identifier" in w for w in warnings)


def test_one_ambiguous_hint_does_not_prevent_a_different_resolvable_hint_in_the_same_call():
    available = ["NC_011818.1_[0:5000]", "NC_011818.1_[5000:16775]", "GOOD.1"]
    resolved = ff.parse_circular_id_set("NC_011818.1,GOOD.1", available)
    assert resolved == {"GOOD.1"}


def test_empty_arg_returns_empty_set_without_calling_log_fn():
    calls = []
    resolved = ff.parse_circular_id_set(None, ["OTHER.1"], log_fn=calls.append)
    assert resolved == set()
    assert calls == []


def test_log_fn_defaults_to_a_silent_noop_when_omitted():
    # Must not raise even though the hint is unresolvable and no log_fn was supplied.
    # (2+ available targets -- with only one candidate overall, resolve_ref_hint_to_target
    # has its own separate "just use the only option" fallback, unrelated to this path.)
    resolved = ff.parse_circular_id_set("MISSING.1", ["OTHER.1", "OTHER.2"])
    assert resolved == set()

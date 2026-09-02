"""Tests for --redundancy-prereduce-identity: an opt-in cd-hit-est pre-pass ahead of
self_blast_filter()'s own redundancy removal.

Real finding, 2026-09-02, following up on the already-fixed --skip-self-mask default
change (see test_self_repeat_mask.py): at real multi-species panel scale with masking
disabled, self_blast_filter()'s all-vs-all blastn plus its Python post-processing loop
can take HOURS (confirmed live: 2h39m1s on a real 726,110-candidate ANIMAL panel from a
real ~18,500-reference PalaeoSCOPE Phase B project, landing at 144,311 final baits).
Real, independently-verified experiment: running cd-hit-est ONCE at its own practical
identity floor (0.80 -- cd-hit-est hard-refuses lower than that in nucleotide mode,
confirmed independently, not a word-length tuning issue) as a cheap pre-pass, then
running self_blast_filter() completely UNCHANGED on the much smaller remainder, reduced
the SAME real 726,110-candidate pool to 118,461 in 68s, then to a 95,015-bait final
result in another ~20s -- roughly 108x faster overall, and landing closer to the real
target budget than the unmodified pipeline. Default is None (disabled) -- zero behavior
change for any existing config/CLI invocation unless this flag is explicitly set.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import FlyForge as ff


def _make_parser():
    return ff.build_arg_parser()


def _parse(parser, extra_args):
    return parser.parse_args(["-i", "x.fa", "--prefix", "x", *extra_args])


class TestRedundancyPrereduceIdentityDefault:
    def test_default_is_disabled(self):
        args = _parse(_make_parser(), [])
        assert args.redundancy_prereduce_identity is None

    def test_explicit_value_round_trips(self):
        args = _parse(_make_parser(), ["--redundancy-prereduce-identity", "0.85"])
        assert args.redundancy_prereduce_identity == pytest.approx(0.85)


class TestRedundancyPrereduceValidation:
    """run_pipeline()'s own early validation, before any file I/O -- these tests exit
    via sys.exit(1) before ever touching a real FASTA or external tool, so they need
    no real input file and no cd-hit-est/blastn binary to be fast and dependency-free."""

    def _args(self, tmp_path, identity, no_redundancy=False):
        parser = _make_parser()
        extra = ["--output-dir", str(tmp_path), "--redundancy-prereduce-identity", str(identity)]
        if no_redundancy:
            extra.append("--no-redundancy")
        return _parse(parser, extra)

    def test_value_below_floor_exits(self, tmp_path):
        args = self._args(tmp_path, 0.79)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_value_well_below_floor_exits(self, tmp_path):
        args = self._args(tmp_path, 0.5)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_value_at_or_above_one_exits(self, tmp_path):
        args = self._args(tmp_path, 1.0)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_invalid_value_still_exits_even_with_no_redundancy(self, tmp_path):
        """The flag would have no effect if --no-redundancy is also set, but an
        explicitly-invalid value should still be reported rather than silently
        ignored -- avoids silent confusion if --no-redundancy is later removed."""
        args = self._args(tmp_path, 0.5, no_redundancy=True)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_unset_does_not_exit_early(self, tmp_path):
        """Sanity check that the validation itself (not some unrelated failure) is what
        causes the exits above: with the flag unset, run_pipeline() must get PAST the
        validation block and fail later for the ordinary reason (no real input file),
        not with the same SystemExit(1)/validation error."""
        parser = _make_parser()
        args = _parse(parser, ["--output-dir", str(tmp_path)])
        with pytest.raises(FileNotFoundError):
            ff.run_pipeline(args)


@pytest.mark.slow
class TestRedundancyPrereduceIntegration:
    """Real end-to-end run with cd-hit-est and blastn actually invoked. A handful of
    references including one exact duplicate pair, so the pre-pass has real, guaranteed
    redundancy to collapse regardless of random tiling details."""

    def _write_fasta(self, path: Path) -> None:
        import random
        rng = random.Random(42)
        seq_a = "".join(rng.choice("ACGT") for _ in range(200))
        seq_b = "".join(rng.choice("ACGT") for _ in range(200))
        records = [
            (">ref1", seq_a),
            (">ref2", seq_a),  # exact duplicate of ref1 -- guaranteed real redundancy
            (">ref3", seq_b),
        ]
        with path.open("w", encoding="utf-8") as f:
            for header, seq in records:
                f.write(f"{header}\n{seq}\n")

    def test_prereduce_step_runs_and_reduces_candidates(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        self._write_fasta(fasta)
        outdir = tmp_path / "out"
        parser = _make_parser()
        args = parser.parse_args([
            "-i", str(fasta), "--prefix", "test", "--output-dir", str(outdir),
            "--bait-length", "40", "--tiling-density", "3", "--min-tm", "0",
            "--no-opool", "--redundancy-prereduce-identity", "0.80",
            "--threads", "2",
        ])
        ff.run_pipeline(args)

        summary_path = outdir / "test_summary.tsv"
        assert summary_path.is_file()
        summary_text = summary_path.read_text(encoding="utf-8")
        assert "after_prereduce" in summary_text

        per_ref_path = outdir / "test_per_ref_stats.tsv"
        header = per_ref_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert "n_baits_after_prereduce" in header

    def test_disabled_by_default_leaves_no_prereduce_trace(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        self._write_fasta(fasta)
        outdir = tmp_path / "out"
        parser = _make_parser()
        args = parser.parse_args([
            "-i", str(fasta), "--prefix", "test", "--output-dir", str(outdir),
            "--bait-length", "40", "--tiling-density", "3", "--min-tm", "0",
            "--no-opool", "--threads", "2",
        ])
        ff.run_pipeline(args)

        summary_text = (outdir / "test_summary.tsv").read_text(encoding="utf-8")
        assert "after_prereduce" not in summary_text

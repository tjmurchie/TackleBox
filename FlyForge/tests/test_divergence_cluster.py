"""Tests for divergence_bounded_cluster() -- a real vsearch-identity binary search that
GUARANTEES a hard maximum bait count, never searching below a real scientific divergence
floor (e.g. 0.75 identity = 25% max divergence).

Real motivating context, 2026-09-02: Tyler pushed back on simply disabling self-repeat
masking (v1.3.0) with nothing filling its actual role -- large panels (PLANT: 506,888
baits) became unaffordably big. The real fix isn't reference-count reduction (rejected,
correctly, as throwing away real biological variation to hit a number) -- it's a
divergence-tunable COLLAPSE mechanism: keep one real representative per group of baits
within an acceptable divergence of each other, with a HARD guarantee on the final count
(unlike --max-baits/--probe-num-cutoff, both confirmed to overshoot at real scale,
because self_blast_filter()'s id_cutoff loop can exhaust its range before crossing an
arbitrary low target). cd-hit-est can't do this -- confirmed hard floor of 0.80 identity
in nucleotide mode, stricter than a real 25%-divergence (0.75 identity) tolerance.
vsearch's `--id` genuinely accepts any value in [0, 1], and real production-scale timing
(15s-400s on real 726K-6.4M-candidate pools) makes an iterative search practical.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import FlyForge as ff


def _make_parser():
    return ff.build_arg_parser()


def _parse(parser, extra_args):
    return parser.parse_args(["-i", "x.fa", "--prefix", "x", *extra_args])


def _make_bait(i: int) -> "ff.Bait":
    return ff.Bait(
        bait_id=f"bait{i}", seq="ACGT" * 20, ref_id="ref1", ref_start=1, ref_end=80,
        gc_frac=0.5, masked_frac=0.0, ambiguous_count=0, tm=65.0,
    )


class TestDivergenceBoundedClusterLogic:
    """Fast tests: cluster_baits_vsearch() is monkeypatched with a deterministic,
    MONOTONIC fake (bait count increases with identity, exactly like real vsearch
    behavior confirmed on real ANIMAL/PLANT data) so these test the SEARCH algorithm's
    own correctness -- convergence, floor respect, infeasibility reporting -- without
    needing the real binary or real timing."""

    def _install_fake(self, monkeypatch, n_at_identity):
        calls = []

        def fake_cluster(baits, identity, threads=1, vsearch_path="vsearch"):
            calls.append(identity)
            n = n_at_identity(identity)
            return baits[:n], len(baits) - n

        monkeypatch.setattr(ff, "cluster_baits_vsearch", fake_cluster)
        return calls

    def test_already_under_cap_short_circuits_without_clustering(self, monkeypatch):
        calls = self._install_fake(monkeypatch, lambda identity: 999)
        baits = [_make_bait(i) for i in range(10)]
        kept, result = ff.divergence_bounded_cluster(
            baits, min_identity=0.75, hard_max_baits=100)
        assert len(kept) == 10
        assert result["feasible"] is True
        assert result["achieved_identity"] is None
        assert calls == []  # never even calls vsearch -- nothing to do

    def test_floor_infeasible_reports_clearly_but_returns_best_effort(self, monkeypatch):
        # Even at the floor (0.75), 500 baits remain -- over the 100 cap.
        self._install_fake(monkeypatch, lambda identity: 500)
        baits = [_make_bait(i) for i in range(1000)]
        kept, result = ff.divergence_bounded_cluster(
            baits, min_identity=0.75, hard_max_baits=100)
        assert result["feasible"] is False
        assert result["achieved_identity"] == 0.75
        assert len(kept) == 500  # best-effort result still returned, not dropped

    def test_finds_tightest_fitting_identity(self, monkeypatch):
        # n(0.75) = 100, n(1.0) = 1000, linear -- n(0.85) = 460 exactly, a clean,
        # checkable crossing point for hard_max_baits=460.
        def n_at(identity):
            frac = (identity - 0.75) / 0.25
            return int(round(100 + frac * 900))

        self._install_fake(monkeypatch, n_at)
        baits = [_make_bait(i) for i in range(1000)]
        kept, result = ff.divergence_bounded_cluster(
            baits, min_identity=0.75, hard_max_baits=460,
            identity_tolerance=0.001, max_search_iterations=20)
        assert result["feasible"] is True
        assert len(kept) <= 460
        assert abs(result["achieved_identity"] - 0.85) < 0.01

    def test_never_searches_below_min_identity(self, monkeypatch):
        seen = []

        def n_at(identity):
            seen.append(identity)
            return 100

        self._install_fake(monkeypatch, n_at)
        baits = [_make_bait(i) for i in range(1000)]
        ff.divergence_bounded_cluster(baits, min_identity=0.75, hard_max_baits=50)
        assert min(seen) >= 0.75  # floor doesn't fit, but never searches looser than it

    def test_respects_max_search_iterations(self, monkeypatch):
        # Floor (identity=0.0) gives 0 baits -- fits immediately, so the search proceeds
        # to hunt for the tightest fitting identity, bounded by max_search_iterations.
        calls = self._install_fake(
            monkeypatch, lambda identity: int(round(identity * 1000)))
        baits = [_make_bait(i) for i in range(1000)]
        ff.divergence_bounded_cluster(
            baits, min_identity=0.0, hard_max_baits=500, max_search_iterations=3)
        assert len(calls) <= 4  # 1 floor check + at most 3 search iterations

    def test_iterations_list_records_every_real_call(self, monkeypatch):
        self._install_fake(monkeypatch, lambda identity: 500)
        baits = [_make_bait(i) for i in range(1000)]
        kept, result = ff.divergence_bounded_cluster(
            baits, min_identity=0.75, hard_max_baits=100)
        assert len(result["iterations"]) == 1
        assert result["iterations"][0]["identity"] == 0.75
        assert result["iterations"][0]["n_baits"] == 500

    def test_result_never_exceeds_hard_cap_across_a_range_of_targets(self, monkeypatch):
        """The one property that actually matters for a 'hard cap': whatever target is
        asked for, the real returned count must never exceed it (when feasible)."""
        def n_at(identity):
            frac = (identity - 0.75) / 0.25
            return int(round(100 + frac * 900))

        self._install_fake(monkeypatch, n_at)
        baits = [_make_bait(i) for i in range(1000)]
        for cap in (150, 300, 500, 700, 999):
            kept, result = ff.divergence_bounded_cluster(
                baits, min_identity=0.75, hard_max_baits=cap,
                identity_tolerance=0.001, max_search_iterations=20)
            assert result["feasible"] is True
            assert len(kept) <= cap


class TestCLIValidation:
    """run_pipeline()'s own early validation, before any file I/O -- exits via
    sys.exit(1) before touching a real FASTA or external tool."""

    def _args(self, tmp_path, min_identity=None, hard_max=None):
        parser = _make_parser()
        extra = ["--output-dir", str(tmp_path)]
        if min_identity is not None:
            extra += ["--min-cluster-identity", str(min_identity)]
        if hard_max is not None:
            extra += ["--hard-max-baits", str(hard_max)]
        return _parse(parser, extra)

    def test_min_identity_without_hard_max_exits(self, tmp_path):
        args = self._args(tmp_path, min_identity=0.75)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_hard_max_without_min_identity_exits(self, tmp_path):
        args = self._args(tmp_path, hard_max=60000)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_min_identity_zero_exits(self, tmp_path):
        args = self._args(tmp_path, min_identity=0.0, hard_max=60000)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_min_identity_one_exits(self, tmp_path):
        args = self._args(tmp_path, min_identity=1.0, hard_max=60000)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_hard_max_zero_exits(self, tmp_path):
        args = self._args(tmp_path, min_identity=0.75, hard_max=0)
        with pytest.raises(SystemExit) as exc_info:
            ff.run_pipeline(args)
        assert exc_info.value.code == 1

    def test_valid_pair_does_not_exit_early(self, tmp_path):
        """Sanity check that a valid pair gets PAST validation (fails later for the
        ordinary reason -- no real input file -- not the same SystemExit(1))."""
        args = self._args(tmp_path, min_identity=0.75, hard_max=60000)
        with pytest.raises(FileNotFoundError):
            ff.run_pipeline(args)

    def test_neither_flag_set_does_not_exit_early(self, tmp_path):
        args = self._args(tmp_path)
        with pytest.raises(FileNotFoundError):
            ff.run_pipeline(args)


@pytest.mark.slow
class TestClusterBaitsVsearchIntegration:
    """Real end-to-end tests with the actual vsearch binary."""

    def _bait(self, bait_id, seq):
        return ff.Bait(bait_id=bait_id, seq=seq, ref_id="r", ref_start=1, ref_end=80,
                       gc_frac=0.5, masked_frac=0.0, ambiguous_count=0, tm=65.0)

    def test_real_vsearch_collapses_exact_duplicate(self):
        import random
        rng = random.Random(1)
        seq_a = "".join(rng.choice("ACGT") for _ in range(80))
        seq_b = "".join(rng.choice("ACGT") for _ in range(80))
        baits = [
            self._bait("a1", seq_a),
            self._bait("a2", seq_a),  # exact duplicate of a1
            self._bait("b1", seq_b),
        ]
        kept, removed = ff.cluster_baits_vsearch(baits, identity=0.95, threads=2)
        assert len(kept) == 2  # a1/a2 collapse to one representative
        assert removed == 1
        kept_ids = {b.bait_id for b in kept}
        assert "b1" in kept_ids
        assert len({"a1", "a2"} & kept_ids) == 1

    def test_real_divergence_bounded_cluster_end_to_end(self):
        import random
        rng = random.Random(2)
        seq_a = "".join(rng.choice("ACGT") for _ in range(80))
        seq_b = "".join(rng.choice("ACGT") for _ in range(80))
        seq_c = "".join(rng.choice("ACGT") for _ in range(80))
        baits = [
            self._bait("a1", seq_a), self._bait("a2", seq_a),  # exact dup pair
            self._bait("b1", seq_b),
            self._bait("c1", seq_c),
        ]
        kept, result = ff.divergence_bounded_cluster(
            baits, min_identity=0.75, hard_max_baits=3, threads=2)
        assert result["feasible"] is True
        assert len(kept) <= 3
        assert len(result["iterations"]) >= 1


@pytest.mark.slow
class TestPipelineIntegration:
    """Real end-to-end run_pipeline() calls with the actual vsearch binary."""

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

    def test_hard_max_baits_replaces_old_trio_and_produces_expected_output(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        self._write_fasta(fasta)
        outdir = tmp_path / "out"
        parser = _make_parser()
        # Real data check (see test_hard_cap_reports_infeasibility_when_unachievable
        # below): this dataset's 42 raw candidates only compress to 12 even at the
        # 0.75 floor -- 20 is a realistic, meaningfully-tighter-than-42 but genuinely
        # achievable cap, so this test exercises real binary-search convergence rather
        # than immediately hitting the floor/infeasible path.
        args = parser.parse_args([
            "-i", str(fasta), "--prefix", "test", "--output-dir", str(outdir),
            "--bait-length", "40", "--tiling-density", "3", "--min-tm", "0",
            "--no-opool", "--min-cluster-identity", "0.75", "--hard-max-baits", "20",
            "--threads", "2",
        ])
        ff.run_pipeline(args)

        summary_text = (outdir / "test_summary.tsv").read_text(encoding="utf-8")
        assert "after_divergence_cluster" in summary_text
        # The old trio must NOT have run when hard_max_baits is set.
        assert "after_prereduce" not in summary_text
        assert "after_self_blast_filter" not in summary_text
        assert "after_clustering" not in summary_text

        per_ref_path = outdir / "test_per_ref_stats.tsv"
        header = per_ref_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert "n_baits_after_divergence_cluster" in header

        final_baits = list(ff.read_fasta(str(outdir / "test_final_baits.fa")).keys())
        assert len(final_baits) <= 20  # the real, guaranteed hard cap

        assert "hard_max_baits_feasible\tTrue" in summary_text
        assert "hard_max_baits_achieved_identity\t0." in summary_text  # a real identity, not N/A

    def test_hard_cap_reports_infeasibility_when_unachievable(self, tmp_path):
        """This exact dataset's 42 raw candidates only compress to 12 even at the 0.75
        floor -- a cap of 5 is genuinely unachievable without violating the divergence
        floor. Must return a real best-effort result (not crash, not silently claim
        success) and log a clear warning -- exactly the real infeasibility signal
        Tyler's design asked for ("suggest reducing targets" rather than forcing it)."""
        fasta = tmp_path / "in.fasta"
        self._write_fasta(fasta)
        outdir = tmp_path / "out"
        progress_log_path = outdir / "test_progress.log"
        parser = _make_parser()
        args = parser.parse_args([
            "-i", str(fasta), "--prefix", "test", "--output-dir", str(outdir),
            "--bait-length", "40", "--tiling-density", "3", "--min-tm", "0",
            "--no-opool", "--min-cluster-identity", "0.75", "--hard-max-baits", "5",
            "--threads", "2", "--progress-log", str(progress_log_path),
        ])
        ff.run_pipeline(args)  # must not raise/crash

        progress_log = progress_log_path.read_text(encoding="utf-8")
        assert "WARNING" in progress_log
        assert "could not be met even at" in progress_log
        assert "INFEASIBLE" in progress_log

        final_baits = list(ff.read_fasta(str(outdir / "test_final_baits.fa")).keys())
        assert len(final_baits) == 12  # real best-effort result at the floor, not 5

        summary_text = (outdir / "test_summary.tsv").read_text(encoding="utf-8")
        assert "hard_max_baits_feasible\tFalse" in summary_text
        assert "hard_max_baits_achieved_identity\t0.75" in summary_text

    def test_disabled_by_default_leaves_no_divergence_cluster_trace(self, tmp_path):
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
        assert "after_divergence_cluster" not in summary_text
        assert "after_self_blast_filter" in summary_text  # old path still runs by default
        assert "hard_max_baits_feasible\tN/A" in summary_text
        assert "hard_max_baits_achieved_identity\tN/A" in summary_text

    @staticmethod
    def _read_per_ref_rows_with_final_baits(per_ref_path: Path) -> list[dict[str, str]]:
        lines = per_ref_path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]
        return [r for r in rows if int(r["n_baits_final"]) > 0]

    def test_hard_cap_carries_forward_skipped_trio_columns(self, tmp_path):
        """Real bug found 2026-09-03 while investigating PLANT's own real
        per_ref_stats.tsv (an ANIMAL --hard-max-baits run's per-ref columns showed
        n_baits_after_divergence_cluster > 0 for references where the preceding
        n_baits_after_prereduce/redundancy/cluster columns all read 0 -- internally
        impossible if those were real sequential counts). Root cause: when
        --hard-max-baits is set, the old prereduce/self-BLAST/clustering trio never
        runs at all, so RefStats' dataclass default of 0 silently stood in for "not
        applicable," indistinguishable from "this reference lost every bait here."
        Fixed by carrying the real post-divergence-cluster count forward into all
        three columns instead."""
        fasta = tmp_path / "in.fasta"
        self._write_fasta(fasta)
        outdir = tmp_path / "out"
        parser = _make_parser()
        args = parser.parse_args([
            "-i", str(fasta), "--prefix", "test", "--output-dir", str(outdir),
            "--bait-length", "40", "--tiling-density", "3", "--min-tm", "0",
            "--no-opool", "--min-cluster-identity", "0.75", "--hard-max-baits", "20",
            "--threads", "2",
        ])
        ff.run_pipeline(args)

        rows = self._read_per_ref_rows_with_final_baits(outdir / "test_per_ref_stats.tsv")
        assert rows  # sanity: the fixture must produce >=1 surviving reference
        for row in rows:
            div_cluster_count = int(row["n_baits_after_divergence_cluster"])
            assert div_cluster_count > 0
            assert int(row["n_baits_after_prereduce"]) == div_cluster_count
            assert int(row["n_baits_after_redundancy"]) == div_cluster_count
            assert int(row["n_baits_after_cluster"]) == div_cluster_count

    def test_soft_cap_carries_forward_divergence_cluster_column(self, tmp_path):
        """Mirror image of the hard-cap case above: without --hard-max-baits,
        divergence_cluster never runs, so ITS column must carry the trio's real
        final count forward instead of staying at the same misleading 0 default."""
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

        rows = self._read_per_ref_rows_with_final_baits(outdir / "test_per_ref_stats.tsv")
        assert rows
        for row in rows:
            cluster_count = int(row["n_baits_after_cluster"])
            assert cluster_count > 0
            assert int(row["n_baits_after_divergence_cluster"]) == cluster_count

    def test_skipped_tm_and_blast_filter_columns_carry_forward_not_zero(self, tmp_path):
        """Same bug class, but for the earlier, always-optional tm_filter (only runs
        if --min-tm > 0) and blast_filter (only runs if --blast-db is set) -- these
        were already silently zeroing their own columns whenever skipped, in EVERY
        run mode (not just --hard-max-baits), since before this fix they were never
        even reached when hard_max_baits was set either."""
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

        rows = self._read_per_ref_rows_with_final_baits(outdir / "test_per_ref_stats.tsv")
        assert rows
        for row in rows:
            mask_count = int(row["n_baits_after_mask"])
            assert mask_count > 0  # sanity on the fixture itself
            # --min-tm 0 skips tm_filter entirely; masked_filter is the immediately
            # preceding step with no other step in between, so this must be exact.
            assert int(row["n_baits_after_tm"]) == mask_count
            # No --blast-db set, so blast_filter is skipped -- must carry forward
            # the real (nonzero) count rather than the old stale-0 default.
            assert int(row["n_baits_after_blast"]) > 0

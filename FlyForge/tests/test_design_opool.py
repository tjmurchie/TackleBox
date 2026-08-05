"""Tests for FlyForge.design_opool() -- constructs the literal oligo
sequence that gets sent for physical synthesis:

    5'-[T7 promoter][probe][reverse_complement(amplification primer)]-3'

This is the single highest-stakes function in FlyForge: a subtle bug here
doesn't just produce a wrong report, it produces a real, expensive
synthesis order. The function does have one good internal runtime check
(asserting the selected primer's reverse complement starts with the
expected 'CGAAGAGC' BspQI/LguI motif before writing anything) but that only
guards THIS run's specific data -- it is not a regression test, and nothing
previously verified that the probe sequence itself survives into the final
oligo unmodified and in the right position.

design_opool() does a real exhaustive amplification-primer search (via
primer3) plus a real blastn subprocess call, so the full end-to-end test
takes roughly 30-90 seconds. It is marked `slow` and run once per test
session via a module-scoped fixture, with several fast assertions checked
against that one real run. Deliberately NOT mocked -- BLAST+ is a required
runtime dependency for this function already (it raises RuntimeError itself
if blastn is missing), and mocking away the primer-selection/BLAST step
would defeat the point of protecting the real assembled sequence.

Run just these tests with: pytest -m slow tests/test_design_opool.py
"""

import shutil

import pytest
from Bio import SeqIO

import FlyForge as ff

pytestmark = pytest.mark.slow


def _make_bait(bait_id: str, seq: str) -> "ff.Bait":
    return ff.Bait(
        bait_id=bait_id, seq=seq, ref_id="ref1", ref_start=1, ref_end=len(seq),
        gc_frac=0.5, masked_frac=0.0, ambiguous_count=0, tm=65.0,
    )


def _rand_seq(n: int, seed: int) -> str:
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture(scope="module")
def opool_result(tmp_path_factory):
    baits = [
        _make_bait("ref1|b1|pos1-80", _rand_seq(80, seed=1)),
        _make_bait("ref1|b2|pos30-110", _rand_seq(80, seed=2)),
        _make_bait("ref2|b1|pos1-80", _rand_seq(80, seed=3)),
    ]
    outdir = tmp_path_factory.mktemp("opool")
    logs = []
    oligo_path, bare_path, primer_fasta_path, primer_seq = ff.design_opool(
        baits, str(outdir), "test", logs.append,
    )
    return {
        "baits": baits,
        "oligo_path": oligo_path,
        "bare_path": bare_path,
        "primer_fasta_path": primer_fasta_path,
        "primer_seq": primer_seq,
        "logs": logs,
    }


class TestOligoStructure:
    def test_one_oligo_record_per_input_bait(self, opool_result):
        recs = {r.id: str(r.seq) for r in SeqIO.parse(opool_result["oligo_path"], "fasta")}
        assert set(recs.keys()) == {b.bait_id for b in opool_result["baits"]}

    def test_every_oligo_starts_with_the_t7_promoter(self, opool_result):
        for rec in SeqIO.parse(opool_result["oligo_path"], "fasta"):
            assert str(rec.seq).startswith(ff.T7_PROMOTER)

    def test_every_oligo_ends_with_the_expected_bspqi_lgui_motif(self, opool_result):
        """This mirrors design_opool's OWN internal runtime assertion, but as
        a regression test against the actual written FASTA output rather
        than trusting the function's internal check never regresses."""
        for rec in SeqIO.parse(opool_result["oligo_path"], "fasta"):
            primer_rc = str(rec.seq)[len(ff.T7_PROMOTER) + 80:]
            assert primer_rc.startswith("CGAAGAGC"), (
                f"{rec.id}: expected the post-probe suffix to start with the "
                f"BspQI/LguI motif CGAAGAGC, got {primer_rc[:12]!r}"
            )

    def test_probe_sequence_survives_unmodified_in_the_middle(self, opool_result):
        """The single most consequential correctness property: the exact
        probe sequence that passed all of FlyForge's filtering stages must
        appear, byte-for-byte, between the T7 prefix and the primer suffix
        -- not truncated, not reverse-complemented, not off-by-one shifted."""
        baits_by_id = {b.bait_id: b.seq.upper() for b in opool_result["baits"]}
        for rec in SeqIO.parse(opool_result["oligo_path"], "fasta"):
            full = str(rec.seq)
            middle = full[len(ff.T7_PROMOTER): len(ff.T7_PROMOTER) + 80]
            assert middle == baits_by_id[rec.id]

    def test_primer_rc_is_consistent_with_returned_primer_seq(self, opool_result):
        from Bio.Seq import Seq
        expected_rc = str(Seq(opool_result["primer_seq"]).reverse_complement())
        rec = next(SeqIO.parse(opool_result["oligo_path"], "fasta"))
        actual_rc = str(rec.seq)[len(ff.T7_PROMOTER) + 80:]
        assert actual_rc == expected_rc


class TestBareProbesOutput:
    def test_bare_probes_have_no_t7_or_primer_tails(self, opool_result):
        baits_by_id = {b.bait_id: b.seq.upper() for b in opool_result["baits"]}
        recs = {r.id: str(r.seq) for r in SeqIO.parse(opool_result["bare_path"], "fasta")}
        assert recs == baits_by_id


class TestPrimerFastaOutput:
    def test_primer_fasta_contains_t7_and_amplification_primer(self, opool_result):
        recs = {r.id: str(r.seq) for r in SeqIO.parse(opool_result["primer_fasta_path"], "fasta")}
        assert recs.get("T7_primer") == ff.T7_PROMOTER
        assert opool_result["primer_seq"] in recs.values()


class TestErrorHandling:
    def test_empty_bait_list_raises_without_needing_blast(self):
        with pytest.raises(RuntimeError, match="No baits"):
            ff.design_opool([], "/tmp", "empty", lambda msg: None)

    def test_missing_blastn_raises_actionable_error(self, monkeypatch, tmp_path):
        monkeypatch.setattr(ff.shutil, "which", lambda name: None)
        baits = [_make_bait("ref1|b1|pos1-80", _rand_seq(80, seed=99))]
        with pytest.raises(RuntimeError, match="blastn"):
            ff.design_opool(baits, str(tmp_path), "noblast", lambda msg: None)

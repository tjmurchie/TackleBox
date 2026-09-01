"""Tests for spinner.clustering.run_vsearch()'s cluster.by grouping.

Real bug, found 2026-09-01 via a live ~10,000-species PalaeoSCOPE Phase B run:
cluster.by (default ["species_guess", "marker_class"] in config.py -- every real
project inherits this default) was declared in the config schema but never actually
read by run_vsearch() -- every run clustered the ENTIRE input FASTA in ONE global
vsearch call regardless, which is fine at smoke-test scale (a few hundred records) but
catastrophic at real production scale (367,014 candidate records for one panel: 26+
CPU-hours and still zero output after 3.5+ wall-clock hours, a ~100 HOUR ETA).

vsearch itself is never invoked in these tests -- ``run()`` is monkeypatched with a
deterministic fake that mimics vsearch's own UC/centroids output shape (first record in
the group becomes the seed/centroid, the rest become members) so these tests exercise
the real grouping/merging/resume logic without any dependency on the real binary.
"""
from __future__ import annotations

import copy

from spinner.annotation import Annotation
from spinner.clustering import run_vsearch
from spinner.config import DEFAULT_CONFIG
from spinner.fasta import parse_fasta, write_fasta


def _ann(key: str, species: str, marker_class: str = "Mito") -> Annotation:
    return Annotation(
        accession=key,
        record_key=key,
        source_file="",
        header=f">{key} {species} {marker_class.lower()}",
        length=200,
        seq_sha256=key,
        species_guess=species,
        genus_guess=species.split()[0] if species else "",
        kingdom="Animal",
        marker_class=marker_class,
    )


def _write_input_fasta(path, ids):
    with path.open("w", encoding="utf-8") as f:
        for i in ids:
            f.write(f">{i} original=irrelevant\n")
            f.write("ACGT" * 20 + "\n")


def _cfg(**cluster_overrides):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["cluster"].update(cluster_overrides)
    return cfg


def _install_fake_vsearch(monkeypatch, calls):
    """Deterministic stand-in for vsearch --cluster_fast: first record in the group's
    own input FASTA becomes the seed/centroid, every other record becomes a member."""
    def fake_run(cmd, log=None, verbose=True):
        calls.append(list(cmd))
        group_fasta = cmd[2]
        cen_path = cmd[cmd.index("--centroids") + 1]
        uc_path = cmd[cmd.index("--uc") + 1]
        recs = parse_fasta(group_fasta)
        with open(uc_path, "w", encoding="utf-8") as f:
            f.write(f"S\t0\t*\t*\t*\t*\t*\t*\t{recs[0].id}\t*\n")
            for r in recs[1:]:
                f.write(f"H\t0\t*\t*\t*\t*\t*\t*\t{r.id}\t{recs[0].id}\n")
        write_fasta([recs[0]], cen_path)
        return ""

    monkeypatch.setattr("spinner.clustering.run", fake_run)


def test_by_unset_falls_back_to_exactly_one_global_vsearch_call(tmp_path, monkeypatch):
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ids = ["A.1", "B.1", "C.1"]
    ann = {i: _ann(i, species=f"Species{i}") for i in ids}
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, ids)
    cfg = _cfg(by=[])

    run_vsearch(str(fasta), str(tmp_path / "out"), cfg, ann, total_queries=len(ids))

    assert len(calls) == 1
    assert calls[0][2] == str(fasta)  # the ORIGINAL full input, not a per-group subset


def test_by_set_groups_by_species_and_marker_class(tmp_path, monkeypatch):
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    # Two records for Species A / Mito (a real clustering group), one singleton for
    # Species B / Mito, one singleton for Species A / Plastid (different marker_class,
    # must NOT be grouped with Species A / Mito despite sharing a species name).
    ann = {
        "A1.1": _ann("A1.1", "Species A", "Mito"),
        "A2.1": _ann("A2.1", "Species A", "Mito"),
        "B1.1": _ann("B1.1", "Species B", "Mito"),
        "A3.1": _ann("A3.1", "Species A", "Plastid"),
    }
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, list(ann.keys()))
    cfg = _cfg(by=["species_guess", "marker_class"])

    run_vsearch(str(fasta), str(tmp_path / "out"), cfg, ann, total_queries=len(ann))

    # Only the genuine 2-member group (Species A / Mito) needs a real vsearch call --
    # the two singletons must never trigger one.
    assert len(calls) == 1
    group_fasta_recs = parse_fasta(calls[0][2])
    assert {r.id for r in group_fasta_recs} == {"A1.1", "A2.1"}


def test_singleton_group_is_its_own_centroid_without_calling_vsearch(tmp_path, monkeypatch):
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ann = {"LONELY.1": _ann("LONELY.1", "Solo species")}
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, ["LONELY.1"])
    cfg = _cfg(by=["species_guess", "marker_class"])

    run_vsearch(str(fasta), str(tmp_path / "out"), cfg, ann, total_queries=1)

    assert calls == []
    assert ann["LONELY.1"].cluster_role == "centroid"
    assert "cluster_representative" in ann["LONELY.1"].reasons


def test_grouped_clustering_annotates_centroid_and_member_correctly(tmp_path, monkeypatch):
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ann = {
        "A1.1": _ann("A1.1", "Species A"),
        "A2.1": _ann("A2.1", "Species A"),
        "A3.1": _ann("A3.1", "Species A"),
    }
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, list(ann.keys()))
    cfg = _cfg(by=["species_guess", "marker_class"])

    run_vsearch(str(fasta), str(tmp_path / "out"), cfg, ann, total_queries=3)

    centroids = [k for k, a in ann.items() if a.cluster_role == "centroid"]
    members = [k for k, a in ann.items() if a.cluster_role == "member"]
    assert len(centroids) == 1
    assert len(members) == 2
    assert ann[members[0]].cluster_id == ann[members[1]].cluster_id == ann[centroids[0]].cluster_id


def test_cluster_ids_do_not_collide_across_independently_numbered_groups(tmp_path, monkeypatch):
    """Each group's own vsearch call independently numbers its clusters starting at 0 --
    merging without remapping would make two DIFFERENT species' unrelated cluster "0"s
    look identical."""
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ann = {
        "A1.1": _ann("A1.1", "Species A"),
        "A2.1": _ann("A2.1", "Species A"),
        "B1.1": _ann("B1.1", "Species B"),
        "B2.1": _ann("B2.1", "Species B"),
    }
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, list(ann.keys()))
    cfg = _cfg(by=["species_guess", "marker_class"])

    run_vsearch(str(fasta), str(tmp_path / "out"), cfg, ann, total_queries=4)

    assert len(calls) == 2
    id_a = ann["A1.1"].cluster_id
    id_b = ann["B1.1"].cluster_id
    assert id_a != id_b
    assert ann["A2.1"].cluster_id == id_a
    assert ann["B2.1"].cluster_id == id_b


def test_resume_skips_group_whose_uc_file_already_exists(tmp_path, monkeypatch):
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ann = {
        "A1.1": _ann("A1.1", "Species A"),
        "A2.1": _ann("A2.1", "Species A"),
    }
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, list(ann.keys()))
    cfg = _cfg(by=["species_guess", "marker_class"])
    outprefix = str(tmp_path / "out")

    # Pre-populate group 0's own UC file, as if a prior partial run already completed it.
    group_dir = tmp_path / "out.vsearch_groups"
    group_dir.mkdir()
    (group_dir / "group_000000.uc").write_text(
        "S\t0\t*\t*\t*\t*\t*\t*\tA1.1\t*\nH\t0\t*\t*\t*\t*\t*\t*\tA2.1\tA1.1\n",
        encoding="utf-8",
    )

    run_vsearch(str(fasta), outprefix, cfg, ann, total_queries=2)

    assert calls == []  # vsearch never invoked -- the cached UC file was reused
    assert ann["A1.1"].cluster_role == "centroid"
    assert ann["A2.1"].cluster_role == "member"


def test_records_missing_from_ann_are_grouped_together_rather_than_crashing(tmp_path, monkeypatch):
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ann = {"KNOWN.1": _ann("KNOWN.1", "Species A")}
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, ["KNOWN.1", "ORPHAN.1"])  # ORPHAN.1 has no ann entry
    cfg = _cfg(by=["species_guess", "marker_class"])

    run_vsearch(str(fasta), str(tmp_path / "out"), cfg, ann, total_queries=2)

    # Must not raise -- ORPHAN.1 is simply absent from any resulting annotation.
    assert "ORPHAN.1" not in ann
    assert ann["KNOWN.1"].cluster_role == "centroid"


def test_merged_output_files_are_written_at_the_usual_paths(tmp_path, monkeypatch):
    """Downstream code reads <outprefix>.clusters.uc / .cluster_centroids.fasta --
    the grouped code path must still produce them at the same paths as the old
    single-global-call behavior."""
    calls = []
    _install_fake_vsearch(monkeypatch, calls)
    ann = {
        "A1.1": _ann("A1.1", "Species A"),
        "A2.1": _ann("A2.1", "Species A"),
        "B1.1": _ann("B1.1", "Species B"),
    }
    fasta = tmp_path / "input.fasta"
    _write_input_fasta(fasta, list(ann.keys()))
    cfg = _cfg(by=["species_guess", "marker_class"])
    outprefix = str(tmp_path / "out")

    run_vsearch(str(fasta), outprefix, cfg, ann, total_queries=3)

    uc_path = tmp_path / "out.clusters.uc"
    cen_path = tmp_path / "out.cluster_centroids.fasta"
    assert uc_path.is_file() and uc_path.stat().st_size > 0
    assert cen_path.is_file() and cen_path.stat().st_size > 0
    centroid_ids = {r.id for r in parse_fasta(str(cen_path))}
    assert centroid_ids == {"A1.1", "B1.1"}  # one centroid per real group

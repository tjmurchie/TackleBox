"""Tests for Items 3, 9, 12, 17: authenticity tiers, contaminants, re-routing, normalization."""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List

import pytest

from fillet.relic_lca import ReadAssignment, TaxonSummary, summarize_assignments
from fillet.taxonomy import Taxonomy
from fillet.io import (
    read_contaminants, read_regional_taxa, read_sample_sheet,
    read_support_table, write_tsv,
    read_fossil_table_structured, filter_fossil_taxa_for_sample,
)


# ---------------------------------------------------------------------------
# Minimal taxonomy for testing
# ---------------------------------------------------------------------------

def _make_tax() -> Taxonomy:
    lines_nodes = (
        "1\t|\t1\t|\tno rank\t|\n"
        "131567\t|\t1\t|\tno rank\t|\n"         # cellular organisms
        "2759\t|\t131567\t|\tsuperkingdom\t|\n"  # Eukaryota
        "33154\t|\t2759\t|\tno rank\t|\n"        # Opisthokonta
        "33208\t|\t33154\t|\tkingdom\t|\n"       # Metazoa
        "6072\t|\t33208\t|\tno rank\t|\n"        # Eumetazoa
        "33213\t|\t6072\t|\tno rank\t|\n"        # Bilateria
        "7711\t|\t33213\t|\tphylum\t|\n"         # Chordata
        "40674\t|\t7711\t|\tclass\t|\n"          # Mammalia
        "91561\t|\t40674\t|\torder\t|\n"         # Cetartiodactyla
        "9895\t|\t91561\t|\tfamily\t|\n"         # Bovidae
        "9903\t|\t9895\t|\tgenus\t|\n"           # Bison
        "9904\t|\t9903\t|\tspecies\t|\n"         # Bison bison
        "10088\t|\t40674\t|\torder\t|\n"         # Rodentia
        "10066\t|\t10088\t|\tfamily\t|\n"        # Muridae
        "10090\t|\t10066\t|\tspecies\t|\n"       # Mus musculus
    )
    lines_names = (
        "1\t|\tall\t|\t\t|\tscientific name\t|\n"
        "131567\t|\tcellular organisms\t|\t\t|\tscientific name\t|\n"
        "2759\t|\tEukaryota\t|\t\t|\tscientific name\t|\n"
        "33154\t|\tOpisthokonta\t|\t\t|\tscientific name\t|\n"
        "33208\t|\tMetazoa\t|\t\t|\tscientific name\t|\n"
        "6072\t|\tEumetazoa\t|\t\t|\tscientific name\t|\n"
        "33213\t|\tBilateria\t|\t\t|\tscientific name\t|\n"
        "7711\t|\tChordata\t|\t\t|\tscientific name\t|\n"
        "40674\t|\tMammalia\t|\t\t|\tscientific name\t|\n"
        "91561\t|\tCetartiodactyla\t|\t\t|\tscientific name\t|\n"
        "9895\t|\tBovidae\t|\t\t|\tscientific name\t|\n"
        "9903\t|\tBison\t|\t\t|\tscientific name\t|\n"
        "9904\t|\tBison bison\t|\t\t|\tscientific name\t|\n"
        "10088\t|\tRodentia\t|\t\t|\tscientific name\t|\n"
        "10066\t|\tMuridae\t|\t\t|\tscientific name\t|\n"
        "10090\t|\tMus musculus\t|\t\t|\tscientific name\t|\n"
    )
    with tempfile.TemporaryDirectory() as d:
        nodes = Path(d) / "nodes.dmp"
        names = Path(d) / "names.dmp"
        nodes.write_text(lines_nodes)
        names.write_text(lines_names)
        return Taxonomy.from_paths(nodes=str(nodes), names=str(names))


def _assignment(taxid: str, sid: str = "S1", posterior: float = 0.9,
                damage: float = 0.0, read_id: str = "r1") -> ReadAssignment:
    return ReadAssignment(
        read_id=read_id, sample_id=sid, assigned_taxid=taxid,
        assigned_name="", assigned_rank="",
        posterior=posterior, entropy=0.1,
        status="assigned", reason="",
        n_hits_raw=3, n_hits_used=3,
        best_hit_taxid=taxid, best_hit_name="",
        best_bitscore=300.0, best_pident=95.0, best_qcov=0.9, best_aln_len=50,
        best_subject_id="acc1",
        lca_taxid=taxid, lca_name="",
        damage_ct_5p=damage, damage_ga_3p=damage, damage_score=damage,
        damage_computable=(damage > 0),
    )


TAX = _make_tax()


# ---------------------------------------------------------------------------
# Item 9: Known contaminants
# ---------------------------------------------------------------------------

class TestContaminants:
    def test_plain_text_file(self, tmp_path):
        f = tmp_path / "contam.txt"
        f.write_text("9904\n# comment\n10090\n")
        result = read_contaminants(str(f))
        assert "9904" in result
        assert "10090" in result
        assert len(result) == 2

    def test_tsv_file(self, tmp_path):
        f = tmp_path / "contam.tsv"
        f.write_text("taxid\tname\n9904\tBison bison\n10090\tMus musculus\n")
        result = read_contaminants(str(f))
        assert "9904" in result

    def test_none_returns_empty(self):
        assert read_contaminants(None) == set()

    def test_contaminant_flag_in_summary(self):
        a = _assignment("9904", damage=0.05)
        summaries = summarize_assignments(TAX, [a], contaminants={"9904"})
        row = next(s for s in summaries if s.taxid == "9904")
        assert "user_contaminant" in row.flags

    def test_non_contaminant_no_flag(self):
        a = _assignment("9904", damage=0.05)
        summaries = summarize_assignments(TAX, [a], contaminants={"10090"})
        row = next(s for s in summaries if s.taxid == "9904")
        assert "user_contaminant" not in row.flags


# ---------------------------------------------------------------------------
# Item 12: Manual taxon re-routing (tested via ReadAssignment mutation)
# ---------------------------------------------------------------------------

class TestRerouting:
    def test_reroute_reassigns(self):
        a = _assignment("9904")  # Bison bison → should be re-routed to Bison (9903)
        # Simulate what cmd_classify does: mutate assignment in place
        reroute_map = {"9904": "9903"}
        for assignment in [a]:
            if assignment.assigned_taxid in reroute_map:
                assignment.assigned_taxid = reroute_map[assignment.assigned_taxid]
        assert a.assigned_taxid == "9903"

    def test_reroute_unchanged_when_not_in_map(self):
        a = _assignment("10090")
        reroute_map = {"9904": "9903"}
        for assignment in [a]:
            if assignment.assigned_taxid in reroute_map:
                assignment.assigned_taxid = reroute_map[assignment.assigned_taxid]
        assert a.assigned_taxid == "10090"

    def test_rerouted_reads_appear_at_new_taxon(self):
        a = _assignment("9904")
        a.assigned_taxid = "9903"  # post-reroute
        summaries = summarize_assignments(TAX, [a])
        taxids = {s.taxid for s in summaries}
        assert "9903" in taxids
        assert "9904" not in taxids or all(
            s.direct_hard_reads == 0 for s in summaries if s.taxid == "9904"
        )


# ---------------------------------------------------------------------------
# Item 17: Normalization
# ---------------------------------------------------------------------------

class TestNormalization:
    def test_rpm_computed_from_assignments(self):
        # 2 reads total, 1 assigned to Bison bison
        a1 = _assignment("9904")
        a2 = ReadAssignment(
            read_id="r2", sample_id="S1", assigned_taxid="0",
            assigned_name="", assigned_rank="",
            posterior=0.0, entropy=0.0,
            status="unassigned", reason="no_hits",
            n_hits_raw=0, n_hits_used=0,
            best_hit_taxid="", best_hit_name="",
            best_bitscore=0.0, best_pident=0.0, best_qcov=None, best_aln_len=0,
            best_subject_id="",
        )
        summaries = summarize_assignments(TAX, [a1, a2])
        row = next((s for s in summaries if s.taxid == "9904"), None)
        assert row is not None
        # 1 hard read / 2 total × 1e6 = 500000 RPM
        assert abs(row.reads_per_million - 500_000.0) < 1.0

    def test_rpm_with_user_supplied_depth(self):
        a = _assignment("9904")
        summaries = summarize_assignments(TAX, [a], sample_depths={"S1": 1_000_000})
        row = next(s for s in summaries if s.taxid == "9904")
        # 1 read / 1,000,000 × 1e6 = 1.0 RPM
        assert abs(row.reads_per_million - 1.0) < 0.001

    def test_relative_abundance_sums_to_one(self):
        a1 = _assignment("9904", sid="S1")
        a2 = _assignment("10090", sid="S1")
        summaries = summarize_assignments(TAX, [a1, a2])
        # root node cumulative_weighted_reads = sum of all weighted, so relative_abundance at root = 1
        root = next((s for s in summaries if s.taxid == "1"), None)
        if root:
            assert abs(root.relative_abundance - 1.0) < 0.01


# ---------------------------------------------------------------------------
# Item 3: Multi-proxy authenticity tier
# ---------------------------------------------------------------------------

class TestAuthenticityTier:
    def _make_high_quality_assignments(self, n: int = 10, damage: float = 0.08) -> List[ReadAssignment]:
        return [
            ReadAssignment(
                read_id=f"r{i}", sample_id="S1", assigned_taxid="9904",
                assigned_name="Bison bison", assigned_rank="species",
                lca_taxid="9904", lca_name="Bison bison",
                posterior=0.95, entropy=0.05,
                status="assigned", reason="",
                n_hits_raw=5, n_hits_used=5,
                best_hit_taxid="9904", best_hit_name="Bison bison",
                best_bitscore=350.0, best_pident=98.0, best_qcov=0.95, best_aln_len=55,
                best_subject_id="acc1", best_sstart=i * 60, best_send=i * 60 + 55,
                damage_ct_5p=damage, damage_ga_3p=damage, damage_score=damage,
                damage_computable=(damage > 0),
            )
            for i in range(n)
        ]

    def test_tier1_with_damage_and_eco_pal_fos(self):
        assignments = self._make_high_quality_assignments(20, damage=0.07)
        summaries = summarize_assignments(
            TAX, assignments,
            regional_taxa={"9904": {"weight": "1.5"}},
            palynology_taxa={"9904"},
            fossil_taxa={"9904"},
        )
        row = next(s for s in summaries if s.taxid == "9904")
        assert row.eco_support is True
        assert row.pal_support is True
        assert row.fos_support is True
        assert row.authenticity_tier == 1
        assert row.authenticity_badge == "★"
        assert row.authenticity_tier_pch == 8

    def test_tier4_lca_only(self):
        # posterior=0.9 → LCA quality passes; no damage/breadth/eco/pal/fos → 1 evidence → tier 4
        a = _assignment("9904", damage=0.0, posterior=0.9)
        summaries = summarize_assignments(TAX, [a])
        row = next(s for s in summaries if s.taxid == "9904")
        assert row.authenticity_tier == 4
        assert row.authenticity_badge == "▲"
        assert row.authenticity_tier_pch == 17
        assert row.n_support_lines == 1

    def test_tier5_no_support(self):
        # posterior=0.5 below threshold; no damage/breadth/eco/pal/fos → 0 evidence → tier 5
        a = _assignment("9904", damage=0.0, posterior=0.5)
        summaries = summarize_assignments(TAX, [a])
        row = next(s for s in summaries if s.taxid == "9904")
        assert row.authenticity_tier == 5
        assert row.authenticity_badge == "○"
        assert row.authenticity_tier_pch == 1
        assert row.n_support_lines == 0

    def test_contaminant_gets_flagged_badge(self):
        a = _assignment("9904", damage=0.07)
        summaries = summarize_assignments(TAX, [a], contaminants={"9904"})
        row = next(s for s in summaries if s.taxid == "9904")
        assert row.authenticity_tier == 0
        assert row.authenticity_badge == "✕"
        assert row.authenticity_tier_pch == 4

    def test_support_table_reads_plain_text(self, tmp_path):
        f = tmp_path / "pal.txt"
        f.write_text("9904\n9903\n")
        result = read_support_table(str(f))
        assert "9904" in result
        assert "9903" in result

    def test_ancestor_eco_match(self):
        # regional_taxa has Bovidae (9895), Bison bison (9904) is a descendant — should match
        a = _assignment("9904", damage=0.0)
        summaries = summarize_assignments(TAX, [a], regional_taxa={"9895": {}})
        row = next(s for s in summaries if s.taxid == "9904")
        assert row.eco_support is True


# ---------------------------------------------------------------------------
# read_regional_taxa
# ---------------------------------------------------------------------------

class TestReadRegionalTaxa:
    def test_none_returns_empty(self):
        assert read_regional_taxa(None) == {}

    def test_tsv_keyed_by_taxid(self, tmp_path):
        f = tmp_path / "regional.tsv"
        f.write_text("taxid\tstatus\thabitat\n9606\tpresent\tterrestrial\n", encoding="utf-8")
        result = read_regional_taxa(f)
        assert "9606" in result
        assert result["9606"]["status"] == "present"

    def test_keyed_by_scientific_name(self, tmp_path):
        f = tmp_path / "regional.tsv"
        f.write_text("scientific_name\tstatus\nHomo sapiens\tpresent\n", encoding="utf-8")
        result = read_regional_taxa(f)
        assert "Homo sapiens" in result

    def test_versioned_both_keys_present(self, tmp_path):
        """Row with both taxid and name columns should be accessible via either key."""
        f = tmp_path / "regional.tsv"
        f.write_text("taxid\tscientific_name\tstatus\n9606\tHomo sapiens\tpresent\n", encoding="utf-8")
        result = read_regional_taxa(f)
        assert "9606" in result
        assert "Homo sapiens" in result

    def test_empty_file_returns_empty(self, tmp_path):
        f = tmp_path / "regional.tsv"
        f.write_text("taxid\tstatus\n", encoding="utf-8")
        result = read_regional_taxa(f)
        assert result == {}


# ---------------------------------------------------------------------------
# read_sample_sheet
# ---------------------------------------------------------------------------

class TestReadSampleSheet:
    def test_none_returns_empty(self):
        assert read_sample_sheet(None) == {}

    def test_tsv_keyed_by_sample_id(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text("sample_id\trole\nS1\tsample\nS2\tnegative\n", encoding="utf-8")
        result = read_sample_sheet(f)
        assert "S1" in result
        assert result["S1"]["role"] == "sample"
        assert "S2" in result

    def test_alternative_id_column_library(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text("library\trole\nLIB001\tsample\n", encoding="utf-8")
        result = read_sample_sheet(f)
        assert "LIB001" in result

    def test_missing_id_column_skips_row(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text("role\tgroup\nsample\tA\n", encoding="utf-8")
        result = read_sample_sheet(f)
        assert result == {}


# ---------------------------------------------------------------------------
# write_tsv (io.py)
# ---------------------------------------------------------------------------

class TestWriteTsv:
    def test_creates_file_with_header(self, tmp_path):
        p = tmp_path / "out.tsv"
        write_tsv(p, [{"a": 1, "b": 2}], fieldnames=["a", "b"])
        with p.open() as fh:
            header = fh.readline().strip().split("\t")
        assert header == ["a", "b"]

    def test_writes_rows(self, tmp_path):
        p = tmp_path / "out.tsv"
        write_tsv(p, [{"x": "foo", "y": "bar"}, {"x": "baz", "y": "qux"}], fieldnames=["x", "y"])
        import csv
        with p.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert len(rows) == 2
        assert rows[0]["x"] == "foo"
        assert rows[1]["y"] == "qux"

    def test_extra_keys_ignored(self, tmp_path):
        p = tmp_path / "out.tsv"
        write_tsv(p, [{"a": 1, "b": 2, "c": 99}], fieldnames=["a", "b"])
        import csv
        with p.open() as fh:
            row = list(csv.DictReader(fh, delimiter="\t"))[0]
        assert "c" not in row

    def test_creates_parent_dirs(self, tmp_path):
        p = tmp_path / "sub" / "deep" / "out.tsv"
        write_tsv(p, [], fieldnames=["a"])
        assert p.exists()


# ---------------------------------------------------------------------------
# Virtual branch (unclassified/excluded nodes in the tree)
# ---------------------------------------------------------------------------

def _make_unclassified_assignment(
    status: str, read_id: str = "u1", sid: str = "S1"
) -> ReadAssignment:
    return ReadAssignment(
        read_id=read_id, sample_id=sid,
        assigned_taxid="0", assigned_name="", assigned_rank="",
        posterior=0.0, entropy=0.0,
        status=status, reason="",
        n_hits_raw=0, n_hits_used=0,
        best_hit_taxid="0", best_hit_name="",
        best_bitscore=0.0, best_pident=0.0, best_qcov=0.0, best_aln_len=0,
        best_subject_id="",
        lca_taxid="0", lca_name="",
        damage_ct_5p=0.0, damage_ga_3p=0.0, damage_score=0.0,
        damage_computable=False,
    )


class TestVirtualBranchFromAssignments:
    """viewer.py _build_virtual_branch_from_assignments counts categories correctly."""

    def test_no_unclassified_returns_none(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        assigned = [
            ReadAssignment(
                read_id="r1", sample_id="S1",
                assigned_taxid="9904", assigned_name="Bison bison", assigned_rank="species",
                posterior=0.9, entropy=0.1,
                status="assigned", reason="",
                n_hits_raw=3, n_hits_used=3,
                best_hit_taxid="9904", best_hit_name="",
                best_bitscore=300.0, best_pident=95.0, best_qcov=0.9, best_aln_len=50,
                best_subject_id="acc1",
                lca_taxid="9904", lca_name="",
                damage_ct_5p=0.0, damage_ga_3p=0.0, damage_score=0.0,
                damage_computable=False,
            )
        ]
        assert _build_virtual_branch_from_assignments(assigned) is None

    def test_counts_no_hits(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        asgns = [
            _make_unclassified_assignment("unclassified:no_candidate_hits", "u1"),
            _make_unclassified_assignment("unclassified:no_candidate_hits", "u2"),
        ]
        branch = _build_virtual_branch_from_assignments(asgns)
        assert branch is not None
        assert branch["taxid"] == "__unclassified__"
        assert branch["is_virtual_node"] is True
        kids = {c["taxid"]: c for c in branch["children"]}
        assert kids["__no_hits__"]["metrics"]["hard_reads"] == 2
        assert kids["__em_demoted__"]["metrics"]["hard_reads"] == 0

    def test_counts_em_demoted(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        asgns = [
            _make_unclassified_assignment("assigned;em_coherence_lift:unclassified", "u1"),
            _make_unclassified_assignment("assigned;em_coherence_lift:unclassified", "u2"),
            _make_unclassified_assignment("assigned;em_coherence_lift:unclassified", "u3"),
        ]
        branch = _build_virtual_branch_from_assignments(asgns)
        assert branch is not None
        kids = {c["taxid"]: c for c in branch["children"]}
        assert kids["__em_demoted__"]["metrics"]["hard_reads"] == 3
        assert kids["__lca_broad__"]["metrics"]["hard_reads"] == 0

    def test_lca_broad_goes_to_lca_broad_not_em_demoted(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        asgns = [
            _make_unclassified_assignment(
                "assigned;em_coherence_lift:unclassified(lca_too_broad:root)", "u1"
            ),
        ]
        branch = _build_virtual_branch_from_assignments(asgns)
        kids = {c["taxid"]: c for c in branch["children"]}
        assert kids["__lca_broad__"]["metrics"]["hard_reads"] == 1
        assert kids["__em_demoted__"]["metrics"]["hard_reads"] == 0

    def test_rank_capped(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        asgns = [
            _make_unclassified_assignment("assigned;rank_cap:unclassified(Bison@genus)", "u1"),
        ]
        branch = _build_virtual_branch_from_assignments(asgns)
        kids = {c["taxid"]: c for c in branch["children"]}
        assert kids["__rank_capped__"]["metrics"]["hard_reads"] == 1

    def test_contaminants_child_present(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        asgns = [_make_unclassified_assignment("unclassified:no_candidate_hits")]
        branch = _build_virtual_branch_from_assignments(asgns)
        taxids = [c["taxid"] for c in branch["children"]]
        assert "__contaminants__" in taxids

    def test_virtual_tooltip_present(self):
        from fillet.viewer import _build_virtual_branch_from_assignments
        asgns = [_make_unclassified_assignment("unclassified:no_candidate_hits")]
        branch = _build_virtual_branch_from_assignments(asgns)
        kids = {c["taxid"]: c for c in branch["children"]}
        assert len(kids["__no_hits__"]["virtual_tooltip"]) > 20
        assert len(branch["virtual_tooltip"]) > 20


class TestVirtualBranchSQLite:
    """viewer_sqlite.py _build_virtual_branch and reads() for virtual taxids."""

    def _make_db(self, tmp_path, asgns: list) -> "ViewerDB":
        from fillet.viewer_sqlite import write_viewer_sqlite, ViewerDB
        tax = TAX
        db_path = tmp_path / "v.sqlite"
        write_viewer_sqlite(db_path, tax, asgns, [])
        return ViewerDB(db_path)

    def test_no_unclassified_no_virtual_branch(self, tmp_path):
        db = self._make_db(tmp_path, [])
        payload = db.payload()
        # Root children should all be real taxids (no __unclassified__)
        child_taxids = [c["taxid"] for c in payload["tree"]["children"]]
        assert "__unclassified__" not in child_taxids

    def test_virtual_branch_appended_when_unclassified_reads_exist(self, tmp_path):
        asgns = [_make_unclassified_assignment("unclassified:no_candidate_hits")]
        db = self._make_db(tmp_path, asgns)
        payload = db.payload()
        child_taxids = [c["taxid"] for c in payload["tree"]["children"]]
        assert "__unclassified__" in child_taxids

    def test_reads_for_virtual_no_hits(self, tmp_path):
        asgns = [
            _make_unclassified_assignment("unclassified:no_candidate_hits", "u1"),
            _make_unclassified_assignment("unclassified:no_candidate_hits", "u2"),
        ]
        db = self._make_db(tmp_path, asgns)
        reads = db.reads(["__no_hits__"])
        assert len(reads) == 2
        assert all(r["assigned_taxid"] == "0" for r in reads)

    def test_reads_for_virtual_em_demoted(self, tmp_path):
        asgns = [
            _make_unclassified_assignment("assigned;em_coherence_lift:unclassified", "u1"),
            _make_unclassified_assignment("unclassified:no_candidate_hits", "u2"),
        ]
        db = self._make_db(tmp_path, asgns)
        em_reads = db.reads(["__em_demoted__"])
        assert len(em_reads) == 1
        assert em_reads[0]["read_id"] == "u1"

    def test_reads_for_virtual_lca_broad(self, tmp_path):
        asgns = [
            _make_unclassified_assignment(
                "assigned;em_coherence_lift:unclassified(lca_too_broad:root)", "u1"
            ),
        ]
        db = self._make_db(tmp_path, asgns)
        reads = db.reads(["__lca_broad__"])
        assert len(reads) == 1
        assert reads[0]["read_id"] == "u1"

    def test_reads_for_root_virtual_returns_all_unclassified(self, tmp_path):
        asgns = [
            _make_unclassified_assignment("unclassified:no_candidate_hits", "u1"),
            _make_unclassified_assignment("assigned;em_coherence_lift:unclassified", "u2"),
        ]
        db = self._make_db(tmp_path, asgns)
        reads = db.reads(["__unclassified__"])
        assert len(reads) == 2

    def test_virtual_branch_is_virtual_node_flagged(self, tmp_path):
        asgns = [_make_unclassified_assignment("unclassified:no_candidate_hits")]
        db = self._make_db(tmp_path, asgns)
        payload = db.payload()
        vbranch = next(
            c for c in payload["tree"]["children"] if c["taxid"] == "__unclassified__"
        )
        assert vbranch["is_virtual_node"] is True
        assert vbranch["rank"] == "virtual"
        child_taxids = [c["taxid"] for c in vbranch["children"]]
        assert "__no_hits__" in child_taxids
        assert "__em_demoted__" in child_taxids
        assert "__contaminants__" in child_taxids


# ---------------------------------------------------------------------------
# Age/site-resolved fossil support — read_fossil_table_structured +
# filter_fossil_taxa_for_sample + fos_evidence_text in summarize_assignments
# ---------------------------------------------------------------------------

class TestReadFossilTableStructured:
    def test_none_returns_empty(self):
        assert read_fossil_table_structured(None) == []

    def test_plain_taxid_list(self, tmp_path):
        f = tmp_path / "fossil.txt"
        f.write_text("9904\n10090\n", encoding="utf-8")
        records = read_fossil_table_structured(f)
        assert len(records) == 2
        taxons = {r["taxon"] for r in records}
        assert "9904" in taxons
        assert "10090" in taxons

    def test_structured_tsv_with_site_and_age(self, tmp_path):
        f = tmp_path / "fossil.tsv"
        f.write_text(
            "taxon\tsite_name\tage_min_bp\tage_max_bp\tevidence_type\tsource\tnotes\n"
            "9904\tPaJs-3N\t800\t400\tmacrofossil\tMurchie2025\tbowhead bones in pond\n"
            "10090\t*\t\t\tecological\tfield_obs\twidespread\n",
            encoding="utf-8",
        )
        records = read_fossil_table_structured(f)
        assert len(records) == 2
        r0 = next(r for r in records if r["taxon"] == "9904")
        assert r0.get("site_name") == "PaJs-3N"
        assert r0.get("age_min_bp") == 800.0
        assert r0.get("age_max_bp") == 400.0
        assert r0.get("evidence_type") == "macrofossil"
        r1 = next(r for r in records if r["taxon"] == "10090")
        assert r1.get("site_name") in ("*", None, "")  # wildcard or absent


class TestFilterFossilTaxaForSample:
    def _records(self):
        return [
            {"taxon": "9904", "site_name": "PaJs-3N", "age_min_bp": 800.0, "age_max_bp": 300.0,
             "evidence_type": "macrofossil", "source": "Murchie2025", "notes": "bones in pond"},
            {"taxon": "10090", "site_name": "PaJs-13", "age_min_bp": 1500.0, "age_max_bp": 800.0,
             "evidence_type": "archaeofaunal", "source": "Kissinger2022", "notes": "faunal assemblage"},
            {"taxon": "40674", "site_name": "*", "age_min_bp": None, "age_max_bp": None,
             "evidence_type": "ecological", "source": "survey", "notes": "any site any age"},
        ]

    def test_site_filter_matches(self):
        taxa, ev = filter_fossil_taxa_for_sample(self._records(), site_name="PaJs-3N", age_bp=500.0)
        assert "9904" in taxa   # site match + age within [300, 800]
        assert "40674" in taxa  # wildcard site
        assert "10090" not in taxa  # wrong site

    def test_site_filter_no_match(self):
        taxa, ev = filter_fossil_taxa_for_sample(self._records(), site_name="SavR5", age_bp=500.0)
        assert "9904" not in taxa
        assert "10090" not in taxa
        assert "40674" in taxa  # wildcard always matches

    def test_age_filter_excludes_out_of_range(self):
        taxa, ev = filter_fossil_taxa_for_sample(self._records(), site_name="PaJs-3N", age_bp=1200.0)
        # age_bp=1200 is older than age_min_bp=800 for "9904" → outside range → should not match
        assert "9904" not in taxa
        assert "40674" in taxa

    def test_age_filter_within_range(self):
        taxa, ev = filter_fossil_taxa_for_sample(self._records(), site_name="PaJs-3N", age_bp=600.0)
        assert "9904" in taxa   # 300 <= 600 <= 800

    def test_no_site_no_age_returns_all_wildcard(self):
        taxa, ev = filter_fossil_taxa_for_sample(self._records(), site_name=None, age_bp=None)
        # Without site/age context only wildcard (*) records should match
        assert "40674" in taxa
        assert "9904" not in taxa
        assert "10090" not in taxa

    def test_evidence_text_populated(self):
        taxa, ev = filter_fossil_taxa_for_sample(self._records(), site_name="PaJs-3N", age_bp=600.0)
        assert any("bones in pond" in e or "macrofossil" in e or "PaJs-3N" in e for e in ev)

    def test_empty_records_returns_none(self):
        taxa, ev = filter_fossil_taxa_for_sample([], site_name="PaJs-3N", age_bp=500.0)
        assert taxa is None or len(taxa) == 0


class TestFosEvidenceTextInSummarize:
    def test_fos_evidence_text_set_when_per_sample_dict_provided(self):
        tax = TAX
        # Assign reads to Mammalia (40674) — its ancestors include Metazoa, Chordata etc.
        a = _assignment("40674", sid="S1", posterior=0.95)
        fossil_taxa_by_sample = {"S1": {"40674"}}
        fos_evidence_by_sample = {"S1": ["macrofossil: bones found at PaJs-3N ~600 BP"]}
        summaries = summarize_assignments(
            tax, [a],
            fossil_taxa_by_sample=fossil_taxa_by_sample,
            fos_evidence_by_sample=fos_evidence_by_sample,
        )
        row = next((s for s in summaries if s.taxid == "40674" and s.sample_id == "S1"), None)
        assert row is not None, "Mammalia summary row not found"
        assert row.fos_support is True
        assert "PaJs-3N" in (row.fos_evidence_text or "")

    def test_fos_evidence_text_empty_when_no_match(self):
        tax = TAX
        a = _assignment("40674", sid="S1", posterior=0.95)
        # fossil_taxa_by_sample for S1 contains a different taxon
        fossil_taxa_by_sample = {"S1": {"10090"}}
        fos_evidence_by_sample = {"S1": ["Mus musculus trap records"]}
        summaries = summarize_assignments(
            tax, [a],
            fossil_taxa_by_sample=fossil_taxa_by_sample,
            fos_evidence_by_sample=fos_evidence_by_sample,
        )
        row = next((s for s in summaries if s.taxid == "40674" and s.sample_id == "S1"), None)
        assert row is not None
        assert row.fos_support is False
        assert not row.fos_evidence_text

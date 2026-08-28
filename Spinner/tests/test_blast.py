"""Tests for BLAST result parsers using small mock TSV files.

These tests do NOT require blastn to be installed.  They parse pre-written
TSV files that mimic BLAST tabular output.
"""
from __future__ import annotations

import copy
import textwrap
from pathlib import Path


from spinner.annotation import Annotation
from spinner.clustering import parse_uchimeout
from spinner.config import DEFAULT_CONFIG
from spinner.taxonomy_blast import (
    parse_tax_blast,
    parse_tax_blast_nt_fallback,
    parse_windowed_blast,
)
from spinner.vector_screen import parse_vector_blast


def _ann(key: str, sp: str = "", gen: str = "", kingdom: str = "Animal",
         length: int = 500) -> Annotation:
    return Annotation(
        accession=key, record_key=key, source_file="", header=f">{key}",
        length=length, seq_sha256=key,
        species_guess=sp, genus_guess=gen, kingdom=kingdom,
    )


def _cfg():
    return copy.deepcopy(DEFAULT_CONFIG)


# ---------------------------------------------------------------------------
# Taxonomy BLAST parser
# ---------------------------------------------------------------------------

# outfmt: qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames
BLAST_TAX_TEMPLATE = (
    "{qid}\tNM_001234.1\t{pident}\t{length}\t{qlen}\t1\t{length}\t1e-50\t500\t9615\t{sciname}\n"
)


def write_tax_blast(tmp_path: Path, lines: list) -> str:
    p = tmp_path / "tax.tsv"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_parse_tax_blast_pass_species(tmp_path):
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="TEST.1", pident=99, length=400, qlen=500,
        sciname="Rangifer tarandus"
    )]
    path = write_tax_blast(tmp_path, lines)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert a.taxonomy_status == "PASS_SPECIES"
    assert "taxonomy_same_species" in a.reasons


def test_parse_tax_blast_pass_genus(tmp_path):
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="TEST.1", pident=95, length=400, qlen=500,
        sciname="Rangifer sp."  # genus matches but not species
    )]
    path = write_tax_blast(tmp_path, lines)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert a.taxonomy_status == "PASS_GENUS"
    assert "taxonomy_same_genus" in a.reasons


def test_parse_tax_blast_no_expected_match(tmp_path):
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="TEST.1", pident=95, length=400, qlen=500,
        sciname="Homo sapiens"
    )]
    path = write_tax_blast(tmp_path, lines)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert a.taxonomy_status == "NO_EXPECTED_MATCH"
    assert "taxonomy_no_expected_match" in a.reasons


def test_parse_tax_blast_below_min_pident(tmp_path):
    """Hits below min_pident threshold should be ignored -> taxonomy_not_checked."""
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="TEST.1", pident=60, length=400, qlen=500,  # below 70 default
        sciname="Rangifer tarandus"
    )]
    path = write_tax_blast(tmp_path, lines)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert "taxonomy_not_checked" in a.reasons


def test_parse_tax_blast_below_min_qcov(tmp_path):
    """Hits below min_qcov (50%) should be ignored."""
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="TEST.1", pident=99, length=100, qlen=500,  # qcov=20% -> below 50 default
        sciname="Rangifer tarandus"
    )]
    path = write_tax_blast(tmp_path, lines)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert "taxonomy_not_checked" in a.reasons


def test_parse_tax_blast_no_hits(tmp_path):
    """Query with no BLAST hits -> taxonomy_not_checked."""
    path = write_tax_blast(tmp_path, [])  # empty file
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert "taxonomy_not_checked" in a.reasons


def test_parse_tax_blast_missing_file_still_applies_not_checked_reason(tmp_path):
    """Real bug, found 2026-08-28 via a live palaeoscope run at Holarctic-beetle scale:
    when the taxonomy_blast subprocess itself fails to run at all (e.g. the mmseqs binary
    is missing from PATH), pipeline.py has no real output file to point parse_tax_blast()
    at. The old implementation `return`ed immediately when the path did not exist,
    skipping the per-record loop entirely -- so NOT ONE record ever got the
    "taxonomy_not_checked" reason, even though every one of them was genuinely never
    checked. Since only the REASON string (not the `taxonomy_status` field) drives
    review/reject scoring, this silently made a broken/missing tool MORE permissive than
    either a working search or an explicitly-disabled one: the review-forcing safety net
    that a length-exempted record already correctly gets (see
    test_parse_tax_blast_no_hits above) vanished along with the tool. A missing file must
    be treated identically to a genuinely-searched-but-empty result set."""
    missing_path = str(tmp_path / "does_not_exist.taxonomy_blast.tsv")
    assert not Path(missing_path).exists()
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(missing_path, ann, _cfg())
    assert "taxonomy_not_checked" in a.reasons
    assert a.taxonomy_status == "NOT_CHECKED"


def test_parse_tax_blast_missing_file_does_not_crash_with_multiple_records(tmp_path):
    """The missing-file path must scale to the real multi-record case, not just one."""
    missing_path = str(tmp_path / "does_not_exist.taxonomy_blast.tsv")
    ann = {
        "A.1": _ann("A.1", sp="Rangifer tarandus", gen="Rangifer"),
        "B.1": _ann("B.1", sp="Amara alpina", gen="Amara"),
    }
    parse_tax_blast(missing_path, ann, _cfg())
    for a in ann.values():
        assert "taxonomy_not_checked" in a.reasons


def test_parse_tax_blast_only_first_hit_used(tmp_path):
    """Only the first valid hit per query should be used."""
    lines = [
        BLAST_TAX_TEMPLATE.format(qid="TEST.1", pident=99, length=400, qlen=500,
                                   sciname="Rangifer tarandus"),
        BLAST_TAX_TEMPLATE.format(qid="TEST.1", pident=90, length=400, qlen=500,
                                   sciname="Homo sapiens"),
    ]
    path = write_tax_blast(tmp_path, lines)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert a.taxonomy_status == "PASS_SPECIES"


# ---------------------------------------------------------------------------
# NT fallback taxonomy parser (regression: Bug C)
#
# Legacy non-coding marker sequences (rRNA/tRNA genes, D-loop, spacers) have
# no ORF for the primary protein/translated search to find, so they are left
# taxonomy_not_checked forever unless a real nucleotide-level check is run.
# ---------------------------------------------------------------------------

def test_nt_fallback_resolves_record_left_not_checked_by_primary_search(tmp_path):
    """A record with no ORF hit from the primary search (still NOT_CHECKED)
    gets a real taxonomy verdict from the nucleotide fallback search."""
    a = _ann("Z54097.1", sp="Mammuthus primigenius", gen="Mammuthus", length=93)
    a.add_reason("taxonomy_not_checked")  # left over from the primary protein-mode search
    ann = {"Z54097.1": a}

    lines = [BLAST_TAX_TEMPLATE.format(
        qid="Z54097.1", pident=99, length=90, qlen=93,
        sciname="Mammuthus primigenius",
    )]
    path = write_tax_blast(tmp_path, lines)

    resolved = parse_tax_blast_nt_fallback(path, ann, _cfg(), {"Z54097.1"})

    assert resolved == 1
    assert a.taxonomy_status == "PASS_SPECIES"
    assert "taxonomy_same_species" in a.reasons
    assert "taxonomy_checked_nt_fallback" in a.reasons
    assert "taxonomy_not_checked" not in a.reasons


def test_nt_fallback_leaves_record_not_checked_when_still_no_hit(tmp_path):
    """A candidate with no hit even in the fallback search stays NOT_CHECKED."""
    a = _ann("Z54097.1", sp="Mammuthus primigenius", gen="Mammuthus", length=93)
    a.add_reason("taxonomy_not_checked")
    ann = {"Z54097.1": a}

    path = write_tax_blast(tmp_path, [])  # no hits at all

    resolved = parse_tax_blast_nt_fallback(path, ann, _cfg(), {"Z54097.1"})

    assert resolved == 0
    assert "taxonomy_not_checked" in a.reasons
    assert "taxonomy_checked_nt_fallback" not in a.reasons


def test_nt_fallback_only_touches_candidate_keys(tmp_path):
    """Records not in candidate_keys must be left completely untouched, even
    if the fallback search output happens to contain a hit for them."""
    candidate = _ann("CAND.1", sp="Mammuthus primigenius", gen="Mammuthus")
    candidate.add_reason("taxonomy_not_checked")
    other = _ann("OTHER.1", sp="Rangifer tarandus", gen="Rangifer")
    other.add_reason("taxonomy_not_checked")
    ann = {"CAND.1": candidate, "OTHER.1": other}

    lines = [
        BLAST_TAX_TEMPLATE.format(qid="CAND.1", pident=99, length=90, qlen=93,
                                   sciname="Mammuthus primigenius"),
        BLAST_TAX_TEMPLATE.format(qid="OTHER.1", pident=99, length=90, qlen=93,
                                   sciname="Rangifer tarandus"),
    ]
    path = write_tax_blast(tmp_path, lines)

    resolved = parse_tax_blast_nt_fallback(path, ann, _cfg(), {"CAND.1"})

    assert resolved == 1
    assert "taxonomy_not_checked" not in candidate.reasons
    assert "taxonomy_not_checked" in other.reasons  # untouched — not a candidate


def test_nt_fallback_below_min_pident_stays_not_checked(tmp_path):
    a = _ann("Z54097.1", sp="Mammuthus primigenius", gen="Mammuthus", length=93)
    a.add_reason("taxonomy_not_checked")
    ann = {"Z54097.1": a}

    lines = [BLAST_TAX_TEMPLATE.format(
        qid="Z54097.1", pident=60, length=90, qlen=93,  # below default 80% NT threshold
        sciname="Mammuthus primigenius",
    )]
    path = write_tax_blast(tmp_path, lines)

    resolved = parse_tax_blast_nt_fallback(path, ann, _cfg(), {"Z54097.1"})

    assert resolved == 0
    assert "taxonomy_not_checked" in a.reasons


def test_nt_fallback_cross_kingdom_hit_still_rejected(tmp_path):
    """Cross-kingdom classification logic must apply identically in the
    fallback path (shared _apply_taxonomy_hit helper)."""
    a = _ann("PLANT.1", sp="Salix arctica", gen="Salix", kingdom="Plant", length=200)
    a.add_reason("taxonomy_not_checked")
    ann = {"PLANT.1": a}

    lines = [BLAST_TAX_TEMPLATE.format(
        qid="PLANT.1", pident=95, length=190, qlen=200,
        sciname="Escherichia coli",
    )]
    # Swap in a bacterial staxid so the taxdump-based cross-kingdom check fires.
    lines = [lines[0].replace("\t9615\t", "\t562\t")]
    path = write_tax_blast(tmp_path, lines)

    class _StubTaxdb:
        def get_kingdom(self, staxid):
            return "Bacteria"

        def get_rank_name(self, staxid, rank):
            return "Escherichia" if rank == "genus" else ""

    resolved = parse_tax_blast_nt_fallback(path, ann, _cfg(), {"PLANT.1"}, taxdb=_StubTaxdb())

    assert resolved == 1
    assert a.taxonomy_status == "REJECT_CROSS_KINGDOM"
    assert "taxonomy_cross_kingdom" in a.reasons
    assert "taxonomy_checked_nt_fallback" in a.reasons


def test_nt_fallback_disabled_by_default_in_config():
    """taxonomy_blast.nt_fallback_db must default to "" (disabled) so this
    new feature has zero effect on any existing config that doesn't opt in."""
    assert DEFAULT_CONFIG["taxonomy_blast"]["nt_fallback_db"] == ""


def test_nt_fallback_no_candidates_is_noop(tmp_path):
    a = _ann("Z54097.1", sp="Mammuthus primigenius", gen="Mammuthus")
    a.add_reason("taxonomy_not_checked")
    ann = {"Z54097.1": a}
    resolved = parse_tax_blast_nt_fallback(str(tmp_path / "missing.tsv"), ann, _cfg(), set())
    assert resolved == 0
    assert "taxonomy_not_checked" in a.reasons


# ---------------------------------------------------------------------------
# Vector BLAST parser
# ---------------------------------------------------------------------------

# outfmt: qseqid sseqid pident length qstart qend sstart send evalue bitscore
VECTOR_TEMPLATE = (
    "{qid}\tUniVec:123\t{pident}\t{length}\t{qstart}\t{qend}\t1\t{length}\t1e-5\t200\n"
)


def write_vec_blast(tmp_path: Path, lines: list) -> str:
    p = tmp_path / "vec.tsv"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_parse_vector_blast_internal(tmp_path):
    # Hit at position 100-120 in a 300 bp seq -> internal (window=25, both ends safe).
    lines = [VECTOR_TEMPLATE.format(qid="TEST.1", pident=95, length=30, qstart=100, qend=130)]
    path = write_vec_blast(tmp_path, lines)
    a = _ann("TEST.1", length=300)
    ann = {"TEST.1": a}
    parse_vector_blast(path, ann, _cfg())
    assert a.vector_hit
    assert a.vector_internal
    assert not a.vector_terminal
    assert "vector_internal" in a.reasons


def test_parse_vector_blast_terminal_start(tmp_path):
    # Hit at position 5-25 in a 200 bp seq -> terminal (within window=25).
    lines = [VECTOR_TEMPLATE.format(qid="TEST.1", pident=95, length=20, qstart=5, qend=25)]
    path = write_vec_blast(tmp_path, lines)
    a = _ann("TEST.1", length=200)
    ann = {"TEST.1": a}
    parse_vector_blast(path, ann, _cfg())
    assert a.vector_terminal
    assert "vector_terminal" in a.reasons


def test_parse_vector_blast_below_min_pident(tmp_path):
    lines = [VECTOR_TEMPLATE.format(qid="TEST.1", pident=80, length=30, qstart=5, qend=35)]
    path = write_vec_blast(tmp_path, lines)
    a = _ann("TEST.1", length=200)
    ann = {"TEST.1": a}
    parse_vector_blast(path, ann, _cfg())
    assert not a.vector_hit


def test_parse_vector_blast_below_min_length(tmp_path):
    lines = [VECTOR_TEMPLATE.format(qid="TEST.1", pident=99, length=5, qstart=5, qend=10)]
    path = write_vec_blast(tmp_path, lines)
    a = _ann("TEST.1", length=200)
    ann = {"TEST.1": a}
    parse_vector_blast(path, ann, _cfg())
    assert not a.vector_hit


# ---------------------------------------------------------------------------
# vsearch UC parser
# ---------------------------------------------------------------------------

def test_parse_uc_centroid_and_member(tmp_path):
    from spinner.clustering import parse_uc

    # Format: type cluster size/id ... query target
    uc_content = textwrap.dedent("""\
        S\t0\t200\t*\t*\t*\t*\t*\tCENT.1\t*
        H\t0\t199\t99.5\t+\t0\t0\t200M\tMEM.1\tCENT.1
        C\t0\t2\t*\t*\t*\t*\t*\tCENT.1\t*
    """)
    p = tmp_path / "clusters.uc"
    p.write_text(uc_content)

    cent = _ann("CENT.1")
    mem = _ann("MEM.1")
    ann = {"CENT.1": cent, "MEM.1": mem}
    parse_uc(str(p), ann)

    assert ann["CENT.1"].cluster_role == "centroid"
    assert "cluster_representative" in ann["CENT.1"].reasons
    assert ann["MEM.1"].cluster_role == "member"
    assert "cluster_nonrepresentative" in ann["MEM.1"].reasons


def test_parse_uc_missing_file():
    from spinner.clustering import parse_uc
    a = _ann("X.1")
    ann = {"X.1": a}
    parse_uc("/no/such/file.uc", ann)  # should not raise
    assert a.cluster_role == ""


# ---------------------------------------------------------------------------
# Windowed BLAST parser
# ---------------------------------------------------------------------------

def test_parse_windowed_blast_conflict(tmp_path):
    # Two windows map to different organisms -> conflict.
    win_blast = textwrap.dedent("""\
        LONG.1|win1|1-500\tACC1\t99\t500\t500\t1\t500\t1e-50\t800\t9615\tRangifer tarandus
        LONG.1|win2|251-750\tACC2\t99\t500\t500\t1\t500\t1e-50\t800\t3702\tArabidopsis thaliana
    """)
    p = tmp_path / "windowed.tsv"
    p.write_text(win_blast)

    a = _ann("LONG.1", length=1000)
    ann = {"LONG.1": a}
    parse_windowed_blast(str(p), ann, _cfg())
    assert a.windowed_status == "WINDOWED_CONFLICT"
    assert "windowed_blast_conflict" in a.reasons


def test_parse_windowed_blast_ok(tmp_path):
    # Both windows map to same genus -> no conflict.
    win_blast = textwrap.dedent("""\
        LONG.1|win1|1-500\tACC1\t99\t500\t500\t1\t500\t1e-50\t800\t9615\tRangifer tarandus
        LONG.1|win2|251-750\tACC2\t98\t500\t500\t1\t500\t1e-50\t780\t9615\tRangifer tarandus
    """)
    p = tmp_path / "windowed.tsv"
    p.write_text(win_blast)

    a = _ann("LONG.1", length=1000)
    ann = {"LONG.1": a}
    parse_windowed_blast(str(p), ann, _cfg())
    assert a.windowed_status == "WINDOWED_OK"
    assert "windowed_blast_conflict" not in a.reasons


class _MockTaxdb:
    """Minimal taxdb stub for windowed BLAST lineage tests."""

    def __init__(self, family_map: dict):
        # taxid (int) -> family name string
        self._family = family_map

    def get_rank_name(self, taxid: int, rank: str) -> str:
        if rank == "family":
            return self._family.get(taxid, "")
        return ""


def test_parse_windowed_blast_taxdb_conflict(tmp_path):
    """With taxdb, windows from different families are flagged WINDOWED_CONFLICT."""
    win_blast = textwrap.dedent("""\
        LONG.1|win1|1-500\tACC1\t99\t500\t500\t1\t500\t1e-50\t800\t9615\tRangifer tarandus
        LONG.1|win2|251-750\tACC2\t99\t500\t500\t1\t500\t1e-50\t800\t3702\tArabidopsis thaliana
    """)
    # taxid 9615 (dog genus) -> Cervidae; taxid 3702 (Arabidopsis) -> Brassicaceae
    taxdb = _MockTaxdb({9615: "Cervidae", 3702: "Brassicaceae"})
    p = tmp_path / "windowed.tsv"
    p.write_text(win_blast)

    a = _ann("LONG.1", length=1000)
    ann = {"LONG.1": a}
    parse_windowed_blast(str(p), ann, _cfg(), taxdb=taxdb)
    assert a.windowed_status == "WINDOWED_CONFLICT"
    assert "windowed_blast_conflict" in a.reasons


def test_parse_windowed_blast_taxdb_ok_same_family(tmp_path):
    """With taxdb, two different species in the same family are NOT a conflict."""
    win_blast = textwrap.dedent("""\
        LONG.1|win1|1-500\tACC1\t99\t500\t500\t1\t500\t1e-50\t800\t9615\tRangifer tarandus
        LONG.1|win2|251-750\tACC2\t98\t500\t500\t1\t500\t1e-50\t780\t9913\tBos taurus
    """)
    # Both taxids map to Cervidae/Bovidae... use same family for test.
    taxdb = _MockTaxdb({9615: "Cervidae", 9913: "Cervidae"})
    p = tmp_path / "windowed.tsv"
    p.write_text(win_blast)

    a = _ann("LONG.1", length=1000)
    ann = {"LONG.1": a}
    parse_windowed_blast(str(p), ann, _cfg(), taxdb=taxdb)
    assert a.windowed_status == "WINDOWED_OK"
    assert "windowed_blast_conflict" not in a.reasons


def test_parse_windowed_blast_taxdb_missing_taxid_falls_back(tmp_path):
    """Windows with no taxid in taxdb fall back to genus-string comparison."""
    win_blast = textwrap.dedent("""\
        LONG.1|win1|1-500\tACC1\t99\t500\t500\t1\t500\t1e-50\t800\t0\tRangifer tarandus
        LONG.1|win2|251-750\tACC2\t99\t500\t500\t1\t500\t1e-50\t800\t0\tRangifer tarandus
    """)
    taxdb = _MockTaxdb({})  # taxid 0 -> no family
    p = tmp_path / "windowed.tsv"
    p.write_text(win_blast)

    a = _ann("LONG.1", length=1000)
    ann = {"LONG.1": a}
    parse_windowed_blast(str(p), ann, _cfg(), taxdb=taxdb)
    # Both fall back to genus "rangifer" — same genus, no conflict.
    assert a.windowed_status == "WINDOWED_OK"


def test_parse_tax_blast_mmseqs_format_compatible(tmp_path):
    """MMSeqs2 easy-search output with the configured format-output string
    is parsed correctly by parse_tax_blast() — same column layout as BLAST outfmt 6.
    """
    # MMSeqs2 columns: query target pident alnlen qlen qstart qend evalue bits taxid taxname
    mmseqs_line = (
        "TEST.1\tNM_001234.1\t98.0\t400\t500\t1\t400\t1e-50\t500\t9615\tRangifer tarandus\n"
    )
    p = tmp_path / "mmseqs.tsv"
    p.write_text(mmseqs_line)
    a = _ann("TEST.1", sp="Rangifer tarandus", gen="Rangifer")
    ann = {"TEST.1": a}
    parse_tax_blast(str(p), ann, _cfg())
    assert a.taxonomy_status == "PASS_SPECIES"
    assert "taxonomy_same_species" in a.reasons


# ---------------------------------------------------------------------------
# parse_uchimeout (chimera screen)
# ---------------------------------------------------------------------------

def _make_uchimeout_line(label: str, verdict: str) -> str:
    """Build a minimal 18-column uchimeout line for testing."""
    # Fields: score, query, parentA, parentB, top, idQM, idQA, idQB, idAB, idQT,
    #         LY, LN, LA, RY, RN, RA, div, verdict
    fields = ["0.0", label, "*", "*", "*", "0", "0", "0", "0", "0",
              "0", "0", "0", "0", "0", "0", "0.0", verdict]
    return "\t".join(fields)


def test_parse_uchimeout_chimera_detected(tmp_path):
    content = _make_uchimeout_line("CHIM.1", "Y") + "\n"
    p = tmp_path / "uchime.tsv"
    p.write_text(content)
    a = _ann("CHIM.1")
    ann = {"CHIM.1": a}
    chimeras, borderline = parse_uchimeout(str(p), ann, reject_chimeras=True, review_borderline=True)
    assert chimeras == 1
    assert borderline == 0
    assert "chimera_detected" in a.reasons


def test_parse_uchimeout_clean(tmp_path):
    content = _make_uchimeout_line("CLEAN.1", "N") + "\n"
    p = tmp_path / "uchime.tsv"
    p.write_text(content)
    a = _ann("CLEAN.1")
    ann = {"CLEAN.1": a}
    chimeras, borderline = parse_uchimeout(str(p), ann, reject_chimeras=True, review_borderline=True)
    assert chimeras == 0
    assert borderline == 0
    assert "chimera_detected" not in a.reasons
    assert "chimera_borderline" not in a.reasons


def test_parse_uchimeout_borderline(tmp_path):
    content = _make_uchimeout_line("BORDER.1", "?") + "\n"
    p = tmp_path / "uchime.tsv"
    p.write_text(content)
    a = _ann("BORDER.1")
    ann = {"BORDER.1": a}
    chimeras, borderline = parse_uchimeout(str(p), ann, reject_chimeras=True, review_borderline=True)
    assert chimeras == 0
    assert borderline == 1
    assert "chimera_borderline" in a.reasons


def test_parse_uchimeout_borderline_suppressed(tmp_path):
    """With review_borderline=False, ? verdicts should be ignored."""
    content = _make_uchimeout_line("BORDER.1", "?") + "\n"
    p = tmp_path / "uchime.tsv"
    p.write_text(content)
    a = _ann("BORDER.1")
    ann = {"BORDER.1": a}
    chimeras, borderline = parse_uchimeout(str(p), ann, reject_chimeras=True, review_borderline=False)
    assert borderline == 0
    assert "chimera_borderline" not in a.reasons


def test_parse_uchimeout_mixed(tmp_path):
    content = (
        _make_uchimeout_line("CHIM.1", "Y") + "\n"
        + _make_uchimeout_line("CLEAN.1", "N") + "\n"
        + _make_uchimeout_line("BORDER.1", "?") + "\n"
    )
    p = tmp_path / "uchime.tsv"
    p.write_text(content)
    a_chim = _ann("CHIM.1")
    a_clean = _ann("CLEAN.1")
    a_border = _ann("BORDER.1")
    ann = {"CHIM.1": a_chim, "CLEAN.1": a_clean, "BORDER.1": a_border}
    chimeras, borderline = parse_uchimeout(str(p), ann)
    assert chimeras == 1
    assert borderline == 1
    assert "chimera_detected" in a_chim.reasons
    assert "chimera_borderline" in a_border.reasons
    assert not a_clean.reasons

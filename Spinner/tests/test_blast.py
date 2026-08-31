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
    mark_length_exempt_records,
    parse_tax_blast,
    parse_tax_blast_escalation_genus_species,
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
    that a genuinely-searched-but-empty result set already correctly gets (see
    test_parse_tax_blast_no_hits above) vanished along with the tool. A missing file must
    be treated identically to a genuinely-searched-but-empty result set (records already
    tagged `taxonomy_exempt_length` by the caller are the one deliberate exception --
    see test_parse_tax_blast_exempt_length_records_skip_not_checked_even_with_no_output_file
    below)."""
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


def test_parse_tax_blast_exempt_length_records_skip_not_checked(tmp_path):
    """Real bug, found 2026-08-29 via forensic analysis of a real complete-mitogenome
    cluster centroid (EU153454.1, 16,480 bp): records excluded from the search entirely
    by `taxonomy_blast.max_query_length` (full organelle genomes -- see the
    long-standing config comment "full organelles don't need kingdom verification") used
    to get the exact same "taxonomy_not_checked" reason a genuinely-searched-but-no-hit
    record gets. Since "taxonomy_not_checked" is unconditionally in
    decision_rules.review_reasons, this permanently blocked otherwise-clean, high-scoring
    complete genomes from ever reaching auto-KEEP. The pipeline.py caller now tags these
    records `taxonomy_exempt_length` BEFORE calling parse_tax_blast(); this function must
    recognize that tag and skip them, never adding "taxonomy_not_checked" on top."""
    path = write_tax_blast(tmp_path, [])  # empty file -- record was never submitted
    a = _ann("EU153454.1", sp="Mammuthus primigenius", gen="Mammuthus", length=16480)
    a.taxonomy_status = "EXEMPT_LENGTH"
    a.add_reason("taxonomy_exempt_length")
    ann = {"EU153454.1": a}
    parse_tax_blast(path, ann, _cfg())
    assert "taxonomy_not_checked" not in a.reasons
    assert a.reasons == ["taxonomy_exempt_length"]
    assert a.taxonomy_status == "EXEMPT_LENGTH"


def test_parse_tax_blast_exempt_length_records_skip_not_checked_even_with_no_output_file(tmp_path):
    """Same guarantee as the test above, but for the tool-failure path (no output file
    at all, e.g. mmseqs missing from PATH) -- must not regress to double-tagging exempt
    records just because the search subprocess itself never ran."""
    missing_path = str(tmp_path / "does_not_exist.taxonomy_blast.tsv")
    assert not Path(missing_path).exists()
    a = _ann("EU153454.1", sp="Mammuthus primigenius", gen="Mammuthus", length=16480)
    a.taxonomy_status = "EXEMPT_LENGTH"
    a.add_reason("taxonomy_exempt_length")
    ann = {"EU153454.1": a}
    parse_tax_blast(missing_path, ann, _cfg())
    assert "taxonomy_not_checked" not in a.reasons
    assert a.reasons == ["taxonomy_exempt_length"]


def test_parse_tax_blast_exempt_length_does_not_affect_sibling_records(tmp_path):
    """The exemption must be strictly per-record: a normal-length sibling record with a
    genuine no-hit result must still correctly get "taxonomy_not_checked"."""
    path = write_tax_blast(tmp_path, [])
    exempt = _ann("EU153454.1", sp="Mammuthus primigenius", gen="Mammuthus", length=16480)
    exempt.taxonomy_status = "EXEMPT_LENGTH"
    exempt.add_reason("taxonomy_exempt_length")
    normal = _ann("SHORT.1", sp="Mammuthus primigenius", gen="Mammuthus", length=400)
    ann = {"EU153454.1": exempt, "SHORT.1": normal}
    parse_tax_blast(path, ann, _cfg())
    assert "taxonomy_not_checked" not in exempt.reasons
    assert "taxonomy_not_checked" in normal.reasons


def test_taxonomy_exempt_length_is_neutral_score_and_not_a_review_reason():
    """Config-level guarantee the fix depends on: `taxonomy_exempt_length` must score 0
    (a confident pass, not a bonus or penalty) and must never appear in
    decision_rules.review_reasons or hard_reject_reasons -- otherwise it would still
    force review/reject exactly like the bug it replaces."""
    cfg = _cfg()
    assert cfg["scoring"]["taxonomy_exempt_length"] == 0
    assert "taxonomy_exempt_length" not in cfg["decision_rules"]["review_reasons"]
    assert "taxonomy_exempt_length" not in cfg["decision_rules"]["hard_reject_reasons"]


class _FakeRecord:
    """Minimal stand-in for a FastaRecord -- mark_length_exempt_records only reads .id."""
    def __init__(self, rid: str):
        self.id = rid


def test_mark_length_exempt_records_tags_mito_and_plastid_only():
    """Pinned real scenario: EU153454.1 (Mito, 16,480 bp) should be tagged; a long
    NucMark/Other-class record that merely happens to cross max_query_length (e.g. a
    large nuclear scaffold) has no "can't plausibly be a contaminant" guarantee and must
    NOT be tagged -- it keeps the normal taxonomy_not_checked review-forcing path."""
    mito = _ann("EU153454.1", sp="Mammuthus primigenius", gen="Mammuthus", length=16480)
    mito.marker_class = "Mito"
    plastid = _ann("PLASTID_BIG.1", length=20000)
    plastid.marker_class = "Plastid"
    nucmark = _ann("NUCMARK_BIG.1", length=15000)
    nucmark.marker_class = "NucMark"
    other = _ann("OTHER_BIG.1", length=15000)
    other.marker_class = "Other"
    ann = {a.record_key: a for a in (mito, plastid, nucmark, other)}
    records = [_FakeRecord(k) for k in ann]

    n_organelle, n_other = mark_length_exempt_records(records, ann)

    assert n_organelle == 2
    assert n_other == 2
    assert "taxonomy_exempt_length" in mito.reasons
    assert mito.taxonomy_status == "EXEMPT_LENGTH"
    assert "taxonomy_exempt_length" in plastid.reasons
    assert "taxonomy_exempt_length" not in nucmark.reasons
    assert "taxonomy_exempt_length" not in other.reasons
    # NucMark/Other records are untouched here -- parse_tax_blast() applies
    # taxonomy_not_checked to them itself, since they're absent from `hits`.
    assert nucmark.taxonomy_status == "NOT_CHECKED"
    assert other.taxonomy_status == "NOT_CHECKED"


def test_mark_length_exempt_records_empty_input():
    assert mark_length_exempt_records([], {}) == (0, 0)


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
# NR-protein escalation for taxonomy_no_expected_match (genus/species win-condition)
# ---------------------------------------------------------------------------

class _GenusStubTaxdb:
    """Minimal taxdb stub: staxid 9615 -> genus 'Canis', anything else -> ''."""

    def get_rank_name(self, staxid, rank):
        if rank == "genus" and staxid == 9615:
            return "Canis"
        return ""


def test_escalation_no_expected_match_disabled_by_default_in_config():
    assert DEFAULT_CONFIG["taxonomy_blast"]["escalate_no_expected_match"] is False


def test_escalation_no_expected_match_rescues_on_species_string_match(tmp_path):
    a = _ann("WOLF.1", sp="Canis lupus", gen="Canis")
    a.add_reason("taxonomy_no_expected_match")
    ann = {"WOLF.1": a}
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="WOLF.1", pident=95, length=190, qlen=200, sciname="Canis lupus",
    )]
    lines = [lines[0].replace("\t9615\t", "\tN/A\t")]  # no taxdump staxid -> string fallback
    path = write_tax_blast(tmp_path, lines)

    rescued = parse_tax_blast_escalation_genus_species(path, ann, _cfg(), None, "nr_no_expected_match")

    assert rescued == 1
    assert a.taxonomy_status == "PASS_SPECIES"
    assert "taxonomy_no_expected_match" not in a.reasons
    assert "taxonomy_same_species" in a.reasons
    assert "taxonomy_rescued_nr_no_expected_match" in a.reasons


def test_escalation_no_expected_match_rescues_on_genus_via_taxdb(tmp_path):
    a = _ann("DOG.1", sp="Canis familiaris", gen="Canis")
    a.add_reason("taxonomy_no_expected_match")
    ann = {"DOG.1": a}
    # Sciname deliberately does NOT contain "Canis" as a substring, so only the
    # taxdb-based genus lookup (not the string fallback) can rescue this record.
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="DOG.1", pident=95, length=190, qlen=200, sciname="Vulpes vulpes",
    )]
    path = write_tax_blast(tmp_path, lines)  # staxid 9615 from the shared template

    rescued = parse_tax_blast_escalation_genus_species(
        path, ann, _cfg(), _GenusStubTaxdb(), "nr_no_expected_match"
    )

    assert rescued == 1
    assert a.taxonomy_status == "PASS_GENUS"
    assert "taxonomy_no_expected_match" not in a.reasons
    assert "taxonomy_same_genus" in a.reasons


def test_escalation_no_expected_match_no_hit_leaves_record_untouched(tmp_path):
    a = _ann("NOHIT.1", sp="Canis lupus", gen="Canis")
    a.add_reason("taxonomy_no_expected_match")
    ann = {"NOHIT.1": a}

    rescued = parse_tax_blast_escalation_genus_species(
        str(tmp_path / "missing.tsv"), ann, _cfg(), None, "nr_no_expected_match"
    )

    assert rescued == 0
    assert "taxonomy_no_expected_match" in a.reasons
    assert a.taxonomy_status == "NOT_CHECKED"


def test_escalation_no_expected_match_wrong_taxon_not_rescued(tmp_path):
    """A hit that shares neither genus nor species must NOT rescue the record --
    this is the whole point of the stricter genus/species win-condition vs. the
    cross-kingdom escalation's same-kingdom-is-enough check."""
    a = _ann("CAT.1", sp="Felis catus", gen="Felis")
    a.add_reason("taxonomy_no_expected_match")
    ann = {"CAT.1": a}
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="CAT.1", pident=95, length=190, qlen=200, sciname="Canis lupus",
    )]
    path = write_tax_blast(tmp_path, lines)

    rescued = parse_tax_blast_escalation_genus_species(
        path, ann, _cfg(), _GenusStubTaxdb(), "nr_no_expected_match"
    )

    assert rescued == 0
    assert "taxonomy_no_expected_match" in a.reasons


def test_escalation_no_expected_match_ignores_records_without_that_reason(tmp_path):
    """Only records currently flagged taxonomy_no_expected_match are candidates --
    a record that already passed (or was never checked) must not be touched even if
    it happens to appear in the escalation search output."""
    a = _ann("PASSED.1", sp="Canis lupus", gen="Canis")
    a.taxonomy_status = "PASS_SPECIES"
    a.add_reason("taxonomy_same_species")
    ann = {"PASSED.1": a}
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="PASSED.1", pident=95, length=190, qlen=200, sciname="Canis lupus",
    )]
    path = write_tax_blast(tmp_path, lines)

    rescued = parse_tax_blast_escalation_genus_species(
        path, ann, _cfg(), None, "nr_no_expected_match"
    )

    assert rescued == 0
    assert a.taxonomy_status == "PASS_SPECIES"


def test_escalation_no_expected_match_below_qcov_threshold_not_rescued(tmp_path):
    a = _ann("LOWCOV.1", sp="Canis lupus", gen="Canis")
    a.add_reason("taxonomy_no_expected_match")
    ann = {"LOWCOV.1": a}
    # length/qlen = 25% query coverage, below the default 50% min_qcov filter.
    lines = [BLAST_TAX_TEMPLATE.format(
        qid="LOWCOV.1", pident=95, length=50, qlen=200, sciname="Canis lupus",
    )]
    path = write_tax_blast(tmp_path, lines)

    rescued = parse_tax_blast_escalation_genus_species(
        path, ann, _cfg(), None, "nr_no_expected_match"
    )

    assert rescued == 0
    assert "taxonomy_no_expected_match" in a.reasons


def test_escalation_rescued_no_expected_match_score_key_is_neutral():
    """The provenance reason itself must score 0 -- the real score comes from
    taxonomy_same_genus/taxonomy_same_species, not double-counted here."""
    assert DEFAULT_CONFIG["scoring"]["taxonomy_rescued_nr_no_expected_match"] == 0


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

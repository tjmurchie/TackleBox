"""Tests for spinner.taxonomy_blast.TaxdumpDB's synonym-aware name resolution and
canonicalize_species_guesses().

Real bug, found 2026-09-01 via a live ~10,000-species PalaeoSCOPE Phase B run's
candidate list: real taxonomic synonym pairs (e.g. Antigone canadensis / Grus
canadensis, a genus reassignment) parse to two DIFFERENT species_guess strings from
real GenBank headers even though they are the same species -- confirmed live against
the real NCBI taxdump that both names already resolve to the same taxid (1977160) via
its "synonym" name class, which TaxdumpDB previously never indexed at all (only
"scientific name" was kept). These fixtures use the real taxids for both on-record
example pairs (Antigone/Grus canadensis = 1977160; Bison bison/Bos bison = 9901) so the
tests are grounded in real data, not fabricated IDs.

No real taxdump files are used -- small synthetic nodes.dmp/names.dmp fixtures are
written per test, matching the real NCBI dump `field|field|...` shape closely enough
for TaxdumpDB._load()'s own parser (it only strips and splits on "|", so exact tab
padding doesn't matter).
"""
from __future__ import annotations

import argparse
import csv

from spinner.annotation import Annotation
from spinner.pipeline import run_pipeline
from spinner.taxonomy_blast import TaxdumpDB, canonicalize_species_guesses


def _write_taxdump(tmp_path, nodes_lines, names_lines):
    d = tmp_path / "taxdump"
    d.mkdir()
    (d / "nodes.dmp").write_text("\n".join(nodes_lines) + "\n", encoding="utf-8")
    (d / "names.dmp").write_text("\n".join(names_lines) + "\n", encoding="utf-8")
    return str(d)


def _real_synonym_taxdump(tmp_path):
    """A tiny synthetic taxdump covering both real on-record synonym pairs plus one
    unrelated species with no synonym, for contrast."""
    nodes = [
        "1977160|8496|species|",   # Antigone canadensis
        "9901|9895|species|",      # Bison bison
        "9606|9605|species|",      # Homo sapiens -- no synonym, contrast case
    ]
    names = [
        "1977160|Antigone canadensis||scientific name|",
        "1977160|Grus canadensis||synonym|",
        "9901|Bison bison||scientific name|",
        "9901|Bos bison||synonym|",
        "9606|Homo sapiens||scientific name|",
    ]
    return _write_taxdump(tmp_path, nodes, names)


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


class TestTaxdumpDBSynonymIndex:
    def test_scientific_name_and_synonym_resolve_to_same_canonical_name(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.canonical_species_name("Antigone canadensis") == "Antigone canadensis"
        assert taxdb.canonical_species_name("Grus canadensis") == "Antigone canadensis"

    def test_second_real_pair_bison_bos(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.canonical_species_name("Bison bison") == "Bison bison"
        assert taxdb.canonical_species_name("Bos bison") == "Bison bison"

    def test_case_insensitive_resolution(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.canonical_species_name("grus CANADENSIS") == "Antigone canadensis"

    def test_unknown_name_returned_unchanged(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.canonical_species_name("Nonexistent species") == "Nonexistent species"

    def test_name_with_no_synonym_returned_unchanged(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.canonical_species_name("Homo sapiens") == "Homo sapiens"

    def test_empty_string_returned_unchanged(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.canonical_species_name("") == ""

    def test_name_to_taxid_covers_both_classes(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        assert taxdb.name_to_taxid["antigone canadensis"] == 1977160
        assert taxdb.name_to_taxid["grus canadensis"] == 1977160


class TestCanonicalizeSpeciesGuesses:
    def test_synonym_pair_collapses_to_one_species_guess(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        ann = {
            "rec1": _ann("rec1", "Antigone canadensis"),
            "rec2": _ann("rec2", "Grus canadensis"),
        }
        changed = canonicalize_species_guesses(ann, taxdb)
        assert changed == 1  # only rec2's guess actually changes
        assert ann["rec1"].species_guess == "Antigone canadensis"
        assert ann["rec2"].species_guess == "Antigone canadensis"

    def test_grouping_key_now_identical_after_canonicalization(self, tmp_path):
        """The actual point of the fix: cluster.by/capping/rescue all group by raw
        species_guess text, so after canonicalization the two records must produce the
        SAME grouping key -- this is what lets them be treated as one species pool
        downstream instead of two independently-clustered/independently-capped ones."""
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        ann = {
            "rec1": _ann("rec1", "Antigone canadensis", marker_class="Mito"),
            "rec2": _ann("rec2", "Grus canadensis", marker_class="Mito"),
        }
        canonicalize_species_guesses(ann, taxdb)
        key1 = (ann["rec1"].species_guess, ann["rec1"].marker_class)
        key2 = (ann["rec2"].species_guess, ann["rec2"].marker_class)
        assert key1 == key2

    def test_unresolvable_species_guess_left_untouched(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        ann = {"rec1": _ann("rec1", "Totally Unknown species")}
        changed = canonicalize_species_guesses(ann, taxdb)
        assert changed == 0
        assert ann["rec1"].species_guess == "Totally Unknown species"

    def test_empty_species_guess_skipped(self, tmp_path):
        taxdb = TaxdumpDB(_real_synonym_taxdump(tmp_path))
        ann = {"rec1": _ann("rec1", "")}
        changed = canonicalize_species_guesses(ann, taxdb)
        assert changed == 0
        assert ann["rec1"].species_guess == ""


class TestPipelineWiring:
    """Real end-to-end proof that pipeline.py's wiring (early taxdb load, gated on
    cluster.canonicalize_species_guess, calling canonicalize_species_guesses before
    decisions/output are written) actually takes effect -- not just the unit-level
    functions in isolation above. No vsearch/blastn involved (cluster/taxonomy_blast
    steps stay disabled): this only proves the canonicalization wiring itself fires."""

    def _run(self, tmp_path, canonicalize: bool):
        import yaml

        fasta = tmp_path / "in.fasta"
        fasta.write_text(
            ">REC1.1 Antigone canadensis mitochondrion, complete genome\n"
            + "ACGT" * 30 + "\n"
            + ">REC2.1 Grus canadensis mitochondrion, complete genome\n"
            + "TGCA" * 30 + "\n",
            encoding="utf-8",
        )
        config = tmp_path / "config.yml"
        config.write_text(
            yaml.safe_dump({
                "cluster": {"canonicalize_species_guess": canonicalize},
                "taxonomy_blast": {"taxdump_dir": _real_synonym_taxdump(tmp_path)},
            }),
            encoding="utf-8",
        )
        outprefix = str(tmp_path / "out")
        args = argparse.Namespace(
            fasta=[str(fasta)], outprefix=outprefix, config=str(config),
            species_kingdom="", regions_config="", adapters="", bad_keywords="",
            keep_temp=False, taxonomy_blast_db="", vector_blast_db="",
            windowed_blast_db="", enable_cluster=False,
        )
        run_pipeline(args, filter_mode=True)
        with open(outprefix + ".decisions.tsv", encoding="utf-8") as f:
            rows = {r["record_key"]: r for r in csv.DictReader(f, delimiter="\t")}
        return rows["REC1.1"]["species_guess"], rows["REC2.1"]["species_guess"]

    def test_canonicalize_enabled_merges_synonym_pair(self, tmp_path):
        sp1, sp2 = self._run(tmp_path, canonicalize=True)
        assert sp1 == "Antigone canadensis"
        assert sp2 == "Antigone canadensis"

    def test_canonicalize_disabled_keeps_original_header_names(self, tmp_path):
        """Default (opt-in flag off) behavior must be completely unchanged."""
        sp1, sp2 = self._run(tmp_path, canonicalize=False)
        assert sp1 == "Antigone canadensis"
        assert sp2 == "Grus canadensis"

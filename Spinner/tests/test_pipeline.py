"""Integration tests: run the pipeline on toy data and verify output."""
from __future__ import annotations

import argparse
import copy
import csv
from pathlib import Path

import pytest

from spinner.config import DEFAULT_CONFIG
from spinner.fasta import parse_fasta
from spinner.pipeline import run_pipeline


def _args(fasta_paths, outprefix, config="", species_kingdom="",
          regions_config="", adapters="", bad_keywords="", keep_temp=False):
    """Create a minimal argparse.Namespace suitable for run_pipeline()."""
    return argparse.Namespace(
        fasta=fasta_paths,
        outprefix=outprefix,
        config=config,
        species_kingdom=species_kingdom,
        regions_config=regions_config,
        adapters=adapters,
        bad_keywords=bad_keywords,
        keep_temp=keep_temp,
        # These are only set by screen-* subcommands:
        taxonomy_blast_db="",
        vector_blast_db="",
        windowed_blast_db="",
        enable_cluster=False,
    )


# ---------------------------------------------------------------------------
# Minimal filter run
# ---------------------------------------------------------------------------

def test_minimal_filter_run(
    tmp_path, fasta_file, adapters_file, bad_keywords_file,
    species_kingdom_file, regions_file
):
    """filter subcommand should write decisions.tsv and split FASTAs."""
    outprefix = str(tmp_path / "out")
    args = _args(
        [str(fasta_file)],
        outprefix,
        species_kingdom=str(species_kingdom_file),
        adapters=str(adapters_file),
        bad_keywords=str(bad_keywords_file),
        regions_config=str(regions_file),
    )
    run_pipeline(args, filter_mode=True)

    assert Path(outprefix + ".decisions.tsv").exists()
    assert Path(outprefix + ".keep.fasta").exists()
    assert Path(outprefix + ".review.fasta").exists()
    assert Path(outprefix + ".reject.fasta").exists()
    assert Path(outprefix + ".summary.tsv").exists()
    assert Path(outprefix + ".summary.html").exists()


def test_command_txt_written(tmp_path, fasta_file):
    outprefix = str(tmp_path / "out")
    args = _args([str(fasta_file)], outprefix)
    run_pipeline(args, filter_mode=False)
    assert Path(outprefix + ".command.txt").exists()


def test_resolved_config_written(tmp_path, fasta_file):
    outprefix = str(tmp_path / "out")
    args = _args([str(fasta_file)], outprefix)
    run_pipeline(args, filter_mode=False)
    assert Path(outprefix + ".run_config.resolved.yml").exists()


# ---------------------------------------------------------------------------
# All records appear in decisions.tsv
# ---------------------------------------------------------------------------

def test_all_records_in_decisions(
    tmp_path, fasta_file, adapters_file, bad_keywords_file,
    species_kingdom_file, regions_file
):
    """Every input record must have exactly one decisions row."""
    input_recs = parse_fasta([str(fasta_file)])
    outprefix = str(tmp_path / "out")
    args = _args(
        [str(fasta_file)], outprefix,
        species_kingdom=str(species_kingdom_file),
        adapters=str(adapters_file),
        bad_keywords=str(bad_keywords_file),
        regions_config=str(regions_file),
    )
    run_pipeline(args, filter_mode=False)

    with open(outprefix + ".decisions.tsv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    assert len(rows) == len(input_recs)


def test_duplicate_accession_both_in_decisions(tmp_path, adapters_file, bad_keywords_file):
    """Both records with the same accession must appear in decisions.tsv."""
    fasta_content = (
        ">DUP.1 Betula papyrifera chloroplast rbcL gene\n" + "ACGT" * 20 + "\n"
        ">DUP.1 Betula papyrifera chloroplast rbcL gene duplicate\n" + "GGGG" * 20 + "\n"
    )
    p = tmp_path / "dup.fasta"
    p.write_text(fasta_content)
    outprefix = str(tmp_path / "out")
    args = _args([str(p)], outprefix, adapters=str(adapters_file),
                 bad_keywords=str(bad_keywords_file))
    run_pipeline(args, filter_mode=False)

    with open(outprefix + ".decisions.tsv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    assert len(rows) == 2
    accessions = [r["accession"] for r in rows]
    assert accessions.count("DUP.1") == 2
    # The second copy should be flagged.
    dup_rows = [r for r in rows if r.get("duplicate_accession") == "True"]
    assert len(dup_rows) == 1


# ---------------------------------------------------------------------------
# Specific detection tests
# ---------------------------------------------------------------------------

def test_adapter_internal_gets_reject(tmp_path):
    """A record with an internal adapter hit should be REJECT."""
    # Internal adapter: position 120, well inside both terminal windows.
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    seq = "A" * 120 + adapter + "A" * 120
    p = tmp_path / "t.fasta"
    p.write_text(f">INT.1 Rangifer tarandus mitochondrion\n{seq}\n")
    ad_file = tmp_path / "ad.tsv"
    ad_file.write_text(f"name\tsequence\tmax_mismatch\taction\nTruSeq\t{adapter}\t0\treject\n")
    outprefix = str(tmp_path / "out")
    args = _args([str(p)], outprefix, adapters=str(ad_file))
    run_pipeline(args, filter_mode=False)

    with open(outprefix + ".decisions.tsv") as f:
        rows = {r["record_key"]: r for r in csv.DictReader(f, delimiter="\t")}

    assert rows["INT.1"]["decision"] == "REJECT"
    assert "adapter_internal" in rows["INT.1"]["reasons"]


def test_high_n_gets_penalty_not_rejected(tmp_path):
    """High-N records receive a score penalty but are NOT hard-rejected.
    A high-N sequence from a rare/ancient species is better than no reference.
    Score: 100 - 20 (n_fraction_high) = 80 >= keep_min(65) -> KEEP.
    """
    p = tmp_path / "t.fasta"
    p.write_text(">HIGHN.1 Salix arctica ITS\n" + "N" * 80 + "\n")
    outprefix = str(tmp_path / "out")
    args = _args([str(p)], outprefix)
    run_pipeline(args, filter_mode=False)

    with open(outprefix + ".decisions.tsv") as f:
        rows = {r["record_key"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert "n_fraction_high" in rows["HIGHN.1"]["reasons"]
    assert rows["HIGHN.1"]["decision"] == "KEEP"


def test_bad_keyword_reject(tmp_path, bad_keywords_file):
    p = tmp_path / "t.fasta"
    p.write_text(">VEC.1 synthetic construct vector pUC19\n" + "ACGT" * 25 + "\n")
    outprefix = str(tmp_path / "out")
    args = _args([str(p)], outprefix, bad_keywords=str(bad_keywords_file))
    run_pipeline(args, filter_mode=False)

    with open(outprefix + ".decisions.tsv") as f:
        rows = {r["record_key"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["VEC.1"]["decision"] == "REJECT"


# ---------------------------------------------------------------------------
# explain subcommand
# ---------------------------------------------------------------------------

def test_explain_found(tmp_path, fasta_file, adapters_file, bad_keywords_file, capsys):
    outprefix = str(tmp_path / "out")
    args = _args([str(fasta_file)], outprefix, adapters=str(adapters_file),
                 bad_keywords=str(bad_keywords_file))
    run_pipeline(args, filter_mode=False)

    from spinner.cli import _cmd_explain
    explain_args = argparse.Namespace(
        decisions=outprefix + ".decisions.tsv",
        accession="NC_CLEAN.1",
    )
    _cmd_explain(explain_args)
    captured = capsys.readouterr()
    assert "NC_CLEAN.1" in captured.out
    assert "Decision" in captured.out


def test_explain_not_found(tmp_path, fasta_file, capsys):
    outprefix = str(tmp_path / "out")
    args = _args([str(fasta_file)], outprefix)
    run_pipeline(args, filter_mode=False)

    from spinner.cli import _cmd_explain
    explain_args = argparse.Namespace(
        decisions=outprefix + ".decisions.tsv",
        accession="DOES_NOT_EXIST.999",
    )
    with pytest.raises(SystemExit):
        _cmd_explain(explain_args)


# ---------------------------------------------------------------------------
# report subcommand
# ---------------------------------------------------------------------------

def test_report_subcommand(tmp_path, fasta_file):
    outprefix = str(tmp_path / "out")
    args = _args([str(fasta_file)], outprefix)
    run_pipeline(args, filter_mode=False)

    from spinner.reporting import report_from_decisions
    report_prefix = str(tmp_path / "report")
    report_from_decisions(outprefix + ".decisions.tsv", report_prefix)
    assert Path(report_prefix + ".summary.tsv").exists()
    assert Path(report_prefix + ".summary.html").exists()

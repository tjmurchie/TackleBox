"""Regression tests for NCBI-NT_Downloader.pl's --enable-regions/--disable-regions
handling of the organelle (MITO/PLASTID) block.

Background: the script's header comment documents --enable-regions/--disable-regions
as working generically for any regions_config.tsv region_id, including MITO and
PLASTID. Before this fix, MITO/PLASTID had an empty ncbi_title_clause column, so
load_marker_clauses_from_regions() always skipped them regardless of
--enable-regions/--disable-regions; the actual mitochondrial/plastid search terms
came from a separate hardcoded $organelle_block string that never consulted the
CLI overrides at all -- i.e. --enable-regions/--disable-regions had zero observable
effect on organelle terms.

NCBI-NT_Downloader.pl has no unit-test framework of its own (it's a standalone Perl
script, not a module), so these tests exercise the real, unmodified script via its
--dry-run flag (added alongside this fix) as a real subprocess invocation -- the
script performs its normal ESearch term construction and prints it instead of
querying NCBI, so this validates the actual code path a live run takes without
requiring network access.
"""
import os
import pathlib
import shutil
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "NCBI-NT_Downloader.pl"


def _find_working_perl():
    """Return a perl binary that actually has Bio::DB::EUtilities available, or
    None if none can be found. NCBI-NT_Downloader.pl `use`s that module at
    compile time (even under --dry-run), and on this host it typically only
    exists in the dedicated `ncbi-perl_flyguide` conda env, not on the bare
    PATH -- see palaeoscope/adapters/flyguide.py's perl_bin handling for the
    same environment-resolution problem from the PalaeoSCOPE side.
    """
    candidates = []
    env_override = os.environ.get("FLYGUIDE_PERL_BIN")
    if env_override:
        candidates.append(env_override)
    candidates.append("perl")

    conda_exe = os.environ.get("CONDA_EXE")
    conda_roots = []
    if conda_exe:
        conda_roots.append(pathlib.Path(conda_exe).resolve().parents[1])
    conda_roots.append(pathlib.Path.home() / "miniforge3")
    conda_roots.append(pathlib.Path.home() / "miniconda3")
    for conda_root in conda_roots:
        candidates.append(str(conda_root / "envs" / "ncbi-perl_flyguide" / "bin" / "perl"))

    for candidate in candidates:
        perl_path = shutil.which(candidate) or (candidate if os.path.isfile(candidate) else None)
        if not perl_path:
            continue
        probe = subprocess.run(
            [perl_path, "-e", "use Bio::DB::EUtilities; print 'ok';"],
            capture_output=True, text=True,
        )
        if probe.returncode == 0 and probe.stdout.strip() == "ok":
            return perl_path
    return None


PERL_BIN = _find_working_perl()

pytestmark = pytest.mark.skipif(
    PERL_BIN is None,
    reason=(
        "No perl with Bio::DB::EUtilities found (expected in the ncbi-perl_flyguide "
        "conda env, or set FLYGUIDE_PERL_BIN to a working perl binary)"
    ),
)


def _dry_run_term(tmp_path, species, extra_args=()):
    species_file = tmp_path / "species.txt"
    species_file.write_text(species + "\n", encoding="utf-8")
    out_prefix = tmp_path / "out"
    cmd = [
        PERL_BIN, str(SCRIPT),
        "--dry-run", "--no-tui",
        *extra_args,
        str(species_file), str(out_prefix), "test@example.com",
    ]
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    line = result.stdout.strip().splitlines()[0]
    _, term = line.split("\t", 1)
    return term


def test_default_organelle_block_matches_pre_fix_hardcoded_string(tmp_path):
    """No --enable-regions/--disable-regions given: must reproduce the exact
    literal string the script always used to emit unconditionally, so this fix
    is behavior-preserving for every existing invocation that doesn't pass
    either flag."""
    term = _dry_run_term(tmp_path, "Bison bison")
    assert "(mitochondrion OR mitochondrial OR chloroplast OR plastid)" in term


def test_disable_regions_plastid_removes_plastid_terms_only(tmp_path):
    term = _dry_run_term(tmp_path, "Bison bison", ["--disable-regions", "PLASTID"])
    assert "chloroplast" not in term
    assert "plastid" not in term
    assert "(mitochondrion OR mitochondrial)" in term


def test_disable_regions_mito_removes_mito_terms_only(tmp_path):
    term = _dry_run_term(tmp_path, "Bison bison", ["--disable-regions", "MITO"])
    assert "mitochondrion" not in term
    assert "mitochondrial" not in term
    assert "(chloroplast OR plastid)" in term


def test_no_organelle_removes_both_mito_and_plastid_terms(tmp_path):
    term = _dry_run_term(tmp_path, "Bison bison", ["--no-organelle"])
    for word in ("mitochondrion", "mitochondrial", "chloroplast", "plastid"):
        assert word not in term


def test_enable_regions_mito_overrides_no_organelle(tmp_path):
    """--enable-regions MITO must force mitochondrial terms back into the query
    even when --no-organelle is also given -- the CLI-override semantics
    documented in the script's own header comment (force_enable takes
    priority over any other setting, including --no-organelle)."""
    term = _dry_run_term(
        tmp_path, "Bison bison", ["--no-organelle", "--enable-regions", "MITO"]
    )
    assert "(mitochondrion OR mitochondrial)" in term
    assert "chloroplast" not in term
    assert "plastid" not in term


def test_disable_regions_both_mito_and_plastid_matches_no_organelle(tmp_path):
    disabled_term = _dry_run_term(
        tmp_path, "Bison bison", ["--disable-regions", "MITO,PLASTID"]
    )
    no_organelle_term = _dry_run_term(tmp_path, "Bison bison", ["--no-organelle"])
    assert disabled_term == no_organelle_term

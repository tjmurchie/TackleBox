#!/usr/bin/env python3
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "neotoma_extinct_to_gbif.py"
DATA = ROOT / "tests" / "data"

def run_cmd(args, cwd):
    result = subprocess.run(args, cwd=cwd, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr + result.stdout
    return result

def test_offline_extinct_animals_binomial(tmp_path):
    out = tmp_path / "out.csv"
    run_cmd([
        sys.executable, str(SCRIPT),
        "--offline-taxa-json", str(DATA / "neotoma_taxa_fixture.json"),
        "--offline-occurrences-json", str(DATA / "neotoma_occurrences_fixture.json"),
        "--region", "northern_hemisphere",
        "--period", "quaternary",
        "--organisms", "animals",
        "--status", "extinct",
        "--out", str(out),
        "--write-flyguide-files",
        "--out-prefix", str(tmp_path / "out")
    ], cwd=ROOT)
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    names = {r["species"] for r in rows}
    assert "Mammuthus primigenius" in names
    # antiquus collapses to the binomial bucket by default for broad NCBI retrievability
    assert "Bison bison" in names
    assert "Canis" not in names
    assert (tmp_path / "out_species_search.txt").exists()
    assert (tmp_path / "out_species_kingdom.tsv").exists()
    summary = json.load(open(tmp_path / "out.summary.json", encoding="utf-8"))
    assert summary["output_taxa"] == 2

def test_offline_status_all_keeps_genus_bucket(tmp_path):
    out = tmp_path / "all.csv"
    run_cmd([
        sys.executable, str(SCRIPT),
        "--offline-taxa-json", str(DATA / "neotoma_taxa_fixture.json"),
        "--offline-occurrences-json", str(DATA / "neotoma_occurrences_fixture.json"),
        "--organisms", "animals",
        "--status", "all",
        "--out", str(out)
    ], cwd=ROOT)
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    names = {r["species"] for r in rows}
    assert "Canis" in names

def test_offline_plants_morphotype_to_genus(tmp_path):
    out = tmp_path / "plants.csv"
    run_cmd([
        sys.executable, str(SCRIPT),
        "--offline-taxa-json", str(DATA / "neotoma_taxa_fixture.json"),
        "--offline-occurrences-json", str(DATA / "neotoma_occurrences_fixture.json"),
        "--organisms", "plants",
        "--status", "all",
        "--out", str(out)
    ], cwd=ROOT)
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    assert {r["species"] for r in rows} == {"Ambrosia"}

import csv
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PBDB = ROOT / "pbdb_to_gbif.py"
NOW = ROOT / "now_to_gbif.py"
MERGE = ROOT / "flyguide_merge_palaeo_sources.py"
DATA = ROOT / "tests" / "data"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def run_cmd(cmd, cwd=None):
    return subprocess.run(cmd, cwd=cwd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def test_pbdb_fixture_exports_flyguide_files(tmp_path):
    out = tmp_path / "pbdb_out.csv"
    run_cmd([
        sys.executable, str(PBDB),
        "--fixture", str(DATA / "pbdb_fixture.json"),
        "--region", "northern_hemisphere",
        "--period", "quaternary",
        "--organisms", "animals",
        "--out", str(out),
        "--write-flyguide-files",
        "--out-prefix", str(tmp_path / "pbdb_test"),
        "--quiet",
    ])
    rows = read_csv(out)
    names = {r["species"] for r in rows}
    assert "Bison bison" in names
    assert "Mammuthus primigenius" in names
    assert (tmp_path / "pbdb_test_species_search.txt").exists()
    assert (tmp_path / "pbdb_test_species_kingdom.tsv").exists()
    text = (tmp_path / "pbdb_test_species_search.txt").read_text()
    assert "Bison bison" in text


def test_now_merged_fixture_exports(tmp_path):
    out = tmp_path / "now_out.csv"
    run_cmd([
        sys.executable, str(NOW),
        "--input", str(DATA / "now_merged_fixture.csv"),
        "--region", "northern_hemisphere",
        "--period", "quaternary",
        "--out", str(out),
        "--write-flyguide-files",
        "--out-prefix", str(tmp_path / "now_test"),
        "--quiet",
    ])
    rows = read_csv(out)
    names = {r["species"] for r in rows}
    assert "Bison bison" in names
    assert "Mammuthus primigenius" in names
    assert "Canis" not in names  # Tanzania is outside northern hemisphere? lat -2.99, so excluded
    assert all(r["kingdom"] == "Animalia" for r in rows)


def test_merge_outputs(tmp_path):
    pbdb_out = tmp_path / "pbdb.csv"
    now_out = tmp_path / "now.csv"
    merged = tmp_path / "merged.csv"
    run_cmd([sys.executable, str(PBDB), "--fixture", str(DATA / "pbdb_fixture.json"), "--out", str(pbdb_out), "--quiet"])
    run_cmd([sys.executable, str(NOW), "--input", str(DATA / "now_merged_fixture.csv"), "--out", str(now_out), "--quiet"])
    run_cmd([
        sys.executable, str(MERGE),
        "--inputs", str(pbdb_out), str(now_out),
        "--out", str(merged),
        "--write-flyguide-files",
        "--out-prefix", str(tmp_path / "merged_test"),
    ])
    rows = read_csv(merged)
    names = {r["species"] for r in rows}
    assert "Bison bison" in names
    assert "Mammuthus primigenius" in names
    assert (tmp_path / "merged_test_species_search.txt").exists()

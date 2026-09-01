#!/usr/bin/env python3
"""Regression coverage for guidecheck.sh's real, live-hit bug (2026-08-31): NCBI EUtils
occasionally returns a non-JSON body under sustained load, and `eutils()` used to pipe
that straight into `jq` with no validation -- crashing the ENTIRE script hard on the
first bad response (`jq`'s own parse error), losing every already-processed row. Uses a
fake `curl` (a tiny script placed first on PATH) to deterministically simulate NCBI
returning garbage without touching the real network.
"""
import os
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "guidecheck.sh"

#: A fake `curl` standing in for NCBI EUtils. The taxon NAME only ever appears in the
#: very first (esearch/taxonomy) call of a species' 4-call sequence -- every follow-up
#: call (esummary, and esearch/nuccore|sra|assembly) only carries the numeric taxid, so
#: each scenario is keyed on BOTH its name (initial lookup) and its own reserved taxid
#: (every later call). `esummary.fcgi` is deliberately matched before the generic
#: `db=taxonomy` pattern below it -- the real esummary endpoint's own query string also
#: contains `db=taxonomy` (a different, unrelated `db=` param than esearch's), so matching
#: order must disambiguate by endpoint filename first.
FAKE_CURL = r"""#!/usr/bin/env bash
set -euo pipefail
url="${@: -1}"
state_dir="${FAKE_CURL_STATE_DIR:?}"

respond_for_taxid() {
  local taxid="$1" name="$2" nuccore="$3" sra="$4" assembly="$5"
  case "$url" in
    *esummary.fcgi*) echo "{\"result\":{\"${taxid}\":{\"scientificname\":\"${name}\",\"rank\":\"species\"}}}" ;;
    *db=taxonomy*)   echo "{\"esearchresult\":{\"idlist\":[\"${taxid}\"]}}" ;;
    *db=nuccore*)    echo "{\"esearchresult\":{\"count\":\"${nuccore}\"}}" ;;
    *db=sra*)        echo "{\"esearchresult\":{\"count\":\"${sra}\"}}" ;;
    *db=assembly*)   echo "{\"esearchresult\":{\"count\":\"${assembly}\"}}" ;;
  esac
}

case "$url" in
  *Permafailus*)
    # Every single call for this species returns a non-JSON body (simulating NCBI's
    # real transient rate-limit/error page) -- always, never recovers.
    echo "Rate limit exceeded, please try again later"
    ;;
  *Retryus*|*999888*)
    counter_file="$state_dir/retryus_calls"
    n=0
    [[ -f "$counter_file" ]] && n="$(cat "$counter_file")"
    n=$((n + 1))
    echo "$n" > "$counter_file"
    if (( n < 2 )); then
      echo "Rate limit exceeded, please try again later"
    else
      respond_for_taxid 999888 "Retryus transientbad" 5 3 0
    fi
    ;;
  *Goodus*|*111222*)
    respond_for_taxid 111222 "Goodus normal" 10 0 0
    ;;
  *)
    echo '{"esearchresult":{"idlist":[]}}'
    ;;
esac
"""


def _install_fake_curl(bin_dir: Path) -> None:
    bin_dir.mkdir(parents=True, exist_ok=True)
    fake_curl = bin_dir / "curl"
    fake_curl.write_text(FAKE_CURL, encoding="utf-8")
    fake_curl.chmod(fake_curl.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def _run_guidecheck(tmp_path: Path, taxa_names: list[str]) -> subprocess.CompletedProcess:
    bin_dir = tmp_path / "bin"
    _install_fake_curl(bin_dir)
    state_dir = tmp_path / "fake_curl_state"
    state_dir.mkdir()
    taxa_file = tmp_path / "taxa.txt"
    taxa_file.write_text("\n".join(taxa_names) + "\n", encoding="utf-8")

    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["FAKE_CURL_STATE_DIR"] = str(state_dir)
    # No real sleeping in tests -- the retry logic itself is what's under test, not timing.
    env["GUIDECHECK_EUTILS_RETRY_DELAY"] = "0"

    return subprocess.run(
        ["bash", str(SCRIPT), "-i", str(taxa_file)],
        cwd=tmp_path, env=env, text=True, capture_output=True,
    )


def test_permanently_bad_ncbi_response_degrades_to_no_taxid_instead_of_crashing(tmp_path):
    result = _run_guidecheck(tmp_path, ["Permafailus alwaysbad"])
    assert result.returncode == 0, result.stderr + result.stdout
    rows = result.stdout.strip().splitlines()
    assert rows[0].startswith("query_name\t")
    assert rows[1] == "Permafailus alwaysbad\tNA\tNA\tNA\t0\t0\t0\tNO_TAXID"
    # The real, resilience-specific behavior: warnings on stderr, not a jq parse-error crash.
    assert "jq: error" not in result.stderr
    assert "non-JSON/empty NCBI response" in result.stderr


def test_transient_bad_response_recovers_via_retry_with_real_data(tmp_path):
    result = _run_guidecheck(tmp_path, ["Retryus transientbad"])
    assert result.returncode == 0, result.stderr + result.stdout
    rows = result.stdout.strip().splitlines()
    assert rows[1] == "Retryus transientbad\t999888\tRetryus transientbad\tspecies\t5\t3\t0\tHAS_NUCCORE_AND_SRA"


def test_binomial_matched_name_survives_intact_not_split_on_its_own_space(tmp_path):
    """Second bug in the same file, found while writing the retry-recovery test above:
    `matched_name`/`rank` used to be parsed via `awk '{print $1, $2}'`, which re-joins the
    real TSV fields with a plain space, destroying the tab boundary -- corrupting any
    `matched_name` with its own internal space (i.e. nearly every real binomial species
    name, not an edge case) into a truncated genus-only value plus a garbled rank."""
    result = _run_guidecheck(tmp_path, ["Goodus normal"])
    assert result.returncode == 0, result.stderr + result.stdout
    rows = result.stdout.strip().splitlines()
    assert rows[1] == "Goodus normal\t111222\tGoodus normal\tspecies\t10\t0\t0\tNUCCORE_ONLY"


def test_one_bad_species_does_not_stop_the_rest_of_the_batch_from_processing(tmp_path):
    """The real bug's most damaging effect: ONE bad response anywhere in a multi-thousand
    species run used to kill the entire batch, losing every already-computed row along
    with it. A species after the permanently-failing one must still be processed."""
    result = _run_guidecheck(tmp_path, ["Goodus normal", "Permafailus alwaysbad", "Goodus normal"])
    assert result.returncode == 0, result.stderr + result.stdout
    rows = result.stdout.strip().splitlines()
    assert len(rows) == 4  # header + 3 species, none dropped
    assert rows[1].startswith("Goodus normal\t111222\t")
    assert rows[2] == "Permafailus alwaysbad\tNA\tNA\tNA\t0\t0\t0\tNO_TAXID"
    assert rows[3].startswith("Goodus normal\t111222\t")


if __name__ == "__main__":
    sys.exit(subprocess.run([sys.executable, "-m", "pytest", __file__, "-v"]).returncode)

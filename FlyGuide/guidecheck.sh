#!/usr/bin/env bash
# guidecheck.sh
#
# TackleBox: FlyGuide - GuideCheck
#
# Quick NCBI "coverage check" for a list of taxa (species or genera).
#
# Examples:
#   bash guidecheck.sh -i taxa.txt > ncbi_counts.tsv
#   bash guidecheck.sh -i taxa.txt -o ncbi_counts.tsv
#   bash guidecheck.sh -i taxa.txt --name-field ALL -o ncbi_counts.tsv
#   bash guidecheck.sh -i taxa.txt --api-key "YOUR_KEY" -o ncbi_counts.tsv
#
# Typical use inside FlyGuide:
#   - Input is OUTPREFIX_species_search.txt
#   - Output is OUTPREFIX_ncbi_guidecheck.tsv
#
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  guidecheck.sh -i taxa.txt [-o results.tsv] [--name-field SCIN|ALL] [--api-key KEY]

Inputs:
  taxa.txt: one taxon name per line.

Options:
  -i, --input       Input taxa file (required)
  -o, --output      Output TSV file (default: stdout)
  --name-field      How to search taxonomy:
                     SCIN = scientific name field (default; stricter, fewer false hits)
                     ALL  = general text search (more hits, more ambiguity)
  --api-key         NCBI API key (optional; higher rate limits)
  -h, --help        Show this help

Output columns:
  query_name        Name from the input file
  taxid             Matched NCBI taxonomy ID (or NA)
  matched_name      Canonical scientific name at taxid
  rank              Taxonomic rank (species, genus, etc.)
  nuccore           Count of nuccore records for txid[Organism:exp]
  sra               Count of SRA records for txid[Organism:exp]
  assembly          Count of assembly records for txid[Organism:exp]
  status            Simple label summarizing availability:
                      NO_TAXID
                      ASSEMBLY
                      HAS_NUCCORE_AND_SRA
                      SRA_ONLY
                      NUCCORE_ONLY
                      NONE
EOF
}

INPUT=""
OUTPUT=""
NAME_FIELD="SCIN"
API_KEY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="${2:-}"; shift 2 ;;
    -o|--output) OUTPUT="${2:-}"; shift 2 ;;
    --name-field) NAME_FIELD="${2:-}"; shift 2 ;;
    --api-key) API_KEY="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -z "$INPUT" ]] && { echo "ERROR: --input is required" >&2; usage; exit 2; }
[[ -r "$INPUT" ]] || { echo "ERROR: Can't read input file: $INPUT" >&2; exit 2; }

command -v curl >/dev/null || { echo "ERROR: curl not found" >&2; exit 2; }
command -v jq   >/dev/null || { echo "ERROR: jq not found" >&2; exit 2; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 2; }

urlenc() {
  python3 - <<'PY' "$1"
import sys, urllib.parse
print(urllib.parse.quote(sys.argv[1]))
PY
}

#: NCBI EUtils occasionally returns a non-JSON body under sustained load (a transient
#: rate-limit/error page rather than the documented JSON response) -- previously this
#: crashed the ENTIRE script hard via `jq`'s own parse error ("Invalid numeric literal at
#: line 1, column 10"), losing every already-processed row along with it, since there was
#: no retry/validation of any kind. Real failure, hit twice in a row on a live ~10,000-
#: species run (species #3, then again at #37 on retry), 2026-08-31. `eutils()` now
#: validates the response is real JSON (`jq empty`) before returning it, retrying a few
#: times on failure, and degrading to a safe empty JSON object (`{}`) after exhausting
#: retries rather than ever propagating bad data downstream -- every caller's own `jq`
#: filter already treats a missing/null field as "no data" (NO_TAXID / NA / a zero count),
#: so `{}` degrades that one query to "nothing found" instead of crashing the whole run.
#: Both knobs are overridable (used by this script's own test suite to avoid real sleeps).
EUTILS_MAX_ATTEMPTS="${GUIDECHECK_EUTILS_MAX_ATTEMPTS:-3}"
EUTILS_RETRY_DELAY="${GUIDECHECK_EUTILS_RETRY_DELAY:-1}"

eutils() {
  local url="$1"
  if [[ -n "$API_KEY" ]]; then
    if [[ "$url" == *"?"* ]]; then
      url="${url}&api_key=$(urlenc "$API_KEY")"
    else
      url="${url}?api_key=$(urlenc "$API_KEY")"
    fi
  fi
  local attempt body
  for (( attempt=1; attempt<=EUTILS_MAX_ATTEMPTS; attempt++ )); do
    body="$(curl -s "$url" || true)"
    if [[ -n "$body" ]] && jq empty <<<"$body" >/dev/null 2>&1; then
      printf '%s' "$body"
      return 0
    fi
    echo "WARNING: guidecheck.sh: non-JSON/empty NCBI response (attempt ${attempt}/${EUTILS_MAX_ATTEMPTS}): ${url}" >&2
    if (( attempt < EUTILS_MAX_ATTEMPTS )); then
      sleep "$EUTILS_RETRY_DELAY"
    fi
  done
  echo "WARNING: guidecheck.sh: giving up after ${EUTILS_MAX_ATTEMPTS} attempts, treating as empty result: ${url}" >&2
  printf '{}'
}

tax_term() {
  local name="$1"
  case "$NAME_FIELD" in
    SCIN) printf "%s[SCIN]" "$name" ;;
    ALL)  printf "%s" "$name" ;;
    *) echo "ERROR: --name-field must be SCIN or ALL (got: $NAME_FIELD)" >&2; exit 2 ;;
  esac
}

get_taxid() {
  local name="$1"
  local term; term="$(urlenc "$(tax_term "$name")")"
  eutils "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmode=json&term=${term}" \
    | jq -r '.esearchresult.idlist[0] // empty'
}

get_tax_summary() {
  local taxid="$1"
  eutils "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&retmode=json&id=${taxid}" \
    | jq -r --arg id "$taxid" '
      .result[$id] as $r
      | [($r.scientificname // "NA"), ($r.rank // "NA")] | @tsv
    '
}

count_db() {
  local db="$1" taxid="$2"
  local term; term="$(urlenc "txid${taxid}[Organism:exp]")"
  eutils "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=${db}&retmode=json&term=${term}" \
    | jq -r '.esearchresult.count // "0"'
}

status_label() {
  local nuccore="$1" sra="$2" assembly="$3" taxid="$4"
  if [[ "$taxid" == "NA" ]]; then echo "NO_TAXID"; return; fi
  local n s a; n=$((nuccore)); s=$((sra)); a=$((assembly))
  if (( a > 0 )); then echo "ASSEMBLY"
  elif (( n > 0 && s > 0 )); then echo "HAS_NUCCORE_AND_SRA"
  elif (( s > 0 && n == 0 )); then echo "SRA_ONLY"
  elif (( n > 0 && s == 0 )); then echo "NUCCORE_ONLY"
  else echo "NONE"
  fi
}

# ---- progress helpers (stderr only) ----
TOTAL="$(awk '!/^[[:space:]]*($|#)/{c++} END{print c+0}' "$INPUT")"
SHOW_PROGRESS=0
[[ -t 2 ]] && SHOW_PROGRESS=1

progress() {
  # progress <i> <name>
  local i="$1" nm="$2"
  (( SHOW_PROGRESS == 1 )) || return 0
  # truncate the displayed name a bit
  local disp="$nm"
  if (( ${#disp} > 60 )); then disp="${disp:0:57}..."; fi
  printf '\r[%d/%d] %s' "$i" "$TOTAL" "$disp" >&2
}

progress_done() {
  (( SHOW_PROGRESS == 1 )) || return 0
  printf '\n' >&2
}

# -------- run --------
{
  echo -e "query_name\ttaxid\tmatched_name\trank\tnuccore\tsra\tassembly\tstatus"

  i=0
  while IFS= read -r raw || [[ -n "$raw" ]]; do
    name="$(echo "$raw" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    [[ -z "$name" ]] && continue
    [[ "$name" == \#* ]] && continue

    i=$((i+1))
    progress "$i" "$name"

    taxid="$(get_taxid "$name")"
    if [[ -z "$taxid" ]]; then
      echo -e "${name}\tNA\tNA\tNA\t0\t0\t0\tNO_TAXID"
      continue
    fi

    # Second real bug, found 2026-08-31 while writing a regression test for the retry fix
    # above: piping the TSV line through `awk '{print $1, $2}'` re-joins its two fields
    # with awk's own OFS (a single space) BEFORE `read` ever sees them -- destroying the
    # original tab boundary. Any real binomial `matched_name` (e.g. "Mammuthus
    # primigenius" -- the overwhelming majority of real successful matches, not an edge
    # case) then has its own internal space misread as the matched_name/rank boundary:
    # `read -r matched_name rank` silently produced `matched_name="Mammuthus"`,
    # `rank="primigenius species"`. Reading the TSV line directly with `IFS=$'\t'`
    # preserves the real tab-delimited field boundary regardless of internal spaces.
    IFS=$'\t' read -r matched_name rank < <(get_tax_summary "$taxid")
    nuccore="$(count_db nuccore "$taxid")"
    sra="$(count_db sra "$taxid")"
    assembly="$(count_db assembly "$taxid")"
    status="$(status_label "$nuccore" "$sra" "$assembly" "$taxid")"

    echo -e "${name}\t${taxid}\t${matched_name}\t${rank}\t${nuccore}\t${sra}\t${assembly}\t${status}"
  done < "$INPUT"

  progress_done
} | {
  if [[ -n "$OUTPUT" ]]; then
    cat > "$OUTPUT"
  else
    cat
  fi
}

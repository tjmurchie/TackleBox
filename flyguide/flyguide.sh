#!/usr/bin/env bash
set -euo pipefail

# TackleBox: FlyGuide
#
# Build region-specific NCBI nucleotide reference panels from a GBIF download:
#   1) Parse GBIF CSV/TSV → species/genus search list + species↔kingdom map.
#   2) Download curated NCBI nucleotide records for those taxa.
#   3) Split combined FASTA by kingdom + genome/marker type.
#   4) Optionally run "GuideCheck" (guidecheck.sh) to summarize NCBI content
#      for each query name.
#
# Usage:
#   flyguide.sh [--no-guidecheck] GBIF_download.csv OUTPREFIX EMAIL [NCBI_API_KEY]
#
# Examples:
#   # Default: run full pipeline + GuideCheck
#   ./flyguide.sh RegionX_GBIF.csv RegionX_refs you@example.org MY_NCBI_KEY
#
#   # Skip GuideCheck (faster, fewer API calls)
#   ./flyguide.sh --no-guidecheck RegionX_GBIF.csv RegionX_refs you@example.org MY_NCBI_KEY
#
# Requirements (on PATH / in env):
#   - bash
#   - python or python3 (for gbif_prep_from_csv.py)
#   - Perl 5.x with Bio::DB::EUtilities (for NCBI-NT_Downloader.pl)
#   - curl, jq, python3 (for guidecheck.sh)
#
# All helper scripts are expected to live next to this wrapper.

usage() {
  cat <<'EOF'
Usage:
  flyguide.sh [--no-guidecheck] GBIF_download.csv OUTPREFIX EMAIL [NCBI_API_KEY]

Positional arguments:
  GBIF_download.csv   GBIF occurrence or checklist CSV/TSV
                      (must have columns: species, genus, kingdom)
  OUTPREFIX           Prefix for generated files (in current working directory)
  EMAIL               Contact email for NCBI E-utilities (required)
  NCBI_API_KEY        Optional NCBI API key (recommended)

Options:
  --no-guidecheck     Do NOT run GuideCheck (guidecheck.sh) on the species list
  -h, --help          Show this help and exit

Typical workflow:
  1) Export a GBIF CSV for your region/taxa.
  2) Run flyguide.sh as above.
  3) Use the split FASTAs and the *_ncbi_guidecheck.tsv file as inputs to
     downstream mapping / classifier tools.

EOF
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

RUN_GUIDECHECK=1

GBIF_CSV=""
OUTPREFIX=""
EMAIL=""
API_KEY=""

# Parse flags + positional args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-guidecheck)
      RUN_GUIDECHECK=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 2
      ;;
    *)
      # positional
      if [[ -z "$GBIF_CSV" ]]; then
        GBIF_CSV="$1"
      elif [[ -z "$OUTPREFIX" ]]; then
        OUTPREFIX="$1"
      elif [[ -z "$EMAIL" ]]; then
        EMAIL="$1"
      elif [[ -z "$API_KEY" ]]; then
        API_KEY="$1"
      else
        echo "ERROR: Too many positional arguments." >&2
        usage
        exit 2
      fi
      shift
      ;;
  esac
done

if [[ -z "$GBIF_CSV" || -z "$OUTPREFIX" || -z "$EMAIL" ]]; then
  echo "ERROR: GBIF_download.csv, OUTPREFIX, and EMAIL are required." >&2
  usage
  exit 2
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GBIF_PREP="${SCRIPT_DIR}/gbif_prep_from_csv.py"
NCBI_DL="${SCRIPT_DIR}/NCBI-NT_Downloader.pl"
SPLIT_FASTA="${SCRIPT_DIR}/split_fasta_by_kingdom_organelle_simple.pl"
GUIDECHECK="${SCRIPT_DIR}/guidecheck.sh"

# Basic checks
[[ -f "$GBIF_CSV" ]] || { echo "ERROR: GBIF CSV not found: $GBIF_CSV" >&2; exit 1; }
[[ -f "$GBIF_PREP" ]] || { echo "ERROR: gbif_prep_from_csv.py not found in $SCRIPT_DIR" >&2; exit 1; }
[[ -f "$NCBI_DL" ]] || { echo "ERROR: NCBI-NT_Downloader.pl not found in $SCRIPT_DIR" >&2; exit 1; }
[[ -f "$SPLIT_FASTA" ]] || { echo "ERROR: split_fasta_by_kingdom_organelle_simple.pl not found in $SCRIPT_DIR" >&2; exit 1; }

if [[ "$RUN_GUIDECHECK" -eq 1 ]]; then
  if [[ ! -f "$GUIDECHECK" ]]; then
    echo "WARNING: guidecheck.sh not found in $SCRIPT_DIR; GuideCheck will be skipped." >&2
    RUN_GUIDECHECK=0
  fi
fi

echo "=== TackleBox: FlyGuide ==="
echo "GBIF CSV      : $GBIF_CSV"
echo "OUTPREFIX     : $OUTPREFIX"
echo "EMAIL         : $EMAIL"
echo "API_KEY       : ${API_KEY:-<none>}"
echo "GuideCheck    : $([[ $RUN_GUIDECHECK -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo

###############################################################################
# Step 1: GBIF prep → species search list + species↔kingdom map
###############################################################################

echo "=== Step 1: Preparing GBIF species/genus search list and species↔kingdom map ==="
"$PYTHON_BIN" "$GBIF_PREP" "$GBIF_CSV" "$OUTPREFIX"

SPECIES_LIST="${OUTPREFIX}_species_search.txt"
SPECIES_KINGDOM="${OUTPREFIX}_species_kingdom.tsv"

if [[ ! -s "$SPECIES_LIST" ]]; then
  echo "ERROR: Species search list missing or empty: $SPECIES_LIST" >&2
  exit 1
fi

if [[ ! -s "$SPECIES_KINGDOM" ]]; then
  echo "ERROR: Species↔kingdom TSV missing or empty: $SPECIES_KINGDOM" >&2
  exit 1
fi

###############################################################################
# Optional Step 1.5: GuideCheck → NCBI coverage summary
###############################################################################

if [[ "$RUN_GUIDECHECK" -eq 1 ]]; then
  echo "=== Step 1.5: Running GuideCheck (NCBI coverage summary) ==="
  GUIDE_OUT="${OUTPREFIX}_ncbi_guidecheck.tsv"

  # Build args for guidecheck.sh
  GUIDECHECK_ARGS=( -i "$SPECIES_LIST" -o "$GUIDE_OUT" )
  if [[ -n "$API_KEY" ]]; then
    GUIDECHECK_ARGS+=( --api-key "$API_KEY" )
  fi

  bash "$GUIDECHECK" "${GUIDECHECK_ARGS[@]}"

  if [[ -s "$GUIDE_OUT" ]]; then
    echo "GuideCheck summary written to: $GUIDE_OUT"
  else
    echo "WARNING: GuideCheck did not produce a non-empty TSV: $GUIDE_OUT" >&2
  fi
fi

###############################################################################
# Step 2: NCBI download
###############################################################################

echo "=== Step 2: Running NCBI-NT_Downloader.pl ==="
if [[ -n "$API_KEY" ]]; then
  perl "$NCBI_DL" "$SPECIES_LIST" "$OUTPREFIX" "$EMAIL" "$API_KEY"
else
  perl "$NCBI_DL" "$SPECIES_LIST" "$OUTPREFIX" "$EMAIL"
fi

FASTA_OUT="${OUTPREFIX}.fasta"
if [[ ! -s "$FASTA_OUT" ]]; then
  echo "ERROR: Downloader FASTA output not found or empty: $FASTA_OUT" >&2
  exit 1
fi

###############################################################################
# Step 3: Split FASTA by kingdom + genome/marker type
###############################################################################

echo "=== Step 3: Splitting FASTA by kingdom + genome/marker type ==="
perl "$SPLIT_FASTA" "$FASTA_OUT" "$OUTPREFIX" "$SPECIES_KINGDOM"

echo "=== FlyGuide pipeline complete ==="
echo "Key outputs (in current working directory):"
echo "  Species search list           : $SPECIES_LIST"
echo "  Species↔kingdom map (TSV)     : $SPECIES_KINGDOM"
echo "  Combined FASTA from NCBI      : $FASTA_OUT"
if [[ "$RUN_GUIDECHECK" -eq 1 ]]; then
  echo "  NCBI coverage summary (TSV)   : ${OUTPREFIX}_ncbi_guidecheck.tsv"
fi
echo "  Split FASTAs (kingdom/type)   : ${OUTPREFIX}.<Kingdom>-<Type>.fasta"

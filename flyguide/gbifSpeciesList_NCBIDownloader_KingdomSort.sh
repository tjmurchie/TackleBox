#!/usr/bin/env bash
set -euo pipefail

# gbifSpeciesList_NCBIDownloader_KingdomSort
#
# Wrapper to:
#   1) Take a GBIF CSV download.
#   2) Build an NCBI species/genus search list and a species↔kingdom map.
#   3) Run NCBI-NT_Downloader.pl on that search list.
#   4) Split the resulting FASTA by kingdom + genome/marker type.
#
# Usage:
#   ./gbifSpeciesList_NCBIDownloader_KingdomSort.sh GBIF_download.csv OUTPREFIX EMAIL [NCBI_API_KEY]
#
# You can run this from the directory containing your GBIF CSV.
# The script will automatically locate the helper scripts based on
# its own install location.
#
# Requirements:
#   - gbif_prep_from_csv.py
#   - NCBI-NT_Downloader.pl
#   - split_fasta_by_kingdom_organelle_simple.pl
#   All three should live in the same directory as this wrapper.

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 GBIF_download.csv OUTPREFIX EMAIL [NCBI_API_KEY]" >&2
  exit 1
fi

GBIF_CSV="$1"
OUTPREFIX="$2"
EMAIL="$3"
API_KEY="${4:-}"

PYTHON_BIN="${PYTHON_BIN:-python}"

# Directory where this script resides (tool installation directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "$GBIF_CSV" ]]; then
  echo "ERROR: GBIF CSV not found: $GBIF_CSV" >&2
  exit 1
fi

if [[ ! -f "$SCRIPT_DIR/gbif_prep_from_csv.py" ]]; then
  echo "ERROR: gbif_prep_from_csv.py not found next to this wrapper (in $SCRIPT_DIR)." >&2
  exit 1
fi

if [[ ! -f "$SCRIPT_DIR/NCBI-NT_Downloader.pl" ]]; then
  echo "ERROR: NCBI-NT_Downloader.pl not found next to this wrapper (in $SCRIPT_DIR)." >&2
  exit 1
fi

if [[ ! -f "$SCRIPT_DIR/split_fasta_by_kingdom_organelle_simple.pl" ]]; then
  echo "ERROR: split_fasta_by_kingdom_organelle_simple.pl not found next to this wrapper (in $SCRIPT_DIR)." >&2
  exit 1
fi

echo "=== Step 1: Preparing GBIF species search list and species↔kingdom map ==="
$PYTHON_BIN "$SCRIPT_DIR/gbif_prep_from_csv.py" "$GBIF_CSV" "$OUTPREFIX"

SPECIES_LIST="${OUTPREFIX}_species_search.txt"
SPECIES_KINGDOM="${OUTPREFIX}_species_kingdom.tsv"

if [[ ! -s "$SPECIES_LIST" ]]; then
  echo "ERROR: Species search list is empty or missing: $SPECIES_LIST" >&2
  exit 1
fi

if [[ ! -s "$SPECIES_KINGDOM" ]]; then
  echo "ERROR: Species↔kingdom TSV is empty or missing: $SPECIES_KINGDOM" >&2
  exit 1
fi

echo "=== Step 2: Running NCBI-NT_Downloader.pl ==="
if [[ -n "$API_KEY" ]]; then
  perl "$SCRIPT_DIR/NCBI-NT_Downloader.pl" "$SPECIES_LIST" "$OUTPREFIX" "$EMAIL" "$API_KEY"
else
  perl "$SCRIPT_DIR/NCBI-NT_Downloader.pl" "$SPECIES_LIST" "$OUTPREFIX" "$EMAIL"
fi

FASTA_OUT="${OUTPREFIX}.fasta"
if [[ ! -s "$FASTA_OUT" ]]; then
  echo "ERROR: Downloader FASTA output not found or empty: $FASTA_OUT" >&2
  exit 1
fi

echo "=== Step 3: Splitting FASTA by kingdom + genome/marker type ==="
perl "$SCRIPT_DIR/split_fasta_by_kingdom_organelle_simple.pl" \
  "$FASTA_OUT" \
  "$OUTPREFIX" \
  "$SPECIES_KINGDOM"

echo "=== Pipeline complete ==="
echo "Key outputs (in current working directory):"
echo "  Species search list          : $SPECIES_LIST"
echo "  Species↔kingdom map (TSV)    : $SPECIES_KINGDOM"
echo "  Combined FASTA from NCBI     : $FASTA_OUT"
echo "  Split FASTAs (kingdom/type)  : ${OUTPREFIX}.<Kingdom>-<Type>.fasta"

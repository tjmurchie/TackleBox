#!/usr/bin/env bash
set -euo pipefail

FLYGUIDE_VERSION="1.3.0"

# TackleBox: FlyGuide
#
# Build region-specific NCBI nucleotide reference panels from a GBIF download:
#   1) Parse GBIF CSV/TSV -> species/genus search list + species<->kingdom map.
#   2) (Optional) run GuideCheck to summarize NCBI coverage per taxon.
#   3) Download curated NCBI nucleotide records for those taxa.
#   4) Split combined FASTA by kingdom (+ optionally phylum) + region class.
#
# Usage:
#   flyguide.sh [OPTIONS] GBIF_download.csv OUTPREFIX EMAIL [NCBI_API_KEY]
#
# Examples:
#   ./flyguide.sh region.csv Region_refs you@example.org MY_KEY
#   ./flyguide.sh --split-phylum region.csv Region_refs you@example.org
#   ./flyguide.sh --enable-regions COI --max-per-taxon 500 region.csv out you@x.org
#   ./flyguide.sh --no-tui region.csv out you@x.org          # plain text output
#
# Requirements:
#   - bash, python3 (or $PYTHON_BIN), perl + Bio::DB::EUtilities
#   - curl, jq, python3 for guidecheck.sh (when --guidecheck is used)

usage() {
  cat <<'EOF'
Usage:
  flyguide.sh [OPTIONS] GBIF_download.csv OUTPREFIX EMAIL [NCBI_API_KEY]

Positional arguments:
  GBIF_download.csv   GBIF occurrence or checklist CSV/TSV
                      (must have columns: species, genus, kingdom;
                       phylum is used automatically when present)
  OUTPREFIX           Prefix for all generated files
  EMAIL               Contact email for NCBI E-utilities (required)
  NCBI_API_KEY        Optional NCBI API key (strongly recommended)

Pipeline options:
  --guidecheck        Enable GuideCheck (NCBI coverage summary TSV).
                      Adds extra NCBI calls but produces per-taxon diagnostics.
  --no-guidecheck     Explicitly disable GuideCheck (default).
  --split-phylum      Split output FASTAs by phylum in addition to kingdom
                      and region type (e.g. Plant.Tracheophyta-Plastid.fasta).
                      Requires phylum column in GBIF file.

Download options (passed to NCBI-NT_Downloader.pl):
  --max-per-taxon N   Max NCBI records per taxon (default: 1000)
  --min-slen N        Min sequence length in bp (default: 50)
  --max-slen N        Max sequence length in bp (default: 400000)
  --no-organelle      Fetch nuclear markers only; skip mito/plastid records
  --enable-regions R  Comma-separated region_id(s) to force-enable
                      (e.g. --enable-regions COI to add COI to every query)
  --disable-regions R Comma-separated region_id(s) to force-disable
                      (e.g. --disable-regions NUC_H3,NUCRDNA_ITS)
  --no-tui            Disable the interactive TUI dashboard; plain line output
  --ncbi-name-mode M  Name normalization mode passed to NCBI downloader:
                      species (legacy default), trinomial, or as-is

General:
  -h, --help          Show this help and exit
  -V, --version       Print version and exit

Query customization:
  The NCBI search query is built from two parts for each taxon:
    1) Organelle block  : mitochondrial / chloroplast / plastid records
       (disable with --no-organelle)
    2) Marker block     : built from regions_config.tsv (ncbi_title_clause column)
       or built-in defaults (18S/28S/ITS/5.8S/H3) if no config clauses found
  To customize which markers are searched, edit regions_config.tsv:
    - Set enabled_default=0 to disable a marker by default
    - Enable per-run with --enable-regions <region_id>
    - Disable per-run with --disable-regions <region_id>
  Sequence length range: --min-slen / --max-slen (default 50–400,000 bp)
EOF
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

# ── Defaults ─────────────────────────────────────────────────────────────
RUN_GUIDECHECK=0
SPLIT_PHYLUM=0

GBIF_CSV=""
OUTPREFIX=""
EMAIL=""
API_KEY=""

# Downloader pass-through options (accumulated as array)
DOWNLOADER_OPTS=()

# ── Parse flags + positional args ────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --guidecheck)
      RUN_GUIDECHECK=1
      shift ;;
    --no-guidecheck)
      RUN_GUIDECHECK=0
      shift ;;
    --split-phylum)
      SPLIT_PHYLUM=1
      shift ;;
    --max-per-taxon)
      DOWNLOADER_OPTS+=("--max-per-taxon" "$2")
      shift 2 ;;
    --min-slen)
      DOWNLOADER_OPTS+=("--min-slen" "$2")
      shift 2 ;;
    --max-slen)
      DOWNLOADER_OPTS+=("--max-slen" "$2")
      shift 2 ;;
    --no-organelle)
      DOWNLOADER_OPTS+=("--no-organelle")
      shift ;;
    --enable-regions)
      DOWNLOADER_OPTS+=("--enable-regions" "$2")
      shift 2 ;;
    --disable-regions)
      DOWNLOADER_OPTS+=("--disable-regions" "$2")
      shift 2 ;;
    --no-tui)
      DOWNLOADER_OPTS+=("--no-tui")
      shift ;;
    --ncbi-name-mode)
      DOWNLOADER_OPTS+=("--name-mode" "$2")
      shift 2 ;;
    -h|--help)
      usage
      exit 0 ;;
    -V|--version)
      echo "FlyGuide version ${FLYGUIDE_VERSION}"
      exit 0 ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 2 ;;
    *)
      if   [[ -z "$GBIF_CSV"   ]]; then GBIF_CSV="$1"
      elif [[ -z "$OUTPREFIX"  ]]; then OUTPREFIX="$1"
      elif [[ -z "$EMAIL"      ]]; then EMAIL="$1"
      elif [[ -z "$API_KEY"    ]]; then API_KEY="$1"
      else
        echo "ERROR: Too many positional arguments." >&2
        usage
        exit 2
      fi
      shift ;;
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
REGIONS_CONF="${SCRIPT_DIR}/regions_config.tsv"

[[ -f "$GBIF_CSV"    ]] || { echo "ERROR: GBIF CSV not found: $GBIF_CSV" >&2; exit 1; }
[[ -f "$GBIF_PREP"   ]] || { echo "ERROR: gbif_prep_from_csv.py not found in $SCRIPT_DIR" >&2; exit 1; }
[[ -f "$NCBI_DL"     ]] || { echo "ERROR: NCBI-NT_Downloader.pl not found in $SCRIPT_DIR" >&2; exit 1; }
[[ -f "$SPLIT_FASTA" ]] || { echo "ERROR: split_fasta_by_kingdom_organelle_simple.pl not found in $SCRIPT_DIR" >&2; exit 1; }

if [[ "$RUN_GUIDECHECK" -eq 1 && ! -f "$GUIDECHECK" ]]; then
  echo "WARNING: guidecheck.sh not found; GuideCheck will be skipped." >&2
  RUN_GUIDECHECK=0
fi

###############################################################################
# Step 1: GBIF prep
###############################################################################
echo "=== Step 1: Preparing GBIF species/genus search list and species<->kingdom map ==="
"$PYTHON_BIN" "$GBIF_PREP" "$GBIF_CSV" "$OUTPREFIX"

SPECIES_LIST="${OUTPREFIX}_species_search.txt"
SPECIES_KINGDOM="${OUTPREFIX}_species_kingdom.tsv"

[[ -s "$SPECIES_LIST"    ]] || { echo "ERROR: Species search list missing or empty: $SPECIES_LIST" >&2; exit 1; }
[[ -s "$SPECIES_KINGDOM" ]] || { echo "ERROR: Species<->kingdom TSV missing or empty: $SPECIES_KINGDOM" >&2; exit 1; }

###############################################################################
# Step 1.5 (optional): GuideCheck
###############################################################################
if [[ "$RUN_GUIDECHECK" -eq 1 ]]; then
  echo "=== Step 1.5: Running GuideCheck (NCBI coverage summary) ==="
  GUIDE_OUT="${OUTPREFIX}_ncbi_guidecheck.tsv"
  GUIDECHECK_ARGS=(-i "$SPECIES_LIST" -o "$GUIDE_OUT")
  [[ -n "$API_KEY" ]] && GUIDECHECK_ARGS+=(--api-key "$API_KEY")
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
DOWNLOADER_ARGS=("${DOWNLOADER_OPTS[@]}" "$SPECIES_LIST" "$OUTPREFIX" "$EMAIL")
[[ -n "$API_KEY" ]] && DOWNLOADER_ARGS+=("$API_KEY")
perl "$NCBI_DL" "${DOWNLOADER_ARGS[@]}"

FASTA_OUT="${OUTPREFIX}.fasta"
[[ -s "$FASTA_OUT" ]] || { echo "ERROR: Downloader FASTA output not found or empty: $FASTA_OUT" >&2; exit 1; }

###############################################################################
# Step 3: Split FASTA
###############################################################################
echo "=== Step 3: Splitting FASTA by kingdom$([ "$SPLIT_PHYLUM" -eq 1 ] && echo " + phylum") + region class ==="

SPLITTER_ARGS=()
[[ "$SPLIT_PHYLUM" -eq 1 ]] && SPLITTER_ARGS+=("--split-phylum")
SPLITTER_ARGS+=("$FASTA_OUT" "$OUTPREFIX" "$SPECIES_KINGDOM")
[[ -f "$REGIONS_CONF" ]] && SPLITTER_ARGS+=("$REGIONS_CONF")

perl "$SPLIT_FASTA" "${SPLITTER_ARGS[@]}"

###############################################################################
# Summary
###############################################################################
echo
echo "=== FlyGuide pipeline complete ==="
echo "  Species search list           : $SPECIES_LIST"
echo "  Species<->kingdom map (TSV)   : $SPECIES_KINGDOM"
echo "  Combined FASTA from NCBI      : $FASTA_OUT"
[[ "$RUN_GUIDECHECK" -eq 1 ]] && echo "  NCBI coverage summary (TSV)   : ${OUTPREFIX}_ncbi_guidecheck.tsv"
if [[ "$SPLIT_PHYLUM" -eq 1 ]]; then
  echo "  Split FASTAs (kingdom/phylum) : ${OUTPREFIX}.<Kingdom>.<Phylum>-<Class>.fasta"
else
  echo "  Split FASTAs (kingdom/class)  : ${OUTPREFIX}.<Kingdom>-<Class>.fasta"
fi

#!/usr/bin/env bash
set -euo pipefail

# Convenience wrapper: Neotoma -> GBIF-like CSV -> FlyGuide.
# For maximum control, run neotoma_extinct_to_gbif.py directly, inspect the CSV,
# then run flyguide.sh manually.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"

usage() {
  cat <<'EOF'
Usage:
  flyguide_neotoma.sh [NEOTOMA OPTIONS] -- OUTPREFIX EMAIL [NCBI_API_KEY] [FLYGUIDE OPTIONS]

Examples:
  ./flyguide_neotoma.sh \
    --region northern_hemisphere --period quaternary --organisms animals --status extinct \
    -- My_NH_Quaternary_refs you@example.org NCBI_KEY --guidecheck --max-per-taxon 500

  ./flyguide_neotoma.sh \
    --region tanzania --period pleistocene --organisms animals --ncbi-name-mode binomial \
    -- Tanzania_Pleistocene_refs you@example.org

Notes:
  - Options before '--' are passed to neotoma_extinct_to_gbif.py.
  - OUTPREFIX/EMAIL/API_KEY and options after '--' are passed to flyguide.sh.
  - The wrapper writes OUTPREFIX_neotoma_gbif.csv and then runs FlyGuide.
  - For very broad regions, first run neotoma_extinct_to_gbif.py alone and inspect
    the CSV/rejected log before launching NCBI downloads.
EOF
}

if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

NEOTOMA_ARGS=()
while [[ $# -gt 0 ]]; do
  if [[ "$1" == "--" ]]; then
    shift
    break
  fi
  NEOTOMA_ARGS+=("$1")
  shift
done

if [[ $# -lt 2 ]]; then
  echo "ERROR: expected OUTPREFIX EMAIL [NCBI_API_KEY] after --" >&2
  usage >&2
  exit 2
fi

OUTPREFIX="$1"; shift
EMAIL="$1"; shift
API_KEY=""
if [[ $# -gt 0 && "${1:-}" != --* ]]; then
  API_KEY="$1"
  shift
fi
FLYGUIDE_OPTS=("$@")

CSV_OUT="${OUTPREFIX}_neotoma_gbif.csv"

# If user supplied a name mode for Neotoma export, keep FlyGuide downloader in
# the matching mode where possible. Otherwise default both to species/binomial.
NCBI_MODE="species"
for ((i=0; i<${#NEOTOMA_ARGS[@]}; i++)); do
  if [[ "${NEOTOMA_ARGS[$i]}" == "--ncbi-name-mode" && $((i+1)) -lt ${#NEOTOMA_ARGS[@]} ]]; then
    case "${NEOTOMA_ARGS[$((i+1))]}" in
      binomial) NCBI_MODE="species" ;;
      trinomial) NCBI_MODE="trinomial" ;;
      as-is) NCBI_MODE="as-is" ;;
    esac
  fi
done

"$PYTHON_BIN" "${SCRIPT_DIR}/neotoma_extinct_to_gbif.py" \
  "${NEOTOMA_ARGS[@]}" \
  --out "$CSV_OUT" \
  --out-prefix "$OUTPREFIX" \
  --write-flyguide-files

FLY_ARGS=("--ncbi-name-mode" "$NCBI_MODE" "${FLYGUIDE_OPTS[@]}" "$CSV_OUT" "$OUTPREFIX" "$EMAIL")
[[ -n "$API_KEY" ]] && FLY_ARGS+=("$API_KEY")

bash "${SCRIPT_DIR}/flyguide.sh" "${FLY_ARGS[@]}"

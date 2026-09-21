#!/usr/bin/env bash
# Download (if needed), extract malignant counts, score, write 4/4 gene lists.
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
RAW="${1:-/tmp/concordant4_raw}"
MALIG="${2:-/tmp/concordant4_malig}"
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
export COMET_HERE="$HERE"
bash "$HERE/scripts/download.sh" "$RAW"
python3 "$HERE/scripts/extract_malignant.py" --raw "$RAW" --out "$MALIG" --cohort all
Rscript "$HERE/scripts/extract_gse205335.R" "$RAW" "$MALIG/GSE205335"
python3 "$HERE/scripts/extract_malignant.py" --raw "$RAW" --out "$MALIG" --cohort check
python3 "$HERE/scripts/score_comet.py" --malig "$MALIG" --results "$HERE/results/tables"
python3 "$HERE/scripts/concordance_enrichr.py" --malig "$MALIG" --tables "$HERE/results/tables"
echo "DONE $HERE"

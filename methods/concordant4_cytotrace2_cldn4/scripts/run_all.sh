#!/usr/bin/env bash
# Concordant-4 malignant CytoTRACE2 vs CLDN4. Data stay outside the repo.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DATA="${DATA:-/tmp/concordant4_ct2}"
PY="${PY:-python3}"
"$PY" "$ROOT/methods/concordant4_cytotrace2_cldn4/scripts/download.py" --out "$DATA"
"$PY" "$ROOT/methods/concordant4_cytotrace2_cldn4/scripts/extract_malignant.py" --data "$DATA" --out "$DATA/h5ad"
"$PY" "$ROOT/methods/concordant4_cytotrace2_cldn4/scripts/analyze.py" \
  --h5ad "$DATA/h5ad" \
  --outdir "$ROOT/methods/concordant4_cytotrace2_cldn4/results" \
  --finding "$ROOT/methods/concordant4_cytotrace2_cldn4/FINDING.md"

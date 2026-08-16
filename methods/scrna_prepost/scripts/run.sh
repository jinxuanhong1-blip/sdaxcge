#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DATA="$ROOT/data/scrna_prepost"
OUT="$ROOT/methods/scrna_prepost/tables"
mkdir -p "$OUT"
python3 "$ROOT/methods/scrna_prepost/scripts/analyze.py" \
  --gse207422-matrix "$DATA/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz" \
  --gse207422-metadata "$DATA/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx" \
  --gse337519-raw "$DATA/GSE337519/raw" \
  --outdir "$OUT"

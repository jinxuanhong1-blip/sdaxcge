#!/usr/bin/env bash
# Reproducible download of the open GEO datasets used by the fable_adc slice.
# Only public, processed files with outcome annotation are fetched. No controlled data.
set -euo pipefail

DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../results/fable_adc/raw" && pwd 2>/dev/null || true)"
if [[ -z "${DEST}" ]]; then
  DEST="$(dirname "${BASH_SOURCE[0]}")/../../results/fable_adc/raw"
  mkdir -p "${DEST}"
  DEST="$(cd "${DEST}" && pwd)"
fi
cd "${DEST}"

FILES=(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz"
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz"
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz"
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz"
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"
)

for url in "${FILES[@]}"; do
  fn="$(basename "$url")"
  if [[ -s "$fn" ]]; then
    echo "[skip] $fn already present"
    continue
  fi
  echo "[get ] $fn"
  for attempt in 1 2 3 4; do
    if curl -fsSL -o "$fn" "$url"; then
      echo "       ok ($(wc -c < "$fn") bytes)"
      break
    fi
    echo "       retry $attempt failed; backing off"
    sleep $((attempt * attempt))
  done
done

echo "Done. Files in: ${DEST}"

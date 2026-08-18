#!/usr/bin/env bash
# Public GSE131907 processed files only. Skip 2.86 GB log2TPM text and EGA FASTQ.
set -euo pipefail
OUT="${1:-/tmp/gse131907}"
mkdir -p "$OUT"
UA="sdaxcge-seurat-gse131907/1.0"
fetch() {
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then
    echo "exists $dest ($(wc -c < "$dest") bytes)"
    return
  fi
  echo "GET $url"
  curl -L --retry 5 --retry-delay 4 -A "$UA" -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "wrote $dest ($(wc -c < "$dest") bytes)"
}
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907"
fetch "$BASE/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz" "$OUT/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
fetch "$BASE/matrix/GSE131907_series_matrix.txt.gz" "$OUT/GSE131907_series_matrix.txt.gz"
fetch "$BASE/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" "$OUT/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
echo "skipped GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
echo "skipped EGA FASTQ EGAD00001005054"
echo "skipped GSE148071"

#!/usr/bin/env bash
# Public GSE131907 processed UMI + cell annotation. No log2TPM text, no EGA FASTQ.
set -euo pipefail
OUT="${1:-/tmp/gse131907}"
mkdir -p "$OUT"
UA="sdaxcge-gse131907-maxeffect/1.0"
fetch() {
  local url="$1" dest="$2" expect="${3:-}"
  if [[ -s "$dest" ]]; then
    echo "exists $dest ($(wc -c < "$dest") bytes)"
    return
  fi
  echo "GET $url"
  curl -L --retry 5 --retry-delay 4 -A "$UA" -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "wrote $dest ($(wc -c < "$dest") bytes)"
  if [[ -n "$expect" ]]; then
    local got
    got=$(wc -c < "$dest")
    if [[ "$got" != "$expect" ]]; then
      echo "size mismatch for $dest: got $got expected $expect" >&2
      exit 1
    fi
  fi
}
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907"
fetch "$BASE/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz" "$OUT/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
fetch "$BASE/matrix/GSE131907_series_matrix.txt.gz" "$OUT/GSE131907_series_matrix.txt.gz"
fetch "$BASE/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" "$OUT/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" 408736818
echo "skipped GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
echo "skipped EGA FASTQ EGAD00001005054"

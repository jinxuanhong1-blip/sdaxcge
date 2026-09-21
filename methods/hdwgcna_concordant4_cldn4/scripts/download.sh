#!/usr/bin/env bash
# Public GEO processed matrices for the concordant-4 cohort. No FASTQ.
set -euo pipefail
OUT="${1:-/tmp/geo_hdwgcna}"
mkdir -p "$OUT"
UA="sdaxcge-hdwgcna/1.0"
fetch() {
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then
    echo "exists $dest ($(wc -c < "$dest") bytes)"
    return
  fi
  echo "GET $url"
  curl -fL --retry 5 --retry-delay 4 -A "$UA" -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "wrote $dest ($(wc -c < "$dest") bytes)"
}
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar" \
  "$OUT/GSE123902_RAW.tar"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar" \
  "$OUT/GSE189357_RAW.tar"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  "$OUT/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" \
  "$OUT/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  "$OUT/GSE205335_Lung_IO_CellIdentity.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  "$OUT/GSE205335_Lung_IO_UMI_matrix.rds.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz" \
  "$OUT/GSE205335_family.soft.gz"

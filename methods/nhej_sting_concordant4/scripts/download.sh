#!/usr/bin/env bash
# Public GEO supplements for the concordant-4 malignant NHEJ / STING check.
# Not GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526.
set -euo pipefail
OUT="${1:-/tmp/geo_nhej}"
mkdir -p "$OUT"

dl() {
  local dest="$1" url="$2"
  if [[ -s "$dest" ]] && [[ "$(stat -c%s "$dest")" -gt 1000 ]]; then
    echo "HAVE $dest ($(stat -c%s "$dest") bytes)"
    return 0
  fi
  echo "GET $url"
  curl -fL --retry 6 --retry-delay 8 --retry-all-errors -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "OK $dest ($(stat -c%s "$dest") bytes)"
}

dl "$OUT/GSE123902_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
dl "$OUT/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
dl "$OUT/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz
dl "$OUT/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
dl "$OUT/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
dl "$OUT/GSE205335_family.soft.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz
dl "$OUT/GSE189357_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar

echo "done $OUT"
ls -lh "$OUT"

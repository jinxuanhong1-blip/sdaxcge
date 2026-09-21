#!/usr/bin/env bash
# Public concordant-4 matrices used by the Mixscape-like signature.
# Nothing here is private 8KL data. Matrices stay outside the git repo.
set -euo pipefail
OUT="${1:-/tmp/geo_mixscape}"
mkdir -p "$OUT/gse123902" "$OUT/gse189357"
UA="sdaxcge-mixscape-like-concordant4/1.0"
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
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar" "$OUT/GSE123902_RAW.tar"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar" "$OUT/GSE189357_RAW.tar"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz" "$OUT/GSE205335_Lung_IO_CellIdentity.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz" "$OUT/GSE205335_Lung_IO_UMI_matrix.rds.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz" "$OUT/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" "$OUT/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
if [[ ! -e "$OUT/gse123902/GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz" ]]; then
  tar -C "$OUT/gse123902" -xf "$OUT/GSE123902_RAW.tar"
fi
if [[ ! -e "$OUT/gse189357/GSM5699777_TD1_matrix.mtx.gz" ]]; then
  tar -C "$OUT/gse189357" -xf "$OUT/GSE189357_RAW.tar"
fi
echo "done $OUT"

#!/usr/bin/env bash
# Public downloads for TACSTD2–CLDN4 same-cell bridge.
set -euo pipefail
GEO="${1:-/tmp/geo_c4}"
COSMX_DIR="${2:-/workspace/data/cosmx_nsclc}"
mkdir -p "$GEO"/gse123902 "$GEO"/gse131907 "$GEO"/gse205335 "$GEO"/gse189357 "$COSMX_DIR"

dl() {
  local dest="$1" url="$2"
  if [[ -s "$dest" ]]; then
    echo "HAVE $dest ($(stat -c%s "$dest") bytes)"
    return 0
  fi
  echo "GET $url"
  curl -fL --retry 6 --retry-delay 8 --retry-all-errors -o "${dest}.partial" "$url"
  mv "${dest}.partial" "$dest"
  echo "OK $dest ($(stat -c%s "$dest") bytes)"
}

dl "$GEO/gse123902/GSE123902_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
dl "$GEO/gse131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
dl "$GEO/gse131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz
dl "$GEO/gse205335/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
dl "$GEO/gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
dl "$GEO/gse189357/GSE189357_RAW.tar" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_RAW.tar

COSMX="$COSMX_DIR/cosmx_human_nsclc_clustered.h5ad"
EXPECTED=2755776882
if [[ -f "$COSMX" && $(stat -c%s "$COSMX") -eq $EXPECTED ]]; then
  echo "HAVE $COSMX"
else
  dl "$COSMX" https://ndownloader.figshare.com/files/46841842
  got=$(stat -c%s "$COSMX")
  [[ "$got" -eq $EXPECTED ]] || { echo "CosMx size $got != $EXPECTED"; exit 1; }
fi
echo "done"

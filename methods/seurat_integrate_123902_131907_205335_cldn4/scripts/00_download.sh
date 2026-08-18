#!/usr/bin/env bash
# Public processed GEO only. Triple GSE123902 + GSE131907 + GSE205335.
# Not GSE148071. Not GSE189357.
set -euo pipefail
DEST="${1:-/tmp/triple_geo}"
mkdir -p "$DEST/GSE123902" "$DEST/GSE131907"

if [[ ! -s "$DEST/GSE123902/GSE123902_RAW.tar" ]]; then
  curl -fL --retry 4 --retry-delay 4 -o "$DEST/GSE123902/GSE123902_RAW.tar.partial" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
  mv "$DEST/GSE123902/GSE123902_RAW.tar.partial" "$DEST/GSE123902/GSE123902_RAW.tar"
fi
mkdir -p "$DEST/GSE123902/csv"
if ! ls "$DEST/GSE123902/csv"/GSM*_dense.csv.gz >/dev/null 2>&1; then
  tar -xf "$DEST/GSE123902/GSE123902_RAW.tar" -C "$DEST/GSE123902/csv"
fi

if [[ ! -s "$DEST/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz" ]]; then
  curl -fL --retry 4 --retry-delay 4 \
    -o "$DEST/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz
fi
if [[ ! -s "$DEST/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" ]]; then
  curl -fL --retry 4 --retry-delay 8 \
    -o "$DEST/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz.partial" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
  mv "$DEST/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz.partial" \
     "$DEST/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
fi

if [[ ! -s "$DEST/GSE205335_Lung_IO_CellIdentity.txt.gz" ]]; then
  curl -fL --retry 4 --retry-delay 4 -o "$DEST/GSE205335_Lung_IO_CellIdentity.txt.gz" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
fi
if [[ ! -s "$DEST/GSE205335_Lung_IO_UMI_matrix.rds.gz" ]]; then
  curl -fL --retry 4 --retry-delay 8 \
    -o "$DEST/GSE205335_Lung_IO_UMI_matrix.rds.gz.partial" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
  mv "$DEST/GSE205335_Lung_IO_UMI_matrix.rds.gz.partial" \
     "$DEST/GSE205335_Lung_IO_UMI_matrix.rds.gz"
fi

python3 - <<PY
import gzip
from pathlib import Path
src = Path("${DEST}/GSE205335_Lung_IO_UMI_matrix.rds.gz")
dst = Path("${DEST}/GSE205335_Lung_IO_UMI_matrix.rds")
if src.exists() and (not dst.exists() or dst.stat().st_size < 1000):
    with gzip.open(src, "rb") as f, dst.open("wb") as o:
        while True:
            chunk = f.read(16 * 1024 * 1024)
            if not chunk:
                break
            o.write(chunk)
    print("wrote", dst, dst.stat().st_size)
PY

ls -lh "$DEST" "$DEST/GSE123902" "$DEST/GSE131907"

#!/usr/bin/env bash
# Public processed GEO only. Pair GSE123902 + GSE205335.
# Do not download GSE148071 or GSE127465.
set -euo pipefail
DEST="${1:-/tmp/geo_pair_123902_205335}"
mkdir -p "$DEST/GSE123902/csv" "$DEST/GSE205335"

if [[ ! -s "$DEST/GSE123902/GSE123902_RAW.tar" ]]; then
  curl -fL --retry 5 --retry-delay 5 --max-time 600 \
    -o "$DEST/GSE123902/GSE123902_RAW.tar.partial" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar
  mv "$DEST/GSE123902/GSE123902_RAW.tar.partial" "$DEST/GSE123902/GSE123902_RAW.tar"
fi
if ! ls "$DEST/GSE123902/csv"/GSM*_dense.csv.gz >/dev/null 2>&1; then
  tar -xf "$DEST/GSE123902/GSE123902_RAW.tar" -C "$DEST/GSE123902/csv"
fi

if [[ ! -s "$DEST/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz" ]]; then
  curl -fL --retry 5 --retry-delay 5 --max-time 300 \
    -o "$DEST/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz.partial" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz
  mv "$DEST/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz.partial" \
     "$DEST/GSE205335/GSE205335_Lung_IO_CellIdentity.txt.gz"
fi
if [[ ! -s "$DEST/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz" ]]; then
  curl -fL --retry 5 --retry-delay 8 --max-time 1800 \
    -o "$DEST/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz.partial" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz
  mv "$DEST/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz.partial" \
     "$DEST/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz"
fi

# GEO ships a gzip-of-gzip RDS; one peel makes it readRDS-able.
python3 - <<PY
import gzip
from pathlib import Path
src = Path("${DEST}/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz")
dst = Path("${DEST}/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds")
if src.exists() and (not dst.exists() or dst.stat().st_size < 1000):
    with gzip.open(src, "rb") as f, dst.open("wb") as o:
        while True:
            chunk = f.read(16 * 1024 * 1024)
            if not chunk:
                break
            o.write(chunk)
    print("wrote", dst, dst.stat().st_size)
else:
    print("rds ready", dst, dst.stat().st_size if dst.exists() else 0)
PY

ls -lh "$DEST/GSE123902" "$DEST/GSE123902/csv" "$DEST/GSE205335" | head -40

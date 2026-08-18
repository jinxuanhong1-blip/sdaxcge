#!/usr/bin/env bash
# Public processed GEO only. Winning pair GSE131907 + GSE205335.
set -euo pipefail
DEST="${1:-/tmp/winpair_geo}"
mkdir -p "$DEST/GSE131907"
curl -fL --max-time 180 -o "$DEST/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
curl -fL --max-time 1200 -o "$DEST/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
curl -fL --max-time 180 -o "$DEST/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
curl -fL --max-time 1500 -o "$DEST/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
# GEO ships a gzip-of-gzip RDS; one peel makes it readRDS-able.
python3 - <<'PY'
import gzip
from pathlib import Path
src = Path("/tmp/winpair_geo/GSE205335_Lung_IO_UMI_matrix.rds.gz")
dst = Path("/tmp/winpair_geo/GSE205335_Lung_IO_UMI_matrix.rds")
if src.exists() and (not dst.exists() or dst.stat().st_size < 1000):
    with gzip.open(src, "rb") as f, dst.open("wb") as o:
        while True:
            chunk = f.read(16 * 1024 * 1024)
            if not chunk:
                break
            o.write(chunk)
    print("wrote", dst, dst.stat().st_size)
PY
ls -lh "$DEST" "$DEST/GSE131907"

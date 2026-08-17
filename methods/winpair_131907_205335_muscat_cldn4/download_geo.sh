#!/usr/bin/env bash
# Public processed GEO only. Winning pair GSE131907 + GSE205335.
set -euo pipefail
DEST="${1:-/tmp/winpair_geo}"
mkdir -p "$DEST/GSE131907"
curl -fL --max-time 180 -o "$DEST/GSE131907/cell_annotation.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
curl -fL --max-time 1200 -o "$DEST/GSE131907/raw_UMI_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
curl -fL --max-time 180 -o "$DEST/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
curl -fL --max-time 180 -o "$DEST/GSE205335_family.soft.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz"
curl -fL --max-time 1500 -o "$DEST/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
ls -lh "$DEST" "$DEST/GSE131907"

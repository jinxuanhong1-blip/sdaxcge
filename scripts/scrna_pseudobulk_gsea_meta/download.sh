#!/usr/bin/env bash
# Public processed GEO files only. GSE253013 9.3 GB RDS is intentionally not fetched.
set -euo pipefail
DEST="${1:-/tmp/scrna_pb_gsea_data}"
mkdir -p "$DEST/GSE131907" "$DEST/GSE241934" "$DEST/GSE291670"
cd "$DEST"

curl -fL --max-time 180 -o GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"
curl -fL --max-time 900 -o GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"

curl -fL --max-time 180 -o GSE205335_Lung_IO_CellIdentity.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
curl -fL --max-time 180 -o GSE205335_family.soft.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz"
curl -fL --max-time 1200 -o GSE205335_Lung_IO_UMI_matrix.rds.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"

curl -fL --max-time 180 -o GSE131907/cell_annotation.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
curl -fL --max-time 900 -o GSE131907/raw_UMI_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"

curl -fL --max-time 180 -o GSE241934/IIT_Meta.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz"
curl -fL --max-time 180 -o GSE241934/Real_Meta.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz"
curl -fL --max-time 60 -o GSE241934/IIT_features.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz"
curl -fL --max-time 60 -o GSE241934/IIT_barcodes.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz"
curl -fL --max-time 60 -o GSE241934/RWC_features.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz"
curl -fL --max-time 60 -o GSE241934/RWC_barcodes.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz"
curl -fL --max-time 1200 -o GSE241934/IIT_Matrix.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz"
curl -fL --max-time 1800 -o GSE241934/Real_Matrix.mtx.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz"

curl -fL --max-time 600 -o GSE291670/GSE291670_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291670/suppl/GSE291670_RAW.tar"
tar -xf GSE291670/GSE291670_RAW.tar -C GSE291670

curl -fL --max-time 60 -o GSE253013_series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/matrix/GSE253013_series_matrix.txt.gz"

echo "NOTE: GSE253013_all_luad_garnett_temp.rds.gz (9.3 GB) is not downloaded (16 GB RAM)."
ls -lh "$DEST" "$DEST/GSE131907" "$DEST/GSE241934" "$DEST/GSE291670"

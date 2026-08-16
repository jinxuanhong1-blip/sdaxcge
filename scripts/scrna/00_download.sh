#!/usr/bin/env bash
# Download processed GEO matrices only (no raw FASTQ). Files staged under /tmp/scrna_data.
set -euo pipefail
DEST=/tmp/scrna_data
mkdir -p "$DEST"
cd "$DEST"

curl -fL --max-time 180 -o GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"
curl -fL --max-time 600 -o GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"

curl -fL --max-time 180 -o GSE205335_Lung_IO_CellIdentity.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz"
curl -fL --max-time 900 -o GSE205335_Lung_IO_UMI_matrix.rds.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz"
curl -fL --max-time 180 -o GSE205335_family.soft.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz"

curl -fL --max-time 180 -o GSE271689_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/GSE271689_RAW.tar"
curl -fL --max-time 60 -o GSE271689_filelist.txt \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/filelist.txt"

ls -lh "$DEST"

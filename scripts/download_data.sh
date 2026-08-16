#!/usr/bin/env bash
# Download all inputs for Claims B3 + B4 into ./data
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data && cd data

echo "[1/3] TCGA-LUAD expression (UCSC Xena, log2 RSEM)..."
curl -sL -o TCGA_LUAD_HiSeqV2.gz \
  "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz"
gunzip -kf TCGA_LUAD_HiSeqV2.gz

echo "[2/3] GSE126044 counts + series matrix (GEO)..."
curl -sL -o GSE126044_counts.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz"
gunzip -kf GSE126044_counts.txt.gz
curl -sL -o GSE126044_series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz"
gunzip -kf GSE126044_series_matrix.txt.gz

echo "[3/3] KEGG tight junction pathway (hsa04530)..."
curl -sL -o kegg_hsa04530.txt "https://rest.kegg.jp/get/hsa04530"

echo "Done. Files in $(pwd):"
ls -la

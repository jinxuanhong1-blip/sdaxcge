#!/usr/bin/env bash
# Download processed GEO data for the fable_blood_ici slice.
#
# Datasets (all human NSCLC, blood compartment, ICI context):
#   GSE285888  scRNA-seq of baseline PBMC from NSCLC patients prior to ICI
#              (response groups CR/DR/PD + irAE severity). PMID 40404203.
#   GSE305086  Whole-blood bulk gene-expression array under immunotherapy
#              as a candidate ICB biomarker in NSCLC (baseline/follow-up/controls).
#
# All downloads are processed matrices; each file is < 2 GB.
# Raw data is written to $GEO_DIR (default /tmp/geo) and is NOT committed.
set -euo pipefail

GEO_DIR="${GEO_DIR:-/tmp/geo}"
mkdir -p "$GEO_DIR"
cd "$GEO_DIR"

base="https://ftp.ncbi.nlm.nih.gov/geo"

echo "[1/4] GSE285888 scRNA-seq cell metadata (~5.5 MB)"
curl -fsS -o GSM8712033_metadata.csv.gz \
  "$base/samples/GSM8712nnn/GSM8712033/suppl/GSM8712033_metadata.csv.gz"

echo "[2/4] GSE285888 scRNA-seq expression matrix (~486 MB, genes x cells)"
curl -fsS -o GSM8712033_matrix.txt.gz \
  "$base/samples/GSM8712nnn/GSM8712033/suppl/GSM8712033_matrix.txt.gz"

echo "[3/4] GSE305086 whole-blood expression matrix (~25 MB)"
curl -fsS -o GSE305086_Expression_matrix_final.csv.gz \
  "$base/series/GSE305nnn/GSE305086/suppl/GSE305086_Expression_matrix_final.csv.gz"

echo "[4/4] GSE305086 series matrix (sample annotation)"
curl -fsS -o GSE305086_series_matrix.txt.gz \
  "$base/series/GSE305nnn/GSE305086/matrix/GSE305086_series_matrix.txt.gz"

echo "Done. Files in $GEO_DIR:"
ls -la "$GEO_DIR"

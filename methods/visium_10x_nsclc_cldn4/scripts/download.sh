#!/usr/bin/env bash
set -euo pipefail
ROOT="${VISIUM_DATA:-/workspace/data}"
mkdir -p "$ROOT/10x_visium" "$ROOT/geo/GSE189487" "$ROOT/geo/GSE273378" "$ROOT/geo/GSE300676"

cd "$ROOT/10x_visium"
curl -fsSL -o CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5 \
  "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
curl -fsSL -o CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz \
  "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz"
curl -fsSL -o CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5 \
  "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
curl -fsSL -o CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz \
  "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz"
mkdir -p LUSC_FFPE NEC_11mm
tar -xzf CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz -C LUSC_FFPE
tar -xzf CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz -C NEC_11mm

curl -fsSL -o "$ROOT/geo/GSE189487/GSE189487_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189487/suppl/GSE189487_RAW.tar"
mkdir -p "$ROOT/geo/GSE189487/raw"
tar -xf "$ROOT/geo/GSE189487/GSE189487_RAW.tar" -C "$ROOT/geo/GSE189487/raw"

# Optional additional open GEO Visium LUAD (large)
curl -fsSL -o "$ROOT/geo/GSE273378/GSE273378_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273378/suppl/GSE273378_RAW.tar" || true
curl -fsSL -o "$ROOT/geo/GSE300676/GSE300676_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300676/suppl/GSE300676_RAW.tar" || true

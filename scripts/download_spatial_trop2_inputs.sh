#!/usr/bin/env bash
# Public inputs for scripts/spatial_trop2_cd8_after_cldn4.py.
# Raw matrices stay outside the repo (default /tmp).
set -euo pipefail
ROOT="${1:-/tmp}"
mkdir -p "$ROOT/cosmx" "$ROOT/visium"

if [[ ! -s "$ROOT/cosmx/cosmx_lung.zip" ]]; then
  curl -fL --retry 3 -o "$ROOT/cosmx/cosmx_lung.zip" \
    "https://zenodo.org/records/15487520/files/cosmx_lung.zip?download=1"
fi

BASE10="https://cf.10xgenomics.com/samples/spatial-exp"
curl -fL --retry 3 -o "$ROOT/visium/lusc.h5" \
  "$BASE10/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
curl -fL --retry 3 -o "$ROOT/visium/lusc_spatial.tar.gz" \
  "$BASE10/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz"
curl -fL --retry 3 -o "$ROOT/visium/nec.h5" \
  "$BASE10/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
curl -fL --retry 3 -o "$ROOT/visium/nec_spatial.tar.gz" \
  "$BASE10/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz"

curl -fL --retry 3 -o "$ROOT/visium/GSE189487_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189487/suppl/GSE189487_RAW.tar"
curl -fL --retry 3 -o "$ROOT/visium/GSE273378_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273378/suppl/GSE273378_RAW.tar"

echo "inputs in $ROOT (not committed)"

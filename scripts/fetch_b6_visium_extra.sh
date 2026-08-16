#!/usr/bin/env bash
# Fetch only the small processed 10x Visium lung SCC files used by B6.
# Deliberately does NOT download huge raw / image / Loupe files.
# Dataset: https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard
# NOT E-MTAB-13530.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${ROOT}/data/CytAssist_FFPE_Human_Lung_SCC"
BASE="https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"
mkdir -p "${DEST}"
cd "${DEST}"

echo "Downloading filtered matrix (~24 MB) and spatial folder (~31 MB)..."
curl -fL -o CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5 \
  "${BASE}_filtered_feature_bc_matrix.h5"
curl -fL -o CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz \
  "${BASE}_spatial.tar.gz"
curl -fL -o CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_metrics_summary.csv \
  "${BASE}_metrics_summary.csv"
tar xzf CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz

echo "SKIPPED (huge / unused raw):"
echo "  tissue_image.tif      ~3.0 GB  (full-res H&E)"
echo "  *.cloupe              ~1.6 GB  (Loupe browser)"
echo "  molecule_info.h5      ~377 MB"
echo "  FASTQ / raw FASTQ     multi-GB (not fetched)"
echo "  raw_feature_bc_matrix.h5  ~31 MB (unfiltered; not needed)"
echo "Done. Layout: ${DEST}"

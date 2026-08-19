#!/usr/bin/env bash
# Fetch compact processed Visium files only (no FASTQ, no full-res H&E, no GSE307534 9.4 GB RAW).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="${ROOT}/data"
mkdir -p "${DATA}"/{GSE277206,GSE189487,GSE273378,tenx_lung_scc,tenx_lung_11mm,zenodo_13337961}

retry() { curl -fL --retry 4 --retry-delay 4 "$@"; }

echo "== GSE277206 (2 CytAssist FFPE LUAD sections, 51 MB) =="
retry -o "${DATA}/GSE277206/GSE277206_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277206/suppl/GSE277206_RAW.tar"
tar xf "${DATA}/GSE277206/GSE277206_RAW.tar" -C "${DATA}/GSE277206"

echo "== GSE189487 (6 early LUAD Visium sections, 210 MB) =="
retry -o "${DATA}/GSE189487/GSE189487_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189487/suppl/GSE189487_RAW.tar"
tar xf "${DATA}/GSE189487/GSE189487_RAW.tar" -C "${DATA}/GSE189487"

echo "== GSE273378 (16 stage-I LUAD Visium sections; matrices+coords only) =="
retry -o "${DATA}/GSE273378/GSE273378_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273378/suppl/GSE273378_RAW.tar"
tar xf "${DATA}/GSE273378/GSE273378_RAW.tar" -C "${DATA}/GSE273378" \
  --wildcards '*matrix.mtx*' '*barcodes.tsv*' '*features.tsv*' '*tissue_positions*' || true

echo "== 10x NSCLC SCC FFPE (6.5 mm CytAssist) =="
BASE="https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"
retry -o "${DATA}/tenx_lung_scc/filtered_feature_bc_matrix.h5" "${BASE}_filtered_feature_bc_matrix.h5"
retry -o "${DATA}/tenx_lung_scc/spatial.tar.gz" "${BASE}_spatial.tar.gz"
tar xzf "${DATA}/tenx_lung_scc/spatial.tar.gz" -C "${DATA}/tenx_lung_scc"

echo "== 10x NSCLC 11 mm neuroendocrine =="
BASE11="https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer"
retry -o "${DATA}/tenx_lung_11mm/filtered_feature_bc_matrix.h5" "${BASE11}_filtered_feature_bc_matrix.h5"
retry -o "${DATA}/tenx_lung_11mm/spatial.tar.gz" "${BASE11}_spatial.tar.gz"
tar xzf "${DATA}/tenx_lung_11mm/spatial.tar.gz" -C "${DATA}/tenx_lung_11mm"

echo "== Zenodo 13337961 (lepidic + solid LUAD matrices; skip pathology images) =="
retry -o "${DATA}/zenodo_13337961/lepidic-filtered_feature_bc_matrix.rar" \
  "https://zenodo.org/records/13337961/files/lepidic-filtered_feature_bc_matrix.rar?download=1"
retry -o "${DATA}/zenodo_13337961/solid-filtered_feature_bc_matrix.rar" \
  "https://zenodo.org/records/13337961/files/solid-filtered_feature_bc_matrix.rar?download=1"
unar -o "${DATA}/zenodo_13337961/lepidic" "${DATA}/zenodo_13337961/lepidic-filtered_feature_bc_matrix.rar"
unar -o "${DATA}/zenodo_13337961/solid" "${DATA}/zenodo_13337961/solid-filtered_feature_bc_matrix.rar"

echo "SKIPPED: GSE307534 series RAW 9.4 GB (per-sample tars 50–350 MB, image-heavy). Not re-downloaded."
echo "Done."

#!/usr/bin/env bash
# Download public GSE271689 + paper source data. Large files stay local (see .gitignore).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"
mkdir -p "$DATA/paper" "$DATA/dcc" "$DATA/repo"

curl -fsSL -o "$DATA/GSE271689_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/matrix/GSE271689_series_matrix.txt.gz"
curl -fsSL -o "$DATA/GSE271689_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271689/suppl/GSE271689_RAW.tar"
tar -xf "$DATA/GSE271689_RAW.tar" -C "$DATA/dcc"
curl -fsSL -o "$DATA/Hs_R_NGS_WTA_v1.0.pkc" \
  "https://zenodo.org/records/12752405/files/Hs_R_NGS_WTA_v1.0.pkc?download=1"
curl -fsSL -A "Mozilla/5.0" -o "$DATA/paper/41588_2025_2351_MOESM10_ESM.xlsx" \
  "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02351-7/MediaObjects/41588_2025_2351_MOESM10_ESM.xlsx"
curl -fsSL -o "$DATA/repo/CK_av_3B.csv" \
  "https://raw.githubusercontent.com/tznaung/NSCLC_SpatialOmics/main/data/CK_av_3B.csv"

echo "next: python3 scripts/hunt_gse271689/01_build_matrix.py"
echo "      python3 scripts/hunt_gse271689/02_assemble_cohorts.py"
echo "      python3 scripts/hunt_gse271689/03_analyze_trop2_os.py"

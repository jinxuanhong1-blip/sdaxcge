#!/usr/bin/env bash
# Re-download the open source data used in this hunt.
# Run from repo root: bash results/hunt_oak/scripts/fetch_data.sh
set -euo pipefail
RAW="results/hunt_oak/data/raw"
mkdir -p "$RAW"

# --- IMvigor210 (atezolizumab, mUC) authoritative CountDataSet -------------
# Original package research-pub.gene.com/IMvigor210CoreBiologies is offline (404
# as of 2026-08); use the community mirror that preserves the CountDataSet with
# real Entrez/symbol annotations.
curl -sSL -o "$RAW/imvigor210_cds.RData" \
  "https://raw.githubusercontent.com/SiYangming/IMvigor210CoreBiologies/master/data/cds.RData"

# --- Bessede 2024 supplements (aggregate tables only; cannot re-test) -----
mkdir -p "$RAW/bessede_supplements"
# figshare file IDs from collection 10.1158/1078-0432.c.7077754
curl -sSL -o "$RAW/bessede_supplements/ccr-23-2566_supplementary_table_s1_suppts1.docx" \
  "https://ndownloader.figshare.com/files/44570659"
curl -sSL -o "$RAW/bessede_supplements/ccr-23-2566_supplementary_table_s2_suppts2.docx" \
  "https://ndownloader.figshare.com/files/44570656"

# --- Gandara 2018 Nat Med ESM (OAK/POPLAR clinical + bTMB; no TACSTD2) ----
mkdir -p "$RAW/gandara2018"
curl -sSL -o "$RAW/gandara2018/41591_2018_134_MOESM3_ESM.xlsx" \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-018-0134-3/MediaObjects/41591_2018_134_MOESM3_ESM.xlsx"

# --- GSE135222 (NSCLC, anti-PD-1/PD-L1) TPM + metadata --------------------
curl -sSL -o "$RAW/GSE135222_exp.tsv.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
curl -sSL -o "$RAW/GSE135222_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz"
gunzip -kf "$RAW/GSE135222_exp.tsv.gz" "$RAW/GSE135222_series_matrix.txt.gz"

echo "Downloaded. Now run:"
echo "  Rscript results/hunt_oak/scripts/00_extract_imvigor210_cds.R"
echo "  python3 results/hunt_oak/scripts/01_imvigor210.py"
echo "  python3 results/hunt_oak/scripts/02_gse135222.py"
echo "  python3 results/hunt_oak/scripts/03_figures.py"

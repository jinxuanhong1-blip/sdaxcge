#!/usr/bin/env bash
# Download open processed data for the fable_china_ici slice.
# Total download size: ~630 MB (well under the 2 GB budget).
# Large raw matrices are cached outside the repo (default /tmp/fable_china_ici_data);
# only derived tables/figures are written into results/fable_china_ici/.
#
# Usage: bash scripts/fable_china_ici/01_download_data.sh [core|all]
#   core (default): bulk/processed matrices + metadata + series matrices (~10 MB)
#   all:            core + GSE207422 scRNA UMI matrix (175 MB) + GSE241934 IIT scRNA (447 MB)
set -euo pipefail

DATA_DIR="${DATA_DIR:-/tmp/fable_china_ici_data}"
MODE="${1:-core}"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

dl () { # dl <url> <out>
  if [ -s "$2" ]; then echo "exists: $2"; else
    echo "downloading: $2"
    curl -fsSL --retry 4 --retry-delay 5 -o "$2" "$1"
  fi
}

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series

# --- GSE207422: Shanghai Pulmonary Hospital, neoadjuvant anti-PD-1 + chemo, resectable NSCLC ---
dl "$GEO/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz"  GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz
dl "$GEO/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx"   GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx
dl "$GEO/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx"      GSE207422_NSCLC_scRNAseq_metadata.xlsx

# --- GSE126044: Yonsei (South Korea), anti-PD-1 monotherapy, advanced NSCLC, R vs NR ---
dl "$GEO/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz"                     GSE126044_counts.txt.gz
dl "$GEO/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz"             GSE126044_series_matrix.txt.gz

# --- GSE135222: South Korea, anti-PD-1/PD-L1, advanced NSCLC, clinical benefit + PFS ---
dl "$GEO/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"   GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz
dl "$GEO/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz"             GSE135222_series_matrix.txt.gz

# --- GSE260770: Guangzhou Medical University, sintilimab on GGO lesions, exosomal mRNA (exploratory) ---
dl "$GEO/GSE260nnn/GSE260770/suppl/GSE260770_mRNA_FPKM.txt.gz"                  GSE260770_mRNA_FPKM.txt.gz
dl "$GEO/GSE260nnn/GSE260770/matrix/GSE260770_series_matrix.txt.gz"             GSE260770_series_matrix.txt.gz

if [ "$MODE" = "all" ]; then
  # GSE207422 scRNA UMI matrix (175 MB, dense text)
  dl "$GEO/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz" GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz
  # GSE241934 (NEOTIDE/CTONG2104, Peking University) scRNA:
  #   IIT cohort  = EGFR-mutant trial arm, neoadjuvant sintilimab + chemo (444 MB mtx)
  #   Real cohort = real-world EGFR-WT neoadjuvant immunochemotherapy (1.2 GB mtx)
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz"    GSE241934_IIT_Matrix.mtx.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz"      GSE241934_IIT_Meta.txt.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz"  GSE241934_IIT_barcodes.tsv.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz"  GSE241934_IIT_features.tsv.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz"   GSE241934_Real_Matrix.mtx.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz"     GSE241934_Real_Meta.txt.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz"  GSE241934_RWC_barcodes.tsv.gz
  dl "$GEO/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz"  GSE241934_RWC_features.tsv.gz
fi

echo "--- files in $DATA_DIR ---"
ls -lh "$DATA_DIR"

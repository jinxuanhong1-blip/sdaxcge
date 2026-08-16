#!/usr/bin/env bash
# Public GEO downloads only. GSE253013 (9.3 GB Seurat RDS) is not fetched
# unless explicitly requested — it is not used for scVI HVG integration.
set -euo pipefail
OUT="${1:-/tmp/scrna_scvi_combo}"
mkdir -p "$OUT"
cd "$OUT"

base207="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl"
base241="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl"

fetch() {
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then
    echo "exists $dest"
    return 0
  fi
  curl -L --retry 5 --retry-delay 4 -o "$dest" "$url"
}

fetch "$base207/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz" \
  GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz
fetch "$base207/GSE207422_NSCLC_scRNAseq_metadata.xlsx" \
  GSE207422_NSCLC_scRNAseq_metadata.xlsx

fetch "$base241/GSE241934_IIT_Matrix.mtx.gz" GSE241934_IIT_Matrix.mtx.gz
fetch "$base241/GSE241934_IIT_Meta.txt.gz" GSE241934_IIT_Meta.txt.gz
fetch "$base241/GSE241934_IIT_barcodes.tsv.gz" GSE241934_IIT_barcodes.tsv.gz
fetch "$base241/GSE241934_IIT_features.tsv.gz" GSE241934_IIT_features.tsv.gz

# Optional real-world slice of the same series (not a second GSE).
if [[ "${FETCH_GSE241934_RWC:-0}" == "1" ]]; then
  fetch "$base241/GSE241934_Real_Matrix.mtx.gz" GSE241934_Real_Matrix.mtx.gz
  fetch "$base241/GSE241934_Real_Meta.txt.gz" GSE241934_Real_Meta.txt.gz
  fetch "$base241/GSE241934_RWC_barcodes.tsv.gz" GSE241934_RWC_barcodes.tsv.gz
  fetch "$base241/GSE241934_RWC_features.tsv.gz" GSE241934_RWC_features.tsv.gz
fi

ls -lh

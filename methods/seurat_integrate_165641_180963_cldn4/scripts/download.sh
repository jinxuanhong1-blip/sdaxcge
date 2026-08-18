#!/usr/bin/env bash
# Public GEO processed 10x MTX only. No SRA. No private 8-KL.
set -euo pipefail
OUT="${1:-/tmp/geo/work}"
mkdir -p "${OUT}/GSE165641" "${OUT}/GSE180963"

dl() {
  local url="$1" dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "exists ${dest}"
    return 0
  fi
  echo "GET ${url}"
  wget -q -O "${dest}" "${url}"
}

# GSE165641 — 2 KL mice, Cell Ranger filtered MTX
dl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5047302&format=file&file=GSM5047302_KL1_count.tar.gz" \
  "${OUT}/GSE165641/GSM5047302_KL1_count.tar.gz"
dl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5047303&format=file&file=GSM5047303_KL2_count.tar.gz" \
  "${OUT}/GSE165641/GSM5047303_KL2_count.tar.gz"

# GSE180963 — 1 K + 1 KL, author-QC 10x MTX
dl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5481386&format=file&file=GSM5481386_K.tar.gz" \
  "${OUT}/GSE180963/GSM5481386_K.tar.gz"
dl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5481387&format=file&file=GSM5481387_KL.tar.gz" \
  "${OUT}/GSE180963/GSM5481387_KL.tar.gz"

if [[ ! -d "${OUT}/GSE165641/KL1_count/filtered_feature_bc_matrix" ]]; then
  tar -xzf "${OUT}/GSE165641/GSM5047302_KL1_count.tar.gz" -C "${OUT}/GSE165641"
fi
if [[ ! -d "${OUT}/GSE165641/KL2_count/filtered_feature_bc_matrix" ]]; then
  tar -xzf "${OUT}/GSE165641/GSM5047303_KL2_count.tar.gz" -C "${OUT}/GSE165641"
fi
if [[ ! -f "${OUT}/GSE180963/K/matrix.mtx" ]]; then
  tar -xzf "${OUT}/GSE180963/GSM5481386_K.tar.gz" -C "${OUT}/GSE180963"
fi
if [[ ! -f "${OUT}/GSE180963/KL/matrix.mtx" ]]; then
  tar -xzf "${OUT}/GSE180963/GSM5481387_KL.tar.gz" -C "${OUT}/GSE180963"
fi

echo "downloaded under ${OUT}"
find "${OUT}" -name 'matrix.mtx*' -o -name 'barcodes.tsv*' -o -name 'genes.tsv*' -o -name 'features.tsv*' | sort

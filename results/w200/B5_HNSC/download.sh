#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${ROOT}/data"

mkdir -p \
  "${DATA}/raw" \
  "${DATA}/candidate_GSE93157" \
  "${DATA}/candidate_GSE190575" \
  "${DATA}/candidate_GSE212549" \
  "${DATA}/candidate_GSE301741"

download() {
  local url="$1"
  local destination="$2"
  if [[ -s "${destination}" ]]; then
    return
  fi
  curl --fail --location --retry 4 --retry-all-errors \
    --output "${destination}.part" "${url}"
  mv "${destination}.part" "${destination}"
}

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159067/matrix/GSE159067_series_matrix.txt.gz" \
  "${DATA}/raw/GSE159067_series_matrix.txt.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159067/suppl/GSE159067_IHN_log2cpm_data.txt.gz" \
  "${DATA}/raw/GSE159067_IHN_log2cpm_data.txt.gz"

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/matrix/GSE93157_series_matrix.txt.gz" \
  "${DATA}/candidate_GSE93157/GSE93157_series_matrix.txt.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/suppl/GSE93157_raw_data_values.txt.gz" \
  "${DATA}/candidate_GSE93157/GSE93157_raw_data_values.txt.gz"

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190575/matrix/GSE190575_series_matrix.txt.gz" \
  "${DATA}/candidate_GSE190575/GSE190575_series_matrix.txt.gz"

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212549/matrix/GSE212549_series_matrix.txt.gz" \
  "${DATA}/candidate_GSE212549/GSE212549_series_matrix.txt.gz"

ANNOTATION_TAR="${DATA}/candidate_GSE212549/clariomdhumantranscriptcluster.db_8.8.0.tar.gz"
ANNOTATION_SQLITE="${DATA}/candidate_GSE212549/clariomdhumantranscriptcluster.db/inst/extdata/clariomdhumantranscriptcluster.sqlite"
if [[ ! -s "${ANNOTATION_SQLITE}" ]]; then
  download \
    "https://bioconductor.org/packages/release/data/annotation/src/contrib/clariomdhumantranscriptcluster.db_8.8.0.tar.gz" \
    "${ANNOTATION_TAR}"
  tar -xzf "${ANNOTATION_TAR}" \
    -C "${DATA}/candidate_GSE212549" \
    "clariomdhumantranscriptcluster.db/inst/extdata/clariomdhumantranscriptcluster.sqlite"
  rm -f "${ANNOTATION_TAR}"
fi

download \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301741/matrix/GSE301741_series_matrix.txt.gz" \
  "${DATA}/candidate_GSE301741/GSE301741_series_matrix.txt.gz"

python3 "${ROOT}/analyze.py"

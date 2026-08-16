#!/usr/bin/env bash
# Fetch the public GSE283829 (NSCLC, ICI-treated) expression matrix and metadata from NCBI GEO.
# Everything written here is gitignored; re-run to reproduce the inputs.
set -euo pipefail

DATA_DIR="${1:-data/GSE283829}"
mkdir -p "$DATA_DIR"

BASE_ACC="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
SUPPL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283829/suppl/GSE283829_raw_express_matrix_all_samples.txt.gz"
HGNC="https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.tsv"

curl -sSL --retry 4 --retry-delay 4 -m 300 \
  "${BASE_ACC}?acc=GSE283829&targ=self&form=text&view=brief" \
  -o "$DATA_DIR/GSE283829_series.txt"

curl -sSL --retry 4 --retry-delay 4 -m 300 \
  "${BASE_ACC}?acc=GSE283829&targ=gsm&form=text&view=brief" \
  -o "$DATA_DIR/GSE283829_samples.txt"

curl -sSL --retry 4 --retry-delay 4 -m 600 \
  "$SUPPL" -o "$DATA_DIR/GSE283829_raw_express_matrix_all_samples.txt.gz"

curl -sSL --retry 4 --retry-delay 4 -m 600 \
  "$HGNC" -o "$DATA_DIR/hgnc_complete_set.tsv"

echo "Downloaded into $DATA_DIR:"
( cd "$DATA_DIR" && md5sum ./* )

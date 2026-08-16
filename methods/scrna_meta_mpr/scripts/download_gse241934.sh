#!/usr/bin/env bash
set -euo pipefail
DEST=/tmp/gse241934
mkdir -p "$DEST"
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl
cd "$DEST"
# small files first
for f in \
  GSE241934_IIT_features.tsv.gz \
  GSE241934_IIT_barcodes.tsv.gz \
  GSE241934_IIT_Meta.txt.gz \
  GSE241934_RWC_features.tsv.gz \
  GSE241934_RWC_barcodes.tsv.gz \
  GSE241934_Real_Meta.txt.gz \
  GSE241934_IIT_Matrix.mtx.gz \
  GSE241934_Real_Matrix.mtx.gz
do
  if [[ -s "$f" ]]; then
    echo "have $f"
    continue
  fi
  echo "GET $f"
  curl -L --retry 5 --retry-delay 4 -o "$f.part" "$BASE/$f"
  mv "$f.part" "$f"
done
ls -lh

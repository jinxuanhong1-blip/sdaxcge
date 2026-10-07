#!/usr/bin/env bash
# Re-download the author-provided PROCESSED tables (small FPKM matrices) used here.
# No raw/multi-GB matrices. Run from the methods/.../ directory.
set -euo pipefail
cd "$(dirname "$0")/../data"
base="https://ftp.ncbi.nlm.nih.gov/geo/series"
for gse in GSE304294 GSE312098 GSE311016; do
  pre="$(echo "$gse" | sed -E 's/[0-9]{3}$//')nnn"
  url="${base}/${pre}/${gse}/suppl/${gse}_gene_fpkm.txt.gz"
  echo "fetching $gse ..."
  curl -fsSL --retry 4 --max-time 120 -o "${gse}_gene_fpkm.txt.gz" "$url"
  ls -la "${gse}_gene_fpkm.txt.gz"
done
echo "done."

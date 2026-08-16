#!/usr/bin/env bash
# Download processed (only) matrices/DE tables for the verified TROP2/TACSTD2
# loss-of-function datasets. All files are << 2 GB. Raw FASTQ / >2GB files skipped.
set -euo pipefail
cd "$(dirname "$0")/../../results/fable_tacstd2_kdko"
mkdir -p data
cd data

dl() { # url outfile
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then echo "[skip] $out exists"; return; fi
  echo "[get ] $out"
  for i in 1 2 3 4; do
    if curl -fsSL -m 300 "$url" -o "$out"; then return; fi
    echo "  retry $i..."; sleep $((2**i))
  done
  echo "FAILED: $url" >&2; return 1
}

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series"

# GSE334497 (mouse 4T1 Trop2 KO vs WT) - normalized counts
dl "$BASE/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz" GSE334497_normalized_counts.csv.gz

# GSE289287 (human T-47D) - author DESeq2 tables (Trop2 KO + DSG2 KO)
for f in \
  GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz \
  GSE289287_DESeq2-DSG2KO_tumors_vs_WT.tsv.gz \
  GSE289287_DESeq2-DSG2KO_cells_vs_DSG2WT.tsv.gz ; do
  dl "$BASE/GSE289nnn/GSE289287/suppl/$f" "$f"
done

# GSE245459 (human ovarian shTACSTD2) - FPKM
dl "$BASE/GSE245nnn/GSE245459/suppl/GSE245459_fpkm.anno.txt.gz" GSE245459_fpkm.anno.txt.gz

echo "[done] downloads in $(pwd)"
ls -la

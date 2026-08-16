#!/usr/bin/env bash
# Download every processed GEO supplement for GSE131907 (Kim et al. 2020).
# Size is not a skip reason. Both the raw UMI text matrix and the authors'
# normalized log2TPM text matrix are fetched in full. RDS copies are the same
# matrices in a different container; the text files are used because this
# environment has no R runtime.
set -euo pipefail

DEST="${1:-$(cd "$(dirname "$0")/.." && pwd)/data}"
mkdir -p "$DEST"

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907"
FILES=(
  "$BASE/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
  "$BASE/suppl/GSE131907_Lung_Cancer_Feature_Summary.xlsx"
  "$BASE/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
  "$BASE/suppl/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
  "$BASE/matrix/GSE131907_series_matrix.txt.gz"
)

for url in "${FILES[@]}"; do
  f="$(basename "$url")"
  echo "=== $f"
  wget -c --tries=5 --waitretry=10 --timeout=60 --progress=dot:giga \
       -O "$DEST/$f" "$url"
done

echo "=== done"
ls -l "$DEST"

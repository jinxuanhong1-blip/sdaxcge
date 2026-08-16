#!/usr/bin/env bash
# Extract header + lineage-marker rows from GSE207422 dense genes-x-cells TSV.gz.
set -eu
SRC="${1:?matrix gz}"
OUT="${2:?out gz}"
zcat "$SRC" | awk -F'\t' 'NR==1 || $1 ~ /^(TACSTD2|EPCAM|KRT19|KRT18|KRT8|KRT7|KRT5|KRT17|NAPSA|NKX2-1|PTPRC|CD3D|CD3E|CD3G|CD2|CD8A|CD4|IL7R|NKG7|GNLY|KLRD1|NCAM1|GZMB|PRF1|MS4A1|CD79A|LYZ|CD68|PECAM1|COL1A1)$/' | gzip -c > "$OUT"
echo "wrote $OUT"

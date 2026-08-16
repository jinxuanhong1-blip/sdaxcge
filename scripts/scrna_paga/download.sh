#!/usr/bin/env bash
# Download processed public GSE131907 files for the PAGA slice.
# Raw FASTQ / EGA and the 2.9 GB log2TPM text matrix are refused.
set -euo pipefail

DEST="${1:-/tmp/scrna_paga_data}"
mkdir -p "$DEST"

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907"
UA="sdaxcge-scrna-paga/1.0 (public GEO extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

dl() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "skip exists: $out ($(wc -c < "$out") bytes)"
    return 0
  fi
  echo "GET $url"
  local tmp="${out}.part"
  local n=0
  while [[ $n -lt 5 ]]; do
    if curl -fL --retry 3 --retry-delay 8 --connect-timeout 60 \
      -A "$UA" -o "$tmp" "$url"; then
      mv -f "$tmp" "$out"
      echo "wrote $out ($(wc -c < "$out") bytes)"
      return 0
    fi
    n=$((n + 1))
    echo "retry $n: $url" >&2
    sleep $((4 * n))
  done
  echo "FAILED $url" >&2
  return 1
}

# Per-cell author annotation (1.8 MB) + series matrix (sample meta).
dl "${BASE}/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz" \
  "${DEST}/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
dl "${BASE}/matrix/GSE131907_series_matrix.txt.gz" \
  "${DEST}/GSE131907_series_matrix.txt.gz"

# Raw UMI text (~390 MB). Refuse the 2.9 GB log2TPM text and both RDS blobs.
UMI="${DEST}/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
if [[ -s "$UMI" ]]; then
  echo "skip exists: $UMI ($(wc -c < "$UMI") bytes)"
else
  dl "${BASE}/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz" "$UMI"
fi

python3 - <<'PY' "$DEST"
import hashlib, json, os, sys
dest = sys.argv[1]
rows = []
for name in sorted(os.listdir(dest)):
    if name.endswith(".part"):
        continue
    path = os.path.join(dest, name)
    if not os.path.isfile(path):
        continue
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    rows.append({
        "file": name,
        "bytes": os.path.getsize(path),
        "sha256": h.hexdigest(),
    })
out = os.path.join(dest, "DOWNLOAD_MANIFEST.json")
with open(out, "w") as f:
    json.dump({"accession": "GSE131907", "files": rows}, f, indent=2)
print("manifest", out)
PY

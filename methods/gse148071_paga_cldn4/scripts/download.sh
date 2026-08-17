#!/usr/bin/env bash
# Public GSE148071 UMI matrices + TISCH2 cell labels (annotation only).
# Raw FASTQ are not on GEO. TISCH expression.h5 is not required.
set -euo pipefail

DEST="${1:-/tmp/gse148071_paga_data}"
mkdir -p "$DEST"

GEO="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071"
TISCH="https://tisch.compbio.cn/static/data/NSCLC_GSE148071"
UA="sdaxcge-gse148071-paga-cldn4/1.0 (public GEO extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

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

dl "${GEO}/suppl/GSE148071_RAW.tar" "${DEST}/GSE148071_RAW.tar"
dl "${GEO}/matrix/GSE148071_series_matrix.txt.gz" "${DEST}/GSE148071_series_matrix.txt.gz"
dl "${TISCH}/NSCLC_GSE148071_CellMetainfo_table.tsv" \
  "${DEST}/NSCLC_GSE148071_CellMetainfo_table.tsv"

python3 - <<'PY' "$DEST"
import hashlib, json, os, sys
dest = sys.argv[1]
rows = []
for name in sorted(os.listdir(dest)):
    if name.endswith(".part") or name == "raw":
        continue
    path = os.path.join(dest, name)
    if not os.path.isfile(path):
        continue
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    rows.append({"file": name, "bytes": os.path.getsize(path), "sha256": h.hexdigest()})
out = os.path.join(dest, "DOWNLOAD_MANIFEST.json")
with open(out, "w") as f:
    json.dump({"accession": "GSE148071", "files": rows}, f, indent=2)
print("manifest", out)
PY

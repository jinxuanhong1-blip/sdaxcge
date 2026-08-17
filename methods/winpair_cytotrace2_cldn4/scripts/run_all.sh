#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MOD="$ROOT/winpair_cytotrace2_cldn4"
DATA="${WINPAIR_CT2_DATA:-/tmp/winpair_cytotrace2_cldn4}"
VENV="${WINPAIR_CT2_VENV:-/tmp/winpair_ct2_venv}"

if [[ ! -x "$VENV/bin/python" ]]; then
  python3 -m venv "$VENV"
  "$VENV/bin/pip" install -U pip
  "$VENV/bin/pip" install -r "$MOD/requirements.txt"
fi

"$VENV/bin/python" "$MOD/scripts/download.py" --out "$DATA"
"$VENV/bin/python" "$MOD/scripts/extract_malignant.py" --data "$DATA" --out "$DATA/malignant.h5ad"
"$VENV/bin/python" "$MOD/scripts/analyze.py" \
  --input "$DATA/malignant.h5ad" \
  --outdir "$MOD/results" \
  --finding "$MOD/FINDING.md" \
  --work "$DATA/work"

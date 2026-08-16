#!/usr/bin/env bash
# Reproduce the opus_microarray slice (discovery → screen → probes → extract → analyze → verify).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="${ROOT}/.venv-opus-microarray/bin/python"
if [[ ! -x "$PY" ]]; then
  PY="${PYTHON:-python3}"
fi
cd "$ROOT"
"$PY" scripts/opus_microarray/01_discover_geo.py
"$PY" scripts/opus_microarray/02_screen_series.py
"$PY" scripts/opus_microarray/03_platform_probes.py
"$PY" scripts/opus_microarray/04_extract_expression.py
"$PY" scripts/opus_microarray/05_analyze.py
"$PY" scripts/opus_microarray/06_verify.py
echo "[done] results under results/opus_microarray/"

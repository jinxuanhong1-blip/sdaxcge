#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${ROOT}/.venv/bin/python"
if [[ ! -x "$PY" ]]; then
  PY=python3
fi
"$PY" "${ROOT}/scripts/01_build_genesets.py"
"$PY" "${ROOT}/scripts/02_bulk_GSE239485.py"
"$PY" "${ROOT}/scripts/03_sc_GSE133604.py"

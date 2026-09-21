#!/usr/bin/env bash
# MEBOCOST (kaifuchenlab) on the system Python. No conda required.
# CINE is not installed: no metabolite-CCC package by that name was found.
set -euo pipefail
ROOT="${MEBOCOST_ROOT:-/tmp/src/MEBOCOST}"
if [[ ! -f "$ROOT/pyproject.toml" && ! -f "$ROOT/setup.py" ]]; then
  mkdir -p "$(dirname "$ROOT")"
  git clone --depth 1 https://github.com/kaifuchenlab/MEBOCOST.git "$ROOT"
fi
python3 -m pip install -r "$ROOT/requirements.txt"
python3 -m pip install "$ROOT"
python3 - << 'PY'
from mebocost import mebocost
print("mebocost import ok", mebocost.__file__)
PY
echo "MEBOCOST_ROOT=$ROOT"

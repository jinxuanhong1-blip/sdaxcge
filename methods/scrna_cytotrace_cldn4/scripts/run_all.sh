#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$ROOT/scripts/download.py"
python3 "$ROOT/scripts/extract.py"
python3 "$ROOT/scripts/analyze.py"

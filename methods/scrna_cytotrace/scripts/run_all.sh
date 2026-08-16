#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python3 "$ROOT/download.py"
python3 "$ROOT/extract_gse207422.py"
python3 "$ROOT/extract_gse241934.py"
python3 "$ROOT/analyze.py"

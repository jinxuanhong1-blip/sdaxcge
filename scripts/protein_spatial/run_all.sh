#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 download_open_files.py
python3 analyze_cptac.py
python3 analyze_pride.py
python3 analyze_spatial.py
python3 write_catalog.py
echo "done"

#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 01_search_geo_2023.py
python3 02_classify_leftovers.py
python3 03_probe_leftovers.py
python3 04_analyze_usable.py
python3 04b_analyze_gse248378.py
python3 05_build_catalog.py
echo "done: results/w200/GEO_2023/"

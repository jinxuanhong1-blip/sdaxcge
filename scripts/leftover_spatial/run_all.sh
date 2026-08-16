#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python3 scripts/leftover_spatial/write_catalog.py
python3 scripts/leftover_spatial/analyze_visium.py
python3 scripts/leftover_spatial/analyze_geomx.py
python3 scripts/leftover_spatial/analyze_cosmx.py
python3 scripts/leftover_spatial/make_figures.py

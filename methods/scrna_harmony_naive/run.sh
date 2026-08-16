#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 download.py
python3 extract.py --datasets GSE131907 GSE148071 GSE127465
python3 extract_gse253013.py
python3 assemble_gse253013.py
python3 integrate.py
echo "done -> /workspace/results/scrna_harmony_naive"

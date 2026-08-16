#!/usr/bin/env bash
# Full replication pipeline. Everything it writes lands under results/align_tismo/.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NULL_GENES="${NULL_GENES:-500}"
WORKERS="${WORKERS:-6}"

python3 "$HERE/01_download.py" --null-genes "$NULL_GENES" --workers "$WORKERS"
python3 "$HERE/02_replicate.py"       | tee "$HERE/../../results/align_tismo/replication_report.txt"
python3 "$HERE/03_null_calibration.py" | tee "$HERE/../../results/align_tismo/null_calibration_report.txt"
python3 "$HERE/04_figures.py"

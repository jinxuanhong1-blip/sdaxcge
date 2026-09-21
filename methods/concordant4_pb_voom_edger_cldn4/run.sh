#!/usr/bin/env bash
# Patient-level malignant pseudo-bulk: limma-voom + edgeR, then cohort meta.
set -euo pipefail
cd "$(dirname "$0")"
python3 prep_units.py
Rscript analyze.R
Rscript analyze_robust.R

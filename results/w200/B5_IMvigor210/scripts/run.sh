#!/usr/bin/env bash
# Re-run B5 IMvigor210 CLDN4 analysis from a local IMvigor210CoreBiologies tarball.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
PKG="${1:-/tmp/dl/IMvigor210CoreBiologies_1.0.0.tar.gz}"
WORKDIR="${TMPDIR:-/tmp}/imvigor210_extract"
mkdir -p "$WORKDIR"
tar -C "$WORKDIR" -xzf "$PKG" IMvigor210CoreBiologies/data/cds.RData
Rscript "$ROOT/results/w200/B5_IMvigor210/scripts/01_extract_cds.R" \
  "$WORKDIR/IMvigor210CoreBiologies/data/cds.RData" \
  "$ROOT/results/w200/B5_IMvigor210/data/sample_level.tsv"
Rscript "$ROOT/results/w200/B5_IMvigor210/scripts/02_analyze_cldn4.R" \
  "$ROOT/results/w200/B5_IMvigor210/data/sample_level.tsv" \
  "$ROOT/results/w200/B5_IMvigor210"
echo "Done. Headline:"
column -t -s $'\t' "$ROOT/results/w200/B5_IMvigor210/tables/headline_OR_n_p.tsv" | head

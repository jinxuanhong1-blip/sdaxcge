#!/usr/bin/env bash
# Drive the official CIBERSORTx Docker image (Newman et al. 2019).
#
# This repository reimplements classic CIBERSORT relative-mode SVR. Absolute
# mode, B-mode / S-mode batch correction, and single-cell signature inference
# are only available from the authors. If you report "CIBERSORTx" numbers in a
# paper they must come from this container (or the web server), not from
# bulkimmune.cibersort.
#
# Prerequisites
#   1. Register at https://cibersortx.stanford.edu and accept the academic licence.
#   2. Request a Docker token (username + token) from the same portal.
#   3. docker login --username <portal-username> docker.synapse.org
#      (the image is distributed via Synapse; follow the portal's current
#      pull instructions -- the registry has moved more than once).
#
# Usage
#   export CSX_TOKEN=...
#   bash scripts/R/run_cibersortx.sh \\
#       --mixture data/GSE126044.counts.tsv.gz \\
#       --sig resources/LM22.txt \\
#       --outdir results/GSE126044/cibersortx
set -euo pipefail

MIXTURE=""
SIG=""
OUTDIR="cibersortx_out"
PERM=100
ABS="FALSE"
BMODE="FALSE"
IMAGE="${CIBERSORTX_IMAGE:-cibersortx/fractions}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mixture) MIXTURE="$2"; shift 2 ;;
    --sig) SIG="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --perm) PERM="$2"; shift 2 ;;
    --absolute) ABS="TRUE"; shift ;;
    --bmode) BMODE="TRUE"; shift ;;
    --image) IMAGE="$2"; shift 2 ;;
    *) echo "unknown arg $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$MIXTURE" || -z "$SIG" ]]; then
  echo "usage: $0 --mixture MIX.tsv --sig LM22.txt [--outdir DIR]" >&2
  exit 2
fi

mkdir -p "$OUTDIR"
# CIBERSORTx wants uncompressed tab-separated files named mixture / signature
python3 - <<'PY' "$MIXTURE" "$OUTDIR/mixture.txt"
import sys, pandas as pd
src, dest = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, sep="\t", index_col=0, compression="infer")
df.to_csv(dest, sep="\t")
PY
cp "$SIG" "$OUTDIR/signature.txt"

echo "Running $IMAGE (absolute=$ABS bmode=$BMODE perm=$PERM)"
# The flag names below match the public CIBERSORTx fractions container circa 2024.
# Re-check the portal if a pull fails -- Stanford has renamed flags before.
docker run --rm \
  -v "$(cd "$OUTDIR" && pwd)":/src/data \
  "$IMAGE" \
  --username "${CIBERSORTX_USERNAME:?set CIBERSORTX_USERNAME}" \
  --token "${CIBERSORTX_TOKEN:?set CIBERSORTX_TOKEN}" \
  --mixture mixture.txt \
  --sigmatrix signature.txt \
  --perm "$PERM" \
  --absolute "$ABS" \
  --rmbatchBmode "$BMODE" \
  --outdir /src/data

echo "CIBERSORTx finished. Outputs in $OUTDIR"
echo "Cite: Newman et al. Nat Biotechnol 2019;37:773-782."

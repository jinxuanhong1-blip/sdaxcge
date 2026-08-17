#!/usr/bin/env bash
# Install scanpy/PAGA and real Slingshot (R/Bioconductor) when possible.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python3 -m pip install --user -r "$ROOT/requirements.txt"
if ! command -v Rscript >/dev/null 2>&1; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -y
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
fi
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("slingshot", quietly=TRUE)) BiocManager::install("slingshot", ask=FALSE, update=FALSE)'
Rscript -e 'cat("slingshot", as.character(packageVersion("slingshot")), "\n")'

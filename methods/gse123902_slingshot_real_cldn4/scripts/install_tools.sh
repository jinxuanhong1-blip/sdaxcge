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
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p "$R_LIBS_USER"
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(Sys.getenv("R_LIBS_USER")); if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'
Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER")); if (!requireNamespace("slingshot", quietly=TRUE)) BiocManager::install(c("slingshot","DelayedMatrixStats","sparseMatrixStats"), ask=FALSE, update=FALSE, lib=Sys.getenv("R_LIBS_USER"))'
Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER")); cat("slingshot", as.character(packageVersion("slingshot")), "\n")'

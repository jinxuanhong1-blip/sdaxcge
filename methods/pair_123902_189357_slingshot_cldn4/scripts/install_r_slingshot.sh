#!/usr/bin/env bash
# Install R + Bioconductor slingshot (Street et al. 2018). Required — not optional.
set -euo pipefail

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Installing r-base via apt..."
  sudo apt-get update -y
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev
fi

echo "Rscript: $(command -v Rscript)"
Rscript --version

# System site-library is often not writable in this environment.
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p "$R_LIBS_USER"

Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(Sys.getenv("R_LIBS_USER")); if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'
Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER")); if (!requireNamespace("slingshot", quietly=TRUE) || !requireNamespace("DelayedMatrixStats", quietly=TRUE)) BiocManager::install(c("slingshot","SingleCellExperiment","DelayedMatrixStats"), ask=FALSE, update=FALSE, lib=Sys.getenv("R_LIBS_USER"))'
Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER")); cat("slingshot", as.character(packageVersion("slingshot")), "\n")'
echo "OK slingshot installed"

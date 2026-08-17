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

Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("slingshot", quietly=TRUE)) BiocManager::install(c("slingshot","SingleCellExperiment"), ask=FALSE, update=FALSE)'
Rscript -e 'cat("slingshot", as.character(packageVersion("slingshot")), "\n")'
echo "OK slingshot installed"

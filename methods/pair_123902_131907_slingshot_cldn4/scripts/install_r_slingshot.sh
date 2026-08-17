#!/usr/bin/env bash
# Install R + Bioconductor slingshot. Do not stop at R-missing.
set -euo pipefail

if ! command -v Rscript >/dev/null 2>&1; then
  sudo apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev \
    libgit2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev gfortran \
    liblapack-dev libblas-dev cmake
fi

LIB="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p "$LIB"
export R_LIBS_USER="$LIB"

Rscript -e "
Sys.setenv(R_LIBS_USER='$LIB')
.libPaths(c('$LIB', .libPaths()))
options(Ncpus = max(1L, parallel::detectCores()), timeout = 600)
options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/latest'))
if (!requireNamespace('BiocManager', quietly=TRUE)) {
  install.packages('BiocManager', lib='$LIB')
}
if (!requireNamespace('slingshot', quietly=TRUE)) {
  BiocManager::install(c('slingshot', 'SingleCellExperiment', 'TrajectoryUtils'),
                       ask=FALSE, update=FALSE, Ncpus=4, lib='$LIB')
}
cat('slingshot ', as.character(packageVersion('slingshot')), '\n')
"

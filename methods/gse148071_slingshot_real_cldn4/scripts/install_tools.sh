#!/usr/bin/env bash
# Install Python scanpy stack + R Bioconductor slingshot (REAL Slingshot).
set -euo pipefail

python3 -m pip install --user -q \
  'numpy>=1.26' 'pandas>=2.2' 'scipy>=1.11' 'matplotlib>=3.8' \
  'anndata>=0.10' 'scanpy>=1.10' 'igraph>=0.11' 'leidenalg>=0.10' \
  'seaborn>=0.13' 'scikit-misc>=0.3'

if ! command -v Rscript >/dev/null 2>&1; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    gfortran liblapack-dev libblas-dev
fi

export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p "$R_LIBS_USER"
Rscript -e '
.libPaths(Sys.getenv("R_LIBS_USER"))
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))
need <- c("slingshot", "SingleCellExperiment", "DelayedMatrixStats")
miss <- need[!vapply(need, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(miss))
  BiocManager::install(miss, ask=FALSE, update=FALSE, Ncpus=4, lib=Sys.getenv("R_LIBS_USER"))
library(slingshot)
cat("SLINGSHOT_OK", as.character(packageVersion("slingshot")), "\n")
'
echo "tools ready"

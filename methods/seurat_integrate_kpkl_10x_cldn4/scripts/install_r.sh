#!/usr/bin/env bash
# R + Seurat + Harmony. Binary packages via r2u on Ubuntu 24.04 when possible.
set -euo pipefail

if command -v Rscript >/dev/null 2>&1; then
  if Rscript -e 'quit(status=if (requireNamespace("Seurat", quietly=TRUE) && requireNamespace("harmony", quietly=TRUE) && requireNamespace("hdf5r", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1; then
    Rscript -e 'cat("R", R.version.string, "\n"); cat("Seurat", as.character(packageVersion("Seurat")), "\n"); cat("harmony", as.character(packageVersion("harmony")), "\n"); cat("hdf5r", as.character(packageVersion("hdf5r")), "\n")'
    exit 0
  fi
fi

export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends wget ca-certificates gnupg
if [[ ! -f /etc/apt/sources.list.d/cranapt.list ]]; then
  wget -q -O /tmp/add_cranapt_noble.sh \
    https://raw.githubusercontent.com/eddelbuettel/r2u/master/inst/scripts/add_cranapt_noble.sh
  sudo bash /tmp/add_cranapt_noble.sh
fi

sudo apt-get install -y --no-install-recommends \
  r-base r-cran-seurat r-cran-seuratobject r-cran-matrix r-cran-ggplot2 \
  r-cran-dplyr r-cran-harmony r-cran-hdf5r r-cran-jsonlite r-cran-patchwork \
  r-cran-uwot

# Harmony is required. If the r2u package name is missing, install from CRAN.
Rscript -e '
ok <- requireNamespace("harmony", quietly=TRUE)
if (!ok) {
  install.packages("harmony", repos="https://cloud.r-project.org")
}
ok2 <- requireNamespace("hdf5r", quietly=TRUE)
if (!ok2) {
  install.packages("hdf5r", repos="https://cloud.r-project.org")
}
cat("R", R.version.string, "\n")
cat("Seurat", as.character(packageVersion("Seurat")), "\n")
cat("harmony", as.character(packageVersion("harmony")), "\n")
cat("hdf5r", as.character(packageVersion("hdf5r")), "\n")
'

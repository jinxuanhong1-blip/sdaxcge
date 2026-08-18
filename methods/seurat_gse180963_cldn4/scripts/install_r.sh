#!/usr/bin/env bash
# R + Seurat (CreateSeuratObject). Binary packages via r2u on Ubuntu when possible.
set -euo pipefail
if command -v Rscript >/dev/null 2>&1; then
  if Rscript -e 'quit(status=if (requireNamespace("Seurat", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1; then
    Rscript -e 'cat("Seurat", as.character(packageVersion("Seurat")), "\n")'
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
  r-base r-cran-seurat r-cran-seuratobject r-cran-matrix r-cran-ggplot2 r-cran-dplyr
Rscript -e 'cat("R", R.version.string, "\n"); cat("Seurat", as.character(packageVersion("Seurat")), "\n")'

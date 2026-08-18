#!/usr/bin/env bash
# R + Seurat + Harmony. Binary packages via r2u on Ubuntu when possible.
# If Seurat cannot install, this script exits non-zero and the analysis MUST stop.
set -euo pipefail

have_seurat() {
  command -v Rscript >/dev/null 2>&1 && \
    Rscript -e 'quit(status=if (requireNamespace("Seurat", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1
}

have_harmony() {
  Rscript -e 'quit(status=if (requireNamespace("harmony", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1
}

if have_seurat && have_harmony; then
  Rscript -e 'cat("R", R.version.string, "\n"); cat("Seurat", as.character(packageVersion("Seurat")), "\n"); cat("harmony", as.character(packageVersion("harmony")), "\n")'
  exit 0
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
  r-cran-dplyr r-cran-hdf5r r-cran-jsonlite r-cran-patchwork

if ! have_seurat; then
  echo "STOP: Seurat could not be installed. No Python-only primary. Analysis halted." >&2
  exit 2
fi

# Harmony: r2u binary if present, else CRAN source via the already-installed Seurat stack.
if ! have_harmony; then
  if apt-cache show r-cran-harmony >/dev/null 2>&1; then
    sudo apt-get install -y --no-install-recommends r-cran-harmony || true
  fi
fi
if ! have_harmony; then
  sudo Rscript -e 'install.packages("harmony", repos="https://cloud.r-project.org")' || true
fi
if ! have_harmony; then
  echo "STOP: Seurat is present but Harmony could not be installed. Pair integration halted." >&2
  exit 3
fi

Rscript -e 'cat("R", R.version.string, "\n"); cat("Seurat", as.character(packageVersion("Seurat")), "\n"); cat("harmony", as.character(packageVersion("harmony")), "\n")'

#!/usr/bin/env bash
# Install Python scanpy/Harmony stack and real Slingshot (R/Bioconductor).
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"
python3 -m pip install --user -r "$(dirname "$0")/../requirements.txt"

if ! command -v Rscript >/dev/null 2>&1; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
fi

mkdir -p "$HOME/R/library"
export R_LIBS_USER="$HOME/R/library"
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))
BiocManager::install("slingshot", ask=FALSE, update=FALSE, lib=Sys.getenv("R_LIBS_USER"))
cat("slingshot", as.character(packageVersion("slingshot")), "\n")'
echo "tools installed"

#!/usr/bin/env bash
# Install REAL Slingshot (R/Bioconductor) and the Python scanpy stack.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${HOME}/.local/bin:${PATH}"
export R_LIBS_USER="${R_LIBS_USER:-/tmp/r-libs}"
mkdir -p "${R_LIBS_USER}"

python3 -m pip install --user -r "${ROOT}/requirements.txt"

if ! command -v Rscript >/dev/null 2>&1; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev gfortran liblapack-dev libblas-dev cmake
fi

Rscript --vanilla -e '
Sys.setenv(R_LIBS_USER="'"${R_LIBS_USER}"'")
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings=FALSE, recursive=TRUE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
options(Ncpus=max(1L, parallel::detectCores()), repos=c(CRAN="https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
if (!requireNamespace("slingshot", quietly=TRUE))
  BiocManager::install(c("slingshot", "SingleCellExperiment"), ask=FALSE, update=FALSE)
cat("slingshot", as.character(packageVersion("slingshot")), "\n")
'

echo "tools ready"
Rscript --version
python3 -c "import scanpy,harmonypy,leidenalg; print('scanpy', scanpy.__version__)"

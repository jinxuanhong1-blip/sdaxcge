#!/usr/bin/env bash
# Install REAL Slingshot (R/Bioconductor) + scanpy/PAGA. Idempotent.
set -euo pipefail
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p "$R_LIBS_USER"

if ! command -v Rscript >/dev/null 2>&1; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libhdf5-dev libfontconfig1-dev gfortran cmake
fi

python3 -m pip install --user -q -r "$(dirname "$0")/../requirements.txt"

Rscript - <<'RS'
dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org",
                   lib=Sys.getenv("R_LIBS_USER"))
}
if (!requireNamespace("slingshot", quietly=TRUE)) {
  BiocManager::install(
    c("slingshot", "SingleCellExperiment", "TrajectoryUtils",
      "DelayedMatrixStats", "sparseMatrixStats"),
    ask=FALSE, update=FALSE, Ncpus=4,
    lib=Sys.getenv("R_LIBS_USER"))
}
cat("slingshot", as.character(packageVersion("slingshot")), "\n")
RS
echo "tools ok"

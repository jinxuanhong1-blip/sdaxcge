#!/usr/bin/env Rscript
# Install Seurat + CellChat. If either fails, stop. Do not invent a Python primary.
options(Ncpus = max(1L, parallel::detectCores() - 1L))
options(repos = c(CRAN = "https://cloud.r-project.org"))
lib <- Sys.getenv("R_LIBS_USER", "/tmp/r_lib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)

cran_min <- c(
  "remotes", "BiocManager", "data.table", "Matrix", "ggplot2", "Rcpp",
  "RcppArmadillo", "RcppEigen", "irlba", "RcppAnnoy", "uwot", "Rtsne",
  "SeuratObject", "sctransform", "patchwork", "RColorBrewer", "igraph",
  "future", "future.apply", "jsonlite", "NMF", "ggalluvial", "svglite",
  "circlize", "FNN", "hdf5r"
)
need <- cran_min[!vapply(cran_min, ok, logical(1))]
if (length(need)) {
  message("CRAN: ", paste(need, collapse = ", "))
  install.packages(need, lib = lib, dependencies = c("Depends", "Imports"))
}
if (!ok("Seurat")) {
  message("Installing Seurat")
  install.packages("Seurat", lib = lib, dependencies = c("Depends", "Imports"))
}
if (!ok("BiocManager")) stop("BiocManager missing")
for (p in c("ComplexHeatmap", "BiocNeighbors")) {
  if (!ok(p)) {
    message("Bioconductor: ", p)
    BiocManager::install(p, lib = lib, ask = FALSE, update = FALSE)
  }
}
if (!ok("CellChat")) {
  message("Installing CellChat from jinworks/CellChat")
  remotes::install_github(
    "jinworks/CellChat",
    lib = lib,
    upgrade = "never",
    dependencies = TRUE
  )
}

failed <- c("Seurat", "CellChat")[!vapply(c("Seurat", "CellChat"), ok, logical(1))]
if (length(failed)) {
  stop("INSTALL FAILED: ", paste(failed, collapse = ", "))
}
cat("OK Seurat", as.character(packageVersion("Seurat")),
    "CellChat", as.character(packageVersion("CellChat")), "\n")

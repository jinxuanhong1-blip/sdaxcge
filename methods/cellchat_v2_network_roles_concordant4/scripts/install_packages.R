#!/usr/bin/env Rscript
# Seurat + CellChat v2. If either fails, stop. Do not invent a Python primary.
# Ubuntu 24.04 also needs libglpk (igraph) and libuv.
options(Ncpus = 4L)
ua <- sprintf(
  "R/%s R (%s)", getRversion(),
  paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
)
options(HTTPUserAgent = ua)
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))
lib <- Sys.getenv("R_LIBS_USER", "/tmp/r_lib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
Sys.setenv(USE_BUNDLED_LIBUV = "1")

ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)

if (!ok("fs")) install.packages("fs", lib = lib)
if (!ok("Seurat")) {
  install.packages("Seurat", lib = lib, dependencies = c("Depends", "Imports"))
}
if (!ok("Seurat")) stop("Seurat is not installed. Stop.")

if (!ok("BiocManager")) install.packages("BiocManager", lib = lib)
for (p in c("Biobase", "BiocGenerics", "ComplexHeatmap", "BiocNeighbors")) {
  if (!ok(p)) BiocManager::install(p, lib = lib, ask = FALSE, update = FALSE)
}
if (!ok("NMF")) install.packages("NMF", lib = lib, dependencies = c("Depends", "Imports"))
for (p in c(
  "ggalluvial", "ggpubr", "ggnetwork", "collapse", "RcppEigen",
  "igraph", "sna", "FNN", "circlize", "remotes", "data.table"
)) {
  if (!ok(p)) {
    install.packages(p, lib = lib, dependencies = c("Depends", "Imports"))
  }
}
if (!ok("RcppArmadillo")) {
  install.packages("RcppArmadillo", lib = lib, dependencies = c("Depends", "Imports"))
}
if (!ok("presto")) {
  remotes::install_github(
    "immunogenomics/presto", lib = lib, upgrade = "never",
    dependencies = FALSE
  )
}
if (!ok("CellChat")) {
  remotes::install_github(
    "jinworks/CellChat", lib = lib, upgrade = "never", dependencies = FALSE
  )
}
if (!ok("CellChat")) stop("CellChat is not installed. Stop.")
cat(
  "OK Seurat", as.character(packageVersion("Seurat")),
  "CellChat", as.character(packageVersion("CellChat")), "\n"
)

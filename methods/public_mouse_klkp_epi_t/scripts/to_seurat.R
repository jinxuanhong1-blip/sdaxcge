#!/usr/bin/env Rscript
# Build Seurat objects from the cleaned h5ad files written by preprocess.py.
# Requires Seurat and reticulate (anndata) on the library path.
# Expression objects are written under objects/ and are gitignored.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "methods/public_mouse_klkp_epi_t"
lib <- if (length(args) >= 2) args[[2]] else "/tmp/rlibs"
.libPaths(c(lib, .libPaths()))

status <- file.path(out, "objects", "SEURAT_STATUS.txt")
dir.create(dirname(status), recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("reticulate", quietly = TRUE)) {
  writeLines(
    c(
      "Seurat RDS was not written.",
      "Cleaned objects are the h5ad files in objects/ (log-normalized scored genes + cell metadata).",
      "Install Seurat and reticulate, then rerun this script."
    ),
    status
  )
  quit(save = "no", status = 0)
}

library(Seurat)
library(reticulate)
use_python(Sys.which("python3"), required = TRUE)
anndata <- import("anndata")
obj_dir <- file.path(out, "objects")
h5ads <- list.files(obj_dir, pattern = "\\.h5ad$", full.names = TRUE)
if (!length(h5ads)) stop("no h5ad files in ", obj_dir)

written <- character()
for (path in h5ads) {
  ad <- anndata$read_h5ad(path)
  # anndata X is cells x genes. Seurat wants genes x cells.
  mat <- t(ad$X)
  if (inherits(mat, "dgRMatrix") || inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "CsparseMatrix")
  }
  meta <- as.data.frame(ad$obs)
  colnames(mat) <- rownames(meta)
  rownames(mat) <- rownames(ad$var)
  obj <- CreateSeuratObject(counts = mat, meta.data = meta, assay = "RNA")
  dest <- sub("\\.h5ad$", ".rds", path)
  saveRDS(obj, dest)
  written <- c(written, dest)
  message("wrote ", dest)
}
writeLines(c("Seurat RDS written:", written), status)

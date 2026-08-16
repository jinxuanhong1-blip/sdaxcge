#!/usr/bin/env Rscript
# Optional numerical check: Python ssGSEA vs GSVA::gsva(method="ssgsea").
# Skipped automatically when GSVA is not installed.
if (!requireNamespace("GSVA", quietly = TRUE)) {
  cat("SKIP: GSVA not installed\n")
  quit(status = 0)
}

suppressPackageStartupMessages(library(GSVA))

set.seed(1)
n_genes <- 80
n_samples <- 6
mat <- matrix(rnorm(n_genes * n_samples), nrow = n_genes,
              dimnames = list(paste0("g", seq_len(n_genes)),
                              paste0("s", seq_len(n_samples))))
sets <- list(A = paste0("g", 1:12), B = paste0("g", 20:35))

# Write inputs for the Python side
tmpdir <- tempdir()
write.table(mat, file.path(tmpdir, "expr.tsv"), sep = "\t", quote = FALSE, col.names = NA)
writeLines(c(
  paste(c("A", "na", sets$A), collapse = "\t"),
  paste(c("B", "na", sets$B), collapse = "\t")
), file.path(tmpdir, "sets.gmt"))

if (packageVersion("GSVA") >= "1.50.0") {
  param <- GSVA::ssgseaParam(exprData = mat, geneSets = sets, normalize = FALSE, alpha = 0.25)
  ref <- GSVA::gsva(param)
} else {
  ref <- GSVA::gsva(mat, sets, method = "ssgsea", ssgsea.norm = FALSE, tau = 0.25)
}

write.table(ref, file.path(tmpdir, "gsva.tsv"), sep = "\t", quote = FALSE, col.names = NA)

py <- paste0("
import sys
from pathlib import Path
sys.path.insert(0, '", normalizePath(".."), "')
import pandas as pd
from bulkimmune.ssgsea import ssgsea
from bulkimmune.signatures import read_gmt
expr = pd.read_csv('", file.path(tmpdir, "expr.tsv"), "', sep='\\t', index_col=0)
sets = read_gmt('", file.path(tmpdir, "sets.gmt"), "')
out = ssgsea(expr, sets, alpha=0.25, normalize='none', tie_method='average_int')
out.to_csv('", file.path(tmpdir, "py.tsv"), "', sep='\\t')
")
status <- system2("python3", c("-c", py))
if (status != 0) stop("python ssGSEA failed")

py_scores <- as.matrix(read.delim(file.path(tmpdir, "py.tsv"), row.names = 1))
ref <- as.matrix(ref)
# align
common_sets <- intersect(rownames(ref), rownames(py_scores))
common_samples <- intersect(colnames(ref), colnames(py_scores))
delta <- max(abs(ref[common_sets, common_samples] - py_scores[common_sets, common_samples]))
cat(sprintf("max|python - GSVA| = %.6g on %d sets x %d samples\n",
            delta, length(common_sets), length(common_samples)))
if (delta > 1e-4) {
  stop("ssGSEA disagreement with GSVA exceeds 1e-4")
}
cat("OK\n")

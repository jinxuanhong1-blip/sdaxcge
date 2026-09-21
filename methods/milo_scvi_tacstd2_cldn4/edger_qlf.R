# edgeR quasi-likelihood F tests, the same GLM miloR::testNhoods calls
# (glmQLFit + glmQLFTest). Graph construction stays in Python.
#
# Usage:
#   Rscript edger_qlf.R <models.json> <out_dir>
#
# models.json is a list of objects:
#   name, counts, meta, formula, coef, norm (TMM|logMS), lib_mode
#   (colsum|column), lib_column, sample_column, keep_column (optional)

.libPaths(c(path.expand("~/R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: edger_qlf.R models.json out_dir")
spec <- jsonlite::fromJSON(args[[1]], simplifyVector = FALSE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# jsonlite may be missing; fall back to a tiny parser if needed.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required")
}

run_one <- function(job) {
  counts <- as.matrix(read.delim(job$counts, row.names = 1, check.names = FALSE))
  storage.mode(counts) <- "double"
  meta <- read.delim(job$meta, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta) <- meta[[job$sample_column]]
  if (!is.null(job$keep_column) && nzchar(job$keep_column)) {
    keep <- meta[[job$keep_column]] %in% c(TRUE, "TRUE", "true", 1, "1")
    meta <- meta[keep, , drop = FALSE]
  }
  missing <- setdiff(rownames(meta), colnames(counts))
  if (length(missing)) stop(job$name, " missing samples in counts: ", paste(head(missing, 5), collapse = ","))
  counts <- counts[, rownames(meta), drop = FALSE]
  if (job$lib_mode == "colsum") {
    lib <- colSums(counts)
  } else if (job$lib_mode == "column") {
    lib <- as.numeric(meta[[job$lib_column]])
  } else {
    stop("bad lib_mode")
  }
  names(lib) <- rownames(meta)
  ok <- is.finite(lib) & lib > 0
  if (!all(ok)) {
    meta <- meta[ok, , drop = FALSE]
    counts <- counts[, ok, drop = FALSE]
    lib <- lib[ok]
  }
  # Drop all-zero neighbourhoods in this sample set.
  keep_rows <- rowSums(counts) > 0
  n_drop <- sum(!keep_rows)
  counts <- counts[keep_rows, , drop = FALSE]
  design <- model.matrix(as.formula(job$formula), data = meta)
  if (nrow(design) != ncol(counts)) {
    stop(job$name, " design rows ", nrow(design), " != samples ", ncol(counts))
  }
  rk <- qr(design)$rank
  info <- data.frame(
    model = job$name,
    n_samples = ncol(counts),
    n_nhoods_in = nrow(counts),
    n_zero_rows_dropped = n_drop,
    n_coef = ncol(design),
    rank = rk,
    singular = rk < ncol(design),
    coef = job$coef,
    formula = job$formula,
    norm = job$norm,
    lib_mode = job$lib_mode,
    stringsAsFactors = FALSE
  )
  write.table(info, file.path(out_dir, paste0(job$name, ".design.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  if (isTRUE(info$singular)) {
    message(job$name, " SINGULAR rank ", rk, "/", ncol(design))
    return(invisible(NULL))
  }
  if (!(job$coef %in% colnames(design))) {
    stop(job$name, " coef ", job$coef, " not in ", paste(colnames(design), collapse = ","))
  }
  dge <- DGEList(counts = counts, lib.size = lib)
  if (job$norm == "TMM") {
    dge <- calcNormFactors(dge, method = "TMM")
  } else if (job$norm != "logMS") {
    stop("bad norm")
  }
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = job$coef)
  tab <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none")$table)
  tab$nhood_id <- rownames(counts)
  tab$model <- job$name
  write.table(tab, file.path(out_dir, paste0(job$name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(job$name, " nhoods ", nrow(tab), " minP ", signif(min(tab$PValue), 3))
}

for (job in spec) run_one(job)

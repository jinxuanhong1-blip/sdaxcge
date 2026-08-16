# Per-dataset / per-tumor-type log2 + unit-variance scaling for EcoTyper recovery.
# Independent of the EcoTyper executable so the step is auditable.
# Spec: playbook.md §4.4 and bayesprism_ecotyper_vs_tacstd2_cldn4.md §B5.1
#
# Input:  genes x samples, linear TPM/FPKM or already-logged values (declared in config).
# Output: the same shape, log2 (if needed), then scale each gene to unit variance
#         *within* each level of the grouping column (tumor type or dataset).

suppressPackageStartupMessages({
  library(yaml)
})

scale_for_ecotyper <- function(expr, group, already_log2 = FALSE, min_sd = 1e-8) {
  stopifnot(is.matrix(expr), length(group) == ncol(expr), !anyNA(group))
  if (!already_log2) {
    if (any(expr < 0, na.rm = TRUE)) {
      stop("negative values in a linear matrix; refuse to log2")
    }
    expr <- log2(expr + 1)
  }
  out <- expr
  for (g in unique(group)) {
    idx <- which(group == g)
    block <- expr[, idx, drop = FALSE]
    sds <- apply(block, 1L, sd, na.rm = TRUE)
    sds[!is.finite(sds) | sds < min_sd] <- NA_real_
    out[, idx] <- sweep(block, 1L, sds, "/")
  }
  out
}

# CLI wrapper when sourced as a script:
if (sys.nframe() == 0L && length(commandArgs(trailingOnly = TRUE))) {
  cfg <- yaml::read_yaml(commandArgs(trailingOnly = TRUE)[[1]])
  expr <- as.matrix(read.delim(cfg$input_expr, row.names = 1, check.names = FALSE))
  anno <- read.delim(cfg$annotation, check.names = FALSE)
  stopifnot(cfg$group_column %in% names(anno))
  sample_id <- if ("sample_id" %in% names(anno)) anno$sample_id else anno[[1]]
  stopifnot(all(colnames(expr) %in% sample_id))
  group <- anno[[cfg$group_column]][match(colnames(expr), sample_id)]
  scaled <- scale_for_ecotyper(expr, group, already_log2 = isTRUE(cfg$already_log2))
  dir.create(dirname(cfg$output_expr), recursive = TRUE, showWarnings = FALSE)
  write.table(scaled, cfg$output_expr, sep = "\t", quote = FALSE, col.names = NA)
}

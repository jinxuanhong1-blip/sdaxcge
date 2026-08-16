#!/usr/bin/env Rscript
# ssGSEA implementation for the hMENA-TGF-beta CAF9 and Pan-F-TBRS sets.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: score_ssgsea.R expression.tsv signatures.tsv output.tsv")
}
if (!requireNamespace("GSVA", quietly = TRUE)) stop("Bioconductor package 'GSVA' is required")

expr <- read.delim(args[[1]], row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
storage.mode(expr) <- "double"
rownames(expr) <- toupper(rownames(expr))
if (anyDuplicated(rownames(expr))) stop("duplicate gene symbols")

defs <- read.delim(args[[2]], stringsAsFactors = FALSE)
wanted <- c("hMENA_TGFB_CAF9_mean", "Pan_F_TBRS19_mean")
sets <- lapply(wanted, function(name) unique(defs$gene[defs$signature == name]))
names(sets) <- sub("_mean$", "_ssGSEA", wanted)
missing <- lapply(sets, function(genes) setdiff(genes, rownames(expr)))
if (any(lengths(missing))) {
  detail <- paste(
    names(missing)[lengths(missing) > 0],
    vapply(missing[lengths(missing) > 0], paste, collapse = ",", character(1)),
    sep = ": ", collapse = "; "
  )
  stop("exact ssGSEA requires complete sets; missing ", detail)
}

if (exists("ssgseaParam", envir = asNamespace("GSVA"), inherits = FALSE)) {
  param <- GSVA::ssgseaParam(expr, sets, minSize = 1, maxSize = Inf, normalize = TRUE)
  scored <- GSVA::gsva(param, verbose = FALSE)
} else {
  scored <- GSVA::gsva(
    expr, sets, method = "ssgsea", kcdf = "Gaussian",
    abs.ranking = FALSE, ssgsea.norm = TRUE, verbose = FALSE
  )
}

result <- data.frame(sample = colnames(expr), t(scored), check.names = FALSE)
if ("TGFB1" %in% rownames(expr)) {
  cutoff <- as.numeric(stats::quantile(expr["TGFB1", ], 0.10, names = FALSE))
  result$hMENA_TGFB1_above_p10 <- as.logical(expr["TGFB1", ] > cutoff)
} else {
  result$hMENA_TGFB1_above_p10 <- NA
  warning("TGFB1 absent: hMENA study eligibility gate cannot be evaluated")
}
write.table(result, args[[3]], sep = "\t", quote = FALSE, row.names = FALSE)

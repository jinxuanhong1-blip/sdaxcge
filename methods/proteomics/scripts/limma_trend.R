#!/usr/bin/env Rscript
# Optional Bioconductor limma / limma-trend wrapper.
# Usage:
#   Rscript limma_trend.R --matrix expr.tsv --group group.tsv --out fit.tsv \
#       --method limma-trend --num Tumor --den NAT
#
# TMT log2-ratio: method=limma (trend optional). Do not voom.
# RNA counts: convert to log-CPM first, then limma-trend.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(matrix=NULL, group=NULL, out=NULL, method="limma",
            sample_col="sample_id", group_col="group", num="Tumor", den="NAT",
            counts=FALSE)

i <- 1
while (i <= length(args)) {
  key <- args[[i]]
  if (key == "--counts") { opt$counts <- TRUE; i <- i + 1; next }
  val <- args[[i + 1]]
  if (key == "--matrix") opt$matrix <- val
  if (key == "--group") opt$group <- val
  if (key == "--out") opt$out <- val
  if (key == "--method") opt$method <- val
  if (key == "--sample-col") opt$sample_col <- val
  if (key == "--group-col") opt$group_col <- val
  if (key == "--num") opt$num <- val
  if (key == "--den") opt$den <- val
  i <- i + 2
}

if (is.null(opt$matrix) || is.null(opt$group) || is.null(opt$out)) {
  stop("required: --matrix --group --out")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Bioconductor limma is not installed. Use the Python limma_trend.py instead.")
}

mat <- as.matrix(read.delim(opt$matrix, row.names = 1, check.names = FALSE))
grp <- read.delim(opt$group, check.names = FALSE)
norm_id <- function(x) {
  x <- gsub("[.]", "-", as.character(x))
  x
}
colnames(mat) <- norm_id(colnames(mat))
grp[[opt$sample_col]] <- norm_id(grp[[opt$sample_col]])
common <- intersect(colnames(mat), grp[[opt$sample_col]])
mat <- mat[, common, drop = FALSE]
g <- grp[[opt$group_col]][match(common, grp[[opt$sample_col]])]
g <- factor(g, levels = c(opt$den, opt$num))
keep <- !is.na(g)
mat <- mat[, keep, drop = FALSE]
g <- droplevels(g[keep])

if (isTRUE(opt$counts)) {
  lib <- colSums(mat, na.rm = TRUE)
  mat <- log2(t(t(mat) / lib * 1e6) + 2)
}

# median-impute remaining NAs (document in playbook)
for (i in seq_len(nrow(mat))) {
  row <- mat[i, ]
  row[is.na(row)] <- median(row, na.rm = TRUE)
  mat[i, ] <- row
}

design <- model.matrix(~ 0 + g)
colnames(design) <- levels(g)
fit <- limma::lmFit(mat, design)
contrast <- limma::makeContrasts(contrasts = paste(opt$num, opt$den, sep = "-"), levels = design)
fit2 <- limma::contrasts.fit(fit, contrast)
trend <- identical(opt$method, "limma-trend")
fit2 <- limma::eBayes(fit2, trend = trend, robust = TRUE)
tab <- limma::topTable(fit2, number = Inf, sort.by = "P")
tab$gene <- rownames(tab)
tab$method <- opt$method
write.table(tab, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", opt$out)

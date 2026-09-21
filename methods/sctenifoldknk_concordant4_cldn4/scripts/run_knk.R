#!/usr/bin/env Rscript
# Official scTenifoldKnk virtual knockout. qc=FALSE because cells and genes
# were already filtered in prepare_grn.py. Network defaults are the package defaults
# except nCores=1 (one process per cohort) and nc_nCells capped by the matrix.

.libPaths(c("/home/ubuntu/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldKnk)
  library(RhpcBLASctl)
})
blas_set_num_threads(1)
omp_set_num_threads(1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: run_knk.R <knk_input_dir> <outdir> <gKO>")
}
indir <- args[[1]]
outdir <- args[[2]]
gko <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- readLines(file.path(indir, "genes.tsv"))
genes <- genes[nzchar(genes)]
m <- readMM(file.path(indir, "counts.mtx"))
m <- as(m, "dgCMatrix")
if (nrow(m) != length(genes)) stop("gene file length ", length(genes), " != mtx rows ", nrow(m))
rownames(m) <- genes
rs <- Matrix::rowSums(m)
if (any(rs <= 0)) {
  message("dropping zero-sum genes ", sum(rs <= 0))
  m <- m[rs > 0, , drop = FALSE]
}
cs <- Matrix::colSums(m)
if (any(cs <= 0)) m <- m[, cs > 0, drop = FALSE]
if (!gko %in% rownames(m)) stop(gko, " not in the count matrix")
if (ncol(m) < 200) stop("too few cells: ", ncol(m))

n_cells <- min(500L, ncol(m))
message("scTenifoldKnk ", gko, " genes ", nrow(m), " cells ", ncol(m), " nc_nCells ", n_cells)
t0 <- proc.time()
set.seed(1)
res <- scTenifoldKnk(
  countMatrix = m,
  gKO = gko,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = n_cells,
  nc_nComp = 3,
  nc_q = 0.9,
  td_K = 3,
  ma_nDim = 2,
  nCores = 1
)
elapsed <- unname((proc.time() - t0)[["elapsed"]])
dr <- res$diffRegulation
dr$Z <- as.numeric(dr$Z)
dr$distance <- as.numeric(dr$distance)
dr$FC <- as.numeric(dr$FC)
dr$p.value <- as.numeric(dr$p.value)
dr$p.adj <- as.numeric(dr$p.adj)
write.table(dr, file.path(outdir, "diffregulation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

W <- as.matrix(res$tensorNetworks$WT)
if (!gko %in% rownames(W)) stop("KO gene missing from WT network rownames")
edges <- data.frame(gene = colnames(W), weight = as.numeric(W[gko, ]), stringsAsFactors = FALSE)
write.table(edges, file.path(outdir, "outgoing_edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

meta <- data.frame(
  gKO = gko,
  n_genes = nrow(m),
  n_cells = ncol(m),
  nc_nNet = 10,
  nc_nCells = n_cells,
  nc_nComp = 3,
  nc_q = 0.9,
  td_K = 3,
  ma_nDim = 2,
  qc = FALSE,
  n_fdr05 = sum(dr$p.adj < 0.05, na.rm = TRUE),
  n_fdr05_excluding_ko = sum(dr$p.adj < 0.05 & dr$gene != gko, na.rm = TRUE),
  elapsed_sec = elapsed,
  scTenifoldKnk = as.character(packageVersion("scTenifoldKnk")),
  scTenifoldNet = as.character(packageVersion("scTenifoldNet")),
  stringsAsFactors = FALSE
)
write.table(meta, file.path(outdir, "run_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("done ", gko, " FDR<0.05 ", meta$n_fdr05, " elapsed ", round(elapsed, 1), "s")

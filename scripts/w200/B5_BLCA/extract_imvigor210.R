#!/usr/bin/env Rscript
# Extract IMvigor210 (Mariathasan 2018, atezolizumab, metastatic urothelial carcinoma)
# counts + phenotype from the archived IMvigor210CoreBiologies `cds` object into
# plain-text tables for downstream analysis.
#
# Output: data/interim/imvigor210_pheno.csv
#         data/interim/imvigor210_gene_tpm.csv   (subset of genes of interest)
#         data/interim/imvigor210_qc.json

suppressWarnings(suppressMessages({
  e <- new.env()
  loaded <- load("/workspace/data/raw/IMvigor210/cds.RData", envir = e)
}))
cds <- get(loaded[1], envir = e)

counts <- get("counts", envir = attr(cds, "assayData"))
pheno <- attr(attr(cds, "phenoData"), "data")
fdat <- attr(attr(cds, "featureData"), "data")

stopifnot(nrow(counts) == nrow(fdat), ncol(counts) == nrow(pheno))

# ---- TPM: length-normalise then scale to 1e6 per sample -----------------------
len_kb <- fdat$length / 1000
stopifnot(all(is.finite(len_kb)), all(len_kb > 0))
rpk <- counts / len_kb
tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6

# ---- DESeq-style size-factor normalised counts (sensitivity) -----------------
sf <- pheno$sizeFactor
norm_counts <- sweep(counts, 2, sf, "/")

# ---- genes of interest -------------------------------------------------------
genes <- list(
  primary          = c("CLDN4"),
  claudin_tj       = c("CLDN1", "CLDN3", "CLDN7", "CLDN18", "CLDN2", "TJP1", "OCLN", "CDH1", "EPCAM"),
  project_context  = c("TACSTD2"),
  # Positive control: CD8 T-effector / IFN-gamma programme, a signature reported to
  # associate with atezolizumab response in this very cohort (Mariathasan 2018).
  cd8_teff         = c("CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21", "PRF1"),
  checkpoint       = c("CD274", "PDCD1", "CTLA4", "LAG3"),
  tgfb             = c("TGFB1", "TGFB2", "TGFB3"),
  # Negative control: stably expressed housekeeping genes, expected to be unrelated
  # to response.
  housekeeping     = c("ACTB", "GAPDH", "TBP", "RPL13A", "PGK1")
)
want <- unique(unlist(genes, use.names = FALSE))

sym <- as.character(fdat$Symbol)
qc <- list()
rows <- list()
for (g in want) {
  idx <- which(sym == g)
  if (length(idx) == 0) {
    qc[[g]] <- "absent"
    next
  }
  qc[[g]] <- paste0("n_probes=", length(idx))
  # If a symbol maps to several rows, sum counts (same-gene fragments) before TPM
  # would be ideal; here TPM rows are already length-normalised, so sum TPM.
  tv <- if (length(idx) == 1) tpm[idx, ] else colSums(tpm[idx, , drop = FALSE])
  nv <- if (length(idx) == 1) norm_counts[idx, ] else colSums(norm_counts[idx, , drop = FALSE])
  rows[[paste0(g, "__tpm")]] <- tv
  rows[[paste0(g, "__normcount")]] <- nv
}

mat <- as.data.frame(rows, check.names = FALSE)
mat <- cbind(sample_id = colnames(counts), mat)

pheno_out <- cbind(sample_id = rownames(pheno), pheno)

dir.create("/workspace/data/interim", recursive = TRUE, showWarnings = FALSE)
write.csv(mat, "/workspace/data/interim/imvigor210_gene_tpm.csv", row.names = FALSE)
write.csv(pheno_out, "/workspace/data/interim/imvigor210_pheno.csv", row.names = FALSE)

# ---- QC / provenance summary -------------------------------------------------
bcor <- table(pheno[["Best Confirmed Overall Response"]], useNA = "ifany")
bin <- table(pheno[["binaryResponse"]], useNA = "ifany")
cat("n_samples:", ncol(counts), "\n")
cat("n_genes:", nrow(counts), "\n")
cat("\nBest Confirmed Overall Response:\n"); print(bcor)
cat("\nbinaryResponse:\n"); print(bin)
cat("\nTissue:\n"); print(table(pheno$Tissue, useNA = "ifany"))
cat("\nImmune phenotype:\n"); print(table(pheno[["Immune phenotype"]], useNA = "ifany"))
cat("\nCLDN4 probe count:", qc[["CLDN4"]], "\n")
cat("\nCLDN4 TPM summary:\n"); print(summary(mat[["CLDN4__tpm"]]))
cat("\ncolSums(tpm) check (should be 1e6):\n"); print(summary(colSums(tpm)))
cat("\nGenes absent from featureData:",
    paste(names(qc)[unlist(qc) == "absent"], collapse = ", "), "\n")

qc_json <- list(
  n_samples = ncol(counts),
  n_genes = nrow(counts),
  best_confirmed_overall_response = as.list(setNames(as.integer(bcor), names(bcor))),
  binary_response = as.list(setNames(as.integer(bin), names(bin))),
  gene_probe_counts = qc,
  cldn4_tpm_median = median(mat[["CLDN4__tpm"]]),
  tpm_colsum_min = min(colSums(tpm)),
  tpm_colsum_max = max(colSums(tpm))
)
writeLines(jsonlite_like <- {
  # minimal hand-rolled JSON to avoid an extra dependency
  esc <- function(x) gsub('"', '\\\\"', as.character(x))
  fmt <- function(v) {
    if (is.null(v)) return("null")
    if (is.list(v)) {
      inner <- paste0('"', esc(names(v)), '": ', vapply(v, fmt, ""), collapse = ", ")
      return(paste0("{", inner, "}"))
    }
    if (is.numeric(v) && length(v) == 1) return(as.character(v))
    paste0('"', esc(v), '"')
  }
  fmt(qc_json)
}, "/workspace/data/interim/imvigor210_qc.json")
cat("\nwrote data/interim/imvigor210_{pheno,gene_tpm}.csv and imvigor210_qc.json\n")

#!/usr/bin/env Rscript
# Extract IMvigor210 sample-level CLDN4 + clinical fields from the official
# IMvigor210CoreBiologies 1.0.0 CountDataSet (cds.RData).
#
# Does not redistribute the full transcriptome. Writes one row per RNA sample.

suppressPackageStartupMessages({
  .libPaths(c("/tmp/Rlib", .libPaths()))
  library(Biobase)
  library(DESeq)
})

args <- commandArgs(trailingOnly = TRUE)
cds_path <- if (length(args) >= 1) args[[1]] else "/tmp/dl/IMvigor210CoreBiologies/data/cds.RData"
out_path <- if (length(args) >= 2) args[[2]] else "results/w200/B5_IMvigor210/data/sample_level.tsv"

e <- new.env()
load(cds_path, envir = e)
if (!exists("cds", envir = e)) stop("cds object not found in ", cds_path)
cds <- e$cds

ad <- slot(cds, "assayData")
counts <- get("counts", envir = ad)
pd <- pData(slot(cds, "phenoData"))
fd <- pData(slot(cds, "featureData"))

if (!identical(colnames(counts), rownames(pd))) {
  stop("Sample names in counts and phenoData do not match")
}
if (!identical(rownames(counts), rownames(fd))) {
  stop("Feature names in counts and featureData do not match")
}

cldn4_idx <- which(fd$symbol == "CLDN4")
if (length(cldn4_idx) != 1) {
  stop("Expected exactly one CLDN4 row, found ", length(cldn4_idx))
}

raw <- as.numeric(counts[cldn4_idx, ])
sf <- as.numeric(pd$sizeFactor)
if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid sizeFactor")

# DESeq size-factor normalized count, then log2(x+1). Gene length is constant
# across samples so this ranking equals CPM ranking up to a monotone transform
# of library size vs sizeFactor.
norm <- raw / sf
log2_norm <- log2(norm + 1)
libsize <- colSums(counts)
cpm <- raw / libsize * 1e6
log2_cpm <- log2(cpm + 1)

out <- data.frame(
  sample_id = rownames(pd),
  ANONPT_ID = as.character(pd$ANONPT_ID),
  CLDN4_entrez = fd$entrez_id[cldn4_idx],
  CLDN4_symbol = fd$symbol[cldn4_idx],
  CLDN4_length = fd$length[cldn4_idx],
  CLDN4_raw_count = raw,
  sizeFactor = sf,
  library_size = as.numeric(libsize),
  CLDN4_deseq_norm = norm,
  CLDN4_log2_deseq = log2_norm,
  CLDN4_cpm = cpm,
  CLDN4_log2_cpm = log2_cpm,
  Best_Confirmed_Overall_Response = as.character(pd[["Best Confirmed Overall Response"]]),
  binaryResponse = as.character(pd$binaryResponse),
  Enrollment_IC = as.character(pd[["Enrollment IC"]]),
  IC_Level = as.character(pd[["IC Level"]]),
  TC_Level = as.character(pd[["TC Level"]]),
  Immune_phenotype = as.character(pd[["Immune phenotype"]]),
  FMOne_TMB_per_MB = as.numeric(pd[["FMOne mutation burden per MB"]]),
  Neoantigen_burden_per_MB = as.numeric(pd[["Neoantigen burden per MB"]]),
  Sex = as.character(pd$Sex),
  Race = as.character(pd$Race),
  Intravesical_BCG = as.character(pd[["Intravesical BCG administered"]]),
  Baseline_ECOG = as.numeric(as.character(pd[["Baseline ECOG Score"]])),
  Tobacco_Use_History = as.character(pd[["Tobacco Use History"]]),
  Met_Disease_Status = as.character(pd[["Met Disease Status"]]),
  Sample_age = as.character(pd[["Sample age"]]),
  Tissue = as.character(pd$Tissue),
  Received_platinum = as.character(pd[["Received platinum"]]),
  Sample_collected_pre_platinum = as.character(pd[["Sample collected pre-platinum"]]),
  os_months = as.numeric(pd$os),
  censOS = as.integer(as.character(pd$censOS)),
  Lund = as.character(pd$Lund),
  Lund2 = as.character(pd$Lund2),
  TCGA_Subtype = as.character(pd[["TCGA Subtype"]]),
  stringsAsFactors = FALSE
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote", nrow(out), "rows to", out_path, "\n")
cat("CLDN4 raw range:", paste(range(raw), collapse = "-"), "\n")
cat("CLDN4 log2 DESeq range:", paste(sprintf("%.3f", range(log2_norm)), collapse = "-"), "\n")
cat("Unique ANONPT_ID:", length(unique(out$ANONPT_ID)), "\n")
dup <- out$ANONPT_ID[duplicated(out$ANONPT_ID)]
if (length(dup)) {
  cat("Duplicate ANONPT_ID values:\n")
  print(out[out$ANONPT_ID %in% unique(dup), c("sample_id", "ANONPT_ID", "Tissue", "Best_Confirmed_Overall_Response", "CLDN4_log2_deseq", "os_months")])
}

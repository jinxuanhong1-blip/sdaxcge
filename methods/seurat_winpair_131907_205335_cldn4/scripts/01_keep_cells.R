#!/usr/bin/env Rscript
# Build capped malignant + T/NK barcode lists. Unit = sample (GSE131907)
# or patient (GSE205335). Locked PR #320 tables gate eligibility.
suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
geo <- if (length(args) >= 1) args[[1]] else "/tmp/winpair_geo"
root <- if (length(args) >= 2) args[[2]] else {
  this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(this)) {
    "/workspace/methods/seurat_winpair_131907_205335_cldn4"
  } else {
    dirname(dirname(normalizePath(this)))
  }
}

cap_mal <- 120L
cap_tnk <- 120L
seed <- 1L
set.seed(seed)

TUMOR_ORIGINS <- c("tLung", "tL/B", "mLN", "PE", "mBrain")
NORMAL_TISSUE <- c("Normal Lung", "Normal LN", "Normal Brain")

data_dir <- file.path(root, "data")
out_dir <- file.path(root, "extract")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s131 <- read.delim(file.path(data_dir, "GSE131907_samples.tsv"), check.names = FALSE)
p205 <- read.delim(file.path(data_dir, "GSE205335_patients.tsv"), check.names = FALSE)
gsm <- read.csv(file.path(data_dir, "GSE205335_gsm_sample_metadata.csv"), check.names = FALSE)

elig131 <- s131$sample[s131$n_malignant >= 20 & s131$origin %in% TUMOR_ORIGINS]
elig205 <- p205$patient[p205$n_malignant >= 20]

ann <- read.delim(
  file.path(geo, "GSE131907", "GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  check.names = FALSE
)
ann$Index <- as.character(ann$Index)
ann$compartment <- NA_character_
ann$compartment[ann$Cell_subtype == "Malignant cells" & ann$Sample %in% elig131] <- "malignant"
ann$compartment[ann$Cell_type %in% c("T lymphocytes", "NK cells") & ann$Sample %in% elig131] <- "TNK"
keep131 <- subset(ann, !is.na(compartment))
keep131$dataset <- "GSE131907"
keep131$unit_id <- as.character(keep131$Sample)
keep131$patient_id <- gsub(".*_(\\d+)$", "\\1", keep131$Sample)
keep131$barcode <- keep131$Index

cap_by_unit <- function(df, unit_col, cap, seed_i) {
  set.seed(seed_i)
  parts <- split(df, df[[unit_col]], drop = TRUE)
  out <- lapply(parts, function(sub) {
    if (nrow(sub) <= cap) return(sub)
    sub[sample.int(nrow(sub), cap), , drop = FALSE]
  })
  do.call(rbind, out)
}

mal131 <- cap_by_unit(subset(keep131, compartment == "malignant"), "unit_id", cap_mal, seed)
tnk131 <- cap_by_unit(subset(keep131, compartment == "TNK"), "unit_id", cap_tnk, seed + 1L)
sel131 <- rbind(mal131, tnk131)
rownames(sel131) <- NULL

cid <- read.delim(
  file.path(geo, "GSE205335_Lung_IO_CellIdentity.txt.gz"),
  check.names = FALSE
)
map <- unique(gsm[, c("orig.ident", "patient", "tissue")])
cid$orig.ident <- as.character(cid$orig.ident)
cid <- merge(cid, map, by = "orig.ident", all.x = TRUE)
cid <- subset(cid, patient %in% elig205 & !tissue %in% NORMAL_TISSUE)
cid$compartment <- NA_character_
cid$compartment[cid$lineage.sub == "Malignant cells"] <- "malignant"
cid$compartment[cid$lineage.total == "T/NK cells"] <- "TNK"
keep205 <- subset(cid, !is.na(compartment))
keep205$dataset <- "GSE205335"
keep205$unit_id <- as.character(keep205$patient)
keep205$patient_id <- as.character(keep205$patient)
keep205$barcode <- as.character(keep205$barcode)

mal205 <- cap_by_unit(subset(keep205, compartment == "malignant"), "unit_id", cap_mal, seed + 2L)
tnk205 <- cap_by_unit(subset(keep205, compartment == "TNK"), "unit_id", cap_tnk, seed + 3L)
sel205 <- rbind(mal205, tnk205)
rownames(sel205) <- NULL

cols131 <- c("barcode", "dataset", "unit_id", "patient_id", "compartment",
             "Sample", "Sample_Origin", "Cell_type", "Cell_subtype")
cols205 <- c("barcode", "dataset", "unit_id", "patient_id", "compartment",
             "orig.ident", "tissue", "lineage.total", "lineage.sub", "celltype")
write.table(sel131[, cols131], file.path(out_dir, "keep_GSE131907.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sel205[, cols205], file.path(out_dir, "keep_GSE205335.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

inv <- rbind(
  data.frame(
    dataset = "GSE131907",
    n_eligible_units = length(elig131),
    n_unique_patient_ids = length(unique(sel131$patient_id)),
    n_mal_cells_kept = sum(sel131$compartment == "malignant"),
    n_tnk_cells_kept = sum(sel131$compartment == "TNK"),
    n_mal_available = sum(keep131$compartment == "malignant"),
    n_tnk_available = sum(keep131$compartment == "TNK")
  ),
  data.frame(
    dataset = "GSE205335",
    n_eligible_units = length(elig205),
    n_unique_patient_ids = length(unique(sel205$patient_id)),
    n_mal_cells_kept = sum(sel205$compartment == "malignant"),
    n_tnk_cells_kept = sum(sel205$compartment == "TNK"),
    n_mal_available = sum(keep205$compartment == "malignant"),
    n_tnk_available = sum(keep205$compartment == "TNK")
  )
)
write.table(inv, file.path(out_dir, "keep_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote keep lists\n")
print(inv)
cat("GSE131907 unique Sample vs patient_id:",
    length(unique(sel131$unit_id)), length(unique(sel131$patient_id)), "\n")

#!/usr/bin/env Rscript
# Reproduce the GSE205335 malignant-cell TACSTD2/CLDN4 analysis.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else file.path(dirname(normalizePath(sys.frame(1)$ofile)), "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache <- tempfile("gse205335_")
dir.create(cache)
on.exit(unlink(cache, recursive = TRUE), add = TRUE)

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/"
matrix_gz <- file.path(cache, "GSE205335_Lung_IO_UMI_matrix.rds.gz")
identity_gz <- file.path(cache, "GSE205335_Lung_IO_CellIdentity.txt.gz")
download.file(paste0(base_url, basename(matrix_gz)), matrix_gz, mode = "wb", quiet = TRUE)
download.file(paste0(base_url, basename(identity_gz)), identity_gz, mode = "wb", quiet = TRUE)

# GEO wrapped an already gzip-compressed RDS in another gzip layer.
matrix_inner <- file.path(cache, "matrix_inner.rds.gz")
input <- gzfile(matrix_gz, "rb")
output <- file(matrix_inner, "wb")
repeat {
  block <- readBin(input, "raw", 1024 * 1024)
  if (!length(block)) break
  writeBin(block, output)
}
close(input)
close(output)

expression <- readRDS(matrix_inner)
annotation <- read.delim(identity_gz, check.names = FALSE)
annotation <- annotation[match(colnames(expression), annotation$barcode), ]
stopifnot(identical(colnames(expression), annotation$barcode))

# Clinical labels are the 14 "Core" rows in eLife Supplementary file 1.
clinical <- read.table(
  text = "
sample patient response recist
EBUS-06-3P P1006 Responder PR
EFFUSION-06-3P P1006 Responder PR
EBUS-27-5P P1027 Responder PR
EBUS-37-3P P1037 Responder PR
EBUS-90-5P P1090 Responder PR
EBUS-17-5P P1017 Non-responder SD
LM-17-5P P1017 Non-responder SD
PCNB-01-5P P4001 Non-responder SD
EBUS-30-5P P1030 Non-responder PD
EBUS-62-5P P1062 Non-responder PD
EBUS-76-3P P1076 Non-responder PD
NECK-05-3P P1076 Non-responder PD
EBUS-89-3P P1089 Non-responder PD
EBUS-119-3P P1119 Non-responder PD
", header = TRUE, stringsAsFactors = FALSE
)
stopifnot(all(clinical$sample %in% annotation$orig.ident))

keep <- (
  !is.na(annotation$lineage.sub) &
    annotation$lineage.sub == "Malignant cells" &
    annotation$orig.ident %in% clinical$sample
)
stopifnot(sum(keep) == 12975)
malignant <- expression[, keep, drop = FALSE]
malignant_annotation <- annotation[keep, ]
cell_library <- Matrix::colSums(malignant)

sample_rows <- list()
for (sample_id in clinical$sample) {
  indices <- which(malignant_annotation$orig.ident == sample_id)
  label <- clinical[clinical$sample == sample_id, ]
  target_values <- sapply(c("TACSTD2", "CLDN4"), function(gene) {
    counts <- as.numeric(malignant[gene, indices])
    c(umi = sum(counts), detection = mean(counts > 0))
  })
  total_umi <- sum(cell_library[indices])
  log2_cpm <- log2(target_values["umi", ] / total_umi * 1e6 + 1)
  both <- mean(
    as.numeric(malignant["TACSTD2", indices]) > 0 &
      as.numeric(malignant["CLDN4", indices]) > 0
  )
  sample_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    patient = label$patient,
    response = label$response,
    recist = label$recist,
    n_malignant_cells = length(indices),
    total_malignant_umi = total_umi,
    TACSTD2_umi = target_values["umi", "TACSTD2"],
    CLDN4_umi = target_values["umi", "CLDN4"],
    TACSTD2_detect = target_values["detection", "TACSTD2"],
    CLDN4_detect = target_values["detection", "CLDN4"],
    both_detect = both,
    TACSTD2_log2cpm = log2_cpm["TACSTD2"],
    CLDN4_log2cpm = log2_cpm["CLDN4"]
  )
}
samples <- do.call(rbind, sample_rows)
write.table(
  samples,
  file.path(out_dir, "GSE205335_core_sample_targets.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

patients <- aggregate(
  cbind(
    TACSTD2_log2cpm, CLDN4_log2cpm,
    TACSTD2_detect, CLDN4_detect, both_detect
  ) ~ patient + response,
  data = samples,
  FUN = mean
)
write.table(
  patients,
  file.path(out_dir, "GSE205335_core_patient_targets.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

correlation <- cor.test(
  patients$TACSTD2_log2cpm,
  patients$CLDN4_log2cpm,
  method = "spearman",
  exact = FALSE
)

permutation_result <- function(column) {
  values <- patients[[column]]
  responder <- patients$response == "Responder"
  observed <- mean(values[responder]) - mean(values[!responder])
  assignments <- combn(length(values), sum(responder))
  null <- apply(assignments, 2, function(index) {
    mean(values[index]) - mean(values[-index])
  })
  data.frame(
    metric = column,
    responder_n = sum(responder),
    nonresponder_n = sum(!responder),
    responder_median = median(values[responder]),
    nonresponder_median = median(values[!responder]),
    responder_minus_nonresponder_mean = observed,
    exact_two_sided_permutation_p = mean(abs(null) >= abs(observed) - 1e-12)
  )
}

tests <- do.call(
  rbind,
  lapply(
    c("TACSTD2_log2cpm", "CLDN4_log2cpm", "TACSTD2_detect", "CLDN4_detect"),
    permutation_result
  )
)
tests <- rbind(
  tests,
  data.frame(
    metric = "TACSTD2_CLDN4_patient_spearman",
    responder_n = nrow(patients),
    nonresponder_n = NA,
    responder_median = unname(correlation$estimate),
    nonresponder_median = NA,
    responder_minus_nonresponder_mean = NA,
    exact_two_sided_permutation_p = correlation$p.value
  )
)
write.table(
  tests,
  file.path(out_dir, "GSE205335_statistics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

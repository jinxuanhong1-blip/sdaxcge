# GSE205335: author malignant CLDN4 (and TACSTD2 check) per patient.
# Uses public RDS + CellIdentity + the already-curated orig.ident → patient/RECIST map.

args <- commandArgs(trailingOnly = TRUE)
rds_gz <- if (length(args) >= 1) args[[1]] else "/tmp/gse205335/GSE205335_Lung_IO_UMI_matrix.rds"
id_path <- if (length(args) >= 2) args[[2]] else "/tmp/gse205335/GSE205335_Lung_IO_CellIdentity.txt.gz"
map_path <- if (length(args) >= 3) args[[3]] else "/workspace/methods/scrna_meta_mpr/inputs/GSE205335_sample_map.csv"
out_dir <- if (length(args) >= 4) args[[4]] else "/workspace/methods/scrna_meta_mpr/inputs"

message("reading ", rds_gz)
mat <- readRDS(rds_gz)
message("dim ", paste(dim(mat), collapse = "x"), " class ", paste(class(mat), collapse = "/"))
stopifnot("CLDN4" %in% rownames(mat), "TACSTD2" %in% rownames(mat))

id <- read.delim(id_path, check.names = FALSE)
id$barcode <- as.character(id$barcode)
map <- read.csv(map_path, stringsAsFactors = FALSE)

lib <- Matrix::colSums(mat)
cldn_umi <- as.numeric(mat["CLDN4", ])
tac_umi <- as.numeric(mat["TACSTD2", ])
cldn <- log1p(cldn_umi / pmax(lib, 1) * 1e4)
tac <- log1p(tac_umi / pmax(lib, 1) * 1e4)

cells <- colnames(mat)
idx <- match(cells, id$barcode)
if (mean(!is.na(idx)) < 0.9) {
  # try stripping prefixes
  message("direct barcode match ", mean(!is.na(idx)), "; trying orig.ident alignment by order")
}
id2 <- id[idx, ]
pc <- data.frame(
  barcode = cells,
  orig.ident = as.character(id2$orig.ident),
  lineage.sub = as.character(id2$lineage.sub),
  lineage.total = as.character(id2$lineage.total),
  CLDN4_umi = cldn_umi,
  TACSTD2_umi = tac_umi,
  CLDN4_log1p_cp10k = cldn,
  TACSTD2_log1p_cp10k = tac,
  stringsAsFactors = FALSE
)
pc$compartment <- ifelse(pc$lineage.sub == "Malignant cells", "Malignant",
                  ifelse(pc$lineage.total == "T/NK cells", "T/NK", "other"))
pc <- merge(pc, unique(map[, c("orig.ident", "patient", "recist", "response")]),
            by = "orig.ident", all.x = TRUE)

agg <- function(gene_log, gene_umi) {
  do.call(rbind, lapply(split(pc, list(pc$patient, pc$compartment), drop = TRUE), function(d) {
    data.frame(
      patient = d$patient[1],
      recist = d$recist[1],
      response = d$response[1],
      compartment = d$compartment[1],
      n_cells = nrow(d),
      pct_pos = 100 * mean(d[[gene_umi]] > 0),
      mean_log1p_cp10k = mean(d[[gene_log]]),
      stringsAsFactors = FALSE
    )
  }))
}

cld_pt <- agg("CLDN4_log1p_cp10k", "CLDN4_umi")
tac_pt <- agg("TACSTD2_log1p_cp10k", "TACSTD2_umi")
write.csv(cld_pt, file.path(out_dir, "GSE205335_per_patient_cldn4.csv"), row.names = FALSE)
write.csv(tac_pt, file.path(out_dir, "GSE205335_per_patient_tacstd2_reextract.csv"), row.names = FALSE)
message("CLDN4 malignant n patients: ", sum(cld_pt$compartment == "Malignant"))
print(table(cld_pt$compartment, cld_pt$response, useNA = "ifany"))

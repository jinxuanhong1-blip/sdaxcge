#!/usr/bin/env Rscript
# Pull a shared gene panel for Q1/Q4 GSE205335 cells. Double-gzip RDS.

suppressPackageStartupMessages(library(Matrix))

geo <- "/tmp/geo_c4"
plain <- file.path(geo, "GSE205335_plain.rds")
if (!file.exists(plain)) {
  message("peeling double gzip")
  system(sprintf(
    "gzip -dc %s | gzip -dc > %s",
    file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds.gz"),
    plain
  ))
}

root <- "/workspace/methods/concordant4_augur_dirichlet_ilr"
locked <- read.delim(file.path(root, "results/tables/composition_counts.tsv"), stringsAsFactors = FALSE)
locked <- locked[locked$dataset == "GSE205335" & locked$cldn4_quartile %in% c("Q1", "Q4"), ]
gsm <- read.delim(file.path(root, "data/GSE205335_gsm_map.tsv"), stringsAsFactors = FALSE)
normal <- gsm$orig.ident[startsWith(gsm$tissue, "Normal")]
patient_of <- setNames(gsm$patient, gsm$orig.ident)

ident <- read.delim(file.path(geo, "GSE205335_Lung_IO_CellIdentity.txt.gz"), stringsAsFactors = FALSE)
ident <- ident[!(ident$orig.ident %in% normal) & ident$orig.ident %in% names(patient_of), ]
ident$patient <- unname(patient_of[ident$orig.ident])
ident <- ident[ident$patient %in% locked$unit_id, ]
q_of <- setNames(locked$cldn4_quartile, locked$unit_id)
ident$quartile <- unname(q_of[ident$patient])

lab <- rep(NA_character_, nrow(ident))
lab[ident$lineage.sub == "Malignant cells"] <- "malignant"
lab[ident$lineage.sub == "NK cells"] <- "NK"
lab[ident$lineage.total == "T/NK cells" & ident$lineage.sub != "NK cells"] <- "T"
lab[ident$lineage.total == "Myeloid cells"] <- "myeloid"
lab[ident$lineage.total == "B/Plasma cells"] <- "B"
ident$cell_type <- lab
ident <- ident[!is.na(ident$cell_type), ]

set.seed(1)
per_unit <- 25L
min_cells <- 8L
picked <- do.call(rbind, lapply(split(ident, paste(ident$patient, ident$cell_type)), function(d) {
  if (nrow(d) < min_cells) return(NULL)
  if (nrow(d) > per_unit) d <- d[sample.int(nrow(d), per_unit), , drop = FALSE]
  d
}))
message("picked cells ", nrow(picked))

panel <- readLines(file.path(root, "data/gene_panel.txt"))
panel <- toupper(trimws(sub("#.*", "", panel)))
panel <- unique(panel[nzchar(panel)])

message("reading RDS")
mat <- readRDS(plain)
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
rn <- toupper(rownames(mat))
keep_gene <- !duplicated(rn)
mat <- mat[keep_gene, , drop = FALSE]
rownames(mat) <- rn[keep_gene]
message("matrix ", nrow(mat), " x ", ncol(mat))

bc <- picked$barcode
common <- intersect(bc, colnames(mat))
if (length(common) < 100) {
  alt <- gsub("_", "-", colnames(mat))
  names(alt) <- colnames(mat)
  common_mat <- colnames(mat)[alt %in% bc]
  if (!length(common_mat)) stop("barcode mismatch: ", colnames(mat)[1], " vs ", bc[1])
  map <- setNames(colnames(mat), alt)
  picked$mat_bc <- unname(map[picked$barcode])
} else {
  picked$mat_bc <- picked$barcode
}
picked <- picked[!is.na(picked$mat_bc) & picked$mat_bc %in% colnames(mat), ]
message("matched ", nrow(picked))

sub <- mat[intersect(panel, rownames(mat)), picked$mat_bc, drop = FALSE]
lib <- Matrix::colSums(mat[, picked$mat_bc, drop = FALSE])
# dense panel in gene-panel order
expr <- matrix(0, nrow = length(panel), ncol = nrow(picked), dimnames = list(panel, picked$mat_bc))
have <- intersect(panel, rownames(sub))
expr[have, ] <- as.matrix(sub[have, , drop = FALSE])
out_dir <- file.path(geo, "augur_extract")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# cells x genes, same layout as the python npz
write.table(
  data.frame(
    dataset = "GSE205335",
    unit_id = picked$patient,
    quartile = picked$quartile,
    cell_type = picked$cell_type,
    lib_size = as.numeric(lib[picked$mat_bc]),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "gse205335_panel.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
# save as a plain csv.gz of the expression (cells x genes) to avoid an R->npy bridge
expr_t <- t(expr)
write.table(
  data.frame(expr_t, check.names = FALSE),
  file = gzfile(file.path(out_dir, "gse205335_expr.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
message("wrote GSE205335 panel")

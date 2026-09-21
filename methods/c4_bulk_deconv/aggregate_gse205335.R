# Pool GSE205335 UMI counts by author lineage. Sparse throughout.
args <- commandArgs(trailingOnly = TRUE)
mat_path <- args[[1]]
id_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(Matrix))

message("reading ", mat_path)
m <- readRDS(mat_path)
message("class: ", paste(class(m), collapse = " | "))
if (is.list(m) && !inherits(m, "dgCMatrix")) {
  message("list names: ", paste(names(m), collapse = ", "))
  stop("RDS is not a sparse matrix; inspect and extend the loader")
}
message("dim: ", paste(dim(m), collapse = " x "))
message("rownames: ", paste(head(rownames(m), 3), collapse = " | "))
message("colnames: ", paste(head(colnames(m), 3), collapse = " | "))

ids <- read.delim(id_path, stringsAsFactors = FALSE, check.names = FALSE)
message("id barcodes: ", paste(head(ids$barcode, 3), collapse = " | "))

class_of <- function(sub, total) {
  if (is.na(sub)) sub <- ""
  if (is.na(total)) total <- ""
  if (sub == "Malignant cells") return("malignant")
  if (sub == "CD8+ T cells") return("cd8")
  if (sub == "CD4+ T cells") return("cd4")
  if (sub == "NK cells") return("nk")
  if (sub == "B/Plasma cells" || total == "B/Plasma cells") return("b")
  if (sub == "Myeloid cells") return("myeloid")
  if (total == "Endothelial cells") return("endothelial")
  if (total == "Fibroblasts") return("fibroblast")
  return(NA_character_)
}
ids$class <- mapply(class_of, ids$lineage.sub, ids$lineage.total)
ids <- ids[!is.na(ids$class), c("barcode", "class")]

cn <- colnames(m)
if (mean(rownames(m) %in% ids$barcode) > 0.5 && mean(cn %in% ids$barcode) < 0.5) {
  message("transposing cells x genes -> genes x cells")
  m <- Matrix::t(m)
  cn <- colnames(m)
}
hit <- mean(cn %in% ids$barcode)
message("fraction of matrix columns matched to a kept barcode: ", signif(hit, 4))
if (hit < 0.2) {
  # try barcode after the first underscore, or before
  bare <- sub("^.*_", "", cn)
  hit2 <- mean(bare %in% ids$barcode)
  message("bare-barcode match: ", signif(hit2, 4))
  if (hit2 > hit) {
    colnames(m) <- bare
    cn <- bare
    hit <- hit2
  }
}
if (hit < 0.2) stop("could not match GSE205335 barcodes to the UMI matrix")

keep_cells <- cn %in% ids$barcode
m <- m[, keep_cells, drop = FALSE]
lab <- ids$class[match(colnames(m), ids$barcode)]
classes <- c("malignant", "cd8", "cd4", "nk", "b", "myeloid", "endothelial", "fibroblast")
message("n cells: ", paste(classes, sapply(classes, function(k) sum(lab == k)), sep = "=", collapse = " "))

genes <- rownames(m)
# Sum symbols if the matrix is already symbols; keep first-row names as-is.
# Aggregate with a cell-class indicator (sparse).
ind <- Matrix(0, nrow = ncol(m), ncol = length(classes), sparse = TRUE)
for (j in seq_along(classes)) ind[lab == classes[[j]], j] <- 1
sums <- m %*% ind
libs <- as.numeric(Matrix::colSums(m) %*% ind)
sums <- as.matrix(sums)
colnames(sums) <- classes
out <- data.frame(gene = genes, sums, check.names = FALSE)
write.table(out, file.path(out_dir, "gse205335_sums.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(class = classes, n_cells = as.integer(Matrix::colSums(ind)),
                        library_umi = libs),
            file.path(out_dir, "gse205335_ncells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote GSE205335 pseudobulk")

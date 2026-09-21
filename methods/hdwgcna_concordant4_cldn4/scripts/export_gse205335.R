#!/usr/bin/env Rscript
# Export author-malignant cells for the locked GSE205335 patients.
# Drops Normal* tissue. Does not load immune expression into the export.

.libPaths(c("/home/ubuntu/R/library", .libPaths()))
suppressPackageStartupMessages(library(Matrix))

geo <- "/tmp/geo_hdwgcna"
out <- "/tmp/hdwgcna_c4/gse205335_export"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rds_path <- "/tmp/hdwgcna_c4/gse205335.rds"
locked_path <- "/workspace/methods/hdwgcna_concordant4_cldn4/data/patient_units_locked.tsv"
universe_path <- "/tmp/hdwgcna_c4/universe_genes.txt"

say <- function(...) {
  msg <- paste(format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " "))
  cat(msg, "\n")
  flush(stdout())
}

norm_code <- function(x) {
  x <- toupper(gsub("-", "_", x, fixed = TRUE))
  sub("_[35]P$", "", x)
}

locked <- read.delim(locked_path, stringsAsFactors = FALSE)
keep_pats <- locked$unit_id[locked$dataset == "GSE205335"]
say("locked patients", length(keep_pats))

soft_lines <- readLines(gzfile(file.path(geo, "GSE205335_family.soft.gz")))
soft <- list()
cur <- list()
for (line in soft_lines) {
  if (startsWith(line, "^SAMPLE")) {
    if (!is.null(cur$title)) soft[[length(soft) + 1L]] <- cur
    cur <- list()
  } else if (startsWith(line, "!Sample_title")) {
    cur$title <- sub("^!Sample_title = ", "", line)
  } else if (startsWith(line, "!Sample_characteristics_ch1")) {
    val <- sub("^!Sample_characteristics_ch1 = ", "", line)
    if (grepl(": ", val, fixed = TRUE)) {
      kv <- strsplit(val, ": ", fixed = TRUE)[[1]]
      cur[[kv[[1]]]] <- paste(kv[-1], collapse = ": ")
    }
  }
}
if (!is.null(cur$title)) soft[[length(soft) + 1L]] <- cur

soft_map <- do.call(rbind, lapply(soft, function(s) {
  code <- if (!is.null(s$title) && grepl(" ", s$title, fixed = TRUE)) {
    sub("^\\S+\\s+", "", s$title)
  } else {
    NA_character_
  }
  data.frame(
    title = if (is.null(s$title)) NA_character_ else s$title,
    patient = if (is.null(s$patient)) NA_character_ else s$patient,
    tissue = if (is.null(s$tissue)) NA_character_ else s$tissue,
    code = toupper(gsub("-", "_", code, fixed = TRUE)),
    stringsAsFactors = FALSE
  )
}))
say("soft rows", nrow(soft_map), "unique codes", length(unique(soft_map$code)))

ident <- read.delim(gzfile(file.path(geo, "GSE205335_Lung_IO_CellIdentity.txt.gz")),
                    stringsAsFactors = FALSE, check.names = FALSE)
ident$orig_code <- norm_code(ident$orig.ident)
ident <- merge(ident, soft_map, by.x = "orig_code", by.y = "code", all.x = TRUE)
ident$tissue[is.na(ident$tissue)] <- ""
ident$is_normal <- grepl("^Normal", ident$tissue)
ident$author_malignant <- ident$lineage.sub == "Malignant cells"
say("ident rows", nrow(ident), "unmatched patient", sum(is.na(ident$patient)),
    "malignant", sum(ident$author_malignant %in% TRUE))

say("reading RDS", rds_path)
mat <- readRDS(rds_path)
say("class", paste(class(mat), collapse = ","), "dim", paste(dim(mat), collapse = "x"))
if (!inherits(mat, "dgCMatrix")) {
  mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
}
say("example colnames", paste(head(colnames(mat), 3), collapse = " | "))
say("example ident barcode", paste(head(ident$barcode, 3), collapse = " | "))

bc_ident <- ident$barcode
common <- intersect(colnames(mat), bc_ident)
say("direct intersect", length(common))
if (length(common) < 1000) {
  alt <- gsub("_", "-", colnames(mat), fixed = TRUE)
  hit <- alt %in% bc_ident
  say("underscore-swap intersect", sum(hit))
  if (sum(hit) < 1000) {
    alt2 <- gsub("-", "_", colnames(mat), fixed = TRUE)
    hit2 <- alt2 %in% bc_ident
    say("dash-swap intersect", sum(hit2))
    stop("barcode mismatch")
  }
  map <- ident
  rownames(map) <- map$barcode
  ident2 <- map[alt[hit], , drop = FALSE]
  ident2$mat_barcode <- colnames(mat)[hit]
} else {
  map <- ident
  rownames(map) <- map$barcode
  ident2 <- map[common, , drop = FALSE]
  ident2$mat_barcode <- common
}
say("matched cells", nrow(ident2))

# Tumor rows of locked patients, author malignant. which() drops NA flags.
keep <- which(ident2$patient %in% keep_pats &
              ident2$is_normal %in% FALSE &
              ident2$author_malignant %in% TRUE)
say("export cells", length(keep), "patients", length(unique(ident2$patient[keep])))
ident_keep <- ident2[keep, , drop = FALSE]
# per-patient counts for the audit
tab <- as.data.frame(table(ident_keep$patient), stringsAsFactors = FALSE)
colnames(tab) <- c("unit_id", "n_malignant_export")
write.table(tab, file.path(out, "patient_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mat <- mat[, ident_keep$mat_barcode, drop = FALSE]
rm(ident, ident2)
gc(verbose = FALSE)
say("subset dim", paste(dim(mat), collapse = "x"))

rn <- toupper(rownames(mat))
if (anyDuplicated(rn)) {
  say("collapsing duplicated symbols", sum(duplicated(rn)))
  grp <- split(seq_along(rn), rn)
  acc <- vector("list", length(grp))
  names(acc) <- names(grp)
  i <- 0L
  for (g in names(grp)) {
    ix <- grp[[g]]
    acc[[g]] <- if (length(ix) == 1L) mat[ix, ] else Matrix::colSums(mat[ix, , drop = FALSE])
    i <- i + 1L
    if (i %% 5000L == 0L) say("  collapsed", i)
  }
  # rebuild sparse from a dense chunk if needed: use rsparsematrix path via cbind of rows is slow.
  # Column-bind is wrong. Build dgCMatrix from row sums stored as a matrix in chunks.
  genes <- names(grp)
  chunk <- 2000L
  pieces <- list()
  for (start in seq(1L, length(genes), by = chunk)) {
    end <- min(length(genes), start + chunk - 1L)
    block <- do.call(rbind, acc[genes[start:end]])
    pieces[[length(pieces) + 1L]] <- as(block, "dgCMatrix")
    say("  block", start, end)
  }
  mat <- do.call(rbind, pieces)
  rownames(mat) <- genes
  rm(acc, pieces)
  gc(verbose = FALSE)
} else {
  rownames(mat) <- rn
}

universe <- readLines(universe_path)
universe <- unique(universe[nzchar(universe)])
in_u <- rownames(mat) %in% universe
say("genes in 123902 universe", sum(in_u), "of", nrow(mat))
mat <- mat[in_u, , drop = FALSE]
say("final dim", paste(dim(mat), collapse = "x"), "nnz", nnzero(mat))

meta <- data.frame(
  barcode = colnames(mat),
  unit_id = ident_keep$patient,
  tissue = ident_keep$tissue,
  stringsAsFactors = FALSE
)
write.table(meta, file.path(out, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(rownames(mat), file.path(out, "genes.txt"))
Matrix::writeMM(mat, file.path(out, "matrix.mtx"))
say("wrote", out)
# free the decompressed RDS; the .gz remains
unlink(rds_path)
say("removed decompressed RDS")

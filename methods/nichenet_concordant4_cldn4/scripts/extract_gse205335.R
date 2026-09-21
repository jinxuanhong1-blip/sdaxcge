#!/usr/bin/env Rscript
# GSE205335 only. Author malignant -> T/NK. Patient is the unit.
# Normal tissue is dropped. Cells from multiple GSMs of one patient are pooled.
# Does not read any other accession.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}
RAW <- parse_opt("--raw", "/tmp/concordant4_raw")
HERE <- parse_opt("--here", "methods/nichenet_concordant4_cldn4")
OUT <- parse_opt("--out", "/tmp/nichenet_work/extract")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ")

is_gzip <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && magic[[1]] == as.raw(0x1f) && magic[[2]] == as.raw(0x8b)
}

read_geo_rds <- function(path) {
  cur <- path
  temps <- character()
  on.exit(unlink(temps[file.exists(temps)]), add = TRUE)
  for (i in seq_len(4)) {
    if (!is_gzip(cur)) break
    dest <- tempfile(pattern = paste0("geo205335_", i, "_"), tmpdir = tempdir())
    temps <- c(temps, dest)
    logmsg("gzip -dc layer", i)
    st <- system2("gzip", c("-dc", cur), stdout = dest)
    if (!identical(st, 0L)) stop("gzip -dc failed")
    cur <- dest
  }
  logmsg("readRDS", file.info(cur)$size)
  readRDS(cur)
}

ALIASES <- c(PVRL2 = "NECTIN2", JAM1 = "F11R", PVRL1 = "NECTIN1", PVRL3 = "NECTIN3")
MIN_MAL <- 40L
MIN_TNK <- 20L

canon <- function(g) {
  g <- toupper(as.character(g))
  ifelse(g %in% names(ALIASES), unname(ALIASES[g]), g)
}

gsm <- read.delim(file.path(HERE, "data", "GSE205335_gsm_map.tsv"), stringsAsFactors = FALSE)
locked <- read.delim(file.path(HERE, "data", "GSE205335_patients.tsv"), stringsAsFactors = FALSE)
keep_pt <- as.character(locked$patient[locked$n_malignant > 0])
gsm <- gsm[gsm$patient %in% keep_pt, , drop = FALSE]
gsm$normal <- grepl("^Normal", gsm$tissue)
logmsg("GSMs", nrow(gsm), "normal dropped", sum(gsm$normal), "patients", length(unique(gsm$patient[!gsm$normal])))

ident <- read.delim(gzfile(file.path(RAW, "GSE205335", "GSE205335_Lung_IO_CellIdentity.txt.gz")),
                    stringsAsFactors = FALSE, check.names = FALSE)
mat <- read_geo_rds(file.path(RAW, "GSE205335", "GSE205335_Lung_IO_UMI_matrix.rds.gz"))
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
logmsg("matrix", nrow(mat), "x", ncol(mat))
rownames(mat) <- canon(rownames(mat))
if (anyDuplicated(rownames(mat))) {
  logmsg("collapsing duplicated gene symbols", sum(duplicated(rownames(mat))))
  mat <- rowsum(mat, rownames(mat))
}

common <- intersect(colnames(mat), ident$barcode)
if (length(common) < 1000) stop("barcode overlap too small: ", length(common))
mat <- mat[, common, drop = FALSE]
ident <- ident[match(common, ident$barcode), , drop = FALSE]
ident <- merge(ident, gsm[, c("orig.ident", "patient", "tissue", "normal")],
               by = "orig.ident", all.x = TRUE, sort = FALSE)
ident <- ident[match(colnames(mat), ident$barcode), , drop = FALSE]
ident$normal[is.na(ident$normal)] <- TRUE
ok <- !is.na(ident$patient) & ident$patient %in% keep_pt & !ident$normal
logmsg("cells kept", sum(ok), "of", ncol(mat))
mat <- mat[, ok, drop = FALSE]
ident <- ident[ok, , drop = FALSE]
lib <- Matrix::colSums(mat)
mal <- !is.na(ident$lineage.sub) & ident$lineage.sub == "Malignant cells"
tnk <- !is.na(ident$lineage.total) & ident$lineage.total == "T/NK cells" & !mal
logmsg("mal", sum(mal), "tnk", sum(tnk))

ligands <- readLines("/tmp/nichenet_prior/ligands.txt")
ligands <- intersect(canon(ligands), rownames(mat))
targets <- readLines("/tmp/nichenet_prior/targets.txt")
targets <- intersect(canon(targets), rownames(mat))
logmsg("ligands in matrix", length(ligands), "targets in matrix", length(targets))

mean_log_pct <- function(sub_mat, lib_sub) {
  n <- ncol(sub_mat)
  ng <- nrow(sub_mat)
  if (n == 0) return(list(mean = rep(NA_real_, ng), pct = rep(NA_real_, ng)))
  scale <- 1e4 / pmax(as.numeric(lib_sub), 1)
  sm <- summary(sub_mat)
  sum_log <- numeric(ng)
  nnz <- numeric(ng)
  if (length(sm$i)) {
    val <- log1p(sm$x * scale[sm$j])
    rs <- rowsum(cbind(val = val, nnz = 1), group = sm$i, reorder = TRUE)
    present <- as.integer(rownames(rs))
    sum_log[present] <- rs[, "val"]
    nnz[present] <- rs[, "nnz"]
  }
  list(mean = sum_log / n, pct = nnz / n)
}

quartile_idx <- function(x) {
  n <- length(x)
  if (n < 4) return(list(high = integer(), low = integer(), ok = FALSE))
  r <- rank(x, ties.method = "first")
  q1 <- floor(n * 0.25)
  q4 <- ceiling(n * 0.75)
  if (q1 < 1 || q4 >= n) return(list(high = integer(), low = integer(), ok = FALSE))
  list(high = which(r > q4), low = which(r <= q1), ok = TRUE)
}

patients <- unique(as.character(ident$patient))
patients <- intersect(keep_pt, patients)

unit_rows <- list()
lig_rows <- list()
# gene x patient matrices filled after we know units
gene_names <- targets
# include CLDN4 and ligands even if somehow not in targets
gene_names <- unique(c(gene_names, ligands, "CLDN4"))
gene_names <- intersect(gene_names, rownames(mat))
gindex <- match(gene_names, rownames(mat))

mean_mat <- matrix(NA_real_, nrow = length(gene_names), ncol = length(patients),
                   dimnames = list(gene_names, patients))
pct_mat <- mean_mat

for (pt in patients) {
  ix <- which(ident$patient == pt)
  mal_ix <- ix[mal[ix]]
  tnk_ix <- ix[tnk[ix]]
  n_mal <- length(mal_ix)
  n_tnk <- length(tnk_ix)
  logmsg(pt, "mal", n_mal, "tnk", n_tnk)
  cldn4_umi <- if ("CLDN4" %in% rownames(mat)) as.numeric(mat["CLDN4", mal_ix]) else rep(0, n_mal)
  cldn4_log <- log1p(1e4 * cldn4_umi / pmax(lib[mal_ix], 1))
  q <- quartile_idx(cldn4_log)
  med <- if (n_mal) stats::median(cldn4_log) else NA_real_
  hi_med <- if (n_mal) which(cldn4_log > med) else integer()
  lo_med <- if (n_mal) which(cldn4_log <= med) else integer()
  hi_pos <- if (n_mal) which(cldn4_umi > 0) else integer()
  lo_pos <- if (n_mal) which(cldn4_umi == 0) else integer()
  eligible <- n_mal >= MIN_MAL && n_tnk >= MIN_TNK && isTRUE(q$ok)
  separated <- eligible && length(q$high) && length(q$low) &&
    mean(cldn4_log[q$high]) > mean(cldn4_log[q$low])
  unit_rows[[pt]] <- data.frame(
    cohort = "GSE205335",
    unit_id = pt,
    n_cells = length(ix),
    n_mal = n_mal,
    n_tnk = n_tnk,
    n_high = length(q$high),
    n_low = length(q$low),
    cldn4_pct = if (n_mal) mean(cldn4_umi > 0) else NA_real_,
    cldn4_mean = if (n_mal) mean(cldn4_log) else NA_real_,
    cldn4_mean_high = if (length(q$high)) mean(cldn4_log[q$high]) else NA_real_,
    cldn4_mean_low = if (length(q$low)) mean(cldn4_log[q$low]) else NA_real_,
    eligible_q4 = eligible,
    cldn4_separated = separated,
    stringsAsFactors = FALSE
  )
  if (n_tnk >= 1) {
    agg <- mean_log_pct(mat[gindex, tnk_ix, drop = FALSE], lib[tnk_ix])
    mean_mat[, pt] <- agg$mean
    pct_mat[, pt] <- agg$pct
  }
  if (!eligible) next
  # ligand sender splits on malignant cells
  mal_mat <- mat[, mal_ix, drop = FALSE]
  lib_mal <- lib[mal_ix]
  arms <- list(
    q4 = q$high, q1 = q$low,
    med_hi = hi_med, med_lo = lo_med,
    pos = hi_pos, neg = lo_pos
  )
  arm_mean <- lapply(arms, function(w) {
    if (!length(w)) return(rep(NA_real_, length(ligands)))
    mean_log_pct(mal_mat[ligands, w, drop = FALSE], lib_mal[w])$mean
  })
  all_mean <- mean_log_pct(mal_mat[ligands, , drop = FALSE], lib_mal)
  for (i in seq_along(ligands)) {
    lig_rows[[length(lig_rows) + 1L]] <- data.frame(
      cohort = "GSE205335",
      unit_id = pt,
      ligand = ligands[[i]],
      delta_q4q1 = arm_mean$q4[[i]] - arm_mean$q1[[i]],
      mean_q4 = arm_mean$q4[[i]],
      mean_q1 = arm_mean$q1[[i]],
      delta_median = arm_mean$med_hi[[i]] - arm_mean$med_lo[[i]],
      delta_pctpos = arm_mean$pos[[i]] - arm_mean$neg[[i]],
      pct_mal = all_mean$pct[[i]],
      mean_mal = all_mean$mean[[i]],
      stringsAsFactors = FALSE
    )
  }
  rm(mal_mat); gc(verbose = FALSE)
}

units <- do.call(rbind, unit_rows)
lig <- do.call(rbind, lig_rows)
write.table(units, file.path(OUT, "units_GSE205335.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lig, file.path(OUT, "ligands_GSE205335.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(genes = gene_names, units = patients, mean = mean_mat, pct = pct_mat),
        file.path(OUT, "tnk_GSE205335.rds"))
logmsg("wrote", nrow(units), "units", nrow(lig), "ligand rows")

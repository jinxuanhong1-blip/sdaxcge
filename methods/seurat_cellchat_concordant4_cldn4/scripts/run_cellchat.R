#!/usr/bin/env Rscript
# Concordant-four CLDN4-high vs CLDN4-low senders.
# MUST: R + Seurat + CellChat (jinworks). No Python primary. No GSE148071.
# Honest unit = patient / locked sample. Primary split = malignant Q4 vs Q1.

suppressPackageStartupMessages({
  lib <- Sys.getenv("R_LIBS_USER", "/tmp/r_lib")
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed. Stop. Do not fall back to a Python-only primary.")
  }
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat is not installed. Stop. Do not fall back to a Python-only primary.")
  }
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(data.table)
})

options(warn = 1)
set.seed(1)

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)[[1]]) a else b

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}

RAW <- parse_opt("--raw", "/tmp/concordant4_raw")
HERE <- tryCatch({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  normalizePath(file.path(dirname(f), ".."))
}, error = function(e) normalizePath("methods/seurat_cellchat_concordant4_cldn4"))
OUT <- parse_opt("--out", HERE)

DIR_RES <- file.path(OUT, "results")
DIR_FIG <- file.path(OUT, "figures")
DIR_TAB <- file.path(DIR_RES, "tables")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

logmsg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ", flush = TRUE)
}

# ---------------------------------------------------------------------------
# Pre-specified pairs (CellChatDB v2 protein names). Thesis already correct.
# ---------------------------------------------------------------------------
PAIRS <- data.frame(
  interaction_name = c(
    "JAM1_ITGAL_ITGB2", "NECTIN2_TIGIT", "CDH1_ITGAE_ITGB7", "CDH1_KLRG1",
    "LGALS9_HAVCR2", "LGALS9_CD44", "LGALS9_CD45",
    "CXCL9_CXCR3", "CXCL10_CXCR3", "CCL5_CCR5", "CCL5_CCR1",
    "HLA-A_CD8A", "HLA-B_CD8A", "HLA-C_CD8A"
  ),
  axis = c(
    "F11R", "NECTIN2-TIGIT", "CDH1", "CDH1",
    "LGALS9", "LGALS9", "LGALS9",
    "CXCL9/10-CXCR3", "CXCL9/10-CXCR3", "CCL5", "CCL5",
    "HLA-CD8", "HLA-CD8", "HLA-CD8"
  ),
  family = c(
    rep("barrier_inhibitory", 7),
    rep("ifn_recruit_mhci", 7)
  ),
  thesis_expect = c(
    rep("high>low", 7),
    rep("low>high", 7)
  ),
  stringsAsFactors = FALSE
)

EPI <- c("EPCAM", "KRT8", "KRT18", "KRT19")
TNK_MARKERS <- c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
ALIASES <- c(PVRL2 = "NECTIN2", JAM1 = "F11R")
MALIG_SUB <- c("Malignant cells", "tS1", "tS2", "tS3")
TNK_TYPES <- c("T lymphocytes", "NK cells")
MIN_ARM <- 10L
MIN_TNK <- 20L
MIN_MAL_Q4 <- 40L

WANTED <- unique(c(
  "CLDN4", EPI, TNK_MARKERS, "PTPRC",
  "F11R", "JAM1", "ITGAL", "ITGB2", "NECTIN2", "PVRL2", "TIGIT",
  "CDH1", "ITGAE", "ITGB7", "KLRG1", "LGALS9", "HAVCR2", "CD44",
  "CXCL9", "CXCL10", "CXCR3", "CCL5", "CCR5", "CCR1",
  "HLA-A", "HLA-B", "HLA-C", "CD8A"
))

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-3) sprintf("%.2e", p) else sprintf("%.3g", p)
}
fmt_num <- function(x, d = 3) {
  if (!is.finite(x)) return("NA")
  sprintf(paste0("%+.", d, "f"), x)
}

is_gzip <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && magic[[1]] == as.raw(0x1f) && magic[[2]] == as.raw(0x8b)
}

# Peel one or more gzip layers (GEO files are sometimes double-gzipped).
# Never write a giant uncompressed object to a durable path.
read_geo_rds <- function(path) {
  if (!file.exists(path)) stop("missing ", path)
  cur <- path
  temps <- character()
  on.exit(unlink(temps[file.exists(temps)]), add = TRUE)
  for (i in seq_len(4)) {
    if (!is_gzip(cur)) break
    dest <- tempfile(pattern = paste0("geo_rds_", i, "_"), tmpdir = tempdir())
    temps <- c(temps, dest)
    logmsg("  gzip -dc layer", i, basename(cur), "->", dest)
    st <- system2("gzip", c("-dc", cur), stdout = dest)
    if (!identical(st, 0L) && !is.null(st) && st != 0) {
      stop("gzip -dc failed on ", cur, " status=", st)
    }
    if (!file.exists(dest) || file.info(dest)$size < 10) {
      stop("gzip -dc produced empty file from ", cur)
    }
    cur <- dest
  }
  logmsg("  readRDS", cur, "bytes", file.info(cur)$size)
  obj <- try(readRDS(cur), silent = TRUE)
  if (inherits(obj, "try-error")) {
    logmsg("  readRDS failed; trying gzfile connection")
    obj <- readRDS(gzfile(cur, open = "rb"))
  }
  obj
}

apply_aliases <- function(mat) {
  rn <- as.character(rownames(mat))
  rn <- ifelse(rn %in% names(ALIASES), unname(ALIASES[rn]), rn)
  # keep first occurrence after alias
  if (anyDuplicated(rn)) {
    keep <- !duplicated(rn)
    mat <- mat[keep, , drop = FALSE]
    rn <- rn[keep]
  }
  rownames(mat) <- rn
  mat
}

gene_row <- function(mat, g) {
  if (!g %in% rownames(mat)) return(rep(0, ncol(mat)))
  as.numeric(mat[g, ])
}

pos_any <- function(mat, genes) {
  hit <- intersect(genes, rownames(mat))
  if (!length(hit)) return(rep(FALSE, ncol(mat)))
  if (length(hit) == 1L) return(as.numeric(mat[hit, ]) > 0)
  as.numeric(Matrix::colSums(mat[hit, , drop = FALSE] > 0)) > 0
}

subset_genes <- function(mat, extra = character()) {
  mat <- apply_aliases(mat)
  keep <- intersect(unique(c(WANTED, extra)), rownames(mat))
  if (!length(keep)) stop("none of the wanted genes are in the matrix")
  mat[keep, , drop = FALSE]
}

get_layer <- function(obj, which = "counts") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    GetAssayData(obj, assay = "RNA", layer = which)
  } else {
    GetAssayData(obj, assay = "RNA", slot = which)
  }
}

quartile_high_low <- function(x) {
  # ties.method="first" keeps four arms when many cells share CLDN4=0.
  r <- rank(as.numeric(x), ties.method = "first")
  n <- length(r)
  if (n < 4) {
    return(list(high = rep(FALSE, n), low = rep(FALSE, n), ok = FALSE))
  }
  q1 <- floor(n * 0.25)
  q4 <- ceiling(n * 0.75)
  if (q1 < 1 || q4 > n || q1 >= q4) {
    return(list(high = rep(FALSE, n), low = rep(FALSE, n), ok = FALSE))
  }
  list(high = r > q4, low = r <= q1, ok = TRUE)
}

median_high_low <- function(x) {
  med <- stats::median(as.numeric(x), na.rm = TRUE)
  list(high = x > med, low = x <= med, ok = is.finite(med))
}

pctpos_high_low <- function(umi) {
  list(high = umi > 0, low = umi == 0, ok = TRUE)
}

subset_cellchat_db <- function() {
  db <- CellChatDB.human
  inter <- db$interaction
  keep <- inter$interaction_name %in% PAIRS$interaction_name
  if (!any(keep)) {
    # some builds store names without hyphen variants
    keep <- toupper(inter$interaction_name) %in% toupper(PAIRS$interaction_name)
  }
  db$interaction <- inter[keep, , drop = FALSE]
  if (!nrow(db$interaction)) {
    stop("CellChatDB.human has none of the pre-specified pairs")
  }
  db
}

run_one_cellchat <- function(counts, group) {
  stopifnot(ncol(counts) == length(group))
  keep <- !is.na(group) & group %in% c("CLDN4_high", "CLDN4_low", "TNK")
  counts <- counts[, keep, drop = FALSE]
  group <- factor(as.character(group[keep]),
                  levels = c("CLDN4_high", "CLDN4_low", "TNK"))
  if (length(unique(group)) < 3) {
    return(list(ok = FALSE, reason = "missing_group", df = NULL))
  }
  tab <- table(group)
  if (any(tab < MIN_ARM)) {
    return(list(ok = FALSE, reason = paste0("n_floor:", paste(tab, collapse = "/")), df = NULL))
  }

  obj <- CreateSeuratObject(counts = counts, project = "unit", min.cells = 0, min.features = 0)
  obj$group <- as.character(group)
  Idents(obj) <- obj$group
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  data.input <- get_layer(obj, "data")
  meta <- data.frame(labels = obj$group,
                     samples = factor("sample1"),
                     row.names = colnames(obj),
                     stringsAsFactors = FALSE)

  cellchat <- createCellChat(object = as.matrix(data.input), meta = meta,
                             group.by = "labels")
  cellchat@DB <- subset_cellchat_db()
  cellchat <- subsetData(cellchat)
  # Keep the pre-specified pairs even if they are not DE among 3 groups.
  cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 1, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  if (!is.null(cellchat@DB$interaction) && nrow(cellchat@DB$interaction)) {
    cellchat@LR$LRsig <- cellchat@DB$interaction
  }
  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    nboot = 1,
    population.size = TRUE
  )
  df <- tryCatch(
    subsetCommunication(cellchat, slot.name = "net"),
    error = function(e) NULL
  )
  list(ok = TRUE, reason = "ok", df = df, n = as.list(tab),
       seurat = paste(as.character(packageVersion("Seurat")), collapse = "."),
       cellchat = paste(as.character(packageVersion("CellChat")), collapse = "."))
}

extract_pair_rows <- function(df, cohort, patient, split, n_high, n_low, n_tnk, n_mal) {
  out <- PAIRS
  out$cohort <- cohort
  out$patient <- patient
  out$split <- split
  out$n_mal <- n_mal
  out$n_high <- n_high
  out$n_low <- n_low
  out$n_tnk <- n_tnk
  out$prob_high <- NA_real_
  out$prob_low <- NA_real_
  out$pval_high <- NA_real_
  out$pval_low <- NA_real_
  out$detected <- FALSE
  out$delta <- NA_real_
  if (is.null(df) || !nrow(df)) return(out)
  nm <- names(df)
  iname <- if ("interaction_name" %in% nm) "interaction_name" else if ("interaction_name_2" %in% nm) "interaction_name_2" else NA
  if (is.na(iname)) return(out)
  for (i in seq_len(nrow(out))) {
    hi <- df[df$source == "CLDN4_high" & df$target == "TNK" &
               df[[iname]] == out$interaction_name[i], , drop = FALSE]
    lo <- df[df$source == "CLDN4_low" & df$target == "TNK" &
               df[[iname]] == out$interaction_name[i], , drop = FALSE]
    if (nrow(hi)) {
      out$prob_high[i] <- as.numeric(hi$prob[[1]])
      if ("pval" %in% names(hi)) out$pval_high[i] <- as.numeric(hi$pval[[1]])
    } else {
      out$prob_high[i] <- 0
    }
    if (nrow(lo)) {
      out$prob_low[i] <- as.numeric(lo$prob[[1]])
      if ("pval" %in% names(lo)) out$pval_low[i] <- as.numeric(lo$pval[[1]])
    } else {
      out$prob_low[i] <- 0
    }
    out$detected[i] <- is.finite(out$prob_high[i]) && is.finite(out$prob_low[i]) &&
      (out$prob_high[i] > 0 || out$prob_low[i] > 0)
  }
  out$delta <- out$prob_high - out$prob_low
  out
}

score_unit <- function(counts, mal, tnk, cohort, patient, splits = c("q4q1", "median", "pctpos")) {
  rows <- list()
  inv <- list()
  mal <- as.logical(mal) %in% TRUE
  tnk <- as.logical(tnk) %in% TRUE
  n_mal <- as.integer(sum(mal))
  n_tnk <- as.integer(sum(tnk))
  cldn4_umi <- gene_row(counts, "CLDN4")
  # log1p CP10k on malignant only for ranking
  lib <- Matrix::colSums(counts)
  cldn4_log <- log1p(cldn4_umi / pmax(lib, 1) * 1e4)
  inv[[1]] <- data.frame(
    cohort = cohort, patient = patient, n_cells = ncol(counts),
    n_mal = n_mal, n_tnk = n_tnk,
    mal_cldn4_mean = if (isTRUE(n_mal > 0)) mean(cldn4_log[mal]) else NA_real_,
    mal_cldn4_pct = if (isTRUE(n_mal > 0)) mean(cldn4_umi[mal] > 0) else NA_real_,
    stringsAsFactors = FALSE
  )
  if (n_tnk < MIN_TNK || n_mal < MIN_ARM * 2) {
    return(list(lr = NULL, inv = inv[[1]]))
  }
  for (sp in splits) {
    if (sp == "q4q1") {
      if (n_mal < MIN_MAL_Q4) next
      hl <- quartile_high_low(cldn4_log[mal])
    } else if (sp == "median") {
      hl <- median_high_low(cldn4_log[mal])
    } else {
      hl <- pctpos_high_low(cldn4_umi[mal])
    }
    if (!isTRUE(hl$ok)) next
    high <- mal
    low <- mal
    high[mal] <- hl$high
    low[mal] <- hl$low
    n_high <- as.integer(sum(high))
    n_low <- as.integer(sum(low))
    if (n_high < MIN_ARM || n_low < MIN_ARM) next
    group <- rep(NA_character_, ncol(counts))
    group[high] <- "CLDN4_high"
    group[low] <- "CLDN4_low"
    group[tnk] <- "TNK"
    keep <- !is.na(group)
    logmsg("  CellChat", cohort, patient, sp,
           "high", n_high, "low", n_low, "tnk", n_tnk)
    cc <- tryCatch(
      run_one_cellchat(counts[, keep, drop = FALSE], group[keep]),
      error = function(e) list(ok = FALSE, reason = conditionMessage(e), df = NULL)
    )
    if (!isTRUE(cc$ok)) {
      logmsg("    skip", cc$reason %||% "fail")
      next
    }
    rows[[length(rows) + 1]] <- extract_pair_rows(
      cc$df, cohort, patient, sp, n_high, n_low, n_tnk, n_mal
    )
  }
  list(lr = if (length(rows)) do.call(rbind, rows) else NULL, inv = inv[[1]])
}

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
load_gse123902 <- function() {
  logmsg("==== GSE123902 CSV ====")
  units <- read.delim(file.path(HERE, "data", "GSE123902_marker_units.tsv"),
                      stringsAsFactors = FALSE)
  tumor <- units[units$tissue %in% c("PRIMARY", "METASTASIS"), , drop = FALSE]
  tumor <- tumor[order(tumor$patient, ifelse(tumor$tissue == "PRIMARY", 0, 1)), ]
  tumor <- tumor[!duplicated(tumor$patient), ]
  keep <- as.character(tumor$patient)
  tar <- file.path(RAW, "GSE123902", "GSE123902_RAW.tar")
  csv_dir <- file.path(RAW, "GSE123902", "csv")
  dir.create(csv_dir, showWarnings = FALSE)
  if (!length(list.files(csv_dir, pattern = "_dense\\.csv\\.gz$"))) {
    if (!file.exists(tar)) stop("missing ", tar)
    logmsg("untar GSE123902")
    untar(tar, exdir = csv_dir)
  }
  files <- list.files(csv_dir, pattern = "_dense\\.csv\\.gz$", full.names = TRUE)
  out_lr <- list(); out_inv <- list()
  for (i in seq_len(nrow(tumor))) {
    patient <- as.character(tumor$patient[i])
    want <- as.character(tumor$file[i])
    fp <- file.path(csv_dir, want)
    if (!file.exists(fp)) {
      hit <- grep(paste0("_", patient, "_"), files, value = TRUE)
      hit <- hit[!grepl("_NORMAL_", hit)]
      hit <- hit[order(!grepl("_PRIMARY_", hit))]
      if (!length(hit)) {
        logmsg("  missing CSV for", patient)
        next
      }
      fp <- hit[[1]]
    }
    bn <- basename(fp)
    logmsg("  fread", bn)
    dt <- data.table::fread(cmd = paste("gzip -dc", shQuote(fp)), sep = ",",
                            header = TRUE, data.table = FALSE, showProgress = FALSE)
    genes <- toupper(colnames(dt)[-1])
    mat <- t(as.matrix(dt[, -1, drop = FALSE]))
    storage.mode(mat) <- "double"
    rownames(mat) <- genes
    colnames(mat) <- paste0(patient, "_", as.character(dt[[1]]))
    rm(dt); gc(verbose = FALSE)
    mat <- subset_genes(Matrix::Matrix(mat, sparse = TRUE))
    mal <- pos_any(mat, EPI) & (gene_row(mat, "PTPRC") == 0)
    tnk <- pos_any(mat, TNK_MARKERS) & (!mal)
    rec <- score_unit(mat, mal, tnk, "GSE123902", patient)
    if (!is.null(rec$inv)) out_inv[[length(out_inv) + 1]] <- rec$inv
    if (!is.null(rec$lr)) out_lr[[length(out_lr) + 1]] <- rec$lr
    rm(mat); gc(verbose = FALSE)
  }
  list(lr = out_lr, inv = out_inv)
}

stream_gene_matrix <- function(path, wanted) {
  cache <- file.path(dirname(path), "wanted_stream_cache.rds")
  if (file.exists(cache) && file.info(cache)$size > 1000) {
    logmsg("  load stream cache", cache)
    return(readRDS(cache))
  }
  wanted_u <- unique(c(wanted, names(ALIASES), unname(ALIASES)))
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  cell_ids <- header[-1]
  n <- length(cell_ids)
  found <- list()
  n_streamed <- 0L
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    n_streamed <- n_streamed + 1L
    tab <- regexpr("\t", line, fixed = TRUE)
    gene <- if (tab[[1]] > 0) substr(line, 1, tab[[1]] - 1) else line
    gene_u <- toupper(gene)
    if (gene_u %in% wanted_u || gene %in% wanted_u) {
      rest <- if (tab[[1]] > 0) substr(line, tab[[1]] + 1, nchar(line)) else ""
      vals <- as.numeric(strsplit(rest, "\t", fixed = TRUE)[[1]])
      if (length(vals) != n) {
        # pad / trim
        if (length(vals) < n) vals <- c(vals, rep(0, n - length(vals)))
        if (length(vals) > n) vals <- vals[seq_len(n)]
      }
      key <- if (gene_u %in% names(ALIASES)) unname(ALIASES[[gene_u]]) else gene_u
      found[[key]] <- vals
      logmsg("    kept", key, "nz", sum(vals > 0))
    }
    if (n_streamed %% 5000L == 0L) logmsg("    streamed", n_streamed, "genes")
  }
  if (!length(found)) stop("no wanted genes in ", path)
  mat <- Matrix::Matrix(do.call(rbind, found), sparse = TRUE)
  rownames(mat) <- names(found)
  colnames(mat) <- cell_ids
  saveRDS(mat, cache)
  logmsg("  wrote cache", cache, "genes", nrow(mat), "cells", ncol(mat))
  mat
}

load_gse131907 <- function() {
  logmsg("==== GSE131907 UMI TXT ====")
  samples <- read.delim(file.path(HERE, "data", "GSE131907_samples.tsv"),
                        stringsAsFactors = FALSE)
  keep <- as.character(samples$sample[samples$n_malignant > 0])
  annot <- data.table::fread(
    cmd = paste("gzip -dc", shQuote(file.path(RAW, "GSE131907",
      "GSE131907_Lung_Cancer_cell_annotation.txt.gz"))),
    sep = "\t", header = TRUE, data.table = FALSE
  )
  # tolerate column name variants
  cn <- names(annot)
  cell_col <- cn[tolower(cn) %in% c("index", "cell", "barcode", "cell_id")][[1]]
  samp_col <- cn[grepl("sample|orig", tolower(cn))][[1]]
  type_col <- cn[grepl("cell_type$", tolower(cn)) | tolower(cn) == "cell_type"][[1]]
  sub_col <- cn[grepl("subtype|cell_subtype", tolower(cn))][[1]]
  if (is.na(cell_col) || is.na(samp_col)) {
    logmsg("  annot columns:", paste(cn, collapse = ", "))
    stop("cannot parse GSE131907 annotation columns")
  }
  cells <- as.character(annot[[cell_col]])
  samp <- as.character(annot[[samp_col]])
  subtype <- if (!is.na(sub_col)) as.character(annot[[sub_col]]) else rep("", nrow(annot))
  ctype <- if (!is.na(type_col)) as.character(annot[[type_col]]) else rep("", nrow(annot))
  names(samp) <- cells
  names(subtype) <- cells
  names(ctype) <- cells

  mat <- stream_gene_matrix(
    file.path(RAW, "GSE131907", "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"),
    WANTED
  )
  common <- intersect(colnames(mat), cells)
  mat <- mat[, common, drop = FALSE]
  samp <- samp[common]; subtype <- subtype[common]; ctype <- ctype[common]
  out_lr <- list(); out_inv <- list()
  for (s in keep) {
    idx <- which(samp == s)
    if (!length(idx)) {
      logmsg("  missing sample", s)
      next
    }
    sub <- mat[, idx, drop = FALSE]
    mal <- subtype[idx] %in% MALIG_SUB
    tnk <- ctype[idx] %in% TNK_TYPES
    rec <- score_unit(sub, mal, tnk, "GSE131907", s)
    if (!is.null(rec$inv)) out_inv[[length(out_inv) + 1]] <- rec$inv
    if (!is.null(rec$lr)) out_lr[[length(out_lr) + 1]] <- rec$lr
    rm(sub); gc(verbose = FALSE)
  }
  rm(mat); gc(verbose = FALSE)
  list(lr = out_lr, inv = out_inv)
}

parse_soft <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  recs <- list(); cur <- NULL; titles <- character(); descs <- character()
  flush_cur <- function() {
    if (is.null(cur)) return(invisible(NULL))
    cur$title <- if (length(titles)) titles[[1]] else ""
    cur$description <- if (length(descs)) descs[[1]] else ""
    recs[[length(recs) + 1]] <<- cur
  }
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "^SAMPLE = ")) {
      flush_cur()
      cur <- list(gsm = sub("^SAMPLE = ", "", line))
      titles <- character(); descs <- character()
    } else if (!is.null(cur) && startsWith(line, "!Sample_title = ")) {
      titles <- c(titles, sub("^!Sample_title = ", "", line))
    } else if (!is.null(cur) && startsWith(line, "!Sample_description = ")) {
      descs <- c(descs, sub("^!Sample_description = ", "", line))
    } else if (!is.null(cur) && startsWith(line, "!Sample_characteristics_ch1 = ")) {
      val <- sub("^!Sample_characteristics_ch1 = ", "", line)
      if (grepl(": ", val, fixed = TRUE)) {
        key <- gsub(" ", "_", sub(": .*", "", val))
        cur[[key]] <- sub("^[^:]+: ", "", val)
      }
    }
  }
  flush_cur()
  md <- do.call(rbind, lapply(recs, function(x) {
    keys <- c("gsm", "patient", "tissue", "tumor_stage", "cancer_subtype",
              "recist", "platform", "description", "title")
    vals <- vapply(keys, function(k) {
      v <- x[[k]]
      if (is.null(v) || !length(v)) "" else as.character(v[[1]])
    }, character(1))
    as.data.frame(as.list(vals), stringsAsFactors = FALSE)
  }))
  read_end <- sub(".*Single Cell ([35])'.*", "\\1", md$platform)
  read_end[!grepl("Single Cell [35]'", md$platform)] <- NA_character_
  md$orig.ident <- paste0(gsub("_", "-", md$description, fixed = TRUE), "-", read_end, "P")
  md
}

load_gse205335 <- function() {
  logmsg("==== GSE205335 RDS (double-gzip safe) ====")
  ident <- utils::read.delim(
    file.path(RAW, "GSE205335", "GSE205335_Lung_IO_CellIdentity.txt.gz"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  gsm_map <- read.delim(file.path(HERE, "data", "GSE205335_gsm_map.tsv"),
                        stringsAsFactors = FALSE)
  rds_gz <- file.path(RAW, "GSE205335", "GSE205335_Lung_IO_UMI_matrix.rds.gz")
  mat <- read_geo_rds(rds_gz)
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!inherits(mat, "dgCMatrix") && !inherits(mat, "Matrix") && !is.matrix(mat)) {
    stop("GSE205335 RDS class=", paste(class(mat), collapse = ","))
  }
  if (!inherits(mat, "Matrix")) mat <- Matrix::Matrix(mat, sparse = TRUE)
  logmsg("  full matrix", nrow(mat), "x", ncol(mat), class(mat)[1])
  mat <- subset_genes(mat)
  logmsg("  subset", nrow(mat), "x", ncol(mat))
  gc(verbose = FALSE)

  bc <- ident$barcode
  if (is.null(bc)) stop("identity missing barcode")
  common <- intersect(colnames(mat), bc)
  mat <- mat[, common, drop = FALSE]
  ident <- ident[match(colnames(mat), ident$barcode), , drop = FALSE]
  ident <- merge(ident, gsm_map[, c("orig.ident", "patient", "tissue"), drop = FALSE],
                 by = "orig.ident", all.x = TRUE, sort = FALSE)
  ident <- ident[match(colnames(mat), ident$barcode), , drop = FALSE]
  if (anyNA(ident$patient) || any(ident$patient == "")) {
    stop("unmapped orig.ident: ",
         paste(unique(ident$orig.ident[is.na(ident$patient) | ident$patient == ""]),
               collapse = ", "))
  }
  is_normal <- grepl("^Normal ", ident$tissue %||% "")
  is_normal[is.na(is_normal)] <- FALSE
  mal <- !is.na(ident$lineage.sub) & ident$lineage.sub == "Malignant cells" & !is_normal
  tnk <- !is.na(ident$lineage.total) & ident$lineage.total == "T/NK cells" & !is_normal
  logmsg("  mapped patients", length(unique(ident$patient)),
         "mal", sum(mal), "tnk", sum(tnk))
  locked <- read.delim(file.path(HERE, "data", "GSE205335_patients.tsv"),
                       stringsAsFactors = FALSE)
  keep <- as.character(locked$patient[locked$n_malignant > 0])
  out_lr <- list(); out_inv <- list()
  for (pt in keep) {
    idx <- which(ident$patient == pt & !is_normal)
    if (!length(idx)) next
    rec <- score_unit(mat[, idx, drop = FALSE], mal[idx], tnk[idx], "GSE205335", pt)
    if (!is.null(rec$inv)) out_inv[[length(out_inv) + 1]] <- rec$inv
    if (!is.null(rec$lr)) out_lr[[length(out_lr) + 1]] <- rec$lr
    gc(verbose = FALSE)
  }
  rm(mat); gc(verbose = FALSE)
  list(lr = out_lr, inv = out_inv)
}

load_gse189357 <- function() {
  logmsg("==== GSE189357 MTX ====")
  meta <- read.delim(file.path(HERE, "data", "GSE189357_sample_metadata.tsv"),
                     stringsAsFactors = FALSE)
  tar <- file.path(RAW, "GSE189357", "GSE189357_RAW.tar")
  ex <- file.path(RAW, "GSE189357", "raw")
  dir.create(ex, showWarnings = FALSE)
  needed <- paste0(meta$gsm, "_", meta$patient, "_matrix.mtx.gz")
  if (!all(file.exists(file.path(ex, needed)))) {
    if (!file.exists(tar)) stop("missing ", tar)
    logmsg("untar GSE189357")
    untar(tar, exdir = ex)
  }
  out_lr <- list(); out_inv <- list()
  for (i in seq_len(nrow(meta))) {
    patient <- meta$patient[i]
    gsm <- meta$gsm[i]
    prefix <- file.path(ex, paste0(gsm, "_", patient))
    mtx <- paste0(prefix, "_matrix.mtx.gz")
    cells <- paste0(prefix, "_barcodes.tsv.gz")
    features <- paste0(prefix, "_features.tsv.gz")
    if (!all(file.exists(c(mtx, cells, features)))) {
      # some tars flatten names
      alt <- list.files(ex, pattern = paste0(patient, "_matrix\\.mtx\\.gz$"), full.names = TRUE)
      if (length(alt)) {
        mtx <- alt[[1]]
        cells <- sub("_matrix\\.mtx\\.gz$", "_barcodes.tsv.gz", mtx)
        features <- sub("_matrix\\.mtx\\.gz$", "_features.tsv.gz", mtx)
      }
    }
    logmsg("  ReadMtx", patient)
    mat <- Seurat::ReadMtx(
      mtx = mtx, cells = cells, features = features,
      cell.column = 1, feature.column = 2, unique.features = TRUE
    )
    colnames(mat) <- paste0(patient, "_", colnames(mat))
    mat <- subset_genes(mat)
    mal <- pos_any(mat, EPI) & (gene_row(mat, "PTPRC") == 0)
    tnk <- pos_any(mat, TNK_MARKERS) & (!mal)
    rec <- score_unit(mat, mal, tnk, "GSE189357", patient)
    if (!is.null(rec$inv)) out_inv[[length(out_inv) + 1]] <- rec$inv
    if (!is.null(rec$lr)) out_lr[[length(out_lr) + 1]] <- rec$lr
    rm(mat); gc(verbose = FALSE)
  }
  list(lr = out_lr, inv = out_inv)
}

summarize_split <- function(lr, split) {
  d <- lr[lr$split == split & lr$detected %in% TRUE, , drop = FALSE]
  if (!nrow(d)) d <- lr[lr$split == split, , drop = FALSE]
  rows <- list()
  for (nm in unique(PAIRS$interaction_name)) {
    sub <- d[d$interaction_name == nm, , drop = FALSE]
    sub <- sub[is.finite(sub$delta), , drop = FALSE]
    n <- nrow(sub)
    n_by <- as.list(table(sub$cohort))
    if (n >= 2) {
      wt <- suppressWarnings(stats::wilcox.test(sub$delta, mu = 0, exact = FALSE))
      p <- wt$p.value
    } else {
      p <- NA_real_
    }
    mean_d <- if (n) mean(sub$delta) else NA_real_
    obs <- if (!n || !is.finite(mean_d)) "NA" else if (mean_d > 0) "high>low" else if (mean_d < 0) "low>high" else "tie"
    expect <- PAIRS$thesis_expect[PAIRS$interaction_name == nm][[1]]
    agrees <- if (obs == "NA") "thin" else if (obs == expect) "yes" else if (sign(mean_d) == 0) "tie" else "opposite"
    rows[[length(rows) + 1]] <- data.frame(
      pair = nm,
      axis = PAIRS$axis[PAIRS$interaction_name == nm][[1]],
      family = PAIRS$family[PAIRS$interaction_name == nm][[1]],
      expect = expect,
      split = split,
      n = n,
      n_gse123902 = n_by$GSE123902 %||% 0,
      n_gse131907 = n_by$GSE131907 %||% 0,
      n_gse205335 = n_by$GSE205335 %||% 0,
      n_gse189357 = n_by$GSE189357 %||% 0,
      mean_delta = mean_d,
      p_W = p,
      observed = obs,
      agrees = agrees,
      stringsAsFactors = FALSE
    )
  }
  # family aggregates: per-patient mean of detected pair deltas
  for (fam in unique(PAIRS$family)) {
    subp <- d[d$family == fam & is.finite(d$delta), , drop = FALSE]
    if (!nrow(subp)) next
    agg <- aggregate(delta ~ cohort + patient, data = subp, FUN = mean)
    n <- nrow(agg)
    n_by <- as.list(table(agg$cohort))
    if (n >= 2) {
      p <- suppressWarnings(stats::wilcox.test(agg$delta, mu = 0, exact = FALSE)$p.value)
    } else p <- NA_real_
    mean_d <- mean(agg$delta)
    expect <- unique(PAIRS$thesis_expect[PAIRS$family == fam])[[1]]
    obs <- if (mean_d > 0) "high>low" else if (mean_d < 0) "low>high" else "tie"
    agrees <- if (obs == expect) "yes" else "opposite"
    rows[[length(rows) + 1]] <- data.frame(
      pair = paste0("FAMILY_", fam),
      axis = "family",
      family = fam,
      expect = expect,
      split = split,
      n = n,
      n_gse123902 = n_by$GSE123902 %||% 0,
      n_gse131907 = n_by$GSE131907 %||% 0,
      n_gse205335 = n_by$GSE205335 %||% 0,
      n_gse189357 = n_by$GSE189357 %||% 0,
      mean_delta = mean_d,
      p_W = p,
      observed = obs,
      agrees = agrees,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

write_finding <- function(lig, inv, per, versions) {
  q <- lig[lig$split == "q4q1", ]
  md <- lig[lig$split == "median", ]
  inv_n <- aggregate(patient ~ cohort, data = inv, FUN = function(x) length(unique(x)))
  both <- inv[inv$n_mal >= 10 & inv$n_tnk >= 20, ]
  q_per <- per[per$split == "q4q1" & is.finite(per$delta), , drop = FALSE]
  q_units <- unique(q_per[, c("cohort", "patient"), drop = FALSE])
  md_path <- file.path(OUT, "FINDING.md")
  line_pair <- function(tab, nm) {
    r <- tab[tab$pair == nm, ]
    if (!nrow(r)) return("| | | | | | | | | | | | |")
    r <- r[1, ]
    sprintf("| %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |",
            r$pair, r$axis, r$expect, r$n,
            r$n_gse123902, r$n_gse131907, r$n_gse205335, r$n_gse189357,
            fmt_num(r$mean_delta), fmt_p(r$p_W), r$observed, r$agrees)
  }
  fam_q <- q[q$pair == "FAMILY_barrier_inhibitory", ]
  con <- file(md_path, open = "wt")
  on.exit(close(con))
  writeLines(c(
    "# FINDING — Seurat + CellChat concordant-four CLDN4-high vs low senders",
    "",
    "ADDITIVE. **Thesis already correct. Ligands stay.**",
    "CLDN4 only. No dual-high. Concordant four only",
    "(GSE123902 + GSE131907 + GSE205335 + GSE189357).",
    "Do **not** add GSE148071 / GSE127465 / CD45-only. This is **not** a full-pool.",
    "",
    "Engine: **R + Seurat + CellChat** (not a Python reimplementation).",
    sprintf("Seurat %s. CellChat %s.", versions$seurat, versions$cellchat),
    "Outgoing communication probability is CellChat `computeCommunProb`",
    "(truncatedMean, trim=0.1, population.size=TRUE).",
    "Senders = malignant CLDN4-high vs CLDN4-low; receivers = T/NK.",
    "Honest n = patient / locked sample.",
    "",
    "Thesis:",
    "",
    "- Barrier/inhibitory outgoing **UP from CLDN4-high** (F11R, NECTIN2–TIGIT, CDH1, LGALS9) is **ON-thesis**.",
    "- IFN/T-recruit outgoing UP from CLDN4-low is the KD-like arm (often **not** detected; CXCL9/10 sparse).",
    "  Do not bury barrier-up-in-high as a recruit-up skip.",
    "",
    "Primary split is malignant **Q4 vs Q1**. Extra: median and %pos.",
    "",
    "## Honest n",
    "",
    "Locked four n=65 (13+21+22+9) is **not** the ligand n.",
    "",
    "| gate | n | note |",
    "|---|---:|---|",
    sprintf("| Locked four | 65 | 13+21+22+9 |"),
    sprintf("| Inventory units loaded | %d | %s |",
            nrow(inv),
            paste(sprintf("%s=%s", inv_n$cohort, inv_n$patient), collapse = ", ")),
    sprintf("| Both compartments (n_mal≥10, n_tnk≥20) | %d | honest inventory |", nrow(both)),
    sprintf("| Q4 vs Q1 CellChat units | **%d** | primary ligand n |", nrow(q_units)),
    "",
    "GSE123902 / GSE189357 malignant labels are thin → epithelium marker-malignant",
    "(EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0) → T/NK.",
    "GSE131907 / GSE205335 use author malignant → T/NK.",
    "TACSTD2 is never a gate. PVRL2 is aliased to NECTIN2 on GSE131907.",
    "",
    "## Primary — barrier / inhibitory family ΔP (high − low)",
    "",
    "ON-thesis. Primary table: `results/tables/ligand_table.tsv` (Q4 vs Q1).",
    "",
    "| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|",
    line_pair(q, "JAM1_ITGAL_ITGB2"),
    line_pair(q, "NECTIN2_TIGIT"),
    line_pair(q, "CDH1_ITGAE_ITGB7"),
    line_pair(q, "CDH1_KLRG1"),
    line_pair(q, "LGALS9_HAVCR2"),
    line_pair(q, "LGALS9_CD44"),
    line_pair(q, "LGALS9_CD45"),
    line_pair(q, "FAMILY_barrier_inhibitory"),
    "",
    if (nrow(fam_q)) {
      sprintf("Family aggregate: n=%d mean ΔP=%s p_W=%s agrees=%s.",
              fam_q$n[[1]], fmt_num(fam_q$mean_delta[[1]]),
              fmt_p(fam_q$p_W[[1]]), fam_q$agrees[[1]])
    } else "Family aggregate: not computed.",
    "",
    "## KD-like arm — IFN / T-recruit / MHC-I (expect low > high)",
    "",
    "Do not file the barrier result as a recruit-up skip.",
    "",
    "| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|",
    line_pair(q, "CXCL9_CXCR3"),
    line_pair(q, "CXCL10_CXCR3"),
    line_pair(q, "CCL5_CCR5"),
    line_pair(q, "CCL5_CCR1"),
    line_pair(q, "HLA-A_CD8A"),
    line_pair(q, "HLA-B_CD8A"),
    line_pair(q, "HLA-C_CD8A"),
    line_pair(q, "FAMILY_ifn_recruit_mhci"),
    "",
    "## Extra: median-split barrier family",
    "",
    "| pair | axis | expect | n | 123902 | 131907 | 205335 | 189357 | mean ΔP | p_W | observed | agrees |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|",
    line_pair(md, "JAM1_ITGAL_ITGB2"),
    line_pair(md, "NECTIN2_TIGIT"),
    line_pair(md, "CDH1_ITGAE_ITGB7"),
    line_pair(md, "CDH1_KLRG1"),
    line_pair(md, "LGALS9_HAVCR2"),
    line_pair(md, "LGALS9_CD44"),
    line_pair(md, "LGALS9_CD45"),
    line_pair(md, "FAMILY_barrier_inhibitory"),
    "",
    "## What is not claimed",
    "",
    "- TACSTD2 is not used to define high/low. This is not dual-high.",
    "- GSE148071, GSE127465, and CD45-only libraries are not added.",
    "- This is not a CellChat discovery screen and not a 7-pool.",
    "- Cell-pooled tests are not reported. Honest n is the patient/sample.",
    "",
    "## Reproduce",
    "",
    "```bash",
    "Rscript methods/seurat_cellchat_concordant4_cldn4/scripts/install_packages.R",
    "bash methods/seurat_cellchat_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw",
    "Rscript methods/seurat_cellchat_concordant4_cldn4/scripts/run_cellchat.R --raw=/tmp/concordant4_raw",
    "```",
    "",
    "CellChat: type=truncatedMean, trim=0.1, population.size=TRUE.",
    "Arm floor ≥10; T/NK ≥20; Q4 vs Q1 also n_mal≥40.",
    ""
  ), con)
  logmsg("wrote", md_path)
}

main <- function() {
  logmsg("Seurat", as.character(packageVersion("Seurat")),
         "CellChat", as.character(packageVersion("CellChat")))
  if (!exists("CreateSeuratObject")) stop("CreateSeuratObject missing")
  if (!exists("createCellChat")) stop("createCellChat missing")
  if (!exists("computeCommunProb")) stop("computeCommunProb missing")

  versions <- list(
    seurat = as.character(packageVersion("Seurat")),
    cellchat = as.character(packageVersion("CellChat"))
  )
  writeLines(
    c(paste("Seurat", versions$seurat),
      paste("CellChat", versions$cellchat),
      capture.output(sessionInfo())),
    file.path(DIR_RES, "sessionInfo.txt")
  )

  parts <- list(
    load_gse123902(),
    load_gse131907(),
    load_gse205335(),
    load_gse189357()
  )
  lr_list <- unlist(lapply(parts, `[[`, "lr"), recursive = FALSE)
  inv_list <- unlist(lapply(parts, `[[`, "inv"), recursive = FALSE)
  if (!length(lr_list)) stop("CellChat produced no per-unit rows")
  per <- do.call(rbind, lr_list)
  inv <- do.call(rbind, inv_list)
  rownames(per) <- NULL
  rownames(inv) <- NULL

  lig_q <- summarize_split(per, "q4q1")
  lig_m <- summarize_split(per, "median")
  lig_p <- summarize_split(per, "pctpos")
  lig <- rbind(lig_q, lig_m, lig_p)

  utils::write.table(per, file.path(DIR_TAB, "per_patient_lr.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(inv, file.path(DIR_TAB, "patient_inventory.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(lig, file.path(DIR_TAB, "ligand_table.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(lig[lig$split == "q4q1", ],
                     file.path(DIR_TAB, "ligand_table_q4q1.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  # also top-level ligand table as requested
  utils::write.table(lig[lig$split == "q4q1", ],
                     file.path(DIR_RES, "ligand_table.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  write_finding(lig, inv, per, versions)
  logmsg("DONE units", length(unique(paste(per$cohort, per$patient))),
         "rows", nrow(per))
}

main()

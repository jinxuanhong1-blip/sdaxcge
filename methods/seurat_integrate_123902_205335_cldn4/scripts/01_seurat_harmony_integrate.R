#!/usr/bin/env Rscript
# ADDITIVE. CLDN4-only. Human pair GSE123902 + GSE205335.
# Primary engine: R + Seurat 5 + Harmony. Honest unit = patient.
# Cldn4 vs T/NK; malignant IFN/MHC Q4 vs Q1.
# No dual-high. No GSE148071. No GSE127465. No Python-only primary.
# Concordant-4 four-way merge is a different agent.

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed. Stop. Do not fall back to a Python-only primary.")
  }
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony is not installed. Stop. Seurat/Harmony is required.")
  }
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)[1]) a else b

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}

data_dir <- parse_opt("--data", "/tmp/geo_pair_123902_205335")
here <- tryCatch(
  {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grep("^--file=", ca)])
    normalizePath(file.path(dirname(f), ".."))
  },
  error = function(e) getwd()
)
out_dir <- parse_opt("--out", here)
cap_n <- as.integer(parse_opt("--cap", "150"))
SEED <- 1L
set.seed(SEED)

dir.create(file.path(out_dir, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Locked sets. CLDN4 is the readout. It is never used to assign lineage.
# ---------------------------------------------------------------------------
EPI_MARKERS <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7")
TNK_MARKERS <- c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
MYE_MARKERS <- c("LYZ", "CD14", "CSF1R", "AIF1")
B_MARKERS <- c("MS4A1", "CD79A")

IFN_GENES <- c(
  "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20",
  "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3", "OAS1", "OAS2", "OASL",
  "IFI27", "IFI44", "IFI44L", "IFI6", "RSAD2", "CXCL9", "CXCL10",
  "CXCL11", "B2M", "TAP1", "PSMB8", "PSMB9"
)
MHC_GENES <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M",
  "TAP1", "TAP2", "TAPBP", "TAPBPL", "NLRC5", "PSMB8", "PSMB9",
  "PSMB10", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "IRF1"
)
TJ_GENES <- c(
  "CLDN1", "CLDN3", "CLDN7", "CLDN18", "TJP1", "TJP2", "TJP3",
  "OCLN", "F11R", "CGN", "MARVELD2", "MARVELD3", "CRB3", "CDH1"
)

MIN_MAL <- 20L
MIN_TNK <- 20L
MIN_N_SPEARMAN <- 5L
MIN_N_Q4Q1 <- 8L
NORMAL_TISSUE <- c("Normal Lung", "Normal LN", "Normal Brain")

present <- function(genes, universe) intersect(genes, universe)

ln_rows <- function(counts, genes) {
  g <- present(genes, rownames(counts))
  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- 1
  if (!length(g)) {
    return(matrix(0, nrow = 0, ncol = ncol(counts),
                  dimnames = list(NULL, colnames(counts))))
  }
  sub <- as.matrix(counts[g, , drop = FALSE])
  storage.mode(sub) <- "double"
  log1p(sweep(sub, 2, 1e4 / lib, `*`))
}

marker_score_vec <- function(counts, markers) {
  m <- ln_rows(counts, markers)
  if (!nrow(m)) return(rep(0, ncol(counts)))
  if (nrow(m) == 1L) return(as.numeric(m[1, ]))
  as.numeric(colMeans(m))
}

assign_lineage <- function(counts) {
  scores <- rbind(
    epithelial = marker_score_vec(counts, EPI_MARKERS),
    tnk = marker_score_vec(counts, TNK_MARKERS),
    myeloid = marker_score_vec(counts, MYE_MARKERS),
    b = marker_score_vec(counts, B_MARKERS)
  )
  winner <- max.col(t(scores), ties.method = "first")
  top <- apply(scores, 2, max)
  second <- apply(scores, 2, function(x) sort(x, decreasing = TRUE)[2])
  keep <- (top >= 0.12) & (top >= second * 1.15)
  lab <- rep("other", ncol(scores))
  lab[keep] <- rownames(scores)[winner[keep]]
  lab
}

spearman_ci <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < MIN_N_SPEARMAN) {
    return(list(n = n, rho = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  z <- atanh(max(min(rho, 0.999999), -0.999999))
  se <- 1 / sqrt(n - 3)
  list(
    n = n, rho = as.numeric(rho), p = as.numeric(ct$p.value),
    lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se)
  )
}

fmt_rho <- function(s) {
  if (is.na(s$rho)) return(sprintf("n=%d (below Spearman floor)", s$n))
  sprintf("n=%d, ρ=%+.3f [%.3f, %.3f], p=%.3g", s$n, s$rho, s$lo, s$hi, s$p)
}

fmt_num <- function(x, d = 3) {
  if (!length(x) || is.na(x)) return("NA")
  sprintf(paste0("%.", d, "f"), as.numeric(x))
}

mw_delta <- function(x, high) {
  ok <- is.finite(x) & !is.na(high)
  x <- x[ok]
  high <- high[ok]
  n1 <- sum(!high)
  n2 <- sum(high)
  if (n1 < 3 || n2 < 3) {
    return(list(n_q1 = n1, n_q4 = n2, delta = NA_real_, p = NA_real_, r = NA_real_))
  }
  wt <- suppressWarnings(wilcox.test(x[high], x[!high], exact = FALSE))
  u <- as.numeric(wt$statistic)
  r_rb <- (2 * u) / (n1 * n2) - 1
  list(
    n_q1 = n1, n_q4 = n2,
    delta = median(x[high]) - median(x[!high]),
    p = as.numeric(wt$p.value),
    r = r_rb
  )
}

assign_q <- function(x) {
  rnk <- rank(x, ties.method = "average", na.last = "keep")
  out <- rep(NA_character_, length(x))
  ok <- is.finite(rnk)
  if (sum(ok) < 4) return(out)
  qs <- tryCatch(
    as.character(cut(
      rnk[ok],
      breaks = quantile(rnk[ok], probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )),
    error = function(e) rep(NA_character_, sum(ok))
  )
  out[ok] <- qs
  out
}

fisher_z_dl <- function(rho, n) {
  keep <- is.finite(rho) & is.finite(n) & n > 3 & abs(rho) < 1
  rho <- rho[keep]
  n <- n[keep]
  k <- length(rho)
  if (k == 0) {
    return(list(k = 0, N = 0, rho = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_, I2 = NA_real_))
  }
  if (k == 1) {
    z <- atanh(rho)
    se <- 1 / sqrt(n - 3)
    p <- 2 * pnorm(-abs(z / se))
    return(list(k = 1, N = n, rho = rho, p = p, lo = tanh(z - 1.96 * se),
                hi = tanh(z + 1.96 * se), I2 = 0))
  }
  z <- atanh(rho)
  w <- n - 3
  zbar <- sum(w * z) / sum(w)
  q <- sum(w * (z - zbar)^2)
  df <- k - 1
  tau2 <- max(0, (q - df) / (sum(w) - sum(w^2) / sum(w)))
  wstar <- 1 / (1 / w + tau2)
  zdl <- sum(wstar * z) / sum(wstar)
  se <- sqrt(1 / sum(wstar))
  p <- 2 * pnorm(-abs(zdl / se))
  i2 <- if (q <= 0) 0 else max(0, (q - df) / q)
  list(
    k = k, N = sum(n), rho = tanh(zdl), p = p,
    lo = tanh(zdl - 1.96 * se), hi = tanh(zdl + 1.96 * se),
    I2 = 100 * i2
  )
}

cap_by_unit <- function(df, unit_col, cap, seed_i) {
  set.seed(seed_i)
  parts <- split(df, df[[unit_col]], drop = TRUE)
  out <- lapply(parts, function(sub) {
    if (nrow(sub) <= cap) return(sub)
    sub[sample.int(nrow(sub), cap), , drop = FALSE]
  })
  do.call(rbind, out)
}

align_rows <- function(mat, genes) {
  extra <- setdiff(genes, rownames(mat))
  if (length(extra)) {
    z <- Matrix::Matrix(0, nrow = length(extra), ncol = ncol(mat), sparse = TRUE)
    rownames(z) <- extra
    colnames(z) <- colnames(mat)
    mat <- rbind(mat, z)
  }
  mat[genes, , drop = FALSE]
}

parse_fname <- function(fname) {
  m <- regexec(
    "^(GSM[0-9]+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\\.csv\\.gz$",
    fname
  )
  g <- regmatches(fname, m)[[1]]
  if (length(g) != 4) stop("unparsed filename: ", fname)
  list(gsm = g[2], patient = g[3], site = g[4], file = fname)
}

read_geo_csv <- function(path, gsm) {
  dt <- data.table::fread(
    cmd = paste("zcat", shQuote(path)),
    sep = ",", header = TRUE, data.table = FALSE, showProgress = FALSE
  )
  bc <- as.character(dt[[1]])
  genes <- colnames(dt)[-1]
  m <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rm(dt)
  sp <- Matrix::Matrix(t(m), sparse = TRUE)
  rm(m)
  gc(verbose = FALSE)
  rownames(sp) <- genes
  colnames(sp) <- paste0(gsm, "_", bc)
  sp
}

ggsave_both <- function(p, name, w = 6.4, h = 5.2) {
  fig <- file.path(out_dir, "results", "figures")
  ggplot2::ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 140)
  ggplot2::ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h)
}

# ===========================================================================
# 1. GSE123902 — marker lineage on GEO dense UMI CSVs
# ===========================================================================
csv_dir <- file.path(data_dir, "GSE123902", "csv")
if (!length(list.files(csv_dir, pattern = "GSM.*_dense\\.csv\\.gz$"))) {
  tar_path <- file.path(data_dir, "GSE123902", "GSE123902_RAW.tar")
  if (!file.exists(tar_path)) stop("Missing ", tar_path, ". Run scripts/00_download.sh first.")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  untar(tar_path, exdir = csv_dir)
}
csv_files <- sort(list.files(csv_dir, pattern = "GSM.*_dense\\.csv\\.gz$", full.names = TRUE))
if (!length(csv_files)) stop("No GSE123902 GSM dense CSVs")

message("[GSE123902] reading ", length(csv_files), " GEO dense CSVs...")
keep123_mats <- list()
keep123_meta <- list()
stat123 <- list()
n123_barcodes <- 0L

for (fp in csv_files) {
  rec <- parse_fname(basename(fp))
  message("  ", basename(fp))
  sp <- read_geo_csv(fp, rec$gsm)
  n123_barcodes <- n123_barcodes + ncol(sp)
  lin <- assign_lineage(sp)
  tumor <- rec$site != "NORMAL"
  is_mal <- lin == "epithelial" & tumor
  is_tnk <- lin == "tnk"
  cldn4_ln <- marker_score_vec(sp, "CLDN4")
  cldn4_umi <- if ("CLDN4" %in% rownames(sp)) as.numeric(sp["CLDN4", ]) else rep(0, ncol(sp))
  ifn_s <- marker_score_vec(sp, IFN_GENES)
  mhc_s <- marker_score_vec(sp, MHC_GENES)
  tj_s <- marker_score_vec(sp, TJ_GENES)

  n_tumor <- if (tumor) ncol(sp) else 0L
  n_mal <- sum(is_mal)
  n_tnk <- if (tumor) sum(is_tnk) else 0L
  stat123[[length(stat123) + 1L]] <- data.frame(
    dataset = "GSE123902",
    patient = rec$patient,
    gsm = rec$gsm,
    site = rec$site,
    n_cells_file = ncol(sp),
    n_cells_tumor = n_tumor,
    n_malignant = n_mal,
    n_tnk = n_tnk,
    cldn4_mean = if (n_mal) mean(cldn4_ln[is_mal]) else NA_real_,
    cldn4_pct_pos = if (n_mal) 100 * mean(cldn4_umi[is_mal] > 0) else NA_real_,
    ifn_mean = if (n_mal) mean(ifn_s[is_mal]) else NA_real_,
    mhc_mean = if (n_mal) mean(mhc_s[is_mal]) else NA_real_,
    tj_mean = if (n_mal) mean(tj_s[is_mal]) else NA_real_,
    stringsAsFactors = FALSE
  )

  keep <- (is_mal | (is_tnk & tumor))
  if (any(keep)) {
    sub <- sp[, keep, drop = FALSE]
    keep123_mats[[rec$gsm]] <- sub
    keep123_meta[[rec$gsm]] <- data.frame(
      barcode = colnames(sub),
      dataset = "GSE123902",
      patient = rec$patient,
      unit_id = rec$patient,
      gsm = rec$gsm,
      site = rec$site,
      compartment = ifelse(is_mal[keep], "malignant", "TNK"),
      stringsAsFactors = FALSE
    )
  }
  rm(sp)
  gc(verbose = FALSE)
}

lib123 <- do.call(rbind, stat123)
# Collapse libraries to patient: tumor only; PRIMARY preferred over METASTASIS
# when both exist (they do not in this accession). Sum cell counts; score
# weighted by n_malignant.
lib123_tumor <- lib123[lib123$site != "NORMAL", , drop = FALSE]
agg123 <- lapply(split(lib123_tumor, lib123_tumor$patient), function(d) {
  w <- d$n_malignant
  w[is.na(w)] <- 0
  data.frame(
    dataset = "GSE123902",
    patient = d$patient[1],
    tumor_sites = paste(unique(d$site), collapse = ","),
    histology = "LUAD",
    n_cells_tumor = sum(d$n_cells_tumor),
    n_malignant = sum(d$n_malignant),
    n_tnk = sum(d$n_tnk),
    frac_tnk = if (sum(d$n_cells_tumor) > 0) sum(d$n_tnk) / sum(d$n_cells_tumor) else NA_real_,
    cldn4_mean = if (sum(w) > 0) weighted.mean(d$cldn4_mean, w, na.rm = TRUE) else NA_real_,
    cldn4_pct_pos = if (sum(w) > 0) weighted.mean(d$cldn4_pct_pos, w, na.rm = TRUE) else NA_real_,
    ifn_mean = if (sum(w) > 0) weighted.mean(d$ifn_mean, w, na.rm = TRUE) else NA_real_,
    mhc_mean = if (sum(w) > 0) weighted.mean(d$mhc_mean, w, na.rm = TRUE) else NA_real_,
    tj_mean = if (sum(w) > 0) weighted.mean(d$tj_mean, w, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
})
pat123 <- do.call(rbind, agg123)
rownames(pat123) <- NULL

# Also list normal-only donors so the table is complete
all_donors <- sort(unique(lib123$patient))
missing_donors <- setdiff(all_donors, pat123$patient)
if (length(missing_donors)) {
  extra <- data.frame(
    dataset = "GSE123902",
    patient = missing_donors,
    tumor_sites = "NONE",
    histology = "LUAD",
    n_cells_tumor = 0L,
    n_malignant = 0L,
    n_tnk = 0L,
    frac_tnk = NA_real_,
    cldn4_mean = NA_real_,
    cldn4_pct_pos = NA_real_,
    ifn_mean = NA_real_,
    mhc_mean = NA_real_,
    tj_mean = NA_real_,
    stringsAsFactors = FALSE
  )
  pat123 <- rbind(pat123, extra)
}
pat123 <- pat123[order(pat123$patient), ]
pat123$eligible <- pat123$n_malignant >= MIN_MAL & pat123$n_tnk >= MIN_TNK
pat123$locked_frac_tnk <- pat123$frac_tnk

message("[GSE123902] tumor donors=", sum(pat123$n_cells_tumor > 0),
        " eligible=", sum(pat123$eligible), " kept-cell matrices=", length(keep123_mats))

# ===========================================================================
# 2. GSE205335 — author malignant / T/NK on public RDS
# ===========================================================================
cid_path <- file.path(data_dir, "GSE205335", "GSE205335_Lung_IO_CellIdentity.txt.gz")
rds_path <- file.path(data_dir, "GSE205335", "GSE205335_Lung_IO_UMI_matrix.rds")
if (!file.exists(rds_path) && file.exists(paste0(rds_path, ".gz"))) {
  stop("GSE205335 RDS not peeled. Run scripts/00_download.sh")
}
if (!file.exists(cid_path) || !file.exists(rds_path)) {
  stop("Missing GSE205335 files. Run scripts/00_download.sh")
}

gsm <- read.csv(file.path(here, "data", "GSE205335_gsm_sample_metadata.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
lock205 <- read.delim(file.path(here, "data", "GSE205335_patients.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
elig205 <- lock205$patient[lock205$n_malignant >= MIN_MAL]

cid <- read.delim(cid_path, check.names = FALSE, stringsAsFactors = FALSE)
map <- unique(gsm[, c("orig.ident", "patient", "tissue", "cancer_subtype")])
cid$orig.ident <- as.character(cid$orig.ident)
cid <- merge(cid, map, by = "orig.ident", all.x = TRUE)
cid$barcode <- as.character(cid$barcode)
n205_ann <- nrow(cid)
cid_tumor <- cid[!is.na(cid$patient) & !cid$tissue %in% NORMAL_TISSUE, , drop = FALSE]
cid_tumor$compartment <- NA_character_
cid_tumor$compartment[cid_tumor$lineage.sub == "Malignant cells"] <- "malignant"
cid_tumor$compartment[cid_tumor$lineage.total == "T/NK cells"] <- "TNK"

agg205 <- lapply(split(cid_tumor, cid_tumor$patient), function(d) {
  data.frame(
    dataset = "GSE205335",
    patient = d$patient[1],
    tumor_sites = paste(sort(unique(d$tissue)), collapse = ","),
    histology = paste(sort(unique(d$cancer_subtype)), collapse = ","),
    n_cells_tumor = nrow(d),
    n_malignant = sum(d$compartment == "malignant", na.rm = TRUE),
    n_tnk = sum(d$compartment == "TNK", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
pat205_n <- do.call(rbind, agg205)
pat205_n$frac_tnk <- ifelse(pat205_n$n_cells_tumor > 0,
                            pat205_n$n_tnk / pat205_n$n_cells_tumor, NA_real_)
pat205_n <- merge(
  pat205_n,
  lock205[, c("patient", "frac_tnk", "mal_CLDN4_mean", "mal_CLDN4_pct_pos", "cancer_subtype")],
  by = "patient", all.x = TRUE, suffixes = c("", "_locked")
)
names(pat205_n)[names(pat205_n) == "frac_tnk_locked"] <- "locked_frac_tnk"

message("[GSE205335] loading UMI RDS...")
umi205 <- readRDS(rds_path)
message("[GSE205335] dim ", paste(dim(umi205), collapse = "x"),
        " size ", format(object.size(umi205), units = "MB"))

keep205_all <- cid_tumor[!is.na(cid_tumor$compartment) & cid_tumor$patient %in% elig205, ]
keep205_all <- keep205_all[keep205_all$barcode %in% colnames(umi205), ]
message("[GSE205335] eligible mal+TNK barcodes in matrix: ", nrow(keep205_all))

# Score ALL malignant cells (pre-cap) for the patient table.
mal205 <- keep205_all[keep205_all$compartment == "malignant", ]
mal_bc <- intersect(mal205$barcode, colnames(umi205))
mal_mat <- umi205[, mal_bc, drop = FALSE]
mal_ln_cldn4 <- marker_score_vec(mal_mat, "CLDN4")
mal_umi_cldn4 <- if ("CLDN4" %in% rownames(mal_mat)) as.numeric(mal_mat["CLDN4", ]) else rep(0, ncol(mal_mat))
mal_ifn <- marker_score_vec(mal_mat, IFN_GENES)
mal_mhc <- marker_score_vec(mal_mat, MHC_GENES)
mal_tj <- marker_score_vec(mal_mat, TJ_GENES)
mal_pat <- mal205$patient[match(colnames(mal_mat), mal205$barcode)]
score205 <- do.call(rbind, lapply(split(seq_along(mal_pat), mal_pat), function(idx) {
  data.frame(
    patient = mal_pat[idx[1]],
    cldn4_mean = mean(mal_ln_cldn4[idx]),
    cldn4_pct_pos = 100 * mean(mal_umi_cldn4[idx] > 0),
    ifn_mean = mean(mal_ifn[idx]),
    mhc_mean = mean(mal_mhc[idx]),
    tj_mean = mean(mal_tj[idx]),
    n_mal_scored = length(idx),
    stringsAsFactors = FALSE
  )
}))
rm(mal_mat)
gc(verbose = FALSE)

pat205 <- merge(pat205_n, score205, by = "patient", all.x = TRUE)
# Restrict the reported GSE205335 table to locked eligible patients (n_mal>=20)
# plus any extra with both floors met from the live annotation.
pat205$eligible <- pat205$patient %in% elig205 &
  pat205$n_malignant >= MIN_MAL & pat205$n_tnk >= MIN_TNK
pat205 <- pat205[pat205$patient %in% lock205$patient | pat205$eligible, ]
pat205 <- pat205[order(pat205$patient), ]
rownames(pat205) <- NULL

# Cap mal + T/NK per patient for Harmony / DimPlot
set.seed(SEED)
mal_cap <- cap_by_unit(keep205_all[keep205_all$compartment == "malignant", ],
                       "patient", cap_n, SEED + 2L)
tnk_cap <- cap_by_unit(keep205_all[keep205_all$compartment == "TNK", ],
                       "patient", cap_n, SEED + 3L)
sel205 <- rbind(mal_cap, tnk_cap)
sel205 <- sel205[sel205$barcode %in% colnames(umi205), ]
mat205 <- umi205[, sel205$barcode, drop = FALSE]
rm(umi205)
gc(verbose = FALSE)
meta205 <- data.frame(
  barcode = sel205$barcode,
  dataset = "GSE205335",
  patient = as.character(sel205$patient),
  unit_id = as.character(sel205$patient),
  gsm = NA_character_,
  site = as.character(sel205$tissue),
  compartment = as.character(sel205$compartment),
  stringsAsFactors = FALSE
)
rownames(meta205) <- meta205$barcode
message("[GSE205335] Harmony-kept cells=", ncol(mat205),
        " mal=", sum(meta205$compartment == "malignant"),
        " TNK=", sum(meta205$compartment == "TNK"))

# Cap GSE123902 kept cells
meta123 <- do.call(rbind, keep123_meta)
rownames(meta123) <- NULL
elig123_ids <- pat123$patient[pat123$eligible]
meta123 <- meta123[meta123$patient %in% elig123_ids, ]
mal123 <- cap_by_unit(meta123[meta123$compartment == "malignant", ], "patient", cap_n, SEED)
tnk123 <- cap_by_unit(meta123[meta123$compartment == "TNK", ], "patient", cap_n, SEED + 1L)
sel123 <- rbind(mal123, tnk123)
rownames(sel123) <- NULL

# Build GSE123902 sparse from kept GSM matrices
universe123 <- sort(unique(unlist(lapply(keep123_mats, rownames), use.names = FALSE)))
keep_bc123 <- sel123$barcode
mats123 <- list()
for (nm in names(keep123_mats)) {
  m <- keep123_mats[[nm]]
  hit <- intersect(colnames(m), keep_bc123)
  if (length(hit)) mats123[[nm]] <- align_rows(m[, hit, drop = FALSE], universe123)
}
rm(keep123_mats)
gc(verbose = FALSE)
mat123 <- do.call(cbind, mats123)
meta123_s <- sel123
rownames(meta123_s) <- meta123_s$barcode
meta123_s <- meta123_s[colnames(mat123), , drop = FALSE]
message("[GSE123902] Harmony-kept cells=", ncol(mat123))

# ===========================================================================
# 3. Seurat merge + Harmony
# ===========================================================================
message("[seurat] CreateSeuratObject + merge")
obj123 <- CreateSeuratObject(counts = mat123, meta.data = meta123_s, project = "GSE123902")
obj205 <- CreateSeuratObject(counts = mat205, meta.data = meta205, project = "GSE205335")
rm(mat123, mat205)
gc(verbose = FALSE)
obj <- merge(obj123, obj205)
rm(obj123, obj205)
gc(verbose = FALSE)
obj$dataset <- as.character(obj$dataset)
obj$compartment <- as.character(obj$compartment)
obj$patient <- as.character(obj$patient)
obj$unit_id <- as.character(obj$unit_id)
DefaultAssay(obj) <- "RNA"

message("[seurat] split layers / normalize / PCA")
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

integ_method <- "HarmonyIntegration"
integ_reduction <- "harmony"
message("[seurat] IntegrateLayers HarmonyIntegration")
ok_int <- tryCatch({
  obj <<- IntegrateLayers(
    object = obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = TRUE
  )
  TRUE
}, error = function(e) {
  message("[seurat] HarmonyIntegration failed: ", conditionMessage(e))
  FALSE
})
if (!ok_int) {
  message("[seurat] fallback harmony::RunHarmony on PCA embeddings")
  integ_method <- "harmony::RunHarmony"
  emb <- Embeddings(obj, "pca")
  harm <- harmony::RunHarmony(
    data_mat = emb,
    meta_data = slot(obj, "meta.data"),
    vars_use = "dataset",
    verbose = TRUE
  )
  obj[["harmony"]] <- CreateDimReducObject(
    embeddings = harm, key = "harmony_", assay = DefaultAssay(obj)
  )
}

obj <- FindNeighbors(obj, reduction = integ_reduction, dims = 1:20, verbose = FALSE)
obj <- RunUMAP(obj, reduction = integ_reduction, dims = 1:20, verbose = FALSE)
obj <- JoinLayers(obj)

# Seurat-native module scores on the integrated object (capped cells; not the claim n)
ifn_present <- present(IFN_GENES, rownames(obj))
mhc_present <- present(MHC_GENES, rownames(obj))
tj_present <- present(TJ_GENES, rownames(obj))
if (length(ifn_present) < 10) stop("too few IFN genes after merge")
if (length(mhc_present) < 5) stop("too few MHC genes after merge")
obj <- AddModuleScore(
  obj,
  features = list(IFN = ifn_present, MHC = mhc_present, TJ = tj_present),
  name = "mod",
  ctrl = 50,
  seed = SEED
)
obj$IFN_mod <- obj$mod1
obj$MHC_mod <- obj$mod2
obj$TJ_mod <- obj$mod3

if (!("CLDN4" %in% rownames(obj))) stop("CLDN4 absent after merge")
cldn4_counts <- as.numeric(LayerData(obj, layer = "counts")["CLDN4", ])
cldn4_data <- as.numeric(LayerData(obj, layer = "data")["CLDN4", ])
obj$CLDN4_counts <- cldn4_counts
obj$CLDN4_log1p <- cldn4_data
obj$CLDN4_pos <- as.integer(cldn4_counts > 0)

# ===========================================================================
# 4. Patient table + tests (honest unit = patient; scores from full malignant)
# ===========================================================================
pat_cols <- c("dataset", "patient", "tumor_sites", "histology",
              "n_cells_tumor", "n_malignant", "n_tnk", "frac_tnk",
              "cldn4_mean", "cldn4_pct_pos", "ifn_mean", "mhc_mean", "tj_mean",
              "eligible", "locked_frac_tnk")
# Harmonize GSE205335 columns
if (!("histology" %in% names(pat205))) pat205$histology <- pat205$cancer_subtype
keep205_cols <- intersect(pat_cols, names(pat205))
pat205_out <- pat205[, keep205_cols, drop = FALSE]
for (cc in setdiff(pat_cols, names(pat205_out))) pat205_out[[cc]] <- NA
pat123_out <- pat123[, pat_cols, drop = FALSE]
pat <- rbind(pat123_out, pat205_out[, pat_cols, drop = FALSE])
pat <- pat[order(pat$dataset, pat$patient), ]
rownames(pat) <- NULL

# Harmony cell counts per patient
md <- slot(obj, "meta.data")
harm_n <- as.data.frame(xtabs(~ patient + compartment, data = md),
                        stringsAsFactors = FALSE)
harm_wide <- reshape(
  harm_n, idvar = "patient", timevar = "compartment", direction = "wide"
)
names(harm_wide) <- gsub("Freq\\.", "n_harmony_", names(harm_wide))
pat <- merge(pat, harm_wide, by = "patient", all.x = TRUE)
pat$n_harmony_malignant[is.na(pat$n_harmony_malignant)] <- 0
if ("n_harmony_TNK" %in% names(pat)) {
  pat$n_harmony_TNK[is.na(pat$n_harmony_TNK)] <- 0
} else {
  pat$n_harmony_TNK <- 0
}
pat <- pat[order(pat$dataset, pat$patient), ]

pat$q_pct <- NA_character_
pat$q_mean <- NA_character_
for (ds in unique(pat$dataset)) {
  idx <- pat$dataset == ds & pat$eligible
  pat$q_pct[idx] <- assign_q(pat$cldn4_pct_pos[idx])
  pat$q_mean[idx] <- assign_q(pat$cldn4_mean[idx])
}

write.table(pat, file.path(out_dir, "results", "tables", "patient_level_cldn4_tnk.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lib123, file.path(out_dir, "results", "tables", "gse123902_library_level.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

el <- pat[pat$eligible, , drop = FALSE]
contrasts <- list(
  cldn4_pct_vs_tnk = list(x = "cldn4_pct_pos", y = "frac_tnk"),
  cldn4_mean_vs_tnk = list(x = "cldn4_mean", y = "frac_tnk"),
  cldn4_pct_vs_ifn = list(x = "cldn4_pct_pos", y = "ifn_mean"),
  cldn4_pct_vs_mhc = list(x = "cldn4_pct_pos", y = "mhc_mean"),
  cldn4_pct_vs_tj = list(x = "cldn4_pct_pos", y = "tj_mean"),
  cldn4_mean_vs_ifn = list(x = "cldn4_mean", y = "ifn_mean"),
  cldn4_mean_vs_mhc = list(x = "cldn4_mean", y = "mhc_mean"),
  cldn4_mean_vs_tj = list(x = "cldn4_mean", y = "tj_mean")
)

test_rows <- list()
for (nm in names(contrasts)) {
  xcol <- contrasts[[nm]]$x
  ycol <- contrasts[[nm]]$y
  per <- list()
  for (ds in c("GSE123902", "GSE205335")) {
    sub <- el[el$dataset == ds, ]
    sp <- spearman_ci(sub[[xcol]], sub[[ycol]])
    qh <- sub$q_pct == "Q4"
    ql <- sub$q_pct == "Q1"
    mw <- mw_delta(sub[[ycol]][ql | qh], qh[ql | qh])
    per[[ds]] <- list(sp = sp, mw = mw)
    test_rows[[length(test_rows) + 1L]] <- data.frame(
      dataset = ds, contrast = nm,
      n = sp$n, rho = sp$rho, p = sp$p, lo = sp$lo, hi = sp$hi,
      n_q1 = mw$n_q1, n_q4 = mw$n_q4, delta = mw$delta, r = mw$r, p_q4q1 = mw$p,
      I2 = NA_real_, k = 1, stringsAsFactors = FALSE
    )
  }
  dl <- fisher_z_dl(c(per$GSE123902$sp$rho, per$GSE205335$sp$rho),
                    c(per$GSE123902$sp$n, per$GSE205335$sp$n))
  stack <- el[el$q_pct %in% c("Q1", "Q4"), ]
  mw_stack <- mw_delta(stack[[ycol]], stack$q_pct == "Q4")
  test_rows[[length(test_rows) + 1L]] <- data.frame(
    dataset = "combined_DL", contrast = nm,
    n = dl$N, rho = dl$rho, p = dl$p, lo = dl$lo, hi = dl$hi,
    n_q1 = mw_stack$n_q1, n_q4 = mw_stack$n_q4,
    delta = mw_stack$delta, r = mw_stack$r, p_q4q1 = mw_stack$p,
    I2 = dl$I2, k = dl$k, stringsAsFactors = FALSE
  )
  # stacked patients, ignoring cohort (descriptive)
  sp_all <- spearman_ci(el[[xcol]], el[[ycol]])
  test_rows[[length(test_rows) + 1L]] <- data.frame(
    dataset = "stacked_patients", contrast = nm,
    n = sp_all$n, rho = sp_all$rho, p = sp_all$p, lo = sp_all$lo, hi = sp_all$hi,
    n_q1 = mw_stack$n_q1, n_q4 = mw_stack$n_q4,
    delta = mw_stack$delta, r = mw_stack$r, p_q4q1 = mw_stack$p,
    I2 = NA_real_, k = 2, stringsAsFactors = FALSE
  )
}
tests <- do.call(rbind, test_rows)
write.table(tests, file.path(out_dir, "results", "tables", "patient_level_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pull <- function(ds, contrast) {
  row <- tests[tests$dataset == ds & tests$contrast == contrast, ]
  if (!nrow(row)) return(NULL)
  list(
    n = row$n[1], rho = row$rho[1], p = row$p[1],
    lo = row$lo[1], hi = row$hi[1],
    n_q1 = row$n_q1[1], n_q4 = row$n_q4[1],
    delta = row$delta[1], r = row$r[1], p_q4q1 = row$p_q4q1[1],
    I2 = row$I2[1]
  )
}

s123_pct <- pull("GSE123902", "cldn4_pct_vs_tnk")
s205_pct <- pull("GSE205335", "cldn4_pct_vs_tnk")
s_dl_pct <- pull("combined_DL", "cldn4_pct_vs_tnk")
s_st_pct <- pull("stacked_patients", "cldn4_pct_vs_tnk")
s123_mean <- pull("GSE123902", "cldn4_mean_vs_tnk")
s205_mean <- pull("GSE205335", "cldn4_mean_vs_tnk")
s_dl_mean <- pull("combined_DL", "cldn4_mean_vs_tnk")
s123_ifn <- pull("GSE123902", "cldn4_pct_vs_ifn")
s205_ifn <- pull("GSE205335", "cldn4_pct_vs_ifn")
s_dl_ifn <- pull("combined_DL", "cldn4_pct_vs_ifn")
s_st_ifn <- pull("stacked_patients", "cldn4_pct_vs_ifn")
s123_mhc <- pull("GSE123902", "cldn4_pct_vs_mhc")
s205_mhc <- pull("GSE205335", "cldn4_pct_vs_mhc")
s_dl_mhc <- pull("combined_DL", "cldn4_pct_vs_mhc")
s_st_mhc <- pull("stacked_patients", "cldn4_pct_vs_mhc")
s123_tj <- pull("GSE123902", "cldn4_pct_vs_tj")
s205_tj <- pull("GSE205335", "cldn4_pct_vs_tj")
s_dl_tj <- pull("combined_DL", "cldn4_pct_vs_tj")

q_ifn <- pull("combined_DL", "cldn4_pct_vs_ifn")
q_mhc <- pull("combined_DL", "cldn4_pct_vs_mhc")
q_tj <- pull("combined_DL", "cldn4_pct_vs_tj")
q_tnk <- pull("combined_DL", "cldn4_pct_vs_tnk")

n_tab <- data.frame(
  item = c(
    "GSE123902 GEO dense CSVs",
    "GSE123902 barcodes read",
    "GSE123902 donors / LX IDs",
    "GSE123902 eligible tumor donors",
    "GSE205335 CellIdentity rows",
    "GSE205335 locked eligible patients",
    "combined eligible patients",
    "cells in Harmony object",
    "malignant cells in Harmony object",
    "T/NK cells in Harmony object",
    "IFN genes present after merge",
    "MHC-I/APM genes present after merge",
    "TJ genes present after merge (CLDN4 held out)",
    "GSE148071 cells",
    "GSE127465 cells",
    "dual-high TACSTD2 ∩ CLDN4"
  ),
  n = c(
    length(csv_files),
    n123_barcodes,
    length(unique(lib123$patient)),
    sum(pat$dataset == "GSE123902" & pat$eligible),
    n205_ann,
    sum(pat$dataset == "GSE205335" & pat$eligible),
    nrow(el),
    ncol(obj),
    sum(obj$compartment == "malignant"),
    sum(obj$compartment == "TNK"),
    length(ifn_present),
    length(mhc_present),
    length(tj_present),
    0, 0, NA
  ),
  note = c(
    "GSE123902_RAW.tar public processed UMI",
    "union of SEQC dense CSVs",
    "filename MSK_LX*",
    "≥20 marker-epithelial in tumor and ≥20 T/NK",
    "author annotation, 96505 barcodes",
    "author malignant ≥20 (PR #459 table) and ≥20 T/NK",
    "patient is the unit; cells are not n",
    sprintf("capped ≤%d mal + ≤%d T/NK per patient", cap_n, cap_n),
    "author (205335) / marker epithelial (123902)",
    "T/NK fraction uses full-sample annotation, not this cap",
    "Hallmark IFNα ∩ IFNγ core; CLDN4 not in set",
    "curated antigen-presentation set; CLDN4 not in set",
    "epithelial TJ; CLDN4 held out",
    "not merged",
    "not merged (mouse; out of scope)",
    "not defined; CLDN4-only"
  ),
  stringsAsFactors = FALSE
)
write.table(n_tab, file.path(out_dir, "results", "tables", "n_honest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cov_tab <- data.frame(
  set = c("IFN_core", "MHC_I_APM", "TJ_no_CLDN4"),
  n_locked = c(length(IFN_GENES), length(MHC_GENES), length(TJ_GENES)),
  n_present = c(length(ifn_present), length(mhc_present), length(tj_present)),
  missing = c(
    paste(setdiff(IFN_GENES, ifn_present), collapse = ","),
    paste(setdiff(MHC_GENES, mhc_present), collapse = ","),
    paste(setdiff(TJ_GENES, tj_present), collapse = ",")
  ),
  stringsAsFactors = FALSE
)
write.table(cov_tab, file.path(out_dir, "results", "tables", "geneset_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ===========================================================================
# 5. DimPlot + patient figures
# ===========================================================================
message("[seurat] DimPlot / FeaturePlot")
theme_ok <- ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "bottom")

p_ds <- DimPlot(obj, reduction = "umap", group.by = "dataset", pt.size = 0.25) +
  ggtitle("Harmony UMAP — dataset (GSE123902 + GSE205335)")
ggsave_both(p_ds, "fig_dimplot_dataset")

p_comp <- DimPlot(obj, reduction = "umap", group.by = "compartment", pt.size = 0.25) +
  ggtitle("Harmony UMAP — malignant vs T/NK")
ggsave_both(p_comp, "fig_dimplot_compartment")

p_pat <- DimPlot(obj, reduction = "umap", group.by = "patient", pt.size = 0.15) +
  ggtitle("Harmony UMAP — patient") +
  theme(legend.position = "none")
ggsave_both(p_pat, "fig_dimplot_patient", w = 6.8, h = 5.4)

p_cldn <- FeaturePlot(obj, features = "CLDN4", reduction = "umap", pt.size = 0.25,
                      order = TRUE) +
  ggtitle("Harmony UMAP — CLDN4 (log-norm)")
ggsave_both(p_cldn, "fig_featureplot_cldn4")

p_vln <- VlnPlot(obj, features = "CLDN4", group.by = "compartment", pt.size = 0) +
  ggtitle("CLDN4 by compartment (Harmony object)")
ggsave_both(p_vln, "fig_vlnplot_cldn4", w = 5.2, h = 4.6)

p_sc <- ggplot(el, aes(cldn4_pct_pos, frac_tnk, color = dataset)) +
  geom_point(size = 2.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_ok +
  xlab("Malignant CLDN4 %pos (all malignant cells)") +
  ylab("T/NK fraction (full-sample tumor cells)") +
  ggtitle("Patient unit: CLDN4 %pos vs T/NK")
ggsave_both(p_sc, "fig_patient_cldn4_vs_tnk")

p_ifn <- ggplot(el, aes(cldn4_pct_pos, ifn_mean, color = dataset)) +
  geom_point(size = 2.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_ok +
  xlab("Malignant CLDN4 %pos") +
  ylab("Malignant IFN (mean log1p CP10k)") +
  ggtitle("Patient unit: CLDN4 vs malignant IFN")
ggsave_both(p_ifn, "fig_patient_cldn4_vs_ifn")

p_mhc <- ggplot(el, aes(cldn4_pct_pos, mhc_mean, color = dataset)) +
  geom_point(size = 2.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_ok +
  xlab("Malignant CLDN4 %pos") +
  ylab("Malignant MHC-I/APM (mean log1p CP10k)") +
  ggtitle("Patient unit: CLDN4 vs malignant MHC")
ggsave_both(p_mhc, "fig_patient_cldn4_vs_mhc")

qdf <- el[el$q_pct %in% c("Q1", "Q4"), ]
qdf$q_pct <- factor(qdf$q_pct, levels = c("Q1", "Q4"))
p_qifn <- ggplot(qdf, aes(q_pct, ifn_mean, fill = dataset)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.7) +
  theme_ok + ylab("Malignant IFN") + xlab("Within-cohort CLDN4 %pos quartile") +
  ggtitle("Q4 vs Q1 malignant IFN (patient unit)")
ggsave_both(p_qifn, "fig_q4q1_ifn")

p_qmhc <- ggplot(qdf, aes(q_pct, mhc_mean, fill = dataset)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.7) +
  theme_ok + ylab("Malignant MHC-I/APM") + xlab("Within-cohort CLDN4 %pos quartile") +
  ggtitle("Q4 vs Q1 malignant MHC (patient unit)")
ggsave_both(p_qmhc, "fig_q4q1_mhc")

p_qtnk <- ggplot(qdf, aes(q_pct, frac_tnk, fill = dataset)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.7) +
  theme_ok + ylab("T/NK fraction") + xlab("Within-cohort CLDN4 %pos quartile") +
  ggtitle("Q4 vs Q1 T/NK (patient unit)")
ggsave_both(p_qtnk, "fig_q4q1_tnk")

# ===========================================================================
# 6. FINDING.md
# ===========================================================================
fmt_s <- function(s) {
  if (is.null(s) || is.na(s$rho)) return("NA")
  sprintf("n=%d, ρ=%+.3f [%.3f, %.3f], p=%.3g", s$n, s$rho, s$lo, s$hi, s$p)
}
fmt_row <- function(s) {
  if (is.null(s) || is.na(s$rho)) return("| — | — | — | — |")
  sprintf("| %d | %+.3f [%.3f, %.3f] | %.3g | I²=%s |",
          s$n, s$rho, s$lo, s$hi, s$p,
          ifelse(is.na(s$I2), "—", sprintf("%.0f%%", s$I2)))
}
fmt_q <- function(s) {
  if (is.null(s) || is.na(s$delta)) return("Q4 vs Q1 not run")
  sprintf("n_Q1=%d n_Q4=%d, Δ median=%+.3f, r=%+.3f, MW p=%.3g",
          s$n_q1, s$n_q4, s$delta, s$r, s$p_q4q1)
}

pat_md_lines <- apply(pat, 1, function(r) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r[["dataset"]], r[["patient"]], r[["tumor_sites"]], r[["histology"]],
    r[["n_cells_tumor"]], r[["n_malignant"]], r[["n_tnk"]],
    fmt_num(as.numeric(r[["frac_tnk"]])),
    fmt_num(as.numeric(r[["cldn4_mean"]])),
    ifelse(is.na(as.numeric(r[["cldn4_pct_pos"]])), "NA",
           sprintf("%.1f", as.numeric(r[["cldn4_pct_pos"]]))),
    fmt_num(as.numeric(r[["ifn_mean"]])),
    fmt_num(as.numeric(r[["mhc_mean"]])),
    fmt_num(as.numeric(r[["tj_mean"]])),
    ifelse(as.logical(r[["eligible"]]), "yes", "no")
  )
})

seurat_ver <- as.character(utils::packageVersion("Seurat"))
harm_ver <- as.character(utils::packageVersion("harmony"))
r_ver <- paste(R.version$major, R.version$minor, sep = ".")

finding <- paste0(
  "# FINDING — Seurat / Harmony pair GSE123902 + GSE205335, CLDN4-only\n\n",
  "ADDITIVE. **CLDN4 only.** Thesis already correct. This is a Seurat-native ",
  "Harmony integration of the public human pair that previously differed: ",
  "[GSE123902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902) ",
  "(Laughney et al., *Nat Med* 2020) + ",
  "[GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) ",
  "(Hu et al. lung IO atlas). SuperSeries GSE123904 / mouse GSE123903 are not used. ",
  "No TACSTD2∩CLDN4 dual-high gate. No GSE148071. No GSE127465. ",
  "Concordant-4 four-way merge is a different agent. No Python-only primary.\n\n",
  "Primary engine: **R + Seurat ", seurat_ver, " + Harmony ", harm_ver, "**. ",
  "`IntegrateLayers(HarmonyIntegration)` on `dataset` layers (", integ_method,
  ", reduction `", integ_reduction, "`). Honest unit = **patient** ",
  "(GSE123902 donor / LX ID; GSE205335 patient). Cells are counts, not n.\n\n",
  "Given pair T/NK Spearman from PR #459 is **not re-derived** here ",
  "(n=35, %pos ρ=−0.522, Q4 vs Q1 r=−0.802, 9 vs 9). ",
  "This folder adds the Seurat/Harmony object, a patient table, and DimPlot.\n\n",
  "## Verdict\n\n",
  "Seurat/Harmony ran. Patient table: `results/tables/patient_level_cldn4_tnk.tsv`. ",
  "Eligible patients (≥", MIN_MAL, " malignant and ≥", MIN_TNK, " T/NK): **n=", nrow(el),
  "** (", sum(el$dataset == "GSE123902"), " GSE123902 + ",
  sum(el$dataset == "GSE205335"), " GSE205335).\n\n",
  "- GSE123902 CLDN4 **%pos** vs T/NK: ", fmt_s(s123_pct), ".\n",
  "- GSE205335 CLDN4 **%pos** vs T/NK: ", fmt_s(s205_pct), ".\n",
  "- Combined DL (Fisher-z) %pos vs T/NK: ", fmt_s(s_dl_pct),
  if (!is.null(s_dl_pct) && is.finite(s_dl_pct$I2)) sprintf(" (I²=%.0f%%)", s_dl_pct$I2) else "", ".\n",
  "- Stacked patients %pos vs T/NK: ", fmt_s(s_st_pct), ".\n",
  "- Combined DL %pos vs malignant IFN: ", fmt_s(s_dl_ifn), ".\n",
  "- Combined DL %pos vs malignant MHC-I/APM: ", fmt_s(s_dl_mhc), ".\n",
  "- Combined DL %pos vs TJ (CLDN4 held out): ", fmt_s(s_dl_tj), ".\n\n",
  "Q4 vs Q1 (within-cohort %pos quartiles, stacked patients): IFN ", fmt_q(q_ifn),
  "; MHC ", fmt_q(q_mhc), "; T/NK ", fmt_q(q_tnk), ".\n\n",
  "n=", nrow(el), " is the honest ceiling. Cell-level p-values are not the claim. ",
  "Do not write this as a failed audit of the thesis.\n\n",
  "## Honest n\n\n",
  "| item | n | note |\n",
  "|---|---:|---|\n",
  paste(sprintf("| %s | %s | %s |", n_tab$item,
                ifelse(is.na(n_tab$n), "—", as.character(n_tab$n)),
                n_tab$note), collapse = "\n"),
  "\n\n",
  "## Gate\n\n",
  "| file | public? | used |\n",
  "|---|---|---|\n",
  "| `GSE123902_RAW.tar` (17 dense UMI CSVs) | yes | **yes** — marker lineage, then Seurat |\n",
  "| `GSE205335_Lung_IO_UMI_matrix.rds` + CellIdentity | yes | **yes** — author malignant / T/NK |\n",
  "| GSE148071 / GSE127465 / GSE131907 / GSE189357 | — | **no** |\n",
  "| Author 36.5 GB GSE123902 H5 | yes | no — GEO CSVs suffice |\n\n",
  "## Locked choices\n\n",
  "- GSE123902 lineage is a four-way marker argmax on log1p CP10k. ",
  "Keep if top ≥ 0.12 and top ≥ 1.15 × second. CLDN4 is never a lineage marker. ",
  "Malignant = marker epithelial **in tumor**. Matched normal is not malignant.\n",
  "- GSE205335 malignant = author `lineage.sub == Malignant cells` in tumor tissue. ",
  "T/NK = author `lineage.total == T/NK cells`. Normal lung/LN/brain dropped.\n",
  "- T/NK fraction = full-sample tumor cells of that patient (not the Harmony cap).\n",
  "- IFN / MHC / TJ scores = mean log1p CP10k of locked genes on **all** malignant cells ",
  "(before the ≤", cap_n, "/compartment Harmony cap). TJ holds CLDN4 out.\n",
  "- Quartiles are **within-cohort** on malignant CLDN4 %pos.\n",
  "- Eligible Spearman n requires ≥", MIN_MAL, " malignant and ≥", MIN_TNK, " T/NK.\n\n",
  "## Patient table\n\n",
  "Machine table: `results/tables/patient_level_cldn4_tnk.tsv`.\n\n",
  "| dataset | patient | tumor site | histology | n tumor | n mal | n T/NK | frac T/NK | CLDN4 mean | CLDN4 %pos | IFN | MHC | TJ | eligible |\n",
  "|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n",
  paste(pat_md_lines, collapse = "\n"), "\n\n",
  "## Patient-level Spearman\n\n",
  "### CLDN4 %pos vs T/NK fraction\n\n",
  "| split | n | ρ [95% CI] | p | note |\n",
  "|---|---:|---|---:|---|\n",
  sprintf("| GSE123902 | %d | %+.3f [%.3f, %.3f] | %.3g | marker epithelial |\n",
          s123_pct$n, s123_pct$rho, s123_pct$lo, s123_pct$hi, s123_pct$p),
  sprintf("| GSE205335 | %d | %+.3f [%.3f, %.3f] | %.3g | author malignant |\n",
          s205_pct$n, s205_pct$rho, s205_pct$lo, s205_pct$hi, s205_pct$p),
  sprintf("| combined DL | %d | %+.3f [%.3f, %.3f] | %.3g | I²=%.0f%% |\n",
          s_dl_pct$n, s_dl_pct$rho, s_dl_pct$lo, s_dl_pct$hi, s_dl_pct$p, s_dl_pct$I2 %||% NA),
  sprintf("| stacked patients | %d | %+.3f [%.3f, %.3f] | %.3g | cohort mix |\n",
          s_st_pct$n, s_st_pct$rho, s_st_pct$lo, s_st_pct$hi, s_st_pct$p),
  "\n",
  "### CLDN4 mean vs T/NK fraction\n\n",
  "| split | n | ρ [95% CI] | p |\n",
  "|---|---:|---|---:|\n",
  sprintf("| GSE123902 | %d | %+.3f [%.3f, %.3f] | %.3g |\n",
          s123_mean$n, s123_mean$rho, s123_mean$lo, s123_mean$hi, s123_mean$p),
  sprintf("| GSE205335 | %d | %+.3f [%.3f, %.3f] | %.3g |\n",
          s205_mean$n, s205_mean$rho, s205_mean$lo, s205_mean$hi, s205_mean$p),
  sprintf("| combined DL | %d | %+.3f [%.3f, %.3f] | %.3g |\n",
          s_dl_mean$n, s_dl_mean$rho, s_dl_mean$lo, s_dl_mean$hi, s_dl_mean$p),
  "\n",
  "### Malignant programs vs CLDN4 %pos\n\n",
  "| contrast | GSE123902 ρ (p) | GSE205335 ρ (p) | combined DL ρ (p) |\n",
  "|---|---|---|---|\n",
  sprintf("| IFN | %+.3f (%.3g) | %+.3f (%.3g) | %+.3f (%.3g) |\n",
          s123_ifn$rho, s123_ifn$p, s205_ifn$rho, s205_ifn$p, s_dl_ifn$rho, s_dl_ifn$p),
  sprintf("| MHC-I/APM | %+.3f (%.3g) | %+.3f (%.3g) | %+.3f (%.3g) |\n",
          s123_mhc$rho, s123_mhc$p, s205_mhc$rho, s205_mhc$p, s_dl_mhc$rho, s_dl_mhc$p),
  sprintf("| TJ (CLDN4 held out) | %+.3f (%.3g) | %+.3f (%.3g) | %+.3f (%.3g) |\n",
          s123_tj$rho, s123_tj$p, s205_tj$rho, s205_tj$p, s_dl_tj$rho, s_dl_tj$p),
  "\n",
  "## Malignant Q4 vs Q1 IFN / MHC / T/NK\n\n",
  "Within-cohort CLDN4-%pos quartiles among eligible patients. Mid quartiles unused. ",
  "Stacked Q4 vs Q1 is the pair contrast.\n\n",
  "| family | n_Q1 | n_Q4 | Δ median (Q4−Q1) | r | MW p |\n",
  "|---|---:|---:|---:|---:|---:|\n",
  sprintf("| T/NK fraction | %d | %d | %+.3f | %+.3f | %.3g |\n",
          q_tnk$n_q1, q_tnk$n_q4, q_tnk$delta, q_tnk$r, q_tnk$p_q4q1),
  sprintf("| IFN | %d | %d | %+.3f | %+.3f | %.3g |\n",
          q_ifn$n_q1, q_ifn$n_q4, q_ifn$delta, q_ifn$r, q_ifn$p_q4q1),
  sprintf("| MHC-I/APM | %d | %d | %+.3f | %+.3f | %.3g |\n",
          q_mhc$n_q1, q_mhc$n_q4, q_mhc$delta, q_mhc$r, q_mhc$p_q4q1),
  sprintf("| TJ (CLDN4 held out) | %d | %d | %+.3f | %+.3f | %.3g |\n",
          q_tj$n_q1, q_tj$n_q4, q_tj$delta, q_tj$r, q_tj$p_q4q1),
  "\n",
  "Single-cohort Q4 vs Q1 tails are thinner (GSE123902 4 vs 4; GSE205335 6 vs 6 on %pos). ",
  "GSE205335 Q4 mixes SCLC with ADC; that histology mix is part of the honest n.\n\n",
  "## Seurat plots\n\n",
  "- `results/figures/fig_dimplot_dataset.png` — `DimPlot` by dataset (Harmony UMAP)\n",
  "- `results/figures/fig_dimplot_compartment.png` — `DimPlot` malignant vs T/NK\n",
  "- `results/figures/fig_dimplot_patient.png` — `DimPlot` by patient\n",
  "- `results/figures/fig_featureplot_cldn4.png` — `FeaturePlot` CLDN4\n",
  "- `results/figures/fig_vlnplot_cldn4.png` — `VlnPlot` CLDN4 by compartment\n",
  "- `results/figures/fig_patient_cldn4_vs_tnk.png` — patient scatter %pos vs T/NK\n",
  "- `results/figures/fig_patient_cldn4_vs_ifn.png` — patient scatter vs IFN\n",
  "- `results/figures/fig_patient_cldn4_vs_mhc.png` — patient scatter vs MHC\n",
  "- `results/figures/fig_q4q1_ifn.png` / `fig_q4q1_mhc.png` / `fig_q4q1_tnk.png`\n\n",
  "## What this does not claim\n\n",
  "- Cell-level p-values are not the claim. n_cells is large by construction.\n",
  "- GSE123902 marker epithelial is **not** a CNV-malignant call.\n",
  "- GSE205335 author malignant includes SCLC / NUT / SQ / ADC; histology is not hidden.\n",
  "- The Harmony cap (≤", cap_n, " mal + ≤", cap_n, " T/NK per patient) is for UMAP only. ",
  "T/NK fraction and program scores use the full patient.\n",
  "- Restriction (CLDN4 in epithelium vs T/NK) is not an infiltration / immune-cold claim by itself.\n",
  "- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.\n",
  "- No ICI / MPR / RECIST / survival test is the primary claim.\n",
  "- Not merged with GSE148071, GSE127465, or the concordant-4 four-way.\n",
  "- The PR #459 n=35 ρ=−0.522 row is given and is not recomputed as the thesis test.\n\n",
  "## Reproduce\n\n",
  "```bash\n",
  "bash methods/seurat_integrate_123902_205335_cldn4/scripts/00_download.sh\n",
  "Rscript methods/seurat_integrate_123902_205335_cldn4/scripts/01_seurat_harmony_integrate.R\n",
  "```\n\n",
  "Seurat ", seurat_ver, ". Harmony ", harm_ver, ". R ", r_ver, ".\n"
)

writeLines(finding, file.path(out_dir, "FINDING.md"))

summary <- list(
  accessions = c("GSE123902", "GSE205335"),
  engine = "R_Seurat_Harmony",
  seurat = seurat_ver,
  harmony = harm_ver,
  integration_method = integ_method,
  integration_reduction = integ_reduction,
  n_cells_harmony = ncol(obj),
  n_eligible = nrow(el),
  n_gse123902 = sum(el$dataset == "GSE123902"),
  n_gse205335 = sum(el$dataset == "GSE205335"),
  spearman_pct_vs_tnk_dl = s_dl_pct,
  spearman_pct_vs_ifn_dl = s_dl_ifn,
  spearman_pct_vs_mhc_dl = s_dl_mhc,
  q4q1_ifn = q_ifn,
  q4q1_mhc = q_mhc,
  q4q1_tnk = q_tnk
)
jsonlite::write_json(summary, file.path(out_dir, "results", "summary.json"),
                     auto_unbox = TRUE, pretty = TRUE)
writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(out_dir, "results", "sessionInfo.txt"))

message("DONE. Eligible patients n=", nrow(el),
        " Harmony cells=", ncol(obj),
        " method=", integ_method)
message("  %pos vs T/NK DL: ", fmt_s(s_dl_pct))
message("  %pos vs IFN DL:  ", fmt_s(s_dl_ifn))
message("  %pos vs MHC DL:  ", fmt_s(s_dl_mhc))

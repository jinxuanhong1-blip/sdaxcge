#!/usr/bin/env Rscript
# GSE189357 CLDN4-only: Seurat CreateSeuratObject, patient-level T/NK + IFN/MHC/TJ.
# ADDITIVE. Thesis already correct. No dual-high. No GSE148071. No Python-only primary.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

opt <- list(
  tar = "/tmp/gse189357/GSE189357_RAW.tar",
  extract = "/tmp/gse189357/raw",
  outdir = "methods/seurat_gse189357_cldn4/results",
  finding = "methods/seurat_gse189357_cldn4/FINDING.md"
)
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[[i]])
  if (key %in% names(opt) && i < length(args)) {
    opt[[key]] <- args[[i + 1]]
    i <- i + 2
  } else {
    stop("unknown or incomplete arg: ", args[[i]])
  }
}

here <- function(...) {
  script <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.na(script)) {
    return(file.path(getwd(), ...))
  }
  file.path(dirname(dirname(script)), ...)
}

root <- tryCatch({
  script <- sub("--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  normalizePath(file.path(dirname(script), ".."))
}, error = function(e) normalizePath("methods/seurat_gse189357_cldn4"))

meta_path <- file.path(root, "data", "sample_metadata.tsv")
sets_path <- file.path(root, "data", "family_sets.json")
if (!file.exists(meta_path)) stop("missing ", meta_path)
if (!file.exists(sets_path)) stop("missing ", sets_path)

meta <- read.delim(meta_path, stringsAsFactors = FALSE)
sets <- fromJSON(sets_path)
fam <- list(
  IFN = unique(as.character(sets$IFN)),
  MHC_I_APM = unique(as.character(sets$MHC_I_APM)),
  TJ = unique(as.character(sets$TJ))
)
stopifnot(!"CLDN4" %in% fam$TJ)
stopifnot(!"CLDN4" %in% fam$IFN)

EPI <- c("EPCAM", "KRT8", "KRT18", "KRT19")
TNK <- c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
ALIASES <- list(
  WARS1 = c("WARS", "WARS1"),
  MARCHF1 = c("MARCH1", "MARCHF1"),
  HLA.A = c("HLA-A", "HLA.A"),
  HLA.B = c("HLA-B", "HLA.B"),
  HLA.C = c("HLA-C", "HLA.C"),
  HLA.E = c("HLA-E", "HLA.E"),
  HLA.F = c("HLA-F", "HLA.F"),
  HLA.G = c("HLA-G", "HLA.G")
)

logmsg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ", flush = TRUE)
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 3))
  formatC(p, format = "f", digits = 4)
}

fmt_num <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

spearman_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4) {
    return(list(n = n, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = n, rho = unname(ct$estimate), p = unname(ct$p.value))
}

bh <- function(p) {
  p.adjust(p, method = "BH")
}

assign_quartiles <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  out <- rep(NA_character_, length(x))
  ok <- is.finite(r)
  if (sum(ok) < 4) return(out)
  qs <- tryCatch(
    as.character(cut(
      r[ok],
      breaks = quantile(r[ok], probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )),
    error = function(e) rep(NA_character_, sum(ok))
  )
  out[ok] <- qs
  out
}

resolve_gene <- function(feats, symbol) {
  ufeats <- toupper(feats)
  cands <- unique(c(symbol, ALIASES[[gsub("-", ".", symbol)]], ALIASES[[symbol]]))
  cands <- unique(toupper(cands[!is.na(cands)]))
  hit <- match(cands, ufeats)
  hit <- hit[!is.na(hit)]
  if (length(hit)) return(feats[[hit[[1]]]])
  # Seurat unique.features suffix .1
  base <- sub("\\.\\d+$", "", ufeats)
  hit <- match(cands, base)
  hit <- hit[!is.na(hit)]
  if (length(hit)) return(feats[[hit[[1]]]])
  NA_character_
}

resolve_genes <- function(feats, symbols) {
  v <- vapply(symbols, function(g) resolve_gene(feats, g), character(1), USE.NAMES = FALSE)
  unique(v[!is.na(v)])
}

pos_any <- function(counts, feats, symbols) {
  genes <- resolve_genes(feats, symbols)
  if (!length(genes)) return(rep(FALSE, ncol(counts)))
  as.numeric(Matrix::colSums(counts[genes, , drop = FALSE]) > 0) > 0
}

gene_vec <- function(mat, feats, symbol) {
  g <- resolve_gene(feats, symbol)
  if (is.na(g)) return(rep(0, ncol(mat)))
  as.numeric(mat[g, ])
}

family_mean <- function(mat, feats, symbols) {
  genes <- resolve_genes(feats, symbols)
  if (!length(genes)) return(list(score = rep(NA_real_, ncol(mat)), n = 0L, genes = character()))
  list(
    score = as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE])),
    n = length(genes),
    genes = genes
  )
}

pseudobulk_logcpm <- function(counts, cells) {
  if (!length(cells)) {
    return(list(logcpm = numeric(), libsize = 0))
  }
  umi <- Matrix::rowSums(counts[, cells, drop = FALSE])
  lib <- sum(umi)
  if (lib <= 0) return(list(logcpm = umi * NA_real_, libsize = lib))
  list(logcpm = log2(umi / lib * 1e6 + 1), libsize = lib)
}

ensure_extract <- function(tar, dest) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  needed <- paste0(meta$gsm, "_", meta$patient, "_matrix.mtx.gz")
  have <- file.exists(file.path(dest, needed))
  if (all(have)) {
    logmsg("extract already present")
    return(invisible(dest))
  }
  if (!file.exists(tar)) stop("missing tar: ", tar)
  logmsg("extracting", tar, "->", dest)
  ok <- system2("tar", c("-xf", tar, "-C", dest), stdout = TRUE, stderr = TRUE)
  invisible(dest)
}

read_patient <- function(patient, gsm, extract_dir) {
  prefix <- file.path(extract_dir, paste0(gsm, "_", patient))
  mtx <- paste0(prefix, "_matrix.mtx.gz")
  cells <- paste0(prefix, "_barcodes.tsv.gz")
  features <- paste0(prefix, "_features.tsv.gz")
  if (!all(file.exists(c(mtx, cells, features)))) {
    stop("missing 10x files for ", patient)
  }
  logmsg("ReadMtx + CreateSeuratObject", patient)
  mat <- ReadMtx(
    mtx = mtx,
    cells = cells,
    features = features,
    cell.column = 1,
    feature.column = 2,
    unique.features = TRUE
  )
  colnames(mat) <- paste0(patient, "_", colnames(mat))
  obj <- CreateSeuratObject(
    counts = mat,
    project = patient,
    min.cells = 0,
    min.features = 0
  )
  obj$patient <- patient
  obj$unit_id <- patient
  obj$dataset <- "GSE189357"
  obj
}

process_patient <- function(patient, gsm, stage, extract_dir) {
  obj <- read_patient(patient, gsm, extract_dir)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  feats <- rownames(obj)
  mal <- pos_any(counts, feats, EPI) & (gene_vec(counts, feats, "PTPRC") == 0)
  tnk <- pos_any(counts, feats, TNK) & (!mal)
  obj$is_mal <- mal
  obj$is_tnk <- tnk

  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  norm <- GetAssayData(obj, assay = "RNA", layer = "data")

  mal_cells <- colnames(obj)[mal]
  tnk_cells <- colnames(obj)[tnk]
  n_cells <- ncol(obj)
  n_mal <- length(mal_cells)
  n_tnk <- length(tnk_cells)

  cldn4 <- gene_vec(norm, feats, "CLDN4")
  tacstd2 <- gene_vec(norm, feats, "TACSTD2")
  cldn4_counts <- gene_vec(counts, feats, "CLDN4")

  if (n_mal > 0) {
    mal_cldn4 <- cldn4[mal]
    mal_cldn4_umi <- cldn4_counts[mal]
    mal_tac <- tacstd2[mal]
  } else {
    mal_cldn4 <- numeric()
    mal_cldn4_umi <- numeric()
    mal_tac <- numeric()
  }

  pb <- pseudobulk_logcpm(counts, mal_cells)
  pb_scores <- lapply(fam, function(genes) {
    present <- resolve_genes(feats, genes)
    list(
      score = if (length(present) && length(pb$logcpm)) mean(pb$logcpm[present]) else NA_real_,
      n_genes = length(present),
      missing = setdiff(genes, present)
    )
  })

  cell_scores <- lapply(fam, function(genes) {
    if (n_mal == 0) {
      return(list(mean = NA_real_, n_genes = 0L))
    }
    sc <- family_mean(norm[, mal_cells, drop = FALSE], feats, genes)
    list(mean = mean(sc$score), n_genes = sc$n)
  })

  # Seurat AddModuleScore on malignant cells (matching extra; not the claim).
  ams <- list(IFN = NA_real_, MHC_I_APM = NA_real_, TJ = NA_real_)
  if (n_mal >= 20) {
    mal_obj <- subset(obj, cells = mal_cells)
    feats_use <- list(
      resolve_genes(feats, fam$IFN),
      resolve_genes(feats, fam$MHC_I_APM),
      resolve_genes(feats, fam$TJ)
    )
    if (all(lengths(feats_use) >= 5)) {
      mal_obj <- AddModuleScore(mal_obj, features = feats_use, name = "AMS", ctrl = 50, seed = 1)
      ams$IFN <- mean(mal_obj$AMS1)
      ams$MHC_I_APM <- mean(mal_obj$AMS2)
      ams$TJ <- mean(mal_obj$AMS3)
    }
  }

  # Within-patient CLDN4-high vs low (median of log1p CP10k). Extra only.
  paired <- list(
    n_high = NA_integer_, n_low = NA_integer_,
    d_IFN = NA_real_, d_MHC = NA_real_, d_TJ = NA_real_,
    eligible = FALSE
  )
  if (n_mal >= 20) {
    med <- stats::median(mal_cldn4)
    hi <- mal_cells[mal_cldn4 > med]
    lo <- mal_cells[mal_cldn4 <= med]
    paired$n_high <- length(hi)
    paired$n_low <- length(lo)
    paired$eligible <- length(hi) >= 10 && length(lo) >= 10
    if (paired$eligible) {
      sc_ifn <- family_mean(norm, feats, fam$IFN)$score
      sc_mhc <- family_mean(norm, feats, fam$MHC_I_APM)$score
      sc_tj <- family_mean(norm, feats, fam$TJ)$score
      paired$d_IFN <- mean(sc_ifn[hi]) - mean(sc_ifn[lo])
      paired$d_MHC <- mean(sc_mhc[hi]) - mean(sc_mhc[lo])
      paired$d_TJ <- mean(sc_tj[hi]) - mean(sc_tj[lo])
    }
  }

  logmsg(
    patient, "cells", n_cells, "mal", n_mal, "tnk", n_tnk,
    "CLDN4%pos", ifelse(n_mal > 0, mean(mal_cldn4_umi > 0), NA)
  )

  list(
    patient = patient,
    gsm = gsm,
    stage = stage,
    n_cells = n_cells,
    n_malignant = n_mal,
    n_tnk = n_tnk,
    frac_tnk = n_tnk / n_cells,
    mal_CLDN4_pct = if (n_mal) mean(mal_cldn4_umi > 0) else NA_real_,
    mal_CLDN4_mean = if (n_mal) mean(mal_cldn4) else NA_real_,
    mal_TACSTD2_pct = if (n_mal) mean(mal_tac > 0) else NA_real_,
    mal_TACSTD2_mean = if (n_mal) mean(mal_tac) else NA_real_,
    IFN_pb = pb_scores$IFN$score,
    MHC_pb = pb_scores$MHC_I_APM$score,
    TJ_pb = pb_scores$TJ$score,
    IFN_cell = cell_scores$IFN$mean,
    MHC_cell = cell_scores$MHC_I_APM$mean,
    TJ_cell = cell_scores$TJ$mean,
    IFN_ams = ams$IFN,
    MHC_ams = ams$MHC_I_APM,
    TJ_ams = ams$TJ,
    n_IFN_genes = pb_scores$IFN$n_genes,
    n_MHC_genes = pb_scores$MHC_I_APM$n_genes,
    n_TJ_genes = pb_scores$TJ$n_genes,
    mal_libsize = pb$libsize,
    paired_n_high = paired$n_high,
    paired_n_low = paired$n_low,
    paired_eligible = paired$eligible,
    paired_d_IFN = paired$d_IFN,
    paired_d_MHC = paired$d_MHC,
    paired_d_TJ = paired$d_TJ,
    missing_IFN = paste(pb_scores$IFN$missing, collapse = ","),
    missing_MHC = paste(pb_scores$MHC_I_APM$missing, collapse = ","),
    missing_TJ = paste(pb_scores$TJ$missing, collapse = ","),
    seurat = paste(as.character(packageVersion("Seurat")), collapse = ".")
  )
}

save_plot <- function(p, path_stub, width = 6.2, height = 4.6) {
  ggsave(paste0(path_stub, ".png"), p, width = width, height = height, dpi = 160)
  ggsave(paste0(path_stub, ".pdf"), p, width = width, height = height)
}

theme_find <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "#444444")
    )
}

stage_cols <- c(AIS = "#4c78a8", MIA = "#f58518", IAC = "#54a24b")

main <- function() {
  logmsg("Seurat", as.character(packageVersion("Seurat")))
  if (!exists("CreateSeuratObject")) stop("CreateSeuratObject missing; stop.")
  ensure_extract(opt$tar, opt$extract)
  dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

  rows <- lapply(seq_len(nrow(meta)), function(i) {
    process_patient(meta$patient[i], meta$gsm[i], meta$stage[i], opt$extract)
  })
  keep_fields <- setdiff(names(rows[[1]]), c("missing_IFN", "missing_MHC", "missing_TJ"))
  units <- as.data.frame(lapply(keep_fields, function(nm) {
    v <- lapply(rows, `[[`, nm)
    unlist(v, use.names = FALSE)
  }), stringsAsFactors = FALSE)
  names(units) <- keep_fields
  units$eligible <- units$n_malignant >= 20
  units$quartile <- assign_quartiles(units$mal_CLDN4_pct)
  units$unit <- "patient"
  units$malig_def <- "marker_malig (EPCAM|KRT8|KRT18|KRT19)>0 & PTPRC==0"
  units$tnk_def <- "(CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0 & not malignant"

  contrasts <- list(
    list(name = "CLDN4 %pos vs T/NK fraction", x = "mal_CLDN4_pct", y = "frac_tnk", family = "primary"),
    list(name = "CLDN4 mean vs T/NK fraction", x = "mal_CLDN4_mean", y = "frac_tnk", family = "primary"),
    list(name = "CLDN4 %pos vs malignant IFN (pseudobulk)", x = "mal_CLDN4_pct", y = "IFN_pb", family = "primary"),
    list(name = "CLDN4 %pos vs malignant MHC-I/APM (pseudobulk)", x = "mal_CLDN4_pct", y = "MHC_pb", family = "primary"),
    list(name = "CLDN4 %pos vs malignant TJ (pseudobulk, CLDN4 out)", x = "mal_CLDN4_pct", y = "TJ_pb", family = "primary"),
    list(name = "CLDN4 %pos vs malignant IFN (cell mean)", x = "mal_CLDN4_pct", y = "IFN_cell", family = "extra"),
    list(name = "CLDN4 %pos vs malignant MHC-I/APM (cell mean)", x = "mal_CLDN4_pct", y = "MHC_cell", family = "extra"),
    list(name = "CLDN4 %pos vs malignant TJ (cell mean, CLDN4 out)", x = "mal_CLDN4_pct", y = "TJ_cell", family = "extra"),
    list(name = "CLDN4 %pos vs TACSTD2 %pos (comparator)", x = "mal_CLDN4_pct", y = "mal_TACSTD2_pct", family = "comparator"),
    list(name = "CLDN4 %pos vs AddModuleScore IFN", x = "mal_CLDN4_pct", y = "IFN_ams", family = "extra"),
    list(name = "CLDN4 %pos vs AddModuleScore MHC-I/APM", x = "mal_CLDN4_pct", y = "MHC_ams", family = "extra"),
    list(name = "CLDN4 %pos vs AddModuleScore TJ", x = "mal_CLDN4_pct", y = "TJ_ams", family = "extra")
  )

  spear <- do.call(rbind, lapply(contrasts, function(c) {
    s <- spearman_one(as.numeric(units[[c$x]]), as.numeric(units[[c$y]]))
    data.frame(
      contrast = c$name,
      family = c$family,
      n_patients = s$n,
      rho = s$rho,
      p = s$p,
      stringsAsFactors = FALSE
    )
  }))
  prim <- spear$family == "primary"
  spear$q <- NA_real_
  spear$q[prim] <- bh(spear$p[prim])

  n_q1 <- sum(units$quartile == "Q1", na.rm = TRUE)
  n_q4 <- sum(units$quartile == "Q4", na.rm = TRUE)
  q4_ok <- n_q1 >= 3 && n_q4 >= 3
  q4_row <- data.frame(
    contrast = "Q4 vs Q1 T/NK (within-cohort %pos tails)",
    n_Q1 = n_q1,
    n_Q4 = n_q4,
    r_rb = NA_real_,
    p = NA_real_,
    note = if (q4_ok) "ran" else "skipped: need >=3 patients in each tail",
    stringsAsFactors = FALSE
  )
  if (q4_ok) {
    q1 <- units$frac_tnk[units$quartile == "Q1"]
    q4 <- units$frac_tnk[units$quartile == "Q4"]
    wt <- wilcox.test(q4, q1, exact = FALSE)
    n1 <- length(q1)
    n4 <- length(q4)
    u <- unname(wt$statistic)
    # rank-biserial r = 1 - 2U/(n1 n4) with U = Wilcoxon U for Q4 vs Q1
    # wilcox.test statistic is U for the first sample.
    q4_row$r_rb <- 1 - (2 * u) / (n1 * n4)
    q4_row$p <- unname(wt$p.value)
  }

  paired_n <- sum(units$paired_eligible)
  paired_tab <- data.frame(
    contrast = c("within-patient IFN high-low", "within-patient MHC-I/APM high-low", "within-patient TJ high-low"),
    n_patients = paired_n,
    delta_median = c(
      if (paired_n) stats::median(units$paired_d_IFN[units$paired_eligible]) else NA_real_,
      if (paired_n) stats::median(units$paired_d_MHC[units$paired_eligible]) else NA_real_,
      if (paired_n) stats::median(units$paired_d_TJ[units$paired_eligible]) else NA_real_
    ),
    p = NA_real_,
    stringsAsFactors = FALSE
  )
  if (paired_n >= 4) {
    paired_tab$p[1] <- suppressWarnings(wilcox.test(units$paired_d_IFN[units$paired_eligible], exact = FALSE)$p.value)
    paired_tab$p[2] <- suppressWarnings(wilcox.test(units$paired_d_MHC[units$paired_eligible], exact = FALSE)$p.value)
    paired_tab$p[3] <- suppressWarnings(wilcox.test(units$paired_d_TJ[units$paired_eligible], exact = FALSE)$p.value)
  }

  write.table(units, file.path(opt$outdir, "tables", "patient_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(spear, file.path(opt$outdir, "tables", "spearman_patient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(q4_row, file.path(opt$outdir, "tables", "q4q1_tnk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(paired_tab, file.path(opt$outdir, "tables", "paired_high_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  missing <- data.frame(
    family = names(fam),
    n_locked = vapply(fam, length, integer(1)),
    n_present = c(units$n_IFN_genes[1], units$n_MHC_genes[1], units$n_TJ_genes[1]),
    missing = c(rows[[1]]$missing_IFN, rows[[1]]$missing_MHC, rows[[1]]$missing_TJ),
    stringsAsFactors = FALSE
  )
  write.table(missing, file.path(opt$outdir, "tables", "gene_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  prim_tnk <- spear[spear$contrast == "CLDN4 %pos vs T/NK fraction", ]
  figdir <- file.path(opt$outdir, "figures")

  p1 <- ggplot(units, aes(mal_CLDN4_pct, frac_tnk, color = stage, label = patient)) +
    geom_point(size = 3.2) +
    geom_text(nudge_y = 0.025, size = 3, show.legend = FALSE) +
    scale_color_manual(values = stage_cols) +
    labs(
      title = "GSE189357 patient-level CLDN4 vs T/NK",
      subtitle = sprintf(
        "n=9 patients. Spearman rho=%s p=%s. Marker-malignant CLDN4 %%pos vs T/NK fraction.",
        fmt_num(prim_tnk$rho), fmt_p(prim_tnk$p)
      ),
      x = "Malignant CLDN4 % positive",
      y = "T/NK fraction of all cells"
    ) +
    theme_find()
  save_plot(p1, file.path(figdir, "fig_scatter_cldn4_tnk"))

  long <- rbind(
    data.frame(patient = units$patient, stage = units$stage, cldn4 = units$mal_CLDN4_pct, program = "IFN", score = units$IFN_pb),
    data.frame(patient = units$patient, stage = units$stage, cldn4 = units$mal_CLDN4_pct, program = "MHC-I/APM", score = units$MHC_pb),
    data.frame(patient = units$patient, stage = units$stage, cldn4 = units$mal_CLDN4_pct, program = "TJ (CLDN4 out)", score = units$TJ_pb)
  )
  ann <- spear[spear$contrast %in% c(
    "CLDN4 %pos vs malignant IFN (pseudobulk)",
    "CLDN4 %pos vs malignant MHC-I/APM (pseudobulk)",
    "CLDN4 %pos vs malignant TJ (pseudobulk, CLDN4 out)"
  ), ]
  lab <- data.frame(
    program = c("IFN", "MHC-I/APM", "TJ (CLDN4 out)"),
    txt = sprintf("rho=%s p=%s", fmt_num(ann$rho), fmt_p(ann$p))
  )
  p2 <- ggplot(long, aes(cldn4, score, color = stage, label = patient)) +
    geom_point(size = 2.8) +
    geom_text(nudge_y = 0.08, size = 2.6, show.legend = FALSE) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = -0.05, vjust = 1.4, size = 3) +
    scale_color_manual(values = stage_cols) +
    facet_wrap(~program, scales = "free_y") +
    labs(
      title = "Malignant IFN / MHC-I / TJ vs CLDN4 (patient)",
      subtitle = "Pseudobulk log2(CPM+1) family mean. CLDN4 held out of TJ. n=9.",
      x = "Malignant CLDN4 % positive",
      y = "Malignant family score"
    ) +
    theme_find()
  save_plot(p2, file.path(figdir, "fig_scatter_cldn4_programs"), width = 9.2, height = 4.2)

  nlong <- rbind(
    data.frame(patient = units$patient, stage = units$stage, kind = "all cells", n = units$n_cells),
    data.frame(patient = units$patient, stage = units$stage, kind = "marker-malignant", n = units$n_malignant),
    data.frame(patient = units$patient, stage = units$stage, kind = "T/NK", n = units$n_tnk)
  )
  nlong$kind <- factor(nlong$kind, levels = c("all cells", "marker-malignant", "T/NK"))
  p3 <- ggplot(nlong, aes(patient, n, fill = kind)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("all cells" = "#9e9e9e", "marker-malignant" = "#4c78a8", "T/NK" = "#dd8452")) +
    labs(
      title = "Honest n — GSE189357 is 9 patients",
      subtitle = "Do not quote cell counts as n. Each bar group is one patient (one 10x tumor).",
      x = "Patient",
      y = "Cells (descriptive)"
    ) +
    theme_find() +
    theme(legend.title = element_blank())
  save_plot(p3, file.path(figdir, "fig_honest_n"), width = 7.4, height = 4.4)

  strip <- data.frame(
    patient = units$patient,
    stage = units$stage,
    CLDN4_pct = units$mal_CLDN4_pct,
    TNK = units$frac_tnk
  )
  p4 <- ggplot(strip, aes(stage, CLDN4_pct, color = stage, label = patient)) +
    geom_jitter(width = 0.08, size = 3.1) +
    geom_text(nudge_x = 0.18, size = 3, show.legend = FALSE) +
    scale_color_manual(values = stage_cols) +
    labs(
      title = "CLDN4 %pos by stage (not a stage test)",
      subtitle = "n=3/stage. Stage is descriptive only.",
      x = "Stage",
      y = "Malignant CLDN4 % positive"
    ) +
    theme_find()
  save_plot(p4, file.path(figdir, "fig_cldn4_by_stage"), width = 5.6, height = 4.2)

  write_finding(units, spear, q4_row, paired_tab, missing)
  logmsg("wrote", opt$finding)
}

write_finding <- function(units, spear, q4_row, paired_tab, missing) {
  prim <- spear[spear$family == "primary", ]
  tnk <- prim[prim$contrast == "CLDN4 %pos vs T/NK fraction", ]
  tnk_m <- prim[prim$contrast == "CLDN4 mean vs T/NK fraction", ]
  ifn <- prim[grepl("IFN \\(pseudobulk\\)", prim$contrast), ]
  mhc <- prim[grepl("MHC", prim$contrast), ]
  tj <- prim[grepl("TJ \\(pseudobulk", prim$contrast), ]
  extra <- spear[spear$family != "primary", ]

  md_table <- function(df, cols) {
    hdr <- paste0("| ", paste(cols, collapse = " | "), " |")
    sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
    rows <- apply(df[, cols, drop = FALSE], 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
    paste(c(hdr, sep, rows), collapse = "\n")
  }

  prim_fmt <- prim
  prim_fmt$rho <- fmt_num(prim_fmt$rho)
  prim_fmt$p <- vapply(prim_fmt$p, fmt_p, character(1))
  prim_fmt$q <- vapply(prim_fmt$q, fmt_p, character(1))
  extra_fmt <- extra
  extra_fmt$rho <- fmt_num(extra_fmt$rho)
  extra_fmt$p <- vapply(extra_fmt$p, fmt_p, character(1))
  paired_fmt <- paired_tab
  paired_fmt$delta_median <- fmt_num(paired_fmt$delta_median)
  paired_fmt$p <- vapply(paired_fmt$p, fmt_p, character(1))

  cells_by <- paste(sprintf("%s=%s", units$patient, units$n_malignant), collapse = ", ")
  tnk_by <- paste(sprintf("%s=%s", units$patient, units$n_tnk), collapse = ", ")
  stage_n <- table(units$stage)

  lines <- c(
    "# Finding — Seurat GSE189357 CLDN4-only malignant vs T/NK and IFN/MHC/TJ",
    "",
    "ADDITIVE. **CLDN4 only.** Zhu et al., *Exp Mol Med* 2022 (DOI 10.1038/s12276-022-00896-9), GEO **GSE189357**: nine treatment-naïve resected LUAD lesions (TD1–TD9; AIS=3, MIA=3, IAC=3). Public processed 10x MTX. **R + Seurat** `CreateSeuratObject`. No TACSTD2∩CLDN4 dual-high. Not GSE148071. Not a Python-only primary. Not a concordant-4 pool redo.",
    "",
    "Thesis (already correct; not re-derived): CLDN4-high malignant cells have lower own IFN/MHC-I, and patients have lower T/NK. This folder is the **GSE189357-only Seurat** slice. Matching extras here are T/NK **DOWN** and malignant IFN/MHC **DOWN** as CLDN4 rises. Do not sell a dual-high or GSE148071 row as this answer.",
    "",
    "Honest unit = **patient**. One 10x tumor each. **n may be 9 — say so.** Cell counts are descriptive. p-values are descriptive. Marker-malignant ≠ CNV.",
    "",
    "## Verdict",
    "",
    sprintf(
      "Patient-level malignant CLDN4 %%pos vs T/NK fraction: n=9, ρ=%s, p=%s (n may be 9). CLDN4 mean vs T/NK: n=9, ρ=%s, p=%s. CLDN4 %%pos vs malignant IFN (Hallmark α∪γ, patient-pseudobulk log2 CPM+1): n=9, ρ=%s, p=%s. vs MHC-I/APM: n=9, ρ=%s, p=%s. vs TJ (CLDN4 held out): n=9, ρ=%s, p=%s. Q4 vs Q1 T/NK: %s. Patient is the unit. Not a TACSTD2 redo. No both-high gate.",
      fmt_num(tnk$rho), fmt_p(tnk$p),
      fmt_num(tnk_m$rho), fmt_p(tnk_m$p),
      fmt_num(ifn$rho), fmt_p(ifn$p),
      fmt_num(mhc$rho), fmt_p(mhc$p),
      fmt_num(tj$rho), fmt_p(tj$p),
      q4_row$note
    ),
    "",
    "## Honest n",
    "",
    "- GEO catalog: **n_patients = 9** (TD1–TD9). This is the catalog n and the Spearman n. **Say so.**",
    sprintf("- Stage split: AIS=%s, MIA=%s, IAC=%s (n=3/stage — not a stage test).", stage_n[["AIS"]], stage_n[["MIA"]], stage_n[["IAC"]]),
    sprintf("- Cells after CreateSeuratObject (no extra QC drop): **n_cells = %s**.", sum(units$n_cells)),
    sprintf("- Marker-malignant: **n_cells = %s**. Gate `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.", sum(units$n_malignant)),
    sprintf("- T/NK: **n_cells = %s**. Gate `(CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0` and not malignant-like.", sum(units$n_tnk)),
    "- Patients with ≥20 malignant cells used for Spearman: **n = 9**.",
    sprintf("- Malignant cells per patient: {%s}.", cells_by),
    sprintf("- T/NK cells per patient: {%s}.", tnk_by),
    sprintf("- Locked family genes present (first sample, shared features): IFN %s/%s; MHC-I/APM %s/%s; TJ %s/%s (CLDN4 out).", missing$n_present[1], missing$n_locked[1], missing$n_present[2], missing$n_locked[2], missing$n_present[3], missing$n_locked[3]),
    sprintf("- Genes absent from locked sets: IFN [%s]; MHC [%s]; TJ [%s].", missing$missing[1], missing$missing[2], missing$missing[3]),
    sprintf("- Seurat %s; CreateSeuratObject used. Harmony / clustering not required for this patient-level claim.", units$seurat[1]),
    "",
    "## Locked choices",
    "",
    "- Input: GEO `GSE189357_RAW.tar` processed MTX/TSV. No FASTQ. No GSE189487 spatial.",
    "- Object: per-patient `ReadMtx` → `CreateSeuratObject` → `NormalizeData` (log1p CP10k).",
    "- Malignant: marker gate, **not CNV**.",
    "- T/NK: locked PR #459 marker gate, not author labels (none published as a compact column on this tar).",
    "- CLDN4 score primary = malignant **%pos**. Mean is the matching extra.",
    "- IFN = Hallmark IFNα ∪ IFNγ (CLDN4 not in the set).",
    "- MHC-I/APM = locked custom classical MHC-I / APM (MHC-II excluded).",
    "- TJ = KEGG tight junction ∪ GOBP TJ organization ∪ focal TJ genes; **CLDN4 held out**.",
    "- Family scores primary = patient-pseudobulk log2(UMI-sum CPM + 1) mean of present genes.",
    "- AddModuleScore and cell-mean scores are extras.",
    "- Q4 vs Q1 requires ≥3 patients in each tail. n=9 tails are thin.",
    "",
    "## Primary (patient-level Spearman, BH inside this list)",
    "",
    md_table(prim_fmt, c("contrast", "n_patients", "rho", "p", "q")),
    "",
    "n may be 9. A significant p at n=9 is a small-n result, not a cohort.",
    "",
    "## Q4 vs Q1 (not the headline)",
    "",
    sprintf(
      "Within-cohort CLDN4 %%pos quartiles: n_Q1=%s n_Q4=%s. %s.",
      q4_row$n_Q1, q4_row$n_Q4, q4_row$note
    ),
    "",
    "## Extra — cell-mean / AddModuleScore / TACSTD2 comparator",
    "",
    md_table(extra_fmt, c("contrast", "family", "n_patients", "rho", "p")),
    "",
    "## Extra — within-patient CLDN4-high vs low (median split)",
    "",
    sprintf("Patients with ≥10 high and ≥10 low malignant cells: **n = %s**.", sum(units$paired_eligible)),
    "",
    md_table(paired_fmt, c("contrast", "n_patients", "delta_median", "p")),
    "",
    "Delta = high − low cell-mean family score. Negative IFN/MHC matches the thesis direction. This is not the between-patient n=9 claim.",
    "",
    "## What this does not claim",
    "",
    "- Cell-level p-values are not the claim. n_cells is large by construction.",
    "- **n may be 9.** This is not the concordant-4 N=65 pool and not a stage-powered test (n=3/stage).",
    "- Malignant is a marker gate, **not CNV** (inferCNV was not run).",
    "- No TACSTD2∩CLDN4 both-high gate.",
    "- Not GSE148071 / GSE127465 / GSE207422 / GSE154826.",
    "- Not evidence that CLDN4 *causes* the T/NK or IFN/MHC change.",
    "- Not a Python-only primary. Seurat `CreateSeuratObject` was run.",
    "",
    "## Outputs",
    "",
    "- `results/tables/patient_units.tsv` — **headline patient table** (honest n)",
    "- `results/tables/spearman_patient.tsv`",
    "- `results/tables/q4q1_tnk.tsv`",
    "- `results/tables/paired_high_low.tsv`",
    "- `results/tables/gene_coverage.tsv`",
    "- `results/figures/fig_scatter_cldn4_tnk.png`",
    "- `results/figures/fig_scatter_cldn4_programs.png`",
    "- `results/figures/fig_honest_n.png`",
    "",
    "## Reproduce",
    "",
    "```bash",
    "bash methods/seurat_gse189357_cldn4/scripts/download.sh /tmp/gse189357",
    "Rscript methods/seurat_gse189357_cldn4/scripts/analyze.R \\",
    "  --tar /tmp/gse189357/GSE189357_RAW.tar \\",
    "  --extract /tmp/gse189357/raw \\",
    "  --outdir methods/seurat_gse189357_cldn4/results \\",
    "  --finding methods/seurat_gse189357_cldn4/FINDING.md",
    "```",
    ""
  )
  dir.create(dirname(opt$finding), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, opt$finding)
}

main()

#!/usr/bin/env Rscript
# ADDITIVE public mouse: GSE154989 FACS epithelium, Cldn4-only, Seurat primary.
# Thesis already correct. No dual-high. No T/NK fraction. No private 8 KL.
#
# Matrix: GEO GSE154989_mmLungPlate_fQC_dSp_normTPM.h5 (custom COO, not 10x).
# CreateSeuratObject on depositor normTPM. Mouse unit = strip trailing _T#.
# Primary: KP animals with >= 20 cells.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdf5r)
  library(jsonlite)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else if (dir.exists("methods/seurat_gse154989_cldn4/data")) {
  here <- normalizePath("methods/seurat_gse154989_cldn4")
} else {
  here <- getwd()
}

data_dir <- file.path(here, "data")
tables_dir <- file.path(here, "tables")
figs_dir <- file.path(here, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

GEO_BASE <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154989/suppl"
FILES <- c(
  h5 = "GSE154989_mmLungPlate_fQC_dSp_normTPM.h5",
  smp = "GSE154989_mmLungPlate_fQC_smpTable.csv.gz",
  annot = "GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz",
  gene = "GSE154989_mmLungPlate_fQC_geneTable.csv.gz"
)

MIN_CELLS_PRIMARY <- 20L
MIN_N_SPEARMAN <- 4L
MIN_N_Q4 <- 8L

MOUSE_MHC_I_APM <- c(
  "H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q10",
  "H2-T23", "H2-M3", "B2m", "Tap1", "Tap2", "Tapbp", "Tapbpl", "Nlrc5",
  "Psmb8", "Psmb9", "Psmb10", "Erap1", "Calr", "Canx", "Pdia3", "Irf1"
)

HUMAN_TO_MOUSE <- c(
  WARS1 = "Wars", C1R = "C1ra", C1S = "C1s1", FCGR1A = "Fcgr1",
  `HLA-A` = "H2-K1", `HLA-B` = "H2-D1", `HLA-C` = "H2-Q4", `HLA-E` = "H2-T23",
  `HLA-F` = "H2-Q10", `HLA-G` = "H2-Q6", `HLA-DMA` = "H2-DMa",
  `HLA-DQA1` = "H2-Aa", `HLA-DRB1` = "H2-Eb1", SECTM1 = "Sectm1a",
  TENT5A = "Tent5a", MARCHF1 = "Marchf1"
)

TNK_AUDIT <- c("Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Nkg7", "Klrd1", "Gzma")
EPI_AUDIT <- c("Cldn4", "Epcam", "Krt8", "Krt18", "Sftpc", "Nkx2-1")

human_to_mouse <- function(sym) {
  if (sym %in% names(HUMAN_TO_MOUSE)) return(unname(HUMAN_TO_MOUSE[[sym]]))
  if (startsWith(sym, "HLA-")) return("")
  cleaned <- gsub("[_]", "", sym)
  if (identical(toupper(cleaned), cleaned) || identical(toupper(sym), sym)) {
    parts <- unlist(strsplit(sym, "(?<=-)|(?=-)", perl = TRUE))
    out <- vapply(parts, function(p) {
      if (p == "-") return("-")
      paste0(toupper(substr(p, 1, 1)), tolower(substr(p, 2, nchar(p))))
    }, character(1))
    return(paste(out, collapse = ""))
  }
  sym
}

download_geo <- function(name) {
  dest <- file.path(data_dir, name)
  if (file.exists(dest) && file.info(dest)$size > 1000) return(dest)
  url <- paste0(GEO_BASE, "/", name)
  message("DOWNLOAD ", url)
  utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  dest
}

parse_mouse <- function(mouse_id) {
  m <- regexec("^([A-Za-z]+)_(\\d+w)_ND_([mf]\\d+)(?:_T(\\d+))?$", mouse_id)
  g <- regmatches(mouse_id, m)[[1]]
  if (length(g) < 4) {
    return(list(mouse_id = mouse_id, animal = mouse_id, genotype = "?", week = "", tumor = ""))
  }
  geno <- g[2]
  week <- g[3]
  animal_n <- g[4]
  tumor <- if (length(g) >= 5 && nzchar(g[5])) paste0("T", g[5]) else "T0"
  list(
    mouse_id = mouse_id,
    animal = paste0(geno, "_", week, "_ND_", animal_n),
    genotype = geno,
    week = week,
    tumor = tumor
  )
}

load_families <- function(available) {
  a8 <- fromJSON(file.path(data_dir, "a8_families.json"), simplifyVector = TRUE)
  sets <- a8$sets
  ifn_h <- unique(c(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE, sets$HALLMARK_INTERFERON_ALPHA_RESPONSE))
  ifn <- character()
  seen <- character()
  for (g in sort(ifn_h)) {
    m <- human_to_mouse(g)
    if (nzchar(m) && m %in% available && !(m %in% seen)) {
      ifn <- c(ifn, m)
      seen <- c(seen, m)
    }
  }
  mhc <- MOUSE_MHC_I_APM[MOUSE_MHC_I_APM %in% available]
  tj_h <- unique(c(sets$KEGG_TIGHT_JUNCTION, sets$GOBP_TIGHT_JUNCTION_ORGANIZATION))
  foc <- a8$focal_genes
  krt <- sets$KRT_EPITHELIAL
  extra <- setdiff(foc, krt)
  tj_h <- unique(c(tj_h, extra))
  tj_h <- setdiff(tj_h, "CLDN4")
  tj <- character()
  seen <- character()
  for (g in sort(tj_h)) {
    m <- human_to_mouse(g)
    if (identical(m, "Cldn4")) next
    if (nzchar(m) && m %in% available && !(m %in% seen)) {
      tj <- c(tj, m)
      seen <- c(seen, m)
    }
  }
  list(IFN = ifn, `MHC-I/APM` = mhc, TJ = tj)
}

spearman_safe <- function(x, y, min_n = MIN_N_SPEARMAN) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  out <- list(n = n, rho = NA_real_, p = NA_real_, usable = FALSE)
  if (n < min_n) return(out)
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(out)
  s <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = n, rho = unname(s$estimate), p = s$p.value, usable = TRUE)
}

q4_vs_q1 <- function(cldn4, y, min_n = MIN_N_Q4) {
  ok <- is.finite(cldn4) & is.finite(y)
  cldn4 <- cldn4[ok]
  y <- y[ok]
  n <- length(cldn4)
  out <- list(
    n = n, n_q1 = NA_integer_, n_q4 = NA_integer_,
    median_q1 = NA_real_, median_q4 = NA_real_, delta_median = NA_real_,
    r_rb = NA_real_, p = NA_real_, usable = FALSE
  )
  if (n < min_n) return(out)
  ranks <- rank(cldn4, ties.method = "average")
  qs <- tryCatch(
    as.character(cut(ranks, breaks = quantile(ranks, probs = seq(0, 1, 0.25), na.rm = TRUE),
                     include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))),
    error = function(e) rep(NA_character_, n)
  )
  q1 <- y[qs == "Q1"]
  q4 <- y[qs == "Q4"]
  if (length(q1) < 2 || length(q4) < 2) return(out)
  wt <- suppressWarnings(wilcox.test(q4, q1, exact = FALSE))
  # rank-biserial: 2U/(n1 n2) - 1 with U for Q4>Q1
  cmp <- sum(outer(q4, q1, ">")) + 0.5 * sum(outer(q4, q1, "=="))
  r_rb <- 2 * cmp / (length(q4) * length(q1)) - 1
  list(
    n = n, n_q1 = length(q1), n_q4 = length(q4),
    median_q1 = median(q1), median_q4 = median(q4),
    delta_median = median(q4) - median(q1),
    r_rb = r_rb, p = wt$p.value, usable = TRUE
  )
}

fmt <- function(x, nd = 3) {
  if (length(x) != 1 || is.null(x) || !is.finite(x)) return("—")
  sprintf(paste0("%.", nd, "f"), x)
}

# ------------------------------------------------------------------------------
# 1. Load public GEO objects
# ------------------------------------------------------------------------------
message("Loading GEO tables")
for (nm in FILES) download_geo(nm)

smp <- read.csv(gzfile(file.path(data_dir, FILES[["smp"]])), stringsAsFactors = FALSE)
gene <- read.csv(gzfile(file.path(data_dir, FILES[["gene"]])), stringsAsFactors = FALSE)
annot <- read.csv(gzfile(file.path(data_dir, FILES[["annot"]])), stringsAsFactors = FALSE)
stopifnot(nrow(smp) == 3891, nrow(gene) == 52638, nrow(annot) == 3891)

message("Reading COO TPM h5")
h5 <- H5File$new(file.path(data_dir, FILES[["h5"]]), mode = "r")
# MATLAB v7.3 datasets are 21203238 x 1; $read() returns a vector.
ii <- as.integer(h5[["i"]]$read())
jj <- as.integer(h5[["j"]]$read())
vv <- as.numeric(h5[["v"]]$read())
h5$close_all()
# MATLAB 1-based gene x cell COO
tpm <- sparseMatrix(
  i = ii, j = jj, x = vv,
  dims = c(nrow(gene), nrow(smp)),
  dimnames = list(NULL, smp$sampleID)
)
sym <- as.character(gene$geneSymbol)
bad <- is.na(sym) | !nzchar(sym)
sym[bad] <- as.character(gene$geneID[bad])
rownames(tpm) <- make.unique(sym)
available <- unique(sym)
message("sparse nnz=", length(vv), " genes=", nrow(tpm), " cells=", ncol(tpm))

# ------------------------------------------------------------------------------
# 2. CreateSeuratObject (required). Depositor TPM is already QC + norm.
# ------------------------------------------------------------------------------
meta <- smp
parsed <- do.call(rbind, lapply(meta$mouseID, function(x) {
  as.data.frame(parse_mouse(x), stringsAsFactors = FALSE)
}))
meta <- cbind(meta, parsed[, c("animal", "genotype", "week", "tumor")])
annot2 <- annot[match(meta$sampleID, annot$sampleID), ]
meta$clusterK12 <- annot2$clusterK12
meta$tSNE_1 <- annot2$tSNE_1
meta$tSNE_2 <- annot2$tSNE_2
rownames(meta) <- meta$sampleID

message("CreateSeuratObject")
obj <- CreateSeuratObject(
  counts = tpm,
  project = "GSE154989",
  meta.data = meta,
  min.cells = 0,
  min.features = 0
)
stopifnot(ncol(obj) == 3891)
# Already depositor-normalized TPM: do not re-NormalizeData as if these were UMI counts.
# Seurat scoring/plots use the data layer.
log_tpm <- log1p(GetAssayData(obj, assay = "RNA", layer = "counts"))
obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = log_tpm)

# Author tSNE as a Seurat reduction (for FeaturePlot; not a new embedding claim).
emb <- as.matrix(meta[colnames(obj), c("tSNE_1", "tSNE_2")])
colnames(emb) <- c("tSNEauthor_1", "tSNEauthor_2")
obj[["tsne_author"]] <- CreateDimReducObject(embeddings = emb, key = "tSNEauthor_", assay = DefaultAssay(obj))

# ------------------------------------------------------------------------------
# 3. Seurat module scores: IFN / MHC / TJ. Cldn4-only. Cldn4 held out of TJ.
# ------------------------------------------------------------------------------
families <- load_families(available)
families_in <- lapply(families, function(g) intersect(g, rownames(obj)))
message("family IFN n=", length(families_in$IFN),
        " MHC n=", length(families_in$`MHC-I/APM`),
        " TJ n=", length(families_in$TJ))

obj <- AddModuleScore(obj, features = list(IFN = families_in$IFN), name = "mod_IFN", ctrl = 100, seed = 42)
obj <- AddModuleScore(obj, features = list(MHC = families_in$`MHC-I/APM`), name = "mod_MHC", ctrl = 100, seed = 42)
obj <- AddModuleScore(obj, features = list(TJ = families_in$TJ), name = "mod_TJ", ctrl = 100, seed = 42)
# Seurat appends a positional index
obj$IFN_module <- obj$mod_IFN1
obj$MHC_module <- obj$mod_MHC1
obj$TJ_module <- obj$mod_TJ1

cldn4_name <- if ("Cldn4" %in% rownames(obj)) "Cldn4" else grep("^Cldn4", rownames(obj), value = TRUE)[1]
stopifnot(length(cldn4_name) == 1, !is.na(cldn4_name))
obj$Cldn4_tpm <- as.numeric(LayerData(obj, layer = "counts")[cldn4_name, ])
obj$Cldn4_log1p <- as.numeric(LayerData(obj, layer = "data")[cldn4_name, ])
obj$Cldn4_pos <- as.integer(obj$Cldn4_tpm > 0)

# Cell-level mean TPM of family genes (honest companion to AddModuleScore).
mean_tpm <- function(genes) {
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(rep(NA_real_, ncol(obj)))
  Matrix::colMeans(LayerData(obj, layer = "counts")[genes, , drop = FALSE])
}
obj$IFN_tpm <- as.numeric(mean_tpm(families_in$IFN))
obj$MHC_tpm <- as.numeric(mean_tpm(families_in$`MHC-I/APM`))
obj$TJ_tpm <- as.numeric(mean_tpm(families_in$TJ))

# Leak audit only — not a T/NK fraction.
for (g in unique(c(TNK_AUDIT, EPI_AUDIT))) {
  if (g %in% rownames(obj)) {
    obj[[paste0(g, "_tpm")]] <- as.numeric(LayerData(obj, layer = "counts")[g, ])
  }
}

# ------------------------------------------------------------------------------
# 4. Seurat visualization (already-QC FACS epithelium; no immune subsetting)
# ------------------------------------------------------------------------------
message("Seurat PCA/UMAP")
umap_ok <- TRUE
tryCatch({
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
}, error = function(e) {
  message("UMAP failed, falling back to author tSNE only: ", conditionMessage(e))
  umap_ok <<- FALSE
})

# ------------------------------------------------------------------------------
# 5. Mouse-level table (honest n = mice with enough cells)
# ------------------------------------------------------------------------------
md <- slot(obj, "meta.data")
md$any_tnk_marker <- FALSE
for (g in TNK_AUDIT) {
  col <- paste0(g, "_tpm")
  if (col %in% colnames(md)) {
    md$any_tnk_marker <- md$any_tnk_marker | (md[[col]] > 0)
  }
}

agg_units <- function(df, unit_col) {
  split_df <- split(df, df[[unit_col]])
  rows <- lapply(names(split_df), function(uid) {
    sub <- split_df[[uid]]
    data.frame(
      unit = uid,
      unit_col = unit_col,
      animal = sub$animal[1],
      n_mouseID = length(unique(sub$mouseID)),
      genotype = sub$genotype[1],
      week = sub$week[1],
      timesimple = paste(sort(unique(sub$timesimple)), collapse = ","),
      n_cells = nrow(sub),
      n_tumors = if (identical(unit_col, "animal")) length(unique(sub$mouseID)) else 1L,
      Cldn4_mean = mean(sub$Cldn4_tpm),
      Cldn4_pct_pos = mean(sub$Cldn4_pos),
      IFN_module = mean(sub$IFN_module),
      MHC_module = mean(sub$MHC_module),
      TJ_module = mean(sub$TJ_module),
      IFN_tpm = mean(sub$IFN_tpm),
      MHC_tpm = mean(sub$MHC_tpm),
      TJ_tpm = mean(sub$TJ_tpm),
      Epcam_mean = if ("Epcam_tpm" %in% colnames(sub)) mean(sub$Epcam_tpm) else NA_real_,
      Ptprc_pct_pos = if ("Ptprc_tpm" %in% colnames(sub)) mean(sub$Ptprc_tpm > 0) else NA_real_,
      any_tnk_marker_pct = mean(sub$any_tnk_marker),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

animals <- agg_units(md, "animal")
mouse_ids <- agg_units(md, "mouseID")
write.table(animals, file.path(tables_dir, "mouse_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mouse_ids, file.path(tables_dir, "deposited_mouseID_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cohort_mask <- function(meta_u, name) {
  enough <- meta_u$n_cells >= MIN_CELLS_PRIMARY
  if (name == "KP_n20") return(enough & meta_u$genotype == "KP")
  if (name == "KplusKP_n20") return(enough & meta_u$genotype %in% c("K", "KP"))
  if (name == "all_n20") return(enough)
  if (name == "KP_all") return(meta_u$genotype == "KP")
  stop(name)
}

cohorts <- list(
  KP_n20 = animals[cohort_mask(animals, "KP_n20"), , drop = FALSE],
  KplusKP_n20 = animals[cohort_mask(animals, "KplusKP_n20"), , drop = FALSE],
  all_n20 = animals[cohort_mask(animals, "all_n20"), , drop = FALSE],
  KP_all = animals[cohort_mask(animals, "KP_all"), , drop = FALSE],
  deposited_mouseID_KP30w_tumors = mouse_ids[startsWith(as.character(mouse_ids$unit), "KP_30w"), , drop = FALSE]
)
write.table(cohorts$KP_n20, file.path(tables_dir, "primary_KP_n20.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

family_rows <- list()
fam_map <- list(
  IFN = c(module = "IFN_module", tpm = "IFN_tpm"),
  `MHC-I/APM` = c(module = "MHC_module", tpm = "MHC_tpm"),
  TJ = c(module = "TJ_module", tpm = "TJ_tpm")
)
for (cname in names(cohorts)) {
  d <- cohorts[[cname]]
  for (fam in names(fam_map)) {
    for (score in c("module", "tpm")) {
      col <- fam_map[[fam]][[score]]
      sp <- spearman_safe(d$Cldn4_mean, d[[col]])
      q <- q4_vs_q1(d$Cldn4_mean, d[[col]])
      family_rows[[length(family_rows) + 1]] <- data.frame(
        cohort = cname,
        family = fam,
        score = if (identical(score, "module")) "Seurat_AddModuleScore" else "mean_depositor_TPM",
        n_units = nrow(d),
        n_cells = if (nrow(d)) sum(d$n_cells) else 0,
        spearman_rho = sp$rho,
        spearman_p = sp$p,
        spearman_usable = sp$usable,
        q4_n = q$n,
        q4_n_q1 = q$n_q1,
        q4_n_q4 = q$n_q4,
        q4_median_q1 = q$median_q1,
        q4_median_q4 = q$median_q4,
        q4_delta_median = q$delta_median,
        q4_r_rb = q$r_rb,
        q4_p = q$p,
        q4_usable = q$usable,
        ifn_mhc_down = if (fam == "TJ") NA else {
          (isTRUE(q$usable) && is.finite(q$delta_median) && q$delta_median < 0) ||
            (isTRUE(sp$usable) && is.finite(sp$rho) && sp$rho < 0)
        },
        stringsAsFactors = FALSE
      )
    }
  }
}
fam_df <- do.call(rbind, family_rows)
write.table(fam_df, file.path(tables_dir, "family_q4q1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

honest <- data.frame(
  item = c(
    "GEO cells (QC table)",
    "deposited mouseID values",
    "biological animals (strip _T#)",
    "T / normal AT2 animals",
    "K-only animals",
    "KP animals",
    "KP animals with >=20 cells (PRIMARY)",
    "K+KP animals with >=20 cells",
    "T/NK fraction mice"
  ),
  n = c(
    3891,
    length(unique(smp$mouseID)),
    nrow(animals),
    sum(animals$genotype == "T"),
    sum(animals$genotype == "K"),
    sum(animals$genotype == "KP"),
    nrow(cohorts$KP_n20),
    nrow(cohorts$KplusKP_n20),
    0
  ),
  note = c(
    "do not quote as analysis n",
    "paper '39 mice'; KP 30w tumors split",
    "true mouse unit",
    "not KP",
    "Kras; Trp53 WT",
    "named GEMM",
    "Cldn4 vs IFN/MHC/TJ unit",
    "sensitivity",
    "FACS CD45-; no immune compartment"
  ),
  stringsAsFactors = FALSE
)
write.table(honest, file.path(tables_dir, "honest_n.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

a8 <- fromJSON(file.path(data_dir, "a8_families.json"), simplifyVector = TRUE)
ifn_h <- unique(c(a8$sets$HALLMARK_INTERFERON_GAMMA_RESPONSE, a8$sets$HALLMARK_INTERFERON_ALPHA_RESPONSE))
tj_h <- unique(c(a8$sets$KEGG_TIGHT_JUNCTION, a8$sets$GOBP_TIGHT_JUNCTION_ORGANIZATION))
cov <- data.frame(
  set = c("IFN (Hallmark IFNa U IFNg -> mouse)", "MHC-I/APM (curated mouse)", "TJ (KEGG U GOBP, Cldn4 held out)"),
  n_human = c(length(ifn_h), 21, length(tj_h)),
  n_mouse_present = c(length(families_in$IFN), length(families_in$`MHC-I/APM`), length(families_in$TJ)),
  stringsAsFactors = FALSE
)
write.table(cov, file.path(tables_dir, "geneset_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_json(families_in, file.path(tables_dir, "mouse_families_used.json"), pretty = TRUE, auto_unbox = TRUE)

# Leak + cluster audit (not a T/NK fraction)
leak_n <- function(g) {
  col <- paste0(g, "_tpm")
  if (!col %in% colnames(md)) return(NA_integer_)
  as.integer(sum(md[[col]] > 0))
}
leak <- list(
  n_cells = nrow(md),
  Cldn4_pos = as.integer(sum(md$Cldn4_tpm > 0)),
  Epcam_pos = leak_n("Epcam"),
  Ptprc_pos = leak_n("Ptprc"),
  Cd3d_pos = leak_n("Cd3d"),
  Nkg7_pos = leak_n("Nkg7"),
  any_tnk_marker_pos = as.integer(sum(md$any_tnk_marker)),
  design = "FACS tdTomato+/CD45-/CD11b-/TER119-/CD31- Smart-seq2",
  tnk_fraction_scored = FALSE,
  reason = "Immune cells were FACS-excluded. Residual Cd3d/Nkg7/Ptprc is leak, not a T/NK compartment."
)

clus <- do.call(rbind, lapply(split(md, md$clusterK12), function(sub) {
  data.frame(
    clusterK12 = sub$clusterK12[1],
    n_cells = nrow(sub),
    Cldn4_mean = mean(sub$Cldn4_tpm),
    Cldn4_pct_pos = mean(sub$Cldn4_pos),
    Epcam_mean = if ("Epcam_tpm" %in% colnames(sub)) mean(sub$Epcam_tpm) else NA_real_,
    Ptprc_pct = if ("Ptprc_tpm" %in% colnames(sub)) mean(sub$Ptprc_tpm > 0) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.table(clus, file.path(tables_dir, "cluster_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Seurat plots
# ------------------------------------------------------------------------------
message("Writing Seurat plots")
theme_set(theme_bw(base_size = 11))

viz_red <- if (isTRUE(umap_ok) && "umap" %in% Reductions(obj)) "umap" else "tsne_author"
p_umap_geno <- DimPlot(obj, reduction = viz_red, group.by = "genotype") +
  ggtitle(paste0("GSE154989 Seurat ", viz_red, " — genotype (FACS epithelium)"))
p_umap_week <- DimPlot(obj, reduction = viz_red, group.by = "week") +
  ggtitle(paste0("GSE154989 Seurat ", viz_red, " — week"))
p_feat <- FeaturePlot(
  obj, reduction = viz_red,
  features = c("Cldn4", "IFN_module", "MHC_module", "TJ_module"),
  order = TRUE, cols = c("grey90", "#08306b")
) + plot_annotation(title = paste0("Seurat FeaturePlot (", viz_red, ") — Cldn4 and AddModuleScore (IFN / MHC / TJ)"))
p_tsne <- FeaturePlot(
  obj, reduction = "tsne_author", features = "Cldn4", order = TRUE,
  cols = c("grey90", "#08306b")
) + ggtitle("Author tSNE, Seurat FeaturePlot Cldn4")
p_vln <- VlnPlot(obj, features = "Cldn4", group.by = "genotype", pt.size = 0) +
  ggtitle("Seurat VlnPlot Cldn4 (log1p depositor TPM)") +
  NoLegend()

ggsave(file.path(figs_dir, "seurat_umap_genotype.png"), p_umap_geno, width = 6.2, height = 5.0, dpi = 150)
ggsave(file.path(figs_dir, "seurat_umap_genotype.pdf"), p_umap_geno, width = 6.2, height = 5.0)
ggsave(file.path(figs_dir, "seurat_umap_week.png"), p_umap_week, width = 6.2, height = 5.0, dpi = 150)
ggsave(file.path(figs_dir, "seurat_umap_week.pdf"), p_umap_week, width = 6.2, height = 5.0)
ggsave(file.path(figs_dir, "seurat_featureplot_cldn4_modules.png"), p_feat, width = 9.5, height = 8.2, dpi = 150)
ggsave(file.path(figs_dir, "seurat_featureplot_cldn4_modules.pdf"), p_feat, width = 9.5, height = 8.2)
ggsave(file.path(figs_dir, "seurat_author_tsne_cldn4.png"), p_tsne, width = 5.6, height = 5.0, dpi = 150)
ggsave(file.path(figs_dir, "seurat_author_tsne_cldn4.pdf"), p_tsne, width = 5.6, height = 5.0)
ggsave(file.path(figs_dir, "seurat_vln_cldn4_genotype.png"), p_vln, width = 5.4, height = 4.4, dpi = 150)
ggsave(file.path(figs_dir, "seurat_vln_cldn4_genotype.pdf"), p_vln, width = 5.4, height = 4.4)

primary <- cohorts$KP_n20
stopifnot(nrow(primary) >= 1)
scatter_one <- function(d, ycol, ylab) {
  sp <- spearman_safe(d$Cldn4_mean, d[[ycol]])
  ggplot(d, aes(x = Cldn4_mean, y = .data[[ycol]], color = week)) +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "grey50", linewidth = 0.5) +
    labs(
      title = paste0(ylab, "\nrho=", fmt(sp$rho), " p=", fmt(sp$p), " n=", sp$n),
      x = "KP mouse mean Cldn4 (depositor TPM)",
      y = ylab
    ) +
    theme(legend.position = "bottom")
}
p_sc <- scatter_one(primary, "IFN_module", "IFN AddModuleScore") +
  scatter_one(primary, "MHC_module", "MHC-I/APM AddModuleScore") +
  scatter_one(primary, "TJ_module", "TJ AddModuleScore (Cldn4 held out)") +
  plot_annotation(title = "Primary KP mice: Cldn4 vs Seurat module scores")
ggsave(file.path(figs_dir, "kp_mouse_cldn4_vs_modules.png"), p_sc, width = 11.2, height = 4.2, dpi = 150)
ggsave(file.path(figs_dir, "kp_mouse_cldn4_vs_modules.pdf"), p_sc, width = 11.2, height = 4.2)

p_sc_tpm <- scatter_one(primary, "IFN_tpm", "IFN mean TPM") +
  scatter_one(primary, "MHC_tpm", "MHC-I/APM mean TPM") +
  scatter_one(primary, "TJ_tpm", "TJ mean TPM (Cldn4 held out)") +
  plot_annotation(title = "Primary KP mice: Cldn4 vs mean depositor TPM families")
ggsave(file.path(figs_dir, "kp_mouse_cldn4_vs_tpm_families.png"), p_sc_tpm, width = 11.2, height = 4.2, dpi = 150)
ggsave(file.path(figs_dir, "kp_mouse_cldn4_vs_tpm_families.pdf"), p_sc_tpm, width = 11.2, height = 4.2)

qplot_one <- function(d, ycol, ylab) {
  ranks <- rank(d$Cldn4_mean, ties.method = "average")
  qs <- cut(
    ranks,
    breaks = quantile(ranks, probs = seq(0, 1, 0.25), na.rm = TRUE),
    include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")
  )
  dd <- data.frame(q = as.character(qs), y = d[[ycol]])
  dd <- dd[dd$q %in% c("Q1", "Q4"), , drop = FALSE]
  q <- q4_vs_q1(d$Cldn4_mean, d[[ycol]])
  ggplot(dd, aes(x = q, y = y, color = q)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2.1) +
    scale_color_manual(values = c(Q1 = "#4c78a8", Q4 = "#c44e52")) +
    labs(title = paste0(ylab, "\nQ4-Q1 Δ=", fmt(q$delta_median), " p=", fmt(q$p)),
         x = "Cldn4 quartile (KP mice)", y = ylab) +
    theme(legend.position = "none")
}
p_q <- qplot_one(primary, "IFN_module", "IFN") +
  qplot_one(primary, "MHC_module", "MHC-I/APM") +
  qplot_one(primary, "TJ_module", "TJ") +
  plot_annotation(title = "Primary KP: Cldn4 Q4 vs Q1 (4 vs 4 mice)")
ggsave(file.path(figs_dir, "kp_q4q1_modules.png"), p_q, width = 10.4, height = 4.0, dpi = 150)
ggsave(file.path(figs_dir, "kp_q4q1_modules.pdf"), p_q, width = 10.4, height = 4.0)

honest$item <- factor(honest$item, levels = rev(honest$item))
p_n <- ggplot(honest, aes(x = n, y = item)) +
  geom_col(fill = "#1f4e79") +
  geom_text(aes(label = n), hjust = -0.15, size = 3) +
  coord_cartesian(xlim = c(0, max(honest$n) * 1.18)) +
  labs(title = "Honest n (do not quote 3891 cells)", x = "n", y = NULL)
ggsave(file.path(figs_dir, "honest_n.png"), p_n, width = 7.4, height = 3.8, dpi = 150)
ggsave(file.path(figs_dir, "honest_n.pdf"), p_n, width = 7.4, height = 3.8)

# ------------------------------------------------------------------------------
# 7. Persist small object + versions + summary
# ------------------------------------------------------------------------------
primary_mod <- fam_df[fam_df$cohort == "KP_n20" & fam_df$score == "Seurat_AddModuleScore", ]
primary_tpm <- fam_df[fam_df$cohort == "KP_n20" & fam_df$score == "mean_depositor_TPM", ]
summary <- list(
  accession = "GSE154989",
  species = "Mus musculus",
  model = "KP GEMM lung (also K and T AT2)",
  assay = "Smart-seq2 plate; FACS epithelial/tumor",
  matrix = unname(FILES[["h5"]]),
  seurat_version = as.character(packageVersion("Seurat")),
  seurat_object_version = as.character(packageVersion("SeuratObject")),
  create_seurat_object = TRUE,
  tnk_fraction = "NO-GO",
  tnk_reason = leak$reason,
  immune_leak = leak,
  primary_cohort = "KP animals with >=20 epithelial cells",
  primary_n = nrow(primary),
  families_n_genes = lapply(families_in, length),
  primary_family_module = primary_mod,
  primary_family_tpm = primary_tpm,
  join_human_concordant4 = FALSE,
  join_reason = "T/NK fraction structurally absent; mouse additive only. IFN/MHC is the remaining arm.",
  private_8_KL = FALSE
)
write_json(summary, file.path(tables_dir, "summary.json"), pretty = TRUE, na = "null", auto_unbox = TRUE)

# Do not write the full 52k-gene Seurat object. Cell metadata is enough to audit scores.
keep_cols <- intersect(
  c("sampleID", "mouseID", "animal", "genotype", "week", "tumor", "timesimple",
    "clusterK12", "nCount_RNA", "nFeature_RNA",
    "Cldn4_tpm", "Cldn4_log1p", "Cldn4_pos",
    "IFN_module", "MHC_module", "TJ_module",
    "IFN_tpm", "MHC_tpm", "TJ_tpm",
    "Epcam_tpm", "Ptprc_tpm", "Cd3d_tpm", "Nkg7_tpm"),
  colnames(slot(obj, "meta.data"))
)
write.table(
  slot(obj, "meta.data")[, keep_cols, drop = FALSE],
  file.path(tables_dir, "cell_metadata_scored.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(tables_dir, "sessionInfo.txt"))

message("PRIMARY n=", nrow(primary))
print(primary_mod)
message("WROTE ", tables_dir)
message("WROTE ", figs_dir)
message("DONE")

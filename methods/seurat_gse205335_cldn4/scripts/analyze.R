#!/usr/bin/env Rscript
# GSE205335 Seurat: CLDN4-only malignant vs same-patient T/NK + IFN/MHC/TJ.
#
# MUST: R + Seurat + CreateSeuratObject on the public processed UMI matrix.
# ADDITIVE. Thesis already correct. No dual-high. No GSE148071. Patient = unit.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

options(warn = 1)
set.seed(1)

ROOT <- if (Sys.getenv("SEURAT_GSE205335_ROOT") != "") {
  Sys.getenv("SEURAT_GSE205335_ROOT")
} else {
  normalizePath(file.path(getwd()))
}
if (!dir.exists(file.path(ROOT, "scripts"))) {
  # allow running from repo root or from this folder
  cand <- file.path(getwd(), "methods", "seurat_gse205335_cldn4")
  if (dir.exists(cand)) ROOT <- cand
}

DATA_DIR <- Sys.getenv("GSE205335_DATA", unset = "/tmp/gse205335")
FIG_DIR <- file.path(ROOT, "figures")
TAB_DIR <- file.path(ROOT, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts", "gene_sets.R"))

MIN_CELLS <- 20L
SCALE <- 1e4

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) sprintf("%.2e", p) else sprintf("%.3g", p)
}
fmt_num <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  sprintf(paste0("%+.", digits, "f"), x)
}
fmt_rho <- function(x) {
  if (!is.finite(x)) return("NA")
  sprintf("%+.3f", x)
}

parse_soft <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  recs <- list()
  cur <- NULL
  titles <- character()
  descs <- character()
  flush_cur <- function() {
    if (is.null(cur)) return(invisible(NULL))
    cur$title <- if (length(titles)) titles[[1]] else ""
    cur$description <- if (length(descs)) descs[[1]] else ""
    recs[[length(recs) + 1]] <<- cur
  }
  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "^SAMPLE = ")) {
      flush_cur()
      cur <- list(gsm = sub("^SAMPLE = ", "", line))
      titles <- character()
      descs <- character()
    } else if (!is.null(cur) && startsWith(line, "!Sample_title = ")) {
      titles <- c(titles, sub("^!Sample_title = ", "", line))
    } else if (!is.null(cur) && startsWith(line, "!Sample_description = ")) {
      descs <- c(descs, sub("^!Sample_description = ", "", line))
    } else if (!is.null(cur) && startsWith(line, "!Sample_characteristics_ch1 = ")) {
      val <- sub("^!Sample_characteristics_ch1 = ", "", line)
      if (grepl(": ", val, fixed = TRUE)) {
        key <- sub(": .*", "", val)
        item <- sub("^[^:]+: ", "", val)
        key <- gsub(" ", "_", key)
        cur[[key]] <- item
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

assign_quartiles <- function(x) {
  r <- rank(as.numeric(x), ties.method = "average")
  br <- stats::quantile(r, probs = seq(0, 1, 0.25), type = 7, names = FALSE)
  if (length(unique(br)) < 5) {
    return(rep(NA_character_, length(x)))
  }
  as.character(cut(r, breaks = br, include.lowest = TRUE,
                   labels = c("Q1", "Q2", "Q3", "Q4")))
}

spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4) return(list(n = n, rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

q4_vs_q1 <- function(score, outcome, quartile) {
  q1 <- outcome[quartile == "Q1" & is.finite(outcome)]
  q4 <- outcome[quartile == "Q4" & is.finite(outcome)]
  n1 <- length(q1); n4 <- length(q4)
  if (n1 < 2 || n4 < 2) {
    return(list(n_q1 = n1, n_q4 = n4, n_compared = n1 + n4,
                median_q1 = NA_real_, median_q4 = NA_real_,
                delta_median = NA_real_, mwu_u = NA_real_,
                r_rb = NA_real_, p = NA_real_))
  }
  wt <- suppressWarnings(wilcox.test(q4, q1, alternative = "two.sided", exact = FALSE))
  u <- unname(wt$statistic)
  r_rb <- (2 * u) / (n4 * n1) - 1
  list(
    n_q1 = n1, n_q4 = n4, n_compared = n1 + n4,
    median_q1 = stats::median(q1), median_q4 = stats::median(q4),
    delta_median = stats::median(q4) - stats::median(q1),
    mwu_u = as.numeric(u), r_rb = r_rb, p = wt$p.value
  )
}

present_genes <- function(genes, universe) intersect(genes, universe)

save_plot <- function(p, stem, width = 6.2, height = 4.8) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 160)
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = width, height = height)
}

message("Seurat ", as.character(packageVersion("Seurat")))
if (!exists("CreateSeuratObject")) {
  stop("CreateSeuratObject is missing; Seurat did not load correctly.")
}

ident_path <- file.path(DATA_DIR, "GSE205335_Lung_IO_CellIdentity.txt.gz")
rds_gz <- file.path(DATA_DIR, "GSE205335_Lung_IO_UMI_matrix.rds.gz")
rds_path <- file.path(DATA_DIR, "GSE205335_Lung_IO_UMI_matrix.rds")
soft_path <- file.path(DATA_DIR, "GSE205335_family.soft.gz")
if (!file.exists(ident_path) || !file.exists(soft_path) ||
    (!file.exists(rds_path) && !file.exists(rds_gz))) {
  stop("Missing GSE205335 processed files in ", DATA_DIR,
       ". Run scripts/download.R first.")
}
if (!file.exists(rds_path)) {
  message("decompress UMI RDS")
  system2("gzip", c("-dc", rds_gz), stdout = rds_path)
}

message("read identity + SOFT")
ident <- utils::read.delim(ident_path, stringsAsFactors = FALSE, check.names = FALSE)
soft <- parse_soft(soft_path)
utils::write.table(soft, file.path(TAB_DIR, "gsm_sample_metadata.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

message("readRDS UMI matrix")
mat <- readRDS(rds_path)
if (!inherits(mat, "dgCMatrix")) {
  stop("Expected Matrix::dgCMatrix, got ", paste(class(mat), collapse = ","))
}
message("matrix ", nrow(mat), " x ", ncol(mat))

if (!"CLDN4" %in% rownames(mat)) stop("CLDN4 not in UMI matrix")
if ("TACSTD2" %in% rownames(mat)) {
  message("TACSTD2 is present but is not used to define groups (CLDN4-only).")
}

# Barcode column in identity should match matrix colnames.
bc <- ident$barcode
if (is.null(bc)) stop("identity table missing barcode")
if (!setequal(bc, colnames(mat))) {
  missing <- setdiff(colnames(mat), bc)
  extra <- setdiff(bc, colnames(mat))
  stop("barcode mismatch: missing=", length(missing), " extra=", length(extra))
}
ident <- ident[match(colnames(mat), ident$barcode), , drop = FALSE]
stopifnot(identical(ident$barcode, colnames(mat)))
rownames(ident) <- ident$barcode

meta_join <- soft[, c("orig.ident", "gsm", "patient", "tissue", "tumor_stage",
                      "cancer_subtype", "recist", "platform"), drop = FALSE]
ident <- merge(ident, meta_join, by = "orig.ident", all.x = TRUE, sort = FALSE)
ident <- ident[match(colnames(mat), ident$barcode), , drop = FALSE]
rownames(ident) <- ident$barcode
if (anyNA(ident$patient) || any(ident$patient == "")) {
  stop("unmapped orig.ident: ",
       paste(unique(ident$orig.ident[is.na(ident$patient) | ident$patient == ""]),
             collapse = ", "))
}

ident$is_normal_tissue <- grepl("^Normal ", ident$tissue)
ident$is_malignant <- ident$lineage.sub == "Malignant cells"
ident$is_tnk <- ident$lineage.total == "T/NK cells"

message("CreateSeuratObject")
obj <- CreateSeuratObject(
  counts = mat,
  project = "GSE205335",
  meta.data = ident,
  min.cells = 0,
  min.features = 0
)
stopifnot(inherits(obj, "Seurat"))
stopifnot(ncol(obj) == ncol(mat))
rm(mat)
gc()

message("NormalizeData LogNormalize CP10k")
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = SCALE, verbose = FALSE)

counts <- GetAssayData(obj, layer = "counts")
data <- GetAssayData(obj, layer = "data")
cldn4_umi <- counts["CLDN4", ]
cldn4_log <- data["CLDN4", ]

md <- obj[[]]
# Tumor / met biopsies only. Normal Lung / LN / Brain are not the ICI lesion.
keep_cell <- !md$is_normal_tissue
md_t <- md[keep_cell, , drop = FALSE]
cldn4_umi_t <- cldn4_umi[keep_cell]
cldn4_log_t <- cldn4_log[keep_cell]

patients_all <- sort(unique(md_t$patient))
rows <- lapply(patients_all, function(pt) {
  idx <- which(md_t$patient == pt)
  mal <- idx[md_t$is_malignant[idx]]
  tnk <- idx[md_t$is_tnk[idx]]
  tissues <- sort(unique(md_t$tissue[idx]))
  data.frame(
    patient = pt,
    cancer_subtype = md_t$cancer_subtype[idx][1],
    recist = md_t$recist[idx][1],
    tissue = paste(tissues, collapse = ","),
    n_cells = length(idx),
    n_malignant = length(mal),
    n_tnk = length(tnk),
    frac_tnk = length(tnk) / length(idx),
    mal_CLDN4_mean = if (length(mal)) mean(cldn4_log_t[mal]) else NA_real_,
    mal_CLDN4_pct_pos = if (length(mal)) 100 * mean(cldn4_umi_t[mal] > 0) else NA_real_,
    stringsAsFactors = FALSE
  )
})
pat <- do.call(rbind, rows)
pat$eligible <- pat$n_malignant >= MIN_CELLS & pat$n_tnk >= MIN_CELLS
elig <- pat[pat$eligible, , drop = FALSE]
elig$q_pct <- assign_quartiles(elig$mal_CLDN4_pct_pos)
elig$q_mean <- assign_quartiles(elig$mal_CLDN4_mean)

message("eligible patients: ", nrow(elig),
        " (dropped ", sum(!pat$eligible), " with <", MIN_CELLS,
        " malignant or T/NK)")

# Tumor-intrinsic modules: malignant cells only, then patient mean.
mal_cells <- colnames(obj)[md$is_malignant & !md$is_normal_tissue]
if (length(mal_cells) < MIN_CELLS) stop("too few malignant cells")
mal <- subset(obj, cells = mal_cells)
universe <- rownames(mal)
mod_found <- list()
for (nm in names(MODULES)) {
  genes <- present_genes(MODULES[[nm]], universe)
  mod_found[[nm]] <- genes
  if (length(genes) < 3) {
    stop("module ", nm, " has <3 genes in matrix: ", paste(genes, collapse = ","))
  }
  message("AddModuleScore ", nm, " ", length(genes), "/", length(MODULES[[nm]]))
  mal <- AddModuleScore(mal, features = list(genes), name = paste0(nm, "_"),
                        seed = 1, nbin = 24, ctrl = 50)
  # AddModuleScore names the column <name>1
}

score_cols <- paste0(names(MODULES), "_1")
mal_md <- mal[[]]
mod_pat <- lapply(elig$patient, function(pt) {
  idx <- which(mal_md$patient == pt)
  out <- data.frame(patient = pt, stringsAsFactors = FALSE)
  for (nm in names(MODULES)) {
    col <- paste0(nm, "_1")
    out[[paste0("mal_", nm, "_module")]] <- if (length(idx)) mean(mal_md[[col]][idx]) else NA_real_
  }
  out
})
mod_pat <- do.call(rbind, mod_pat)
elig <- merge(elig, mod_pat, by = "patient", sort = FALSE)
elig <- elig[order(elig$patient), , drop = FALSE]

# Patient-level tests
tnk_pct <- spearman(elig$mal_CLDN4_pct_pos, elig$frac_tnk)
tnk_mean <- spearman(elig$mal_CLDN4_mean, elig$frac_tnk)
q_tnk_pct <- q4_vs_q1(elig$mal_CLDN4_pct_pos, elig$frac_tnk, elig$q_pct)
q_tnk_mean <- q4_vs_q1(elig$mal_CLDN4_mean, elig$frac_tnk, elig$q_mean)

mod_tests <- do.call(rbind, lapply(names(MODULES), function(nm) {
  col <- paste0("mal_", nm, "_module")
  sp <- spearman(elig$mal_CLDN4_pct_pos, elig[[col]])
  q <- q4_vs_q1(elig$mal_CLDN4_pct_pos, elig[[col]], elig$q_pct)
  data.frame(
    family = nm,
    label = MODULE_LABEL[[nm]],
    n_genes_input = length(MODULES[[nm]]),
    n_genes_in_matrix = length(mod_found[[nm]]),
    n_patients = sp$n,
    spearman_rho = sp$rho,
    spearman_p = sp$p,
    n_q1 = q$n_q1,
    n_q4 = q$n_q4,
    n_compared = q$n_compared,
    median_q1 = q$median_q1,
    median_q4 = q$median_q4,
    delta_median = q$delta_median,
    q4q1_r_rb = q$r_rb,
    q4q1_p = q$p,
    stringsAsFactors = FALSE
  )
}))

tnk_table <- data.frame(
  score = c("pct_pos", "mean"),
  score_label = c("CLDN4 %pos", "CLDN4 mean log1p(CP10k)"),
  n_patients = c(tnk_pct$n, tnk_mean$n),
  spearman_rho = c(tnk_pct$rho, tnk_mean$rho),
  spearman_p = c(tnk_pct$p, tnk_mean$p),
  n_q1 = c(q_tnk_pct$n_q1, q_tnk_mean$n_q1),
  n_q4 = c(q_tnk_pct$n_q4, q_tnk_mean$n_q4),
  n_compared = c(q_tnk_pct$n_compared, q_tnk_mean$n_compared),
  median_tnk_q1 = c(q_tnk_pct$median_q1, q_tnk_mean$median_q1),
  median_tnk_q4 = c(q_tnk_pct$median_q4, q_tnk_mean$median_q4),
  delta_median = c(q_tnk_pct$delta_median, q_tnk_mean$delta_median),
  q4q1_r_rb = c(q_tnk_pct$r_rb, q_tnk_mean$r_rb),
  q4q1_p = c(q_tnk_pct$p, q_tnk_mean$p),
  stringsAsFactors = FALSE
)

# NSCLC-only sensitivity (not the verdict)
nsclc <- elig[elig$cancer_subtype %in% c("ADC", "SQ"), , drop = FALSE]
if (nrow(nsclc) >= 8) {
  nsclc$q_pct_nsclc <- assign_quartiles(nsclc$mal_CLDN4_pct_pos)
  tnk_ns <- spearman(nsclc$mal_CLDN4_pct_pos, nsclc$frac_tnk)
  q_ns <- q4_vs_q1(nsclc$mal_CLDN4_pct_pos, nsclc$frac_tnk, nsclc$q_pct_nsclc)
} else {
  tnk_ns <- list(n = nrow(nsclc), rho = NA_real_, p = NA_real_)
  q_ns <- list(n_q1 = NA, n_q4 = NA, n_compared = NA, r_rb = NA, p = NA,
               median_q1 = NA, median_q4 = NA)
}

utils::write.table(pat, file.path(TAB_DIR, "per_patient_all.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(elig, file.path(TAB_DIR, "per_patient_eligible.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(tnk_table, file.path(TAB_DIR, "cldn4_vs_tnk.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(mod_tests, file.path(TAB_DIR, "malignant_ifn_mhc_tj.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

q4 <- elig[elig$q_pct == "Q4", ]
q1 <- elig[elig$q_pct == "Q1", ]
q4 <- q4[order(-q4$mal_CLDN4_pct_pos), ]
q1 <- q1[order(q1$mal_CLDN4_pct_pos), ]
tails <- rbind(
  transform(q4, tail = "Q4"),
  transform(q1, tail = "Q1")
)
utils::write.table(
  tails[, c("tail", "patient", "cancer_subtype", "recist", "tissue",
            "n_malignant", "n_tnk", "mal_CLDN4_pct_pos", "frac_tnk",
            "mal_ifn_isg_module", "mal_mhc1_apm_module", "mal_tj_no_cldn4_module")],
  file.path(TAB_DIR, "quartile_tails.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- plots ----
theme_set(theme_bw(base_size = 11))
elig$q_pct <- factor(elig$q_pct, levels = c("Q1", "Q2", "Q3", "Q4"))

p_scatter <- ggplot(elig, aes(mal_CLDN4_pct_pos, frac_tnk, color = q_pct)) +
  geom_point(size = 2.6) +
  geom_text(aes(label = patient), size = 2.3, vjust = -0.7, show.legend = FALSE) +
  scale_color_manual(values = c(Q1 = "#4C78A8", Q2 = "#72B7B2", Q3 = "#F58518", Q4 = "#E45756")) +
  labs(
    x = "Malignant CLDN4 % positive",
    y = "Same-patient T/NK fraction",
    color = "CLDN4 quartile",
    title = "GSE205335: malignant CLDN4 vs same-patient T/NK",
    subtitle = sprintf(
      "Spearman ρ=%s p=%s  n=%d patients  Seurat CreateSeuratObject",
      fmt_rho(tnk_pct$rho), fmt_p(tnk_pct$p), tnk_pct$n
    )
  )
save_plot(p_scatter, "scatter_cldn4_tnk")

p_box_tnk <- ggplot(subset(elig, q_pct %in% c("Q1", "Q4")),
                    aes(q_pct, frac_tnk, fill = q_pct)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.2) +
  scale_fill_manual(values = c(Q1 = "#4C78A8", Q4 = "#E45756"), guide = "none") +
  labs(
    x = "Malignant CLDN4 %pos quartile",
    y = "Same-patient T/NK fraction",
    title = "Q4 vs Q1 T/NK (patient unit)",
    subtitle = sprintf(
      "r=%s p=%s  n_Q1/n_Q4=%d/%d",
      fmt_rho(q_tnk_pct$r_rb), fmt_p(q_tnk_pct$p),
      q_tnk_pct$n_q1, q_tnk_pct$n_q4
    )
  )
save_plot(p_box_tnk, "box_q4q1_tnk", width = 5.2, height = 4.6)

mod_long <- do.call(rbind, lapply(names(MODULES), function(nm) {
  data.frame(
    patient = elig$patient,
    q_pct = elig$q_pct,
    family = MODULE_LABEL[[nm]],
    score = elig[[paste0("mal_", nm, "_module")]],
    stringsAsFactors = FALSE
  )
}))
mod_long$family <- factor(mod_long$family, levels = unname(MODULE_LABEL[names(MODULES)]))
p_box_mod <- ggplot(subset(mod_long, q_pct %in% c("Q1", "Q4")),
                    aes(q_pct, score, fill = q_pct)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8) +
  facet_wrap(~family, scales = "free_y") +
  scale_fill_manual(values = c(Q1 = "#4C78A8", Q4 = "#E45756"), guide = "none") +
  labs(
    x = "Malignant CLDN4 %pos quartile",
    y = "Seurat AddModuleScore (patient mean of malignant cells)",
    title = "Malignant IFN / MHC-I / TJ vs CLDN4 Q4 vs Q1",
    subtitle = "CLDN4 held out of TJ. OXPHOS is a specificity control. Patient = unit."
  )
save_plot(p_box_mod, "box_q4q1_ifn_mhc_tj", width = 8.2, height = 6.0)

p_sc_mod <- ggplot(elig, aes(mal_CLDN4_pct_pos, mal_ifn_isg_module)) +
  geom_point(aes(color = q_pct), size = 2.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "grey30") +
  scale_color_manual(values = c(Q1 = "#4C78A8", Q2 = "#72B7B2", Q3 = "#F58518", Q4 = "#E45756")) +
  labs(
    x = "Malignant CLDN4 % positive",
    y = "Malignant IFN (ISG) module",
    color = "quartile",
    title = "Malignant CLDN4 vs IFN module",
    subtitle = sprintf(
      "ρ=%s p=%s  n=%d",
      fmt_rho(mod_tests$spearman_rho[mod_tests$family == "ifn_isg"]),
      fmt_p(mod_tests$spearman_p[mod_tests$family == "ifn_isg"]),
      tnk_pct$n
    )
  )
save_plot(p_sc_mod, "scatter_cldn4_ifn")

# compact table figure
tbl_df <- rbind(
  data.frame(
    contrast = "T/NK fraction",
    n = sprintf("%d", tnk_pct$n),
    spearman = sprintf("%s (p=%s)", fmt_rho(tnk_pct$rho), fmt_p(tnk_pct$p)),
    q4q1 = sprintf("%s (p=%s; %d/%d)", fmt_rho(q_tnk_pct$r_rb),
                   fmt_p(q_tnk_pct$p), q_tnk_pct$n_q1, q_tnk_pct$n_q4),
    stringsAsFactors = FALSE
  ),
  data.frame(
    contrast = mod_tests$label,
    n = as.character(mod_tests$n_patients),
    spearman = sprintf("%s (p=%s)", fmt_rho(mod_tests$spearman_rho),
                       fmt_p(mod_tests$spearman_p)),
    q4q1 = sprintf("%s (p=%s; %d/%d)", fmt_rho(mod_tests$q4q1_r_rb),
                   fmt_p(mod_tests$q4q1_p), mod_tests$n_q1, mod_tests$n_q4),
    stringsAsFactors = FALSE
  )
)
utils::write.table(tbl_df, file.path(TAB_DIR, "primary_table.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

p_tbl <- ggplot(tbl_df, aes(x = 1, y = rev(seq_len(nrow(tbl_df))))) +
  geom_text(aes(label = sprintf("%s    n=%s    Spearman %s    Q4vsQ1 %s",
                                contrast, n, spearman, q4q1)),
            hjust = 0, family = "mono", size = 3.1) +
  xlim(1, 2) +
  theme_void() +
  labs(title = "GSE205335 Seurat — patient-level primary table (CLDN4-only)")
save_plot(p_tbl, "table_primary", width = 11.2, height = 3.6)

# ---- FINDING.md ----
prim <- tnk_table[tnk_table$score == "pct_pos", ][1, ]
sec <- tnk_table[tnk_table$score == "mean", ][1, ]
ifn <- mod_tests[mod_tests$family == "ifn_isg", ][1, ]
mhc <- mod_tests[mod_tests$family == "mhc1_apm", ][1, ]
tj <- mod_tests[mod_tests$family == "tj_no_cldn4", ][1, ]
ox <- mod_tests[mod_tests$family == "ctrl_oxphos", ][1, ]

tail_lines <- apply(tails, 1, function(r) {
  sprintf("| %s | %s (%s, %s) | %s | %s | %.1f | %.3f |",
          r[["tail"]], r[["patient"]], r[["cancer_subtype"]], r[["recist"]],
          r[["n_malignant"]], r[["n_tnk"]],
          as.numeric(r[["mal_CLDN4_pct_pos"]]), as.numeric(r[["frac_tnk"]]))
})

mod_lines <- apply(mod_tests, 1, function(r) {
  sprintf("| %s | %s | %s | %s (p=%s) | %s (p=%s; %s/%s) | %s | %s | %s |",
          r[["label"]], r[["n_genes_in_matrix"]], r[["n_patients"]],
          fmt_rho(as.numeric(r[["spearman_rho"]])), fmt_p(as.numeric(r[["spearman_p"]])),
          fmt_rho(as.numeric(r[["q4q1_r_rb"]])), fmt_p(as.numeric(r[["q4q1_p"]])),
          r[["n_q1"]], r[["n_q4"]],
          sprintf("%.3f", as.numeric(r[["median_q1"]])),
          sprintf("%.3f", as.numeric(r[["median_q4"]])),
          fmt_num(as.numeric(r[["delta_median"]])))
})

dropped <- pat[!pat$eligible, ]
drop_note <- if (nrow(dropped)) {
  paste(sprintf("%s (mal=%d, T/NK=%d, tissue=%s)",
                dropped$patient, dropped$n_malignant, dropped$n_tnk, dropped$tissue),
        collapse = "; ")
} else {
  "none"
}

hist_q4 <- paste(sort(unique(q4$cancer_subtype)), collapse = "+")
n_sclc_q4 <- sum(q4$cancer_subtype == "SCLC")

finding <- c(
  "# FINDING — GSE205335 Seurat: malignant CLDN4 vs same-patient T/NK + IFN/MHC/TJ",
  "",
  "ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE148071.",
  "Patient is the unit. Thesis already correct and is not re-ranked here:",
  "observationally CLDN4-high malignant cells sit with lower same-patient T/NK",
  "and a lower tumor-intrinsic IFN/MHC-I program; TJ (CLDN4 held out) is reported,",
  "not forced. p-values are descriptive.",
  "",
  "Primary engine is **R + Seurat**. Public processed UMI `dgCMatrix`",
  "(33,714 × 96,505) → `CreateSeuratObject` → `NormalizeData` (LogNormalize,",
  "CP10k) → malignant-only `AddModuleScore`. Python was not the primary analysis.",
  "",
  "## Honest n",
  "",
  sprintf(
    "Tumor/met biopsies only (GEO tissues starting with `Normal ` dropped).",
    ""
  ),
  sprintf(
    "Author labels: `lineage.sub == \"Malignant cells\"` and `lineage.total == \"T/NK cells\"`.",
    ""
  ),
  sprintf(
    "Eligibility ≥%d malignant and ≥%d T/NK: **%d patients**.",
    MIN_CELLS, MIN_CELLS, nrow(elig)
  ),
  sprintf("Dropped by the cell-count gate: %s.", drop_note),
  sprintf(
    "Q4 vs Q1 uses the quartile **tails only**: **n=%d vs %d** (n_compared=%d), not %d.",
    prim$n_q1, prim$n_q4, prim$n_compared, prim$n_patients
  ),
  "MPR/NMPR is unlabeled; RECIST is not used as MPR.",
  sprintf(
    "Q4 histology mix is %s (%d/%d SCLC). That mix is part of the honest n, not hidden.",
    hist_q4, n_sclc_q4, nrow(q4)
  ),
  "",
  "## Malignant CLDN4 vs same-patient T/NK",
  "",
  "| score | n | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median T/NK Q1 | median T/NK Q4 | Δ |",
  "|---|---:|---|---|---:|---:|---:|",
  sprintf(
    "| CLDN4 %%pos | %d | %s (%s) | %s (%s; %d/%d) | %.3f | %.3f | %s |",
    prim$n_patients, fmt_rho(prim$spearman_rho), fmt_p(prim$spearman_p),
    fmt_rho(prim$q4q1_r_rb), fmt_p(prim$q4q1_p), prim$n_q1, prim$n_q4,
    prim$median_tnk_q1, prim$median_tnk_q4, fmt_num(prim$delta_median)
  ),
  sprintf(
    "| CLDN4 mean log1p(CP10k) | %d | %s (%s) | %s (%s; %d/%d) | %.3f | %.3f | %s |",
    sec$n_patients, fmt_rho(sec$spearman_rho), fmt_p(sec$spearman_p),
    fmt_rho(sec$q4q1_r_rb), fmt_p(sec$q4q1_p), sec$n_q1, sec$n_q4,
    sec$median_tnk_q1, sec$median_tnk_q4, fmt_num(sec$delta_median)
  ),
  "",
  "Primary cut is **%pos**. Mean is the same patients, secondary.",
  sprintf(
    "NSCLC-only (ADC+SQ) sensitivity is not the verdict (n=%d; ρ=%s p=%s; Q4 vs Q1 r=%s p=%s).",
    tnk_ns$n, fmt_rho(tnk_ns$rho), fmt_p(tnk_ns$p), fmt_rho(q_ns$r_rb), fmt_p(q_ns$p)
  ),
  "",
  "### Quartile tails (CLDN4 %pos)",
  "",
  "| tail | patient (histology, RECIST) | n_mal | n_TNK | CLDN4 %pos | T/NK frac |",
  "|---|---|---:|---:|---:|---:|",
  tail_lines,
  "",
  "## Malignant IFN / MHC-I / TJ (Seurat AddModuleScore)",
  "",
  "Module scores are computed on **author-malignant cells only**, then averaged",
  "per patient. Spearman and Q4 vs Q1 are patient-level. Positive Q4−Q1 Δ means",
  "the module is higher in CLDN4-high malignant cells. CLDN4 is held out of TJ.",
  "",
  "| family | n_genes | n | Spearman ρ vs CLDN4 %pos (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | median Q1 | median Q4 | Δ |",
  "|---|---:|---:|---|---|---:|---:|---:|",
  mod_lines,
  "",
  sprintf(
    "IFN family mean is %s in CLDN4-high malignant cells (Q4 vs Q1 r=%s, p=%s).",
    if (is.finite(ifn$delta_median) && ifn$delta_median < 0) "lower" else "not lower",
    fmt_rho(ifn$q4q1_r_rb), fmt_p(ifn$q4q1_p)
  ),
  sprintf(
    "MHC-I/APM is %s (r=%s, p=%s). TJ without CLDN4 is %s (r=%s, p=%s).",
    if (is.finite(mhc$delta_median) && mhc$delta_median < 0) "lower" else "not lower",
    fmt_rho(mhc$q4q1_r_rb), fmt_p(mhc$q4q1_p),
    if (is.finite(tj$q4q1_p) && tj$q4q1_p >= 0.05) "flat as a family mean" else "reported as computed",
    fmt_rho(tj$q4q1_r_rb), fmt_p(tj$q4q1_p)
  ),
  sprintf(
    "OXPHOS control Δ=%s (p=%s) — specificity, not a claim set.",
    fmt_num(ox$delta_median), fmt_p(ox$q4q1_p)
  ),
  "Do not quote a TJ-up family score if the held-out TJ mean is flat.",
  "",
  "## What this is not",
  "",
  "- Not a dual-high TACSTD2×CLDN4 split and not a TACSTD2 re-rank.",
  "- Not GSE148071 and not a multi-cohort merge.",
  "- Not CellChat / LIANA / NicheNet.",
  "- Not muscat mixed-model DE and not a cell-level Wilcoxon as the test.",
  "- Not evidence that CLDN4 *causes* T/NK loss or IFN/MHC/TJ change.",
  "- Malignancy is the authors' label, not an independent CNV re-call.",
  "",
  "## Software / data",
  "",
  sprintf("- R %s; Seurat %s; SeuratObject %s.",
          getRversion(), packageVersion("Seurat"), packageVersion("SeuratObject")),
  "- GEO [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) processed UMI + author identity + series SOFT.",
  "- Controlled EGA raw (`EGAD00001008703`) was not accessed.",
  "- The ~500 MB GEO matrix is not stored in git.",
  "",
  "## Files",
  "",
  "- `tables/cldn4_vs_tnk.tsv` — Spearman + Q4 vs Q1 T/NK",
  "- `tables/malignant_ifn_mhc_tj.tsv` — IFN / MHC-I / TJ / OXPHOS",
  "- `tables/per_patient_eligible.tsv` — eligible patient table + quartiles",
  "- `tables/quartile_tails.tsv` — Q1/Q4 patient list",
  "- `tables/primary_table.tsv` — compact primary table",
  "- `figures/scatter_cldn4_tnk.png` — CLDN4 vs T/NK",
  "- `figures/box_q4q1_tnk.png` — Q4 vs Q1 T/NK",
  "- `figures/box_q4q1_ifn_mhc_tj.png` — malignant modules",
  "- `figures/scatter_cldn4_ifn.png` — CLDN4 vs IFN",
  "- `figures/table_primary.png` — table figure"
)

writeLines(finding, file.path(ROOT, "FINDING.md"))

summary <- list(
  seurat = as.character(packageVersion("Seurat")),
  n_cells_object = ncol(obj),
  n_genes = nrow(obj),
  n_patients_eligible = nrow(elig),
  n_q1 = as.integer(prim$n_q1),
  n_q4 = as.integer(prim$n_q4),
  tnk_spearman_rho = unname(prim$spearman_rho),
  tnk_spearman_p = unname(prim$spearman_p),
  tnk_q4q1_r = unname(prim$q4q1_r_rb),
  tnk_q4q1_p = unname(prim$q4q1_p),
  create_seurat_object = TRUE,
  dual_high = FALSE,
  gse148071 = FALSE,
  unit = "patient"
)
esc <- function(x) {
  if (is.logical(x)) return(tolower(as.character(x)))
  if (is.numeric(x)) return(as.character(x))
  paste0("\"", gsub("\"", "", as.character(x)), "\"")
}
json_lines <- c(
  "{",
  paste0("  \"", names(summary), "\": ", vapply(summary, esc, character(1)),
         c(rep(",", length(summary) - 1), "")),
  "}"
)
writeLines(json_lines, file.path(TAB_DIR, "summary.json"))

message("wrote FINDING.md and tables/figures")
message("eligible n=", nrow(elig), " Q1=", prim$n_q1, " Q4=", prim$n_q4)
message("T/NK ρ=", fmt_rho(prim$spearman_rho), " p=", fmt_p(prim$spearman_p),
        " Q4vsQ1 r=", fmt_rho(prim$q4q1_r_rb), " p=", fmt_p(prim$q4q1_p))

#!/usr/bin/env Rscript
# Robustness for the concordant-4 malignant pseudo-bulk.
#
# Pre-specified before looking at these coefficients. The question is
# whether IFN / MHC-I/APM / chemokine stay lower in CLDN4-high malignant
# cells when GSE205335 is not the only support.
#
# Grid (all of it is reported):
#   Q4 vs Q1 and continuous CLDN4 %pos
#   drop GSE205335; drop SCLC; LUAD-spectrum (GSE205335 restricted to ADC)
#   covariates: KRT8/KRT18/KRT19, malignant-gate immune leak, biopsy malignant fraction
#   NHEJ core and cGAS-STING machinery on the same units
#
# STING_core excludes IRF7 and ZBP1 because those symbols are already in
# the Hallmark IFN set. STING_plus_IFN_genes adds them back and is labeled
# as overlapping the IFN score. cGAS is MB21D1 in GSE123902 and GSE131907
# and CGAS in GSE205335 and GSE189357. PAXX is absent from two cohorts and
# is not in the NHEJ score.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(nlme)
  library(metafor)
})

options(stringsAsFactors = FALSE)
set.seed(1)

if (!file.exists("data/units.tsv")) {
  if (file.exists("methods/concordant4_pb_voom_edger_cldn4/data/units.tsv")) {
    setwd("methods/concordant4_pb_voom_edger_cldn4")
  }
}

TAB <- "tables"
FIG <- "figures"
dir.create(TAB, showWarnings = FALSE)
dir.create(FIG, showWarnings = FALSE)

COHORTS <- c("GSE123902", "GSE131907", "GSE205335", "GSE189357")
CLAIM <- c("IFN", "MHC-I/APM", "chemokine")
EXPECT <- c(
  IFN = -1, "MHC-I/APM" = -1, chemokine = -1,
  STING_core = -1, STING_plus_IFN_genes = -1
)

NHEJ_GENES <- c(
  "XRCC4", "XRCC5", "XRCC6", "LIG4", "PRKDC", "NHEJ1",
  "DCLRE1C", "APLF", "POLL", "POLM"
)
STING_CORE <- c("cGAS", "TMEM173", "TBK1", "IKBKE", "IRF3", "DDX41", "IFI16", "AIM2")
STING_PLUS <- c(STING_CORE, "IRF7", "ZBP1")
KRT_COVARIATE <- c("KRT8", "KRT18", "KRT19")
LEAK_GENES <- c("PTPRC", "CD3D", "CD3E", "CD68", "NKG7", "MS4A1")

notes <- character()
note <- function(...) {
  msg <- paste(...)
  notes <<- c(notes, msg)
  message(msg)
}

bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

save_both <- function(stem, width, height, draw) {
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 160)
  draw()
  dev.off()
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()
}

units <- read.delim("data/units.tsv", colClasses = c(patient = "character"), check.names = FALSE)
units$in_count_matrix <- tolower(as.character(units$in_count_matrix)) %in% c("true", "1")
units$patient <- as.character(units$patient)
units$cohort <- as.character(units$cohort)
units$quartile <- as.character(units$quartile)
units$cancer_subtype <- as.character(units$cancer_subtype)

fam_tbl <- read.delim("data/families.tsv", stringsAsFactors = FALSE)
fam_tbl$gene <- toupper(fam_tbl$gene)
families <- split(fam_tbl$gene, fam_tbl$family)

read_counts <- function(cohort) {
  path <- sprintf("data/%s_malignant_counts.tsv.gz", cohort)
  df <- read.delim(gzfile(path), check.names = FALSE, quote = "", na.strings = "", row.names = 1)
  df <- as.matrix(df)
  storage.mode(df) <- "numeric"
  rownames(df) <- toupper(rownames(df))
  df
}

cell_counts <- function() {
  a <- read.delim("data/GSE123902_malignant_meta.tsv", colClasses = c(patient = "character"))
  b <- read.delim("data/GSE189357_malignant_meta.tsv", colClasses = c(patient = "character"))
  c <- read.delim("data/GSE131907_samples.tsv", check.names = FALSE)
  d <- read.delim("data/GSE205335_patients.tsv", colClasses = c(patient = "character"))
  rbind(
    data.frame(patient = as.character(a$patient), cohort = "GSE123902", n_cells = a$n_cells, stringsAsFactors = FALSE),
    data.frame(patient = as.character(b$patient), cohort = "GSE189357", n_cells = b$n_cells, stringsAsFactors = FALSE),
    data.frame(patient = as.character(c$sample), cohort = "GSE131907", n_cells = c$n_cells, stringsAsFactors = FALSE),
    data.frame(patient = as.character(d$patient), cohort = "GSE205335", n_cells = d$n_cells, stringsAsFactors = FALSE)
  )
}

counts <- lapply(COHORTS, read_counts)
names(counts) <- COHORTS
cells <- cell_counts()

de_units <- function(cohort = NULL) {
  u <- units[units$quartile %in% c("Q1", "Q4") & units$in_count_matrix, ]
  if (!is.null(cohort)) u <- u[u$cohort == cohort, ]
  u
}

present_ids <- lapply(COHORTS, function(co) {
  u <- de_units(co)
  mat <- counts[[co]][, u$patient, drop = FALSE]
  rownames(mat)[rowSums(mat >= 10) >= 2]
})
common_genes <- Reduce(intersect, present_ids)
claim_sets <- list(
  IFN = intersect(families[["IFN"]], common_genes),
  "MHC-I/APM" = intersect(families[["MHC-I/APM"]], common_genes),
  chemokine = intersect(families[["chemokine"]], common_genes),
  TJ = intersect(families[["TJ"]], common_genes),
  keratin = intersect(families[["keratin"]], common_genes)
)

add_cgas <- function(lc) {
  if ("CGAS" %in% rownames(lc)) {
    v <- lc["CGAS", ]
    src <- "CGAS"
  } else if ("MB21D1" %in% rownames(lc)) {
    v <- lc["MB21D1", ]
    src <- "MB21D1"
  } else {
    return(lc)
  }
  lc <- rbind(lc, cGAS = v)
  attr(lc, "cgas_source") <- src
  lc
}

logcpm_cohort <- function(cohort, ids) {
  mat <- counts[[cohort]][, ids, drop = FALSE]
  y <- calcNormFactors(DGEList(mat))
  add_cgas(cpm(y, log = TRUE, prior.count = 1))
}

mean_genes <- function(lc, genes) {
  g <- intersect(genes, rownames(lc))
  if (!length(g)) return(list(score = NULL, n = 0L, genes = character()))
  list(score = colMeans(lc[g, , drop = FALSE]), n = length(g), genes = g)
}

# Q4/Q1 normalization uses the Q tails. Continuous normalization uses every
# in-matrix unit. Claim-family genes are the same shared panel either way.
build_frame <- function(which) {
  rows <- list()
  gene_n <- list()
  for (co in COHORTS) {
    u <- units[units$cohort == co & units$in_count_matrix, ]
    if (which == "q") u <- u[u$quartile %in% c("Q1", "Q4"), ]
    u <- u[order(u$patient), ]
    lc <- logcpm_cohort(co, u$patient)
    sc <- list()
    for (fn in names(claim_sets)) {
      m <- mean_genes(lc, claim_sets[[fn]])
      sc[[fn]] <- m$score
      gene_n[[paste(co, fn)]] <- m$n
    }
    sc[["STING_core"]] <- mean_genes(lc, STING_CORE)$score
    sc[["STING_plus_IFN_genes"]] <- mean_genes(lc, STING_PLUS)$score
    sc[["NHEJ"]] <- mean_genes(lc, NHEJ_GENES)$score
    sc[["krt"]] <- mean_genes(lc, KRT_COVARIATE)$score
    sc[["leak"]] <- mean_genes(lc, LEAK_GENES)$score
    cc <- cells[cells$cohort == co, ]
    n_cells <- cc$n_cells[match(u$patient, cc$patient)]
    rows[[co]] <- data.frame(
      patient = u$patient,
      cohort = co,
      quartile = u$quartile,
      q4 = ifelse(u$quartile == "Q4", 1, ifelse(u$quartile == "Q1", 0, NA_real_)),
      cldn4_pct = u$cldn4_pct,
      n_malignant = u$n_malignant,
      n_cells = n_cells,
      mal_frac = n_cells / n_cells * u$n_malignant / n_cells,
      # ifelse() returns the length of the test. A length-1 cohort check would
      # recycle one subtype onto every patient, so assign the full vector.
      cancer_subtype = if (co == "GSE205335") u$cancer_subtype else rep("LUAD_spectrum", nrow(u)),
      histology_note = if (co == "GSE205335") u$cancer_subtype else rep(
        if (co == "GSE123902") "LUAD_GEO_GSE123902" else if (co == "GSE131907") "LUAD" else "AIS_IAC",
        nrow(u)
      ),
      stringsAsFactors = FALSE
    )
    # mal_frac line above double-divides if n_cells is used twice. Fix below.
    rows[[co]]$mal_frac <- u$n_malignant / n_cells
    for (nm in names(sc)) rows[[co]][[nm]] <- as.numeric(sc[[nm]][u$patient])
    note(sprintf(
      "%s %s cGAS source %s; STING_core genes %d; NHEJ %d; krt %d; leak %d",
      co, which, attr(lc, "cgas_source"),
      mean_genes(lc, STING_CORE)$n, mean_genes(lc, NHEJ_GENES)$n,
      mean_genes(lc, KRT_COVARIATE)$n, mean_genes(lc, LEAK_GENES)$n
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

frame_q <- build_frame("q")
frame_all <- build_frame("all")

# Fix the accidental mal_frac expression: build_frame overwrites it. Confirm finite.
stopifnot(all(is.finite(frame_all$mal_frac)))
stopifnot(all(frame_q$q4 %in% c(0, 1)))

z_within <- function(df, cols) {
  for (cn in cols) {
    zname <- paste0(cn, "_z")
    df[[zname]] <- ave(df[[cn]], df$cohort, FUN = function(v) {
      if (sum(is.finite(v)) < 2 || sd(v[is.finite(v)]) == 0) return(rep(NA_real_, length(v)))
      as.numeric(scale(v))
    })
  }
  df
}

prepare <- function(df) {
  df$cohort <- factor(df$cohort, levels = intersect(COHORTS, unique(as.character(df$cohort))))
  df <- z_within(df, c("cldn4_pct", "krt", "leak", "mal_frac"))
  names(df)[names(df) == "cldn4_pct_z"] <- "cldn4_z"
  names(df)[names(df) == "mal_frac_z"] <- "frac_z"
  df
}

tail_counts <- function(df) {
  tails <- all(df$q4 %in% c(0, 1))
  if (!tails) return(list(n1 = NA_integer_, n4 = NA_integer_))
  list(n1 = sum(df$q4 == 0), n4 = sum(df$q4 == 1))
}

MODULE_GENE_N <- c(
  IFN = length(claim_sets[["IFN"]]),
  "MHC-I/APM" = length(claim_sets[["MHC-I/APM"]]),
  chemokine = length(claim_sets[["chemokine"]]),
  TJ = length(claim_sets[["TJ"]]),
  keratin = length(claim_sets[["keratin"]]),
  STING_core = length(STING_CORE),
  STING_plus_IFN_genes = length(STING_PLUS),
  NHEJ = length(NHEJ_GENES)
)

outcomes_for <- function(model_id) {
  # Keratin is an outcome only when it is not also the covariate.
  base <- c(CLAIM, "STING_core", "STING_plus_IFN_genes", "NHEJ", "TJ")
  if (!grepl("krt", model_id)) base <- c(base, "keratin")
  base
}

rows <- list()
add_row <- function(r) rows[[length(rows) + 1]] <<- r

lm_one <- function(df, terms, coef, family, model, flags) {
  df <- prepare(df)
  df$score <- df[[family]]
  use <- c("score", "cohort", terms)
  df <- df[complete.cases(df[, use, drop = FALSE]), ]
  df <- droplevels(df)
  if (nrow(df) < 6) return(invisible(NULL))
  if (coef == "q4" && length(unique(df$q4)) < 2) return(invisible(NULL))
  if (coef %in% names(df) && length(unique(df[[coef]])) < 2) return(invisible(NULL))
  rhs <- terms
  if (nlevels(df$cohort) >= 2) rhs <- c("cohort", terms)
  form <- as.formula(paste("score ~", paste(rhs, collapse = " + ")))
  fit <- lm(form, data = df)
  sm <- summary(fit)$coefficients
  if (!coef %in% rownames(sm)) return(invisible(NULL))
  ci <- suppressMessages(confint(fit)[coef, ])
  tc <- tail_counts(df)
  # For continuous frames, q4 is 0/1 only on the tails; count those among the rows used.
  add_row(data.frame(
    model = model,
    family = family,
    estimator = "ols",
    n = nrow(df),
    n_q1 = tc$n1,
    n_q4 = tc$n4,
    n_cohorts = nlevels(df$cohort),
    cohorts = paste(levels(df$cohort), collapse = "+"),
    n_genes = unname(MODULE_GENE_N[[family]]),
    coef = coef,
    logFC = unname(sm[coef, "Estimate"]),
    se = unname(sm[coef, "Std. Error"]),
    df = unname(fit$df.residual),
    p = unname(sm[coef, "Pr(>|t|)"]),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    excludes_gse205335 = flags$ex205,
    drops_sclc = flags$sclc,
    luad_matched = flags$luad,
    adjusts_keratin = flags$krt,
    adjusts_leak = flags$leak,
    adjusts_malfrac = flags$frac,
    I2 = NA_real_,
    stringsAsFactors = FALSE
  ))
  invisible(fit)
}

lme_one <- function(df, terms, coef, family, model, flags) {
  df <- prepare(df)
  df$score <- df[[family]]
  df <- df[complete.cases(df[, c("score", "cohort", terms), drop = FALSE]), ]
  df <- droplevels(df)
  if (nlevels(df$cohort) < 3 || nrow(df) < 12) return(invisible(NULL))
  form <- as.formula(paste("score ~", paste(terms, collapse = " + ")))
  fit <- tryCatch(
    lme(form, random = ~ 1 | cohort, data = df, method = "REML",
        control = lmeControl(msMaxIter = 200, opt = "optim", returnObject = TRUE)),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    note("LMM failed", model, family, fit$message)
    return(invisible(NULL))
  }
  sm <- summary(fit)$tTable
  if (!coef %in% rownames(sm)) return(invisible(NULL))
  df_l <- unname(sm[coef, "DF"])
  est <- unname(sm[coef, "Value"])
  se <- unname(sm[coef, "Std.Error"])
  tc <- tail_counts(df)
  add_row(data.frame(
    model = model,
    family = family,
    estimator = "lmm",
    n = nrow(df),
    n_q1 = tc$n1,
    n_q4 = tc$n4,
    n_cohorts = nlevels(df$cohort),
    cohorts = paste(levels(df$cohort), collapse = "+"),
    n_genes = unname(MODULE_GENE_N[[family]]),
    coef = coef,
    logFC = est,
    se = se,
    df = df_l,
    p = unname(sm[coef, "p-value"]),
    ci_low = est - qt(0.975, df_l) * se,
    ci_high = est + qt(0.975, df_l) * se,
    excludes_gse205335 = flags$ex205,
    drops_sclc = flags$sclc,
    luad_matched = flags$luad,
    adjusts_keratin = flags$krt,
    adjusts_leak = flags$leak,
    adjusts_malfrac = flags$frac,
    I2 = NA_real_,
    stringsAsFactors = FALSE
  ))
}

within_meta <- function(df, terms, coef, family, model, flags) {
  df0 <- prepare(df)
  ys <- ss <- numeric()
  used <- character()
  ns <- n1s <- n4s <- integer()
  for (co in COHORTS) {
    d <- df0[df0$cohort == co, ]
    if (!nrow(d)) next
    d$score <- d[[family]]
    d <- d[complete.cases(d[, c("score", terms), drop = FALSE]), ]
    if (nrow(d) < 5) next
    if (coef == "q4" && length(unique(d$q4)) < 2) next
    form <- as.formula(paste("score ~", paste(terms, collapse = " + ")))
    fit <- tryCatch(lm(form, data = d), error = function(e) e)
    if (inherits(fit, "error")) next
    sm <- summary(fit)$coefficients
    if (!coef %in% rownames(sm)) next
    ci <- suppressMessages(confint(fit)[coef, ])
    n1_here <- if (coef == "q4") sum(d$q4 == 0) else NA_integer_
    n4_here <- if (coef == "q4") sum(d$q4 == 1) else NA_integer_
    add_row(data.frame(
      model = paste0(model, "_", co),
      family = family,
      estimator = "within",
      n = nrow(d),
      n_q1 = n1_here,
      n_q4 = n4_here,
      n_cohorts = 1L,
      cohorts = co,
      n_genes = unname(MODULE_GENE_N[[family]]),
      coef = coef,
      logFC = unname(sm[coef, "Estimate"]),
      se = unname(sm[coef, "Std. Error"]),
      df = unname(fit$df.residual),
      p = unname(sm[coef, "Pr(>|t|)"]),
      ci_low = unname(ci[1]),
      ci_high = unname(ci[2]),
      excludes_gse205335 = co != "GSE205335" && flags$ex205,
      drops_sclc = flags$sclc,
      luad_matched = flags$luad,
      adjusts_keratin = flags$krt,
      adjusts_leak = flags$leak,
      adjusts_malfrac = flags$frac,
      I2 = NA_real_,
      stringsAsFactors = FALSE
    ))
    # This row is one cohort. The flag records whether that cohort is GSE205335.
    rows[[length(rows)]]$excludes_gse205335 <<- (co != "GSE205335")
    ys <- c(ys, unname(sm[coef, "Estimate"]))
    ss <- c(ss, unname(sm[coef, "Std. Error"]))
    used <- c(used, co)
    ns <- c(ns, nrow(d))
    n1s <- c(n1s, n1_here)
    n4s <- c(n4s, n4_here)
  }
  if (length(ys) < 2) return(invisible(NULL))
  mk <- tryCatch(metafor::rma(yi = ys, sei = ss, method = "DL", test = "knha"), error = function(e) e)
  if (inherits(mk, "error")) {
    note("meta failed", model, family, mk$message)
    return(invisible(NULL))
  }
  add_row(data.frame(
    model = paste0(model, "_meta"),
    family = family,
    estimator = "re_meta",
    n = sum(ns),
    n_q1 = if (coef == "q4") sum(n1s) else NA_integer_,
    n_q4 = if (coef == "q4") sum(n4s) else NA_integer_,
    n_cohorts = mk$k,
    cohorts = paste(used, collapse = "+"),
    n_genes = unname(MODULE_GENE_N[[family]]),
    coef = coef,
    logFC = as.numeric(mk$beta),
    se = as.numeric(mk$se),
    df = mk$k - 1,
    p = mk$pval,
    ci_low = mk$ci.lb,
    ci_high = mk$ci.ub,
    excludes_gse205335 = !"GSE205335" %in% used,
    drops_sclc = flags$sclc,
    luad_matched = flags$luad,
    adjusts_keratin = flags$krt,
    adjusts_leak = flags$leak,
    adjusts_malfrac = flags$frac,
    I2 = mk$I2,
    stringsAsFactors = FALSE
  ))
}

spearman_meta <- function(df, family, model, flags) {
  df <- df[is.finite(df[[family]]) & is.finite(df$cldn4_pct), ]
  rhos <- ns <- character()
  ys <- ss <- numeric()
  used <- character()
  for (co in COHORTS) {
    d <- df[df$cohort == co, ]
    if (nrow(d) < 5) next
    ct <- suppressWarnings(cor.test(d$cldn4_pct, d[[family]], method = "spearman", exact = FALSE))
    rho <- unname(ct$estimate)
    n <- nrow(d)
    add_row(data.frame(
      model = paste0(model, "_", co),
      family = family,
      estimator = "spearman",
      n = n,
      n_q1 = NA_integer_,
      n_q4 = NA_integer_,
      n_cohorts = 1L,
      cohorts = co,
      n_genes = unname(MODULE_GENE_N[[family]]),
      coef = "spearman_rho",
      logFC = rho,
      se = NA_real_,
      df = n - 2,
      p = ct$p.value,
      ci_low = NA_real_,
      ci_high = NA_real_,
      excludes_gse205335 = co != "GSE205335",
      drops_sclc = flags$sclc,
      luad_matched = flags$luad,
      adjusts_keratin = FALSE,
      adjusts_leak = FALSE,
      adjusts_malfrac = FALSE,
      I2 = NA_real_,
      stringsAsFactors = FALSE
    ))
    ys <- c(ys, atanh(max(min(rho, 0.999), -0.999)))
    ss <- c(ss, 1 / sqrt(n - 3))
    used <- c(used, co)
  }
  if (length(ys) < 2) return(invisible(NULL))
  mk <- metafor::rma(yi = ys, sei = ss, method = "DL", test = "knha")
  add_row(data.frame(
    model = paste0(model, "_meta"),
    family = family,
    estimator = "spearman_meta",
    n = sum(df$cohort %in% used),
    n_q1 = NA_integer_,
    n_q4 = NA_integer_,
    n_cohorts = mk$k,
    cohorts = paste(used, collapse = "+"),
    n_genes = unname(MODULE_GENE_N[[family]]),
    coef = "spearman_rho",
    logFC = tanh(as.numeric(mk$beta)),
    # mk$se is the SE of atanh(rho). Leave se blank so it is not read as an SE of rho.
    # ci_low and ci_high are the Knapp-Hartung interval on the Fisher-z scale, transformed back.
    se = NA_real_,
    df = mk$k - 1,
    p = mk$pval,
    ci_low = tanh(mk$ci.lb),
    ci_high = tanh(mk$ci.ub),
    excludes_gse205335 = !"GSE205335" %in% used,
    drops_sclc = flags$sclc,
    luad_matched = flags$luad,
    adjusts_keratin = FALSE,
    adjusts_leak = FALSE,
    adjusts_malfrac = FALSE,
    I2 = mk$I2,
    stringsAsFactors = FALSE
  ))
}

flag <- function(ex205 = FALSE, sclc = FALSE, luad = FALSE, krt = FALSE, leak = FALSE, frac = FALSE) {
  list(ex205 = ex205, sclc = sclc, luad = luad, krt = krt, leak = leak, frac = frac)
}

run_block <- function(df, terms, coef, model, flags, families, do_lme = FALSE, do_within = FALSE) {
  for (fn in families) {
    lm_one(df, terms, coef, fn, model, flags)
    if (do_lme) lme_one(df, terms, coef, fn, model, flags)
    if (do_within) within_meta(df, terms, coef, fn, paste0(model, "_within"), flags)
  }
}

# Sample subsets.
q_all <- frame_q
q_drop205 <- frame_q[frame_q$cohort != "GSE205335", ]
q_drop_sclc <- frame_q[frame_q$cancer_subtype != "SCLC", ]
q_luad <- frame_q[frame_q$cohort != "GSE205335" | frame_q$cancer_subtype == "ADC", ]
all_u <- frame_all
all_drop205 <- frame_all[frame_all$cohort != "GSE205335", ]
all_luad <- frame_all[frame_all$cohort != "GSE205335" | frame_all$cancer_subtype == "ADC", ]
adc_205 <- frame_all[frame_all$cohort == "GSE205335" & frame_all$cancer_subtype == "ADC", ]

note(sprintf(
  "n q4q1=%d drop205=%d dropSCLC=%d luad_q=%d cont=%d cont_drop205=%d cont_luad=%d adc205=%d",
  nrow(q_all), nrow(q_drop205), nrow(q_drop_sclc), nrow(q_luad),
  nrow(all_u), nrow(all_drop205), nrow(all_luad), nrow(adc_205)
))
note(sprintf(
  "q4q1 Q1/Q4 %d/%d; drop205 %d/%d; dropSCLC %d/%d; luad %d/%d",
  sum(q_all$q4 == 0), sum(q_all$q4 == 1),
  sum(q_drop205$q4 == 0), sum(q_drop205$q4 == 1),
  sum(q_drop_sclc$q4 == 0), sum(q_drop_sclc$q4 == 1),
  sum(q_luad$q4 == 0), sum(q_luad$q4 == 1)
))
stopifnot(nrow(q_all) == 34, sum(q_all$q4 == 0) == 18, sum(q_all$q4 == 1) == 16)
stopifnot(nrow(q_drop205) == 23, sum(q_drop205$q4 == 0) == 13, sum(q_drop205$q4 == 1) == 10)
stopifnot(nrow(q_drop_sclc) == 31, sum(q_drop_sclc$q4 == 0) == 18, sum(q_drop_sclc$q4 == 1) == 13)
stopifnot(nrow(q_luad) == 29, sum(q_luad$q4 == 0) == 17, sum(q_luad$q4 == 1) == 12)
stopifnot(nrow(all_u) == 64, nrow(all_drop205) == 43, nrow(all_luad) == 56, nrow(adc_205) == 13)
stopifnot(!any(q_drop_sclc$cancer_subtype == "SCLC"))
stopifnot(all(q_luad$cancer_subtype %in% c("LUAD_spectrum", "ADC")))

FAM_FULL <- c(CLAIM, "STING_core", "STING_plus_IFN_genes", "NHEJ", "TJ", "keratin")

run_block(q_all, "q4", "q4", "q4q1", flag(), FAM_FULL, do_lme = TRUE, do_within = TRUE)
run_block(q_drop205, "q4", "q4", "q4q1_drop205335", flag(ex205 = TRUE), FAM_FULL, do_lme = TRUE, do_within = TRUE)
run_block(q_drop_sclc, "q4", "q4", "q4q1_dropSCLC", flag(sclc = TRUE), FAM_FULL, do_lme = TRUE)
run_block(q_luad, "q4", "q4", "q4q1_luad", flag(luad = TRUE), FAM_FULL, do_lme = TRUE)

run_block(q_all, c("q4", "krt_z"), "q4", "q4q1_krt", flag(krt = TRUE), c(CLAIM, "STING_core", "NHEJ", "TJ"), do_lme = TRUE)
run_block(q_drop205, c("q4", "krt_z"), "q4", "q4q1_krt_drop205335", flag(ex205 = TRUE, krt = TRUE), c(CLAIM, "STING_core", "NHEJ"), do_lme = TRUE)
run_block(q_all, c("q4", "leak_z"), "q4", "q4q1_leak", flag(leak = TRUE), c(CLAIM, "STING_core", "NHEJ"))
run_block(q_drop205, c("q4", "leak_z"), "q4", "q4q1_leak_drop205335", flag(ex205 = TRUE, leak = TRUE), CLAIM)
run_block(q_all, c("q4", "krt_z", "leak_z"), "q4", "q4q1_krt_leak", flag(krt = TRUE, leak = TRUE), c(CLAIM, "STING_core", "NHEJ"), do_lme = TRUE)
run_block(q_drop205, c("q4", "krt_z", "leak_z"), "q4", "q4q1_krt_leak_drop205335", flag(ex205 = TRUE, krt = TRUE, leak = TRUE), c(CLAIM, "STING_core", "NHEJ"), do_lme = TRUE)
run_block(q_all, c("q4", "frac_z"), "q4", "q4q1_malfrac", flag(frac = TRUE), CLAIM)
run_block(q_drop205, c("q4", "frac_z"), "q4", "q4q1_malfrac_drop205335", flag(ex205 = TRUE, frac = TRUE), CLAIM)
run_block(q_luad, c("q4", "krt_z", "leak_z"), "q4", "q4q1_luad_krt_leak", flag(luad = TRUE, krt = TRUE, leak = TRUE), CLAIM)

run_block(all_u, "cldn4_z", "cldn4_z", "cont", flag(), FAM_FULL, do_lme = TRUE, do_within = TRUE)
run_block(all_drop205, "cldn4_z", "cldn4_z", "cont_drop205335", flag(ex205 = TRUE), FAM_FULL, do_lme = TRUE, do_within = TRUE)
run_block(all_u, c("cldn4_z", "krt_z"), "cldn4_z", "cont_krt", flag(krt = TRUE), c(CLAIM, "STING_core", "NHEJ", "TJ"))
run_block(all_drop205, c("cldn4_z", "krt_z"), "cldn4_z", "cont_krt_drop205335", flag(ex205 = TRUE, krt = TRUE), c(CLAIM, "STING_core", "NHEJ"), do_lme = TRUE, do_within = TRUE)
run_block(all_u, c("cldn4_z", "leak_z"), "cldn4_z", "cont_leak", flag(leak = TRUE), CLAIM)
run_block(all_drop205, c("cldn4_z", "leak_z"), "cldn4_z", "cont_leak_drop205335", flag(ex205 = TRUE, leak = TRUE), CLAIM)
run_block(all_u, c("cldn4_z", "krt_z", "leak_z"), "cldn4_z", "cont_krt_leak", flag(krt = TRUE, leak = TRUE), c(CLAIM, "STING_core", "NHEJ"))
run_block(all_drop205, c("cldn4_z", "krt_z", "leak_z"), "cldn4_z", "cont_krt_leak_drop205335", flag(ex205 = TRUE, krt = TRUE, leak = TRUE), c(CLAIM, "STING_core", "NHEJ"), do_lme = TRUE)
run_block(all_luad, "cldn4_z", "cldn4_z", "cont_luad", flag(luad = TRUE), FAM_FULL, do_lme = TRUE)
run_block(all_luad, c("cldn4_z", "krt_z", "leak_z"), "cldn4_z", "cont_luad_krt_leak", flag(luad = TRUE, krt = TRUE, leak = TRUE), CLAIM)
run_block(adc_205, "cldn4_z", "cldn4_z", "cont_GSE205335_ADC", flag(luad = TRUE), FAM_FULL)
run_block(all_u, c("cldn4_z", "frac_z"), "cldn4_z", "cont_malfrac", flag(frac = TRUE), CLAIM)
run_block(all_drop205, c("cldn4_z", "frac_z"), "cldn4_z", "cont_malfrac_drop205335", flag(ex205 = TRUE, frac = TRUE), CLAIM)

for (fn in c(CLAIM, "STING_core", "NHEJ")) {
  spearman_meta(all_u, fn, "cont_spearman", flag())
  spearman_meta(all_drop205, fn, "cont_spearman_drop205335", flag(ex205 = TRUE))
  spearman_meta(all_luad, fn, "cont_spearman_luad", flag(luad = TRUE))
}

# --- voom / edgeR gene-set checks on the same contrasts ----------------------

harmonize_cgas <- function(mat) {
  if ("CGAS" %in% rownames(mat)) {
    v <- mat["CGAS", , drop = FALSE]
  } else if ("MB21D1" %in% rownames(mat)) {
    v <- mat["MB21D1", , drop = FALSE]
  } else {
    return(mat)
  }
  rownames(v) <- "cGAS"
  rbind(mat, v)
}

bind_mat <- function(df) {
  pieces <- lapply(seq_len(nrow(df)), function(i) {
    # as.character: a factor index uses the level code, which is the wrong
    # cohort once a middle cohort has been dropped from the levels.
    co <- as.character(df$cohort[i])
    id <- as.character(df$patient[i])
    harmonize_cgas(counts[[co]][, id, drop = FALSE])
  })
  genes <- Reduce(intersect, lapply(pieces, rownames))
  mat <- do.call(cbind, lapply(pieces, function(m) m[genes, , drop = FALSE]))
  colnames(mat) <- paste(df$cohort, df$patient, sep = "|")
  stopifnot(ncol(mat) == nrow(df))
  list(mat = mat, meta = df)
}

fry_one <- function(mat, design, coef, sets, model, method) {
  y <- DGEList(mat)
  y <- calcNormFactors(y)
  keep <- filterByExpr(y, design = design, min.count = 10)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  if (method == "voom") {
    v <- voom(y, design, plot = FALSE)
    fit <- eBayes(lmFit(v, design), robust = TRUE)
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    idx <- lapply(sets, function(g) intersect(g, rownames(v)))
    idx <- idx[lengths(idx) >= 3]
    fr <- fry(v, idx, design, contrast = coef)
    list(fit_tab = tt, fry = fr, genes = rownames(v))
  } else {
    y <- estimateDisp(y, design, robust = TRUE)
    qlf <- glmQLFit(y, design, robust = TRUE)
    et <- glmQLFTest(qlf, coef = coef)$table
    # fry needs a voom-like object; report gene signs from edgeR only.
    list(fit_tab = et, fry = NULL, genes = rownames(et))
  }
}

design_from <- function(df, terms) {
  df <- prepare(df)
  rhs <- terms
  if (nlevels(droplevels(df$cohort)) >= 2) rhs <- c("cohort", terms)
  # q4 and z columns exist after prepare.
  mf <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(mf, data = df)
  list(design = design, df = df)
}

GENE_SETS <- c(claim_sets, list(STING_core = STING_CORE, NHEJ = NHEJ_GENES))
# STING_core symbols in the count matrix do not include the cGAS alias row.
# Fry uses raw symbols: MB21D1 and CGAS both listed; a cohort contributes the one it has.
GENE_SETS[["STING_core"]] <- c("cGAS", "TMEM173", "TBK1", "IKBKE", "IRF3", "DDX41", "IFI16", "AIM2")

run_contrast <- function(job, method) {
  df <- prepare(job$df)
  use <- c("cohort", job$terms)
  df <- df[complete.cases(df[, use, drop = FALSE]), ]
  df$cohort <- factor(df$cohort, levels = intersect(COHORTS, unique(as.character(df$cohort))))
  rhs <- job$terms
  if (nlevels(df$cohort) >= 2) rhs <- c("cohort", job$terms)
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = df)
  built <- bind_mat(df)
  if (nrow(design) != ncol(built$mat)) stop("design/matrix mismatch in ", job$id)
  if (!job$coef %in% colnames(design)) stop("missing coef ", job$coef, " in ", job$id, ": ", paste(colnames(design), collapse = ","))
  note(method, job$id, "n", nrow(df), "coef", job$coef, "genes", nrow(built$mat))
  out <- fry_one(built$mat, design, job$coef, GENE_SETS, job$id, method)
  tc <- tail_counts(df)
  out$n <- nrow(df)
  out$n_q1 <- tc$n1
  out$n_q4 <- tc$n4
  out
}

voom_jobs <- list(
  list(id = "q4q1_drop205335", df = q_drop205, terms = "q4", coef = "q4"),
  list(id = "q4q1_dropSCLC", df = q_drop_sclc, terms = "q4", coef = "q4"),
  list(id = "q4q1_luad", df = q_luad, terms = "q4", coef = "q4"),
  list(id = "q4q1_krt", df = q_all, terms = c("q4", "krt_z"), coef = "q4"),
  list(id = "q4q1_krt_leak", df = q_all, terms = c("q4", "krt_z", "leak_z"), coef = "q4"),
  list(id = "q4q1_krt_leak_drop205335", df = q_drop205, terms = c("q4", "krt_z", "leak_z"), coef = "q4"),
  list(id = "cont", df = all_u, terms = "cldn4_z", coef = "cldn4_z"),
  list(id = "cont_drop205335", df = all_drop205, terms = "cldn4_z", coef = "cldn4_z"),
  list(id = "cont_krt_drop205335", df = all_drop205, terms = c("cldn4_z", "krt_z"), coef = "cldn4_z"),
  list(id = "cont_krt_leak_drop205335", df = all_drop205, terms = c("cldn4_z", "krt_z", "leak_z"), coef = "cldn4_z"),
  list(id = "cont_luad", df = all_luad, terms = "cldn4_z", coef = "cldn4_z")
)

gst_rows <- list()
gene_sign_rows <- list()
for (job in voom_jobs) {
  vo <- run_contrast(job, "voom")
  fr <- vo$fry
  for (fn in rownames(fr)) {
    genes <- intersect(GENE_SETS[[fn]], vo$genes)
    lfc <- vo$fit_tab$logFC[match(genes, rownames(vo$fit_tab))]
    gst_rows[[length(gst_rows) + 1]] <- data.frame(
      model = job$id, method = "limma-voom", family = fn,
      n = vo$n,
      n_q1 = vo$n_q1,
      n_q4 = vo$n_q4,
      n_genes = length(genes),
      median_logFC = median(lfc, na.rm = TRUE),
      n_down = sum(lfc < 0, na.rm = TRUE),
      n_up = sum(lfc > 0, na.rm = TRUE),
      fry_direction = fr[fn, "Direction"],
      fry_p = fr[fn, "PValue"],
      stringsAsFactors = FALSE
    )
  }
}

# edgeR signs for the three contrasts that answer "not only GSE205335 / not SCLC".
edger_jobs <- list(
  list(id = "q4q1_dropSCLC", df = q_drop_sclc, terms = "q4", coef = "q4"),
  list(id = "cont_drop205335", df = all_drop205, terms = "cldn4_z", coef = "cldn4_z"),
  list(id = "q4q1_krt_leak_drop205335", df = q_drop205, terms = c("q4", "krt_z", "leak_z"), coef = "q4")
)
for (job in edger_jobs) {
  er <- run_contrast(job, "edger")
  for (fn in names(GENE_SETS)) {
    genes <- intersect(GENE_SETS[[fn]], rownames(er$fit_tab))
    if (length(genes) < 3) next
    lfc <- er$fit_tab[genes, "logFC"]
    p <- er$fit_tab[genes, "PValue"]
    gene_sign_rows[[length(gene_sign_rows) + 1]] <- data.frame(
      model = job$id, method = "edgeR-QL", family = fn,
      n = er$n, n_genes = length(genes),
      median_logFC = median(lfc, na.rm = TRUE),
      n_down = sum(lfc < 0, na.rm = TRUE),
      n_up = sum(lfc > 0, na.rm = TRUE),
      n_p05_expected = {
        # NHEJ, TJ, and keratin have no pre-specified sign in this grid.
        exp <- if (fn %in% names(EXPECT)) unname(EXPECT[[fn]]) else NA_real_
        if (is.na(exp)) NA_integer_
        else if (exp < 0) sum(p < 0.05 & lfc < 0, na.rm = TRUE)
        else sum(p < 0.05 & lfc > 0, na.rm = TRUE)
      },
      stringsAsFactors = FALSE
    )
  }
}

res <- do.call(rbind, rows)
rownames(res) <- NULL
if (!"I2" %in% names(res)) res$I2 <- NA_real_

# Knapp-Hartung can shrink the SE below the usual random-effects SE when Q is
# tiny. test="adhoc" floors that scale factor at 1 (Jackson et al. 2017).
# p stays the pre-specified Knapp-Hartung p-value. p_adhoc is the floored one.
res$p_adhoc <- NA_real_
res$se_adhoc <- NA_real_
meta_idx <- which(res$estimator %in% c("re_meta", "spearman_meta"))
for (i in meta_idx) {
  prefix <- sub("_meta$", "", res$model[i])
  child_est <- if (res$estimator[i] == "spearman_meta") "spearman" else "within"
  ch <- res[res$estimator == child_est & res$family == res$family[i] & startsWith(res$model, paste0(prefix, "_")), ]
  if (nrow(ch) < 2 || any(!is.finite(ch$se))) next
  ad <- tryCatch(
    metafor::rma(yi = ch$logFC, sei = ch$se, method = "DL", test = "adhoc"),
    error = function(e) e
  )
  if (inherits(ad, "error")) next
  # Spearman rows store rho in logFC and a blank SE, so this branch is the
  # within-cohort logFC meta only. Spearman meta keeps p_adhoc blank.
  if (res$estimator[i] == "spearman_meta") next
  res$p_adhoc[i] <- ad$pval
  res$se_adhoc[i] <- as.numeric(ad$se)
}
res$expect <- EXPECT[res$family]
res$sign_matches <- ifelse(is.na(res$expect), NA, sign(res$logFC) == res$expect)
res$fdr_claim <- NA_real_
claim_hit <- res$family %in% CLAIM & res$estimator %in% c("ols", "lmm", "re_meta")
res$fdr_claim[claim_hit] <- ave(res$p[claim_hit], res$model[claim_hit], res$estimator[claim_hit], FUN = bh)

gst <- do.call(rbind, gst_rows)
rownames(gst) <- NULL
edger_sign <- do.call(rbind, gene_sign_rows)
rownames(edger_sign) <- NULL

write.table(res, file.path(TAB, "robust_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gst, file.path(TAB, "robust_fry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(edger_sign, file.path(TAB, "robust_edger_signs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mod_genes <- rbind(
  data.frame(module = "NHEJ", gene = NHEJ_GENES, in_IFN_hallmark = FALSE),
  data.frame(module = "STING_core", gene = STING_CORE, in_IFN_hallmark = FALSE),
  data.frame(module = "STING_plus_IFN_genes_only", gene = c("IRF7", "ZBP1"), in_IFN_hallmark = TRUE),
  data.frame(module = "keratin_covariate", gene = KRT_COVARIATE, in_IFN_hallmark = FALSE),
  data.frame(module = "immune_leak", gene = LEAK_GENES, in_IFN_hallmark = FALSE)
)
write.table(mod_genes, file.path(TAB, "robust_module_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# n audit
n_audit <- data.frame(
  subset = c("q4q1", "q4q1_drop205335", "q4q1_dropSCLC", "q4q1_luad", "cont", "cont_drop205335", "cont_luad", "cont_GSE205335_ADC"),
  n = c(nrow(q_all), nrow(q_drop205), nrow(q_drop_sclc), nrow(q_luad), nrow(all_u), nrow(all_drop205), nrow(all_luad), nrow(adc_205)),
  n_q1 = c(sum(q_all$q4 == 0), sum(q_drop205$q4 == 0), sum(q_drop_sclc$q4 == 0), sum(q_luad$q4 == 0), NA, NA, NA, NA),
  n_q4 = c(sum(q_all$q4 == 1), sum(q_drop205$q4 == 1), sum(q_drop_sclc$q4 == 1), sum(q_luad$q4 == 1), NA, NA, NA, NA),
  stringsAsFactors = FALSE
)
write.table(n_audit, file.path(TAB, "robust_n.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(notes, file.path(TAB, "robust_notes.txt"))

# Bridge check against the main stacked OLS IFN number (-0.643 on the shared panel).
bridge <- res[res$model == "q4q1" & res$estimator == "ols" & res$family == "IFN", ]
note(sprintf(
  "bridge q4q1 OLS IFN logFC=%.6f p=%.6g n=%d n_q1=%s n_q4=%s",
  bridge$logFC, bridge$p, bridge$n, bridge$n_q1, bridge$n_q4
))
# Same estimand as tables/family_effects.tsv stacked_ols IFN (-0.643005, n=34, Q1=18, Q4=16).
stopifnot(nrow(bridge) == 1)
stopifnot(bridge$n == 34, bridge$n_q1 == 18, bridge$n_q4 == 16)
stopifnot(abs(bridge$logFC - (-0.643005243383009)) < 0.02)

# --- figure: IFN coefficient across the pre-specified grid -------------------

plot_ids <- c(
  "q4q1", "q4q1_drop205335", "q4q1_dropSCLC", "q4q1_luad",
  "q4q1_krt", "q4q1_leak", "q4q1_malfrac", "q4q1_krt_leak",
  "q4q1_krt_leak_drop205335", "q4q1_krt_drop205335",
  "cont", "cont_drop205335", "cont_krt_drop205335", "cont_krt_leak_drop205335",
  "cont_luad", "cont_malfrac_drop205335",
  "cont_drop205335_within_meta", "cont_krt_drop205335_within_meta"
)
labels <- c(
  q4q1 = "Q4 vs Q1, 4 cohorts",
  q4q1_drop205335 = "Q4 vs Q1, drop GSE205335",
  q4q1_dropSCLC = "Q4 vs Q1, drop SCLC",
  q4q1_luad = "Q4 vs Q1, LUAD / ADC only",
  q4q1_krt = "Q4 vs Q1 + KRT8/18/19",
  q4q1_leak = "Q4 vs Q1 + immune leak",
  q4q1_malfrac = "Q4 vs Q1 + malignant fraction",
  q4q1_krt_leak = "Q4 vs Q1 + KRT + immune leak",
  q4q1_krt_leak_drop205335 = "Q4 vs Q1 + KRT + leak, drop GSE205335",
  q4q1_krt_drop205335 = "Q4 vs Q1 + KRT, drop GSE205335",
  cont = "Continuous %pos, 4 cohorts",
  cont_drop205335 = "Continuous %pos, drop GSE205335",
  cont_krt_drop205335 = "Continuous + KRT, drop GSE205335",
  cont_krt_leak_drop205335 = "Continuous + KRT + leak, drop GSE205335",
  cont_luad = "Continuous, LUAD / ADC only",
  cont_malfrac_drop205335 = "Continuous + malignant fraction, drop GSE205335",
  cont_drop205335_within_meta = "RE meta of slopes, drop GSE205335",
  cont_krt_drop205335_within_meta = "RE meta of KRT-adjusted slopes, drop GSE205335"
)

draw_ifn <- function() {
  d <- res[res$family == "IFN" & res$estimator == "ols" & res$model %in% plot_ids, ]
  # within meta rows use estimator re_meta and model suffix _meta. The plot ids
  # above use _within_meta; the stored model is cont_drop205335_within_meta.
  extra <- res[res$family == "IFN" & res$estimator == "re_meta" & res$model %in% paste0(c("cont_drop205335_within", "cont_krt_drop205335_within"), "_meta"), ]
  if (nrow(extra)) {
    extra$model <- sub("_within_meta$", "_within_meta", extra$model)
    # stored name is cont_drop205335_within_meta. Map to plot id.
    extra$model <- sub("_within_meta$", "_within_meta", extra$model)
  }
  # Actual stored ids:
  meta_map <- c(
    cont_drop205335_within_meta = "cont_drop205335_within_meta",
    cont_krt_drop205335_within_meta = "cont_krt_drop205335_within_meta"
  )
  d <- res[res$family == "IFN" & (
    (res$estimator == "ols" & res$model %in% plot_ids) |
    (res$estimator == "re_meta" & res$model %in% c("cont_drop205335_within_meta", "cont_krt_drop205335_within_meta"))
  ), ]
  d$label <- labels[d$model]
  d$label[is.na(d$label)] <- d$model[is.na(d$label)]
  d <- d[match(plot_ids, d$model), ]
  d <- d[!is.na(d$model), ]
  d <- d[nrow(d):1, ]
  par(mar = c(5, 22, 4, 2))
  xlim <- range(c(d$ci_low, d$ci_high, -1.5, 0.8), finite = TRUE)
  plot(NA, xlim = xlim, ylim = c(0.5, nrow(d) + 0.5), yaxt = "n",
       xlab = "IFN score coefficient (negative = lower in CLDN4-high)",
       ylab = "", main = "IFN score without relying on one cohort")
  abline(v = 0, lty = 2, col = "grey40")
  cols <- ifelse(d$excludes_gse205335, "#08519c", "grey35")
  points(d$logFC, seq_len(nrow(d)), pch = ifelse(d$estimator == "re_meta", 18, 16), col = cols)
  segments(d$ci_low, seq_len(nrow(d)), d$ci_high, seq_len(nrow(d)), col = cols)
  axis(2, at = seq_len(nrow(d)), labels = d$label, las = 1, cex.axis = 0.62)
  legend("bottomleft", c("Includes GSE205335", "GSE205335 left out"), col = c("grey35", "#08519c"), pch = 16, bty = "n", cex = 0.75)
  mtext("Continuous coefficients are per within-cohort SD of CLDN4 %pos. Q4 vs Q1 coefficients are the group difference.",
        side = 1, line = 3.6, cex = 0.65)
}
save_both(file.path(FIG, "forest_ifn_robust"), 9.4, 7.2, draw_ifn)

# Heat map of claim-family coefficients for the same models.
draw_heat <- function() {
  ids <- c(
    "q4q1", "q4q1_drop205335", "q4q1_dropSCLC", "q4q1_luad",
    "q4q1_krt_leak", "q4q1_krt_leak_drop205335",
    "cont", "cont_drop205335", "cont_krt_leak_drop205335", "cont_luad"
  )
  fams <- c("IFN", "MHC-I/APM", "chemokine", "STING_core", "NHEJ")
  mat <- matrix(NA_real_, nrow = length(fams), ncol = length(ids), dimnames = list(fams, ids))
  ann <- matrix("", nrow = length(fams), ncol = length(ids))
  for (i in seq_along(fams)) {
    for (j in seq_along(ids)) {
      hit <- res$family == fams[i] & res$model == ids[j] & res$estimator == "ols"
      if (!any(hit)) next
      r <- res[hit, ][1, ]
      mat[i, j] <- r$logFC
      ann[i, j] <- sprintf("%.2f\np=%.3g", r$logFC, r$p)
    }
  }
  par(mar = c(12, 8, 3, 2))
  lim <- max(abs(mat), na.rm = TRUE)
  cols_r <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
  plot(NA, xlim = c(0.5, ncol(mat) + 0.5), ylim = c(0.5, nrow(mat) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       main = "Family-score coefficients across the robustness grid")
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      v <- mat[i, j]
      ci <- if (!is.finite(v)) 51 else max(1, min(101, 1 + round(100 * (v + lim) / (2 * lim))))
      rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = cols_r[ci], border = "white")
      text(j, i, ann[i, j], cex = 0.55)
    }
  }
  axis(1, at = seq_len(ncol(mat)), labels = ids, las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(nrow(mat)), labels = rownames(mat), las = 1)
  mtext("Blue = lower in CLDN4-high. Continuous models are per SD.", side = 3, line = 0.2, cex = 0.7)
}
save_both(file.path(FIG, "heatmap_robust_grid"), 10, 5.4, draw_heat)

# Print the rows that decide the writeup.
show <- res[res$estimator %in% c("ols", "lmm", "re_meta", "spearman_meta") & res$family %in% c(CLAIM, "STING_core", "NHEJ"),
            c("model", "family", "estimator", "n", "n_q1", "n_q4", "n_cohorts", "logFC", "p", "excludes_gse205335")]
show <- show[order(show$family, show$model, show$estimator), ]
write.table(show, file.path(TAB, "robust_headline.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("DONE robust")

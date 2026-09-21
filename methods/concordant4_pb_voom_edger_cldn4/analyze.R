#!/usr/bin/env Rscript
# Malignant CLDN4 Q4 vs Q1, patient-level pseudo-bulk.
#
# The count matrices are already malignant UMI sums (one column per
# donor / sample / patient). That is the muscat::pbDS collapse. This
# script does not re-sum cells. It fits the two engines pbDS uses:
#   limma-voom + robust eBayes
#   edgeR quasi-likelihood (robust glmQLFit)
# Positive logFC means higher in CLDN4-high (Q4).
#
# Quartile labels are locked in data/units.tsv (PR #503 %pos split).
# GSE189357 Q4 has 2 units. Those units stay in the fit and are flagged
# thin. A pre-specified sensitivity drops any cohort with an arm < 3.
#
# Pooled layers:
#   1. Stacked limma-voom / edgeR with a cohort intercept (common logFC).
#   2. Two-stage random-effects meta of the four cohort logFCs (DL).
#   3. Family-score linear mixed model, random intercept for cohort.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(nlme)
  library(metafor)
})

options(stringsAsFactors = FALSE)
set.seed(1)

ROOT <- Sys.getenv("PB_ROOT", unset = ".")
if (!file.exists(file.path(ROOT, "data", "units.tsv"))) {
  # Allow running from the repo root or from this directory.
  if (file.exists("methods/concordant4_pb_voom_edger_cldn4/data/units.tsv")) {
    ROOT <- "methods/concordant4_pb_voom_edger_cldn4"
  }
}
setwd(ROOT)

TAB <- "tables"
FIG <- "figures"
dir.create(TAB, showWarnings = FALSE)
dir.create(FIG, showWarnings = FALSE)

COHORTS <- c("GSE123902", "GSE131907", "GSE205335", "GSE189357")
FAM_ORDER <- c("IFN", "MHC-I/APM", "chemokine", "TJ", "keratin")
EXPECT <- c(IFN = -1, "MHC-I/APM" = -1, chemokine = -1, TJ = 1, keratin = 1)
KEY <- c(
  "CLDN4", "STAT1", "IRF1", "IFIT1", "ISG15", "MX1",
  "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9",
  "CXCL9", "CXCL10", "CCL5", "CXCL8",
  "OCLN", "TJP1", "CLDN3", "CLDN7", "CDH1", "F11R",
  "KRT7", "KRT8", "KRT18", "KRT19"
)
FAM_COLS <- c(
  IFN = "#d62728", "MHC-I/APM" = "#1f77b4", chemokine = "#ff7f0e",
  TJ = "#2ca02c", keratin = "#9467bd"
)
COH_COLS <- c(
  GSE123902 = "#4c78a8", GSE131907 = "#f58518",
  GSE205335 = "#54a24b", GSE189357 = "#b279a2"
)

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

# --- locked inputs ----------------------------------------------------------

units <- read.delim(
  "data/units.tsv",
  colClasses = c(patient = "character"),
  check.names = FALSE
)
units$in_count_matrix <- tolower(as.character(units$in_count_matrix)) %in% c("true", "1")
units$cohort <- as.character(units$cohort)
units$quartile <- as.character(units$quartile)

fam_tbl <- read.delim("data/families.tsv", stringsAsFactors = FALSE)
fam_tbl$gene <- toupper(fam_tbl$gene)
families <- split(fam_tbl$gene, fam_tbl$family)
families <- families[FAM_ORDER]
stopifnot(all(FAM_ORDER %in% names(families)))
stopifnot(!"CLDN4" %in% families[["TJ"]])

read_counts <- function(cohort) {
  path <- sprintf("data/%s_malignant_counts.tsv.gz", cohort)
  df <- read.delim(gzfile(path), check.names = FALSE, quote = "", na.strings = "", row.names = 1)
  df <- as.matrix(df)
  storage.mode(df) <- "numeric"
  if (anyNA(df)) stop("NA counts in ", cohort)
  rownames(df) <- toupper(rownames(df))
  if (anyDuplicated(rownames(df))) {
    spl <- split(seq_len(nrow(df)), rownames(df))
    df <- do.call(rbind, lapply(spl, function(i) colSums(df[i, , drop = FALSE])))
  }
  df
}

de_units <- function(cohort = NULL) {
  u <- units[units$quartile %in% c("Q1", "Q4") & units$in_count_matrix, ]
  if (!is.null(cohort)) u <- u[u$cohort == cohort, ]
  u <- u[order(match(u$cohort, COHORTS), u$patient), ]
  u$group <- factor(u$quartile, levels = c("Q1", "Q4"))
  u
}

# Hard stop if the DE tails drift away from PR #503.
expect_de <- list(
  GSE123902 = list(Q1 = c("LX675", "LX682", "LX699", "LX701"), Q4 = c("LX653", "LX680", "LX684")),
  GSE131907 = list(
    Q1 = c("EBUS_13", "EBUS_15", "EBUS_49", "NS_02", "NS_06", "NS_16"),
    Q4 = c("EBUS_19", "EBUS_28", "NS_03", "NS_04", "NS_07")
  ),
  GSE205335 = list(
    Q1 = c("P1015", "P1062", "P1063", "P1090", "P1119"),
    Q4 = c("P1016", "P1025", "P1037", "P1084", "P1089", "P1115")
  ),
  GSE189357 = list(Q1 = c("TD2", "TD4", "TD7"), Q4 = c("TD6", "TD9"))
)
for (co in COHORTS) {
  u <- de_units(co)
  for (arm in c("Q1", "Q4")) {
    got <- sort(u$patient[u$quartile == arm])
    exp <- sort(expect_de[[co]][[arm]])
    if (!identical(got, exp)) stop("DE units drifted for ", co, " ", arm)
  }
}

# --- engines ----------------------------------------------------------------

fit_engines <- function(counts, group, tag) {
  group <- factor(group, levels = c("Q1", "Q4"))
  if (length(group) != ncol(counts)) stop("group length mismatch (", tag, ")")
  y0 <- DGEList(counts)
  y0 <- calcNormFactors(y0)
  keep <- filterByExpr(y0, group = group, min.count = 10)
  if (sum(keep) < 500) note(tag, "filterByExpr kept only", sum(keep), "genes")
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", "groupQ4")

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = "groupQ4", number = Inf, sort.by = "none")
  se_v <- abs(tt$logFC / tt$t)
  se_v[!is.finite(tt$t) | tt$t == 0] <- NA_real_
  voom_df <- data.frame(
    gene = rownames(tt),
    logFC = tt$logFC,
    se = se_v,
    t = tt$t,
    p = tt$P.Value,
    fdr = tt$adj.P.Val,
    AveExpr = tt$AveExpr,
    method = "limma-voom",
    contrast = tag,
    stringsAsFactors = FALSE
  )

  y <- estimateDisp(y, design, robust = TRUE)
  qlf <- glmQLFit(y, design, robust = TRUE)
  qt <- glmQLFTest(qlf, coef = "groupQ4")
  et <- qt$table
  se_e <- rep(NA_real_, nrow(et))
  ok_f <- is.finite(et$F) & et$F > 0
  se_e[ok_f] <- abs(et$logFC[ok_f]) / sqrt(et$F[ok_f])
  edger_df <- data.frame(
    gene = rownames(et),
    logFC = et$logFC,
    se = se_e,
    t = et$logFC / se_e,
    p = et$PValue,
    fdr = bh(et$PValue),
    AveExpr = et$logCPM,
    method = "edgeR-QL",
    contrast = tag,
    stringsAsFactors = FALSE
  )

  # log2 CPM with the filtered-matrix TMM factors, prior.count = 1.
  lc <- cpm(y, log = TRUE, prior.count = 1)

  list(voom = voom_df, edger = edger_df, v = v, design = design, y = y, logcpm = lc, keep = keep)
}

fit_stacked <- function(counts, meta, tag) {
  group <- factor(meta$group, levels = c("Q1", "Q4"))
  cohort <- factor(meta$cohort, levels = COHORTS)
  cohort <- droplevels(cohort)
  y0 <- DGEList(counts)
  y0 <- calcNormFactors(y0)
  keep <- filterByExpr(y0, group = group, min.count = 10)
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  if (nlevels(cohort) >= 2) {
    design <- model.matrix(~ cohort + group)
  } else {
    design <- model.matrix(~ group)
  }
  # The Q4 coefficient is always the last column.
  coef <- colnames(design)[ncol(design)]
  if (!grepl("groupQ4", coef)) stop("unexpected coef ", coef, " in ", tag)

  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  se_v <- abs(tt$logFC / tt$t)
  se_v[!is.finite(tt$t) | tt$t == 0] <- NA_real_
  voom_df <- data.frame(
    gene = rownames(tt), logFC = tt$logFC, se = se_v, t = tt$t,
    p = tt$P.Value, fdr = tt$adj.P.Val, AveExpr = tt$AveExpr,
    method = "limma-voom", contrast = tag, stringsAsFactors = FALSE
  )

  y <- estimateDisp(y, design, robust = TRUE)
  qlf <- glmQLFit(y, design, robust = TRUE)
  et <- glmQLFTest(qlf, coef = coef)$table
  se_e <- rep(NA_real_, nrow(et))
  ok_f <- is.finite(et$F) & et$F > 0
  se_e[ok_f] <- abs(et$logFC[ok_f]) / sqrt(et$F[ok_f])
  edger_df <- data.frame(
    gene = rownames(et), logFC = et$logFC, se = se_e, t = et$logFC / se_e,
    p = et$PValue, fdr = bh(et$PValue), AveExpr = et$logCPM,
    method = "edgeR-QL", contrast = tag, stringsAsFactors = FALSE
  )
  list(
    voom = voom_df, edger = edger_df, v = v, design = design, y = y,
    coef = coef, logcpm = cpm(y, log = TRUE, prior.count = 1)
  )
}

# DL random-effects meta. Matches metafor::rma(method="DL", test="z")
# on the point estimate, tau2, and z p-value. p_t uses t with k-1 df.
dl_one <- function(yi, sei) {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]
  sei <- sei[ok]
  k <- length(yi)
  empty <- list(k = k, b = NA_real_, se = NA_real_, p_z = NA_real_, p_t = NA_real_,
                tau2 = NA_real_, I2 = NA_real_, Q = NA_real_)
  if (k < 2) return(empty)
  vi <- sei^2
  wi <- 1 / vi
  bfe <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - bfe)^2)
  C <- sum(wi) - sum(wi^2) / sum(wi)
  df <- k - 1
  tau2 <- if (is.finite(C) && C > 0) max(0, (Q - df) / C) else 0
  wre <- 1 / (vi + tau2)
  b <- sum(wre * yi) / sum(wre)
  se <- sqrt(1 / sum(wre))
  z <- b / se
  I2 <- if (Q > 0) max(0, (Q - df) / Q) else 0
  list(
    k = k, b = b, se = se,
    p_z = 2 * pnorm(-abs(z)),
    p_t = 2 * pt(-abs(z), df = df),
    tau2 = tau2, I2 = I2, Q = Q, w = wre / sum(wre)
  )
}

dl_mat <- function(yi, sei) {
  vi <- sei^2
  bad <- !is.finite(yi) | !is.finite(sei) | sei <= 0
  yi[bad] <- NA
  vi[bad] <- NA
  wi <- 1 / vi
  k <- rowSums(is.finite(yi))
  sw <- rowSums(wi, na.rm = TRUE)
  bfe <- rowSums(wi * yi, na.rm = TRUE) / sw
  Q <- rowSums(wi * (yi - bfe)^2, na.rm = TRUE)
  sw2 <- rowSums(wi^2, na.rm = TRUE)
  C <- sw - sw2 / sw
  df <- k - 1
  tau2 <- pmax(0, (Q - df) / C)
  tau2[k < 2 | !is.finite(C) | C <= 0] <- 0
  tau2[k < 2] <- NA
  wre <- 1 / (vi + tau2)
  wre[!is.finite(yi)] <- NA
  swre <- rowSums(wre, na.rm = TRUE)
  b <- rowSums(replace(wre * yi, !is.finite(yi), 0), na.rm = TRUE) / swre
  se <- sqrt(1 / swre)
  z <- b / se
  I2 <- ifelse(Q > 0, pmax(0, (Q - df) / Q), 0)
  p_z <- 2 * pnorm(-abs(z))
  p_t <- 2 * pt(-abs(z), df = df)
  b[k < 2] <- NA
  se[k < 2] <- NA
  p_z[k < 2] <- NA
  p_t[k < 2] <- NA
  I2[k < 2] <- NA
  # Dominant cohort, for an honesty check on whether one study carries the gene.
  share <- wre / swre
  share[!is.finite(share)] <- NA
  max_w <- apply(share, 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else max(x)
  })
  share[!is.finite(share)] <- -Inf
  dom <- max.col(share, ties.method = "first")
  cn <- colnames(yi)
  if (is.null(cn)) cn <- paste0("c", seq_len(ncol(yi)))
  data.frame(
    k = as.integer(k), logFC = b, se = se, p_z = p_z, p_t = p_t,
    tau2 = tau2, I2 = I2, Q = Q,
    max_weight_share = as.numeric(max_w),
    dominant_cohort = cn[dom],
    stringsAsFactors = FALSE
  )
}

# Toy check against metafor before any biology is touched.
toy_y <- c(0.2, -0.4, 0.15, 0.05)
toy_s <- c(0.1, 0.3, 0.12, 0.25)
toy <- dl_one(toy_y, toy_s)
toy_m <- metafor::rma(yi = toy_y, sei = toy_s, method = "DL")
if (abs(toy$b - as.numeric(toy_m$beta)) > 1e-8) stop("DL point estimate != metafor")
if (abs(toy$se - as.numeric(toy_m$se)) > 1e-8) stop("DL SE != metafor")
if (abs(toy$tau2 - toy_m$tau2) > 1e-8) stop("DL tau2 != metafor")
toy_mat <- dl_mat(rbind(toy_y), rbind(toy_s))
if (abs(toy_mat$logFC[1] - toy$b) > 1e-8) stop("vectorized DL != dl_one")
note("DL meta matches metafor on the toy example")

# --- load counts and fit each cohort ---------------------------------------

counts <- lapply(COHORTS, read_counts)
names(counts) <- COHORTS

cohort_fit <- list()
gene_rows <- list()
for (co in COHORTS) {
  u <- de_units(co)
  note(sprintf(
    "%s DE n_Q1=%d n_Q4=%d (unit=%s). Thin arm=%s",
    co, sum(u$quartile == "Q1"), sum(u$quartile == "Q4"), u$unit[1],
    ifelse(min(sum(u$quartile == "Q1"), sum(u$quartile == "Q4")) < 3, "YES", "no")
  ))
  mat <- counts[[co]][, u$patient, drop = FALSE]
  fit <- fit_engines(mat, u$group, co)
  cohort_fit[[co]] <- fit
  gene_rows[[paste(co, "voom")]] <- fit$voom
  gene_rows[[paste(co, "edger")]] <- fit$edger
  note(sprintf(
    "%s tested genes voom=%d edger=%d",
    co, nrow(fit$voom), nrow(fit$edger)
  ))
}

# Family-score gene panel: count >= 10 in >= 2 DE samples of every cohort.
present_ids <- lapply(COHORTS, function(co) {
  u <- de_units(co)
  mat <- counts[[co]][, u$patient, drop = FALSE]
  rownames(mat)[rowSums(mat >= 10) >= 2]
})
common_genes <- Reduce(intersect, present_ids)
note("genes with count>=10 in >=2 DE samples of every cohort:", length(common_genes))

score_genes <- lapply(families, function(g) intersect(g, common_genes))
note("family score panel:", paste(sprintf("%s=%d", names(score_genes), lengths(score_genes)), collapse = ", "))

# Per-unit family scores. One shared gene panel per family, TMM logCPM
# from the full Q1/Q4 matrix (not the DE filter), so cohorts are scored
# on the same genes.
score_long <- list()
for (co in COHORTS) {
  u <- de_units(co)
  mat <- counts[[co]][, u$patient, drop = FALSE]
  y_sc <- calcNormFactors(DGEList(mat))
  lc <- cpm(y_sc, log = TRUE, prior.count = 1)
  for (fn in FAM_ORDER) {
    genes <- intersect(score_genes[[fn]], rownames(lc))
    if (length(genes) < 3) {
      note("skip score", fn, co, "n_genes", length(genes))
      next
    }
    if (length(genes) != length(score_genes[[fn]])) {
      note("score panel shrunk", fn, co, length(genes), "of", length(score_genes[[fn]]))
    }
    sc <- colMeans(lc[genes, u$patient, drop = FALSE])
    score_long[[paste(co, fn)]] <- data.frame(
      patient = u$patient,
      cohort = co,
      quartile = u$quartile,
      q4 = as.integer(u$quartile == "Q4"),
      n_malignant = u$n_malignant,
      unit = u$unit,
      family = fn,
      n_genes = length(genes),
      score = as.numeric(sc),
      stringsAsFactors = FALSE
    )
  }
}
scores <- do.call(rbind, score_long)
rownames(scores) <- NULL

# Within-cohort family-score lm, then metafor KNHA meta.
fam_effects <- list()
add_effect <- function(row) fam_effects[[length(fam_effects) + 1]] <<- row

for (fn in FAM_ORDER) {
  ys <- numeric()
  ss <- numeric()
  used <- character()
  for (co in COHORTS) {
    d <- scores[scores$family == fn & scores$cohort == co, ]
    n1 <- sum(d$quartile == "Q1")
    n4 <- sum(d$quartile == "Q4")
    fit <- lm(score ~ q4, data = d)
    sm <- summary(fit)$coefficients
    est <- unname(sm["q4", "Estimate"])
    se <- unname(sm["q4", "Std. Error"])
    df <- unname(fit$df.residual)
    p <- unname(sm["q4", "Pr(>|t|)"])
    ci <- suppressMessages(confint(fit)["q4", ])
    add_effect(data.frame(
      family = fn, estimator = paste0("cohort_", co),
      n_q1 = n1, n_q4 = n4, n = n1 + n4, k = 1,
      n_genes = d$n_genes[1], logFC = est, se = se, df = df, p = p,
      ci_low = ci[1], ci_high = ci[2], I2 = NA_real_, tau2 = NA_real_,
      thin = n1 < 3 || n4 < 3,
      stringsAsFactors = FALSE
    ))
    ys <- c(ys, est)
    ss <- c(ss, se)
    used <- c(used, co)
  }
  # All cohorts, including the thin GSE189357 arm.
  mk <- metafor::rma(yi = ys, sei = ss, method = "DL", test = "knha")
  add_effect(data.frame(
    family = fn, estimator = "re_meta",
    n_q1 = sum(de_units()$quartile == "Q1"),
    n_q4 = sum(de_units()$quartile == "Q4"),
    n = nrow(de_units()), k = mk$k,
    n_genes = unique(scores$n_genes[scores$family == fn])[1],
    logFC = as.numeric(mk$beta), se = as.numeric(mk$se), df = mk$k - 1,
    p = mk$pval, ci_low = mk$ci.lb, ci_high = mk$ci.ub,
    I2 = mk$I2, tau2 = mk$tau2, thin = TRUE,
    stringsAsFactors = FALSE
  ))
  # Pre-specified: drop cohorts with an arm smaller than 3.
  keep_co <- used[sapply(used, function(co) {
    d <- scores[scores$family == fn & scores$cohort == co, ]
    min(sum(d$quartile == "Q1"), sum(d$quartile == "Q4")) >= 3
  })]
  mk3 <- metafor::rma(
    yi = ys[used %in% keep_co], sei = ss[used %in% keep_co],
    method = "DL", test = "knha"
  )
  n_u <- de_units()
  n_u <- n_u[n_u$cohort %in% keep_co, ]
  add_effect(data.frame(
    family = fn, estimator = "re_meta_drop_thin",
    n_q1 = sum(n_u$quartile == "Q1"), n_q4 = sum(n_u$quartile == "Q4"),
    n = nrow(n_u), k = mk3$k,
    n_genes = unique(scores$n_genes[scores$family == fn])[1],
    logFC = as.numeric(mk3$beta), se = as.numeric(mk3$se), df = mk3$k - 1,
    p = mk3$pval, ci_low = mk3$ci.lb, ci_high = mk3$ci.ub,
    I2 = mk3$I2, tau2 = mk3$tau2, thin = FALSE,
    stringsAsFactors = FALSE
  ))

  d <- scores[scores$family == fn, ]
  ols <- lm(score ~ cohort + q4, data = d)
  sm <- summary(ols)$coefficients
  ci <- suppressMessages(confint(ols)["q4", ])
  add_effect(data.frame(
    family = fn, estimator = "stacked_ols",
    n_q1 = sum(d$quartile == "Q1"), n_q4 = sum(d$quartile == "Q4"),
    n = nrow(d), k = length(unique(d$cohort)),
    n_genes = d$n_genes[1],
    logFC = unname(sm["q4", "Estimate"]), se = unname(sm["q4", "Std. Error"]),
    df = unname(ols$df.residual), p = unname(sm["q4", "Pr(>|t|)"]),
    ci_low = ci[1], ci_high = ci[2], I2 = NA_real_, tau2 = NA_real_, thin = FALSE,
    stringsAsFactors = FALSE
  ))

  lmm_ok <- TRUE
  lmm <- tryCatch(
    lme(score ~ q4, random = ~ 1 | cohort, data = d, method = "REML",
        control = lmeControl(msMaxIter = 200, opt = "optim", returnObject = TRUE)),
    error = function(e) e
  )
  if (inherits(lmm, "error")) {
    note("LMM failed for", fn, ":", lmm$message)
    lmm_ok <- FALSE
  }
  if (lmm_ok) {
    sm <- summary(lmm)$tTable
    vc <- VarCorr(lmm)
    tau <- as.numeric(vc[1, "StdDev"])
    add_effect(data.frame(
      family = fn, estimator = "lmm_cohort",
      n_q1 = sum(d$quartile == "Q1"), n_q4 = sum(d$quartile == "Q4"),
      n = nrow(d), k = length(unique(d$cohort)),
      n_genes = d$n_genes[1],
      logFC = unname(sm["q4", "Value"]), se = unname(sm["q4", "Std.Error"]),
      df = unname(sm["q4", "DF"]), p = unname(sm["q4", "p-value"]),
      ci_low = unname(sm["q4", "Value"] - qt(0.975, sm["q4", "DF"]) * sm["q4", "Std.Error"]),
      ci_high = unname(sm["q4", "Value"] + qt(0.975, sm["q4", "DF"]) * sm["q4", "Std.Error"]),
      I2 = NA_real_, tau2 = tau^2, thin = FALSE,
      stringsAsFactors = FALSE
    ))
    note(sprintf("LMM %s logFC=%.3f SE=%.3f p=%.3g cohort_SD=%.3f",
                 fn, sm["q4", "Value"], sm["q4", "Std.Error"], sm["q4", "p-value"], tau))
  }
}
effects <- do.call(rbind, fam_effects)
rownames(effects) <- NULL
effects$fdr <- ave(effects$p, effects$estimator, FUN = bh)
effects$expect <- EXPECT[effects$family]
effects$sign_matches <- sign(effects$logFC) == effects$expect

# --- stacked voom / edgeR ---------------------------------------------------

u_all <- de_units()
# Column-bind genes by intersection so a missing symbol cannot shift columns.
gene_universe <- Reduce(intersect, lapply(counts, rownames))
mat_all <- do.call(cbind, lapply(COHORTS, function(co) {
  u <- de_units(co)
  counts[[co]][gene_universe, u$patient, drop = FALSE]
}))
mat_all <- mat_all[, u_all$patient, drop = FALSE]
stopifnot(identical(colnames(mat_all), u_all$patient))
stacked <- fit_stacked(mat_all, u_all, "stacked")
note(sprintf("stacked voom genes=%d edger genes=%d coef=%s",
             nrow(stacked$voom), nrow(stacked$edger), stacked$coef))

# Gene-set tests on the stacked voom fit (self-contained fry, competitive camera).
idx <- lapply(families, function(g) intersect(g, rownames(stacked$v)))
idx <- idx[lengths(idx) >= 3]
fry_tab <- fry(stacked$v, idx, stacked$design, contrast = stacked$coef)
fry_tab$family <- rownames(fry_tab)
fry_tab$test <- "fry"
fry_tab$contrast <- "stacked"
cam_tab <- camera(stacked$v, idx, stacked$design, contrast = stacked$coef, inter.gene.cor = 0.01)
cam_tab$family <- rownames(cam_tab)
cam_tab$test <- "camera"
cam_tab$contrast <- "stacked"

# Within-cohort fry as well. df is small; p-values are coarse.
fry_co <- list(fry_tab)
cam_co <- list(cam_tab)
for (co in COHORTS) {
  v <- cohort_fit[[co]]$v
  design <- cohort_fit[[co]]$design
  ix <- lapply(families, function(g) intersect(g, rownames(v)))
  ix <- ix[lengths(ix) >= 3]
  fr <- fry(v, ix, design, contrast = "groupQ4")
  fr$family <- rownames(fr)
  fr$test <- "fry"
  fr$contrast <- co
  ca <- camera(v, ix, design, contrast = "groupQ4", inter.gene.cor = 0.01)
  ca$family <- rownames(ca)
  ca$test <- "camera"
  ca$contrast <- co
  fry_co[[length(fry_co) + 1]] <- fr
  cam_co[[length(cam_co) + 1]] <- ca
}
fry_all <- do.call(rbind, lapply(fry_co, function(x) {
  x$family <- rownames(x)
  x
}))
# rbind of fry may duplicate family column; rebuild cleanly.
clean_gst <- function(lst) {
  do.call(rbind, lapply(lst, function(x) {
    data.frame(
      contrast = x$contrast,
      test = x$test,
      family = x$family,
      NGenes = x$NGenes,
      Direction = x$Direction,
      PValue = x$PValue,
      FDR = x$FDR,
      stringsAsFactors = FALSE
    )
  }))
}
gst <- rbind(clean_gst(fry_co), clean_gst(cam_co))
rownames(gst) <- NULL

# --- gene-level RE meta of cohort voom and edgeR logFCs --------------------

meta_genes <- function(method) {
  mats_y <- list()
  mats_s <- list()
  genes <- NULL
  for (co in COHORTS) {
    df <- if (method == "limma-voom") cohort_fit[[co]]$voom else cohort_fit[[co]]$edger
    if (is.null(genes)) genes <- df$gene
    # Align to the first cohort, then add genes only in later cohorts below.
    mats_y[[co]] <- df$logFC
    mats_s[[co]] <- df$se
    names(mats_y[[co]]) <- df$gene
    names(mats_s[[co]]) <- df$gene
  }
  genes <- sort(unique(unlist(lapply(mats_y, names))))
  yi <- sapply(COHORTS, function(co) mats_y[[co]][genes])
  sei <- sapply(COHORTS, function(co) mats_s[[co]][genes])
  rownames(yi) <- genes
  rownames(sei) <- genes
  colnames(yi) <- COHORTS
  colnames(sei) <- COHORTS
  out <- dl_mat(yi, sei)
  out$gene <- genes
  out$method <- method
  out$fdr_kge2 <- NA_real_
  out$fdr_kge3 <- NA_real_
  out$fdr_kge2[out$k >= 2] <- bh(out$p_t[out$k >= 2])
  out$fdr_kge3[out$k >= 3] <- bh(out$p_t[out$k >= 3])
  # Spot-check 25 genes against dl_one.
  set.seed(2)
  take <- sample(which(out$k >= 2), 25)
  for (i in take) {
    one <- dl_one(yi[i, ], sei[i, ])
    if (abs(one$b - out$logFC[i]) > 1e-6) stop("meta row mismatch for ", genes[i])
  }
  out
}

meta_voom <- meta_genes("limma-voom")
meta_edger <- meta_genes("edgeR-QL")
note(sprintf(
  "voom meta genes k>=3: %d; k=4: %d",
  sum(meta_voom$k >= 3), sum(meta_voom$k == 4)
))

# Spot-check one real gene against metafor.
chk_gene <- "CLDN4"
if (chk_gene %in% meta_voom$gene) {
  i <- match(chk_gene, meta_voom$gene)
  yi <- sapply(COHORTS, function(co) {
    df <- cohort_fit[[co]]$voom
    df$logFC[match(chk_gene, df$gene)]
  })
  sei <- sapply(COHORTS, function(co) {
    df <- cohort_fit[[co]]$voom
    df$se[match(chk_gene, df$gene)]
  })
  mk <- metafor::rma(yi = yi, sei = sei, method = "DL")
  if (abs(as.numeric(mk$beta) - meta_voom$logFC[i]) > 1e-6) {
    stop("CLDN4 meta != metafor")
  }
  note(sprintf("CLDN4 voom RE logFC=%.3f (metafor match) p_t=%.3g k=%d",
               meta_voom$logFC[i], meta_voom$p_t[i], meta_voom$k[i]))
}

# --- LOO and n_malignant<50 sensitivity on the stacked model ---------------

loo_rows <- list()
run_loo <- function(drop_cohort, drop_patients, tag) {
  u <- u_all
  if (!is.null(drop_cohort)) u <- u[u$cohort != drop_cohort, ]
  if (!is.null(drop_patients)) u <- u[!u$patient %in% drop_patients, ]
  mat <- mat_all[, u$patient, drop = FALSE]
  fit <- fit_stacked(mat, u, tag)
  fr_idx <- lapply(families, function(g) intersect(g, rownames(fit$v)))
  fr_idx <- fr_idx[lengths(fr_idx) >= 3]
  fr <- fry(fit$v, fr_idx, fit$design, contrast = fit$coef)
  for (fn in rownames(fr)) {
    genes <- intersect(families[[fn]], fit$voom$gene)
    lfc <- fit$voom$logFC[match(genes, fit$voom$gene)]
    loo_rows[[length(loo_rows) + 1]] <<- data.frame(
      sensitivity = tag,
      family = fn,
      n_q1 = sum(u$quartile == "Q1"),
      n_q4 = sum(u$quartile == "Q4"),
      n = nrow(u),
      n_genes = length(genes),
      median_logFC = median(lfc, na.rm = TRUE),
      mean_logFC = mean(lfc, na.rm = TRUE),
      n_up = sum(lfc > 0, na.rm = TRUE),
      n_down = sum(lfc < 0, na.rm = TRUE),
      fry_direction = fr[fn, "Direction"],
      fry_p = fr[fn, "PValue"],
      stringsAsFactors = FALSE
    )
  }
  cl <- fit$voom[fit$voom$gene == "CLDN4", ]
  loo_rows[[length(loo_rows) + 1]] <<- data.frame(
    sensitivity = tag, family = "CLDN4",
    n_q1 = sum(u$quartile == "Q1"), n_q4 = sum(u$quartile == "Q4"), n = nrow(u),
    n_genes = nrow(cl),
    median_logFC = if (nrow(cl)) cl$logFC else NA_real_,
    mean_logFC = if (nrow(cl)) cl$logFC else NA_real_,
    n_up = NA_integer_, n_down = NA_integer_,
    fry_direction = NA_character_, fry_p = if (nrow(cl)) cl$p else NA_real_,
    stringsAsFactors = FALSE
  )
}

for (co in COHORTS) run_loo(co, NULL, paste0("drop_", co))
small <- u_all$patient[u_all$n_malignant < 50]
note("DE units with n_malignant<50:", if (length(small)) paste(small, collapse = ",") else "none")
if (length(small)) run_loo(NULL, small, "drop_n_malignant_lt50")
loo <- do.call(rbind, loo_rows)
rownames(loo) <- NULL

# --- gene direction summaries ----------------------------------------------

gene_family <- function(gene) {
  hit <- names(families)[vapply(families, function(g) gene %in% g, logical(1))]
  if (!length(hit)) return("other")
  paste(hit, collapse = "|")
}

as_gene_table <- function(df) {
  if (!"p" %in% names(df)) df$p <- df$p_t
  if (!"fdr" %in% names(df)) df$fdr <- df$fdr_kge3
  df
}

direction_block <- function(df, contrast, method) {
  df <- as_gene_table(df)
  rows <- lapply(FAM_ORDER, function(fn) {
    genes <- intersect(families[[fn]], df$gene)
    sub <- df[match(genes, df$gene), ]
    sub <- sub[is.finite(sub$logFC), ]
    n_up <- sum(sub$logFC > 0)
    n_down <- sum(sub$logFC < 0)
    exp <- EXPECT[[fn]]
    succ <- if (exp > 0) n_up else n_down
    bt <- binom.test(succ, n_up + n_down, p = 0.5, alternative = "greater")
    data.frame(
      contrast = contrast, method = method, family = fn,
      n_genes = nrow(sub),
      n_up = n_up, n_down = n_down,
      n_p05 = sum(sub$p < 0.05),
      n_p05_up = sum(sub$p < 0.05 & sub$logFC > 0),
      n_p05_down = sum(sub$p < 0.05 & sub$logFC < 0),
      n_fdr05 = sum(sub$fdr < 0.05),
      median_logFC = median(sub$logFC),
      mean_logFC = mean(sub$logFC),
      sign_p_expected = bt$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

dir_rows <- rbind(
  direction_block(stacked$voom, "stacked", "limma-voom"),
  direction_block(stacked$edger, "stacked", "edgeR-QL"),
  direction_block(meta_voom[meta_voom$k >= 3, ], "re_meta_kge3", "limma-voom"),
  direction_block(meta_edger[meta_edger$k >= 3, ], "re_meta_kge3", "edgeR-QL")
)
for (co in COHORTS) {
  dir_rows <- rbind(
    dir_rows,
    direction_block(cohort_fit[[co]]$voom, co, "limma-voom"),
    direction_block(cohort_fit[[co]]$edger, co, "edgeR-QL")
  )
}
dir_rows$fdr_sign <- ave(dir_rows$sign_p_expected, dir_rows$contrast, dir_rows$method, FUN = bh)

# Agreement of the two engines.
agree_one <- function(a, b, label) {
  m <- merge(a[, c("gene", "logFC")], b[, c("gene", "logFC")], by = "gene", suffixes = c("_voom", "_edger"))
  same <- sign(m$logFC_voom) == sign(m$logFC_edger)
  data.frame(
    contrast = label,
    n = nrow(m),
    cor = cor(m$logFC_voom, m$logFC_edger, use = "complete.obs"),
    frac_same_sign = mean(same, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
agree <- rbind(
  agree_one(stacked$voom, stacked$edger, "stacked"),
  agree_one(meta_voom, meta_edger, "re_meta")
)
for (co in COHORTS) {
  agree <- rbind(agree, agree_one(cohort_fit[[co]]$voom, cohort_fit[[co]]$edger, co))
}

# --- honest n --------------------------------------------------------------

n_rows <- lapply(COHORTS, function(co) {
  u <- units[units$cohort == co, ]
  d <- de_units(co)
  dropped <- sort(u$patient[u$quartile %in% c("Q1", "Q4") & !u$in_count_matrix])
  data.frame(
    cohort = co,
    unit = u$unit[1],
    malig_def = u$malig_def[1],
    n_quartile_vector = nrow(u),
    n_q1_vector = sum(u$quartile == "Q1"),
    n_q4_vector = sum(u$quartile == "Q4"),
    n_q1_de = sum(d$quartile == "Q1"),
    n_q4_de = sum(d$quartile == "Q4"),
    n_de = nrow(d),
    thin_arm = min(sum(d$quartile == "Q1"), sum(d$quartile == "Q4")) < 3,
    patients_q1_de = paste(sort(d$patient[d$quartile == "Q1"]), collapse = ","),
    patients_q4_de = paste(sort(d$patient[d$quartile == "Q4"]), collapse = ","),
    dropped_q_tail = paste(dropped, collapse = ","),
    min_n_malignant_de = min(d$n_malignant),
    median_n_malignant_de = median(d$n_malignant),
    stringsAsFactors = FALSE
  )
})
n_honest <- do.call(rbind, n_rows)
# Composition that can confound a malignant program (GSE205335 histology).
comp <- units[units$cohort == "GSE205335" & units$quartile %in% c("Q1", "Q4"), ]
comp$in_de <- comp$in_count_matrix
comp_tab <- as.data.frame(table(comp$quartile, comp$cancer_subtype, comp$in_de), stringsAsFactors = FALSE)
colnames(comp_tab) <- c("quartile", "cancer_subtype", "in_de", "n")
comp_tab <- comp_tab[comp_tab$n > 0, ]

# --- key genes -------------------------------------------------------------

pull_key <- function(df, contrast, method) {
  df <- as_gene_table(df)
  sub <- df[df$gene %in% KEY, ]
  if (!nrow(sub)) return(NULL)
  data.frame(
    contrast = contrast, method = method, gene = sub$gene,
    family = vapply(sub$gene, gene_family, character(1)),
    logFC = sub$logFC, se = sub$se, p = sub$p, fdr = sub$fdr,
    stringsAsFactors = FALSE
  )
}
key_df <- rbind(
  pull_key(stacked$voom, "stacked", "limma-voom"),
  pull_key(stacked$edger, "stacked", "edgeR-QL"),
  pull_key(meta_voom, "re_meta", "limma-voom"),
  pull_key(meta_edger, "re_meta", "edgeR-QL")
)
for (co in COHORTS) {
  key_df <- rbind(
    key_df,
    pull_key(cohort_fit[[co]]$voom, co, "limma-voom"),
    pull_key(cohort_fit[[co]]$edger, co, "edgeR-QL")
  )
}

# Family genes, stacked + meta, both engines. Full gene tables are gzipped.
fam_genes <- unique(unlist(families))
fam_de <- rbind(
  stacked$voom[stacked$voom$gene %in% fam_genes | stacked$voom$gene == "CLDN4", ],
  stacked$edger[stacked$edger$gene %in% fam_genes | stacked$edger$gene == "CLDN4", ]
)
fam_de$family <- vapply(fam_de$gene, gene_family, character(1))

write_gz <- function(df, path) {
  con <- gzfile(path, "wt")
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}

# --- write tables ----------------------------------------------------------

write.table(n_honest, file.path(TAB, "n_honest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comp_tab, file.path(TAB, "gse205335_histology_qtails.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(effects, file.path(TAB, "family_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gst, file.path(TAB, "geneset_fry_camera.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dir_rows, file.path(TAB, "family_direction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(agree, file.path(TAB, "engine_agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(loo, file.path(TAB, "sensitivity_loo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(key_df, file.path(TAB, "key_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(scores, file.path(TAB, "family_scores_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fam_de, file.path(TAB, "de_family_genes_stacked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(units, file.path(TAB, "units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

meta_out <- function(df) {
  df$family <- vapply(df$gene, gene_family, character(1))
  df[, c("gene", "family", "method", "k", "logFC", "se", "p_z", "p_t", "fdr_kge2", "fdr_kge3",
         "tau2", "I2", "max_weight_share", "dominant_cohort")]
}
write_gz(meta_out(meta_voom), file.path(TAB, "de_meta_voom.tsv.gz"))
write_gz(meta_out(meta_edger), file.path(TAB, "de_meta_edger.tsv.gz"))
write_gz(rbind(stacked$voom, stacked$edger), file.path(TAB, "de_stacked.tsv.gz"))
per_co <- do.call(rbind, gene_rows)
write_gz(per_co, file.path(TAB, "de_by_cohort.tsv.gz"))

# Plain extract of meta rows for family genes + CLDN4, so the PR is readable without gzip.
meta_fam <- rbind(meta_out(meta_voom), meta_out(meta_edger))
meta_fam <- meta_fam[meta_fam$family != "other" | meta_fam$gene == "CLDN4", ]
write.table(meta_fam, file.path(TAB, "de_meta_family_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(TAB, "session_info.txt"))
cat(R.version.string, "\n")
cat("limma", as.character(packageVersion("limma")),
    "edgeR", as.character(packageVersion("edgeR")),
    "metafor", as.character(packageVersion("metafor")),
    "nlme", as.character(packageVersion("nlme")), "\n")
sessionInfo()
sink()

writeLines(notes, file.path(TAB, "run_notes.txt"))

# --- figures ----------------------------------------------------------------

draw_n <- function() {
  par(mar = c(7, 5, 4, 1))
  labs <- sprintf("%s\n%s", n_honest$cohort, n_honest$unit)
  m <- rbind(n_honest$n_q1_de, n_honest$n_q4_de)
  bp <- barplot(
    m, beside = TRUE, names.arg = labs, col = c("#9ecae1", "#e6550d"),
    ylim = c(0, max(n_honest$n_quartile_vector) + 4),
    ylab = "Units in the Q4 vs Q1 contrast",
    main = "Honest n: quartile vector vs units that enter DE"
  )
  legend("topright", c("Q1 in DE", "Q4 in DE"), fill = c("#9ecae1", "#e6550d"), bty = "n")
  text(colMeans(bp), n_honest$n_quartile_vector + 0.6,
       sprintf("vector n=%d", n_honest$n_quartile_vector), cex = 0.75)
  mtext("GSE189357 Q4 n=2 is below 3. P4001 is Q1 on the 22-patient vector and is not in the UMI-sum. GSE131907 units are samples.",
        side = 1, line = 5.5, cex = 0.7)
}
save_both(file.path(FIG, "n_honest"), 8, 5.2, draw_n)

draw_strip <- function() {
  par(mfrow = c(1, 4), mar = c(4, 3, 3, 1), oma = c(0, 2, 2, 0))
  for (co in COHORTS) {
    u <- units[units$cohort == co, ]
    u <- u[order(u$cldn4_pct), ]
    cols <- ifelse(u$quartile == "Q4", "#e6550d", ifelse(u$quartile == "Q1", "#3182bd", "grey70"))
    pch <- ifelse(u$in_count_matrix, 16, 1)
    plot(seq_len(nrow(u)), u$cldn4_pct, col = cols, pch = pch,
         xlab = "Units ordered by CLDN4 %pos", ylab = "",
         main = co, las = 1)
    if (any(!u$in_count_matrix & u$quartile %in% c("Q1", "Q4"))) {
      miss <- which(!u$in_count_matrix & u$quartile %in% c("Q1", "Q4"))
      text(miss, u$cldn4_pct[miss], u$patient[miss], pos = 3, cex = 0.6)
    }
  }
  mtext("Malignant CLDN4 % positive. Open circle = Q tail absent from the UMI-sum.",
        outer = TRUE, cex = 0.8)
}
save_both(file.path(FIG, "cldn4_quartile_strip"), 10, 3.6, draw_strip)

draw_forest <- function() {
  # Cohort lm + RE meta + LMM + stacked OLS.
  keep_est <- c(
    paste0("cohort_", COHORTS),
    "re_meta", "re_meta_drop_thin", "lmm_cohort", "stacked_ols"
  )
  d <- effects[effects$estimator %in% keep_est, ]
  d$family <- factor(d$family, levels = FAM_ORDER)
  d$estimator <- factor(d$estimator, levels = rev(keep_est))
  d <- d[order(d$family, d$estimator), ]
  d$row <- seq_len(nrow(d))
  par(mar = c(5, 16, 4, 2))
  xlim <- range(c(d$ci_low, d$ci_high), finite = TRUE)
  xlim <- c(min(xlim[1], -1.5), max(xlim[2], 1.2))
  plot(NA, xlim = xlim, ylim = c(0.5, nrow(d) + 0.5), yaxt = "n",
       xlab = "Family-score log2 difference (Q4 - Q1)", ylab = "",
       main = "Malignant family score, CLDN4 Q4 vs Q1")
  abline(v = 0, lty = 2, col = "grey40")
  cols <- ifelse(grepl("^cohort_", d$estimator), COH_COLS[sub("^cohort_", "", d$estimator)], "black")
  pch <- ifelse(grepl("^cohort_", d$estimator), 16, 18)
  segments(d$ci_low, d$row, d$ci_high, d$row, col = cols)
  points(d$logFC, d$row, pch = pch, col = cols, cex = ifelse(pch == 18, 1.4, 0.9))
  labs <- sprintf("%s  %s", d$family, gsub("_", " ", d$estimator))
  axis(2, at = d$row, labels = labs, las = 1, cex.axis = 0.62)
  mtext("Positive = higher in CLDN4-high. Diamonds: RE meta (Knapp-Hartung), mixed model, stacked OLS.",
        side = 1, line = 3.6, cex = 0.7)
}
save_both(file.path(FIG, "forest_family_score"), 9, 11, draw_forest)

draw_volcano <- function(df, pcol, title, path) {
  draw <- function() {
    p <- df[[pcol]]
    x <- df$logFC
    y <- -log10(pmax(p, 1e-300))
    fam <- df$family
    if (is.null(fam)) fam <- vapply(df$gene, gene_family, character(1))
    col <- rep("#d0d0d0", nrow(df))
    for (fn in FAM_ORDER) {
      hit <- grepl(fn, fam, fixed = TRUE)
      col[hit] <- FAM_COLS[[fn]]
    }
    par(mar = c(5, 5, 4, 1))
    plot(x, y, pch = 16, cex = 0.35, col = col,
         xlab = "logFC (Q4 - Q1)", ylab = "-log10 p", main = title)
    abline(v = 0, lty = 2, col = "grey40")
    # Label CLDN4 and the strongest gene in each family.
    lab_genes <- "CLDN4"
    for (fn in FAM_ORDER) {
      hit <- grepl(fn, fam, fixed = TRUE) & is.finite(p)
      if (!any(hit)) next
      lab_genes <- c(lab_genes, df$gene[which(hit)[which.min(p[hit])]])
    }
    lab_genes <- unique(lab_genes)
    for (g in lab_genes) {
      i <- match(g, df$gene)
      if (is.na(i) || !is.finite(y[i])) next
      text(x[i], y[i], g, cex = 0.6, pos = 3)
    }
    legend("topright", c(FAM_ORDER, "other"),
           col = c(FAM_COLS[FAM_ORDER], "#d0d0d0"), pch = 16, bty = "n", cex = 0.75)
  }
  save_both(path, 7.2, 5.6, draw)
}

mv <- meta_voom[meta_voom$k >= 3, ]
mv$family <- vapply(mv$gene, gene_family, character(1))
draw_volcano(stacked$voom, "p", "Stacked limma-voom, cohort-adjusted (n = 18 vs 16)", file.path(FIG, "volcano_stacked_voom"))
draw_volcano(mv, "p_t", "Random-effects meta of cohort limma-voom logFCs (k >= 3)", file.path(FIG, "volcano_re_meta"))

draw_heat <- function() {
  # Median limma-voom logFC. Columns: 4 cohorts, stacked, meta.
  cols <- c(COHORTS, "stacked", "re_meta_kge3")
  mat <- matrix(NA_real_, nrow = length(FAM_ORDER), ncol = length(cols),
                dimnames = list(FAM_ORDER, cols))
  ann <- matrix("", nrow = length(FAM_ORDER), ncol = length(cols))
  for (i in seq_along(FAM_ORDER)) {
    for (j in seq_along(cols)) {
      hit <- dir_rows$family == FAM_ORDER[i] & dir_rows$contrast == cols[j] & dir_rows$method == "limma-voom"
      if (!any(hit)) next
      r <- dir_rows[hit, ][1, ]
      mat[i, j] <- r$median_logFC
      ann[i, j] <- sprintf("%.2f\n%d down / %d up", r$median_logFC, r$n_down, r$n_up)
    }
  }
  draw <- function() {
    par(mar = c(8, 8, 4, 6))
    lim <- max(abs(mat), na.rm = TRUE)
    lim <- max(lim, 0.5)
    cols_r <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
    z <- mat
    plot(NA, xlim = c(0.5, ncol(z) + 0.5), ylim = c(0.5, nrow(z) + 0.5),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         main = "Median limma-voom logFC by family")
    for (i in seq_len(nrow(z))) {
      for (j in seq_len(ncol(z))) {
        v <- z[i, j]
        ci <- 1 + round(100 * (v + lim) / (2 * lim))
        ci <- max(1, min(101, ci))
        rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = cols_r[ci], border = "white")
        text(j, i, ann[i, j], cex = 0.62)
      }
    }
    axis(1, at = seq_len(ncol(z)), labels = colnames(z), las = 2, cex.axis = 0.75)
    axis(2, at = seq_len(nrow(z)), labels = rownames(z), las = 1)
    mtext("Red = higher in CLDN4 Q4. Cell text is median logFC and genes down / up.", side = 1, line = 6.5, cex = 0.75)
  }
  save_both(file.path(FIG, "heatmap_family_median_logfc"), 8.5, 5.2, draw)
}
draw_heat()

draw_key <- function() {
  d <- key_df[key_df$contrast == "stacked" & key_df$method == "limma-voom", ]
  d <- d[is.finite(d$logFC) & is.finite(d$se), ]
  d <- d[order(d$logFC), ]
  d$row <- seq_len(nrow(d))
  draw <- function() {
    par(mar = c(5, 8, 3, 1))
    ci_l <- d$logFC - 1.96 * d$se
    ci_h <- d$logFC + 1.96 * d$se
    plot(NA, xlim = range(c(ci_l, ci_h)), ylim = c(0.5, nrow(d) + 0.5),
         yaxt = "n", xlab = "Stacked voom logFC (Q4 - Q1)", ylab = "",
         main = "Key genes, cohort-adjusted limma-voom")
    abline(v = 0, lty = 2, col = "grey40")
    primary <- sub("\\|.*", "", d$family)
    cols <- ifelse(primary %in% names(FAM_COLS), FAM_COLS[primary], "grey30")
    cols[d$gene == "CLDN4"] <- "black"
    segments(ci_l, d$row, ci_h, d$row, col = cols)
    points(d$logFC, d$row, pch = 16, col = cols)
    axis(2, at = d$row, labels = d$gene, las = 1, cex.axis = 0.7)
  }
  save_both(file.path(FIG, "forest_key_genes"), 7, 7, draw)
}
draw_key()

draw_agree <- function() {
  m <- merge(
    stacked$voom[, c("gene", "logFC")],
    stacked$edger[, c("gene", "logFC")],
    by = "gene", suffixes = c("_voom", "_edger")
  )
  fam <- vapply(m$gene, gene_family, character(1))
  draw <- function() {
    par(mar = c(5, 5, 3, 1))
    plot(m$logFC_voom, m$logFC_edger, pch = 16, cex = 0.25, col = "#d0d0d0",
         xlab = "limma-voom logFC", ylab = "edgeR QL logFC",
         main = "Stacked engines agree on the Q4 coefficient")
    abline(0, 1, lty = 2)
    for (fn in FAM_ORDER) {
      hit <- grepl(fn, fam, fixed = TRUE)
      points(m$logFC_voom[hit], m$logFC_edger[hit], pch = 16, cex = 0.45, col = FAM_COLS[[fn]])
    }
    legend("topleft", FAM_ORDER, col = FAM_COLS[FAM_ORDER], pch = 16, bty = "n", cex = 0.75)
  }
  save_both(file.path(FIG, "agreement_voom_edger"), 6.2, 6, draw)
}
draw_agree()

note("wrote tables to ", TAB, " and figures to ", FIG)
message("DONE")

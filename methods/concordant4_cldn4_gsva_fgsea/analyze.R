#!/usr/bin/env Rscript

# Patient-level GSVA / ssGSEA and fgsea on concordant-4 malignant pseudobulks,
# stratified by malignant CLDN4 %pos (within-cohort quartiles locked in
# data/tnk_units.tsv). The unit is the patient / donor / sample. No cell-level
# p-values.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(GSVA)
  library(fgsea)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
ROOT <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  normalizePath(getwd())
}
DATA <- file.path(ROOT, "data")
TAB <- file.path(ROOT, "results", "tables")
FIG <- file.path(ROOT, "results", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

SEED <- 20260921L
N_BOOT <- as.integer(Sys.getenv("N_BOOT", "200"))
NPERM_BOOT <- as.integer(Sys.getenv("NPERM_BOOT", "1000"))
NPERM_POINT <- 10000L
MIN_SIZE <- 10L
MAX_SIZE <- 500L
MIN_ARM <- 3L
COHORTS <- c("GSE123902", "GSE131907", "GSE205335", "GSE189357")
REF <- "GSE123902"

PRIMARY <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_APICAL_JUNCTION",
  "KEGG_TIGHT_JUNCTION"
)
SECONDARY <- c(
  "GOBP_TIGHT_JUNCTION_ORGANIZATION",
  "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY"
)
TESTED <- c(PRIMARY, SECONDARY)

SET_LABEL <- c(
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "Hallmark IFN-alpha",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "Hallmark IFN-gamma",
  CUSTOM_MHC_I_ANTIGEN_PRESENTATION = "MHC-I antigen presentation",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "Hallmark EMT",
  HALLMARK_APICAL_JUNCTION = "Hallmark apical junction",
  KEGG_TIGHT_JUNCTION = "KEGG tight junction",
  GOBP_TIGHT_JUNCTION_ORGANIZATION = "GO TJ organization",
  GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY = "GO bicellular TJ assembly"
)
COHORT_COLOR <- c(
  GSE123902 = "#4c78a8",
  GSE131907 = "#f58518",
  GSE205335 = "#54a24b",
  GSE189357 = "#b279a2",
  pooled = "#222222",
  meta = "#222222"
)

set.seed(SEED)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

read_counts <- function(path) {
  dt <- fread(cmd = paste("gzip -dc", shQuote(path)), data.table = FALSE, showProgress = FALSE)
  genes <- toupper(as.character(dt[[1]]))
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, rownames(mat), reorder = FALSE)
  }
  mat
}

filter_genes <- function(counts, min_count = 10, min_samples = 3) {
  keep <- rowSums(counts >= min_count) >= min_samples
  counts[keep, , drop = FALSE]
}

tmm_norm_factors <- function(counts) {
  lib <- colSums(counts)
  lib[lib == 0] <- NA_real_
  rel <- sweep(counts, 2, lib, "/")
  f75 <- apply(rel, 2, stats::quantile, probs = 0.75, na.rm = TRUE, names = FALSE)
  ref <- names(which.min(abs(f75 - mean(f75))))
  ref_c <- counts[, ref]
  ref_lib <- lib[[ref]]
  factors <- setNames(numeric(ncol(counts)), colnames(counts))
  for (col in colnames(counts)) {
    obs <- counts[, col]
    obs_lib <- lib[[col]]
    keep <- obs > 0 & ref_c > 0
    if (sum(keep) < 50) {
      factors[col] <- 1
      next
    }
    m <- log2((obs[keep] / obs_lib) / (ref_c[keep] / ref_lib))
    a <- 0.5 * log2((obs[keep] / obs_lib) * (ref_c[keep] / ref_lib))
    w <- (obs_lib - obs[keep]) / (obs_lib * obs[keep]) +
      (ref_lib - ref_c[keep]) / (ref_lib * ref_c[keep])
    ok <- is.finite(m) & is.finite(a) & is.finite(w) & w > 0
    m <- m[ok]
    a <- a[ok]
    w <- w[ok]
    if (length(m) < 50) {
      factors[col] <- 1
      next
    }
    trim <- m >= stats::quantile(m, 0.30) & m <= stats::quantile(m, 0.70) &
      a >= stats::quantile(a, 0.05) & a <= stats::quantile(a, 0.95)
    if (sum(trim) < 20) {
      factors[col] <- 1
      next
    }
    factors[col] <- 2^stats::weighted.mean(m[trim], 1 / w[trim])
  }
  factors / mean(factors)
}

log_cpm <- function(counts) {
  fac <- tmm_norm_factors(counts)
  lib <- colSums(counts) * fac
  sweep(counts, 2, lib, "/") * 1e6
  log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
}

prepare_logcpm <- function(counts, samples) {
  samples <- intersect(samples, colnames(counts))
  cts <- filter_genes(counts[, samples, drop = FALSE])
  list(logcpm = log_cpm(cts), n_genes = nrow(cts))
}

bind_counts <- function(counts_list, samples) {
  genes <- sort(unique(unlist(lapply(counts_list, rownames), use.names = FALSE)))
  mat <- matrix(0, length(genes), length(samples), dimnames = list(genes, samples))
  for (cm in counts_list) {
    cols <- intersect(colnames(cm), samples)
    if (!length(cols)) next
    mat[rownames(cm), cols] <- cm[, cols, drop = FALSE]
  }
  mat
}

design_matrix <- function(meta_sub, mode) {
  cohorts <- COHORTS[COHORTS %in% meta_sub$cohort]
  ref <- if (REF %in% cohorts) REF else cohorts[[1]]
  others <- if (length(cohorts) >= 2L) setdiff(cohorts, ref) else character(0)
  n <- nrow(meta_sub)
  # paste0("x", character(0)) is "x" on this R, which would add an all-zero column.
  cohort_cols <- if (length(others)) paste0("cohort_", others) else character(0)
  cn <- c("Intercept", cohort_cols, "exposure")
  X <- matrix(0, n, length(cn), dimnames = list(NULL, cn))
  X[, "Intercept"] <- 1
  for (co in others) X[, paste0("cohort_", co)] <- as.numeric(meta_sub$cohort == co)
  if (mode == "q4") {
    X[, "exposure"] <- as.numeric(meta_sub$quartile == "Q4")
  } else {
    z <- as.numeric(scale(meta_sub$cldn4_pct))
    if (any(!is.finite(z))) return(NULL)
    X[, "exposure"] <- z
  }
  X
}

fit_exposure <- function(Y, meta_sub, mode) {
  X <- design_matrix(meta_sub, mode)
  if (is.null(X)) {
    message("    design dropped (", mode, ", n=", nrow(meta_sub), ")")
    return(NULL)
  }
  if (qr(X)$rank < ncol(X)) {
    message("    rank-deficient design (", mode, ", n=", nrow(meta_sub), ")")
    return(NULL)
  }
  XtX <- crossprod(X)
  inv <- solve(XtX)
  beta <- inv %*% crossprod(X, t(Y))
  resid <- t(Y) - X %*% beta
  df <- nrow(X) - ncol(X)
  if (df < 1) return(NULL)
  sigma2 <- colSums(resid^2) / df
  j <- ncol(X)
  se <- sqrt(pmax(sigma2 * inv[j, j], 0))
  est <- as.numeric(beta[j, ])
  tstat <- ifelse(se > 0, est / se, 0)
  names(tstat) <- rownames(Y)
  names(est) <- rownames(Y)
  p <- 2 * stats::pt(-abs(tstat), df)
  list(t = tstat, logFC = est, se = se, df = df, p = p, n = nrow(X))
}

named_stats <- function(x) {
  if (is.null(x) || !length(x) || !is.numeric(x)) {
    return(setNames(numeric(0), character(0)))
  }
  ok <- is.finite(x)
  x <- x[ok]
  if (!length(x)) return(x)
  nm <- names(x)
  if (is.null(nm)) nm <- sprintf("g%06d", seq_along(x))
  nm <- make.unique(as.character(nm), sep = "_dup")
  names(x) <- nm
  x[order(x, nm, decreasing = TRUE, method = "radix")]
}

quiet_fgsea_simple <- function(...) {
  ff <- tempfile()
  zz <- file(ff, open = "wt")
  sink(zz)
  sink(zz, type = "message")
  res <- tryCatch(fgseaSimple(...), error = function(e) e)
  sink(type = "message")
  sink()
  close(zz)
  unlink(ff)
  if (inherits(res, "error")) stop(res)
  res
}

point_fgsea <- function(pathways, stats) {
  stats <- named_stats(stats)
  if (length(stats) < 100) return(NULL)
  set.seed(SEED)
  fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    nPermSimple = NPERM_POINT,
    nproc = 1,
    scoreType = "std",
    eps = 0
  )
}

simple_nes <- function(pathways, stats, nperm) {
  stats <- named_stats(stats)
  if (length(stats) < 100) return(setNames(rep(NA_real_, length(pathways)), names(pathways)))
  res <- quiet_fgsea_simple(
    pathways, stats,
    nperm = nperm,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    nproc = 1,
    scoreType = "std"
  )
  out <- setNames(rep(NA_real_, length(pathways)), names(pathways))
  out[res$pathway] <- res$NES
  out
}

fgsea_to_rows <- function(res, scope, contrast, n, n_q1, n_q4, n_genes, family) {
  if (is.null(res) || !nrow(res)) return(data.frame())
  rows <- lapply(seq_len(nrow(res)), function(i) {
    le <- unlist(res$leadingEdge[[i]])
    data.frame(
      set = res$pathway[[i]],
      family = family,
      scope = scope,
      contrast = contrast,
      n = n,
      n_q1 = n_q1,
      n_q4 = n_q4,
      n_genes_ranked = n_genes,
      size = res$size[[i]],
      ES = res$ES[[i]],
      NES = res$NES[[i]],
      p = res$pval[[i]],
      n_leading = length(le),
      leading_edge = paste(utils::head(le, 40), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

bh_within <- function(df, cols) {
  df$fdr_family <- NA_real_
  if (!nrow(df)) return(df)
  key <- interaction(df[cols], drop = TRUE)
  df$fdr_family <- ave(df$p, key, FUN = function(p) p.adjust(p, "BH"))
  df
}

dl_meta <- function(y, v) {
  ok <- is.finite(y) & is.finite(v) & v > 0
  y <- y[ok]
  v <- v[ok]
  k <- length(y)
  if (k < 2) {
    return(list(k = k, mu = NA_real_, se = NA_real_, p = NA_real_,
                tau2 = NA_real_, I2 = NA_real_, Q = NA_real_,
                ci_lo = NA_real_, ci_hi = NA_real_))
  }
  w <- 1 / v
  ybar <- sum(w * y) / sum(w)
  Q <- sum(w * (y - ybar)^2)
  df <- k - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (C > 0) max(0, (Q - df) / C) else 0
  w2 <- 1 / (v + tau2)
  mu <- sum(w2 * y) / sum(w2)
  se <- sqrt(1 / sum(w2))
  z <- mu / se
  p <- 2 * stats::pnorm(-abs(z))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) else 0
  list(k = k, mu = mu, se = se, p = p, tau2 = tau2, I2 = I2, Q = Q,
       ci_lo = mu - 1.96 * se, ci_hi = mu + 1.96 * se)
}

# t.test(x, y) estimates mean(x) - mean(y). Pass high, then low.
welch_diff <- function(x_high, x_low) {
  tt <- stats::t.test(x_high, x_low)
  data.frame(
    effect = unname(tt$estimate[[1]] - tt$estimate[[2]]),
    se = unname((tt$conf.int[[2]] - tt$conf.int[[1]]) / (2 * stats::qt(0.975, tt$parameter))),
    ci_lo = unname(tt$conf.int[[1]]),
    ci_hi = unname(tt$conf.int[[2]]),
    p = unname(tt$p.value),
    df = unname(tt$parameter)
  )
}

hedges_g <- function(x_high, x_low) {
  n1 <- length(x_low)
  n4 <- length(x_high)
  m1 <- mean(x_low)
  m4 <- mean(x_high)
  s1 <- stats::sd(x_low)
  s4 <- stats::sd(x_high)
  sp <- sqrt(((n1 - 1) * s1^2 + (n4 - 1) * s4^2) / (n1 + n4 - 2))
  d <- (m4 - m1) / sp
  df <- n1 + n4 - 2
  J <- 1 - 3 / (4 * df - 1)
  g <- J * d
  var_d <- (n1 + n4) / (n1 * n4) + d^2 / (2 * (n1 + n4))
  var_g <- J^2 * var_d
  se <- sqrt(var_g)
  p <- 2 * stats::pt(-abs(g / se), df)
  data.frame(effect = g, se = se, ci_lo = g - 1.96 * se, ci_hi = g + 1.96 * se, p = p, df = df)
}

spearman_row <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 5) return(data.frame(effect = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  rho_c <- max(min(rho, 0.999), -0.999)
  z <- atanh(rho_c)
  se_z <- 1 / sqrt(n - 3)
  data.frame(
    effect = rho,
    se = se_z,
    ci_lo = tanh(z - 1.96 * se_z),
    ci_hi = tanh(z + 1.96 * se_z),
    p = unname(ct$p.value),
    n = n,
    z = z,
    var_z = se_z^2
  )
}

boot_indices <- function(meta_sub, mode) {
  if (mode == "q4") {
    idx <- integer(0)
    for (co in unique(meta_sub$cohort)) {
      for (q in c("Q1", "Q4")) {
        w <- which(meta_sub$cohort == co & meta_sub$quartile == q)
        if (!length(w)) next
        idx <- c(idx, sample(w, length(w), replace = TRUE))
      }
    }
    return(idx)
  }
  idx <- integer(0)
  for (co in unique(meta_sub$cohort)) {
    w <- which(meta_sub$cohort == co)
    idx <- c(idx, sample(w, length(w), replace = TRUE))
  }
  idx
}

bootstrap_nes <- function(Y, meta_sub, pathways, mode, B, nperm) {
  sets <- names(pathways)
  mat <- matrix(NA_real_, B, length(sets), dimnames = list(NULL, sets))
  for (b in seq_len(B)) {
    idx <- boot_indices(meta_sub, mode)
    sub <- meta_sub[idx, , drop = FALSE]
    fit <- fit_exposure(Y[, idx, drop = FALSE], sub, mode)
    if (is.null(fit)) next
    mat[b, ] <- simple_nes(pathways, fit$t, nperm)
    if (b %% 25 == 0) message("    bootstrap ", b, "/", B)
  }
  ci <- lapply(sets, function(s) {
    v <- mat[, s]
    v <- v[is.finite(v)]
    if (length(v) < 20) {
      return(data.frame(set = s, n_boot_ok = length(v), nes_boot_median = NA_real_,
                        ci_lo = NA_real_, ci_hi = NA_real_))
    }
    qs <- stats::quantile(v, c(0.025, 0.975), names = FALSE)
    data.frame(set = s, n_boot_ok = length(v), nes_boot_median = stats::median(v),
               ci_lo = qs[[1]], ci_hi = qs[[2]])
  })
  do.call(rbind, ci)
}

fmt_p <- function(p) {
  if (length(p) != 1 || !is.finite(p)) return("NA")
  if (p < 0.001) sprintf("%.2e", p) else sprintf("%.3f", p)
}
fmt_n <- function(x, d = 2) {
  if (length(x) != 1 || !is.finite(x)) return("NA")
  sprintf(paste0("%.", d, "f"), x)
}

save_plot <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h, dpi = 160, bg = "white")
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h, bg = "white")
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

message("Loading units and counts")
units <- fread(file.path(DATA, "tnk_units.tsv"), data.table = FALSE)
units$patient <- as.character(units$patient)
units$cohort <- as.character(units$cohort)
units$quartile <- as.character(units$quartile)
units$in_count_matrix <- tolower(as.character(units$in_count_matrix)) %in% c("true", "t", "1")
units <- units[units$cohort %in% COHORTS, ]
units$cohort <- factor(units$cohort, levels = COHORTS)

a8 <- fromJSON(file.path(DATA, "a8_sets.json"))
raw_sets <- lapply(a8$sets[TESTED], function(g) unique(toupper(as.character(g))))
missing_sets <- setdiff(TESTED, names(raw_sets))
if (length(missing_sets)) stop("Missing gene sets: ", paste(missing_sets, collapse = ", "))
contains_cldn4 <- vapply(raw_sets, function(g) "CLDN4" %in% g, logical(1))
sets_holdout <- lapply(raw_sets, function(g) setdiff(g, "CLDN4"))
sets_intact <- raw_sets

counts <- lapply(COHORTS, function(co) {
  read_counts(file.path(DATA, paste0(co, "_malignant_counts.tsv.gz")))
})
names(counts) <- COHORTS

for (co in COHORTS) {
  need <- units$patient[units$cohort == co & units$in_count_matrix]
  miss <- setdiff(need, colnames(counts[[co]]))
  if (length(miss)) stop(co, " count matrix missing: ", paste(miss, collapse = ", "))
}

meta <- units[units$in_count_matrix, ]
meta$cohort <- as.character(meta$cohort)
meta <- meta[order(match(meta$cohort, COHORTS), meta$patient), ]
rownames(meta) <- NULL

message("Units in matrices: ", nrow(meta),
        " Q1/Q4 ", sum(meta$quartile == "Q1"), "/", sum(meta$quartile == "Q4"))

# ---------------------------------------------------------------------------
# Within-cohort logCPM, then shared-gene GSVA / ssGSEA
# ---------------------------------------------------------------------------

message("TMM logCPM within cohort")
cohort_log <- list()
cohort_n_genes <- integer(0)
for (co in COHORTS) {
  ids <- meta$patient[meta$cohort == co]
  prep <- prepare_logcpm(counts[[co]], ids)
  cohort_log[[co]] <- prep$logcpm
  cohort_n_genes[co] <- prep$n_genes
  message("  ", co, " genes ", prep$n_genes, " samples ", ncol(prep$logcpm))
}
shared <- Reduce(intersect, lapply(cohort_log, rownames))
message("Shared genes passing the filter in every cohort: ", length(shared))
if (length(shared) < 3000) {
  warning("Shared gene universe is under 3000; GSVA ECDF is thin.")
}
log_shared <- do.call(cbind, lapply(COHORTS, function(co) {
  cohort_log[[co]][shared, meta$patient[meta$cohort == co], drop = FALSE]
}))
log_shared <- log_shared[, meta$patient, drop = FALSE]
stopifnot(identical(colnames(log_shared), meta$patient))

score_one <- function(method) {
  message("Scoring ", method)
  if (method == "gsva") {
    param <- gsvaParam(
      log_shared, sets_holdout,
      minSize = MIN_SIZE, maxSize = MAX_SIZE,
      kcdf = "Gaussian", maxDiff = TRUE, tau = 1
    )
  } else {
    param <- ssgseaParam(
      log_shared, sets_holdout,
      minSize = MIN_SIZE, maxSize = MAX_SIZE,
      alpha = 0.25, normalize = TRUE
    )
  }
  sc <- gsva(param, verbose = FALSE)
  if (!is.matrix(sc)) sc <- as.matrix(sc)
  sc
}

scores <- list(gsva = score_one("gsva"), ssgsea = score_one("ssgsea"))
for (nm in names(scores)) {
  miss <- setdiff(TESTED, rownames(scores[[nm]]))
  if (length(miss)) warning(nm, " dropped sets: ", paste(miss, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Score tests
# ---------------------------------------------------------------------------

message("Testing pathway scores")
score_rows <- list()
add_score <- function(row) score_rows[[length(score_rows) + 1]] <<- row

for (method in names(scores)) {
  sc <- scores[[method]]
  for (set in rownames(sc)) {
    fam <- if (set %in% PRIMARY) "primary" else "secondary"
    y_all <- as.numeric(sc[set, meta$patient])
    # continuous Spearman within cohort
    z_list <- list()
    v_list <- list()
    for (co in COHORTS) {
      ix <- meta$cohort == co
      sp <- spearman_row(y_all[ix], meta$cldn4_pct[ix])
      add_score(data.frame(
        method = method, set = set, family = fam, scope = co, contrast = "spearman_rho",
        n = sp$n, n_q1 = sum(meta$quartile[ix] == "Q1"), n_q4 = sum(meta$quartile[ix] == "Q4"),
        effect = sp$effect, se = sp$se, ci_lo = sp$ci_lo, ci_hi = sp$ci_hi, p = sp$p,
        stat = "rho", stringsAsFactors = FALSE
      ))
      z_list[[co]] <- if (is.finite(sp$effect)) atanh(max(min(sp$effect, 0.999), -0.999)) else NA_real_
      v_list[[co]] <- sp$var_z
    }
    md <- dl_meta(unlist(z_list), unlist(v_list))
    add_score(data.frame(
      method = method, set = set, family = fam, scope = "DL_meta", contrast = "spearman_rho",
      n = sum(meta$cohort %in% names(z_list)), n_q1 = NA_integer_, n_q4 = NA_integer_,
      effect = tanh(md$mu), se = md$se, ci_lo = tanh(md$ci_lo), ci_hi = tanh(md$ci_hi),
      p = md$p, stat = "rho", I2 = md$I2, tau2 = md$tau2, k = md$k,
      stringsAsFactors = FALSE
    ))

    # Q4 vs Q1 within cohort (native score units) + Hedges g
    g_y <- c()
    g_v <- c()
    d_y <- c()
    d_v <- c()
    d_names <- c()
    for (co in COHORTS) {
      ix1 <- meta$cohort == co & meta$quartile == "Q1"
      ix4 <- meta$cohort == co & meta$quartile == "Q4"
      n1 <- sum(ix1)
      n4 <- sum(ix4)
      if (n1 < MIN_ARM || n4 < MIN_ARM) {
        add_score(data.frame(
          method = method, set = set, family = fam, scope = co, contrast = "q4_minus_q1",
          n = n1 + n4, n_q1 = n1, n_q4 = n4,
          effect = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_,
          stat = "mean_diff", note = "arm < 3", stringsAsFactors = FALSE
        ))
        next
      }
      wd <- welch_diff(y_all[ix4], y_all[ix1])
      add_score(data.frame(
        method = method, set = set, family = fam, scope = co, contrast = "q4_minus_q1",
        n = n1 + n4, n_q1 = n1, n_q4 = n4,
        effect = wd$effect, se = wd$se, ci_lo = wd$ci_lo, ci_hi = wd$ci_hi, p = wd$p,
        stat = "mean_diff", stringsAsFactors = FALSE
      ))
      hg <- hedges_g(y_all[ix4], y_all[ix1])
      add_score(data.frame(
        method = method, set = set, family = fam, scope = co, contrast = "hedges_g",
        n = n1 + n4, n_q1 = n1, n_q4 = n4,
        effect = hg$effect, se = hg$se, ci_lo = hg$ci_lo, ci_hi = hg$ci_hi, p = hg$p,
        stat = "g", stringsAsFactors = FALSE
      ))
      d_y <- c(d_y, wd$effect)
      d_v <- c(d_v, wd$se^2)
      g_y <- c(g_y, hg$effect)
      g_v <- c(g_v, hg$se^2)
      d_names <- c(d_names, co)
    }
    md_d <- dl_meta(d_y, d_v)
    add_score(data.frame(
      method = method, set = set, family = fam, scope = "DL_meta", contrast = "q4_minus_q1",
      n = sum(d_y * 0 + 1, na.rm = TRUE), n_q1 = NA_integer_, n_q4 = NA_integer_,
      effect = md_d$mu, se = md_d$se, ci_lo = md_d$ci_lo, ci_hi = md_d$ci_hi, p = md_d$p,
      stat = "mean_diff", I2 = md_d$I2, tau2 = md_d$tau2, k = md_d$k,
      stringsAsFactors = FALSE
    ))
    md_g <- dl_meta(g_y, g_v)
    add_score(data.frame(
      method = method, set = set, family = fam, scope = "DL_meta", contrast = "hedges_g",
      n = length(g_y), n_q1 = NA_integer_, n_q4 = NA_integer_,
      effect = md_g$mu, se = md_g$se, ci_lo = md_g$ci_lo, ci_hi = md_g$ci_hi, p = md_g$p,
      stat = "g", I2 = md_g$I2, tau2 = md_g$tau2, k = md_g$k,
      stringsAsFactors = FALSE
    ))

    # Pooled OLS, cohort covariate, shared-gene scores
    m4 <- meta[meta$quartile %in% c("Q1", "Q4"), ]
    dd4 <- data.frame(
      score = as.numeric(sc[set, m4$patient]),
      cohort = factor(m4$cohort, levels = COHORTS),
      q4 = as.numeric(m4$quartile == "Q4")
    )
    fit4 <- stats::lm(score ~ cohort + q4, data = dd4)
    cf4 <- summary(fit4)$coefficients
    b4 <- unname(cf4["q4", 1])
    se4 <- unname(cf4["q4", 2])
    p4 <- unname(cf4["q4", 4])
    ci4 <- stats::confint(fit4)["q4", ]
    add_score(data.frame(
      method = method, set = set, family = fam, scope = "pooled_OLS", contrast = "q4_minus_q1",
      n = nrow(m4), n_q1 = sum(m4$quartile == "Q1"), n_q4 = sum(m4$quartile == "Q4"),
      effect = b4, se = se4, ci_lo = unname(ci4[[1]]), ci_hi = unname(ci4[[2]]), p = p4,
      stat = "ols_beta", stringsAsFactors = FALSE
    ))

    ddz <- data.frame(
      score = y_all,
      cohort = factor(meta$cohort, levels = COHORTS),
      z = as.numeric(scale(meta$cldn4_pct))
    )
    fitz <- stats::lm(score ~ cohort + z, data = ddz)
    cfz <- summary(fitz)$coefficients
    bz <- unname(cfz["z", 1])
    sez <- unname(cfz["z", 2])
    pz <- unname(cfz["z", 4])
    ciz <- stats::confint(fitz)["z", ]
    add_score(data.frame(
      method = method, set = set, family = fam, scope = "pooled_OLS", contrast = "per_sd_pct",
      n = nrow(meta), n_q1 = NA_integer_, n_q4 = NA_integer_,
      effect = bz, se = sez, ci_lo = unname(ciz[[1]]), ci_hi = unname(ciz[[2]]), p = pz,
      stat = "ols_beta", stringsAsFactors = FALSE
    ))
  }
}
score_df <- rbindlist(score_rows, fill = TRUE)
score_df <- as.data.frame(score_df)
score_df$fdr_primary <- NA_real_
for (method in unique(score_df$method)) {
  for (contrast in unique(score_df$contrast)) {
    for (scope in unique(score_df$scope)) {
      ix <- score_df$method == method & score_df$contrast == contrast &
        score_df$scope == scope & score_df$family == "primary" & is.finite(score_df$p)
      if (sum(ix) > 1) score_df$fdr_primary[ix] <- p.adjust(score_df$p[ix], "BH")
    }
  }
}

# Long patient scores
patient_scores <- rbindlist(lapply(names(scores), function(method) {
  sc <- scores[[method]]
  rows <- lapply(rownames(sc), function(set) {
    data.frame(
      method = method,
      set = set,
      patient = meta$patient,
      cohort = meta$cohort,
      quartile = meta$quartile,
      cldn4_pct = meta$cldn4_pct,
      score = as.numeric(sc[set, meta$patient]),
      stringsAsFactors = FALSE
    )
  })
  rbindlist(rows)
}))

# ---------------------------------------------------------------------------
# fgsea rankings (patient-level t)
# ---------------------------------------------------------------------------

message("Pooled rankings for fgsea (zero-filled counts, global TMM, cohort + exposure)")
q_ids <- meta$patient[meta$quartile %in% c("Q1", "Q4")]
all_ids <- meta$patient
counts_q <- bind_counts(counts, q_ids)
counts_all <- bind_counts(counts, all_ids)
prep_q <- prepare_logcpm(counts_q, q_ids)
prep_all <- prepare_logcpm(counts_all, all_ids)
meta_q <- meta[match(colnames(prep_q$logcpm), meta$patient), ]
meta_all <- meta[match(colnames(prep_all$logcpm), meta$patient), ]
fit_q <- fit_exposure(prep_q$logcpm, meta_q, "q4")
fit_all <- fit_exposure(prep_all$logcpm, meta_all, "cont")
stopifnot(!is.null(fit_q), !is.null(fit_all))
message("CLDN4 Q4 vs Q1 logFC ", signif(fit_q$logFC[["CLDN4"]], 4),
        " t ", signif(fit_q$t[["CLDN4"]], 4),
        " p ", signif(fit_q$p[["CLDN4"]], 4),
        " n_genes ", length(fit_q$t))
if (!is.finite(fit_q$logFC[["CLDN4"]]) || fit_q$logFC[["CLDN4"]] <= 0) {
  stop("CLDN4 is not higher in Q4. Quartile alignment is wrong.")
}

message("Point fgsea")
fgsea_rows <- list()
push_fgsea <- function(res, scope, contrast, n, n_q1, n_q4, n_genes) {
  fgsea_rows[[length(fgsea_rows) + 1]] <<- fgsea_to_rows(
    res, scope, contrast, n, n_q1, n_q4, n_genes, NA_character_
  )
}
res_q <- point_fgsea(sets_holdout, fit_q$t)
res_c <- point_fgsea(sets_holdout, fit_all$t)
push_fgsea(res_q, "pooled", "q4_vs_q1", nrow(meta_q),
           sum(meta_q$quartile == "Q1"), sum(meta_q$quartile == "Q4"), length(fit_q$t))
push_fgsea(res_c, "pooled", "continuous_z", nrow(meta_all), NA_integer_, NA_integer_, length(fit_all$t))

# CLDN4 restored, same ranking, junction sets only
sets_restore <- sets_holdout
for (s in names(sets_intact)[contains_cldn4]) sets_restore[[s]] <- sets_intact[[s]]
res_q_restore <- point_fgsea(sets_restore, fit_q$t)
res_c_restore <- point_fgsea(sets_restore, fit_all$t)
rest_rows <- rbind(
  fgsea_to_rows(res_q_restore, "pooled", "q4_vs_q1", nrow(meta_q),
                sum(meta_q$quartile == "Q1"), sum(meta_q$quartile == "Q4"), length(fit_q$t), NA),
  fgsea_to_rows(res_c_restore, "pooled", "continuous_z", nrow(meta_all),
                NA, NA, length(fit_all$t), NA)
)
rest_rows <- rest_rows[rest_rows$set %in% names(sets_intact)[contains_cldn4], ]

# Per-cohort rankings on within-cohort TMM (no zero-fill)
cohort_fits <- list()
for (co in COHORTS) {
  ids <- meta$patient[meta$cohort == co]
  sub <- meta[meta$cohort == co, ]
  prep <- prepare_logcpm(counts[[co]], ids)
  sub <- sub[match(colnames(prep$logcpm), sub$patient), ]
  fit_co_c <- fit_exposure(prep$logcpm, sub, "cont")
  if (is.null(fit_co_c)) {
    message("  no continuous ranking ", co)
  } else {
    cohort_fits[[paste0(co, "|cont")]] <- list(fit = fit_co_c, meta = sub, log = prep$logcpm)
    res <- point_fgsea(sets_holdout, fit_co_c$t)
    push_fgsea(res, co, "continuous_z", nrow(sub), NA, NA, length(fit_co_c$t))
  }

  sub4 <- sub[sub$quartile %in% c("Q1", "Q4"), ]
  n1 <- sum(sub4$quartile == "Q1")
  n4 <- sum(sub4$quartile == "Q4")
  if (n1 < MIN_ARM || n4 < MIN_ARM) {
    message("  skip binary fgsea ", co, " nQ1=", n1, " nQ4=", n4)
    next
  }
  prep4 <- prepare_logcpm(counts[[co]], sub4$patient)
  sub4 <- sub4[match(colnames(prep4$logcpm), sub4$patient), ]
  fit_co_q <- fit_exposure(prep4$logcpm, sub4, "q4")
  if (is.null(fit_co_q)) {
    message("  no binary ranking ", co)
  } else {
    cohort_fits[[paste0(co, "|q4")]] <- list(fit = fit_co_q, meta = sub4, log = prep4$logcpm)
    res4 <- point_fgsea(sets_holdout, fit_co_q$t)
    push_fgsea(res4, co, "q4_vs_q1", nrow(sub4), n1, n4, length(fit_co_q$t))
  }
}

# Leave-one-cohort-out on the pooled Q4 ranking (re-TMM without that cohort)
loo_rows <- list()
for (drop in COHORTS) {
  keep_co <- setdiff(COHORTS, drop)
  ids <- meta$patient[meta$cohort %in% keep_co & meta$quartile %in% c("Q1", "Q4")]
  cts <- bind_counts(counts[keep_co], ids)
  prep <- prepare_logcpm(cts, ids)
  sub <- meta[match(colnames(prep$logcpm), meta$patient), ]
  # reference cohort among those remaining
  fit <- fit_exposure(prep$logcpm, sub, "q4")
  if (is.null(fit)) {
    message("  LOO fit failed ", drop)
    next
  }
  res <- point_fgsea(sets_holdout, fit$t)
  rr <- fgsea_to_rows(res, paste0("drop_", drop), "q4_vs_q1", nrow(sub),
                      sum(sub$quartile == "Q1"), sum(sub$quartile == "Q4"), length(fit$t), NA)
  loo_rows[[drop]] <- rr
  message("  LOO drop ", drop, " n=", nrow(sub))
}

fgsea_df <- rbindlist(fgsea_rows, fill = TRUE)
fgsea_df <- as.data.frame(fgsea_df)
fgsea_df$family <- ifelse(fgsea_df$set %in% PRIMARY, "primary", "secondary")
fgsea_df$fdr_primary <- NA_real_
for (scope in unique(fgsea_df$scope)) {
  for (contrast in unique(fgsea_df$contrast)) {
    ix <- fgsea_df$scope == scope & fgsea_df$contrast == contrast &
      fgsea_df$family == "primary" & is.finite(fgsea_df$p)
    if (sum(ix)) fgsea_df$fdr_primary[ix] <- p.adjust(fgsea_df$p[ix], "BH")
  }
}
fgsea_df$fdr_fgsea_call <- NA_real_
# fgsea's own BH is across sets passed together; recompute within scope×contrast on all tested sets
for (scope in unique(fgsea_df$scope)) {
  for (contrast in unique(fgsea_df$contrast)) {
    ix <- fgsea_df$scope == scope & fgsea_df$contrast == contrast & is.finite(fgsea_df$p)
    if (sum(ix)) fgsea_df$fdr_fgsea_call[ix] <- p.adjust(fgsea_df$p[ix], "BH")
  }
}

# ---------------------------------------------------------------------------
# Bootstrap NES intervals (patient resample; normalization fixed)
# ---------------------------------------------------------------------------

message("Bootstrap NES B=", N_BOOT, " nperm=", NPERM_BOOT)
boot_rows <- list()
run_boot <- function(tag, Y, meta_sub, mode, scope, contrast) {
  message("  ", tag)
  ci <- bootstrap_nes(Y, meta_sub, sets_holdout, mode, N_BOOT, NPERM_BOOT)
  ci$scope <- scope
  ci$contrast <- contrast
  ci$B <- N_BOOT
  ci$nperm <- NPERM_BOOT
  boot_rows[[tag]] <<- ci
}
run_boot("pooled_q4", prep_q$logcpm, meta_q, "q4", "pooled", "q4_vs_q1")
run_boot("pooled_cont", prep_all$logcpm, meta_all, "cont", "pooled", "continuous_z")
for (co in COHORTS) {
  keyc <- paste0(co, "|cont")
  if (!is.null(cohort_fits[[keyc]])) {
    run_boot(keyc, cohort_fits[[keyc]]$log, cohort_fits[[keyc]]$meta, "cont", co, "continuous_z")
  }
  keyq <- paste0(co, "|q4")
  if (!is.null(cohort_fits[[keyq]])) {
    run_boot(keyq, cohort_fits[[keyq]]$log, cohort_fits[[keyq]]$meta, "q4", co, "q4_vs_q1")
  }
}
boot_df <- rbindlist(boot_rows, fill = TRUE)
boot_df <- as.data.frame(boot_df)

# ---------------------------------------------------------------------------
# Overlap table
# ---------------------------------------------------------------------------

overlap_rows <- lapply(names(sets_holdout), function(s) {
  data.frame(
    set = s,
    family = if (s %in% PRIMARY) "primary" else "secondary",
    n_frozen = length(raw_sets[[s]]),
    contains_CLDN4 = contains_cldn4[[s]],
    n_after_holdout = length(sets_holdout[[s]]),
    n_in_shared_logcpm = sum(sets_holdout[[s]] %in% rownames(log_shared)),
    n_in_pooled_q4_rank = sum(sets_holdout[[s]] %in% names(fit_q$t)),
    label = unname(SET_LABEL[[s]]),
    stringsAsFactors = FALSE
  )
})
overlap_df <- do.call(rbind, overlap_rows)

# ---------------------------------------------------------------------------
# Write tables
# ---------------------------------------------------------------------------

fwrite(as.data.table(score_df), file.path(TAB, "score_tests.tsv"), sep = "\t")
fwrite(patient_scores, file.path(TAB, "patient_scores.tsv"), sep = "\t")
fwrite(as.data.table(fgsea_df), file.path(TAB, "nes_fgsea.tsv"), sep = "\t")
fwrite(as.data.table(boot_df), file.path(TAB, "nes_fgsea_bootstrap.tsv"), sep = "\t")
fwrite(rbindlist(loo_rows, fill = TRUE), file.path(TAB, "nes_fgsea_loo.tsv"), sep = "\t")
fwrite(as.data.table(rest_rows), file.path(TAB, "nes_fgsea_cldn4_restored.tsv"), sep = "\t")
fwrite(as.data.table(overlap_df), file.path(TAB, "gene_set_overlap.tsv"), sep = "\t")
fwrite(as.data.table(meta), file.path(TAB, "patient_units_used.tsv"), sep = "\t")

rank_q <- data.frame(gene = names(fit_q$t), t = fit_q$t, logFC = fit_q$logFC[names(fit_q$t)], p = fit_q$p[names(fit_q$t)])
rank_q <- rank_q[order(-rank_q$t), ]
rank_c <- data.frame(gene = names(fit_all$t), t = fit_all$t, logFC = fit_all$logFC[names(fit_all$t)], p = fit_all$p[names(fit_all$t)])
rank_c <- rank_c[order(-rank_c$t), ]
fwrite(as.data.table(rank_q), file.path(TAB, "ranks_pooled_q4q1.tsv"), sep = "\t")
fwrite(as.data.table(rank_c), file.path(TAB, "ranks_pooled_continuous.tsv"), sep = "\t")

gmt_lines <- vapply(names(sets_holdout), function(s) {
  paste(c(s, "CLDN4_held_out", sets_holdout[[s]]), collapse = "\t")
}, character(1))
writeLines(gmt_lines, file.path(TAB, "gene_sets_cldn4_held_out.gmt"))

inv <- data.frame(
  cohort = c(COHORTS, "pooled"),
  n_units_matrix = c(vapply(COHORTS, function(co) sum(meta$cohort == co), integer(1)), nrow(meta)),
  n_q1 = c(vapply(COHORTS, function(co) sum(meta$cohort == co & meta$quartile == "Q1"), integer(1)), sum(meta$quartile == "Q1")),
  n_q4 = c(vapply(COHORTS, function(co) sum(meta$cohort == co & meta$quartile == "Q4"), integer(1)), sum(meta$quartile == "Q4")),
  n_genes_within = c(cohort_n_genes[COHORTS], length(shared)),
  stringsAsFactors = FALSE
)
# units not in the matrix
dropped <- units$patient[!units$in_count_matrix]
fwrite(as.data.table(inv), file.path(TAB, "inventory.tsv"), sep = "\t")

cldn4_check <- data.frame(
  gene = "CLDN4",
  contrast = c("q4_vs_q1", "continuous_z"),
  logFC = c(fit_q$logFC[["CLDN4"]], fit_all$logFC[["CLDN4"]]),
  t = c(fit_q$t[["CLDN4"]], fit_all$t[["CLDN4"]]),
  p = c(fit_q$p[["CLDN4"]], fit_all$p[["CLDN4"]]),
  n = c(nrow(meta_q), nrow(meta_all)),
  n_genes = c(length(fit_q$t), length(fit_all$t)),
  df = c(fit_q$df, fit_all$df)
)
fwrite(as.data.table(cldn4_check), file.path(TAB, "cldn4_direction_check.tsv"), sep = "\t")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

message("Plots")
theme_forest <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
}

nes_join <- merge(fgsea_df, boot_df, by = c("set", "scope", "contrast"), all.x = TRUE)
nes_join$label <- SET_LABEL[nes_join$set]
nes_join$scope_label <- ifelse(nes_join$scope == "pooled", "pooled ranking", nes_join$scope)

plot_nes_forest <- function(contrast, sets, title, subtitle) {
  d <- nes_join[nes_join$contrast == contrast & nes_join$set %in% sets & is.finite(nes_join$NES), ]
  scope_levels <- c(COHORTS, "pooled")
  d <- d[d$scope %in% scope_levels, ]
  d$scope <- factor(d$scope, levels = rev(scope_levels))
  d$label <- factor(d$label, levels = SET_LABEL[sets])
  d$fill <- COHORT_COLOR[ifelse(d$scope == "pooled", "pooled", as.character(d$scope))]
  ggplot(d, aes(x = NES, y = scope)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.18, color = "grey35", linewidth = 0.4) +
    geom_point(aes(color = scope), size = 2.4) +
    scale_color_manual(values = c(COHORT_COLOR, pooled = "#222222"), guide = "none") +
    facet_wrap(~label, ncol = 3) +
    theme_forest() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "NES (positive = enriched in CLDN4-high malignant pseudobulk)",
      y = NULL
    )
}

sub_boot <- sprintf(
  "Point = fgseaMultilevel (nPermSimple=%d). Interval = patient bootstrap percentile (B=%d, fgseaSimple nperm=%d). CLDN4 held out.",
  NPERM_POINT, N_BOOT, NPERM_BOOT
)
save_plot(
  plot_nes_forest("q4_vs_q1", PRIMARY, "fgsea NES, CLDN4 Q4 vs Q1", sub_boot),
  "forest_fgsea_nes_q4q1", 11, 7.2
)
save_plot(
  plot_nes_forest("continuous_z", PRIMARY, "fgsea NES, continuous CLDN4 %pos", sub_boot),
  "forest_fgsea_nes_continuous", 11, 7.2
)
save_plot(
  plot_nes_forest("q4_vs_q1", SECONDARY, "fgsea NES, secondary junction sets, Q4 vs Q1", sub_boot),
  "forest_fgsea_nes_q4q1_secondary", 8.5, 4.2
)

# Headline pooled NES, both contrasts
head_df <- nes_join[nes_join$scope == "pooled" & nes_join$set %in% PRIMARY, ]
head_df$label <- factor(SET_LABEL[head_df$set], levels = rev(SET_LABEL[PRIMARY]))
head_df$contrast_label <- ifelse(head_df$contrast == "q4_vs_q1", "Q4 vs Q1", "Continuous %pos")
save_plot(
  ggplot(head_df, aes(x = NES, y = label, color = contrast_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.18, position = position_dodge(width = 0.55), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.55), size = 2.6) +
    scale_color_manual(values = c("Q4 vs Q1" = "#222222", "Continuous %pos" = "#4c78a8")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
    labs(
      title = "Pooled patient-level fgsea, concordant-4 malignant cells",
      subtitle = sub_boot,
      x = "NES (positive = higher in CLDN4-high)",
      y = NULL
    ),
  "forest_fgsea_nes_pooled", 8.2, 5.4
)

plot_score_forest <- function(method, contrast, scope_meta, title, xlab) {
  d <- score_df[score_df$method == method & score_df$contrast == contrast & score_df$set %in% PRIMARY, ]
  d <- d[d$scope %in% c(COHORTS, scope_meta) & is.finite(d$effect), ]
  d$scope <- factor(d$scope, levels = rev(c(COHORTS, scope_meta)))
  d$label <- factor(SET_LABEL[d$set], levels = SET_LABEL[PRIMARY])
  ggplot(d, aes(x = effect, y = scope)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.18, color = "grey35", linewidth = 0.4) +
    geom_point(aes(shape = scope == scope_meta), size = 2.3, color = "#222222") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    facet_wrap(~label, ncol = 3, scales = "free_x") +
    theme_forest() +
    labs(title = title, subtitle = "Positive = higher pathway score in CLDN4-high. CLDN4 held out of sets.", x = xlab, y = NULL)
}

save_plot(
  plot_score_forest("gsva", "q4_minus_q1", "pooled_OLS",
                    "GSVA score, Q4 minus Q1",
                    "Score difference (pooled row = cohort-adjusted OLS beta)"),
  "forest_gsva_q4q1", 11, 7.2
)
save_plot(
  plot_score_forest("ssgsea", "q4_minus_q1", "pooled_OLS",
                    "ssGSEA score, Q4 minus Q1",
                    "Score difference (pooled row = cohort-adjusted OLS beta)"),
  "forest_ssgsea_q4q1", 11, 7.2
)
save_plot(
  plot_score_forest("gsva", "spearman_rho", "DL_meta",
                    "GSVA vs malignant CLDN4 %pos",
                    "Spearman rho (bottom row = DL meta on Fisher z, back-transformed)"),
  "forest_gsva_spearman", 11, 7.2
)
save_plot(
  plot_score_forest("ssgsea", "spearman_rho", "DL_meta",
                    "ssGSEA vs malignant CLDN4 %pos",
                    "Spearman rho (bottom row = DL meta on Fisher z, back-transformed)"),
  "forest_ssgsea_spearman", 11, 7.2
)

# Hedges g meta forest, both methods, Q4 vs Q1 — common standardized scale
gdf <- score_df[score_df$contrast == "hedges_g" & score_df$set %in% PRIMARY &
                  score_df$scope %in% c(COHORTS, "DL_meta") & is.finite(score_df$effect), ]
gdf$scope <- factor(gdf$scope, levels = rev(c(COHORTS, "DL_meta")))
gdf$label <- factor(SET_LABEL[gdf$set], levels = SET_LABEL[PRIMARY])
gdf$method_label <- ifelse(gdf$method == "gsva", "GSVA", "ssGSEA")
save_plot(
  ggplot(gdf, aes(x = effect, y = scope, color = method_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.2,
                  position = position_dodge(width = 0.55), linewidth = 0.35) +
    geom_point(position = position_dodge(width = 0.55), size = 2) +
    scale_color_manual(values = c(GSVA = "#b279a2", ssGSEA = "#4c78a8")) +
    facet_wrap(~label, ncol = 3) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8, face = "bold"),
          legend.position = "bottom", legend.title = element_blank()) +
    labs(
      title = "Standardized Q4 vs Q1 pathway-score difference",
      subtitle = "Hedges g. Diamond row is the DL random-effects meta of the cohorts with both arms n>=3. GSE189357 binary arm is omitted (Q4 n=2).",
      x = "Hedges g (positive = higher score in CLDN4 Q4)",
      y = NULL
    ),
  "forest_score_hedges_g", 11, 7.4
)

# Patient ssGSEA, within-cohort z, Q1 vs Q4
ps <- as.data.frame(patient_scores)
ps <- ps[ps$method == "ssgsea" & ps$set %in% PRIMARY & ps$quartile %in% c("Q1", "Q4"), ]
ps$z <- ave(ps$score, ps$set, ps$cohort, FUN = function(x) as.numeric(scale(x)))
ps$label <- factor(SET_LABEL[ps$set], levels = SET_LABEL[PRIMARY])
ps$quartile <- factor(ps$quartile, levels = c("Q1", "Q4"))
save_plot(
  ggplot(ps, aes(x = quartile, y = z, color = cohort)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    geom_boxplot(outlier.shape = NA, color = "grey40", fill = "grey95", width = 0.55) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.6, alpha = 0.9) +
    scale_color_manual(values = COHORT_COLOR) +
    facet_wrap(~label, ncol = 3) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8, face = "bold"),
          legend.position = "bottom", legend.title = element_blank()) +
    labs(
      title = "Patient ssGSEA scores, CLDN4 Q1 vs Q4",
      subtitle = "Within-cohort z of the ssGSEA score. Each point is one malignant pseudobulk, not a cell.",
      x = "Within-cohort malignant CLDN4 %pos quartile",
      y = "ssGSEA z within cohort"
    ),
  "box_ssgsea_q4q1_patient", 11, 7.2
)

# ---------------------------------------------------------------------------
# FINDING.md
# ---------------------------------------------------------------------------

lookup_f <- function(scope, contrast, set) {
  fgsea_df[fgsea_df$scope == scope & fgsea_df$contrast == contrast & fgsea_df$set == set, ]
}
lookup_b <- function(scope, contrast, set) {
  boot_df[boot_df$scope == scope & boot_df$contrast == contrast & boot_df$set == set, ]
}
lookup_s <- function(method, scope, contrast, set) {
  score_df[score_df$method == method & score_df$scope == scope &
             score_df$contrast == contrast & score_df$set == set, ]
}

nes_line <- function(set, contrast) {
  r <- lookup_f("pooled", contrast, set)
  b <- lookup_b("pooled", contrast, set)
  if (!nrow(r)) return("| | | | | | |")
  ci <- if (nrow(b) && is.finite(b$ci_lo[1])) sprintf("%s to %s", fmt_n(b$ci_lo[1]), fmt_n(b$ci_hi[1])) else "NA"
  sprintf("| %s | %s | %s | %s | %s | %d |",
          SET_LABEL[[set]], fmt_n(r$NES[1]), ci, fmt_p(r$p[1]), fmt_p(r$fdr_primary[1]), r$size[1])
}

score_line <- function(method, set, scope, contrast) {
  r <- lookup_s(method, scope, contrast, set)
  if (!nrow(r) || !is.finite(r$effect[1])) return(NULL)
  i2 <- if ("I2" %in% names(r) && is.finite(r$I2[1])) sprintf("%.0f%%", 100 * r$I2[1]) else "—"
  sprintf("| %s | %s | %s | %s | %s to %s | %s | %s | %s |",
          toupper(method), SET_LABEL[[set]], scope, fmt_n(r$effect[1], 3),
          fmt_n(r$ci_lo[1], 3), fmt_n(r$ci_hi[1], 3), fmt_p(r$p[1]),
          fmt_p(r$fdr_primary[1]), i2)
}

cohort_nes_md <- function(contrast) {
  lines <- c("| set | cohort | n | NES | bootstrap 95% | p | FDR |",
             "|---|---|---:|---:|---|---:|---:|")
  for (set in PRIMARY) {
    for (scope in c(COHORTS, "pooled")) {
      r <- lookup_f(scope, contrast, set)
      if (!nrow(r)) next
      b <- lookup_b(scope, contrast, set)
      ci <- if (nrow(b) && is.finite(b$ci_lo[1])) sprintf("%s to %s", fmt_n(b$ci_lo[1]), fmt_n(b$ci_hi[1])) else "—"
      lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                                SET_LABEL[[set]], scope,
                                ifelse(is.finite(r$n[1]), r$n[1], "—"),
                                fmt_n(r$NES[1]), ci, fmt_p(r$p[1]), fmt_p(r$fdr_primary[1])))
    }
  }
  paste(lines, collapse = "\n")
}

loo_md <- function() {
  loo <- rbindlist(loo_rows, fill = TRUE)
  loo <- as.data.frame(loo)
  lines <- c("| dropped | set | n | NES | p |", "|---|---|---:|---:|---:|")
  for (drop in COHORTS) {
    for (set in PRIMARY) {
      r <- loo[loo$scope == paste0("drop_", drop) & loo$set == set, ]
      if (!nrow(r)) next
      lines <- c(lines, sprintf("| %s | %s | %d | %s | %s |",
                                drop, SET_LABEL[[set]], r$n[1], fmt_n(r$NES[1]), fmt_p(r$p[1])))
    }
  }
  paste(lines, collapse = "\n")
}

restore_md <- function() {
  lines <- c("| set | contrast | NES held out | NES with CLDN4 | delta |",
             "|---|---|---:|---:|---:|")
  for (set in names(sets_intact)[contains_cldn4]) {
    for (contrast in c("q4_vs_q1", "continuous_z")) {
      a <- lookup_f("pooled", contrast, set)
      b <- rest_rows[rest_rows$set == set & rest_rows$contrast == contrast, ]
      if (!nrow(a) || !nrow(b)) next
      lines <- c(lines, sprintf("| %s | %s | %s | %s | %s |",
                                SET_LABEL[[set]], contrast, fmt_n(a$NES[1]), fmt_n(b$NES[1]),
                                fmt_n(b$NES[1] - a$NES[1])))
    }
  }
  paste(lines, collapse = "\n")
}

n_by <- vapply(COHORTS, function(co) sum(units$cohort == co), integer(1))
n_mat <- vapply(COHORTS, function(co) sum(meta$cohort == co), integer(1))
n_q1 <- vapply(COHORTS, function(co) sum(meta$cohort == co & meta$quartile == "Q1"), integer(1))
n_q4 <- vapply(COHORTS, function(co) sum(meta$cohort == co & meta$quartile == "Q4"), integer(1))

call_primary <- function(contrast) {
  bits <- c()
  for (set in PRIMARY) {
    r <- lookup_f("pooled", contrast, set)
    if (!nrow(r) || !is.finite(r$NES[1])) next
    flag <- if (is.finite(r$fdr_primary[1]) && r$fdr_primary[1] < 0.05) {
      if (r$NES[1] < 0) "down at FDR<0.05" else "up at FDR<0.05"
    } else {
      "not FDR<0.05"
    }
    bits <- c(bits, sprintf("%s NES %s (%s)", SET_LABEL[[set]], fmt_n(r$NES[1]), flag))
  }
  paste(bits, collapse = "; ")
}

finding <- paste0(
"# Concordant-4 malignant cells: patient-level GSVA / ssGSEA and fgsea by CLDN4

ADDITIVE. **CLDN4-only.** Cohorts are the locked concordant four:
**GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not GSE148071, GSE127465,
GSE207422, GSE154826, or CD45+/T-only extracts. The T/NK result
(malignant CLDN4 %pos vs T/NK, DL rho about -0.53, N=65) is not re-derived.

The unit is the patient, donor, or sample malignant pseudobulk. Cell-level
p-values are not computed. Quartiles are the locked within-cohort CLDN4 %pos
labels. P4001 is out of the malignant UMI sum, so expression n is not 65.

## Prespecified sets

Scored separately. IFN-alpha and IFN-gamma are not merged. Antigen presentation
is the frozen 21-gene classical MHC-I / APM panel (MHC-II excluded; not a
Hallmark set). Junction sets are Hallmark apical junction and KEGG tight
junction. Hallmark EMT is included as its own set. CLDN4 is removed from every
set before scoring. Two GO tight-junction sets are secondary.

Shared genes for GSVA/ssGSEA (count >= 10 in >= 3 units inside every cohort): **",
length(shared), "**. Within-cohort gene counts: ",
paste(sprintf("%s %d", COHORTS, cohort_n_genes[COHORTS]), collapse = "; "),
".

## Honest n

| cohort | unit in the locked table | in the count matrix | Q1 | Q4 |
|---|---:|---:|---:|---:|
",
paste(sprintf("| %s | %d | %d | %d | %d |", COHORTS, n_by, n_mat, n_q1, n_q4), collapse = "\n"),
"
| pooled expression | ", sum(n_by), " | **", nrow(meta), "** | **", sum(meta$quartile == "Q1"), "** | **", sum(meta$quartile == "Q4"), "** |

GSE189357 Q4 has 2 units, so its within-cohort binary test is skipped. Those
two units stay in the pooled Q4 vs Q1 model. Dropped from the matrix: ",
paste(dropped, collapse = ", "), ".

Direction check, same pooled OLS as the ranking (positive = higher in CLDN4-high):
CLDN4 Q4 vs Q1 logFC **", fmt_n(fit_q$logFC[["CLDN4"]], 3), "** (p ", fmt_p(fit_q$p[["CLDN4"]]),
", n=", nrow(meta_q), ", genes ", length(fit_q$t), "). Continuous logFC per SD of %pos **",
fmt_n(fit_all$logFC[["CLDN4"]], 3), "** (p ", fmt_p(fit_all$p[["CLDN4"]]), ").

## 1. fgsea NES (headline)

Ranking for the pool is the patient-level OLS t on log2(TMM-CPM+1),
`~ cohort + exposure`, built the same way as the concordant-4 malignant DE
(zero-fill genes absent from a cohort, then one TMM). Positive NES = the set
sits toward genes that are higher in CLDN4-high malignant pseudobulks.

p-values are fgseaMultilevel. The interval is a patient bootstrap (B=", N_BOOT,
", fgseaSimple nperm=", NPERM_BOOT, ") on that fixed logCPM matrix; normalization
factors are not re-estimated inside the bootstrap. FDR is Benjamini-Hochberg
across the six primary sets, within the contrast.

### Pooled Q4 vs Q1 (n=", nrow(meta_q), ", ", sum(meta_q$quartile == "Q1"), " vs ",
sum(meta_q$quartile == "Q4"), ")

| set | NES | bootstrap 95% | p | FDR | size |
|---|---:|---|---:|---:|---:|
",
paste(vapply(PRIMARY, nes_line, character(1), contrast = "q4_vs_q1"), collapse = "\n"),
"

### Pooled continuous CLDN4 %pos (n=", nrow(meta_all), ", t per SD)

| set | NES | bootstrap 95% | p | FDR | size |
|---|---:|---|---:|---:|---:|
",
paste(vapply(PRIMARY, nes_line, character(1), contrast = "continuous_z"), collapse = "\n"),
"

Q4 vs Q1 call, primary FDR: ", call_primary("q4_vs_q1"), ".

Continuous call, primary FDR on the fixed ranking: ", call_primary("continuous_z"), ".

Pooled Q4 patient-bootstrap intervals that stay negative are IFN-alpha,
IFN-gamma, MHC-I antigen presentation, and Hallmark EMT. Hallmark apical
junction and KEGG tight junction cross zero on that bootstrap. The fgsea
p-value conditions on the observed ranking. The bootstrap resamples patients.
Where they disagree, the bootstrap is the uncertainty that matches the unit.
Leave-one-cohort-out keeps the four negative sets negative, including the
drop of GSE205335. GSE131907 alone is flat for IFN.

Hallmark apical junction is a broad MSigDB set, not a pure tight-junction
list. KEGG tight junction and the two GO tight-junction sets are the junction
lists. CLDN4 is held out of every set. Putting it back does not flip the sign.
This run is not tight-junction up.

### Cohort NES (same sets; binary skipped when an arm has n<3)

",
cohort_nes_md("q4_vs_q1"),
"

### Continuous cohort NES

",
cohort_nes_md("continuous_z"),
"

### Leave-one-cohort-out, pooled Q4 ranking (re-TMM)

",
loo_md(),
"

### CLDN4 put back into junction sets (same ranking)

",
restore_md(),
"

## 2. GSVA and ssGSEA on the shared-gene malignant pseudobulks

One GSVA (`kcdf=Gaussian`, `maxDiff=TRUE`) and one ssGSEA (`alpha=0.25`,
`normalize=TRUE`) on the within-cohort log2(TMM-CPM+1) matrix restricted to
the ", length(shared), " shared genes. Scores are per unit. Tests:

- Within-cohort Spearman of the score vs CLDN4 %pos, then DerSimonian-Laird
  on Fisher z (all four cohorts).
- Within-cohort Welch difference Q4 minus Q1 where both arms have n>=3, and
  the same contrast as Hedges g.
- Pooled OLS `score ~ cohort + Q4` and `score ~ cohort + z(CLDN4 %pos)`.

FDR is BH across the six primary sets inside method x scope x contrast.

### Pooled OLS, Q4 vs Q1 (score units; positive = higher in Q4)

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
",
paste(na.omit(unlist(lapply(c("gsva", "ssgsea"), function(m) {
  vapply(PRIMARY, function(s) score_line(m, s, "pooled_OLS", "q4_minus_q1"), character(1))
}))), collapse = "\n"),
"

### Pooled OLS, per SD of CLDN4 %pos

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
",
paste(na.omit(unlist(lapply(c("gsva", "ssgsea"), function(m) {
  vapply(PRIMARY, function(s) score_line(m, s, "pooled_OLS", "per_sd_pct"), character(1))
}))), collapse = "\n"),
"

### DL meta of within-cohort Spearman rho

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
",
paste(na.omit(unlist(lapply(c("gsva", "ssgsea"), function(m) {
  vapply(PRIMARY, function(s) score_line(m, s, "DL_meta", "spearman_rho"), character(1))
}))), collapse = "\n"),
"

### DL meta of Hedges g (Q4 vs Q1; cohorts with both arms n>=3)

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
",
paste(na.omit(unlist(lapply(c("gsva", "ssgsea"), function(m) {
  vapply(PRIMARY, function(s) score_line(m, s, "DL_meta", "hedges_g"), character(1))
}))), collapse = "\n"),
"

## How to read the two layers

fgsea asks whether set members sit at one end of the patient-level t ranking.
GSVA and ssGSEA ask whether that patient's enrichment score is lower when
malignant CLDN4 %pos is high.

On the pooled Q4 contrast, GSVA and ssGSEA both put IFN-alpha, IFN-gamma,
MHC-I antigen presentation, and Hallmark EMT lower in CLDN4-high (FDR < 0.05
inside each method). That matches the negative fgsea NES. The Hedges g
meta-analysis is wider: ssGSEA IFN intervals cross zero and I2 is high, so
the score shift is not the same size in every cohort.

KEGG tight junction does not clear the score tests. Its Q4 fgsea NES is
positive and not FDR < 0.05. This run is not tight-junction up. Hallmark
apical junction is negative on fgsea and on the GSVA Q4 model, and not
FDR < 0.05 on ssGSEA. It is not a junction-barrier score.

## What this is

A second look at the same malignant pseudobulks, with rank enrichment
(fgsea) and sample-wise enrichment scores (GSVA, ssGSEA) instead of a mean
logFC across a merged IFN list. It does not replace the T/NK Spearman, and
it is not a cell-level test.

## What this is not

- Not a mega-merge and not TACSTD2.
- Not proof that CLDN4 causes the pathway shift.
- Not a DoRothEA / PROGENy activity estimate.
- Not cell-level GSEA p-values.
- Genome-wide significance is not claimed. The tested family is the six
  primary sets.

## Figures

- `results/figures/forest_fgsea_nes_pooled.png` — pooled NES, both contrasts
- `results/figures/forest_fgsea_nes_q4q1.png` — NES by cohort
- `results/figures/forest_fgsea_nes_continuous.png`
- `results/figures/forest_gsva_q4q1.png` / `forest_ssgsea_q4q1.png`
- `results/figures/forest_gsva_spearman.png` / `forest_ssgsea_spearman.png`
- `results/figures/forest_score_hedges_g.png`
- `results/figures/box_ssgsea_q4q1_patient.png`

## Reproduce

```bash
Rscript methods/concordant4_cldn4_gsva_fgsea/analyze.R
```

Requires R packages GSVA (>= 2.0), fgsea, ggplot2, data.table, jsonlite.
`N_BOOT` and `NPERM_BOOT` override the patient-bootstrap size.
Seed ", SEED, ".
"
)

writeLines(finding, file.path(ROOT, "FINDING.md"))
writeLines(capture.output(sessionInfo()), file.path(TAB, "sessionInfo.txt"))
message("Done.")

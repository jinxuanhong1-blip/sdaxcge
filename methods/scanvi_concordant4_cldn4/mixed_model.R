#!/usr/bin/env Rscript
# Patient-unit mixed models for continuous CLDN4 vs T/NK fraction.
# Grouping factor unit_id is the locked donor / sample / patient (n = 65).
# The binomial GLMM uses an observation-level patient random intercept so cells
# inside a unit are not treated as independent trials.

suppressPackageStartupMessages(library(lme4))

args <- commandArgs(trailingOnly = TRUE)
infile <- args[[1]]
outfile <- args[[2]]
df <- read.delim(infile, stringsAsFactors = FALSE)
df$dataset <- factor(df$dataset)
df$unit_id <- factor(df$unit_id)
df$n_other <- df$n_cells - df$n_tnk
if (any(df$n_other < 0)) stop("n_tnk exceeds n_cells")

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
rows <- list()

add_coef <- function(model_name, term, estimate, se, stat, p, n_units, note) {
  rows[[length(rows) + 1]] <<- data.frame(
    model = model_name, term = term, estimate = estimate, se = se, stat = stat,
    p = p, n_units = n_units, note = note, stringsAsFactors = FALSE
  )
}

fit_glmer <- function(score, label) {
  d <- df[is.finite(df[[score]]) & is.finite(df$frac_tnk), , drop = FALSE]
  d$x <- as.numeric(scale(d[[score]]))
  note <- ""
  m <- tryCatch(
    glmer(cbind(n_tnk, n_other) ~ x + dataset + (1 | unit_id), data = d, family = binomial, control = ctrl),
    error = function(e) e
  )
  if (inherits(m, "error")) {
    add_coef(label, "x", NA, NA, NA, NA, nrow(d), paste("glmer failed:", conditionMessage(m)))
    return()
  }
  singular <- isSingular(m, tol = 1e-4)
  vc <- as.data.frame(VarCorr(m))
  re_sd <- if (nrow(vc)) vc$sdcor[1] else NA_real_
  co <- summary(m)$coefficients
  # term name is "x"
  est <- co["x", "Estimate"]
  se <- co["x", "Std. Error"]
  z <- co["x", "z value"]
  p <- co["x", "Pr(>|z|)"]
  note <- sprintf("patient RE sd=%.4f singular=%s family=binomial Wald", re_sd, singular)
  add_coef(label, "cldn4_z", est, se, z, p, nrow(d), note)
  m0 <- tryCatch(
    glmer(cbind(n_tnk, n_other) ~ dataset + (1 | unit_id), data = d, family = binomial, control = ctrl),
    error = function(e) e
  )
  if (!inherits(m0, "error")) {
    ao <- anova(m0, m)
    # second row is the larger model
    lrt_p <- ao$`Pr(>Chisq)`[2]
    lrt_chi <- ao$Chisq[2]
    add_coef(paste0(label, "_LRT"), "cldn4_z", est, se, lrt_chi, lrt_p, nrow(d),
             paste(note, "LRT vs dataset + (1|unit)"))
  } else {
    add_coef(paste0(label, "_LRT"), "cldn4_z", NA, NA, NA, NA, nrow(d), paste("LRT failed:", conditionMessage(m0)))
  }
}

fit_lmer <- function(score, label) {
  d <- df[is.finite(df[[score]]) & is.finite(df$frac_tnk), , drop = FALSE]
  d$x <- as.numeric(scale(d[[score]]))
  m <- tryCatch(lmer(frac_tnk ~ x + (1 | dataset), data = d, REML = FALSE), error = function(e) e)
  if (inherits(m, "error")) {
    add_coef(label, "cldn4_z", NA, NA, NA, NA, nrow(d), paste("lmer failed:", conditionMessage(m)))
    return()
  }
  co <- summary(m)$coefficients
  est <- co["x", "Estimate"]
  se <- co["x", "Std. Error"]
  tval <- co["x", "t value"]
  df_res <- nrow(d) - 2
  p <- 2 * pt(-abs(tval), df = df_res)
  singular <- isSingular(m, tol = 1e-4)
  vc <- as.data.frame(VarCorr(m))
  re_sd <- if (nrow(vc)) vc$sdcor[1] else NA_real_
  add_coef(label, "cldn4_z", est, se, tval, p, nrow(d),
           sprintf("dataset RE sd=%.4f singular=%s gaussian t df=%d", re_sd, singular, df_res))
}

fit_lm <- function(score, label) {
  d <- df[is.finite(df[[score]]) & is.finite(df$frac_tnk), , drop = FALSE]
  d$x <- as.numeric(scale(d[[score]]))
  m <- lm(frac_tnk ~ x + dataset, data = d)
  co <- summary(m)$coefficients
  est <- co["x", "Estimate"]
  se <- co["x", "Std. Error"]
  tval <- co["x", "t value"]
  p <- co["x", "Pr(>|t|)"]
  add_coef(label, "cldn4_z", est, se, tval, p, nrow(d), "OLS equal-unit weight, dataset fixed effect")
  m2 <- lm(frac_tnk ~ x, data = d)
  co2 <- summary(m2)$coefficients
  add_coef(paste0(label, "_marginal"), "cldn4_z", co2["x", "Estimate"], co2["x", "Std. Error"],
           co2["x", "t value"], co2["x", "Pr(>|t|)"], nrow(d), "OLS equal-unit weight, no dataset term")
}

for (score in c("mal_CLDN4_pct", "mal_CLDN4_mean_log1p", "latent_CLDN4_seed", "latent_CLDN4")) {
  if (!score %in% names(df)) next
  fit_glmer(score, paste0("glmer_patientRE_", score))
  fit_lmer(score, paste0("lmer_datasetRE_", score))
  fit_lm(score, paste0("lm_datasetFE_", score))
}

out <- do.call(rbind, rows)
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
print(out)

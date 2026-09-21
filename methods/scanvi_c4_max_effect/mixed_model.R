#!/usr/bin/env Rscript
# Patient-level binomial GLMMs. One row per locked unit.
# Likelihood-ratio statistic is the patient-model LR effect:
# cbind(success, fail) ~ cl_z + dataset + (1 | unit_id)
# versus the same model without cl_z.

suppressPackageStartupMessages(library(lme4))

args <- commandArgs(trailingOnly = TRUE)
infile <- args[[1]]
outfile <- args[[2]]
df <- read.delim(infile, stringsAsFactors = FALSE)
df$dataset <- factor(df$dataset)
df$unit_id <- factor(df$unit_id)

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), calc.derivs = FALSE)

score_cols <- strsplit(args[[3]], ",", fixed = TRUE)[[1]]
outcomes <- list(
  all_cells = c("n_tnk", "n_fail_all"),
  compartment = c("n_tnk", "n_fail_comp")
)

z_global <- function(x) as.numeric(scale(x))
z_within <- function(x, dataset) {
  out <- rep(NA_real_, length(x))
  for (ds in unique(dataset)) {
    m <- dataset == ds & is.finite(x)
    s <- sd(x[m])
    if (!is.finite(s) || s == 0) {
      out[m] <- 0
    } else {
      out[m] <- (x[m] - mean(x[m])) / s
    }
  }
  out
}
z_rank <- function(x, dataset) {
  ranked <- rep(NA_real_, length(x))
  for (ds in unique(dataset)) {
    m <- dataset == ds & is.finite(x)
    ranked[m] <- rank(x[m], ties.method = "average")
  }
  z_within(ranked, dataset)
}

rows <- list()
add <- function(...) {
  rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
}

for (score in score_cols) {
  for (outcome in names(outcomes)) {
    succ <- outcomes[[outcome]][[1]]
    fail <- outcomes[[outcome]][[2]]
    for (coding in c("global_z", "within_z", "within_rank")) {
      d <- df[is.finite(df[[score]]) & is.finite(df[[succ]]) & is.finite(df[[fail]]) & df[[fail]] >= 0, , drop = FALSE]
      d$x <- switch(
        coding,
        global_z = z_global(d[[score]]),
        within_z = z_within(d[[score]], d$dataset),
        within_rank = z_rank(d[[score]], d$dataset)
      )
      d <- d[is.finite(d$x), , drop = FALSE]
      label <- paste(score, outcome, coding, sep = "|")
      m <- tryCatch(
        glmer(
          cbind(d[[succ]], d[[fail]]) ~ x + dataset + (1 | unit_id),
          data = d, family = binomial, control = ctrl
        ),
        error = function(e) e
      )
      if (inherits(m, "error")) {
        add(spec = label, score = score, outcome = outcome, coding = coding,
            estimate = NA, se = NA, wald_z = NA, wald_p = NA, lr_chi2 = NA, lr_p = NA,
            n_units = nrow(d), re_sd = NA, singular = NA, note = conditionMessage(m))
        next
      }
      singular <- isSingular(m, tol = 1e-4)
      vc <- as.data.frame(VarCorr(m))
      re_sd <- if (nrow(vc)) vc$sdcor[[1]] else NA_real_
      co <- summary(m)$coefficients
      est <- co["x", "Estimate"]
      se <- co["x", "Std. Error"]
      z <- co["x", "z value"]
      p <- co["x", "Pr(>|z|)"]
      m0 <- tryCatch(
        glmer(
          cbind(d[[succ]], d[[fail]]) ~ dataset + (1 | unit_id),
          data = d, family = binomial, control = ctrl
        ),
        error = function(e) e
      )
      if (inherits(m0, "error")) {
        add(spec = label, score = score, outcome = outcome, coding = coding,
            estimate = est, se = se, wald_z = z, wald_p = p, lr_chi2 = NA, lr_p = NA,
            n_units = nrow(d), re_sd = re_sd, singular = singular, note = conditionMessage(m0))
        next
      }
      ao <- anova(m0, m)
      add(spec = label, score = score, outcome = outcome, coding = coding,
          estimate = est, se = se, wald_z = z, wald_p = p,
          lr_chi2 = ao$Chisq[[2]], lr_p = ao$`Pr(>Chisq)`[[2]],
          n_units = nrow(d), re_sd = re_sd, singular = singular, note = "ok")
      message(sprintf("%s  beta=%.3f  chi2=%.2f  p=%.3g", label, est, ao$Chisq[[2]], ao$`Pr(>Chisq)`[[2]]))
    }
  }
}

out <- do.call(rbind, rows)
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

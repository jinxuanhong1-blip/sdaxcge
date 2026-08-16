#!/usr/bin/env Rscript
# Template: median cut vs maximally selected log-rank (survival / maxstat).
# The naive p-value at the 'best' cut is not a valid type-I error.

suppressPackageStartupMessages({
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: cutpoint.R input.csv [time_col event_col marker_col min_frac n_perm]")
}
infile <- args[[1]]
time_col <- if (length(args) >= 2) args[[2]] else "time"
event_col <- if (length(args) >= 3) args[[3]] else "event"
marker_col <- if (length(args) >= 4) args[[4]] else "marker"
min_frac <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.25
n_perm <- if (length(args) >= 6) as.integer(args[[6]]) else 1000L

df <- read.csv(infile, stringsAsFactors = FALSE)
surv <- Surv(df[[time_col]], df[[event_col]])
x <- df[[marker_col]]

med <- median(x, na.rm = TRUE)
high_med <- as.integer(x > med)
fit_med <- survdiff(surv ~ high_med)
p_med <- pchisq(fit_med$chisq, 1, lower.tail = FALSE)
cat(sprintf("median cut: chi2=%.2f p=%.4f n_high=%d n_low=%d\n",
            fit_med$chisq, p_med, sum(high_med), sum(1 - high_med)))

if (requireNamespace("maxstat", quietly = TRUE)) {
  ms <- maxstat::maxstat.test(
    surv ~ x, data = data.frame(surv = surv, x = x),
    smethod = "LogRank", pmethod = "Lau92", minprop = min_frac, maxprop = 1 - min_frac
  )
  cat(sprintf(
    "maxstat cut=%.4g statistic=%.2f p_corrected=%.4f\n",
    unname(ms$estimate), unname(ms$statistic), ms$p.value
  ))
} else {
  cat("maxstat not installed; falling back to a permutation of a grid search\n")
  qs <- quantile(x, probs = seq(min_frac, 1 - min_frac, length.out = 20), na.rm = TRUE)
  best <- list(stat = -Inf, cut = NA_real_)
  for (cut in unique(qs)) {
    g <- as.integer(x > cut)
    if (min(sum(g), sum(1 - g)) < 2) next
    st <- survdiff(surv ~ g)$chisq
    if (st > best$stat) best <- list(stat = st, cut = cut)
  }
  null <- replicate(n_perm, {
    xp <- sample(x)
    mx <- -Inf
    for (cut in unique(qs)) {
      g <- as.integer(xp > cut)
      if (min(sum(g), sum(1 - g)) < 2) next
      mx <- max(mx, survdiff(surv ~ g)$chisq)
    }
    mx
  })
  p_perm <- (1 + sum(null >= best$stat)) / (1 + n_perm)
  cat(sprintf("grid-search cut=%.4g chi2=%.2f permutation p=%.4f\n",
              best$cut, best$stat, p_perm))
}

cat("Do not report the naive p-value as the result.\n")

#!/usr/bin/env Rscript
# Template: Cox PH for a continuous ICI biomarker (survival / survminer).
# Primary analysis is the continuous marker (HR per 1 SD). See playbook.md §§3–5.

suppressPackageStartupMessages({
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: cox_ph.R input.csv [time_col event_col marker_col min_epv]")
}
infile <- args[[1]]
time_col <- if (length(args) >= 2) args[[2]] else "time"
event_col <- if (length(args) >= 3) args[[3]] else "event"
marker_col <- if (length(args) >= 4) args[[4]] else "marker"
min_epv <- if (length(args) >= 5) as.numeric(args[[5]]) else 10

df <- read.csv(infile, stringsAsFactors = FALSE)
df <- df[complete.cases(df[, c(time_col, event_col, marker_col)]), ]
n_events <- sum(df[[event_col]] == 1)
epv <- n_events / 1
if (epv < min_epv) {
  stop(sprintf(
    "refuse: events per variable = %.1f (%d events) < %.1f",
    epv, n_events, min_epv
  ))
}

df$z <- as.numeric(scale(df[[marker_col]]))
surv <- Surv(df[[time_col]], df[[event_col]])
fit <- coxph(surv ~ z, data = df, ties = "efron")
print(summary(fit))
cat(sprintf("Harrell C (apparent) = %.3f\n", concordance(fit)$concordance))
cat(sprintf("events per variable = %.1f\n", epv))

zph <- cox.zph(fit)
print(zph)
if (any(zph$table[, "p"] < 0.05, na.rm = TRUE)) {
  cat("PH is questionable. Report RMST or a time-varying coefficient.\n")
}

if (requireNamespace("survminer", quietly = TRUE)) {
  print(survminer::ggcoxzph(zph))
}

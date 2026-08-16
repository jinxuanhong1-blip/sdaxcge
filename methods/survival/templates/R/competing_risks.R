#!/usr/bin/env Rscript
# Template: Aalen–Johansen / Fine–Gray vs 1−KM that censors the competing event.
# Requires cause: 0 = censored, 1 = event of interest, 2+ = competing.
# Most GEO ICI series do not deposit this field.

suppressPackageStartupMessages({
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: competing_risks.R input.csv horizon [time_col cause_col]")
}
infile <- args[[1]]
horizon <- as.numeric(args[[2]])
time_col <- if (length(args) >= 3) args[[3]] else "time"
cause_col <- if (length(args) >= 4) args[[4]] else "cause"

df <- read.csv(infile, stringsAsFactors = FALSE)
df <- df[complete.cases(df[, c(time_col, cause_col)]), ]
causes <- sort(unique(df[[cause_col]]))
if (all(causes %in% c(0, 1))) {
  stop("no competing event is coded. Do not run a competing-risks model just because the playbook mentions one.")
}

# 1-KM that incorrectly censors competing events
km <- survfit(Surv(df[[time_col]], as.integer(df[[cause_col]] == 1)) ~ 1)
km_at <- summary(km, times = horizon, extend = TRUE)
naive <- 1 - km_at$surv

# Aalen-Johansen via survfit(Surv(..., type = "mstate"))
# status must be a factor whose first level is censoring
st <- factor(df[[cause_col]], levels = sort(unique(c(0, df[[cause_col]]))))
aj <- survfit(Surv(df[[time_col]], st, type = "mstate") ~ 1)
aj_sum <- summary(aj, times = horizon, extend = TRUE)
# pstate columns follow the factor levels excluding the censoring level
cif <- NA_real_
if (!is.null(aj_sum$pstate)) {
  # column for cause 1
  cause_levels <- levels(st)[-1]
  idx <- match("1", cause_levels)
  if (!is.na(idx)) cif <- aj_sum$pstate[1, idx]
}

cat(sprintf("Aalen-Johansen CIF(cause=1, t=%g) = %.3f\n", horizon, cif))
cat(sprintf("1-KM treating competing events as censoring     = %.3f\n", naive))
cat(sprintf("absolute overestimate                           = %.3f\n", naive - cif))

if (requireNamespace("cmprsk", quietly = TRUE) && "marker" %in% names(df)) {
  cat("\nFine-Gray on a binary high/low marker (illustrative):\n")
  fg <- cmprsk::crr(
    ftime = df[[time_col]],
    fstatus = df[[cause_col]],
    cov1 = data.frame(high = as.integer(df$marker > median(df$marker, na.rm = TRUE))),
    failcode = 1,
    cencode = 0
  )
  print(summary(fg))
}

#!/usr/bin/env Rscript
# Template: time-dependent (cumulative/dynamic) AUC with IPCW (timeROC).
# Pass a RISK score. For a protective marker use the Cox linear predictor.

suppressPackageStartupMessages({
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: time_dependent_auc.R input.csv [time_col event_col marker_col times_csv]")
}
infile <- args[[1]]
time_col <- if (length(args) >= 2) args[[2]] else "time"
event_col <- if (length(args) >= 3) args[[3]] else "event"
marker_col <- if (length(args) >= 4) args[[4]] else "marker"
times <- if (length(args) >= 5) {
  as.numeric(strsplit(args[[5]], ",", fixed = TRUE)[[1]])
} else {
  c(3, 6, 12)
}

df <- read.csv(infile, stringsAsFactors = FALSE)
df <- df[complete.cases(df[, c(time_col, event_col, marker_col)]), ]
fit <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ df[[marker_col]], ties = "efron")
risk <- as.numeric(predict(fit, type = "lp"))
last_event <- max(df[[time_col]][df[[event_col]] == 1])
times <- times[times < last_event]
if (!length(times)) {
  stop(sprintf("no evaluation time is earlier than the last event (%.2f)", last_event))
}

if (!requireNamespace("timeROC", quietly = TRUE)) {
  stop("install.packages('timeROC') is required for this template")
}

tr <- timeROC::timeROC(
  T = df[[time_col]],
  delta = df[[event_col]],
  marker = risk,
  cause = 1,
  times = times,
  iid = TRUE
)
print(tr)
cat(sprintf("Cox log-HR = %.3f\n", unname(coef(fit))))
cat("Apparent AUC; nest inside CV or apply to a locked external cohort.\n")

#!/usr/bin/env Rscript
# Continuous TACSTD2 survival and treatment-interaction analysis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: bessede_survival.R clinical.tsv results.tsv")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("R package 'survival' is required")
}

dat <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = TRUE)
required <- c(
  "pfs_time", "pfs_event", "os_time", "os_event", "treatment",
  "TACSTD2", "histology", "CD274", "TLS"
)
missing <- setdiff(required, names(dat))
if (length(missing)) stop("missing columns: ", paste(missing, collapse = ", "))
dat <- dat[complete.cases(dat[, required]), required]
if (nrow(dat) < 20L) stop("fewer than 20 complete observations")

dat$TACSTD2_z <- as.numeric(scale(dat$TACSTD2))
dat$CD274_z <- as.numeric(scale(dat$CD274))
dat$TLS_z <- as.numeric(scale(dat$TLS))
dat$TACSTD2_median <- factor(
  ifelse(dat$TACSTD2 > median(dat$TACSTD2), "high", "low"),
  levels = c("low", "high")
)

extract <- function(fit, model_name, term_pattern) {
  table <- summary(fit)$coefficients
  rows <- grep(term_pattern, rownames(table))
  if (!length(rows)) return(NULL)
  data.frame(
    model = model_name,
    term = rownames(table)[rows],
    hazard_ratio = exp(table[rows, "coef"]),
    ci_low = exp(table[rows, "coef"] - 1.96 * table[rows, "se(coef)"]),
    ci_high = exp(table[rows, "coef"] + 1.96 * table[rows, "se(coef)"]),
    p_value = table[rows, "Pr(>|z|)"],
    n = fit$n,
    row.names = NULL
  )
}

results <- list()
for (endpoint in c("pfs", "os")) {
  outcome <- survival::Surv(
    dat[[paste0(endpoint, "_time")]],
    dat[[paste0(endpoint, "_event")]]
  )
  adjusted <- survival::coxph(
    outcome ~ TACSTD2_z + treatment + histology + CD274_z + TLS_z,
    data = dat, ties = "efron"
  )
  interaction <- survival::coxph(
    outcome ~ TACSTD2_z * treatment + histology + CD274_z + TLS_z,
    data = dat, ties = "efron"
  )
  median_sensitivity <- survival::coxph(
    outcome ~ TACSTD2_median + treatment + histology + CD274_z + TLS_z,
    data = dat, ties = "efron"
  )
  results[[length(results) + 1L]] <- extract(
    adjusted, paste0(endpoint, "_adjusted_continuous"), "^TACSTD2_z$"
  )
  results[[length(results) + 1L]] <- extract(
    interaction, paste0(endpoint, "_treatment_interaction"), "TACSTD2_z:treatment"
  )
  results[[length(results) + 1L]] <- extract(
    median_sensitivity, paste0(endpoint, "_median_sensitivity"), "^TACSTD2_median"
  )
}

write.table(
  do.call(rbind, results), args[[2]], sep = "\t", quote = FALSE, row.names = FALSE
)

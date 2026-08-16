#!/usr/bin/env Rscript
# Prespecified CLDN4-high vs ORR / OS in IMvigor210 (atezolizumab, urothelial).
#
# Primary exposure: median split of CLDN4_log2_deseq (high = >= median).
# Primary ORR: package binaryResponse CR/PR vs SD/PD (NE excluded).
# Primary OS: os_months + censOS (1 = death).
# One RNA sample per ANONPT_ID (keep first in file order if duplicated).

suppressPackageStartupMessages({
  library(survival)
  library(stats)
})

args <- commandArgs(trailingOnly = TRUE)
in_path <- if (length(args) >= 1) args[[1]] else "results/w200/B5_IMvigor210/data/sample_level.tsv"
out_dir <- if (length(args) >= 2) args[[2]] else "results/w200/B5_IMvigor210"

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

d0 <- read.delim(in_path, stringsAsFactors = FALSE, check.names = FALSE)
n_rna_samples <- nrow(d0)
n_patients_raw <- length(unique(d0$ANONPT_ID))

# One row per patient. File order is CountDataSet column order.
d <- d0[!duplicated(d0$ANONPT_ID), ]
n_dropped_dup <- n_rna_samples - nrow(d)

med <- median(d$CLDN4_log2_deseq)
d$CLDN4_z <- as.numeric(scale(d$CLDN4_log2_deseq))
d$CLDN4_high <- ifelse(d$CLDN4_log2_deseq >= med, "high", "low")
d$CLDN4_high <- factor(d$CLDN4_high, levels = c("low", "high"))

ter <- quantile(d$CLDN4_log2_deseq, probs = c(1 / 3, 2 / 3), names = FALSE)
d$CLDN4_tertile <- cut(
  d$CLDN4_log2_deseq,
  breaks = c(-Inf, ter[1], ter[2], Inf),
  labels = c("T1_low", "T2_mid", "T3_high"),
  right = FALSE
)
# cut with right=FALSE can leave the exact max in NA if == Inf; force
d$CLDN4_tertile[is.na(d$CLDN4_tertile)] <- "T3_high"

qcuts <- quantile(d$CLDN4_log2_deseq, probs = c(0.25, 0.75), names = FALSE)
d$CLDN4_quartile <- NA_character_
d$CLDN4_quartile[d$CLDN4_log2_deseq <= qcuts[1]] <- "Q1"
d$CLDN4_quartile[d$CLDN4_log2_deseq >= qcuts[2]] <- "Q4"

d$responder <- NA_integer_
d$responder[d$binaryResponse == "CR/PR"] <- 1L
d$responder[d$binaryResponse == "SD/PD"] <- 0L
d$orr_eval <- !is.na(d$responder)

fmt_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 1e-4) return(format(p, digits = 2, scientific = TRUE))
  sprintf("%.4f", p)
}
fmt_num <- function(x, d = 3) {
  if (length(x) == 0 || !is.finite(x)) return(NA_character_)
  sprintf(paste0("%.", d, "f"), x)
}

extract_glm <- function(fit, term) {
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) return(NULL)
  est <- s[term, "Estimate"]
  se <- s[term, "Std. Error"]
  p <- s[term, "Pr(>|z|)"]
  list(
    or = exp(est),
    or_lo = exp(est - 1.96 * se),
    or_hi = exp(est + 1.96 * se),
    p = p,
    n = nobs(fit)
  )
}

extract_cox <- function(fit, term) {
  s <- summary(fit)
  if (!term %in% rownames(s$coefficients)) return(NULL)
  list(
    hr = s$conf.int[term, "exp(coef)"],
    hr_lo = s$conf.int[term, "lower .95"],
    hr_hi = s$conf.int[term, "upper .95"],
    p = s$coefficients[term, "Pr(>|z|)"],
    n = fit$n,
    n_event = as.integer(fit$nevent)
  )
}

orr_2x2 <- function(df, group_col, high_level, low_level) {
  x <- df[df$orr_eval & df[[group_col]] %in% c(high_level, low_level), ]
  x$g <- factor(x[[group_col]], levels = c(low_level, high_level))
  tab <- table(g = x$g, resp = factor(x$responder, levels = c(0, 1)))
  # rows: low, high; cols: nonresp, resp
  n_high <- sum(x$g == high_level)
  n_low <- sum(x$g == low_level)
  ev_high <- sum(x$g == high_level & x$responder == 1)
  ev_low <- sum(x$g == low_level & x$responder == 1)
  ft <- fisher.test(tab)
  fit <- glm(responder ~ g, data = x, family = binomial())
  lg <- extract_glm(fit, "ghigh")
  if (is.null(lg)) lg <- extract_glm(fit, paste0("g", high_level))
  list(
    n_eval = nrow(x),
    n_high = n_high,
    n_low = n_low,
    n_resp_high = ev_high,
    n_resp_low = ev_low,
    orr_high = ev_high / n_high,
    orr_low = ev_low / n_low,
    fisher_or = unname(ft$estimate),
    fisher_or_lo = ft$conf.int[1],
    fisher_or_hi = ft$conf.int[2],
    fisher_p = ft$p.value,
    logit_or = lg$or,
    logit_or_lo = lg$or_lo,
    logit_or_hi = lg$or_hi,
    logit_p = lg$p,
    table = tab
  )
}

os_binary <- function(df, group_col, high_level, low_level) {
  x <- df[df[[group_col]] %in% c(high_level, low_level) & is.finite(df$os_months) & df$censOS %in% c(0, 1), ]
  x$g <- factor(x[[group_col]], levels = c(low_level, high_level))
  sdif <- survdiff(Surv(os_months, censOS) ~ g, data = x)
  p_lr <- 1 - pchisq(sdif$chisq, length(sdif$n) - 1)
  fit <- coxph(Surv(os_months, censOS) ~ g, data = x)
  cx <- extract_cox(fit, "ghigh")
  if (is.null(cx)) cx <- extract_cox(fit, paste0("g", high_level))
  fit_km <- survfit(Surv(os_months, censOS) ~ g, data = x)
  list(
    n = nrow(x),
    n_event = sum(x$censOS == 1),
    n_high = sum(x$g == high_level),
    n_low = sum(x$g == low_level),
    logrank_p = p_lr,
    hr = cx$hr,
    hr_lo = cx$hr_lo,
    hr_hi = cx$hr_hi,
    cox_p = cx$p,
    km = fit_km,
    data = x
  )
}

# ---- primary ORR ----
prim_orr <- orr_2x2(d, "CLDN4_high", "high", "low")

# continuous ORR
d_orr <- d[d$orr_eval, ]
fit_orr_z <- glm(responder ~ CLDN4_z, data = d_orr, family = binomial())
orr_z <- extract_glm(fit_orr_z, "CLDN4_z")
wt <- wilcox.test(CLDN4_log2_deseq ~ factor(responder), data = d_orr)

# tertile / quartile ORR
orr_ter <- orr_2x2(d, "CLDN4_tertile", "T3_high", "T1_low")
orr_q <- orr_2x2(d, "CLDN4_quartile", "Q4", "Q1")

# CPM-median sensitivity
med_cpm <- median(d$CLDN4_log2_cpm)
d$CLDN4_high_cpm <- factor(ifelse(d$CLDN4_log2_cpm >= med_cpm, "high", "low"), levels = c("low", "high"))
orr_cpm <- orr_2x2(d, "CLDN4_high_cpm", "high", "low")

# ---- primary OS ----
prim_os <- os_binary(d, "CLDN4_high", "high", "low")
fit_os_z <- coxph(Surv(os_months, censOS) ~ CLDN4_z, data = d)
os_z <- extract_cox(fit_os_z, "CLDN4_z")
os_ter <- os_binary(d, "CLDN4_tertile", "T3_high", "T1_low")
os_q <- os_binary(d, "CLDN4_quartile", "Q4", "Q1")
os_cpm <- os_binary(d, "CLDN4_high_cpm", "high", "low")

# ---- multivariable (complete case) ----
# Covariates chosen a priori from trial-relevant clinical fields present in cds.
# TMB log1p; IC level; ECOG; sex; liver mets; prior platinum.
d$logTMB <- log1p(d$FMOne_TMB_per_MB)
d$IC_Level_f <- factor(d$IC_Level, levels = c("IC0", "IC1", "IC2+"))
d$Sex_f <- factor(d$Sex, levels = c("M", "F"))
d$ECOG_f <- factor(d$Baseline_ECOG)
d$Liver <- factor(ifelse(is.na(d$Met_Disease_Status), NA, ifelse(d$Met_Disease_Status == "Liver", "Liver", "NoLiver")), levels = c("NoLiver", "Liver"))
d$Plat_f <- factor(d$Received_platinum, levels = c("N", "Y"))

mv_vars <- c("CLDN4_high", "IC_Level_f", "ECOG_f", "Sex_f", "Liver", "Plat_f", "logTMB")
d_mv <- d[complete.cases(d[, c(mv_vars, "responder", "orr_eval", "os_months", "censOS")]), ]

fit_orr_mv <- glm(responder ~ CLDN4_high + IC_Level_f + ECOG_f + Sex_f + Liver + Plat_f + logTMB,
                  data = d_mv[d_mv$orr_eval, ], family = binomial())
orr_mv <- extract_glm(fit_orr_mv, "CLDN4_highhigh")

fit_os_mv <- coxph(Surv(os_months, censOS) ~ CLDN4_high + IC_Level_f + ECOG_f + Sex_f + Liver + Plat_f + logTMB,
                   data = d_mv)
os_mv <- extract_cox(fit_os_mv, "CLDN4_highhigh")

# continuous MV
fit_orr_mv_z <- glm(responder ~ CLDN4_z + IC_Level_f + ECOG_f + Sex_f + Liver + Plat_f + logTMB,
                    data = d_mv[d_mv$orr_eval, ], family = binomial())
orr_mv_z <- extract_glm(fit_orr_mv_z, "CLDN4_z")
fit_os_mv_z <- coxph(Surv(os_months, censOS) ~ CLDN4_z + IC_Level_f + ECOG_f + Sex_f + Liver + Plat_f + logTMB,
                     data = d_mv)
os_mv_z <- extract_cox(fit_os_mv_z, "CLDN4_z")

# ---- subgroups (unadjusted, primary median) ----
subgroup_defs <- list(
  all = function(x) rep(TRUE, nrow(x)),
  IC0 = function(x) x$IC_Level == "IC0",
  IC1 = function(x) x$IC_Level == "IC1",
  `IC2+` = function(x) x$IC_Level == "IC2+",
  platinum_Y = function(x) x$Received_platinum == "Y",
  platinum_N = function(x) x$Received_platinum == "N",
  liver = function(x) x$Met_Disease_Status == "Liver",
  no_liver = function(x) !is.na(x$Met_Disease_Status) & x$Met_Disease_Status != "Liver",
  TCGA_I = function(x) x$TCGA_Subtype == "I",
  TCGA_II = function(x) x$TCGA_Subtype == "II",
  TCGA_III = function(x) x$TCGA_Subtype == "III",
  TCGA_IV = function(x) x$TCGA_Subtype == "IV",
  inflamed = function(x) x$Immune_phenotype == "inflamed",
  excluded = function(x) x$Immune_phenotype == "excluded",
  desert = function(x) x$Immune_phenotype == "desert"
)

sub_rows <- list()
for (nm in names(subgroup_defs)) {
  keep <- subgroup_defs[[nm]](d)
  keep[is.na(keep)] <- FALSE
  xs <- d[keep, ]
  if (nrow(xs) < 20) next
  o <- try(orr_2x2(xs, "CLDN4_high", "high", "low"), silent = TRUE)
  s <- try(os_binary(xs, "CLDN4_high", "high", "low"), silent = TRUE)
  if (inherits(o, "try-error") || inherits(s, "try-error")) next
  sub_rows[[nm]] <- data.frame(
    subgroup = nm,
    n_patients = nrow(xs),
    n_orr = o$n_eval,
    n_resp_high = o$n_resp_high,
    n_high_orr = o$n_high,
    n_resp_low = o$n_resp_low,
    n_low_orr = o$n_low,
    orr_high = o$orr_high,
    orr_low = o$orr_low,
    fisher_OR = o$fisher_or,
    fisher_OR_lo = o$fisher_or_lo,
    fisher_OR_hi = o$fisher_or_hi,
    fisher_p = o$fisher_p,
    n_os = s$n,
    n_event = s$n_event,
    HR = s$hr,
    HR_lo = s$hr_lo,
    HR_hi = s$hr_hi,
    cox_p = s$cox_p,
    logrank_p = s$logrank_p,
    stringsAsFactors = FALSE
  )
}
sub_tab <- do.call(rbind, sub_rows)

# ---- 4-level response vs CLDN4 ----
d4 <- d[d$Best_Confirmed_Overall_Response %in% c("CR", "PR", "SD", "PD"), ]
resp_tab <- table(
  CLDN4 = d4$CLDN4_high,
  RECIST = factor(d4$Best_Confirmed_Overall_Response, levels = c("CR", "PR", "SD", "PD"))
)
resp_fisher <- fisher.test(resp_tab, simulate.p.value = TRUE, B = 20000)

# ---- write tables ----
headline <- data.frame(
  item = c(
    "primary_ORR_fisher_OR_high_vs_low",
    "primary_ORR_fisher_OR_95CI",
    "primary_ORR_fisher_p",
    "primary_ORR_n_eval",
    "primary_ORR_n_high_resp_over_n_high",
    "primary_ORR_n_low_resp_over_n_low",
    "primary_ORR_rate_high",
    "primary_ORR_rate_low",
    "primary_OS_HR_high_vs_low",
    "primary_OS_HR_95CI",
    "primary_OS_cox_p",
    "primary_OS_logrank_p",
    "primary_OS_n",
    "primary_OS_n_event",
    "CLDN4_median_log2_deseq",
    "n_rna_samples_in_cds",
    "n_patients_after_dedup",
    "n_duplicate_samples_dropped",
    "n_NE_excluded_from_ORR"
  ),
  value = c(
    fmt_num(prim_orr$fisher_or),
    sprintf("%.3f-%.3f", prim_orr$fisher_or_lo, prim_orr$fisher_or_hi),
    fmt_p(prim_orr$fisher_p),
    as.character(prim_orr$n_eval),
    sprintf("%d/%d", prim_orr$n_resp_high, prim_orr$n_high),
    sprintf("%d/%d", prim_orr$n_resp_low, prim_orr$n_low),
    fmt_num(prim_orr$orr_high, 3),
    fmt_num(prim_orr$orr_low, 3),
    fmt_num(prim_os$hr),
    sprintf("%.3f-%.3f", prim_os$hr_lo, prim_os$hr_hi),
    fmt_p(prim_os$cox_p),
    fmt_p(prim_os$logrank_p),
    as.character(prim_os$n),
    as.character(prim_os$n_event),
    fmt_num(med, 4),
    as.character(n_rna_samples),
    as.character(nrow(d)),
    as.character(n_dropped_dup),
    as.character(sum(d$binaryResponse == "" | is.na(d$binaryResponse) | d$Best_Confirmed_Overall_Response == "NE"))
  ),
  stringsAsFactors = FALSE
)
write.table(headline, file.path(tab_dir, "headline_OR_n_p.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

orr_rows <- rbind(
  data.frame(analysis = "primary_median_DESeq", contrast = "high_vs_low",
             n = prim_orr$n_eval, n_high = prim_orr$n_high, n_low = prim_orr$n_low,
             resp_high = prim_orr$n_resp_high, resp_low = prim_orr$n_resp_low,
             orr_high = prim_orr$orr_high, orr_low = prim_orr$orr_low,
             fisher_OR = prim_orr$fisher_or, fisher_OR_lo = prim_orr$fisher_or_lo,
             fisher_OR_hi = prim_orr$fisher_or_hi, fisher_p = prim_orr$fisher_p,
             logit_OR = prim_orr$logit_or, logit_OR_lo = prim_orr$logit_or_lo,
             logit_OR_hi = prim_orr$logit_or_hi, logit_p = prim_orr$logit_p),
  data.frame(analysis = "sensitivity_tertile_T3_vs_T1", contrast = "T3_vs_T1",
             n = orr_ter$n_eval, n_high = orr_ter$n_high, n_low = orr_ter$n_low,
             resp_high = orr_ter$n_resp_high, resp_low = orr_ter$n_resp_low,
             orr_high = orr_ter$orr_high, orr_low = orr_ter$orr_low,
             fisher_OR = orr_ter$fisher_or, fisher_OR_lo = orr_ter$fisher_or_lo,
             fisher_OR_hi = orr_ter$fisher_or_hi, fisher_p = orr_ter$fisher_p,
             logit_OR = orr_ter$logit_or, logit_OR_lo = orr_ter$logit_or_lo,
             logit_OR_hi = orr_ter$logit_or_hi, logit_p = orr_ter$logit_p),
  data.frame(analysis = "sensitivity_Q4_vs_Q1", contrast = "Q4_vs_Q1",
             n = orr_q$n_eval, n_high = orr_q$n_high, n_low = orr_q$n_low,
             resp_high = orr_q$n_resp_high, resp_low = orr_q$n_resp_low,
             orr_high = orr_q$orr_high, orr_low = orr_q$orr_low,
             fisher_OR = orr_q$fisher_or, fisher_OR_lo = orr_q$fisher_or_lo,
             fisher_OR_hi = orr_q$fisher_or_hi, fisher_p = orr_q$fisher_p,
             logit_OR = orr_q$logit_or, logit_OR_lo = orr_q$logit_or_lo,
             logit_OR_hi = orr_q$logit_or_hi, logit_p = orr_q$logit_p),
  data.frame(analysis = "sensitivity_median_log2CPM", contrast = "high_vs_low",
             n = orr_cpm$n_eval, n_high = orr_cpm$n_high, n_low = orr_cpm$n_low,
             resp_high = orr_cpm$n_resp_high, resp_low = orr_cpm$n_resp_low,
             orr_high = orr_cpm$orr_high, orr_low = orr_cpm$orr_low,
             fisher_OR = orr_cpm$fisher_or, fisher_OR_lo = orr_cpm$fisher_or_lo,
             fisher_OR_hi = orr_cpm$fisher_or_hi, fisher_p = orr_cpm$fisher_p,
             logit_OR = orr_cpm$logit_or, logit_OR_lo = orr_cpm$logit_or_lo,
             logit_OR_hi = orr_cpm$logit_or_hi, logit_p = orr_cpm$logit_p)
)
write.table(orr_rows, file.path(tab_dir, "orr_binary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cont_orr <- data.frame(
  analysis = c("continuous_per_SD_unadjusted", "continuous_per_SD_multivariable", "wilcoxon_CLDN4_by_binaryResponse"),
  n = c(orr_z$n, orr_mv_z$n, nrow(d_orr)),
  OR_or_stat = c(orr_z$or, orr_mv_z$or, unname(wt$statistic)),
  lo = c(orr_z$or_lo, orr_mv_z$or_lo, NA),
  hi = c(orr_z$or_hi, orr_mv_z$or_hi, NA),
  p = c(orr_z$p, orr_mv_z$p, wt$p.value),
  stringsAsFactors = FALSE
)
write.table(cont_orr, file.path(tab_dir, "orr_continuous.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

os_rows <- rbind(
  data.frame(analysis = "primary_median_DESeq", contrast = "high_vs_low",
             n = prim_os$n, n_event = prim_os$n_event, n_high = prim_os$n_high, n_low = prim_os$n_low,
             HR = prim_os$hr, HR_lo = prim_os$hr_lo, HR_hi = prim_os$hr_hi,
             cox_p = prim_os$cox_p, logrank_p = prim_os$logrank_p),
  data.frame(analysis = "sensitivity_tertile_T3_vs_T1", contrast = "T3_vs_T1",
             n = os_ter$n, n_event = os_ter$n_event, n_high = os_ter$n_high, n_low = os_ter$n_low,
             HR = os_ter$hr, HR_lo = os_ter$hr_lo, HR_hi = os_ter$hr_hi,
             cox_p = os_ter$cox_p, logrank_p = os_ter$logrank_p),
  data.frame(analysis = "sensitivity_Q4_vs_Q1", contrast = "Q4_vs_Q1",
             n = os_q$n, n_event = os_q$n_event, n_high = os_q$n_high, n_low = os_q$n_low,
             HR = os_q$hr, HR_lo = os_q$hr_lo, HR_hi = os_q$hr_hi,
             cox_p = os_q$cox_p, logrank_p = os_q$logrank_p),
  data.frame(analysis = "sensitivity_median_log2CPM", contrast = "high_vs_low",
             n = os_cpm$n, n_event = os_cpm$n_event, n_high = os_cpm$n_high, n_low = os_cpm$n_low,
             HR = os_cpm$hr, HR_lo = os_cpm$hr_lo, HR_hi = os_cpm$hr_hi,
             cox_p = os_cpm$cox_p, logrank_p = os_cpm$logrank_p),
  data.frame(analysis = "continuous_per_SD_unadjusted", contrast = "per_SD",
             n = os_z$n, n_event = os_z$n_event, n_high = NA, n_low = NA,
             HR = os_z$hr, HR_lo = os_z$hr_lo, HR_hi = os_z$hr_hi,
             cox_p = os_z$p, logrank_p = NA)
)
write.table(os_rows, file.path(tab_dir, "os_cox.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mv_tab <- data.frame(
  endpoint = c("ORR", "ORR", "OS", "OS"),
  exposure = c("CLDN4_high_vs_low", "CLDN4_per_SD", "CLDN4_high_vs_low", "CLDN4_per_SD"),
  n = c(orr_mv$n, orr_mv_z$n, os_mv$n, os_mv_z$n),
  n_event = c(NA, NA, os_mv$n_event, os_mv_z$n_event),
  effect = c("OR", "OR", "HR", "HR"),
  estimate = c(orr_mv$or, orr_mv_z$or, os_mv$hr, os_mv_z$hr),
  lo = c(orr_mv$or_lo, orr_mv_z$or_lo, os_mv$hr_lo, os_mv_z$hr_lo),
  hi = c(orr_mv$or_hi, orr_mv_z$or_hi, os_mv$hr_hi, os_mv_z$hr_hi),
  p = c(orr_mv$p, orr_mv_z$p, os_mv$p, os_mv_z$p),
  covariates = "IC_Level + ECOG + Sex + Liver + prior_platinum + log1p(TMB)",
  stringsAsFactors = FALSE
)
write.table(mv_tab, file.path(tab_dir, "multivariable.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sub_tab, file.path(tab_dir, "subgroups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

two_by_two <- as.data.frame.matrix(prim_orr$table)
two_by_two$CLDN4 <- rownames(two_by_two)
write.table(two_by_two, file.path(tab_dir, "orr_2x2_primary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

recist4 <- as.data.frame.matrix(resp_tab)
recist4$CLDN4 <- rownames(recist4)
write.table(recist4, file.path(tab_dir, "recist4_by_cldn4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# annotated analysis table (no full transcriptome)
d$orr_eval <- as.integer(d$orr_eval)
write.table(d, file.path(out_dir, "data", "analysis_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- figures ----
png(file.path(fig_dir, "km_os_median.png"), width = 1800, height = 1400, res = 200)
par(mar = c(5, 5, 3, 1))
plot(prim_os$km, col = c("#2c7bb6", "#d7191c"), lwd = 2, xlab = "Overall survival (months)",
     ylab = "Survival probability", mark.time = TRUE)
legend("topright",
       legend = c(
         sprintf("CLDN4-low  n=%d", prim_os$n_low),
         sprintf("CLDN4-high n=%d", prim_os$n_high),
         sprintf("HR %.2f (%.2f-%.2f)", prim_os$hr, prim_os$hr_lo, prim_os$hr_hi),
         sprintf("log-rank p=%s", fmt_p(prim_os$logrank_p))
       ),
       col = c("#2c7bb6", "#d7191c", NA, NA), lwd = c(2, 2, NA, NA), bty = "n")
dev.off()

png(file.path(fig_dir, "orr_bar_median.png"), width = 1400, height = 1400, res = 200)
par(mar = c(5, 5, 3, 1))
rates <- c(prim_orr$orr_low, prim_orr$orr_high)
bp <- barplot(rates, names.arg = c("CLDN4-low", "CLDN4-high"), ylim = c(0, max(rates) * 1.35),
              col = c("#2c7bb6", "#d7191c"), ylab = "ORR (CR/PR)", border = NA)
text(bp, rates, sprintf("%.1f%%\n%d/%d", 100 * rates,
                        c(prim_orr$n_resp_low, prim_orr$n_resp_high),
                        c(prim_orr$n_low, prim_orr$n_high)), pos = 3, xpd = TRUE)
title(sprintf("Fisher OR %.2f (%.2f-%.2f), p=%s, n=%d",
              prim_orr$fisher_or, prim_orr$fisher_or_lo, prim_orr$fisher_or_hi,
              fmt_p(prim_orr$fisher_p), prim_orr$n_eval))
dev.off()

png(file.path(fig_dir, "cldn4_by_binary_response.png"), width = 1400, height = 1400, res = 200)
par(mar = c(5, 5, 3, 1))
boxplot(CLDN4_log2_deseq ~ factor(binaryResponse, levels = c("SD/PD", "CR/PR")),
        data = d_orr, ylab = "CLDN4 log2(DESeq-norm + 1)", xlab = "binaryResponse",
        col = c("#bbbbbb", "#abdda4"), outline = TRUE)
stripchart(CLDN4_log2_deseq ~ factor(binaryResponse, levels = c("SD/PD", "CR/PR")),
           data = d_orr, method = "jitter", vertical = TRUE, add = TRUE, pch = 16, cex = 0.5, col = "#33333388")
title(sprintf("Wilcoxon p=%s, n=%d", fmt_p(wt$p.value), nrow(d_orr)))
dev.off()

png(file.path(fig_dir, "cldn4_by_recist4.png"), width = 1600, height = 1400, res = 200)
par(mar = c(5, 5, 3, 1))
boxplot(CLDN4_log2_deseq ~ factor(Best_Confirmed_Overall_Response, levels = c("CR", "PR", "SD", "PD")),
        data = d4, ylab = "CLDN4 log2(DESeq-norm + 1)", xlab = "Best confirmed overall response",
        col = c("#1a9641", "#a6d96a", "#fdae61", "#d7191c"))
title(sprintf("4-level Fisher (simulated) p=%s, n=%d", fmt_p(resp_fisher$p.value), nrow(d4)))
dev.off()

# subgroup forest for OS
if (nrow(sub_tab)) {
  png(file.path(fig_dir, "forest_os_subgroups.png"), width = 1800, height = 1600, res = 200)
  st <- sub_tab[is.finite(sub_tab$HR), ]
  st <- st[order(st$HR), ]
  par(mar = c(5, 12, 3, 3))
  ys <- seq_len(nrow(st))
  plot(st$HR, ys, pch = 16, xlim = c(min(0.4, min(st$HR_lo, na.rm = TRUE)), max(2.5, max(st$HR_hi, na.rm = TRUE))),
       yaxt = "n", xlab = "OS HR (CLDN4-high vs low)", ylab = "", log = "x")
  abline(v = 1, lty = 2, col = "grey50")
  segments(st$HR_lo, ys, st$HR_hi, ys)
  axis(2, at = ys, labels = sprintf("%s (n=%d)", st$subgroup, st$n_os), las = 1, cex.axis = 0.8)
  title("Exploratory OS subgroups (unadjusted)")
  dev.off()

  png(file.path(fig_dir, "forest_orr_subgroups.png"), width = 1800, height = 1600, res = 200)
  st <- sub_tab[is.finite(sub_tab$fisher_OR), ]
  st <- st[order(st$fisher_OR), ]
  par(mar = c(5, 12, 3, 3))
  ys <- seq_len(nrow(st))
  plot(st$fisher_OR, ys, pch = 16, xlim = c(min(0.15, min(st$fisher_OR_lo, na.rm = TRUE)), max(4, max(st$fisher_OR_hi, na.rm = TRUE))),
       yaxt = "n", xlab = "ORR Fisher OR (CLDN4-high vs low)", ylab = "", log = "x")
  abline(v = 1, lty = 2, col = "grey50")
  segments(st$fisher_OR_lo, ys, st$fisher_OR_hi, ys)
  axis(2, at = ys, labels = sprintf("%s (n=%d)", st$subgroup, st$n_orr), las = 1, cex.axis = 0.8)
  title("Exploratory ORR subgroups (unadjusted)")
  dev.off()
}

# machine-readable JSON-ish key results via dput
sink(file.path(log_dir, "primary_console.txt"))
cat("PRIMARY ORR\n")
print(prim_orr[setdiff(names(prim_orr), "table")])
cat("\n2x2\n"); print(prim_orr$table)
cat("\nPRIMARY OS\n")
print(prim_os[setdiff(names(prim_os), c("km", "data"))])
cat("\nCONTINUOUS ORR per SD\n"); print(orr_z)
cat("\nCONTINUOUS OS per SD\n"); print(os_z)
cat("\nMV ORR high\n"); print(orr_mv)
cat("\nMV OS high\n"); print(os_mv)
cat("\nWilcoxon\n"); print(wt)
cat("\n4-level simulated Fisher p=", resp_fisher$p.value, "\n")
sink()

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

cat("PRIMARY ORR: OR=", fmt_num(prim_orr$fisher_or),
    " n=", prim_orr$n_eval,
    " p=", fmt_p(prim_orr$fisher_p), "\n", sep = "")
cat("PRIMARY OS:  HR=", fmt_num(prim_os$hr),
    " n=", prim_os$n,
    " p=", fmt_p(prim_os$cox_p), "\n", sep = "")

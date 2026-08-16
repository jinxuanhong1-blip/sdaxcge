#!/usr/bin/env Rscript
# 09 · ICI overall-survival (OS) modelling. ICI 总生存（OS）建模。
# ---------------------------------------------------------------------------
# RUNS ONLY IF SURVIVAL LABELS EXIST. If the required time/event columns are not
# present in the input table, the script prints a notice and exits 0 WITHOUT
# producing any statistics. Nothing is imputed or invented.
# 仅当存在生存标签时运行。若输入表缺少所需的时间/事件列，脚本仅提示并以 0 退出，
# 不产生任何统计量，不做任何插补或编造。
#
# Typical inputs: an AOI/patient-level table that joins a spatial feature (e.g. a
# compartment DE signature score, a deconvolved cell-type fraction, or a
# TACSTD2/CLDN4 niche fraction) to clinical OS. Suits the GSE271689 design
# (GeoMx-WTA NSCLC, PD-1 immunotherapy, per-patient outcomes).
# 典型输入：将空间特征（如区室 DE 信号评分、去卷积细胞比例、TACSTD2/CLDN4 生态位
# 比例）与临床 OS 关联的 AOI/患者级表格。契合 GSE271689 设计。
#
# Kaplan-Meier + log-rank by group, univariable Cox, and multivariable Cox with
# covariates. Analysis is at the level of the 'id' column (aggregate AOIs to
# patient first to avoid pseudo-replication).
# KM + log-rank；单因素 Cox；含协变量的多因素 Cox。分析以 id 列为单位（先把 AOI
# 聚合到患者，避免伪重复）。
#
# Tools: survival, survminer, dplyr. No results are fabricated.
#
# Usage:
#   Rscript 09_ici_os_survival.R --table <clinical_features.csv> \
#       --time os_months --event os_event --group signature_high \
#       --covariates age,sex,stage --id patient --out methods/spatial/demo/out/survival
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--table", type = "character"),
  make_option("--time", type = "character", default = "os_months"),
  make_option("--event", type = "character", default = "os_event"),
  make_option("--group", type = "character", default = NA,
              help = "categorical grouping for KM/log-rank (optional)"),
  make_option("--covariates", type = "character", default = "",
              help = "comma-separated covariates for multivariable Cox"),
  make_option("--id", type = "character", default = NA,
              help = "unit of analysis; rows are aggregated to this id if given"),
  make_option("--out", type = "character", default = "methods/spatial/demo/out/survival")
)))

stopifnot(!is.null(opt$table))
df <- read.csv(opt$table, check.names = FALSE, stringsAsFactors = FALSE)

# ---- Label gate. 标签门控。 ----
required <- c(opt$time, opt$event)
missing <- setdiff(required, colnames(df))
if (length(missing) > 0) {
  message(sprintf(
    "[09_survival] survival labels not found (missing: %s). Skipping OS analysis; nothing fabricated.",
    paste(missing, collapse = ", ")))
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({ library(survival); library(survminer); library(dplyr) })
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

# Aggregate to the analysis unit if requested (mean of numeric features).
# 若指定分析单位，按其聚合（数值特征取均值），避免同一患者多 AOI 的伪重复。
if (!is.na(opt$id) && opt$id %in% colnames(df)) {
  num_cols <- names(df)[sapply(df, is.numeric)]
  df <- df %>% group_by(across(all_of(opt$id))) %>%
    summarise(across(all_of(setdiff(num_cols, opt$id)), ~ mean(.x, na.rm = TRUE)),
              across(where(is.character), ~ dplyr::first(.x)), .groups = "drop") %>%
    as.data.frame()
}

df <- df[!is.na(df[[opt$time]]) & !is.na(df[[opt$event]]), ]
surv_obj <- Surv(time = df[[opt$time]], event = df[[opt$event]])

# ---- Kaplan-Meier + log-rank. ----
if (!is.na(opt$group) && opt$group %in% colnames(df)) {
  fit <- survfit(surv_obj ~ df[[opt$group]])
  lr <- survdiff(surv_obj ~ df[[opt$group]])
  pval <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  message(sprintf("[09_survival] log-rank p = %.4g (group = %s)", pval, opt$group))
  ggsave(file.path(opt$out, "km_curve.png"),
         ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE)$plot,
         width = 7, height = 6, dpi = 150)
  writeLines(sprintf("logrank_p\t%.6g", pval), file.path(opt$out, "logrank.txt"))
}

# ---- Univariable Cox for each candidate feature (numeric, non-clinical). ----
covs <- if (nchar(opt$covariates)) strsplit(opt$covariates, ",")[[1]] else character(0)
reserved <- c(opt$time, opt$event, opt$group, opt$id, covs)
feature_cols <- names(df)[sapply(df, is.numeric)]
feature_cols <- setdiff(feature_cols, reserved)

uni <- lapply(feature_cols, function(f) {
  m <- tryCatch(coxph(surv_obj ~ df[[f]]), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)
  data.frame(feature = f,
             HR = s$coefficients[1, "exp(coef)"],
             lower95 = s$conf.int[1, "lower .95"],
             upper95 = s$conf.int[1, "upper .95"],
             pval = s$coefficients[1, "Pr(>|z|)"])
})
uni <- do.call(rbind, Filter(Negate(is.null), uni))
if (!is.null(uni)) {
  uni <- uni[order(uni$pval), ]
  write.csv(uni, file.path(opt$out, "cox_univariable.csv"), row.names = FALSE)
  message(sprintf("[09_survival] univariable Cox for %d features -> cox_univariable.csv",
                  nrow(uni)))
}

# ---- Multivariable Cox: leading feature + covariates. ----
if (!is.null(uni) && nrow(uni) > 0 && length(covs) > 0 && all(covs %in% colnames(df))) {
  lead <- uni$feature[1]
  fml <- as.formula(paste("surv_obj ~", paste(c(lead, covs), collapse = " + ")))
  m <- coxph(fml, data = df)
  capture.output(summary(m), file = file.path(opt$out, "cox_multivariable.txt"))
  message(sprintf("[09_survival] multivariable Cox (%s + %s) -> cox_multivariable.txt",
                  lead, paste(covs, collapse = "+")))
}
message("[09_survival] done.")

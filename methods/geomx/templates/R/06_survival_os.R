## 06_survival_os.R
## Patient-level Cox is PRIMARY. Clustered / frailty AOI-level models are
## sensitivity only. Continuous exposures, scaled per SD. Cut points off by
## default; if enabled, lock in the discovery cohort and apply unchanged.
##
## 主分析：患者层面 Cox。聚类/脆弱的 AOI 层面模型仅作敏感性分析。
## 暴露保持连续并按 SD 标准化。切点默认关闭；若启用，仅在发现集锁定后原样用于验证集。
##
## Inputs:  results/05_patient_exposures.tsv
##          results/04_targets_aoi.tsv          (for sensitivity models)
## Outputs: results/06_cox_patient.tsv
##          results/06_cox_sensitivity.tsv
##          results/06_validation.tsv
##          results/figures/06_km.pdf

source("templates/R/00_setup.R")
need(c("survival"), hard = TRUE)
need(c("coxme", "rms", "maxstat"), hard = FALSE)

log_step("06  overall survival")

exp_path <- out_path("05_patient_exposures.tsv")
if (!file.exists(exp_path)) stop("Run 05_patient_aggregation.R first.", call. = FALSE)
pdat <- utils::read.delim(exp_path, stringsAsFactors = FALSE)

time_col  <- cfg$survival$time
event_col <- cfg$survival$event
if (!all(c(time_col, event_col) %in% names(pdat))) {
  stop("Clinical columns ", time_col, " / ", event_col, " not in patient table. ",
       "Join them in the annotation (one value per patient).", call. = FALSE)
}

## Landmark (immortal time): drop patients with events before the landmark and
## subtract the landmark from time. Leave NULL unless tissue was on-treatment.
## Landmark（永生时间）：剔除 landmark 前发生事件的患者，并从时间中减去 landmark。
## 仅当组织在治疗中采集时才设置。
lm_m <- cfg$survival$landmark_months
if (!is.null(lm_m) && is.finite(lm_m) && lm_m > 0) {
  keep <- is.finite(pdat[[time_col]]) & pdat[[time_col]] >= lm_m
  pdat[[time_col]] <- pdat[[time_col]] - lm_m
  pdat <- pdat[keep, , drop = FALSE]
  log_step("landmark ", lm_m, " ", cfg$survival$time_unit, ": kept ", nrow(pdat), " rows")
}

covars <- intersect(unlist(cfg$survival$covariates), names(pdat))
strata <- intersect(unlist(cfg$survival$strata), names(pdat))
expos  <- intersect(unlist(cfg$survival$exposures), names(pdat))
if (!length(expos)) stop("No configured exposures found in the patient table.", call. = FALSE)

scale_if <- function(x) {
  if (isTRUE(cfg$survival$scale_exposures)) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(x)
    as.numeric(scale(x))
  } else x
}

fit_cox <- function(df, exposure, label) {
  df <- df[is.finite(df[[time_col]]) & is.finite(df[[event_col]]) & is.finite(df[[exposure]]), ]
  if (nrow(df) < 10 || sum(df[[event_col]]) < 5) {
    return(data.frame(label = label, exposure = exposure, n = nrow(df),
                      events = sum(df[[event_col]]), note = "too few events",
                      stringsAsFactors = FALSE))
  }
  df$x <- scale_if(df[[exposure]])
  rhs <- "x"
  if (length(covars)) rhs <- paste(c(rhs, covars), collapse = " + ")
  if (length(strata)) rhs <- paste(rhs, "+", paste(sprintf("strata(%s)", strata), collapse = " + "))
  form <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time_col, event_col, rhs))
  fit <- try(survival::coxph(form, data = df, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(label = label, exposure = exposure, note = as.character(fit),
                      stringsAsFactors = FALSE))
  }
  sm <- summary(fit)
  ph <- if (isTRUE(cfg$survival$ph_check)) {
    z <- try(survival::cox.zph(fit), silent = TRUE)
    if (inherits(z, "try-error")) NA_real_ else z$table["x", "p"]
  } else NA_real_
  ev <- sum(df[[event_col]])
  nvar <- length(stats::coef(fit))
  data.frame(
    label = label,
    exposure = exposure,
    n = sm$n,
    events = ev,
    events_per_variable = ev / max(nvar, 1),
    hr_per_sd = sm$conf.int["x", 1],
    hr_lo = sm$conf.int["x", 3],
    hr_hi = sm$conf.int["x", 4],
    p = sm$coefficients["x", "Pr(>|z|)"],
    c_index = sm$concordance[1],
    ph_p_exposure = ph,
    formula = rhs,
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## Primary: one row per patient, per gene, per exposure, discovery cohort
## 主分析：发现集中每位患者一行，逐基因逐暴露
## ---------------------------------------------------------------------------
disc <- cfg$study$discovery_cohort
val  <- unlist(cfg$study$validation_cohorts)
if ("cohort" %in% names(pdat) && !is.null(disc)) {
  p_disc <- pdat[pdat$cohort == disc, , drop = FALSE]
  if (!nrow(p_disc)) {
    warning("Discovery cohort '", disc, "' not found; using all patients for the primary fit. ",
            "This is not a validated analysis.")
    p_disc <- pdat
  }
} else {
  p_disc <- pdat
}

primary <- list()
for (g in unique(p_disc$gene)) {
  dg <- p_disc[p_disc$gene == g, ]
  for (ex in expos) {
    row <- fit_cox(dg, ex, label = paste("primary", g, sep = ":"))
    row$gene <- g
    row$cohort <- disc %||% "all"
    primary[[paste(g, ex, sep = "_")]] <- row
  }
}
prim_tab <- do.call(rbind, primary)
rownames(prim_tab) <- NULL
write_tsv(prim_tab, out_path("06_cox_patient.tsv"))

## Spline vs linear (functional form), discovery, tumor_score only.
## 函数形式：限制性立方样条 vs 线性，仅发现集 tumor_score。
if (isTRUE(cfg$survival$spline_check) && requireNamespace("rms", quietly = TRUE)) {
  sink(out_path("06_spline_check.txt"))
  for (g in unique(p_disc$gene)) {
    dg <- p_disc[p_disc$gene == g & is.finite(p_disc$tumor_score), ]
    if (nrow(dg) < 20) next
    dd <- rms::datadist(dg); options(datadist = "dd")
    f_lin <- rms::cph(stats::as.formula(sprintf("survival::Surv(%s, %s) ~ tumor_score",
                                                time_col, event_col)), data = dg)
    f_rcs <- try(rms::cph(stats::as.formula(sprintf("survival::Surv(%s, %s) ~ rms::rcs(tumor_score, 3)",
                                                    time_col, event_col)), data = dg), silent = TRUE)
    cat("==== ", g, " tumor_score linear vs rcs(3) ====\n", sep = "")
    if (!inherits(f_rcs, "try-error")) print(stats::anova(f_lin, f_rcs))
    else cat("rcs fit failed\n")
  }
  sink()
}

## ---------------------------------------------------------------------------
## Validation: apply the discovery model direction, do not refit cut points
## 验证：沿用发现集方向，不重新拟合切点
## ---------------------------------------------------------------------------
val_rows <- list()
if ("cohort" %in% names(pdat) && length(val)) {
  for (vc in val) {
    pv <- pdat[pdat$cohort == vc, , drop = FALSE]
    if (!nrow(pv)) next
    for (g in unique(pv$gene)) {
      dg <- pv[pv$gene == g, ]
      for (ex in expos) {
        row <- fit_cox(dg, ex, label = paste("validation", vc, g, sep = ":"))
        row$gene <- g
        row$cohort <- vc
        val_rows[[paste(vc, g, ex, sep = "_")]] <- row
      }
    }
  }
}
if (length(val_rows)) {
  val_tab <- do.call(rbind, val_rows)
  rownames(val_tab) <- NULL
  write_tsv(val_tab, out_path("06_validation.tsv"))
}

## ---------------------------------------------------------------------------
## Optional cut point, locked in discovery
## 可选切点，仅在发现集锁定
## ---------------------------------------------------------------------------
if (isTRUE(cfg$survival$cutpoint$use)) {
  method <- cfg$survival$cutpoint$method %||% "prespecified_median"
  cuts <- list()
  for (g in unique(p_disc$gene)) {
    dg <- p_disc[p_disc$gene == g, ]
    x <- dg$tumor_score
    cut <- switch(method,
      prespecified_median = stats::median(x, na.rm = TRUE),
      tertile = stats::quantile(x, 2/3, na.rm = TRUE),
      maxstat = {
        if (!requireNamespace("maxstat", quietly = TRUE)) {
          warning("maxstat not installed; falling back to median.")
          stats::median(x, na.rm = TRUE)
        } else {
          ms <- maxstat::maxstat.test(
            stats::as.formula(sprintf("survival::Surv(%s, %s) ~ tumor_score",
                                      time_col, event_col)),
            data = dg, smethod = "LogRank", pmethod = "Lau92")
          unname(ms$estimate)
        }
      },
      stats::median(x, na.rm = TRUE)
    )
    cuts[[g]] <- data.frame(gene = g, method = method, cut = as.numeric(cut),
                            locked_in = disc %||% "all",
                            stringsAsFactors = FALSE)
  }
  cut_tab <- do.call(rbind, cuts)
  write_tsv(cut_tab, out_path("06_cutpoints_locked.tsv"))
  ## Apply the locked cut to every cohort, including validation — do not re-estimate.
  ## 把锁定切点用于所有队列（含验证集）——不得重新估计。
  km_df <- merge(pdat, cut_tab[, c("gene", "cut")], by = "gene", all.x = TRUE)
  km_df$hi <- as.integer(km_df$tumor_score >= km_df$cut)
  pdf(fig_path("06_km.pdf"), width = 7, height = 6)
  for (g in unique(km_df$gene)) {
    for (coh in unique(km_df$cohort)) {
      d <- km_df[km_df$gene == g & km_df$cohort == coh & is.finite(km_df$hi), ]
      if (nrow(d) < 10) next
      sf <- survival::survfit(
        stats::as.formula(sprintf("survival::Surv(%s, %s) ~ hi", time_col, event_col)),
        data = d)
      plot(sf, col = c("#4C78A8", "#E45756"), lwd = 2,
           xlab = paste("time (", cfg$survival$time_unit, ")", sep = ""),
           ylab = "overall survival",
           main = sprintf("%s  %s  cut locked in %s", g, coh, disc %||% "all"))
      legend("topright", c("low", "high"), col = c("#4C78A8", "#E45756"), lwd = 2, bty = "n")
    }
  }
  dev.off()
}

## ---------------------------------------------------------------------------
## Sensitivity: AOI-level clustered Cox and coxme frailty
## 敏感性：AOI 层面聚类 Cox 与 coxme 脆弱模型
## These DUPLICATE the event; they do not create power (playbook §8.2).
## 这些模型复制事件，不带来效能（手册 §8.2）。
## ---------------------------------------------------------------------------
sens <- list()
aoi_path <- out_path("04_targets_aoi.tsv")
if (file.exists(aoi_path) && length(cfg$survival$sensitivity_models)) {
  aoi <- utils::read.delim(aoi_path, stringsAsFactors = FALSE)
  ## Join patient-level clinical onto AOIs.
  clin <- pdat[, intersect(c("patient_id", time_col, event_col, covars, "cohort"),
                           names(pdat)), drop = FALSE]
  clin <- clin[!duplicated(clin$patient_id), ]
  aoi <- merge(aoi, clin, by.x = "patient", by.y = "patient_id", all.x = TRUE)
  aoi <- aoi[aoi$segment == TUMOR_SEGMENT, ]
  for (g in unique(aoi$gene)) {
    dg <- aoi[aoi$gene == g, ]
    dg <- dg[is.finite(dg$expr) & is.finite(dg[[time_col]]) & is.finite(dg[[event_col]]), ]
    if (nrow(dg) < 20) next
    dg$x <- scale_if(dg$expr)
    form <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ x", time_col, event_col))
    if ("cluster_robust" %in% unlist(cfg$survival$sensitivity_models)) {
      fit <- try(survival::coxph(form, data = dg, cluster = dg$patient), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        sm <- summary(fit)
        sens[[paste(g, "cluster")]] <- data.frame(
          gene = g, model = "cluster_robust",
          hr_per_sd = sm$conf.int["x", 1],
          hr_lo = sm$conf.int["x", 3],
          hr_hi = sm$conf.int["x", 4],
          p = sm$coefficients["x", "Pr(>|z|)"],
          n_aoi = sm$n,
          note = "sensitivity only; event is duplicated across AOIs",
          stringsAsFactors = FALSE
        )
      }
    }
    if ("coxme_frailty" %in% unlist(cfg$survival$sensitivity_models) &&
        requireNamespace("coxme", quietly = TRUE)) {
      fit <- try(coxme::coxme(
        stats::as.formula(sprintf("survival::Surv(%s, %s) ~ x + (1 | patient)",
                                  time_col, event_col)),
        data = dg), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        cf <- stats::coef(fit)
        sens[[paste(g, "coxme")]] <- data.frame(
          gene = g, model = "coxme_frailty",
          hr_per_sd = exp(cf["x"]),
          hr_lo = NA_real_, hr_hi = NA_real_,
          p = NA_real_,
          n_aoi = nrow(dg),
          note = "sensitivity only; conditional HR given patient frailty",
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(sens)) {
  sens_tab <- do.call(rbind, sens)
  rownames(sens_tab) <- NULL
  write_tsv(sens_tab, out_path("06_cox_sensitivity.tsv"))
}

save_session_info("06")
log_step("06 done. Primary rows: ", nrow(prim_tab))

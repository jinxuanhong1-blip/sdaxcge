## 07_power_simulation.R
## Nested simulation: patients -> ROIs -> segments -> expression, then
## (a) segment contrast, (b) between-group contrast, (c) OS hazard ratio
## at a given event rate and exposure reliability.
##
## 嵌套模拟：患者 -> ROI -> 分区 -> 表达，然后给出
## (a) 分区对比、(b) 组间对比、(c) 给定事件率与暴露可靠性下的 OS 风险比 的经验效能。
##
## No real data required. Uses only base R + survival (optional).
## 不需要真实数据，仅依赖 base R 与可选的 survival。
##
## Outputs: results/07_power.tsv
##          results/figures/07_power.pdf

source("templates/R/00_setup.R")
need(c("survival"), hard = FALSE)

log_step("07  power simulation (design planning)")

pw <- cfg$power
set.seed(pw$seed %||% 1)

## Variance decomposition used in every draw:
## y_ij = mu + a_i + e_ij,  a_i ~ N(0, sigma_b^2), e_ij ~ N(0, sigma_w^2)
## ICC = sigma_b^2 / (sigma_b^2 + sigma_w^2). We fix total variance = 1.
## 每次抽样的方差分解。总方差固定为 1。

draw_y <- function(n, m, icc, delta = 0) {
  ## n patients, m AOIs each, patient random intercept, optional mean shift delta.
  ## n 例患者、每人 m 个 AOI、患者随机截距、可选均值平移 delta。
  sigma_b <- sqrt(icc)
  sigma_w <- sqrt(1 - icc)
  a <- stats::rnorm(n, 0, sigma_b)
  y <- a[rep(seq_len(n), each = m)] + stats::rnorm(n * m, 0, sigma_w) + delta
  list(y = y, patient = factor(rep(seq_len(n), each = m)))
}

## (a) Paired segment contrast: each patient has m tumor and m immune AOIs.
##     True tumor-immune difference = log2FC (on the already-log scale).
## (a) 配对分区对比：每位患者 m 个肿瘤 AOI 与 m 个免疫 AOI。
power_segment <- function(n, m, icc, lfc, nsim) {
  sig <- logical(nsim)
  for (i in seq_len(nsim)) {
    tum <- draw_y(n, m, icc, delta = lfc)
    imm <- draw_y(n, m, icc, delta = 0)
    ## Paired patient means (the honest test).
    pt <- tapply(tum$y, tum$patient, mean)
    pi <- tapply(imm$y, imm$patient, mean)
    p <- try(stats::t.test(pt, pi, paired = TRUE)$p.value, silent = TRUE)
    sig[i] <- is.numeric(p) && p < 0.05
  }
  mean(sig)
}

## (b) Between-group contrast within one segment (n/2 vs n/2 patients).
## (b) 同一分区内的组间对比（n/2 vs n/2）。
power_group <- function(n, m, icc, lfc, nsim) {
  n1 <- floor(n / 2); n2 <- n - n1
  sig <- logical(nsim)
  for (i in seq_len(nsim)) {
    g1 <- draw_y(n1, m, icc, delta = lfc)
    g2 <- draw_y(n2, m, icc, delta = 0)
    p1 <- tapply(g1$y, g1$patient, mean)
    p2 <- tapply(g2$y, g2$patient, mean)
    p <- try(stats::t.test(p1, p2)$p.value, silent = TRUE)
    sig[i] <- is.numeric(p) && p < 0.05
  }
  mean(sig)
}

## (c) OS: Cox on a patient-level exposure with reliability R, true HR.
##     Event times ~ Exp, admin censoring to hit the configured event rate.
## (c) OS：在可靠性为 R 的患者层面暴露上拟合 Cox，真实 HR 给定。
power_os <- function(n, m, icc, hr, event_rate, nsim) {
  R <- spearman_brown(icc, m)
  sig <- logical(nsim)
  if (!requireNamespace("survival", quietly = TRUE)) return(NA_real_)
  for (i in seq_len(nsim)) {
    z_true <- stats::rnorm(n)
    z_obs  <- sqrt(R) * z_true + sqrt(1 - R) * stats::rnorm(n)
    ## Exponential PH: lambda = lambda0 * hr^{z_true / sd}.
    beta <- log(hr)
    lam <- 0.05 * exp(beta * as.numeric(scale(z_true)))
    t_evt <- stats::rexp(n, rate = lam)
    ## Calibrate censoring so that P(event) ~ event_rate.
    t_cen <- stats::rexp(n, rate = 0.05 * (1 - event_rate) / max(event_rate, 0.05))
    time <- pmin(t_evt, t_cen)
    event <- as.integer(t_evt <= t_cen)
    fit <- try(survival::coxph(survival::Surv(time, event) ~ z_obs), silent = TRUE)
    if (inherits(fit, "try-error")) {
      sig[i] <- FALSE
    } else {
      p <- summary(fit)$coefficients[1, "Pr(>|z|)"]
      sig[i] <- is.finite(p) && p < 0.05
    }
  }
  mean(sig)
}

grid <- expand.grid(
  n_patients = unlist(pw$n_patients),
  n_roi = unlist(pw$n_roi_per_patient),
  icc = unlist(pw$icc_patient),
  lfc = unlist(pw$effect_log2fc),
  hr = unlist(pw$hazard_ratio),
  stringsAsFactors = FALSE
)
nsim <- as.integer(pw$n_sim %||% 200)
log_step("grid rows: ", nrow(grid), "  nsim: ", nsim)

out <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  out[[i]] <- data.frame(
    n_patients = g$n_patients,
    n_roi = g$n_roi,
    icc = g$icc,
    reliability = spearman_brown(g$icc, g$n_roi),
    n_eff = n_eff(g$n_patients, g$n_roi, g$icc),
    lfc = g$lfc,
    power_segment = power_segment(g$n_patients, g$n_roi, g$icc, g$lfc, nsim),
    power_group = power_group(g$n_patients, g$n_roi, g$icc, g$lfc, nsim),
    hr = g$hr,
    event_rate = pw$event_rate,
    power_os = power_os(g$n_patients, g$n_roi, g$icc, g$hr, pw$event_rate, nsim),
    nsim = nsim,
    stringsAsFactors = FALSE
  )
  if (i %% 10 == 0) log_step("sim ", i, " / ", nrow(grid))
}
tab <- do.call(rbind, out)
write_tsv(tab, out_path("07_power.tsv"))

pdf(fig_path("07_power.pdf"), width = 8, height = 6)
for (metric in c("power_segment", "power_group", "power_os")) {
  ## One panel per ICC, x = n_patients, lines = n_roi, using the first lfc/hr.
  lfc0 <- unlist(pw$effect_log2fc)[1]
  hr0  <- unlist(pw$hazard_ratio)[1]
  sub <- tab[tab$lfc == lfc0 & tab$hr == hr0, ]
  iccs <- sort(unique(sub$icc))
  op <- par(mfrow = c(1, length(iccs)), mar = c(4, 4, 3, 1))
  for (ic in iccs) {
    s <- sub[sub$icc == ic, ]
    plot(NA, xlim = range(s$n_patients), ylim = c(0, 1),
         xlab = "n patients", ylab = "empirical power",
         main = sprintf("%s  ICC=%.2f", metric, ic))
    abline(h = 0.8, lty = 2, col = "grey")
    for (m in sort(unique(s$n_roi))) {
      ss <- s[s$n_roi == m, ]
      ss <- ss[order(ss$n_patients), ]
      lines(ss$n_patients, ss[[metric]], type = "b", pch = 16)
    }
    legend("bottomright", legend = paste("m =", sort(unique(s$n_roi))),
           bty = "n", cex = 0.7)
  }
  par(op)
}
dev.off()

save_session_info("07")
log_step("07 done -> ", out_path("07_power.tsv"))

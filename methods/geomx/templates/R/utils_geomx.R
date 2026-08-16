## utils_geomx.R -- shared helpers for the GeoMx WTA playbook
## No Bioconductor dependency. Safe to source from 00_setup.R before optional packages load.
## 不依赖 Bioconductor，可在可选包加载前由 00_setup.R source。

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

ann_col <- function(cfg, role) {
  ## Map a role (patient, segment, ...) to the user's annotation column name.
  ## 把角色映射到用户注释表中的列名。
  val <- cfg$annotation[[role]]
  if (is.null(val) || !nzchar(val)) return(NA_character_)
  val
}

require_cols <- function(df, cols, context = "data") {
  cols <- cols[!is.na(cols) & nzchar(cols)]
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(sprintf("%s is missing required columns: %s", context, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

factor_segment <- function(x, levels = SEGMENT_LEVELS) {
  x <- as.character(x)
  unknown <- setdiff(unique(x[!is.na(x)]), levels)
  if (length(unknown)) {
    stop("Unknown segment levels: ", paste(unknown, collapse = ", "),
         ". Expected: ", paste(levels, collapse = ", "), call. = FALSE)
  }
  factor(x, levels = levels)
}

is_tumor <- function(segment, tumor = TUMOR_SEGMENT) {
  as.character(segment) == tumor
}

is_immune <- function(segment, immune = IMMUNE_SEGMENTS) {
  as.character(segment) %in% immune
}

## ---------------------------------------------------------------------------
## LOQ and detection / 定量下限与检出
## ---------------------------------------------------------------------------
loq_from_neg <- function(neg_geomean, neg_geosd, cutoff = 2, min_loq = 2) {
  pmax(min_loq, as.numeric(neg_geomean) * as.numeric(neg_geosd)^cutoff)
}

detection_rate <- function(counts, loq, na.rm = TRUE) {
  ## counts: genes x AOIs; loq: length = ncol(counts) or genes x AOIs
  if (is.null(dim(loq))) {
    above <- sweep(counts, 2, loq, FUN = ">")
  } else {
    above <- counts > loq
  }
  colMeans(above, na.rm = na.rm)
}

gene_detection_by_segment <- function(counts, loq, segment) {
  ## Returns a genes x segments matrix of detection rates.
  ## 返回 基因 x 分区 的检出率矩阵。
  segment <- as.character(segment)
  if (is.null(dim(loq))) loq <- matrix(loq, nrow = nrow(counts), ncol = ncol(counts), byrow = TRUE)
  above <- counts > loq
  segs <- unique(segment)
  out <- sapply(segs, function(s) {
    idx <- which(segment == s)
    if (!length(idx)) return(rep(NA_real_, nrow(counts)))
    rowMeans(above[, idx, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1, dimnames = list(rownames(counts), segs))
  out
}

keep_genes_per_segment <- function(det_by_seg, min_rate = 0.05, how = c("union", "intersection")) {
  how <- match.arg(how)
  ok <- det_by_seg >= min_rate
  if (how == "union") rownames(det_by_seg)[rowSums(ok, na.rm = TRUE) > 0]
  else rownames(det_by_seg)[rowSums(ok, na.rm = TRUE) == ncol(ok)]
}

## ---------------------------------------------------------------------------
## Variance components / ICC / reliability / 方差成分、组内相关、可靠性
## ---------------------------------------------------------------------------
icc_oneway <- function(y, group) {
  ## Unbiased one-way ICC from a random-intercept ANOVA.
  ## 单因素随机截距 ANOVA 的无偏 ICC。
  ok <- is.finite(y) & !is.na(group)
  y <- y[ok]; group <- droplevels(as.factor(group[ok]))
  if (nlevels(group) < 2 || length(y) < 4) return(NA_real_)
  fit <- try(stats::aov(y ~ group), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  ms <- summary(fit)[[1]][["Mean Sq"]]
  if (length(ms) < 2) return(NA_real_)
  n_bar <- mean(table(group))
  sigma_b <- max(ms[1] - ms[2], 0) / n_bar
  sigma_w <- ms[2]
  sigma_b / (sigma_b + sigma_w)
}

spearman_brown <- function(rho, m) {
  ## Reliability of the mean of m exchangeable observations with ICC rho.
  ## m 个可交换观测均值的可靠性。
  ifelse(is.finite(rho) & is.finite(m) & m > 0,
         (m * rho) / (1 + (m - 1) * rho),
         NA_real_)
}

design_effect <- function(rho, m) 1 + (m - 1) * rho

n_eff <- function(n, m, rho) (n * m) / design_effect(rho, m)

n_roi_for_reliability <- function(rho, R = 0.8) {
  ## Smallest integer m such that Spearman-Brown reliability >= R.
  ## 使可靠性达到 R 所需的最少 ROI 数。
  if (!is.finite(rho) || rho <= 0) return(Inf)
  if (rho >= R) return(1L)
  m <- (R * (1 - rho)) / (rho * (1 - R))
  as.integer(ceiling(m - 1e-10))
}

## ---------------------------------------------------------------------------
## Aggregation AOI -> patient / AOI 到患者的聚合
## ---------------------------------------------------------------------------
aggregate_aoi <- function(expr, patient, segment = NULL, nuclei = NULL,
                          method = c("mean", "median", "weighted_mean", "blup"),
                          min_n = 1L) {
  method <- match.arg(method)
  stopifnot(length(expr) == length(patient))
  df <- data.frame(
    expr = as.numeric(expr),
    patient = as.character(patient),
    segment = if (is.null(segment)) NA_character_ else as.character(segment),
    nuclei = if (is.null(nuclei)) 1 else as.numeric(nuclei),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$expr), , drop = FALSE]
  split_key <- if (all(is.na(df$segment))) df$patient else paste(df$patient, df$segment, sep = "\t")
  pieces <- split(df, split_key)
  rows <- lapply(pieces, function(d) {
    n <- nrow(d)
    val <- switch(method,
      mean = mean(d$expr),
      median = stats::median(d$expr),
      weighted_mean = {
        w <- d$nuclei
        w[!is.finite(w) | w <= 0] <- 1
        stats::weighted.mean(d$expr, w)
      },
      blup = {
        if (n < 2) mean(d$expr) else mean(d$expr)  # placeholder; BLUP is computed across patients
      }
    )
    data.frame(
      patient_id = d$patient[1],
      segment = d$segment[1],
      n_aoi = n,
      score = val,
      sd_within = if (n >= 2) stats::sd(d$expr) else NA_real_,
      flagged_low_n = n < min_n,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (method == "blup") {
    ## Shrink patient means toward the grand mean, separately per segment.
    ## 按分区分别把患者均值向总均值收缩。
    out$score <- blup_shrink(out$score, out$patient_id, out$segment, out$n_aoi)
  }
  out
}

blup_shrink <- function(y, patient, segment, n_aoi) {
  ## Closed-form random-intercept BLUP using method-of-moments variance components.
  ## 用矩估计方差成分给出的随机截距 BLUP（不依赖 lme4，便于无包环境）。
  y_out <- y
  for (s in unique(segment)) {
    idx <- which(segment == s & is.finite(y))
    if (length(idx) < 3) next
    rho <- icc_oneway(y[idx], patient[idx])
    if (!is.finite(rho) || rho <= 0) next
    mu <- mean(y[idx])
    ## reliability of each patient's mean; shrink factor = R
    R <- spearman_brown(rho, n_aoi[idx])
    y_out[idx] <- mu + R * (y[idx] - mu)
  }
  y_out
}

## ---------------------------------------------------------------------------
## Spillover index / 渗漏指数
## ---------------------------------------------------------------------------
spillover_index <- function(expr_mat, genes, transform = c("mean_z", "mean")) {
  ## expr_mat: genes x AOIs, already on the log2-normalized scale.
  ## 输入为 log2 归一化后的 基因 x AOI 矩阵。
  transform <- match.arg(transform)
  present <- intersect(genes, rownames(expr_mat))
  if (!length(present)) {
    warning("None of the spillover genes are in the matrix: ", paste(genes, collapse = ", "))
    return(rep(NA_real_, ncol(expr_mat)))
  }
  sub <- expr_mat[present, , drop = FALSE]
  if (transform == "mean_z") {
    sub <- t(scale(t(sub)))
    sub[!is.finite(sub)] <- 0
  }
  colMeans(sub, na.rm = TRUE)
}

## ---------------------------------------------------------------------------
## Contrast helpers / 对比辅助
## ---------------------------------------------------------------------------
sanitize_level <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", as.character(x))
}

make_segment_group <- function(segment, group) {
  factor(paste(sanitize_level(segment), sanitize_level(group), sep = "."))
}

## ---------------------------------------------------------------------------
## I/O helpers / 读写辅助
## ---------------------------------------------------------------------------
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ", path, "  (", nrow(x), " rows)")
  invisible(path)
}

read_annotation <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Reading .xlsx requires the readxl package.", call. = FALSE)
    }
    as.data.frame(readxl::read_excel(path, sheet = sheet %||% 1), stringsAsFactors = FALSE)
  } else {
    stop("Unsupported annotation format: ", ext, call. = FALSE)
  }
}

confound_check <- function(df, a, b, context = "") {
  ## Perfect confounding: every level of a maps to exactly one level of b (or vice versa).
  ## 完全混杂：a 的每个水平只对应 b 的一个水平（或反过来）。
  tab <- table(df[[a]], df[[b]], useNA = "no")
  a_pure <- all(rowSums(tab > 0) <= 1)
  b_pure <- all(colSums(tab > 0) <= 1)
  list(
    table = tab,
    a_confounded_with_b = a_pure && ncol(tab) > 1 && nrow(tab) > 1,
    b_confounded_with_a = b_pure && ncol(tab) > 1 && nrow(tab) > 1,
    context = context
  )
}

report_confound <- function(chk) {
  if (isTRUE(chk$a_confounded_with_b) || isTRUE(chk$b_confounded_with_a)) {
    warning(chk$context, " appears perfectly confounded. The contrast is not identifiable. ",
            "Declare this in the paper; do not 'adjust' and hope.",
            immediate. = TRUE)
    print(chk$table)
  }
  invisible(chk)
}

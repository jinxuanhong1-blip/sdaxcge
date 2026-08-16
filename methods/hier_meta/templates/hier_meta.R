#!/usr/bin/env Rscript
# hier_meta.R -- frequentist companion to templates/hier_meta.py
#
# OPEN COHORTS ONLY. This script never downloads EGA/dbGaP. It reads the same
# results/*.csv contract as the Python template and fits metafor::rma.uni
# (REML + Hartung-Knapp). The Bayesian quadrature model lives in hier_meta.py;
# a brms formula is printed so a DAC-approved re-run can be copied, not executed
# against restricted files.
#
# Usage:
#   Rscript methods/hier_meta/templates/hier_meta.R \
#       --results-dir results --gene TACSTD2 --endpoint DCB --outdir methods/hier_meta/out/TACSTD2
#   Rscript methods/hier_meta/templates/hier_meta.R --demo

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Install metafor: install.packages('metafor')", call. = FALSE)
  }
})

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) if (is.null(a) || is.na(a) || identical(a, "")) b else a

parse_args <- function(args) {
  out <- list(
    results_dir = "results",
    gene = "TACSTD2",
    endpoint = "DCB",
    outdir = file.path("methods", "hier_meta", "out"),
    demo = FALSE,
    access_policy = "open_only"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--demo") {
      out$demo <- TRUE
    } else if (i < length(args) && startsWith(key, "--")) {
      name <- gsub("-", "_", sub("^--", "", key))
      out[[name]] <- args[[i + 1L]]
      i <- i + 1L
    }
    i <- i + 1L
  }
  out
}

opt <- parse_args(args)
if (isTRUE(opt$demo) || identical(opt$demo, "TRUE")) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(file_arg)) {
    dirname(sub("^--file=", "", file_arg[[1]]))
  } else {
    file.path("methods", "hier_meta", "templates")
  }
  opt$results_dir <- file.path(dirname(here), "example", "results")
  message("[hier_meta.R] DEMO MODE: synthetic inputs. Numbers mean nothing.")
}

csvs <- list.files(opt$results_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(csvs)) {
  message("[hier_meta.R] no CSV under ", opt$results_dir,
          ". This template reads results/*.csv if present and otherwise exits.")
  quit(status = 0)
}

norm <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
alias <- c(
  cohortid = "cohort_id", cohort = "cohort_id", gse = "cohort_id",
  accession = "cohort_id", gene = "gene", genesymbol = "gene",
  endpoint = "endpoint", outcome = "endpoint",
  unit = "unit", yi = "yi", hedgesg = "yi", smd = "yi", estimate = "yi",
  sei = "sei", se = "sei", stderr = "sei",
  ni = "ni", n = "ni", access = "access", availability = "access"
)

read_one <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- vapply(names(d), function(nm) {
    key <- norm(nm)
    if (key %in% names(alias)) unname(alias[[key]]) else nm
  }, character(1))
  d$`.file` <- path
  d
}

raw <- do.call(rbind, lapply(csvs, read_one))
need <- c("cohort_id", "gene", "endpoint", "yi", "sei")
if (!all(need %in% names(raw))) {
  stop("CSV contract missing columns: ", paste(setdiff(need, names(raw)), collapse = ", "))
}

open_tok <- c("open", "public", "geo", "arrayexpress", "ae", "zenodo", "pride")
restr_tok <- c("ega", "egad", "egas", "dbgap", "phs", "controlled", "restricted", "dac")
acc <- if ("access" %in% names(raw)) norm(raw$access) else rep("", nrow(raw))
is_restr <- acc %in% restr_tok | grepl("^(ega|phs)", acc)
cid <- toupper(as.character(raw$cohort_id))
is_open <- acc %in% open_tok |
  grepl("^GSE|^GDS|^E-MTAB|^E-GEOD|^SYNTH_OPEN", cid)
keep_access <- if (identical(opt$access_policy, "open_only")) {
  is_open & !is_restr
} else {
  !is_restr
}

keep <- keep_access &
  norm(raw$gene) == norm(opt$gene) &
  norm(raw$endpoint) == norm(opt$endpoint) &
  is.finite(as.numeric(raw$yi)) &
  is.finite(as.numeric(raw$sei)) &
  as.numeric(raw$sei) > 0

dat <- raw[keep, , drop = FALSE]
# Default overlap policy matches Python: one row per sample_group, prefer pretreatment.
if ("sample_group" %in% names(dat)) {
  sg <- as.character(dat$sample_group)
  timed <- if ("tissue_timing" %in% names(dat)) norm(dat$tissue_timing) else rep("", nrow(dat))
  prec <- c("pretreatment", "baseline", "ontreatment", "posttreatment")
  rank <- match(timed, prec, nomatch = length(prec) + 1L)
  drop <- logical(nrow(dat))
  for (g in unique(sg[nzchar(sg)])) {
    idx <- which(sg == g)
    if (length(idx) > 1L) {
      keep_i <- idx[order(rank[idx], as.numeric(dat$sei[idx]))][1]
      drop[setdiff(idx, keep_i)] <- TRUE
    }
  }
  dat <- dat[!drop, , drop = FALSE]
}
message("[hier_meta.R] ", nrow(raw), " rows in, ", nrow(dat), " retained for ",
        opt$gene, " / ", opt$endpoint, " (open_only)")

if (nrow(dat) < 2) {
  message("[hier_meta.R] fewer than 2 open cohorts; not pooling.")
  quit(status = 0)
}

yi <- as.numeric(dat$yi)
sei <- as.numeric(dat$sei)
# Harm orientation for a response SMD stored as responder-minus-nonresponder.
yi <- -yi

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

fe <- metafor::rma.uni(yi = yi, sei = sei, method = "FE")
dl <- metafor::rma.uni(yi = yi, sei = sei, method = "DL")
pm <- metafor::rma.uni(yi = yi, sei = sei, method = "PM")
reml <- metafor::rma.uni(yi = yi, sei = sei, method = "REML", test = "knha")

pack <- function(fit, method) {
  data.frame(
    method = method,
    k = fit$k,
    estimate = as.numeric(fit$b),
    se = fit$se,
    ci_low = fit$ci.lb,
    ci_high = fit$ci.ub,
    pvalue = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2,
    Q = fit$QE,
    Q_pvalue = fit$QEp,
    stringsAsFactors = FALSE
  )
}

pooled <- rbind(
  pack(fe, "fixed_effect"),
  pack(dl, "random_DL"),
  pack(pm, "random_PM"),
  pack(reml, "random_REML_HKSJ")
)
utils::write.csv(pooled, file.path(opt$outdir, "r_pooled_estimates.csv"), row.names = FALSE)

egger <- tryCatch(metafor::regtest(reml, model = "lm"), error = function(e) NULL)
if (!is.null(egger)) {
  utils::write.csv(
    data.frame(test = "egger_metafor", statistic = unname(egger$zval),
               pvalue = egger$pval),
    file.path(opt$outdir, "r_egger.csv"), row.names = FALSE
  )
}

tf <- tryCatch(metafor::trimfill(reml), error = function(e) NULL)
if (!is.null(tf)) {
  utils::write.csv(
    data.frame(k0 = tf$k0, estimate = as.numeric(tf$b),
               ci_low = tf$ci.lb, ci_high = tf$ci.ub),
    file.path(opt$outdir, "r_trimfill.csv"), row.names = FALSE
  )
}

message("[hier_meta.R] REML+HKSJ estimate = ", signif(as.numeric(reml$b), 3),
        " [", signif(reml$ci.lb, 3), ", ", signif(reml$ci.ub, 3), "]")
message("[hier_meta.R] wrote frequentist companions under ", opt$outdir)
message("[hier_meta.R] Bayesian primary estimator is hier_meta.py (tau quadrature).")
message("[hier_meta.R] Optional brms formula (do not run on EGA files):")
message("  brms::brm(yi | se(sei) ~ 1 + (1 | cohort_id),")
message("            prior = c(prior(normal(0, 1), class = Intercept),")
message("                      prior(normal(0, 0.5), class = sd)),")
message("            data = dat, family = gaussian())")

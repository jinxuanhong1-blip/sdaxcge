#!/usr/bin/env Rscript
# Sample-size gates for methods/wgcna/playbook.md
# 本脚本只做门槛判定，不构建网络、不计算疗效 p 值。
#
# Decisions:
#   forbid  — do not construct a network / do not run the association
#   explore — allowed but must be labeled exploratory
#   allow   — method may run as a primary analysis of that class
#
# Passing a gate is not a sample-size justification.

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

.norm_decision <- function(decision, reason) {
  data.frame(
    decision = decision,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

#' Soft-threshold suggestion from the WGCNA FAQ (Dec 2017 table).
#' @param n unique donors after QC
#' @param network_type "signed", "signed hybrid", or "unsigned"
faq_soft_power <- function(n, network_type = "signed") {
  n <- as.integer(n)
  signed <- tolower(network_type) %in% c("signed")
  if (is.na(n) || n < 1) {
    return(NA_integer_)
  }
  if (n < 20) {
    return(if (signed) 18L else 9L)
  }
  if (n < 30) {
    return(if (signed) 16L else 8L)
  }
  if (n < 40) {
    return(if (signed) 14L else 7L)
  }
  if (signed) 12L else 6L
}

#' Gate for constructing a bulk WGCNA network.
#' @param n_donors unique patients after QC (not aliquots)
gate_bulk_network <- function(n_donors) {
  n <- as.integer(n_donors)
  if (is.na(n) || n < 15) {
    return(.norm_decision(
      "forbid",
      "n<15: WGCNA FAQ — correlations too noisy; use locked junction score (playbook §7)"
    ))
  }
  if (n < 20) {
    return(.norm_decision(
      "explore",
      "15<=n<20: exploratory signed network only; FAQ power 18; no hub signature; no covariate ME-ICI model"
    ))
  }
  if (n < 30) {
    return(.norm_decision(
      "explore",
      "20<=n<30: conservative network allowed; ME-trait unadjusted and labeled exploratory"
    ))
  }
  if (n < 40) {
    return(.norm_decision(
      "allow",
      "30<=n<40: preferred band for module-trait; still do not train a new ICI-supervised signature"
    ))
  }
  .norm_decision(
    "allow",
    "n>=40: standard discovery cohort; lock the module before applying to ICI sets"
  )
}

#' Gate for donor-level ME (or locked score) versus ICI labels.
#' @param n_donors donors with both ME and the endpoint
#' @param n_arm_a,n_arm_b sizes of the two response groups (binary endpoint)
#' @param n_events events for Cox (PFS/OS); leave NA if not survival
#' @param endpoint "binary", "continuous", or "survival"
gate_me_ici <- function(n_donors,
                        n_arm_a = NA_integer_,
                        n_arm_b = NA_integer_,
                        n_events = NA_integer_,
                        endpoint = c("binary", "continuous", "survival")) {
  endpoint <- match.arg(endpoint)
  n <- as.integer(n_donors)

  if (endpoint == "continuous") {
    if (is.na(n) || n < 4) {
      return(.norm_decision("forbid", "n<4: project convention — Spearman/bicor is NA"))
    }
    if (n < 15) {
      return(.norm_decision("explore", "4<=n<15: report rho with n; do not over-interpret"))
    }
    return(.norm_decision("allow", "n>=15: Spearman/bicor of a locked ME vs continuous trait"))
  }

  if (endpoint == "survival") {
    ev <- as.integer(n_events %||% NA_integer_)
    if (is.na(ev) || ev < 10) {
      return(.norm_decision(
        "forbid",
        "survival events<10: do not fit Cox as a primary ME-ICI test"
      ))
    }
    return(.norm_decision(
      "explore",
      "events>=10: Cox per SD of locked ME; median KM is display only; no cutoff search"
    ))
  }

  a <- as.integer(n_arm_a)
  b <- as.integer(n_arm_b)
  if (is.na(a) || is.na(b) || a < 5 || b < 5) {
    return(.norm_decision(
      "forbid",
      "binary: <5 donors in an arm — do not run Mann-Whitney / logistic"
    ))
  }
  if (a < 10 || b < 10) {
    return(.norm_decision(
      "explore",
      "binary: 5-9 per arm — Mann-Whitney + rank-biserial as exploratory supplement only"
    ))
  }
  .norm_decision(
    "allow",
    "binary: >=10 per arm — Mann-Whitney of locked ME; still not a trained signature"
  )
}

#' Gate for hdWGCNA network construction in one compartment (e.g. epithelial).
#' @param n_donors donors contributing cells to this compartment
#' @param n_metacells metacells actually built in the analysis group
#' @param min_cells_failed_frac fraction of sample×type groups dropped by min_cells
gate_hdwgcna_network <- function(n_donors,
                                 n_metacells,
                                 min_cells_failed_frac = 0) {
  n_mc <- as.integer(n_metacells)
  n <- as.integer(n_donors)
  fail <- as.numeric(min_cells_failed_frac %||% 0)

  if (is.na(n_mc) || n_mc < 20) {
    return(.norm_decision(
      "forbid",
      "metacells<20: do not ConstructNetwork; project a locked module or use §7 score"
    ))
  }
  if (!is.na(fail) && fail > 0.5) {
    return(.norm_decision(
      "forbid",
      ">50% sample×type groups failed min_cells; do not lower min_cells to keep everyone"
    ))
  }
  if (is.na(n) || n < 15) {
    return(.norm_decision(
      "explore",
      "metacells ok but donor n<15: network is cell-state description only; ME-ICI still gated by donors"
    ))
  }
  if (n_mc < 30) {
    return(.norm_decision(
      "explore",
      "20-29 metacells: borderline hdWGCNA; label exploratory"
    ))
  }
  .norm_decision(
    "allow",
    ">=30 metacells and donor n>=15: hdWGCNA construction allowed in this compartment"
  )
}

#' One-row report for a cohort (TSV-friendly).
gate_report_row <- function(cohort,
                            n_donors,
                            n_arm_a = NA_integer_,
                            n_arm_b = NA_integer_,
                            n_events = NA_integer_,
                            n_metacells = NA_integer_,
                            assay = c("bulk", "scrna"),
                            endpoint = c("binary", "continuous", "survival")) {
  assay <- match.arg(assay)
  endpoint <- match.arg(endpoint)

  net <- if (assay == "bulk") {
    gate_bulk_network(n_donors)
  } else {
    gate_hdwgcna_network(n_donors, n_metacells)
  }
  assoc <- gate_me_ici(n_donors, n_arm_a, n_arm_b, n_events, endpoint)

  data.frame(
    cohort = cohort,
    assay = assay,
    endpoint = endpoint,
    n_donors = as.integer(n_donors),
    n_arm_a = as.integer(n_arm_a %||% NA_integer_),
    n_arm_b = as.integer(n_arm_b %||% NA_integer_),
    n_events = as.integer(n_events %||% NA_integer_),
    n_metacells = as.integer(n_metacells %||% NA_integer_),
    faq_power_signed = faq_soft_power(n_donors, "signed"),
    network_decision = net$decision,
    network_reason = net$reason,
    association_decision = assoc$decision,
    association_reason = assoc$reason,
    stringsAsFactors = FALSE
  )
}

# Self-check (does not use real GEO labels).
if (identical(Sys.getenv("WGCNA_GATES_SELFTEST"), "1")) {
  stopifnot(gate_bulk_network(14)$decision == "forbid")
  stopifnot(gate_bulk_network(16)$decision == "explore")
  stopifnot(gate_bulk_network(35)$decision == "allow")
  stopifnot(faq_soft_power(16, "signed") == 18L)
  stopifnot(gate_me_ici(16, 8, 8, endpoint = "binary")$decision == "explore")
  stopifnot(gate_me_ici(16, 4, 12, endpoint = "binary")$decision == "forbid")
  stopifnot(gate_me_ici(3, endpoint = "continuous")$decision == "forbid")
  stopifnot(gate_hdwgcna_network(20, 12)$decision == "forbid")
  stopifnot(gate_hdwgcna_network(20, 40)$decision == "allow")
  message("sample_size_gates.R self-test OK")
}

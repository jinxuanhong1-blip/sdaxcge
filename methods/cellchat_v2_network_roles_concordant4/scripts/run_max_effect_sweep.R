#!/usr/bin/env Rscript
# Max-effect sweep of barrier/inhibitory outgoing probabilities.
# Same concordant-4 patients, same full-library LogNormalize, same cap.
# Does not replace the primary Q4/Q1 table in FINDING.md.
# Probabilities follow the CellChat Hill model and are checked against computeCommunProb.

Sys.setenv(CELLCHAT_SOURCE_ONLY = "1")
ca <- commandArgs(trailingOnly = FALSE)
src_file <- sub("^--file=", "", ca[grep("^--file=", ca)])
source(file.path(dirname(normalizePath(src_file)), "run_cellchat_v2.R"))
logmsg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ")
  flush.console()
}

SWEEP_VERSION <- "sweep1"
BENCH <- parse_int("--benchmark", 0L)
SWEEP_CACHE <- file.path(parse_opt("--cache", "/tmp/cellchat_v2_cache"), SWEEP_VERSION)
dir.create(SWEEP_CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SWEEP_CACHE, "mat"), recursive = TRUE, showWarnings = FALSE)

DIR_MAX <- file.path(OUT, "results", "max_effect")
DIR_MAX_TAB <- file.path(DIR_MAX, "tables")
DIR_MAX_FIG <- file.path(DIR_MAX, "figures")
dir.create(DIR_MAX_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MAX_FIG, recursive = TRUE, showWarnings = FALSE)
options(digits = 15, scipen = 999)

SPLITS <- c("q4q1", "top10", "top20", "top25", "top30", "detected")
RECVS <- c("TNK", "CD8")
UNIVERSES <- c("edges14", "barrier_pathway")
MEANS <- data.frame(
  mean_id = c("trunc10", "trunc05", "trimean", "thresh10", "thresh05", "median"),
  type = c("truncatedMean", "truncatedMean", "triMean", "thresholdedMean", "thresholdedMean", "median"),
  trim = c(0.1, 0.05, 0.1, 0.1, 0.05, 0.1),
  stringsAsFactors = FALSE
)
FILTERS <- c(
  "all", "drop_zero_delta", "detected_ge3", "min_arm_30", "min_arm_50",
  "min_recv_50", "min_arm_30_detected", "both_arms_positive"
)
PROB_METRICS <- c("edge_sum", "edge_mean_detected", "pathway_sum", "pathway_mean_detected")
FOLD_METRICS <- c("log2_fold_edge", "log2_fold_pathway", "logodds_edge")
CD8_131907 <- c("CD8 low T", "Cytotoxic CD8+ T", "Exhausted CD8+ T", "Naive CD8+ T")
N_SLOTS <- length(SPLITS) * length(RECVS) * nrow(MEANS) * length(UNIVERSES)
VALID <- new.env(parent = emptyenv())
VALID$formula <- FALSE
VALID$td1 <- FALSE

logmsg(
  "SWEEP", SWEEP_VERSION, "slots", N_SLOTS,
  "cohorts", paste(COHORTS, collapse = ","),
  "limit", LIMIT, "benchmark", BENCH
)

# ---------------------------------------------------------------------------
# CellChat-equivalent probabilities (population.size applied later)
# ---------------------------------------------------------------------------

order_lr <- function(lr) {
  lr <- lr[seq_len(nrow(lr)), , drop = FALSE]
  if ("annotation" %in% names(lr) && length(unique(stats::na.omit(lr$annotation))) > 1) {
    lr$annotation <- factor(lr$annotation, levels = c(
      "Secreted Signaling", "ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"
    ))
    lr <- lr[order(lr$annotation), , drop = FALSE]
    lr$annotation <- as.character(lr$annotation)
  }
  lr
}

quietly <- function(expr) {
  con <- file(nullfile(), open = "wt")
  sink(con, type = "output")
  on.exit({
    sink(type = "output")
    close(con)
  }, add = TRUE)
  force(expr)
}

mean_fun <- function(type, trim) {
  ns <- asNamespace("CellChat")
  switch(type,
    triMean = ns$triMean,
    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
    thresholdedMean = function(x) ns$thresholdedMean(x, trim = trim, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    stop("unknown mean type ", type)
  )
}

fast_prob <- function(object, type, trim, lr.use) {
  ns <- asNamespace("CellChat")
  lr.use <- order_lr(lr.use)
  data <- as.matrix(object@data.signaling)
  mx <- suppressWarnings(max(data))
  n_group <- nlevels(object@idents)
  empty <- array(
    0,
    dim = c(n_group, n_group, nrow(lr.use)),
    dimnames = list(levels(object@idents), levels(object@idents), lr.use$interaction_name)
  )
  if (!is.finite(mx) || mx <= 0) return(list(prob = empty, expr_max = 0))
  group <- object@idents
  data.use <- data / mx
  data.use.avg <- stats::aggregate(t(data.use), list(group), FUN = mean_fun(type, trim))
  data.use.avg <- t(data.use.avg[, -1, drop = FALSE])
  colnames(data.use.avg) <- levels(group)
  if (is.null(rownames(data.use.avg)) || any(rownames(data.use.avg) == "")) {
    rownames(data.use.avg) <- rownames(data)
  }
  geneL <- as.character(lr.use$ligand)
  geneR <- as.character(lr.use$receptor)
  dataLavg <- ns$computeExpr_LR(geneL, data.use.avg, object@DB$complex)
  dataRavg <- ns$computeExpr_LR(geneR, data.use.avg, object@DB$complex)
  coA <- ns$computeExpr_coreceptor(object@DB$cofactor, data.use.avg, lr.use, type = "A")
  coI <- ns$computeExpr_coreceptor(object@DB$cofactor, data.use.avg, lr.use, type = "I")
  dataRavg <- dataRavg * coA / coI
  Kh <- 0.5
  nHill <- 1
  agonist_idx <- which(!is.na(lr.use$agonist) & lr.use$agonist != "")
  antagonist_idx <- which(!is.na(lr.use$antagonist) & lr.use$antagonist != "")
  prob <- empty
  for (i in seq_len(nrow(lr.use))) {
    dataLR <- outer(as.numeric(dataLavg[i, ]), as.numeric(dataRavg[i, ]))
    P1 <- dataLR^nHill / (Kh^nHill + dataLR^nHill)
    if (isTRUE(sum(P1) == 0)) next
    P2 <- matrix(1, n_group, n_group)
    P3 <- P2
    if (i %in% agonist_idx) {
      ag <- ns$computeExpr_agonist(
        data.use = data.use.avg, pairLRsig = lr.use, cofactor_input = object@DB$cofactor,
        index.agonist = i, Kh = Kh, n = nHill
      )
      P2 <- outer(as.numeric(ag), as.numeric(ag))
    }
    if (i %in% antagonist_idx) {
      an <- ns$computeExpr_antagonist(
        data.use = data.use.avg, pairLRsig = lr.use, cofactor_input = object@DB$cofactor,
        index.antagonist = i, Kh = Kh, n = nHill
      )
      P3 <- outer(as.numeric(an), as.numeric(an))
    }
    sl <- P1 * P2 * P3
    sl[!is.finite(sl)] <- 0
    prob[, , i] <- sl
  }
  list(prob = prob, expr_max = mx)
}

slice_named <- function(prob, sender, target) {
  nms <- dimnames(prob)[[3]]
  setNames(as.numeric(prob[sender, target, ]), nms)
}

max_abs_named <- function(a, b) {
  nms <- intersect(names(a), names(b))
  if (!length(nms)) return(Inf)
  max(abs(a[nms] - b[nms]))
}

make_cc <- function(norm, group, db) {
  meta <- data.frame(
    labels = group,
    samples = factor("sample1"),
    row.names = colnames(norm),
    stringsAsFactors = FALSE
  )
  cc <- quietly(createCellChat(object = as.matrix(norm), meta = meta, group.by = "labels"))
  cc@DB <- db
  cc@DB$interaction <- db$interaction[seq_len(nrow(db$interaction)), , drop = FALSE]
  quietly(subsetData(cc))
}

percentile_split <- function(x, frac) {
  r <- rank(as.numeric(x), ties.method = "first")
  n <- length(r)
  k <- floor(n * frac)
  if (k < MIN_ARM) return(list(high = rep(FALSE, n), low = rep(FALSE, n), ok = FALSE))
  list(high = r > (n - k), low = r <= k, ok = TRUE)
}

assign_groups <- function(mal, recv, cldn4_log, cldn4_count, split) {
  group <- rep(NA_character_, length(mal))
  idx <- which(mal)
  if (split == "detected") {
    sp <- list(
      high = cldn4_count[idx] > 0,
      low = cldn4_count[idx] == 0,
      ok = TRUE
    )
  } else if (split == "q4q1") {
    sp <- quartile_high_low(cldn4_log[idx])
  } else {
    frac <- switch(split, top10 = 0.10, top20 = 0.20, top25 = 0.25, top30 = 0.30, NA_real_)
    if (!is.finite(frac)) stop("bad split ", split)
    sp <- percentile_split(cldn4_log[idx], frac)
  }
  n_high <- sum(sp$high)
  n_low <- sum(sp$low)
  if (!isTRUE(sp$ok) || n_high < MIN_ARM || n_low < MIN_ARM) {
    return(list(ok = FALSE, reason = "arm_floor", n_high = n_high, n_low = n_low, n_recv = NA_integer_))
  }
  group[idx[sp$high]] <- "CLDN4_high"
  group[idx[sp$low]] <- "CLDN4_low"
  group[recv & !mal] <- "TNK"
  n_recv <- sum(group == "TNK", na.rm = TRUE)
  if (n_recv < MIN_TNK) {
    return(list(ok = FALSE, reason = "recv_floor", n_high = n_high, n_low = n_low, n_recv = n_recv))
  }
  list(ok = TRUE, reason = "ok", group = group, n_high = n_high, n_low = n_low, n_recv = n_recv)
}

cap_keep <- function(group, seed) {
  set.seed(seed)
  keep <- rep(TRUE, length(group))
  for (g in GROUPS) {
    idx <- which(group == g)
    if (length(idx) > MAX_CELLS) {
      drop <- sample(idx, length(idx) - MAX_CELLS)
      keep[drop] <- FALSE
    }
  }
  keep[is.na(group)] <- FALSE
  keep
}

safe_id <- function(cohort, patient) {
  paste0(cohort, "__", gsub("[^A-Za-z0-9._-]", "_", patient))
}

spec_path <- function(cohort, patient) file.path(SWEEP_CACHE, paste0(safe_id(cohort, patient), ".rds"))
mat_path <- function(cohort, patient) file.path(SWEEP_CACHE, "mat", paste0(safe_id(cohort, patient), ".rds"))

empty_spec <- function(universe, split, receiver, mean_id, type, trim, status, n_high, n_low, n_recv) {
  data.frame(
    universe = universe, split = split, receiver = receiver, mean_id = mean_id,
    type = type, trim = trim, status = status,
    n_high = n_high, n_low = n_low, n_recv = n_recv,
    n_high_used = NA_integer_, n_low_used = NA_integer_, n_recv_used = NA_integer_,
    expr_max = NA_real_, stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Score one prepared patient
# ---------------------------------------------------------------------------

score_prepared <- function(prep, db, lr_map, gene_map, fam_of) {
  cohort <- prep$cohort
  patient <- prep$patient
  fp <- spec_path(cohort, patient)
  if (file.exists(fp)) {
    old <- readRDS(fp)
    if (identical(old$version, SWEEP_VERSION) && nrow(old$specs) == N_SLOTS) {
      logmsg("  cache", cohort, patient)
      if (cohort == "GSE189357" && patient == "TD1") validate_td1_edges(old)
      return(old)
    }
  }
  t0 <- Sys.time()
  counts <- prep$counts
  lib <- as.numeric(prep$lib)
  mal <- as.logical(prep$mal)
  cldn4 <- gene_row(counts, "CLDN4")
  cldn4_log <- log1p(cldn4 / pmax(lib, 1) * 1e4)
  spec_rows <- vector("list", N_SLOTS)
  edge_rows <- vector("list", N_SLOTS)
  slot <- 0L
  for (split in SPLITS) {
    for (recv_name in RECVS) {
      recv <- if (recv_name == "CD8") prep$cd8 else prep$tnk
      built <- assign_groups(mal, as.logical(recv), cldn4_log, cldn4, split)
      if (!isTRUE(built$ok)) {
        for (univ in UNIVERSES) {
          for (mi in seq_len(nrow(MEANS))) {
            slot <- slot + 1L
            spec_rows[[slot]] <- empty_spec(
              univ, split, recv_name, MEANS$mean_id[mi], MEANS$type[mi], MEANS$trim[mi],
              built$reason, built$n_high, built$n_low, built$n_recv
            )
          }
        }
        next
      }
      seed <- seed_of(paste(cohort, patient, "cap"))
      keep <- cap_keep(built$group, seed)
      norm <- lognorm_full(counts[, keep, drop = FALSE], lib[keep])
      present <- Matrix::rowSums(norm) > 0
      if (sum(present) < 2) {
        for (univ in UNIVERSES) {
          for (mi in seq_len(nrow(MEANS))) {
            slot <- slot + 1L
            spec_rows[[slot]] <- empty_spec(
              univ, split, recv_name, MEANS$mean_id[mi], MEANS$type[mi], MEANS$trim[mi],
              "no_expr", built$n_high, built$n_low, built$n_recv
            )
          }
        }
        next
      }
      norm <- norm[present, , drop = FALSE]
      group <- factor(built$group[keep], levels = GROUPS)
      n_high_used <- sum(group == "CLDN4_high")
      n_low_used <- sum(group == "CLDN4_low")
      n_recv_used <- sum(group == "TNK")
      cc <- make_cc(norm, group, db)
      for (univ in UNIVERSES) {
        cc_u <- tryCatch(subset_signaling(cc, gene_map[[univ]]), error = function(e) NULL)
        lr <- lr_map[[univ]]
        for (mi in seq_len(nrow(MEANS))) {
          slot <- slot + 1L
          if (is.null(cc_u)) {
            spec_rows[[slot]] <- empty_spec(
              univ, split, recv_name, MEANS$mean_id[mi], MEANS$type[mi], MEANS$trim[mi],
              "no_genes", built$n_high, built$n_low, built$n_recv
            )
            next
          }
          fp_prob <- fast_prob(cc_u, MEANS$type[mi], MEANS$trim[mi], lr)
          if (!VALID$formula && split == "q4q1" && recv_name == "TNK" &&
              univ == "edges14" && MEANS$mean_id[mi] == "trunc10") {
            validate_formula(cc_u, lr, fp_prob$prob, n_high_used, n_low_used, n_recv_used)
          }
          if (!VALID$td1 && cohort == "GSE189357" && patient == "TD1" &&
              split == "q4q1" && recv_name == "TNK" &&
              univ == "edges14" && MEANS$mean_id[mi] == "trunc10") {
            validate_td1(fp_prob$prob, n_high_used, n_low_used, n_recv_used)
          }
          ph <- slice_named(fp_prob$prob, "CLDN4_high", "TNK")
          pl <- slice_named(fp_prob$prob, "CLDN4_low", "TNK")
          role <- unname(fam_of[lr$interaction_name])
          role[is.na(role)] <- "extra"
          edge_rows[[slot]] <- data.frame(
            universe = univ, split = split, receiver = recv_name, mean_id = MEANS$mean_id[mi],
            interaction_name = lr$interaction_name,
            pathway_name = as.character(lr$pathway_name),
            role = role,
            prob_high = unname(ph[lr$interaction_name]),
            prob_low = unname(pl[lr$interaction_name]),
            stringsAsFactors = FALSE
          )
          edge_rows[[slot]]$prob_high[!is.finite(edge_rows[[slot]]$prob_high)] <- 0
          edge_rows[[slot]]$prob_low[!is.finite(edge_rows[[slot]]$prob_low)] <- 0
          spec_rows[[slot]] <- data.frame(
            universe = univ, split = split, receiver = recv_name, mean_id = MEANS$mean_id[mi],
            type = MEANS$type[mi], trim = MEANS$trim[mi], status = "ok",
            n_high = built$n_high, n_low = built$n_low, n_recv = built$n_recv,
            n_high_used = n_high_used, n_low_used = n_low_used, n_recv_used = n_recv_used,
            expr_max = fp_prob$expr_max, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (slot != N_SLOTS) stop("slot count ", slot, " != ", N_SLOTS, " for ", cohort, " ", patient)
  specs <- do.call(rbind, spec_rows)
  edges <- do.call(rbind, edge_rows[!vapply(edge_rows, is.null, logical(1))])
  out <- list(
    version = SWEEP_VERSION, cohort = cohort, patient = patient,
    n_mal = prep$n_mal, n_tnk = prep$n_tnk, n_cd8 = prep$n_cd8,
    specs = specs, edges = edges
  )
  saveRDS(out, fp, compress = FALSE)
  logmsg(
    "  scored", cohort, patient,
    "mal", prep$n_mal, "tnk", prep$n_tnk, "cd8", prep$n_cd8,
    "ok_slots", sum(specs$status == "ok"),
    "sec", round(as.numeric(Sys.time() - t0, units = "secs"), 1)
  )
  out
}

validate_formula <- function(cc_u, lr, fast, n_high_used, n_low_used, n_recv_used) {
  ref <- quietly(computeCommunProb(
    cc_u, type = "truncatedMean", trim = 0.1, LR.use = lr,
    raw.use = TRUE, population.size = FALSE, nboot = 1L, seed.use = 1L
  ))
  d_raw <- max(c(
    max_abs_named(slice_named(fast, "CLDN4_high", "TNK"), slice_named(ref@net$prob, "CLDN4_high", "TNK")),
    max_abs_named(slice_named(fast, "CLDN4_low", "TNK"), slice_named(ref@net$prob, "CLDN4_low", "TNK"))
  ))
  ref_pop <- quietly(computeCommunProb(
    cc_u, type = "truncatedMean", trim = 0.1, LR.use = lr,
    raw.use = TRUE, population.size = TRUE, nboot = 1L, seed.use = 1L
  ))
  Ntot <- n_high_used + n_low_used + n_recv_used
  scale_h <- (n_high_used / Ntot) * (n_recv_used / Ntot)
  scale_l <- (n_low_used / Ntot) * (n_recv_used / Ntot)
  raw_h <- slice_named(fast, "CLDN4_high", "TNK")
  raw_l <- slice_named(fast, "CLDN4_low", "TNK")
  d_pop <- max(c(
    max_abs_named(raw_h * scale_h, slice_named(ref_pop@net$prob, "CLDN4_high", "TNK")),
    max_abs_named(raw_l * scale_l, slice_named(ref_pop@net$prob, "CLDN4_low", "TNK"))
  ))
  logmsg("  formula check max|Δ| raw", signif(d_raw, 4), "pop", signif(d_pop, 4))
  if (!is.finite(d_raw) || d_raw > 1e-6 || !is.finite(d_pop) || d_pop > 1e-6) {
    stop("fast_prob does not match computeCommunProb (raw ", d_raw, ", pop ", d_pop, ")")
  }
  VALID$formula <- TRUE
}

validate_td1_edges <- function(rec) {
  sp <- rec$specs
  hit <- which(
    sp$universe == "edges14" & sp$split == "q4q1" & sp$receiver == "TNK" &
      sp$mean_id == "trunc10" & sp$status == "ok"
  )
  if (!length(hit)) stop("TD1 reference spec missing")
  sp1 <- sp[hit[[1]], , drop = FALSE]
  ed <- rec$edges
  ed <- ed[
    ed$universe == "edges14" & ed$split == "q4q1" & ed$receiver == "TNK" & ed$mean_id == "trunc10",
    , drop = FALSE
  ]
  pub <- read.delim(file.path(OUT, "results", "tables", "per_patient_pairs.tsv"), stringsAsFactors = FALSE)
  pub <- pub[pub$cohort == "GSE189357" & pub$patient == "TD1", , drop = FALSE]
  Ntot <- sp1$n_high_used + sp1$n_low_used + sp1$n_recv_used
  scale_h <- (sp1$n_high_used / Ntot) * (sp1$n_recv_used / Ntot)
  scale_l <- (sp1$n_low_used / Ntot) * (sp1$n_recv_used / Ntot)
  m <- match(pub$interaction_name, ed$interaction_name)
  if (anyNA(m)) stop("TD1 cache is missing a published pair")
  d_h <- max(abs(ed$prob_high[m] * scale_h - pub$prob_high))
  d_l <- max(abs(ed$prob_low[m] * scale_l - pub$prob_low))
  logmsg("  TD1 cache vs published max|Δ| high", signif(d_h, 4), "low", signif(d_l, 4))
  if (!is.finite(d_h) || d_h > 1e-6 || !is.finite(d_l) || d_l > 1e-6) {
    stop("TD1 cache does not match the published pair table")
  }
  VALID$td1 <- TRUE
}

validate_td1 <- function(fast, n_high_used, n_low_used, n_recv_used) {
  pub_path <- file.path(OUT, "results", "tables", "per_patient_pairs.tsv")
  pub <- read.delim(pub_path, stringsAsFactors = FALSE)
  pub <- pub[pub$cohort == "GSE189357" & pub$patient == "TD1", , drop = FALSE]
  Ntot <- n_high_used + n_low_used + n_recv_used
  scale_h <- (n_high_used / Ntot) * (n_recv_used / Ntot)
  scale_l <- (n_low_used / Ntot) * (n_recv_used / Ntot)
  raw_h <- slice_named(fast, "CLDN4_high", "TNK")
  raw_l <- slice_named(fast, "CLDN4_low", "TNK")
  d_h <- max(abs(unname(raw_h[pub$interaction_name]) * scale_h - pub$prob_high))
  d_l <- max(abs(unname(raw_l[pub$interaction_name]) * scale_l - pub$prob_low))
  logmsg("  TD1 vs published max|Δ| high", signif(d_h, 4), "low", signif(d_l, 4))
  if (!is.finite(d_h) || d_h > 1e-6 || !is.finite(d_l) || d_l > 1e-6) {
    stop("TD1 probabilities do not match the published pair table")
  }
  VALID$td1 <- TRUE
}

# ---------------------------------------------------------------------------
# Loaders — return prepared counts, do not run the primary score_unit
# ---------------------------------------------------------------------------

marker_cd8 <- function(mat, mal) {
  gene_row(mat, "CD8A") > 0 & pos_any(mat, c("CD3D", "CD3E")) & !mal
}

save_prep <- function(prep) {
  saveRDS(
    prep[c("cohort", "patient", "counts", "lib", "mal", "tnk", "cd8", "n_mal", "n_tnk", "n_cd8")],
    mat_path(prep$cohort, prep$patient),
    compress = FALSE
  )
  invisible(prep)
}

prep_from_counts <- function(counts, lib, mal, tnk, cd8, cohort, patient, needed) {
  counts <- subset_needed(as_sparse(counts), needed)
  mal <- as.logical(mal) & !is.na(mal)
  tnk <- as.logical(tnk) & !is.na(tnk) & !mal
  cd8 <- as.logical(cd8) & !is.na(cd8) & !mal
  if (is.null(lib)) lib <- Matrix::colSums(counts)
  lib <- as.numeric(lib)
  prep <- list(
    cohort = cohort, patient = patient, counts = counts, lib = lib,
    mal = mal, tnk = tnk, cd8 = cd8,
    n_mal = as.integer(sum(mal)), n_tnk = as.integer(sum(tnk)), n_cd8 = as.integer(sum(cd8))
  )
  save_prep(prep)
}

load_prepared_gse123902 <- function(needed) {
  logmsg("==== GSE123902 ====")
  units <- read.delim(file.path(HERE, "data", "GSE123902_marker_units.tsv"), stringsAsFactors = FALSE)
  tumor <- units[units$tissue %in% c("PRIMARY", "METASTASIS"), , drop = FALSE]
  tumor <- tumor[order(tumor$patient, ifelse(tumor$tissue == "PRIMARY", 0, 1)), ]
  tumor <- tumor[!duplicated(tumor$patient), ]
  csv_dir <- file.path(RAW, "GSE123902", "csv")
  files <- list.files(csv_dir, pattern = "_dense\\.csv\\.gz$", full.names = TRUE)
  if (LIMIT > 0) tumor <- head(tumor, LIMIT)
  out <- list()
  for (i in seq_len(nrow(tumor))) {
    patient <- as.character(tumor$patient[i])
    if (file.exists(spec_path("GSE123902", patient)) && file.exists(mat_path("GSE123902", patient))) {
      out[[patient]] <- readRDS(mat_path("GSE123902", patient))
      next
    }
    fp <- file.path(csv_dir, as.character(tumor$file[i]))
    if (!file.exists(fp)) {
      hit <- grep(paste0("_", patient, "_"), files, value = TRUE)
      hit <- hit[!grepl("_NORMAL_", hit)]
      hit <- hit[order(!grepl("_PRIMARY_", hit))]
      if (!length(hit)) stop("missing csv for ", patient)
      fp <- hit[[1]]
    }
    logmsg("  fread", basename(fp))
    dt <- data.table::fread(
      cmd = paste("gzip -dc", shQuote(fp)), sep = ",",
      header = TRUE, data.table = FALSE, showProgress = FALSE
    )
    genes <- colnames(dt)[-1]
    mat <- t(as.matrix(dt[, -1, drop = FALSE]))
    storage.mode(mat) <- "double"
    rownames(mat) <- genes
    colnames(mat) <- paste0(patient, "_", as.character(dt[[1]]))
    lib <- colSums(mat)
    rm(dt)
    mat <- harmonize_genes(as_sparse(mat), needed)
    mal <- pos_any(mat, EPI) & (gene_row(mat, "PTPRC") == 0)
    tnk <- pos_any(mat, TNK_MARKERS) & !mal
    cd8 <- marker_cd8(mat, mal)
    out[[patient]] <- prep_from_counts(mat, lib, mal, tnk, cd8, "GSE123902", patient, needed)
    gc(verbose = FALSE)
  }
  out
}

load_prepared_gse131907 <- function(needed) {
  logmsg("==== GSE131907 ====")
  slim <- file.path(RAW, "GSE131907", "slim")
  mat <- Matrix::readMM(file.path(slim, "matrix.mtx"))
  genes <- readLines(file.path(slim, "genes.tsv"))
  cells <- read.delim(file.path(slim, "cells.tsv"), stringsAsFactors = FALSE)
  rownames(mat) <- genes
  colnames(mat) <- cells$cell
  mat <- harmonize_genes(as_sparse(mat), needed)
  samples <- unique(as.character(cells$sample))
  if (LIMIT > 0) samples <- head(samples, LIMIT)
  out <- list()
  for (s in samples) {
    if (file.exists(spec_path("GSE131907", s)) && file.exists(mat_path("GSE131907", s))) {
      out[[s]] <- readRDS(mat_path("GSE131907", s))
      next
    }
    idx <- which(cells$sample == s)
    mal <- cells$malignant[idx] == 1
    tnk <- cells$tnk[idx] == 1 & !mal
    cd8 <- cells$cell_subtype[idx] %in% CD8_131907 & !mal
    out[[s]] <- prep_from_counts(
      mat[, idx, drop = FALSE], cells$libsize[idx], mal, tnk, cd8,
      "GSE131907", s, needed
    )
  }
  rm(mat)
  gc(verbose = FALSE)
  out
}

load_prepared_gse205335 <- function(needed) {
  logmsg("==== GSE205335 ====")
  blob <- readRDS(file.path(RAW, "GSE205335", "slim_counts.rds"))
  mat <- harmonize_genes(as_sparse(blob$mat), needed)
  ident <- blob$ident
  rm(blob)
  is_normal <- grepl("^Normal ", ident$tissue %||% "")
  is_normal[is.na(is_normal)] <- FALSE
  mal_all <- !is.na(ident$lineage.sub) & ident$lineage.sub == "Malignant cells" & !is_normal
  tnk_all <- !is.na(ident$lineage.total) & ident$lineage.total == "T/NK cells" & !is_normal
  cd8_all <- !is.na(ident$lineage.sub) & ident$lineage.sub == "CD8+ T cells" & !is_normal
  locked <- read.delim(file.path(HERE, "data", "GSE205335_patients.tsv"), stringsAsFactors = FALSE)
  keep <- as.character(locked$patient[locked$n_malignant > 0])
  if (LIMIT > 0) keep <- head(keep, LIMIT)
  out <- list()
  for (pt in keep) {
    if (file.exists(spec_path("GSE205335", pt)) && file.exists(mat_path("GSE205335", pt))) {
      out[[pt]] <- readRDS(mat_path("GSE205335", pt))
      next
    }
    idx <- which(ident$patient == pt & !is_normal)
    if (!length(idx)) stop("no cells for ", pt)
    out[[pt]] <- prep_from_counts(
      mat[, idx, drop = FALSE], ident$libsize[idx],
      mal_all[idx], tnk_all[idx], cd8_all[idx],
      "GSE205335", pt, needed
    )
  }
  rm(mat)
  gc(verbose = FALSE)
  out
}

load_prepared_gse189357 <- function(needed) {
  logmsg("==== GSE189357 ====")
  meta <- read.delim(file.path(HERE, "data", "GSE189357_sample_metadata.tsv"), stringsAsFactors = FALSE)
  ex <- file.path(RAW, "GSE189357", "raw")
  if (LIMIT > 0) meta <- head(meta, LIMIT)
  out <- list()
  for (i in seq_len(nrow(meta))) {
    patient <- meta$patient[i]
    if (file.exists(spec_path("GSE189357", patient)) && file.exists(mat_path("GSE189357", patient))) {
      out[[patient]] <- readRDS(mat_path("GSE189357", patient))
      next
    }
    gsm <- meta$gsm[i]
    prefix <- file.path(ex, paste0(gsm, "_", patient))
    mtx <- paste0(prefix, "_matrix.mtx.gz")
    cells <- paste0(prefix, "_barcodes.tsv.gz")
    features <- paste0(prefix, "_features.tsv.gz")
    logmsg("  ReadMtx", patient)
    mat <- Seurat::ReadMtx(
      mtx = mtx, cells = cells, features = features,
      cell.column = 1, feature.column = 2, unique.features = TRUE
    )
    colnames(mat) <- paste0(patient, "_", colnames(mat))
    lib <- Matrix::colSums(mat)
    mat <- harmonize_genes(as_sparse(mat), needed)
    mal <- pos_any(mat, EPI) & (gene_row(mat, "PTPRC") == 0)
    tnk <- pos_any(mat, TNK_MARKERS) & !mal
    cd8 <- marker_cd8(mat, mal)
    out[[patient]] <- prep_from_counts(mat, lib, mal, tnk, cd8, "GSE189357", patient, needed)
    gc(verbose = FALSE)
  }
  out
}

# ---------------------------------------------------------------------------
# Patient-level metrics and the sweep summary
# ---------------------------------------------------------------------------

logit_diff <- function(ph, pl) {
  clip <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  both <- ph > 0 & pl > 0
  if (!any(both)) return(NA_real_)
  lp <- function(p) log(clip(p) / (1 - clip(p)))
  mean(lp(ph[both]) - lp(pl[both]))
}

apply_pop <- function(ph, pl, n_high_used, n_low_used, n_recv_used, pop) {
  if (!pop) return(list(ph = ph, pl = pl))
  Ntot <- n_high_used + n_low_used + n_recv_used
  if (!is.finite(Ntot) || Ntot <= 0) return(list(ph = ph * NA, pl = pl * NA))
  list(
    ph = ph * (n_high_used / Ntot) * (n_recv_used / Ntot),
    pl = pl * (n_low_used / Ntot) * (n_recv_used / Ntot)
  )
}

metrics_one <- function(ed, n_high_used, n_low_used, n_recv_used, pop, barrier_pw) {
  bucket <- new.env(parent = emptyenv())
  bucket$rows <- list()
  add <- function(metric, family, sum_high, sum_low, n_detected, delta) {
    log2_fold <- if (is.finite(sum_high) && is.finite(sum_low) && sum_high > 0 && sum_low > 0) {
      log2(sum_high / sum_low)
    } else NA_real_
    bucket$rows[[length(bucket$rows) + 1]] <- data.frame(
      metric = metric, family = family,
      sum_high = sum_high, sum_low = sum_low, delta = delta,
      n_detected = n_detected, log2_fold = log2_fold,
      stringsAsFactors = FALSE
    )
  }
  for (fam in c("barrier_inhibitory", "ifn_recruit")) {
    sub <- ed[ed$role == fam, , drop = FALSE]
    if (!nrow(sub)) next
    sc <- apply_pop(sub$prob_high, sub$prob_low, n_high_used, n_low_used, n_recv_used, pop)
    det <- sc$ph > 0 | sc$pl > 0
    add("edge_sum", fam, sum(sc$ph), sum(sc$pl), sum(det), sum(sc$ph) - sum(sc$pl))
    if (any(det)) {
      add(
        "edge_mean_detected", fam,
        mean(sc$ph[det]), mean(sc$pl[det]), sum(det),
        mean(sc$ph[det] - sc$pl[det])
      )
    }
    add(
      "logodds_edge", fam,
      sum(sc$ph), sum(sc$pl), sum(sc$ph > 0 & sc$pl > 0),
      logit_diff(sc$ph, sc$pl)
    )
    add(
      "log2_fold_edge", fam,
      sum(sc$ph), sum(sc$pl), sum(det),
      if (sum(sc$ph) > 0 && sum(sc$pl) > 0) log2(sum(sc$ph) / sum(sc$pl)) else NA_real_
    )
  }
  pw_ed <- ed[ed$pathway_name %in% barrier_pw, , drop = FALSE]
  if (nrow(pw_ed) && any(ed$universe[1] == "barrier_pathway")) {
    sc <- apply_pop(pw_ed$prob_high, pw_ed$prob_low, n_high_used, n_low_used, n_recv_used, pop)
    pw_ed$ph <- sc$ph
    pw_ed$pl <- sc$pl
    pw_sum_h <- tapply(pw_ed$ph, pw_ed$pathway_name, sum)
    pw_sum_l <- tapply(pw_ed$pl, pw_ed$pathway_name, sum)
    pw_delta <- pw_sum_h - pw_sum_l
    det_pw <- (pw_sum_h + pw_sum_l) > 0
    add(
      "pathway_sum", "barrier_inhibitory",
      sum(pw_sum_h), sum(pw_sum_l), sum(det_pw), sum(pw_delta)
    )
    if (any(det_pw)) {
      add(
        "pathway_mean_detected", "barrier_inhibitory",
        mean(pw_sum_h[det_pw]), mean(pw_sum_l[det_pw]), sum(det_pw),
        mean(pw_delta[det_pw])
      )
    }
    add(
      "log2_fold_pathway", "barrier_inhibitory",
      sum(pw_sum_h), sum(pw_sum_l), sum(det_pw),
      if (sum(pw_sum_h) > 0 && sum(pw_sum_l) > 0) log2(sum(pw_sum_h) / sum(pw_sum_l)) else NA_real_
    )
  }
  if (!length(bucket$rows)) return(NULL)
  do.call(rbind, bucket$rows)
}

build_patient_metrics <- function(recs, barrier_pw) {
  rows <- list()
  for (rec in recs) {
    ok <- rec$specs[rec$specs$status == "ok", , drop = FALSE]
    if (!nrow(ok) || is.null(rec$edges) || !nrow(rec$edges)) next
    edges <- rec$edges
    key <- paste(ok$universe, ok$split, ok$receiver, ok$mean_id, sep = "|")
    ekey <- paste(edges$universe, edges$split, edges$receiver, edges$mean_id, sep = "|")
    for (i in seq_len(nrow(ok))) {
      ed <- edges[ekey == key[i], , drop = FALSE]
      if (!nrow(ed)) next
      for (pop in c(FALSE, TRUE)) {
        met <- metrics_one(
          ed, ok$n_high_used[i], ok$n_low_used[i], ok$n_recv_used[i], pop, barrier_pw
        )
        if (is.null(met)) next
        met$cohort <- rec$cohort
        met$patient <- rec$patient
        met$universe <- ok$universe[i]
        met$split <- ok$split[i]
        met$receiver <- ok$receiver[i]
        met$mean_id <- ok$mean_id[i]
        met$pop <- pop
        met$n_high <- ok$n_high[i]
        met$n_low <- ok$n_low[i]
        met$n_recv <- ok$n_recv[i]
        met$n_high_used <- ok$n_high_used[i]
        met$n_low_used <- ok$n_low_used[i]
        met$n_recv_used <- ok$n_recv_used[i]
        met$expr_max <- ok$expr_max[i]
        met$n_mal <- rec$n_mal
        rows[[length(rows) + 1]] <- met
      }
    }
  }
  do.call(rbind, rows)
}

pass_filter <- function(df, filter) {
  if (filter == "all") return(df)
  if (filter == "drop_zero_delta") return(df[is.finite(df$delta) & df$delta != 0, , drop = FALSE])
  if (filter == "detected_ge3") return(df[is.finite(df$n_detected) & df$n_detected >= 3, , drop = FALSE])
  if (filter == "min_arm_30") return(df[df$n_high >= 30 & df$n_low >= 30, , drop = FALSE])
  if (filter == "min_arm_50") return(df[df$n_high >= 50 & df$n_low >= 50, , drop = FALSE])
  if (filter == "min_recv_50") return(df[df$n_recv >= 50, , drop = FALSE])
  if (filter == "min_arm_30_detected") {
    return(df[df$n_high >= 30 & df$n_low >= 30 & is.finite(df$n_detected) & df$n_detected >= 1, , drop = FALSE])
  }
  if (filter == "both_arms_positive") return(df[df$sum_high > 0 & df$sum_low > 0, , drop = FALSE])
  stop("unknown filter ", filter)
}

cliff_delta <- function(d) {
  d <- d[is.finite(d)]
  if (!length(d)) return(NA_real_)
  (sum(d > 0) - sum(d < 0)) / length(d)
}

cliff_notie <- function(d) {
  d <- d[is.finite(d)]
  den <- sum(d > 0) + sum(d < 0)
  if (!den) return(NA_real_)
  (sum(d > 0) - sum(d < 0)) / den
}

cohen_d <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) < 2) return(NA_real_)
  s <- stats::sd(d)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(d) / s
}

signflip_trio <- function(d, B = SIGNFLIP, seed = 3979L) {
  d <- d[is.finite(d)]
  n <- length(d)
  empty <- list(greater = NA_real_, less = NA_real_, two = NA_real_)
  if (n < 2) return(empty)
  obs <- mean(d)
  set.seed(seed)
  signs <- matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = B, ncol = n)
  means <- as.numeric(signs %*% d / n)
  list(
    greater = (1 + sum(means >= obs - 1e-15)) / (B + 1),
    less = (1 + sum(means <= obs + 1e-15)) / (B + 1),
    two = (1 + sum(abs(means) >= abs(obs) - 1e-15)) / (B + 1)
  )
}

summarize_metric_rows <- function(pat, filters) {
  keys <- unique(pat[, c(
    "universe", "split", "receiver", "mean_id", "pop", "metric", "family"
  ), drop = FALSE])
  rows <- vector("list", nrow(keys) * length(filters))
  k <- 0L
  for (i in seq_len(nrow(keys))) {
    base <- pat[
      pat$universe == keys$universe[i] & pat$split == keys$split[i] &
        pat$receiver == keys$receiver[i] & pat$mean_id == keys$mean_id[i] &
        pat$pop == keys$pop[i] & pat$metric == keys$metric[i] &
        pat$family == keys$family[i],
      , drop = FALSE
    ]
    for (filter in filters) {
      sub <- pass_filter(base, filter)
      sub <- sub[is.finite(sub$delta), , drop = FALSE]
      k <- k + 1L
      if (k %% 500L == 0L) logmsg("  summary rows", k)
      nb <- n_by_cohort(sub$cohort)
      means <- tapply(sub$delta, sub$cohort, mean)
      cohort_mean <- function(nm) if (nm %in% names(means)) unname(means[[nm]]) else NA_real_
      cm <- c(cohort_mean("GSE123902"), cohort_mean("GSE131907"), cohort_mean("GSE205335"), cohort_mean("GSE189357"))
      n_pos_cohorts <- sum(is.finite(cm) & cm > 0)
      all4 <- all(nb >= 1) && n_pos_cohorts == 4L
      d <- sub$delta
      pflip <- signflip_trio(d)
      p_g <- pflip$greater
      p_l <- pflip$less
      p_t <- pflip$two
      rows[[k]] <- data.frame(
        universe = keys$universe[i], split = keys$split[i], receiver = keys$receiver[i],
        mean_id = keys$mean_id[i], pop = keys$pop[i], metric = keys$metric[i],
        family = keys$family[i], filter = filter,
        n = length(d),
        n_gse123902 = unname(nb["GSE123902"]), n_gse131907 = unname(nb["GSE131907"]),
        n_gse205335 = unname(nb["GSE205335"]), n_gse189357 = unname(nb["GSE189357"]),
        mean_delta = if (length(d)) mean(d) else NA_real_,
        median_delta = if (length(d)) stats::median(d) else NA_real_,
        cohort_mean_123902 = cm[1], cohort_mean_131907 = cm[2],
        cohort_mean_205335 = cm[3], cohort_mean_189357 = cm[4],
        n_cohorts_positive = n_pos_cohorts, all4_positive = all4,
        p_signflip_greater = p_g, p_signflip_less = p_l, p_signflip_two = p_t, p_wilcox = wilcox_p(d),
        cliff = cliff_delta(d), cliff_notie = cliff_notie(d), cohen_d = cohen_d(d),
        n_pos = sum(d > 0), n_neg = sum(d < 0), n_zero = sum(d == 0),
        median_expr_max = if (nrow(sub)) stats::median(sub$expr_max) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$eligible <- out$family == "barrier_inhibitory" &
    out$metric %in% PROB_METRICS &
    out$all4_positive %in% TRUE &
    is.finite(out$p_signflip_greater) &
    out$p_signflip_greater <= 0.05
  out$eligible_fold <- out$family == "barrier_inhibitory" &
    out$metric %in% FOLD_METRICS &
    out$all4_positive %in% TRUE &
    is.finite(out$p_signflip_greater) &
    out$p_signflip_greater <= 0.05
  out
}

assert_reference <- function(pat) {
  pub <- read.delim(file.path(OUT, "results", "tables", "per_patient_family.tsv"), stringsAsFactors = FALSE)
  pub <- pub[pub$family == "barrier_inhibitory", c("cohort", "patient", "delta_sum"), drop = FALSE]
  mine <- pat[
    pat$universe == "edges14" & pat$split == "q4q1" & pat$receiver == "TNK" &
      pat$mean_id == "trunc10" & pat$pop %in% TRUE & pat$metric == "edge_sum" &
      pat$family == "barrier_inhibitory",
    c("cohort", "patient", "delta"), drop = FALSE
  ]
  m <- merge(pub, mine, by = c("cohort", "patient"), all = TRUE)
  if (anyNA(m$delta_sum) || anyNA(m$delta)) {
    stop("reference patient set does not match the published family table")
  }
  diff <- max(abs(m$delta_sum - m$delta))
  logmsg("reference vs published max|Δ|", signif(diff, 4), "mean", signif(mean(m$delta), 6))
  if (!is.finite(diff) || diff > 1e-6) stop("reference edge_sum does not match published delta_sum")
  invisible(mean(m$delta))
}

# ---------------------------------------------------------------------------
# Full-network residual on the winning spec
# ---------------------------------------------------------------------------

rebuild_norm <- function(prep, split, recv_name) {
  mal <- as.logical(prep$mal)
  recv <- if (recv_name == "CD8") as.logical(prep$cd8) else as.logical(prep$tnk)
  cldn4 <- gene_row(prep$counts, "CLDN4")
  lib <- as.numeric(prep$lib)
  cldn4_log <- log1p(cldn4 / pmax(lib, 1) * 1e4)
  built <- assign_groups(mal, recv, cldn4_log, cldn4, split)
  if (!isTRUE(built$ok)) return(NULL)
  keep <- cap_keep(built$group, seed_of(paste(prep$cohort, prep$patient, "cap")))
  norm <- lognorm_full(prep$counts[, keep, drop = FALSE], lib[keep])
  group <- factor(built$group[keep], levels = GROUPS)
  list(
    norm = norm, group = group, built = built,
    n_high_used = sum(group == "CLDN4_high"),
    n_low_used = sum(group == "CLDN4_low"),
    n_recv_used = sum(group == "TNK")
  )
}

residual_one <- function(prep, db, pair_df, winner) {
  built <- rebuild_norm(prep, winner$split, winner$receiver)
  if (is.null(built)) return(NULL)
  cc <- make_cc(built$norm, built$group, db)
  cc <- tryCatch(
    identifyOverExpressedGenes(
      cc, thresh.p = 0.05, thresh.fc = 0, thresh.pc = 0,
      do.fast = requireNamespace("presto", quietly = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(cc)) return(NULL)
  cc <- tryCatch(identifyOverExpressedInteractions(cc), error = function(e) NULL)
  if (is.null(cc) || is.null(cc@LR$LRsig) || !nrow(cc@LR$LRsig)) return(NULL)
  lrsig <- cc@LR$LRsig
  cols <- intersect(names(lrsig), names(pair_df))
  extra <- pair_df[!pair_df$interaction_name %in% lrsig$interaction_name, cols, drop = FALSE]
  lr <- rbind(lrsig[, cols, drop = FALSE], extra)
  mean_row <- MEANS[MEANS$mean_id == winner$mean_id, , drop = FALSE]
  fp <- fast_prob(cc, mean_row$type, mean_row$trim, lr)
  sc_h <- slice_named(fp$prob, "CLDN4_high", "TNK")
  sc_l <- slice_named(fp$prob, "CLDN4_low", "TNK")
  scaled <- apply_pop(
    sc_h, sc_l, built$n_high_used, built$n_low_used, built$n_recv_used,
    isTRUE(winner$pop)
  )
  ph <- scaled$ph
  pl <- scaled$pl
  oe_names <- intersect(lrsig$interaction_name, names(ph))
  barrier_names <- intersect(pair_df$interaction_name[pair_df$family == "barrier_inhibitory"], names(ph))
  ifn_names <- intersect(pair_df$interaction_name[pair_df$family == "ifn_recruit"], names(ph))
  pw_names <- intersect(lr$interaction_name[lr$pathway_name %in% unique(pair_df$pathway_name[pair_df$family == "barrier_inhibitory"])], names(ph))
  sum_names <- function(nms, v) sum(v[nms], na.rm = TRUE)
  data.frame(
    cohort = prep$cohort, patient = prep$patient,
    total_high = sum_names(oe_names, ph), total_low = sum_names(oe_names, pl),
    barrier_high = sum_names(barrier_names, ph), barrier_low = sum_names(barrier_names, pl),
    pathway_high = sum_names(pw_names, ph), pathway_low = sum_names(pw_names, pl),
    ifn_high = sum_names(ifn_names, ph), ifn_low = sum_names(ifn_names, pl),
    n_oe = length(oe_names), expr_max = fp$expr_max,
    stringsAsFactors = FALSE
  )
}

run_residual <- function(preps, db, pair_df, winner, pat) {
  logmsg("RESIDUAL", winner$split, winner$receiver, winner$mean_id, "pop", winner$pop)
  rows <- list()
  for (prep in preps) {
    hit <- tryCatch(
      residual_one(prep, db, pair_df, winner),
      error = function(e) {
        logmsg("  residual failed", prep$cohort, prep$patient, conditionMessage(e))
        NULL
      }
    )
    if (!is.null(hit)) rows[[length(rows) + 1]] <- hit
  }
  if (!length(rows)) return(NULL)
  res <- do.call(rbind, rows)
  res$total_delta <- res$total_high - res$total_low
  res$barrier_delta <- res$barrier_high - res$barrier_low
  res$pathway_delta <- res$pathway_high - res$pathway_low
  res$ifn_delta <- res$ifn_high - res$ifn_low
  res$excess <- ifelse(
    res$total_low > 0,
    res$barrier_high - (res$barrier_low / res$total_low) * res$total_high,
    NA_real_
  )
  res$excess_pathway <- ifelse(
    res$total_low > 0,
    res$pathway_high - (res$pathway_low / res$total_low) * res$total_high,
    NA_real_
  )
  res$composition <- ifelse(
    res$total_high > 0 & res$total_low > 0,
    res$barrier_high / res$total_high - res$barrier_low / res$total_low,
    NA_real_
  )
  head_pat <- pat[
    pat$universe == winner$universe & pat$split == winner$split &
      pat$receiver == winner$receiver & pat$mean_id == winner$mean_id &
      pat$pop %in% winner$pop & pat$metric == winner$metric &
      pat$family == "barrier_inhibitory",
    , drop = FALSE
  ]
  head_pat <- pass_filter(head_pat, winner$filter)
  res$key <- paste(res$cohort, res$patient, sep = "|")
  head_pat$key <- paste(head_pat$cohort, head_pat$patient, sep = "|")
  res$in_winner_set <- res$key %in% head_pat$key
  res$headline_delta <- head_pat$delta[match(res$key, head_pat$key)]
  res
}

effect_line <- function(d) {
  d <- d[is.finite(d)]
  data.frame(
    n = length(d),
    mean_delta = if (length(d)) mean(d) else NA_real_,
    median_delta = if (length(d)) stats::median(d) else NA_real_,
    p_signflip_greater = signflip_p(d, "greater", seed = 3979L),
    p_signflip_two = signflip_p(d, "two.sided", seed = 3979L),
    p_wilcox = wilcox_p(d),
    cliff = cliff_delta(d),
    cohen_d = cohen_d(d),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Finding and figures
# ---------------------------------------------------------------------------

fmt_delta <- function(x) {
  if (length(x) != 1 || !is.finite(x)) return("NA")
  sprintf("%+.6g", x)
}

spec_label <- function(r) {
  sprintf(
    "%s | %s | %s | %s | pop=%s | %s | filter=%s",
    r$universe, r$split, r$receiver, r$mean_id, r$pop, r$metric, r$filter
  )
}

pick_best <- function(df) {
  df <- df[order(-df$mean_delta, -df$n, -df$cliff), , drop = FALSE]
  df[1, , drop = FALSE]
}

write_figures_max <- function(summary_df, winner_pts, residual) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  library(ggplot2)
  elig <- summary_df[summary_df$eligible %in% TRUE, , drop = FALSE]
  if (nrow(elig)) {
    top <- head(elig[order(-elig$mean_delta), , drop = FALSE], 12)
    top$label <- factor(vapply(seq_len(nrow(top)), function(i) spec_label(top[i, ]), character(1)),
                        levels = rev(vapply(seq_len(nrow(top)), function(i) spec_label(top[i, ]), character(1))))
    p <- ggplot(top, aes(x = label, y = mean_delta)) +
      geom_col(fill = "#2c7fb8") +
      coord_flip() +
      labs(
        x = NULL, y = "Mean patient delta (probability)",
        title = "Largest eligible barrier outgoing specs"
      ) +
      theme_bw(base_size = 9)
    ggsave(file.path(DIR_MAX_FIG, "top_eligible_specs.png"), p, width = 11, height = 6, dpi = 120)
    ggsave(file.path(DIR_MAX_FIG, "top_eligible_specs.pdf"), p, width = 11, height = 6)
  }
  if (!is.null(winner_pts) && nrow(winner_pts)) {
    p2 <- ggplot(winner_pts, aes(x = cohort, y = delta, color = cohort)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_jitter(width = 0.15, height = 0, size = 1.6) +
      labs(x = NULL, y = "Patient delta", title = "Winning spec, one point per patient") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
    ggsave(file.path(DIR_MAX_FIG, "winner_patient_delta.png"), p2, width = 7, height = 4.5, dpi = 120)
    ggsave(file.path(DIR_MAX_FIG, "winner_patient_delta.pdf"), p2, width = 7, height = 4.5)
  }
  if (!is.null(residual) && nrow(residual)) {
    use <- residual[residual$in_winner_set %in% TRUE & is.finite(residual$barrier_delta) & is.finite(residual$total_delta), , drop = FALSE]
    if (nrow(use)) {
      p3 <- ggplot(use, aes(x = total_delta, y = barrier_delta, color = cohort)) +
        geom_hline(yintercept = 0, linewidth = 0.3) +
        geom_vline(xintercept = 0, linewidth = 0.3) +
        geom_point(size = 1.6) +
        labs(
          x = "Total overexpressed sender delta",
          y = "Barrier-edge delta (same network)",
          title = "Winning split, shared max-normalization"
        ) +
        theme_bw(base_size = 11)
      ggsave(file.path(DIR_MAX_FIG, "residual_barrier_vs_total.png"), p3, width = 7, height = 5, dpi = 120)
      ggsave(file.path(DIR_MAX_FIG, "residual_barrier_vs_total.pdf"), p3, width = 7, height = 5)
    }
  }
}

write_max_finding <- function(summary_df, pat, residual, pair_tab, ref_mean, n_patients) {
  prob <- summary_df[summary_df$family == "barrier_inhibitory" & summary_df$metric %in% PROB_METRICS, , drop = FALSE]
  elig <- prob[prob$eligible %in% TRUE, , drop = FALSE]
  fold <- summary_df[summary_df$eligible_fold %in% TRUE, , drop = FALSE]
  lines <- c(
    "# FINDING — max-effect barrier outgoing, concordant-four",
    "",
    "ADDITIVE. The primary Q4/Q1 CellChat table in `FINDING.md` is unchanged.",
    "CLDN4 only. No dual-high. No TACSTD2 gate. Concordant four only.",
    "Honest unit = patient / locked sample.",
    "",
    "Barrier/inhibitory outgoing is the JAM, NECTIN2–TIGIT, CDH1, and LGALS9 family.",
    "IFN/recruit is a separate arm and is not added into the barrier score.",
    "",
    sprintf(
      "Search space for the probability-scale gate: %d barrier specs (split × receiver × mean × population.size × universe × aggregation × patient filter).",
      nrow(prob)
    ),
    "A spec clears the gate when all four cohort means are positive and the one-sided sign-flip p is ≤ 0.05 (B = 10000, seed 3979).",
    "That p is the p of the reported spec. It is not multiplied by the number of specs.",
    "",
    sprintf(
      "Reference spec (14-pair universe, Q4/Q1, T/NK, truncatedMean trim 0.1, population.size TRUE, edge sum, all units): mean Δ = %s. This matches the primary table (max absolute patient-level difference below 1e-6).",
      fmt_delta(ref_mean)
    ),
    ""
  )
  if (!nrow(elig)) {
    lines <- c(lines, "No probability-scale spec cleared the four-cohort sign-flip gate.", "")
  } else {
    w <- pick_best(elig)
    lines <- c(
      lines,
      "## Largest probability-scale mean Δ that clears the gate",
      "",
      sprintf("**Mean Δ = %s** on %d patients.", fmt_delta(w$mean_delta), w$n),
      sprintf("Spec: `%s`.", spec_label(w)),
      sprintf(
        "Cliff's delta (zeros stay in the denominator) = %s. Cliff's delta among nonzero pairs = %s. Cohen's d (mean / sd) = %s.",
        fmt_delta(w$cliff), fmt_delta(w$cliff_notie), fmt_delta(w$cohen_d)
      ),
      sprintf(
        "Median patient Δ = %s. Sign-flip one-sided p = %s. Sign-flip two-sided p = %s. Wilcoxon p = %s.",
        fmt_delta(w$median_delta), fmt_p(w$p_signflip_greater), fmt_p(w$p_signflip_two), fmt_p(w$p_wilcox)
      ),
      sprintf(
        "Cohort n (123902 / 131907 / 205335 / 189357) = %d / %d / %d / %d.",
        w$n_gse123902, w$n_gse131907, w$n_gse205335, w$n_gse189357
      ),
      sprintf(
        "Cohort mean Δ = %s / %s / %s / %s.",
        fmt_delta(w$cohort_mean_123902), fmt_delta(w$cohort_mean_131907),
        fmt_delta(w$cohort_mean_205335), fmt_delta(w$cohort_mean_189357)
      ),
      "A cohort mean is the mean of the patients who remain after the filter. One patient makes that cohort mean equal to that patient's delta.",
      sprintf(
        "Signs: %d positive, %d negative, %d zero. Median CellChat expression max in this gene universe = %s.",
        w$n_pos, w$n_neg, w$n_zero, fmt_delta(w$median_expr_max)
      ),
      ""
    )
    same <- elig[elig$universe == "edges14" & elig$metric %in% c("edge_sum", "edge_mean_detected"), , drop = FALSE]
    if (nrow(same)) {
      s <- pick_best(same)
      lines <- c(
        lines,
        sprintf(
          "Largest eligible Δ that keeps the primary 14-pair gene universe (edge sum or mean of detected edges): mean Δ = %s, Cliff = %s, Cohen's d = %s, n = %d, spec `%s`.",
          fmt_delta(s$mean_delta), fmt_delta(s$cliff), fmt_delta(s$cohen_d), s$n, spec_label(s)
        ),
        ""
      )
    }
    stable <- elig[
      elig$n_gse123902 >= 3 & elig$n_gse131907 >= 3 &
        elig$n_gse205335 >= 3 & elig$n_gse189357 >= 3,
      , drop = FALSE
    ]
    if (nrow(stable)) {
      st <- pick_best(stable)
      lines <- c(
        lines,
        sprintf(
          "Largest eligible Δ with at least 3 patients in every cohort: mean Δ = %s, Cliff = %s, Cohen's d = %s, n = %d (%d/%d/%d/%d), cohort means %s / %s / %s / %s, spec `%s`.",
          fmt_delta(st$mean_delta), fmt_delta(st$cliff), fmt_delta(st$cohen_d), st$n,
          st$n_gse123902, st$n_gse131907, st$n_gse205335, st$n_gse189357,
          fmt_delta(st$cohort_mean_123902), fmt_delta(st$cohort_mean_131907),
          fmt_delta(st$cohort_mean_205335), fmt_delta(st$cohort_mean_189357),
          spec_label(st)
        ),
        ""
      )
    }
    all_units <- elig[elig$filter == "all", , drop = FALSE]
    if (nrow(all_units)) {
      au <- pick_best(all_units)
      lines <- c(
        lines,
        sprintf(
          "Largest eligible Δ that keeps every unit passing the cell-count floors (filter = all): mean Δ = %s, Cliff = %s, Cohen's d = %s, n = %d (%d/%d/%d/%d), spec `%s`.",
          fmt_delta(au$mean_delta), fmt_delta(au$cliff), fmt_delta(au$cohen_d), au$n,
          au$n_gse123902, au$n_gse131907, au$n_gse205335, au$n_gse189357,
          spec_label(au)
        ),
        ""
      )
    }
    top <- head(elig[order(-elig$mean_delta), , drop = FALSE], 8)
    lines <- c(lines, "### Next eligible probability specs", "", "| mean Δ | Cliff | Cohen d | n | spec |", "|---:|---:|---:|---:|---|")
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      lines <- c(lines, sprintf(
        "| %s | %s | %s | %d | `%s` |",
        fmt_delta(r$mean_delta), fmt_delta(r$cliff), fmt_delta(r$cohen_d), r$n, spec_label(r)
      ))
    }
    lines <- c(lines, "")
    lines <- c(lines, "## IFN/recruit on the same patients", "")
    barrier_pts <- pat[
      pat$family == "barrier_inhibitory" & pat$metric == w$metric &
        pat$universe == w$universe & pat$split == w$split &
        pat$receiver == w$receiver & pat$mean_id == w$mean_id &
        pat$pop %in% w$pop,
      , drop = FALSE
    ]
    barrier_pts <- pass_filter(barrier_pts, w$filter)
    barrier_pts <- barrier_pts[is.finite(barrier_pts$delta), , drop = FALSE]
    ifn_pts <- pat[
      pat$family == "ifn_recruit" & pat$metric == "edge_sum" & pat$universe == "edges14" &
        pat$split == w$split & pat$receiver == w$receiver & pat$mean_id == w$mean_id &
        pat$pop %in% w$pop,
      , drop = FALSE
    ]
    ifn_keys <- paste(barrier_pts$cohort, barrier_pts$patient, sep = "|")
    ifn_pts <- ifn_pts[paste(ifn_pts$cohort, ifn_pts$patient, sep = "|") %in% ifn_keys, , drop = FALSE]
    ifn_pts <- ifn_pts[is.finite(ifn_pts$delta), , drop = FALSE]
    if (!nrow(ifn_pts)) {
      lines <- c(lines, "IFN/recruit edge sum was not available for these patients.", "")
    } else {
      fe <- effect_line(ifn_pts$delta)
      direction <- if (!is.finite(fe$mean_delta)) "NA" else if (fe$mean_delta > 0) "high>low" else if (fe$mean_delta < 0) "low>high" else "tie"
      cm <- tapply(ifn_pts$delta, ifn_pts$cohort, mean)
      pickc <- function(nm) if (nm %in% names(cm)) unname(cm[[nm]]) else NA_real_
      p_less <- signflip_p(ifn_pts$delta, alternative = "less", seed = 3979L)
      lines <- c(
        lines,
        "Same patients as the barrier winner. The seven IFN/recruit edges are scored in the 14-pair universe, so an HLA-scale max-normalization is the one used for this arm even when the barrier winner uses the pathway universe.",
        "",
        sprintf(
          "IFN/recruit edge sum: mean Δ = %s (observed %s; thesis expect low>high), n = %d, Cliff = %s, Cohen's d = %s.",
          fmt_delta(fe$mean_delta), direction, fe$n, fmt_delta(fe$cliff), fmt_delta(fe$cohen_d)
        ),
        sprintf(
          "Thesis-direction sign-flip (low>high) p = %s. Two-sided sign-flip p = %s. Wilcoxon p = %s.",
          fmt_p(p_less), fmt_p(fe$p_signflip_two), fmt_p(fe$p_wilcox)
        ),
        sprintf(
          "Cohort mean Δ (123902 / 131907 / 205335 / 189357) = %s / %s / %s / %s.",
          fmt_delta(pickc("GSE123902")), fmt_delta(pickc("GSE131907")),
          fmt_delta(pickc("GSE205335")), fmt_delta(pickc("GSE189357"))
        ),
        ""
      )
    }
    if (!is.null(pair_tab) && nrow(pair_tab)) {
      lines <- c(
        lines,
        "### Pairs inside the winning probability",
        "",
        "| pair | n | mean Δ | cohort means 123902 / 131907 / 205335 / 189357 | sign-flip greater p |",
        "|---|---:|---:|---|---:|"
      )
      for (i in seq_len(nrow(pair_tab))) {
        r <- pair_tab[i, ]
        lines <- c(lines, sprintf(
          "| %s | %d | %s | %s / %s / %s / %s | %s |",
          r$interaction_name, r$n, fmt_delta(r$mean_delta),
          fmt_delta(r$c1), fmt_delta(r$c2), fmt_delta(r$c3), fmt_delta(r$c4),
          fmt_p(r$p_signflip_greater)
        ))
      }
      if (identical(w$metric, "pathway_sum")) {
        lines <- c(
          lines,
          "",
          sprintf(
            "Sum of these seven pair-mean deltas = %s. The winning pathway sum is %s. Other interactions in the JAM, NECTIN, CDH1, and GALECTIN pathways make up the difference.",
            fmt_delta(sum(pair_tab$mean_delta)), fmt_delta(w$mean_delta)
          ),
          ""
        )
      }
    }
  }
  if (nrow(fold)) {
    fw <- pick_best(fold)
    lines <- c(
      lines,
      "## Largest fold / log-odds that clears the same gate",
      "",
      "These are not probability differences. Log2 fold uses units with both arm sums positive. Log-odds uses edges positive on both arms, clipped to [1e-6, 1−1e-6].",
      "",
      sprintf(
        "**Mean = %s** (%s). Cliff = %s. Cohen's d = %s. Median = %s. n = %d. Sign-flip p = %s.",
        fmt_delta(fw$mean_delta), fw$metric, fmt_delta(fw$cliff), fmt_delta(fw$cohen_d),
        fmt_delta(fw$median_delta), fw$n, fmt_p(fw$p_signflip_greater)
      ),
      sprintf("Spec: `%s`.", spec_label(fw)),
      ""
    )
  }
  lines <- c(lines, "## Residual after total sender strength", "")
  if (is.null(residual) || !nrow(residual)) {
    lines <- c(lines, "Residual was not computed.", "")
  } else {
    use <- residual[residual$in_winner_set %in% TRUE, , drop = FALSE]
    b <- effect_line(use$barrier_delta)
    tline <- effect_line(use$total_delta)
    e <- effect_line(use$excess)
    comp <- effect_line(use$composition)
    ratio <- if (is.finite(tline$mean_delta) && tline$mean_delta != 0) b$mean_delta / tline$mean_delta else NA_real_
    r2 <- NA_real_
    intercept <- NA_real_
    slope <- NA_real_
    ok_lm <- use[is.finite(use$barrier_delta) & is.finite(use$total_delta), , drop = FALSE]
    if (nrow(ok_lm) >= 5 && stats::sd(ok_lm$total_delta) > 0) {
      fit <- stats::lm(barrier_delta ~ total_delta, data = ok_lm)
      r2 <- summary(fit)$r.squared
      intercept <- unname(coef(fit)[[1]])
      slope <- unname(coef(fit)[[2]])
    }
    r2_head <- NA_real_
    ok_h <- use[is.finite(use$headline_delta) & is.finite(use$total_delta), , drop = FALSE]
    if (nrow(ok_h) >= 5 && stats::sd(ok_h$total_delta) > 0) {
      r2_head <- summary(stats::lm(headline_delta ~ total_delta, data = ok_h))$r.squared
    }
    lines <- c(
      lines,
      "On the winning split, receiver, mean, and population.size, one overexpressed full network (thresh.p = 0.05) supplies both the barrier edges and the total sender sum. They share one max-normalization. The seven barrier pairs are kept even when the overexpression filter would drop them. Total sender strength sums the overexpressed interactions only.",
      "",
      sprintf(
        "Within that network, on the winning patient set (n = %d): barrier-edge mean Δ = %s (sign-flip p %s, Cliff %s, Cohen d %s). Total sender mean Δ = %s (sign-flip p %s, Cliff %s, Cohen d %s).",
        b$n, fmt_delta(b$mean_delta), fmt_p(b$p_signflip_greater), fmt_delta(b$cliff), fmt_delta(b$cohen_d),
        fmt_delta(tline$mean_delta), fmt_p(tline$p_signflip_greater), fmt_delta(tline$cliff), fmt_delta(tline$cohen_d)
      ),
      sprintf("Mean barrier Δ / mean total Δ = %s.", fmt_delta(ratio)),
      sprintf(
        "Proportional residual (barrier_high − barrier_low × total_high / total_low): mean = %s, n = %d, Cliff = %s, Cohen's d = %s, sign-flip p = %s.",
        fmt_delta(e$mean_delta), e$n, fmt_delta(e$cliff), fmt_delta(e$cohen_d), fmt_p(e$p_signflip_greater)
      ),
      sprintf(
        "Composition residual (barrier/total on high minus barrier/total on low): mean = %s, n = %d, sign-flip p = %s, Cliff = %s, Cohen's d = %s.",
        fmt_delta(comp$mean_delta), comp$n, fmt_p(comp$p_signflip_greater), fmt_delta(comp$cliff), fmt_delta(comp$cohen_d)
      ),
      sprintf(
        "OLS of within-network barrier Δ on total Δ: intercept = %s, slope = %s, R² = %s. OLS R² of the headline sweep Δ on this total Δ = %s.",
        fmt_delta(intercept), fmt_delta(slope), fmt_delta(r2), fmt_delta(r2_head)
      ),
      ""
    )
  }
  lines <- c(
    lines,
    "## What this sweep is",
    "",
    sprintf("Prepared patients with a matrix: %d.", n_patients),
    "Sender arms need at least 10 cells. Receivers need at least 20. Each group is capped at 200.",
    "population.size TRUE multiplies by (n_sender / N) × (n_receiver / N) after the Hill probability.",
    "The 14-pair universe is the normalization used in the primary table. The barrier-pathway universe renormalizes inside JAM, NECTIN, CDH1, and GALECTIN, so its absolute probabilities are a different scale.",
    "Top 25% vs bottom 25% is the same rank rule as Q4/Q1.",
    "",
    "## Reproduce",
    "",
    "```bash",
    "Rscript methods/cellchat_v2_network_roles_concordant4/scripts/run_max_effect_sweep.R --raw=/tmp/concordant4_raw",
    "```",
    ""
  )
  writeLines(lines, file.path(OUT, "FINDING_MAX_EFFECT.md"))
  logmsg("wrote FINDING_MAX_EFFECT.md")
}

pair_breakdown <- function(recs, winner, barrier_names, keep_keys) {
  rows <- list()
  for (rec in recs) {
    ed <- rec$edges
    if (is.null(ed) || !nrow(ed)) next
    sp <- rec$specs
    hit <- sp$universe == winner$universe & sp$split == winner$split &
      sp$receiver == winner$receiver & sp$mean_id == winner$mean_id & sp$status == "ok"
    if (!any(hit)) next
    sp1 <- sp[hit, , drop = FALSE][1, ]
    sub <- ed[
      ed$universe == winner$universe & ed$split == winner$split &
        ed$receiver == winner$receiver & ed$mean_id == winner$mean_id &
        ed$interaction_name %in% barrier_names,
      , drop = FALSE
    ]
    if (!nrow(sub)) next
    # Patient filter uses family-level n_detected / arms. Approximate with this pair table's patient
    # only if the patient survives the same arm filter; pair rows are per interaction.
    sc <- apply_pop(sub$prob_high, sub$prob_low, sp1$n_high_used, sp1$n_low_used, sp1$n_recv_used, isTRUE(winner$pop))
    part <- data.frame(
      cohort = rec$cohort, patient = rec$patient,
      interaction_name = sub$interaction_name,
      delta = sc$ph - sc$pl,
      n_high = sp1$n_high, n_low = sp1$n_low, n_recv = sp1$n_recv,
      stringsAsFactors = FALSE
    )
    part <- part[paste(part$cohort, part$patient, sep = "|") %in% keep_keys, , drop = FALSE]
    if (nrow(part)) rows[[length(rows) + 1]] <- part
  }
  if (!length(rows)) return(NULL)
  allp <- do.call(rbind, rows)
  out <- list()
  for (nm in barrier_names) {
    d <- allp$delta[allp$interaction_name == nm]
    cohorts <- allp$cohort[allp$interaction_name == nm]
    means <- tapply(d, cohorts, mean)
    pick <- function(k) if (k %in% names(means)) unname(means[[k]]) else NA_real_
    out[[length(out) + 1]] <- data.frame(
      interaction_name = nm, n = length(d), mean_delta = if (length(d)) mean(d) else NA_real_,
      c1 = pick("GSE123902"), c2 = pick("GSE131907"), c3 = pick("GSE205335"), c4 = pick("GSE189357"),
      p_signflip_greater = signflip_p(d, "greater", seed = 3979L),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

startup_checks <- function() {
  set.seed(1)
  x <- rnorm(2186)
  a <- quartile_high_low(x)
  b <- percentile_split(x, 0.25)
  if (!identical(a$high, b$high) || !identical(a$low, b$low)) {
    stop("top25 split does not match Q4/Q1 on n=2186")
  }
}

main_sweep <- function() {
  startup_checks()
  logmsg("Seurat", as.character(packageVersion("Seurat")), "CellChat", as.character(packageVersion("CellChat")))
  db <- CellChatDB.human
  db$interaction <- db$interaction[seq_len(nrow(db$interaction)), , drop = FALSE]
  if ("annotation" %in% names(db$interaction)) {
    keep_ann <- db$interaction$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
    # Keep every thesis pair even if its annotation were outside those three.
    pair_preview <- resolve_pairs(db)
    keep_ann <- keep_ann | db$interaction$interaction_name %in% pair_preview$interaction_name
    if (any(keep_ann)) db$interaction <- db$interaction[keep_ann, , drop = FALSE]
  }
  pair_df <- resolve_pairs(db)
  barrier_pw <- unique(pair_df$pathway_name[pair_df$family == "barrier_inhibitory"])
  barrier_lr <- db$interaction[db$interaction$pathway_name %in% barrier_pw, , drop = FALSE]
  logmsg(
    "barrier pathways", paste(barrier_pw, collapse = ","),
    "pathway LRs", nrow(barrier_lr),
    "pairs", nrow(pair_df)
  )
  lr_map <- list(edges14 = order_lr(pair_df), barrier_pathway = order_lr(barrier_lr))
  gene_map <- list(
    edges14 = genes_for_pairs(db, pair_df),
    barrier_pathway = genes_for_pairs(db, barrier_lr)
  )
  fam_of <- setNames(pair_df$family, pair_df$interaction_name)
  needed <- unique(c(
    "CLDN4", "PTPRC", EPI, TNK_MARKERS, "CD3D", "CD3E", "CD8A",
    collect_db_genes(db), names(CANON_FROM), unname(CANON_FROM),
    unlist(gene_map, use.names = FALSE)
  ))
  needed <- needed[!is.na(needed) & nzchar(needed)]
  db_head <- paste(utils::head(db$interaction$interaction_name, 3), collapse = "|")
  preps <- list()
  if ("GSE189357" %in% COHORTS) preps <- c(preps, load_prepared_gse189357(needed))
  if ("GSE123902" %in% COHORTS) preps <- c(preps, load_prepared_gse123902(needed))
  if ("GSE205335" %in% COHORTS) preps <- c(preps, load_prepared_gse205335(needed))
  if ("GSE131907" %in% COHORTS) preps <- c(preps, load_prepared_gse131907(needed))
  if (!identical(paste(utils::head(db$interaction$interaction_name, 3), collapse = "|"), db_head)) {
    stop("CellChatDB.interaction was mutated during loading")
  }
  recs <- list()
  for (prep in preps) {
    recs[[length(recs) + 1]] <- score_prepared(prep, db, lr_map, gene_map, fam_of)
  }
  logmsg("scored patients", length(recs), "formula", VALID$formula, "td1", VALID$td1)
  if (!VALID$formula) {
    logmsg("formula check not repeated: every unit was read from cache. The published-table match is the CellChat check.")
  }
  if (BENCH == 1L) {
    logmsg("benchmark stop before the full summary")
    return(invisible(recs))
  }
  logmsg("building patient metrics")
  pat <- build_patient_metrics(recs, barrier_pw)
  write_tsv(pat, file.path(DIR_MAX_TAB, "per_patient_metrics.tsv"))
  full_cohorts <- all(c("GSE123902", "GSE131907", "GSE205335", "GSE189357") %in% COHORTS) && LIMIT == 0L
  if (full_cohorts) {
    if (!VALID$td1) {
      td1_rec <- Filter(function(r) identical(r$cohort, "GSE189357") && identical(r$patient, "TD1"), recs)
      if (!length(td1_rec)) stop("TD1 result missing")
      validate_td1_edges(td1_rec[[1]])
    }
    assert_reference(pat)
  }
  # Sign-flip the barrier search and the fold search. IFN is summarized for the winner only,
  # plus every IFN edge_sum row is still in `pat` for that later step. To keep the gate honest
  # for IFN cohort means we do summarize IFN edge_sum under the same filters, which is one
  # metric × edges14 × the grid. That is part of the reported IFN line, not the barrier search.
  barrier_pat <- pat[pat$family == "barrier_inhibitory" & pat$metric %in% c(PROB_METRICS, FOLD_METRICS), , drop = FALSE]
  ifn_pat <- pat[pat$family == "ifn_recruit" & pat$metric == "edge_sum" & pat$universe == "edges14", , drop = FALSE]
  logmsg("summarizing barrier rows", nrow(barrier_pat), "ifn rows", nrow(ifn_pat))
  summary_df <- rbind(
    summarize_metric_rows(barrier_pat, FILTERS),
    summarize_metric_rows(ifn_pat, FILTERS)
  )
  write_tsv(summary_df, file.path(DIR_MAX_TAB, "summary_sweep.tsv"))
  logmsg(
    "eligible probability specs",
    sum(summary_df$eligible %in% TRUE),
    "eligible fold specs",
    sum(summary_df$eligible_fold %in% TRUE)
  )
  if (!full_cohorts) {
    logmsg("partial cohort run: tables written, finding withheld")
    return(invisible(summary_df))
  }
  prob <- summary_df[summary_df$eligible %in% TRUE, , drop = FALSE]
  winner <- if (nrow(prob)) pick_best(prob) else NULL
  residual <- NULL
  pairs <- NULL
  winner_pts <- NULL
  if (!is.null(winner)) {
    winner_pts <- pat[
      pat$universe == winner$universe & pat$split == winner$split &
        pat$receiver == winner$receiver & pat$mean_id == winner$mean_id &
        pat$pop %in% winner$pop & pat$metric == winner$metric &
        pat$family == "barrier_inhibitory",
      , drop = FALSE
    ]
    winner_pts <- pass_filter(winner_pts, winner$filter)
    winner_pts <- winner_pts[is.finite(winner_pts$delta), , drop = FALSE]
    write_tsv(winner_pts, file.path(DIR_MAX_TAB, "winner_patients.tsv"))
    pairs <- pair_breakdown(
      recs, winner, pair_df$interaction_name[pair_df$family == "barrier_inhibitory"],
      paste(winner_pts$cohort, winner_pts$patient, sep = "|")
    )
    if (!is.null(pairs)) write_tsv(pairs, file.path(DIR_MAX_TAB, "winner_pairs.tsv"))
    residual <- run_residual(preps, db, pair_df, winner, pat)
    if (!is.null(residual)) write_tsv(residual, file.path(DIR_MAX_TAB, "residual_patients.tsv"))
  }
  write_figures_max(summary_df, winner_pts, residual)
  ref_mean <- mean(pat$delta[
    pat$universe == "edges14" & pat$split == "q4q1" & pat$receiver == "TNK" &
      pat$mean_id == "trunc10" & pat$pop %in% TRUE & pat$metric == "edge_sum" &
      pat$family == "barrier_inhibitory"
  ])
  write_max_finding(summary_df, pat, residual, pairs, ref_mean, length(recs))
  writeLines(
    c(
      paste("Seurat", packageVersion("Seurat")),
      paste("CellChat", packageVersion("CellChat")),
      paste("sweep", SWEEP_VERSION),
      paste("slots", N_SLOTS),
      capture.output(sessionInfo())
    ),
    file.path(DIR_MAX, "sessionInfo.txt")
  )
  logmsg("SWEEP DONE")
  invisible(summary_df)
}

main_sweep()

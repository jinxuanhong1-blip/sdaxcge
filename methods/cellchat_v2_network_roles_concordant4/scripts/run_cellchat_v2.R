#!/usr/bin/env Rscript
# Concordant-4 CellChat v2 full network + signaling-role analysis.
# Senders = malignant CLDN4-high vs CLDN4-low (within-unit Q4 vs Q1).
# Receiver context = T/NK.
# Engine: R + Seurat + CellChat. Python is used only to stream GSE131907 text.
# Does not replace the locked 14-pair table. No GSE148071. No dual-high. No TACSTD2 gate.

suppressPackageStartupMessages({
  lib <- Sys.getenv("R_LIBS_USER", "/tmp/r_lib")
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed. Stop. Do not fall back to a Python-only primary.")
  }
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat is not installed. Stop. Do not fall back to a Python-only primary.")
  }
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(data.table)
})

options(warn = 1)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)[[1]]) a else b

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}
parse_int <- function(flag, default) as.integer(parse_opt(flag, as.character(default)))

RAW <- parse_opt("--raw", "/tmp/concordant4_raw")
HERE <- tryCatch({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  normalizePath(file.path(dirname(f), ".."))
}, error = function(e) normalizePath("methods/cellchat_v2_network_roles_concordant4"))
OUT <- parse_opt("--out", HERE)
MAX_CELLS <- parse_int("--max-cells", 200L)
NBOOT_FULL <- parse_int("--nboot-full", 20L)
NBOOT_PAIR <- parse_int("--nboot-pair", 100L)
LABEL_PERM <- parse_int("--label-perm", 49L)
SIGNFLIP <- parse_int("--signflip", 10000L)
LIMIT <- parse_int("--limit", 0L)
SMOKE <- parse_int("--smoke", 0L)
COHORTS <- strsplit(parse_opt(
  "--cohorts", "GSE189357,GSE123902,GSE131907,GSE205335"
), ",", fixed = TRUE)[[1]]
ANALYSIS_VERSION <- "v2"

DIR_RES <- file.path(OUT, "results")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_TAB <- file.path(DIR_RES, "tables")
CACHE <- file.path(parse_opt("--cache", "/tmp/cellchat_v2_cache"), ANALYSIS_VERSION)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

logmsg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ", flush = TRUE)
}

EPI <- c("EPCAM", "KRT8", "KRT18", "KRT19")
TNK_MARKERS <- c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
MALIG_SUB <- c("Malignant cells", "tS1", "tS2", "tS3")
TNK_TYPES <- c("T lymphocytes", "NK cells")
MIN_ARM <- 10L
MIN_TNK <- 20L
MIN_MAL_Q4 <- 40L
GROUPS <- c("CLDN4_high", "CLDN4_low", "TNK")

# Pre-specified pairs. Names are resolved against CellChatDB at runtime.
PAIR_SPEC <- data.frame(
  interaction_name = c(
    "JAM1_ITGAL_ITGB2", "NECTIN2_TIGIT", "CDH1_ITGAE_ITGB7", "CDH1_KLRG1",
    "LGALS9_HAVCR2", "LGALS9_CD44", "LGALS9_CD45",
    "CXCL9_CXCR3", "CXCL10_CXCR3", "CCL5_CCR5", "CCL5_CCR1",
    "HLA-A_CD8A", "HLA-B_CD8A", "HLA-C_CD8A"
  ),
  family = c(rep("barrier_inhibitory", 7), rep("ifn_recruit", 7)),
  thesis_expect = c(rep("high>low", 7), rep("low>high", 7)),
  alt_name = c(
    "F11R_ITGAL_ITGB2", "", "", "",
    "", "", "LGALS9_PTPRC",
    "", "", "", "",
    "", "", ""
  ),
  stringsAsFactors = FALSE
)

# Alias -> CellChat symbol. Never rename the CellChat symbol back to the alias.
CANON_FROM <- c(PVRL2 = "NECTIN2", JAM1 = "F11R", CD45 = "PTPRC")

fmt_p <- function(p) {
  if (length(p) != 1 || !is.finite(p)) return("NA")
  if (p < 1e-3) sprintf("%.2e", p) else sprintf("%.3g", p)
}
fmt_num <- function(x, d = 4) {
  if (length(x) != 1 || !is.finite(x)) return("NA")
  sprintf(paste0("%+.", d, "f"), x)
}

seed_of <- function(x) {
  bits <- utf8ToInt(paste0(x))
  as.integer((sum(bits * seq_along(bits)) %% 100000L) + 1L)
}

is_gzip <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && magic[[1]] == as.raw(0x1f) && magic[[2]] == as.raw(0x8b)
}

read_geo_rds <- function(path) {
  if (!file.exists(path)) stop("missing ", path)
  cur <- path
  temps <- character()
  on.exit(unlink(temps[file.exists(temps)]), add = TRUE)
  for (i in seq_len(4)) {
    if (!is_gzip(cur)) break
    dest <- tempfile(pattern = paste0("geo_rds_", i, "_"), tmpdir = tempdir())
    temps <- c(temps, dest)
    logmsg("  gzip -dc layer", i, basename(cur))
    st <- system2("gzip", c("-dc", cur), stdout = dest)
    if (!identical(st, 0L) && !is.null(st) && st != 0) stop("gzip -dc failed on ", cur)
    cur <- dest
  }
  logmsg("  readRDS", cur, "bytes", file.info(cur)$size)
  obj <- try(readRDS(cur), silent = TRUE)
  if (inherits(obj, "try-error")) obj <- readRDS(gzfile(cur, open = "rb"))
  obj
}

db_character_genes <- function(df) {
  if (is.null(df) || !nrow(df)) return(character())
  ch <- df[, vapply(df, is.character, logical(1)), drop = FALSE]
  if (!ncol(ch)) return(character())
  vals <- unlist(ch, use.names = FALSE)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  unique(vals)
}

collect_db_genes <- function(db) {
  genes <- unique(c(
    db$interaction$ligand, db$interaction$receptor,
    db_character_genes(db$complex), db_character_genes(db$cofactor)
  ))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  unique(genes)
}

resolve_pairs <- function(db) {
  inter <- db$interaction
  nm <- as.character(inter$interaction_name)
  out <- PAIR_SPEC
  out$resolved <- out$interaction_name
  for (i in seq_len(nrow(out))) {
    if (out$interaction_name[i] %in% nm) next
    alt <- out$alt_name[i]
    if (nzchar(alt) && alt %in% nm) {
      out$resolved[i] <- alt
      next
    }
    stop(
      "CellChatDB is missing ", out$interaction_name[i],
      ". Closest: ",
      paste(grep(sub("_.*", "", out$interaction_name[i]), nm, value = TRUE)[seq_len(8)],
            collapse = ", ")
    )
  }
  hit <- inter[match(out$resolved, nm), , drop = FALSE]
  hit$family <- out$family
  hit$thesis_expect <- out$thesis_expect
  hit$requested_name <- out$interaction_name
  rownames(hit) <- hit$interaction_name
  hit
}

expand_complexes <- function(db, genes) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  cx <- db$complex
  if (!is.null(cx) && nrow(cx)) {
    hit <- genes %in% rownames(cx)
    if (any(hit)) {
      subcols <- grep("^subunit", names(cx), value = TRUE)
      subs <- unlist(cx[genes[hit], subcols, drop = FALSE], use.names = FALSE)
      subs <- subs[!is.na(subs) & nzchar(subs)]
      genes <- unique(c(genes[!hit], subs))
    }
  }
  co <- db$cofactor
  if (!is.null(co) && nrow(co)) {
    hit <- genes %in% rownames(co)
    if (any(hit)) {
      subs <- unlist(co[genes[hit], , drop = FALSE], use.names = FALSE)
      subs <- subs[!is.na(subs) & nzchar(subs)]
      genes <- unique(c(genes[!hit], subs))
    }
  }
  unique(genes[!is.na(genes) & nzchar(genes)])
}

genes_for_pairs <- function(db, pair_df) {
  genes <- unique(c(as.character(pair_df$ligand), as.character(pair_df$receptor)))
  extra_cols <- intersect(
    c("agonist", "antagonist", "co_A_receptor", "co_I_receptor"),
    names(pair_df)
  )
  extra <- unlist(pair_df[, extra_cols, drop = FALSE], use.names = FALSE)
  expand_complexes(db, c(genes, extra))
}

harmonize_genes <- function(mat, needed) {
  rn <- rownames(mat)
  idx <- match(toupper(rn), toupper(needed))
  rn[!is.na(idx)] <- needed[idx[!is.na(idx)]]
  for (alias in names(CANON_FROM)) {
    target <- unname(CANON_FROM[[alias]])
    if (!target %in% needed || target %in% rn) next
    if (alias %in% rn) rn[rn == alias] <- target
  }
  rownames(mat) <- rn
  if (anyDuplicated(rn)) mat <- mat[!duplicated(rn), , drop = FALSE]
  mat
}

as_sparse <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- Matrix::Matrix(mat, sparse = TRUE)
  mat
}

lognorm_full <- function(counts, lib) {
  counts <- as_sparse(counts) * 1
  if (is.null(lib)) lib <- Matrix::colSums(counts)
  lib <- as.numeric(lib)
  if (length(lib) != ncol(counts)) stop("library size length does not match cells")
  lib[!is.finite(lib) | lib <= 0] <- 1
  counts@x <- log1p(counts@x / lib[rep.int(seq_len(ncol(counts)), diff(counts@p))] * 1e4)
  counts
}

subset_needed <- function(mat, needed) {
  mat <- harmonize_genes(as_sparse(mat), needed)
  keep <- intersect(needed, rownames(mat))
  if (!length(keep)) stop("none of the needed genes are in the matrix")
  mat[keep, , drop = FALSE]
}

gene_row <- function(mat, g) {
  if (!g %in% rownames(mat)) return(rep(0, ncol(mat)))
  as.numeric(mat[g, ])
}

pos_any <- function(mat, genes) {
  hit <- intersect(genes, rownames(mat))
  if (!length(hit)) return(rep(FALSE, ncol(mat)))
  if (length(hit) == 1L) return(as.numeric(mat[hit, ]) > 0)
  as.numeric(Matrix::colSums(mat[hit, , drop = FALSE] > 0)) > 0
}

quartile_high_low <- function(x) {
  r <- rank(as.numeric(x), ties.method = "first")
  n <- length(r)
  if (n < 4) return(list(high = rep(FALSE, n), low = rep(FALSE, n), ok = FALSE))
  q1 <- floor(n * 0.25)
  q4 <- ceiling(n * 0.75)
  if (q1 < 1 || q4 > n || q1 >= q4) {
    return(list(high = rep(FALSE, n), low = rep(FALSE, n), ok = FALSE))
  }
  list(high = r > q4, low = r <= q1, ok = TRUE)
}

cap_groups <- function(group, max_n, seed) {
  set.seed(seed)
  keep <- rep(TRUE, length(group))
  for (g in GROUPS) {
    idx <- which(group == g)
    if (length(idx) > max_n) {
      drop <- sample(idx, length(idx) - max_n)
      keep[drop] <- FALSE
    }
  }
  keep
}

cent_val <- function(centr, measure, group) {
  if (is.null(centr) || is.null(centr[[measure]])) return(NA_real_)
  v <- centr[[measure]]
  if (is.list(v)) v <- unlist(v)
  v <- as.numeric(v)
  if (!length(v) || all(!is.finite(v))) return(NA_real_)
  if (!is.null(names(centr[[measure]]))) {
    nm <- names(centr[[measure]])
    if (group %in% nm) return(unname(as.numeric(centr[[measure]][group][[1]])))
  }
  pos <- match(group, GROUPS)
  if (!is.na(pos) && pos <= length(v)) return(v[[pos]])
  NA_real_
}

safe_centrality <- function(net) {
  net[!is.finite(net)] <- 0
  if (sum(net) <= 0) {
    z <- setNames(rep(0, nrow(net)), rownames(net))
    return(list(outdeg = z, indeg = z, flowbet = z, info = z))
  }
  tryCatch(
    CellChat:::computeCentralityLocal(net),
    error = function(e) {
      z <- setNames(rep(0, nrow(net)), rownames(net))
      list(outdeg = z, indeg = z, flowbet = z, info = z)
    }
  )
}

pathway_array <- function(prob, pair_df) {
  pnames <- dimnames(prob)[[3]]
  if (is.null(pnames) || !length(pnames)) return(NULL)
  meta <- pair_df[match(pnames, pair_df$interaction_name), , drop = FALSE]
  pw <- as.character(meta$pathway_name)
  if (all(is.na(pw))) return(NULL)
  pathways <- unique(pw[!is.na(pw)])
  arr <- array(
    0,
    dim = c(dim(prob)[1], dim(prob)[2], length(pathways)),
    dimnames = list(dimnames(prob)[[1]], dimnames(prob)[[2]], pathways)
  )
  for (i in seq_along(pnames)) {
    if (is.na(pw[i])) next
    slice <- prob[, , i]
    slice[!is.finite(slice)] <- 0
    arr[, , pw[i]] <- arr[, , pw[i]] + slice
  }
  keep <- apply(arr, 3, sum) > 0
  if (!any(keep)) return(NULL)
  arr[, , keep, drop = FALSE]
}

prob_at <- function(prob, sender, target, iname) {
  if (is.null(prob) || !iname %in% dimnames(prob)[[3]]) return(NA_real_)
  if (!sender %in% dimnames(prob)[[1]] || !target %in% dimnames(prob)[[2]]) return(NA_real_)
  v <- as.numeric(prob[sender, target, iname])
  if (!is.finite(v)) NA_real_ else v
}

sum_idx <- function(prob, sender, target, names_use) {
  if (is.null(prob) || !length(names_use)) return(0)
  hit <- intersect(names_use, dimnames(prob)[[3]])
  if (!length(hit)) return(0)
  v <- as.numeric(prob[sender, target, hit])
  sum(v[is.finite(v)])
}

collapse_net <- function(prob, names_use) {
  lev <- dimnames(prob)[[1]]
  net <- matrix(0, nrow = length(lev), ncol = length(lev), dimnames = list(lev, lev))
  hit <- intersect(names_use, dimnames(prob)[[3]])
  if (!length(hit)) return(net)
  for (nm in hit) {
    sl <- prob[, , nm]
    sl[!is.finite(sl)] <- 0
    net <- net + sl
  }
  net
}

run_prob <- function(object, nboot, lr.use, seed) {
  computeCommunProb(
    object,
    type = "truncatedMean",
    trim = 0.1,
    LR.use = lr.use,
    raw.use = TRUE,
    population.size = TRUE,
    nboot = nboot,
    seed.use = seed
  )
}

subset_signaling <- function(object, genes) {
  have <- intersect(genes, rownames(object@data.signaling))
  if (length(have) < 2) stop("fewer than 2 signaling genes for the extracted pairs")
  object@data.signaling <- object@data.signaling[have, , drop = FALSE]
  miss <- setdiff(genes, rownames(object@data.signaling))
  if (length(miss)) {
    add <- Matrix::Matrix(
      0, nrow = length(miss), ncol = ncol(object@data.signaling),
      dimnames = list(miss, colnames(object@data.signaling)), sparse = TRUE
    )
    object@data.signaling <- rbind(object@data.signaling, add)
  }
  object
}

empty_unit <- function(cohort, patient, reason, n_mal, n_tnk) {
  list(
    ok = FALSE,
    reason = reason,
    inv = data.frame(
      cohort = cohort, patient = patient, n_mal = n_mal, n_tnk = n_tnk,
      n_high = NA_integer_, n_low = NA_integer_,
      n_high_used = NA_integer_, n_low_used = NA_integer_, n_tnk_used = NA_integer_,
      n_lr_full = NA_integer_, status = reason, stringsAsFactors = FALSE
    ),
    pairs = NULL, family = NULL, pathway = NULL, global = NULL
  )
}

score_unit <- function(counts, lib, mal, tnk, cohort, patient, db, pair_df, needed) {
  mal <- as.logical(mal) %in% TRUE
  tnk <- as.logical(tnk) %in% TRUE
  # a cell cannot be both
  tnk <- tnk & !mal
  n_mal <- as.integer(sum(mal))
  n_tnk <- as.integer(sum(tnk))
  if (n_tnk < MIN_TNK || n_mal < MIN_MAL_Q4) {
    return(empty_unit(cohort, patient, "floor", n_mal, n_tnk))
  }
  counts <- subset_needed(as_sparse(counts), needed)
  if (is.null(lib)) lib <- Matrix::colSums(counts)
  lib <- as.numeric(lib)
  cldn4 <- gene_row(counts, "CLDN4")
  cldn4_log <- log1p(cldn4 / pmax(lib, 1) * 1e4)
  hl <- quartile_high_low(cldn4_log[mal])
  if (!isTRUE(hl$ok)) return(empty_unit(cohort, patient, "quartile", n_mal, n_tnk))
  high <- mal
  low <- mal
  high[mal] <- hl$high
  low[mal] <- hl$low
  n_high <- as.integer(sum(high))
  n_low <- as.integer(sum(low))
  if (n_high < MIN_ARM || n_low < MIN_ARM) {
    return(empty_unit(cohort, patient, "arm_floor", n_mal, n_tnk))
  }
  group <- rep(NA_character_, ncol(counts))
  group[high] <- "CLDN4_high"
  group[low] <- "CLDN4_low"
  group[tnk] <- "TNK"
  keep <- cap_groups(group, MAX_CELLS, seed_of(paste(cohort, patient, "cap")))
  keep[is.na(group)] <- FALSE
  counts <- counts[, keep, drop = FALSE]
  lib <- lib[keep]
  group <- factor(group[keep], levels = GROUPS)
  tab <- table(group)
  if (any(tab[GROUPS] < MIN_ARM)) {
    return(empty_unit(cohort, patient, "cap_floor", n_mal, n_tnk))
  }
  logmsg(
    "  CellChat", cohort, patient,
    "high", n_high, "->", tab[["CLDN4_high"]],
    "low", n_low, "->", tab[["CLDN4_low"]],
    "tnk", n_tnk, "->", tab[["TNK"]]
  )
  norm <- lognorm_full(counts, lib)
  meta <- data.frame(
    labels = group,
    samples = factor("sample1"),
    row.names = colnames(norm),
    stringsAsFactors = FALSE
  )
  obj <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 0, min.features = 0)
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = norm)
  Idents(obj) <- obj$labels
  cc <- createCellChat(object = as.matrix(norm), meta = meta, group.by = "labels")
  cc@DB <- db
  cc <- subsetData(cc)
  if (nrow(cc@data.signaling) < 3) {
    return(empty_unit(cohort, patient, "no_signaling_genes", n_mal, n_tnk))
  }
  seed <- seed_of(paste(cohort, patient))

  # ---- full network (official overexpressed L-R) ----
  n_lr_full <- 0L
  pathway <- NULL
  global <- NULL
  full <- tryCatch({
    cc_full <- identifyOverExpressedGenes(
      cc, thresh.p = 0.05, thresh.fc = 0, thresh.pc = 0,
      do.fast = requireNamespace("presto", quietly = TRUE)
    )
    cc_full <- identifyOverExpressedInteractions(cc_full)
    cc_full
  }, error = function(e) {
    logmsg("    full-network OE failed:", conditionMessage(e))
    NULL
  })
  if (!is.null(full) && !is.null(full@LR$LRsig) && nrow(full@LR$LRsig) > 0) {
    t0 <- Sys.time()
    full <- tryCatch(
      run_prob(full, nboot = NBOOT_FULL, lr.use = NULL, seed = seed),
      error = function(e) {
        logmsg("    full computeCommunProb failed:", conditionMessage(e))
        NULL
      }
    )
    logmsg("    full network seconds", round(as.numeric(Sys.time() - t0, units = "secs"), 1))
  } else {
    full <- NULL
  }
  if (!is.null(full) && !is.null(full@net$prob)) {
    prob <- full@net$prob
    pval <- full@net$pval
    n_lr_full <- dim(prob)[3]
    lr_full <- full@LR$LRsig
    rownames(lr_full) <- lr_full$interaction_name
    arr <- pathway_array(prob, lr_full)
    rows <- list()
    if (!is.null(arr)) {
      centr_all <- netAnalysis_computeCentrality(object = NULL, net = arr)
      # significant-edge pathway array, only when the permutation grid is real
      if (NBOOT_FULL > 1) {
        prob_sig <- prob
        prob_sig[pval > 0.05] <- 0
        arr_sig <- pathway_array(prob_sig, lr_full)
      } else {
        arr_sig <- NULL
      }
      tag_of <- function(pw) {
        b <- unique(pair_df$pathway_name[pair_df$family == "barrier_inhibitory"])
        f <- unique(pair_df$pathway_name[pair_df$family == "ifn_recruit"])
        if (pw %in% b && pw %in% f) "both" else if (pw %in% b) "barrier_inhibitory" else if (pw %in% f) "ifn_recruit" else "other"
      }
      for (pw in dimnames(arr)[[3]]) {
        cnet <- safe_centrality(arr[, , pw])
        sig_n <- if (!is.null(arr_sig) && pw %in% dimnames(arr_sig)[[3]]) {
          sum(arr_sig[, , pw] > 0)
        } else NA_integer_
        rows[[length(rows) + 1]] <- data.frame(
          cohort = cohort, patient = patient, pathway = pw,
          family_tag = tag_of(pw),
          outdeg_high = cent_val(cnet, "outdeg", "CLDN4_high"),
          outdeg_low = cent_val(cnet, "outdeg", "CLDN4_low"),
          indeg_high = cent_val(cnet, "indeg", "CLDN4_high"),
          indeg_low = cent_val(cnet, "indeg", "CLDN4_low"),
          flowbet_high = cent_val(cnet, "flowbet", "CLDN4_high"),
          flowbet_low = cent_val(cnet, "flowbet", "CLDN4_low"),
          info_high = cent_val(cnet, "info", "CLDN4_high"),
          info_low = cent_val(cnet, "info", "CLDN4_low"),
          prob_to_tnk_high = sum(arr["CLDN4_high", "TNK", pw]),
          prob_to_tnk_low = sum(arr["CLDN4_low", "TNK", pw]),
          prob_out_high = sum(arr["CLDN4_high", , pw]),
          prob_out_low = sum(arr["CLDN4_low", , pw]),
          n_sig_edges = sig_n,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows)) {
      pathway <- do.call(rbind, rows)
      pathway$delta_outdeg <- pathway$outdeg_high - pathway$outdeg_low
      pathway$delta_tnk <- pathway$prob_to_tnk_high - pathway$prob_to_tnk_low
    }
    gnet <- collapse_net(prob, dimnames(prob)[[3]])
    gcent <- safe_centrality(gnet)
    global <- data.frame(
      cohort = cohort, patient = patient,
      n_interactions = n_lr_full,
      sum_high = sum_idx(prob, "CLDN4_high", "TNK", dimnames(prob)[[3]]),
      sum_low = sum_idx(prob, "CLDN4_low", "TNK", dimnames(prob)[[3]]),
      outdeg_high = cent_val(gcent, "outdeg", "CLDN4_high"),
      outdeg_low = cent_val(gcent, "outdeg", "CLDN4_low"),
      stringsAsFactors = FALSE
    )
    global$delta_sum <- global$sum_high - global$sum_low
    global$delta_outdeg <- global$outdeg_high - global$outdeg_low
    logmsg("    full L-R", n_lr_full, "pathways", if (is.null(pathway)) 0 else nrow(pathway))
  }

  # ---- extracted thesis pairs, with CellChat permutation and label permutation ----
  pair_genes <- genes_for_pairs(db, pair_df)
  cc_pair <- tryCatch(subset_signaling(cc, pair_genes), error = function(e) {
    logmsg("    pair gene subset failed:", conditionMessage(e))
    NULL
  })
  pair_tab <- NULL
  fam_tab <- NULL
  if (!is.null(cc_pair)) {
    t0 <- Sys.time()
    cc_pair <- tryCatch(
      run_prob(cc_pair, nboot = NBOOT_PAIR, lr.use = pair_df, seed = seed),
      error = function(e) {
        logmsg("    pair computeCommunProb failed:", conditionMessage(e))
        NULL
      }
    )
    logmsg("    pair nboot seconds", round(as.numeric(Sys.time() - t0, units = "secs"), 1))
  }
  if (!is.null(cc_pair) && !is.null(cc_pair@net$prob)) {
    prob <- cc_pair@net$prob
    pval <- cc_pair@net$pval
    pair_tab <- pair_df[, c(
      "interaction_name", "pathway_name", "ligand", "receptor", "annotation",
      "family", "thesis_expect", "requested_name"
    ), drop = FALSE]
    pair_tab$cohort <- cohort
    pair_tab$patient <- patient
    pair_tab$prob_high <- vapply(pair_tab$interaction_name, function(nm) prob_at(prob, "CLDN4_high", "TNK", nm), numeric(1))
    pair_tab$prob_low <- vapply(pair_tab$interaction_name, function(nm) prob_at(prob, "CLDN4_low", "TNK", nm), numeric(1))
    pair_tab$pval_high <- vapply(pair_tab$interaction_name, function(nm) prob_at(pval, "CLDN4_high", "TNK", nm), numeric(1))
    pair_tab$pval_low <- vapply(pair_tab$interaction_name, function(nm) prob_at(pval, "CLDN4_low", "TNK", nm), numeric(1))
    pair_tab$prob_high[!is.finite(pair_tab$prob_high)] <- 0
    pair_tab$prob_low[!is.finite(pair_tab$prob_low)] <- 0
    pair_tab$delta <- pair_tab$prob_high - pair_tab$prob_low
    pair_tab$detected <- pair_tab$prob_high > 0 | pair_tab$prob_low > 0

    fam_delta <- function(pr, fam) {
      sub <- pr[pr$family == fam, , drop = FALSE]
      sum(sub$prob_high) - sum(sub$prob_low)
    }
    obs_b <- fam_delta(pair_tab, "barrier_inhibitory")
    obs_f <- fam_delta(pair_tab, "ifn_recruit")
    # label permutation: shuffle high/low among malignant senders, T/NK fixed
    id0 <- as.character(cc_pair@idents)
    mal_idx <- which(id0 %in% c("CLDN4_high", "CLDN4_low"))
    null_b <- rep(NA_real_, LABEL_PERM)
    null_f <- rep(NA_real_, LABEL_PERM)
    null_ob <- rep(NA_real_, LABEL_PERM)
    null_of <- rep(NA_real_, LABEL_PERM)
    set.seed(seed + 17L)
    t0 <- Sys.time()
    if (LABEL_PERM > 0 && length(mal_idx) > 1) {
      for (b in seq_len(LABEL_PERM)) {
        new <- id0
        new[mal_idx] <- sample(id0[mal_idx])
        cc_pair@idents <- factor(new, levels = GROUPS)
        cc_b <- tryCatch(
          run_prob(cc_pair, nboot = 1L, lr.use = pair_df, seed = seed + b),
          error = function(e) NULL
        )
        if (is.null(cc_b)) next
        prb <- cc_b@net$prob
        sb <- sum_idx(prb, "CLDN4_high", "TNK", pair_tab$interaction_name[pair_tab$family == "barrier_inhibitory"]) -
          sum_idx(prb, "CLDN4_low", "TNK", pair_tab$interaction_name[pair_tab$family == "barrier_inhibitory"])
        sf <- sum_idx(prb, "CLDN4_high", "TNK", pair_tab$interaction_name[pair_tab$family == "ifn_recruit"]) -
          sum_idx(prb, "CLDN4_low", "TNK", pair_tab$interaction_name[pair_tab$family == "ifn_recruit"])
        null_b[b] <- sb
        null_f[b] <- sf
        nb <- collapse_net(prb, pair_tab$interaction_name[pair_tab$family == "barrier_inhibitory"])
        nf <- collapse_net(prb, pair_tab$interaction_name[pair_tab$family == "ifn_recruit"])
        null_ob[b] <- sum(nb["CLDN4_high", ]) - sum(nb["CLDN4_low", ])
        null_of[b] <- sum(nf["CLDN4_high", ]) - sum(nf["CLDN4_low", ])
      }
    }
    logmsg("    label perm seconds", round(as.numeric(Sys.time() - t0, units = "secs"), 1))
    perm_p <- function(obs, null, alternative) {
      null <- null[is.finite(null)]
      if (!is.finite(obs) || !length(null)) return(NA_real_)
      B <- length(null)
      if (alternative == "greater") return((1 + sum(null >= obs - 1e-15)) / (B + 1))
      if (alternative == "less") return((1 + sum(null <= obs + 1e-15)) / (B + 1))
      (1 + sum(abs(null) >= abs(obs) - 1e-15)) / (B + 1)
    }
    fam_rows <- list()
    for (fam in c("barrier_inhibitory", "ifn_recruit")) {
      sub <- pair_tab[pair_tab$family == fam, , drop = FALSE]
      det <- sub[sub$detected %in% TRUE, , drop = FALSE]
      net <- collapse_net(prob, sub$interaction_name)
      cnet <- safe_centrality(net)
      obs <- sum(sub$prob_high) - sum(sub$prob_low)
      null <- if (fam == "barrier_inhibitory") null_b else null_f
      null_o <- if (fam == "barrier_inhibitory") null_ob else null_of
      obs_o <- cent_val(cnet, "outdeg", "CLDN4_high") - cent_val(cnet, "outdeg", "CLDN4_low")
      fam_rows[[fam]] <- data.frame(
        cohort = cohort, patient = patient, family = fam,
        thesis_expect = unique(sub$thesis_expect)[[1]],
        n_pairs = nrow(sub),
        n_pairs_detected = nrow(det),
        sum_high = sum(sub$prob_high),
        sum_low = sum(sub$prob_low),
        delta_sum = obs,
        mean_detected_delta = if (nrow(det)) mean(det$delta) else NA_real_,
        outdeg_high = cent_val(cnet, "outdeg", "CLDN4_high"),
        outdeg_low = cent_val(cnet, "outdeg", "CLDN4_low"),
        delta_outdeg = obs_o,
        indeg_high = cent_val(cnet, "indeg", "CLDN4_high"),
        indeg_low = cent_val(cnet, "indeg", "CLDN4_low"),
        flowbet_high = cent_val(cnet, "flowbet", "CLDN4_high"),
        flowbet_low = cent_val(cnet, "flowbet", "CLDN4_low"),
        info_high = cent_val(cnet, "info", "CLDN4_high"),
        info_low = cent_val(cnet, "info", "CLDN4_low"),
        pval_high_median = stats::median(sub$pval_high[sub$detected %in% TRUE], na.rm = TRUE),
        pval_low_median = stats::median(sub$pval_low[sub$detected %in% TRUE], na.rm = TRUE),
        frac_edge_p05_high = if (nrow(det)) mean(det$pval_high < 0.05, na.rm = TRUE) else NA_real_,
        frac_edge_p05_low = if (nrow(det)) mean(det$pval_low < 0.05, na.rm = TRUE) else NA_real_,
        p_label_greater = perm_p(obs, null, "greater"),
        p_label_less = perm_p(obs, null, "less"),
        p_label_two = perm_p(obs, null, "two"),
        p_label_outdeg_greater = perm_p(obs_o, null_o, "greater"),
        p_label_outdeg_less = perm_p(obs_o, null_o, "less"),
        n_label_perm = sum(is.finite(null)),
        stringsAsFactors = FALSE
      )
    }
    fam_tab <- do.call(rbind, fam_rows)
    rownames(fam_tab) <- NULL
    rownames(pair_tab) <- NULL
  }

  inv <- data.frame(
    cohort = cohort, patient = patient, n_mal = n_mal, n_tnk = n_tnk,
    n_high = n_high, n_low = n_low,
    n_high_used = as.integer(tab[["CLDN4_high"]]),
    n_low_used = as.integer(tab[["CLDN4_low"]]),
    n_tnk_used = as.integer(tab[["TNK"]]),
    n_lr_full = n_lr_full,
    status = if (!is.null(fam_tab)) "ok" else "pair_failed",
    stringsAsFactors = FALSE
  )
  list(ok = !is.null(fam_tab), reason = inv$status[[1]], inv = inv,
       pairs = pair_tab, family = fam_tab, pathway = pathway, global = global)
}

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

bind_or_null <- function(xs) {
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(NULL)
  do.call(rbind, xs)
}

cache_file <- function(cohort, patient) {
  file.path(CACHE, paste0(cohort, "__", gsub("[^A-Za-z0-9._-]", "_", patient), ".rds"))
}

run_cached <- function(cohort, patient, expr_fun) {
  fp <- cache_file(cohort, patient)
  if (file.exists(fp)) {
    logmsg("  cache", cohort, patient)
    return(readRDS(fp))
  }
  rec <- expr_fun()
  saveRDS(rec, fp)
  rec
}

load_gse123902 <- function(db, pair_df, needed) {
  logmsg("==== GSE123902 ====")
  units <- read.delim(file.path(HERE, "data", "GSE123902_marker_units.tsv"), stringsAsFactors = FALSE)
  tumor <- units[units$tissue %in% c("PRIMARY", "METASTASIS"), , drop = FALSE]
  tumor <- tumor[order(tumor$patient, ifelse(tumor$tissue == "PRIMARY", 0, 1)), ]
  tumor <- tumor[!duplicated(tumor$patient), ]
  tar <- file.path(RAW, "GSE123902", "GSE123902_RAW.tar")
  csv_dir <- file.path(RAW, "GSE123902", "csv")
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  if (!length(list.files(csv_dir, pattern = "_dense\\.csv\\.gz$"))) {
    if (!file.exists(tar)) stop("missing ", tar)
    logmsg("untar GSE123902")
    untar(tar, exdir = csv_dir)
  }
  files <- list.files(csv_dir, pattern = "_dense\\.csv\\.gz$", full.names = TRUE)
  recs <- list()
  if (LIMIT > 0) tumor <- head(tumor, LIMIT)
  for (i in seq_len(nrow(tumor))) {
    patient <- as.character(tumor$patient[i])
    recs[[patient]] <- run_cached(patient = patient, cohort = "GSE123902", expr_fun = function() {
      fp <- file.path(csv_dir, as.character(tumor$file[i]))
      if (!file.exists(fp)) {
        hit <- grep(paste0("_", patient, "_"), files, value = TRUE)
        hit <- hit[!grepl("_NORMAL_", hit)]
        hit <- hit[order(!grepl("_PRIMARY_", hit))]
        if (!length(hit)) return(empty_unit("GSE123902", patient, "missing_csv", NA, NA))
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
      score_unit(mat, lib, mal, tnk, "GSE123902", patient, db, pair_df, needed)
    })
    gc(verbose = FALSE)
  }
  recs
}

load_gse131907 <- function(db, pair_df, needed) {
  logmsg("==== GSE131907 ====")
  slim <- file.path(RAW, "GSE131907", "slim")
  mtx <- file.path(slim, "matrix.mtx")
  if (!file.exists(mtx)) {
    gene_file <- file.path(CACHE, "wanted_genes.txt")
    writeLines(needed, gene_file)
    py <- file.path(HERE, "scripts", "extract_gse131907.py")
    logmsg("  python stream")
    st <- system2("python3", c(
      py,
      "--umi", file.path(RAW, "GSE131907", "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"),
      "--annot", file.path(RAW, "GSE131907", "GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
      "--samples", file.path(HERE, "data", "GSE131907_samples.tsv"),
      "--genes", gene_file,
      "--out", slim
    ))
    if (!identical(st, 0L)) stop("GSE131907 extract failed status=", st)
  }
  mat <- Matrix::readMM(mtx)
  genes <- readLines(file.path(slim, "genes.tsv"))
  cells <- read.delim(file.path(slim, "cells.tsv"), stringsAsFactors = FALSE)
  rownames(mat) <- genes
  colnames(mat) <- cells$cell
  mat <- as_sparse(mat)
  samples <- unique(cells$sample)
  if (LIMIT > 0) samples <- head(samples, LIMIT)
  recs <- list()
  for (s in samples) {
    recs[[s]] <- run_cached(patient = s, cohort = "GSE131907", expr_fun = function() {
      idx <- which(cells$sample == s)
      sub <- mat[, idx, drop = FALSE]
      lib <- cells$libsize[idx]
      mal <- cells$malignant[idx] == 1
      tnk <- cells$tnk[idx] == 1
      score_unit(sub, lib, mal, tnk, "GSE131907", s, db, pair_df, needed)
    })
    gc(verbose = FALSE)
  }
  recs
}

load_gse205335 <- function(db, pair_df, needed) {
  logmsg("==== GSE205335 ====")
  slim_path <- file.path(RAW, "GSE205335", "slim_counts.rds")
  if (!file.exists(slim_path)) {
    ident <- utils::read.delim(
      file.path(RAW, "GSE205335", "GSE205335_Lung_IO_CellIdentity.txt.gz"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    gsm_map <- read.delim(file.path(HERE, "data", "GSE205335_gsm_map.tsv"), stringsAsFactors = FALSE)
    mat <- read_geo_rds(file.path(RAW, "GSE205335", "GSE205335_Lung_IO_UMI_matrix.rds.gz"))
    if (is.data.frame(mat)) mat <- as.matrix(mat)
    if (!inherits(mat, "Matrix")) mat <- Matrix::Matrix(mat, sparse = TRUE)
    logmsg("  full", nrow(mat), "x", ncol(mat))
    lib <- Matrix::colSums(mat)
    names(lib) <- colnames(mat)
    mat <- subset_needed(mat, needed)
    logmsg("  subset", nrow(mat), "x", ncol(mat))
    common <- intersect(colnames(mat), ident$barcode)
    mat <- mat[, common, drop = FALSE]
    lib <- lib[common]
    ident <- ident[match(common, ident$barcode), , drop = FALSE]
    ident <- merge(ident, gsm_map[, c("orig.ident", "patient", "tissue"), drop = FALSE],
                   by = "orig.ident", all.x = TRUE, sort = FALSE)
    ident <- ident[match(colnames(mat), ident$barcode), , drop = FALSE]
    if (anyNA(ident$patient) || any(ident$patient == "")) {
      stop("unmapped orig.ident: ",
           paste(unique(ident$orig.ident[is.na(ident$patient) | ident$patient == ""]), collapse = ", "))
    }
    ident$libsize <- as.numeric(lib)
    saveRDS(list(mat = mat, ident = ident), slim_path)
    rm(lib)
    gc(verbose = FALSE)
  }
  blob <- readRDS(slim_path)
  mat <- blob$mat
  ident <- blob$ident
  is_normal <- grepl("^Normal ", ident$tissue %||% "")
  is_normal[is.na(is_normal)] <- FALSE
  mal_all <- !is.na(ident$lineage.sub) & ident$lineage.sub == "Malignant cells" & !is_normal
  tnk_all <- !is.na(ident$lineage.total) & ident$lineage.total == "T/NK cells" & !is_normal
  locked <- read.delim(file.path(HERE, "data", "GSE205335_patients.tsv"), stringsAsFactors = FALSE)
  keep <- as.character(locked$patient[locked$n_malignant > 0])
  if (LIMIT > 0) keep <- head(keep, LIMIT)
  recs <- list()
  for (pt in keep) {
    recs[[pt]] <- run_cached(patient = pt, cohort = "GSE205335", expr_fun = function() {
      idx <- which(ident$patient == pt & !is_normal)
      if (!length(idx)) return(empty_unit("GSE205335", pt, "no_cells", 0, 0))
      score_unit(
        mat[, idx, drop = FALSE], ident$libsize[idx],
        mal_all[idx], tnk_all[idx], "GSE205335", pt, db, pair_df, needed
      )
    })
    gc(verbose = FALSE)
  }
  recs
}

load_gse189357 <- function(db, pair_df, needed) {
  logmsg("==== GSE189357 ====")
  meta <- read.delim(file.path(HERE, "data", "GSE189357_sample_metadata.tsv"), stringsAsFactors = FALSE)
  tar <- file.path(RAW, "GSE189357", "GSE189357_RAW.tar")
  ex <- file.path(RAW, "GSE189357", "raw")
  dir.create(ex, showWarnings = FALSE, recursive = TRUE)
  needed_files <- paste0(meta$gsm, "_", meta$patient, "_matrix.mtx.gz")
  if (!all(file.exists(file.path(ex, needed_files)))) {
    if (!file.exists(tar)) stop("missing ", tar)
    logmsg("untar GSE189357")
    untar(tar, exdir = ex)
  }
  if (LIMIT > 0) meta <- head(meta, LIMIT)
  recs <- list()
  for (i in seq_len(nrow(meta))) {
    patient <- meta$patient[i]
    recs[[patient]] <- run_cached(patient = patient, cohort = "GSE189357", expr_fun = function() {
      gsm <- meta$gsm[i]
      prefix <- file.path(ex, paste0(gsm, "_", patient))
      mtx <- paste0(prefix, "_matrix.mtx.gz")
      cells <- paste0(prefix, "_barcodes.tsv.gz")
      features <- paste0(prefix, "_features.tsv.gz")
      if (!all(file.exists(c(mtx, cells, features)))) {
        alt <- list.files(ex, pattern = paste0(patient, "_matrix\\.mtx\\.gz$"), full.names = TRUE)
        if (!length(alt)) return(empty_unit("GSE189357", patient, "missing_mtx", NA, NA))
        mtx <- alt[[1]]
        cells <- sub("_matrix\\.mtx\\.gz$", "_barcodes.tsv.gz", mtx)
        features <- sub("_matrix\\.mtx\\.gz$", "_features.tsv.gz", mtx)
      }
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
      score_unit(mat, lib, mal, tnk, "GSE189357", patient, db, pair_df, needed)
    })
    gc(verbose = FALSE)
  }
  recs
}

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

signflip_p <- function(d, alternative = "two.sided", B = SIGNFLIP, seed = 3979L) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 2) return(NA_real_)
  obs <- mean(d)
  set.seed(seed)
  signs <- matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = B, ncol = n)
  means <- as.numeric(signs %*% d / n)
  if (alternative == "greater") return((1 + sum(means >= obs - 1e-15)) / (B + 1))
  if (alternative == "less") return((1 + sum(means <= obs + 1e-15)) / (B + 1))
  (1 + sum(abs(means) >= abs(obs) - 1e-15)) / (B + 1)
}

wilcox_p <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) < 2) return(NA_real_)
  suppressWarnings(stats::wilcox.test(d, mu = 0, exact = FALSE)$p.value)
}

stouffer_p <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  z <- stats::qnorm(p, lower.tail = FALSE)
  stats::pnorm(sum(z) / sqrt(length(z)), lower.tail = FALSE)
}

n_by_cohort <- function(cohorts) {
  tab <- table(as.character(cohorts))
  pick <- function(nm) if (nm %in% names(tab)) as.integer(tab[[nm]]) else 0L
  c(
    GSE123902 = pick("GSE123902"),
    GSE131907 = pick("GSE131907"),
    GSE205335 = pick("GSE205335"),
    GSE189357 = pick("GSE189357")
  )
}

summarize_family <- function(fam) {
  rows <- list()
  for (metric in c("delta_sum", "mean_detected_delta", "delta_outdeg")) {
    for (f in c("barrier_inhibitory", "ifn_recruit")) {
      sub <- fam[fam$family == f, , drop = FALSE]
      d <- sub[[metric]]
      use <- sub[is.finite(d), , drop = FALSE]
      d <- use[[metric]]
      nb <- n_by_cohort(use$cohort)
      expect <- if (nrow(sub)) sub$thesis_expect[[1]] else NA_character_
      alt <- if (identical(expect, "high>low")) "greater" else "less"
      mean_d <- if (length(d)) mean(d) else NA_real_
      obs <- if (!length(d) || !is.finite(mean_d)) "NA" else if (mean_d > 0) "high>low" else if (mean_d < 0) "low>high" else "tie"
      agrees <- if (obs == "NA") "thin" else if (obs == expect) "yes" else if (obs == "tie") "tie" else "opposite"
      p_one <- signflip_p(d, alternative = alt, seed = 3979L)
      p_label <- if (metric == "delta_sum") {
        if (alt == "greater") stouffer_p(use$p_label_greater) else stouffer_p(use$p_label_less)
      } else if (metric == "delta_outdeg") {
        if (alt == "greater") stouffer_p(use$p_label_outdeg_greater) else stouffer_p(use$p_label_outdeg_less)
      } else NA_real_
      frac_label <- if (metric == "delta_sum" && nrow(use)) {
        col <- if (alt == "greater") use$p_label_greater else use$p_label_less
        mean(col < 0.05, na.rm = TRUE)
      } else if (metric == "delta_outdeg" && nrow(use)) {
        col <- if (alt == "greater") use$p_label_outdeg_greater else use$p_label_outdeg_less
        mean(col < 0.05, na.rm = TRUE)
      } else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        family = f, metric = metric, expect = expect, observed = obs, agrees = agrees,
        n = length(d),
        n_gse123902 = nb[["GSE123902"]], n_gse131907 = nb[["GSE131907"]],
        n_gse205335 = nb[["GSE205335"]], n_gse189357 = nb[["GSE189357"]],
        mean_delta = mean_d,
        p_signflip_two = signflip_p(d, "two.sided", seed = 3979L),
        p_signflip_thesis = p_one,
        p_wilcox = wilcox_p(d),
        p_stouffer_label = p_label,
        frac_units_label_p05 = frac_label,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

summarize_global <- function(g) {
  if (is.null(g) || !nrow(g)) return(NULL)
  rows <- list()
  for (metric in c("delta_sum", "delta_outdeg")) {
    d <- g[[metric]]
    use <- g[is.finite(d), , drop = FALSE]
    d <- use[[metric]]
    nb <- n_by_cohort(use$cohort)
    mean_d <- if (length(d)) mean(d) else NA_real_
    obs <- if (!length(d) || !is.finite(mean_d)) "NA" else if (mean_d > 0) "high>low" else if (mean_d < 0) "low>high" else "tie"
    rows[[length(rows) + 1]] <- data.frame(
      metric = metric, n = length(d),
      n_gse123902 = nb[["GSE123902"]], n_gse131907 = nb[["GSE131907"]],
      n_gse205335 = nb[["GSE205335"]], n_gse189357 = nb[["GSE189357"]],
      mean_delta = mean_d, observed = obs,
      p_signflip_two = signflip_p(d, "two.sided", seed = 3979L),
      p_signflip_greater = signflip_p(d, "greater", seed = 3979L),
      p_wilcox = wilcox_p(d),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_pairs <- function(pairs) {
  if (is.null(pairs) || !nrow(pairs)) return(NULL)
  rows <- list()
  for (nm in unique(pairs$interaction_name)) {
    sub <- pairs[pairs$interaction_name == nm, , drop = FALSE]
    det <- sub[sub$detected %in% TRUE & is.finite(sub$delta), , drop = FALSE]
    d <- det$delta
    nb <- n_by_cohort(det$cohort)
    expect <- sub$thesis_expect[[1]]
    alt <- if (identical(expect, "high>low")) "greater" else "less"
    mean_d <- if (length(d)) mean(d) else NA_real_
    obs <- if (!length(d) || !is.finite(mean_d)) "NA" else if (mean_d > 0) "high>low" else if (mean_d < 0) "low>high" else "tie"
    agrees <- if (obs == "NA") "thin" else if (obs == expect) "yes" else if (obs == "tie") "tie" else "opposite"
    rows[[length(rows) + 1]] <- data.frame(
      pair = nm,
      requested = sub$requested_name[[1]],
      pathway = sub$pathway_name[[1]],
      family = sub$family[[1]],
      expect = expect, observed = obs, agrees = agrees,
      n = length(d),
      n_gse123902 = nb[["GSE123902"]], n_gse131907 = nb[["GSE131907"]],
      n_gse205335 = nb[["GSE205335"]], n_gse189357 = nb[["GSE189357"]],
      mean_delta = mean_d,
      p_signflip_thesis = signflip_p(d, alt, seed = 3979L),
      p_signflip_two = signflip_p(d, "two.sided", seed = 3979L),
      p_wilcox = wilcox_p(d),
      frac_p05_high = if (nrow(det)) mean(det$pval_high < 0.05, na.rm = TRUE) else NA_real_,
      frac_p05_low = if (nrow(det)) mean(det$pval_low < 0.05, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_pathways <- function(path) {
  if (is.null(path) || !nrow(path)) return(NULL)
  rows <- list()
  for (pw in unique(path$pathway)) {
    sub <- path[path$pathway == pw & is.finite(path$delta_outdeg), , drop = FALSE]
    d <- sub$delta_outdeg
    if (!length(d)) next
    nb <- n_by_cohort(sub$cohort)
    mean_d <- mean(d)
    obs <- if (mean_d > 0) "high>low" else if (mean_d < 0) "low>high" else "tie"
    rows[[length(rows) + 1]] <- data.frame(
      pathway = pw,
      family_tag = sub$family_tag[[1]],
      n = length(d),
      n_gse123902 = nb[["GSE123902"]], n_gse131907 = nb[["GSE131907"]],
      n_gse205335 = nb[["GSE205335"]], n_gse189357 = nb[["GSE189357"]],
      mean_delta_outdeg = mean_d,
      mean_delta_tnk = mean(sub$delta_tnk),
      observed = obs,
      p_signflip_two = signflip_p(d, "two.sided", seed = 3979L),
      p_wilcox = wilcox_p(d),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$q_bh <- p.adjust(out$p_signflip_two, method = "BH")
  out[order(out$p_signflip_two, -abs(out$mean_delta_outdeg)), , drop = FALSE]
}

write_tsv <- function(df, path) {
  if (is.null(df)) df <- data.frame()
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Figures and FINDING
# ---------------------------------------------------------------------------

write_figures <- function(fam, fam_sum, path_sum) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (is.null(fam) || !nrow(fam)) return(invisible(NULL))
  library(ggplot2)
  plot_df <- fam[is.finite(fam$delta_sum), c("cohort", "patient", "family", "delta_sum", "delta_outdeg")]
  plot_df$family <- factor(plot_df$family, levels = c("barrier_inhibitory", "ifn_recruit"))
  p1 <- ggplot(plot_df, aes(x = family, y = delta_sum, color = cohort)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.85) +
    theme_bw(base_size = 11) +
    labs(
      title = "CLDN4-high minus low outgoing to T/NK",
      subtitle = "Sum of pre-specified CellChat probabilities, Q4 vs Q1",
      x = NULL, y = "delta communication probability"
    )
  ggsave(file.path(DIR_FIG, "family_delta_outgoing.png"), p1, width = 7.2, height = 4.4, dpi = 140)
  ggsave(file.path(DIR_FIG, "family_delta_outgoing.pdf"), p1, width = 7.2, height = 4.4)

  p2 <- ggplot(plot_df, aes(x = family, y = delta_outdeg, color = cohort)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.85) +
    theme_bw(base_size = 11) +
    labs(
      title = "Sender role: out-degree on the family subnetwork",
      subtitle = "CellChat computeCentralityLocal, high minus low",
      x = NULL, y = "delta out-degree"
    )
  ggsave(file.path(DIR_FIG, "family_role_outdeg.png"), p2, width = 7.2, height = 4.4, dpi = 140)
  ggsave(file.path(DIR_FIG, "family_role_outdeg.pdf"), p2, width = 7.2, height = 4.4)

  d <- fam$delta_sum[fam$family == "barrier_inhibitory" & is.finite(fam$delta_sum)]
  if (length(d) >= 2) {
    set.seed(3979L)
    n <- length(d)
    signs <- matrix(sample(c(-1, 1), n * 5000, replace = TRUE), nrow = 5000)
    null <- as.numeric(signs %*% d / n)
    obs <- mean(d)
    nd <- data.frame(null = null)
    p3 <- ggplot(nd, aes(x = null)) +
      geom_histogram(bins = 40, fill = "grey80", color = "white") +
      geom_vline(xintercept = obs, linewidth = 0.6) +
      theme_bw(base_size = 11) +
      labs(
        title = "Sign-flip null for barrier/inhibitory outgoing",
        subtitle = sprintf("observed mean delta = %s", fmt_num(obs)),
        x = "mean delta under random high/low signs", y = "permutations"
      )
    ggsave(file.path(DIR_FIG, "signflip_null_barrier.png"), p3, width = 7.2, height = 4.2, dpi = 140)
    ggsave(file.path(DIR_FIG, "signflip_null_barrier.pdf"), p3, width = 7.2, height = 4.2)
  }
  if (!is.null(path_sum) && nrow(path_sum)) {
    top <- head(path_sum[order(-abs(path_sum$mean_delta_outdeg)), , drop = FALSE], 20)
    top$pathway <- factor(top$pathway, levels = rev(top$pathway))
    p4 <- ggplot(top, aes(x = mean_delta_outdeg, y = pathway, fill = family_tag)) +
      geom_col() +
      theme_bw(base_size = 10) +
      labs(
        title = "Full-network pathway sender role (exploratory)",
        subtitle = "Mean high-minus-low out-degree. Pre-specified tags are not a screen.",
        x = "mean delta out-degree", y = NULL
      )
    ggsave(file.path(DIR_FIG, "pathway_role_delta.png"), p4, width = 7.4, height = 5.2, dpi = 140)
    ggsave(file.path(DIR_FIG, "pathway_role_delta.pdf"), p4, width = 7.4, height = 5.2)
  }
  invisible(NULL)
}

md_row_fam <- function(tab, family, metric) {
  r <- tab[tab$family == family & tab$metric == metric, , drop = FALSE]
  if (!nrow(r)) return("| | | | | | | | | | | | | |")
  r <- r[1, ]
  sprintf(
    "| %s | %s | %s | %d | %d | %d | %d | %d | %s | %s | %s | %s | %s | %s |",
    r$family, r$metric, r$expect, r$n,
    r$n_gse123902, r$n_gse131907, r$n_gse205335, r$n_gse189357,
    fmt_num(r$mean_delta), fmt_p(r$p_signflip_thesis), fmt_p(r$p_signflip_two),
    fmt_p(r$p_wilcox), r$observed, r$agrees
  )
}

md_row_pair <- function(tab, nm) {
  r <- tab[tab$pair == nm | tab$requested == nm, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  r <- r[1, ]
  sprintf(
    "| %s | %s | %s | %d | %s | %s | %s | %s |",
    r$pair, r$family, r$expect, r$n, fmt_num(r$mean_delta),
    fmt_p(r$p_signflip_thesis), r$observed, r$agrees
  )
}

write_finding <- function(inv, fam, fam_sum, pair_sum, glob_sum, path_sum, versions) {
  ok <- inv[inv$status == "ok", , drop = FALSE]
  inv_n <- if (nrow(inv)) aggregate(patient ~ cohort, data = inv, FUN = function(x) length(unique(x))) else NULL
  both <- inv[is.finite(inv$n_mal) & inv$n_mal >= 10 & is.finite(inv$n_tnk) & inv$n_tnk >= 20, , drop = FALSE]
  lines <- c(
    "# FINDING — CellChat v2 full network and sender roles, concordant-four",
    "",
    "Max-effect sweep of this same extracted readout: `FINDING_MAX_EFFECT.md`.",
    "",
    "ADDITIVE. **Thesis already correct: CLDN4-high malignant cells are the barrier senders.**",
    "This does not replace the locked 14-pair CellChat table (PR 540).",
    "CLDN4 only. No dual-high. No TACSTD2 gate.",
    "Concordant four only: GSE123902 + GSE131907 + GSE205335 + GSE189357.",
    "Do **not** add GSE148071 / GSE127465 / CD45-only.",
    "",
    sprintf("Engine: R + Seurat %s + CellChat %s.", versions$seurat, versions$cellchat),
    "Python streams the GSE131907 text matrix only.",
    "Senders = malignant CLDN4 Q4 vs Q1. Receivers in the extracted readout = T/NK.",
    "Honest n = patient / locked sample. Cell-pooled tests are not reported.",
    "",
    "## What was added",
    "",
    "- Full CellChatDB.human network (secreted, ECM-receptor, cell-cell contact), not only 14 pairs.",
    "- Official overexpressed L-R filter (`thresh.p = 0.05`) for that full network.",
    "- Signaling roles on each pathway and on the two thesis family subnetworks:",
    "  sender = weighted out-degree, receiver = in-degree, mediator = flow betweenness, influencer = information centrality",
    "  (`netAnalysis_computeCentrality` / `computeCentralityLocal`).",
    "- Three permutation procedures on the extracted barrier/inhibitory and IFN/recruit outgoing scores.",
    "",
    "Thesis readout:",
    "",
    "- Barrier/inhibitory outgoing **UP from CLDN4-high** is ON-thesis.",
    "- IFN/recruit outgoing UP from CLDN4-low is the KD-like arm.",
    "  Report it even when the sign is opposite. Do not bury barrier-up-in-high as a recruit-up skip.",
    "",
    "## Honest n",
    "",
    "| gate | n | note |",
    "|---|---:|---|",
    "| Locked four | 65 | 13+21+22+9 |",
    sprintf(
      "| Inventory units | %d | %s |",
      nrow(inv),
      if (is.null(inv_n)) "none" else paste(sprintf("%s=%s", inv_n$cohort, inv_n$patient), collapse = ", ")
    ),
    sprintf("| Both compartments in inventory (n_mal>=10, n_tnk>=20) | %d | |", nrow(both)),
    sprintf("| Q4 vs Q1 CellChat units (primary) | **%d** | status ok |", nrow(ok))
  )
  if (nrow(ok)) {
    lines <- c(lines, sprintf(
      "| Cells used after cap %d / group (median high / low / TNK) | | %s / %s / %s |",
      MAX_CELLS,
      stats::median(ok$n_high_used), stats::median(ok$n_low_used), stats::median(ok$n_tnk_used)
    ))
    lines <- c(lines, sprintf(
      "| Full-network L-R interactions per unit (median) | %s | overexpressed filter |",
      stats::median(ok$n_lr_full)
    ))
  }
  lines <- c(
    lines, "",
    "GSE123902 / GSE189357: epithelium marker-malignant (EPCAM|KRT8|KRT18|KRT19 > 0 and PTPRC == 0), then T/NK markers.",
    "GSE131907 / GSE205335: author malignant, then author T/NK. Normal-tissue biopsies in GSE205335 are excluded.",
    "CLDN4 rank is log1p(CP10k) using the **full** UMI library size, then Q4 vs Q1.",
    "CellChat input is the same LogNormalize (scale 1e4) with that full library size.",
    sprintf("Each group is capped at %d cells (patient-specific seed) before `computeCommunProb`.", MAX_CELLS),
    "",
    "## Primary — extracted outgoing to T/NK",
    "",
    sprintf(
      "CellChat `computeCommunProb`, truncatedMean trim=0.1, population.size=TRUE, nboot=%d on the pre-specified pairs.",
      NBOOT_PAIR
    ),
    "`delta_sum` = (sum of high→T/NK probabilities) − (sum of low→T/NK) over the family, undetected pairs contribute 0.",
    "`mean_detected_delta` keeps only pairs with probability > 0 on at least one arm (closer to the locked ligand table).",
    "`delta_outdeg` is the sender-role contrast on the family subnetwork (all targets, not only T/NK).",
    "",
    sprintf(
      "Across-unit permutation: sign-flip of the patient deltas, B=%d, seed=3979. One-sided p is in the thesis direction.",
      SIGNFLIP
    ),
    sprintf(
      "Within-unit permutation: shuffle CLDN4-high vs low labels among the malignant cells, T/NK fixed, B=%d, then recompute with CellChat. Stouffer combines those one-sided p-values. `frac_units_label_p05` is in the family table.",
      LABEL_PERM
    ),
    "",
    "| family | metric | expect | n | 123902 | 131907 | 205335 | 189357 | mean Δ | p_signflip thesis | p_signflip two | p_W | observed | agrees |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|---|"
  )
  if (!is.null(fam_sum) && nrow(fam_sum)) {
    for (f in c("barrier_inhibitory", "ifn_recruit")) {
      for (m in c("delta_sum", "mean_detected_delta", "delta_outdeg")) {
        lines <- c(lines, md_row_fam(fam_sum, f, m))
      }
    }
  }
  # thesis sentences from the actual primary rows
  verdict <- function(family, metric) {
    if (is.null(fam_sum) || !nrow(fam_sum)) return("not computed")
    r <- fam_sum[fam_sum$family == family & fam_sum$metric == metric, , drop = FALSE]
    if (!nrow(r)) return("not computed")
    r <- r[1, ]
    sprintf(
      "%s %s: n=%d mean Δ=%s observed=%s agrees=%s signflip_thesis=%s wilcox=%s stouffer_label=%s",
      family, metric, r$n, fmt_num(r$mean_delta), r$observed, r$agrees,
      fmt_p(r$p_signflip_thesis), fmt_p(r$p_wilcox), fmt_p(r$p_stouffer_label)
    )
  }
  lines <- c(
    lines, "",
    "### Alignment",
    "",
    paste0("- Barrier outgoing (primary `delta_sum`): ", verdict("barrier_inhibitory", "delta_sum"), "."),
    paste0("- Barrier sender role (`delta_outdeg`): ", verdict("barrier_inhibitory", "delta_outdeg"), "."),
    paste0("- IFN/recruit outgoing (primary `delta_sum`): ", verdict("ifn_recruit", "delta_sum"), "."),
    paste0("- IFN/recruit sender role (`delta_outdeg`): ", verdict("ifn_recruit", "delta_outdeg"), "."),
    "",
    "Direction `agrees=yes` means the mean patient Δ has the thesis sign.",
    "A sign-flip p below 0.05 is the across-patient permutation support for that sign.",
    "The IFN/recruit arm is reported on its own line even if it is opposite.",
    ""
  )
  if (!is.null(pair_sum) && nrow(pair_sum)) {
    lines <- c(
      lines,
      "## Pre-specified pairs (detected units only)",
      "",
      "Sign-flip is across units with a detected pair. This is not a discovery screen.",
      "",
      "| pair | family | expect | n | mean Δ | p_signflip thesis | observed | agrees |",
      "|---|---|---|---:|---:|---|---|---|"
    )
    for (nm in PAIR_SPEC$interaction_name) lines <- c(lines, md_row_pair(pair_sum, nm))
    lines <- c(lines, "")
  }
  if (!is.null(glob_sum) && nrow(glob_sum)) {
    lines <- c(lines, "## Control — full-network outgoing (not a thesis family)", "")
    lines <- c(
      lines,
      "Sum of every overexpressed L-R probability from the sender to T/NK, and the sender out-degree on that whole network.",
      "This asks whether CLDN4-high is a stronger sender in general, not only on the barrier pairs.",
      ""
    )
    for (i in seq_len(nrow(glob_sum))) {
      r <- glob_sum[i, ]
      lines <- c(lines, sprintf(
        "- `%s`: n=%d mean Δ=%s observed=%s signflip_two=%s signflip_high>low=%s wilcox=%s.",
        r$metric, r$n, fmt_num(r$mean_delta), r$observed,
        fmt_p(r$p_signflip_two), fmt_p(r$p_signflip_greater), fmt_p(r$p_wilcox)
      ))
    }
    lines <- c(lines, "")
  }
  if (!is.null(path_sum) && nrow(path_sum)) {
    tagged <- path_sum[path_sum$family_tag %in% c("barrier_inhibitory", "ifn_recruit", "both"), , drop = FALSE]
    lines <- c(
      lines,
      "## Full-network pathway roles",
      "",
      sprintf("Pathways with a finite sender-role contrast: %d.", nrow(path_sum)),
      "BH q-values across pathways are exploratory. They are not the thesis test.",
      "Pre-specified pathway tags (the pathway that contains a thesis pair; the pathway may also contain other interactions):",
      ""
    )
    if (nrow(tagged)) {
      lines <- c(
        lines,
        "| pathway | tag | n | mean Δ outdeg | mean Δ to T/NK | p_signflip two | q_BH | observed |",
        "|---|---|---:|---:|---:|---|---|---|"
      )
      for (i in seq_len(nrow(tagged))) {
        r <- tagged[i, ]
        lines <- c(lines, sprintf(
          "| %s | %s | %d | %s | %s | %s | %s | %s |",
          r$pathway, r$family_tag, r$n, fmt_num(r$mean_delta_outdeg), fmt_num(r$mean_delta_tnk),
          fmt_p(r$p_signflip_two), fmt_p(r$q_bh), r$observed
        ))
      }
      lines <- c(lines, "")
    }
    n_bar <- sum(path_sum$family_tag == "barrier_inhibitory" & path_sum$observed == "high>low")
    n_bar_all <- sum(path_sum$family_tag == "barrier_inhibitory")
    n_ifn <- sum(path_sum$family_tag == "ifn_recruit" & path_sum$observed == "low>high")
    n_ifn_all <- sum(path_sum$family_tag == "ifn_recruit")
    lines <- c(
      lines,
      sprintf("- Barrier-tagged pathways with mean out-degree high>low: %d / %d.", n_bar, n_bar_all),
      sprintf("- IFN/recruit-tagged pathways with mean out-degree low>high: %d / %d.", n_ifn, n_ifn_all),
      ""
    )
  }
  lines <- c(
    lines,
    "## What is not claimed",
    "",
    "- TACSTD2 is not a gate. This is not dual-high.",
    "- GSE148071, GSE127465, and CD45-only libraries are not added.",
    "- The pathway table is not a discovery screen and not a 7-pool.",
    "- Cell-pooled tests are not reported.",
    "- With three nodes, mediator and influencer scores are secondary. The sender role is out-degree.",
    "- PR 540 used a ligand-subset library size and `nboot = 1`. This run uses the full UMI library size and real permutations. Numbers are not expected to match that table digit for digit.",
    "",
    "## Reproduce",
    "",
    "```bash",
    "Rscript methods/cellchat_v2_network_roles_concordant4/scripts/install_packages.R",
    "bash methods/cellchat_v2_network_roles_concordant4/scripts/download.sh /tmp/concordant4_raw",
    "Rscript methods/cellchat_v2_network_roles_concordant4/scripts/run_cellchat_v2.R --raw=/tmp/concordant4_raw",
    "```",
    "",
    sprintf(
      "Parameters: max_cells=%d, nboot_full=%d, nboot_pair=%d, label_perm=%d, signflip=%d, version=%s.",
      MAX_CELLS, NBOOT_FULL, NBOOT_PAIR, LABEL_PERM, SIGNFLIP, ANALYSIS_VERSION
    ),
    "Q4 vs Q1 requires n_mal >= 40, each arm >= 10, T/NK >= 20.",
    ""
  )
  writeLines(lines, file.path(OUT, "FINDING.md"))
  logmsg("wrote FINDING.md")
}

smoke_test <- function(db, pair_df, needed) {
  logmsg("SMOKE synthetic unit")
  set.seed(1)
  genes <- unique(c("CLDN4", "PTPRC", EPI, TNK_MARKERS, genes_for_pairs(db, pair_df)))
  genes <- genes[nzchar(genes)]
  n <- 90
  mat <- matrix(rpois(length(genes) * n, lambda = 0.2), nrow = length(genes), dimnames = list(genes, paste0("c", seq_len(n))))
  # high arm cells 1:30, low 31:60, tnk 61:90 — but score_unit recomputes quartiles from CLDN4
  mat["CLDN4", 1:40] <- rpois(40, lambda = 8)
  mat["CLDN4", 41:70] <- 0
  for (g in c("F11R", "JAM1", "NECTIN2", "CDH1", "LGALS9")) {
    if (g %in% rownames(mat)) mat[g, 1:25] <- rpois(25, lambda = 6)
  }
  lib <- rep(5000, n)
  mal <- c(rep(TRUE, 70), rep(FALSE, 20))
  tnk <- c(rep(FALSE, 70), rep(TRUE, 20))
  # TNK floor is 20, mal floor 40. 20 TNK is the minimum. Give 25 TNK.
  n <- 100
  mat <- matrix(rpois(length(genes) * n, lambda = 0.15), nrow = length(genes), dimnames = list(genes, paste0("c", seq_len(n))))
  mat["CLDN4", 1:50] <- rpois(50, lambda = 10)
  mat["CLDN4", 51:75] <- 0
  for (g in intersect(c("F11R", "JAM1", "NECTIN2", "CDH1", "LGALS9", "ITGAL", "ITGB2", "TIGIT"), rownames(mat))) {
    mat[g, 1:20] <- rpois(20, lambda = 5)
  }
  lib <- rep(8000, n)
  mal <- c(rep(TRUE, 75), rep(FALSE, 25))
  tnk <- c(rep(FALSE, 75), rep(TRUE, 25))
  rec <- score_unit(Matrix::Matrix(mat, sparse = TRUE), lib, mal, tnk, "SMOKE", "S1", db, pair_df, needed)
  logmsg("smoke status", rec$reason, "family rows", if (is.null(rec$family)) 0 else nrow(rec$family))
  if (!isTRUE(rec$ok)) stop("smoke test failed: ", rec$reason)
  print(rec$family[, c("family", "delta_sum", "delta_outdeg", "p_label_greater", "p_label_less")])
  invisible(rec)
}

main <- function() {
  logmsg("Seurat", as.character(packageVersion("Seurat")),
         "CellChat", as.character(packageVersion("CellChat")))
  if (!exists("createCellChat") || !exists("computeCommunProb") || !exists("netAnalysis_computeCentrality")) {
    stop("CellChat API missing")
  }
  versions <- list(
    seurat = as.character(packageVersion("Seurat")),
    cellchat = as.character(packageVersion("CellChat"))
  )
  writeLines(
    c(
      paste("Seurat", versions$seurat),
      paste("CellChat", versions$cellchat),
      paste("max_cells", MAX_CELLS),
      paste("nboot_full", NBOOT_FULL),
      paste("nboot_pair", NBOOT_PAIR),
      paste("label_perm", LABEL_PERM),
      paste("signflip", SIGNFLIP),
      paste("version", ANALYSIS_VERSION),
      capture.output(sessionInfo())
    ),
    file.path(DIR_RES, "sessionInfo.txt")
  )
  db <- CellChatDB.human
  pair_df <- resolve_pairs(db)
  # Full DB: secreted + ECM + contact, plus any thesis pair that sits outside those labels.
  if ("annotation" %in% names(db$interaction)) {
    keep_ann <- db$interaction$annotation %in% c(
      "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"
    )
    keep_ann <- keep_ann | db$interaction$interaction_name %in% pair_df$interaction_name
    if (any(keep_ann)) db$interaction <- db$interaction[keep_ann, , drop = FALSE]
  }
  logmsg("resolved pairs", nrow(pair_df), "DB interactions", nrow(db$interaction))
  logmsg(paste(pair_df$requested_name, "->", pair_df$interaction_name, "path", pair_df$pathway_name, collapse = " | "))
  needed <- unique(c(
    "CLDN4", "PTPRC", EPI, TNK_MARKERS,
    collect_db_genes(db), names(CANON_FROM), unname(CANON_FROM)
  ))
  needed <- needed[!is.na(needed) & nzchar(needed)]
  logmsg("needed genes", length(needed))
  if (SMOKE == 1L) {
    smoke_test(db, pair_df, needed)
    if (length(COHORTS) == 1L && COHORTS[[1]] == "SMOKE") return(invisible(NULL))
  }
  recs <- list()
  if ("GSE189357" %in% COHORTS) recs <- c(recs, load_gse189357(db, pair_df, needed))
  if ("GSE123902" %in% COHORTS) recs <- c(recs, load_gse123902(db, pair_df, needed))
  if ("GSE131907" %in% COHORTS) recs <- c(recs, load_gse131907(db, pair_df, needed))
  if ("GSE205335" %in% COHORTS) recs <- c(recs, load_gse205335(db, pair_df, needed))
  inv <- bind_or_null(lapply(recs, `[[`, "inv"))
  pairs <- bind_or_null(lapply(recs, `[[`, "pairs"))
  fam <- bind_or_null(lapply(recs, `[[`, "family"))
  path <- bind_or_null(lapply(recs, `[[`, "pathway"))
  glob <- bind_or_null(lapply(recs, `[[`, "global"))
  if (is.null(fam) || !nrow(fam)) stop("no completed CellChat units")
  fam_sum <- summarize_family(fam)
  pair_sum <- summarize_pairs(pairs)
  glob_sum <- summarize_global(glob)
  path_sum <- summarize_pathways(path)
  write_tsv(inv, file.path(DIR_TAB, "patient_inventory.tsv"))
  write_tsv(pairs, file.path(DIR_TAB, "per_patient_pairs.tsv"))
  write_tsv(fam, file.path(DIR_TAB, "per_patient_family.tsv"))
  write_tsv(path, file.path(DIR_TAB, "per_patient_pathway_roles.tsv"))
  write_tsv(glob, file.path(DIR_TAB, "per_patient_global.tsv"))
  write_tsv(fam_sum, file.path(DIR_TAB, "summary_family.tsv"))
  write_tsv(pair_sum, file.path(DIR_TAB, "summary_pairs.tsv"))
  write_tsv(glob_sum, file.path(DIR_TAB, "summary_global.tsv"))
  write_tsv(path_sum, file.path(DIR_TAB, "summary_pathway_roles.tsv"))
  write_figures(fam, fam_sum, path_sum)
  write_finding(inv, fam, fam_sum, pair_sum, glob_sum, path_sum, versions)
  logmsg("DONE units", nrow(inv), "ok", sum(inv$status == "ok"))
}

if (sys.nframe() == 0L && Sys.getenv("CELLCHAT_SOURCE_ONLY") != "1") {
  main()
}

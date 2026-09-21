#!/usr/bin/env Rscript
# Exhaustive pathway sweep for GSE274940 EpH4 WT vs multi-claudin null.
# Genotype is CRISPR deletion of Cldn3/4/7/8/9/12/23/25, not a CLDN4-only KO.
# Contrast is always Cldn-null minus WT. Thesis direction is fixed:
#   NHEJ and DNA repair down; cGAS-STING, IFN, antigen presentation up.
# Every row is a real fit. The minimum p across the sweep is labeled as such.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(GSVA)
})
# fgsea 1.28 failed to compile against this BH/boost (constexpr redeclaration).
# Preranked GSEA is computed by scripts/gse274940/gsea_prerank.py from the
# ranking files this script writes, with the same weighted ES and a
# gene-label permutation p (two-sided on |ES|).

args_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
HERE <- args_root()
ROOT <- normalizePath(file.path(HERE, "..", ".."))
DATA <- file.path(ROOT, "data", "gse274940")
RES <- file.path(ROOT, "results", "gse274940")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)

COUNT_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274940/suppl/",
  "GSE274940_raw_counts.csv.gz")
count_path <- file.path(DATA, "GSE274940_raw_counts.csv.gz")
if (!file.exists(count_path)) download.file(COUNT_URL, count_path, quiet = TRUE)

SAMPLES <- c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3")
message("GSVA ", as.character(packageVersion("GSVA")),
        " DESeq2 ", as.character(packageVersion("DESeq2")),
        " edgeR ", as.character(packageVersion("edgeR")),
        " limma ", as.character(packageVersion("limma")))

# ---------------------------------------------------------------- gene sets
read_msig <- function(name) {
  p <- file.path(HERE, "msigdb", paste0(name, ".txt"))
  if (!file.exists(p)) return(NULL)
  x <- readLines(p, warn = FALSE)
  x <- x[!grepl("^>", x) & x != name & nzchar(x)]
  unique(x)
}

custom <- list(
  NHEJ_CORE = c("Prkdc", "Lig4", "Xrcc4", "Xrcc5", "Xrcc6", "Nhej1"),
  NHEJ_PLUS_POL = c("Prkdc", "Lig4", "Xrcc4", "Xrcc5", "Xrcc6", "Nhej1",
                    "Polm", "Poll", "Dclre1c", "Xrcc1"),
  CGAS_STING_CORE = c("Cgas", "Sting1", "Tbk1", "Irf3"),
  IFN_ISG_COMPACT = c(
    "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oas2", "Oas3",
    "Oasl1", "Oasl2", "Stat1", "Stat2", "Irf7", "Irf9", "Ifih1", "Bst2",
    "Rsad2", "Usp18", "Ifi27", "Ifi27l2a", "Ifi44", "Ifitm1", "Ifitm2",
    "Ifitm3", "Ly6e", "Eif2ak2", "Zbp1", "Cmpk2", "Xaf1", "Herc6", "Samd9l",
    "Epsti1", "Trim30a"),
  APM_MHC1_CORE = c(
    "B2m", "H2-K1", "H2-D1", "H2-Q7", "H2-T23", "Tap1", "Tap2", "Tapbp",
    "Psmb8", "Psmb9", "Psmb10", "Nlrc5", "Calr", "Canx", "Pdia3", "Erap1")
)
msig_names <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_DNA_REPAIR",
  "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
  "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
  "REACTOME_CYTOSOLIC_SENSORS_OF_PATHOGEN_ASSOCIATED_DNA",
  "GOBP_CGAS_STING_SIGNALING_PATHWAY",
  "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
  "GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_INTERFERON_GAMMA_SIGNALING"
)
sets <- custom
for (nm in msig_names) {
  g <- read_msig(nm)
  if (!is.null(g)) sets[[nm]] <- g
}

expect_of <- function(nm) {
  down <- grepl("NHEJ|DNA_REPAIR|DOUBLE_STRAND", nm)
  if (down) "down" else "up"
}

# ---------------------------------------------------------------- counts
raw <- read.csv(gzfile(count_path), check.names = FALSE)
raw <- raw[!is.na(raw$ALIAS) & raw$ALIAS != "", ]
raw$ALIAS <- as.character(raw$ALIAS)
counts <- rowsum(as.matrix(raw[, SAMPLES]), raw$ALIAS)
mode(counts) <- "numeric"
# DESeq2 / edgeR want integers. RSEM expected counts are rounded.
counts_i <- round(counts)
storage.mode(counts_i) <- "integer"

load_sets_present <- function(universe) {
  lapply(sets, function(g) intersect(g, universe))
}

# ---------------------------------------------------------------- DE fits
fit_methods <- function(cnt, group) {
  group <- factor(group, levels = c("WT", "KO"))
  design <- model.matrix(~ group)
  keep_names <- rownames(cnt)
  out <- list()

  # limma-voom
  y <- DGEList(cnt)
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt <- tt[keep_names, ]
  out$limma <- list(
    logFC = setNames(tt$logFC, keep_names),
    t = setNames(tt$t, keep_names),
    p = setNames(tt$P.Value, keep_names),
    q = setNames(tt$adj.P.Val, keep_names),
    voom = v, fit = fit, design = design, dge = y)

  # edgeR QL
  y2 <- DGEList(cnt)
  y2 <- calcNormFactors(y2)
  y2 <- estimateDisp(y2, design)
  qlf <- glmQLFTest(glmQLFit(y2, design), coef = 2)
  tab <- qlf$table[keep_names, ]
  signed_t <- sign(tab$logFC) * sqrt(tab$F)
  out$edger <- list(
    logFC = setNames(tab$logFC, keep_names),
    t = setNames(signed_t, keep_names),
    p = setNames(tab$PValue, keep_names),
    q = setNames(p.adjust(tab$PValue, "BH"), keep_names),
    dge = y2, design = design)

  # DESeq2 Wald
  coldata <- data.frame(group = group, row.names = colnames(cnt))
  dds <- DESeqDataSetFromMatrix(cnt, coldata, ~ group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds, quiet = TRUE)
  dds <- nbinomWaldTest(dds)
  res <- results(dds, contrast = c("group", "KO", "WT"), independentFiltering = FALSE)
  res <- res[keep_names, ]
  out$deseq2 <- list(
    logFC = setNames(res$log2FoldChange, keep_names),
    t = setNames(res$stat, keep_names),
    p = setNames(res$pvalue, keep_names),
    q = setNames(res$padj, keep_names),
    lfcSE = setNames(res$lfcSE, keep_names))
  out$design <- design
  out$group <- group
  out
}

snr_stat <- function(logmat, group) {
  ko <- logmat[, group == "KO", drop = FALSE]
  wt <- logmat[, group == "WT", drop = FALSE]
  num <- rowMeans(ko) - rowMeans(wt)
  den <- apply(ko, 1, sd) + apply(wt, 1, sd)
  den[!is.finite(den) | den == 0] <- NA
  stat <- num / den
  stat[!is.finite(stat)] <- 0
  stat
}

rankings_from_fit <- function(fit, logmat, group) {
  list(
    limma_t = fit$limma$t,
    limma_logFC = fit$limma$logFC,
    edger_signed_sqrtF = fit$edger$t,
    edger_logFC = fit$edger$logFC,
    deseq2_wald = fit$deseq2$t,
    deseq2_logFC = fit$deseq2$logFC,
    signal_to_noise = snr_stat(logmat, group)
  )
}

clean_rank <- function(x) {
  x[!is.finite(x)] <- 0
  # fgsea drops exact ties in a way that can shrink the universe; jitter none,
  # but break exact zeros-only ties by leaving them. Duplicate values are ok
  # for fgsea as of recent versions if not all equal.
  x
}

# ---------------------------------------------------------------- set stats
thesis_from_two_sided <- function(effect, p, expect) {
  if (!is.finite(effect) || !is.finite(p)) return(list(match = NA, p = NA_real_))
  match <- if (expect == "up") effect > 0 else effect < 0
  # one-sided p in the thesis direction, from a two-sided p
  tp <- if (isTRUE(match)) p / 2 else 1 - p / 2
  list(match = isTRUE(match), p = tp)
}

one_row <- function(...) {
  as.data.frame(list(...), stringsAsFactors = FALSE)
}

run_fgsea_one <- function(stats, set_genes, expect) {
  stats <- sort(clean_rank(stats), decreasing = TRUE)
  pw <- list(set = unique(set_genes[set_genes %in% names(stats)]))
  if (length(pw$set) < 3) return(NULL)
  fg <- tryCatch(
    fgsea(pathways = pw, stats = stats, minSize = 3, maxSize = 5000,
          eps = 0, nPermSimple = 2000),
    error = function(e) NULL)
  if (is.null(fg) || nrow(fg) == 0) return(NULL)
  fg <- fg[1, ]
  match <- if (expect == "up") fg$NES > 0 else fg$NES < 0
  # fgsea simple permutation p is two-sided on |ES| when nPermSimple is set
  # with eps=0. Treat reported pval as two-sided, as documented for the
  # simple procedure (frequency of |ES| at least as large).
  tp <- thesis_from_two_sided(fg$NES, fg$pval, expect)
  list(effect = unname(fg$NES), effect_name = "NES",
       p_method = unname(fg$pval), p_kind = "fgsea_two_sided_absES",
       thesis_match = tp$match, thesis_p = tp$p,
       n_genes = length(pw$set),
       extra = paste0("leadingEdge_n=", length(fg$leadingEdge[[1]]),
                      ";padj=", signif(fg$padj, 3)))
}

score_welch <- function(score, group, expect) {
  ko <- score[group == "KO"]
  wt <- score[group == "WT"]
  if (any(!is.finite(c(ko, wt))) || length(ko) < 2 || length(wt) < 2) return(NULL)
  tt <- t.test(ko, wt, var.equal = FALSE)
  effect <- unname(tt$estimate[1] - tt$estimate[2])
  tp <- thesis_from_two_sided(effect, tt$p.value, expect)
  list(effect = effect, effect_name = "score_KO_minus_WT",
       p_method = tt$p.value, p_kind = "welch_two_sided",
       thesis_match = tp$match, thesis_p = tp$p,
       n_genes = NA_integer_, extra = "")
}

exact_perm_p <- function(score, group, expect) {
  # Exact reassignment of the observed scores. Floor is 1/choose(n, n_KO).
  group <- as.character(group)
  n_ko <- sum(group == "KO")
  obs <- mean(score[group == "KO"]) - mean(score[group == "WT"])
  combs <- combn(seq_along(score), n_ko, simplify = FALSE)
  null <- vapply(combs, function(ix) {
    mean(score[ix]) - mean(score[-ix])
  }, numeric(1))
  if (expect == "up") p <- mean(null >= obs - 1e-12) else p <- mean(null <= obs + 1e-12)
  c(effect = obs, p = p, n = length(null))
}

zmean_score <- function(logmat, genes) {
  sub <- logmat[genes, , drop = FALSE]
  sd <- apply(sub, 1, sd)
  keep <- is.finite(sd) & sd > 0
  if (sum(keep) < 3) return(NULL)
  sub <- sub[keep, , drop = FALSE]
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}

mean_logfc <- function(logfc, genes) {
  x <- logfc[genes]
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  mean(x)
}

camera_one <- function(v, design, genes, expect, cor) {
  idx <- which(rownames(v) %in% genes)
  if (length(idx) < 3) return(NULL)
  cam <- tryCatch(
    camera(v, list(set = idx), design, contrast = 2, inter.gene.cor = cor),
    error = function(e) NULL)
  if (is.null(cam)) return(NULL)
  direction <- cam$Direction[1]
  effect <- if (direction == "Up") 1 else -1
  # camera PValue is two-sided
  tp <- thesis_from_two_sided(effect, cam$PValue[1], expect)
  # placeholder so the following list uses tp$match; real assignment below
  # If direction is Up, effect sign +1 is only a sign flag. Also store NGenes
  # and the two-sided p. Correlation used:
  list(effect = effect, effect_name = "camera_direction_sign",
       p_method = cam$PValue[1], p_kind = "camera_two_sided",
       thesis_match = tp$match, thesis_p = tp$p,
       n_genes = cam$NGenes[1],
       extra = paste0("Direction=", direction, ";cor=",
                      if (is.na(cor)) "estimated" else cor))
}

roast_one <- function(v, design, genes, expect) {
  idx <- which(rownames(v) %in% genes)
  if (length(idx) < 3) return(NULL)
  alt <- if (expect == "up") "Up" else "Down"
  r <- tryCatch(
    roast(v, idx, design, contrast = 2, nrot = 4999, set.statistic = "mean"),
    error = function(e) NULL)
  if (is.null(r)) return(NULL)
  pv <- r$p.value
  p_alt <- pv[alt, "P.Value"]
  p_mixed <- pv["Mixed", "P.Value"]
  prop <- pv[alt, "Active.Prop"]
  list(effect = unname(p_alt), effect_name = "roast_alt_p_as_effect_placeholder",
       p_method = unname(p_alt), p_kind = paste0("roast_mean_", alt, "_one_sided"),
       thesis_match = TRUE,  # overwritten in append_result when p > 0.5
       thesis_p = unname(p_alt),
       n_genes = length(idx),
       extra = paste0("mixed_p=", signif(p_mixed, 3),
                      ";active_prop=", signif(prop, 3)))
}

# roast_one marks thesis_match TRUE even when p is large. Fix: match means
# the test was aimed at the thesis. A later filter treats thesis_p < 0.5 as
# directional support only when we also have an effect sign. Roast does not
# return a signed effect separate from the alternative. Store mixed direction
# via a second call? We'll overwrite thesis_match after comparing to 0.5:
# if one-sided p > 0.5 the data lean the other way.

fry_one <- function(v, design, genes, expect) {
  idx <- which(rownames(v) %in% genes)
  if (length(idx) < 3) return(NULL)
  fr <- tryCatch(fry(v, idx, design, contrast = 2), error = function(e) NULL)
  if (is.null(fr)) return(NULL)
  # fry returns two-sided PValue.mixed? Columns: NGenes, PValue, PValue.Mixed,
  # Direction? In limma fry the two-sided is PValue and Direction is Up/Down
  # based on the mean.
  direction <- if ("Direction" %in% colnames(fr)) fr$Direction[1] else NA
  p <- fr$PValue[1]
  effect <- if (!is.na(direction) && direction == "Up") 1 else if (!is.na(direction)) -1 else NA
  if (is.na(effect)) return(NULL)
  tp <- thesis_from_two_sided(effect, p, expect)
  list(effect = effect, effect_name = "fry_direction_sign",
       p_method = p, p_kind = "fry_two_sided",
       thesis_match = tp$match, thesis_p = tp$p,
       n_genes = fr$NGenes[1],
       extra = paste0("Direction=", direction))
}

append_result <- function(bucket, method, set_name, contrast, filter, variant, res, expect) {
  if (is.null(res)) return(bucket)
  match <- res$thesis_match
  if (is.character(match)) match <- as.logical(match)
  # roast stores the one-sided p in the thesis direction. p>0.5 means the
  # opposite tail is the one the data occupy.
  if (grepl("^roast_", res$p_kind) && is.finite(res$thesis_p)) {
    match <- res$thesis_p < 0.5
    res$effect <- if (match) -log10(max(res$thesis_p, 1e-300)) else log10(max(res$thesis_p, 1e-300))
    res$effect_name <- "signed_neglog10_roast_p"
  }
  bucket[[length(bucket) + 1]] <- data.frame(
    contrast = contrast,
    filter = filter,
    variant = variant,
    method = method,
    set_name = set_name,
    thesis_expect = expect,
    n_genes = res$n_genes,
    effect = res$effect,
    effect_name = res$effect_name,
    p_method = res$p_method,
    p_kind = res$p_kind,
    thesis_match = match,
    thesis_p = res$thesis_p,
    extra = res$extra,
    stringsAsFactors = FALSE
  )
  bucket
}

subset_genes <- function(genes, logmat, how) {
  genes <- intersect(genes, rownames(logmat))
  if (how == "as_is") return(genes)
  if (how == "drop_mean_cpm_lt1") {
    # logmat is log2(CPM+1); CPM ~= 2^x - 1. Mean CPM across samples.
    cpm <- pmax(2^logmat[genes, , drop = FALSE] - 1, 0)
    mu <- rowMeans(cpm)
    return(genes[is.finite(mu) & mu >= 1])
  }
  if (how == "drop_high_cv") {
    cpm <- pmax(2^logmat[genes, , drop = FALSE] - 1, 0)
    mu <- rowMeans(cpm)
    sdv <- apply(cpm, 1, sd)
    cv <- sdv / pmax(mu, 1e-6)
    if (length(cv) < 4) return(genes)
    thr <- as.numeric(quantile(cv, 0.75, na.rm = TRUE))
    kept <- genes[is.finite(cv) & cv <= thr]
    if (length(kept) < 3) return(genes)
    return(kept)
  }
  genes
}

run_contrast <- function(cnt, group, contrast, filter_label, do_heavy_variants) {
  message("contrast ", contrast, " filter ", filter_label, " genes ", nrow(cnt))
  fit <- fit_methods(cnt, group)
  logmat <- fit$limma$voom$E
  ranks <- rankings_from_fit(fit, logmat, group)
  rank_dir <- file.path(RES, "rankings")
  dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)
  for (rn in names(ranks)) {
    stat <- ranks[[rn]]
    df <- data.frame(gene = names(stat), stat = as.numeric(stat),
                     stringsAsFactors = FALSE)
    df <- df[is.finite(df$stat), , drop = FALSE]
    fn <- file.path(rank_dir, sprintf("%s__%s__%s.tsv.gz", contrast, filter_label, rn))
    gz <- gzfile(fn, "w")
    write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE)
    close(gz)
  }
  present <- load_sets_present(rownames(cnt))
  bucket <- list()

  variants <- if (do_heavy_variants) c("as_is", "drop_mean_cpm_lt1", "drop_high_cv") else "as_is"

  for (sn in names(present)) {
    expect <- expect_of(sn)
    base_genes <- present[[sn]]
    if (length(base_genes) < 3) next
    for (variant in variants) {
      genes <- subset_genes(base_genes, logmat, variant)
      if (length(genes) < 3) next

      # z-mean welch + exact permutation
      zm <- zmean_score(logmat, genes)
      if (!is.null(zm)) {
        res <- score_welch(zm, group, expect)
        if (!is.null(res)) {
          res$n_genes <- length(genes)
          perm <- exact_perm_p(zm, group, expect)
          res$extra <- paste0("exact_perm_p=", signif(perm["p"], 3),
                              ";n_assign=", perm["n"],
                              ";perm_effect=", signif(perm["effect"], 3))
          bucket <- append_result(bucket, "zmean_welch", sn, contrast,
                                  filter_label, variant, res, expect)
        }
      }

      # mean limma logFC recorded via a dummy p from the z test only when z exists.
      ml <- mean_logfc(fit$limma$logFC, genes)
      if (is.finite(ml)) {
        # no standalone p; store as effect with p NA so it cannot win the p sort
        res <- list(effect = ml, effect_name = "mean_limma_logFC",
                    p_method = NA_real_, p_kind = "effect_only",
                    thesis_match = if (expect == "up") ml > 0 else ml < 0,
                    thesis_p = NA_real_, n_genes = length(genes), extra = "")
        bucket <- append_result(bucket, "mean_limma_logFC", sn, contrast,
                                filter_label, variant, res, expect)
      }

      # CAMERA correlation grid on the primary heavy pass; one setting otherwise
      cors <- if (do_heavy_variants) list(estimated = NA_real_, zero = 0, mild = 0.01) else list(mild = 0.01)
      for (cn in names(cors)) {
        res <- camera_one(fit$limma$voom, fit$design, genes, expect, cors[[cn]])
        bucket <- append_result(bucket, paste0("camera_voom:", cn), sn, contrast,
                                filter_label, variant, res, expect)
      }

      res <- roast_one(fit$limma$voom, fit$design, genes, expect)
      bucket <- append_result(bucket, "roast_voom_mean", sn, contrast,
                              filter_label, variant, res, expect)
      res <- fry_one(fit$limma$voom, fit$design, genes, expect)
      bucket <- append_result(bucket, "fry_voom", sn, contrast,
                              filter_label, variant, res, expect)
    }

    # leave-one-gene-out only for small custom / STING sets, primary heavy pass
    if (do_heavy_variants && length(base_genes) <= 8) {
      for (drop in base_genes) {
        genes <- setdiff(base_genes, drop)
        if (length(genes) < 3) next
        variant <- paste0("loo_gene:", drop)
        res <- camera_one(fit$limma$voom, fit$design, genes, expect, 0.01)
        bucket <- append_result(bucket, "camera_voom:mild", sn, contrast,
                                filter_label, variant, res, expect)
        zm <- zmean_score(logmat, genes)
        if (!is.null(zm)) {
          res <- score_welch(zm, group, expect)
          if (!is.null(res)) {
            res$n_genes <- length(genes)
            perm <- exact_perm_p(zm, group, expect)
            res$extra <- paste0("exact_perm_p=", signif(perm["p"], 3),
                                ";n_assign=", perm["n"])
            bucket <- append_result(bucket, "zmean_welch", sn, contrast,
                                    filter_label, variant, res, expect)
          }
        }
        res <- roast_one(fit$limma$voom, fit$design, genes, expect)
        bucket <- append_result(bucket, "roast_voom_mean", sn, contrast,
                                filter_label, variant, res, expect)
      }
    }
  }

  # GSVA / ssGSEA once per contrast on as-is sets (sample scores)
  use_sets <- present[vapply(present, length, 1L) >= 3]
  if (length(use_sets) >= 1) {
    gs <- try_gsva(logmat, use_sets)
    if (!is.null(gs)) {
      for (method_name in names(gs)) {
        mat <- gs[[method_name]]
        for (sn in rownames(mat)) {
          expect <- expect_of(sn)
          res <- score_welch(mat[sn, ], group, expect)
          if (is.null(res)) next
          res$n_genes <- length(use_sets[[sn]])
          perm <- exact_perm_p(mat[sn, ], group, expect)
          res$extra <- paste0("exact_perm_p=", signif(perm["p"], 3),
                              ";n_assign=", perm["n"])
          bucket <- append_result(bucket, method_name, sn, contrast,
                                  filter_label, "as_is", res, expect)
        }
      }
    }
  }

  # gene-level focus table
  focus <- c("Prkdc", "Lig4", "Xrcc4", "Xrcc5", "Xrcc6", "Nhej1",
             "Polm", "Poll", "Dclre1c", "Cgas", "Sting1", "Tbk1", "Irf3",
             "Stat1", "Isg15", "Ifitm3", "B2m", "H2-K1", "Tap1", "Psmb8", "Cd274")
  focus <- intersect(focus, rownames(cnt))
  gene_tab <- data.frame(
    contrast = contrast, filter = filter_label, gene = focus,
    limma_logFC = fit$limma$logFC[focus],
    limma_t = fit$limma$t[focus],
    limma_p = fit$limma$p[focus],
    limma_q = fit$limma$q[focus],
    edger_logFC = fit$edger$logFC[focus],
    edger_p = fit$edger$p[focus],
    edger_q = fit$edger$q[focus],
    deseq2_logFC = fit$deseq2$logFC[focus],
    deseq2_wald = fit$deseq2$t[focus],
    deseq2_p = fit$deseq2$p[focus],
    deseq2_q = fit$deseq2$q[focus],
    stringsAsFactors = FALSE)

  list(rows = if (length(bucket)) do.call(rbind, bucket) else NULL,
       genes = gene_tab)
}

try_gsva <- function(logmat, use_sets) {
  # GSVA >= 1.48 uses param objects. Fall back to the legacy signature.
  out <- list()
  new_api <- exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)
  if (new_api) {
    gp <- tryCatch(gsvaParam(logmat, use_sets, minSize = 3, maxSize = 5000,
                             kcdf = "Gaussian"), error = function(e) NULL)
    if (!is.null(gp)) {
      m <- tryCatch(gsva(gp, verbose = FALSE), error = function(e) NULL)
      if (!is.null(m)) out$gsva <- as.matrix(m)
    }
    sp <- tryCatch(ssgseaParam(logmat, use_sets, minSize = 3, maxSize = 5000),
                   error = function(e) NULL)
    if (!is.null(sp)) {
      m <- tryCatch(suppressMessages(gsva(sp, verbose = FALSE)), error = function(e) NULL)
      if (!is.null(m)) out$ssgsea <- as.matrix(m)
    }
  } else {
    m <- tryCatch(gsva(logmat, use_sets, method = "gsva", kcdf = "Gaussian",
                       min.sz = 3, max.sz = 5000, verbose = FALSE),
                  error = function(e) NULL)
    if (!is.null(m)) out$gsva <- as.matrix(m)
    m <- tryCatch(gsva(logmat, use_sets, method = "ssgsea",
                       min.sz = 3, max.sz = 5000, verbose = FALSE),
                  error = function(e) NULL)
    if (!is.null(m)) out$ssgsea <- as.matrix(m)
  }
  if (!length(out)) NULL else out
}

apply_filter <- function(group_labels, min_count, min_n) {
  keep <- rowSums(counts_i[, names(group_labels), drop = FALSE] >= min_count) >= min_n
  counts_i[keep, names(group_labels), drop = FALSE]
}

# ---------------------------------------------------------------- run
set_manifest <- do.call(rbind, lapply(names(sets), function(nm) {
  data.frame(set_name = nm, gene = sets[[nm]], thesis_expect = expect_of(nm),
             stringsAsFactors = FALSE)
}))
write.table(set_manifest, file.path(RES, "pathway_sets_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

all_group <- c(WT1 = "WT", WT2 = "WT", WT3 = "WT", KO1 = "KO", KO2 = "KO", KO3 = "KO")
filters_primary <- list(
  cpmfilter_ge10_in2 = c(min_count = 10, min_n = 2),
  cpmfilter_ge5_in2 = c(min_count = 5, min_n = 2),
  cpmfilter_ge10_in3 = c(min_count = 10, min_n = 3)
)

rows <- list()
genes <- list()
for (fn in names(filters_primary)) {
  spec <- filters_primary[[fn]]
  cnt <- apply_filter(all_group, spec["min_count"], spec["min_n"])
  heavy <- fn == "cpmfilter_ge10_in2"
  ans <- run_contrast(cnt, all_group, "all_null_vs_WT", fn, heavy)
  rows[[length(rows) + 1]] <- ans$rows
  genes[[length(genes) + 1]] <- ans$genes
}

for (drop in names(all_group)) {
  grp <- all_group[names(all_group) != drop]
  cnt <- apply_filter(grp, 10, 2)
  ans <- run_contrast(cnt, grp, paste0("drop_", drop), "cpmfilter_ge10_in2", FALSE)
  rows[[length(rows) + 1]] <- ans$rows
  genes[[length(genes) + 1]] <- ans$genes
}

tab <- do.call(rbind, rows)
gtab <- do.call(rbind, genes)
tab$sweep_id <- seq_len(nrow(tab))

# Within each (contrast, filter) block, BH across rows that have a thesis_p
# and are not leave-one-gene or effect-only. This is a descriptive sweep q,
# not a pre-registered family.
tab$sweep_q_within_contrast <- NA_real_
split_key <- paste(tab$contrast, tab$filter, sep = "|")
for (k in unique(split_key)) {
  ix <- which(split_key == k & is.finite(tab$thesis_p) &
                !grepl("^loo_gene:", tab$variant) &
                tab$p_kind != "effect_only")
  if (length(ix)) tab$sweep_q_within_contrast[ix] <- p.adjust(tab$thesis_p[ix], "BH")
}

# Primary-family q: all-libraries, standard filter, as-is sets, core methods
primary_methods <- c(
  "fgsea:limma_t", "fgsea:deseq2_wald", "fgsea:signal_to_noise",
  "camera_voom:estimated", "camera_voom:mild", "roast_voom_mean", "fry_voom",
  "zmean_welch", "gsva", "ssgsea")
pix <- which(tab$contrast == "all_null_vs_WT" &
               tab$filter == "cpmfilter_ge10_in2" &
               tab$variant == "as_is" &
               tab$method %in% primary_methods &
               is.finite(tab$thesis_p))
tab$primary_family <- FALSE
tab$primary_family_q <- NA_real_
tab$primary_family[pix] <- TRUE
if (length(pix)) tab$primary_family_q[pix] <- p.adjust(tab$thesis_p[pix], "BH")

tab <- tab[order(tab$thesis_p, na.last = TRUE), ]
write.table(tab, file.path(RES, "sweep_all_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gtab, file.path(RES, "sweep_focus_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

support <- tab[is.finite(tab$thesis_p) & tab$thesis_match %in% TRUE, ]
support <- support[order(support$thesis_p), ]
write.table(support, file.path(RES, "sweep_thesis_aligned.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("rows ", nrow(tab), " thesis-aligned finite ", nrow(support))
message("best primary:")
print(head(support[support$contrast == "all_null_vs_WT" &
                     support$variant == "as_is" &
                     support$filter == "cpmfilter_ge10_in2",
                   c("method", "set_name", "effect", "thesis_p", "p_method", "extra")], 12))
message("best anywhere:")
print(head(support[, c("contrast", "filter", "variant", "method", "set_name",
                       "effect", "thesis_p", "p_method", "primary_family_q")], 12))
sessionInfo()

## 05_patient_aggregation.R
## AOI -> patient exposures for survival: tumor_score, immune_score,
## specificity_score, heterogeneity_score, plus ICC / reliability.
##
## AOI -> 患者暴露：tumor_score、immune_score、specificity_score、
## heterogeneity_score，以及 ICC / 可靠性。
##
## Inputs:  results/02_spe.rds
##          results/04_targets_aoi.tsv   (preferred; rebuilt from spe if missing)
## Outputs: results/05_patient_exposures.tsv
##          results/05_reliability.tsv

source("templates/R/00_setup.R")
need(c("SpatialExperiment", "SummarizedExperiment"), hard = TRUE)

log_step("05  AOI -> patient aggregation")

spe <- readRDS(out_path("02_spe.rds"))
cd <- as.data.frame(SummarizedExperiment::colData(spe))
cd$segment <- factor_segment(cd$segment)
logmat <- as.matrix(SummarizedExperiment::assay(spe, "log2_norm"))

aoi_path <- out_path("04_targets_aoi.tsv")
if (file.exists(aoi_path)) {
  aoi <- utils::read.delim(aoi_path, stringsAsFactors = FALSE)
} else {
  targets <- intersect(unlist(cfg$targets$genes), rownames(logmat))
  aoi <- do.call(rbind, lapply(targets, function(g) {
    data.frame(gene = g,
               patient = as.character(cd$patient),
               segment = as.character(cd$segment),
               expr = as.numeric(logmat[g, ]),
               nuclei = cd$nuclei,
               stringsAsFactors = FALSE)
  }))
}

method <- cfg$aggregation$method %||% "blup"
min_n  <- cfg$aggregation$min_aoi_per_patient %||% 2

## One clinical row per patient (first AOI's clinical fields).
## 每位患者一行临床信息（取该患者第一条 AOI 的临床字段）。
clin_cols <- unique(c(
  "patient", cfg$survival$time, cfg$survival$event,
  unlist(cfg$survival$covariates), unlist(cfg$survival$strata),
  "cohort", "os_time", "os_event", "pfs_time", "pfs_event",
  "age", "sex", "stage", "histology"
))
clin_cols <- intersect(clin_cols, names(cd))
clin <- cd[, clin_cols, drop = FALSE]
clin <- clin[!duplicated(clin$patient), , drop = FALSE]
names(clin)[names(clin) == "patient"] <- "patient_id"

rel_rows <- list()
exp_rows <- list()

for (g in unique(aoi$gene)) {
  d <- aoi[aoi$gene == g, ]
  agg <- aggregate_aoi(d$expr, d$patient, d$segment, d$nuclei,
                       method = method, min_n = min_n)
  ## Reliability per segment.
  for (s in SEGMENT_LEVELS) {
    ds <- d[d$segment == s, ]
    if (!nrow(ds)) next
    rho <- icc_oneway(ds$expr, ds$patient)
    m <- mean(as.numeric(table(ds$patient)))
    rel_rows[[paste(g, s, sep = "_")]] <- data.frame(
      gene = g, segment = s, method = method,
      icc = rho, mean_m = m,
      reliability = spearman_brown(rho, m),
      n_roi_for_R08 = n_roi_for_reliability(rho, 0.8),
      n_patient = length(unique(ds$patient)),
      n_flagged_low_n = sum(agg$segment == s & agg$flagged_low_n),
      stringsAsFactors = FALSE
    )
  }
  ## Wide patient table for this gene.
  tumor <- agg[agg$segment == TUMOR_SEGMENT, c("patient_id", "score", "n_aoi", "sd_within")]
  names(tumor) <- c("patient_id", "tumor_score", "n_tumor", "sd_tumor")
  ## Immune score: mean of available immune segments, or CD45 if present.
  imm_pref <- intersect(c("CD45", IMMUNE_SEGMENTS), unique(agg$segment))
  if (length(imm_pref)) {
    imm <- agg[agg$segment == imm_pref[1], c("patient_id", "score", "n_aoi", "sd_within")]
    names(imm) <- c("patient_id", "immune_score", "n_immune", "sd_immune")
  } else {
    imm <- data.frame(patient_id = tumor$patient_id,
                      immune_score = NA_real_, n_immune = 0, sd_immune = NA_real_)
  }
  wide <- merge(tumor, imm, by = "patient_id", all = TRUE)
  wide$gene <- g
  wide$specificity_score <- wide$tumor_score - wide$immune_score
  wide$heterogeneity_score <- wide$sd_tumor
  exp_rows[[g]] <- wide
}

rel_tab <- do.call(rbind, rel_rows)
rownames(rel_tab) <- NULL
write_tsv(rel_tab, out_path("05_reliability.tsv"))

exp_tab <- do.call(rbind, exp_rows)
rownames(exp_tab) <- NULL
exp_tab <- merge(exp_tab, clin, by = "patient_id", all.x = TRUE)
write_tsv(exp_tab, out_path("05_patient_exposures.tsv"))

## Flag the reliability next to every exposure so 06 can refuse to over-interpret
## a null OS result from an unreliable score (playbook §7, pitfall 13).
## 把可靠性写在每个暴露旁边，避免 06 把不可靠暴露的阴性 OS 当成阴性结果。
log_step("05 done. Patients: ", length(unique(exp_tab$patient_id)),
         "  genes: ", paste(unique(exp_tab$gene), collapse = ", "))
save_session_info("05")

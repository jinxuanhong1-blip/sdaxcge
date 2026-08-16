## 01_qc_dcc_to_spe.R
## Public GeoMx DCC/PKC -> GeomxSet -> segment/probe/LOQ QC -> SpatialExperiment
## 公开 GeoMx DCC/PKC -> 质控 -> SpatialExperiment
##
## Inputs  (from study_config.yml):
##   paths.dcc_dir, paths.pkc_files, paths.annotation
## Outputs:
##   results/01_spe.rds
##   results/01_qc_aoi.tsv
##   results/01_qc_genes.tsv
##   results/01_confound_checks.txt
##
## This script talks to local files only. Fetch public GSE271689 DCCs first with
## 08_fetch_public_geo.R. No private / restricted access data.
## 本脚本只读本地文件。请先用 08_fetch_public_geo.R 拉取公开的 GSE271689 DCC。

source("templates/R/00_setup.R")
need(c("GeomxTools", "NanoStringNCTools", "SpatialExperiment"), hard = TRUE)
need(c("standR"), hard = FALSE)

log_step("01  QC: DCC -> SpatialExperiment")

## ---------------------------------------------------------------------------
## Annotation / 注释
## ---------------------------------------------------------------------------
ann_path <- cfg$paths$annotation
if (!file.exists(ann_path)) {
  stop("Annotation file not found: ", ann_path,
       "\nCopy templates/config/annotation_template.csv and fill it, ",
       "or run 08_fetch_public_geo.R for the public GSE271689 scaffold.",
       call. = FALSE)
}
ann <- read_annotation(ann_path, sheet = cfg$paths$annotation_sheet)
map <- cfg$annotation
require_cols(ann, c(map$dcc_file, map$slide, map$patient, map$segment, map$cohort),
             context = "annotation")

ann$`.dcc`     <- as.character(ann[[map$dcc_file]])
ann$`.slide`   <- as.character(ann[[map$slide]])
ann$`.patient` <- as.character(ann[[map$patient]])
ann$`.segment` <- factor_segment(ann[[map$segment]])
ann$`.cohort`  <- as.character(ann[[map$cohort]])
ann$`.roi`     <- if (!is.na(ann_col(cfg, "roi"))) as.character(ann[[map$roi]]) else NA_character_
ann$`.area`    <- if (!is.na(ann_col(cfg, "area"))) as.numeric(ann[[map$area]]) else NA_real_
ann$`.nuclei`  <- if (!is.na(ann_col(cfg, "nuclei"))) as.numeric(ann[[map$nuclei]]) else NA_real_

## Structural checks before any count is touched.
## 碰计数之前先做结构检查。
sink(out_path("01_confound_checks.txt"))
cat("AOIs per patient x segment\n")
print(table(ann$`.patient`, ann$`.segment`))
cat("\nslide x cohort\n")
print(table(ann$`.slide`, ann$`.cohort`))
cat("\nslide x segment\n")
print(table(ann$`.slide`, ann$`.segment`))
sink()
report_confound(confound_check(ann, ".slide", ".cohort", "slide vs cohort"))

## ---------------------------------------------------------------------------
## Read DCC / 读取 DCC
## ---------------------------------------------------------------------------
dcc_dir <- cfg$paths$dcc_dir
pkc <- unlist(cfg$paths$pkc_files)
if (!dir.exists(dcc_dir)) stop("DCC directory not found: ", dcc_dir, call. = FALSE)
dcc_files <- list.files(dcc_dir, pattern = "\\.dcc$", full.names = TRUE)
if (!length(dcc_files)) stop("No .dcc files in ", dcc_dir, call. = FALSE)
log_step("found ", length(dcc_files), " DCC files")

## GeomxTools::readNanoStringGeoMxSet expects a protocolData annotation with
## Sample_ID matching the DCC basename. We rewrite a temporary annotation so
## the user's column names do not have to match the vendor template exactly.
## 临时重写注释，使 Sample_ID 与 DCC 文件名对齐，用户列名不必与厂商模板完全一致。
tmp_ann <- ann
tmp_ann$Sample_ID <- tools::file_path_sans_ext(basename(tmp_ann$`.dcc`))
tmp_csv <- tempfile(fileext = ".csv")
utils::write.csv(tmp_ann, tmp_csv, row.names = FALSE)

gx <- GeomxTools::readNanoStringGeoMxSet(
  dccFiles = dcc_files,
  pkcFiles = pkc,
  phenoDataFile = tmp_csv,
  phenoDataSheet = NULL,
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = intersect(c("aoi", "roi", "slide name", "scan name"),
                                   names(tmp_ann)),
  experimentDataColNames = character(0)
)

## Harmonize the role columns onto pData so later scripts never re-map.
## 把角色列写回 pData，后续脚本不再做列名映射。
pd <- Biobase::pData(gx)
pd$slide   <- pd$`.slide`   %||% pd[[map$slide]]
pd$patient <- pd$`.patient` %||% pd[[map$patient]]
pd$segment <- factor_segment(pd$`.segment` %||% pd[[map$segment]])
pd$cohort  <- pd$`.cohort`  %||% pd[[map$cohort]]
pd$roi_id  <- pd$`.roi`
pd$area_um2 <- pd$`.area`
pd$nuclei   <- pd$`.nuclei`
Biobase::pData(gx) <- pd

## ---------------------------------------------------------------------------
## Segment QC / AOI 质控
## ---------------------------------------------------------------------------
qc <- cfg$qc
gx <- GeomxTools::setSegmentQCFlags(
  gx,
  qcCutoffs = list(
    minSegmentReads = qc$min_segment_reads,
    percentTrimmed = qc$percent_trimmed,
    percentStitched = qc$percent_stitched,
    percentAligned = qc$percent_aligned,
    percentSaturation = qc$percent_saturation,
    minNegativeCount = qc$min_negative_count,
    maxNTCCount = qc$max_ntc_count,
    minNuclei = min(vapply(cfg$segments$qc, `[[`, numeric(1), "min_nuclei")),
    minArea = min(vapply(cfg$segments$qc, `[[`, numeric(1), "min_area"))
  )
)

## Per-segment nuclei/area floors (the vendor call above used the *minimum*
## across segments so it does not silently drop CD68). Apply the real floors here.
## 上面的厂商调用用的是各分区阈值的最小值，以免悄悄删掉 CD68；真正的分区阈值在这里施加。
pd <- Biobase::pData(gx)
seg_qc <- cfg$segments$qc
drop_seg <- rep(FALSE, nrow(pd))
for (s in names(seg_qc)) {
  idx <- which(as.character(pd$segment) == s)
  if (!length(idx)) next
  thr <- seg_qc[[s]]
  fail_n <- !is.na(pd$nuclei[idx]) & pd$nuclei[idx] < thr$min_nuclei
  fail_a <- !is.na(pd$area_um2[idx]) & pd$area_um2[idx] < thr$min_area
  drop_seg[idx] <- fail_n | fail_a
}
flag_cols <- grep("QCFlags|qcFlags|LowReads|LowSaturation|LowAligned",
                  names(pd), value = TRUE, ignore.case = TRUE)
vendor_fail <- if (length(flag_cols)) {
  rowSums(as.data.frame(pd[, flag_cols, drop = FALSE]) == TRUE, na.rm = TRUE) > 0
} else {
  rep(FALSE, nrow(pd))
}
keep_aoi <- !drop_seg & !vendor_fail
log_step("AOI QC: keep ", sum(keep_aoi), " / ", length(keep_aoi),
         "  dropped per segment: ",
         paste(names(table(pd$segment[!keep_aoi])),
               table(pd$segment[!keep_aoi]), sep = "=", collapse = ", "))
gx <- gx[, keep_aoi]

## ---------------------------------------------------------------------------
## Probe QC + collapse / 探针质控并聚合到基因
## ---------------------------------------------------------------------------
gx <- GeomxTools::setBioProbeQCFlags(gx)
probe_fail <- if ("QCFlags" %in% names(Biobase::fData(gx))) {
  rowSums(as.data.frame(Biobase::fData(gx)$QCFlags) == TRUE, na.rm = TRUE) > 0
} else {
  rep(FALSE, nrow(gx))
}
gx <- gx[!probe_fail, ]
gx <- GeomxTools::aggregateCounts(gx)
log_step("after probe collapse: ", nrow(gx), " targets x ", ncol(gx), " AOIs")

## ---------------------------------------------------------------------------
## LOQ and gene/segment filters / LOQ 与基因、AOI 过滤
## ---------------------------------------------------------------------------
neg_mean_col <- grep("NegGeoMean", names(Biobase::pData(gx)), value = TRUE)[1]
neg_sd_col   <- grep("NegGeoSD",   names(Biobase::pData(gx)), value = TRUE)[1]
if (is.na(neg_mean_col) || is.na(neg_sd_col)) {
  warning("NegGeoMean / NegGeoSD not in pData; computing from negative probes if present.")
  ## Fallback: geometric mean/SD of any feature tagged as Negative.
  fd <- Biobase::fData(gx)
  neg_idx <- grepl("Neg", fd$TargetName %||% rownames(fd), ignore.case = TRUE) |
    grepl("Negative", fd$CodeClass %||% "", ignore.case = TRUE)
  if (any(neg_idx)) {
    neg <- as.matrix(Biobase::exprs(gx)[neg_idx, , drop = FALSE])
    Biobase::pData(gx)$NegGeoMean_WTA <- exp(colMeans(log(pmax(neg, 1))))
    Biobase::pData(gx)$NegGeoSD_WTA <- exp(apply(log(pmax(neg, 1)), 2, stats::sd))
    neg_mean_col <- "NegGeoMean_WTA"
    neg_sd_col <- "NegGeoSD_WTA"
  } else {
    stop("Cannot compute LOQ: no negative-probe summary in pData and no Negative features.",
         call. = FALSE)
  }
}

loq <- loq_from_neg(Biobase::pData(gx)[[neg_mean_col]],
                    Biobase::pData(gx)[[neg_sd_col]],
                    cutoff = qc$loq_sd_cutoff,
                    min_loq = qc$min_loq)
Biobase::pData(gx)$LOQ <- loq
counts <- as.matrix(Biobase::exprs(gx))
Biobase::pData(gx)$GeneDetectionRate <- detection_rate(counts, loq)

keep_aoi2 <- Biobase::pData(gx)$GeneDetectionRate >= qc$min_segment_detection_rate
log_step("LOQ AOI filter: keep ", sum(keep_aoi2), " / ", length(keep_aoi2))
gx <- gx[, keep_aoi2]
counts <- as.matrix(Biobase::exprs(gx))
loq <- Biobase::pData(gx)$LOQ

det <- gene_detection_by_segment(counts, loq, Biobase::pData(gx)$segment)
## Intersection for cross-compartment work; union is written out for within-compartment work.
## 跨分区用交集，分区内用并集，两者都写出。
genes_union <- keep_genes_per_segment(det, qc$min_gene_detection_rate, "union")
genes_inter <- keep_genes_per_segment(det, qc$min_gene_detection_rate, "intersection")
log_step("gene filter: union ", length(genes_union),
         "  intersection ", length(genes_inter),
         "  of ", nrow(gx))

## Default object keeps the UNION (within-compartment analyses need it).
## Cross-compartment scripts subset to the intersection themselves.
## 默认对象保留并集；跨分区脚本自行子集到交集。
gx <- gx[genes_union, ]

## ---------------------------------------------------------------------------
## SpatialExperiment / 转为 SpatialExperiment
## ---------------------------------------------------------------------------
spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = as.matrix(Biobase::exprs(gx))),
  colData = S4Vectors::DataFrame(Biobase::pData(gx)),
  rowData = S4Vectors::DataFrame(Biobase::fData(gx))
)
S4Vectors::metadata(spe)$genes_union <- genes_union
S4Vectors::metadata(spe)$genes_intersection <- genes_inter
S4Vectors::metadata(spe)$detection_by_segment <- det
S4Vectors::metadata(spe)$config_path <- CONFIG_PATH
S4Vectors::metadata(spe)$qc_summary <- list(
  n_dcc = length(dcc_files),
  n_aoi_in = nrow(ann),
  n_aoi_out = ncol(spe),
  n_gene_in = nrow(det),
  n_gene_union = length(genes_union),
  n_gene_intersection = length(genes_inter)
)

saveRDS(spe, out_path("01_spe.rds"))

aoi_qc <- as.data.frame(SummarizedExperiment::colData(spe))
aoi_qc$sample_id <- colnames(spe)
write_tsv(aoi_qc, out_path("01_qc_aoi.tsv"))

gene_qc <- data.frame(
  gene = rownames(det),
  det,
  keep_union = rownames(det) %in% genes_union,
  keep_intersection = rownames(det) %in% genes_inter,
  stringsAsFactors = FALSE
)
write_tsv(gene_qc, out_path("01_qc_genes.tsv"))

save_session_info("01")
log_step("01 done -> ", out_path("01_spe.rds"))

#!/usr/bin/env Rscript
# 06 · GeoMx DSP-WTA QC + normalization. GeoMx DSP-WTA 质控 + 归一化。
# ---------------------------------------------------------------------------
# Builds a NanoStringGeoMxSet from DCC + PKC + annotation, applies segment (AOI)
# and probe QC, aggregates probes to genes, computes the limit of quantification
# (LOQ), filters genes by detection rate, and applies Q3 (upper-quartile)
# normalization -- the GeoMx standard.
# 由 DCC + PKC + 注释构建 GeoMxSet，进行 AOI 与探针质控，探针聚合到基因，计算
# 定量下限（LOQ），按检出率过滤基因，并做 Q3（上四分位）归一化——GeoMx 标准流程。
#
# Tools: GeomxTools, GeoMxWorkflows, NanoStringNCTools (Bioconductor 3.18+, 2024+).
# No results are fabricated; every number comes from the assay data.
# 工具：GeomxTools 等。不伪造结果；所有数值来自实测数据。
#
# Inputs (see config.yaml and demo/README.md):
#   dcc_dir      : directory of *.dcc files (GSE271689_RAW.tar extracts here)
#   pkc          : Hs_R_NGS_WTA_v1.0.pkc  (from NanoString; NOT in GEO RAW.tar)
#   annotation   : AOI-level metadata sheet (slide, ROI, segment/compartment,
#                  patient, plus clinical columns if survival is planned)
#
# Usage:
#   Rscript 06_geomx_qc_normalization.R \
#       --dcc <dir> --pkc <file.pkc> --annotation <sheet.xlsx|csv> --out <dir>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(GeomxTools)
  library(GeoMxWorkflows)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dcc", type = "character"),
  make_option("--pkc", type = "character"),
  make_option("--annotation", type = "character"),
  make_option("--out", type = "character", default = "methods/spatial/demo/out/geomx"),
  # QC gates mirror config.yaml:qc.geomx; override on the CLI if needed.
  make_option("--min_reads", type = "double", default = 1000),
  make_option("--pct_trimmed", type = "double", default = 80),
  make_option("--pct_aligned", type = "double", default = 75),
  make_option("--pct_saturation", type = "double", default = 50),
  make_option("--min_nuclei", type = "double", default = 20),
  make_option("--min_area", type = "double", default = 1600),
  make_option("--loq_sd", type = "double", default = 2),
  make_option("--detection_rate", type = "double", default = 0.05)
)))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

dcc_files <- dir(opt$dcc, pattern = "\\.dcc$", full.names = TRUE, recursive = TRUE)
if (length(dcc_files) == 0) {
  # GEO ships .dcc.gz; decompress once if that is what we have.
  gz <- dir(opt$dcc, pattern = "\\.dcc\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(gz) > 0) {
    for (f in gz) R.utils::gunzip(f, overwrite = TRUE, remove = FALSE)
    dcc_files <- dir(opt$dcc, pattern = "\\.dcc$", full.names = TRUE, recursive = TRUE)
  }
}
stopifnot(length(dcc_files) > 0)

geomx <- readNanoStringGeoMxSet(
  dccFiles = dcc_files,
  pkcFiles = opt$pkc,
  phenoDataFile = opt$annotation,
  phenoDataSheet = "Template",           # adjust to your sheet name
  phenoDataDccColName = "Sample_ID",     # column linking annotation rows to DCC ids
  protocolDataColNames = c("aoi", "roi"),
  experimentDataColNames = c("panel")
)

# --- Shift zero counts by 1 so log/geomean operations are defined. ---
geomx <- shiftCountsOne(geomx, useDALogic = TRUE)

# --- Segment (AOI) QC. AOI 质控。 ---
geomx <- setSegmentQCFlags(geomx, qcCutoffs = list(
  minSegmentReads = opt$min_reads,
  percentTrimmed  = opt$pct_trimmed,
  percentStitched = 80,
  percentAligned  = opt$pct_aligned,
  percentSaturation = opt$pct_saturation,
  minNegativeCount = 1,
  maxNTCCount = 9000,
  minNuclei = opt$min_nuclei,
  minArea   = opt$min_area
))
qc <- sData(geomx)
write.csv(qc, file.path(opt$out, "segment_qc.csv"))
geomx <- geomx[, qc$QCFlags$LowReads == FALSE &
                 qc$QCFlags$LowSaturation == FALSE &
                 qc$QCFlags$LowNuclei == FALSE]

# --- Probe QC + aggregate probes to gene targets. 探针质控并聚合到基因。 ---
geomx <- setBioProbeQCFlags(geomx)
probe_qc <- fData(geomx)
geomx <- subset(geomx, subset = !(probe_qc$QCFlags$LowProbeRatio |
                                   probe_qc$QCFlags$GlobalGrubbsOutlier))
target_geomx <- aggregateCounts(geomx)

# --- Limit of quantification (LOQ) per segment from negative probes. ---
# 每个 AOI 由阴性探针几何均值 * SD^n 计算 LOQ。
neg_geomeans <- esBy(negativeControlSubset(target_geomx), GROUP = "Module",
                     FUN = function(x) assayDataApply(x, 2, ngeoMean, elt = "exprs"))
loq <- pmax(1, neg_geomeans * (opt$loq_sd))
pData(target_geomx)$LOQ <- loq

exprs_mat <- exprs(target_geomx)
above_loq <- exprs_mat > matrix(loq, nrow = nrow(exprs_mat), ncol = ncol(exprs_mat), byrow = TRUE)
gene_detect_rate <- rowSums(above_loq) / ncol(above_loq)
keep_genes <- gene_detect_rate >= opt$detection_rate
message(sprintf("[06_geomx] keeping %d/%d genes at detection rate >= %.2f",
                sum(keep_genes), length(keep_genes), opt$detection_rate))
target_geomx <- target_geomx[keep_genes, ]

# --- Q3 normalization (upper-quartile). Q3 归一化。 ---
target_geomx <- normalize(target_geomx, norm_method = "quant",
                          desiredQuantile = 0.75, toElt = "q_norm")

saveRDS(target_geomx, file.path(opt$out, "geomx_target_qnorm.rds"))
write.csv(assayDataElement(target_geomx, "q_norm"),
          file.path(opt$out, "geomx_q3_normalized.csv"))
message(sprintf("[06_geomx] wrote normalized set: %d genes x %d segments -> %s",
                nrow(target_geomx), ncol(target_geomx), opt$out))

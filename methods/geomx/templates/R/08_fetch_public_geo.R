## 08_fetch_public_geo.R
## Download PUBLIC GeoMx WTA files for GSE271689-class analyses.
## 下载用于 GSE271689 类分析的【公开】GeoMx WTA 文件。
##
## Scope (public only / 仅公开数据):
##   * GEO series GSE271689  -- Yale NSCLC GeoMx WTA DCC (Aung et al., Nat Genet 2025)
##   * GEO series GSE292098  -- companion GeoMx DSP accession cited by the same paper
##   * Human WTA PKC from the Bruker/NanoString public PKC bundle (user-supplied path
##     if the automated URL is unavailable)
##
## This script never touches dbGaP, restricted access, or unpublished counts.
## 本脚本不接触 dbGaP、受限访问或未发表计数。
##
## Usage:
##   Rscript templates/R/08_fetch_public_geo.R
##   GEOMX_GEO=GSE271689 Rscript templates/R/08_fetch_public_geo.R
##
## Outputs (under data/public/<GSE>/):
##   RAW.tar            original GEO supplementary archive
##   dcc/*.dcc          extracted digital count conversion files
##   annotation_scaffold.csv   DCC filenames only; YOU must join public phenotype
##                             from the GEO SOFT / series matrix / the paper
##   README_PUBLIC.txt  provenance

source("templates/R/00_setup.R")

gse <- Sys.getenv("GEOMX_GEO", unset = "GSE271689")
allowed <- c("GSE271689", "GSE292098")
if (!gse %in% allowed) {
  stop("Refusing to fetch '", gse, "'. This template only downloads the public ",
       "GSE271689-class accessions: ", paste(allowed, collapse = ", "),
       call. = FALSE)
}

## Keep all downloaded public files under methods/geomx/ (playbook output rule).
## 所有公开下载文件都放在 methods/geomx/ 下（手册输出范围）。
dest_root <- file.path(.GEOMX_DIR, "data", "public", gse)
dcc_dir   <- file.path(dest_root, "dcc")
dir.create(dcc_dir, recursive = TRUE, showWarnings = FALSE)

## GEO supplementary file naming is stable: GSE#####_RAW.tar
## GEO 补充文件命名稳定：GSE#####_RAW.tar
ftp <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s_RAW.tar",
               sub("\\d{3}$", "nnn", gse), gse, gse)
tar_path <- file.path(dest_root, paste0(gse, "_RAW.tar"))

log_step("public fetch: ", gse)
log_step("URL: ", ftp)

if (!file.exists(tar_path) || file.info(tar_path)$size < 1000) {
  ok <- try(utils::download.file(ftp, tar_path, mode = "wb", quiet = FALSE), silent = TRUE)
  if (inherits(ok, "try-error") || !file.exists(tar_path)) {
    stop("Download failed. Confirm the accession is public and the network allows ",
         "ftp.ncbi.nlm.nih.gov. Manual fallback:\n  ", ftp, "\n",
         call. = FALSE)
  }
} else {
  log_step("using cached archive: ", tar_path)
}

## Extract. GEO RAW.tar typically contains gzipped DCC files, one per GSM.
## 解压。GEO RAW.tar 通常每个 GSM 一个 gzip 压缩的 DCC。
untar_dir <- file.path(dest_root, "untar")
dir.create(untar_dir, showWarnings = FALSE)
utils::untar(tar_path, exdir = untar_dir)

gz <- list.files(untar_dir, pattern = "\\.(dcc|dcc\\.gz)$", recursive = TRUE, full.names = TRUE)
if (!length(gz)) {
  ## Sometimes the tar holds a nested tar or zip of DCCs.
  nested <- list.files(untar_dir, pattern = "\\.(tar|tar\\.gz|tgz|zip)$",
                       recursive = TRUE, full.names = TRUE)
  for (n in nested) {
    if (grepl("\\.zip$", n)) utils::unzip(n, exdir = untar_dir)
    else utils::untar(n, exdir = untar_dir)
  }
  gz <- list.files(untar_dir, pattern = "\\.(dcc|dcc\\.gz)$", recursive = TRUE, full.names = TRUE)
}
if (!length(gz)) {
  stop("No .dcc files inside ", tar_path,
       ". Inspect ", untar_dir, " and extract manually.", call. = FALSE)
}

copied <- 0L
for (f in gz) {
  dest <- file.path(dcc_dir, sub("\\.gz$", "", basename(f)))
  if (grepl("\\.gz$", f)) {
    con <- gzfile(f, "rb")
    writeBin(readBin(con, "raw", n = file.info(f)$size %||% 1e8), dest)
    close(con)
  } else {
    file.copy(f, dest, overwrite = TRUE)
  }
  copied <- copied + 1L
}
log_step("extracted ", copied, " DCC files -> ", dcc_dir)

dcc_names <- list.files(dcc_dir, pattern = "\\.dcc$")
scaffold <- data.frame(
  Sample_ID = dcc_names,
  slide_name = NA_character_,
  roi_id = NA_character_,
  patient_id = NA_character_,
  segment = NA_character_,
  area_um2 = NA_real_,
  nuclei = NA_real_,
  cohort = NA_character_,
  os_time = NA_real_,
  os_event = NA_integer_,
  stringsAsFactors = FALSE
)
utils::write.csv(scaffold, file.path(dest_root, "annotation_scaffold.csv"), row.names = FALSE)

readme <- c(
  paste("Public GeoMx fetch for", gse),
  paste("Downloaded:", format(Sys.time(), tz = "UTC"), "UTC"),
  paste("URL:", ftp),
  paste("DCC count:", copied),
  "",
  "PROVENANCE / 来源",
  "  Aung et al. Nat Genet 2025;57:2482-2493. PMID 41073787. doi:10.1038/s41588-025-02351-7",
  "  GEO GSE271689 / GSE292098 (public supplementary DCC).",
  "",
  "WHAT THIS DOES NOT INCLUDE / 本下载不包含",
  "  * Patient-level OS/PFS tables (not always in the GEO series matrix).",
  "  * PKC probe kit (download the public Human WTA PKC from Bruker/NanoString).",
  "  * Restricted or unpublished counts. Do not add any.",
  "",
  "NEXT / 下一步",
  "  1. Join public phenotype from the GEO SOFT / series matrix / the paper supplement.",
  "  2. Map morphology-marker AOIs to segment = Tumor | CD45 | CD68.",
  "     NOTE: the GSE271689 GEO protocol text lists S100B as a tumor marker;",
  "     the Nature Genetics methods describe PanCK for the NSCLC cohorts.",
  "     Reconcile against the published methods before modelling.",
  "  3. Point study_config.yml paths.dcc_dir at methods/geomx/data/public/<GSE>/dcc",
  "     and paths.annotation at methods/geomx/data/public/<GSE>/annotation.csv.",
  "  4. Run 00 -> 01 -> 02 -> 03 -> 04 -> 05 -> 06."
)
writeLines(readme, file.path(dest_root, "README_PUBLIC.txt"))
writeLines(readme, out_path("08_public_fetch.txt"))

log_step("08 done. Scaffold: ", file.path(dest_root, "annotation_scaffold.csv"))
log_step("Join PUBLIC phenotype only. Do not add restricted clinical tables.")

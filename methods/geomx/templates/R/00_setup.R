## 00_setup.R -- packages, configuration, output scaffolding
## GeoMx WTA compartment/mixed-effects/survival playbook
##
## Source this at the top of every other script:
##   source("templates/R/00_setup.R")
##
## 在其他脚本开头 source 本文件。

## ---------------------------------------------------------------------------
## Packages / 依赖包
## ---------------------------------------------------------------------------
## Bioconductor: NanoStringNCTools GeomxTools GeoMxWorkflows standR SpatialExperiment
##               SpatialDecon limma edgeR GeoDiff
## CRAN:         yaml dplyr tidyr ggplot2 lme4 lmerTest emmeans variancePartition
##               survival coxme rms timeROC maxstat broom
##
## install.packages("BiocManager")
## BiocManager::install(c("NanoStringNCTools", "GeomxTools", "GeoMxWorkflows", "standR",
##                        "SpatialExperiment", "SpatialDecon", "limma", "edgeR", "GeoDiff"))
## install.packages(c("yaml", "dplyr", "tidyr", "ggplot2", "lme4", "lmerTest", "emmeans",
##                    "variancePartition", "survival", "coxme", "rms", "timeROC",
##                    "maxstat", "broom"))

need <- function(pkgs, hard = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    msg <- sprintf("Missing packages: %s", paste(missing, collapse = ", "))
    if (hard) stop(msg, call. = FALSE) else message(msg, " (optional, skipping related steps)")
  }
  invisible(length(missing) == 0)
}

need(c("yaml"))

`%||%` <- function(x, y) if (is.null(x)) y else x

## ---------------------------------------------------------------------------
## Configuration / 配置
## ---------------------------------------------------------------------------
## Resolve 00_setup.R's own directory so utils_geomx.R is found regardless of cwd.
## 解析 00_setup.R 自身所在目录，使无论工作目录在哪都能找到 utils_geomx.R。
.SETUP_FILE <- (function() {
  of <- sys.frames()
  for (i in rev(seq_along(of))) {
    f <- of[[i]]$ofile
    if (!is.null(f)) return(normalizePath(f))
  }
  ca <- grep("[.]R$", commandArgs(trailingOnly = FALSE), value = TRUE)
  ca <- sub("^--file=", "", ca)
  if (length(ca) && file.exists(ca[1])) return(normalizePath(ca[1]))
  if (file.exists("templates/R/00_setup.R")) return(normalizePath("templates/R/00_setup.R"))
  if (file.exists("methods/geomx/templates/R/00_setup.R"))
    return(normalizePath("methods/geomx/templates/R/00_setup.R"))
  stop("Cannot locate 00_setup.R; set cwd to the repo root.", call. = FALSE)
})()
.SETUP_DIR <- dirname(.SETUP_FILE)
.TMPL_DIR  <- dirname(.SETUP_DIR)   # .../templates
.GEOMX_DIR <- dirname(.TMPL_DIR)    # .../geomx

CONFIG_PATH <- Sys.getenv("GEOMX_CONFIG", unset = "")
if (!nzchar(CONFIG_PATH) || !file.exists(CONFIG_PATH)) {
  candidates <- c(
    file.path(.TMPL_DIR, "config", "study_config.yml"),
    file.path(.TMPL_DIR, "config", "study_config.example.yml"),
    "templates/config/study_config.yml",
    "templates/config/study_config.example.yml"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("No config found. Copy study_config.example.yml to study_config.yml.", call. = FALSE)
  }
  if (grepl("example", hit)) {
    message("Using the EXAMPLE config: ", hit, ". Copy and edit it before a real run.")
  }
  CONFIG_PATH <- hit
}
cfg <- yaml::read_yaml(CONFIG_PATH)

## Fail loudly on the settings that silently corrupt an analysis if wrong.
## 对那些一旦设错就会悄悄毁掉整个分析的配置项，直接报错。
stopifnot(
  length(cfg$segments$levels) >= 2,
  cfg$segments$tumor %in% cfg$segments$levels,
  all(cfg$segments$immune %in% cfg$segments$levels),
  cfg$normalization$primary %in% c("q3", "tmm", "quantile", "geodiff"),
  cfg$differential_expression$block == "patient"
)

SEGMENT_LEVELS <- unlist(cfg$segments$levels)
TUMOR_SEGMENT  <- cfg$segments$tumor
IMMUNE_SEGMENTS <- unlist(cfg$segments$immune)

set.seed(cfg$reproducibility$seed %||% 1)

## ---------------------------------------------------------------------------
## Output directories / 输出目录
## ---------------------------------------------------------------------------
RESULTS_DIR <- cfg$paths$results_dir
FIGURES_DIR <- cfg$paths$figures_dir
for (d in c(RESULTS_DIR, FIGURES_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

out_path <- function(...) file.path(RESULTS_DIR, ...)
fig_path <- function(...) file.path(FIGURES_DIR, ...)

## ---------------------------------------------------------------------------
## Helpers / 辅助函数
## ---------------------------------------------------------------------------
source(file.path(.SETUP_DIR, "utils_geomx.R"))

log_step <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

## Write a session record next to the results. Reviewers ask for this; write it once, automatically.
## 把会话信息与结果放在一起自动保存；审稿人一定会要。
save_session_info <- function(tag = "session") {
  if (isTRUE(cfg$reproducibility$save_session_info)) {
    f <- out_path(sprintf("%s_sessionInfo.txt", tag))
    writeLines(capture.output(utils::sessionInfo()), f)
    log_step("session info -> ", f)
  }
}

log_step("config: ", CONFIG_PATH)
log_step("segments: ", paste(SEGMENT_LEVELS, collapse = " | "),
         "  (tumor = ", TUMOR_SEGMENT, ")")
log_step("primary normalization: ", cfg$normalization$primary,
         "  |  primary DE method: ", cfg$differential_expression$primary_method)

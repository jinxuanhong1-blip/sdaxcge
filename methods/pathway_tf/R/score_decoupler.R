#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript score_decoupler.R path/to/config.yml", call. = FALSE)
}

required <- c("yaml", "decoupleR")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

cfg <- yaml::read_yaml(args[[1]])
expr <- readRDS(cfg$inputs$expression_rds)

if (!is.matrix(expr) || !is.numeric(expr) || is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("expression_rds must contain a numeric genes-by-samples matrix with dimnames.", call. = FALSE)
}
if (anyNA(expr) || anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) {
  stop("Expression matrix must have no missing values or duplicated dimnames.", call. = FALSE)
}
if (!identical(cfg$analysis$assay_scale, "log")) {
  stop("PROGENy/DoRothEA activity inference requires normalized log-scale expression.", call. = FALSE)
}

metadata <- read.delim(
  cfg$inputs$metadata_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(metadata$sample_id, colnames(expr))) {
  stop("metadata sample_id must exactly match expression columns in the same order.", call. = FALSE)
}

organism <- cfg$analysis$organism
if (!organism %in% c("human", "mouse")) {
  stop("analysis.organism must be human or mouse.", call. = FALSE)
}

progeny <- decoupleR::get_progeny(
  organism = organism,
  top = cfg$decoupler$progeny_top
)
dorothea <- decoupleR::get_dorothea(
  organism = organism,
  levels = unlist(cfg$decoupler$dorothea_confidence)
)

run_method <- function(method, network, network_name) {
  runner <- switch(
    method,
    ulm = decoupleR::run_ulm,
    mlm = decoupleR::run_mlm,
    stop("Unsupported decoupleR method: ", method, call. = FALSE)
  )
  result <- runner(
    mat = expr,
    net = network,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = cfg$decoupler$minsize
  )
  result$method <- method
  result$network <- network_name
  result
}

# DoRothEA stores signed regulator-target weights in `mor`; standardize its name.
dorothea$weight <- dorothea$mor
methods <- unlist(cfg$decoupler$methods)
scores <- do.call(
  rbind,
  c(
    lapply(methods, run_method, network = progeny, network_name = "PROGENy"),
    lapply(methods, run_method, network = dorothea, network_name = "DoRothEA")
  )
)

out_dir <- cfg$analysis$output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  scores,
  file.path(out_dir, "decoupler_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(scores, file.path(out_dir, "decoupler_scores.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "decoupler_sessionInfo.txt"))

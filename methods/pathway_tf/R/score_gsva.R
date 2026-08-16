#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript score_gsva.R path/to/config.yml", call. = FALSE)
}

required <- c("yaml", "GSVA", "msigdbr")
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
  stop("This template accepts normalized log-scale expression only.", call. = FALSE)
}

metadata <- read.delim(
  cfg$inputs$metadata_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(metadata$sample_id, colnames(expr))) {
  stop("metadata sample_id must exactly match expression columns in the same order.", call. = FALSE)
}

read_gmt <- function(path) {
  if (is.null(path) || !nzchar(path) || identical(path, "CHANGE_ME") || !file.exists(path)) {
    stop("inputs.custom_gmt must point to a curated GMT containing junction signatures.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines[nzchar(lines)], "\t", fixed = TRUE)
  if (!length(fields) || any(lengths(fields) < 3L)) {
    stop("Each GMT row needs a set name, description, and at least one gene.", call. = FALSE)
  }
  sets <- lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  names(sets) <- vapply(fields, `[[`, character(1), 1L)
  if (anyDuplicated(names(sets))) stop("GMT set names must be unique.", call. = FALSE)
  if (any(grepl("^REPLACE_", unlist(sets)))) {
    stop("Replace all REPLACE_GENE_* tokens in the custom GMT before scoring.", call. = FALSE)
  }
  sets
}

hallmark_names <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)
hallmark_table <- msigdbr::msigdbr(
  species = cfg$gsva$hallmark_species,
  collection = "H"
)
symbol_column <- if ("gene_symbol" %in% names(hallmark_table)) {
  "gene_symbol"
} else {
  "db_gene_symbol"
}
hallmark_table <- hallmark_table[hallmark_table$gs_name %in% hallmark_names, , drop = FALSE]
hallmark_sets <- split(hallmark_table[[symbol_column]], hallmark_table$gs_name)
hallmark_sets <- lapply(hallmark_sets, unique)

custom_sets <- read_gmt(cfg$inputs$custom_gmt)
if (!any(grepl("JUNCTION", names(custom_sets), ignore.case = TRUE))) {
  stop("Custom GMT must include at least one set with JUNCTION in its name.", call. = FALSE)
}
gene_sets <- c(hallmark_sets, custom_sets)
overlap <- vapply(gene_sets, function(x) sum(x %in% rownames(expr)), integer(1))
if (any(overlap < cfg$gsva$min_size)) {
  message(
    "Sets below min_size will be excluded: ",
    paste(names(overlap)[overlap < cfg$gsva$min_size], collapse = ", ")
  )
}

method <- cfg$gsva$method
param <- switch(
  method,
  ssgsea = GSVA::ssgseaParam(
    exprData = expr,
    geneSets = gene_sets,
    minSize = cfg$gsva$min_size,
    maxSize = cfg$gsva$max_size,
    normalize = cfg$gsva$ssgsea_normalize
  ),
  gsva = GSVA::gsvaParam(
    exprData = expr,
    geneSets = gene_sets,
    minSize = cfg$gsva$min_size,
    maxSize = cfg$gsva$max_size,
    kcdf = cfg$gsva$kcdf
  ),
  stop("gsva.method must be ssgsea or gsva.", call. = FALSE)
)
scores <- GSVA::gsva(param, verbose = FALSE)

out_dir <- cfg$analysis$output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  cbind(signature = rownames(scores), as.data.frame(scores, check.names = FALSE)),
  file.path(out_dir, paste0(method, "_scores.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(scores, file.path(out_dir, paste0(method, "_scores.rds")))
write.table(
  data.frame(signature = names(overlap), detected_genes = overlap),
  file.path(out_dir, paste0(method, "_gene_set_overlap.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "gsva_sessionInfo.txt"))

#!/usr/bin/env Rscript
# Extract selected genes from the open GSE205335 UMI RDS.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: extract_gse205335_genes.R <rds.gz> <out.tsv>")
}
rds_path <- args[[1]]
out_path <- args[[2]]

genes <- c(
  "TACSTD2", "CLDN4",
  "HLA-A", "HLA-B", "HLA-C", "B2M", "STAT1", "IRF1", "IFITM1", "CXCL10"
)

message("Reading ", rds_path)
if (grepl("\\.gz$", rds_path, ignore.case = TRUE)) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  status <- system2("gzip", c("-dc", rds_path), stdout = tmp)
  if (!identical(status, 0L)) {
    stop("Failed to decompress ", rds_path)
  }
  obj <- readRDS(tmp)
} else {
  obj <- readRDS(rds_path)
}
message("Class: ", paste(class(obj), collapse = "/"))

as_matrix <- function(x) {
  if (inherits(x, "dgCMatrix") || inherits(x, "Matrix")) {
    return(x)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    return(as.matrix(x))
  }
  NULL
}

counts <- NULL
if (!is.null(as_matrix(obj))) {
  counts <- as_matrix(obj)
} else if (isS4(obj) && "assays" %in% slotNames(obj)) {
  assays <- slot(obj, "assays")
  if (is.list(assays) || is(assays, "SimpleList")) {
    assay <- assays[[1]]
    if (isS4(assay) && "counts" %in% slotNames(assay)) {
      counts <- slot(assay, "counts")
    } else if (is.list(assay) && "counts" %in% names(assay)) {
      counts <- assay$counts
    }
  }
} else if (is.list(obj)) {
  for (name in c("counts", "data", "RNA")) {
    if (name %in% names(obj) && !is.null(as_matrix(obj[[name]]))) {
      counts <- as_matrix(obj[[name]])
      break
    }
  }
}

if (is.null(counts)) {
  stop("Could not locate a count matrix in the RDS object")
}

message("Matrix dim: ", paste(dim(counts), collapse = " x "))
gene_names <- rownames(counts)
if (is.null(gene_names)) {
  stop("Count matrix has no rownames")
}

present <- intersect(genes, gene_names)
missing <- setdiff(genes, gene_names)
if (length(missing)) {
  message("Missing genes: ", paste(missing, collapse = ", "))
}
if (!length(present)) {
  stop("None of the requested genes are present")
}

selected <- as.matrix(counts[present, , drop = FALSE])
out <- data.frame(barcode = colnames(counts), t(selected), check.names = FALSE)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out_path, " genes=", paste(present, collapse = ","))

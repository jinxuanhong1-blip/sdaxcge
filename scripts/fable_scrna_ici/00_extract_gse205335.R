#!/usr/bin/env Rscript
# Extract per-cell panel-gene UMI counts + per-cell total UMI from the
# GSE205335 processed UMI matrix (dgCMatrix in RDS).
# Output: /tmp/data_scrna_ici/gse205335_panel_cells.tsv.gz

suppressPackageStartupMessages(library(Matrix))

data_dir <- "/tmp/data_scrna_ici"
# NOTE: the GEO file is double-gzipped; decompress once first:
#   zcat GSE205335_Lung_IO_UMI_matrix.rds.gz > GSE205335_Lung_IO_UMI_matrix.rds
rds_file <- file.path(data_dir, "GSE205335_Lung_IO_UMI_matrix.rds")
panel_file <- "/workspace/scripts/fable_scrna_ici/gene_panel.tsv"
out_file <- file.path(data_dir, "gse205335_panel_cells.tsv.gz")

m <- readRDS(rds_file)
cat("class:", class(m), "dim:", dim(m)[1], "x", dim(m)[2], "\n")

panel <- read.delim(panel_file, stringsAsFactors = FALSE)$gene
present <- intersect(panel, rownames(m))
missing <- setdiff(panel, rownames(m))
cat("panel genes present:", length(present), "missing:", paste(missing, collapse = ","), "\n")

total_umi <- Matrix::colSums(m)
sub <- as.matrix(m[present, , drop = FALSE])

df <- data.frame(barcode = colnames(m), total_umi = total_umi,
                 t(sub), check.names = FALSE)
gz <- gzfile(out_file, "w")
write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
cat("wrote", out_file, "cells:", nrow(df), "\n")

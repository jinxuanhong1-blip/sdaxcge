#!/usr/bin/env Rscript
# Extract raw counts (for genes of interest) + full clinical pData from the
# authoritative IMvigor210CoreBiologies CountDataSet `cds`.
#
# Input : data/raw/imvigor210_cds.RData  (mirror: github.com/SiYangming/
#         IMvigor210CoreBiologies, file data/cds.RData; original package:
#         Mariathasan et al., Nature 2018)
# Output: data/raw/imvigor210_counts_selected.csv  (samples x genes, raw counts)
#         data/raw/imvigor210_pdata_full.csv        (samples x clinical vars)
#
# Slots are read directly (base R) so DESeq/Biobase need not be installed; the
# `namespace 'DESeq' not available` warning on load is expected and harmless.

RAW <- "results/hunt_oak/data/raw"
load(file.path(RAW, "imvigor210_cds.RData"))          # -> object `cds`

counts <- get("counts", envir = attr(cds, "assayData"))  # Entrez rownames
pdat   <- attr(attr(cds, "phenoData"),   "data")
fdat   <- attr(attr(cds, "featureData"), "data")
stopifnot(all(rownames(counts) == rownames(fdat)))

genes <- c(TACSTD2 = 4070, CD274 = 29126, EPCAM = 4072,
           ACTB = 60, PTPRC = 5788, CD8A = 925)          # Entrez ids
sel <- t(counts[as.character(genes), , drop = FALSE])
colnames(sel) <- names(genes)

write.csv(data.frame(sample = rownames(sel), sel, check.names = FALSE),
          file.path(RAW, "imvigor210_counts_selected.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(pdat), pdat, check.names = FALSE),
          file.path(RAW, "imvigor210_pdata_full.csv"), row.names = FALSE)

cat("Extracted", nrow(sel), "samples;",
    "TACSTD2 median raw count =", median(sel[, "TACSTD2"]), "\n")

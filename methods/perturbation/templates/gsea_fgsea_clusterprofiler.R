#!/usr/bin/env Rscript
# =============================================================================
# gsea_fgsea_clusterprofiler.R  -  pre-ranked GSEA in R
# -----------------------------------------------------------------------------
# INPUT : a .rnk file (gene <tab> rank_metric) produced by the DESeq2/edgeR
#         templates. Rank by shrunken LFC or signed -log10(p).
# SETS  : MSigDB Hallmark (H) via msigdbr, focusing on:
#           HALLMARK_INTERFERON_ALPHA_RESPONSE
#           HALLMARK_INTERFERON_GAMMA_RESPONSE
#           HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#         plus a tight-junction set from C2:CP:KEGG / GO:CC.
# USAGE : Rscript gsea_fgsea_clusterprofiler.R ranked.rnk human out_dir
#         species: "human" (Homo sapiens) or "mouse" (Mus musculus)
# =============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
rnk_file <- ifelse(length(args) >= 1, args[1], "ranked.rnk")
species  <- ifelse(length(args) >= 2, args[2], "human")
out_dir  <- ifelse(length(args) >= 3, args[3], "gsea_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp <- if (species == "mouse") "Mus musculus" else "Homo sapiens"

# ---- ranked vector ---------------------------------------------------------
rnk <- fread(rnk_file, header = FALSE, col.names = c("gene", "metric"))
rnk <- rnk[is.finite(metric)]
rnk <- rnk[order(-metric)]
ranks <- setNames(rnk$metric, rnk$gene)
# collapse duplicate symbols (keep the most extreme metric)
ranks <- tapply(ranks, names(ranks), function(x) x[which.max(abs(x))])
ranks <- sort(ranks, decreasing = TRUE)

# ---- gene sets -------------------------------------------------------------
hall <- msigdbr(species = sp, category = "H")
want_h <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
tj <- msigdbr(species = sp, category = "C2", subcategory = "CP:KEGG")
tj <- tj[grepl("TIGHT_JUNCTION", tj$gs_name), ]

msig <- rbind(hall[hall$gs_name %in% want_h, ], tj)
pathways <- split(msig$gene_symbol, msig$gs_name)
pathways <- lapply(pathways, unique)
message("Gene sets tested: ", paste(names(pathways), collapse = ", "))

# ---- fgsea -----------------------------------------------------------------
set.seed(42)
fg <- fgsea(pathways = pathways, stats = ranks, minSize = 5, maxSize = 1000, eps = 0)
fg <- fg[order(fg$padj)]
fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ",")
fwrite(fg, file.path(out_dir, "fgsea_results.tsv"), sep = "\t")
print(fg[, c("pathway", "NES", "pval", "padj", "size")])
message("Done. Outputs in: ", normalizePath(out_dir))

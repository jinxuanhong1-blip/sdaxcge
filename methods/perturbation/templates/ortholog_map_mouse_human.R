#!/usr/bin/env Rscript
# =============================================================================
# ortholog_map_mouse_human.R  -  mouse <-> human symbol / ortholog harmonization
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   * Gene SYMBOLS differ by case and history between species:
#       human CLDN4 / TACSTD2 / EPCAM   <->   mouse Cldn4 / Tacstd2 / Epcam
#     A naive toupper() works for MANY genes but SILENTLY fails for the rest
#     (e.g. human C1orf* vs mouse, human ORFs, and true 1:many orthologs).
#   * Most curated gene sets (MSigDB Hallmark, KEGG) are HUMAN-symbol based, so
#     to run GSEA on a mouse experiment you must map mouse -> human orthologs.
#   * DEPRECATED/alias symbols: harmonize to current symbols before merging
#     tables from different studies (e.g. TROP2 -> TACSTD2, older claudin names).
#
# TWO SUPPORTED ROUTES (pick one):
#   A) babelgene  (offline, bundled ortholog tables; fastest & reproducible)
#   B) biomaRt    (queries Ensembl; needs network but always current)
#
#   Rscript ortholog_map_mouse_human.R genes.txt mouse2human out.tsv
#     direction: mouse2human | human2mouse
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
genes_file <- ifelse(length(args) >= 1, args[1], "genes.txt")
direction  <- ifelse(length(args) >= 2, args[2], "mouse2human")
out_file   <- ifelse(length(args) >= 3, args[3], "orthologs.tsv")

genes <- readLines(genes_file)
genes <- trimws(genes[nzchar(genes)])

# ---- ROUTE A: babelgene (recommended default) -----------------------------
route_babelgene <- function(genes, direction) {
  library(babelgene)
  if (direction == "mouse2human") {
    o <- orthologs(genes = genes, species = "mouse", human = FALSE)
    data.frame(input = o$symbol, mapped = o$human_symbol,
               human_entrez = o$human_entrez, support = o$support_n)
  } else {
    o <- orthologs(genes = genes, species = "mouse", human = TRUE)
    data.frame(input = o$human_symbol, mapped = o$symbol,
               mouse_entrez = o$entrez, support = o$support_n)
  }
}

# ---- ROUTE B: biomaRt (network) -------------------------------------------
route_biomart <- function(genes, direction) {
  library(biomaRt)
  hs <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
  mm <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl")
  if (direction == "mouse2human") {
    getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = genes,
           mart = mm, attributesL = "hgnc_symbol", martL = hs)
  } else {
    getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = genes,
           mart = hs, attributesL = "mgi_symbol", martL = mm)
  }
}

res <- tryCatch(route_babelgene(genes, direction),
                error = function(e) {
                  message("babelgene failed (", conditionMessage(e),
                          "); falling back to biomaRt.")
                  route_biomart(genes, direction)
                })

write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Mapped %d input genes -> %d rows. Wrote %s",
                length(genes), nrow(res), out_file))

# ---- NOTE on the quick-and-dirty path -------------------------------------
# For a fast eyeball check only (NOT for publication tables):
#   human_like <- toupper(mouse_symbols)   # Cldn4 -> CLDN4
# This is acceptable for a sanity glance at a handful of focus genes, but use a
# real ortholog table for anything that feeds statistics or figures.

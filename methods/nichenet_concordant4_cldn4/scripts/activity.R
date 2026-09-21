#!/usr/bin/env Rscript
# Official NicheNet ligand activity.
# Sources saeyslab/nichenetr @ 66f90d5eeafef280b2b2f339b3fd70ffec1781dd
# predict_ligand_activities -> get_single_ligand_importances -> evaluate_target_prediction
# -> classification_evaluation_continuous_pred (ROCR + caTools::trapz).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(magrittr)
  library(ROCR)
  library(caTools)
})

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}
IN <- parse_opt("--in", "/tmp/nichenet_work/activity_in")
OUT <- parse_opt("--out", "/tmp/nichenet_work/activity_out")
SRC <- parse_opt("--src", "/tmp")
PRIOR <- parse_opt("--prior", "/tmp/nichenet_prior/ligand_target_matrix_nsga2r_final.rds")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(SRC, "supporting_functions.R"))
source(file.path(SRC, "evaluate_model_target_prediction.R"))
source(file.path(SRC, "evaluate_model_ligand_prediction.R"))

predict_ligand_activities <- function(geneset, background_expressed_genes, ligand_target_matrix, potential_ligands) {
  setting <- list(geneset) %>%
    lapply(convert_gene_list_settings_evaluation, name = "gene set", ligands_oi = potential_ligands, background = background_expressed_genes)
  settings_ligand_prediction <- setting %>%
    convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = TRUE)
  ligand_importances <- settings_ligand_prediction %>%
    lapply(get_single_ligand_importances, ligand_target_matrix = ligand_target_matrix, known = FALSE) %>%
    bind_rows()
  ligand_importances %>% select(test_ligand, auroc, aupr, aupr_corrected, pearson)
}

logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", sep = " ")

logmsg("load ligand-target matrix")
lt <- readRDS(PRIOR)
stopifnot(is.matrix(lt))
storage.mode(lt) <- "double"

background <- readLines(file.path(IN, "background_genes.txt"))
potential <- readLines(file.path(IN, "potential_ligands.txt"))
background <- intersect(background, rownames(lt))
potential <- intersect(potential, colnames(lt))
logmsg("background", length(background), "potential", length(potential))
lt_sub <- lt[background, potential, drop = FALSE]
rm(lt); gc(verbose = FALSE)

sets <- list.files(IN, pattern = "^geneset_.*\\.txt$", full.names = TRUE)
if (!length(sets)) stop("no genesets")

for (fp in sets) {
  nm <- sub("^geneset_", "", sub("\\.txt$", "", basename(fp)))
  geneset <- intersect(readLines(fp), rownames(lt_sub))
  logmsg("geneset", nm, "n", length(geneset))
  if (length(geneset) < 5) {
    logmsg("  skip, fewer than 5 genes in matrix")
    write.table(
      data.frame(geneset = nm, test_ligand = character(), auroc = numeric(),
                 aupr = numeric(), aupr_corrected = numeric(), pearson = numeric()),
      file.path(OUT, paste0("activity_", nm, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    next
  }
  pred <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = lt_sub,
    potential_ligands = potential
  )
  pred$geneset <- nm
  pred$rank_aupr <- rank(-pred$aupr_corrected, ties.method = "min")
  pred$rank_pearson <- rank(-pred$pearson, ties.method = "min")
  pred$n_geneset <- length(geneset)
  pred$n_background <- nrow(lt_sub)
  pred$n_potential <- length(potential)
  write.table(pred, file.path(OUT, paste0("activity_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  logmsg("  wrote", nrow(pred))
}

# Continuous companion: Pearson of regulatory potential vs a named numeric vector
cont <- file.path(IN, "continuous_response.tsv")
if (file.exists(cont)) {
  resp <- read.delim(cont, stringsAsFactors = FALSE)
  resp <- resp[resp$gene %in% rownames(lt_sub), , drop = FALSE]
  y <- resp$z
  names(y) <- resp$gene
  y <- y[rownames(lt_sub)]
  logmsg("continuous response genes", sum(is.finite(y)))
  rows <- lapply(potential, function(lig) {
    x <- lt_sub[, lig]
    ok <- is.finite(x) & is.finite(y)
    data.frame(
      test_ligand = lig,
      pearson = suppressWarnings(cor(x[ok], y[ok], method = "pearson")),
      spearman = suppressWarnings(cor(x[ok], y[ok], method = "spearman")),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows)
  out$geneset <- "continuous_meta_z"
  write.table(out, file.path(OUT, "activity_continuous_meta_z.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  logmsg("continuous wrote", nrow(out))
}
logmsg("activity done")

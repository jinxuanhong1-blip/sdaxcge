#!/usr/bin/env Rscript
# %% [markdown]
# 03 - NicheNet: which epithelial ligands best explain a T-cell response program
#
# READ THIS FIRST (scope). NicheNet links a sender LIGAND to a *transcriptional*
# response in the receiver via a generic prior model (ligand -> signaling -> target
# genes). Chemokine-driven "recruitment" is chemotaxis/migration, which is largely
# NOT a transcriptional target program -- so NicheNet is the WRONG primary tool for
# "reduced recruitment" and will happily rank chemokines low or high for the wrong
# reason. Use it for a *different, legitimate* sub-question: do TACSTD2/CLDN4-high
# epithelia send ligands that reshape the transcriptional STATE of infiltrating
# T cells (exhaustion, dysfunction, activation)? Keep recruitment claims to
# LIANA/CellChat + the pseudobulk guardrail. See playbook §NicheNet overclaims.
#
# What NicheNet does / does not do. Ligand "activity" = how well a ligand's
# prior-predicted targets match the receiver's DE gene set (AUPR vs background).
# By default it does NOT require the ligand or its receptor to be expressed --
# you must add expression filters. The prior model is tissue/context-agnostic.
# Results are correlational hypotheses, not causal effects.

# %% ---- setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(nichenetr)   # v2
  library(Seurat)
  library(tidyverse)
  library(yaml)
})
cfg <- tryCatch(yaml::read_yaml("../config/gene_sets.yaml"),
                error = function(e) yaml::read_yaml("methods/grn_ccc/config/gene_sets.yaml"))
P <- cfg$params

IN_RDS <- "data/lung_ici_preprocessed.rds"
OUTDIR <- "results/nichenet"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# NicheNet v2 human prior model (download once from Zenodo 7074291):
#   ligand_target_matrix_nsga2r_final.rds, lr_network_human_21122021.rds,
#   weighted_networks_nsga2r_final.rds
lt_path <- Sys.getenv("NICHENET_LT", "networks/ligand_target_matrix_nsga2r_final.rds")
lr_path <- Sys.getenv("NICHENET_LR", "networks/lr_network_human_21122021.rds")
wn_path <- Sys.getenv("NICHENET_WN", "networks/weighted_networks_nsga2r_final.rds")
ligand_target_matrix <- readRDS(lt_path)
lr_network  <- readRDS(lr_path) %>% distinct(from, to)
weighted_networks <- readRDS(wn_path)

obj <- readRDS(IN_RDS)
Idents(obj) <- P$celltype_key

# %% [markdown]
# ## 1. Define receiver, senders, and the response gene set (the crux)
# receiver = T cells. geneset_oi = genes DE in T cells that reside in a
# TACSTD2/CLDN4-HIGH niche vs a LOW niche (compute per patient, then pool).
# GARBAGE IN / GARBAGE OUT: a sloppy DE set makes every downstream number
# meaningless. Use a proper contrast with `min.pct`/logFC and, ideally,
# patient-level pseudobulk DE for the gene set too.

# %%
receiver <- "T_cell"
senders  <- c("Epithelial")   # split further to high/low below via subsets

expressed_genes_receiver <- get_expressed_genes(receiver, obj, pct = 0.10)
background_expressed_genes <- expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target_matrix)]

# DE in receiver between niches (requires obj$epi_niche in {high,low} per T cell,
# derived from the dominant neighboring epithelial state; or use condition_key).
Tsub <- subset(obj, idents = receiver)
Idents(Tsub) <- Tsub@meta.data[["epi_niche"]]   # "high"/"low"; EDIT to your design
de <- FindMarkers(Tsub, ident.1 = "high", ident.2 = "low",
                  min.pct = 0.10, logfc.threshold = 0.25) %>% rownames_to_column("gene")
geneset_oi <- de %>% filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]
message(length(geneset_oi), " response genes define the T-cell program.")

# %% [markdown]
# ## 2. Candidate ligands = expressed by sender AND receptor expressed by receiver
# This expression filter is what keeps NicheNet honest. We compute it separately
# for the HIGH and LOW epithelial states to see whether the high state loses
# credible ligands.

# %%
ligands   <- unique(lr_network$from)
receptors <- unique(lr_network$to)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)
potential_from_receptor <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

run_state <- function(state_label) {
  Esub <- subset(obj, subset = epi_state == paste0("TACSTD2_CLDN4_", state_label))
  expr_lig <- get_expressed_genes("Epithelial", Esub, pct = 0.10)
  potential_ligands <- intersect(intersect(ligands, expr_lig), potential_from_receptor)
  act <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands) %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(state = state_label,
           is_recruitment_chemokine = test_ligand %in%
             unlist(lapply(cfg$t_cell_recruitment, function(x) x$ligands)))
  act
}

act.high <- run_state("high")
act.low  <- run_state("low")
bind_rows(act.high, act.low) %>%
  write_csv(file.path(OUTDIR, "nichenet_ligand_activities_by_state.csv"))
print(head(act.high, 15))
print(head(act.low, 15))

# %% [markdown]
# ## 3. Ligand-target links for the top ligands (hypothesis map, not proof)

# %%
top_ligands <- act.low %>% top_n(20, aupr_corrected) %>% pull(test_ligand)
lt_links <- top_ligands %>%
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix, n = 100) %>%
  bind_rows()
write_csv(lt_links, file.path(OUTDIR, "nichenet_ligand_target_links.csv"))

message("[done] NicheNet outputs in ", OUTDIR)

# %% [markdown]
# ## Interpretation guardrails
# - A top-ranked ligand is a *hypothesis*: verify it is actually expressed in the
#   sender state (we filtered, but re-check fraction + level) and that its receptor
#   is on T cells.
# - Ligand activity depends entirely on `geneset_oi`; re-run with a perturbed gene
#   set to check stability. If the ranking swings wildly, do not report it.
# - Do NOT phrase NicheNet output as "epithelium causes T-cell state X". It is a
#   prior-model prediction consistent with the data.

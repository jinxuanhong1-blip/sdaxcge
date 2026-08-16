#!/usr/bin/env Rscript
# %% [markdown]
# 02 - CellChat: pathway-level cell-cell communication
#
# Question: is the CXCR3 / CCR5 / CXCR6 recruitment signaling from
# TACSTD2/CLDN4-high epithelium to T/NK cells weaker than from the low state?
#
# What CellChat does / does not do. CellChat scores communication with a
# law-of-mass-action + Hill "communication probability" using the curated
# CellChatDB, then aggregates ligand-receptor pairs into signaling PATHWAYS.
# The "probability" is a heuristic score, NOT a physical probability, and is
# sensitive to `type` (triMean vs truncatedMean), `population.size`, and the
# `trim` fraction. Its permutation test is per-object label-shuffling, so
# comparing two CellChat objects is descriptive, not a cross-patient test.
# See playbook §CellChat overclaims.
#
# Run: Rscript 02_cellchat.R   (needs a Seurat/SCE export of 00_preprocess.py)

# %% ---- setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(yaml)
  library(Matrix)
})
options(stringsAsFactors = FALSE)

cfg <- yaml::read_yaml(file.path(dirname(sys.frame(1)$ofile %||% "."),
                                 "..", "config", "gene_sets.yaml"))
# Fallback if run interactively:
if (is.null(cfg)) cfg <- yaml::read_yaml("../config/gene_sets.yaml")
P <- cfg$params

IN_RDS  <- "data/lung_ici_preprocessed.rds"   # Seurat obj with cc_group + patient_id
OUTDIR  <- "results/cellchat"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## 1. Load and build the per-condition input
# We build ONE CellChat object per state (high / low) so we can compare them.
# `data.input` must be log-normalized data (NOT batch-corrected), `meta` must
# carry `cc_group` (senders/receivers, with Epi split into high/low).

# %%
obj <- readRDS(IN_RDS)                          # e.g. SeuratObject
data.input <- as(GetAssayData(obj, slot = "data"), "dgCMatrix")  # log-norm
meta <- obj@meta.data
meta$cc_group <- as.character(meta$cc_group)

make_cellchat <- function(cells) {
  cc <- createCellChat(object = data.input[, cells],
                       meta = meta[cells, , drop = FALSE],
                       group.by = "cc_group")
  cc@DB <- CellChatDB.human                     # curated, literature-biased
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # KEY KNOBS (report them): triMean is conservative (~>=25% expressing);
  # population.size=FALSE avoids abundance inflating probabilities.
  cc <- computeCommunProb(cc, type = "triMean", population.size = FALSE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# Compare high vs low senders in the SAME tissue context:
cc.high <- make_cellchat(rownames(meta)[meta$cc_group != "Epi_TACSTD2low"])
cc.low  <- make_cellchat(rownames(meta)[meta$cc_group != "Epi_TACSTD2high"])

# %% [markdown]
# ## 2. Extract LR-level table and keep the recruitment axes only

# %%
recr_ligands <- unique(unlist(lapply(cfg$t_cell_recruitment,
                                     function(x) x$ligands)))

get_focus <- function(cc, sender) {
  df <- subsetCommunication(cc)                 # data.frame of significant LR
  df <- df[df$source == sender &
             grepl("T_cell|NK|Tcell", df$target) &
             df$ligand %in% recr_ligands, ]
  df[order(-df$prob), ]
}

focus.high <- get_focus(cc.high, "Epi_TACSTD2high")
focus.low  <- get_focus(cc.low,  "Epi_TACSTD2low")
write.csv(focus.high, file.path(OUTDIR, "cellchat_focus_high.csv"), row.names = FALSE)
write.csv(focus.low,  file.path(OUTDIR, "cellchat_focus_low.csv"),  row.names = FALSE)
print(focus.high[, c("source","target","ligand","receptor","prob","pval")])
print(focus.low [, c("source","target","ligand","receptor","prob","pval")])

# %% [markdown]
# ## 3. Formal comparison of the two objects (descriptive only)
# mergeCellChat + rankNet contrasts pathway information flow. This shuffles
# labels within each object; it does NOT test across patients. Frame as
# "the high object shows lower CXCL/CCL information flow", then confirm with
# the patient-level pseudobulk guardrail (05_pseudobulk_guardrail.py).

# %%
object.list <- list(low = cc.low, high = cc.high)
merged <- mergeCellChat(object.list, add.names = names(object.list))

pdf(file.path(OUTDIR, "cellchat_rankNet_information_flow.pdf"), width = 6, height = 8)
tryCatch(
  print(rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)),
  error = function(e) message("rankNet failed: ", conditionMessage(e)))
dev.off()

# Pathway-level flow for the recruitment chemokine families (CXCL, CCL):
saveRDS(object.list, file.path(OUTDIR, "cellchat_object_list.rds"))
message("[done] CellChat outputs in ", OUTDIR)

# %% [markdown]
# ## Interpretation guardrails
# - Report `type`, `population.size`, `trim`, `min.cells`; results move with them.
# - A pathway "present in low, absent in high" can reflect detection/dropout at
#   low chemokine expression, not biology. Cross-check raw fraction-expressing.
# - CellChatDB omits many interactions; absence is not evidence of absence.

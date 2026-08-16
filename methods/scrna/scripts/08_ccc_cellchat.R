## SECONDARY analysis: cell-cell communication with CellChat v2.
## Hypothesis-generating only. For group comparison, build one CellChat object per
## condition (e.g. MPR / NMPR) and use compareInteractions / netVisual_diffInteraction.
## 仅作次要/假设生成；分组比较需每组各建对象后对比。见 playbook.md 第8节。
suppressPackageStartupMessages({library(CellChat); library(Seurat)})

args <- commandArgs(trailingOnly = TRUE)
in_rds   <- args[[1]]     # Seurat with 'lineage' (or finer) in meta.data
out_rds  <- args[[2]]
group_by <- ifelse(length(args) >= 3, args[[3]], "lineage")

seu <- readRDS(in_rds)
data_input <- GetAssayData(seu, assay = "RNA", slot = "data")  # log-normalized
meta <- data.frame(labels = seu@meta.data[[group_by]], row.names = colnames(seu))

cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
cc@DB <- CellChatDB.human
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, type = "triMean")
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

saveRDS(cc, out_rds)
message("CellChat object -> ", out_rds,
        " (use netVisual_* / compareInteractions for group contrasts)")

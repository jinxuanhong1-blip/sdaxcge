## Malignant vs normal epithelial via inferCNV. Reference = immune/stromal cells.
## Malignant cells show broad CNV and (correctly) do NOT integrate across patients.
## 用免疫/基质细胞作参考，inferCNV 判定恶性上皮；恶性细胞的患者特异性来自CNV而非批次。
suppressPackageStartupMessages({library(infercnv)})

args <- commandArgs(trailingOnly = TRUE)
counts_mat  <- args[[1]]   # genes x cells raw counts (matrix or 10x dir)
annot_file  <- args[[2]]   # 2-col: cell, group (e.g. "epithelial","Tcell","myeloid"...)
gene_order  <- args[[3]]   # gene ordering file (chr, start, end) for the reference genome
out_dir     <- args[[4]]

ref_groups <- c("Tcell", "NK", "Bcell", "myeloid", "fibroblast", "endothelial")

obj <- CreateInfercnvObject(raw_counts_matrix = counts_mat,
                            annotations_file  = annot_file,
                            gene_order_file   = gene_order,
                            ref_group_names   = ref_groups)

obj <- infercnv::run(obj,
                     cutoff = 0.1,                 # 0.1 for 10x, 1 for smart-seq
                     out_dir = out_dir,
                     cluster_by_groups = TRUE,
                     denoise = TRUE,
                     HMM = TRUE)
## Downstream: cells in epithelial clusters with high CNV burden = malignant.
## Combine HMM state / CNV score with a threshold, or use CopyKAT/numbat as cross-check.
message("inferCNV done -> ", out_dir)

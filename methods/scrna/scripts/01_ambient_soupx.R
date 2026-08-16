## Ambient-RNA removal with SoupX (fallback when CellBender GPU unavailable).
## Needs BOTH the filtered and raw 10x matrices for one lane. Run per lane.
## 需要同一lane的 filtered 与 raw 矩阵；每个lane单独运行。
suppressPackageStartupMessages({library(SoupX); library(Seurat); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
# args: <cellranger_outs_dir> <out_rds>
outs_dir <- args[[1]]      # contains filtered_feature_bc_matrix/ and raw_feature_bc_matrix/
out_rds  <- args[[2]]

sc <- load10X(outs_dir)                 # auto-loads filtered + raw
sc <- autoEstCont(sc)                   # estimate contamination fraction (inspect the plot)
adj <- adjustCounts(sc, roundToInt = TRUE)

## Sanity check: contamination fraction should be modest (e.g. < 0.2 typically).
message("Estimated rho (contamination): ",
        round(mean(sc$metaData$rho, na.rm = TRUE), 3))

seu <- CreateSeuratObject(counts = adj)
saveRDS(seu, out_rds)
message("SoupX-corrected Seurat object -> ", out_rds)

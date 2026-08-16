# Environment notes / 环境说明

Two environments: a **Python** stack (Scanpy/scvi-tools/decoupler/LIANA) and an **R**
stack (Seurat/edgeR/DESeq2/UCell/CellChat/scDblFinder). Pin versions and set seeds.

固定版本、设定随机种子。Python 与 R 两套环境分开管理。

## Python (conda) / Python 环境

```bash
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate scrna-ici
python -c "import scanpy, scvi, decoupler; print(scanpy.__version__)"
```

Key packages: `scanpy`, `anndata`, `scvi-tools`, `harmonypy`, `celltypist`,
`decoupler`, `liana`, `scikit-misc`, `leidenalg`, `igraph`.
GPU strongly recommended for scVI/scANVI (install the CUDA build of `jax`/`torch`).

> The `demo/` in this repo only needs the lighter subset:
> `scanpy anndata pandas scipy matplotlib openpyxl harmonypy igraph leidenalg scikit-misc`.
> Harmony (not scVI) is used in the demo to keep it CPU-friendly.

## R (renv / conda-forge) / R 环境

```r
# option A: renv
renv::init()
renv::restore()      # from renv.lock (generate on your machine)

# option B: install.packages + Bioconductor
install.packages(c("Seurat", "harmony", "Matrix", "dplyr", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scDblFinder", "edgeR", "DESeq2", "UCell", "AUCell",
                       "CellChat", "SingleCellExperiment", "scran", "scater",
                       "infercnv", "muscat", "decoupleR", "liana"))
```

Key packages: `Seurat` (v5), `harmony`, `scDblFinder`, `edgeR`, `DESeq2`, `muscat`,
`UCell`, `AUCell`, `CellChat` (v2), `liana`, `infercnv`/`copykat`, `scran`/`scater`.

## Reproducibility / 可复现

- Record versions into every output (`sessionInfo()` / `sc.logging.print_versions()`).
- Set seeds: numpy/scanpy `random_state=0`, `set.seed(1)` in R, `scvi.settings.seed=0`.
- Save `.h5ad`/`.rds`, the pseudobulk matrix, and the sample design table as artifacts.

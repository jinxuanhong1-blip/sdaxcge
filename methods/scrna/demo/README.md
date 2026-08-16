# Demo: GSE207422 NSCLC neoadjuvant ICI scRNA-seq

**Dataset / 数据集:** [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422)
(Hu et al., *Genome Medicine* 15:14, 2023). 15 stage-IIIA NSCLC patients,
neoadjuvant PD-1 + chemo. Labels: **pCR / MPR / NMPR** and **RECIST**.

**What this is / 这是什么:** a *minimal, honest* run of the downstream half of
`../playbook.md` on the authors' **post-QC** UMI matrix (92,330 cells × 24,292
genes; 175 MB gzipped, < 2 GB). Ambient RNA and doublet removal were already
applied upstream; we do **not** re-invent those numbers. All figures/tables
are computed from the downloaded matrix. **No fabricated statistics.**
Small-n MPR vs NMPR comparisons are labeled **exploratory**.

**GSE205335** (Kim et al., *eLife* 2024; 499 MB processed `.rds`) is documented
in the playbook as an external replication cohort and was not re-run here.

## How to run / 如何运行

```bash
# 1. download processed files from GEO (< 2 GB)
mkdir -p /tmp/geo && cd /tmp/geo
curl -fsSLO https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz
curl -fsSLO https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx

# 2. Python deps (see ../env/)
pip install scanpy anndata pandas scipy matplotlib openpyxl harmonypy igraph leidenalg scikit-misc pydeseq2

# 3. downstream playbook (Harmony, annotation, module scores, pseudobulk export)
cd /path/to/methods/scrna/demo
python run_demo.py \
    --matrix /tmp/geo/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
    --meta   /tmp/geo/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
    --outdir .
# default: stratified subsample, max 1,200 cells/sample (CPU-friendly).
# all cells: add  --max-per-sample 0

# 4. sample-level DESeq2 on epithelial pseudobulk (MPR-like vs NMPR)
python pseudobulk_de.py
# optional R equivalent (needs edgeR + DESeq2):  Rscript pseudobulk_de.R
```

## Outputs / 输出

| File | Content |
|---|---|
| `figures/01_qc_distributions.png` | genes / UMI / %MT (author-filtered matrix) |
| `figures/umap_02_lineage.png` | Leiden → marker-based lineage |
| `figures/umap_03_sample.png` | Harmony-integrated UMAP by sample |
| `figures/umap_04_*.png` | TACSTD2/CLDN4/junction, CD8, TLS, CXCL13 scores |
| `figures/05_epi_module_by_response.png` | epithelial module vs MPR/NMPR (descriptive) |
| `figures/06_de_goi_points.png` | per-sample TACSTD2 / CLDN4 CPM (exploratory DE) |
| `demo_summary.json` | n, lineage counts, integration method |
| `pseudobulk_epithelial_counts.csv` + `pseudobulk_design.csv` | sample × gene counts + MPR/RECIST |
| `de_epithelial_MPRlike_vs_NMPR.csv` | DESeq2 table (exploratory) |
| `de_epithelial_genes_of_interest.csv` | TACSTD2/CLDN4/junction genes only |

MPR and RECIST are kept as **separate** columns. pCR is collapsed into
`MPR_like` only for the two-group DE contrast. Pre-treatment / `NE` samples
are excluded from that contrast.

## Caveats / 注意

- This is a **methods demo**, not a discovery paper. n ≈ 4 MPR-like vs ≈ 8 NMPR.
- Wilcoxon-on-cells is **not** used as the primary DE. See `../playbook.md` §5.
- Cell–cell communication (LIANA/CellChat) is **not** run here (secondary only;
  templates live in `../scripts/08_*`).

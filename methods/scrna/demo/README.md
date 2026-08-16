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

## Results from this run (computed, not fabricated) / 本次运行的真实结果

Logged in `demo_summary.json` and the DE CSVs. Snapshot of the committed run:

- Matrix: **92,330 cells × 24,292 genes**, 15 samples, 100% matched to GEO metadata.
- Analyzed: **18,000 cells** (stratified, max 1,200 / sample). Integration: **Harmony**
  (`X_pca_harmony`, 18,000 × 30). **22** Leiden clusters.
- Lineage (marker argmax): T/NK 7,622; Myeloid 3,570; Epithelial 2,262; B/Plasma 2,117;
  Endothelial 2,077; Mast 177; Fibroblast 175.
- Junction module genes present: TACSTD2, CLDN4, CLDN3, CLDN7, CDH1, TJP1, OCLN, F11R
  (**CLDN18 absent** from this matrix). Score is highest on the epithelial cluster
  (`figures/umap_04_TACSTD2_CLDN4_junction.png`).
- **Epithelial pseudobulk DESeq2**, post-treatment only, ≥20 epithelial cells/sample:
  **4 MPR-like (pCR+MPR) vs 7 NMPR**. Exploratory; no GOI survives BH.
  | gene | log2FC (MPR-like vs NMPR) | p | padj |
  |---|---:|---:|---:|
  | CLDN7 | −0.94 | 0.042 | 0.76 |
  | TACSTD2 | −0.85 | 0.077 | 0.79 |
  | F11R | +0.69 | 0.082 | 0.79 |
  | CLDN4 | −0.45 | 0.22 | 0.89 |
  | CXCL13 | +1.68 | 0.30 | 0.93 |
- Per-sample epithelial module scores overlap between MPR and NMPR
  (`figures/05_epi_module_by_response.png`). Do **not** overclaim.

## Caveats / 注意

- This is a **methods demo**, not a discovery paper. n = 4 MPR-like vs 7 NMPR.
- Wilcoxon-on-cells is **not** used as the primary DE. See `../playbook.md` §5.
- Cell–cell communication (LIANA/CellChat) is **not** run here (secondary only;
  templates live in `../scripts/08_*`).
- Pre-treatment biopsies are excluded from the DE contrast (timepoint confound).

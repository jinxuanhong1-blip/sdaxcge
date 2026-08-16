# Target hunt: Tacstd2 (Trop2) in "cold" KL vs K lung tumors — GSE180963

**Question.** In the LKB1-deficient ("cold", immune-desert) KL lung tumor
microenvironment, is `Tacstd2`/Trop2 expressed selectively by the tumor
(epithelial) compartment relative to immune cells? Tumor-selective, immune-sparing
surface expression is the property that would make Trop2 attractive as a
tumor-directed target (e.g. the Trop2 antibody-drug conjugate sacituzumab
govitecan).

**Bottom line (honest).** Yes, within this dataset `Tacstd2`/Trop2 is essentially
restricted to the tumor epithelial compartment and near-absent from immune cells,
and this selectivity holds inside the cold KL tumor. But the dataset is weak for
any *KL-vs-K* claim, and it does **not** reproduce a lower total immune fraction in
KL. Treat the tumor-vs-immune selectivity as the real, robust result; treat
everything comparing the two genotypes as suggestive at best.

---

## Dataset

GSE180963 — "Single cell RNA sequencing of tumor sections from GEMM harboring
KrasG12D/+ or KrasG12D/+Lkb1fl/fl (KL) mutation" (Southern Medical University,
public 2022). Mouse (FVB), 10x Chromium 3' v3.1, aligned to mm10 by the authors.

| GSM | label | genotype | cells (author-filtered) |
|---|---|---|---|
| GSM5481386 | **K** | KrasG12D/+ | 6,696 |
| GSM5481387 | **KL** | KrasG12D/+ ; Lkb1fl/fl | 7,564 |

`Lkb1` = `Stk11`. **n = 1 mouse per genotype.** The two samples were pooled into a
single 10x library and demultiplexed by label after sequencing. The GEO
supplementary matrices are already QC-filtered by the authors (500–6000 genes/cell,
mito < 20%, genes in ≥3 cells; 20,304 genes); no raw/unfiltered matrices are
available.

## What was done

`analysis/run_hunt_gse180963.py`: load both matrices → QC metrics → CP10k
normalize + log1p → HVG(2000) → PCA(50) → neighbors → UMAP → Leiden(res=1.0) →
marker-based lineage annotation (per-cluster majority vote of per-cell lineage
scores) → `Tacstd2` quantification by cell type / compartment / genotype →
annotation-independent robustness checks. Deterministic (seed 0).

**No batch integration was performed, on purpose:** genotype is fully confounded
with sample (one K, one KL), so "correcting" batch would erase the very biology in
question. Annotation was therefore kept at the lineage level, which is robust to
this.

## Key result — KL (cold tumor): epithelial/tumor vs immune

Two independent definitions of "tumor cell" agree:

| definition | epithelial cells | % Trop2+ (epi) | % Trop2+ (immune) | mean log-norm epi | mean log-norm immune | log2FC | Mann–Whitney p |
|---|---|---|---|---|---|---|---|
| Leiden cluster = Epithelial/Tumor | 231 | **62.8%** | 3.6% | 0.577 | 0.039 | 3.91 | ≈0 |
| Marker rule (Epcam+ & keratin/Cldn18+ & Ptprc−) | 170 | **64.7%** | 4.3% | 0.614 | 0.044 | 3.79 | ≈0 |

`Tacstd2`/Trop2 is ~15–18× higher (mean) in KL tumor epithelial cells than in KL
immune cells, expressed in roughly two-thirds of tumor cells vs a few percent of
immune cells. The low immune "positivity" (~4%) is low-level and consistent with
ambient contamination / doublets rather than genuine immune expression (see the
UMAP feature plot — signal is confined to the epithelial cluster).

The same pattern holds in K (control): epithelial 40% (cluster) / 17.9% (marker)
Trop2+ vs ~2% in immune cells (log2FC ~3.0–3.3, p ≪ 1e-7), just with far fewer
epithelial cells to work with.

## KL vs K (read with heavy skepticism)

Trop2 positivity/level in the epithelial compartment was nominally higher in KL
than K (62.8% vs 40% by cluster; 64.7% vs 17.9% by marker rule; mean-expression
log2FC ~1.1). **This is not a trustworthy KL>K claim:** n = 1 mouse per genotype,
genotype fully confounded with sample and the single pooled library, and only
20–28 confidently-epithelial cells were recovered from K. Any p-value here is
cell-level pseudoreplication, not a test across biological replicates.

## Honesty caveats (please read)

1. **No biological replication.** One mouse per genotype. All cross-genotype
   comparisons are confounded with mouse and batch; reported p-values treat cells
   as independent replicates (pseudoreplication) and should be read as descriptive,
   not inferential.
2. **The "cold" phenotype is not visible as a lower immune fraction here.** Both
   samples are ~82% CD45/`Ptprc`+ (K 82.8%, KL 82.3%); overall immune fraction is
   essentially identical (K 81.8%, KL 81.7% by lineage). The immune-desert
   phenotype reported for LKB1 loss is compositional/functional (T-cell exclusion,
   suppressive myeloid skew), which this quick pass did **not** attempt to
   reproduce. We are not independently confirming that this KL sample is "cold".
3. **Tumor cells are rare and inconsistently recovered.** The dissociation is
   immune-dominated; confidently-epithelial cells were ~20–231 (K vs KL). K is
   especially thin (~20–28 cells), so K epithelial statistics are fragile.
4. **Ambient surfactant contamination.** `Sftpc` is "positive" in ~77% of all
   cells (including immune) — clearly ambient. We therefore did **not** rely on
   surfactant genes for epithelial calling and used `Epcam`/keratins/`Cldn18`
   instead. Ambient RNA likely also explains the small non-zero Trop2 signal in
   immune/stromal cells.
5. **Automated marker-based annotation**, not expert-reviewed. Cluster-to-lineage
   assignments are auditable in `tables/cluster_lineage_scores.csv`.
6. **Author-prefiltered input only** — we cannot redo cell calling or ambient
   correction (e.g. SoupX/CellBender) without the raw matrices.

## Files

Tables (`tables/`):
- `Tacstd2_contrasts.csv` — primary tumor-vs-immune contrasts (per genotype) + KL-vs-K epithelial.
- `Tacstd2_contrasts_annotation_independent.csv` — same conclusions using clustering-free definitions (Ptprc− and marker-rule epithelial).
- `Tacstd2_by_celltype_genotype.csv`, `Tacstd2_by_compartment_genotype.csv` — % expressing and mean expression.
- `celltype_counts_by_genotype.csv`, `celltype_fraction_by_genotype.csv`, `compartment_fraction_by_genotype.csv` — composition.
- `immune_dominance_ambient_qc.csv` — %Ptprc+, %Sftpc+ (ambient flag), epithelial counts.
- `cluster_lineage_scores.csv` — per-cluster marker scores + assignment (audit).
- `qc_summary.csv` — per-genotype QC distributions.

Figures (`figures/`):
- `umap_celltype.png`, `umap_compartment.png`, `umap_genotype.png`
- `umap_Tacstd2.png` — Trop2 confined to the epithelial cluster.
- `violin_Tacstd2_compartment_genotype.png`, `violin_Tacstd2_by_celltype.png`
- `dotplot__markers.png` — lineage marker sanity check.

Other: `metrics.json` (machine-readable summary). The processed AnnData
(`adata_gse180963_processed.h5ad`, ~0.8 GB) and the downloaded matrices under
`data/` are git-ignored; regenerate them by running the analysis script.

## Reproduce

```bash
pip install scanpy anndata leidenalg python-igraph pandas scipy numpy matplotlib seaborn
# download K/KL matrices from GEO into data/gse180963/{K,KL}/{K,KL}/{matrix.mtx,genes.tsv,barcodes.tsv}
python3 analysis/run_hunt_gse180963.py
```

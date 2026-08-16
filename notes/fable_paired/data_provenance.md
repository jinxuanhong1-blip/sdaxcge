# Data provenance & methods — fable_paired (ICI lung, TACSTD2/CLDN4)

All raw files are public GEO records, downloaded by `scripts/fable_paired/download_data.py`
into `$FABLE_RAW_DIR` (default `/tmp/fable_raw`, **not** committed). Only compact
processed tables/figures under `results/fable_paired/` are version controlled
(total ~3 MB, far under the 2 GB budget).

## Datasets

| Accession | Cancer | Assay | Timepoints used | Response label | Role |
|---|---|---|---|---|---|
| **GSE248249** | NSCLC (PD-(L)1 blockade) | Clariom D microarray (SST-RMA) | 13 same-patient pre + acquired-resistance pairs | all post = acquired resistance | **A7 analog** — only public same-patient lung ICI tumor RNA |
| **GSE246922** | Mouse KP / LLC1 lung | RNA-seq VST | parental / IFNγ / ICB-resistant | cell-line phenotype | TISMO stand-in (portal not downloadable) |
| **GSE207422** | NSCLC (neoadjuvant anti-PD-1 + platinum chemo) | scRNA-seq (BD Rhapsody WTA) | 3 pre-treatment biopsies + 12 post-treatment resections | pathologic MPR vs NMPR/pCR | POST axis — pre→post, epithelial compartment (0 pairs) |
| **GSE207422** | NSCLC | bulk RNA-seq (log2 TPM) | 24 pre-treatment biopsies | MPR vs NMPR (+ RECIST) | Baseline association |
| **GSE126044** | NSCLC (anti-PD-1) | bulk RNA-seq (counts) | 16 pre-treatment tumors | responder vs non-responder | Baseline verification |
| **GSE135222** | NSCLC (anti-PD-1/PD-L1) | bulk RNA-seq (TPM) | 27 pre-treatment tumors | DCB (PFS ≥ 6 mo) vs NDB | Baseline verification |
| **GSE91061** (Riaz 2017) | **Melanoma** (nivolumab) | bulk RNA-seq (FPKM) | same-patient Pre + On (43 pairs) | PRCR / SD / PD | **Orthogonal** — within-patient pre→on dynamics |

Collectively these cover the **PRE** (all NSCLC baselines), **ON** (GSE91061 early
on-treatment) and **POST** (GSE207422 post-neoadjuvant resection) axes of ICI therapy.

## Gene identifiers
- TACSTD2 (TROP2): symbol `TACSTD2`, Ensembl `ENSG00000184292`, Entrez `4070`.
- CLDN4 (Claudin-4): symbol `CLDN4`, Ensembl `ENSG00000189143`, Entrez `1364`.

## Key methodological choices
- **Epithelial compartment (scRNA).** TACSTD2/CLDN4 are epithelial genes, so per-cell
  values are meaningful only in epithelial cells. Cells are called epithelial when
  `EPCAM>0 & PTPRC==0` (raw UMI). Expression is CP10k-normalised then `log1p`, and
  summarised per sample as the mean over epithelial cells (pseudobulk) so that signals
  are not driven by changing immune infiltration. Samples with <20 epithelial cells are
  dropped from group tests.
- **Bulk normalisation.** GSE207422 already log2 TPM; GSE135222 TPM → `log2(TPM+1)`;
  GSE126044 counts → `log2(CPM+1)`; GSE91061 FPKM → `log2(FPKM+1)`.
- **Statistics.** Two-group contrasts: Mann-Whitney U with AUC = P(responder > non-responder)
  as a nonparametric effect size and mean log2 fold-change. Paired pre→on: Wilcoxon
  signed-rank on within-patient deltas. Multiple testing: Benjamini-Hochberg FDR across
  the 12 response contrasts (`results/fable_paired/tables/master_summary.csv`, `p_fdr_bh`).
- **Response definitions.** GSE207422: pathologic MPR/pCR = responder, NMPR = non-responder.
  GSE126044: deposited responder/non-responder. GSE135222: durable clinical benefit proxied
  by PFS ≥ 6 months (no RECIST in the record). GSE91061: PRCR = responder, PD = non-responder,
  SD reported separately.

## Verification (positive controls) — `05_verify.py`
V1 epithelial ≫ immune for both genes; V2 EPCAM/PTPRC compartment specificity; V3 TACSTD2~CLDN4
co-expression across bulk tumors (Spearman r=0.86); V4 sample counts match metadata; V5 pipeline
determinism. All 7/7 passed (`results/fable_paired/tables/verification_report.csv`).

## Limitations
- GSE207422 pre-vs-post is **cross-patient** (only 3 pre, 12 post), and post-treatment tumors are
  confounded by concurrent chemotherapy; pathologic response ≠ RECIST.
- GSE91061 is **melanoma**, used only to probe the *direction* of on-treatment change, not as a
  lung result.
- Epithelial calling is marker-based, not CNV-refined malignant-cell identification.
- Cohorts are small and use heterogeneous ICI regimens; findings are **hypothesis-generating**.

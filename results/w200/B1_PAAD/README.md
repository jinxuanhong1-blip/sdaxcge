# B1 analog — TCGA-PAAD: TACSTD2 (TROP2) × CLDN4 co-expression

Analog of the B1 co-expression readout from the LUAD/LUSC slice
(`results/fable_tcga`: rho = 0.53 LUAD, 0.39 LUSC), run in TCGA-PAAD.

## TL;DR

**中文**：在 TCGA-PAAD 原发肿瘤（n=178，每患者一个 -01 样本）中，TACSTD2 与
CLDN4 呈强正共表达：Spearman rho = **0.71**（bootstrap 95% CI 0.61–0.79，
p = 1.8e-28；Pearson r = 0.71）。剔除 8 例神经内分泌癌（8246/3，它们全部是
双低离群点，会拉高相关系数）后 rho = **0.68**；以 ABSOLUTE 纯度做秩偏相关
校正后 rho = **0.69**（两基因与 DNA 纯度本身几乎不相关：rho = −0.05 / +0.02）。
在 20,282 个表达基因中，CLDN4 与 TACSTD2 的共表达排名第 **13** 位
（99.94 百分位）。PAAD 的共表达强度**高于** LUAD（0.53）和 LUSC（0.39）。

**English**: In TCGA-PAAD primary tumors (n=178, one -01 sample per patient),
TACSTD2 and CLDN4 are strongly co-expressed: Spearman rho = **0.71**
(bootstrap 95% CI 0.61–0.79, p = 1.8e-28; Pearson r = 0.71). Excluding the
8 neuroendocrine carcinomas (morphology 8246/3) — all of which are double-low
outliers that stretch the correlation — leaves rho = **0.68**. A rank-based
partial correlation given ABSOLUTE tumor purity gives rho = **0.69**
(neither gene tracks DNA purity itself: rho = −0.05 / +0.02). Among 20,282
expressed genes, CLDN4 ranks **13th** (99.94th percentile) as a TACSTD2
co-expression partner, transcriptome-wide. Co-expression in PAAD is
**stronger** than in LUAD (0.53) or LUSC (0.39).

## Data (all open access; nothing controlled)

| File | Source | md5 |
|---|---|---|
| TCGA-PAAD.star_tpm.tsv.gz (log2(TPM+1), GENCODE v36) | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.star_tpm.tsv.gz) | bd179dbd6c3c4ba98dda0770c2aa266f |
| TCGA-PAAD.clinical.tsv.gz (histology) | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.clinical.tsv.gz) | 6d50fc498064228a7d25f15160325ec0 |
| gencode.v36 probemap | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap) | 59d24b459af04b543cf1d5d4161a98fc |
| ABSOLUTE purity (open PanCanAtlas supplement) | [GDC open API](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) | 8ea2ca92c8ae58350538999dfa1174da |

## Methods

- Primary tumors only (barcode sample code -01), one sample per patient →
  n = 178. Matrix is Xena GDC STAR TPM, log2(TPM+1); Ensembl IDs mapped to
  symbols via the GENCODE v36 probemap.
- Pearson and Spearman correlations; 95% CI on Spearman rho by 2,000-fold
  case-resampling bootstrap (seed 20260816).
- Transcriptome-wide ranking: Spearman rho vs TACSTD2 for every gene
  expressed at log2(TPM+1) > 1 in ≥ 20% of tumors (20,282 genes after
  collapsing duplicate symbols to the highest-mean row).
- Sensitivity 1: exclude neuroendocrine carcinomas (morphology 8246/3,
  n = 8), the well-known non-PDAC contaminant of TCGA-PAAD.
- Sensitivity 2: rank-based partial correlation adjusting for ABSOLUTE
  purity (available for 158/178 tumors).

## Results

| Analysis | n | Spearman rho | p |
|---|---|---|---|
| All primary tumors | 178 | **0.709** [0.61, 0.79] | 1.8e-28 |
| Excl. neuroendocrine (8246/3) | 170 | 0.675 | 5.5e-24 |
| Purity subset, unadjusted | 158 | 0.688 | 1.6e-23 |
| Partial rho given ABSOLUTE purity | 158 | 0.690 | 1.1e-23 |

Pearson r (all tumors) = 0.705 (p = 4.5e-28). Median expression:
TACSTD2 9.11, CLDN4 8.14 log2(TPM+1) — both highly expressed in PAAD.

Transcriptome-wide, CLDN4 is the **13th strongest** TACSTD2 partner of
20,282 expressed genes (99.94th percentile). The genes ahead of it (GJB3,
TINAGL1, GNA15, PPP1R13L, SPINT1, RHOC, ST14, PRSS22, PLEKHN1, SLC2A1, …;
rho 0.71–0.76) are themselves epithelial/junctional or squamous-program
genes — consistent with one shared malignant-epithelial program rather
than a CLDN4-specific link. See `top25_tacstd2_partners.csv`.

Files: `coexpression_main.csv`, `purity_adjusted.csv`,
`top25_tacstd2_partners.csv`, `cohort_summary.json`,
`fig1_scatter_tacstd2_cldn4.png`.

## Honest caveats

1. **The 8 neuroendocrine carcinomas are all double-low outliers**
   (visible in `fig1`); they stretch the correlation, and excluding them
   drops rho from 0.71 to 0.68. The signal does not depend on them, but
   the headline number should be quoted with this in mind.
2. **Bulk co-expression in a stroma-rich cancer partly reflects shared
   epithelial content.** ABSOLUTE (DNA-based) purity adjustment leaves rho
   essentially unchanged and neither gene correlates with DNA purity, but
   DNA purity is an imperfect proxy for the epithelial mRNA fraction, so
   compartment confounding is reduced, not eliminated. Single-cell/spatial
   data would be needed to place both genes in the same cells.
3. mRNA ≠ protein; ADC-target relevance rests on IHC, which this cannot
   replace.
4. CLDN4 ranks in the top 0.1% of TACSTD2 partners but is not uniquely
   privileged — it sits inside a broader epithelial co-expression program
   (12 genes rank higher).
5. TCGA-PAAD is a resected, mostly stage-II cohort; no ICI or ADC
   treatment. Nothing here speaks to therapy response.

## Reproduce

```bash
pip install pandas numpy scipy matplotlib
python3 scripts/w200/B1_PAAD/download_data.py   # → /tmp/b1_paad_data (override: B1_PAAD_DATA)
python3 scripts/w200/B1_PAAD/run_analysis.py    # → results/w200/B1_PAAD/
```

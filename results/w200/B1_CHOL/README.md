# B1 analog — TCGA-CHOL: TACSTD2 (TROP2) × CLDN4 surface rank + co-expression

Analog of the B1 surface-rank / co-expression readout
(`results/w200/B1_BRCA`: CLDN4 is Spearman #4 of 2618 surface genes vs TACSTD2;
`results/w200/B1_PAAD`: rho = 0.71; LUAD/LUSC slice rho = 0.53 / 0.39),
run in TCGA-CHOL.

## TL;DR

**中文**：TCGA-CHOL 原发肿瘤只有 **n=36**（每患者一个 -01 样本；另有 9 对癌旁）。
TACSTD2 与 CLDN4 呈**中等、不确定的**正共表达：Spearman rho = **0.41**
（bootstrap 95% CI **0.09–0.68**，p = 0.012；Pearson r = 0.42）。
区间几乎碰到 0，**不能**当成稳定的强共表达。在 2,394 个可排序的
surfaceome 基因中，CLDN4 只排第 **153**（93.7 百分位），FDR q = 0.16，
**不是** TACSTD2 的顶级表面共表达伙伴，也过不了多重检验。
CLDN4 本身在 CHOL 里表达很高（surface 丰度第 **49**/2,619，中位
log2 = 12.50）；TACSTD2 只是中上（第 **299**，中位 10.15）。
ABSOLUTE 纯度校正几乎不动（偏相关 0.40）；32 例肝内胆管癌子集
rho = 0.53，但那是事后切分。

**English**: TCGA-CHOL primary tumors are a **small** cohort (n=36, one -01
sample per patient; 9 matched adjacent normals). TACSTD2 and CLDN4 show
**moderate, uncertain** positive co-expression: Spearman rho = **0.41**
(bootstrap 95% CI **0.09–0.68**, p = 0.012; Pearson r = 0.42). The interval
nearly includes 0; this is **not** a stable strong co-expression claim.
Among 2,394 rankable surfaceome genes, CLDN4 is only **#153** (93.7th
percentile), FDR q = 0.16 — **not** a top TACSTD2 surface partner, and not
significant after multiple testing. CLDN4 itself is highly expressed
(surface-abundance rank **#49**/2,619, median log2 = 12.50); TACSTD2 is only
upper-quartile (#299, median 10.15). ABSOLUTE purity adjustment barely
moves the correlation (partial rho = 0.40). The intrahepatic-only subset
(n=32) is stronger (rho = 0.53) but that cut is post hoc.

**Question this slice answers:** is CLDN4 the top (or a privileged) surface
co-expression partner of TACSTD2 in cholangiocarcinoma?

**Answer: no.**

## Data (all open access; nothing controlled)

| File | Source | md5 |
|---|---|---|
| `CHOL.HiSeqV2.gz` (log2(norm_count+1), HGNC symbols) | [UCSC Xena TCGA hub](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.CHOL.sampleMap%2FHiSeqV2.gz) `TCGA.CHOL.sampleMap/HiSeqV2` | 1d939ba72b7e0d241908d854957ecf6a |
| `table_S3_surfaceome.xlsx` (in-silico surfaceome) | Bausch-Fluck et al., *PNAS* 2018; copy from [steveneschrich/surfaceome](https://github.com/steveneschrich/surfaceome) (`data-raw/surfy/table_S3_surfaceome.xlsx`) | bf95ed25ffc933fbe0c0383ab19b8652 |
| `tcga_absolute_purity.txt` | [PanCanAtlas ABSOLUTE, GDC open API](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) | 8ea2ca92c8ae58350538999dfa1174da |
| `TCGA-CHOL.clinical.tsv.gz` | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CHOL.clinical.tsv.gz) | a3a649719273fa58d9a8ead0f01d5dc5 |

The GDC clinical table lists 51 CHOL primary tumors; HiSeqV2 covers **36**.
The remaining 15 are the known DNA-only / RNA-missing extension cases and
are excluded because they have no expression.

## Methods

- Primary tumors only (barcode sample-type `01`), one sample per patient →
  n = 36. Matrix is Xena HiSeqV2, log2(norm_count+1), already on HGNC symbols
  (same family as B1_BRCA).
- Surface-gene universe: Bausch-Fluck 2018 table S3 “in silico surfaceome
  only” sheet, `UniProt gene` column. 2,886 protein rows → 2,799 unique
  symbols → 2,619 present in HiSeqV2. TACSTD2 itself is removed from the
  partner ranking. 224 surface genes with zero variance in these 36 tumors
  yield undefined correlations and are dropped → **2,394 rankable** partners.
- Primary metric: Spearman correlation. Pearson is reported as a companion.
  95% CI on Spearman rho by 2,000-fold case-resampling bootstrap (seed
  20260816). BH-FDR across the surfaceome ranking.
- Surface-abundance rank: median tumor expression among the 2,619 surface
  genes present in the matrix (anchor **kept** — this is an abundance rank,
  not a partner rank).
- Transcriptome-wide rank: Spearman vs TACSTD2 for every gene with
  log2(norm_count+1) > 1 in ≥ 20% of tumors (17,676 genes after collapsing
  duplicate symbols to the highest-mean row).
- Sensitivity 1: rank-based partial correlation given ABSOLUTE purity
  (available and called for all 36 tumors).
- Sensitivity 2: intrahepatic bile-duct origin only (n=32). The other four
  expression cases are extrahepatic (2), liver (1), gallbladder (1). All 36
  are morphology 8160/3 cholangiocarcinoma — there is no neuroendocrine
  contaminant subset analogous to TCGA-PAAD.
- Tumor vs adjacent normal: 9 matched pairs (Wilcoxon signed-rank) plus
  unpaired Mann–Whitney on 36 vs 9.
- `w200` is the reporting window (top 200 partners written to `top200.csv`).
  It does not change the full ranking.

## Results

### TACSTD2–CLDN4 co-expression

| Analysis | n | Spearman rho | 95% CI | p |
|---|---|---|---|---|
| All primary tumors | 36 | **0.412** | 0.09–0.68 | 0.012 |
| Intrahepatic only | 32 | 0.525 | 0.25–0.74 | 0.002 |
| Partial rho given ABSOLUTE purity | 36 | 0.404 | — | 0.015 |

Pearson r (all tumors) = 0.418 (p = 0.011). Neither gene tracks DNA purity
(TACSTD2 rho = −0.13, CLDN4 rho = −0.11).

Median expression in tumors: TACSTD2 10.15, CLDN4 12.50
log2(norm_count+1). Both are on the high side of the transcriptome; CLDN4
especially so.

### Surface rank (the B1 question)

| metric | TACSTD2 | CLDN4 |
|---|---|---|
| Surface-abundance rank (median expr) | **#299 / 2,619** (88.6th pct) | **#49 / 2,619** (98.2th pct) |
| Surface co-expression rank vs TACSTD2 | — (anchor) | **#153 / 2,394** (93.7th pct) |
| Surface co-expression FDR q | — | **0.16** (not significant) |
| Transcriptome-wide partner rank | — | #1,031 / 17,676 (94.2th pct) |

Actual #1 surface partner of TACSTD2 is **FZD6** (rho = 0.74, q = 4.7e-4).
Only **43 / 2,394** surface genes have BH-FDR q < 0.05. CLDN4 is not among
them. The genes ahead of CLDN4 are a mixed epithelial / adhesion set
(FZD6, MYADM, QSOX1, EPHB3, ADAM8, ENTPD3, SORT1, ITGA2, …) — not a
tight-junction-specific neighborhood that would privilege CLDN4.

### Tumor vs adjacent normal (n_normal = 9)

| gene | median T | median N | unpaired Δ | MWU p | paired Δ | Wilcoxon p |
|---|---|---|---|---|---|---|
| TACSTD2 | 10.15 | 8.28 | +1.87 | 0.067 | +3.58 | 0.027 |
| CLDN4 | 12.50 | 7.58 | +4.92 | 7.5e-5 | +4.12 | 0.0039 |

CLDN4 is clearly higher in tumor than adjacent tissue. TACSTD2 is only
marginally higher in the unpaired test; the 9 paired deltas are larger and
Wilcoxon-significant, but n=9 is small.

Files: `coexpression_main.csv`, `purity_adjusted.csv`,
`coexpression_TACSTD2_surfaceome.csv`, `top200.csv`,
`top25_tacstd2_partners.csv`, `surface_abundance_rank.csv`,
`tumor_vs_normal.csv`, `sample_level.csv`, `cohort_summary.json`,
`summary.json`, `summary.md`,
`fig1_scatter_tacstd2_cldn4.png`, `fig2_top_surface_partners.png`,
`fig3_surface_abundance.png`.

## Honest caveats

1. **n=36 is the binding constraint.** A Spearman CI of 0.09–0.68 is
   compatible with “almost nothing” and with “moderately strong.” Point
   estimates should not be quoted without the interval. Rank among 2,394
   surface genes at this n is noisy; #153 vs #50 would swap under small
   perturbations.
2. **CLDN4 does not survive FDR** among surfaceome partners (q = 0.16).
   The pairwise p = 0.012 is the unadjusted two-gene test only. Claiming
   CLDN4 as a privileged TACSTD2 surface partner in CHOL is not supported.
3. **Bulk co-expression in a stroma-variable cancer partly reflects shared
   epithelial content.** ABSOLUTE (DNA) purity adjustment leaves rho
   essentially unchanged and neither gene correlates with DNA purity, but
   DNA purity is an imperfect proxy for the epithelial mRNA fraction.
   Single-cell / spatial data would be needed to place both genes in the
   same cells.
4. **The intrahepatic-only rho = 0.53 is post hoc.** It is reported as a
   sensitivity, not a primary claim. Four non-intrahepatic cases are too
   few to interpret on their own.
5. **mRNA ≠ protein.** ADC-target relevance rests on IHC / protein, which
   this cannot replace. High CLDN4 mRNA (abundance #49) is consistent with
   a surface-target hypothesis but is not evidence of protein or of
   internalization.
6. **No ICI or TROP2-ADC treatment in TCGA-CHOL.** Nothing here speaks to
   therapy response.
7. **Cross-cohort comparison is qualitative only.** Different B1 slices
   used HiSeqV2 (BRCA, this slice) vs GDC STAR TPM (PAAD). Directionally,
   CHOL (0.41) sits near LUSC (0.39) and BRCA (0.35) and below LUAD (0.53)
   and PAAD (0.71) — all of those CIs/n differ.

## Reproduce

```bash
pip install pandas numpy scipy matplotlib openpyxl
python3 scripts/w200/B1_CHOL/download_data.py   # → /tmp/b1_chol_data (override: B1_CHOL_DATA)
python3 scripts/w200/B1_CHOL/run_analysis.py    # → results/w200/B1_CHOL/
```

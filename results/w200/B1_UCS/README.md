# B1 analog — TCGA-UCS: TACSTD2 surface-gene rank of CLDN4

Analog of claim B1 (“in TCGA, CLDN4 is the *top* cell-surface gene
co-expressed with TACSTD2 / TROP2”), run in uterine carcinosarcoma.

Outputs live in this directory. Scripts: `scripts/w200/B1_UCS/`.

## TL;DR

**中文**：TCGA-UCS 57 例原发肿瘤（每患者一个 -01 样本；全部为混合性 Müllerian
瘤 / 癌肉瘤）。在 2,671 个可计算的 surfaceome 基因中，CLDN4 与 TACSTD2 的
Spearman 共表达排名第 **3**（ρ = **0.835**，bootstrap 95% CI 0.70–0.91；
Pearson 第 1，r = 0.863）。第 1 名是 **CRB3**（ρ = 0.851），第 2 名是 EPHA1。
n=57 很小：对样本做 1,000 次重抽样，CLDN4 的秩中位数是 **8**（95% 1–28），
只有 **8.5%** 的重抽样把它排到第 1。因此 B1 原话“CLDN4 是 *唯一* 最强
surface 共表达伙伴”在 UCS **不成立**；它是 top-tier（约 top 0.1%），但不是
唯一的，并且点估计秩不稳定。HiSeqV2 稳健性：CLDN4 Spearman 第 **6**，CRB3
仍是第 1。

**English:** In TCGA-UCS primary tumors (n=57, one -01 sample per patient;
all mixed Müllerian / carcinosarcoma), CLDN4 is **not** the single top
surfaceome co-expression partner of TACSTD2. On Xena GDC STAR TPM
(log2(TPM+1), GENCODE v36) it ranks **#3 of 2,671** by Spearman
(ρ = **0.835**, bootstrap 95% CI 0.70–0.91, FDR q = 6.4e-13) and **#1**
by Pearson (r = 0.863). The actual Spearman #1 is **CRB3** (ρ = 0.851);
#2 is EPHA1. Because n=57 is small, the point-estimate rank is noisy:
1,000 case-resampling bootstraps put CLDN4 at median rank **8**
(95% 1–28) and at rank 1 in only **8.5%** of draws. The same ranking on
legacy HiSeqV2 puts CLDN4 at Spearman **#6** (CRB3 still #1). Honest
statement: CLDN4 is a genuine top-tier (≈ top 0.1%) TACSTD2 surface
partner in UCS, sitting inside a tight-junction / epithelial program —
not uniquely privileged, and not stably #1.

## What "honest ranking" means here

Every surfaceome gene present in the matrix is scored against TACSTD2.
The focus gene’s true rank is reported; we do not restrict to a
hand-picked neighbourhood. Spearman is the **pre-specified primary**
metric (robust on log-TPM and to outliers). Pearson is reported
alongside so a single favourable metric cannot be cherry-picked — and
in this cohort the two metrics disagree on #1 (Pearson would make the
B1 claim look true; Spearman does not).

`w200` is the reporting window (size of `top200.csv`). It does not
change the full-universe ranking.

## Data (all open access; nothing controlled)

| File | Source | md5 |
|---|---|---|
| TCGA-UCS.star_tpm.tsv.gz (log2(TPM+1), GENCODE v36) | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-UCS.star_tpm.tsv.gz) | a1671b47ce8eb9a6b6ad946162e0ca53 |
| TCGA-UCS.clinical.tsv.gz | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-UCS.clinical.tsv.gz) | 2655410a901045b500adcd2215a65b3a |
| gencode.v36 gene probemap | [Xena GDC hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap) | 59d24b459af04b543cf1d5d4161a98fc |
| TCGA-UCS.HiSeqV2.gz (log2(norm_count+1), gene symbols) | [Xena TCGA hub](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.UCS.sampleMap%2FHiSeqV2.gz) | 7a8334eb7f6350edac2a80eb78f9e699 |
| table_S3_surfaceome.xlsx | [steveneschrich/surfaceome mirror](https://raw.githubusercontent.com/steveneschrich/surfaceome/main/data-raw/surfy/table_S3_surfaceome.xlsx) of Bausch-Fluck et al. 2018 | bf95ed25ffc933fbe0c0383ab19b8652 |
| ABSOLUTE purity (open PanCanAtlas supplement) | [GDC open API](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) | 8ea2ca92c8ae58350538999dfa1174da |

The ETH host (`wlab.ethz.ch/surfaceome/`) now serves a single-page app
for the xlsx path; the GitHub workbook is the same table S3
(“in silico surfaceome only”, 2,886 proteins). Raw matrices are
git-ignored and re-fetched by the downloader.

## Methods

- Primary tumors only (barcode sample-type `01`), one sample per
  patient → n = 57 / 57 patients. Clinical histology: Müllerian mixed
  tumor 45, carcinosarcoma NOS 11, mesodermal mixed tumor 1.
- Expression: Xena GDC STAR TPM, log2(TPM+1). Ensembl IDs mapped via
  the GENCODE v36 probemap. Surfaceome genes matched first by Ensembl
  (2,541) then by UniProt gene symbol (301); 40 unmatched; 87
  zero-variance genes dropped (undefined correlation). Universe after
  removing TACSTD2: **2,671**.
- Spearman + Pearson; two-sided p from the t approximation; BH-FDR
  across the surface universe. 95% CI on TACSTD2–CLDN4 Spearman rho by
  2,000-fold case-resampling bootstrap (seed 20260816).
- CLDN4 rank uncertainty: 1,000 sample-resampling bootstraps of the
  full surface ranking.
- Sensitivity 1: rank-based partial Spearman given ABSOLUTE purity
  (available for 56/57 tumors).
- Sensitivity 2: the same surface ranking on legacy HiSeqV2
  (symbol-indexed, n=57, universe 2,618).
- Transcriptome-wide context (not the B1 claim): Spearman vs TACSTD2
  for every gene with log2(TPM+1) > 1 in ≥ 20% of tumors (21,774 genes
  after collapsing duplicate symbols).

## Results

| Analysis | n | CLDN4 vs TACSTD2 | note |
|---|---|---|---|
| STAR TPM, Spearman (primary) | 57 | ρ = **0.835** [0.70, 0.91] | surface rank **#3 / 2,671** |
| STAR TPM, Pearson | 57 | r = 0.863 | surface rank **#1 / 2,671** |
| Partial Spearman \| ABSOLUTE purity | 56 | ρ = 0.819 | unadjusted on same 56: 0.829 |
| HiSeqV2, Spearman | 57 | ρ = 0.823 | surface rank **#6 / 2,618** |
| Transcriptome-wide Spearman | 57 | ρ = 0.835 | rank **#20 / 21,774** |

Median expression: TACSTD2 5.91, CLDN4 6.34 log2(TPM+1) — both clearly
expressed (TPM ≈ 60 and 80).

### Top surface partners of TACSTD2 (STAR TPM, Spearman)

| rank | gene | Spearman ρ | Pearson r | note |
|---|---|---|---|---|
| 1 | **CRB3** | 0.851 | 0.860 | crumbs cell-polarity complex |
| 2 | EPHA1 | 0.835 | 0.843 | |
| 3 | **CLDN4** | 0.835 | 0.863 | focus; Pearson #1 |
| 4 | SCNN1A | 0.833 | 0.832 | |
| 5 | SLC44A4 | 0.831 | 0.832 | |
| 6 | MFSD6L | 0.829 | 0.805 | |
| 7 | CLDN7 | 0.827 | 0.825 | tight junction |
| 8 | F11R | 0.826 | 0.810 | JAM-A, tight junction |
| 9 | VTCN1 | 0.816 | 0.826 | B7-H4 |
| 10 | PRSS8 | 0.815 | 0.856 | |

The genes around CLDN4 are themselves epithelial / junctional
(CRB3, CLDN7, CLDN3, F11R, CDH1, EPCAM, PVRL4/NECTIN4).
Transcriptome-wide the pattern is the same: TMEM125, TJP3, TMPRSS4,
ST14, KRT23, RAB25, CRB3, EHF, OVOL1 sit above CLDN4. This is one
shared malignant-epithelial / tight-junction program, not a
CLDN4-unique link.

HiSeqV2 top 3: CRB3, PROM2, EPHA1 — CRB3 is #1 on both matrices.

## Honest caveats

1. **The B1 wording is false in UCS on the primary metric.** CLDN4 is
   #3 (STAR) / #6 (HiSeqV2) by Spearman. It is #1 by Pearson on STAR
   TPM. Quoting only Pearson would manufacture a “yes”. Spearman was
   pre-specified.
2. **n=57 is small and the rank is unstable.** Bootstrap median rank
   is 8, 95% interval 1–28; rank 1 in 8.5% of resamples, top-5 in 39%,
   top-20 in 91%. “#3” is a point estimate, not a settled ordering
   versus CRB3 / EPHA1 (ρ differs by < 0.02).
3. **UCS is biphasic (carcinoma + sarcoma).** TACSTD2 and CLDN4 are
   epithelial. Bulk co-expression can track the carcinoma-component
   fraction, which ABSOLUTE (DNA tumor purity) does not separate from
   sarcoma. Partial rho given ABSOLUTE purity is essentially unchanged
   (0.819 vs 0.829), and neither gene correlates positively with DNA
   purity (ρ = −0.26 / −0.23, p ≈ 0.05 / 0.09) — so this is not
   simple “more tumor DNA → both genes up”. It still does not prove
   the two proteins sit on the same cells. Single-cell or spatial
   data would be required.
4. mRNA ≠ protein. ADC-target relevance needs IHC / surface protein,
   which this cannot replace.
5. TCGA-UCS is a resected cohort with no ICI or ADC treatment.
   Nothing here speaks to therapy response.
6. Surfaceome symbols follow the 2018 table (legacy HGNC; e.g. PVRL4
   = NECTIN4).

## Outputs

| file | contents |
|---|---|
| `coexpression_TACSTD2_surfaceome.csv` | full ranked surface table (Spearman/Pearson ρ, p, BH-FDR q, ranks, mean expression, % expressed) |
| `top200.csv` | the w200 window |
| `coexpression_main.csv` | TACSTD2 × CLDN4 rho + bootstrap CI |
| `purity_adjusted.csv` | partial Spearman given ABSOLUTE purity |
| `cldn4_rank_bootstrap.csv` | 1,000 bootstrap ranks of CLDN4 |
| `top25_tacstd2_partners_transcriptome.csv` | transcriptome-wide top 25 |
| `hiseqv2_top20.csv` | HiSeqV2 surface top 20 |
| `summary.json` | machine-readable run record |
| `summary.md` | short verdict + top-15 |
| `scatter_TACSTD2_vs_CLDN4.png` | expression scatter |
| `top_partners_TACSTD2.png` | top-20 bar chart, CLDN4 highlighted |

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/w200/B1_UCS/download_data.py   # → data/ (override: B1_UCS_DATA)
python3 scripts/w200/B1_UCS/run_analysis.py    # → results/w200/B1_UCS/
```

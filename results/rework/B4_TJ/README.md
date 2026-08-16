# B4 rework · CLDN4 and a documented TJ score, R vs NR

**Honest verdict: DIRECTION_ONLY_NS.** User claim *GSE126044 NR higher TJ p=0.019* is **not reproduced** for CLDN4 or for the documented Reactome tight-junction score. Median direction is NR > R, but two-sided p is 0.090 (Reactome) / 0.115 (CLDN4), not 0.019. The same features are null as a binary DCB split in GSE135222.

| Source | Cohort | Feature | n | Claimed | This recompute |
|---|---|---|---|---|---|
| User B4 | GSE126044 | TJ score | 5 R / 11 NR | NR higher, **p=0.019** | Reactome TJ mean-z: median NR 0.106 vs R -0.269; two-sided MWU **p=0.090**; one-sided NR>R p=0.045 |
| This rework | GSE126044 | **CLDN4** log2(CPM+1) | 5 R / 11 NR | (not stated) | 3.771 vs 2.542 (NR vs R); MWU two-sided p=0.115; one-sided NR>R p=0.057; AUC=0.76; Cliff δ=+0.53 |
| This rework | GSE126044 | **Reactome TJ** (R-HSA-420029) mean-z | 5 R / 11 NR | — | 0.106 vs -0.269 (NR vs R); MWU two-sided p=0.090; one-sided NR>R p=0.045; AUC=0.78; Cliff δ=+0.56 |
| This rework | GSE126044 | KEGG hsa04530 mean-z | 5 R / 11 NR | — | 0.092 vs 0.008 (NR vs R); MWU two-sided p=0.510; one-sided NR>R p=0.255; AUC=0.62; Cliff δ=+0.24 |
| This rework | GSE135222 | CLDN4 log2(TPM+1) | 7 DCB / 20 NDB | — | 7.126 vs 7.689 (NDB vs DCB); MWU two-sided p=0.685; one-sided NDB>DCB p=0.677; AUC=0.44; Cliff δ=-0.11 |
| This rework | GSE135222 | Reactome TJ mean-z | 7 DCB / 20 NDB | — | -0.010 vs 0.041 (NDB vs DCB); MWU two-sided p=0.646; one-sided NDB>DCB p=0.323; AUC=0.56; Cliff δ=+0.13 |

`p=0.019` recovered in any pre-specified GSE126044 all-sample CLDN4/TJ test: **no**.

CD8A (positive control, expect R > NR / DCB > NDB; table uses NR/NDB as the “high” arm so AUC<0.5 is the immune-hot direction):

- GSE126044 CD8A: 3.259 vs 5.807 (NR vs R); MWU two-sided p=0.001; one-sided NR>R p=1.000; AUC=0.02; Cliff δ=-0.96
- GSE135222 CD8A: 2.441 vs 3.676 (NDB vs DCB); MWU two-sided p=0.072; one-sided NDB>DCB p=0.969; AUC=0.26; Cliff δ=-0.47

---

## 中文一句话

**GSE126044 上复算不到 NR 更高 TJ、p=0.019。** CLDN4 与 Reactome 紧密连接评分都是不显著；GSE135222（DCB vs NDB）同样是空结果。CD8A 在 GSE126044 能分开 R/NR，所以不是队列完全没有信号。

---

## What was tested (pre-specified)

Two public NSCLC ICI RNA-seq series. Two primary features. All rows are reported.

1. **CLDN4** — single gene, not a signature.
2. **Documented TJ score** — Reactome **R-HSA-420029 Tight junction interactions** (MSigDB `REACTOME_TIGHT_JUNCTION_INTERACTIONS`, 30 genes). Score = mean of per-gene z-scores across the cohort. This is a published pathway list, not a custom claudin mash-up.

Secondary, still reported:

- KEGG **hsa04530 Tight junction** (standard named TJ pathway; actin/myosin-heavy).
- Reactome TJ **without CLDN4**.
- CD8A and CD8A+CD8B (sanity check that the labels work).
- TACSTD2.
- Mean expression (no z) for every set.
- GSE126044 **fresh-only** (all 5 responders are fresh; all 5 FFPE samples are NR).
- GSE135222 Spearman vs PFS days.

No gene was added or dropped to chase p=0.019. No one-sided p-value is treated as the primary claim match.

## Data

| Item | Value |
|---|---|
| GSE126044 paper | Cho JW et al. *Exp Mol Med* 2020;52:155–165. DOI 10.1038/s12276-020-0384-2 |
| GSE126044 counts | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz |
| GSE126044 labels | GEO series matrix `patient response` / `sample` (fresh vs FFPE) |
| GSE126044 transform | log2(CPM+1) from deposited raw counts |
| GSE126044 n | 5 responder / 11 non-responder (author RECIST: PR or SD>6 mo = R) |
| GSE126044 FFPE | 5 samples, **all NR**; all 5 R are fresh |
| GSE135222 paper | Jung H et al. *Nat Commun* 2019;10:2244. DOI 10.1038/s41467-019-10126-y |
| GSE135222 matrix | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz (RSEM TPM, hg19) |
| GSE135222 labels | GEO `pfs.time` + `pfs` event. **DCB = PFS ≥ 183 days**. Not author-tagged DCB in the series matrix. All censored patients have follow-up ≥ 205 d, so the 6-month cut is unambiguous. |
| GSE135222 n | 7 DCB / 20 NDB |
| GSE135222 transform | log2(TPM+1) |
| Reactome TJ | R-HSA-420029 · https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_TIGHT_JUNCTION_INTERACTIONS · https://reactome.org/content/detail/R-HSA-420029 |
| KEGG TJ | hsa04530 · https://www.kegg.jp/pathway/hsa04530 |
| Reactome genes found | GSE126044 29/30 (missing PATJ); GSE135222 30/30 |
| KEGG genes found | GSE126044 166/171; GSE135222 171/171 |

Sample IDs are matched **by name**, not by column order. The GSE126044 series-matrix order and the count-table order differ (Dis_06 / Dis_07 / Dis_10).

## Statistics

- Primary test: two-sided Mann–Whitney U.
- Claim-direction test: one-sided MWU, NR (or NDB) > R (or DCB).
- Effect: median difference, Cliff’s δ, AUC for “higher value in the NR/NDB arm”, Welch t, Cohen’s d.
- Score: mean of per-gene z-scores (population SD, ddof=0) computed **inside each cohort**. Single genes use the log-expression value itself (`mean_expr`).
- No multiple-testing correction is applied to hide or promote a hit. The user claim is a single p-value; the table above is the honest comparison.

## Results

### GSE126044 (claim cohort)

| Feature | Scoring | Median NR | Median R | Δ med | MWU p two-sided | MWU p NR>R | AUC (NR high) | Cliff δ | Fresh-only p two-sided |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | log2CPM | 3.771 | 2.542 | +1.229 | 0.115 | 0.057 | 0.76 | +0.53 | 0.329 |
| Reactome TJ | mean-z | 0.106 | -0.269 | +0.375 | 0.090 | 0.045 | 0.78 | +0.56 | 0.177 |
| Reactome TJ no CLDN4 | mean-z | 0.093 | -0.279 | +0.372 | 0.090 | 0.045 | 0.78 | +0.56 | — |
| KEGG hsa04530 | mean-z | 0.092 | 0.008 | +0.084 | 0.510 | 0.255 | 0.62 | +0.24 | — |
| CD8A | log2CPM | 3.259 | 5.807 | -2.548 | 0.001 | 1.000 | 0.02 | -0.96 | — |

### GSE135222 (second ICI NSCLC series)

Jung et al. do not deposit a binary R/NR column. DCB vs NDB at PFS ≥ 6 months is the documented clinical-benefit split used for this series.

| Feature | Scoring | Median NDB | Median DCB | Δ med | MWU p two-sided | MWU p NDB>DCB | AUC (NDB high) | Cliff δ |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | log2TPM | 7.126 | 7.689 | -0.563 | 0.685 | 0.677 | 0.44 | -0.11 |
| Reactome TJ | mean-z | -0.010 | 0.041 | -0.050 | 0.646 | 0.323 | 0.56 | +0.13 |
| KEGG hsa04530 | mean-z | 0.122 | 0.050 | +0.072 | 0.314 | 0.157 | 0.64 | +0.27 |
| CD8A | log2TPM | 2.441 | 3.676 | -1.235 | 0.072 | 0.969 | 0.26 | -0.47 |

Secondary (not the B4 claim): Spearman vs continuous PFS days. Higher Reactome/KEGG TJ tracks **shorter** PFS at the margin; CLDN4 does not. This is a trend, not a confirmation of p=0.019, and it does **not** appear in the binary DCB split.

| Feature | Spearman ρ vs PFS | p |
|---|---:|---:|
| CLDN4 | -0.143 | 0.476 |
| Reactome TJ mean-z | -0.382 | 0.049 |
| KEGG hsa04530 mean-z | -0.367 | 0.060 |
| CD8A | +0.379 | 0.051 |

## Honest caveats

- **n=16 (5 vs 11) is tiny.** Only a large effect would reliably give p=0.019. Failure to hit 0.019 is not proof that TJ is irrelevant; it is proof that **this public matrix does not support that p-value** under a documented TJ definition.
- **FFPE is confounded with NR.** Restricting to fresh tissue does not create a p=0.019 hit (CLDN4 fresh p=0.329; Reactome TJ fresh p=0.177).
- **The original TJ gene list was not provided.** A different unpublished list could in principle give p=0.019. That would be a custom score, not a documented TJ score. We did not search lists to recover 0.019.
- **KEGG hsa04530 is a poor “barrier” score** (myosins, ARP2/3, tubulin, kinases). It is included because it is the named KEGG Tight junction pathway.
- **GSE135222 DCB is derived** from PFS ≥ 183 days. 7 vs 20 is also underpowered.
- CD8A **does** separate GSE126044 R vs NR, so the response labels and the matrix are not globally inert.

## Files

- `summary.json` — machine-readable claim, methods, all tests
- `verdict.txt` — one-screen verdict
- `tables/stats.tsv` — every test
- `tables/GSE126044_sample_scores.tsv` / `GSE135222_sample_scores.tsv`
- `tables/gene_coverage.tsv`
- `genesets/` — frozen Reactome + KEGG lists
- `figures/` — boxplots
- `scripts/rework/B4_TJ/` — download + analyze

## Rerun

```bash
python3 scripts/rework/B4_TJ/download.py
python3 scripts/rework/B4_TJ/analyze.py
```

Requires: numpy, pandas, scipy, matplotlib (see `scripts/rework/B4_TJ/requirements.txt`).

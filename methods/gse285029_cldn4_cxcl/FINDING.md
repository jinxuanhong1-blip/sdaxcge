# FINDING — GSE285029 CLDN4 vs CXCL9 / CXCL10 / CXCL13 and GEP-like

**Additive only. CLDN4 only. Chemokine extra.** Public pre-ICI NSCLC WTS, **n = 234** (Koh et al., *JITC* 2025; GEO [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029)). This folder does **not** audit or retract any slide. TACSTD2 is **not** an anchor here.

IFN-compact Spearman ρ vs CLDN4 is **already known** from the sibling GSE285029 score folder (`results/w200/A11_GSE285029`; ρ = +0.204, p = 0.002, n = 234) and is restated only as context. The extra cut here is the single-gene chemokines **CXCL9 / CXCL10 / CXCL13** plus **GEP-like** (Ayers 2017 GEP18, unweighted z-mean; not NanoString TIS weights). CXCL9 and CXCL10 sit inside IFN-compact and GEP18; **CXCL13 does not**.

GEO has no RECIST, histology, PD-L1 IHC, TMB, or purity. This is bulk tumour WTS, not a cell-intrinsic call.

Numbers below are written from `tables/one_row.tsv` and `tables/spearman.tsv`.

---

## 一句话 / TL;DR

| Contrast | n | CLDN4–CXCL9 ρ (p, BH q) | CLDN4–CXCL10 ρ (p, BH q) | CLDN4–CXCL13 ρ (p, BH q) | CLDN4–GEP18 ρ (p, BH q) |
|---|---:|---|---|---|---|
| GSE285029 CLDN4 vs chemokine extra | 234 | +0.100 (0.127, q=0.127) | +0.112 (0.088, q=0.127) | +0.102 (0.118, q=0.127) | +0.190 (0.004, q=0.014) |

**Verdict:** CXCL9 / CXCL10 / CXCL13 are **NULL** on continuous Spearman (all |ρ| ≈ 0.10, all BH q = 0.127). GEP-like is the already-known **WEAK_POSITIVE** (ρ = +0.190, q = 0.014). IFN-compact (already known, not in BH): ρ = +0.204 (p = 0.002, n = 234). The IFN program association is **not** a strong single-gene CXCL9/10/13 correlation. Q4 vs Q1 is a secondary cut and is slightly more positive for CXCL9 (MWU p = 0.042; Q4∩Q4 Fisher p = 0.009) — do not upgrade the primary Spearman NULL on that tail.

---

## Honest n

| item | n | rule |
|---|---:|---|
| GEO series | 1 | [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029) |
| author ICI–RNA-seq WTS columns | **234** | Case1–Case234 on `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| genes after all-NA drop | 20069 | deposited symbols |
| CLDN4 finite | 234 | complete-case Spearman n |
| CXCL9 / CXCL10 / CXCL13 finite | 234 / 234 / 234 | all present |
| GEP18 genes present | 18/18 | Ayers 2017 symbols |
| IFN-compact genes present | 8/8 | context only |
| RECIST / histology / PD-L1 IHC / purity | **0** | not on GEO; skipped |

Primary tests use **n = 234**. Do not write a larger n. The filename says 235 columns (gene id + 234 samples). Author paper n is 234.

---

## Locked design (before ρ)

| Item | Choice |
|---|---|
| Matrix | Public author WTS only. No FASTQ / SRA. |
| File | `GSE285029_WTS_expr_count_235_032820.txt.gz` |
| Samples | **n = 234** (Case1–Case234). Author ICI–RNA-seq cohort. |
| Transform | `log2(pmax(x,0)+1)` (same clip as the sibling GSE285029 score / GSEA folders) |
| Anchor | **CLDN4 only** |
| Extra primary | CXCL9, CXCL10, CXCL13, GEP-like (Ayers GEP18 unweighted z-mean) |
| Already known | IFN-compact (`IFNG STAT1 IRF1 CXCL9 CXCL10 CXCL11 IDO1 GBP1`) — restated, not in BH |
| Sensitivity | CXCL11; CHEMO3 = z-mean of CXCL9/CXCL10/CXCL13; CD8A; CD274; raw-author Spearman; epithelial partial |
| Score | unweighted mean of per-gene z-scores on the log2(clip0+1) matrix |
| Statistic | Spearman ρ, two-sided, Fisher-z 95% CI, n = complete pairs |
| FDR | BH inside the **four extra primary** partners only |
| Epithelial z-mean | `EPCAM KRT8 KRT18 KRT19` — purity-like **sensitivity** covariate, not ABSOLUTE / ESTIMATE |
| Out of scope | TACSTD2 ranks; response labels; FASTQ; claiming a new IFN ρ |

Series: [GSE285029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285029). Title on GEO: pre-ICI NSCLC tumour WTS (PD-1 or PD-L1 blockade). Paper: Koh et al., *J Immunother Cancer* 2025; PMID [40050048](https://pubmed.ncbi.nlm.nih.gov/40050048/).

The filename says “count”. Values are continuous with negatives; they are the author-processed matrix, not raw integer counts and not TPM. Rank correlations on the clipped log matrix are the claim. Sensitivity Spearman on the raw author values (negatives kept) is reported in `tables/spearman.tsv`.

---

## Extra primary — CLDN4 vs chemokines / GEP-like

| feature | role | ρ | 95% CI | p | BH q | n | label | epi-partial ρ | epi-partial p |
|---|---|---:|---|---:|---:|---:|---|---:|---:|
| CXCL9 | extra | +0.100 | -0.029 to +0.225 | 0.127 | 0.127 | 234 | NULL | +0.044 | 0.506 |
| CXCL10 | extra | +0.112 | -0.017 to +0.237 | 0.088 | 0.127 | 234 | NULL | +0.030 | 0.650 |
| CXCL13 | extra | +0.102 | -0.026 to +0.228 | 0.118 | 0.127 | 234 | NULL | +0.089 | 0.176 |
| GEP18 | GEP-like | +0.190 | +0.063 to +0.310 | 0.004 | 0.014 | 234 | WEAK_POSITIVE | +0.153 | 0.020 |
| IFN_compact | already known | +0.204 | +0.077 to +0.323 | 0.002 | — | 234 | ASSOCIATED_POSITIVE | +0.102 | 0.120 |
| CXCL11 | sensitivity | +0.142 | +0.014 to +0.266 | 0.030 | — | 234 | WEAK_POSITIVE | +0.012 | 0.853 |
| CHEMO3 | sensitivity | +0.136 | +0.008 to +0.260 | 0.038 | — | 234 | WEAK_POSITIVE | +0.083 | 0.207 |
| CD8A | sensitivity | +0.067 | -0.061 to +0.194 | 0.305 | — | 234 | NULL | +0.093 | 0.159 |
| CD274 | sensitivity | +0.167 | +0.040 to +0.289 | 0.011 | — | 234 | WEAK_POSITIVE | +0.101 | 0.124 |

GEP18 on this same matrix was already in the sibling score folder (ρ = +0.190, p = 0.004). It is restated here because the request is CXCL + GEP-like. The new single-gene extra is CXCL9 / CXCL10 / CXCL13; CXCL13 is not a member of IFN-compact or GEP18.

CHEMO3 is the unweighted z-mean of CXCL9, CXCL10, and CXCL13 (3/3 present). It is sensitivity, not a published GEP.

---

## CLDN4-high vs CLDN4-low (Q4 vs Q1)

Ties stay in the tail. n Q4 / Q1 = 59 / 59.

| feature | median Q4 | median Q1 | Δmedian | MWU p | n Q4 | n Q1 |
|---|---:|---:|---:|---:|---:|---:|
| CXCL9 | +3.018 | +1.923 | +1.095 | 0.042 | 59 | 59 |
| CXCL10 | +4.322 | +1.557 | +2.765 | 0.049 | 59 | 59 |
| CXCL13 | +2.072 | +1.081 | +0.991 | 0.087 | 59 | 59 |
| GEP18 | +0.236 | -0.182 | +0.418 | 0.002 | 59 | 59 |
| IFN_compact | +0.256 | -0.539 | +0.796 | 8.63e-04 | 59 | 59 |
| CHEMO3 | +0.208 | -0.253 | +0.461 | 0.023 | 59 | 59 |

Independence expectation for two Q4 calls is 25%.

| feature | both Q4 | anchor Q4 | frac in feature Q4 | expected if independent | OR | Fisher p |
|---|---:|---:|---:|---:|---:|---:|
| CXCL9 | 23 | 59 | 0.390 | 0.250 | 2.47 | 0.009 |
| CXCL10 | 19 | 59 | 0.322 | 0.250 | 1.60 | 0.167 |
| CXCL13 | 20 | 59 | 0.339 | 0.250 | 1.79 | 0.085 |
| GEP18 | 21 | 59 | 0.356 | 0.250 | 1.99 | 0.039 |
| IFN_compact | 25 | 59 | 0.424 | 0.250 | 3.05 | 8.54e-04 |

---

## Extra figures

- `methods/gse285029_cldn4_cxcl/figures/fig1_spearman_forest.png`
- `methods/gse285029_cldn4_cxcl/figures/fig2_scatter.png`

Tables: `methods/gse285029_cldn4_cxcl/tables/one_row.tsv` (the headline table), `spearman.tsv`, `q4_vs_q1.tsv`, `q4_overlap.tsv`, `inventory.tsv`, `gene_coverage.tsv`, `sample_scores.tsv`.

---

## What this is not

- Not a TACSTD2 analysis. CLDN4 only.
- Not a new IFN-compact ρ. That number is already known and is restated as context.
- Not a re-run of the sibling GSE285029 score / GSEA folders. Those stay as written.
- Not FASTQ / salmon / DESeq2 from SRA (author matrix is used as deposited).
- Not a cell-intrinsic chemokine call. This is bulk pre-ICI NSCLC WTS.
- Not a response / PD-L1 IHC / purity analysis. GEO does not release those labels.
- GEP-like is unweighted Ayers GEP18 z-mean, not NanoString TIS weights.

---

## 中文摘要

只补公开 **GSE285029**（n=234）上 **CLDN4** 对 **CXCL9 / CXCL10 / CXCL13** 和 **GEP-like（Ayers GEP18）** 的 Spearman，不审不撤已有页。IFN-compact ρ 已知（ρ = +0.204, p = 0.002），本文件夹只做 chemokine extra。只做 CLDN4，不做 TACSTD2。FDR 只在四个 extra primary 内做 BH。

- 234 例 ICI 前 NSCLC 肿瘤 WTS。完整配对 n = 234。
- CLDN4–CXCL9 ρ = +0.100 (p = 0.127, q = 0.127)。
- CLDN4–CXCL10 ρ = +0.112 (p = 0.088, q = 0.127)。
- CLDN4–CXCL13 ρ = +0.102 (p = 0.118, q = 0.127)。
- CLDN4–GEP18 ρ = +0.190 (p = 0.004, q = 0.014)。
- 单基因 CXCL9/10/13 连续 Spearman 为 NULL；IFN 程序相关不是强单基因 chemokine 相关。

---

## Files

- `download.py` — GEO author matrix
- `analyze.py` — complete-case n, Spearman, Q4, extra scatters, this FINDING.md
- `tables/one_row.tsv`, `spearman.tsv`, `inventory.tsv`, `summary.json`

```bash
python3 methods/gse285029_cldn4_cxcl/download.py
python3 methods/gse285029_cldn4_cxcl/analyze.py
```

# FINDING — GSE283829 leftover extra: CLDN4 vs CD274/HLA and Q4 vs ImmuneScore

**Additive only.** Public GEO [GSE283829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283829) tumor RNA-seq counts (Lindberg et al., *J Thorac Oncol* 2025, PMID 39743139). Patient is the unit. **n = 27.**

**Already known (not re-claimed):** leftover OPEN lung ICI bulk, CLDN4 vs ESTIMATE ImmuneScore (Yoshihara immune-list **mean-z** on log2(CPM+1)) Spearman **ρ = −0.57, p = 0.0020, n = 27**. This folder reproduces that number exactly and does not treat it as a new row.

**This page adds extra cuts only:**

1. CLDN4 vs **CD274** and **HLA** (HLA-A / HLA-B / HLA-C / HLA-DRA, plus HLA-A/B/C mean-z).
2. CLDN4 **Q4 vs Q1** vs ImmuneScore (and the same extra genes).

Continuous ImmuneScore is the powered test. Q4 vs Q1 is **7 vs 7**, not 27.

Reproduce: `python3 methods/gse283829_cldn4_extra/analyze.py` (GEO files cached under `/tmp/gse283829_cldn4_extra/`).

---

## Verdict

| Extra cut | n | Effect | p | Read |
|---|---:|---:|---:|---|
| CLDN4 vs CD274 (Spearman) | **27** | ρ = **−0.37** | **0.054** | same sign as immune-low; not <0.05 |
| CLDN4 vs HLA-A / B / C / DRA | **27** | ρ = −0.16 / −0.10 / −0.32 / −0.19 | 0.42 / 0.62 / 0.10 / 0.35 | no HLA gene is a hit |
| CLDN4 vs MHC-I mean-z (A/B/C) | **27** | ρ = **−0.24** | **0.22** | matches leftover MHC mean-z (ρ=−0.21, p=0.29) |
| same six, partial \| ImmuneScore | **27** | all \|ρ\| ≤ 0.33 | all p ≥ 0.098 | infiltrate accounts for the raw CD274 tilt |
| **CLDN4 Q4 vs Q1 ImmuneScore** | **7 vs 7** | δ_median = **−0.82**; r = **−0.67** | **0.038** | same sign as the known continuous test; **do not write n=27** |
| Q4 vs Q1 CD274 / HLA / MHC-I | **7 vs 7** | all r ∈ [−0.47, −0.14] | all p ≥ 0.16 | extra Q4 cuts on CD274/HLA are null |

**What holds:** the already-known continuous ImmuneScore anti-correlation. The Q4 vs Q1 ImmuneScore cut is the same direction at honest **7 vs 7** (p=0.038). CD274 is the closest extra gene (ρ=−0.37, p=0.054) and goes to null after ImmuneScore.

**What does not hold:** a CLDN4-high PD-L1 / HLA program in this leftover. No extra HLA gene, and no Q4 vs Q1 CD274/HLA contrast, is significant. Do not upgrade p=0.054 CD274 to a claim.

---

## Honest n

| Item | n | Note |
|---|---|---|
| GEO RNA-seq samples (matrix n) | **27** | do not write a larger n |
| CR / SD / PD | **7 / 10 / 10** | GEO `disease stage` (not TNM); no PR label |
| Continuous extra tests (CLDN4 vs CD274/HLA) | **27** | pairwise-complete; all genes present |
| ESTIMATE ImmuneScore genes used | **141 / 141** | Yoshihara immune list; same definition as leftover bulk |
| MHC-I mean-z genes | **3 / 3** | HLA-A, HLA-B, HLA-C |
| CLDN4 Q1 / Q2 / Q3 / Q4 | **7 / 7 / 6 / 7** | `rank(method='first')` then `qcut` |
| **Q4 vs Q1 used** | **7 vs 7** | middle **13** samples excluded |

Do not write n=27 on a Q4 vs Q1 row. Do not write n=7 vs 10 (CR vs PD) on these extra immune-axis rows.

---

## Extra continuous: CLDN4 vs CD274 / HLA (n=27)

Expression is log2(CPM+1) from the public supplementary raw-count matrix. MHC-I = mean of per-gene z-scores for HLA-A/B/C (not the leftover 6-gene MHC that also includes B2M/TAP1/TAP2). Partial Spearman residualizes ranks on ImmuneScore.

| Endpoint | n | ρ | p | partial ρ \| ImmuneScore | partial p |
|---|---:|---:|---:|---:|---:|
| CD274 | 27 | **−0.375** | **0.054** | −0.233 | 0.253 |
| HLA-A | 27 | −0.161 | 0.424 | +0.054 | 0.793 |
| HLA-B | 27 | −0.101 | 0.617 | +0.124 | 0.547 |
| HLA-C | 27 | −0.322 | 0.102 | −0.252 | 0.215 |
| HLA-DRA | 27 | −0.186 | 0.352 | +0.331 | 0.098 |
| MHC-I (A/B/C mean-z) | 27 | −0.242 | 0.223 | −0.053 | 0.798 |
| ImmuneScore *(already known)* | 27 | **−0.567** | **0.0020** | — | — |

CD274 is the only extra gene near 0.05. After ImmuneScore it is not. HLA-DRA’s partial p=0.098 is a sign flip (raw negative, residual positive) and is not a claim.

---

## Extra cut: CLDN4 Q4 vs Q1

Quartiles on CLDN4 log2(CPM+1) among all 27 tumors. Q1 = lowest 7, Q4 = highest 7. Two-sided Mann–Whitney U. **r** = rank-biserial \(2U/(n_4 n_1)-1\) (positive = higher in Q4). 95% CI on r = percentile bootstrap, 2000 resamples, seed 20260817. δ = median(Q4) − median(Q1).

### ImmuneScore (the requested extra cut)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | δ | U | r (95% CI) | p |
|---|---:|---:|---:|---:|---:|---|---:|
| ImmuneScore | **7 / 7** | −0.398 | +0.421 | **−0.819** | 8 | **−0.67 (−1.00 to −0.14)** | **0.038** |

Same sign as the known continuous ρ=−0.57 (n=27). The quartile cut throws away Q2+Q3 (13 tumors). That is why p moves from 0.0020 to 0.038. Report **7 vs 7**.

### CD274 / HLA (same Q4 vs Q1 gate)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | δ | r | p |
|---|---:|---:|---:|---:|---:|---:|
| CD274 | 7 / 7 | 3.47 | 4.76 | −1.29 | −0.43 | 0.209 |
| HLA-A | 7 / 7 | 9.46 | 10.04 | −0.58 | −0.27 | 0.456 |
| HLA-B | 7 / 7 | 10.41 | 10.70 | −0.29 | −0.22 | 0.535 |
| HLA-C | 7 / 7 | 9.47 | 10.12 | −0.65 | −0.47 | 0.165 |
| HLA-DRA | 7 / 7 | 8.64 | 8.71 | −0.07 | −0.14 | 0.710 |
| MHC-I mean-z | 7 / 7 | 0.028 | 0.867 | −0.84 | −0.35 | 0.318 |

All extra Q4 vs Q1 CD274/HLA intervals include chance. Do not use the ImmuneScore Q4 p as if it applied to PD-L1 or HLA.

---

## Methods (fixed before looking at p)

- Source: GEO supplementary `GSE283829_raw_express_matrix_all_samples.txt.gz` + series matrix. Series matrix has **no** expression table.
- Genes: CLDN4 = ENSG00000189143; CD274 = ENSG00000120217; HLA-A/B/C/DRA from the leftover HGNC map. Expression = log2(CPM+1).
- ImmuneScore = Yoshihara 2013 immune-signature **mean z** on the author count matrix (same leftover definition; **not** the Affymetrix-calibrated ESTIMATE purity transform). 141/141 immune genes mapped.
- Quartiles: `rank(method='first')` then `qcut` into 4 equal-count bins. Q1 = lowest 25%, Q4 = highest 25%.
- Continuous: two-sided Spearman. Partial Spearman = rank residual after ImmuneScore, t-test on n−3 df.
- Q4 vs Q1: two-sided Mann–Whitney U. No cutpoint search, no in-sample ROC, no multivariable classifier.
- Continuous ImmuneScore is cited from leftover bulk (`methods/lung_ici_bulk_cldn4_leftover`) and reproduced here. It is not a new claim.

---

## Honest caveats

- **n = 27** for continuous extra genes; **7 vs 7** for Q4 vs Q1. Only a large quartile effect is detectable. A null on CD274/HLA does not prove no association; it also does not support one.
- Bulk diagnostic-biopsy RNA. CLDN4 is epithelial; ImmuneScore / CD274 / HLA move with infiltrate and purity. The partials say the raw CD274 tilt is mostly that.
- No PD-L1 IHC percent, no HLA typing, no purity from ESTIMATE’s calibrated transform. Single-arm ICI series: prognostic vs predictive cannot be separated.
- Leftover CR vs PD for CLDN4 remains null (7 vs 10, p=0.96). These extra cuts are immune-axis rows, not a response recut.
- No meta-analysis with other leftover ICI series on this page.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `estimate_gene_sets.csv` — Yoshihara stromal/immune lists (same as leftover bulk)
- `hgnc_symbol_ensembl.tsv` — symbol ↔ Ensembl (CD274 added)
- `tables/extra_cuts.tsv` — all extra tests
- `tables/honest_n.tsv` — n audit
- `tables/summary.json`
- `processed/sample_scores.tsv`
- `figures/CLDN4_vs_CD274.png`, `CLDN4_vs_HLA_A.png`, `CLDN4_vs_MHC_I.png`, `CLDN4_Q4_vs_Q1_ImmuneScore.png`

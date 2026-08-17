# Finding — TCGA-LUAD leftover: CLDN4 vs CXCL9 / CXCL10 / CXCL13 / CXCR3

**Additive single-gene extra.** Predictor is **CLDN4 RNA only**. No TACSTD2 gate.

This folder does **not** re-audit Q4 vs CD8 / ImmuneScore / CD274 (PR #352) and does **not** recompute TJ7 vs CD8 (PR #69 / #90). Those pages stay as written. The leftover axis here is **CXCL9 / CXCL10 / CXCL13 / CXCR3**: continuous Spearman, CLDN4 Q4 vs Q1, and ESTIMATE residual.

Primary matrix: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2`, log2(RSEM norm_count+1), primary tumors (`-01`). ESTIMATE = official MD Anderson RNAseqV2 (`Immune_score`, `Stromal_score`, `ESTIMATE_score`). **Honest n = 515** (HiSeqV2 ∩ ESTIMATE). Do not write n = 502 (ABSOLUTE drop, not used). Do not write n = 515 for the MWU tests.

TCGA-LUAD is treatment-naive surgical RNA. **No ICI response labels.** Do not read CXCL / CXCR3 as immunotherapy outcomes.

Numbers below are written from `tables/one_row.tsv`, `tables/spearman.tsv`, `tables/q4_vs_q1.tsv`, and `tables/estimate_residual.tsv`.

---

## One-row table

| dataset | n | Q4 vs Q1 | CLDN4–CXCL9 ρ (p, BH q) | CLDN4–CXCL10 ρ (p, BH q) | CLDN4–CXCL13 ρ (p, BH q) | CLDN4–CXCR3 ρ (p, BH q) | CXCL9 ESTIMATE residual ρ (p, q) | CXCL10 residual ρ (p, q) | CXCL13 residual ρ (p, q) | CXCR3 residual ρ (p, q) |
|---|---:|---|---|---|---|---|---|---|---|---|
| TCGA-LUAD HiSeqV2 ∩ MDACC ESTIMATE | **515** | **129 vs 129** | −0.111 (0.011, q=0.023) | −0.059 (0.181, q=0.228) | −0.129 (0.0034, q=0.013) | +0.053 (0.228, q=0.228) | −0.080 (0.071, q=0.095) | −0.015 (0.727, q=0.727) | −0.102 (0.021, q=0.042) | +0.133 (0.0026, q=0.010) |

**Verdict:** Unadjusted, CLDN4-high is **weakly CXCL9-low and CXCL13-low**. CXCL10 and CXCR3 are **NULL**. After ESTIMATE residual, CXCL13 leftover stays weakly negative; CXCL9 residual is NS; CXCL10 stays null; CXCR3 leftover flips to a small **positive** (suppressor: CXCR3 tracks ImmuneScore ρ = +0.72 while CLDN4 is nearly orthogonal). Do not quote the CXCR3 residual without the unadjusted null. Do not upgrade this to “CLDN4-high LUAD is chemokine-cold.”

---

## Honest n

| | n |
|---|---:|
| HiSeqV2 primary tumors (`-01`) | **515** |
| MDACC ESTIMATE primaries (`-01`) | **515** |
| HiSeqV2 ∩ ESTIMATE (analysis set) | **515** |
| CLDN4 / CXCL9 / CXCL10 / CXCL13 / CXCR3 / ESTIMATE non-NA | **515 / 515 / 515 / 515 / 515 / 515** |
| Q1 / Q2 / Q3 / Q4 | **129 / 129 / 128 / 129** |
| Q4 vs Q1 used (all four CXCL endpoints) | **129 vs 129** |
| ABSOLUTE-restricted set used elsewhere | 502 (13 dropped; **not used**) |
| TJ7 genes in this test | **0** |
| CD8 / CD274 re-audit | **0** (PR #352) |
| TACSTD2 filter | **none** |

Quartiles are assigned on all 515 CLDN4 values. Do not write n=515 for the MWU tests — those are **129 vs 129**. Do not write n=502 for this slice.

---

## Continuous Spearman (primary, n=515)

Two-sided Spearman. BH inside the four leftover partners only. Fisher-z 95% CI.

| Endpoint | n | ρ | 95% CI | p | BH q | label |
|---|---:|---:|---|---:|---:|---|
| CXCL9 | 515 | −0.111 | −0.196 to −0.025 | **0.011** | **0.023** | WEAK_NEGATIVE |
| CXCL10 | 515 | −0.059 | −0.145 to +0.028 | 0.181 | 0.228 | NULL |
| CXCL13 | 515 | −0.129 | −0.213 to −0.043 | **0.0034** | **0.013** | WEAK_NEGATIVE |
| CXCR3 | 515 | +0.053 | −0.033 to +0.139 | 0.228 | 0.228 | NULL |

---

## Q4 vs Q1 (same n, supporting)

CLDN4 cuts among the 515: Q1 ≤ 12.52, Q4 ≥ 13.66 (HiSeqV2 log2(norm_count+1)). Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | r | p |
|---|---:|---:|---:|---:|---:|---:|
| CXCL9 | 129 / 129 | 10.06 | 10.52 | −0.458 | −0.18 | **0.012** |
| CXCL10 | 129 / 129 | 8.52 | 8.90 | −0.374 | −0.11 | 0.121 |
| CXCL13 | 129 / 129 | 8.45 | 9.15 | −0.698 | −0.18 | **0.015** |
| CXCR3 | 129 / 129 | 6.80 | 6.54 | +0.265 | +0.03 | 0.713 |

Q4 vs Q1 and Spearman agree: CXCL9 / CXCL13 yes (small), CXCL10 / CXCR3 no.

---

## ESTIMATE residual

Partial Spearman of CLDN4 vs each chemokine **| ESTIMATE_score** (algebraic first-order partial). Matching OLS residual of the chemokine on ESTIMATE_score is in `tables/estimate_residual.tsv` and agrees in sign and magnitude.

These four genes **track ImmuneScore** (ρ = +0.65 / +0.64 / +0.60 / +0.72 for CXCL9 / CXCL10 / CXCL13 / CXCR3). CLDN4 itself is nearly orthogonal to ImmuneScore (ρ = −0.065, p = 0.14) and only weakly related to ESTIMATE_score (ρ = −0.082, p = 0.064). Residualising on ESTIMATE_score is therefore a **within-infiltrate leftover**, not an independent purity control. ImmuneScore is a term in ESTIMATE_score; residualising a chemokine on ImmuneScore is the more circular of the two and is reported only as a sensitivity column.

| Endpoint | n | partial ρ \| ESTIMATE_score | p | BH q | label | vs ImmuneScore ρ |
|---|---:|---:|---:|---:|---|---:|
| CXCL9 | 515 | −0.080 | 0.071 | 0.095 | NULL | +0.65 |
| CXCL10 | 515 | −0.015 | 0.727 | 0.727 | NULL | +0.64 |
| CXCL13 | 515 | −0.102 | **0.021** | **0.042** | WEAK_NEGATIVE | +0.60 |
| CXCR3 | 515 | +0.133 | **0.0026** | **0.010** | WEAK_POSITIVE | +0.72 |

CXCR3 residual is a **suppressor**: unadjusted null, leftover positive after removing the shared infiltrate axis. Do not quote +0.133 without the unadjusted +0.053 (p = 0.228). CXCL9 residual loses the unadjusted BH hit. CXCL13 leftover stays weakly negative.

Q4 vs Q1 of the OLS ESTIMATE residual is NS for all four (CXCL9 p=0.070; CXCL10 p=0.51; CXCL13 p=0.11; CXCR3 p=0.13). Do not upgrade the continuous residual on that tail.

---

## Methods (this slice)

- Expression: UCSC Xena HiSeqV2 `https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz`. Primary `-01` only; one column per barcode already.
- ESTIMATE: MD Anderson RNAseqV2 `https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt`. We do not recompute ESTIMATE.
- Predictor: CLDN4 only. TJ7 is not scored. TACSTD2 is not a gate or covariate. CD8A / CD274 are not endpoints here.
- Quartiles: `pd.qcut(rank(method="first"), 4)` so ties still fill four bins.
- Primary tests: two-sided Spearman on all 515; two-sided MWU Q4 vs Q1 (129 vs 129). FDR = BH inside the four unadjusted Spearmans, and separately inside the four ESTIMATE partials.
- ESTIMATE residual: first-order partial Spearman | `ESTIMATE_score`. Sensitivity: partial | ImmuneScore; OLS residual of Y on ESTIMATE_score; OLS residual of both CLDN4 and Y.

---

## Already known (not this extra)

| Claim | Result | Source |
|---|---|---|
| CLDN4 Q4 vs Q1 ImmuneScore / CD8 / CD274 | CD8 small negative; ImmuneScore NS; CD274 null; n=515 / 129 vs 129 | PR #352 |
| TJ7 vs CD8 in TCGA-LUAD | ρ = −0.31, p = 3.3×10⁻¹³, n=515 | PR #90 |

This folder asks whether **CLDN4 alone** also moves the CXCL9 / CXCL10 / CXCL13 / CXCR3 axis. Weak CXCL9 / CXCL13 yes. CXCL10 / unadjusted CXCR3 no.

---

## Extra figures

- `methods/tcga_luad_cldn4_cxcl/figures/fig1_spearman_forest.png` — unadjusted ρ vs ESTIMATE partial
- `methods/tcga_luad_cldn4_cxcl/figures/fig2_scatter.png` — CLDN4 vs each chemokine
- `methods/tcga_luad_cldn4_cxcl/figures/fig3_q4_vs_q1.png` — Q4 vs Q1 violins
- `methods/tcga_luad_cldn4_cxcl/figures/fig4_estimate_residual.png` — CLDN4 vs ESTIMATE residual

Headline table: `methods/tcga_luad_cldn4_cxcl/tables/one_row.tsv`.

---

## What this does not claim

- It does not re-test Q4 vs CD8 / ImmuneScore / CD274 (PR #352).
- It does not re-test TJ7 vs CD8 (PR #69 / #90).
- It does not use a TACSTD2-high / TACSTD2-low gate.
- It does not test ICI response, PFS, or OS.
- It does not treat n=502 as the analysis n.
- It does not treat n=515 as the Q4 vs Q1 n (that is 129 vs 129).
- It does not call CLDN4-high LUAD chemokine-excluded or CXCR3-high.
- ESTIMATE residual is not an independent infiltrate control for genes that track ImmuneScore at |ρ| ≥ 0.60.

---

## 中文摘要

只补公开 **TCGA-LUAD** RNA 上 **CLDN4 单基因** 对 **CXCL9 / CXCL10 / CXCL13 / CXCR3** 的连续 Spearman、Q4 vs Q1 和 ESTIMATE residual。不重审 PR #352 的 Q4 vs CD8，也不重算 TJ7 vs CD8。不做 TACSTD2 门控。

- 诚实 n = **515**（HiSeqV2 ∩ 官方 ESTIMATE 原发瘤）。Q4 vs Q1 = **129 vs 129**。不要写 n=502。
- 未校正：CLDN4–CXCL9 ρ = −0.111（p=0.011, q=0.023）；CLDN4–CXCL13 ρ = −0.129（p=0.0034, q=0.013）。CXCL10 / CXCR3 为 NULL。
- ESTIMATE residual：CXCL13 仍弱负；CXCL9 residual 不再显著；CXCR3 leftover 转为弱正（suppressor，勿单独引用）。
- 效应都小。不能升级成 “CLDN4-high LUAD 趋化因子冷”。不是 ICI 队列。

---

## Outputs

- `tables/one_row.tsv` — headline one-row table
- `tables/n_table.tsv` — honest n
- `tables/spearman.tsv` — continuous Spearman
- `tables/q4_vs_q1.tsv` — MWU Q4 vs Q1 (raw + residual)
- `tables/estimate_residual.tsv` — partial / OLS residual
- `tables/estimate_context.tsv` — CLDN4 and CXCL vs ImmuneScore / ESTIMATE_score
- `tables/sample_scores.tsv`
- `tables/summary.json`
- `figures/fig1_spearman_forest.png` … `fig4_estimate_residual.png`

```bash
python3 methods/tcga_luad_cldn4_cxcl/download.py
pip install -r methods/tcga_luad_cldn4_cxcl/requirements.txt
python3 methods/tcga_luad_cldn4_cxcl/analyze.py
```

The slim gene matrix `harvested/hiseqv2_cldn4_cxcl.tsv` plus the committed ESTIMATE table is enough to reproduce without re-downloading HiSeqV2.

# GSE207422 — malignant CLDN4 only (not dual-high)

**Dual-high was stopped.** This slice scores **malignant CLDN4** only. TACSTD2 is a companion gene and is **never a gate**. The given A3 TACSTD2 analysis (NMPR vs MPR and vs T/NK on this same public UMI) is not re-argued.

**Verdict (n=12 honest):** malignant CLDN4 does not significantly separate NMPR vs MPR and does not significantly anti-correlate with T/NK, CXCL13+, or cyto-high T. Immune correlations are weakly negative (A3-malignant |ρ|≤0.33, n=10, all p≥0.35; epithelial |ρ|≤0.44, n=12, all p≥0.15). NMPR>MPR on A3-malignant mean is Δ=+0.97, exact p=0.18, **n=8 vs 2** (P11/P14 empty; P06 is one cell with CLDN4 UMI=0). Complete-case epithelial CLDN4 is null (Δ=+0.005, p=0.93, n=8 vs 4). Q4 vs Q1 epithelial CXCL13+ hits the 3-vs-3 floor (exact p=0.10, Q4 lower); that is not a claim.

## Data and n

- Public GEO UMI only: **92,330** cells × **24,292** genes. Raw GSA-Human HRA001033 was not used. Author CopyKAT / epithelium RDS barcodes are not on GEO.
- **Unit of every test is the 12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). The three pre-treatment biopsies are excluded. Cell-level p-values are not reported.
- Marker lineages (Hu canonical argmax): epithelial 11,019; A3-malignant-like 6,627; T/NK 36,323; CXCL13+ T/NK 8,109; cyto-high T 7,129.
- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the given A3 TACSTD2 slice; **not** CopyKAT). Post patients with 0 A3-malignant cells (NaN, dropped from A3-malignant tests): **P11, P14**. Post patients with 1–9 A3-malignant cells (kept, noisy): **P06, P13, P15**. A3-malignant CLDN4 n=10/12. Epithelial CLDN4 is the complete-case n=12/12 score.

## Definitions

| Item | Rule |
|---|---|
| CLDN4 mean | patient mean `log1p(CP10k)` inside the named compartment |
| CLDN4 %pos | fraction of cells in the compartment with CLDN4 UMI ≥ 1 |
| T/NK | lineage T or NK / all cells |
| CXCL13+ | (T or NK) AND CXCL13 UMI ≥ 1 / all cells |
| cyto-high T | lineage T AND mean log1p(GZMB, GZMA, PRF1, IFNG, NKG7) ≥ pooled T 75th percentile / all cells |
| Q4 vs Q1 | among patients with a finite rank score, Q4 = highest n//4, Q1 = lowest n//4 (n=12 → 3 vs 3) |
| NMPR vs MPR | exact two-sided Wilcoxon by enumerating label assignments; pCR P06 = MPR |
| TACSTD2 | companion only; not used to define dual-high or to filter patients |

## Primary — A3-malignant CLDN4 vs immune (honest n after empty MPR drop)

| CLDN4 score | vs T/NK | vs CXCL13+ | vs cyto-high T |
|---|---|---|---|
| mean log1p(CP10k) | ρ=-0.13, p=0.73, n=10 | ρ=-0.22, p=0.53, n=10 | ρ=-0.30, p=0.4, n=10 |
| %pos | ρ=-0.33, p=0.35, n=10 | ρ=-0.22, p=0.53, n=10 | ρ=-0.33, p=0.35, n=10 |

## Primary — NMPR vs MPR

| Score | Result |
|---|---|
| A3-malignant CLDN4 mean | mean 1.486 vs 0.512 (Δ=+0.974); exact p=0.18; n=8 vs 2 |
| A3-malignant CLDN4 %pos | mean 72.311 vs 38.070 (Δ=+34.241); exact p=0.53; n=8 vs 2 |
| epithelial CLDN4 mean (n=12 complete) | mean 1.513 vs 1.508 (Δ=+0.005); exact p=0.93; n=8 vs 4 |
| epithelial CLDN4 %pos | mean 81.372 vs 89.043 (Δ=-7.672); exact p=0.57; n=8 vs 4 |
| T/NK fraction | mean 0.378 vs 0.513 (Δ=-0.135); exact p=0.57; n=8 vs 4 |
| CXCL13+ fraction | mean 0.096 vs 0.116 (Δ=-0.020); exact p=0.57; n=8 vs 4 |
| cyto-high T fraction | mean 0.080 vs 0.111 (Δ=-0.031); exact p=0.37; n=8 vs 4 |

## Extra — Q4 vs Q1 immune fractions

n=12 honest: epithelial ranks keep all 12 patients (k=3 vs 3). Exact two-sided p at 3 vs 3 cannot go below 0.10. A3-malignant ranks drop patients with 0 malignant cells.

### Epithelial CLDN4 ranks (n=12 complete, k=3)

- Mean-rank Q4: P13(NMPR), P11(MPR), P02(NMPR)
- Mean-rank Q1: P12(NMPR), P03(MPR), P10(NMPR)
- %pos-rank Q4: P07(NMPR), P13(NMPR), P11(MPR)
- %pos-rank Q1: P12(NMPR), P10(NMPR), P04(NMPR)

| Rank | vs T/NK | vs CXCL13+ | vs cyto-high T |
|---|---|---|---|
| epithelial mean | mean 0.370 vs 0.439 (Δ=-0.070); exact p=0.4; n=3 vs 3 | mean 0.026 vs 0.166 (Δ=-0.140); exact p=0.1; n=3 vs 3 | mean 0.048 vs 0.118 (Δ=-0.070); exact p=0.2; n=3 vs 3 |
| epithelial %pos | mean 0.245 vs 0.475 (Δ=-0.230); exact p=0.2; n=3 vs 3 | mean 0.019 vs 0.127 (Δ=-0.108); exact p=0.1; n=3 vs 3 | mean 0.046 vs 0.118 (Δ=-0.072); exact p=0.2; n=3 vs 3 |

### A3-malignant CLDN4 ranks (empty malignant dropped)

- Mean-rank Q4: P09(NMPR), P13(NMPR)
- Mean-rank Q1: P06(MPR), P12(NMPR)

| Rank | vs T/NK | vs CXCL13+ | vs cyto-high T |
|---|---|---|---|
| A3-malignant mean | mean 0.418 vs 0.404 (Δ=+0.014); exact p=1; n=2 vs 2 | mean 0.169 vs 0.108 (Δ=+0.060); exact p=1; n=2 vs 2 | mean 0.123 vs 0.110 (Δ=+0.013); exact p=1; n=2 vs 2 |
| A3-malignant %pos | mean 0.149 vs 0.404 (Δ=-0.255); exact p=0.33; n=2 vs 2 | mean 0.013 vs 0.108 (Δ=-0.095); exact p=0.33; n=2 vs 2 | mean 0.028 vs 0.110 (Δ=-0.082); exact p=0.33; n=2 vs 2 |

## Complete-case epithelial CLDN4 vs immune (n=12)

| CLDN4 score | vs T/NK | vs CXCL13+ | vs cyto-high T |
|---|---|---|---|
| mean log1p(CP10k) | ρ=-0.06, p=0.86, n=12 | ρ=-0.31, p=0.32, n=12 | ρ=-0.18, p=0.57, n=12 |
| %pos | ρ=-0.29, p=0.35, n=12 | ρ=-0.44, p=0.15, n=12 | ρ=-0.21, p=0.51, n=12 |

## Companion TACSTD2 (not a gate; A3 given)

These numbers are reported so CLDN4 can be read next to TACSTD2. They are **not** a dual-high filter and they do not replace the given A3 TACSTD2 write-up.

- A3-malignant CLDN4 mean vs TACSTD2 mean: ρ=0.21, p=0.56, n=10
- Companion TACSTD2 mean vs T/NK: ρ=-0.30, p=0.4, n=10
- Companion TACSTD2 NMPR vs MPR: mean 1.610 vs 1.466 (Δ=+0.143); exact p=1; n=8 vs 2

## Honest limits

1. n=12 (4 vs 8) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. Q4 vs Q1 is 3 vs 3; the smallest exact two-sided p is 0.10.
2. A3-malignant-like empties some MPR residual tumors (normal-lung program). That is reported, not patched by dual-high or by dropping MPR from the header n.
3. This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak into A3-malignant-like.
4. Dual-high (TACSTD2 AND CLDN4) was not run.
5. Expression is log1p(CP10k) from the public UMI, not author Seurat-normalized values.

## Extra figures

- `figures/fig_cldn4_vs_immune.png` — mean and %pos vs T/NK, CXCL13+, cyto-high T
- `figures/fig_cldn4_nmpr_mpr.png` — NMPR vs MPR
- `figures/fig_extra_q4q1.png` — Q4 vs Q1 on epithelial ranks (n=12)
- `figures/fig_extra_q4q1_a3malig.png` — Q4 vs Q1 on A3-malignant ranks
- `figures/fig_extra_companion_tacstd2.png` — TACSTD2 companion, not a gate
- `figures/fig_extra_compartment.png` — per-patient malignant n (empty MPR visible)

## Files

- `per_patient.tsv` — 15 samples; tests use the 12 post rows
- `spearman.tsv` / `nmpr_vs_mpr.tsv` / `q4_vs_q1.tsv` / `quartile_assignment.tsv`
- `sample_metadata.tsv` / `lineage_counts.tsv` / `sanity.json` / `summary.json`
- Scripts: `scripts/download.py`, `scripts/extract.py`, `scripts/analyze.py`

## Reproduce

```bash
python3 methods/gse207422_dualhigh/scripts/download.py
python3 methods/gse207422_dualhigh/scripts/extract.py
python3 methods/gse207422_dualhigh/scripts/analyze.py
```

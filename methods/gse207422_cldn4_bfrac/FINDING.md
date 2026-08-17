# GSE207422 — B-cell / TLS fraction vs malignant CLDN4

ADDITIVE public-UMI slice. **Honest n is the 12 post-treatment patients.** The three pre-treatment biopsies are excluded. Cell-level p-values are not reported. The prior TLS/B meta that listed GSE207422 as n=15 is not reused. TACSTD2 is a companion gene and is **never a gate**. Dual-high was not run. This is not histologic TLS.

**Verdict (n=12 honest):** malignant CLDN4 does not significantly anti-correlate with B-cell fraction, B+plasma, TLS-like cellular fraction, or TLS12 z. A3-malignant CLDN4 mean vs B ρ=-0.10, p=0.78, n=10; vs TLS-like ρ=-0.08, p=0.83, n=10; vs TLS12 z ρ=0.22, p=0.53, n=10. Complete-case epithelial CLDN4 mean vs B ρ=-0.06, p=0.85, n=12. NMPR vs MPR B/TLS fractions are not significant (n=8 vs 4). Q4 vs Q1 is 3 vs 3 on epithelial ranks; the smallest exact two-sided p is 0.10.

## Data and n

- Public GEO UMI only: **92,330** cells × **24,292** genes. Raw GSA-Human HRA001033 was not used. Author CopyKAT / epithelium RDS barcodes are not on GEO.
- **Unit of every test is the 12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). P01/P05/P08 pre-biopsies are written in `per_patient.tsv` and are not tested.
- Marker lineages (Hu canonical argmax): epithelial 11,019; A3-malignant-like 6,627; B 8,082; plasma 3,096; TLS-like (B + CXCL13+ T/NK) 16,191; CXCL13+ T/NK 8,109.
- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the given A3 TACSTD2 slice; **not** CopyKAT). Post patients with 0 A3-malignant cells (NaN, dropped from A3-malignant tests): **P11, P14**. Post patients with 1–9 A3-malignant cells (kept, noisy): **P06, P13, P15**. A3-malignant CLDN4 n=10/12. Epithelial CLDN4 is the complete-case n=12/12 score.
- TLS12 genes present: 12/12 (CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13).

## Definitions

| Item | Rule |
|---|---|
| CLDN4 mean | patient mean `log1p(CP10k)` inside the named compartment |
| CLDN4 %pos | fraction of cells in the compartment with CLDN4 UMI ≥ 1 |
| B fraction | lineage B / all cells |
| B+plasma | (lineage B or plasma) / all cells |
| TLS-like fraction | (lineage B) OR ((T or NK) AND CXCL13 UMI ≥ 1) / all cells |
| TLS12 z | mean of within-post-n=12 z-scores of patient-mean log1p(CP10k) for the Coppola 12-chemokine set. Not a cell fraction. Not histologic TLS. |
| Q4 vs Q1 | among patients with a finite rank score, Q4 = highest n//4, Q1 = lowest n//4 (n=12 → 3 vs 3) |
| NMPR vs MPR | exact two-sided Wilcoxon; pCR P06 = MPR |
| TACSTD2 | companion only; not a gate |

## Primary — A3-malignant CLDN4 vs B / TLS (honest n after empty MPR drop)

| CLDN4 score | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |
|---|---|---|---|---|
| mean log1p(CP10k) | ρ=-0.10, p=0.78, n=10 | ρ=-0.03, p=0.93, n=10 | ρ=-0.08, p=0.83, n=10 | ρ=0.22, p=0.53, n=10 |
| %pos | ρ=-0.35, p=0.33, n=10 | ρ=-0.31, p=0.38, n=10 | ρ=-0.37, p=0.29, n=10 | ρ=0.31, p=0.38, n=10 |

## Complete-case epithelial CLDN4 vs B / TLS (n=12)

| CLDN4 score | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |
|---|---|---|---|---|
| mean log1p(CP10k) | ρ=-0.06, p=0.85, n=12 | ρ=-0.13, p=0.7, n=12 | ρ=-0.08, p=0.81, n=12 | ρ=0.13, p=0.7, n=12 |
| %pos | ρ=-0.32, p=0.31, n=12 | ρ=-0.50, p=0.1, n=12 | ρ=-0.42, p=0.17, n=12 | ρ=0.20, p=0.53, n=12 |

## Primary — NMPR vs MPR

| Score | Result |
|---|---|
| A3-malignant CLDN4 mean | mean 1.486 vs 0.512 (Δ=+0.974); exact p=0.18; n=8 vs 2 |
| epithelial CLDN4 mean (n=12 complete) | mean 1.513 vs 1.508 (Δ=+0.005); exact p=0.93; n=8 vs 4 |
| B fraction | mean 0.092 vs 0.128 (Δ=-0.036); exact p=0.37; n=8 vs 4 |
| B+plasma fraction | mean 0.110 vs 0.157 (Δ=-0.046); exact p=0.68; n=8 vs 4 |
| TLS-like fraction | mean 0.188 vs 0.243 (Δ=-0.055); exact p=0.68; n=8 vs 4 |
| TLS12 z | mean 0.077 vs -0.154 (Δ=+0.232); exact p=0.93; n=8 vs 4 |

## Extra — Q4 vs Q1 B / TLS

n=12 honest: epithelial ranks keep all 12 patients (k=3 vs 3). Exact two-sided p at 3 vs 3 cannot go below 0.10. A3-malignant ranks drop patients with 0 malignant cells.

### Epithelial CLDN4 ranks (n=12 complete, k=3)

- Mean-rank Q4: P13(NMPR), P11(MPR), P02(NMPR)
- Mean-rank Q1: P12(NMPR), P03(MPR), P10(NMPR)
- %pos-rank Q4: P07(NMPR), P13(NMPR), P11(MPR)
- %pos-rank Q1: P12(NMPR), P10(NMPR), P04(NMPR)

| Rank | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |
|---|---|---|---|---|
| epithelial mean | mean 0.151 vs 0.108 (Δ=+0.043); exact p=1; n=3 vs 3 | mean 0.157 vs 0.166 (Δ=-0.009); exact p=1; n=3 vs 3 | mean 0.177 vs 0.274 (Δ=-0.097); exact p=0.7; n=3 vs 3 | mean 0.088 vs 0.045 (Δ=+0.043); exact p=1; n=3 vs 3 |
| epithelial %pos | mean 0.031 vs 0.089 (Δ=-0.058); exact p=0.7; n=3 vs 3 | mean 0.036 vs 0.121 (Δ=-0.085); exact p=0.2; n=3 vs 3 | mean 0.050 vs 0.216 (Δ=-0.166); exact p=0.1; n=3 vs 3 | mean 0.219 vs 0.121 (Δ=+0.098); exact p=1; n=3 vs 3 |

### A3-malignant CLDN4 ranks (empty malignant dropped)

- Mean-rank Q4: P09(NMPR), P13(NMPR)
- Mean-rank Q1: P06(MPR), P12(NMPR)

| Rank | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |
|---|---|---|---|---|
| A3-malignant mean | mean 0.023 vs 0.109 (Δ=-0.086); exact p=0.67; n=2 vs 2 | mean 0.040 vs 0.128 (Δ=-0.088); exact p=1; n=2 vs 2 | mean 0.192 vs 0.217 (Δ=-0.026); exact p=1; n=2 vs 2 | mean 0.932 vs -0.405 (Δ=+1.337); exact p=0.33; n=2 vs 2 |
| A3-malignant %pos | mean 0.013 vs 0.109 (Δ=-0.096); exact p=0.33; n=2 vs 2 | mean 0.019 vs 0.128 (Δ=-0.109); exact p=0.67; n=2 vs 2 | mean 0.026 vs 0.217 (Δ=-0.192); exact p=0.33; n=2 vs 2 | mean 0.207 vs -0.405 (Δ=+0.612); exact p=0.67; n=2 vs 2 |

## Companion (not a gate)

- A3-malignant CLDN4 mean vs TACSTD2 mean: ρ=0.21, p=0.56, n=10
- A3-malignant CLDN4 mean vs MS4A1 mean: ρ=-0.19, p=0.6, n=10
- A3-malignant CLDN4 mean vs CXCL13 mean: ρ=-0.22, p=0.53, n=10
- Companion TACSTD2 mean vs B: ρ=-0.78, p=0.0075, n=10
- Companion TACSTD2 mean vs TLS-like: ρ=-0.88, p=0.00081, n=10

## Honest limits

1. n=12 (4 vs 8) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. Q4 vs Q1 is 3 vs 3; the smallest exact two-sided p is 0.10.
2. A3-malignant-like empties some MPR residual tumors (normal-lung program). That is reported, not patched by falling back to epithelial inside the malignant column or by counting the 3 pre-biopsies.
3. This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak into A3-malignant-like.
4. TLS-like is a marker proxy (B + CXCL13+ T/NK). TLS12 is a 12-chemokine score. Neither is a pathologist TLS call.
5. Dual-high (TACSTD2 AND CLDN4) was not run. Expression is log1p(CP10k) from the public UMI.

## Extra figures

- `figures/fig_cldn4_vs_b_tls.png` — A3-malignant CLDN4 mean/%pos vs B, B+plasma, TLS-like, TLS12
- `figures/fig_b_tls_nmpr_mpr.png` — NMPR vs MPR B / TLS
- `figures/fig_extra_q4q1.png` — Q4 vs Q1 on epithelial ranks (n=12)
- `figures/fig_extra_compartment.png` — per-patient malignant n and B n (empty MPR visible)

## Files

- `per_patient.tsv` — 15 samples; tests use the 12 post rows
- `spearman.tsv` / `nmpr_vs_mpr.tsv` / `q4_vs_q1.tsv` / `quartile_assignment.tsv`
- `sample_metadata.tsv` / `lineage_counts.tsv` / `sanity.json` / `summary.json`
- Scripts: `scripts/download.py`, `scripts/extract.py`, `scripts/analyze.py`

## Reproduce

```bash
python3 methods/gse207422_cldn4_bfrac/scripts/download.py
python3 methods/gse207422_cldn4_bfrac/scripts/extract.py
python3 methods/gse207422_cldn4_bfrac/scripts/analyze.py
```

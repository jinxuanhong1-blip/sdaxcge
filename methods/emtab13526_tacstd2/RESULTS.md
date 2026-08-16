# Honest n / ρ / p — E-MTAB-13526 extra scRNA (not ICI)

Unit = **patient**. CD235a− tumor lanes only (15 lanes, 12 patients, 276,018
QC cells). Malignant-like = marker epithelial minus normal-lung markers.
Eligible: ≥10 malignant-like and ≥20 T/NK cells.

P16 (LUAD) had 2 epithelial cells; P22 (LUAD) had 969 epithelial but only
1 malignant-like cell. Both are excluded from the eligible malignant contrast.

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| Malignant-like TACSTD2 mean log1p UMI vs T/NK fraction | 10 | +0.45 | 0.19 |
| Malignant-like TACSTD2 mean log1p(CP10k) vs T/NK | 10 | +0.16 | 0.65 |
| Malignant-like TACSTD2 %pos vs T/NK | 10 | +0.27 | 0.45 |
| All-epithelial TACSTD2 mean log1p vs T/NK | 11 | +0.091 | 0.79 |
| All 12 CD235a− patients, malig-like TACSTD2 vs T/NK | 12 | −0.025 | 0.94 |
| Malignant-like TACSTD2 vs T/NK, n_malig ≥ 50 | 9 | +0.50 | 0.17 |
| Malignant-like TACSTD2 vs T/NK, LUSC eligible | 5 | +0.70 | 0.19 |
| Malignant-like TACSTD2 vs T/NK, LUAD eligible | 3 | — | too few |
| Malignant-like CLDN4 mean log1p UMI vs T/NK | 10 | +0.42 | 0.23 |
| Malignant-like CLDN4 mean log1p(CP10k) vs T/NK | 10 | +0.14 | 0.70 |
| All-epithelial CLDN4 mean log1p vs T/NK | 11 | +0.055 | 0.87 |
| All 12 CD235a− patients, malig-like CLDN4 vs T/NK | 12 | −0.046 | 0.89 |
| Response / MPR / R vs TACSTD2 or CLDN4 | 0 | — | no public labels |

Epithelial restriction (eligible n=10, paired Wilcoxon greater): median %pos
malignant-like TACSTD2 45.5 vs T/NK 2.71 (p=0.00098); CLDN4 32.6 vs 1.79
(p=0.00098). P24 T/NK cells are 17.8% TACSTD2+ and 37.1% CLDN4+ — likely
ambient RNA from a large CLDN4-high malignant compartment (9,719 cells).

**Read-out.** This treatment-naive atlas does **not** support an inverse
patient-level association of malignant TACSTD2 or CLDN4 with T/NK fraction.
Point estimates are weakly positive and non-significant. It is extra n, not
an ICI-response replication.

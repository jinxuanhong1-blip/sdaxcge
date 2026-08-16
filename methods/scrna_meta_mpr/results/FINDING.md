# Public scRNA MPR/RECIST pool — TACSTD2 and CLDN4

Patient is the unit. Public processed GEO objects only.

## MPR pool (primary)

Four independent rows, GSE241934 IIT and Real kept separate.

| Cohort | Gene | n NMPR vs MPR | mean NMPR vs MPR | Hedges *g* | MWU *p* |
|---|---|---|---|---|---|
| GSE207422 (A3, epithelial) | TACSTD2 | 8 vs 4 | 1.79 vs 1.59 | +0.32 | 0.68 |
| GSE207422 (A3, epithelial) | CLDN4 | 8 vs 4 | 1.64 vs 1.76 | −0.25 | 0.57 |
| GSE241934 IIT | TACSTD2 | 7 vs 4 | 1.58 vs 1.58 | −0.01 | 1.00 |
| GSE241934 IIT | CLDN4 | 7 vs 4 | 1.57 vs 1.71 | −0.48 | 0.53 |
| GSE241934 Real | TACSTD2 | 18 vs 6 | 1.40 vs 1.48 | −0.18 | 1.00 |
| GSE241934 Real | CLDN4 | 18 vs 6 | 1.59 vs 1.40 | +0.41 | 0.45 |
| GSE291670 | TACSTD2 | 3 vs 3 | 0.10 vs 0.23 | −1.03 | 0.20 |
| GSE291670 | CLDN4 | 3 vs 3 | 0.12 vs 0.21 | −0.78 | 0.40 |

**Total n_patients = 53** (36 NMPR, 17 MPR).

Inverse-variance (Hedges *g*, NMPR − MPR):

- TACSTD2: *g* = −0.12 (95% CI −0.71 to 0.47), *p* = 0.68, *I²* = 0, *k* = 4
- CLDN4: *g* = −0.10 (95% CI −0.69 to 0.49), *p* = 0.74, *I²* = 0, *k* = 4

Collapsing IIT+Real to one row leaves the same 53 patients (*k* = 3): TACSTD2 *g* = −0.16 (*p* = 0.59); CLDN4 *g* = −0.06 (*p* = 0.84).

GSE207422 malignant-like sensitivity (8 vs 2; two MPR have 0 malignant-like cells) is not in the primary pool: TACSTD2 *g* = +0.09, *p* = 0.89; CLDN4 *g* = +0.74, *p* = 0.53.

## Combinations (leave-one-in / pairs / leave-one-out / full)

All 15 non-empty subsets of {GSE207422, IIT, Real, GSE291670}. *k*=1 *p* is Mann–Whitney; *k*≥2 *p* is inverse-variance fixed-effect. NMPR>MPR means Hedges *g* > 0. GSE243013 is not in any subset.

**TACSTD2 — 4 / 15 subsets have *g* > 0; 0 / 15 have *g* > 0 and *p* < 0.05.**

| Kind | Subset | n (NMPR vs MPR) | *g* | *p* | NMPR>MPR |
|---|---|---|---|---|---|
| leave-one-in | GSE207422 | 12 (8 vs 4) | +0.32 | 0.68 | yes |
| leave-one-in | IIT | 11 (7 vs 4) | −0.01 | 1.00 | no |
| leave-one-in | Real | 24 (18 vs 6) | −0.18 | 1.00 | no |
| leave-one-in | GSE291670 | 6 (3 vs 3) | −1.03 | 0.20 | no |
| pair | GSE207422 + IIT | 23 (15 vs 8) | +0.16 | 0.72 | yes |
| pair | GSE207422 + Real | 36 (26 vs 10) | +0.005 | 0.99 | yes |
| pair | GSE207422 + GSE291670 | 18 (11 vs 7) | −0.13 | 0.80 | no |
| pair | IIT + Real | 35 (25 vs 10) | −0.12 | 0.75 | no |
| pair | IIT + GSE291670 | 17 (10 vs 7) | −0.36 | 0.48 | no |
| pair | Real + GSE291670 | 30 (21 vs 9) | −0.38 | 0.37 | no |
| leave-one-out | drop GSE291670 | 47 (33 vs 14) | +0.001 | 1.00 | yes |
| leave-one-out | drop Real | 29 (18 vs 11) | −0.08 | 0.83 | no |
| leave-one-out | drop IIT | 42 (29 vs 13) | −0.16 | 0.65 | no |
| leave-one-out | drop GSE207422 | 41 (28 vs 13) | −0.26 | 0.44 | no |
| full | all four | 53 (36 vs 17) | −0.12 | 0.68 | no |

**CLDN4 — 5 / 15 subsets have *g* > 0; 0 / 15 have *g* > 0 and *p* < 0.05.**

| Kind | Subset | n (NMPR vs MPR) | *g* | *p* | NMPR>MPR |
|---|---|---|---|---|---|
| leave-one-in | GSE207422 | 12 (8 vs 4) | −0.25 | 0.57 | no |
| leave-one-in | IIT | 11 (7 vs 4) | −0.48 | 0.53 | no |
| leave-one-in | Real | 24 (18 vs 6) | +0.41 | 0.45 | yes |
| leave-one-in | GSE291670 | 6 (3 vs 3) | −0.78 | 0.40 | no |
| pair | GSE207422 + IIT | 23 (15 vs 8) | −0.36 | 0.42 | no |
| pair | GSE207422 + Real | 36 (26 vs 10) | +0.16 | 0.66 | yes |
| pair | GSE207422 + GSE291670 | 18 (11 vs 7) | −0.43 | 0.38 | no |
| pair | IIT + Real | 35 (25 vs 10) | +0.09 | 0.81 | yes |
| pair | IIT + GSE291670 | 17 (10 vs 7) | −0.59 | 0.25 | no |
| pair | Real + GSE291670 | 30 (21 vs 9) | +0.13 | 0.76 | yes |
| leave-one-out | drop GSE291670 | 47 (33 vs 14) | −0.002 | 0.99 | no |
| leave-one-out | drop Real | 29 (18 vs 11) | −0.45 | 0.25 | no |
| leave-one-out | drop IIT | 42 (29 vs 13) | +0.009 | 0.98 | yes |
| leave-one-out | drop GSE207422 | 41 (28 vs 13) | −0.06 | 0.87 | no |
| full | all four | 53 (36 vs 17) | −0.10 | 0.74 | no |

TACSTD2 NMPR>MPR subsets all include GSE207422 and exclude GSE291670 (except the GSE207422-alone row). CLDN4 NMPR>MPR subsets all include Real. Full four-cohort *g* is negative for both genes. No subset reaches *p* < 0.05.

## RECIST (separate)

GSE205335 author malignant cells, R = PR, NR = SD/PD, NE dropped.

- TACSTD2: 10 vs 6, *g* = +0.08, *p* = 0.96
- CLDN4: 10 vs 6, *g* = −0.35, *p* = 0.56

## Footnote (not in the malignant table)

GSE243013 has public MPR on 243 chemo-IO patients. The public matrix is CD45-sorted immune cells (T/NK, B, myeloid). There are no malignant/epithelial cells to score TACSTD2 or CLDN4.

PRJNA1068179 (Molecular Cancer 2025; 3 NMPR / 2 MPR / 1 pCR) is raw SRA only.

Other 2023–2026 GEO title hits (GSE229353, GSE280232, GSE299111, GSE300685, GSE337519, GSE218989, GSE302284) fail the malignant + public MPR/RECIST rule; see `leftover_inventory.tsv`.

# FINDING — pair GSE189357 + GSE205335: CLDN4-only high-end Mal→T/NK

ADDITIVE **CLDN4 only** on the pairwise combo that already **differs**
(PR #459: malignant CLDN4 %pos vs same-patient T/NK **n=31 ρ=−0.478** p=0.009 I²=0%;
Q4 vs Q1 **r=−0.750** p=0.010, 8/8). **No TACSTD2∩CLDN4 dual-high.**
The given Spearman / Q4 row is not re-ranked. p-values are descriptive.
CellChat R and LIANA were not run.

## Given combo (PR #459; not re-audited)

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p; n_Q1/n_Q4) |
|---|---|---:|---|---|
| GSE189357+GSE205335 | %pos | 31 | −0.478 (0.009, 0%) | −0.750 (0.010; 8/8) |

Members: GSE189357 marker-malignant n=9 + GSE205335 author-malignant n=22.
Malignant definitions differ and are stated on every row.

## Verdict

Outgoing CellChat-style *P* from CLDN4-high malignant cells to same-patient T/NK
is **higher**, not lower, for inhibitory / barrier and several ECM pairs
(NECTIN2–TIGIT n=19 ΔP +0.137 p=3.8e-6; MDK–NCL n=27 ΔP +0.155 p=3.8e-6;
CDH1–KLRG1; LAMA/LAMB–CD44). CXCL16–CXCR6 is also up (n=9, ΔP +0.098, p=0.0039).
CXCL9/10 and CCL4/5 are **not** detected at scale. This is **not** a recruit-down
copy of the given T/NK Spearman. Within-patient Q4 vs Q1 is thin (n=8) because
CLDN4 zeros collapse `qcut` bins; median is the primary high-end split.

## Method (one paragraph)

Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs: 10% truncated
mean of `log1p(CP10k)`, complexes = geometric mean of subunits, detected if
each complex has expressing-cell fraction ≥ 0.10,
P = (L·R)/(0.5+L·R). Outgoing = malignant → **same-patient** T/NK.
Malignant cells are split **within each patient** by CLDN4 (median and Q4 vs Q1).
A patient enters only with both bins ≥10 cells and T/NK ≥20.
The test is Wilcoxon signed-rank on per-patient P (pair detected on both arms,
n≥6). Unit = patient. GSE189357 malignant = marker-malignant
(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0; GSE205335 malignant = author
`Malignant cells`.

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| Given combo (PR #459 %pos) | **31** | 9 + 22; do not re-audit |
| Given Q4 vs Q1 tails | **8 vs 8** | within-cohort, stacked |
| High-end median paired | **31** | all 9 + 22 given patients pass; 4 extra GEO patients have 0 author-malignant cells (already out of n=22) |
| High-end Q4 vs Q1 paired | **8** | GSE189357 1 (TD9) + GSE205335 7; 23 given patients fail `qcut` (CLDN4 zeros / ties) |

Cells are not n. The given n=31 is the combo Spearman n, not the communication n.
Patients that lack a high bin, a low bin, or T/NK are out of the paired test.
Q4 vs Q1 within-patient is thinner than median and is not the primary claim.

## Primary ligand table — within-patient paired (outgoing)

Detect-gated pairs (median): 68 (45 with p<0.05). Q4 vs Q1: 7 (3 with p<0.05).

### Median split

| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | CD55–ADGRE5 ** | 26 | 9/17 | 0.587 | 0.343 | +0.244 | 2.98e-08 |
| other | APP–CD74 ** | 27 | 8/19 | 0.634 | 0.392 | +0.242 | 3.73e-07 |
| other | MDK–NCL ** | 27 | 7/20 | 0.688 | 0.533 | +0.155 | 3.77e-06 |
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 19 | 6/13 | 0.153 | 0.015 | +0.137 | 3.81e-06 |
| other | ICAM1–SPN ** | 21 | 9/12 | 0.152 | 0.060 | +0.092 | 4.77e-06 |
| other | ICAM1–ITGAL ** | 19 | 9/10 | 0.131 | 0.067 | +0.064 | 1.14e-05 |
| other | LAMB3–CD44 ** | 17 | 3/14 | 0.382 | 0.038 | +0.344 | 1.53e-05 |
| other | LAMA5–CD44 ** | 16 | 4/12 | 0.374 | 0.099 | +0.275 | 3.05e-05 |
| other | APP–SORL1 ** | 22 | 7/15 | 0.167 | 0.043 | +0.124 | 4.20e-05 |
| other | CDH1–KLRG1 ** | 14 | 2/12 | 0.065 | 0.019 | +0.046 | 1.22e-04 |
| other | LAMC1–CD44 ** | 13 | 2/11 | 0.187 | 0.046 | +0.141 | 4.88e-04 |
| other | NECTIN2–CD226 ** | 11 | 3/8 | 0.040 | 0.005 | +0.035 | 9.77e-04 |

### Q4 vs Q1 split

| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | APP–CD74 ** | 8 | 1/7 | 0.589 | 0.473 | +0.116 | 0.00781 |
| other | MDK–NCL ** | 8 | 1/7 | 0.794 | 0.727 | +0.067 | 0.0156 |
| other | APP–SORL1 ** | 6 | 1/5 | 0.153 | 0.058 | +0.095 | 0.0312 |
| other | PPIA–BSG | 8 | 1/7 | 0.696 | 0.666 | +0.030 | 0.195 |
| other | CD99–CD99 | 6 | 0/6 | 0.387 | 0.400 | -0.013 | 0.219 |
| other | MIF–CD74_CD44 | 8 | 1/7 | 0.787 | 0.789 | -0.002 | 0.312 |
| other | MIF–CD74_CXCR4 | 8 | 1/7 | 0.849 | 0.851 | -0.001 | 0.461 |

## Focused pairs (always shown if subunits exist)

MHC-I (HLA-A/B/C–CD8), T-recruit (CXCL9/10/16, CCL4/5), CD274–PDCD1, NECTIN2–TIGIT.
These rows are the same test; they are not a second discovery pass.

### Median

| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 19 | 6/13 | 0.153 | 0.015 | +0.137 | 3.81e-06 |
| T-recruit | CXCL16–CXCR6 ** | 9 | 1/8 | 0.117 | 0.018 | +0.098 | 0.00391 |
| MHC-I | HLA-E–CD94:NKG2C ** | 12 | 1/11 | 0.276 | 0.131 | +0.145 | 0.00342 |
| MHC-I | HLA-E–KLRC2 ** | 12 | 1/11 | 0.170 | 0.045 | +0.125 | 0.00342 |
| MHC-I | HLA-F–CD8B ** | 22 | 7/15 | 0.079 | 0.022 | +0.057 | 0.00927 |
| MHC-I | HLA-E–CD8B ** | 25 | 7/18 | 0.222 | 0.142 | +0.080 | 0.00964 |
| MHC-I | HLA-E–KLRK1 ** | 20 | 0/20 | 0.403 | 0.290 | +0.113 | 0.0107 |
| MHC-I | HLA-F–CD8A ** | 23 | 8/15 | 0.120 | 0.074 | +0.046 | 0.0135 |
| MHC-I | HLA-E–CD8A ** | 27 | 8/19 | 0.342 | 0.293 | +0.049 | 0.0299 |
| MHC-I | HLA-A–CD8B | 25 | 7/18 | 0.325 | 0.289 | +0.035 | 0.0551 |
| MHC-I | HLA-A–CD8A | 27 | 8/19 | 0.475 | 0.413 | +0.062 | 0.0906 |
| CD274–PDCD1 | CD274–PDCD1 | 2 | 0/2 | 0.032 | 0.000 | +0.032 | 1 (below gate) |
| T-recruit | CXCL9–CXCR3 | 0 | 0/0 | — | — | — | not detected |
| T-recruit | CXCL10–CXCR3 | 1 | 0/1 | 0.153 | 0.211 | −0.057 | 1 (below gate) |
| T-recruit | CCL5–CCR5 | 2 | 0/2 | 0.001 | 0.004 | −0.003 | 1 (below gate) |

### Q4 vs Q1

| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| CD274–PDCD1 | CD274–PDCD1 | 0 | 0/0 | nan | nan | +nan | NA |
| MHC-I | HLA-E–KLRK1 | 5 | 0/5 | 0.244 | 0.141 | +0.104 | 0.0625 |
| MHC-I | HLA-A–CD8A | 4 | 0/4 | 0.807 | 0.770 | +0.037 | 0.125 |
| MHC-I | HLA-A–CD8B | 4 | 1/3 | 0.671 | 0.621 | +0.050 | 0.125 |
| MHC-I | HLA-B–CD8A | 4 | 0/4 | 0.769 | 0.733 | +0.036 | 0.125 |
| MHC-I | HLA-C–CD8A | 4 | 0/4 | 0.797 | 0.773 | +0.024 | 0.125 |
| MHC-I | HLA-C–CD8B | 4 | 1/3 | 0.658 | 0.626 | +0.032 | 0.125 |
| MHC-I | HLA-E–CD8A | 4 | 0/4 | 0.530 | 0.439 | +0.091 | 0.125 |
| MHC-I | HLA-E–CD8B | 4 | 1/3 | 0.382 | 0.310 | +0.072 | 0.125 |
| MHC-I | HLA-B–CD8B | 4 | 1/3 | 0.623 | 0.581 | +0.043 | 0.25 |
| MHC-I | HLA-E–CD94:NKG2C | 3 | 0/3 | 0.140 | 0.114 | +0.026 | 0.25 |
| MHC-I | HLA-E–KLRC2 | 3 | 0/3 | 0.037 | 0.030 | +0.008 | 0.25 |
| MHC-I | HLA-F–CD8A | 3 | 0/3 | 0.568 | 0.423 | +0.145 | 0.25 |
| MHC-I | HLA-F–CD8B | 4 | 1/3 | 0.293 | 0.224 | +0.069 | 0.25 |
| MHC-I | HLA-B–KIR3DL2 | 1 | 0/1 | 0.015 | 0.012 | +0.003 | 1 |
| MHC-I | HLA-E–CD94:NKG2A | 1 | 0/1 | 0.025 | 0.006 | +0.020 | 1 |
| MHC-I | HLA-E–KLRC1 | 1 | 0/1 | 0.006 | 0.001 | +0.005 | 1 |
| MHC-I | HLA-F–KIR3DL2 | 1 | 0/1 | 0.003 | 0.002 | +0.001 | 1 |
| MHC-I | HLA-G–CD8A | 1 | 0/1 | 0.182 | 0.026 | +0.156 | 1 |
| MHC-I | HLA-G–CD8B | 1 | 0/1 | 0.101 | 0.013 | +0.088 | 1 |

## Companion — between-patient on the given 31-patient combo

All-malignant → same-patient T/NK. Quartiles / median are **within cohort**,
then stacked. This is the slice whose T/NK association is already given
(not re-ranked). Test = Mann–Whitney on per-patient Mal→T/NK P
(detected ≥3 per arm).

### Between-patient median

| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | LAMA5–CD44 ** | 27 | 8/19 | 0.194 | 0.069 | +0.125 | 0.0186 |
| other | ICAM1–ITGAL ** | 23 | 9/14 | 0.135 | 0.065 | +0.070 | 0.0247 |
| other | LAMB3–CD44 ** | 25 | 8/17 | 0.221 | 0.141 | +0.080 | 0.0306 |
| other | MDK–NCL ** | 31 | 9/22 | 0.690 | 0.343 | +0.347 | 0.038 |
| other | ICAM1–SPN ** | 25 | 9/16 | 0.141 | 0.088 | +0.053 | 0.0471 |
| other | LAMB2–CD44 | 21 | 8/13 | 0.180 | 0.037 | +0.143 | 0.0528 |
| other | LAMA3–CD44 | 15 | 3/12 | 0.191 | 0.059 | +0.132 | 0.0541 |
| other | APP–CD74 | 30 | 9/21 | 0.569 | 0.347 | +0.223 | 0.106 |
| other | HLA-DQA1–CD4 | 16 | 9/7 | 0.036 | 0.084 | -0.048 | 0.13 |
| other | LGALS9–P4HB | 23 | 5/18 | 0.247 | 0.031 | +0.216 | 0.131 |
| other | CLEC2D–KLRB1 | 15 | 1/14 | 0.003 | 0.012 | -0.009 | 0.145 |
| MHC-I | HLA-E–KLRC2 | 12 | 1/11 | 0.042 | 0.253 | -0.210 | 0.149 |

### Between-patient Q4 vs Q1

| focus | pair | n | n_189357/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | CD99–CD99 ** | 16 | 5/11 | 0.285 | 0.525 | -0.240 | 0.0164 |
| other | MDK–NCL ** | 17 | 5/12 | 0.754 | 0.343 | +0.411 | 0.0274 |
| other | FN1–CD44 ** | 9 | 4/5 | 0.052 | 0.116 | -0.065 | 0.0317 |
| other | JAG1–NOTCH1 | 7 | 1/6 | 0.016 | 0.001 | +0.016 | 0.0571 |
| other | GDF15–TGFBR2 | 9 | 5/4 | 0.110 | 0.024 | +0.086 | 0.0952 |
| other | ICAM1–ITGAL | 12 | 5/7 | 0.135 | 0.056 | +0.079 | 0.1 |
| other | THBS3–CD47 | 6 | 0/6 | 0.072 | 0.001 | +0.070 | 0.1 |
| other | CLEC2B–KLRB1 | 11 | 5/6 | 0.009 | 0.055 | -0.046 | 0.133 |
| other | LAMA3–CD44 | 9 | 2/7 | 0.314 | 0.047 | +0.267 | 0.167 |
| other | LAMB2–CD44 | 10 | 5/5 | 0.154 | 0.022 | +0.132 | 0.171 |
| MHC-I | HLA-E–KLRK1 | 10 | 0/10 | 0.168 | 0.406 | -0.238 | 0.171 |
| other | LAMA5–CD44 | 13 | 4/9 | 0.158 | 0.083 | +0.074 | 0.181 |

## What this does not say

- It does not build a TACSTD2∩CLDN4 dual-high score.
- It does not re-audit the given PR #459 Spearman / Q4 row.
- It does not treat this pair as GSE131907+GSE205335 (that is a different folder).
- Cell-pooled stacked means are not the test.
- GSE205335 Q4 remains SCLC-heavy on the given between-patient tails.
- No FASTQ. No CellChat R. No LIANA.

## Files

- `results/ligand_table.tsv` — patient-level paired LR (n / Δ / p)
- `results/ligand_table_key.tsv` — MHC-I / recruit / CD274–PDCD1 / NECTIN2–TIGIT
- `results/ligand_table_between.tsv` — given-combo between-patient companion
- `results/eligibility.tsv` — honest paired n
- `figures/fig_extra_ligand_table.png` — extra median ΔP bars
- `figures/fig_honest_n.png` — extra paired-n scatter
- `figures/fig_n_per_patient.png` — extra per-patient cell floors
- `figures/fig_given_combo_cldn4_tnk.png` — extra given-combo scatter (not re-ranked)
- `figures/fig_key_pairs_paired_median.png` — extra paired P lines
- `METHODS.md` — labels, Hill probability, floors


# FINDING — merged GSE131907 + GSE205335 CellChat-style CLDN4-high → T/NK

ADDITIVE **CLDN4 only** on the merged public UMI slice that already differs
(PR #320: malignant CLDN4 Q4 vs Q1 vs same-patient T/NK n=23 r=−0.705 p=0.0003;
continuous n=43 ρ=−0.479). **No TACSTD2∩CLDN4 dual-high.** GSE207422 was not run.
p-values are descriptive. CellChat R was not run. Probability is Jin et al. 2021.

**English.** Within the same patient, CLDN4-high malignant cells send higher
outgoing Hill *P* to that patient's T/NK than CLDN4-low malignant cells for
**NECTIN2–TIGIT** (median n=39, Δ=+0.094, p=3.6e-12) and classical **MHC-I**
(HLA-E–CD8A n=46, Δ=+0.075, p=5.1e-7; HLA-A/B/C–CD8 same direction). The only
T-recruit pair that clears the detect gate is **CXCL16–CXCR6** (n=25, Δ=+0.034,
p=2.1e-5) — higher from the high arm, not lower. **CXCL9/10 and CCL4/5** are
almost never detected (n≤3). **CD274–PDCD1** is thin (n=6, Δ=+0.009, p=0.094).
Q4 vs Q1 eligible n=27 (20 GSE131907 / 7 GSE205335; GSE205335 tails are thin
because of CLDN4 ties). NECTIN2–TIGIT repeats (n=22, Δ=+0.089, p=1.4e-6). The given
PR #320 T/NK association is not re-ranked. GSE131907 `PVRL2` was aliased to
`NECTIN2`.

## Method (one paragraph)

Jin et al. 2021 Hill probability on CellChatDB v2 protein pairs: 10% truncated
mean of `log1p(CP10k)`, complexes = geometric mean of subunits, detected if
each complex has expressing-cell fraction ≥ 0.10,
P = (L·R)/(0.5+L·R). Outgoing = malignant → **same-patient** T/NK.
Malignant cells are split **within each patient** by CLDN4 (median and Q4 vs Q1).
A patient enters only with both bins ≥10 cells and T/NK ≥20.
The test is Wilcoxon signed-rank on per-patient P (pair detected on both arms,
n≥6). Unit = GEO patient (GSE131907 tumor origins pooled;
PE / nLung / nLN dropped). GSE131907 malignant = author `Malignant cells` + tS1/tS2/tS3;
GSE205335 malignant = author `Malignant cells`.

## Honest paired n

| split | considered (mal>0 and T/NK>0) | eligible | GSE131907 | GSE205335 | out (floor / thin bins) |
|---|---:|---:|---:|---:|---:|
| median | 54 | 51 | 29 | 22 | 3 |
| Q4 vs Q1 | 54 | 27 | 20 | 7 | 27 |

Cells are not n. nLung / nLN / PE-only patients (0 tumor malignant) are not
in the considered column. Q4 vs Q1 loses patients when `qcut` cannot form
four ranks (CLDN4 ties, often a large zero mass).

## Primary ligand table — within-patient paired (outgoing)

Detect-gated pairs (median): 77 (58 with p<0.05). Q4 vs Q1: 58 (42 with p<0.05).

### Median split

| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 39 | 26/13 | 0.181 | 0.087 | +0.094 | 3.64e-12 |
| other | CD55–ADGRE5 ** | 44 | 27/17 | 0.356 | 0.258 | +0.098 | 2.35e-11 |
| other | LAMB3–CD44 ** | 40 | 26/14 | 0.430 | 0.258 | +0.172 | 7.82e-11 |
| other | LAMA5–CD44 ** | 37 | 25/12 | 0.293 | 0.121 | +0.172 | 1.02e-10 |
| other | ICAM1–SPN ** | 36 | 24/12 | 0.193 | 0.123 | +0.069 | 7.28e-10 |
| other | APP–CD74 ** | 47 | 28/19 | 0.682 | 0.600 | +0.082 | 4.24e-09 |
| other | LAMC2–CD44 ** | 28 | 22/6 | 0.191 | 0.057 | +0.135 | 1.04e-07 |
| other | CDH1–KLRG1 ** | 24 | 12/12 | 0.064 | 0.026 | +0.038 | 1.19e-07 |
| other | IGFBP3–TMEM219 ** | 25 | 16/9 | 0.239 | 0.115 | +0.124 | 2.98e-07 |
| MHC-I | HLA-E–CD8A ** | 46 | 27/19 | 0.339 | 0.264 | +0.075 | 5.10e-07 |
| other | LAMC1–CD44 ** | 36 | 25/11 | 0.149 | 0.086 | +0.063 | 7.11e-07 |
| MHC-I | HLA-E–CD8B ** | 41 | 23/18 | 0.238 | 0.211 | +0.027 | 9.40e-07 |

### Q4 vs Q1 split

| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | LAMB3–CD44 ** | 21 | 18/3 | 0.462 | 0.345 | +0.117 | 9.54e-07 |
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 22 | 18/4 | 0.177 | 0.088 | +0.089 | 1.43e-06 |
| other | CD55–ADGRE5 ** | 23 | 19/4 | 0.359 | 0.291 | +0.068 | 1.67e-06 |
| other | LAMA5–CD44 ** | 24 | 19/5 | 0.214 | 0.093 | +0.121 | 2.26e-06 |
| MHC-I | HLA-E–CD8A ** | 22 | 18/4 | 0.339 | 0.274 | +0.065 | 4.77e-06 |
| MHC-I | HLA-A–CD8A ** | 22 | 18/4 | 0.559 | 0.526 | +0.032 | 9.06e-06 |
| other | APP–CD74 ** | 27 | 20/7 | 0.729 | 0.619 | +0.110 | 9.54e-06 |
| MHC-I | HLA-F–CD8A ** | 20 | 17/3 | 0.105 | 0.052 | +0.053 | 2.67e-05 |
| MHC-I | HLA-A–CD8B ** | 20 | 17/3 | 0.476 | 0.437 | +0.039 | 3.62e-05 |
| MHC-I | HLA-E–CD8B ** | 20 | 17/3 | 0.224 | 0.188 | +0.036 | 3.62e-05 |
| other | LGALS9–P4HB ** | 20 | 16/4 | 0.237 | 0.147 | +0.090 | 8.20e-05 |
| MHC-I | HLA-F–CD8B ** | 19 | 16/3 | 0.056 | 0.028 | +0.028 | 9.54e-05 |

## Focused pairs (same test; not a second discovery pass)

MHC-I (HLA-A/B/C/E–CD8), T-recruit (CXCL9/10/16, CCL4/5), CD274–PDCD1,
NECTIN2–TIGIT. Rows with n_paired=0 (never detected on both arms) are omitted
here and kept in `ligand_table_key.tsv`. CXCL9–CXCR3 was not detected.

### Median

| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| MHC-I | HLA-E–CD8A ** | 46 | 27/19 | 0.339 | 0.264 | +0.075 | 5.10e-07 |
| MHC-I | HLA-E–CD8B ** | 41 | 23/18 | 0.238 | 0.211 | +0.027 | 9.40e-07 |
| MHC-I | HLA-A–CD8A ** | 46 | 27/19 | 0.551 | 0.529 | +0.023 | 2.55e-05 |
| MHC-I | HLA-C–CD8A ** | 46 | 27/19 | 0.538 | 0.488 | +0.050 | 2.86e-05 |
| MHC-I | HLA-E–KLRC2 ** | 15 | 4/11 | 0.179 | 0.102 | +0.078 | 6.10e-05 |
| MHC-I | HLA-C–CD8B ** | 41 | 23/18 | 0.419 | 0.415 | +0.004 | 7.86e-05 |
| MHC-I | HLA-A–CD8B ** | 41 | 23/18 | 0.465 | 0.431 | +0.034 | 1.02e-04 |
| MHC-I | HLA-E–KLRC1 ** | 14 | 9/5 | 0.123 | 0.109 | +0.015 | 1.22e-04 |
| MHC-I | HLA-B–CD8A ** | 46 | 27/19 | 0.542 | 0.512 | +0.030 | 3.10e-04 |
| MHC-I | HLA-F–CD8B ** | 36 | 21/15 | 0.072 | 0.044 | +0.028 | 0.00106 |
| MHC-I | HLA-B–CD8B ** | 41 | 23/18 | 0.433 | 0.460 | -0.026 | 0.00113 |
| MHC-I | HLA-F–CD8A ** | 39 | 24/15 | 0.120 | 0.073 | +0.047 | 0.00132 |
| MHC-I | HLA-E–KLRK1 ** | 20 | 0/20 | 0.403 | 0.290 | +0.113 | 0.0107 |
| MHC-I | HLA-G–CD8B | 4 | 2/2 | 0.054 | 0.024 | +0.030 | 0.125 |
| MHC-I | HLA-B–KIR3DL2 | 3 | 1/2 | 0.050 | 0.037 | +0.013 | 0.25 |
| MHC-I | HLA-G–CD8A | 5 | 3/2 | 0.153 | 0.037 | +0.116 | 0.312 |
| MHC-I | HLA-C–KIR2DL3 | 1 | 0/1 | 0.200 | 0.108 | +0.091 | 1 |
| MHC-I | HLA-F–KIR3DL2 | 2 | 1/1 | 0.005 | 0.002 | +0.003 | 1 |
| MHC-I | ULBP2–KLRK1 | 2 | 0/2 | 0.043 | 0.025 | +0.018 | 1 |
| T-recruit | CXCL16–CXCR6 ** | 25 | 17/8 | 0.049 | 0.015 | +0.034 | 2.07e-05 |
| T-recruit | CCL4–CCR5 | 3 | 3/0 | 0.014 | 0.013 | +0.001 | 0.25 |
| T-recruit | CXCL10–CXCR3 | 3 | 2/1 | 0.009 | 0.010 | -0.001 | 0.25 |
| T-recruit | CCL5–CCR1 | 3 | 2/1 | 0.004 | 0.001 | +0.002 | 0.5 |
| T-recruit | CCL5–CCR4 | 1 | 0/1 | 0.002 | 0.001 | +0.001 | 1 |
| CD274–PDCD1 | CD274–PDCD1 | 6 | 4/2 | 0.010 | 0.002 | +0.009 | 0.0938 |
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 39 | 26/13 | 0.181 | 0.087 | +0.094 | 3.64e-12 |

### Q4 vs Q1

| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| MHC-I | HLA-E–CD8A ** | 22 | 18/4 | 0.339 | 0.274 | +0.065 | 4.77e-06 |
| MHC-I | HLA-A–CD8A ** | 22 | 18/4 | 0.559 | 0.526 | +0.032 | 9.06e-06 |
| MHC-I | HLA-F–CD8A ** | 20 | 17/3 | 0.105 | 0.052 | +0.053 | 2.67e-05 |
| MHC-I | HLA-A–CD8B ** | 20 | 17/3 | 0.476 | 0.437 | +0.039 | 3.62e-05 |
| MHC-I | HLA-E–CD8B ** | 20 | 17/3 | 0.224 | 0.188 | +0.036 | 3.62e-05 |
| MHC-I | HLA-F–CD8B ** | 19 | 16/3 | 0.056 | 0.028 | +0.028 | 9.54e-05 |
| MHC-I | HLA-B–CD8A ** | 22 | 18/4 | 0.531 | 0.489 | +0.042 | 1.77e-04 |
| MHC-I | HLA-C–CD8B ** | 20 | 17/3 | 0.441 | 0.436 | +0.006 | 2.10e-04 |
| MHC-I | HLA-C–CD8A ** | 22 | 18/4 | 0.531 | 0.466 | +0.064 | 2.13e-04 |
| MHC-I | HLA-B–CD8B ** | 20 | 17/3 | 0.432 | 0.430 | +0.002 | 8.51e-04 |
| MHC-I | HLA-E–KLRC1 ** | 7 | 6/1 | 0.108 | 0.092 | +0.017 | 0.0156 |
| MHC-I | HLA-E–KLRC2 ** | 6 | 3/3 | 0.115 | 0.088 | +0.027 | 0.0312 |
| MHC-I | HLA-E–KLRK1 | 5 | 0/5 | 0.244 | 0.141 | +0.104 | 0.0625 |
| MHC-I | HLA-G–CD8B | 3 | 2/1 | 0.015 | 0.005 | +0.011 | 0.5 |
| MHC-I | HLA-G–CD8A | 4 | 3/1 | 0.093 | 0.016 | +0.077 | 0.875 |
| MHC-I | HLA-B–KIR3DL2 | 1 | 0/1 | 0.015 | 0.012 | +0.003 | 1 |
| MHC-I | HLA-F–KIR3DL2 | 1 | 0/1 | 0.003 | 0.002 | +0.001 | 1 |
| MHC-I | ULBP2–KLRK1 | 1 | 0/1 | 0.074 | 0.038 | +0.036 | 1 |
| T-recruit | CXCL16–CXCR6 ** | 13 | 12/1 | 0.048 | 0.020 | +0.027 | 4.88e-04 |
| T-recruit | CCL4–CCR5 | 2 | 2/0 | 0.009 | 0.007 | +0.003 | 1 |
| T-recruit | CCL5–CCR1 | 1 | 1/0 | 0.002 | 0.000 | +0.002 | 1 |
| T-recruit | CCL5–CCR5 | 2 | 1/1 | 0.001 | 0.018 | -0.017 | 1 |
| T-recruit | CXCL10–CXCR3 | 1 | 1/0 | 0.004 | 0.011 | -0.007 | 1 |
| CD274–PDCD1 | CD274–PDCD1 | 2 | 2/0 | 0.013 | 0.006 | +0.006 | 1 |
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 22 | 18/4 | 0.177 | 0.088 | +0.089 | 1.43e-06 |

## Companion — between-unit on the locked PR #320 slice

Same 21 GSE131907 author-`Malignant cells` samples (n_mal≥20) + 22 GSE205335
patients. Quartiles / median are **within cohort**, then stacked. This is the
slice whose T/NK association is already given (not re-ranked). Test =
Mann–Whitney on per-unit Mal→T/NK P (detected ≥3 per arm).

### Between-unit median

| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | MDK–NCL ** | 43 | 22/21 | 0.707 | 0.492 | +0.214 | 8.71e-05 |
| other | APP–CD74 ** | 41 | 21/20 | 0.683 | 0.457 | +0.225 | 0.00307 |
| other | LAMB3–CD44 ** | 36 | 18/18 | 0.412 | 0.190 | +0.222 | 0.0184 |
| other | CD69–KLRB1 ** | 13 | 4/9 | 0.005 | 0.123 | -0.118 | 0.0196 |
| other | CLEC2D–KLRB1 ** | 28 | 14/14 | 0.003 | 0.017 | -0.014 | 0.0258 |
| other | LAMA3–CD44 ** | 18 | 8/10 | 0.133 | 0.067 | +0.066 | 0.0343 |
| other | LAMA5–CD44 ** | 40 | 20/20 | 0.208 | 0.087 | +0.121 | 0.0385 |
| other | ICAM1–SPN ** | 31 | 16/15 | 0.151 | 0.100 | +0.051 | 0.0418 |
| other | LGALS9–P4HB ** | 34 | 17/17 | 0.171 | 0.055 | +0.116 | 0.0421 |
| NECTIN2–TIGIT | NECTIN2–TIGIT ** | 37 | 17/20 | 0.154 | 0.095 | +0.059 | 0.0427 |
| other | CLEC2B–KLRB1 | 27 | 12/15 | 0.105 | 0.031 | +0.074 | 0.0923 |
| other | LGALS9–CD44 | 34 | 17/17 | 0.329 | 0.104 | +0.226 | 0.0983 |

### Between-unit Q4 vs Q1

| focus | pair | n | n_131907/n_205335 | median P high | median P low | Δ | p |
|---|---|---:|---|---:|---:|---:|---|
| other | MDK–NCL ** | 23 | 11/12 | 0.724 | 0.460 | +0.264 | 7.96e-04 |
| other | MIF–CD74_CD44 ** | 23 | 11/12 | 0.777 | 0.859 | -0.082 | 0.0247 |
| other | LAMA3–CD44 ** | 11 | 3/8 | 0.314 | 0.047 | +0.267 | 0.0485 |
| other | LAMB3–CD44 | 17 | 7/10 | 0.409 | 0.170 | +0.239 | 0.0553 |
| other | MIF–CD74_CXCR4 | 23 | 11/12 | 0.841 | 0.891 | -0.050 | 0.0694 |
| other | CLEC2B–KLRB1 | 15 | 4/11 | 0.106 | 0.039 | +0.067 | 0.0777 |
| other | LAMA5–CD44 | 20 | 9/11 | 0.201 | 0.092 | +0.109 | 0.0806 |
| other | JAG1–NOTCH1 | 6 | 3/3 | 0.016 | 0.001 | +0.015 | 0.1 |
| MHC-I | HLA-E–KLRC2 | 6 | 3/3 | 0.036 | 0.343 | -0.307 | 0.1 |
| other | LGALS9–P4HB | 17 | 7/10 | 0.339 | 0.031 | +0.307 | 0.109 |
| other | ICAM1–SPN | 15 | 6/9 | 0.136 | 0.095 | +0.041 | 0.145 |
| other | COL6A1–CD44 | 11 | 4/7 | 0.356 | 0.179 | +0.177 | 0.164 |

## What this does not say

- It does not build a TACSTD2∩CLDN4 dual-high score.
- It does not call the thesis a failure. The given T/NK association on this
  merged slice stands (PR #320).
- Cell-pooled stacked means are not the test.
- GSE207422 patient-level CLDN4 vs T/NK is flat and was not run.
- No FASTQ. The 2.86 GB GSE131907 log2TPM text matrix was skipped (UMI exists).

## Files

- `results/ligand_table.tsv` — patient-level paired LR (n / Δ / p)
- `results/ligand_table_key.tsv` — MHC-I / recruit / CD274–PDCD1 / NECTIN2–TIGIT
- `results/ligand_table_between.tsv` — PR #320-slice between-unit companion
- `results/eligibility.tsv` — honest paired n
- `figures/fig_extra_ligand_table.png` — paired median ΔP
- `figures/fig_q4q1_ligand_table.png` — paired Q4 vs Q1 ΔP
- `figures/fig_key_pairs_median.png` — focused-pair bars
- `figures/fig_key_pairs_paired_median.png` / `_q4q1.png` — per-patient high vs low
- `figures/fig_honest_n.png` — both-bin floors
- `figures/fig_pr320_slice_cldn4_tnk.png` — given slice (not re-ranked)
- `figures/fig_between_q4q1.png` — between-unit companion
- `METHODS.md` — labels, Hill probability, floors


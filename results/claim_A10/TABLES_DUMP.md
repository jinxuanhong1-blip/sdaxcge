# Claim A10 — results (honest)

This file is generated from `results/claim_A10/tables/` by `scripts/05_write_report.py`.
The claim, predictions, and decision rules were pre-registered in
`claims/claim_A10.md` before any of these numbers were inspected.

**Claim (operationalised):** ELF3, GRHL1, KLF4 and TFAP2A form a coherent
transcription-factor module that is positively coupled to TACSTD2 and CLDN4,
in a lung-epithelial state where NKX2-1 is down. Expected in human and mouse.

## Overall verdict

_Filled after the per-prediction sections below._

## Data actually used

| dataset | role | n | notes |
|---|---|---|---|
| GTEx-lung (recount3) | human bulk observational | 655 | log2CPM of coverage sums |
| TCGA-LUAD (recount3) | human bulk observational | 601 | log2CPM of coverage sums |
| CELLxGENE Census human lung epithelium | P6 / cell type | 895369 cells | primary data, normal + LUAD |
| CELLxGENE Census mouse lung | P5/P6 | 220674 cells | all primary lung cells in Census |
| GSE129340_H441 | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE129340_H209 | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H358_sg | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H2087_sg | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H441_sg | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H441_sh | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H2087_sh | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H441_sh_set2 | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H2087_CRISPRi | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H358_CRISPRi | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_HCC78_CRISPRi | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_PC9_OE | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE229541_H1975_OE | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE129583_AT1_P5 | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE129583_AT2_P8P9 | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE115899 | perturbation / model | — | see `perturbation_contrasts.csv` |
| GSE145152 | perturbation / model | — | see `perturbation_contrasts.csv` |

### Claim-gene detection (human bulk)

| dataset | gene | mean log2CPM | fraction samples > 1 |
|---|---|---:|---:|
| GTEx-lung | ELF3 | 6.63 | 1.00 |
| GTEx-lung | GRHL1 | 3.19 | 1.00 |
| GTEx-lung | KLF4 | 7.22 | 1.00 |
| GTEx-lung | TFAP2A | 0.86 | 0.28 |
| GTEx-lung | TACSTD2 | 6.20 | 1.00 |
| GTEx-lung | CLDN4 | 5.86 | 1.00 |
| GTEx-lung | NKX2-1 | 6.28 | 1.00 |
| TCGA-LUAD | ELF3 | 8.04 | 1.00 |
| TCGA-LUAD | GRHL1 | 3.93 | 0.99 |
| TCGA-LUAD | KLF4 | 4.61 | 1.00 |
| TCGA-LUAD | TFAP2A | 3.01 | 0.87 |
| TCGA-LUAD | TACSTD2 | 7.71 | 1.00 |
| TCGA-LUAD | CLDN4 | 8.22 | 1.00 |
| TCGA-LUAD | NKX2-1 | 6.63 | 0.97 |

## P1 — the four TFs are a module

Pre-registered: mean pairwise Spearman among ELF3/GRHL1/KLF4/TFAP2A exceeds
expression-matched random 4-gene sets, and |ρ| ≥ 0.3.

| dataset | mean pairwise ρ | matched-null mean | matched-null p95 | empirical p | random-pair percentile |
|---|---:|---:|---:|---:|---:|
| GTEx-lung | 0.157 | 0.014 | 0.200 | 0.084 | 68.0 |
| TCGA-LUAD | 0.099 | 0.031 | 0.172 | 0.1705 | 64.2 |

Pairwise TF–TF (human bulk):

| dataset | pair | ρ | 95% CI | FDR | bg percentile |
|---|---|---:|---|---:|---:|
| GTEx-lung | ELF3–GRHL1 | 0.325 | [0.250, 0.395] | 4.09e-17 | 86.7 |
| GTEx-lung | ELF3–KLF4 | 0.393 | [0.322, 0.460] | 4.85e-25 | 91.7 |
| GTEx-lung | ELF3–TFAP2A | 0.119 | [0.038, 0.199] | 3.29e-03 | 62.4 |
| GTEx-lung | GRHL1–KLF4 | -0.077 | [-0.157, 0.004] | 6.07e-02 | 30.3 |
| GTEx-lung | GRHL1–TFAP2A | 0.198 | [0.118, 0.275] | 6.14e-07 | 73.5 |
| GTEx-lung | KLF4–TFAP2A | -0.014 | [-0.095, 0.067] | 7.30e-01 | 40.2 |
| TCGA-LUAD | ELF3–GRHL1 | 0.285 | [0.205, 0.361] | 3.63e-12 | 90.9 |
| TCGA-LUAD | ELF3–KLF4 | -0.047 | [-0.131, 0.038] | 2.91e-01 | 32.7 |
| TCGA-LUAD | ELF3–TFAP2A | 0.229 | [0.147, 0.308] | 3.51e-08 | 85.3 |
| TCGA-LUAD | GRHL1–KLF4 | -0.146 | [-0.228, -0.062] | 5.73e-04 | 15.6 |
| TCGA-LUAD | GRHL1–TFAP2A | 0.384 | [0.309, 0.454] | 1.14e-21 | 96.5 |
| TCGA-LUAD | KLF4–TFAP2A | -0.112 | [-0.195, -0.028] | 8.67e-03 | 20.6 |

## P2 — module / each TF correlates positively with TACSTD2 and CLDN4

| dataset | pair | ρ | 95% CI | FDR | bg percentile |
|---|---|---:|---|---:|---:|
| GTEx-lung | ELF3–TACSTD2 | 0.409 | [0.339, 0.475] | 3.33e-27 | 92.6 |
| GTEx-lung | ELF3–CLDN4 | 0.651 | [0.602, 0.696] | 2.74e-79 | 99.4 |
| GTEx-lung | GRHL1–TACSTD2 | -0.111 | [-0.190, -0.030] | 6.43e-03 | 25.5 |
| GTEx-lung | GRHL1–CLDN4 | -0.101 | [-0.180, -0.020] | 1.33e-02 | 26.8 |
| GTEx-lung | KLF4–TACSTD2 | 0.320 | [0.245, 0.391] | 1.17e-16 | 86.3 |
| GTEx-lung | KLF4–CLDN4 | 0.385 | [0.313, 0.452] | 5.86e-24 | 91.1 |
| GTEx-lung | TFAP2A–TACSTD2 | -0.017 | [-0.098, 0.064] | 6.80e-01 | 39.7 |
| GTEx-lung | TFAP2A–CLDN4 | 0.076 | [-0.006, 0.156] | 6.47e-02 | 55.4 |
| TCGA-LUAD | ELF3–TACSTD2 | 0.376 | [0.301, 0.447] | 6.59e-21 | 96.2 |
| TCGA-LUAD | ELF3–CLDN4 | 0.567 | [0.507, 0.622] | 5.56e-51 | 99.6 |
| TCGA-LUAD | GRHL1–TACSTD2 | 0.409 | [0.335, 0.477] | 1.10e-24 | 97.4 |
| TCGA-LUAD | GRHL1–CLDN4 | 0.438 | [0.367, 0.504] | 1.16e-28 | 98.1 |
| TCGA-LUAD | KLF4–TACSTD2 | 0.094 | [0.009, 0.177] | 2.89e-02 | 63.3 |
| TCGA-LUAD | KLF4–CLDN4 | -0.219 | [-0.298, -0.137] | 1.39e-07 | 7.7 |
| TCGA-LUAD | TFAP2A–TACSTD2 | 0.278 | [0.198, 0.355] | 1.27e-11 | 90.3 |
| TCGA-LUAD | TFAP2A–CLDN4 | 0.304 | [0.225, 0.379] | 8.50e-14 | 92.4 |

Module score vs targets:

| dataset | y | ρ | 95% CI | FDR |
|---|---|---:|---|---:|
| GTEx-lung | TACSTD2 | 0.284 | [0.208, 0.357] | 1.87e-13 |
| GTEx-lung | CLDN4 | 0.441 | [0.373, 0.504] | 5.12e-32 |
| TCGA-LUAD | TACSTD2 | 0.507 | [0.441, 0.568] | 1.86e-39 |
| TCGA-LUAD | CLDN4 | 0.475 | [0.407, 0.538] | 2.27e-34 |

## P3 — NKX2-1 is anti-correlated with the module and the targets

| dataset | pair | ρ | 95% CI | FDR | bg percentile |
|---|---|---:|---|---:|---:|
| GTEx-lung | ELF3–NKX2-1 | 0.351 | [0.277, 0.420] | 6.65e-20 | 88.8 |
| GTEx-lung | GRHL1–NKX2-1 | -0.278 | [-0.351, -0.201] | 1.01e-12 | 9.0 |
| GTEx-lung | KLF4–NKX2-1 | 0.325 | [0.250, 0.396] | 3.97e-17 | 86.7 |
| GTEx-lung | TFAP2A–NKX2-1 | -0.061 | [-0.142, 0.020] | 1.38e-01 | 32.6 |
| GTEx-lung | TACSTD2–NKX2-1 | 0.558 | [0.500, 0.612] | 3.50e-54 | 98.0 |
| GTEx-lung | CLDN4–NKX2-1 | 0.715 | [0.672, 0.752] | 3.07e-102 | 99.8 |
| TCGA-LUAD | ELF3–NKX2-1 | 0.201 | [0.118, 0.281] | 1.48e-06 | 81.7 |
| TCGA-LUAD | GRHL1–NKX2-1 | 0.117 | [0.033, 0.200] | 6.27e-03 | 67.8 |
| TCGA-LUAD | KLF4–NKX2-1 | -0.383 | [-0.453, -0.308] | 1.43e-21 | 0.9 |
| TCGA-LUAD | TFAP2A–NKX2-1 | -0.218 | [-0.297, -0.136] | 1.57e-07 | 7.8 |
| TCGA-LUAD | TACSTD2–NKX2-1 | -0.034 | [-0.118, 0.051] | 4.46e-01 | 35.5 |
| TCGA-LUAD | CLDN4–NKX2-1 | 0.242 | [0.161, 0.320] | 4.74e-09 | 86.8 |

TCGA-LUAD NKX2-1-low vs NKX2-1-high tumours (quartile split):

| gene | n_low | n_high | log2FC (low−high) | Cliff's δ | FDR |
|---|---:|---:|---:|---:|---:|
| ELF3 | 136 | 136 | -0.612 | -0.276 | 1.58e-04 |
| GRHL1 | 136 | 136 | -0.351 | -0.196 | 7.34e-03 |
| KLF4 | 136 | 136 | 1.465 | 0.683 | 1.83e-21 |
| TFAP2A | 136 | 136 | 0.964 | 0.376 | 3.55e-07 |
| TACSTD2 | 136 | 136 | -0.077 | 0.095 | 1.85e-01 |
| CLDN4 | 136 | 136 | -0.726 | -0.360 | 8.02e-07 |
| EPCAM | 136 | 136 | -0.609 | -0.361 | 8.02e-07 |
| SFTPC | 136 | 136 | -1.511 | -0.215 | 3.45e-03 |
| SCGB1A1 | 136 | 136 | -0.190 | -0.057 | 4.14e-01 |
| KRT5 | 136 | 136 | 1.332 | 0.328 | 6.98e-06 |

## Composition check (bulk stand-in for P6)

If the bulk associations are only 'this biopsy has more airway epithelium',
they should collapse after residualising on EPCAM / SFTPC / SCGB1A1 / lineage markers.

| dataset | covariates | pair | partial ρ | FDR |
|---|---|---|---:|---:|
| GTEx-lung | none | ELF3–TACSTD2 | 0.409 | 6.44e-27 |
| GTEx-lung | none | ELF3–CLDN4 | 0.651 | 7.01e-79 |
| GTEx-lung | none | GRHL1–TACSTD2 | -0.111 | 7.85e-03 |
| GTEx-lung | none | GRHL1–CLDN4 | -0.101 | 1.63e-02 |
| GTEx-lung | none | KLF4–TACSTD2 | 0.320 | 2.07e-16 |
| GTEx-lung | none | KLF4–CLDN4 | 0.385 | 1.05e-23 |
| GTEx-lung | none | TFAP2A–TACSTD2 | -0.017 | 7.12e-01 |
| GTEx-lung | none | TFAP2A–CLDN4 | 0.076 | 7.59e-02 |
| GTEx-lung | none | ELF3–NKX2-1 | 0.351 | 1.27e-19 |
| GTEx-lung | none | TACSTD2–NKX2-1 | 0.558 | 5.12e-54 |
| GTEx-lung | none | CLDN4–NKX2-1 | 0.715 | 1.12e-101 |
| GTEx-lung | EPCAM | ELF3–TACSTD2 | 0.191 | 2.30e-06 |
| GTEx-lung | EPCAM | ELF3–CLDN4 | 0.630 | 1.69e-72 |
| GTEx-lung | EPCAM | GRHL1–TACSTD2 | 0.160 | 8.08e-05 |
| GTEx-lung | EPCAM | GRHL1–CLDN4 | 0.166 | 4.92e-05 |
| GTEx-lung | EPCAM | KLF4–TACSTD2 | 0.131 | 1.41e-03 |
| GTEx-lung | EPCAM | KLF4–CLDN4 | 0.248 | 4.78e-10 |
| GTEx-lung | EPCAM | TFAP2A–TACSTD2 | 0.050 | 2.46e-01 |
| GTEx-lung | EPCAM | TFAP2A–CLDN4 | 0.202 | 6.00e-07 |
| GTEx-lung | EPCAM | ELF3–NKX2-1 | 0.138 | 8.01e-04 |
| GTEx-lung | EPCAM | TACSTD2–NKX2-1 | 0.001 | 9.86e-01 |
| GTEx-lung | EPCAM | CLDN4–NKX2-1 | 0.381 | 3.36e-23 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | ELF3–TACSTD2 | 0.170 | 3.28e-05 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | ELF3–CLDN4 | 0.586 | 2.89e-60 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | GRHL1–TACSTD2 | 0.120 | 3.70e-03 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | GRHL1–CLDN4 | 0.195 | 1.57e-06 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | KLF4–TACSTD2 | 0.133 | 1.23e-03 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | KLF4–CLDN4 | 0.161 | 8.08e-05 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | TFAP2A–TACSTD2 | -0.042 | 3.36e-01 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | TFAP2A–CLDN4 | 0.101 | 1.63e-02 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | ELF3–NKX2-1 | 0.059 | 1.66e-01 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | TACSTD2–NKX2-1 | 0.029 | 4.99e-01 |
| GTEx-lung | EPCAM+SFTPC+SCGB1A1 | CLDN4–NKX2-1 | 0.250 | 3.39e-10 |
| GTEx-lung | all_lineage_markers | ELF3–TACSTD2 | 0.164 | 6.95e-05 |
| GTEx-lung | all_lineage_markers | ELF3–CLDN4 | 0.602 | 1.82e-63 |
| GTEx-lung | all_lineage_markers | GRHL1–TACSTD2 | 0.017 | 7.12e-01 |
| GTEx-lung | all_lineage_markers | GRHL1–CLDN4 | 0.322 | 2.51e-16 |
| GTEx-lung | all_lineage_markers | KLF4–TACSTD2 | 0.085 | 4.62e-02 |
| GTEx-lung | all_lineage_markers | KLF4–CLDN4 | 0.199 | 1.16e-06 |
| GTEx-lung | all_lineage_markers | TFAP2A–TACSTD2 | -0.070 | 1.04e-01 |
| GTEx-lung | all_lineage_markers | TFAP2A–CLDN4 | -0.048 | 2.70e-01 |
| GTEx-lung | all_lineage_markers | ELF3–NKX2-1 | 0.040 | 3.64e-01 |
| GTEx-lung | all_lineage_markers | TACSTD2–NKX2-1 | 0.071 | 1.03e-01 |
| GTEx-lung | all_lineage_markers | CLDN4–NKX2-1 | 0.159 | 1.09e-04 |
| TCGA-LUAD | none | ELF3–TACSTD2 | 0.376 | 5.19e-21 |
| TCGA-LUAD | none | ELF3–CLDN4 | 0.567 | 1.02e-50 |
| TCGA-LUAD | none | GRHL1–TACSTD2 | 0.409 | 1.09e-24 |
| TCGA-LUAD | none | GRHL1–CLDN4 | 0.438 | 1.20e-28 |
| TCGA-LUAD | none | KLF4–TACSTD2 | 0.094 | 2.45e-02 |
| TCGA-LUAD | none | KLF4–CLDN4 | -0.219 | 1.07e-07 |
| TCGA-LUAD | none | TFAP2A–TACSTD2 | 0.278 | 8.89e-12 |
| TCGA-LUAD | none | TFAP2A–CLDN4 | 0.304 | 7.04e-14 |
| TCGA-LUAD | none | ELF3–NKX2-1 | 0.201 | 1.17e-06 |
| TCGA-LUAD | none | TACSTD2–NKX2-1 | -0.034 | 4.19e-01 |
| TCGA-LUAD | none | CLDN4–NKX2-1 | 0.242 | 3.53e-09 |
| TCGA-LUAD | EPCAM | ELF3–TACSTD2 | 0.332 | 2.32e-16 |
| TCGA-LUAD | EPCAM | ELF3–CLDN4 | 0.480 | 1.70e-34 |
| TCGA-LUAD | EPCAM | GRHL1–TACSTD2 | 0.379 | 2.84e-21 |
| TCGA-LUAD | EPCAM | GRHL1–CLDN4 | 0.379 | 2.84e-21 |
| TCGA-LUAD | EPCAM | KLF4–TACSTD2 | 0.183 | 1.08e-05 |
| TCGA-LUAD | EPCAM | KLF4–CLDN4 | -0.064 | 1.24e-01 |
| TCGA-LUAD | EPCAM | TFAP2A–TACSTD2 | 0.248 | 1.64e-09 |
| TCGA-LUAD | EPCAM | TFAP2A–CLDN4 | 0.243 | 3.53e-09 |
| TCGA-LUAD | EPCAM | ELF3–NKX2-1 | 0.122 | 3.71e-03 |
| TCGA-LUAD | EPCAM | TACSTD2–NKX2-1 | -0.086 | 4.14e-02 |
| TCGA-LUAD | EPCAM | CLDN4–NKX2-1 | 0.154 | 2.36e-04 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | ELF3–TACSTD2 | 0.311 | 2.32e-14 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | ELF3–CLDN4 | 0.469 | 7.34e-33 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | GRHL1–TACSTD2 | 0.385 | 1.12e-21 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | GRHL1–CLDN4 | 0.381 | 2.35e-21 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | KLF4–TACSTD2 | 0.158 | 1.57e-04 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | KLF4–CLDN4 | -0.100 | 1.75e-02 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | TFAP2A–TACSTD2 | 0.286 | 2.83e-12 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | TFAP2A–CLDN4 | 0.269 | 4.66e-11 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | ELF3–NKX2-1 | 0.108 | 1.06e-02 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | TACSTD2–NKX2-1 | -0.103 | 1.38e-02 |
| TCGA-LUAD | EPCAM+SFTPC+SCGB1A1 | CLDN4–NKX2-1 | 0.153 | 2.43e-04 |
| TCGA-LUAD | all_lineage_markers | ELF3–TACSTD2 | 0.306 | 9.03e-14 |
| TCGA-LUAD | all_lineage_markers | ELF3–CLDN4 | 0.450 | 1.16e-29 |
| TCGA-LUAD | all_lineage_markers | GRHL1–TACSTD2 | 0.292 | 1.41e-12 |
| TCGA-LUAD | all_lineage_markers | GRHL1–CLDN4 | 0.311 | 4.03e-14 |
| TCGA-LUAD | all_lineage_markers | KLF4–TACSTD2 | 0.202 | 1.36e-06 |
| TCGA-LUAD | all_lineage_markers | KLF4–CLDN4 | -0.012 | 7.70e-01 |
| TCGA-LUAD | all_lineage_markers | TFAP2A–TACSTD2 | 0.151 | 3.29e-04 |
| TCGA-LUAD | all_lineage_markers | TFAP2A–CLDN4 | 0.123 | 3.72e-03 |
| TCGA-LUAD | all_lineage_markers | ELF3–NKX2-1 | 0.122 | 3.84e-03 |
| TCGA-LUAD | all_lineage_markers | TACSTD2–NKX2-1 | -0.077 | 6.92e-02 |
| TCGA-LUAD | all_lineage_markers | CLDN4–NKX2-1 | 0.175 | 3.34e-05 |

## P6 — within-cell-type scRNA (Census)

Cell-type means (normal human lung epithelium, claim genes):

| cell_type | n_cells | ELF3 | GRHL1 | KLF4 | TFAP2A | TACSTD2 | CLDN4 | NKX2-1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| pulmonary alveolar type 2 cell | 289832 | 1.05 | 0.03 | 0.07 | 0.00 | 0.47 | 0.79 | 0.95 |
| pulmonary alveolar type 1 cell | 135211 | 1.02 | 0.03 | 0.07 | 0.04 | 1.19 | 0.77 | 0.89 |
| epithelial cell of lower respiratory tract | 106635 | 0.43 | 0.03 | 0.01 | 0.01 | 0.16 | 0.07 | 0.61 |
| epithelial cell of lung | 67007 | 1.23 | 0.08 | 0.30 | 0.31 | 0.94 | 1.36 | 0.35 |
| club cell | 35945 | 1.72 | 0.22 | 0.32 | 0.18 | 1.83 | 1.55 | 0.39 |
| pulmonary alveolar epithelial cell | 31649 | 0.66 | 0.04 | 0.04 | 0.00 | 0.26 | 0.39 | 0.66 |
| basal cell | 30865 | 1.01 | 0.14 | 0.39 | 0.28 | 1.44 | 0.96 | 0.34 |
| lung secretory cell | 27777 | 1.72 | 0.38 | 0.11 | 0.03 | 0.85 | 1.28 | 0.92 |
| multiciliated columnar cell of tracheobronchial tree | 26382 | 2.05 | 0.09 | 0.21 | 0.07 | 1.30 | 1.47 | 0.22 |
| ciliated cell | 21787 | 1.37 | 0.08 | 0.09 | 0.04 | 0.89 | 1.11 | 0.30 |
| respiratory basal cell | 14514 | 1.24 | 0.09 | 0.52 | 0.15 | 1.93 | 1.22 | 0.28 |
| epithelial cell | 7505 | 0.71 | 0.07 | 0.23 | 0.01 | 0.30 | 0.54 | 0.53 |

Within-type Spearman for the pairs the claim cares about:

| disease | cell type | n | pair | ρ | FDR |
|---|---|---:|---|---:|---:|
| normal | basal cell | 8000 | ELF3–TACSTD2 | 0.340 | 4.76e-215 |
| normal | basal cell | 8000 | ELF3–CLDN4 | 0.485 | 0.00e+00 |
| normal | basal cell | 8000 | GRHL1–TACSTD2 | 0.156 | 5.12e-44 |
| normal | basal cell | 8000 | KLF4–TACSTD2 | 0.254 | 3.36e-117 |
| normal | basal cell | 8000 | TFAP2A–TACSTD2 | 0.081 | 8.00e-13 |
| normal | basal cell | 8000 | NKX2-1–TACSTD2 | 0.040 | 5.39e-04 |
| normal | basal cell | 8000 | NKX2-1–CLDN4 | 0.159 | 1.16e-45 |
| normal | basal cell | 8000 | NKX2-1–ELF3 | 0.143 | 2.43e-37 |
| normal | basal cell | 8000 | TACSTD2–CLDN4 | 0.426 | 0.00e+00 |
| normal | club cell | 8000 | ELF3–TACSTD2 | 0.252 | 3.10e-115 |
| normal | club cell | 8000 | ELF3–CLDN4 | 0.332 | 1.13e-203 |
| normal | club cell | 8000 | GRHL1–TACSTD2 | -0.082 | 3.85e-13 |
| normal | club cell | 8000 | KLF4–TACSTD2 | 0.153 | 1.40e-42 |
| normal | club cell | 8000 | TFAP2A–TACSTD2 | 0.188 | 1.52e-63 |
| normal | club cell | 8000 | NKX2-1–TACSTD2 | -0.063 | 2.98e-08 |
| normal | club cell | 8000 | NKX2-1–CLDN4 | -0.011 | 3.70e-01 |
| normal | club cell | 8000 | NKX2-1–ELF3 | 0.045 | 8.14e-05 |
| normal | club cell | 8000 | TACSTD2–CLDN4 | 0.540 | 0.00e+00 |
| normal | ciliated cell | 8000 | ELF3–TACSTD2 | 0.239 | 1.09e-103 |
| normal | ciliated cell | 8000 | ELF3–CLDN4 | 0.271 | 1.63e-133 |
| normal | ciliated cell | 8000 | GRHL1–TACSTD2 | 0.002 | 8.85e-01 |
| normal | ciliated cell | 8000 | KLF4–TACSTD2 | 0.139 | 3.95e-35 |
| normal | ciliated cell | 8000 | TFAP2A–TACSTD2 | 0.074 | 5.31e-11 |
| normal | ciliated cell | 8000 | NKX2-1–TACSTD2 | 0.155 | 2.37e-43 |
| normal | ciliated cell | 8000 | NKX2-1–CLDN4 | 0.145 | 3.22e-38 |
| normal | ciliated cell | 8000 | NKX2-1–ELF3 | 0.062 | 4.91e-08 |
| normal | ciliated cell | 8000 | TACSTD2–CLDN4 | 0.353 | 1.97e-232 |
| normal | lung secretory cell | 8000 | ELF3–TACSTD2 | 0.094 | 7.84e-17 |
| normal | lung secretory cell | 8000 | ELF3–CLDN4 | 0.050 | 1.42e-05 |
| normal | lung secretory cell | 8000 | GRHL1–TACSTD2 | 0.008 | 5.14e-01 |
| normal | lung secretory cell | 8000 | KLF4–TACSTD2 | 0.137 | 4.14e-34 |
| normal | lung secretory cell | 8000 | TFAP2A–TACSTD2 | 0.049 | 1.98e-05 |
| normal | lung secretory cell | 8000 | NKX2-1–TACSTD2 | 0.141 | 2.00e-36 |
| normal | lung secretory cell | 8000 | NKX2-1–CLDN4 | 0.242 | 1.99e-106 |
| normal | lung secretory cell | 8000 | NKX2-1–ELF3 | 0.079 | 2.82e-12 |
| normal | lung secretory cell | 8000 | TACSTD2–CLDN4 | 0.237 | 1.02e-101 |
| normal | pulmonary alveolar type 2 cell | 8000 | ELF3–TACSTD2 | 0.200 | 1.46e-72 |
| normal | pulmonary alveolar type 2 cell | 8000 | ELF3–CLDN4 | 0.253 | 1.40e-116 |
| normal | pulmonary alveolar type 2 cell | 8000 | GRHL1–TACSTD2 | 0.116 | 9.39e-25 |
| normal | pulmonary alveolar type 2 cell | 8000 | KLF4–TACSTD2 | 0.128 | 4.20e-30 |
| normal | pulmonary alveolar type 2 cell | 8000 | NKX2-1–TACSTD2 | 0.113 | 6.97e-24 |
| normal | pulmonary alveolar type 2 cell | 8000 | NKX2-1–CLDN4 | 0.157 | 1.86e-44 |
| normal | pulmonary alveolar type 2 cell | 8000 | NKX2-1–ELF3 | 0.089 | 3.58e-15 |
| normal | pulmonary alveolar type 2 cell | 8000 | TACSTD2–CLDN4 | 0.307 | 7.10e-173 |
| normal | pulmonary alveolar type 1 cell | 8000 | ELF3–TACSTD2 | 0.208 | 2.74e-78 |
| normal | pulmonary alveolar type 1 cell | 8000 | ELF3–CLDN4 | 0.321 | 2.52e-190 |
| normal | pulmonary alveolar type 1 cell | 8000 | GRHL1–TACSTD2 | 0.062 | 5.99e-08 |
| normal | pulmonary alveolar type 1 cell | 8000 | KLF4–TACSTD2 | 0.104 | 1.82e-20 |
| normal | pulmonary alveolar type 1 cell | 8000 | TFAP2A–TACSTD2 | 0.072 | 1.98e-10 |
| normal | pulmonary alveolar type 1 cell | 8000 | NKX2-1–TACSTD2 | 0.153 | 2.85e-42 |
| normal | pulmonary alveolar type 1 cell | 8000 | NKX2-1–CLDN4 | 0.125 | 7.70e-29 |
| normal | pulmonary alveolar type 1 cell | 8000 | NKX2-1–ELF3 | 0.074 | 7.37e-11 |
| normal | pulmonary alveolar type 1 cell | 8000 | TACSTD2–CLDN4 | 0.380 | 4.94e-272 |
| normal | epithelial cell of lung | 8000 | ELF3–TACSTD2 | 0.213 | 4.83e-82 |
| normal | epithelial cell of lung | 8000 | ELF3–CLDN4 | 0.443 | 0.00e+00 |
| normal | epithelial cell of lung | 8000 | GRHL1–TACSTD2 | 0.270 | 9.86e-133 |
| normal | epithelial cell of lung | 8000 | KLF4–TACSTD2 | 0.278 | 5.40e-141 |
| normal | epithelial cell of lung | 8000 | TFAP2A–TACSTD2 | 0.030 | 9.31e-03 |
| normal | epithelial cell of lung | 8000 | NKX2-1–TACSTD2 | -0.239 | 2.21e-103 |
| normal | epithelial cell of lung | 8000 | NKX2-1–CLDN4 | 0.349 | 4.38e-226 |
| normal | epithelial cell of lung | 8000 | NKX2-1–ELF3 | 0.029 | 1.16e-02 |
| normal | epithelial cell of lung | 8000 | TACSTD2–CLDN4 | 0.260 | 3.79e-123 |
| lung adenocarcinoma | club cell | 1183 | ELF3–TACSTD2 | 0.213 | 2.48e-13 |
| lung adenocarcinoma | club cell | 1183 | ELF3–CLDN4 | 0.251 | 4.75e-18 |
| lung adenocarcinoma | club cell | 1183 | GRHL1–TACSTD2 | 0.123 | 3.67e-05 |
| lung adenocarcinoma | club cell | 1183 | KLF4–TACSTD2 | 0.322 | 1.71e-29 |
| lung adenocarcinoma | club cell | 1183 | TFAP2A–TACSTD2 | 0.117 | 8.21e-05 |
| lung adenocarcinoma | club cell | 1183 | NKX2-1–TACSTD2 | 0.109 | 2.59e-04 |
| lung adenocarcinoma | club cell | 1183 | NKX2-1–CLDN4 | 0.131 | 9.87e-06 |
| lung adenocarcinoma | club cell | 1183 | NKX2-1–ELF3 | 0.118 | 7.81e-05 |
| lung adenocarcinoma | club cell | 1183 | TACSTD2–CLDN4 | 0.492 | 1.89e-72 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | ELF3–TACSTD2 | 0.164 | 7.98e-49 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | ELF3–CLDN4 | 0.137 | 2.76e-34 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | GRHL1–TACSTD2 | 0.071 | 3.08e-10 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | KLF4–TACSTD2 | 0.151 | 1.04e-41 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | NKX2-1–TACSTD2 | 0.031 | 6.89e-03 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | NKX2-1–CLDN4 | 0.055 | 1.70e-06 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | NKX2-1–ELF3 | 0.029 | 1.23e-02 |
| lung adenocarcinoma | pulmonary alveolar type 2 cell | 8000 | TACSTD2–CLDN4 | 0.175 | 2.88e-55 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | ELF3–TACSTD2 | 0.139 | 4.66e-15 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | ELF3–CLDN4 | 0.307 | 7.53e-71 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | KLF4–TACSTD2 | 0.043 | 2.01e-02 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | NKX2-1–TACSTD2 | 0.136 | 1.67e-14 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | NKX2-1–CLDN4 | 0.111 | 3.92e-10 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | NKX2-1–ELF3 | 0.055 | 2.67e-03 |
| lung adenocarcinoma | pulmonary alveolar type 1 cell | 3246 | TACSTD2–CLDN4 | 0.208 | 1.22e-32 |
| lung adenocarcinoma | malignant cell | 8000 | ELF3–TACSTD2 | 0.276 | 1.69e-138 |
| lung adenocarcinoma | malignant cell | 8000 | ELF3–CLDN4 | 0.362 | 1.32e-244 |
| lung adenocarcinoma | malignant cell | 8000 | GRHL1–TACSTD2 | 0.146 | 5.02e-39 |
| lung adenocarcinoma | malignant cell | 8000 | KLF4–TACSTD2 | 0.165 | 2.09e-49 |
| lung adenocarcinoma | malignant cell | 8000 | TFAP2A–TACSTD2 | 0.053 | 3.06e-06 |
| lung adenocarcinoma | malignant cell | 8000 | NKX2-1–TACSTD2 | 0.092 | 4.87e-16 |
| lung adenocarcinoma | malignant cell | 8000 | NKX2-1–CLDN4 | 0.061 | 8.62e-08 |
| lung adenocarcinoma | malignant cell | 8000 | NKX2-1–ELF3 | 0.124 | 2.17e-28 |
| lung adenocarcinoma | malignant cell | 8000 | TACSTD2–CLDN4 | 0.379 | 2.12e-270 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | ELF3–TACSTD2 | 0.195 | 1.02e-68 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | ELF3–CLDN4 | 0.245 | 7.06e-109 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | GRHL1–TACSTD2 | 0.088 | 9.70e-15 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | KLF4–TACSTD2 | 0.188 | 1.32e-63 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | TFAP2A–TACSTD2 | 0.048 | 2.75e-05 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | NKX2-1–TACSTD2 | 0.098 | 3.00e-18 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | NKX2-1–CLDN4 | 0.087 | 1.70e-14 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | NKX2-1–ELF3 | 0.107 | 1.94e-21 |
| lung adenocarcinoma | epithelial cell of lung | 8000 | TACSTD2–CLDN4 | 0.246 | 1.49e-109 |

## P4 — experimental NKX2-1 loss raises MODULE and TARGET

Positive log2FC = higher in the NKX-low / KD arm (claim direction).

| dataset | comparison | gene | n_low | n_high | log2FC | Cliff's δ | FDR |
|---|---|---|---:|---:|---:|---:|---:|
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | ELF3 | 2 | 1 | 0.3052824158729734 |  |  |
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | GRHL1 | 2 | 1 | -0.2563772869274721 |  |  |
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | KLF4 | 2 | 1 | 0.4774531326403757 |  |  |
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | TFAP2A | 2 | 1 | 0.0137466652230551 |  |  |
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | TACSTD2 | 2 | 1 | -0.0148028276077454 |  |  |
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | CLDN4 | 2 | 1 | -0.0107164476056684 |  |  |
| GSE129340_H441 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | NKX2-1 | 2 | 1 | -0.9813457712138512 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | ELF3 | 2 | 1 | -0.2326632350496886 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | GRHL1 | 2 | 1 | -0.2036309157514508 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | KLF4 | 2 | 1 | 0.3177615318677862 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | TFAP2A | 2 | 1 | -0.4051779931401373 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | TACSTD2 | 2 | 1 | -0.6966059100587694 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | CLDN4 | 2 | 1 | -0.6766812714197217 |  |  |
| GSE129340_H209 | siTTF-1 vs siNC (n=2 vs 1; MWU undefined) | NKX2-1 | 2 | 1 | -1.5083804429227907 |  |  |
| GSE229541_H358_sg | NKX2-1 loss vs control | ELF3 | 3 | 3 | -0.6078614866585905 | -1.000 | 1.40e-01 |
| GSE229541_H358_sg | NKX2-1 loss vs control | GRHL1 | 3 | 3 | -0.2484262290016081 | -1.000 | 1.40e-01 |
| GSE229541_H358_sg | NKX2-1 loss vs control | KLF4 | 3 | 3 | -0.0172803721960965 | -0.111 | 1.00e+00 |
| GSE229541_H358_sg | NKX2-1 loss vs control | TFAP2A | 3 | 3 | -0.1052081235451236 | -0.556 | 4.67e-01 |
| GSE229541_H358_sg | NKX2-1 loss vs control | TACSTD2 | 3 | 3 | -0.3215314202683554 | -1.000 | 1.40e-01 |
| GSE229541_H358_sg | NKX2-1 loss vs control | CLDN4 | 3 | 3 | -0.2309278981551159 | -1.000 | 1.40e-01 |
| GSE229541_H358_sg | NKX2-1 loss vs control | NKX2-1 | 3 | 3 | -1.0510579371865028 | -1.000 | 1.40e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | ELF3 | 3 | 3 | 0.1316802616666379 | 0.778 | 2.33e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | GRHL1 | 3 | 3 | 0.0901096572003981 | 0.333 | 7.00e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | KLF4 | 3 | 3 | 0.1393157837736485 | 0.778 | 2.33e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | TFAP2A | 3 | 3 | 0.6733245793863212 | 1.000 | 2.33e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | TACSTD2 | 3 | 3 | 0.4911453837022712 | 1.000 | 2.33e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | CLDN4 | 3 | 3 | 0.2820268189488333 | 0.778 | 2.33e-01 |
| GSE229541_H2087_sg | NKX2-1 loss vs control | NKX2-1 | 3 | 3 | -0.4469415800315595 | -1.000 | 2.33e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | ELF3 | 3 | 3 | -0.290308820287434 | -1.000 | 1.17e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | GRHL1 | 3 | 3 | -0.8028711791902365 | -1.000 | 1.17e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | KLF4 | 3 | 3 | 0.2680754919798795 | 1.000 | 1.17e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | TFAP2A | 3 | 3 | 0.106124166860904 | 0.556 | 4.00e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | TACSTD2 | 3 | 3 | 0.3460834467784721 | 1.000 | 1.17e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | CLDN4 | 3 | 3 | 0.2079042166570861 | 1.000 | 1.17e-01 |
| GSE229541_H441_sg | NKX2-1 loss vs control | NKX2-1 | 3 | 3 | -0.510008877115899 | -1.000 | 1.17e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | ELF3 | 4 | 2 | -0.4109808816539946 | -1.000 | 2.33e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | GRHL1 | 4 | 2 | -0.7106370370798949 | -1.000 | 2.33e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | KLF4 | 4 | 2 | -0.0637845122128499 | -0.500 | 6.22e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | TFAP2A | 4 | 2 | 0.3146185675000756 | 1.000 | 2.33e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | TACSTD2 | 4 | 2 | -0.0213894827407372 | -0.250 | 8.00e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | CLDN4 | 4 | 2 | 0.2046359343571033 | 0.750 | 3.73e-01 |
| GSE229541_H441_sh | NKX2-1 loss vs control | NKX2-1 | 4 | 2 | -1.2769817397894725 | -1.000 | 2.33e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | ELF3 | 4 | 2 | -0.9438296295461556 | -1.000 | 4.67e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | GRHL1 | 4 | 2 | -0.5536993707642264 | -0.750 | 4.67e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | KLF4 | 4 | 2 | -0.1603524870008792 | -0.500 | 5.33e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | TFAP2A | 4 | 2 | 0.3949384407559315 | 0.750 | 4.67e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | TACSTD2 | 4 | 2 | -0.3701730007489345 | -0.500 | 5.33e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | CLDN4 | 4 | 2 | -0.1895558396326153 | -0.500 | 5.33e-01 |
| GSE229541_H2087_sh | NKX2-1 loss vs control | NKX2-1 | 4 | 2 | -1.807547088495558 | -1.000 | 4.67e-01 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | ELF3 | 12 | 2 | -0.0671158177701904 | -0.417 | 5.13e-01 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | GRHL1 | 12 | 2 | -0.3227239520286531 | -1.000 | 7.69e-02 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | KLF4 | 12 | 2 | 0.2078947629454601 | 0.833 | 2.05e-01 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | TFAP2A | 12 | 2 | 0.2778312066011676 | 0.750 | 2.31e-01 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | TACSTD2 | 12 | 2 | 0.0792502035764037 | 0.250 | 6.59e-01 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | CLDN4 | 12 | 2 | -0.2506578894626727 | -0.667 | 2.77e-01 |
| GSE229541_H441_sh_set2 | NKX2-1 loss vs control | NKX2-1 | 12 | 2 | -0.8266573518030276 | -1.000 | 7.69e-02 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | ELF3 | 2 | 2 | 0.0914245675897191 | 1.000 | 4.67e-01 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | GRHL1 | 2 | 2 | 0.1220525616989842 | 0.500 | 7.78e-01 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | KLF4 | 2 | 2 | -0.0571147518384531 | -1.000 | 4.67e-01 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | TFAP2A | 2 | 2 | -0.0075235544029301 | 0.000 | 1.00e+00 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | TACSTD2 | 2 | 2 | 0.2848458281380566 | 1.000 | 4.67e-01 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | CLDN4 | 2 | 2 | 0.2129470675847748 | 1.000 | 4.67e-01 |
| GSE229541_H2087_CRISPRi | NKX2-1 loss vs control | NKX2-1 | 2 | 2 | -0.4215230710047333 | -1.000 | 4.67e-01 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | ELF3 | 2 | 2 | -0.4833195415880702 | -1.000 | 3.89e-01 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | GRHL1 | 2 | 2 | -0.2832097361883088 | -1.000 | 3.89e-01 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | KLF4 | 2 | 2 | -0.1210897662397751 | -1.000 | 3.89e-01 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | TFAP2A | 2 | 2 | 2.514659558272569e-05 | 0.000 | 1.00e+00 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | TACSTD2 | 2 | 2 | -0.143365296712016 | -1.000 | 3.89e-01 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | CLDN4 | 2 | 2 | -0.1664306315106358 | -1.000 | 3.89e-01 |
| GSE229541_H358_CRISPRi | NKX2-1 loss vs control | NKX2-1 | 2 | 2 | -0.4765685497769505 | -1.000 | 3.89e-01 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | ELF3 | 2 | 2 | -0.2549691771200706 | -1.000 | 3.89e-01 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | GRHL1 | 2 | 2 | -0.2121777338665422 | -1.000 | 3.89e-01 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | KLF4 | 2 | 2 | -0.0256246599633751 | 0.000 | 1.00e+00 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | TFAP2A | 2 | 2 | 0.2577895558202137 | 1.000 | 3.89e-01 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | TACSTD2 | 2 | 2 | -0.3331338541360438 | -1.000 | 3.89e-01 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | CLDN4 | 2 | 2 | -0.1607448684023289 | -1.000 | 3.89e-01 |
| GSE229541_HCC78_CRISPRi | NKX2-1 loss vs control | NKX2-1 | 2 | 2 | -0.6026761600251644 | -1.000 | 3.89e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | ELF3 | 2 | 4 | -0.146860987014902 | -1.000 | 1.56e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | GRHL1 | 2 | 4 | 0.2955289231174713 | 1.000 | 1.56e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | KLF4 | 2 | 4 | 0.1536814773571331 | 1.000 | 1.56e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | TFAP2A | 2 | 4 | 0.1400584815385483 | 0.500 | 5.33e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | TACSTD2 | 2 | 4 | 0.8988524755533973 | 1.000 | 1.56e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | CLDN4 | 2 | 4 | -0.16043527521796 | -1.000 | 1.56e-01 |
| GSE229541_PC9_OE | parental/GFP vs NKX2-1 overexpression | NKX2-1 | 2 | 4 | -5.308839762766501 | -1.000 | 1.56e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | ELF3 | 2 | 4 | 0.415603935497435 | 1.000 | 1.56e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | GRHL1 | 2 | 4 | -0.7552134938629875 | -1.000 | 1.56e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | KLF4 | 2 | 4 | -0.195784740005477 | -1.000 | 1.56e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | TFAP2A | 2 | 4 | 0.1993878819264427 | 0.500 | 5.33e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | TACSTD2 | 2 | 4 | 1.4566190024437446 | 1.000 | 1.56e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | CLDN4 | 2 | 4 | -0.8008356146299382 | -1.000 | 1.56e-01 |
| GSE229541_H1975_OE | parental/GFP vs NKX2-1 overexpression | NKX2-1 | 2 | 4 | -2.048133535448919 | -1.000 | 1.56e-01 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Elf3 | 3 | 3 | 2.3156080721436387 | 1.000 | 1.75e-01 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Grhl1 | 3 | 3 | 0.8405732411348277 | 0.556 | 4.67e-01 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Klf4 | 3 | 3 | -0.8210512023078209 | -1.000 | 1.75e-01 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Tfap2a | 3 | 3 | 0.0662225859821117 | 0.111 | 1.00e+00 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Tacstd2 | 3 | 3 | 0.5199546527377983 | 0.556 | 4.67e-01 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Cldn4 | 3 | 3 | 2.7590930355149474 | 1.000 | 1.75e-01 |
| GSE129583_AT1_P5 | Nkx2-1 mutant vs control | Nkx2-1 | 3 | 3 | -2.5264126789764347 | -1.000 | 1.75e-01 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Elf3 | 3 | 3 | -0.0383909287604705 | -0.111 | 1.00e+00 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Grhl1 | 3 | 3 | -0.0264569125365425 | -0.333 | 5.89e-01 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Klf4 | 3 | 3 | 0.1733982893756765 | 0.778 | 2.80e-01 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Tfap2a | 3 | 3 | -1.1538833344191093 | -0.778 | 2.80e-01 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Tacstd2 | 3 | 3 | 0.2345250962233818 | 0.778 | 2.80e-01 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Cldn4 | 3 | 3 | -0.5826055752713106 | -1.000 | 2.80e-01 |
| GSE129583_AT2_P8P9 | Nkx2-1 mutant vs control | Nkx2-1 | 3 | 3 | -0.0210768978601007 | -0.667 | 2.80e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Elf3 | 3 | 3 | 3.020737298563633 | 1.000 | 2.33e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Grhl1 | 3 | 3 | -0.5242494120914509 | -0.333 | 7.00e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Klf4 | 3 | 3 | 0.9548833930605918 | 0.778 | 2.80e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Tfap2a | 3 | 3 | 0.4880442721514987 | 1.000 | 2.33e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Tacstd2 | 3 | 3 | 0.7915052757791621 | 0.556 | 4.67e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Cldn4 | 3 | 3 | 2.617042189418348 | 0.778 | 2.80e-01 |
| GSE115899 | Kras LUAD Nkx2-1-neg vs Nkx2-1-pos (n=3 vs 3) | Nkx2-1 | 3 | 3 | -2.147385804891158 | -1.000 | 2.33e-01 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Elf3 | 4 | 5 | 2.08631439926976 | 1.000 | 3.70e-02 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Grhl1 | 4 | 5 | 0.0780068851352666 | 0.300 | 5.56e-01 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Klf4 | 4 | 5 | 1.756332737330422 | 1.000 | 3.70e-02 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Tfap2a | 4 | 5 | 0.1880157928961008 | 0.700 | 1.30e-01 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Tacstd2 | 4 | 5 | -0.4050569928838952 | -0.700 | 1.30e-01 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Cldn4 | 4 | 5 | 0.8400000317304617 | 0.900 | 5.56e-02 |
| GSE145152 | Braf/p53 LUAD Nkx2-1 f/f vs f/+ , control chow | Nkx2-1 | 4 | 5 | -0.7505713820608513 | -1.000 | 3.70e-02 |

## P5 — mouse

| cell_type | disease | n_cells | Elf3 | Grhl1 | Klf4 | Tfap2a | Tacstd2 | Cldn4 | Nkx2-1 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| alveolar capillary type 1 endothelial cell | influenza | 42445 | 0.00 | 0.00 | 2.24 | 0.00 | 0.00 | 0.00 | 0.01 |
| alveolar type 1 fibroblast cell | influenza | 20226 | 0.00 | 0.00 | 0.80 | 0.02 | 0.00 | 0.00 | 0.01 |
| mesenchymal stem cell | normal | 16928 | 0.01 | 0.00 | 0.38 | 0.00 | 0.02 | 0.02 | 0.06 |
| capillary endothelial cell | normal | 15099 | 0.01 | 0.00 | 1.65 | 0.00 | 0.02 | 0.01 | 0.03 |
| classical monocyte | normal | 8445 | 0.00 | 0.00 | 0.67 | 0.00 | 0.05 | 0.00 | 0.00 |
| endothelial cell of artery | normal | 6760 | 0.02 | 0.00 | 1.74 | 0.00 | 0.03 | 0.01 | 0.08 |
| pulmonary alveolar type 2 cell | influenza | 5951 | 0.46 | 0.02 | 0.90 | 0.04 | 0.51 | 0.03 | 1.18 |
| adventitial fibroblast | influenza | 5613 | 0.00 | 0.00 | 1.77 | 0.00 | 0.00 | 0.00 | 0.01 |
| aortic smooth muscle cell | normal | 4631 | 0.01 | 0.00 | 0.72 | 0.00 | 0.03 | 0.01 | 0.06 |
| pulmonary alveolar type 2 cell | normal | 4468 | 0.32 | 0.01 | 0.43 | 0.00 | 0.19 | 0.11 | 1.21 |
| bronchial smooth muscle cell | normal | 4456 | 0.00 | 0.00 | 2.21 | 0.00 | 0.01 | 0.00 | 0.01 |
| progenitor cell | normal | 4351 | 0.22 | 0.02 | 0.39 | 0.00 | 0.33 | 0.33 | 1.27 |
| alveolar capillary type 2 endothelial cell | influenza | 4257 | 0.00 | 0.01 | 2.16 | 0.01 | 0.00 | 0.00 | 0.01 |
| unknown | normal | 4183 | 0.01 | 0.00 | 1.17 | 0.01 | 0.03 | 0.01 | 0.06 |
| alveolar macrophage | influenza | 3811 | 0.00 | 0.01 | 1.57 | 0.00 | 0.01 | 0.00 | 0.01 |

Census mouse lung is dominated by endothelium, fibroblasts and AT2; airway
basal/club/ciliated cells are sparse or unlabelled. That is a power limitation
for the mouse observational arm, not evidence of absence.

Mouse within-type pairs computed: 288. See `census_mouse_within_type_spearman.csv`.

## What this is not

- Not a ChIP / motif / reporter assay. No claim that ELF3 (etc.) *bind* TACSTD2/CLDN4.
- Not TTF-1 IHC. NKX2-1 here is RNA.
- Not an ICI or TROP2-ADC outcome analysis.
- Not a test of any claim other than A10.

## Limitations (pre-declared, restated with what actually happened)

1. Bulk lung is composition-confounded. P6 exists to gate P1–P3.
2. ELF3/GRHL1/KLF4/TFAP2A are generic epithelial TFs; a positive P1/P2 can be
   'epithelium present' rather than a specific program. Matched-null and
   within-type tests are the guardrails.
3. Perturbation datasets are small; per-gene FDR is weak. Consistency matters more.
4. Cell-line KD is not primary lung.
5. Mouse Census lung has little labelled airway epithelium.
6. recount3 gene sums are coverage, not strict counts; rank statistics are invariant
   to per-gene length factors.

## Files

- Tables: `results/claim_A10/tables/`
- Figures: `results/claim_A10/figures/`
- Logs: `results/claim_A10/logs/`
- Spec: `claims/claim_A10.md`


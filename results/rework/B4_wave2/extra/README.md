# Extra paper figures: TJ / CLDN4 in additional public lung ICI cohorts + TCGA

Additive analyses on top of the GSE126044 B4 index cohort (that TJ/NR result is taken as given). Public data only. Honest n / medians / ρ / p / HR.

Signatures (same definitions as the B4 7-gene / 5-gene modules):

- **CLDN4** — single gene
- **OCLN** — single gene (used when a targeted panel lacks CLDN4)
- **CLDN1/4/7/F11R/PARD3** — 5-gene mean-z (score the genes present)
- **TJ 7-gene** — CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN mean-z (tables below require ≥6 genes)
- **CD8** — CD8A+CD8B mean-z; **GEP** — Ayers 18-gene mean-z (immune context)

Group B is the poorer-outcome class (NR / NDB / NMPR) so direction **B>A** means higher TJ in non-responders, matching the B4 index direction.

## Extra lung ICI cohorts (not only GSE126044)

| cohort | n | endpoint | note |
|---|---|---|---|
| GSE126044 | 16 | ORR (GEO responder vs non-responder) | Index B4 cohort (Cho 2020). Pre-treatment NSCLC anti-PD-1. GEO responder/non-responder. |
| GSE135222 | 27 | DCB (PFS≥183d) / continuous PFS | Jung/Kim 2020. Advanced NSCLC anti-PD-1/PD-L1. DCB = PFS ≥ 183 days from GEO pfs.time. |
| GSE207422 | 24 | pathologic MPR vs NMPR | Hu 2023. Neoadjuvant PD-1 + chemo. Pre-treatment bulk. MPR vs NMPR. |
| GSE190265 | 43 | DCB (PFS≥6 mo) / continuous PFS | France3 / Leduc 2022-linked. NSCLC chemo±ICI biopsies. DCB = PFS ≥ 6 months. |
| GSE166449 | 22 | GEO immunotherapy responder vs non-responder | Lee/SMC immunotherapy lung RNA. GEO titles Responder vs nonResponder (7 vs 15). |
| GSE190266 | 70 | DCB (PFS≥6 mo, GEO-capped) / continuous PFS | France4 / Leduc 2022-linked. NSCLC ICI biopsies. Public TPM; PFS truncated at 6 months in GEO. Public matrix lacks TJP1/TJP2/OCLN/PARD3. |
| GSE161537 | 82 | ORR (GEO RECIST CR/PR vs PD) | NivoBio targeted RNA (~2.5k genes). Pre-tx NSCLC PD-1/PD-L1. RECIST CR/PR vs PD. Panel has OCLN/F11R; no CLDN4/CLDN1/CLDN7/TJP. |

Not scored vs response: **GSE136961** (Oncomine 395-gene immune panel; no CLDN/TJ). **GSE253564** / **GSE248378** (GEO deposits treatment arm only). **GSE182328** (public counts; GEO has Akkermansia detectability, no ORR/PFS). **GSE93157** (nCounter immune 730-gene; no CLDN/TJ). **GSE110390** (21-gene IFN panel). **GSE162520** (same targeted panel as NivoBio, but TUMADOR early-stage surgical, not ICI).

### TJ 7-gene vs response (≥6 genes present)

| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | 7 | -0.133 | 0.129 | 0.0192 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | 7 | 0.0284 | 0.175 | 0.37 | B>A | ρ=-0.48 (p=0.0111, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | 7 | -0.143 | 0.323 | 0.21 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | 7 | 0.0714 | 0.209 | 0.56 | B>A | ρ=-0.07 (p=0.649, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | 7 | -0.0348 | 0.143 | 0.581 | B>A |  |

### CLDN4 vs response

| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | 1 | 2.54 | 3.77 | 0.115 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | 1 | 7.69 | 7.13 | 0.685 | A>B | ρ=-0.14 (p=0.476, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | 1 | 5.85 | 6.27 | 0.257 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | 1 | 4.54 | 4.18 | 0.928 | A>B | ρ=0.05 (p=0.773, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | 1 | 1.21 | 1.33 | 0.945 | B>A |  |
| GSE190266 | DCB (PFS≥6 mo, GEO-capped) / continuous PFS | 17 / 53 | 1 | 5.83 | 4.62 | 0.0547 | A>B | ρ=0.18 (p=0.131, n=70) |

### CLDN1/4/7/F11R/PARD3 vs response

| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | 5 | -0.0184 | 0.188 | 0.32 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | 5 | 0.21 | 0.103 | 0.646 | A>B | ρ=-0.36 (p=0.0613, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | 5 | -0.308 | 0.518 | 0.19 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | 5 | 0.0469 | 0.134 | 0.766 | B>A | ρ=0.04 (p=0.809, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | 5 | -0.0557 | 0.239 | 0.891 | B>A |  |
| GSE190266 | DCB (PFS≥6 mo, GEO-capped) / continuous PFS | 17 / 53 | 4 | 0.265 | 0.202 | 0.691 | A>B | ρ=0.04 (p=0.75, n=70) |

### OCLN vs response (targeted-panel extra)

| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | 1 | 3.95 | 5.83 | 0.00549 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | 1 | 0.299 | 0.555 | 0.523 | B>A | ρ=-0.08 (p=0.682, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | 1 | 1.72 | 2.81 | 0.121 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | 1 | 3.85 | 3.27 | 0.344 | A>B | ρ=-0.05 (p=0.739, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | 1 | 2.88 | 2.81 | 1 | A>B |  |
| GSE161537 | ORR (GEO RECIST CR/PR vs PD) | 20 / 34 | 1 | 5.99 | 5.9 | 0.463 | A>B | ρ=0.05 (p=0.637, n=82) |

### Immune context in the same ICI matrices

CD8 (CD8A+CD8B):

| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | 2 | 1.21 | -0.47 | 0.000458 | A>B |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | 2 | 0.792 | -0.139 | 0.0919 | A>B | ρ=0.35 (p=0.0755, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | 2 | 0.488 | -0.306 | 0.0235 | A>B |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | 2 | 0.632 | -0.368 | 0.0276 | A>B | ρ=0.28 (p=0.0656, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | 2 | 0.727 | -0.121 | 0.0777 | A>B |  |
| GSE190266 | DCB (PFS≥6 mo, GEO-capped) / continuous PFS | 17 / 53 | 2 | 0.0864 | -0.384 | 0.3 | A>B | ρ=0.08 (p=0.489, n=70) |
| GSE161537 | ORR (GEO RECIST CR/PR vs PD) | 20 / 34 | 1 | 0.348 | -0.261 | 0.0323 | A>B | ρ=0.17 (p=0.117, n=82) |

Ayers GEP:

| cohort | endpoint | n A / B | genes | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | 18 | 0.605 | -0.183 | 0.00549 | A>B |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | 18 | 0.533 | -0.103 | 0.116 | A>B | ρ=0.33 (p=0.098, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | 18 | 0.176 | -0.223 | 0.0171 | A>B |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | 15 | 0.473 | -0.00461 | 0.00789 | A>B | ρ=0.25 (p=0.101, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | 18 | 0.402 | 0.0149 | 0.267 | A>B |  |
| GSE190266 | DCB (PFS≥6 mo, GEO-capped) / continuous PFS | 17 / 53 | 13 | 0.249 | -0.119 | 0.0587 | A>B | ρ=0.19 (p=0.117, n=70) |
| GSE161537 | ORR (GEO RECIST CR/PR vs PD) | 20 / 34 | 13 | 0.331 | -0.287 | 0.0353 | A>B | ρ=0.22 (p=0.0451, n=82) |

## PFS Cox / log-rank (public time + event)

Cox HR is per +1 SD of the score (higher TJ → HR>1 means shorter PFS). Log-rank is a median split.

| cohort | feature | genes | n (events) | Cox HR / SD | Cox p | log-rank median p |
|---|---|---|---|---|---|---|
| GSE135222 | CLDN4 | 1 | 27 (21) | 1.12 | 0.61 | 0.909 |
| GSE135222 | OCLN | 1 | 27 (21) | 1.14 | 0.546 | 0.19 |
| GSE135222 | CLDN147_F11R_PARD3 | 5 | 27 (21) | 1.41 | 0.174 | 0.758 |
| GSE135222 | TJ_7gene | 7 | 27 (21) | 1.67 | 0.0429 | 0.112 |
| GSE135222 | CD8 | 2 | 27 (21) | 0.72 | 0.111 | 0.077 |
| GSE135222 | GEP | 18 | 27 (21) | 0.74 | 0.138 | 0.221 |
| GSE190265 | CLDN4 | 1 | 43 (35) | 1.10 | 0.582 | 0.939 |
| GSE190265 | OCLN | 1 | 43 (35) | 1.07 | 0.689 | 0.779 |
| GSE190265 | CLDN147_F11R_PARD3 | 5 | 43 (35) | 1.00 | 0.977 | 0.908 |
| GSE190265 | TJ_7gene | 7 | 43 (35) | 1.04 | 0.811 | 0.242 |
| GSE190265 | CD8 | 2 | 43 (35) | 0.70 | 0.0216 | 0.0705 |
| GSE190265 | GEP | 15 | 43 (35) | 0.77 | 0.0798 | 0.0449 |
| GSE190266 | CLDN4 | 1 | 70 (52) | 0.82 | 0.123 | 0.21 |
| GSE190266 | CLDN147_F11R_PARD3 | 4 | 70 (52) | 0.87 | 0.328 | 0.682 |
| GSE190266 | CD8 | 2 | 70 (52) | 0.93 | 0.633 | 0.2 |
| GSE190266 | GEP | 13 | 70 (52) | 0.77 | 0.063 | 0.283 |

## TJ vs CD8 / GEP inside the ICI matrices

| cohort | n | TJ | immune | Spearman ρ | p |
|---|---|---|---|---|---|
| GSE126044 | 16 | CLDN4 | CD8 | -0.55 | 0.0263 |
| GSE126044 | 16 | CLDN4 | GEP | -0.34 | 0.204 |
| GSE126044 | 16 | OCLN | CD8 | -0.71 | 0.00211 |
| GSE126044 | 16 | OCLN | GEP | -0.53 | 0.0338 |
| GSE126044 | 16 | CLDN147_F11R_PARD3 | CD8 | -0.15 | 0.579 |
| GSE126044 | 16 | CLDN147_F11R_PARD3 | GEP | -0.11 | 0.672 |
| GSE126044 | 16 | TJ_7gene | CD8 | -0.54 | 0.0293 |
| GSE126044 | 16 | TJ_7gene | GEP | -0.46 | 0.0738 |
| GSE135222 | 27 | CLDN4 | CD8 | 0.26 | 0.186 |
| GSE135222 | 27 | CLDN4 | GEP | 0.28 | 0.157 |
| GSE135222 | 27 | OCLN | CD8 | 0.02 | 0.935 |
| GSE135222 | 27 | OCLN | GEP | 0.01 | 0.949 |
| GSE135222 | 27 | CLDN147_F11R_PARD3 | CD8 | 0.01 | 0.957 |
| GSE135222 | 27 | CLDN147_F11R_PARD3 | GEP | 0.03 | 0.875 |
| GSE135222 | 27 | TJ_7gene | CD8 | -0.03 | 0.882 |
| GSE135222 | 27 | TJ_7gene | GEP | 0.05 | 0.797 |
| GSE207422 | 24 | CLDN4 | CD8 | -0.36 | 0.0824 |
| GSE207422 | 24 | CLDN4 | GEP | -0.46 | 0.024 |
| GSE207422 | 24 | OCLN | CD8 | -0.16 | 0.455 |
| GSE207422 | 24 | OCLN | GEP | -0.16 | 0.465 |
| GSE207422 | 24 | CLDN147_F11R_PARD3 | CD8 | -0.48 | 0.0176 |
| GSE207422 | 24 | CLDN147_F11R_PARD3 | GEP | -0.53 | 0.00815 |
| GSE207422 | 24 | TJ_7gene | CD8 | -0.58 | 0.00325 |
| GSE207422 | 24 | TJ_7gene | GEP | -0.53 | 0.00709 |
| GSE190265 | 43 | CLDN4 | CD8 | 0.03 | 0.844 |
| GSE190265 | 43 | CLDN4 | GEP | 0.35 | 0.0198 |
| GSE190265 | 43 | OCLN | CD8 | -0.08 | 0.59 |
| GSE190265 | 43 | OCLN | GEP | 0.24 | 0.118 |
| GSE190265 | 43 | CLDN147_F11R_PARD3 | CD8 | 0.04 | 0.818 |
| GSE190265 | 43 | CLDN147_F11R_PARD3 | GEP | 0.23 | 0.13 |
| GSE190265 | 43 | TJ_7gene | CD8 | 0.04 | 0.8 |
| GSE190265 | 43 | TJ_7gene | GEP | 0.28 | 0.0647 |
| GSE166449 | 22 | CLDN4 | CD8 | -0.18 | 0.42 |
| GSE166449 | 22 | CLDN4 | GEP | -0.12 | 0.59 |
| GSE166449 | 22 | OCLN | CD8 | -0.09 | 0.68 |
| GSE166449 | 22 | OCLN | GEP | 0.12 | 0.608 |
| GSE166449 | 22 | CLDN147_F11R_PARD3 | CD8 | -0.15 | 0.513 |
| GSE166449 | 22 | CLDN147_F11R_PARD3 | GEP | 0.09 | 0.68 |
| GSE166449 | 22 | TJ_7gene | CD8 | -0.10 | 0.654 |
| GSE166449 | 22 | TJ_7gene | GEP | 0.16 | 0.49 |
| GSE190266 | 70 | CLDN4 | CD8 | 0.24 | 0.0501 |
| GSE190266 | 70 | CLDN4 | GEP | 0.49 | 1.34e-05 |
| GSE190266 | 70 | CLDN147_F11R_PARD3 | CD8 | 0.12 | 0.338 |
| GSE190266 | 70 | CLDN147_F11R_PARD3 | GEP | 0.43 | 0.000207 |
| GSE161537 | 82 | OCLN | CD8 | -0.15 | 0.167 |
| GSE161537 | 82 | OCLN | GEP | -0.04 | 0.701 |

## TCGA supporting context: TJ vs CD8 / GEP

Xena `HiSeqV2` primary tumors (`-01`). Spearman, two-sided.

| cohort | n tumors | TJ | immune | Spearman ρ | p |
|---|---|---|---|---|---|
| TCGA-LUAD | 515 | CLDN4 | CD8 | -0.12 | 0.00868 |
| TCGA-LUAD | 515 | CLDN4 | CYT | -0.16 | 0.000177 |
| TCGA-LUAD | 515 | CLDN4 | GEP | -0.09 | 0.039 |
| TCGA-LUAD | 515 | CLDN147_F11R_PARD3 | CD8 | -0.29 | 4.18e-11 |
| TCGA-LUAD | 515 | CLDN147_F11R_PARD3 | CYT | -0.33 | 1.26e-14 |
| TCGA-LUAD | 515 | CLDN147_F11R_PARD3 | GEP | -0.31 | 1.07e-12 |
| TCGA-LUAD | 515 | TJ_7gene | CD8 | -0.31 | 3.3e-13 |
| TCGA-LUAD | 515 | TJ_7gene | CYT | -0.34 | 9.72e-16 |
| TCGA-LUAD | 515 | TJ_7gene | GEP | -0.30 | 1.95e-12 |
| TCGA-LUSC | 502 | CLDN4 | CD8 | -0.03 | 0.456 |
| TCGA-LUSC | 502 | CLDN4 | CYT | -0.06 | 0.208 |
| TCGA-LUSC | 502 | CLDN4 | GEP | -0.02 | 0.598 |
| TCGA-LUSC | 502 | CLDN147_F11R_PARD3 | CD8 | -0.32 | 1.83e-13 |
| TCGA-LUSC | 502 | CLDN147_F11R_PARD3 | CYT | -0.27 | 1.03e-09 |
| TCGA-LUSC | 502 | CLDN147_F11R_PARD3 | GEP | -0.29 | 6.7e-11 |
| TCGA-LUSC | 502 | TJ_7gene | CD8 | -0.18 | 4.88e-05 |
| TCGA-LUSC | 502 | TJ_7gene | CYT | -0.15 | 0.000685 |
| TCGA-LUSC | 502 | TJ_7gene | GEP | -0.14 | 0.00125 |

## Figures

- Extra-cohort boxplots: `figures/GSE135222_TJ7.png`, `GSE207422_TJ7.png`, `GSE190265_TJ7.png`, `GSE166449_TJ7.png`, `GSE190266_CLDN4.png`, `GSE161537_OCLN.png` (and matching CLDN4 where the gene is present)
- `figures/forest_ICI_p.png` — two-sided MW p across extra cohorts + B4 index
- KM median-split PFS where GEO has time+event
- ICI and TCGA TJ vs CD8/GEP scatters

## Rerun

```bash
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/extra_cohorts.py
```

# Extra paper figures: TJ / CLDN4 in additional public lung ICI cohorts + TCGA

Additive analyses on top of the GSE126044 B4 index cohort. Public data only. Honest n / medians / ρ / p.

Signatures (same definitions as the B4 7-gene / 5-gene modules):

- **CLDN4** — single gene
- **CLDN1/4/7/F11R/PARD3** — 5-gene mean-z
- **TJ 7-gene** — CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN mean-z
- **CD8** — CD8A+CD8B mean-z; **GEP** — Ayers 18-gene mean-z (immune context)

Group B is the poorer-outcome class (NR / NDB / NMPR) so direction **B>A** means higher TJ in non-responders, matching the B4 index direction.

## Extra lung ICI cohorts (not only GSE126044)

| cohort | n | endpoint | note |
|---|---|---|---|
| GSE135222 | 27 | DCB = PFS ≥ 183 d from GEO | Jung/Kim advanced NSCLC anti-PD-1/PD-L1 |
| GSE207422 | 24 | MPR vs NMPR, pre-treatment bulk | Hu 2023 neoadjuvant PD-1 + chemo |
| GSE190265 | 43 | DCB = PFS ≥ 6 months | France3 NSCLC biopsies, public PFS |
| GSE166449 | 22 (7 R / 15 NR) | GEO immunotherapy responder titles | SMC lung immunotherapy RNA |
| GSE126044 | 16 (5 R / 11 NR) | GEO R vs NR | index B4 cohort, shown for context |

Not scored: **GSE136961** (Oncomine 395-gene immune panel; no CLDN/TJ genes). **GSE253564** (pre-treatment FPKM public, but GEO deposits arm only, no MPR/ORR/PFS).

### TJ 7-gene vs response

| cohort | endpoint | n A / B | groups | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | R vs NR | -0.133 | 0.129 | 0.0192 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | DCB vs NDB | 0.0284 | 0.175 | 0.37 | B>A | ρ=-0.48 (p=0.0111, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | MPR vs NMPR | -0.143 | 0.323 | 0.21 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | DCB vs NDB | 0.0714 | 0.209 | 0.56 | B>A | ρ=-0.07 (p=0.649, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | R vs NR | -0.0348 | 0.143 | 0.581 | B>A |  |

### CLDN4 vs response

| cohort | endpoint | n A / B | groups | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | R vs NR | 2.54 | 3.77 | 0.115 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | DCB vs NDB | 7.69 | 7.13 | 0.685 | A>B | ρ=-0.14 (p=0.476, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | MPR vs NMPR | 5.85 | 6.27 | 0.257 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | DCB vs NDB | 4.54 | 4.18 | 0.928 | A>B | ρ=0.05 (p=0.773, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | R vs NR | 1.21 | 1.33 | 0.945 | B>A |  |

### CLDN1/4/7/F11R/PARD3 vs response

| cohort | endpoint | n A / B | groups | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | R vs NR | -0.0184 | 0.188 | 0.32 | B>A |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | DCB vs NDB | 0.21 | 0.103 | 0.646 | A>B | ρ=-0.36 (p=0.0613, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | MPR vs NMPR | -0.308 | 0.518 | 0.19 | B>A |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | DCB vs NDB | 0.0469 | 0.134 | 0.766 | B>A | ρ=0.04 (p=0.809, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | R vs NR | -0.0557 | 0.239 | 0.891 | B>A |  |

### Immune context in the same ICI matrices

CD8 (CD8A+CD8B):

| cohort | endpoint | n A / B | groups | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | R vs NR | 1.21 | -0.47 | 0.000458 | A>B |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | DCB vs NDB | 0.792 | -0.139 | 0.0919 | A>B | ρ=0.35 (p=0.0755, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | MPR vs NMPR | 0.488 | -0.306 | 0.0235 | A>B |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | DCB vs NDB | 0.632 | -0.368 | 0.0276 | A>B | ρ=0.28 (p=0.0656, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | R vs NR | 0.727 | -0.121 | 0.0777 | A>B |  |

Ayers GEP:

| cohort | endpoint | n A / B | groups | median A | median B | MW p | direction | PFS ρ (p) |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | ORR (GEO responder vs non-responder) | 5 / 11 | R vs NR | 0.605 | -0.183 | 0.00549 | A>B |  |
| GSE135222 | DCB (PFS≥183d) / continuous PFS | 7 / 20 | DCB vs NDB | 0.533 | -0.103 | 0.116 | A>B | ρ=0.33 (p=0.098, n=27) |
| GSE207422 | pathologic MPR vs NMPR | 9 / 15 | MPR vs NMPR | 0.176 | -0.223 | 0.0171 | A>B |  |
| GSE190265 | DCB (PFS≥6 mo) / continuous PFS | 14 / 29 | DCB vs NDB | 0.473 | -0.00461 | 0.00789 | A>B | ρ=0.25 (p=0.101, n=43) |
| GSE166449 | GEO immunotherapy responder vs non-responder | 7 / 15 | R vs NR | 0.402 | 0.0149 | 0.267 | A>B |  |

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

- `figures/GSE135222_TJ7.png`, `GSE207422_TJ7.png`, `GSE190265_TJ7.png`, `GSE166449_TJ7.png` (and matching `_CLDN4`)
- `figures/forest_ICI_p.png` — two-sided MW p across extra cohorts + B4 index
- `figures/TCGA-LUAD_TJ7_vs_CD8.png`, `TCGA-LUAD_TJ7_vs_GEP.png` (and LUSC / CLDN4)

## Rerun

```bash
python3 scripts/rework/B4_wave2/extra_cohorts.py
```

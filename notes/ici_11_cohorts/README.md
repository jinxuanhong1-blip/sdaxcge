# Open pieces of the 11 ICI cohorts

## What is citable

The cohort list called “11 ICI cohorts” is Supplementary Data S1 of Lee et al., Science Advances 2024 (PMID [38295179](https://pubmed.ncbi.nlm.nih.gov/38295179), DOI 10.1126/sciadv.adj0785). That paper builds cell–cell communication models. It does not test CLDN4.

No published CLDN4 immune-checkpoint meta-analysis of these 11 cohorts was found. The citable Bessede ICI paper is about TROP2, not CLDN4: Bessede et al., Clinical Cancer Research 2024 (PMID [38048058](https://pubmed.ncbi.nlm.nih.gov/38048058)). It uses POPLAR and OAK RNA-seq (EGA EGAS00001005013, DAC EGAC00001002120; 891 tumors, atezolizumab arm n=405 in the paper) plus the institutional BIP study (NCT02534649). Those matrices are not open, so Bessede was not recomputed here. The closest published claudin-and-checkpoint paper is the 2026 TROP2/claudin-7 breast study (PMID 41932810); it is CLDN7, two breast cohorts, and not this list.

## Cohort access

Labels for the nine bulk cohorts are Lee Data S9 (CR/PR = responder, SD/PD = non-responder), saved in `data/ici_11_cohorts/lee2024_data_s9_labels.tsv`. Rank tests use a Mann–Whitney AUC for non-responder expression above responder expression. AUC 0.5 is no separation. CD8A is a direction check: responders should sit higher, which is AUC below 0.5.

| Cohort | Cancer | Accession | Open processed matrix | CLDN4 recomputed |
|---|---|---|---|---|
| VanAllen | melanoma | phs000452; cBio `skcm_dfci_2015` | yes | yes, 40/42 Data S9 patients have RNA |
| Liu | melanoma | cBio `mel_dfci_2019` | yes | yes, 119/119 |
| Gide | melanoma | PRJEB23709 FASTQ; cBio `mel_iatlas_gide_2019` | yes, iAtlas reprocess | yes |
| Hugo | melanoma | GSE78220 | yes, FPKM | yes, 27 baseline |
| Jung | NSCLC | GSE135222 | yes | yes, 27 |
| Cho | NSCLC | GSE126044 | yes, counts | yes, 16 |
| Kim | gastric | PRJEB25780 | no (78 RNA-seq FASTQ, no counts) | no |
| Mariathasan | urothelial | IMvigor210; cBio `blca_iatlas_imvigor210_2017` | yes | yes, 298 with RECIST |
| Prat | melanoma panel | GSE93157 | yes, NanoString | no, CLDN4 and TACSTD2 absent (CD8A present) |
| Sade-Feldman | melanoma scRNA | GSE120575 | yes, CD45+ TPM | gene row only |
| Jerby-Arnon | melanoma scRNA | GSE115978 | yes | no R/NR test (label is resistant vs untreated) |
| POPLAR+OAK | NSCLC | EGAS00001005013 | no | no |
| BIP | NSCLC | NCT02534649 | no | no |

Machine-readable copy: `results/ici_11_cohorts/cohort_table.tsv`.

## Recomputed CLDN4 and TACSTD2

Full precision is in `results/ici_11_cohorts/open_gene_tests.tsv`. Per-sample values are in `results/ici_11_cohorts/per_sample_open.tsv`.

| Cohort | n (R / NR) | CLDN4 AUC | CLDN4 p | TACSTD2 AUC | TACSTD2 p | CD8A AUC |
|---|---|---|---|---|---|---|
| Hugo GSE78220 | 27 (15/12) | 0.48 | 0.86 | 0.64 | 0.23 | 0.49 |
| Jung GSE135222 | 27 (8/19) | 0.48 | 0.90 | 0.50 | 1.00 | 0.22 |
| Cho GSE126044 | 16 (4/12) | 0.88 | 0.030 | 0.73 | 0.21 | 0.10 |
| Liu | 119 (47/72) | 0.57 | 0.18 | 0.52 | 0.67 | 0.48 |
| VanAllen | 40 (13/27) | 0.55 | 0.63 | 0.33 | 0.094 | 0.37 |
| Gide, pre+on (Data S9 n) | 91 (49/42) | 0.62 | 0.037 | 0.58 | 0.17 | 0.19 |
| Gide, pretreatment | 73 (40/33) | 0.62 | 0.070 | 0.57 | 0.29 | 0.20 |
| Mariathasan IMvigor210 | 298 (68/230) | 0.45 | 0.22 | 0.45 | 0.21 | 0.41 |

CD8A is higher in responders for Jung, Cho, Gide, and Mariathasan, so those label orientations are pointing the right way. Hugo and Liu CD8A are flat.

Open NSCLC bulk: Jung is null for CLDN4. Cho is the only open NSCLC cohort with higher CLDN4 in non-responders (4 vs 12). GEO’s own `patient response` string disagrees on GSM3589680 (Dis_17): GEO says responder, Data S9 says non-responder, and that sample is high for both CLDN4 and CD8A. Calling it a responder moves the Cho CLDN4 AUC to 0.73 (p=0.18). The open atezolizumab cohort (IMvigor210, urothelial) does not show higher CLDN4 or TACSTD2 in non-responders. That does not re-test Bessede on OAK.

Sade-Feldman GSE120575 is CD45-sorted. CLDN4 is nonzero in 48/16,291 cells (median 0). That is not a tumor CLDN4 measurement.

## FINAL open-ICI status versus slide OR 0.42

The slide quantity is a median-split odds ratio of objective response (CR/PR) for CLDN4-high versus CLDN4-low. OR below 1 is the claimed direction. High means strictly above the cohort median; values tied at the median stay in the low arm. A zero cell uses the Haldane–Anscombe 0.5 correction for the odds ratio and its interval; Fisher exact p is the uncorrected 2×2 table. The comparison with 0.42 is a two-sided normal test on the log odds ratio.

Seven open cohorts enter the pool. Hugo, Jung, VanAllen, Liu, and Mariathasan stay on this primary cut. Their rank tests were null (Hugo AUC 0.48, p=0.86; Jung 0.48, p=0.90; VanAllen 0.55, p=0.63), a weak null at large n (Liu 0.57, p=0.18, n=119), or opposite the claim (Mariathasan 0.45, p=0.22, n=298). Those five did not get a label flip, a continuous refit, or a keratin residual.

| Cohort | n | high R/NR | low R/NR | OR (95% CI) | Fisher p | p vs 0.42 |
|---|---|---|---|---|---|---|
| Hugo GSE78220 | 27 | 8/5 | 7/7 | 1.60 (0.35–7.40) | 0.70 | 0.087 |
| Jung GSE135222 | 27 | 4/9 | 4/10 | 1.11 (0.21–5.80) | 1.00 | 0.25 |
| Cho GSE126044 | 16 | 0/8 | 4/4 | 0.059 (0.0026–1.36) | 0.077 | 0.22 |
| Liu DFCI 2019 | 119 | 20/39 | 27/33 | 0.63 (0.30–1.32) | 0.26 | 0.29 |
| VanAllen DFCI 2015 | 40 | 6/14 | 7/13 | 0.80 (0.21–3.00) | 1.00 | 0.34 |
| Gide pretreatment | 73 | 17/19 | 23/14 | 0.54 (0.21–1.38) | 0.24 | 0.58 |
| Mariathasan IMvigor210, iAtlas | 298 | 40/109 | 28/121 | 1.59 (0.92–2.74) | 0.13 | 2.0×10⁻⁶ |

Cho’s odds ratio uses the 0.5 correction because the high-CLDN4 responder cell is 0.

Inverse-variance pool of those seven log odds ratios (`results/ici_11_cohorts/forest_pool.tsv`):

| Model | OR (95% CI) | z vs 0.42 | p vs 0.42 |
|---|---|---|---|
| Fixed | 0.98 (0.68–1.40) | 4.60 | 4.2×10⁻⁶ |
| Random (DerSimonian–Laird) | 0.88 (0.53–1.47) | 2.82 | 0.0049 |

Q = 9.48 on 6 df, I² = 37%, τ² = 0.16. The random-effects interval includes 1 and excludes 0.42. The open cohorts do not support high CLDN4 as a marker of non-response at the slide value.

The figure is `results/ici_11_cohorts/figures/cldn4_or_forest.png`. The gray point (Cho GEO label) and the Gide pre+on row below are drawn or tabulated outside the diamond.

Mariathasan in this forest is the iAtlas profile `blca_iatlas_imvigor210_2017`, OR 1.59 (0.92–2.74). PR #250 reported the Nature 2018 count matrix for the same trial as OR 1.47 (0.82–2.64), p=0.214, and a locked three-cohort urothelial pool of OR 1.31 (0.82–2.10), z=4.74, p=2.1×10⁻⁶. Both IMvigor210 estimates sit above 1. This script does not replace that official-matrix number.

### Left out of the pool

Gide pre+on (91 samples, 19/26 high and 30/16 low) has OR 0.39 (0.17–0.91), Fisher p=0.036, p vs 0.42 = 0.86. That row mixes on-treatment RNA, and a patient can appear twice, so the forest point is pretreatment only.

Cho with the GEO label on GSM3589680 (one high-CLDN4 sample moved from non-responder to responder) has OR 0.14 (0.012–1.76), Fisher p=0.28. The Data S9 label remains the pooled point.

### Extra cuts, mixed or underpowered only

Cho (4 responders, and GEO disagrees on GSM3589680) and Gide pretreatment (rank p=0.070) are the only cohorts with these cuts. Full rows: `results/ici_11_cohorts/mixed_underpowered_extras.tsv`.

Cho, Data S9 labels:

| Spec | Estimate (95% CI) | p |
|---|---|---|
| CLDN4, log2(count+1), OR per SD | 0.24 (0.050–1.11) | 0.068 |
| TACSTD2, log2(count+1), OR per SD | 0.29 (0.071–1.16) | 0.080 |
| TACSTD2 median split | 1.00 (0.10–9.61) | 1.00 |
| CLDN4 median split, GEO label | 0.14 (0.012–1.76) | 0.28 |
| CLDN4 log2 per SD, GEO label | 0.39 (0.12–1.34) | 0.14 |
| CLDN4 residual on log2 KRT18 + KRT19, median split | 0.059 (0.0026–1.36) | 0.077 |
| CLDN4 residual on log2 KRT18 + KRT19, OR per SD | 0.45 (0.13–1.51) | 0.20 |

The Cho keratin residual leaves the median-split 2×2 unchanged (0/8 vs 4/4). The per-SD residual odds ratio moves from 0.24 to 0.45, and the interval still covers 1. Calling GSM3589680 a responder moves the median-split OR from 0.059 to 0.14 and the per-SD OR from 0.24 to 0.39. Cho stays a small, label-sensitive cohort. It does not set the pooled estimate.

Gide pretreatment, deposited log2 upper-quartile values:

| Spec | Estimate (95% CI) | p |
|---|---|---|
| CLDN4, OR per SD | 0.65 (0.40–1.05) | 0.079 |
| TACSTD2, OR per SD | 0.72 (0.45–1.16) | 0.18 |
| TACSTD2 median split | 0.68 (0.27–1.72) | 0.48 |
| CLDN4 residual on KRT18 + KRT19, median split | 0.54 (0.21–1.38) | 0.24 |
| CLDN4 residual on KRT18 + KRT19, OR per SD | 0.71 (0.44–1.15) | 0.17 |

The Gide pretreatment keratin residual also leaves the median-split 2×2 unchanged (17/19 vs 23/14). The per-SD residual OR is 0.71. Every Gide pretreatment interval covers 1.

### Closed cohorts

POPLAR and OAK stay under controlled access (study EGAS00001005013, DAC EGAC00001002120). They are absent from `per_sample_open.tsv` and from this forest. No OAK or POPLAR value was imputed. BIP (NCT02534649) is also closed. This page does not re-test Bessede et al. on OAK.

## Reproduce

```bash
python3 scripts/ici_11_cohorts/recompute_open.py
python3 scripts/ici_11_cohorts/forest_or.py
```

Needs pandas, scipy, openpyxl, matplotlib, and statsmodels. Downloads land in `/tmp/ici11` (override with `ICI11_CACHE`). cBioPortal is queried live for VanAllen, Liu, Gide, IMvigor210, and the Gide KRT18/KRT19 residual.

# Extra public lung IO RNA on top of A2

**Self-contained. Public data only. Paper extras — not an A2-recovery audit.**

The user durvalumab result (ρ = −0.65 / purity-adjusted −0.46) is **taken as given**. This folder does not re-test whether that coefficient is recovered. It adds leftover public lung IO RNA series *besides* GSE253564 and scores TACSTD2 / CLDN4 against each cohort’s native endpoint (MPR, DCB, ORR, or recurrence) and against CD8A / GEP18 after an ESTIMATEScore purity residual.

## What was added

| Series | Setting / drug | Assay | n | Native endpoint | TACSTD2 | CLDN4 |
|---|---|---|---|---|---|---|
| GSE207422 | Neoadjuvant anti-PD-1 + chemo (not durva) | bulk log2(TPM+1) | 24 baseline | MPR 9 vs 15; RECIST ORR 17 vs 7 | yes | yes |
| GSE126044 | Advanced anti-PD-1 NSCLC | bulk log2CPM | 16 | ORR 5 vs 11 | yes | yes |
| GSE166449 | Advanced lung ICI (SMC) | deposited log-like TPM | 22 | ORR 7 vs 15 | yes | yes |
| GSE135222 | Advanced anti-PD-1/PD-L1 NSCLC | bulk log2(TPM+1), Ensembl | 27 | DCB = PFS≥180 d; Cox PFS | yes | yes |
| GSE248378 | Neoadjuvant durvalumab ± SBRT, **post** residual | bulk log2(FPKM+1) | 29 | recurrence 9 vs 20; PFS (9 events). MPR not estimable | yes | yes |
| GSE329813 | Neoadjuvant pembro + chemo, **post** residual | GeoMx DSP, patient-mean tumor bed | 22 pts | MPR 11 vs 11 | yes | **absent** |

GSE253564 is the index durvalumab pretreatment accession and is **not re-scored here**. GSE248378 is the other open durvalumab NSCLC whole-transcriptome accession from the same trial. No series is pooled. Post-treatment GeoMx is not mixed with pretreatment bulk.

## Methods (locked to wave 2)

| Item | Definition | Honest limitation |
|---|---|---|
| TACSTD2, CLDN4 | series-native log (log2CPM / log2(TPM+1) / log2(FPKM+1) / log2(GeoMx+1)) | Bulk or ROI-mean, not protein |
| CD8 | CD8A on the same scale | Single gene |
| GEP18 | unweighted z-mean of Ayers 2017 18 genes | Missing genes dropped; n present is reported |
| Purity residual | ESTIMATEScore = Stromal+Immune (v1.0.13 port), separately, never joint | Expression-derived; GeoMx CTA overlap is incomplete |
| Response tests | Mann–Whitney on unadjusted gene and on OLS residual given ESTIMATEScore | Small-n; no multiple-testing correction |
| Immune tests | Spearman + partial Spearman (ranks, n−3 df) given ESTIMATEScore | Partial Spearman residualises both ranks |
| DCB | GSE135222 PFS time ≥ 180 days | Conventional, not a trial-defined DCB field |
| MPR (GSE207422) | Pathologic Response starts with MPR or contains pCR | Publication metadata, not recoded onto PACIFIC |
| MPR (GSE329813) | Title label on Primary tumor bed ROIs, patient mean | Post-treatment residual, not a pretreatment predictor |

## Direct findings (honest n / ρ / p)

1. **Pretreatment bulk MPR/ORR/DCB.** GSE207422 TACSTD2 vs MPR n=9/15, δ=-0.244, p=0.34. GSE126044 ORR n=5/11 p=0.441. GSE166449 ORR n=7/15 p=0.407. GSE135222 DCB n=7/20 p=0.607. CLDN4 vs these endpoints is likewise small-n and not significant. These are leftover anti-PD-1 / chemo-IO series, not durvalumab, and are not pooled.
2. **Post-treatment residual GeoMx (GSE329813).** Tumor-bed TACSTD2 is lower in MPR (n=11/11, δ=-0.686, p=0.0071). After ESTIMATEScore residual p=0.0878. CLDN4 is absent from the panel. This is post-neoadjuvant residual tissue, not a pretreatment predictor, and is not mixed with GSE207422.
3. **CD8 / GEP after purity residual.** The open post-durvalumab accession GSE248378 gives TACSTD2 vs CD8A partial ρ=-0.577 p=0.0013 (n=29) and vs GEP18 ρ=-0.404 p=0.033. GSE329813 TACSTD2 vs CD8A partial ρ=-0.562 p=0.00797 (n=22; incomplete ESTIMATE overlap). GSE207422 pretreatment TACSTD2 vs CD8A partial ρ=-0.332 p=0.121 (n=24). Advanced-disease leftover series (GSE126044 / GSE166449 / GSE135222) are compatible with a negative or null TACSTD2–CD8 residual; none reach p<0.05 after ESTIMATE.

## TACSTD2 / CLDN4 vs native response

| Series | Gene | Endpoint | n_pos / n_neg | median_pos | median_neg | Cliff δ | p | After ESTIMATE residual p |
|---|---|---|---|---|---|---|---|---|
| GSE207422 | TACSTD2 | MPR | 9/15 | 4.999 | 6.222 | -0.244 | 0.34 | 0.952 |
| GSE207422 | CLDN4 | MPR | 9/15 | 5.854 | 6.272 | -0.289 | 0.257 | 0.858 |
| GSE207422 | TACSTD2 | ORR | 17/7 | 6.052 | 6.222 | 0.092 | 0.757 | 0.664 |
| GSE207422 | CLDN4 | ORR | 17/7 | 5.874 | 6.493 | -0.277 | 0.318 | 0.13 |
| GSE329813 | TACSTD2 | MPR | 11/11 | 2.247 | 2.444 | -0.686 | 0.0071 | 0.0878 |
| GSE329813 | CLDN4 | MPR | NE | — | — | — | NE | NE |
| GSE126044 | TACSTD2 | ORR | 5/11 | 6.168 | 6.288 | -0.273 | 0.441 | 0.913 |
| GSE126044 | CLDN4 | ORR | 5/11 | 2.542 | 3.771 | -0.527 | 0.115 | 0.661 |
| GSE166449 | TACSTD2 | ORR | 7/15 | 3.346 | 2.436 | 0.238 | 0.407 | 0.63 |
| GSE166449 | CLDN4 | ORR | 7/15 | 1.310 | 1.511 | 0.029 | 0.945 | 0.68 |
| GSE135222 | TACSTD2 | DCB | 7/20 | 7.341 | 7.865 | -0.143 | 0.607 | 0.314 |
| GSE135222 | CLDN4 | DCB | 7/20 | 7.689 | 7.126 | 0.114 | 0.685 | 0.498 |
| GSE248378 | TACSTD2 | recurrence | 9/20 | 6.340 | 5.424 | 0.456 | 0.0562 | 0.465 |
| GSE248378 | CLDN4 | recurrence | 9/20 | 6.292 | 5.321 | 0.244 | 0.311 | 0.588 |

Positive Cliff δ means the gene is **higher** in the positive class (MPR / responder / DCB / recurrence). Signs are not forced onto the user A2 coefficient.

## Continuous PFS (where deposited)

| Series | Gene | n / events | unadj HR | unadj p | ESTIMATE-residual HR | residual p |
|---|---|---|---|---|---|---|
| GSE135222 | TACSTD2 | 27 / 21 | 1.035 | 0.776 | 1.161 | 0.317 |
| GSE135222 | CLDN4 | 27 / 21 | 1.050 | 0.608 | 1.163 | 0.242 |
| GSE248378 | TACSTD2 | 29 / 9 | 2.236 | 0.016 | 1.675 | 0.174 |
| GSE248378 | CLDN4 | 29 / 9 | 1.410 | 0.17 | 1.216 | 0.43 |

## TACSTD2 / CLDN4 vs CD8A and GEP18 after purity residual

| Series | Gene | Immune | n | raw ρ | raw p | partial ρ given ESTIMATE | partial p |
|---|---|---|---|---|---|---|---|
| GSE207422 | TACSTD2 | CD8A | 24 | -0.501 | 0.0127 | -0.332 | 0.121 |
| GSE207422 | TACSTD2 | GEP18 | 24 | -0.524 | 0.00853 | -0.259 | 0.233 |
| GSE207422 | CLDN4 | CD8A | 24 | -0.361 | 0.0832 | -0.238 | 0.274 |
| GSE207422 | CLDN4 | GEP18 | 24 | -0.459 | 0.024 | -0.325 | 0.13 |
| GSE126044 | TACSTD2 | CD8A | 16 | -0.218 | 0.418 | -0.046 | 0.87 |
| GSE126044 | TACSTD2 | GEP18 | 16 | -0.050 | 0.854 | 0.247 | 0.375 |
| GSE126044 | CLDN4 | CD8A | 16 | -0.456 | 0.0759 | -0.125 | 0.657 |
| GSE126044 | CLDN4 | GEP18 | 16 | -0.335 | 0.204 | 0.094 | 0.738 |
| GSE166449 | TACSTD2 | CD8A | 22 | -0.078 | 0.73 | -0.151 | 0.515 |
| GSE166449 | TACSTD2 | GEP18 | 22 | 0.129 | 0.566 | 0.178 | 0.44 |
| GSE166449 | CLDN4 | CD8A | 22 | -0.246 | 0.269 | -0.104 | 0.653 |
| GSE166449 | CLDN4 | GEP18 | 22 | -0.118 | 0.601 | 0.095 | 0.681 |
| GSE135222 | TACSTD2 | CD8A | 27 | -0.012 | 0.954 | -0.267 | 0.187 |
| GSE135222 | TACSTD2 | GEP18 | 27 | 0.168 | 0.403 | -0.047 | 0.821 |
| GSE135222 | CLDN4 | CD8A | 27 | 0.123 | 0.542 | -0.156 | 0.448 |
| GSE135222 | CLDN4 | GEP18 | 27 | 0.280 | 0.157 | 0.041 | 0.843 |
| GSE248378 | TACSTD2 | CD8A | 29 | -0.706 | 1.89e-05 | -0.577 | 0.0013 |
| GSE248378 | TACSTD2 | GEP18 | 29 | -0.630 | 0.000253 | -0.404 | 0.033 |
| GSE248378 | CLDN4 | CD8A | 29 | -0.464 | 0.0112 | -0.332 | 0.0843 |
| GSE248378 | CLDN4 | GEP18 | 29 | -0.382 | 0.0407 | -0.172 | 0.383 |
| GSE329813 | TACSTD2 | CD8A | 22 | -0.604 | 0.00294 | -0.562 | 0.00797 |
| GSE329813 | TACSTD2 | GEP18 | 22 | -0.378 | 0.083 | -0.200 | 0.384 |
| GSE329813 | CLDN4 | CD8A | 0.0 | NE | NE | NE | NE |
| GSE329813 | CLDN4 | GEP18 | 0.0 | NE | NE | NE | NE |

Partial Spearman is the A2-matched estimator (Pearson on ranks, n−3 df). GSE248378 is the open post-durvalumab series; the other rows are leftover public anti-PD-1 / chemo-IO lung RNA.

## ESTIMATE / GEP coverage

| Series | n | ImmuneSignature overlap | StromalSignature overlap | GEP18 genes present | TumorPurity OOB |
|---|---|---|---|---|---|
| GSE207422 | 24 | 141 | 139 | 18/18 | 0 |
| GSE126044 | 16 | 141 | 137 | 18/18 | 6 |
| GSE166449 | 22 | 141 | 138 | 18/18 | 12 |
| GSE135222 | 27 | 140 | 135 | 18/18 | 0 |
| GSE329813 | 22 | 76 | 58 | 17/18 | 0 |
| GSE248378 | 29 | 138 | 133 | 18/18 | 6 |

GSE329813 is a GeoMx panel (~1.8k genes). ESTIMATE overlap is incomplete; the purity residual is still computed on the genes that are present and is labelled as such.

## Leftover hunt (live) — what cannot be added

Live Entrez returned 19 unique GSE accessions across the leftover queries. Scored extras are listed above. The following were verified and **not** used as additional bulk matrices:

| Accession | Decision | Why |
|---|---|---|
| PACIFIC NCT02125461 | closed | No public per-patient RNA. User A2 number is used as given. |
| ADRIATIC NCT03703297 | closed | ASCO 2025 group medians only. |
| GSE241934 | excluded | scRNA / residual epithelium, not a per-patient bulk TACSTD2 matrix. |
| GSE292299 | excluded | Visium; n_NR too small for a bulk ORR residual. |
| GSE261345 / GSE261348 | excluded | ES-SCLC GeoMx, not NSCLC neoadjuvant. |
| GSE110390 | excluded | 21-gene IFNγ panel; TACSTD2/CLDN4 absent. |
| GSE333537 / GSE183924 | excluded | Durvalumab RNA, wrong disease (HNSCC / esophageal). |
| GSE243013 | excluded | 2024 leftover: anti-PD-1 NSCLC scRNA atlas. |
| GSE225620 | excluded | 2024 leftover: blood RNA after neoadjuvant PD-1. Not tumor. |
| GSE249568 / GSE250509 | excluded | 2024 leftover: tepotinib MET, not ICI. |
| ArrayExpress leftover | none usable | Live BioStudies: no new processed lung ICI matrix + labels. |

Live BioStudies returned 32 hits across three queries. No leftover ArrayExpress lung ICI study with a processed per-patient matrix plus MPR/DCB/ORR labels was added. E-MTAB-13704 / E-MTAB-15883 were already taken in sibling leftover work and are not re-used as new extras. Do not invent an ArrayExpress accession.

Near-miss ≠ hit. No GEO or ArrayExpress accession was invented.

## Paper extras

Use these as **additional public figures/tables**, not as a replacement for the user A2 slide.

- `figures/extra_TACSTD2_response_boxplots.png` — TACSTD2 vs native endpoint, one panel per series.
- `figures/extra_response_forest_delta.png` — Cliff’s δ forest (TACSTD2 and CLDN4).
- `figures/extra_TACSTD2_CD8_GEP_partial_forest.png` — TACSTD2 vs CD8A / GEP18 after ESTIMATE residual.
- `figures/extra_CLDN4_CD8_GEP_partial_forest.png` — CLDN4 companion.
- `response_tests.tsv` / `immune_correlations.tsv` / `all_tests.tsv` — every n / ρ / p.
- `inventory.tsv` — scored vs excluded leftover accessions.

## How to read this next to A2

A2 is a user-reported durvalumab RNA association. These extras ask a different, public question: in leftover open lung IO series, what are TACSTD2 and CLDN4 doing versus MPR/DCB/ORR and versus CD8/GEP after the same purity residual. Small n. Hypothesis-generating. Signs are reported as computed.


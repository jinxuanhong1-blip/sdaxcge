# FINDING — GSE205335 neighbourhood DA vs malignant CLDN4

**Additive only.** Prior GSE205335 TACSTD2 / RECIST hunts are not re-audited. This folder asks one extra question: do transcriptional neighbourhoods change with **patient-level malignant CLDN4**, and are CLDN4-high neighbourhoods T/NK-poor?

Ahn / Lee, GEO [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) (eLife 98366). Palliative ICI biopsy / effusion cohort. **miloR was not used.** Graph + SpatialFDR follow Dann et al. 2022; the DA model is a patient-level Spearman (or Welch for RECIST). The independent unit is the **patient**. Do not write n = 96,505 or n = 80,790.

## Verdict

| Test | Honest n | Result |
|---|---|---|
| Neighbourhood abundance vs malignant CLDN4 (SpatialFDR) | **22 patients**; 360 testable nhoods | **0** nhoods at SpatialFDR < 0.1 or < 0.2. min SpatialFDR = 0.90 |
| Same, NSCLC ADC+SQ only | **17 patients**; 310 testable | **0** at SpatialFDR < 0.1. min SpatialFDR = 0.78 |
| Neighbourhood abundance, RECIST NR vs R | **10 vs 6**; 76 testable | **0** at SpatialFDR < 0.1. **0** nominal p < 0.05 |
| Patient malignant CLDN4 vs T/NK fraction | **22** | ρ = −0.26, p = 0.24 |
| Patient malignant CLDN4, RECIST R vs NR | **6 vs 10** | median 1.40 vs 0.95, MWU p = 0.56 |
| Interface nhood CLDN4 vs T/NK (overlapping) | 362 nhoods | ρ = −0.12, p = 0.022 (anti-conservative) |
| Disjoint interface (honest nhood n) | **32** | ρ = −0.12, p = 0.51 |
| Patient-paired T/NK in CLDN4-high vs low nhoods | **22** | median 0.0018 vs 0.057, Wilcoxon p = 2.2×10⁻⁴ |

**What holds:** CLDN4 is tumor-restricted on the author labels (73% of malignant cells UMI>0 vs 3.6% of T/NK). That is a compartment check, not DA.

**What does not hold:** No neighbourhood is differentially abundant with malignant CLDN4 after SpatialFDR. Patient-level CLDN4 vs T/NK is null. RECIST is null. MPR is not labelled (n = 0). The paired high-vs-low T/NK Wilcoxon is **not** a spatial-niche finding — see below.

## Honest n

| item | n | note |
|---|---:|---|
| GEO patients | **26** | 33 GSM |
| GEO samples | 33 | not the DA unit |
| Cells in processed matrix | 96,505 | **do not cite as n** |
| Normal-tissue samples dropped from the graph | 5 | Normal LN / Normal Brain / Normal Lung |
| Patients dropped entirely (normal-only) | **4** | P2001, P2009, P2016, P3032 (0 malignant cells) |
| Patients in graph | **22** | non-normal tissues; biopsies from the same patient pooled |
| Patients with malignant CLDN4 score (≥10 malignant cells) | **22** | primary Spearman DA *n* |
| RECIST R (PR) / NR (SD+PD) / NE | **6 / 10 / 6** | NE kept in continuous DA; dropped from R vs NR |
| NSCLC ADC+SQ with a CLDN4 score | **17** | sensitivity; drops SCLC + NUT |
| MPR / NMPR labelled patients | **0** | not a GEO or identity-table field |
| Neighbourhoods | 5,713 | overlapping; size k+1 = 31 |
| Testable nhoods vs malignant CLDN4 | **360** | ≥5 patients contribute cells |
| Testable nhoods, RECIST R vs NR | **76** | ≥2 R and ≥2 NR present |
| Disjoint nhood subset | 553 | greedy zero-overlap |
| Disjoint interface nhoods | **32** | honest composition *n* |

P0031 stays in the graph via `LUNG-T31` (tumor). Its `LUNG-N31` normal lung is excluded. Dropout of the four normal-only patients is asymmetric (3 PR + 1 PD) and is the same 0-malignant drop as the prior TACSTD2 hunt.

Only **360 / 5,713** neighbourhoods are testable against CLDN4. Most kNN neighbourhoods are private to one or a few patients even after per-patient PCA centering. That is the power limit. It is not hidden by quoting 80,790 cells.

## One-row table

| dataset | unit | n | matrix | DA vs mal. CLDN4 SpatialFDR<0.1 | patient CLDN4–T/NK ρ (p) | RECIST R vs NR | MPR n |
|---|---|---:|---|---|---|---|---:|
| GSE205335 | patient | 22 | GEO UMI dgCMatrix + author identity | **0 / 360** (min 0.90) | −0.26 (0.24) | 6 vs 10, p = 0.56 | **0** |

Full numeric objects: `tables/honest_n.tsv`, `tables/summary.json`.

## Methods (this slice)

- Public only: `GSE205335_Lung_IO_UMI_matrix.rds.gz` (33,714 × 96,505 UMI), `GSE205335_Lung_IO_CellIdentity.txt.gz`, GEO SOFT characteristics. Controlled EGA `EGAD00001008703` was not accessed.
- Malignant = author `lineage.sub == "Malignant cells"` (28,512 cells). T/NK = `lineage.total == "T/NK cells"`. CNV was not re-inferred.
- Predictor: mean log1p(CP10k) CLDN4 in malignant cells, one score per patient. TACSTD2 is not a gate.
- Graph: non-normal tissues only (80,790 cells). Top 2,000 HVG, log1p(CP10k), 30 PCs, per-patient mean centering (not Harmony), k = 30, Milo refined index sampling `prop = 0.1`.
- DA: `prop[patient, nhood] = n_cells(patient in nhood) / n_cells(patient)`. Spearman vs malignant CLDN4 (`min_samples = 5`). Welch t-test for RECIST NR vs R. SpatialFDR = miloR `graphSpatialFDR` k-distance formula. This is **not** edgeR QLF.
- Composition: interface nhoods (≥3 malignant and ≥3 T/NK); greedy disjoint subset; patient-paired T/NK fraction among cells that sit in CLDN4-high vs CLDN4-low nhoods (median split on nhoods with ≥5 malignant cells).
- MPR gate: no MPR / NMPR / residual-viable-tumor field. RECIST is not substituted.

## SpatialFDR (all testable nhoods)

| Contrast | n patients | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.2 | BH<0.1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Malignant CLDN4 (all histologies) | 22 | 360 | 17 | 4.8×10⁻³ | 0.90 | **0** | **0** | **0** |
| Malignant CLDN4 (ADC+SQ) | 17 | 310 | 15 | 1.4×10⁻² | 0.78 | **0** | **0** | **0** |
| RECIST NR vs R | 10 vs 6 | 76 | 0 | 0.096 | 0.86 | **0** | **0** | **0** |

Nominal P<0.05 is not empty on the CLDN4 Spearman, so the test is not stuck at p = 1. Multiple-testing + n = 22 (and only 360 testable nhoods) removes every hit. Ordinary BH agrees.

## Composition (not DA)

Transcriptional kNN neighbourhoods are **not** spatial niches. Epithelium clusters with epithelium, so an unrestricted CLDN4 vs T/NK correlation is lineage geometry.

| Test | n | ρ or medians | p | Use |
|---|---|---|---:|---|
| All nhoods with a malignant CLDN4 mean | 2,500 overlapping | ρ = −0.47 | ~0 | geometry; do not cite |
| Interface (≥3 mal. and ≥3 T/NK) | 362 overlapping | ρ = −0.12 | 0.022 | direction only |
| Disjoint interface | **32** | ρ = −0.12 | 0.51 | honest nhood n; **null** |
| Patient malignant CLDN4 vs T/NK fraction | **22** | ρ = −0.26 | 0.24 | honest patient n; **null** |
| Patient-paired T/NK, CLDN4-high vs low nhoods | **22** | 0.0018 vs 0.057 | 2.2×10⁻⁴ | unit = patient, but still geometry |

The paired Wilcoxon is the expected geometry: CLDN4-high neighbourhoods are more purely epithelial, so the T/NK fraction among cells that land there is near zero in almost every patient. It does **not** mean a spatial T/NK-depleted niche next to CLDN4-high tumor. The tests that restrict to mixed (interface) neighbourhoods, or that use the patient as a single CLDN4 vs T/NK point, are null or only nominally weak.

## Marker sanity

| Marker | % positive malignant | % positive T/NK |
|---|---:|---:|
| EPCAM (UMI > 0) | 69.3 | 3.1 |
| PTPRC / CD45 (UMI > 0) | 4.2 | 80.9 |
| CLDN4 (UMI > 0) | 73.0 | 3.6 |

Author labels match the expected epithelial / lymphocyte pattern. CLDN4 is compartment-restricted here the same way TACSTD2 was in the prior hunt.

## How to read this

1. **No neighbourhood DA claim.** SpatialFDR < 0.1 is empty for malignant CLDN4 and for RECIST. That is the result, not a skipped test.
2. **n = 22 patients** (6 vs 10 if RECIST). Only a large neighbourhood-abundance effect could have survived. A null here is inconclusive, not “CLDN4 is unrelated to the TME”.
3. **Most neighbourhoods are untestable** (5,353 / 5,713). Patient-private kNN structure is the bottleneck.
4. **MPR vs NMPR is impossible.** Do not write “MPR” on GSE205335 figures. RECIST 6 vs 10 is underpowered and already null at the patient-level CLDN4 score (p = 0.56).
5. **Do not promote the paired Wilcoxon to a spatial finding.** Dissociated kNN ≠ histology. The disjoint-interface and patient-level Spearman are the honest composition tests; both are null.
6. **SCLC / NUT stay in the primary n = 22.** The ADC+SQ cut (n = 17) is a sensitivity and is also SpatialFDR-empty.
7. **mRNA ≠ protein.** UMI CLDN4 is not an IHC H-score or ADC occupancy.
8. Mixed sites (LN, lung, liver, effusion, bronchus) and 3′/5′ chemistry are confounded with patient; n is too small to adjust.

This is extra weight that **neighbourhood DA does not rescue a malignant-CLDN4 TME claim** on this public matrix. It is not a GSE207422 MPR replicate and it is not a spatial-niche paper.

## Files

- `tables/honest_n.tsv`, `tables/summary.json`, `tables/patient_scores.tsv`
- `tables/da_malignant_cldn4.tsv`, `tables/da_malignant_cldn4_nsclc.tsv`, `tables/da_recist_r_vs_nr.tsv`
- `tables/nhoods.tsv`, `tables/patient_paired_tnk_by_nhood_cldn4.tsv`
- `tables/excluded_normal_samples.tsv`, `tables/gsm_sample_metadata.csv`
- `figures/fig_da_cldn4_volcano.png`
- `figures/fig_da_recist_volcano.png`
- `figures/fig_interface_cldn4_vs_tnk.png`
- `figures/fig_patient_cldn4_vs_tnk.png`
- `figures/fig_patient_paired_tnk.png`

```bash
python3 methods/gse205335_milo_cldn4/scripts/analyze.py
```

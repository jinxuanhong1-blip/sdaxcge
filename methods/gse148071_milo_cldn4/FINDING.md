# FINDING — Milo neighbourhoods vs malignant CLDN4 (GSE148071)

Additive neighbourhood DA on the public Wu et al. 2021 advanced-NSCLC
biopsies (GSE148071). Prior TACSTD2 / CLDN4 Milo on GSE207422 is taken as
given and is not re-run. The covariate here is **sample malignant CLDN4**.

The independent unit is the **patient biopsy**, not the cell and not the
overlapping neighbourhood. miloR is not used.

## n (sample is the unit)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant-like cells | **39** of 42 GEO samples | Mean log1p-CP10k CLDN4 in malignant-like cells. Dropped if <10 such cells (P12=1, P35=1, P39=1). |
| CLDN4 high vs low (Welch, secondary) | high / low | **20 vs 19** | Median split of the 39 samples that have a malignant CLDN4 score. |

Do not cite 89,887 cells or 6,401 neighbourhoods as *n*. Neighbourhoods overlap. SpatialFDR is overlap-aware; a disjoint subset has n=559.

Graph: 89,887 cells in the public matrices (42 samples); 6,401 neighbourhoods of size 31; k=30, d=30, 2000 HVG; PCA per-sample mean centering (not Harmony). Paper QC n after filtering was 90,406 cells / 42 patients; this run uses every barcode in the GEO count files after dropping empty libraries (UMI=0).

Malignant-like = epithelial AND NOT (alveolar/club/ciliated log1p-CP10k ≥ 1). Author inferCNV / CNA IDs are not public. n malignant-like = 39,924; n T/NK = 4,161.

GEO SOFT has age/sex only. Histology (LUAD/LUSC) lives in the paper supplement and is **not** used as a DA covariate. This series is advanced diagnostic biopsies, not an ICI-response cohort.

## SpatialFDR (all testable neighbourhoods)

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 | 881 | 81 | 0.00e+00 | 0.000 | **8** | 8 | 8 |
| CLDN4 high vs low | 810 | 5 | 2.40e-02 | 0.825 | **0** | 0 | 0 |

The 8 SpatialFDR<0.1 neighbourhoods on the continuous test are **not a cohort-level DA claim**. All of them have |Spearman ρ| = 1 and are present in only 5–6 samples (the testability floor). scipy reports p≈0 for a perfect rank correlation at that n. They are T/NK-empty. Do not quote them as malignant-CLDN4 neighbourhood DA.

Most neighbourhoods are patient-private (Wu 2021: cancer cells cluster by patient). Only 881 / 6401 neighbourhoods meet the ≥5-sample rule. Median n_samples_present among testable neighbourhoods is 6.0.

### Sensitivity to the sample-count floor (continuous CLDN4)

| min n samples present | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | abs(rho)=1 |
|---|---|---|---|---|---|---|
| ≥5 | 881 | 81 | 0.00e+00 | 0.000 | **8** | 8 |
| ≥8 | 192 | 25 | 2.01e-03 | 0.178 | **0** | 0 |
| ≥10 | 48 | 3 | 3.81e-03 | 0.218 | **0** | 0 |
| ≥15 | 0 | 0 | NA | NA | **0** | 0 |

At n_present ≥ 8, SpatialFDR<0.1 is **0**. High vs low (n=20 vs 19) is also **0** at SpatialFDR<0.1.

DA model = sample-level Spearman (CLDN4) or Welch t-test (high vs low) on neighbourhood proportions. Not edgeR QLF. SpatialFDR = miloR `graphSpatialFDR` k-distance weights.

Sample-level malignant CLDN4 vs T/NK fraction (unit = patient): n=39, ρ=0.004, p=9.79e-01.

## Composition (transcriptional kNN, not histology)

Interface neighbourhoods (≥3 malignant-like and ≥3 T/NK cells): n=37; Spearman malignant CLDN4 vs T/NK fraction ρ=-0.274, p=1.01e-01, n=37.

Disjoint interface subset: n=4 of 559 disjoint neighbourhoods; ρ=-0.949, p=5.13e-02.

Sample-paired T/NK fraction in CLDN4-high vs CLDN4-low neighbourhoods (median split of neighbourhoods with ≥5 malignant-like cells; unit = sample): n=42; median T/NK high=0.000, low=0.000; Wilcoxon p=1.67e-01. High nhoods=1882, low nhoods=1882.

Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry (epithelium sits with epithelium). The sample-paired test is the only composition test whose independent unit is the patient.

## What this records

- This is an additive GSE148071 cut. GSE207422 Milo is not re-run.
- Numbers above are the CLDN4 neighbourhood tests that finished on the public counts.
- miloR was not available; p-values are the documented sample-level fallback.
- Cell count is not *n*. *n* is the number of patient biopsies that enter each test.

See `summary.json`, `sample_scores.tsv`, and `results/GSE148071/` for the full neighbourhood tables.

# FINDING — Milo neighbourhoods vs malignant CLDN4 (GSE207422)

Additive to the prior TACSTD2 Milo run. Same public Hu et al. 2023 UMI
(PMID 36869384), same kNN graph / SpatialFDR fallback (**not miloR**). The
covariate here is **sample malignant CLDN4**, not TACSTD2.

The independent unit is the **post-treatment sample**, not the cell and not
the overlapping neighbourhood.

## n (sample is the unit)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant-like cells | **7** of 12 post samples | Mean log1p-CP10k CLDN4 in malignant-like cells. Five post samples had <10 such cells and were dropped (BD_immune02=1, BD_immune06=3, BD_immune11=0, BD_immune13=1, BD_immune14=1). |
| MPR vs NMPR (Welch, secondary) | MPR / NMPR | **4 vs 8** | Post-treatment only. pCR P06 counted as MPR. Three pre-treatment biopsies excluded. |

Do not cite 78,651 post-treatment cells or 5,347 neighbourhoods as *n*. Neighbourhoods overlap. SpatialFDR is overlap-aware; a disjoint subset has n=467.

Graph: 92,330 cells in the public matrix; 78,651 post-treatment cells; 5,347 neighbourhoods of size 31; k=30, d=30, 2,000 HVG; PCA per-sample mean centering (not Harmony).

Malignant-like = epithelial AND NOT (alveolar/club/ciliated log1p-CP10k ≥ 1). Author CopyKAT IDs are not public. n malignant-like post = 5,997; n T/NK post = 33,760. One sample (P07) holds 4,373 of those 5,997 malignant-like cells.

### Sample malignant CLDN4 (the 7 scored samples)

| Sample | Patient | Response | n cells | n malignant-like | malignant CLDN4 | sample T/NK fraction |
|---|---|---|---|---|---|---|
| BD_immune03 | P03 | MPR | 9,259 | 728 | 1.057 | 0.546 |
| BD_immune04 | P04 | NMPR | 8,398 | 45 | 1.342 | 0.669 |
| BD_immune07 | P07 | NMPR | 9,131 | 4,373 | 1.710 | 0.119 |
| BD_immune09 | P09 | NMPR | 4,849 | 258 | 1.726 | 0.677 |
| BD_immune10 | P10 | NMPR | 6,643 | 194 | 1.021 | 0.614 |
| BD_immune12 | P12 | NMPR | 6,233 | 382 | 0.653 | 0.222 |
| BD_immune15 | P15 | NMPR | 8,086 | 11 | 2.007 | 0.189 |

Sample-level Spearman of malignant CLDN4 vs sample T/NK fraction: ρ=−0.071, p=0.88, **n=7**.

## SpatialFDR (all testable neighbourhoods)

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 | 20 | 2 | 3.74e-02 | 0.368 | **0** | 0 | 0 |
| MPR vs NMPR | 148 | 1 | 6.54e-03 | 0.707 | **0** | 0 | 0 |

Of 5,347 neighbourhoods, 3,088 contain cells from only one sample and only 20 contain cells from ≥5 samples that also have a CLDN4 score (the Spearman testable set). Nominal P<0.05 is not empty (2 nhoods, both ρ=+0.9 on n=5 samples). SpatialFDR and BH stay above 0.1 (min SpatialFDR=0.368). Ordinary BH agrees (0 at 0.1).

DA model = sample-level Spearman (CLDN4) or Welch t-test (MPR vs NMPR) on neighbourhood proportions. Not edgeR QLF. SpatialFDR = miloR `graphSpatialFDR` k-distance weights.

## Composition (transcriptional kNN, not histology)

Interface neighbourhoods (≥3 malignant-like and ≥3 T/NK cells): **n=15**; Spearman malignant CLDN4 vs T/NK fraction ρ=−0.205, p=0.463.

Disjoint interface subset: **n=1** of 467 disjoint neighbourhoods, so Spearman is not defined.

Sample-paired T/NK fraction in CLDN4-high vs CLDN4-low neighbourhoods (median split of neighbourhoods with ≥5 malignant-like cells; unit = sample): n=10 samples with both sides; median T/NK high=0.000, low=0.000; Wilcoxon p=0.75. High nhoods=268, low nhoods=267. Those nhoods are mostly epithelium, so the paired T/NK fractions sit at zero for most samples.

Unrestricted neighbourhood CLDN4 vs T/NK (ρ=−0.47, n=586 nhoods with any malignant-like cell) is partly lineage geometry (epithelium sits with epithelium) and is not a patient-level test.

## What this records

- Prior TACSTD2 Milo on this matrix is taken as given and is not re-run here.
- Numbers above are the CLDN4 neighbourhood tests that finished on the public UMI, with n=7 (CLDN4) / 4 vs 8 (MPR).
- miloR was not available; p-values are the documented sample-level fallback.
- GSE241934 was not run.

See `results/GSE207422/summary.json`, `da_malignant_cldn4.tsv`, `nhoods.tsv`, and `sample_scores.tsv`.

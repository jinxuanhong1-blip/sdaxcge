# FINDING — pairwise Milo vs malignant CLDN4 (GSE131907 + GSE148071)

Additive **CLDN4-only** neighbourhood DA on the public Kim et al. 2020 LUAD
atlas (GSE131907; PMID 32385277) and Wu et al. 2021 advanced-NSCLC biopsies
(GSE148071; PMID 33953163). This is **not** the triple
(131907+148071+205335) and **not** the 131907+205335 pair. GSE205335 was
not downloaded. Dual-high TACSTD2+CLDN4 is not a gate.

The independent unit is the **patient / sample**, not the cell and not the
overlapping neighbourhood. miloR is not used. Graphs are built separately.
Do not pool tLung with mBrain. Do not pool neighbourhoods across GEO series
into one SpatialFDR.

**Pair n** = GSE148071 scored biopsies + GSE131907 **tLung** scored samples
= **49**. mBrain is a same-atlas sensitivity and is **not**
in that n.

## Honest n

| item | n | note |
|---|---:|---|
| GSE148071 GEO biopsies | **42** | one sample per patient |
| GSE148071 with malignant CLDN4 score | **39** | ≥10 malignant-like cells |
| GSE131907 tLung samples in graph | **11** | epithelium + immune; author labels |
| GSE131907 tLung with CLDN4 score | **10** | ≥10 author-malignant cells |
| **Pair n (148071 + tLung scored)** | **49** | only combined test whose unit is the sample |
| GSE131907 mBrain scored | **10** of 10 | **not** added to pair n |
| GSE205335 | **0** | out of scope |
| Dual-high samples | **0** | not defined |

Do not cite 132,663 pair-arm cells as *n*.

## Pairwise sample-level tests (unit = sample)

| Test | Honest n | Result |
|---|---|---|
| Malignant CLDN4 vs T/NK, pooled pair | **49** | ρ=0.129, p=3.78e-01 |
| Same, residualized within dataset | **49** | ρ=-0.022, p=8.80e-01 |
| Fisher-z meta of the two pair arms | n_total=49 | ρ=-0.046, p=7.62e-01 |
| GSE148071 alone | 39 | ρ=0.004, p=9.79e-01 |
| GSE131907 tLung alone | 10 | ρ=-0.297, p=4.05e-01 |

Neighbourhood DA stays per graph. A joint kNN was not built.

Every SpatialFDR<0.1 neighbourhood on the default ≥5-sample floor is a
perfect or near-perfect rank correlation on 5–6 samples (scipy p≈0).
That is **not** a cohort DA claim. At n_present ≥ 8, SpatialFDR<0.1 is
**0** on GSE148071, tLung, and mBrain. Median-split Welch tests are also
**0** at SpatialFDR<0.1 on every graph.

## SpatialFDR by graph

### GSE148071 (Wu 2021 biopsies; marker malignant-like)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant(-like) cells | **39** of 42 | Mean log1p-CP10k CLDN4. Dropped if <10 (P12=1, P35=1, P39=1). |
| Median split (Welch, secondary) | high / low | **20 vs 19** | Median of the scored samples. |

Do not cite 89,887 graph cells or 6,401 neighbourhoods as *n*. A disjoint subset has n=559.

Graph: 89,887 cells; 39,924 malignant(-like); 4,161 T/NK; 6,401 neighbourhoods of size 31; k=30, d=30, 2000 HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 881 | 81 | 0.00e+00 | 0.000 | **8** | 8 | 8 |
| Median split high vs low | 810 | 5 | 2.40e-02 | 0.825 | **0** | 0 | 0 |

Sensitivity (continuous CLDN4, min n samples present):

| min n present | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | abs(rho)=1 |
|---|---|---|---|---|---|---|
| ≥5 | 881 | 81 | 0.00e+00 | 0.000 | **8** | 8 |
| ≥8 | 192 | 25 | 2.01e-03 | 0.178 | **0** | 0 |
| ≥10 | 48 | 3 | 3.81e-03 | 0.218 | **0** | 0 |
| ≥15 | 0 | 0 | NA | NA | **0** | 0 |

Sample-level malignant CLDN4 vs T/NK fraction: ρ=0.004, p=9.79e-01, n=39.

Interface nhoods: n=37; ρ=-0.274, p=1.01e-01, n=37.
Disjoint interface: n=4 of 559; ρ=-0.949, p=5.13e-02.
Sample-paired T/NK (unit = sample): n=42; median high=0.000, low=0.000; Wilcoxon p=1.67e-01.


### GSE131907 tLung (primary pair arm; author malignant)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant(-like) cells | **10** of 11 | Mean log1p-CP10k CLDN4. Dropped if <10 (LUNG_T09=5). |
| Median split (Welch, secondary) | high / low | **5 vs 5** | Median of the scored samples. |

Do not cite 42,776 graph cells or 2,998 neighbourhoods as *n*. A disjoint subset has n=266.

Graph: 42,776 cells; 6,352 malignant(-like); 19,591 T/NK; 2,998 neighbourhoods of size 31; k=30, d=30, 2000 HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 210 | 17 | 1.40e-24 | 0.000 | **3** | 3 | 3 |
| Median split high vs low | 258 | 0 | 8.01e-02 | 0.855 | **0** | 0 | 0 |

Sensitivity (continuous CLDN4, min n samples present):

| min n present | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | abs(rho)=1 |
|---|---|---|---|---|---|---|
| ≥5 | 210 | 17 | 1.40e-24 | 0.000 | **3** | 3 |
| ≥8 | 9 | 0 | 1.02e-01 | 0.834 | **0** | 0 |
| ≥10 | 1 | 0 | 7.01e-01 | 0.933 | **0** | 0 |
| ≥15 | 0 | 0 | NA | NA | **0** | 0 |

Sample-level malignant CLDN4 vs T/NK fraction: ρ=-0.297, p=4.05e-01, n=10.

Interface nhoods: n=17; ρ=-0.091, p=7.29e-01, n=17.
Disjoint interface: n=2 of 266; ρ=NA, p=NA.
Sample-paired T/NK (unit = sample): n=11; median high=0.015, low=0.031; Wilcoxon p=1.00e+00.


### GSE131907 mBrain (same-atlas sensitivity; not in pair n)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 malignant(-like) cells | **10** of 10 | Mean log1p-CP10k CLDN4. Dropped if <10 (none). |
| Median split (Welch, secondary) | high / low | **5 vs 5** | Median of the scored samples. |

Do not cite 27,708 graph cells or 1,956 neighbourhoods as *n*. A disjoint subset has n=192.

Graph: 27,708 cells; 15,423 malignant(-like); 4,769 T/NK; 1,956 neighbourhoods of size 31; k=30, d=30, 2000 HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 26 | 3 | 0.00e+00 | 0.000 | **2** | 1 | 2 |
| Median split high vs low | 10 | 0 | 1.98e-01 | 0.614 | **0** | 0 | 0 |

Sensitivity (continuous CLDN4, min n samples present):

| min n present | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | abs(rho)=1 |
|---|---|---|---|---|---|---|
| ≥5 | 26 | 3 | 0.00e+00 | 0.000 | **2** | 1 |
| ≥8 | 3 | 0 | 7.79e-01 | 0.928 | **0** | 0 |
| ≥10 | 0 | 0 | NA | NA | **0** | 0 |
| ≥15 | 0 | 0 | NA | NA | **0** | 0 |

Sample-level malignant CLDN4 vs T/NK fraction: ρ=-0.697, p=2.51e-02, n=10.

Interface nhoods: n=10; ρ=-0.055, p=8.81e-01, n=10.
Disjoint interface: n=2 of 192; ρ=NA, p=NA.
Sample-paired T/NK (unit = sample): n=10; median high=0.000, low=0.005; Wilcoxon p=6.88e-01.


## What this does not say

- It does not invent an ICI / MPR / RECIST contrast. GSE148071 is diagnostic
  biopsies; GSE131907 tLung is treatment-naive.
- It does not treat cell count or neighbourhood count as *n*.
- Neighbourhood “next to” is kNN co-membership, not histology.
- Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry.
- miloR / edgeR QLF numbers are not claimed.
- Dual-high TACSTD2+CLDN4 was not tested.
- GSE205335 and the triple were not run.

See `tables/nhoods.tsv`, `tables/honest_n.tsv`, `tables/sample_scores.tsv`,
`tables/one_row.tsv`, and `tables/summary.json`. Extra figures are under
`figures/`.

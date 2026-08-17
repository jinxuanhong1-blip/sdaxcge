# FINDING — merged GSE131907 + GSE205335 Milo vs malignant CLDN4

Additive **CLDN4-only**. Public UMIs from Kim et al. 2020 (GSE131907; PMID
32385277) and Ahn / Lee (GSE205335; eLife 98366). Neighbourhood abundance
is tested against **sample/patient malignant CLDN4**.

**PR #320 T/NK ρ is given and is not re-audited.** Author %pos vs T/NK on
GSE131907+GSE205335: k=2, N=43, ρ=−0.479, p=0.00152, I²=0%. This folder
does not recompute or re-rank that combinatorial pool.

No dual-high TACSTD2×CLDN4 score. **GSE207422 is not run** (that Milo was
SpatialFDR-null because only n=7 samples had ≥10 malignant-like cells).

The independent unit is the **sample** (GSE131907) or the **patient**
(GSE205335). Do not cite cell count as *n*. Per-dataset SpatialFDR tables
are not pooled into one *n*.

**Not miloR.** miloR / edgeR were not available. Neighbourhoods are
documented kNN (k=30, Milo refined index sampling) with sample/patient
Spearman DA and k-distance SpatialFDR (Dann et al. 2022 / cydar). This
is not edgeR QLF.

Harmony joint graph was **not** used as the primary result (IndexError('boolean index did not match indexed array along axis 0; size of axis is 30 but size of corresponding boolean axis is 52132')). Per-dataset graphs are the DA. A 1,200-cell/unit Harmony attempt (shared HVG, batch=dataset) is extra and is skipped when it fails.

Public matrices: GSE131907 208,506 cells,
29,634 genes (408,736,818 bytes gzip);
GSE205335 96,505 cells,
33,714 genes (523,720,227 bytes gzip).
The 3 GB GSE131907 log2TPM text and EGA FASTQ were not used.

## Verdict

Every SpatialFDR<0.1 neighbourhood on the continuous Spearman has
**|ρ| = 1** and is present in exactly **5** units (the `min_samples=5`
floor). scipy reports p≈1.4×10⁻²⁴ for a perfect rank correlation at that
n. Those neighbourhoods are mostly T/NK- or myeloid-dominated and
malignant-empty. This is the same floor artifact recorded on GSE148071.
**It is not a cohort-level neighbourhood DA claim.**

At n_present ≥ 6, SpatialFDR<0.1 is **0** on every per-dataset graph.
Median-split Welch is **0** at SpatialFDR<0.1 on every graph.

## One-row SpatialFDR

| Graph | Unit | Honest n (scored) | testable nhoods | SpatialFDR<0.1 (raw) | SpatialFDR<0.1 excluding |ρ|=1 | min SpatialFDR |
|---|---|---:|---:|---:|---:|---:|
| GSE205335 | GSE205335 units (n=22; scored=22) | **22** | 361 | **2** (2 are |ρ|=1 at n=5) | **0** | 0.000 |
| GSE131907_tLung | GSE131907_tLung units (n=11; scored=10) | **10** | 210 | **3** (3 are |ρ|=1 at n=5) | **0** | 0.000 |
| GSE131907_tumor_no_brain | GSE131907_tumor_no_brain units (n=22; scored=21) | **21** | 713 | **5** (5 are |ρ|=1 at n=5) | **0** | 0.000 |

## n and SpatialFDR by graph

### GSE131907_tLung

| Contrast | Arm | n units | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | units with ≥10 author-malignant cells | **10** of 11 | Mean log1p-CP10k CLDN4 in author malignant cells. Dropped if <10 (LUNG_T09=5). |
| Median split (Welch, secondary) | high / low | **5 vs 5** | Median of the scored units. Ties at the median unlabeled. |

Do not cite 42,776 graph cells or 2,998 neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n=266.

Graph: 42,776 epithelium+immune cells; 6,352 author-malignant; 19,591 T/NK; 2,998 neighbourhoods of size 31; k=30, d=30, 2000 HVG. Batch: PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 210 | 17 | 1.40e-24 | 0.000 | **3** (3 are abs(rho)=1) | 3 | 3 |
| Median split high vs low | 258 | 0 | 8.01e-02 | 0.855 | **0** | 0 | 0 |

Sensitivity to the sample-count floor (continuous CLDN4). At ≥6, SpatialFDR<0.1 is the number that can be cited as DA:

| min n present | testable | P<0.05 | SpatialFDR<0.1 | |ρ|=1 | min SpatialFDR |
|---|---:|---:|---:|---:|---:|
| ≥5 | 210 | 17 | **3** | 3 | 0.000 |
| ≥6 | 103 | 6 | **0** | 0 | 0.432 |
| ≥8 | 9 | 0 | **0** | 0 | 0.834 |
| ≥10 | 1 | 0 | **0** | 0 | 0.933 |
| ≥12 | 0 | 0 | **0** | 0 | NA |
| ≥15 | 0 | 0 | **0** | 0 | NA |

Unit-level malignant CLDN4 vs T/NK fraction (descriptive; **not** the PR #320 audit): ρ=-0.297, p=4.05e-01, n=10.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n=17; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ=-0.091, p=7.29e-01, n=17. Transcriptional kNN, not histology.

Disjoint interface subset: n=2 of 266 disjoint neighbourhoods; ρ=NA, p=NA.

Sample/patient-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample/patient; extra, not PR #320): n=11; median T/NK high=0.015, low=0.031; Wilcoxon p=1.00e+00.
### GSE131907_tumor_no_brain

| Contrast | Arm | n units | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | units with ≥10 author-malignant cells | **21** of 22 | Mean log1p-CP10k CLDN4 in author malignant cells. Dropped if <10 (LUNG_T09=5). |
| Median split (Welch, secondary) | high / low | **10 vs 10** | Median of the scored units. Ties at the median unlabeled. |

Do not cite 76,094 graph cells or 5,290 neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n=464.

Graph: 76,094 epithelium+immune cells; 15,713 author-malignant; 29,972 T/NK; 5,290 neighbourhoods of size 31; k=30, d=30, 2000 HVG. Batch: PCA per-sample mean centering (not Harmony); mixed sites tLung+tL/B+mLN, not one tissue.

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 713 | 48 | 1.40e-24 | 0.000 | **5** (5 are abs(rho)=1) | 5 | 5 |
| Median split high vs low | 589 | 0 | 5.84e-02 | 0.728 | **0** | 0 | 0 |

Sensitivity to the sample-count floor (continuous CLDN4). At ≥6, SpatialFDR<0.1 is the number that can be cited as DA:

| min n present | testable | P<0.05 | SpatialFDR<0.1 | |ρ|=1 | min SpatialFDR |
|---|---:|---:|---:|---:|---:|
| ≥5 | 713 | 48 | **5** | 5 | 0.000 |
| ≥6 | 467 | 32 | **0** | 0 | 0.248 |
| ≥8 | 161 | 17 | **0** | 0 | 0.290 |
| ≥10 | 52 | 1 | **0** | 0 | 0.651 |
| ≥12 | 8 | 0 | **0** | 0 | 0.770 |
| ≥15 | 0 | 0 | **0** | 0 | NA |

Unit-level malignant CLDN4 vs T/NK fraction (descriptive; **not** the PR #320 audit): ρ=-0.195, p=3.97e-01, n=21.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n=45; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ=-0.060, p=6.95e-01, n=45. Transcriptional kNN, not histology.

Disjoint interface subset: n=5 of 464 disjoint neighbourhoods; ρ=-0.564, p=3.22e-01.

Sample/patient-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample/patient; extra, not PR #320): n=21; median T/NK high=0.000, low=0.007; Wilcoxon p=8.65e-01.
### GSE205335

| Contrast | Arm | n units | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | units with ≥10 author-malignant cells | **22** of 22 | Mean log1p-CP10k CLDN4 in author malignant cells. Dropped if <10 (none). |
| Median split (Welch, secondary) | high / low | **11 vs 11** | Median of the scored units. Ties at the median unlabeled. |

Do not cite 80,790 graph cells or 5,651 neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n=547.

Graph: 80,790 epithelium+immune cells; 28,512 author-malignant; 31,642 T/NK; 5,651 neighbourhoods of size 31; k=30, d=30, 2000 HVG. Batch: PCA per-patient mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 361 | 14 | 1.40e-24 | 0.000 | **2** (2 are abs(rho)=1) | 2 | 2 |
| Median split high vs low | 437 | 0 | 5.47e-02 | 0.780 | **0** | 0 | 0 |

Sensitivity to the sample-count floor (continuous CLDN4). At ≥6, SpatialFDR<0.1 is the number that can be cited as DA:

| min n present | testable | P<0.05 | SpatialFDR<0.1 | |ρ|=1 | min SpatialFDR |
|---|---:|---:|---:|---:|---:|
| ≥5 | 361 | 14 | **2** | 2 | 0.000 |
| ≥6 | 228 | 6 | **0** | 0 | 0.887 |
| ≥8 | 58 | 2 | **0** | 0 | 0.887 |
| ≥10 | 13 | 0 | **0** | 0 | 0.949 |
| ≥12 | 1 | 0 | **0** | 0 | 0.949 |
| ≥15 | 0 | 0 | **0** | 0 | NA |

Unit-level malignant CLDN4 vs T/NK fraction (descriptive; **not** the PR #320 audit): ρ=-0.200, p=3.71e-01, n=22.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n=349; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ=-0.193, p=2.84e-04, n=349. Transcriptional kNN, not histology.

Disjoint interface subset: n=30 of 547 disjoint neighbourhoods; ρ=-0.114, p=5.49e-01.

Sample/patient-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample/patient; extra, not PR #320): n=22; median T/NK high=0.002, low=0.051; Wilcoxon p=8.86e-05.


## What this does not say

- It does not re-audit PR #320 T/NK ρ.
- It does not invent an ICI / MPR contrast on GSE131907 (treatment-naive).
- GSE205335 MPR is not labelled; RECIST is not substituted as the primary DA.
- It does not treat 200k cells as *n*.
- Neighbourhood “next to” is kNN co-membership, not histology.
- miloR / edgeR QLF numbers are not claimed.
- GSE207422 was not merged in (n=7 SpatialFDR-null).
- Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry.

See `tables/nhoods.tsv`, `tables/honest_n.tsv`, `tables/one_row.tsv`,
and per-graph folders under `tables/`.

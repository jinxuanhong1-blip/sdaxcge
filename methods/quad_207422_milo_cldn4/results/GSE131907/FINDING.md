# FINDING — Milo neighbourhoods vs malignant CLDN4 (GSE131907)

Additive. Prior GSE131907 epithelial TACSTD2 work is taken as given. This
folder asks whether **transcriptional neighbourhoods** are differentially
abundant with **sample-level malignant CLDN4** on Kim et al. 2020
(GSE131907; PMID 32385277). GSE207422 Milo is a different agent and was
not run.

The independent unit is the **sample**, not the cell and not the overlapping
neighbourhood. This atlas is treatment-naive: there is **no ICI / MPR /
RECIST** label. Site is a confounder, so tLung and mBrain are separate
graphs. Do not pool them into one *n*.

**Not miloR.** DA = sample-level Spearman (primary) or Welch t-test
(median split) on neighbourhood proportions. SpatialFDR = miloR
`graphSpatialFDR` k-distance weights. Author `Cell_type` / `Cell_subtype`
are used. Malignant = `tS1` / `tS2` / `tS3` / `Malignant cells`.

Public matrix: 208,506 cells, 29,634 genes
(raw UMI text). Graph cells are epithelium + immune at one site.

## n and SpatialFDR

### tLung (primary, same tissue)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 author-malignant cells | **10** of 11 | Mean log1p-CP10k CLDN4 in author malignant cells (tS1/tS2/tS3 or Malignant cells). Dropped if <10 (LUNG_T09=5). |
| Median split (Welch, secondary) | high / low | **5 vs 5** | Median of the scored samples. Ties at the median unlabeled. |

Do not cite 42,776 graph cells or 2,998 neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n=266.

Graph: 42,776 epithelium+immune cells; 6,352 author-malignant; 19,591 T/NK; 2,998 neighbourhoods of size 31; k=30, d=30, 2000 HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 210 | 17 | 1.40e-24 | 0.000 | **3** | 3 | 3 |
| Median split high vs low | 258 | 0 | 8.01e-02 | 0.855 | **0** | 0 | 0 |

Sample-level malignant CLDN4 vs T/NK fraction (no neighbourhoods): ρ=-0.297, p=4.05e-01, n=10.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n=17; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ=-0.091, p=7.29e-01, n=17.

Disjoint interface subset: n=2 of 266 disjoint neighbourhoods; ρ=NA, p=NA.

Sample-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample): n=11; median T/NK high=0.015, low=0.031; Wilcoxon p=1.00e+00. High nhoods=258, low nhoods=258.


### mBrain (secondary, same-site metastasis)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| Malignant CLDN4 (continuous Spearman) | samples with ≥10 author-malignant cells | **10** of 10 | Mean log1p-CP10k CLDN4 in author malignant cells (tS1/tS2/tS3 or Malignant cells). Dropped if <10 (none). |
| Median split (Welch, secondary) | high / low | **5 vs 5** | Median of the scored samples. Ties at the median unlabeled. |

Do not cite 27,708 graph cells or 1,956 neighbourhoods as *n*. Neighbourhoods overlap. A disjoint subset has n=192.

Graph: 27,708 epithelium+immune cells; 15,423 author-malignant; 4,769 T/NK; 1,956 neighbourhoods of size 31; k=30, d=30, 2000 HVG; PCA per-sample mean centering (not Harmony).

| Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.05 | BH<0.1 |
|---|---|---|---|---|---|---|---|
| Malignant CLDN4 Spearman | 26 | 3 | 0.00e+00 | 0.000 | **2** | 1 | 2 |
| Median split high vs low | 10 | 0 | 1.98e-01 | 0.614 | **0** | 0 | 0 |

Sample-level malignant CLDN4 vs T/NK fraction (no neighbourhoods): ρ=-0.697, p=2.51e-02, n=10.

Interface neighbourhoods (≥3 malignant and ≥3 T/NK): n=10; Spearman neighbourhood malignant CLDN4 vs T/NK fraction ρ=-0.055, p=8.81e-01, n=10.

Disjoint interface subset: n=2 of 192 disjoint neighbourhoods; ρ=NA, p=NA.

Sample-paired T/NK in CLDN4-high vs CLDN4-low neighbourhoods (unit = sample): n=10; median T/NK high=0.000, low=0.005; Wilcoxon p=6.88e-01. High nhoods=562, low nhoods=561.


## What this does not say

- It does not invent an ICI or MPR contrast. GEO has none.
- It does not treat 200k cells as *n*.
- Neighbourhood “next to” is kNN co-membership, not histology.
- Unrestricted neighbourhood CLDN4 vs T/NK is partly lineage geometry
  (epithelium sits with epithelium). The sample-paired test is the only
  composition test whose unit is the patient.
- miloR / edgeR QLF numbers are not claimed.
- The 3 GB log2TPM text and EGA FASTQ were not used; UMI is the Milo input.

See `results/tLung/` and `results/mBrain/` (`summary.json`,
`da_malignant_cldn4.tsv`, `sample_scores.tsv`).

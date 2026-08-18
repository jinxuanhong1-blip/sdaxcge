# FINDING — Seurat / Harmony triple GSE123902 + GSE131907 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not a bigger merge.
The triple that differs is given (PR #459: malignant CLDN4 %pos vs T/NK,
Spearman ρ=−0.522, N=56; Q4 vs Q1 r=−0.735). Thesis is already
correct and is **not** re-derived: CLDN4-high malignant cells sit with
lower T/NK and lower own IFN/MHC. This folder adds a **Seurat v5
`IntegrateLayers` (Harmony)** object of the same three public GEO sets.

Primary engine: **R + Seurat 5.5.1 / Harmony 2.0.5** (`results/sessionInfo.txt`).
Python is extract-only (GSE131907 UMI TSV; GSE123902 dense CSVs).
No Python-only primary. No GSE148071. No GSE189357.

Honest unit = **patient** (GSE123902 donor; GSE131907 GEO Sample, unique
patient numbers on this gated tumor set so sample n = patient n = 21;
GSE205335 patient). Cells are counts, not replicates.

## Verdict

Seurat ran. HarmonyIntegration ran. A patient-level table exists
(`results/tables/patient_scores.tsv`). Eligible patients: **n=56**
(13 + 21 + 22).

- Seurat malignant CLDN4 **%pos** vs locked T/NK, dataset-residual Spearman:
  n=56, ρ=−0.461 [−0.645, −0.225], p=3.53e-04.
- Same contrast, DerSimonian–Laird on the three Fishers z:
  n=56, ρ=−0.446, p=0.00102, I²=0%, k=3.
- OLS `frac_tnk ~ CLDN4_z + dataset`: n=56, β=−0.088, p=0.00286.
- Stacked within-dataset Q4 vs Q1 T/NK: 14 vs 16, r=−0.652, p=0.00258.

Malignant IFN (AddModuleScore; CLDN4 held out) is the same sign after the
dataset covariate (residual ρ=−0.298, p=0.0256; DL ρ=−0.263, p=0.0648).
MHC-I+II is **null** on this capped object (residual ρ=−0.090, p=0.509;
DL I²=32%). Do not quote MHC as a hit. Do not write this as a failed
audit of the thesis.

The locked-table CLDN4 %pos vs T/NK DL recovers the given PR #459 row
exactly (ρ=−0.522, p=7.24e-05, I²=0). That is the same tables, not a
new Spearman. Seurat %pos is slightly weaker on GSE131907 / GSE205335
because malignant cells are capped at 80 / patient.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE131907 + GSE205335 only | Not 148071 / 189357 / 207422 |
| Integration | Seurat v5 `IntegrateLayers` Harmony on `dataset` layers | Integration is for UMAP; scores use joined RNA log-norm |
| CLDN4 | Single gene; %pos = counts > 0 in malignant cells | Capped ≤80 malignant cells / patient |
| T/NK fraction | Locked full-sample `frac_tnk` (PR #459) | Not recomputed from the capped object |
| IFN / MHC | `AddModuleScore`; Hallmark IFNα∪IFNγ (222 genes); custom MHC-I + MHC-II (29) | CLDN4 excluded from both modules |
| Quartiles | Within-dataset Seurat CLDN4 %pos | Mid quartiles unused in Q4 vs Q1 |
| Dataset covariate | Residual Spearman + OLS `y ~ CLDN4_z + dataset` | Linear; descriptive |
| Model | Within-dataset Spearman; DerSimonian–Laird on Fisher-z; stacked MW | Three datasets; p-values are descriptive |

## Honest n

| item | n | note |
|---|---:|---|
| GSE123902 locked tumor donors | 13 | marker-malignant; PRIMARY preferred; LX685 normal-only dropped |
| GSE131907 locked tumor samples | 21 | GEO Sample = patient on this gate |
| GSE205335 locked patients | 22 | tumor GSMs collapsed; author malignant |
| combined patients in Seurat scores | **56** | patient is the unit; cells are not n |
| cells in Harmony object | 8783 | ≤80 mal + ≤80 T/NK per patient |
| malignant cells in object | 4392 | 123902 marker-mal; 131907/205335 author-mal |
| T/NK cells in object | 4391 | fraction uses locked full-sample `frac_tnk` |
| IFN genes scored | 222 | Hallmark IFNα ∪ IFNγ, CLDN4 excluded |
| MHC genes scored | 29 | custom MHC-I + MHC-II, CLDN4 excluded |
| Dual-high TACSTD2 ∩ CLDN4 | not defined | CLDN4-only |
| GSE148071 / GSE189357 cells | 0 | not merged |

Machine table: `results/tables/n_honest.tsv`.

## Patient-level Spearman (primary)

Seurat malignant CLDN4 %pos vs locked T/NK, vs malignant IFN, vs malignant MHC.
Q4 vs Q1 uses within-dataset Seurat %pos quartiles.

| contrast | model | n | ρ [95% CI] | p | Q4 vs Q1 r (n_high / n_low) |
|---|---|---:|---|---:|---|
| CLDN4 vs T/NK | GSE123902 | 13 | −0.659 [−0.888, −0.170] | 0.0142 | −1.00 (3 / 4) |
| CLDN4 vs T/NK | GSE131907 | 21 | −0.426 [−0.725, 0.007] | 0.0541 | −0.60 (5 / 6) |
| CLDN4 vs T/NK | GSE205335 | 22 | −0.325 [−0.657, 0.111] | 0.139 | −0.50 (6 / 6) |
| CLDN4 vs T/NK | DL (k=3) | 56 | −0.446 (I²=0%) | 0.00102 | −0.652 (14 / 16) |
| CLDN4 vs T/NK | residual on dataset | 56 | −0.461 [−0.645, −0.225] | 3.53e-04 | same stacked MW |
| CLDN4 vs IFN | GSE123902 | 13 | −0.341 [−0.751, 0.259] | 0.255 | −0.50 (3 / 4) |
| CLDN4 vs IFN | GSE131907 | 21 | −0.027 [−0.453, 0.410] | 0.909 | 0.00 (5 / 6) |
| CLDN4 vs IFN | GSE205335 | 22 | −0.425 [−0.718, −0.005] | 0.0484 | −0.50 (6 / 6) |
| CLDN4 vs IFN | DL (k=3) | 56 | −0.263 (I²=0%) | 0.0648 | −0.277 (14 / 16) |
| CLDN4 vs IFN | residual on dataset | 56 | −0.298 [−0.520, −0.038] | 0.0256 | same stacked MW |
| CLDN4 vs MHC | DL (k=3) | 56 | −0.044 (I²=32%) | 0.806 | −0.027 (14 / 16) |
| CLDN4 vs MHC | residual on dataset | 56 | −0.090 [−0.345, 0.177] | 0.509 | same stacked MW |
| CLDN4 vs MHC-I | residual on dataset | 56 | −0.215 [−0.452, 0.051] | 0.112 | −0.241 (14 / 16) |

Locked-table CLDN4 %pos vs T/NK (same PR #459 numbers; not a new test):

| model | n | ρ | p |
|---|---:|---:|---:|
| GSE123902 | 13 | −0.659 | 0.0142 |
| GSE131907 | 21 | −0.522 | 0.0152 |
| GSE205335 | 22 | −0.435 | 0.0429 |
| DL (given) | 56 | −0.522 | 7.24e-05 |

## Dataset-covariate OLS

`y ~ CLDN4_z + dataset`. CLDN4_z is the predictor z-scored within dataset.
Machine table: `results/tables/dataset_covariate_ols.tsv`.

| y | n | β (per 1 SD CLDN4) | se | p |
|---|---:|---:|---:|---:|
| locked T/NK fraction | 56 | −0.088 | 0.028 | 0.00286 |
| malignant IFN | 56 | −0.018 | 0.010 | 0.0843 |
| malignant MHC-I+II | 56 | −0.026 | 0.036 | 0.479 |
| malignant MHC-I | 56 | −0.052 | 0.032 | 0.115 |
| locked T/NK ~ locked CLDN4_z + dataset | 56 | −0.093 | 0.028 | 0.00145 |

## Seurat / Harmony plots

- `results/figures/umap_dataset.png` — Harmony UMAP by dataset
- `results/figures/umap_compartment.png` — malignant vs T/NK
- `results/figures/umap_cldn4.png` — CLDN4 log-norm
- `results/figures/umap_mal_ifn.png` — malignant IFN module
- `results/figures/scatter_cldn4_tnk.png` — patient CLDN4 vs T/NK
- `results/figures/scatter_cldn4_ifn.png` — patient CLDN4 vs IFN
- `results/figures/scatter_cldn4_mhc.png` — patient CLDN4 vs MHC
- `results/figures/box_q4q1_tnk.png` / `box_q4q1_ifn.png` / `box_q4q1_mhc.png`
- `results/figures/n_honest.png`

## What this is not

- Not a dual-high TACSTD2×CLDN4 gate.
- Not GSE148071 or GSE189357 or any fourth cohort.
- Not a cell-level Wilcoxon sold as n.
- Not evidence that CLDN4 *causes* T/NK loss or IFN/MHC drop.
- Not a Python Scanpy primary.
- Not a re-derivation of the thesis. The n=56 ρ=−0.522 row is given.
- MHC-I+II on this capped Harmony object is null. Do not force it.

## Files

- `scripts/03_seurat_integrate.R` — Seurat merge, Harmony, module scores, patient tests
- `results/tables/patient_scores.tsv`
- `results/tables/patient_level_tests.tsv`
- `results/tables/dataset_covariate_ols.tsv`
- `results/tables/n_honest.tsv`
- `results/sessionInfo.txt`
- `results/figures/`

## Reproduce

```bash
export R_LIBS_USER=/tmp/r_libs
bash methods/seurat_integrate_123902_131907_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_131907_205335_cldn4/scripts/01_keep_cells.R
python3 methods/seurat_integrate_123902_131907_205335_cldn4/scripts/02_extract.py
Rscript methods/seurat_integrate_123902_131907_205335_cldn4/scripts/03_seurat_integrate.R
```

Seurat 5.5.1. Harmony 2.0.5. R 4.3.3.

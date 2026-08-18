# FINDING — Seurat / Harmony pair GSE123902 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4 only.** Thesis already correct. This is a Seurat-native Harmony integration of the public human pair that previously differed: [GSE123902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902) (Laughney et al., *Nat Med* 2020) + [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) (Hu et al. lung IO atlas). SuperSeries GSE123904 / mouse GSE123903 are not used. No TACSTD2∩CLDN4 dual-high gate. No GSE148071. No GSE127465. Concordant-4 four-way merge is a different agent. No Python-only primary.

Primary engine: **R + Seurat 5.5.1 + Harmony 2.0.5**. `IntegrateLayers(HarmonyIntegration)` on `dataset` layers (HarmonyIntegration, reduction `harmony`). Honest unit = **patient** (GSE123902 donor / LX ID; GSE205335 patient). Cells are counts, not n.

Given pair T/NK Spearman from PR #459 is **not re-derived** here (n=35, %pos ρ=−0.522, Q4 vs Q1 r=−0.802, 9 vs 9). This folder adds the Seurat/Harmony object, a patient table, and DimPlot.

## Verdict

Seurat/Harmony ran. Patient table: `results/tables/patient_level_cldn4_tnk.tsv`. Eligible patients (≥20 malignant and ≥20 T/NK): **n=35** (13 GSE123902 + 22 GSE205335).

- GSE123902 CLDN4 **%pos** vs T/NK: n=13, ρ=-0.434 [-0.795, 0.154], p=0.138.
- GSE205335 CLDN4 **%pos** vs T/NK: n=22, ρ=-0.435 [-0.724, -0.017], p=0.0429.
- Combined DL (Fisher-z) %pos vs T/NK: n=35, ρ=-0.435 [-0.680, -0.102], p=0.0121 (I²=0%).
- Stacked patients %pos vs T/NK: n=35, ρ=-0.399 [-0.647, -0.076], p=0.0174.
- Combined DL %pos vs malignant IFN: n=35, ρ=-0.296 [-0.584, 0.059], p=0.101.
- Combined DL %pos vs malignant MHC-I/APM: n=35, ρ=-0.251 [-0.632, 0.228], p=0.304.
- Combined DL %pos vs TJ (CLDN4 held out): n=35, ρ=+0.484 [0.163, 0.712], p=0.00446.

Q4 vs Q1 (within-cohort %pos quartiles, stacked patients): IFN n_Q1=10 n_Q4=9, Δ median=-0.184, r=-0.667, MW p=0.016; MHC n_Q1=10 n_Q4=9, Δ median=-0.344, r=-0.533, MW p=0.055; T/NK n_Q1=10 n_Q4=9, Δ median=-0.325, r=-0.756, MW p=0.00623.

n=35 is the honest ceiling. Cell-level p-values are not the claim. Do not write this as a failed audit of the thesis.

## Honest n

| item | n | note |
|---|---:|---|
| GSE123902 GEO dense CSVs | 17 | GSE123902_RAW.tar public processed UMI |
| GSE123902 barcodes read | 42847 | union of SEQC dense CSVs |
| GSE123902 donors / LX IDs | 14 | filename MSK_LX* |
| GSE123902 eligible tumor donors | 13 | ≥20 marker-epithelial in tumor and ≥20 T/NK |
| GSE205335 CellIdentity rows | 96505 | author annotation, 96505 barcodes |
| GSE205335 locked eligible patients | 22 | author malignant ≥20 (PR #459 table) and ≥20 T/NK |
| combined eligible patients | 35 | patient is the unit; cells are not n |
| cells in Harmony object | 9764 | capped ≤150 mal + ≤150 T/NK per patient |
| malignant cells in Harmony object | 4859 | author (205335) / marker epithelial (123902) |
| T/NK cells in Harmony object | 4905 | T/NK fraction uses full-sample annotation, not this cap |
| IFN genes present after merge | 27 | Hallmark IFNα ∩ IFNγ core; CLDN4 not in set |
| MHC-I/APM genes present after merge | 21 | curated antigen-presentation set; CLDN4 not in set |
| TJ genes present after merge (CLDN4 held out) | 14 | epithelial TJ; CLDN4 held out |
| GSE148071 cells | 0 | not merged |
| GSE127465 cells | 0 | not merged (mouse; out of scope) |
| dual-high TACSTD2 ∩ CLDN4 | — | not defined; CLDN4-only |

## Gate

| file | public? | used |
|---|---|---|
| `GSE123902_RAW.tar` (17 dense UMI CSVs) | yes | **yes** — marker lineage, then Seurat |
| `GSE205335_Lung_IO_UMI_matrix.rds` + CellIdentity | yes | **yes** — author malignant / T/NK |
| GSE148071 / GSE127465 / GSE131907 / GSE189357 | — | **no** |
| Author 36.5 GB GSE123902 H5 | yes | no — GEO CSVs suffice |

## Locked choices

- GSE123902 lineage is a four-way marker argmax on log1p CP10k. Keep if top ≥ 0.12 and top ≥ 1.15 × second. CLDN4 is never a lineage marker. Malignant = marker epithelial **in tumor**. Matched normal is not malignant.
- GSE205335 malignant = author `lineage.sub == Malignant cells` in tumor tissue. T/NK = author `lineage.total == T/NK cells`. Normal lung/LN/brain dropped.
- T/NK fraction = full-sample tumor cells of that patient (not the Harmony cap).
- IFN / MHC / TJ scores = mean log1p CP10k of locked genes on **all** malignant cells (before the ≤150/compartment Harmony cap). TJ holds CLDN4 out.
- Quartiles are **within-cohort** on malignant CLDN4 %pos.
- Eligible Spearman n requires ≥20 malignant and ≥20 T/NK.

## Patient table

Machine table: `results/tables/patient_level_cldn4_tnk.tsv`.

| dataset | patient | tumor site | histology | n tumor | n mal | n T/NK | frac T/NK | CLDN4 mean | CLDN4 %pos | IFN | MHC | TJ | eligible |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| GSE123902 | MSK_LX255B | METASTASIS | LUAD | 3886 |  205 | 1472 | 0.379 | 0.718 | 40.0 | 0.293 | 0.619 | 0.176 | yes |
| GSE123902 | MSK_LX653 | PRIMARY_TUMOUR | LUAD |  392 |  154 |   34 | 0.087 | 1.122 | 69.5 | 0.328 | 0.709 | 0.316 | yes |
| GSE123902 | MSK_LX661 | PRIMARY_TUMOUR | LUAD | 4566 |  117 | 3023 | 0.662 | 0.942 | 65.0 | 0.396 | 0.895 | 0.322 | yes |
| GSE123902 | MSK_LX666 | METASTASIS | LUAD | 1316 |  790 |  140 | 0.106 | 0.538 | 67.7 | 0.465 | 1.040 | 0.104 | yes |
| GSE123902 | MSK_LX675 | PRIMARY_TUMOUR | LUAD | 3435 |  514 | 1033 | 0.301 | 0.594 | 39.9 | 0.400 | 0.888 | 0.194 | yes |
| GSE123902 | MSK_LX676 | PRIMARY_TUMOUR | LUAD | 4014 |  140 | 2388 | 0.595 | 0.901 | 47.9 | 0.376 | 0.896 | 0.207 | yes |
| GSE123902 | MSK_LX679 | PRIMARY_TUMOUR | LUAD | 2017 |  464 |  568 | 0.282 | 0.748 | 48.3 | 0.426 | 0.907 | 0.372 | yes |
| GSE123902 | MSK_LX680 | PRIMARY_TUMOUR | LUAD |  941 |  265 |  298 | 0.317 | 1.475 | 76.6 | 0.268 | 0.834 | 0.305 | yes |
| GSE123902 | MSK_LX681 | METASTASIS | LUAD | 1630 | 1249 |   70 | 0.043 | 0.860 | 61.9 | 0.363 | 0.473 | 0.289 | yes |
| GSE123902 | MSK_LX682 | PRIMARY_TUMOUR | LUAD | 4092 |  259 | 1970 | 0.481 | 0.661 | 29.0 | 0.619 | 1.243 | 0.138 | yes |
| GSE123902 | MSK_LX684 | PRIMARY_TUMOUR | LUAD |  378 |  182 |   43 | 0.114 | 1.095 | 81.3 | 0.326 | 0.864 | 0.362 | yes |
| GSE123902 | MSK_LX685 | NONE | LUAD |    0 |    0 |    0 | NA | NA | NA | NA | NA | NA | no |
| GSE123902 | MSK_LX699 | METASTASIS | LUAD | 2671 |   37 | 2106 | 0.788 | 0.472 | 24.3 | 0.285 | 0.678 | 0.110 | yes |
| GSE123902 | MSK_LX701 | METASTASIS | LUAD |  788 |   57 |  160 | 0.203 | 0.458 | 22.8 | 0.399 | 0.715 | 0.265 | yes |
| GSE205335 | P0031 | Lung | ADC | 6428 |  319 | 3824 | 0.595 | 0.949 | 59.2 | 0.346 | 0.782 | 0.304 | yes |
| GSE205335 | P1006 | Effusion,Lung | ADC | 6090 | 1282 | 2998 | 0.492 | 1.338 | 74.2 | 0.197 | 0.370 | 0.502 | yes |
| GSE205335 | P1015 | LN | ADC | 1195 |  291 |  473 | 0.396 | 0.768 | 41.6 | 0.479 | 0.706 | 0.186 | yes |
| GSE205335 | P1016 | Liver,LN | SCLC | 6924 | 5336 |  763 | 0.110 | 1.455 | 82.6 | 0.263 | 0.847 | 0.182 | yes |
| GSE205335 | P1017 | Liver,LN | SQ | 6266 | 2186 | 2370 | 0.378 | 1.447 | 64.3 | 0.555 | 0.766 | 0.356 | yes |
| GSE205335 | P1018 | Bronchus | ADC | 5029 | 1180 | 2546 | 0.506 | 1.457 | 67.0 | 0.625 | 1.051 | 0.416 | yes |
| GSE205335 | P1025 | LN | SCLC | 1222 |  980 |  128 | 0.105 | 1.694 | 92.9 | 0.514 | 0.614 | 0.360 | yes |
| GSE205335 | P1027 | LN | ADC | 4600 | 1318 | 2619 | 0.569 | 2.120 | 78.1 | 0.612 | 1.015 | 0.443 | yes |
| GSE205335 | P1030 | LN | ADC |  861 |  603 |  140 | 0.163 | 1.617 | 61.7 | 0.298 | 0.717 | 0.368 | yes |
| GSE205335 | P1037 | LN | SQ | 4983 | 3671 |  611 | 0.123 | 0.954 | 86.1 | 0.193 | 0.550 | 0.096 | yes |
| GSE205335 | P1056 | Liver,LN | NUT | 5103 | 2386 | 1354 | 0.265 | 0.810 | 49.7 | 0.249 | 0.532 | 0.162 | yes |
| GSE205335 | P1062 | LN | ADC | 3057 |  188 | 1934 | 0.633 | 0.986 | 45.2 | 0.542 | 1.131 | 0.268 | yes |
| GSE205335 | P1063 | Bronchus,LN | ADC | 4297 |  131 | 3144 | 0.732 | 0.737 | 35.1 | 0.424 | 1.036 | 0.207 | yes |
| GSE205335 | P1072 | LN | SCLC | 1273 |  406 |  586 | 0.460 | 0.936 | 58.9 | 0.105 | 0.392 | 0.276 | yes |
| GSE205335 | P1076 | LN | ADC | 6975 |  849 | 2089 | 0.299 | 0.924 | 59.4 | 0.623 | 0.729 | 0.258 | yes |
| GSE205335 | P1079 | LN | ADC | 1277 |  543 |  501 | 0.392 | 0.967 | 51.4 | 0.528 | 1.249 | 0.311 | yes |
| GSE205335 | P1084 | LN | ADC | 3534 |  196 | 2350 | 0.665 | 2.256 | 86.7 | 0.550 | 1.335 | 0.601 | yes |
| GSE205335 | P1089 | LN | ADC | 1609 | 1063 |  272 | 0.169 | 2.151 | 87.8 | 0.124 | 0.367 | 0.407 | yes |
| GSE205335 | P1090 | LN | SQ |  976 |  566 |  247 | 0.253 | 0.688 | 38.2 | 0.703 | 1.182 | 0.310 | yes |
| GSE205335 | P1115 | Liver | SCLC | 4397 | 3769 |  271 | 0.062 | 1.623 | 88.9 | 0.182 | 0.385 | 0.287 | yes |
| GSE205335 | P1119 | LN | ADC | 2386 | 1222 |  621 | 0.260 | 0.561 | 36.1 | 0.597 | 1.070 | 0.220 | yes |
| GSE205335 | P4001 | Lung | ADC | 2308 |   27 | 1801 | 0.780 | 0.733 | 48.1 | 0.371 | 1.086 | 0.154 | yes |

## Patient-level Spearman

### CLDN4 %pos vs T/NK fraction

| split | n | ρ [95% CI] | p | note |
|---|---:|---|---:|---|
| GSE123902 | 13 | -0.434 [-0.795, 0.154] | 0.138 | marker epithelial |
| GSE205335 | 22 | -0.435 [-0.724, -0.017] | 0.0429 | author malignant |
| combined DL | 35 | -0.435 [-0.680, -0.102] | 0.0121 | I²=0% |
| stacked patients | 35 | -0.399 [-0.647, -0.076] | 0.0174 | cohort mix |

### CLDN4 mean vs T/NK fraction

| split | n | ρ [95% CI] | p |
|---|---:|---|---:|
| GSE123902 | 13 | -0.132 [-0.637, 0.452] | 0.668 |
| GSE205335 | 22 | -0.200 [-0.574, 0.242] | 0.371 |
| combined DL | 35 | -0.177 [-0.495, 0.183] | 0.335 |

### Malignant programs vs CLDN4 %pos

| contrast | GSE123902 ρ (p) | GSE205335 ρ (p) | combined DL ρ (p) |
|---|---|---|---|
| IFN | -0.286 (0.344) | -0.301 (0.174) | -0.296 (0.101) |
| MHC-I/APM | +0.044 (0.887) | -0.433 (0.0441) | -0.251 (0.304) |
| TJ (CLDN4 held out) | +0.533 (0.0607) | +0.457 (0.0326) | +0.484 (0.00446) |

## Malignant Q4 vs Q1 IFN / MHC / T/NK

Within-cohort CLDN4-%pos quartiles among eligible patients. Mid quartiles unused. Stacked Q4 vs Q1 is the pair contrast.

| family | n_Q1 | n_Q4 | Δ median (Q4−Q1) | r | MW p |
|---|---:|---:|---:|---:|---:|
| T/NK fraction | 10 | 9 | -0.325 | -0.756 | 0.00623 |
| IFN | 10 | 9 | -0.184 | -0.667 | 0.016 |
| MHC-I/APM | 10 | 9 | -0.344 | -0.533 | 0.055 |
| TJ (CLDN4 held out) | 10 | 9 | +0.115 | +0.578 | 0.0373 |

Single-cohort Q4 vs Q1 tails are thinner (GSE123902 4 vs 3 on %pos; GSE205335 6 vs 6). GSE205335 Q4 mixes SCLC with ADC; that histology mix is part of the honest n.

## Seurat plots

- `results/figures/fig_dimplot_dataset.png` — `DimPlot` by dataset (Harmony UMAP)
- `results/figures/fig_dimplot_compartment.png` — `DimPlot` malignant vs T/NK
- `results/figures/fig_dimplot_patient.png` — `DimPlot` by patient
- `results/figures/fig_featureplot_cldn4.png` — `FeaturePlot` CLDN4
- `results/figures/fig_vlnplot_cldn4.png` — `VlnPlot` CLDN4 by compartment
- `results/figures/fig_patient_cldn4_vs_tnk.png` — patient scatter %pos vs T/NK
- `results/figures/fig_patient_cldn4_vs_ifn.png` — patient scatter vs IFN
- `results/figures/fig_patient_cldn4_vs_mhc.png` — patient scatter vs MHC
- `results/figures/fig_q4q1_ifn.png` / `fig_q4q1_mhc.png` / `fig_q4q1_tnk.png`

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- GSE123902 marker epithelial is **not** a CNV-malignant call.
- GSE205335 author malignant includes SCLC / NUT / SQ / ADC; histology is not hidden.
- The Harmony cap (≤150 mal + ≤150 T/NK per patient) is for UMAP only. T/NK fraction and program scores use the full patient.
- Restriction (CLDN4 in epithelium vs T/NK) is not an infiltration / immune-cold claim by itself.
- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.
- No ICI / MPR / RECIST / survival test is the primary claim.
- Not merged with GSE148071, GSE127465, or the concordant-4 four-way.
- The PR #459 n=35 ρ=−0.522 row is given and is not recomputed as the thesis test.

## Reproduce

```bash
bash methods/seurat_integrate_123902_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_integrate_123902_205335_cldn4/scripts/01_seurat_harmony_integrate.R
```

Seurat 5.5.1. Harmony 2.0.5. R 4.6.1.


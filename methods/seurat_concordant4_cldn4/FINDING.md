# Seurat concordant-4 CLDN4-only

ADDITIVE. **CLDN4-only.** No TACSTD2∩CLDN4 dual-high. Not a mega-merge.
The four sets that already point the same way: **GSE123902 + GSE131907 +
GSE205335 + GSE189357**. Not GSE148071, GSE127465, GSE207422, GSE154826,
or CD45+/T-only extracts.

Thesis (already correct; not re-derived): CLDN4-high malignant cells have
lower own IFN/MHC-I, and patients have lower T/NK. CLDN4-low/KD opens IFN/MHC.
Matching extras here are IFN/MHC **DOWN** in CLDN4-high.

Primary stack is **R + Seurat 5.5.1 + Harmony 2.0.5** (`RunHarmony` on the Seurat object, `group.by.vars = dataset`).
Python was used only to stream the GSE131907 genes×cells text into a sparse
atlas subset. Patient-level numbers use **full-unit** counts (not the UMAP cap).

Patient / donor / sample is the unit. Do not quote cell counts as n.
GSE123902 = donor. GSE131907 = sample (tumor-bearing). GSE205335 = patient
(RECIST not required). GSE189357 = patient. p-values are descriptive.

## Honest n

- **n_units = 65** (13 donors + 21 samples + 22 patients + 9 patients).
- **n_cells in the Seurat object (after QC + ≤350/unit cap) = 22653**.
- Do not replace the unit n with this cell count.

- GSE123902: 13 units, 4453 cells used
- GSE131907: 21 units, 7350 cells used
- GSE189357:  9 units, 3150 cells used
- GSE205335: 22 units, 7700 cells used

UMAP cell class (used cells, not the test n): B=1552; malignant=7424; myeloid=4220; NK=1814; other=1285; T=6358.

## 1. Patient-level malignant CLDN4 vs T/NK fraction

Primary score = malignant CLDN4 **%pos** (full unit, not the cap).
Pooling = DerSimonian–Laird on Fisher-z of the four cohort Spearmans.
Q4 vs Q1 = **within-cohort quartiles stacked**, then Mann–Whitney on T/NK
(rank-biserial r). Marker gate (GSE123902, GSE189357): (EPCAM|KRT8|KRT18|KRT19)>0 AND PTPRC==0 vs (CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0.
Author labels (GSE131907, GSE205335). `frac_tnk = n_tnk / n_cells`.

| score | k | N | ρ (p, I², 95% CI) | stacked Q4 vs Q1 r (n_Q1/n_Q4, p) |
|---|---:|---:|---|---|
| %pos | 4 | 65 | -0.531 (1.65e-05, I²=0.0%, -0.697 to -0.312) | -0.724 (19/16, 0.0002879) |

### Singles (context; not the new number)

| cohort | unit | n | ρ | p |
|---|---|---:|---:|---:|
| GSE123902 | donor | 13 | -0.659 | 0.01423 |
| GSE131907 | sample | 21 | -0.522 | 0.0152 |
| GSE205335 | patient | 22 | -0.435 | 0.04286 |
| GSE189357 | patient |  9 | -0.600 | 0.08762 |

## 2. Malignant Q4 vs Q1 IFN / MHC / TJ

Patient-pseudobulk OLS on log2(TMM-CPM+1) family means, `~ cohort + CLDN4_Q4`.
Within-cohort %pos quartiles. Units with n_malignant < 30 are out of this DE
(P4001-style noise gate). Positive logFC = higher in CLDN4-high.
IFN = Hallmark IFNα ∪ IFNγ. MHC = custom MHC-I/APM. TJ = KEGG ∪ GOBP organization
plus CDH1/VIM/ZEB1; **CLDN4 held out**.

| family | n_Q1 / n_Q4 | n_genes | logFC | p |
|---|---|---:|---:|---:|
| IFN | 18 / 16 | 222 | -0.609 | 0.05161 |
| MHC-I/APM | 18 / 16 |  21 | -0.803 | 0.05062 |
| TJ | 18 / 16 | 211 | 0.011 | 0.965 |

Expected under the thesis: IFN down, MHC-I/APM down, TJ up or held.
Do not quote N=65 as the DE n. Do not quote 22653 cells as n.

## 3. Seurat / Harmony UMAP

- QC: n_genes ≥ 200, n_UMI ≥ 500, mitochondrial % < 20.
- Cap: ≤350 cells / unit (memory).
- Inner gene join → NormalizeData → VST 2000 HVG → ScaleData → PCA 30.
- Harmony: `group.by.vars = dataset` only (sample was not a second key).
- Neighbors / UMAP / FindClusters resolution 0.6 on the Harmony embedding.
- Cluster labels: author-malignant fraction ≥ 0.50 → malignant; else highest
  mean log-normalized marker score (EPCAM/KRT/CLDN4, CD3, NK, myeloid, B).
  Endothelium / fibroblast / mixed / low-score → **other**.

DimPlots: `results/figures/DimPlot_dataset.png`, `DimPlot_cell_class.png`.
Seurat object (local, 218 MB): `results/objects/seurat_harmony_integrated.rds`.
Embeddings + meta (committed): `results/objects/seurat_harmony_embeddings.rds`.
Patient table: `results/tables/patient_units.tsv`.

## Reproduce

```
Rscript methods/seurat_concordant4_cldn4/analyze.R
```


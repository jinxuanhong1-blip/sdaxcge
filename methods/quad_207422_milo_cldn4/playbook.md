# QUAD Milo vs malignant CLDN4

**Scope.** Additive CLDN4-only neighbourhood DA on four public NSCLC
scRNA matrices. Include GSE207422. Do not redo the 207422-only Milo
(PR #323, n=7, SpatialFDR empty).

**Not miloR.** Same fallback as the single-dataset folders.

## Graph

Per dataset. GSE131907 is further split tLung / mBrain. PCA per-sample
mean centering. Harmony is the allowed alternative and was not used
(memory + chemistry/label heterogeneity + treatment confounding).

## DA

`prop[s, i] = n_cells(s in i) / n_cells(s)`

Spearman of `prop` vs sample mean log1p-CP10k CLDN4 in malignant cells
(sample dropped if <10 such cells). SpatialFDR = k-distance weights.

## Unit

Patient / sample. Do not cite cell count as *n*.

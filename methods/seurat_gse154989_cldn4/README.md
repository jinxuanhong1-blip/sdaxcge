# Seurat — GSE154989 KP epithelium, Cldn4-only

ADDITIVE public mouse. Cldn4-only. Thesis already correct. No dual-high.

Primary analysis is **R + Seurat** (`CreateSeuratObject`) on the public GEO
processed matrix `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5`.

This is FACS `tdTomato+ / CD45− / CD11b− / TER119− / CD31−` Smart-seq2
epithelium (Marjanovic et al., *Cancer Cell* 2020). There is **no T/NK
compartment**. Residual `Cd3d` / `Nkg7` / `Ptprc` is leak. T/NK fraction is
not scored.

Unit of analysis = biological mouse (`mouseID` with trailing `_T#` stripped)
with enough cells (primary: KP, ≥20 cells).

See `FINDING.md` for the locked result. Mouse-level table:
`tables/mouse_units.tsv` and `tables/primary_KP_n20.tsv`.

```bash
Rscript methods/seurat_gse154989_cldn4/analyze.R
```

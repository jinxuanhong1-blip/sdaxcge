# Concordant-4 malignant CLDN4: patient-level GSVA, ssGSEA, and fgsea

Additive. CLDN4 only. The four cohorts are GSE123902, GSE131907, GSE205335,
and GSE189357. Not GSE148071, GSE127465, GSE207422, or GSE154826.

The unit is the malignant-cell pseudobulk (patient, donor, or sample).
Quartiles are the locked within-cohort CLDN4 % positive labels in
`data/tnk_units.tsv`. Cell-level p-values are not computed.

## Prespecified sets

Scored separately, with CLDN4 removed:

- Hallmark IFN-alpha response
- Hallmark IFN-gamma response
- Custom classical MHC-I / antigen-presentation panel (21 genes, MHC-II out)
- Hallmark EMT
- Hallmark apical junction
- KEGG tight junction

Secondary: GO tight-junction organization and GO bicellular tight-junction
assembly.

## Reproduce

```bash
Rscript methods/concordant4_cldn4_gsva_fgsea/analyze.R
```

Needs GSVA (>= 2.0), fgsea, ggplot2, data.table, and jsonlite.
`N_BOOT` (default 200) and `NPERM_BOOT` (default 1000) set the patient
bootstrap behind the NES intervals. Seed 20260921.

Writeup: `FINDING.md`. Tables: `results/tables/`. Forests: `results/figures/`.

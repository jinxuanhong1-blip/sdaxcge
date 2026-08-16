# xCell reference files (Aran et al., Genome Biology 2017)

Used only to compute the **CD8+ T-cells** xCell score on GSE253564 / GSE248378.

| File | Source |
|---|---|
| `Xcell_signatures.json` | Official xCell 489 signatures, exported as JSON (same gene lists as `xCell.data$signatures` in [dviraran/xCell](https://github.com/dviraran/xCell)). Copy via [sysbio-curie/tipit_benchmark_RNA](https://github.com/sysbio-curie/tipit_benchmark_RNA) `benchmark_RNA/data/deconv_signatures/Xcell_signatures.json`. |
| `xcell_rnaseq_spill_and_calibration.csv` | RNA-seq spillover matrix K plus power / calibration columns from the same export (`Xcell_coef.xlsx` sheet "Spill - Seq" = `xCell.data$spill` for RNA-seq). |

Do not use the microarray ("Spill - array") sheet; these GEO series are RNA-seq FPKM.

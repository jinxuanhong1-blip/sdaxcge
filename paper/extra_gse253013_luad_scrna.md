# Extra scRNA figure · GSE253013 treatment-naïve LUAD

**Placement.** Extra neoadjuvant/ICI-adjacent scRNA panel. Do **not** replace or re-analyze user A3 (GSE207422 malignant TACSTD2 vs T/NK), which is taken as given.

## Honest cohort identity

GSE253013 (Sze, Xiang et al., *Cancer Res* 2024, PMID 38335304) is public LUAD scRNA-seq of **9 treatment-naïve** patients (256,379 cells after author QC; **89** GSM 10x lanes: tumor and adjacent non-tumor lung). The paper discusses immunotherapy resistance as motivation; the sequenced tissues were **not** on-treatment ICI and carry **no public MPR / pCR / R / NR labels** in GEO. This panel is therefore an extra LUAD TME scRNA figure, not a response-stratified neoadjuvant ICI replication.

The series is lung tissue, not blood. After a real gene-name check of the GEO processed object (`GSE253013_all_luad_garnett_temp.rds.gz`), **TACSTD2 and CLDN4 are both present**.

## Suggested results text

In an independent public LUAD scRNA-seq series (GSE253013; 9 treatment-naïve patients, 89 10x lanes), we scored *TACSTD2* and *CLDN4* in marker-defined malignant-like epithelial cells and compared those scores with the per-patient tumor T/NK fraction. No major-pathologic-response or radiographic-response labels are released with this series, so a response contrast was not performed. Patient-level Spearman statistics are reported in `results/gse253013/association_statistics.tsv` (unit = patient, not cell).

## Methods (for the paper methods section)

Processed counts were taken from the GEO supplementary RDS (series matrix has no expression). Because the object is larger than a 16 GB session, a streaming XDR parser extracted a lineage/target gene panel and cell metadata (scripts: `methods/gse253013_*.py`). Lineages were assigned by argmax of mean log1p marker scores. Malignant-like epithelium was epithelial cells with near-zero normal-lung markers (*SFTPA2*, *AGER*, *SCGB1A1*, *SCGB3A1*, *TPPP3*). Per-patient tumor T/NK fraction used T (*CD3D*/*CD3E*/*CD2*) plus NK (*NKG7*/*GNLY*/*FGFBP2*) cells. Associations used Spearman rank correlation on patients with ≥10 malignant-like and ≥20 T/NK cells. GSE207422 was not re-analyzed.

## Figure

`results/gse253013/fig_gse253013_tacstd2_cldn4_vs_tnk.png`

Caption template (fill n/ρ/p from `summary.json` after the run):

**Extra scRNA.** Malignant-like *TACSTD2* and *CLDN4* versus per-patient T/NK fraction in GSE253013 treatment-naïve LUAD (tumor lanes). No public ICI response / MPR labels. Statistics: Spearman, patient unit.

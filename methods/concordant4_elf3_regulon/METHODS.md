# Methods — concordant-4 ELF3 / GRHL module

ADDITIVE. The locked concordant-4 result is not re-derived as a new claim:
malignant CLDN4 % positive versus T/NK fraction, DerSimonian–Laird
ρ = −0.531, N = 65. That number is a pipeline check. If the check misses,
the run stops.

Datasets, and only these: GSE123902 (donor), GSE131907 (sample),
GSE205335 (patient), GSE189357 (patient). Not GSE148071, GSE127465,
GSE207422, GSE154826, or GSE200563. TACSTD2 is not a gate. No dual-high split.
No Visium. No private 8KL matrix.

The unit is the patient / donor / sample already locked in
`data/locked_units.tsv` (n = 65). Cell counts are not n. p-values are
descriptive.

## Compartments (same definitions as the locked Seurat run)

| Cohort | Malignant | T/NK | Unit |
| --- | --- | --- | --- |
| GSE123902, GSE189357 | (EPCAM or KRT8 or KRT18 or KRT19) > 0 and PTPRC = 0 | (CD3D or CD3E or CD8A or NKG7 or GNLY or KLRD1) > 0, and not malignant | donor or patient |
| GSE131907 | author `Cell_subtype == Malignant cells` | author `Cell_type` in {T lymphocytes, NK cells} | sample |
| GSE205335 | author `lineage.sub == Malignant cells` | author `lineage.total == T/NK cells` | patient |

GSE123902 drops NORMAL captures. GSE131907 keeps the 21 tumor-bearing
samples with at least 20 author-malignant cells. GSE205335 drops captures
whose tissue starts with "Normal". `frac_tnk = n_tnk / n_cells` on that
same denominator. Counts are full-unit UMIs, not a 350-cell UMAP cap.

## What is quantified

Pre-specified genes: ELF3, GRHL2, TACSTD2, CLDN4, CLDN7.
GRHL1 and GRHL3 are secondary context. They are not part of the four-gene module.

1. Within malignant cells. For each unit with at least 30 malignant cells,
   Spearman correlation of UMI counts for each pre-specified pair (rank-equivalent
   to log1p). Summary = cohort median, and the median across all such units.
   P4001 (27 malignant cells) is out of this layer only.
2. Across units. Spearman of the malignant mean of log1p(UMI), within cohort,
   then DerSimonian–Laird on Fisher z. A % positive companion uses the same pairs.
3. Co-detection. Per unit, fraction of malignant cells with ELF3, TACSTD2,
   CLDN4, and CLDN7 all > 0, and that fraction divided by the product of the
   four marginal detection rates.
4. Regulon. DoRothEA A+B+C weighted mean (wmean) of ELF3 targets and of GRHL2
   targets on the malignant mean of log1p(CP10k). GRHL2 is also scored with
   CLDN4 removed, because that edge is in the GRHL2 set. A TF is scored only
   if at least 5 targets are present. This is not VIPER, not pySCENIC, and not
   lung ChIP.
5. Specificity. Within each cohort, Spearman of malignant mean log1p(ELF3)
   (and GRHL2) against every detected gene. Percentile of the pre-specified
   genes, plus EPCAM, KRT8/18/19, KRT5, CDH1, and PTPRC.
6. Patient-level T/NK. Within-cohort Spearman versus `frac_tnk`, then the
   same DerSimonian–Laird pool. Primary scores are malignant % positive of
   ELF3, TACSTD2, CLDN4, CLDN7, and GRHL2; the mean of within-cohort z-scores
   of the four % positive values; the same mean without CLDN4; ELF3 wmean;
   GRHL2 wmean; GRHL2 wmean without CLDN4. CLDN4 % positive is the pipeline
   check, not a new result.

## Software

Python (numpy, scipy, pandas, matplotlib). GSE205335 is read in R with the
Matrix package because GEO ships only a dgCMatrix RDS; the R step writes
summaries and does not run the tests. No Seurat, no Harmony, no decoupleR.

# Methods — winning-pair CytoTRACE2 potency, CLDN4 only

ADDITIVE. Lives under `methods/winpair_cytotrace2_cldn4/`. **CLDN4 only.**
Does **not** use GSE148071. GSE207422 is not added. No TACSTD2∩CLDN4 dual-high gate.

**Question.** On the winning public pair GSE131907+GSE205335 **author-malignant cells**, are CLDN4-high cells more differentiated / barrier-locked (**lower potency**) than CLDN4-low cells? Patient is the unit.

## What was (and was not) run

**Primary:** CytoTRACE2 (`cytotrace2-py`, Kang et al. *Nat Methods* 2025). Human model. Score ∈ [0, 1] (0 = differentiated, 1 = totipotent). Plotting disabled. `max_cores=1` on a 16 GB machine.

This is **not** a residual-`n_genes` dump. If CytoTRACE2 cannot import or score, the documented fallback is the original CytoTRACE algorithm (Gulati et al. *Science* 2020): Pearson of each gene vs per-cell gene counts → mean of the top 200 genes → KNN smooth → rank-scale to [0, 1]. FINDING.md states which score was primary.

Sensitivity also always stores the Gulati 2020 score so the two methods can be compared.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Malignant kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text + author annotation | `Cell_subtype == Malignant cells` | tLung (0 author-malignant; cells are tS1/tS2/tS3); nLung; PE unlabeled; 2.86 GB log2TPM; EGA FASTQ |
| GSE205335 | Ahn / Lee, *eLife* 2024 | UMI dgCMatrix + author identity + SOFT | `lineage.sub == Malignant cells` on non-normal tissues | Normal Lung / LN / Brain; EGA FASTQ |

GSE148071 is not downloaded and not analyzed.

Winning-pair membership is taken as given from the CLDN4-first combinatorial search (author malignant CLDN4 %pos vs T/NK). This folder asks a **different** question (potency on malignant cells), not a T/NK re-rank.

## Locked choices

| Item | Choice |
| --- | --- |
| Universe | Author malignant only. No CopyKAT / inferCNV. |
| Cap | ≤200 cells / unit after QC (honest n reports catalog vs analysis) |
| QC | UMI≥200 and n_genes≥200 |
| CLDN4 | `log1p(CP10k)` from raw UMI |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | GSE131907 `Sample` + GSE205335 `patient` |
| Primary test | Spearman of patient-mean CLDN4 vs patient-mean CytoTRACE2. BH inside that list only. |
| Extra | Within-unit CLDN4-high vs low potency and barrier (min 8 cells/arm) |
| Unused | GSE148071; GSE207422; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a potency label |

## Tests

Unit = patient/sample. Cell-level Spearman is stored and labeled exploratory.

1. Spearman of patient-mean malignant CLDN4 vs CytoTRACE2.
2. Spearman of patient-mean CLDN4 vs CLDN4-excluded barrier/keratin.
3. Extra: paired Wilcoxon of potency (and barrier) in within-unit CLDN4-high vs CLDN4-low.

Two-sided p. Honest n. Every test is in `results/tables/stats.tsv`.
The done criterion is `results/tables/patient_cldn4_vs_potency.tsv`.

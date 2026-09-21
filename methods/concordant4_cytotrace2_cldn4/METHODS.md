# Methods — concordant-4 CytoTRACE2 potency and stemness, CLDN4 only

ADDITIVE. Lives under `methods/concordant4_cytotrace2_cldn4/`.

**Question.** On concordant-4 malignant cells, is CLDN4-high more differentiated (lower potency, lower stemness) and more barrier-like than CLDN4-low?

## Datasets

| Cohort | Malignant rule | Unit | Notes |
| --- | --- | --- | --- |
| GSE123902 Laughney 2020 | `(EPCAM\|KRT8\|KRT18\|KRT19)>0` and `PTPRC==0` | donor | PRIMARY + METASTASIS. NORMAL dropped. Marker gate, not CNV. 36.5 GB H5 not used. |
| GSE131907 Kim 2020 | `Cell_subtype==Malignant cells` | sample | Author label. tLung has 0 author-malignant cells. log2TPM text not used. |
| GSE205335 Ahn/Lee | `lineage.sub==Malignant cells`, non-normal tissue | patient | Author label. Histology is mixed. |
| GSE189357 Zhu 2022 | same marker gate as GSE123902 | patient | TD1–TD9. Not CNV. No spatial GSE189487. |

Not used: GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, dual-high TACSTD2∩CLDN4.

## Scores

| Score | Definition |
| --- | --- |
| CytoTRACE2 | `cytotrace2-py` 1.1.0.4, human, raw UMI, **per dataset**. 0 = differentiated, 1 = totipotent. Smoothed score is primary. preKNN score is sensitivity. |
| Gulati 2020 | Top-200 gene-count correlates, KNN smooth, rank to [0,1], per dataset. Sensitivity only. |
| Ben-Porath ES1 | MSigDB c2.cgp v2023.2.Hs `BENPORATH_ES_1`. Mean z of log1p(CP10k) inside the dataset. |
| Wong ESC | `WONG_EMBRYONIC_STEM_CELL_CORE`, same z-score. |
| Core nine | `BENPORATH_ES_CORE_NINE`. Sensitivity. |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR. **CLDN4 is not in this score.** |
| CLDN4 | Patient-mean log1p(CP10k). % positive is sensitivity. |

Stemness sets drop CLDN4, TACSTD2, EPCAM, and the barrier/keratin genes before scoring, so stemness is not the epithelial gate and not the barrier score.

## Tests

Unit = patient / donor / sample with ≥20 malignant cells after QC (UMI≥200, genes≥200) and a cap of ≤200 cells/unit.

1. Within each cohort, Spearman of patient-mean CLDN4 vs patient-mean score.
2. DerSimonian–Laird random-effects meta on Fisher-z of those four ρ values. This is the primary number.
3. BH across the seven primary metas only.
4. Extra: within-unit CLDN4-high vs low (tertile; median split if the unit is small), Wilcoxon signed-rank on the patient arm means. ≥8 cells per arm.

Leave-one-cohort-out, pooled (non-meta) Spearman, %pos, Gulati, preKNN, ADC-only, and GSE123902 primary vs metastasis are sensitivity and are not in the BH family.

Cell counts are not the inferential n.

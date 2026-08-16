# Claim A8/A9 — TROP2-high keratin / tight-junction GSEA

Computed numbers only. No fabricated statistics.

## Question

Are keratin and tight-junction (TJ) programs enriched among genes
positively associated with **TACSTD2 (TROP2)** in lung cancer, and do
the TJ genes **CLDN1, CLDN4, CLDN7, F11R, PARD3** themselves track TROP2-high?

## Method (identical in every dataset)

- Ranking metric for pre-ranked GSEA: Spearman correlation of each gene with TACSTD2.
- TROP2-high vs TROP2-low: top vs bottom tertile of TACSTD2.
- A panel gene is called `up_in_trop2_high` only if **all** of: rho > 0,
  log2FC(high vs low) > 0, Welch t-test BH-FDR < 0.05.
- Gene sets: Enrichr GO BP/CC 2021, KEGG 2021 Human, Reactome 2022
  (local GMT in `data/genesets/`; see provenance JSON).
- GSEA: gseapy prerank, 1000 permutations, min_size=3.

## Datasets

| id | description | role |
|---|---|---|
| tcga_luad | TCGA-LUAD STAR log2(TPM+1), primary tumours, 1/patient | primary bulk |
| gse207422 | GSE207422 NSCLC pre-treatment bulk log2TPM (n=24) | primary bulk |
| tcga_lusc | TCGA-LUSC, same processing | squamous histology check |
| gse131907_tlung_pb | GSE131907 tLung epithelial sample pseudobulk | supporting, small n |
| gse131907_alltumor_pb | + tL/B, mLN, mBrain, PE epithelial pseudobulk | sensitivity |
| gse131907_tlung_cell | tLung epithelial cells, cell-level | exploratory (cells not independent) |

## Panel gene calls

| dataset | role | gene | n_total | log2FC_high_vs_low | spearman_rho | welch_fdr | up_in_trop2_high | same_direction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TCGA-LUAD bulk | primary | CLDN1 | 516 | 1.56 | 0.4225 | 2.135e-16 | True | True |
| TCGA-LUAD bulk | primary | CLDN4 | 516 | 1.199 | 0.5356 | 9.482e-23 | True | True |
| TCGA-LUAD bulk | primary | CLDN7 | 516 | 0.2924 | 0.1254 | 0.01417 | True | True |
| TCGA-LUAD bulk | primary | F11R | 516 | 0.4886 | 0.3127 | 6.049e-10 | True | True |
| TCGA-LUAD bulk | primary | PARD3 | 516 | 0.1665 | 0.1141 | 0.1045 | False | True |
| GSE207422 NSCLC bulk | primary | CLDN1 | 24 | 5.062 | 0.8783 | 0.006833 | True | True |
| GSE207422 NSCLC bulk | primary | CLDN4 | 24 | 3.886 | 0.8583 | 0.09796 | False | True |
| GSE207422 NSCLC bulk | primary | CLDN7 | 24 | 2.917 | 0.8348 | 0.0503 | False | True |
| GSE207422 NSCLC bulk | primary | F11R | 24 | 2.559 | 0.8704 | 0.002037 | True | True |
| GSE207422 NSCLC bulk | primary | PARD3 | 24 | 2.029 | 0.6626 | 0.2495 | False | True |
| TCGA-LUSC bulk | histology_check | CLDN1 | 501 | 1.919 | 0.4258 | 1.567e-16 | True | True |
| TCGA-LUSC bulk | histology_check | CLDN4 | 501 | 1.045 | 0.3914 | 2.452e-12 | True | True |
| TCGA-LUSC bulk | histology_check | CLDN7 | 501 | 0.9462 | 0.3721 | 1.225e-11 | True | True |
| TCGA-LUSC bulk | histology_check | F11R | 501 | 0.4712 | 0.2336 | 2.052e-07 | True | True |
| TCGA-LUSC bulk | histology_check | PARD3 | 501 | 0.3046 | 0.1365 | 0.003199 | True | True |

## Honest intersection of the five genes

Intersection is **not** a genome-wide overlap hunt. It is the subset of
{CLDN1, CLDN4, CLDN7, F11R, PARD3} that pass the same call in the named datasets.

| gene | pass_TCGA_LUAD | pass_GSE207422 | pass_TCGA_LUSC | pass_GSE131907_tLung_pseudobulk | same_direction_GSE131907_tLung_pseudobulk | pass_GSE131907_tLung_cell_exploratory | intersection_LUAD_and_GSE207422 | intersection_LUAD_GSE207422_and_GSE131907_pb | intersection_LUAD_GSE207422_and_GSE131907_pb_direction | intersection_LUAD_LUSC_GSE207422 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CLDN1 | True | True | True | None | None | None | True | None | None | True |
| CLDN4 | True | False | True | None | None | None | False | None | None | False |
| CLDN7 | True | False | True | None | None | None | False | None | None | False |
| F11R | True | True | True | None | None | None | True | None | None | True |
| PARD3 | False | False | True | None | None | None | False | None | None | False |

- **LUAD ∩ GSE207422 (primary):** ['CLDN1', 'F11R']
- **LUAD ∩ LUSC ∩ GSE207422:** ['CLDN1', 'F11R']
- **LUAD ∩ GSE207422 ∩ GSE131907 tLung pseudobulk (strict):** not yet computed
- **LUAD ∩ GSE207422 ∩ GSE131907 tLung pseudobulk (direction only):** not yet computed

PARD3 is listed because it is in the a priori panel, not because it passed.
A False in the table is a real negative under this criterion, not a missing value.

## GSEA

| dataset | role | Term | NES | NOM_p | FDR_q | significant_fdr05 |
| --- | --- | --- | --- | --- | --- | --- |
| TCGA-LUAD bulk | primary | KEGG_TIGHT_JUNCTION | 2.213 | 0 | 0 | True |
| TCGA-LUAD bulk | primary | GOCC_BICELLULAR_TIGHT_JUNCTION | 2.089 | 0 | 0 | True |
| TCGA-LUAD bulk | primary | REACTOME_KERATINIZATION | 2.049 | 0 | 0 | True |
| TCGA-LUAD bulk | primary | GOBP_KERATINOCYTE_DIFFERENTIATION | 2.038 | 0 | 0 | True |
| TCGA-LUAD bulk | primary | GOCC_CORNIFIED_ENVELOPE | 2.029 | 0 | 0 | True |
| TCGA-LUAD bulk | primary | REACTOME_TIGHT_JUNCTION_INTERACTIONS | 1.954 | 0 | 0 | True |
| TCGA-LUAD bulk | primary | GOCC_KERATIN_FILAMENT | 1.776 | 0.001264 | 0.0001622 | True |
| TCGA-LUAD bulk | primary | GOBP_TIGHT_JUNCTION_ORGANIZATION | 1.456 | 0.04713 | 0.02129 | True |
| GSE207422 NSCLC bulk | primary | GOCC_CORNIFIED_ENVELOPE | 2.578 | 0 | 0 | True |
| GSE207422 NSCLC bulk | primary | REACTOME_KERATINIZATION | 2.525 | 0 | 0 | True |
| GSE207422 NSCLC bulk | primary | GOBP_KERATINOCYTE_DIFFERENTIATION | 2.307 | 0 | 0 | True |
| GSE207422 NSCLC bulk | primary | REACTOME_TIGHT_JUNCTION_INTERACTIONS | 2.16 | 0 | 0 | True |
| GSE207422 NSCLC bulk | primary | GOCC_KERATIN_FILAMENT | 2.092 | 0 | 0 | True |
| GSE207422 NSCLC bulk | primary | GOCC_BICELLULAR_TIGHT_JUNCTION | 2.045 | 0 | 0.0001892 | True |
| GSE207422 NSCLC bulk | primary | KEGG_TIGHT_JUNCTION | 1.683 | 0 | 0.00227 | True |
| GSE207422 NSCLC bulk | primary | GOBP_TIGHT_JUNCTION_ORGANIZATION | 1.379 | 0.08359 | 0.05448 | False |
| TCGA-LUSC bulk | histology_check | REACTOME_KERATINIZATION | 3.842 | 0 | 0 | True |
| TCGA-LUSC bulk | histology_check | GOCC_CORNIFIED_ENVELOPE | 3.479 | 0 | 0 | True |
| TCGA-LUSC bulk | histology_check | GOBP_KERATINOCYTE_DIFFERENTIATION | 3.136 | 0 | 0 | True |
| TCGA-LUSC bulk | histology_check | GOCC_KERATIN_FILAMENT | 2.641 | 0 | 0 | True |
| TCGA-LUSC bulk | histology_check | REACTOME_TIGHT_JUNCTION_INTERACTIONS | 2.392 | 0 | 0 | True |
| TCGA-LUSC bulk | histology_check | GOCC_BICELLULAR_TIGHT_JUNCTION | 2.173 | 0 | 0 | True |
| TCGA-LUSC bulk | histology_check | KEGG_TIGHT_JUNCTION | 2.107 | 0 | 0.0006228 | True |
| TCGA-LUSC bulk | histology_check | GOBP_TIGHT_JUNCTION_ORGANIZATION | 1.044 | 0.3852 | 0.3689 | False |

## What this does **not** show

- Causality or that TROP2 drives TJ/keratin transcription.
- Protein-level TROP2 or claudin status.
- A genome-wide intersecting signature beyond the five a priori genes.
- Independent samples in the GSE131907 cell-level view.
- Adequately powered sample-level tests in GSE131907 tLung (n is small).


# ADDITIVE extra — malignant-cell patient pseudobulk GSEA (TACSTD2-high vs low)

**Public data only. Additive to the TCGA A8 slide. That TCGA audit was not re-run.**

User A8 (TROP2-high enrich keratin/TJ, Hallmark EMT down) is **taken as given** for TCGA-LUAD/LUSC. This file asks the same prerank question after collapsing **malignant cells to one patient** in public lung tumor scRNA.

## 一句话结论 / TL;DR

In 4/4 meta-eligible cohorts, TACSTD2-high malignant-cell patient pseudobulk is **keratin/TJ-up and Hallmark EMT-up** (`keratin_TJ_up_Hallmark_EMT_opposite`). Hallmark EMT-down is **not** seen. IFN and MHC-I NES flip by cohort (2/4 up). **No** gated cohort is `supportive`. Per-cohort median verdicts: GSE207422=keratin_TJ_up_Hallmark_EMT_opposite (n_high=5, n_low=5); GSE241934=keratin_TJ_up_Hallmark_EMT_opposite (n_high=17, n_low=17); GSE131907=keratin_TJ_up_Hallmark_EMT_opposite (n_high=10, n_low=10); GSE205335=keratin_TJ_up_Hallmark_EMT_opposite (n_high=10, n_low=10); GSE291670=n_too_small (n_high=3, n_low=3); GSE253013=not_run (n_high=0, n_low=0). Meta (arms ≥5): KEGG_TIGHT_JUNCTION median NES +1.655 (4/4 up; Stouffer p=0.046); GOBP_KERATINIZATION median NES +1.438 (4/4 up; Stouffer p=0.046); HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION median NES +1.687 (4/4 up; Stouffer p=0.046); HALLMARK_INTERFERON_GAMMA_RESPONSE median NES +0.452 (2/4 up; Stouffer p=1.000); CUSTOM_MHC_I_ANTIGEN_PRESENTATION median NES -0.209 (2/4 up; Stouffer p=1.000). This does not retract the TCGA A8 slide. Do not write “A8 conserved in scRNA” unless a gated cohort is `supportive`. n is hypothesis-generating. Hallmark EMT, not GOBP EMT, decides the EMT arm.

## Why this extra exists

A8 is TCGA bulk. Conservation of keratin/TJ-high (and Hallmark EMT-low) inside **malignant epithelium**, after removing stroma/immune mixture, is a different question. Patient is the unit. Cell-level tests are not reported as GSEA.

## Cohorts

| Cohort | n patients (gated) | n high vs low | Meta-eligible | Note |
|---|---:|---|---|---|
| GSE207422 | 10 | 5 vs 5 | yes | Hu 2023 neoadj PD-1+chemo; marker malignant-like (no public CopyKAT table) |
| GSE241934 | 35 | 17 vs 17 | yes | NEOTIDE + real-world neoadj IO; author tumor Epi (IIT+RWC) |
| GSE131907 | 21 | 10 vs 10 | yes | Kim 2020 atlas; author Malignant cells; drop nLung/nLN |
| GSE205335 | 21 | 10 vs 10 | yes | Park/Ahn/Lee ICI atlas; author Malignant cells; patient from SOFT; drop normal tissues |
| GSE291670 | 6 | 3 vs 3 | no | Neoadj anlotinib+camrelizumab; 6 tumors; marker malignant-like |
| GSE253013 | 0 | — | no | Sze/Xiang 2024 treatment-naive LUAD; 9.3 GB RDS not loaded; n would be < meta gate |

## Pre-specified design

| Piece | Choice | Honest limitation |
|---|---|---|
| Split | TACSTD2 median, Welch t | Quartile only if both arms ≥6 |
| Unit | Patient malignant UMI-sum, log2(CPM+1) | Not cell-level; not mixed-cell bulk |
| Min cells | 30 malignant/Epi per patient | Marker gates ≠ CopyKAT |
| GSEA | A8 engine, 1000 gene-set perm, seed=42 | Gene-set permutation |
| FDR | BH within 15 primary sets | Not nested Broad FDR |
| EMT rule | Hallmark EMT decides | GOBP EMT is recorded, not Hallmark |
| IFN / MHC-I | Reported, undirected | Do not enter `supportive` |
| Meta | Median NES + Stouffer of NES signs; both arms ≥5 | Sign Stouffer ignores magnitude; not a TCGA NES |

Positive NES = enriched in TACSTD2-high.

## Primary NES — TACSTD2 median split

| Gene set | GSE207422 NES | GSE207422 FDR | GSE241934 NES | GSE241934 FDR | GSE131907 NES | GSE131907 FDR | GSE205335 NES | GSE205335 FDR | GSE291670 NES | GSE291670 FDR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +1.588 | 0.002 | +1.779 | 0.003 | +1.594 | 0.002 | +2.894 | 0.001 | -1.118 | 0.198 |
| KEGG_TIGHT_JUNCTION | +1.390 | 0.022 | +1.549 | 0.014 | +1.760 | 0.002 | +2.152 | 0.001 | -1.156 | 0.195 |
| HALLMARK_APICAL_JUNCTION | +1.533 | 0.002 | +1.139 | 0.166 | +1.283 | 0.015 | +2.400 | 0.001 | -1.193 | 0.195 |
| GOBP_KERATINIZATION | +1.491 | 0.045 | +0.953 | 0.346 | +1.385 | 0.019 | +2.476 | 0.001 | +0.740 | 0.551 |
| KRT_EPITHELIAL | +1.140 | 0.224 | +1.664 | 0.021 | +1.811 | 0.007 | +2.301 | 0.001 | +1.018 | 0.398 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | -2.012 | 0.002 | -1.685 | 0.003 | +2.777 | 0.002 | +2.590 | 0.001 | -1.678 | 0.015 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | -2.084 | 0.002 | -2.165 | 0.003 | +2.593 | 0.002 | +2.270 | 0.001 | -0.817 | 0.398 |
| CUSTOM_MHC_I_ANTIGEN_PRESENTATION | -1.960 | 0.002 | -2.024 | 0.003 | +1.589 | 0.019 | +1.542 | 0.025 | +1.006 | 0.398 |

| Cohort | n_high | n_low | Verdict | Underpowered |
|---|---:|---:|---|---|
| GSE207422 | 5 | 5 | keratin_TJ_up_Hallmark_EMT_opposite | False |
| GSE241934 | 17 | 17 | keratin_TJ_up_Hallmark_EMT_opposite | False |
| GSE131907 | 10 | 10 | keratin_TJ_up_Hallmark_EMT_opposite | False |
| GSE205335 | 10 | 10 | keratin_TJ_up_Hallmark_EMT_opposite | False |
| GSE291670 | 3 | 3 | n_too_small | True |
| GSE253013 | 0 | 0 | not_run | nan |

## Meta NES (cohorts with both arms ≥ 5)

| Gene set | k | median NES | n up / k | Stouffer z (signs) | Stouffer p | BH-FDR |
|---|---:|---:|---|---:|---:|---:|
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 4 | +1.687 | 4/4 | +2.000 | 0.046 | 0.068 |
| KEGG_TIGHT_JUNCTION | 4 | +1.655 | 4/4 | +2.000 | 0.046 | 0.068 |
| HALLMARK_APICAL_JUNCTION | 4 | +1.408 | 4/4 | +2.000 | 0.046 | 0.068 |
| GOBP_KERATINIZATION | 4 | +1.438 | 4/4 | +2.000 | 0.046 | 0.068 |
| KRT_EPITHELIAL | 4 | +1.738 | 4/4 | +2.000 | 0.046 | 0.068 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 4 | +0.452 | 2/4 | +0.000 | 1.000 | 1.000 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | 4 | +0.093 | 2/4 | +0.000 | 1.000 | 1.000 |
| CUSTOM_MHC_I_ANTIGEN_PRESENTATION | 4 | -0.209 | 2/4 | +0.000 | 1.000 | 1.000 |

Stouffer here is z = sum(sign(NES_i)) / sqrt(k). It tests sign concordance, not a weighted pooled NES. Do not quote it as a TCGA number.

## What this is not

- Not a re-audit of TCGA A8.
- Not evidence that TACSTD2 *causes* keratin/TJ programs.
- Not a TROP2 × ICI interaction test.
- GSE207422 / GSE291670 malignant calls are marker proxies.
- GSE241934 uses author Epi on resected tumors, not a CNV malignant call.
- GSE253013 was not GSEA'd (file size / RAM / n).

## Files

- `tables/gsea_prerank_all.tsv` — all contrasts
- `tables/gsea_tacstd2_median.tsv` — locked split
- `tables/meta_nes.tsv`
- `tables/verdicts.tsv`
- `figures/nes_heatmap_median.png`
- `figures/nes_meta_and_headline.png`
- `pseudobulk/*_counts.tsv.gz`

Reproduce: `bash scripts/scrna_pseudobulk_gsea_meta/download.sh` then the two Python scripts. Methods: `methods/scrna_pseudobulk_gsea_meta/playbook.md`.

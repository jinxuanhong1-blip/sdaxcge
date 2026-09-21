# SWEEP — GSE334497 Trop2 KO, thesis-aligned signals

**This is a Trop2 (Tacstd2) knockout, not a CLDN4 knockdown.** Wu *et al.*, JITC 2026, GEO GSE334497, 4T1 tumors, frozen sections. The pre-specified 5-vs-5 match (Cldn4 Welch and epithelial-ISG sample score) is unchanged in FINDING.md. This file is the exploratory grid run after that call.

Delta for genotype tests is KO minus WT. Thesis direction: barrier and NHEJ **down**, IFN / APM / immune / STING **up**. Cldn4-split delta is Cldn4-low minus Cldn4-high, with *Cldn4* removed from every set. Human GMT sets were case-folded to mouse symbols (`Cldn4` from `CLDN4`); genes absent from the matrix are dropped.

## Grid

| Item | Count |
|---|---:|
| Rows in `tables/sweep_all.tsv` | 1195 |
| Thesis-directed tests (the BH family) | 1083 |
| Full 5 vs 5 genotype × thesis tests | 76 |
| Of those, raw p < 0.05 | 50 |
| Of those, BH-FDR < 0.05 inside the 5 vs 5 grid | 43 |

Methods: exact permutation of the sample mean z-score, exact permutation of the sample median z-score, competitive Mann–Whitney of gene-level Welch *t* versus the rest of the genome, and the same competitive test on limma-style moderated *t*. Filters: all 10 samples, each leave-one-out, and one named pair that widens the *Cldn4* gap (drop KO `KO165`, the KO with the highest *Cldn4*, and WT `RESUB-170R`, the WT with the lowest *Cldn4*).

BH-FDR in the full 5 vs 5 grid is a descriptor inside that grid. The 76 tests reuse the same 10 tumors and overlapping genes, so they are not 76 independent experiments. Gene-rank p-values treat genes as exchangeable; they are smaller than sample-permutation p-values because the null is not “these 10 tumors could have swapped labels.” Sample permutation is the test that matches n = 5 vs 5. Leave-one-out minima are a search, reported below the full-matrix table.

## Headline on all 10 tumors (KO − WT)

| Program | Set | Mean log2FC | Sample mean-z perm p | Gene-rank mod-t p | Call |
|---|---|---:|---:|---:|---|
| Barrier | KEGG_TIGHT_JUNCTION | -0.249 | 0.004 | 5.07e-08 | both |
| Claudin/TJ core | BARRIER_STRUCT | -0.675 | 0.016 | 8.34e-08 | both |
| Keratinization | GOBP_KERATINIZATION | -2.066 | 0.004 | 1.92e-16 | both |
| Bulk immune | BULK_IMMUNE | +0.648 | 0.008 | 4.99e-07 | both |
| Allograft / immune | HALLMARK_ALLOGRAFT | +0.156 | 0.048 | 6.41e-13 | both |
| Hallmark IFN-γ | HALLMARK_IFNG | +0.127 | 0.107 | 2.37e-12 | gene-rank only |
| Hallmark IFN-α | HALLMARK_IFNA | +0.112 | 0.214 | 6.76e-06 | gene-rank only |
| MHC-I / APM | APM_MHCI | +0.156 | 0.175 | 0.000334 | gene-rank only |
| STING core | STING_CORE | +0.036 | 0.048 | 0.045 | both |
| Epithelial ISG | EPITHELIAL_ISG | +0.035 | 0.452 | 0.140 | no |
| NHEJ core | NHEJ_CORE | +0.032 | 0.520 | 0.647 | no |
| DNA repair hallmark | HALLMARK_DNA_REPAIR | -0.041 | 0.329 | 0.215 | no |

“Both” means raw *p* < 0.05 on the sample permutation **and** on the gene-rank test. It is not a claim that the within-grid BH-FDR is < 0.05. STING is the borderline case: sample perm *p* = 0.048, gene-rank Welch *p* = 0.034, and the BH-FDR inside the 76-test grid is **0.059–0.072**. Its mean log2FC is only **+0.036**. *Cgas* moves the other way (log2FC −0.381). *Sting1* (Welch up *p* = 0.035) and *Irf3* (Welch up *p* = 0.009) are the STING genes that rise. *Ifnb1* is absent from the matrix. NHEJ core does not fall (mean log2FC +0.032). Epithelial ISG does not rise. *Cldn4* as one gene does not clear 0.05 on the full 5 vs 5 (moderated-*t* down *p* = 0.094; Welch 0.126). Dropping the KO with the highest *Cldn4* and the WT with the lowest *Cldn4* (`KO165` and `RESUB-170R`) is what brings the *Cldn4* moderated-*t* down *p* to 0.016. That filter was chosen because it widens the *Cldn4* gap.

On the *Cldn4* median split, 3 of the 5 low-*Cldn4* tumors are Trop2 KO, so that contrast is mixed with genotype. *Cldn4* was removed from the sets before testing.

## Best single method per arm on all 10 tumors

| Arm | Set | Method | Filter / contrast | Effect | p | FDR in full 5v5 grid | FDR in whole sweep |
|---|---|---|---|---:|---:|---:|---:|
| epithelial_program | GOBP_KERATINIZATION | gene_rank_modt | all10 | -2.066 | 1.92e-16 | 1.46e-14 | 2.31e-14 |
| immune | HALLMARK_ALLOGRAFT | gene_rank_modt | all10 | +0.156 | 6.41e-13 | 1.62e-11 | 1.88e-11 |
| ifn_gamma | HALLMARK_IFNG | gene_rank_modt | all10 | +0.127 | 2.37e-12 | 3.61e-11 | 6.42e-11 |
| inflammatory | HALLMARK_INFLAMMATORY | gene_rank_modt | all10 | +0.113 | 1.44e-08 | 1.51e-07 | 1.61e-07 |
| barrier | KEGG_TIGHT_JUNCTION | gene_rank_modt | all10 | -0.249 | 5.07e-08 | 3.51e-07 | 4.58e-07 |
| ifn_alpha | HALLMARK_IFNA | gene_rank_modt | all10 | +0.112 | 6.76e-06 | 2.45e-05 | 2.99e-05 |
| apm | APM_MHCI | gene_rank_modt | all10 | +0.156 | 0.000334 | 0.000975 | 0.001 |
| sting | STING_CORE | gene_rank_welch | all10 | +0.036 | 0.034 | 0.059 | 0.065 |
| ifn_isg | EPITHELIAL_ISG | gene_rank_modt | all10 | +0.035 | 0.140 | 0.180 | 0.199 |
| nhej | HALLMARK_DNA_REPAIR | gene_rank_welch | all10 | -0.041 | 0.193 | 0.231 | 0.260 |

## Best row per arm after sample filters

| Arm | Set | Method | Filter / contrast | Effect | p | FDR in full 5v5 grid | FDR in whole sweep |
|---|---|---|---|---:|---:|---:|---:|
| ifn_gamma | HALLMARK_IFNG | gene_rank_welch | drop_RESUB-169R | +0.180 | 3.05e-20 | — | 1.93e-17 |
| epithelial_program | GOBP_KERATINIZATION | gene_rank_modt | drop_RESUB-171R | -2.547 | 4.34e-19 | — | 1.56e-16 |
| immune | HALLMARK_ALLOGRAFT | gene_rank_modt | drop_RESUB-168R | +0.182 | 4.95e-16 | — | 5.36e-14 |
| inflammatory | HALLMARK_INFLAMMATORY | gene_rank_welch | drop_RESUB-169R | +0.141 | 2.39e-11 | — | 4.97e-10 |
| ifn_alpha | HALLMARK_IFNA | gene_rank_modt | drop_RESUB-169R | +0.184 | 4.14e-11 | — | 8.3e-10 |
| barrier | KEGG_TIGHT_JUNCTION | gene_rank_modt | drop_RESUB-170R | -0.264 | 4.59e-10 | — | 7.78e-09 |
| apm | APM_MHCI | gene_rank_welch | drop_RESUB-168R | +0.263 | 8.05e-08 | — | 6.46e-07 |
| ifn_isg | EPITHELIAL_ISG | gene_rank_welch | drop_RESUB-169R | +0.235 | 0.000185 | — | 0.000674 |
| nhej | HALLMARK_DNA_REPAIR | gene_rank_welch | drop_RESUB-KO163R | -0.093 | 0.002 | — | 0.007 |
| sting | STING_CORE | sample_mean_z_perm | drop_RESUB-168R | +0.398 | 0.008 | — | 0.019 |

## Best *Cldn4*-defined contrast per arm

These contrasts are not the knockout. On the median split, the number of Trop2-KO samples inside the *Cldn4*-low half is in `sweep_all.tsv` (`n_KO_in_low`). *Cldn4* itself is removed from the sets.

| Arm | Set | Method | Filter / contrast | Effect | p | FDR in full 5v5 grid | FDR in whole sweep |
|---|---|---|---|---:|---:|---:|---:|
| barrier | GOBP_TJ_ORGANIZATION | gene_rank_modt | Cldn4low_minus_Cldn4high_median | -0.311 | 2.42e-07 | — | 1.65e-06 |
| immune | HALLMARK_ALLOGRAFT | gene_rank_welch | Cldn4low_minus_Cldn4high_median | +0.170 | 1.77e-05 | — | 7.12e-05 |
| epithelial_program | GOBP_KERATINIZATION | gene_rank_modt | Cldn4low_minus_Cldn4high_top3 | -1.319 | 6.65e-05 | — | 0.000251 |
| inflammatory | HALLMARK_INFLAMMATORY | gene_rank_modt | Cldn4low_minus_Cldn4high_median | +0.132 | 0.007 | — | 0.019 |
| ifn_gamma | HALLMARK_IFNG | gene_rank_welch | Cldn4low_minus_Cldn4high_median | +0.080 | 0.018 | — | 0.038 |
| apm | APM_MHCI | gene_rank_welch | Cldn4low_minus_Cldn4high_median | +0.026 | 0.271 | — | 0.339 |
| sting | STING_CORE | sample_median_z_perm | Cldn4low_minus_Cldn4high_median | +0.112 | 0.353 | — | 0.413 |
| ifn_alpha | HALLMARK_IFNA | sample_mean_z_perm | Cldn4low_minus_Cldn4high_median | +0.024 | 0.456 | — | 0.512 |
| ifn_isg | EPITHELIAL_ISG | sample_median_z_perm | Cldn4low_minus_Cldn4high_median | -0.393 | 0.762 | — | 0.793 |
| nhej | HALLMARK_DNA_REPAIR | spearman_score | Spearman_Cldn4 | -0.333 | 0.827 | — | 0.850 |

## Same specification, barrier down and an up-arm together

| Filter | Method | Barrier set | Barrier p | Up set | Up p | Both < 0.05 |
|---|---|---|---:|---|---:|---|
| drop_RESUB-170R | gene_rank_modt | KEGG_TIGHT_JUNCTION | 4.59e-10 | HALLMARK_ALLOGRAFT | 3.25e-13 | yes |
| drop_RESUB-170R | gene_rank_welch | KEGG_TIGHT_JUNCTION | 1.35e-09 | HALLMARK_ALLOGRAFT | 1.45e-12 | yes |
| drop_KO165_and_RESUB-170R_Cldn4gap | gene_rank_modt | BARRIER_STRUCT | 3.14e-09 | HALLMARK_ALLOGRAFT | 1.4e-13 | yes |
| drop_KO165_and_RESUB-170R_Cldn4gap | gene_rank_welch | BARRIER_STRUCT | 6.13e-09 | HALLMARK_ALLOGRAFT | 5.03e-13 | yes |
| drop_KO162 | gene_rank_modt | BARRIER_STRUCT | 9.22e-09 | HALLMARK_IFNG | 7.83e-16 | yes |
| drop_RESUB-171R | gene_rank_modt | KEGG_TIGHT_JUNCTION | 1.06e-08 | HALLMARK_ALLOGRAFT | 1.73e-09 | yes |
| drop_KO164 | gene_rank_modt | KEGG_TIGHT_JUNCTION | 2.03e-08 | HALLMARK_ALLOGRAFT | 6.57e-12 | yes |
| drop_KO162 | gene_rank_welch | BARRIER_STRUCT | 3.38e-08 | HALLMARK_IFNG | 8.8e-16 | yes |

## *Cldn4* gene, three tests

On all 10 tumors the smallest down-sided p is **mod_t p = 0.094** at log2FC -0.817. Welch down p = 0.126, moderated-t down p = 0.094, MWU down p = 0.210.

Across filters and those three tests, the smallest *Cldn4* down p is **drop_KO165_and_RESUB-170R_Cldn4gap / mod_t / p = 0.016** (log2FC -1.414). That row is a sensitivity, not a replacement for the 5 vs 5 test.

## STING and NHEJ genes (all 10, KO − WT)

STING, up side: Cgas log2FC -0.381, mod-t up p 0.955, Welch up p 0.949; Sting1 log2FC +0.307, mod-t up p 0.049, Welch up p 0.035; Tbk1 log2FC +0.028, mod-t up p 0.441, Welch up p 0.437; Irf3 log2FC +0.323, mod-t up p 0.025, Welch up p 0.009; Ifnb1 absent.

NHEJ, down side: Xrcc5 log2FC -0.134, mod-t down p 0.210, Welch down p 0.171; Xrcc6 log2FC -0.205, mod-t down p 0.108, Welch down p 0.078; Lig4 log2FC -0.249, mod-t down p 0.120, Welch down p 0.124; Prkdc log2FC +0.181, mod-t down p 0.738, Welch down p 0.722; Nhej1 log2FC -0.275, mod-t down p 0.227, Welch down p 0.250.

Figures: `figures/fig_sweep_grid.png` (full 5 vs 5 grid), `figures/fig_sweep_best.png` (sample scores of the winning full-matrix set in each arm).

Reproduce: `python3 scripts/gse334497_cldn4_match/sweep.py`

# hdWGCNA-like CLDN4-stratified modules in concordant-4 malignant cells

ADDITIVE. **CLDN4-only.** Same four datasets as the locked concordant-4 patient result (GSE123902, GSE131907, GSE205335, GSE189357). Not GSE148071, GSE127465, GSE207422, GSE154826, or E-MTAB-13526. This does not replace the locked malignant CLDN4 %pos vs T/NK result (ρ=−0.531, P=1.65×10⁻⁵, n=65). That number was recomputed from the locked patient table as a pipeline check and matched.

The R package hdWGCNA was not installed. The network is hdWGCNA-like: metacells follow `ConstructMetacells` (k=25, max_shared=10, average of the raw UMI in each kNN, then log-normalize), grouped by patient × CLDN4 status. A group with fewer than 50 malignant cells contributes no metacells. Each group is capped at 40 metacells so a large sample does not dominate (the package default is 1000). Neighbors are Euclidean kNN on 15 within-group PCs, not a global Harmony reduction. Before the correlation, metacell expression is z-scored within dataset. Soft-threshold QC: a power is not used when scale-free R² first exceeds 0.80 only after mean connectivity falls below 5. The primary power is then the negative-slope power with mean connectivity between 10 and 100 and the highest R², and the scale-free fit is reported as not reached. The network is signed biweight midcorrelation, signed TOM, `dynamicTreeCut` hybrid (deepSplit=2, minClusterSize=30, pamRespectsDendro=FALSE), then eigengene merge at r>0.75. Genes were the intersection of the four datasets, detected in ≥5% of malignant cells in ≥3 datasets, top 2500 by mean within-dataset variance rank of log1p(UMI). CLDN4 was forced in if it passed detection. TJ and IFN gene sets were not used to choose genes.

CLDN4-positive means malignant UMI > 0. CLDN4-negative means UMI = 0. Two networks were built, one in each stratum. A module is called the TJ/barrier module or the IFN module by one-sided Fisher enrichment (background = network genes; overlap ≥ 5; odds ratio > 1; smallest p), not by its correlation with T/NK. TJ = KEGG tight junction ∪ GOBP tight-junction organization ∪ CDH1/VIM/ZEB1, CLDN4 removed. IFN = Hallmark IFNα ∪ IFNγ. The patient-level score is the mean log2(CPM+1) of module genes in the malignant pseudobulk, with CLDN4 held out of the score. The test is a within-cohort Spearman, then DerSimonian–Laird on Fisher z. n is patients, not cells and not metacells. p-values are descriptive.

## Metacells and soft threshold

| network | metacells | 123902 | 131907 | 205335 | 189357 | power | scale-free R² | mean k |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CLDN4pos | 1667 | 169 | 552 | 634 | 312 | 10 | 0.337 | 13.9 |
| CLDN4neg | 1445 | 175 | 343 | 586 | 341 | 10 | 0.456 | 10.6 |

## Which module is TJ/barrier, which is IFN

| network | call | module | genes | overlap | odds ratio | p | q | CLDN4 module (kME) | TACSTD2 module (kME) |
|---|---|---|---:|---:|---:|---:|---:|---|---|
| CLDN4pos | TJ/barrier: called (q<0.05) | ME4 | 177 | 14 | 4.16 | 5.61e-05 | 0.000561 | ME1 (0.24) | ME4 (0.48) |
| CLDN4pos | IFN: called (q<0.05) | ME3 | 311 | 43 | 17.40 | 1.15e-25 | 3.44e-24 | ME1 (0.24) | ME4 (0.48) |
| CLDN4neg | TJ/barrier: best overlap, not called (q≥0.05) | ME3 | 329 | 15 | 2.21 | 0.010 | 0.062 | absent (NA) | ME5 (0.65) |
| CLDN4neg | IFN: called (q<0.05) | ME7 | 74 | 26 | 34.96 | 2.44e-25 | 5.87e-24 | absent (NA) | ME5 (0.65) |

Hub genes are the highest kME genes inside the module in the table above. A module is a TJ or IFN module only when q<0.05.

- CLDN4pos TJ/barrier ME4: SPTBN1 (0.73), AGRN (0.72), KIAA1522 (0.71), UBE2H (0.70), NCKAP1 (0.70), RBM47 (0.69), TJP1 (0.69), CD46 (0.68), SEC31A (0.68), PARP14 (0.67), CD2AP (0.67), RRBP1 (0.67)
- CLDN4pos IFN ME3: HLA-B (0.80), HLA-C (0.79), B2M (0.75), HLA-A (0.74), S100A6 (0.73), PSMB8 (0.71), PSMB9 (0.68), HLA-F (0.68), CTSS (0.68), TIMP1 (0.67), CIB1 (0.67), MYL12A (0.66)
- CLDN4neg TJ/barrier ME3: S100A6 (0.74), NDUFB2 (0.73), UQCRQ (0.73), S100A11 (0.72), COX6C (0.72), S100A13 (0.71), PPIA (0.71), SRI (0.70), NDUFC1 (0.70), TPI1 (0.68), GUK1 (0.68), S100A10 (0.68)
- CLDN4neg IFN ME7: HLA-B (0.83), HLA-C (0.81), HLA-A (0.80), CD74 (0.76), HLA-DRB1 (0.75), B2M (0.74), HLA-DRA (0.72), HLA-DPA1 (0.72), CTSS (0.70), HLA-E (0.70), HLA-DQB1 (0.69), HLA-DPB1 (0.67)

## Patient-level module vs T/NK

Every non-grey module is listed. Role is TJ or IFN only when that family's enrichment q is below 0.05. Cohort rows and the full table are in `results/tables/module_immune_correlation.tsv`. q is BH across non-grey modules inside that network for the T/NK meta-analysis p. Partial ρ is Spearman of the module score vs T/NK given malignant CLDN4 %pos.

| network | module | role | N | ρ vs T/NK (95% CI) | p | q | I² | ρ vs T/NK \| CLDN4 | p | I² partial | ρ vs CLDN4 %pos | p |
|---|---|---|---:|---|---:|---:|---:|---|---:|---:|---|---:|
| CLDN4pos | ME1 | other | 65 | 0.247 (-0.172 to 0.590) | 0.245 | 0.483 | 56.6% | 0.208 | 0.352 | 56.3% | -0.113 | 0.410 |
| CLDN4pos | ME2 | other | 65 | 0.217 (-0.048 to 0.454) | 0.108 | 0.483 | 0.0% | 0.171 | 0.228 | 0.0% | -0.143 | 0.296 |
| CLDN4pos | ME3 | IFN | 65 | 0.298 (-0.074 to 0.598) | 0.114 | 0.483 | 46.5% | 0.439 | 0.000964 | 0.0% | 0.072 | 0.748 |
| CLDN4pos | ME4 | TJ | 65 | 0.123 (-0.316 to 0.519) | 0.591 | 0.694 | 61.3% | 0.265 | 0.171 | 43.5% | 0.185 | 0.172 |
| CLDN4pos | ME5 | other | 65 | 0.130 (-0.365 to 0.569) | 0.617 | 0.694 | 70.2% | 0.130 | 0.605 | 65.0% | 0.040 | 0.777 |
| CLDN4pos | ME6 | other | 65 | 0.060 (-0.419 to 0.513) | 0.816 | 0.816 | 69.3% | 0.154 | 0.490 | 55.7% | 0.189 | 0.443 |
| CLDN4pos | ME7 | other | 65 | 0.107 (-0.220 to 0.412) | 0.524 | 0.694 | 29.9% | 0.260 | 0.062 | 0.0% | 0.199 | 0.342 |
| CLDN4pos | ME8 | other | 65 | 0.249 (-0.156 to 0.582) | 0.226 | 0.483 | 53.6% | 0.161 | 0.268 | 3.9% | -0.086 | 0.543 |
| CLDN4pos | ME9 | other | 65 | 0.269 (-0.209 to 0.643) | 0.268 | 0.483 | 66.9% | 0.225 | 0.178 | 25.2% | -0.110 | 0.582 |
| CLDN4neg | ME1 | other | 65 | 0.229 (-0.036 to 0.464) | 0.090 | 0.235 | 0.0% | 0.247 | 0.146 | 27.8% | -0.036 | 0.858 |
| CLDN4neg | ME2 | other | 65 | 0.322 (0.022 to 0.569) | 0.036 | 0.235 | 22.0% | 0.270 | 0.053 | 0.0% | -0.148 | 0.277 |
| CLDN4neg | ME3 | TJ-best-overlap-not-significant | 65 | 0.214 (-0.179 to 0.548) | 0.286 | 0.400 | 50.6% | 0.335 | 0.019 | 6.1% | 0.117 | 0.426 |
| CLDN4neg | ME4 | other | 65 | 0.064 (-0.484 to 0.576) | 0.833 | 0.833 | 77.6% | 0.079 | 0.805 | 77.9% | 0.065 | 0.636 |
| CLDN4neg | ME5 | other | 65 | 0.146 (-0.362 to 0.588) | 0.584 | 0.681 | 71.6% | 0.247 | 0.359 | 70.1% | 0.188 | 0.315 |
| CLDN4neg | ME6 | other | 65 | 0.311 (-0.242 to 0.711) | 0.267 | 0.400 | 75.6% | 0.272 | 0.227 | 57.9% | -0.110 | 0.494 |
| CLDN4neg | ME7 | IFN | 65 | 0.248 (-0.049 to 0.505) | 0.101 | 0.235 | 17.9% | 0.329 | 0.017 | 0.0% | 0.036 | 0.828 |

Cohort Spearmans for the modules named in the identification table:

| network | module | cohort | n | ρ vs T/NK | p | partial ρ \| CLDN4 | p |
|---|---|---|---:|---:|---:|---:|---:|
| CLDN4pos | ME3 | GSE123902 | 13 | 0.066 | 0.831 | 0.370 | 0.263 |
| CLDN4pos | ME3 | GSE131907 | 21 | 0.022 | 0.924 | 0.354 | 0.137 |
| CLDN4pos | ME3 | GSE205335 | 22 | 0.630 | 0.002 | 0.616 | 0.004 |
| CLDN4pos | ME3 | GSE189357 | 9 | 0.333 | 0.381 | 0.075 | 0.872 |
| CLDN4pos | ME4 | GSE123902 | 13 | -0.286 | 0.344 | -0.121 | 0.724 |
| CLDN4pos | ME4 | GSE131907 | 21 | -0.188 | 0.414 | 0.016 | 0.949 |
| CLDN4pos | ME4 | GSE205335 | 22 | 0.504 | 0.017 | 0.539 | 0.014 |
| CLDN4pos | ME4 | GSE189357 | 9 | 0.433 | 0.244 | 0.593 | 0.161 |
| CLDN4neg | ME3 | GSE123902 | 13 | 0.341 | 0.255 | 0.434 | 0.182 |
| CLDN4neg | ME3 | GSE131907 | 21 | -0.148 | 0.522 | 0.099 | 0.686 |
| CLDN4neg | ME3 | GSE205335 | 22 | 0.549 | 0.008 | 0.557 | 0.011 |
| CLDN4neg | ME3 | GSE189357 | 9 | -0.050 | 0.898 | 8.7e-18 | 1.000 |
| CLDN4neg | ME7 | GSE123902 | 13 | -0.011 | 0.972 | 0.183 | 0.591 |
| CLDN4neg | ME7 | GSE131907 | 21 | 0.026 | 0.911 | 0.238 | 0.327 |
| CLDN4neg | ME7 | GSE205335 | 22 | 0.512 | 0.015 | 0.508 | 0.022 |
| CLDN4neg | ME7 | GSE189357 | 9 | 0.400 | 0.286 | 0.170 | 0.716 |

## CLDN4-positive vs CLDN4-negative cells, same patients

For each patient with at least 20 malignant cells in both strata, the module score was computed separately in each stratum. The table is the within-patient difference (CLDN4-positive minus CLDN4-negative). This is not the T/NK test.

| network | module | family | cohort | n | median Δ | n with pos>neg | Wilcoxon p |
|---|---|---|---|---:|---:|---:|---:|
| CLDN4pos | ME4 | TJ | GSE123902 | 11 | 0.065 | 7 | 0.175 |
| CLDN4pos | ME4 | TJ | GSE131907 | 18 | 0.150 | 16 | 0.000328 |
| CLDN4pos | ME4 | TJ | GSE205335 | 21 | 0.173 | 16 | 0.002 |
| CLDN4pos | ME4 | TJ | GSE189357 | 9 | -0.086 | 3 | 0.129 |
| CLDN4pos | ME3 | IFN | GSE123902 | 11 | 0.115 | 9 | 0.147 |
| CLDN4pos | ME3 | IFN | GSE131907 | 18 | 0.309 | 16 | 5.34e-05 |
| CLDN4pos | ME3 | IFN | GSE205335 | 21 | 0.145 | 17 | 0.000607 |
| CLDN4pos | ME3 | IFN | GSE189357 | 9 | 0.341 | 9 | 0.004 |
| CLDN4neg | ME3 | TJ-overlap-ns | GSE123902 | 11 | 0.251 | 10 | 0.005 |
| CLDN4neg | ME3 | TJ-overlap-ns | GSE131907 | 18 | 0.132 | 17 | 0.000328 |
| CLDN4neg | ME3 | TJ-overlap-ns | GSE205335 | 21 | 0.118 | 17 | 0.000852 |
| CLDN4neg | ME3 | TJ-overlap-ns | GSE189357 | 9 | 0.375 | 9 | 0.004 |
| CLDN4neg | ME7 | IFN | GSE123902 | 11 | -0.042 | 5 | 0.765 |
| CLDN4neg | ME7 | IFN | GSE131907 | 18 | 0.209 | 16 | 0.000252 |
| CLDN4neg | ME7 | IFN | GSE205335 | 21 | 0.085 | 15 | 0.272 |
| CLDN4neg | ME7 | IFN | GSE189357 | 9 | 0.236 | 8 | 0.074 |

## Check against the locked malignant IFN and TJ scores

Spearman of the new module score vs the locked pseudobulk family score on the same patients (DL across cohorts). This asks whether the called module tracks the previous IFN or TJ summary. It is not an immune test.

| network | module | family | ρ vs locked score | p | N |
|---|---|---|---:|---:|---:|
| CLDN4pos | ME4 | TJ | 0.428 | 0.00087 | 65 |
| CLDN4pos | ME3 | IFN | 0.601 | 4.3e-07 | 65 |
| CLDN4neg | ME3 | TJ | 0.352 | 0.007 | 65 |
| CLDN4neg | ME7 | IFN | 0.373 | 0.004 | 65 |

## Reading

The two calls are allowed to land on the same module. They do not. In both strata the significant IFN module is an MHC-I / immunoproteasome block, not the TJ block.

Junction genes in the network, and the module that contains each one (kME for that module). CLDN4 is held out of the TJ enrichment set, so it cannot by itself create the TJ call.

| network | gene | module | kME |
|---|---|---|---:|
| CLDN4pos | TJP1 | ME4 | 0.69 |
| CLDN4pos | TJP2 | not in network |  |
| CLDN4pos | TJP3 | not in network |  |
| CLDN4pos | OCLN | ME4 | 0.57 |
| CLDN4pos | F11R | ME4 | 0.64 |
| CLDN4pos | CDH1 | ME4 | 0.64 |
| CLDN4pos | CLDN1 | not in network |  |
| CLDN4pos | CLDN3 | ME1 | 0.56 |
| CLDN4pos | CLDN4 | ME1 | 0.24 |
| CLDN4pos | CLDN7 | ME7 | 0.55 |
| CLDN4pos | TACSTD2 | ME4 | 0.48 |
| CLDN4pos | EPCAM | ME1 | 0.56 |
| CLDN4pos | CRB3 | ME7 | 0.33 |
| CLDN4neg | TJP1 | ME1 | 0.51 |
| CLDN4neg | TJP2 | not in network |  |
| CLDN4neg | TJP3 | not in network |  |
| CLDN4neg | OCLN | ME1 | 0.34 |
| CLDN4neg | F11R | ME5 | 0.51 |
| CLDN4neg | CDH1 | ME5 | 0.71 |
| CLDN4neg | CLDN1 | not in network |  |
| CLDN4neg | CLDN3 | ME5 | 0.35 |
| CLDN4neg | CLDN4 | not in network |  |
| CLDN4neg | CLDN7 | ME3 | 0.49 |
| CLDN4neg | TACSTD2 | ME5 | 0.65 |
| CLDN4neg | EPCAM | ME3 | 0.33 |
| CLDN4neg | CRB3 | ME3 | 0.45 |

CLDN4pos TJ set inside ME4 (called): TJP1 (0.69), CTTN (0.66), PLEC (0.64), F11R (0.64), CDH1 (0.64), LLGL2 (0.62), MYH9 (0.62), ACTR2 (0.61), OCLN (0.57), ITGB1 (0.56), RUNX1 (0.53), EZR (0.51), CCND1 (0.40), VASP (0.36).

CLDN4neg TJ set inside ME3 (best overlap, not called): RAC1 (0.67), ARPC1A (0.67), MYL12A (0.64), ARPC1B (0.64), ACTB (0.63), MYL12B (0.61), ARPC3 (0.59), RHOA (0.54), ACTG1 (0.52), CLDN7 (0.49), RAB13 (0.46), CRB3 (0.45), TUBA1C (0.41), ARL2 (0.40), ARPC2 (0.33).

Of the locked sets, 61 TJ genes and 63 IFN genes entered the 2500-gene network. Enrichment is against that background.

Patient-level T/NK tests use the full malignant pseudobulk. A 95% CI that crosses 0 is not a concordant module–immune correlation. The locked result is malignant CLDN4 percent-positive versus T/NK fraction. It is not replaced by a module score. The within-patient difference (CLDN4-positive cells minus CLDN4-negative cells) is a different contrast from that between-patient result, and from the locked Q4-versus-Q1 malignant IFN comparison.

## What this does not say

Metacells are not extra patients. A metacell-level p-value is not reported. The module was discovered and scored in the same 65 units; the gene list was not chosen using T/NK, but it is not an external locked signature. Dataset z-scoring removes mean shifts. It is not Harmony and it is not a batch-free claim. Visium co-localization is not in this analysis. Public mouse KL matrices were not added.

Figures: `results/figures/fig_module_enrichment.png`, `fig_module_vs_tnk.png`, `fig_forest_tnk.png`.


# DepMap lung TACSTD2 vs interferon analysis

## Bottom line

- TACSTD2 RNA was associated with **higher IFN-alpha score** (rho=0.515, 95% bootstrap CI 0.400 to 0.620, BH q=3.66e-15, n=208).
- TACSTD2 RNA was associated with **higher IFN-gamma score** (rho=0.490, 95% bootstrap CI 0.368 to 0.593, BH q=5.55e-14, n=208).
- In the candidate CRISPR screen, 4 of 222 IFN-gene associations with TACSTD2 RNA had q<0.05; 1 remained at q<0.05 after model-type adjustment. This is an association screen, not evidence that TACSTD2 controls those genes.
- No IFN-gene effect correlated with TACSTD2 gene effect at q<0.05 (0/222).
- The pooled adjusted candidate hit was not assumed to be histology-general; exploratory subtype estimates are reported below.
- The analysis is observational and cross-sectional. Histology, lineage state, culture conditions, and screen quality can create correlations. No causal claim is supported.

## Data and cohort

- DepMap Public 24Q4 (stable Figshare archive; release DOI [10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1)).
- Primary cohort: `OncotreeLineage == Lung`, `ModelType == Cell Line`, excluding `OncotreePrimaryDisease == Non-Cancerous`.
- Metadata lung cancer cell lines: 254; RNA models analyzed: 208; CRISPR-overlap models: 123.
- RNA values are DepMap log2(TPM+1). CRISPR values are Chronos gene effects; more negative values indicate stronger dependency.

## Prespecified expression tests

| signature   |   n |   rho |   ci_low |   ci_high |        p |        q |
|:------------|----:|------:|---------:|----------:|---------:|---------:|
| ifna_z      | 208 | 0.515 |    0.4   |     0.62  | 1.83e-15 | 3.66e-15 |
| ifng_z      | 208 | 0.49  |    0.368 |     0.593 | 5.55e-14 | 5.55e-14 |

Each Hallmark score is the mean of gene-wise z-scores calculated within the lung cohort. This measures relative IFN-like transcriptional state, not IFN protein or pathway activation. The alpha and gamma sets overlap by 73 genes, so the tests are not independent.

## TACSTD2 dependency vs IFN scores

| signature   |   n |    rho |   ci_low |   ci_high |     p |     q |
|:------------|----:|-------:|---------:|----------:|------:|------:|
| ifna_z      | 123 | -0.043 |   -0.225 |     0.146 | 0.633 | 0.752 |
| ifng_z      | 123 | -0.029 |   -0.216 |     0.163 | 0.752 | 0.752 |

A positive rho here means higher IFN score accompanies a less-negative (weaker) TACSTD2 dependency; a negative rho means stronger TACSTD2 dependency.

## Histology sensitivity analyses

Partial Spearman correlations residualize ranked variables on DepMap model type (rare types with <10 RNA-profiled models pooled as `Other`).

| signature   |   n |   rho |        p |        q |
|:------------|----:|------:|---------:|---------:|
| ifna_z      | 208 | 0.371 | 3.45e-08 | 6.89e-08 |
| ifng_z      | 208 | 0.334 | 8.33e-07 | 8.33e-07 |

Unadjusted major-subtype estimates:

| cohort   | signature   |   n |   rho |     p |     q |
|:---------|:------------|----:|------:|------:|------:|
| LUAD     | ifna_z      |  80 | 0.345 | 0.002 | 0.009 |
| LUAD     | ifng_z      |  80 | 0.259 | 0.02  | 0.031 |
| SCLC     | ifna_z      |  59 | 0.379 | 0.003 | 0.009 |
| SCLC     | ifng_z      |  59 | 0.351 | 0.006 | 0.013 |
| LUSC     | ifna_z      |  27 | 0.291 | 0.141 | 0.169 |
| LUSC     | ifng_z      |  27 | 0.242 | 0.223 | 0.223 |

These are sensitivity analyses, not extra discovery tests. Differences between subtypes may reflect small samples; no formal interaction test was prespecified.

## IFN-gene CRISPR association screen

For every available unique gene in the union of the two Hallmark sets, its Chronos gene effect was correlated with TACSTD2 RNA. BH correction is across all candidate genes for that predictor. Positive rho means high-TACSTD2 lines are less dependent on the gene; negative rho means they are more dependent. Adjusted columns are a sensitivity analysis residualizing ranks on model type.

| gene     | sets      |   n |    rho |     q |   rho_adjusted |   q_adjusted |
|:---------|:----------|----:|-------:|------:|---------------:|-------------:|
| STAT3    | IFNG      | 123 | -0.389 | 0.002 |         -0.357 |        0.011 |
| ITGB7    | IFNG      | 123 | -0.311 | 0.048 |         -0.185 |        0.464 |
| SLAMF7   | IFNG      | 123 |  0.303 | 0.048 |          0.247 |        0.264 |
| PROCR    | IFNA      | 123 | -0.297 | 0.048 |         -0.212 |        0.32  |
| AUTS2    | IFNG      | 123 |  0.27  | 0.097 |          0.254 |        0.264 |
| NFKBIA   | IFNG      | 123 | -0.269 | 0.097 |         -0.251 |        0.264 |
| LGALS3BP | IFNA+IFNG | 123 |  0.26  | 0.108 |          0.218 |        0.316 |
| RNF213   | IFNG      | 123 |  0.258 | 0.108 |          0.121 |        0.761 |
| WARS1    | IFNA+IFNG | 123 |  0.249 | 0.133 |          0.217 |        0.316 |
| PNPT1    | IFNA+IFNG | 123 |  0.245 | 0.142 |          0.199 |        0.377 |
| PTGS2    | IFNG      | 123 |  0.239 | 0.157 |          0.239 |        0.284 |
| DDX60    | IFNA+IFNG | 123 | -0.223 | 0.245 |         -0.187 |        0.464 |
| IL6      | IFNG      | 123 |  0.219 | 0.252 |          0.171 |        0.537 |
| RIPK1    | IFNG      | 123 |  0.217 | 0.254 |          0.18  |        0.464 |
| PTPN1    | IFNG      | 123 | -0.209 | 0.299 |         -0.075 |        0.892 |

Because pooled adjustment does not establish consistency across histologies, the following post-screen estimates show each adjusted q<0.05 hit within major subtypes. These estimates are exploratory and their p-values are intentionally not presented as confirmatory tests:

| gene   | cohort   |   n |    rho |
|:-------|:---------|----:|-------:|
| STAT3  | LUAD     |  51 | -0.604 |
| STAT3  | SCLC     |  25 |  0.002 |
| STAT3  | LUSC     |  20 |  0.385 |

The full table also contains a secondary screen using TACSTD2 gene effect as the predictor, corrected as a separate family.

## Robustness and limitations

- Rank-average versions of both signatures are provided per model for score sensitivity, but the z-score definition is primary.
- Candidate genes were fixed from MSigDB Hallmark IFN-alpha (97 genes) and IFN-gamma (200 genes); they were not selected from these data.
- Hallmark sets include broad antigen-presentation, proteasome, apoptosis, and signaling genes, not only IFN-specific effectors.
- Gene-effect co-variation can arise from screen/library effects and general fitness biology. Correlation does not establish synthetic lethality.
- Cell lines omit immune/stromal compartments and do not model drug response. The results should not be extrapolated to patients or ADC efficacy.
- No mutation, copy-number, media, growth-rate, or batch covariates were modeled. Histology adjustment is limited and cannot remove all confounding.
- Current DepMap releases are newer than 24Q4, but 24Q4 was used because it has stable public URLs and file checksums. Release-specific replication is needed.

## Reproduction

See `README.md`; all source checksums are verified before analysis. Random bootstrap seed: 20260816.

## Files

- `tables/model_scores.csv`: model-level metadata, TACSTD2, and signatures
- `tables/primary_correlations.csv`: primary expression tests
- `tables/tacstd2_dependency_correlations.csv`: TACSTD2 dependency tests
- `tables/histology_adjusted_correlations.csv`: partial Spearman sensitivity
- `tables/subgroup_correlations.csv`: major subtype sensitivity
- `tables/ifn_gene_crispr_correlations.csv`: full raw and model-type-adjusted candidate CRISPR screen
- `tables/top_candidate_subgroup_correlations.csv`: exploratory subtype estimates for adjusted candidate hits
- `figures/`: scatter and candidate-screen plots
- `qc_summary.json`: counts, coverage, checksums, and software-independent inputs

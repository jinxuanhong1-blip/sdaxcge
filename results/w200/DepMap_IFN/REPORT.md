# DepMap lung TACSTD2 vs interferon analysis

## Bottom line

- TACSTD2 RNA was associated with **higher IFN-alpha score** (rho=0.515, 95% bootstrap CI 0.400 to 0.620, BH q=3.66e-15, n=208).
- TACSTD2 RNA was associated with **higher IFN-gamma score** (rho=0.490, 95% bootstrap CI 0.368 to 0.593, BH q=5.55e-14, n=208).
- In the candidate CRISPR screen, 4 of 222 IFN-gene associations with TACSTD2 RNA had q<0.05; 1 remained at q<0.05 after model-type adjustment. This is an association screen, not evidence that TACSTD2 controls those genes.
- No IFN-gene effect correlated with TACSTD2 gene effect at q<0.05 (0/222).
- The pooled adjusted candidate hit was not assumed to be histology-general; exploratory subtype estimates are reported below.
- DepMap has no IFN-treated profiles. In public GEO lung-line IFN experiments, no BH-significant TACSTD2 induction or repression.
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

## IFN-treated lung lines (GEO; not DepMap)

DepMap Public 24Q4 has **no IFN-treated RNA-seq or microarray profiles**. PRISM also lacks recombinant IFN biologics. The treatment test therefore uses public GEO processed matrices of lung lines treated with IFN protein.

- Primary TACSTD2 contrasts with ISG QC pass: 10. no BH-significant TACSTD2 induction or repression.
- Primary contrasts that failed ISG QC (MX1/ISG15/CXCL10 log2FC all < 1): 3.
- Most usable public IFN time courses are A549. In DepMap 24Q4 A549 TACSTD2 is very low (log2(TPM+1)=0.20), so these experiments test induction from a low basal state, not regulation in TROP2-high NSCLC.
- A549-heavy GEO data cannot confirm or refute the basal DepMap TACSTD2–IFN correlation. They only ask whether IFN itself moves TACSTD2.

Primary TACSTD2 contrasts (BH across TACSTD2 primary tests):

| dataset   | cell_line   | ifn_type   | timepoint       |   n_treated |   n_control |   log2fc |     p | q     | isg_qc_pass   |
|:----------|:------------|:-----------|:----------------|------------:|------------:|---------:|------:|:------|:--------------|
| GSE109720 | EBC-1       | IFNG       | unspecified     |           3 |           3 |    1.031 | 0.013 | 0.130 | True          |
| GSE5542   | A549        | IFNA       | 24hr            |           4 |           4 |   -0.77  | 0.394 | 0.780 | True          |
| GSE178640 | A549        | IFNA2      | DMSO_pretreated |           3 |           3 |   -0.023 | 0.423 | 0.780 | True          |
| GSE215771 | A549_WT     | IFNG       | 24h             |           3 |           3 |   -0.05  | 0.423 | 0.780 | True          |
| GSE5542   | A549        | IFNG       | 6hr             |           4 |           4 |    0.547 | 0.436 | 0.780 | True          |
| GSE156295 | A549        | IFNI       | unspecified     |           2 |           2 |   -0.035 | 0.5   | 0.780 | True          |
| GSE60572  | A549        | IFNB       | 6h_IFN_2b_1000  |           3 |           3 |    0.04  | 0.546 | 0.780 | True          |
| GSE60572  | A549        | IFNB       | 6h_IFN_2b_100   |           3 |           3 |   -0.025 | 0.667 | 0.815 | True          |
| GSE5542   | A549        | IFNA       | 6hr             |           4 |           4 |    0.203 | 0.734 | 0.815 | True          |
| GSE5542   | A549        | IFNG       | 24hr            |           4 |           4 |   -0.118 | 0.872 | 0.872 | True          |
| GSE109720 | NCI-H596    | IFNG       | unspecified     |           3 |           3 |   -1.654 | 0.003 | NA    | False         |
| GSE109720 | NCI-H1573   | IFNG       | unspecified     |           3 |           3 |   -0.629 | 0.039 | NA    | False         |
| GSE109720 | NCI-H1993   | IFNG       | unspecified     |           3 |           3 |   -0.06  | 0.63  | NA    | False         |

ISG QC requires MX1 log2FC >= 1, or ISG15/CXCL10 log2FC >= 1 without MX1 collapse. A contrast with MX1 log2FC <= -1 is treated as a failed IFN response. GSE109720 sample columns were aligned to series-matrix order; failed-QC lines are not interpreted as IFN-responsive. BH q-values are only among ISG-passing primary TACSTD2 tests.

Inventory:

| dataset     | included   | reason                                                                            |
|:------------|:-----------|:----------------------------------------------------------------------------------|
| DepMap 24Q4 | False      | No IFN-treated transcriptome; basal RNA/CRISPR only                               |
| GSE5542     | True       | A549 IFNA-con1 / IFNG 6h and 24h, n=4, processed series matrix                    |
| GSE60572    | True       | A549 IFNB 100/1000 IU; only 6h has matched untreated samples                      |
| GSE156295   | True       | A549 and HTBE type I IFN vs mock, n=2, processed counts                           |
| GSE178640   | True       | A549 IFN-alpha2 after DMSO, n=3, processed Ensembl counts                         |
| GSE215771   | True       | A549 IFNG 24h WT/IRF1KO/NF2KO, n=3, processed TPM                                 |
| GSE109720   | True       | EBC-1/H1573/H1993/H596 IFNG; ISG QC applied because some lines do not induce ISGs |
| GSE241977   | False      | A549 AhR KO vs control; no IFN-versus-mock contrast in the deposited samples      |
| GSE261920   | False      | lncRNA-focused A549 IFN series; not used as a coding-gene TACSTD2 test            |

## Robustness and limitations

- Rank-average versions of both signatures are provided per model for score sensitivity, but the z-score definition is primary.
- Candidate genes were fixed from MSigDB Hallmark IFN-alpha (97 genes) and IFN-gamma (200 genes); they were not selected from these data.
- Hallmark sets include broad antigen-presentation, proteasome, apoptosis, and signaling genes, not only IFN-specific effectors.
- Gene-effect co-variation can arise from screen/library effects and general fitness biology. Correlation does not establish synthetic lethality.
- Cell lines omit immune/stromal compartments and do not model drug response. The results should not be extrapolated to patients or ADC efficacy.
- No mutation, copy-number, media, growth-rate, or batch covariates were modeled. Histology adjustment is limited and cannot remove all confounding.
- Current DepMap releases are newer than 24Q4, but 24Q4 was used because it has stable public URLs and file checksums. Release-specific replication is needed.
- DepMap has no IFN-treated profiles. GEO IFN experiments are mostly A549, which is TACSTD2-low in DepMap, and cannot explain the basal cross-line correlation.

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
- `tables/ifn_treated_contrasts.csv`: GEO IFN-versus-control TACSTD2 and ISG tests
- `tables/ifn_treated_inventory.csv`: included and excluded treatment datasets
- `figures/`: scatter, candidate-screen, and IFN-treatment plots
- `qc_summary.json`: counts, coverage, checksums, and software-independent inputs

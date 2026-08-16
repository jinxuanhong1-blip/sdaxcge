# 04 · TACSTD2 / CLDN4 versus IFN, MHC-I and immune genes

Release DepMap Public 24Q4. Lung lines with expression: **208** (NSCLC 143, SCLC 59). Lung lines with both CRISPR and expression: **123** (NSCLC 95).

MHC-I composite uses 16 genes; IFN-response composite uses 14 genes. Missing panel genes are listed in `immune_panel_inventory.csv`.

## Composite scores (headline)

| test | query | partner | cohort | n | r | 95% CI | p | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| expr_vs_expr | TACSTD2 | MHC-I score | lung | 208 | 0.390 | [0.268, 0.500] | 5.66e-09 | 1.91e-08 |
| expr_vs_expr | CLDN4 | MHC-I score | lung | 208 | 0.197 | [0.063, 0.325] | 0.00429 | 0.00947 |
| expr_vs_expr | TACSTD2 | IFN-response score | lung | 208 | 0.489 | [0.378, 0.586] | 6.54e-14 | 5.23e-13 |
| expr_vs_expr | CLDN4 | IFN-response score | lung | 208 | 0.246 | [0.114, 0.370] | 0.000341 | 0.00103 |
| expr_vs_expr | TACSTD2 | MHC-I score | NSCLC | 143 | 0.145 | [-0.020, 0.302] | 0.084 | 0.128 |
| expr_vs_expr | CLDN4 | MHC-I score | NSCLC | 143 | 0.093 | [-0.072, 0.253] | 0.27 | 0.382 |
| expr_vs_expr_partial_NSCLC | TACSTD2 | MHC-I score | lung_partial_NSCLC | 208 | 0.180 | n/a (partial) | 0.00929 | 0.0154 |
| expr_vs_expr_partial_NSCLC | CLDN4 | MHC-I score | lung_partial_NSCLC | 208 | 0.107 | n/a (partial) | 0.126 | 0.201 |
| geneEffect_vs_expr | TACSTD2 | MHC-I score | lung | 123 | -0.078 | [-0.252, 0.100] | 0.389 | 0.867 |
| geneEffect_vs_expr | CLDN4 | MHC-I score | lung | 123 | -0.000 | [-0.177, 0.177] | 0.999 | 0.999 |
| geneEffect_vs_expr | TACSTD2 | IFN-response score | lung | 123 | -0.052 | [-0.227, 0.126] | 0.566 | 0.929 |
| geneEffect_vs_expr | CLDN4 | IFN-response score | lung | 123 | 0.129 | [-0.049, 0.300] | 0.154 | 0.979 |
| geneEffect_vs_expr | TACSTD2 | MHC-I score | NSCLC | 95 | 0.021 | [-0.182, 0.221] | 0.844 | 0.98 |
| geneEffect_vs_expr | CLDN4 | MHC-I score | NSCLC | 95 | -0.048 | [-0.247, 0.155] | 0.643 | 0.931 |
| geneEffect_vs_expr_partial_NSCLC | TACSTD2 | MHC-I score | lung_partial_NSCLC | 123 | -0.024 | n/a (partial) | 0.796 | 0.981 |
| geneEffect_vs_expr_partial_NSCLC | CLDN4 | MHC-I score | lung_partial_NSCLC | 123 | 0.040 | n/a (partial) | 0.659 | 0.917 |

Full per-gene tables: `immune_expr_vs_expr.csv`, `immune_geneEffect_vs_expr.csv`, `immune_geneEffect_vs_geneEffect.csv`, `immune_headline.csv`. Figure: `fig3_immune.png`.

A negative r in test B means higher immune-gene expression accompanies a more negative (more dependent) Chronos score. Given that neither TACSTD2 nor CLDN4 is a CRISPR dependency in these lines, large effects in B are not expected.

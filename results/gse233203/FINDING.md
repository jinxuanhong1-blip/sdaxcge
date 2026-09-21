# GSE233203 effusion CLDN4 vs T/NK

Not merged into concordant-4. Not solid tumor.

**Series:** [GSE233203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233203). Seven stage-IV lung adenocarcinoma pleural fluids after atezolizumab + bevacizumab + carboplatin/paclitaxel. GEO `therapeutic response`: 3 Response, 4 Non-response. Open 10x matrices.

Malignant = marker epithelium with normal-lung score at or below the epithelial 75th percentile in that sample. In four samples that cutoff is 0, so malignant equals epithelium. T/NK = marker T or NK. QC cells: nGene ≥200 and mitochondrial percent <20.

## CLDN4 vs T/NK

| Test | n | Result | p |
|---|---:|---|---:|
| Malignant CLDN4 mean vs T/NK fraction | 7 | ρ=+0.36 | 0.43 |
| Malignant CLDN4 percent vs T/NK fraction | 7 | ρ=+0.21 | 0.64 |
| CLDN4 mean, response vs non-response | 3 vs 4 | Response 0.62 vs non-response 0.058 (response higher) | exact 0.057 |
| CLDN4 percent, response vs non-response | 3 vs 4 | Response 0.41 vs non-response 0.066 | exact 0.11 |

Leave-one-out of the mean correlation stays positive (ρ +0.14 to +0.66, all p≥0.16). There is no inverse to be fragile.

Two non-response samples have almost no T/NK (NCCLu_397 n=2, NCCLu_185 n=25), so the T/NK fraction is unstable. CLDN4 itself is low in five of seven samples (mean log1p(CP10k) 0.02–0.13). The two higher samples are responders (NCCLu_376 1.20, NCCLu_383 0.51). That is the opposite of a resistance-up pattern, and 3 vs 4 cannot go below exact p=0.029. p=0.057 is not a response claim.

Figure: `fig_cldn4_vs_tnk.png`. Table: `per_sample.tsv`.

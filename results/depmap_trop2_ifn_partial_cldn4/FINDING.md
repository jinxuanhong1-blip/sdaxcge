# DepMap/CCLE: does CLDN4 attenuate TROP2 vs IFN / MHC / immune scores?

Numbers below are written by `analyze.py` from `tables/key_stats.json` and `tables/partials.csv`. They are not typed by hand.

## Question

TROP2 (TACSTD2) and CLDN4 protein are already correlated in Gygi CCLE lung columns. This page tests whether the association of **TROP2 RNA or protein** with pre-specified IFN, MHC-I, and immune-ligand scores shrinks after adjustment for CLDN4, and whether CLDN4's own association shrinks after adjustment for TROP2. CRISPR Chronos is reported for the same scores when the screen exists.

Partial correlation is the Pearson correlation of rank residuals after linear adjustment for the covariate ranks. It is an association contrast, not a mediation or causal estimate. Cultured lines have no T-cell infiltrate and no IFN treatment.

## Answer

In NSCLC RNA, CLDN4 does not attenuate the TROP2 associations that clear p < 0.05. Hallmark IFN-γ is +0.339 (n=143, p=3.38e-05, q=0.0001) unadjusted and +0.251 (n=143, p=0.0026, q=0.0079, |partial|/|unadj|=0.74; retained) after CLDN4. The immune-ligand score is +0.405 (n=143, p=5.09e-07, q=3.05e-06) unadjusted and +0.347 (n=143, p=2.31e-05, q=0.0001, |partial|/|unadj|=0.86; retained) after CLDN4. MHC-I has no unadjusted TROP2 association (+0.154 (n=143, p=0.0671, q=0.0671)), so there is nothing there for CLDN4 to attenuate.

On the same NSCLC RNA rows, CLDN4's own positive associations are the ones that shrink. CLDN4 vs Hallmark IFN-γ given TROP2 is +0.015 (n=143, p=0.8637, q=0.8637, |partial|/|unadj|=0.06; attenuated). CLDN4 vs the immune-ligand score given TROP2 is -0.059 (n=143, p=0.4871, q=0.5898, |partial|/|unadj|=0.25; attenuated).

All-lung RNA (n includes SCLC/NET) agrees that TROP2's IFN-γ and immune-ligand associations are retained after CLDN4. All-lung TROP2 vs MHC-I is retained after CLDN4 alone, and is attenuated once a binary NSCLC indicator is added (+0.113 (n=214, p=0.1000, |partial|/|unadj|=0.36; attenuated)). CLDN4 vs the immune-ligand score reverses sign after TROP2 adjustment: -0.203 (n=214, p=0.0029, q=0.0044, |partial|/|unadj|=1.02; sign_reversed). Adding a binary NSCLC indicator leaves that residual -0.188 (n=214, p=0.0062, |partial|/|unadj|=0.94; sign_reversed).

LUAD-only RNA is exploratory and is outside the FDR family. TROP2 vs Hallmark IFN-γ is +0.268 (n=80, p=0.0161) unadjusted and +0.126 (n=80, p=0.2678, |partial|/|unadj|=0.47; attenuated) after CLDN4. TROP2 vs the immune-ligand score after CLDN4 is +0.254 (n=80, p=0.0236, |partial|/|unadj|=0.72; retained). The IFN-γ retention above is an NSCLC result, not a LUAD-only result.

Protein partials do not decide the question. Primary protein partials with BH q < 0.05: 0. All-lung TROP2 vs Hallmark IFN-γ protein is +0.408 (n=44, p=0.0060, q=0.0359) unadjusted and +0.225 (n=44, p=0.1477, q=0.5000, |partial|/|unadj|=0.55; ns_after_adjustment_point_estimate_still_half_or_more) after CLDN4. The protein immune-ligand score uses only CD274, HLA-E (CXCL9/CXCL10/CXCL11 are absent). Its unadjusted TROP2 association is +0.356 (n=38, p=0.0283, q=0.0424); the bootstrap interval is in the table.

CRISPR is available. Lung lines with Chronos: n=126. TACSTD2 mean gene effect -0.0287, dependent 0/126. CLDN4 mean +0.0342, dependent 0/126. On the same lines, KRAS is dependent in 69/126 and EEF2 in 125/126. TROP2 Chronos is not associated with the three RNA scores, before or after CLDN4 Chronos, so there is no CRISPR association for CLDN4 to attenuate.

## Already-known covariate correlation (recomputed)

Gygi TenPx columns with `_LUNG_` in the name, TACSTD2 and CLDN4 both quantified, replicates not collapsed: Spearman **0.693**, n=45, p=1.31e-07. Rounds to 0.69: **True**. This is the same column definition as the existing Gygi lung complete-case result. It is a matrix check, not a new claim. Partials below use **one row per ModelID** (TenPx replicates averaged), so their n is a cell-line count, not this column count.

On the rows used for scoring: RNA lung cell lines TACSTD2–CLDN4 Spearman 0.607 (n=214, p=5.64e-23). Protein, Model.csv lung cell lines, both proteins present: 0.704 (n=44, p=9.81e-08).

## Scores

- **IFNG**: mean of within-cohort z-scores of MSigDB 2024.1 `HALLMARK_INTERFERON_GAMMA_RESPONSE` genes present in the matrix.
- **MHC1**: HLA-A, HLA-B, HLA-C, B2M. RNA requires all 4. Protein requires at least 2 quantified members.
- **IMMUNE**: immune-ligand score, mean z of CD274, PDCD1LG2, CXCL9, CXCL10, CXCL11, HLA-E, IDO1 (genes present). This is cancer-cell RNA or protein, not an immune-infiltrate score.
- CD274 alone and a compact ISG cassette are reported as components or sensitivity and are outside the primary FDR family.

Attenuation label, applied only when the unadjusted Spearman p < 0.05: **attenuated** if the partial p ≥ 0.05 and |partial| ≤ half |unadjusted|; **retained** if the partial p < 0.05 and |partial| > half |unadjusted|; **reduced_but_still_associated** if the partial p < 0.05 but |partial| ≤ half |unadjusted|; **ns_after_adjustment_point_estimate_still_half_or_more** if significance is lost but the point estimate has not halved; **sign_reversed** if the partial stays significant (p < 0.05) with the opposite sign. Retention requires the same sign. BH q for partials is within the six primary partials of that cohort (TROP2|CLDN4 and CLDN4|TROP2 × three scores). Unadjusted q is within the six matching Spearman tests.

## RNA (DepMap 24Q4 log2(TPM+1))

Lung cell lines (`OncotreeLineage == Lung`, `ModelType == Cell Line`) with TACSTD2 and CLDN4: **n = 214** (NSCLC 143, LUAD 80, LUSC 27, SCLC/NET 60). `tables/rna_lung_lines.csv` stores scores z-scored on all lung lines. NSCLC and LUAD rows in `partials.csv` recompute those z-scores inside the cohort.

| cohort | score | TROP2 unadjusted | TROP2 given CLDN4 | CLDN4 unadjusted | CLDN4 given TROP2 |
|---|---|---|---|---|---|
| all lung | Hallmark IFN-γ | +0.482 (n=214, p=7.23e-14, q=2.17e-13); boot CI [+0.36, +0.58] | +0.407 (n=214, p=6.74e-10, q=2.02e-09); boot CI [+0.28, +0.52]; |partial|/|unadj|=0.84; retained | +0.284 (n=214, p=2.50e-05, q=3.75e-05); boot CI [+0.16, +0.41] | -0.013 (n=214, p=0.8485, q=0.8485); boot CI [-0.14, +0.13]; |partial|/|unadj|=0.05; attenuated |
| all lung | MHC-I (HLA-A/B/C+B2M) | +0.314 (n=214, p=2.88e-06, q=5.77e-06); boot CI [+0.18, +0.44] | +0.234 (n=214, p=0.0006, q=0.0011); boot CI [+0.10, +0.36]; |partial|/|unadj|=0.75; retained | +0.217 (n=214, p=0.0014, q=0.0017); boot CI [+0.08, +0.35] | +0.035 (n=214, p=0.6077, q=0.7293); boot CI [-0.10, +0.17]; |partial|/|unadj|=0.16; attenuated |
| all lung | Immune-ligand score | +0.550 (n=214, p=2.71e-18, q=1.63e-17); boot CI [+0.44, +0.64] | +0.551 (n=214, p=2.76e-18, q=1.66e-17); boot CI [+0.45, +0.64]; |partial|/|unadj|=1.00; retained | +0.199 (n=214, p=0.0034, q=0.0034); boot CI [+0.07, +0.32] | -0.203 (n=214, p=0.0029, q=0.0044); boot CI [-0.33, -0.07]; |partial|/|unadj|=1.02; sign_reversed |
| NSCLC | Hallmark IFN-γ | +0.339 (n=143, p=3.38e-05, q=0.0001); boot CI [+0.18, +0.48] | +0.251 (n=143, p=0.0026, q=0.0079); boot CI [+0.08, +0.41]; |partial|/|unadj|=0.74; retained | +0.237 (n=143, p=0.0044, q=0.0084); boot CI [+0.07, +0.39] | +0.015 (n=143, p=0.8637, q=0.8637); boot CI [-0.14, +0.18]; |partial|/|unadj|=0.06; attenuated |
| NSCLC | MHC-I (HLA-A/B/C+B2M) | +0.154 (n=143, p=0.0671, q=0.0671); boot CI [-0.02, +0.32] | +0.058 (n=143, p=0.4915, q=0.5898); boot CI [-0.10, +0.22]; |partial|/|unadj|=0.38; no_unadjusted_association | +0.166 (n=143, p=0.0475, q=0.0570); boot CI [-0.01, +0.33] | +0.086 (n=143, p=0.3067, q=0.5898); boot CI [-0.08, +0.25]; |partial|/|unadj|=0.52; ns_after_adjustment_point_estimate_still_half_or_more |
| NSCLC | Immune-ligand score | +0.405 (n=143, p=5.09e-07, q=3.05e-06); boot CI [+0.26, +0.53] | +0.347 (n=143, p=2.31e-05, q=0.0001); boot CI [+0.19, +0.49]; |partial|/|unadj|=0.86; retained | +0.231 (n=143, p=0.0056, q=0.0084); boot CI [+0.06, +0.39] | -0.059 (n=143, p=0.4871, q=0.5898); boot CI [-0.23, +0.11]; |partial|/|unadj|=0.25; attenuated |
| LUAD (exploratory) | Hallmark IFN-γ | +0.268 (n=80, p=0.0161) | +0.126 (n=80, p=0.2678); |partial|/|unadj|=0.47; attenuated | +0.249 (n=80, p=0.0258) | +0.073 (n=80, p=0.5205); |partial|/|unadj|=0.29; attenuated |
| LUAD (exploratory) | MHC-I (HLA-A/B/C+B2M) | +0.078 (n=80, p=0.4928) | -0.007 (n=80, p=0.9509); |partial|/|unadj|=0.09; no_unadjusted_association | +0.109 (n=80, p=0.3361) | +0.077 (n=80, p=0.5011); |partial|/|unadj|=0.71; no_unadjusted_association |
| LUAD (exploratory) | Immune-ligand score | +0.352 (n=80, p=0.0014) | +0.254 (n=80, p=0.0236); |partial|/|unadj|=0.72; retained | +0.252 (n=80, p=0.0242) | -0.023 (n=80, p=0.8417); |partial|/|unadj|=0.09; attenuated |

All-lung sensitivity, adjusting for CLDN4 or TROP2 plus a binary NSCLC indicator (not in the six-test FDR):

- Hallmark IFN-γ: TROP2|CLDN4+NSCLC +0.295 (n=214, p=1.22e-05); boot CI [+0.16, +0.42]; |partial|/|unadj|=0.61; retained; CLDN4|TROP2+NSCLC +0.004 (n=214, p=0.9506); boot CI [-0.13, +0.14]; |partial|/|unadj|=0.02; attenuated
- MHC-I (HLA-A/B/C+B2M): TROP2|CLDN4+NSCLC +0.113 (n=214, p=0.1000); boot CI [-0.02, +0.25]; |partial|/|unadj|=0.36; attenuated; CLDN4|TROP2+NSCLC +0.056 (n=214, p=0.4176); boot CI [-0.07, +0.18]; |partial|/|unadj|=0.26; attenuated
- Immune-ligand score: TROP2|CLDN4+NSCLC +0.438 (n=214, p=2.30e-11); boot CI [+0.31, +0.54]; |partial|/|unadj|=0.80; retained; CLDN4|TROP2+NSCLC -0.188 (n=214, p=0.0062); boot CI [-0.31, -0.05]; |partial|/|unadj|=0.94; sign_reversed

CD274 and compact ISG are outside the primary FDR. NSCLC RNA:

- CD274_score: TROP2 +0.400 (n=143, p=7.22e-07); TROP2|CLDN4 +0.303 (n=143, p=0.0003); |partial|/|unadj|=0.76; retained; CLDN4|TROP2 +0.012 (n=143, p=0.8878); |partial|/|unadj|=0.04; attenuated
- ISG: TROP2 +0.275 (n=143, p=0.0009); TROP2|CLDN4 +0.228 (n=143, p=0.0063); |partial|/|unadj|=0.83; retained; CLDN4|TROP2 -0.032 (n=143, p=0.7066); |partial|/|unadj|=0.20; no_unadjusted_association

NSCLC RNA expression of immune-ligand members (log2(TPM+1)): CD274 median 2.084, frac>0 1.00; PDCD1LG2 median 0.465, frac>0 0.88; CXCL9 median 0.000, frac>0 0.30; CXCL10 median 0.057, frac>0 0.76; CXCL11 median 0.070, frac>0 0.79; HLA-E median 6.137, frac>0 1.00; IDO1 median 0.084, frac>0 0.94.

## Protein (Gygi/Nusinow CCLE TMT, one row per lung cell line)

Lung cell lines with any mapped Gygi protein: n = 76. With both TACSTD2 and CLDN4 protein: **n = 44** (NSCLC 35). RPPA antibody table from this run: mentions CLDN4 = False, mentions TACSTD2/TROP2 = False. RPPA is not used for the partials.

| cohort | score | TROP2 unadjusted | TROP2 given CLDN4 | CLDN4 unadjusted | CLDN4 given TROP2 |
|---|---|---|---|---|---|
| all lung | Hallmark IFN-γ | +0.408 (n=44, p=0.0060, q=0.0359); boot CI [+0.12, +0.66] | +0.225 (n=44, p=0.1477, q=0.5000); boot CI [-0.14, +0.53]; |partial|/|unadj|=0.55; ns_after_adjustment_point_estimate_still_half_or_more | +0.369 (n=44, p=0.0137, q=0.0412); boot CI [+0.07, +0.62] | +0.126 (n=44, p=0.4206, q=0.5000); boot CI [-0.20, +0.45]; |partial|/|unadj|=0.34; attenuated |
| all lung | MHC-I (HLA-A/B/C+B2M) | +0.260 (n=44, p=0.0877, q=0.0937); boot CI [-0.08, +0.56] | +0.117 (n=44, p=0.4546, q=0.5000); boot CI [-0.25, +0.46]; |partial|/|unadj|=0.45; no_unadjusted_association | +0.256 (n=44, p=0.0937, q=0.0937); boot CI [-0.07, +0.54] | +0.106 (n=44, p=0.5000, q=0.5000); boot CI [-0.20, +0.44]; |partial|/|unadj|=0.41; no_unadjusted_association |
| all lung | Immune-ligand score | +0.356 (n=38, p=0.0283, q=0.0424); boot CI [-0.01, +0.65] | +0.132 (n=38, p=0.4347, q=0.5000); boot CI [-0.28, +0.51]; |partial|/|unadj|=0.37; attenuated | +0.365 (n=38, p=0.0241, q=0.0424); boot CI [+0.07, +0.60] | +0.159 (n=38, p=0.3484, q=0.5000); boot CI [-0.19, +0.48]; |partial|/|unadj|=0.43; attenuated |
| NSCLC | Hallmark IFN-γ | +0.352 (n=35, p=0.0380, q=0.1141); boot CI [+0.00, +0.63] | +0.111 (n=35, p=0.5336, q=0.8004); boot CI [-0.29, +0.48]; |partial|/|unadj|=0.31; attenuated | +0.388 (n=35, p=0.0212, q=0.1141); boot CI [+0.08, +0.62] | +0.206 (n=35, p=0.2427, q=0.8004); boot CI [-0.12, +0.51]; |partial|/|unadj|=0.53; ns_after_adjustment_point_estimate_still_half_or_more |
| NSCLC | MHC-I (HLA-A/B/C+B2M) | +0.155 (n=35, p=0.3725, q=0.3725); boot CI [-0.22, +0.50] | -0.028 (n=35, p=0.8772, q=0.8772); boot CI [-0.41, +0.36]; |partial|/|unadj|=0.18; no_unadjusted_association | +0.239 (n=35, p=0.1664, q=0.2495); boot CI [-0.13, +0.54] | +0.186 (n=35, p=0.2922, q=0.8004); boot CI [-0.15, +0.51]; |partial|/|unadj|=0.78; no_unadjusted_association |
| NSCLC | Immune-ligand score | +0.229 (n=32, p=0.2079, q=0.2495); boot CI [-0.17, +0.59] | +0.077 (n=32, p=0.6801, q=0.8161); boot CI [-0.35, +0.49]; |partial|/|unadj|=0.34; no_unadjusted_association | +0.250 (n=32, p=0.1676, q=0.2495); boot CI [-0.12, +0.57] | +0.129 (n=32, p=0.4894, q=0.8004); boot CI [-0.26, +0.49]; |partial|/|unadj|=0.52; no_unadjusted_association |

Protein immune-ligand score members with enough quantified values: CD274, HLA-E. CXCL9, CXCL10, and CXCL11 are not in the Gygi matrix, so this protein score is not the RNA immune-ligand score. Hallmark IFN-γ proteins used: 147 of 200. MHC-I proteins used: HLA-A, HLA-B, HLA-C, B2M.

No primary protein partial in this table has BH q < 0.05. Protein attenuation labels use the nominal p < 0.05 rule on n = 32–44 and are not FDR findings.

## CRISPR (DepMap 24Q4 Chronos)

Lung cell lines with Chronos for TACSTD2 and CLDN4: **n = 126** (NSCLC 98). Dependent = gene-dependency probability > 0.5.

| gene | n | mean Chronos (SD) | n dependent (prob>0.5) |
|---|---:|---|---:|
| TACSTD2 | 126 | -0.0287 (0.1214) | 0 / 126 |
| CLDN4 | 126 | +0.0342 (0.1316) | 0 / 126 |
| KRAS | 126 | -0.7039 (0.6033) | 69 / 126 |
| EEF2 | 126 | -2.2764 (0.5497) | 125 / 126 |

Chronos versus the same RNA scores, on lines with both CRISPR and RNA. A null unadjusted correlation means there is no association for CLDN4 to attenuate.

| cohort | score | TROP2 Chronos | TROP2 Chronos given CLDN4 Chronos | CLDN4 Chronos given TROP2 Chronos | TROP2 Chronos given CLDN4 RNA |
|---|---|---|---|---|---|
| lung overlap | Hallmark IFN-γ | -0.034 (n=123, p=0.7109, q=0.8531); boot CI [-0.22, +0.16] | -0.037 (n=123, p=0.6858, q=0.8230); boot CI [-0.23, +0.14]; |partial|/|unadj|=1.10; no_unadjusted_association | +0.079 (n=123, p=0.3898, q=0.7120); boot CI [-0.10, +0.25]; |partial|/|unadj|=1.02; no_unadjusted_association | +0.009 (n=123, p=0.9189); boot CI [-0.18, +0.19]; |partial|/|unadj|=0.28; no_unadjusted_association |
| lung overlap | MHC-I (HLA-A/B/C+B2M) | -0.061 (n=123, p=0.5026, q=0.7539); boot CI [-0.24, +0.12] | -0.065 (n=123, p=0.4746, q=0.7120); boot CI [-0.24, +0.12]; |partial|/|unadj|=1.07; no_unadjusted_association | +0.101 (n=123, p=0.2682, q=0.7120); boot CI [-0.08, +0.27]; |partial|/|unadj|=1.03; no_unadjusted_association | -0.029 (n=123, p=0.7498); boot CI [-0.21, +0.15]; |partial|/|unadj|=0.48; no_unadjusted_association |
| lung overlap | Immune-ligand score | -0.099 (n=123, p=0.2739, q=0.7539); boot CI [-0.27, +0.09] | -0.100 (n=123, p=0.2747, q=0.7120); boot CI [-0.28, +0.08]; |partial|/|unadj|=1.00; no_unadjusted_association | +0.008 (n=123, p=0.9296, q=0.9296); boot CI [-0.17, +0.18]; |partial|/|unadj|=2.00; no_unadjusted_association | -0.068 (n=123, p=0.4542); boot CI [-0.25, +0.11]; |partial|/|unadj|=0.69; no_unadjusted_association |
| NSCLC overlap | Hallmark IFN-γ | +0.077 (n=95, p=0.4593, q=0.7846); boot CI [-0.14, +0.28] | +0.062 (n=95, p=0.5556, q=0.8240); boot CI [-0.16, +0.28]; |partial|/|unadj|=0.80; no_unadjusted_association | +0.096 (n=95, p=0.3577, q=0.8240); boot CI [-0.12, +0.30]; |partial|/|unadj|=0.90; no_unadjusted_association | +0.118 (n=95, p=0.2592); boot CI [-0.10, +0.32]; |partial|/|unadj|=1.53; no_unadjusted_association |
| NSCLC overlap | MHC-I (HLA-A/B/C+B2M) | +0.029 (n=95, p=0.7833, q=0.7846); boot CI [-0.18, +0.23] | +0.023 (n=95, p=0.8240, q=0.8240); boot CI [-0.19, +0.22]; |partial|/|unadj|=0.81; no_unadjusted_association | +0.033 (n=95, p=0.7553, q=0.8240); boot CI [-0.17, +0.22]; |partial|/|unadj|=0.89; no_unadjusted_association | +0.055 (n=95, p=0.5980); boot CI [-0.15, +0.25]; |partial|/|unadj|=1.93; no_unadjusted_association |
| NSCLC overlap | Immune-ligand score | -0.028 (n=95, p=0.7846, q=0.7846); boot CI [-0.23, +0.17] | -0.040 (n=95, p=0.7007, q=0.8240); boot CI [-0.24, +0.16]; |partial|/|unadj|=1.41; no_unadjusted_association | +0.078 (n=95, p=0.4557, q=0.8240); boot CI [-0.13, +0.28]; |partial|/|unadj|=1.07; no_unadjusted_association | +0.004 (n=95, p=0.9717); boot CI [-0.20, +0.20]; |partial|/|unadj|=0.13; no_unadjusted_association |

## What the calls say

- rna NSCLC Hallmark IFN-γ: TROP2|CLDN4 call=retained, partial r=0.251
- rna NSCLC MHC-I (HLA-A/B/C+B2M): TROP2|CLDN4 call=no_unadjusted_association, partial r=0.058
- rna NSCLC Immune-ligand score: TROP2|CLDN4 call=retained, partial r=0.347
- rna lung Hallmark IFN-γ: TROP2|CLDN4 call=retained, partial r=0.407
- rna lung MHC-I (HLA-A/B/C+B2M): TROP2|CLDN4 call=retained, partial r=0.234
- rna lung Immune-ligand score: TROP2|CLDN4 call=retained, partial r=0.551
- protein lung Hallmark IFN-γ: TROP2|CLDN4 call=ns_after_adjustment_point_estimate_still_half_or_more, partial r=0.225
- protein lung MHC-I (HLA-A/B/C+B2M): TROP2|CLDN4 call=no_unadjusted_association, partial r=0.117
- protein lung Immune-ligand score: TROP2|CLDN4 call=attenuated, partial r=0.132
- protein NSCLC Hallmark IFN-γ: TROP2|CLDN4 call=attenuated, partial r=0.111
- protein NSCLC MHC-I (HLA-A/B/C+B2M): TROP2|CLDN4 call=no_unadjusted_association, partial r=-0.028
- protein NSCLC Immune-ligand score: TROP2|CLDN4 call=no_unadjusted_association, partial r=0.077
- crispr lung Hallmark IFN-γ: TROP2|CLDN4 call=no_unadjusted_association, partial r=-0.037
- crispr lung MHC-I (HLA-A/B/C+B2M): TROP2|CLDN4 call=no_unadjusted_association, partial r=-0.065
- crispr lung Immune-ligand score: TROP2|CLDN4 call=no_unadjusted_association, partial r=-0.100
- NSCLC RNA Hallmark IFN-γ residual comparison: |TROP2|CLDN4|=0.251 (retained), |CLDN4|TROP2|=0.015 (attenuated); larger absolute partial is TROP2.
- NSCLC RNA MHC-I (HLA-A/B/C+B2M) residual comparison: |TROP2|CLDN4|=0.058 (no_unadjusted_association), |CLDN4|TROP2|=0.086 (ns_after_adjustment_point_estimate_still_half_or_more); larger absolute partial is CLDN4.
- NSCLC RNA Immune-ligand score residual comparison: |TROP2|CLDN4|=0.347 (retained), |CLDN4|TROP2|=0.059 (attenuated); larger absolute partial is TROP2.

## What this is not

- Not a new estimate of the TROP2–CLDN4 protein correlation. That pair is the covariate, recomputed above.
- Not immune exclusion, ICI response, or IFN stimulation. The dish has no T cells.
- Not surface MHC or PD-L1. CD274 protein, where present, is TMT abundance.
- Not a claim that CLDN4 causes or blocks the TROP2 association. Partials remove shared rank variation only.
- Not RPPA n=118. That antibody panel does not contain TROP2 or CLDN4.

## Rerun

```bash
python3 -m pip install -r scripts/depmap_trop2_ifn_partial_cldn4/requirements.txt
python3 scripts/depmap_trop2_ifn_partial_cldn4/download.py --outdir data/depmap_trop2_ifn_partial_cldn4
python3 scripts/depmap_trop2_ifn_partial_cldn4/analyze.py --data data/depmap_trop2_ifn_partial_cldn4 --outdir results/depmap_trop2_ifn_partial_cldn4
```


# TACSTD2–T/NK attenuation: CLDN4 versus negative controls

Partial Spearman of TACSTD2 versus a T/NK measure, conditioning on one gene at a time. The negative controls are CLDN3, CLDN7, EPCAM, MUC1, and KRT19. Attenuation is partial ρ minus unadjusted ρ. A positive value means the correlation moved up. Percent attenuated is reported only as `100 × (ρ_unadj − ρ_partial) / ρ_unadj`. Head-to-head p-values swap the CLDN4 and control labels within cohort (10,000 draws). Benjamini–Hochberg q-values are across those five controls.

**CLDN7 attenuates about as much as CLDN4 where a negative TACSTD2–T/NK association is actually present.** The pre-specified TCGA seven-cohort pool does separate them, but that pool’s unadjusted correlation is weak and mixes opposite signs. Concordant-4 has no TACSTD2–T/NK association to attenuate.

## Pipeline check

Same concordant-4 units as the locked CLDN4 result (n = 65). Maximum absolute difference versus that patient table is 5×10⁻⁵ percentage points for CLDN4 % positive and 5×10⁻⁷ for the T/NK fraction. Malignant CLDN4 % positive versus T/NK reproduces the locked DerSimonian–Laird result: ρ = −0.531, p = 1.65×10⁻⁵, I² = 0%. TACSTD2 is not part of the malignant gate. Conditioning on EPCAM or KRT19 in GSE123902 and GSE189357 is entangled with that gate; the author-label subset (GSE131907 + GSE205335, n = 43) is the check that avoids it.

## Concordant-4

Primary score: malignant percent positive. Outcome: patient T/NK fraction.

| conditioner | partial ρ | p | I² | change in ρ | label-swap p vs CLDN4 | q |
|---|---:|---:|---:|---:|---:|---:|
| none | −0.112 | 0.46 | 15% | 0 |  |  |
| CLDN4 | +0.216 | 0.12 | 0% | +0.328 |  |  |
| CLDN7 | +0.010 | 0.94 | 0% | +0.122 | 0.079 | 0.13 |
| KRT19 | +0.007 | 0.96 | 0% | +0.118 | 0.063 | 0.13 |
| MUC1 | +0.048 | 0.73 | 0% | +0.160 | 0.039 | 0.13 |
| EPCAM | +0.066 | 0.65 | 0% | +0.177 | 0.14 | 0.17 |
| CLDN3 | +0.098 | 0.49 | 0% | +0.210 | 0.29 | 0.29 |

The unadjusted interval is −0.39 to +0.18. This is the same null already reported for TACSTD2 % positive (ρ = −0.112, p = 0.46, I² = 15%). There is no negative association for CLDN4 to explain.

Conditioning on CLDN4 moves the point estimate further than conditioning on CLDN7 (+0.328 versus +0.122) and lands at +0.216 rather than at zero. The label-swap test does not separate those two shifts (p = 0.079, q = 0.13). No control survives Benjamini–Hochberg at 0.05. TACSTD2 is about equally correlated with CLDN4, CLDN7, and EPCAM (DL ρ 0.53–0.55), so the larger CLDN4 shift is not a collinearity artifact that the exchangeability test can pin down at this n.

The author-label subset agrees: unadjusted ρ = +0.034 (p = 0.85); partial given CLDN4 is +0.331 (p = 0.042); versus CLDN7 the label-swap p is 0.13. Mean log1p is not negative either (unadjusted ρ = +0.041, p = 0.78). Partialling CLDN4 makes that positive association larger (ρ = +0.310, p = 0.025). That is shared variance with CLDN4, not attenuation of an immune-cold TACSTD2 effect.

The reverse partial does not move. CLDN4 % positive versus T/NK, adjusted for TACSTD2, stays ρ = −0.543 (p = 2.0×10⁻⁵, I² = 0%). TACSTD2 does not account for the locked CLDN4 association.

## TCGA, pre-specified seven cohorts

LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD. Primary tumors, Xena GDC STAR log2(TPM+1), sample type 01, replicate aliquots averaged. Outcome: mean of CD3D, CD3E, CD3G, CD8A, NKG7, GNLY, and KLRD1. No purity term and no extra keratin term, because KRT19 is one of the controls. N = 3,444.

| conditioner | partial ρ | p | I² | change in ρ | label-swap p vs CLDN4 | q |
|---|---:|---:|---:|---:|---:|---:|
| none | −0.067 | 0.081 | 79% | 0 |  |  |
| CLDN4 | +0.001 | 0.98 | 79% | +0.068 |  |  |
| CLDN7 | −0.024 | 0.49 | 74% | +0.043 | 0.0057 | 0.0057 |
| KRT19 | −0.036 | 0.29 | 73% | +0.031 | 0.0007 | 0.0009 |
| EPCAM | −0.046 | 0.20 | 75% | +0.021 | 0.0001 | 0.0002 |
| MUC1 | −0.052 | 0.24 | 84% | +0.015 | 0.0001 | 0.0002 |
| CLDN3 | −0.052 | 0.12 | 72% | +0.015 | 0.0001 | 0.0002 |

In this pre-specified pool, CLDN4 moves the correlation to zero and CLDN7 does not move it as far. The label-swap test separates them. CLDN3, EPCAM, MUC1, and KRT19 attenuate less than CLDN4. CD3-only and CD8A give the same order (CLDN4 change +0.060 and +0.057; CLDN7 +0.045 and +0.038).

That pool is a poor summary of a single negative association. I² is 79%. The cohort correlations are LUAD −0.097, BRCA −0.068, CESC **+0.146**, KIRC +0.022, STAD −0.174, BLCA −0.165, PAAD −0.133. CESC is positive (p = 0.011). KIRC is flat. BRCA is the cohort where CLDN7 attenuates and CLDN4 does not (changes +0.045 versus −0.002). CLDN4 is the largest shift in the other five cohorts, including the positive CESC cohort, where “attenuation of a negative association” is the wrong description: the partial given CLDN4 is +0.219.

TACSTD2’s correlation with CLDN4 (DL ρ = 0.43) is similar to its correlation with KRT19 (0.41) and higher than with CLDN7 (0.26). KRT19 is therefore a real collinearity control, and CLDN4 still shifts the seven-cohort partial more than KRT19 does.

## Where the unadjusted correlation is negative

LUAD, BRCA, STAD, BLCA, and PAAD (n = 2,607). This set was chosen after seeing the signs. It is not the pre-specified test.

| conditioner | partial ρ | p | I² | change in ρ |
|---|---:|---:|---:|---:|
| none | −0.115 | 5.0×10⁻⁷ | 21% | 0 |
| CLDN4 | −0.058 | 0.012 | 20% | +0.058 |
| CLDN7 | −0.068 | 0.023 | 51% | +0.047 |
| KRT19 | −0.079 | 0.0005 | 19% | +0.037 |
| EPCAM | −0.082 | 0.009 | 56% | +0.033 |
| CLDN3 | −0.094 | 0.0001 | 30% | +0.021 |
| MUC1 | −0.100 | 0.029 | 79% | +0.015 |

Here the unadjusted association is negative and the cohorts agree better (I² = 21%). CLDN4 removes about half of it and a residual remains (ρ = −0.058, p = 0.012, CI −0.102 to −0.013). CLDN7 removes a similar share (residual ρ = −0.068, p = 0.023). The label-swap difference is +0.010, p = 0.27, q = 0.27. **CLDN7 attenuates as much as CLDN4 on this contrast.**

CLDN4 still exceeds CLDN3 (p = 0.0001), MUC1 (p = 0.0006), and EPCAM (p = 0.0056) in this set. KRT19 is close (Δ = +0.021, p = 0.045, q = 0.057).

LUSC is outside the seven. Its unadjusted ρ is −0.248 (p = 1.9×10⁻⁸, n = 501). CLDN4 leaves −0.225 and CLDN7 leaves −0.222. KRT19 leaves −0.146. In squamous lung, CLDN4 does not attenuate the TACSTD2–T/NK score, and KRT19 attenuates more.

LUAD alone, inside both summaries: unadjusted ρ = −0.097 (p = 0.027, n = 516). Partial given CLDN4 is −0.030 (p = 0.50). Partial given CLDN7 is −0.055 (p = 0.21). Partial given KRT19 is −0.096 (p = 0.030). In LUAD, KRT19 removes essentially nothing and CLDN4 removes more of the point estimate than CLDN7 does. That single-cohort gap was not given its own label-swap test.

## What this does not say

Partial correlation is not a mediation test and not a spatial exclusion test. The locked CLDN4–T/NK result in concordant-4 is unchanged, and it is not explained by TACSTD2. Public bulk and single-cell tables do not show a CLDN4-specific attenuation of a TACSTD2–T/NK association once the contrast is restricted to cohorts where that association is negative.

[Partial ρ, concordant-4 and TCGA funnel-7](results/figures/partial_rho_forest.png)

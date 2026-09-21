# CPTAC protein: larger CLDN4 contrasts with the interval kept off zero

Treatment-naive CPTAC TMT freeze v1.2 (LSCC Satpathy *Cell* 2021; LUAD Gillette *Cell* 2020). LSCC and LUAD are not pooled. No ICI column. Endpoint proteins are never filled. Undetected CLDN4 stays missing in the primary rows. `floor` ties those tumors at the lowest measured CLDN4 value; `below_min` places them one log2 unit under that value. Those are rank rules, not imputed abundances.

The grid is fixed in `analyze.py` before the correlations are ranked. MHC-I/CD8 panels are every subset of size 1–4 of HLA-A/B/C/E/F and CD8A (B2M is absent; HLA-G and CD8B are quantified in 20 and 10 tumors and sit under the 40-tumor floor). The score is the mean of cohort z-scores (LSCC ddof=0, the same scale as the locked HLA-A/B/C score). Completeness is strict. Sample cuts are grade G2/G3, stage I/II, a WES-purity median split, and LSCC basaloid versus not. Stage III has too few measured CLDN4 tumors to test. Quantile cuts are Q1/Q4, 20/80, 15/85, 10/90, and tertiles. Spearman rows need n≥20. Group rows need ≥8 tumors on a side. An inverse row is eligible only when its 2,000-draw bootstrap 95% interval lies entirely below 0. A STING row is eligible only when the Q1−Q4 interval lies entirely above 0. Nominal p-values are not a multiplicity correction. The permutation p is for the maximum inside the stated slice.

## Locked LSCC anchors, recomputed

CLDN4 versus the HLA-A/B/C mean z: ρ=-0.410, n=78, bootstrap 95% CI -0.585 to -0.204. WES partial ρ=-0.320 (CI -0.535 to -0.075, n=77). The published point estimate was −0.410, and the published partial was −0.320.

CLDN4 versus CD8A: ρ=-0.444, n=78, bootstrap 95% CI -0.631 to -0.219. WES partial ρ=-0.376 (CI -0.569 to -0.129, n=77). The published point estimate was −0.444, and the published partial was −0.376.

## Largest full-cohort Spearman, measured CLDN4

The largest eligible |ρ| is HLA-A+HLA-C+HLA-F+CD8A (MHC1_CD8, strict, k=4); filter=all; missingness=observed; mode=full; n=78; ρ=-0.531; 95% CI -0.680 to -0.333; p=5.75e-07; imputed CLDN4 in the pair=0.

That is a larger inverse than the locked HLA-A/B/C ρ (-0.410) and the locked CD8A ρ (-0.444), on the same 78 tumors. Every panel in the top of this list contains CD8A and HLA-A. HLA-F is in the winning set; it is nonclassical MHC-I, and the panel is labeled that way.

WES-purity partial ρ=-0.471 (CI -0.644 to -0.273, n=77). The partial interval stays below 0.

| Panel | n | ρ | 95% CI | partial ρ | partial CI |
| --- | --- | --- | --- | --- | --- |
| HLA-A+HLA-C+HLA-F+CD8A | 78 | -0.531 | -0.680 to -0.333 | -0.471 | -0.644 to -0.273 |
| HLA-A+HLA-F+CD8A | 78 | -0.523 | -0.698 to -0.305 | -0.467 | -0.655 to -0.242 |
| HLA-A+HLA-C+CD8A | 78 | -0.515 | -0.671 to -0.323 | -0.451 | -0.617 to -0.243 |
| HLA-A+HLA-B+HLA-F+CD8A | 78 | -0.513 | -0.675 to -0.301 | -0.451 | -0.644 to -0.219 |
| HLA-A+HLA-B+HLA-C+CD8A | 78 | -0.499 | -0.653 to -0.307 | -0.428 | -0.610 to -0.214 |

Search permutation on the full-cohort continuous slice (every MHC/CD8 panel and all three missingness rules; Fisher interval entirely below 0): max |ρ|=0.531, permutation p=9.99e-04 (1000 shuffles; null median 0.000, null 95th 0.309). The quoted interval is the bootstrap. On this slice the Fisher maximum and the bootstrap maximum are the same panel.

## Largest full-cohort Q4−Q1 gap, measured CLDN4

Largest median-z gap: HLA-F+CD8A (MHC1_CD8, strict, k=2); filter=all; missingness=observed; cut=q25_75; n_low=20 vs n_high=20; Δ median z Q4−Q1=-1.567; 95% CI -1.977 to -0.536; Cliff(Q4 vs Q1)=-0.630; MWU p=6.87e-04; imputed in Q1=0.

The strongest rank separation in the same eligible set is HLA-A+HLA-C+CD8A: Cliff(Q4 vs Q1)=-0.755, Δ median z Q4−Q1=-1.006 (CI -1.389 to -0.504), MWU p=4.68e-05. CD8A alone has Cliff=-0.625 and Δz=-1.188. The pre-specified pick is the median-z gap. A mean of z-scores is not on the same scale as one gene, and the HLA-F+CD8A z-gap is not a larger rank separation than HLA-A+HLA-C+CD8A.

Single-protein median log2 (Q4−Q1), same 20 vs 20 tumors: HLA-A -1.081; CD8A -0.764; HLA-F -0.632; HLA-C -0.524; HLA-B -0.501; HLA-E -0.329. HLA-A is the largest single-protein gap. CD8A Q4 is lower than Q1 by 0.764 log2 TMT units.

Search permutation for the median-z rule (measured CLDN4, full cohort, Q1 vs Q4): max (median Q1 − median Q4)=1.567, permutation p=0.00249 (400 shuffles, 400 inner bootstrap draws; null 95th 0.736).

WES residual of the median-z winner (Q4−Q1): -0.781 (CI -1.515 to +0.092, n=20 vs 20). That residual interval is not entirely below 0.

WES residual of the Cliff leader HLA-A+HLA-C+CD8A (Q4−Q1): -0.683 (CI -1.166 to -0.077, n=20 vs 20). That residual interval stays below 0.

CD8A alone, same Q4−Q1 residual: -0.701 (CI -1.394 to -0.064). That interval stays below 0.

HLA-A/B/C alone, same residual: -0.424 (CI -0.898 to +0.204). That interval crosses 0.

## Larger |ρ| once cuts and strata are opened

Largest |ρ| anywhere in the measured-CLDN4 grid: CD8A (CD8, strict, k=1); filter=all; missingness=observed; mode=tails; cut=q15_85; n=24; ρ=-0.807; 95% CI -0.897 to -0.624; p=1.90e-06; imputed CLDN4 in the pair=0.

Tail-only Spearman uses only the tumors outside the middle of the CLDN4 distribution. It is not the full-cohort ρ of −0.41.

The strongest full Spearman inside a stratum, rather than a tail cut, is HLA-A+HLA-F+CD8A (MHC1_CD8, strict, k=3); filter=purity_high; missingness=observed; mode=full; n=41; ρ=-0.666; 95% CI -0.829 to -0.423; p=2.06e-06; imputed CLDN4 in the pair=0.

Search permutation over strata and tail cuts as well as the full cohort: Fisher-gate max |ρ|=0.812, permutation p=0.00599 (1000 shuffles; null median 0.574, null 95th 0.736). The null 95th is already large because the search includes n≈20 slices.

## Missingness rules

Tying the 30 undetected CLDN4 tumors at the floor, or placing them one log2 unit below it, does not increase |ρ|. The strongest filled continuous row is HLA-A+HLA-B+HLA-F+CD8A (MHC1_CD8, strict, k=4); filter=all; missingness=floor; mode=full; n=108; ρ=-0.474; 95% CI -0.629 to -0.283; p=2.21e-07; imputed CLDN4 in the pair=30. That is weaker than the measured-CLDN4 winner. The primary estimand stays the 78 quantified tumors.

## Antigen processing, separate from MHC-I

PSMB10 (APM, strict, k=1); filter=all; missingness=observed; mode=full; n=78; ρ=-0.424; 95% CI -0.597 to -0.208; p=1.09e-04; imputed CLDN4 in the pair=0.

TAP1/2, TAPBP, PSMB8/9/10, and NLRC5 were scored as their own family. The strongest of those rows does not beat the MHC-I/CD8 winner, and it is not called MHC-I.

## LUAD quartiles: DNA-PK and STING

Z-scores use ddof=1, matching the earlier LUAD search. Q1-versus-rest anchors on measured CLDN4 were recomputed before any new panel was ranked. Bootstrap intervals were not part of that earlier report.

| Prior panel, Q1 vs rest | Δ median z | 95% CI | n |
| --- | --- | --- | --- |
| DNA-PKcs+LIG4 | -0.262 | -0.839 to +0.042 | 20 vs 59 |
| DNA-PKcs+Ku80+LIG4 | -0.519 | -0.763 to +0.062 | 20 vs 59 |
| STING+TBK1+IRF3 | +0.418 | +0.086 to +0.653 | 20 vs 59 |
| STING+TBK1+IRF3+IRF7 | +0.500 | +0.226 to +0.728 | 20 vs 59 |

The two DNA-PK Q1-versus-rest gaps stay negative, and both bootstrap intervals cross 0. The two STING Q1-versus-rest gaps stay positive, and both intervals stay above 0.

Same panels, Q1 versus Q4:

| Panel, Q1 vs Q4 | Δ median z | 95% CI | n |
| --- | --- | --- | --- |
| PRKDC+LIG4 | -0.094 | -0.786 to +0.390 | 20 vs 20 |
| PRKDC+XRCC5+LIG4 | -0.357 | -0.713 to +0.219 | 20 vs 20 |
| STING1+TBK1+IRF3 | +0.345 | -0.012 to +0.717 | 20 vs 20 |
| STING1+TBK1+IRF3+IRF7 | +0.467 | +0.023 to +0.788 | 20 vs 20 |

Moving from Q1-versus-rest to Q1-versus-Q4 shrinks these four panels. DNA-PKcs+LIG4 falls from −0.262 to −0.094. DNA-PKcs+Ku80+LIG4 falls from −0.519 to −0.357. STING+TBK1+IRF3+IRF7 falls from +0.500 to +0.467, and that smaller gap still has an interval above 0. STING+TBK1+IRF3’s Q1-versus-Q4 interval crosses 0.

**DNA-PK, full cohort, Q1 vs Q4, measured CLDN4.**

No panel keeps the Q1−Q4 interval entirely below 0. The most negative point estimate is XRCC5+LIG4: Δ=-0.448 (CI -0.850 to +0.154, n=20 vs 20, MWU p=0.126). The interval crosses 0.

Search permutation over every nonempty DNA-PK subset on this slice: max thesis-direction gap=0.000, permutation p=1 (400 shuffles). The observed maximum is 0 because no panel cleared the interval.

**STING, full cohort, Q1 vs Q4, measured CLDN4.**

IRF3+IRF7 (STING, strict, k=2); filter=all; missingness=observed; cut=q25_75; n_low=20 vs n_high=20; Δ median z Q1−Q4=+0.734; 95% CI +0.076 to +1.262; MWU p=0.0499; imputed in Q1=0.

That gap is larger than the locked four-gene Q1-versus-Q4 gap (+0.467) and larger than the locked Q1-versus-rest gap (+0.500). It is an IRF3+IRF7 panel, not the cGAS–STING–TBK1 set.

Search permutation over every nonempty STING-pathway subset on this slice: max thesis-direction gap=0.734, permutation p=0.0623 (400 shuffles, 300 inner bootstrap draws; null 95th 0.760).

WES residual of that STING winner, scored as Q4−Q1 so the thesis direction is negative: -0.390 (CI -1.038 to +0.205, n=20 vs 20). The residual interval crosses 0.

**Strata and other cuts.**

DNA-PK: XRCC5+NHEJ1 (DNA-PK, strict, k=2); filter=purity_low; missingness=observed; cut=q20_80; n_low=9 vs n_high=9; Δ median z Q1−Q4=-0.862; 95% CI -1.181 to -0.042; MWU p=0.0217; imputed in Q1=0.

STING: STING1+IRF3+IRF7 (STING, strict, k=3); filter=non_acinar; missingness=observed; cut=tertile; n_low=9 vs n_high=9; Δ median z Q1−Q4=+1.508; 95% CI +0.129 to +1.847; MWU p=0.00807; imputed in Q1=0.

Those two rows are 9 versus 9. They are not the full cohort, and the permutation above does not cover purity halves, histology, or the 20/80 and tertile cuts. HLA proteins were left out of this LUAD list; adding them previously removed the STING contrast.

## What this does not claim

- The full-cohort LSCC Spearman moves from −0.410 / −0.444 to −0.531 on a declared MHC-I/CD8 grid, with the bootstrap and the WES partial both entirely below 0. It does not replace the locked HLA-A/B/C number, and it is not an ICI result.
- The LSCC Q4−Q1 median-z gap that maximizes |Δz| does not survive a WES residual. The larger rank separation is HLA-A+HLA-C+CD8A, not the z-gap winner.
- Tail-only ρ and grade or purity slices are different estimands. The tail search has a high null.
- LUAD DNA-PK does not gain a larger quartile gap whose interval stays below 0. The earlier Q1-versus-rest point estimates themselves have bootstrap intervals that cross 0.
- LUAD STING Q1 versus Q4 is larger for IRF3+IRF7 than for the locked four-gene panel. The search permutation is 0.062, and the WES residual interval crosses 0.
- Floor and below-min do not increase the LSCC |ρ|.

## Reproduce

```bash
python3 methods/cptac_max_effect/download.py --outdir data/cptac_max_effect
python3 methods/cptac_max_effect/analyze.py --data data/cptac_max_effect --outdir methods/cptac_max_effect
```


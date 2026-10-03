# Concordant-4 malignant CLDN4: NHEJ, cGAS-STING, IFN/APM, mediation

ADDITIVE. **CLDN4-only.** Not a re-derivation of the locked T/NK result
(malignant CLDN4 %pos vs T/NK ρ = −0.53, N = 65). Not a dual-high
TACSTD2×CLDN4 score. Cohorts are the four that already agree:
**GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not GSE148071,
GSE127465, GSE207422, GSE154826, or CD45+/T-only extracts.

The locked tumor-cell result this sits next to: malignant CLDN4-high
patients have lower own IFN and MHC-I/APM. IFN and APM are recomputed
here on the same matrices as a calibration, then NHEJ and cGAS-STING
are added, and the path CLDN4 → NHEJ → IFN is tested. This is an
observational patient-level path. It does not show that CLDN4 causes
the IFN change.

## Unit and split

Patient / donor / sample is the unit. Do not quote cell counts as n.

| Piece | Choice |
|---|---|
| Exposure | malignant CLDN4 %pos, within-cohort quartiles (same labels as PR #503) |
| High vs low | Q4 vs Q1, mid quartiles unused in the binary contrast |
| Scores | mean log2(TMM-CPM+1) of module genes present after the locked count filter |
| Aliases | MB21D1→CGAS, TMEM173→STING1, MRE11A→MRE11, placed on the locked TMM factors. IFN/APM symbols are not relabeled |
| Model | OLS ~ cohort + CLDN4_Q4 (reference cohort GSE123902) |
| Matrix | marker-malignant UMI-sum (GSE123902, GSE189357); author-malignant UMI-sum (GSE131907, GSE205335) |
| In the count matrix | 64 of 65 units (P4001 is Q1 on the %pos vector and is absent from the UMI-sum) |
| Stacked Q4 vs Q1 | n = 34 (18 / 16) |

GSE189357 Q4 has 2 patients, so its within-cohort binary contrast is skipped.
Those two patients stay in the stacked model.

IFN calibration against the locked family score (PR #503, logFC −0.584, p = 0.00243):
this run logFC = -0.584, p = 0.002434. APM locked −0.779; this run -0.779.

## Modules

| module | definition | genes in score | held out of the log matrix |
|---|---|---:|---|
| NHEJ | KEGG hsa03450 `KEGG_NON_HOMOLOGOUS_END_JOINING` | 12 | DNTT |
| NHEJ-core | catalytic c-NHEJ (sensitivity) | 11 | PAXX |
| cGAS-STING | sensors/adapters; IFN-hallmark and NHEJ genes removed | 11 | — |
| IFN | Hallmark IFNα ∪ IFNγ (same lists as PR #503) | 221 | 3 genes below the count filter |
| APM | custom MHC-I / antigen-processing panel (PR #503) | 21 | — |
| proliferation | MKI67, PCNA, TOP2A, MCM2, CDK1, CCNB1, BIRC5, UBE2C | 8 | covariate only |

NHEJ genes used: DCLRE1C, FEN1, LIG4, MRE11, NHEJ1, POLL, POLM, PRKDC, RAD50, XRCC4, XRCC5, XRCC6.
MRE11 is the sum of the MRE11 / MRE11A symbol, whichever the GEO build used.
DNTT is absent (not expressed in these lung malignant sums).

cGAS-STING genes used: CGAS, STING1, TBK1, IKBKE, IRF3, IFI16, DDX41, TRIM56, DHX9, DHX36, PQBP1.
CGAS is MB21D1 in GSE123902 and GSE131907, and CGAS in GSE205335 and GSE189357.
STING1 is TMEM173 in all four matrices.
Held out of cGAS-STING because they sit in Hallmark IFN: IRF7, ZBP1, TRIM21, SAMHD1.
Held out because they sit in KEGG NHEJ: PRKDC, XRCC5, XRCC6, MRE11.
Negative regulators TREX1 and NLRP4 are not shared across the four builds, so that pair was not scored.

## 1. CLDN4-high vs low (stacked Q4 vs Q1)

Positive logFC = higher in CLDN4-high. FDR is Benjamini–Hochberg across the seven scores fit in that model.

| module | n | n_Q1 / n_Q4 | logFC (95% CI) | p | FDR | rank-biserial r (MW p) |
|---|---:|---|---|---:|---:|---|
| NHEJ | 34 | 18 / 16 | +0.013 (-0.166, +0.192) | 0.8838 | 0.9339 | +0.132 (0.5233) |
| cGAS-STING | 34 | 18 / 16 | -0.266 (-0.531, -0.002) | 0.04877 | 0.1138 | -0.354 (0.08144) |
| IFN | 34 | 18 / 16 | -0.584 (-0.944, -0.224) | 0.002434 | 0.01704 | -0.674 (0.0008698) |
| APM | 34 | 18 / 16 | -0.779 (-1.356, -0.202) | 0.009883 | 0.03459 | -0.514 (0.01121) |
| NHEJ-core | 34 | 18 / 16 | +0.007 (-0.154, +0.167) | 0.9339 | 0.9339 | +0.125 (0.546) |
| cGAS-STING-core | 34 | 18 / 16 | -0.243 (-0.591, +0.106) | 0.1651 | 0.2889 | -0.222 (0.2771) |
| proliferation | 34 | 18 / 16 | +0.237 (-0.663, +1.136) | 0.5944 | 0.8322 | +0.132 (0.5233) |

Continuous CLDN4 %pos (global z, same n as the locked continuous family model):

| module | n | logFC per SD CLDN4 %pos (95% CI) | p | FDR |
|---|---:|---|---:|---:|
| NHEJ | 64 | -0.013 (-0.083, +0.057) | 0.7111 | 0.8297 |
| cGAS-STING | 64 | -0.054 (-0.157, +0.050) | 0.3025 | 0.7059 |
| IFN | 64 | -0.161 (-0.316, -0.006) | 0.04215 | 0.1475 |
| APM | 64 | -0.273 (-0.499, -0.046) | 0.01926 | 0.1348 |
| NHEJ-core | 64 | -0.017 (-0.079, +0.044) | 0.5763 | 0.8068 |
| cGAS-STING-core | 64 | +0.008 (-0.139, +0.156) | 0.9086 | 0.9086 |
| proliferation | 64 | +0.102 (-0.226, +0.430) | 0.5357 | 0.8068 |

### Within cohort

| cohort | module | n_Q1 / n_Q4 | logFC | p |
|---|---|---|---:|---:|
| GSE123902 | NHEJ | 4 / 3 | -0.060 | 0.8192 |
| GSE123902 | cGAS-STING | 4 / 3 | -0.586 | 0.01312 |
| GSE123902 | IFN | 4 / 3 | -0.679 | 0.1022 |
| GSE123902 | APM | 4 / 3 | -0.751 | 0.1082 |
| GSE131907 | NHEJ | 6 / 5 | -0.037 | 0.7682 |
| GSE131907 | cGAS-STING | 6 / 5 | +0.289 | 0.2281 |
| GSE131907 | IFN | 6 / 5 | -0.011 | 0.9708 |
| GSE131907 | APM | 6 / 5 | -0.168 | 0.7643 |
| GSE205335 | NHEJ | 5 / 6 | +0.046 | 0.8123 |
| GSE205335 | cGAS-STING | 5 / 6 | -0.708 | 0.007102 |
| GSE205335 | IFN | 5 / 6 | -1.243 | 0.001759 |
| GSE205335 | APM | 5 / 6 | -1.605 | 0.01607 |
| GSE189357 | NHEJ | — | — | skipped (Q tail < 3) |
| GSE189357 | cGAS-STING | — | — | skipped (Q tail < 3) |
| GSE189357 | IFN | — | — | skipped (Q tail < 3) |
| GSE189357 | APM | — | — | skipped (Q tail < 3) |

### Spearman, DerSimonian–Laird on Fisher-z

| module | k | N | ρ (95% CI) | p | I² |
|---|---:|---:|---|---:|---:|
| NHEJ | 4 | 64 | -0.004 (-0.289, +0.282) | 0.9808 | 12.1% |
| cGAS-STING | 4 | 64 | -0.189 (-0.588, +0.283) | 0.4366 | 65.7% |
| IFN | 4 | 64 | -0.324 (-0.565, -0.032) | 0.03027 | 17.2% |
| APM | 4 | 64 | -0.303 (-0.587, +0.048) | 0.08929 | 39.3% |

| cohort | module | n | ρ | p |
|---|---|---:|---:|---:|
| GSE123902 | NHEJ | 13 | -0.093 | 0.7615 |
| GSE131907 | NHEJ | 21 | -0.188 | 0.4137 |
| GSE205335 | NHEJ | 21 | -0.014 | 0.951 |
| GSE189357 | NHEJ | 9 | +0.583 | 0.09919 |
| GSE123902 | cGAS-STING | 13 | -0.522 | 0.06729 |
| GSE131907 | cGAS-STING | 21 | +0.279 | 0.2203 |
| GSE205335 | cGAS-STING | 21 | -0.510 | 0.01808 |
| GSE189357 | cGAS-STING | 9 | +0.133 | 0.7324 |
| GSE123902 | IFN | 13 | -0.549 | 0.05177 |
| GSE131907 | IFN | 21 | +0.017 | 0.9421 |
| GSE205335 | IFN | 21 | -0.475 | 0.02943 |
| GSE189357 | IFN | 9 | -0.317 | 0.4064 |
| GSE123902 | APM | 13 | -0.516 | 0.07075 |
| GSE131907 | APM | 21 | +0.001 | 0.9955 |
| GSE205335 | APM | 21 | -0.552 | 0.009482 |
| GSE189357 | APM | 9 | +0.050 | 0.8984 |

### NHEJ genes (stacked Q4 vs Q1)

| gene | logFC | p | FDR |
|---|---:|---:|---:|
| NHEJ1 | -0.239 | 0.1291 | 0.9242 |
| XRCC6 | -0.143 | 0.3916 | 0.9242 |
| MRE11 | +0.165 | 0.4147 | 0.9242 |
| LIG4 | +0.155 | 0.5075 | 0.9242 |
| RAD50 | +0.097 | 0.615 | 0.9242 |
| FEN1 | +0.156 | 0.663 | 0.9242 |
| PRKDC | -0.111 | 0.7295 | 0.9242 |
| XRCC5 | +0.053 | 0.7755 | 0.9242 |
| POLM | +0.049 | 0.8374 | 0.9242 |
| POLL | +0.024 | 0.8874 | 0.9242 |
| XRCC4 | -0.029 | 0.89 | 0.9242 |
| DCLRE1C | -0.022 | 0.9242 | 0.9242 |

### cGAS-STING genes (stacked Q4 vs Q1)

| gene | logFC | p | FDR |
|---|---:|---:|---:|
| CGAS | -1.187 | 0.001438 | 0.01582 |
| IFI16 | -1.685 | 0.01352 | 0.07438 |
| TBK1 | -0.366 | 0.03785 | 0.1388 |
| IKBKE | -0.577 | 0.1423 | 0.3618 |
| STING1 | +0.639 | 0.1772 | 0.3618 |
| DHX9 | +0.188 | 0.1988 | 0.3618 |
| PQBP1 | +0.292 | 0.2302 | 0.3618 |
| DDX41 | -0.210 | 0.2961 | 0.4072 |
| IRF3 | -0.057 | 0.7647 | 0.8906 |
| TRIM56 | +0.046 | 0.8096 | 0.8906 |
| DHX36 | -0.010 | 0.9628 | 0.9628 |

## 2. Mediation CLDN4 → NHEJ → IFN

Exposure = malignant CLDN4 %pos, global z on the units in the model.
Mediator = NHEJ module score. Outcome = IFN module score.
Cohort fixed effects in every equation. Indirect effect = a×b.
Uncertainty is a cohort-stratified bootstrap (B = 5000, seed 20260921);
the percentile interval is primary. Sobel is the normal approximation.
Proportion mediated = a×b / c, reported only as a descriptive ratio.

Paths: a is CLDN4 → NHEJ, b is NHEJ → IFN given CLDN4, c is the total
CLDN4 → IFN effect, c′ is the direct effect given NHEJ.

| model | n | a (p) | b (p) | total c (p) | direct c′ (p) | a×b (boot 95% CI) | boot p | Sobel p | prop. |
|---|---:|---|---|---|---|---|---:|---:|---:|
| primary: CLDN4% → NHEJ → IFN | 64 | -0.013 (0.7111) | -0.725 (0.01094) | -0.161 (0.04215) | -0.170 (0.02486) | +0.009 (-0.060, +0.078) | 0.7744 | 0.7125 | -0.06 |
| sensitivity: CLDN4% → NHEJ-core → IFN | 64 | -0.017 (0.5763) | -0.400 (0.2271) | -0.161 (0.04215) | -0.168 (0.03418) | +0.007 (-0.029, +0.044) | 0.7068 | 0.6098 | -0.04 |
| sensitivity: CLDN4% → NHEJ → IFN, given proliferation | 64 | -0.022 (0.4894) | -0.730 (0.02097) | -0.154 (0.05155) | -0.171 (0.02679) | +0.016 (-0.043, +0.071) | 0.5676 | 0.5043 | -0.11 |
| parallel: CLDN4% → cGAS-STING → IFN | 64 | -0.054 (0.3025) | +1.143 (1.39e-12) | -0.161 (0.04215) | -0.099 (0.05592) | -0.061 (-0.185, +0.048) | 0.2648 | 0.3015 | 0.38 |
| sensitivity: CLDN4% → cGAS-STING-core → IFN | 64 | +0.008 (0.9086) | +0.693 (1.08e-08) | -0.161 (0.04215) | -0.167 (0.006285) | +0.006 (-0.105, +0.102) | 0.9244 | 0.9083 | -0.04 |
| secondary outcome: CLDN4% → NHEJ → APM | 64 | -0.013 (0.7111) | -0.656 (0.1216) | -0.273 (0.01926) | -0.281 (0.01491) | +0.009 (-0.063, +0.077) | 0.8328 | 0.7173 | -0.03 |
| sensitivity: CLDN4 expr → NHEJ → IFN | 64 | -0.079 (0.0138) | -0.715 (0.02171) | +0.040 (0.5962) | -0.016 (0.8303) | +0.056 (-0.009, +0.163) | 0.0948 | 0.08398 | 1.41 |
| LOO drop GSE123902: CLDN4% → NHEJ → IFN | 51 | -0.125 (0.5183) | -9.608 (0.119) | -0.262 (0.2215) | +0.751 (0.624) | +1.203 (-71.244, +75.933) | 0.9936 | 0.5469 | -4.59 |
| LOO drop GSE131907: CLDN4% → NHEJ → IFN | 43 | +0.005 (0.9174) | -0.643 (0.02902) | -0.317 (0.0007555) | -0.314 (0.0004936) | -0.003 (-0.092, +0.069) | 0.9128 | 0.9169 | 0.01 |
| LOO drop GSE205335: CLDN4% → NHEJ → IFN | 43 | -0.019 (0.6187) | -0.716 (0.03914) | -0.043 (0.6123) | -0.056 (0.4844) | +0.014 (-0.073, +0.105) | 0.732 | 0.6252 | -0.33 |
| LOO drop GSE189357: CLDN4% → NHEJ → IFN | 55 | -0.021 (0.5865) | -0.724 (0.02113) | -0.169 (0.06124) | -0.184 (0.03407) | +0.016 (-0.055, +0.096) | 0.6452 | 0.5936 | -0.09 |

Primary read. KEGG NHEJ does not differ by CLDN4 quartile
(logFC +0.013, p = 0.8838). No NHEJ gene has
FDR < 0.05. The indirect effect is a×b = +0.009
(95% CI -0.060 to +0.078,
boot p = 0.7744). Path a is flat, so these units
do not support NHEJ as the mediator of the CLDN4–IFN association.
Path b is a partial association of the NHEJ score with IFN given CLDN4;
it is not the mediation result. Catalytic NHEJ-core and the
proliferation-adjusted model give the same null indirect effect.
The leave-one-out row that drops GSE123902 has an indirect-effect
interval spanning tens of log units. That row is unstable and is not
a result. The other three leave-one-out intervals still cover 0.

cGAS-STING module, stacked Q4 vs Q1: logFC -0.266
(p = 0.04877, FDR = 0.1138). The continuous
model and the four-cohort Spearman pool do not separate this module
from 0. The four-gene core (CGAS, STING1, TBK1, IRF3) is
logFC -0.243 (p = 0.1651).
Gene-level, the module dip is CGAS
(logFC -1.187, p = 0.001438, FDR = 0.01582)
and IFI16. STING1 is not in that direction
(logFC +0.639, p = 0.1772).
Parallel mediation through the 11-gene module has a bootstrap interval
that includes 0.

IFN and APM reproduce the locked direction (IFN logFC -0.584,
p = 0.002434; APM logFC -0.779, p = 0.009883).

Sensitivities in the same table: catalytic NHEJ-core instead of KEGG NHEJ;
proliferation score as a covariate; cGAS-STING and the four-gene core as
parallel mediators; APM as a parallel outcome; pseudobulk CLDN4 expression
instead of %pos; leave-one-cohort-out.

## 3. Specification search (cutoffs, modules, covariates, estimators)

The pre-specified KEGG model is not the end of the test. `search_mediation.py` refit the same patient table under 16,560 specifications: CLDN4 as %pos, within-cohort rank, pseudobulk expression, or the locked per-cell mean; Q4/Q1, median, and tertile cuts; KEGG, catalytic, Ku, ligase, Reactome (histones removed), equal-weight, rank, PC1, and each KEGG gene; IFN, IFNα, IFNγ, a compact ISG list, and an IFN rank score; covariates cohort, proliferation, malignant-cell count, library size, and KRT8/18/19/EPCAM; subsets all four cohorts, Q4 vs Q1, ADC-only, and drop-GSE131907.

Estimators on the shortlist: Sobel, cohort-stratified case bootstrap (B = 5,000), Imai normal quasi-Bayesian draws, and pingouin 0.6.1 `mediation_analysis` (B = 5,000). No specification was simulated.

Support for the path means the indirect effect is negative (higher CLDN4 pulls IFN down through the mediator) and the total CLDN4 → IFN effect is negative. **None of the 16,560 specifications had Sobel p < 0.05 under that rule.**

The strongest NHEJ-module specification, and the single-gene and DNA-repair neighbors, case-bootstrap 95% intervals:

| specification | n | a (p) | b (p) | total c (p) | a×b | case-bootstrap 95% CI | boot p | pingouin p | fraction of c |
|---|---:|---|---|---|---:|---|---:|---:|---:|
| KEGG NHEJ, pre-specified | 64 | −0.013 (0.71) | −0.725 (0.011) | −0.161 (0.042) | +0.009 | −0.060 to +0.078 | 0.77 | 0.76 | — |
| KEGG NHEJ + keratin | 64 | +0.028 (0.42) | −0.340 (0.23) | −0.266 (5.8e-4) | −0.009 | −0.056 to +0.024 | 0.61 | 0.63 | 0.04 |
| Reactome NHEJ + keratin | 64 | +0.045 (0.14) | −0.542 (0.091) | −0.266 (5.8e-4) | −0.024 | −0.083 to +0.010 | 0.25 | 0.26 | 0.09 |
| MRE11, Q4 vs Q1 + keratin | 34 | +0.169 (0.13) | −0.368 (0.020) | −0.352 (6.3e-4) | −0.062 | −0.144 to +0.047 | 0.29 | 0.30 | 0.18 |
| Hallmark DNA repair (not NHEJ) | 64 | +0.037 (0.036) | −1.028 (0.080) | −0.161 (0.042) | −0.038 | −0.110 to +0.001 | 0.062 | 0.077 | 0.24 |

Reactome NHEJ plus a keratin covariate is the strongest NHEJ-module result. The point estimate is about 9% of the keratin-adjusted CLDN4 → IFN effect. Both the case bootstrap and pingouin put the interval across 0. MRE11 is the strongest single KEGG gene in the same direction (about 18% of the Q4 vs Q1 IFN effect); path b is p = 0.020 and path a is not, and the indirect-effect interval includes 0. Within-cohort fits do not concentrate this in one study.

Hallmark DNA repair is the closest indirect effect in the whole grid (boot p = 0.062; 97% of bootstrap draws negative; upper confidence limit +0.001). It is not an NHEJ module. Of the KEGG NHEJ genes, only FEN1 and POLL sit in that hallmark set. The CLDN4 correlation inside the hallmark is carried by CANT1, AK1, SURF1, and RNA-polymerase subunits, not by Ku, DNA-PKcs, or LIG4.

The smallest case-bootstrap p in the NHEJ search was 0.046, for XRCC6 → IFN-rank in the ADC subset with a library-size covariate. The total CLDN4 effect in that specification is +0.002 (p = 0.84), so the indirect effect is not carrying an IFN decrease. Pingouin p = 0.056. That row is the minimum of the search, not evidence for the path.

Keratin adjustment does change the total effect: CLDN4 %pos → IFN goes from −0.161 (p = 0.042) to −0.266 (p = 5.8×10⁻⁴) after KRT8/18/19/EPCAM. The NHEJ piece of that larger total effect stays small.

## What this is not

- Not a new T/NK estimate and not a re-audit of ρ = −0.53.
- Not cell-level high vs low inside a tumor. The public inputs are
  malignant UMI-sums, one column per donor/sample/patient.
- Not evidence that CLDN4, NHEJ, or cGAS-STING causes the IFN change.
- Not a serial CLDN4 → NHEJ → cGAS-STING → IFN structural model.
- Not GSE148071 / GSE127465 / GSE207422 / GSE154826.
- Not a dual-high TACSTD2 score.
- Genome-wide FDR is not the claim. The pre-specified mediation is the
  KEGG model. The 16,560-specification search is a sensitivity analysis;
  its smallest p-value is not a confirmatory test.
- Not a rebranding of Hallmark DNA repair as NHEJ.

## Files

- `tables/module_q4q1.tsv` — high vs low and continuous module contrasts
- `tables/module_spearman.tsv` — cohort ρ and DL pool
- `tables/mediation.tsv` — pre-specified CLDN4 → NHEJ → IFN and the sensitivities
- `tables/mediation_search.tsv` — all 16,560 specifications (Sobel)
- `tables/mediation_strongest.tsv` — case bootstrap and pingouin on the strongest rows
- `figures/mediation_search_forest.png`
- `tables/module_genes.tsv` — gene-level Q4 vs Q1
- `tables/patient_module_scores.tsv` — one row per unit
- `tables/module_membership.tsv` — which genes entered each score
- `figures/` — forests, boxes, scatters, path, bootstrap

Reproduce:

```bash
python3 methods/concordant4_cldn4_nhej_sting_ifn/analyze.py
python3 methods/concordant4_cldn4_nhej_sting_ifn/search_mediation.py
```

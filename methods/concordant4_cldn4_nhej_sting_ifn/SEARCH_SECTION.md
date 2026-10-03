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

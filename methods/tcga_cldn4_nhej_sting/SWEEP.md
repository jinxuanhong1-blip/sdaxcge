# Sweep — does any pre-declared setting put CLDN4-low toward NHEJ down and IFN/STING/APM up?

Same tumors as `analyze.py` (TCGA-LUAD n=515, TCGA-LUSC n=501, plus the joint NSCLC matrix). Every number below is computed in `sweep.py`. The pre-specified ssGSEA result in `FINDING.md` is unchanged: CLDN4-low has a **higher** KEGG NHEJ score, and IFN-γ / STING / APM are not higher.

**No quadruple in this grid qualifies.** Qualifying was defined before the run as: at least 3 of 4 scores in the thesis direction, at least 2 of those with raw p<0.05, the smallest p in the quadruple also in that direction, and more significant thesis hits than significant opposite hits. Grid size: 5400 score contrasts, 14400 quadruples, 0 qualifying. Bonferroni 0.05/5400 = 9.259e-06.

## Best quadruple on the declared sort

The row below is the maximum sign agreement (then thesis-direction p<0.05 counts, then fewer opposite-direction p<0.05 hits). It is the closest grid point, not a confirmed result.

- Setting: **LUSC**, stratum **krt_high**, method **gsva**, contrast **mw_d1d10**, n=50
- Sets: NHEJ `NHEJ_KEGG`, IFN `IFNg`, STING `STING_reactome`, APM `APM_nonIFN`
- Thesis signs: 4/4. Significant in the thesis direction: 0. Significant the other way: 0. Smallest p is APM `APM_nonIFN` p=0.461 (thesis direction).
- Effects (Spearman ρ, or median(low)−median(high) for a cutoff test): NHEJ -0.069 (p=0.655, match=True); IFN +0.070 (p=0.497, match=True); STING +0.013 (p=0.574, match=True); APM +0.112 (p=0.461, match=True).
- None of the four arms in that closest setting has p<0.05. Matching signs with null tests are not a thesis hit.

## Significant hits by family (raw p<0.05, whole grid)

- **NHEJ:** 0 contrasts in the thesis direction, 1449 in the opposite direction, out of 1800.
- **IFN:** 0 contrasts in the thesis direction, 457 in the opposite direction, out of 1200.
- **STING:** 79 contrasts in the thesis direction, 421 in the opposite direction, out of 1200.
- **APM:** 210 contrasts in the thesis direction, 115 in the opposite direction, out of 1200.

NHEJ has no significant thesis-direction contrast. IFN has none either. STING and APM have some, described below, and they do not co-occur with a significant NHEJ decrease.

Pooled NSCLC can mix the two histologies. Within LUAD, STING-core AUCell vs CLDN4 is ρ=+0.047 (not the thesis direction). The very small pooled AUCell p-values are not a within-LUAD result. Inside one histology, the STING-core AUCell association that does go up in CLDN4-low is LUSC keratin-high, ρ=−0.202, p=0.001, and the NHEJ arm in that same slice is not a significant decrease.

The APM signal that is real inside one histology is mostly the 8 genes outside Hallmark IFN-γ, and it is strongest in immune-high LUAD (z-mean partial on keratin ρ=−0.314, p=2.8×10⁻⁷). The full 19-gene APM score in that same slice, GSVA partial on keratin, is ρ=−0.200, p=0.001. IFN-γ itself has no significant thesis-direction contrast anywhere in the grid, so this is not an IFN-high CLDN4-low state.

## Best core quadruple (stratum = all, unadjusted Spearman, KEGG or ligation-core NHEJ, Hallmark IFN-γ, full APM)

LUSC / aucell / `NHEJ_core` + `IFNg` + `STING_core` + `APM`: signs 2/4, thesis p<0.05 0, opposite p<0.05 1. NHEJ ρ=-0.208 (p=2.67e-06); IFN ρ=+0.004 (p=0.923); STING ρ=-0.061 (p=0.171); APM ρ=-0.018 (p=0.694).

## How often the thesis sign appears (unadjusted Spearman, all tumors)

| Gene set | LUAD matches / tests | LUSC matches / tests | NSCLC matches / tests |
|---|---:|---:|---:|
| NHEJ_KEGG | 0/4 | 0/4 | 0/4 |
| NHEJ_core | 0/4 | 0/4 | 0/4 |
| NHEJ_reactome | 0/4 | 0/4 | 0/4 |
| IFNg | 0/4 | 0/4 | 0/4 |
| IFNa | 0/4 | 0/4 | 0/4 |
| STING_reactome | 0/4 | 0/4 | 1/4 |
| STING_core | 0/4 | 1/4 | 1/4 |
| APM | 4/4 | 3/4 | 0/4 |
| APM_nonIFN | 4/4 | 4/4 | 1/4 |

Four methods are the denominator in each cell. A cell of 0/4 means every method has the opposite sign.

## Named genes (method-invariant)

PRKDC, LIG4, and STING1 (TMEM173) are single rows on the matrix. ssGSEA, GSVA, and AUCell do not change their correlations. Thesis direction for PRKDC/LIG4 is positive ρ (CLDN4-low, lower NHEJ gene). Thesis direction for TMEM173 and HLA is negative ρ.

| Cohort | Gene | Spearman ρ | p | Any contrast with p<0.05 in the thesis direction? |
|---|---|---:|---:|---|
| LUAD | PRKDC | -0.126 | 0.004 | no |
| LUAD | LIG4 | -0.168 | 1.28e-04 | no |
| LUAD | TMEM173 | +0.193 | 1.04e-05 | no |
| LUAD | HLA-A | +0.059 | 0.180 | no |
| LUAD | HLA-B | +0.035 | 0.424 | no |
| LUAD | HLA-C | +0.023 | 0.596 | no |
| LUSC | PRKDC | -0.229 | 2.09e-07 | no |
| LUSC | LIG4 | -0.082 | 0.067 | no |
| LUSC | TMEM173 | -0.001 | 0.978 | yes |
| LUSC | HLA-A | +0.032 | 0.474 | no |
| LUSC | HLA-B | +0.007 | 0.880 | no |
| LUSC | HLA-C | +0.023 | 0.610 | no |

PRKDC does not match the thesis at any cutoff or stratum in either histology. LIG4 does not in LUAD. STING1 (TMEM173) does not in LUAD (ρ stays positive). The LUSC STING1 “yes” is only the immune-low half (ρ=−0.140, p=0.027); the full LUSC cohort is null (ρ=−0.001). One adjacent kinase that was not a primary endpoint does: **TBK1** in LUAD, ρ=−0.228, p=1.8×10⁻⁷, including every CLDN4 cutoff tested. That is higher TBK1 in CLDN4-low LUAD. It does not pull the STING gene-set score with it, because TMEM173 in the same tumors goes the other way (ρ=+0.193).

## What was not done

- Gene sets were not rebuilt by picking members that already correlated with CLDN4 in the thesis direction.
- Samples were not dropped one at a time to chase a p-value.
- Signs were not flipped after the fact. A positive `thesis_stat` is the only match flag.
- KEGG NHEJ unadjusted Spearman matches the thesis in 0/12 method×cohort tests (stratum = all).

## Set sizes on the LUAD matrix

- `APM`: 19 genes (APM)
- `APM_nonIFN`: 8 genes (APM)
- `IFNa`: 96 genes (IFN)
- `IFNg`: 200 genes (IFN)
- `NHEJ_KEGG`: 13 genes (NHEJ)
- `NHEJ_core`: 7 genes (NHEJ)
- `NHEJ_reactome`: 29 genes (NHEJ)
- `STING_core`: 5 genes (STING)
- `STING_reactome`: 11 genes (STING)

## Files

- `tables/sweep_contrasts.tsv` — every score contrast
- `tables/sweep_quadruples.tsv` — every four-score setting, best first
- `tables/gene_cutoffs.tsv` — PRKDC, LIG4, STING1, HLA and the other ligation genes
- `figures/sweep_spearman.png` — unadjusted ρ by method
- Reproduce: `python3 methods/tcga_cldn4_nhej_sting/sweep.py`


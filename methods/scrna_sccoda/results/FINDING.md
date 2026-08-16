# FINDING — compositional analysis (scCODA fallback)

Public only. Patient is the unit. scCODA was **not installable**
(pip `sccoda` requires rpy2 → system R). Fallback: CLR + ILR on
6-part lineage fractions, plus Dirichlet-multinomial two-group LRT
with permutation p.

User trend tested

1. **H1** T/NK fraction drops when patient malignant TACSTD2 is high.
2. **H2** NMPR has higher epithelial / lower T/NK.

## Combination that matches the trend

**GSE241934_REAL** is the only slice whose *directions* match both H1 and H2.

- H1 (CLR T/NK vs TACSTD2): n=24, ρ=-0.018, p=0.933
- H2: NMPR vs MPR 21 vs 13; epi one-sided p=0.101; T/NK one-sided p=0.262
- All direction-matching slices: GSE241934_REAL

Matching means **sign of the effect**, not p<0.05. On this winning
slice H1 is a null correlation (ρ ≈ 0). Fraction Spearman is also
null (ρ = −0.042, p = 0.85, n = 24). Median T/NK fraction:
TACSTD2-high 0.395 vs low 0.476 (CLR high vs low one-sided p = 0.33).
Dirichlet-multinomial high vs low: LRT p_perm = 0.20.

H2 on the same slice: epithelial fraction median NMPR 0.016 vs MPR
0.0065 (one-sided p = 0.10; two-sided p = 0.20). T/NK median 0.471
vs 0.482 (one-sided p = 0.26). Median leftover epithelial *cells*
are 96 (NMPR) vs 17 (MPR) — the residual-tumor definition of MPR.

No two- or three-series combination keeps both directions: adding
GSE207422 flips H1 to positive. The three-series pool is the opposite
of H1 (CLR ρ = +0.36, p = 0.0086, n = 53).

## Per-cohort honest n / p (primary)

| Cohort | H1 n | CLR ρ (TACSTD2 vs T/NK) | ρ p | H1 dir | H2 n (NMPR vs MPR) | epi p (NMPR>) | T/NK p (NMPR<) | H2 dir |
|---|---:|---:|---:|---|---:|---:|---:|---|
| GSE207422_post | 12 | 0.287 | 0.366 | no | 8 vs 4 | 0.659 | 0.467 | no |
| GSE241934_IIT | 11 | -0.064 | 0.853 | yes | 7 vs 4 | 0.885 | 0.606 | no |
| GSE241934_REAL | 24 | -0.018 | 0.933 | yes | 21 vs 13 | 0.101 | 0.262 | yes |
| GSE241934 | 35 | 0.011 | 0.949 | no | 28 vs 17 | 0.203 | 0.367 | no |
| GSE291670 | 6 | -0.429 | 0.397 | yes | 3 vs 3 | 0.650 | 1.000 | no |
| GSE207422_post+GSE241934 | 47 | 0.082 | 0.586 | no | 36 vs 21 | 0.226 | 0.319 | yes |
| GSE207422_post+GSE241934_REAL | 36 | 0.061 | 0.724 | no | 29 vs 17 | 0.128 | 0.255 | yes |
| GSE207422_post+GSE291670 | 18 | 0.738 | 0.000473 | no | 11 vs 7 | 0.761 | 0.813 | no |
| GSE241934+GSE291670 | 41 | 0.381 | 0.014 | no | 31 vs 20 | 0.318 | 0.550 | no |
| all_three_post | 53 | 0.358 | 0.00858 | no | 39 vs 24 | 0.338 | 0.492 | yes |

## Separate-cohort takeaways

- **GSE207422 post-tx (n=12; H2 8 vs 4).** H1 opposite (CLR ρ = +0.29,
  p = 0.37). H2 opposite for epithelium (NMPR median fraction 0.055 vs
  MPR 0.076). DM NMPR vs MPR p_perm = 0.62.
- **GSE241934 IIT (n=11; 7 vs 4).** H1 direction only (ρ = −0.064,
  p = 0.85). H2 opposite: NMPR has *lower* epithelial fraction and
  *higher* T/NK. IIT MPR leftovers still have more epi cells (median
  202 vs 63).
- **GSE291670 (n=6; 3 vs 3).** Strongest H1 *direction* (ρ = −0.43,
  p = 0.40). H2 is the other way: NMPR epithelial 0.25 vs MPR 0.28;
  T/NK 0.098 vs 0.027 (NMPR *higher* T/NK; two-sided p = 0.10, the
  smallest exact p at 3 vs 3). TACSTD2 itself is higher in MPR.

## Confounds that stay in the write-up

- Post-neoadjuvant **MPR/pCR residuals have fewer epithelial cells by definition**.
  H2 epithelial enrichment in NMPR is partly that residual-tumor fact,
  especially in GSE241934.
- GSE207422 / GSE291670 use a marker hierarchy (CopyKAT / author barcodes
  are not on GEO). GSE241934 uses author `major.cell.type`.
- GSE291670 is n=3 vs 3. Exact Wilcoxon cannot go below p=0.10 two-sided.
- Cell-level p-values are not reported (pseudoreplication).
- scCODA credible effects are **not** claimed.

## Files

- `composition_all.tsv` — patient-level counts + TACSTD2
- `tests_summary.tsv` / `tests_full.json`
- `figures/`

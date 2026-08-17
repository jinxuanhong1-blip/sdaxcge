# Methods — winning-pair propeller + CellPhoneDB-style (CLDN4 only)

ADDITIVE. **CLDN4 only.** TACSTD2 does not define groups. Dual-high is not
built. GSE207422 is not run. LIANA and CellChat are not run.

## Winning pair (given)

PR #320 author-malignant **%pos** vs T/NK:

- Spearman: GSE131907 + GSE205335, k=2, N=43, ρ=−0.479, p=0.00152, I²=0%
- Q4 vs Q1: n=23 compared (12/11), r=−0.705, p=0.000301, I²=0%

This recut does not re-rank that combo. It asks two new questions on the
**same locked units**.

## Locked units

| Cohort | Source | Eligibility | Unit |
|---|---|---|---|
| GSE131907 | PR #279 / #320 `GSE131907_samples.tsv` | tumor origins `{tLung, tL/B, mLN, PE, mBrain}` and `n_malignant ≥ 20` | **sample** |
| GSE205335 | PR #279 / #320 `GSE205335_patients.tsv` | ≥20 author-malignant and ≥20 T/NK | **patient** |

GSE131907 unique patient n is reported beside sample n. Quartiles are cut
**within cohort** on malignant CLDN4 %pos with `pd.qcut` on average ranks
(same rule as `methods/cldn4_malig_q4_tnk`). Q4 vs Q1 uses the tails only.

Author labels:

- GSE131907 malignant = `Cell_subtype ∈ {Malignant cells, tS1, tS2, tS3}` on tumor-origin samples. T/NK = `T lymphocytes` + `NK cells`. B = `B lymphocytes`.
- GSE205335 malignant = `lineage.sub == Malignant cells`. T/NK = `lineage.total == T/NK cells`. B = B+plasma (`n_b_plasma` in the locked extract). Normal-tissue GEO samples are dropped.

## 1. Propeller / speckle-style composition

Phipson et al., *Bioinformatics* 2022 (`speckle::propeller`):

1. Cell-type fraction per unit (T/NK, B, CD8).
2. Empirical logit with a 0.5 count offset: `log((c+0.5)/(n−c+0.5))`.
3. Two-group test of Q4 vs Q1 on the transformed values.
4. BH on the pre-specified primary family: T/NK + B × 2 cohorts (4 tests).

R / limma is not available. The engine is therefore:

- **Primary p:** Welch two-sample t on logit (per cohort).
- **Pair companion:** OLS `logit ~ Q4 + cohort` on the stacked tails.
- **eBayes:** Smyth 2004 inverse-χ² squeeze across the 3 compartments in a scope. With k=3 the prior is weakly identified; moderated p is a companion, not the primary p.
- **asin-sqrt** and raw-fraction Mann–Whitney (PR #320 rank-biserial) are companions.

Recovery: logit Δ (Q4 − Q1) < 0 and BH q < 0.10 on the primary 4-test family.

## 2. CellPhoneDB-style documented mean score

CellPhoneDB v5 pair list (`ventolab/cellphonedb-data`) plus a curated MHC-I /
T-recruit overlay. Hyphenated symbols (HLA-A) are split only when both sides
match the gene list. Complexes use `+` subunits.

For sender group S (malignant) and receiver group R (same-unit T/NK):

1. Expression = log1p(CP10k).
2. Partner mean = **minimum of subunit means** (CellPhoneDB complex rule).
3. Pair score = **arithmetic mean of the two partner means** (Efremova 2020; Garcia-Alonso 2022).
4. `pass_expr_prop` if both partners are detected in ≥10% of cells in their group.

Inferential contrast is **between-unit Q4 vs Q1** (Mann–Whitney on the
per-unit score). This is not the within-patient CLDN4-high vs low split used
by the LIANA agents. A pair is listed if it is detected in ≥3 Q1 and ≥3 Q4
units. BH is within each scope × cohort focus family.

Focus pathways: MHC-I and T-recruit. Other CPDB pairs are not the table.

This is **not** a CellChat Hill probability and **not** LIANA `mt.cellphonedb`.

## Software

Python 3 (numpy, pandas, scipy, matplotlib, statsmodels, PyYAML, rdata).
No R. Matrices are not stored in git.

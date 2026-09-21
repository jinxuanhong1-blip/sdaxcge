# Concordant-4 CD8/NK state vs malignant CLDN4

Patient-level test of effector intensity and exhaustion inside CD8 cells and inside NK cells, in the locked concordant-4 units (GSE123902, GSE131907, GSE205335, GSE189357).

This is the scRNA composition analog of the CosMx neighbor-state result. CosMx (He 2022) already says CLDN4-high tumor cells have fewer cytotoxic neighbors, and the neighbors that remain are not lower for GZMB, PRF1, NKG7, or IFNG (hi/lo 1.11–1.22, 0/8 down). scRNA has no coordinates. The question here is whether the CD8 and NK cells that are present in CLDN4-high patients have a different effector or exhaustion intensity than those in CLDN4-low patients.

The locked T/NK fraction result is not re-derived. Malignant CLDN4 % positive and the within-cohort quartile are taken from that patient table. A lower T/NK fraction is exclusion. It is not, by itself, a change in the state of the cells that remain.

## Exposure QC

Recomputed malignant CLDN4 % positive versus the locked value (same gates: author malignant in GSE131907 and GSE205335; marker malignant in GSE123902 and GSE189357). The tests use the locked value.

- GSE123902: Pearson 1.0000, median absolute difference 0.000 percentage points
- GSE131907: Pearson 1.0000, median absolute difference 0.000 percentage points
- GSE189357: Pearson 1.0000, median absolute difference 0.000 percentage points
- GSE205335: Pearson 1.0000, median absolute difference 0.000 percentage points

## Primary result

Unit = patient / donor / sample. Score = mean log1p(CP10k) of the named genes, averaged over cells in the compartment, then averaged over genes. Effector genes are the CosMx set GZMB, PRF1, NKG7, IFNG. Exhaustion genes are PDCD1, HAVCR2, LAG3, TIGIT, TOX. A unit enters a test when that compartment has at least 10 cells. Meta-analysis is DerSimonian–Laird on Fisher z of within-cohort Spearman ρ, cohorts with n≥5. Q4 versus Q1 uses the locked within-cohort quartiles. OLS is score ~ cohort + Q4 among Q1 and Q4 units only (cohort-adjusted difference, Q4 minus Q1). Stacked rank-biserial pools those units across cohorts.

Muzzling, if it were operating at the patient level, would be a negative effector association and a positive exhaustion association. Four tests, two-sided, descriptive.

| compartment | score | k | N | ρ | 95% CI | p | I² | cohorts ρ<0 | OLS Δ | OLS p | Q1/Q4 n | stacked r | stacked p |
|---|---|---:|---:|---:|---|---:|---:|---|---:|---:|---|---:|---:|
| cd8 effector_cosmx | 4 | 61 | 0.052 | -0.395 to 0.480 | 0.828 | 61.2% | 2/4 | -0.141 | 0.374 | 19/14 | -0.165 | 0.434 |
| nk effector_cosmx | 4 | 58 | 0.139 | -0.148 to 0.404 | 0.343 | 0.0% | 1/4 | 0.012 | 0.948 | 17/12 | 0.265 | 0.241 |
| cd8 exhaustion | 4 | 61 | -0.197 | -0.446 to 0.081 | 0.163 | 0.0% | 4/4 | -0.152 | 0.090 | 19/14 | -0.301 | 0.150 |
| nk exhaustion | 4 | 58 | -0.095 | -0.367 to 0.191 | 0.517 | 0.0% | 3/4 | -0.021 | 0.433 | 17/12 | -0.078 | 0.740 |

Cohort rows:

- cd8|ALL|effector_cosmx|min10: GSE123902 n=11 ρ=0.018 (p=0.958) ΔQ=-0.225; GSE131907 n=19 ρ=-0.067 (p=0.786) ΔQ=-0.087; GSE205335 n=22 ρ=-0.347 (p=0.113) ΔQ=-0.424; GSE189357 n=9 ρ=0.733 (p=0.025) ΔQ=0.557
- nk|ALL|effector_cosmx|min10: GSE123902 n=11 ρ=0.427 (p=0.190) ΔQ=0.684; GSE131907 n=18 ρ=0.032 (p=0.900) ΔQ=-0.258; GSE205335 n=20 ρ=0.140 (p=0.556) ΔQ=0.112; GSE189357 n=9 ρ=-0.017 (p=0.966) ΔQ=-0.146
- cd8|ALL|exhaustion|min10: GSE123902 n=11 ρ=-0.045 (p=0.894) ΔQ=-0.091; GSE131907 n=19 ρ=-0.056 (p=0.819) ΔQ=-0.104; GSE205335 n=22 ρ=-0.346 (p=0.115) ΔQ=-0.262; GSE189357 n=9 ρ=-0.267 (p=0.488) ΔQ=-0.041
- nk|ALL|exhaustion|min10: GSE123902 n=11 ρ=0.009 (p=0.979) ΔQ=-0.014; GSE131907 n=18 ρ=-0.020 (p=0.938) ΔQ=-0.027; GSE205335 n=20 ρ=-0.202 (p=0.394) ΔQ=-0.025; GSE189357 n=9 ρ=-0.117 (p=0.765) ΔQ=-0.007

## Verdict

Patient-level scRNA does not show a drop in effector intensity in CLDN4-high patients, in CD8 cells or in NK cells. It also does not show higher exhaustion. The same meta code, run on the locked T/NK fraction, recovers ρ = −0.531 (p = 1.65×10⁻⁵) and stacked rank-biserial r = −0.724 (n = 19/16). A fraction-sized association would have been visible. These state tests are nulls.

CD8 effector ρ = 0.052 (95% CI -0.395 to 0.480), p = 0.828, N = 61, I² = 61%. Cohort-adjusted Q4 minus Q1 = -0.141 (p = 0.374). The I² is GSE189357 (n = 9, ρ = 0.733, p = 0.025) against GSE205335 (n = 22, ρ = −0.347, p = 0.11). Leave-one-out stays non-significant, including after dropping GSE189357 (ρ = −0.180, p = 0.23). Mean Q4/Q1 ratios for this score are 0.83, 0.95, 0.72, and 1.69. That is not a consistent decrease.

NK effector ρ = 0.139 (95% CI -0.148 to 0.404), p = 0.343, N = 58, I² = 0%. The point estimate is slightly higher in CLDN4-high patients, not lower. OLS Δ = 0.012 (p = 0.948).

CD8 exhaustion ρ = -0.197 (95% CI -0.446 to 0.081), p = 0.163, I² = 0%, and 4/4 cohorts have ρ < 0. The direction is less exhaustion, which is the opposite of a muzzling increase. Dropping GSE205335, the cohort with the largest negative ρ (−0.346, p = 0.11), leaves ρ = −0.096 (p = 0.60). NK exhaustion ρ = -0.095, p = 0.517.

Units with fewer than 10 cells in the compartment are out of that test. CD8 keeps 61 of 65. NK keeps 58 of 65. Dropped CD8 units: LX653, LX681, NS_07, NS_12. Dropped NK units: LX653, LX684, EBUS_13, EBUS_28, NS_12, P1090, P1115.

This does not confirm the CosMx neighbor ratios, and it does not contradict them. CosMx scored cells next to a CLDN4-high tumor cell. This scores the CD8 or NK compartment of a CLDN4-high patient. A null here means the patient-level intensity test does not add a muzzling claim on top of the locked exclusion result (fewer T/NK).

Do not quote cell counts as n. Do not add GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, or GSE207422. Do not read a Visium same-spot correlation into this table.

## Sensitivities

Same meta machinery. `effector_residual_on_naive` is the within-cohort residual of the CD8 effector score on the CD8 naive/memory score (IL7R, TCF7, CCR7, SELL, LEF1), so a naive-versus-effector mix is not mistaken for a change inside a fixed program. `GNLY_POS` restricts CD8 to GNLY-positive cells (GNLY is not in the four-gene score) and is the closer analog of intensity inside cells that already look cytotoxic. `frac_nk_among_cd8nk` is lineage composition, not intensity. The Spearman metas all include 0. The closest quartile contrast is the naive-residualized CD8 effector score (OLS p = 0.053); its Spearman meta is ρ = −0.14, p = 0.55.

Gene-level meta-analyses of GZMB, PRF1, NKG7, and IFNG, in CD8 and in NK, all have confidence intervals that include 0. CD8 TOX is the strongest single-gene lean (ρ = -0.261, p = 0.061; Q4−Q1 OLS Δ = -0.113, p = 0.012). The direction is lower TOX, not higher, and it is one gene in a long sensitivity list. It is not a primary result.

Author cell-state names are not pooled across cohorts. In GSE205335 the fraction of CD8 cells labeled CD8+ MT high tracks malignant CLDN4 (ρ = 0.609, p = 0.003, n = 21). Effector intensity inside that label does not (ρ = −0.025). That fraction is one label among many author slices, and it is composition inside CD8, not the intensity test. A transitional-NK intensity Spearman in the same cohort (ρ = 0.83, n = 12) is a small subset whose Q4 arm is a single patient; the NK-wide effector score there is ρ = 0.14 (p = 0.56).

| test | k | N | ρ | p | 95% CI | OLS Δ (p) |
|---|---:|---:|---:|---:|---|---|
| cd8|ALL|naive_memory|min10 | 4 | 61 | -0.106 | 0.690 | -0.558 to 0.394 | -0.064 (0.482) |
| cd8|ALL|effector_residual_on_naive|min10 | 4 | 61 | -0.137 | 0.547 | -0.528 to 0.301 | -0.247 (0.053) |
| cd8nk|ALL|effector_cosmx|min10 | 4 | 65 | -0.027 | 0.864 | -0.325 to 0.275 | -0.166 (0.297) |
| cd8nk|ALL|effector_residual_on_frac_nk|min10 | 4 | 65 | -0.042 | 0.782 | -0.325 to 0.249 | -0.195 (0.176) |
| composition|frac_nk_among_cd8nk|min10 | 4 | 65 | 0.195 | 0.150 | -0.071 to 0.436 | 0.116 (0.149) |
| cd8|GNLY_POS|effector_cosmx|min10 | 4 | 51 | -0.042 | 0.833 | -0.408 to 0.336 | 0.029 (0.811) |
| cd8|UMI500|effector_cosmx|min10 | 4 | 61 | 0.052 | 0.828 | -0.395 to 0.480 | -0.137 (0.384) |
| nk|UMI500|effector_cosmx|min10 | 4 | 58 | 0.147 | 0.317 | -0.140 to 0.411 | 0.036 (0.845) |
| cd8|ALL|effector_no_nkg7|min10 | 4 | 61 | 0.067 | 0.804 | -0.430 to 0.532 | -0.143 (0.336) |
| nk|ALL|effector_no_nkg7|min10 | 4 | 58 | 0.130 | 0.375 | -0.157 to 0.397 | 0.046 (0.748) |
| cd8|ALL|effector_cosmx|min30 | 4 | 54 | 0.061 | 0.788 | -0.368 to 0.469 | -0.018 (0.919) |
| nk|ALL|effector_cosmx|min30 | 4 | 48 | 0.110 | 0.509 | -0.213 to 0.411 | -0.019 (0.924) |

## Reproduce

```bash
bash methods/concordant4_cd8nk_state/scripts/download.sh /tmp/concordant4_raw
# GSE205335 RDS is double-gzipped
gzip -dc /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz | gzip -dc > /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds
Rscript methods/concordant4_cd8nk_state/scripts/export_gse205335.R /tmp/concordant4_raw /tmp/concordant4_state/gse205335_cells.tsv.gz
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE123902
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE189357
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE131907
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE205335
python3 methods/concordant4_cd8nk_state/scripts/analyze_state.py
```

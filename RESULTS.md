# TISMO ICB paired pre/post: study-level REML meta-analysis

Public TISMO Gene-module tables (Zeng et al., *Nucleic Acids Research* 2022, PMID 34534350). Expression is the served log2(TPM+1) scale. Every count below is recomputed in `results/tismo_icb_study_meta/tables/`.

## Locked reference (not replaced)

The pairing unit stays the TISMO slice: one cell-line × study × condition × ICB-regimen label. Baseline is `Baseline == 1`. ICB is `Baseline == 0` (responders and non-responders pooled). Delta is mean(ICB) − mean(baseline). Cldn4 and the TJ score are locked to the same 64 Tacstd2 slices. The TJ score is the per-sample mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, and Ocln (mean 6.99 of 7 genes present; minimum 5).

| Marker | Up | Down | Tie | Wilcoxon *p* | Mean Δ |
|---|---:|---:|---:|---:|---:|
| Tacstd2 | **49/64** | 15 | 0 | 5.84e-05 | +0.265 |
| Cldn4 | **34/64** | 28 | 2 | 0.0767 | +0.129 |
| TJ | **37/64** | 27 | 0 | 0.0869 | +0.048 |

**Tacstd2 49/64 up, Wilcoxon *p* = 5.84e-05**, is the reference sentence. Cldn4 is weaker on that same lock: **34/64 up, 28 down, 2 ties, *p* = 0.0767**. TJ is **37/64**, *p* = 0.0869. This meta-analysis does not revise those three rows.

## Why pool by study

The 64 slices sit in **22 studies**. GSE124821 alone is 19/64 slices, so a slice-level Wilcoxon treats many contrasts from one mammary experiment as independent. Fourteen sample IDs in GSE130472 are shared baselines across anti-PD-L1 and anti-CTLA4 arms (old and young 4T1). No other study reuses a sample ID across slices.

Primary study effect: unweighted mean of the locked slice deltas in that accession, so each slice still counts once inside its study. The mouse-level sampling variance treats mice as independent, uses a within-arm pooled residual variance, and adds a covariance when two slices share mice (GSE130472 only). That technical variance is not the primary weight. Pooled slice-replication variances are Tacstd2 0.460, Cldn4 0.526, TJ 0.109 (studies with at least two slices, denominator Σ(k−1)). Primary vᵢ is the larger of the technical variance and (slice-replication variance / kᵢ). An unexpressed gene has a technical variance near zero and a delta near zero; without the floor it would dominate the meta-analysis. That is what happens to Cldn4 under technical weights alone: GSE109485 carries 99.8% of the weight. That fit is reported as a sensitivity and is not the estimate used below.

Across studies the model is normal-normal, yᵢ ~ N(μ, vᵢ + τ²), with τ² fit by REML (Viechtbauer 2005). The reported interval is a modified Knapp–Hartung interval: t critical value on k−1 degrees of freedom, and the Hartung scale is not allowed to shrink below the Wald standard error. Raw Knapp–Hartung, DerSimonian–Laird, technical-SE weights, slice-variance-only weights, within-study GLS, and ignoring shared mice are sensitivities. Leave-one-study-out refits the primary REML model 22 times.

## Study-level sign count

Each study contributes the sign of its mean slice delta. This is the 49/64 question asked once per study.

| Marker | Studies up | Down | Tie | Binomial *p* | Wilcoxon *p* | Median Δ | Mean Δ (Student 95% CI, *p*) |
|---|---:|---:|---:|---:|---:|---:|---|
| Tacstd2 | **18/22** | 4 | 0 | 0.00434 | 0.00415 | +0.117 | +0.225 (+0.042 to +0.408, *p* = 0.0186) |
| Cldn4 | **12/22** | 10 | 0 | 0.832 | 0.0984 | +0.014 | +0.111 (+0.004 to +0.219, *p* = 0.0434) |
| TJ | **11/22** | 11 | 0 | 1 | 0.276 | +0.019 | +0.045 (-0.014 to +0.105, *p* = 0.129) |

Tacstd2 is up in **18/22 studies** (binomial *p* = 0.00434; Wilcoxon *p* = 0.00415). Cldn4 is up in **12/22** (binomial *p* = 0.832; Wilcoxon *p* = 0.0984). TJ is up in **11/22** (*p* = 1). Cldn4’s equal-weight mean (+0.111, Student *p* = 0.0434) sits above its median (+0.014) because a few studies have large positive deltas. The rank test and the sign count do not call that a consistent increase. Tacstd2 is positive on the sign count, the rank test, and the mean.

## REML mean difference

| Marker | μ | Wald 95% CI | Wald *p* | mKH 95% CI | mKH *p* | τ² | I² | Q *p* | 95% PI | max weight |
|---|---:|---|---:|---|---:|---:|---:|---:|---|---|
| Tacstd2 | +0.267 | +0.047 to +0.487 | 0.0174 | +0.034 to +0.501 | 0.0269 | 0.075 | 16.1% | 0.246 | -0.348 to +0.882 | GSE124821 12.7% |
| Cldn4 | +0.138 | -0.051 to +0.326 | 0.153 | -0.062 to +0.338 | 0.167 | 0.000 | 0.0% | 1 | -0.062 to +0.338 | GSE124821 20.8% |
| TJ | +0.042 | -0.049 to +0.133 | 0.364 | -0.054 to +0.139 | 0.374 | 0.000 | 0.0% | 0.997 | -0.054 to +0.139 | GSE124821 11.0% |

Tacstd2 REML μ = 0.267 (Wald 95% CI 0.047 to 0.487; modified Knapp–Hartung 0.034 to 0.501). Cldn4 REML μ = 0.138 (Wald 95% CI -0.051 to 0.326; modified Knapp–Hartung -0.062 to 0.338). TJ REML μ = 0.042 (Wald 95% CI -0.049 to 0.133; modified Knapp–Hartung -0.054 to 0.139).

Tacstd2’s modified Knapp–Hartung interval lies above 0 (*p* = 0.0269). I² is 16%. The prediction interval (-0.348 to +0.882) includes decreases, so a new study need not be positive. GSE124821 carries 12.7% of the Tacstd2 random-effects weight.

The Tacstd2 mean is pulled by large study effects. GSE146027 (YTN16) is +1.410 (slice deltas +0.133 to +3.210). GSE149825 (B16) is +1.256, the average of slice deltas +0.003 and +2.508 (plain dual checkpoint blockade versus birinapant plus dual blockade). Both slices stay inside the locked 64.

Cldn4’s pooled mean is smaller (+0.138 versus Tacstd2 +0.267) and its Wald and modified Knapp–Hartung intervals both include 0 (*p* = 0.153 and 0.167). The locked slice test remains 34/64, *p* = 0.0767, and the study sign count remains 12/22. Cldn4 is weaker than Tacstd2 on the locked slice test, the 22-study sign count, the study-level Wilcoxon test, and the primary REML interval. The equal-weight Student interval for the Cldn4 mean just excludes 0; that is the summary the right skew can move, and it is not the locked claim.

TJ follows Cldn4 rather than Tacstd2: μ = +0.042, interval includes 0, study signs 11/22, slice signs 37/64.

Forest plots: `results/tismo_icb_study_meta/figures/forest_Tacstd2.png`, `forest_Cldn4.png`, `forest_TJ.png`, and `forest_panel.png`. Pooled comparison: `pooled_comparison.png`.

## Leave-one-study-out

Refitting REML after dropping each study moves Tacstd2 μ across +0.167 to +0.308, Cldn4 across +0.107 to +0.157, and TJ across +0.030 to +0.056. Full table: `results/tismo_icb_study_meta/tables/loso.tsv`. Figure: `figures/loso_panel.png`.

Tacstd2: leaving out GSE146027 changes whether the Wald interval lies entirely above 0.
Cldn4: no single study flips whether the Wald interval lies entirely above 0 (full interval does not).
TJ: no single study flips whether the Wald interval lies entirely above 0 (full interval does not).

Tacstd2 moves most when GSE146027 is removed (μ +0.267 → +0.167).
Cldn4 moves most when GSE146027 is removed (μ +0.138 → +0.107).
TJ moves most when GSE130472 is removed (μ +0.042 → +0.056).

## Sensitivities

Same 22 studies. μ is the REML pooled mean except the one-vote row, which is the unweighted mean of study effects with a Student interval.

| Marker | Primary | Slice variance only | Technical SE | DL on primary v | One-vote mean | Technical max weight |
|---|---:|---:|---:|---:|---:|---|
| Tacstd2 | +0.267 | +0.268 | +0.162 | +0.271 | +0.225 | GSE109485 7.2% |
| Cldn4 | +0.138 | +0.129 | -0.000 | +0.138 | +0.111 | GSE109485 99.8% |
| TJ | +0.042 | +0.048 | +0.018 | +0.042 | +0.045 | GSE150401 8.7% |

DerSimonian–Laird on the primary variances agrees with REML at the reported precision. Slice-variance-only weights are close to the primary fit. Technical-SE weights are not usable for Cldn4: GSE109485 (B16, Cldn4 at the expression floor, study Δ -0.0000) takes 99.8% of the weight and pulls μ to -0.0000. That number is a measurement-floor artifact. It is not evidence that Cldn4 is unchanged in the studies where it is expressed, and it is not used as the result.

### Cancer groups with at least 3 studies

Lung is one study (LLC) and is not pooled as its own meta-analysis. There is no KL subgroup.

| Marker | Group | Studies up | μ | Wald 95% CI | Wald *p* | I² |
|---|---|---:|---:|---|---:|---:|
| Tacstd2 | Colorectal | 4/5 | +0.152 | -0.275 to +0.580 | 0.485 | 0% |
| Tacstd2 | Mammary | 4/6 | +0.089 | -0.158 to +0.336 | 0.481 | 0% |
| Tacstd2 | Melanoma | 7/7 | +0.322 | -0.033 to +0.677 | 0.0758 | 0% |
| Cldn4 | Colorectal | 2/5 | +0.128 | -0.322 to +0.577 | 0.578 | 0% |
| Cldn4 | Mammary | 2/6 | +0.056 | -0.248 to +0.360 | 0.719 | 0% |
| Cldn4 | Melanoma | 4/7 | +0.133 | -0.246 to +0.513 | 0.491 | 0% |
| TJ | Colorectal | 3/5 | +0.070 | -0.134 to +0.275 | 0.501 | 0% |
| TJ | Mammary | 3/6 | +0.039 | -0.125 to +0.203 | 0.639 | 0% |
| TJ | Melanoma | 3/7 | +0.045 | -0.128 to +0.218 | 0.607 | 0% |

Subgroup Wald intervals entirely above 0: 0. Melanoma Tacstd2 is up in 7/7 studies (μ +0.322, Wald -0.033 to +0.677). The overall Tacstd2 mean is not a within-histology certainty.

## Lung is LLC, not KL

TISMO ICB in this 64-slice lock contains one lung study: **GSE155972, line LLC** (2 slices, 33 mice). Study mean Δ: Tacstd2 +0.429, Cldn4 +0.207, TJ -0.145. LLC is Lewis lung carcinoma. It is not a Kras/Stk11 (KL) or Kras/Trp53 (KP) model.

Catalog check: cellLineMeta has 92 lines; vivoMeta has 345 rows. STK11 mentioned in cellLineMeta: False. LKB1 mentioned: False. Same two strings in vivoMeta: False, False. Lung lines on the catalog: LLC (Lung carcinoma), CMT-167 (Lung carcinoma), MLE12 (Lung adenocarcinoma). KPB25L is Mammary cancer, NOS and stays in the mammary group. KPC is Pancreatic ductal adenocarcinoma and is not in this ICB lock.

No KL or KP row was added. The older label `kl_kp_llc` is not used here, because that stratum was LLC twice.

## What this does not say

- It does not replace Tacstd2 49/64.
- It does not say Cldn4 rises after ICB in the same way Tacstd2 does. The locked slice test (34/64, *p* = 0.077), the 22-study sign count (12/22), and the primary REML interval (includes 0) all leave Cldn4 weaker.
- It does not treat a technical standard error of an unexpressed gene as infinite information. That sensitivity is shown and set aside.
- It does not say the next study will be positive. Prediction intervals cover negative deltas.
- It does not split responders from non-responders. The lock pools them.
- It does not use private 8-KL matrices, and it does not call LLC a KL line.

## Reproduce

```bash
pip install -r methods/tismo_icb_study_meta/requirements.txt
python3 methods/tismo_icb_study_meta/analyze.py
python3 -m unittest tests/test_tismo_icb_study_meta.py
```

Inputs are the vendored TISMO Gene-module CSVs in `data/tismo/vivo/` (same export as the 49/64 lock) plus `cellLineMeta.json` and `vivoMeta.json` for the KL audit.

## Files

| Path | Content |
|---|---|
| `results/tismo_icb_study_meta/tables/locked_reference.tsv` | 49/64, 34/64, 37/64 recomputed |
| `results/tismo_icb_study_meta/tables/slice_effects.tsv` | 64 paired deltas per marker |
| `results/tismo_icb_study_meta/tables/study_effects.tsv` | 22 study effects, variances, weights |
| `results/tismo_icb_study_meta/tables/reml_summary.tsv` | REML and sensitivity fits |
| `results/tismo_icb_study_meta/tables/loso.tsv` | Leave-one-study-out |
| `results/tismo_icb_study_meta/tables/subgroup_reml.tsv` | Cancer groups with ≥3 studies |
| `results/tismo_icb_study_meta/figures/` | Forests, LOSO, pooled comparison |

## Cldn4 specification sweep

Tacstd2 stays at 49/64. This section changes only Cldn4 weights, slice filters, expression-floor rules, and genotype subsets. No slice was added. LLC is still LLC. The full grid is `results/tismo_icb_study_meta/tables/cldn4_sweep.tsv`.

Floor rules use the TISMO log2(TPM+1) arm means. `max_arm` is the larger of the baseline and ICB means. `baseline_ge_1` keeps slices whose naive mean is at least 1 (TPM about 1). ICB-only drops slices whose label also names another modality (radiation, birinapant, regorafenib, TGF-β or GARP blockade, IL-17 blockade, exercise, high-fat diet, BRAF inhibitor, swainsonine). Parental drops KO, knockdown, Setdb1, Brca, Aire, Man2a1, MDK, Apobec, and UV labels. Perturbed is the complement.

Biology grid: **54** specifications (6 floors × 3 regimen rules × 3 genotype rules). Two-sided slice Wilcoxon p < 0.05 in **6** of them. Bonferroni threshold for this grid: 0.0009259. Among specs with at least 10 slices and a positive mean, one-sided greater p < 0.05 in **5**.

Locked cell (no filter): none|all|all: 34/64 up, mean Δ +0.129, median +0.040, Wilcoxon p = 0.0767, one-sided greater p = 0.0383.

Smallest two-sided slice p in that biology grid, requiring a positive mean and at least 10 slices: none|icb_plus_other|all: 8/12 up, mean Δ +0.333, median +0.276, Wilcoxon p = 0.0269, one-sided greater p = 0.0134.

Rows with two-sided slice p < 0.05: 6. Every one of them has regimen `icb_plus_other`. They are not six different findings. Parental and all genotypes are the same slices in this subset, and the max-arm floors of 0.25 and 0.5 keep the same 8 slices. The two unique sets are the 12 ICB-plus-other slices (8 up, p = 0.027) and those 12 after dropping 4 slices whose higher arm mean is below 0.25 (7/8 up, p = 0.016). Those 12 labels are radiation, birinapant, swainsonine, regorafenib, GARP:TGF-β (two Fc designs), TGF-β, high-fat diet, TGF-β trap, exercise, IL-17 blockade, and a BRAF inhibitor, each on top of a checkpoint antibody. REML on the 8-slice set is μ = +0.426 (Wald p = 2.6×10⁻⁴; modified Knapp–Hartung p = 0.019). Neither slice p clears the biology-grid Bonferroni line (0.0009259).

The thesis-aligned regimen is checkpoint blockade without a second modality. That row is none|icb_only|all: 26/52 up, mean Δ +0.082, median +0.000, Wilcoxon p = 0.421, one-sided greater p = 0.21. Study-level tests on those slices stay above 0.05 as well. Adding a floor or a parental-only cut does not push the ICB-only two-sided p under 0.05 in any cell with at least 10 slices.

Where Cldn4 is already on at baseline (baseline mean ≥ 1, all regimens, all genotypes): baseline_ge_1|all|all: 8/22 up, mean Δ +0.053, median -0.088, Wilcoxon p = 0.61, one-sided greater p = 0.706. The increase does not get stronger in the slices where a tight-junction reading is possible. It gets weaker. High-baseline models include 4T1 and KPB25L, and those lines are not a consistent rise.

### Weights on the unfiltered 64

| Weight | Estimate | p |
|---|---|---:|
| Equal slice, Wilcoxon two-sided | 34/64 up, mean +0.129 | 0.0767 |
| Equal slice, Wilcoxon greater | same counts | 0.0383 |
| Equal study, sign test | 12/22 studies up | 0.832 |
| Equal study, Wilcoxon | mean of study means +0.111 | 0.0984 |
| Equal study, Student t on study means | +0.111 | 0.0434 |
| REML, primary variance | μ +0.138 (-0.051 to +0.326) | Wald 0.153 |
| REML, technical SE | μ -0.000 | 0.979 |

The one-sided greater p on the locked 64 is 0.038. That is the same signed-rank statistic as the two-sided lock (p = 0.077), read in the thesis direction. It is not a new sample. The two-sided lock stays the reference, because Tacstd2’s two-sided p is 5.84×10⁻⁵ and Cldn4’s is not in that range. The Student interval on 22 study means (p = 0.043) is the equal-study mean, already reported above, and it is pulled by a right skew (median study Δ near 0).

### What was not adopted

Leaving out one cell line at a time, the smallest two-sided slice p is drop_line:KPB25L: 33/58 up, mean Δ +0.200, median +0.074, Wilcoxon p = 0.00842, one-sided greater p = 0.00421. That specification deletes a negative line. KPB25L is a mammary TISMO model inside the 64-slice lock. Removing it is influence, not a genotype subset and not an expression-floor rule. Dropping 4T1, the other high-Cldn4 mammary line, also moves the p-value down because 4T1’s deltas are mostly negative. Those rows are in the table under `influence_drop_line` and are not the result.

Histology splits are in the same table (`family = histology`). They are small and are not a substitute for the 64-slice lock.

**Sweep result.** The locked Cldn4 sentence stays 34/64, Wilcoxon p = 0.077. Checkpoint-only slices do not show a two-sided increase. The grid’s smallest p-values are checkpoint blockade plus a second named modality, on 8 or 12 slices, and they do not survive a Bonferroni cut of the 54-cell grid. Restricting to slices with baseline Cldn4 already on (mean ≥ 1) makes the paired change weaker (8/22, p = 0.61). Tacstd2 remains 49/64.

Figure: `results/tismo_icb_study_meta/figures/cldn4_sweep.png`.

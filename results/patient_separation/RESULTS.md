# Patient bootstrap of CosMx exclusion and concordant-4 CLDN4–T/NK

Additive layer on the locked public results. The CosMx cytotoxic ratios 0.36 / 0.52 (50 / 100 µm, 8/8 sections, 5/5 patients, sign P = 0.031) are not recomputed as a replacement, and the concordant-4 Spearman is not replaced. Inference here uses patients only.

## How a spec was allowed to win

The grid below was fixed before either outcome was ranked. A spec is eligible only when every patient has a usable contrast. The patient-bootstrap 95% interval of the median contrast must lie entirely below zero (CLDN4-high is the lower arm). Among those specs, the largest absolute median wins. If two medians match, the narrower interval wins.

With five CosMx patients the percentile bootstrap of the median is exactly the range of the five patient values: the upper bound is the weakest patient. That is the whole-interval clearance, not a large-sample standard error.

CosMx log2 ratios use a fixed continuity correction ε (count 0.05, contact rate 0.02, neighbor fraction 0.005). The reference arm must be at least 0.05 / 0.02 / 0.01 for those three metrics. Sections need at least 30 tumor cells in each arm. Permutations: 999. Bootstrap draws for concordant-4 and for the five-patient identity check: 4999. Seed 25976224.

## CosMx

Object: figshare 25976224, 302313 author tumor cells, 164656 with CLDN4 count > 0, 8 sections, 5 patients. Neighbors are counted inside the same FOV. Coordinates are CosMx pixels at 0.18 µm/pixel. Cytotoxic cells are author NK plus CD8 T cells. T/NK adds CD4 and Treg. Immune adds B, myeloid, mast, and plasmablast. CLDN4 cuts, inside each section: detected versus undetected, top half versus bottom half, top quartile versus bottom quartile. Ties break on the cell index.

Patient value = unweighted mean of that patient's section means, then log2((high + ε) / (low + ε)). The reported median is the median of the five patients. Within-FOV permutation swaps CLDN4 arms inside each section's FOV and leaves coordinates and neighbor counts where they are.

Winning spec: **immune fraction**, **10 µm**, **detected** CLDN4 cut.

Median corrected log2 ratio **-1.307** (95% patient interval **[-2.041, -0.477]**, width 1.564). Equal-patient-weight ratio of means **0.413**.

Selection-adjusted permutation p = **0.0010** (999 within-FOV shuffles; the null statistic is the best |median| in the same 270-spec grid). Nominal permutation p for this spec alone = **0.0010**. The nominal p does not pay for the search.

| Patient | High-arm mean | Low-arm mean | High − low | Corrected log2 |
|---|---:|---:|---:|---:|
| Lung12 | 0.108 | 0.276 | -0.167 | -1.307 |
| Lung13 | 0.097 | 0.168 | -0.072 | -0.770 |
| Lung5 | 0.025 | 0.070 | -0.045 | -1.342 |
| Lung6 | 0.012 | 0.018 | -0.007 | -0.477 |
| Lung9 | 0.022 | 0.105 | -0.083 | -2.041 |

Patient means are the unweighted average of section means. For a fraction, the mean is taken over tumor cells that have at least one neighbor inside the radius.

Section means for the same spec (8/8 with high < low):

| Section | Patient | High | Low | High/low |
|---|---|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 0.023 | 0.070 | 0.325 |
| LUAD-5 R2 | Lung5 | 0.026 | 0.075 | 0.354 |
| LUAD-5 R3 | Lung5 | 0.024 | 0.065 | 0.377 |
| LUSC-6 | Lung6 | 0.012 | 0.018 | 0.642 |
| LUAD-9 R1 | Lung9 | 0.023 | 0.110 | 0.207 |
| LUAD-9 R2 | Lung9 | 0.021 | 0.100 | 0.207 |
| LUAD-12 | Lung12 | 0.108 | 0.276 | 0.393 |
| LUAD-13 | Lung13 | 0.097 | 0.168 | 0.574 |

### Locked-family reference (cytotoxic neighbor counts)

These rows are the 50 µm and 100 µm cytotoxic counts under every CLDN4 cut. The cell-weighted ratio pools tumor cells and is not a patient-level estimate. It is printed so it can be compared with the locked 0.36 / 0.52 ratios without substituting for them.

| Radius | Cut | Cell-weighted ratio | Median section ratio | Equal-patient ratio | Patient median log2 | Patient interval |
|---:|---|---:|---:|---:|---:|---|
| 50 | detected | 0.935 | 0.760 | 0.782 | -0.273 | [-0.968, 0.099] |
| 50 | median | 0.705 | 0.816 | 0.840 | -0.375 | [-0.795, -0.019] |
| 50 | quartile | 0.669 | 0.867 | 0.788 | -0.059 | [-1.076, 0.103] |
| 100 | detected | 1.121 | 0.947 | 0.911 | -0.026 | [-0.608, 0.128] |
| 100 | median | 0.836 | 0.994 | 0.954 | -0.120 | [-0.500, 0.170] |
| 100 | quartile | 0.861 | 1.078 | 0.963 | 0.129 | [-0.652, 0.277] |

Under this author-label, within-FOV, equal-patient definition, the 50 µm cytotoxic median split has equal-patient ratio 0.840 and a patient interval that just clears zero. The 100 µm intervals include zero. Those ratios are not the locked 0.36 / 0.52 summary.

### Next exclusion specs by |median|

| Neighbor | Metric | µm | Cut | Median log2 | Interval | Width | Equal-patient ratio |
|---|---|---:|---|---:|---|---:|---:|
| immune | fraction | 10 | detected | -1.307 | [-2.041, -0.477] | 1.564 | 0.413 |
| immune | fraction | 20 | detected | -0.830 | [-1.792, -0.308] | 1.484 | 0.507 |
| immune | count | 20 | detected | -0.667 | [-1.770, -0.251] | 1.519 | 0.516 |
| immune | fraction | 10 | median | -0.644 | [-1.961, -0.173] | 1.788 | 0.563 |
| tnk | contact | 20 | detected | -0.556 | [-1.227, -0.094] | 1.134 | 0.617 |
| immune | fraction | 30 | detected | -0.530 | [-1.534, -0.193] | 1.341 | 0.605 |
| tnk | count | 30 | detected | -0.508 | [-1.198, -0.038] | 1.160 | 0.698 |
| immune | fraction | 20 | median | -0.450 | [-1.687, -0.060] | 1.628 | 0.629 |

## Concordant-4

Patient table from the four-cohort integration (65 units: GSE123902 donors, GSE131907 samples, GSE205335 patients, GSE189357 patients). Outcome is `frac_tnk = n_tnk / n_cells`. CLDN4 scores are malignant % positive and malignant mean. Cuts are formed inside each cohort (median halves, tertile extremes, quartile extremes). Estimands: median of within-cohort pairwise high−low differences; pooled difference of medians; median of the four cohort differences. The bootstrap resamples patients inside each cohort and rebuilds the cut.

Locked estimand, recomputed from the same table: DerSimonian–Laird Spearman **ρ = -0.531** (two-sided p = 1.646e-05, I² = 0.0%). Analytical 95% interval on the Fisher-z scale: **[-0.697, -0.312]**. Stratified patient-bootstrap interval: [-0.975, -0.299]. The bootstrap lower tail is wider because a resample of the n=9 cohort can sit on the correlation boundary. Within-cohort permutation p (one-sided, more negative) = **0.0010**.

Winning median-Δ spec: **mal_CLDN4_pct**, **quartile** cut, **median_of_cohort_diffs**.

Median Δ (CLDN4-high − CLDN4-low T/NK fraction) **-0.339** (95% stratified patient-bootstrap interval **[-0.444, -0.093]**, width 0.352).

Selection-adjusted permutation p = **0.0030** (within-cohort CLDN4 shuffles; null statistic = best |median Δ| among the 18 specs). Nominal permutation p for this spec alone = **0.0010**.

| Cohort | n high | n low | Median T/NK high | Median T/NK low | Δ |
|---|---:|---:|---:|---:|---:|
| GSE123902 | 3 | 3 | 0.243 | 0.625 | -0.382 |
| GSE131907 | 5 | 5 | 0.048 | 0.422 | -0.374 |
| GSE189357 | 2 | 2 | 0.312 | 0.616 | -0.303 |
| GSE205335 | 5 | 5 | 0.123 | 0.396 | -0.273 |

The winning Δ is the median of those four cohort differences. GSE189357 has nine patients, so the quartile arms are two and two.

| Score | Cut | Estimand | Median Δ | 95% interval | Width | Eligible |
|---|---|---|---:|---|---:|---|
| mal_CLDN4_pct | median | hl_within_cohort | -0.219 | [-0.281, -0.085] | 0.196 | True |
| mal_CLDN4_pct | median | pooled_diff_of_medians | -0.252 | [-0.359, -0.042] | 0.317 | True |
| mal_CLDN4_pct | median | median_of_cohort_diffs | -0.251 | [-0.330, -0.051] | 0.279 | True |
| mal_CLDN4_pct | tertile | hl_within_cohort | -0.252 | [-0.353, -0.126] | 0.227 | True |
| mal_CLDN4_pct | tertile | pooled_diff_of_medians | -0.288 | [-0.479, -0.057] | 0.422 | True |
| mal_CLDN4_pct | tertile | median_of_cohort_diffs | -0.281 | [-0.407, -0.092] | 0.315 | True |
| mal_CLDN4_pct | quartile | hl_within_cohort | -0.279 | [-0.385, -0.096] | 0.289 | True |
| mal_CLDN4_pct | quartile | pooled_diff_of_medians | -0.308 | [-0.505, -0.057] | 0.448 | True |
| mal_CLDN4_pct | quartile | median_of_cohort_diffs | -0.339 | [-0.444, -0.093] | 0.352 | True |
| mal_CLDN4_mean | median | hl_within_cohort | -0.148 | [-0.255, -0.020] | 0.235 | True |
| mal_CLDN4_mean | median | pooled_diff_of_medians | -0.149 | [-0.297, -0.002] | 0.296 | True |
| mal_CLDN4_mean | median | median_of_cohort_diffs | -0.182 | [-0.296, 0.005] | 0.301 | True |
| mal_CLDN4_mean | tertile | hl_within_cohort | -0.219 | [-0.282, -0.043] | 0.239 | True |
| mal_CLDN4_mean | tertile | pooled_diff_of_medians | -0.252 | [-0.409, -0.012] | 0.397 | True |
| mal_CLDN4_mean | tertile | median_of_cohort_diffs | -0.244 | [-0.356, -0.020] | 0.335 | True |
| mal_CLDN4_mean | quartile | hl_within_cohort | -0.162 | [-0.317, -0.027] | 0.291 | True |
| mal_CLDN4_mean | quartile | pooled_diff_of_medians | -0.237 | [-0.465, 0.042] | 0.507 | True |
| mal_CLDN4_mean | quartile | median_of_cohort_diffs | -0.236 | [-0.416, -0.004] | 0.412 | True |

## What this does not say

Sections are not five extra patients. Cells are not the sample size. The search picked the contrast, so the interval on the winner is the interval for that selected contrast; the permutation p-value is the one that includes the search. A Visium spot correlation is not in this file. Public mouse Harmony objects are not merged with the private KL cohort. Nearby effector cells are not claimed to be muzzled.

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/run_patient_separation.py
```


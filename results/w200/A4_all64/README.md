# A4 — TISMO all 64 ICB slices: Tacstd2 after ICB (honest recompute)

User claim: mouse **Tacstd2 is up in 49/64** TISMO ICI-treated models,
**p = 5.8×10⁻⁵** (described as a sign test).

This folder recomputes that claim from the public TISMO in-vivo ICB table
for `Tacstd2` (all listed ICB treatments × all listed tumors). TISMO returns
exactly **64** paired slices. Cancer types are TISMO's official cell-line
labels.

## Verdict

- **49/64 is real, but it is a mean-based count, not a sign test on medians.**
  `mean(ICB) > mean(baseline)` in 49 of 64 slices (15 down, 0 ties).
- **p = 5.8e-5 is Wilcoxon signed-rank on those 64 mean differences**
  (p = 5.84×10⁻⁵), not a sign / binomial test. The actual two-sided sign
  test of 49/64 is p = 2.44×10⁻⁵.
- **Honest sign test (what a sign test is): median(ICB) vs median(baseline).**
  **40 up / 18 down / 6 ties**, n_tested = 58, two-sided exact binomial
  **p = 5.35×10⁻³**. Still an up-bias, much weaker than the user p-value.
- **Cancer-type split of the median sign test: no group is significant.**
  Mammary, the largest set (29/64 slices, 7 lines), is **16 up / 12 down /
  1 tie, p = 0.57**.
- The 64 “models” are TISMO display slices, not 64 independent cell lines
  (17 lines, 22 studies; GSE124821 alone is 19 slices). Unique cell line
  (pooled unique SRX): **12 up / 4 down / 1 tie, p = 0.077**.

## Pairing

One pair per TISMO `cell_line` label (the 64-slice universe):

- baseline = rows with `Baseline==1` (isotype / vehicle / no_treatment /
  matched co-treatment control)
- ICB = rows with `Baseline==0` (responders + non-responders pooled)

Direction = ICB location − baseline location. Ties (delta == 0) are dropped
from the sign test, as required.

## Primary numbers (64 slices)

| Test | Up | Down | Tie | N tested | p (two-sided) |
|---|---:|---:|---:|---:|---|
| Sign test, **median** (primary) | 40 | 18 | 6 | 58 | **5.35e-3** |
| Sign test, **mean** (user count) | 49 | 15 | 0 | 64 | 2.44e-5 |
| Wilcoxon signed-rank, median deltas | 40 | 18 | 6 | 58 | 7.14e-3 |
| Wilcoxon signed-rank, mean deltas (user p) | 49 | 15 | 0 | 64 | **5.84e-5** |

## Cancer-type split (median sign test)

Cancer type is TISMO `cellLineMeta.cancerType`. Mammary carcinoma /
adenocarcinoma / NOS are one group.

| Group | Slices | Lines | Up | Down | Tie | p (two-sided) |
|---|---:|---:|---:|---:|---:|---|
| Mammary | 29 | 7 | 16 | 12 | 1 | 0.57 |
| Melanoma | 14 | 4 | 8 | 2 | 4 | 0.11 |
| Colorectal | 10 | 2 | 7 | 2 | 1 | 0.18 |
| Gastric | 5 | 1 | 4 | 1 | 0 | 0.38 |
| Lung | 2 | 1 | 2 | 0 | 0 | 0.50 |
| Sarcoma | 2 | 1 | 2 | 0 | 0 | 0.50 |
| Liver | 2 | 1 | 1 | 1 | 0 | 1.00 |

Mean-based split (the method behind 49/64) looks stronger only in
**melanoma (13/1, p = 0.0018)**. That is not a median effect: several B16
slices have median Tacstd2 = 0 in both arms, and a single treated mouse
with tiny expression flips the mean up. Those four B16 slices are the
median ties.

Lung is only LLC (2 slices, both up on the median). This is not a lung-ICI
result.

## Why 64 is not 64 independent models

- 17 unique cell lines, 22 GSE/ERP studies.
- **GSE124821 = 19/64** (KPB25L / T11 / p53 lines × UV/Apobec × day3/day7/end).
- 14 baseline SRX IDs are reused across two ICB arms (4T1 GSE130472 old/young
  isotype controls shared by antiCTLA4 and antiPDL1).
- Several “baselines” are not untreated: BRAFi, radiation, birinapant,
  exercise, high-fat diet, swainsonine, anti-TGFβ, anti-GARP, anti-IL17.

Sensitivity (median sign test):

| Collapse | Up | Down | Tie | p |
|---|---:|---:|---:|---|
| Last timepoint only (49 slices) | 31 | 13 | 5 | 9.6e-3 |
| Drop GSE124821 (45 slices) | 29 | 11 | 5 | 6.4e-3 |
| Unique (line × study × subtype), unique SRX (46 units) | 30 | 11 | 5 | 4.3e-3 |
| Unique GSE (median of slice deltas) | 16 | 5 | 1 | 0.027 |
| Unique cell line, unique SRX | 12 | 4 | 1 | **0.077** |

The pooled up-bias survives some collapsing but is **not significant at
α = 0.05 once the unit is the cell line**.

## What this does and does not support

Supports: in TISMO's 64 ICB-vs-control slices, Tacstd2 more often has a
higher median after ICB than before (40 vs 18, p = 0.005). The user 49/64
count is the same table with means.

Does not support: a 5.8×10⁻⁵ sign-test claim; a general all-cancer ICB
induction of Tacstd2 that holds inside mammary or lung; independence of
the 64 tests; a responder-vs-nonresponder effect (almost no slices have
both R and NR labeled).

## Outputs

- `per_model.tsv` — 64 slices, medians, means, directions, cancer type
- `cancer_type_split.tsv` — median and mean sign tests by cancer group
- `sensitivity.tsv` — sign / Wilcoxon / collapse variants
- `independent_units.tsv` — unique-SRX pairs after dropping ICB agent + timepoint
- `tacstd2_vivo_icb.csv` — TISMO per-sample table
- `summary.json` / `audit.json`
- `figures/fig_waterfall_median_delta.png`
- `figures/fig_cancer_type_sign.png`

Pipeline: `scripts/w200/A4_all64/`.

# Figure 1l panel map

Three separate concordant-4 relations. Data visualization only. Every drawn estimate is recomputed from `source_data/units.tsv` and checked against `source_data/effects.tsv` by `plot_fig1l.py`. The sample size is tumour units (donor, sample, or patient), never cells.

Regenerate with:

```bash
python3 compute_scores.py
python3 plot_fig1l.py
```

Dependencies: numpy, pandas, scipy, matplotlib. Text is set in Arimo. `pubstyle.py` holds the palette and formatters.

Part 2 is not in this figure. Malignant CLDN4 percent-positive versus T/NK (Fig. 5 / PR #780, n = 65, pooled ρ = −0.531) is a different estimand. The Fisher-*z* pooler reproduces that published ρ when given that figure's cohort coefficients; the check is in `compute_scores.py` and is not plotted.

## What is drawn

| File | Panels | Content |
| --- | --- | --- |
| `Fig1l.pdf` / `.png` / `.svg` | a–c | Primary forests. Hallmark apical junction. |
| `Fig1l_ED_scatters.pdf` / `.png` / `.svg` | a–c | One point per included unit. |
| `Fig1l_ED_sensitivity.pdf` / `.png` / `.svg` | a–d | Hallmark without CLDN4, and KEGG tight junction. |

Shared forest geometry: cohort order top to bottom GSE123902, GSE131907, GSE205335, GSE189357, then the pooled diamond. x-axis is Spearman ρ from −1 to 1, with a line at 0. Cohort colour: GSE123902 `#0F4D92`, GSE131907 `#42949E`, GSE205335 `#9A4D8E`, GSE189357 `#8BCF8B`. The diamond is `#272727`. The number to the right of each row is ρ at two decimals. The pooled row also prints I². Marker size is constant. Interval half-width is 1.96 / sqrt(n − 3) on the Fisher-*z* scale.

## Panel a — TACSTD2 vs junction score

- x from `effects.tsv` rows with `relation == tacstd2_vs_junction` and `score == junction_hallmark`.
- y is cohort, then the four-cohort pool.
- n = 13, 21, 21, 9 and pooled 64.
- ρ = 0.483516, 0.266234, 0.614286, 0.116667; pooled 0.427782.
- P = 0.094135, 0.243399, 0.003050, 0.765008; pooled 0.0009781.
- I² = 0.

## Panel b — TACSTD2 vs T/NK fraction

- `relation == tacstd2_vs_tnk`, `score == junction_hallmark` (the score column records which junction definition this relation block belongs to; the y variable is `frac_tnk`).
- ρ = 0.032967, 0.257143, 0.205195, −0.116667; pooled 0.154674.
- P = 0.914856, 0.260473, 0.372237, 0.765008; pooled 0.260846.
- I² = 0.
- Figure label 0.15 is two-decimal rounding of 0.154674, not the panel c pool.

## Panel c — junction score vs T/NK fraction

- `relation == junction_vs_tnk`, `score == junction_hallmark`.
- ρ = −0.065934, −0.166234, 0.437662, 0.416667; pooled 0.146115.
- P = 0.830542, 0.471422, 0.047234, 0.264586; pooled 0.407291.
- I² = 35.272685.
- Figure label 0.15 is two-decimal rounding of 0.146115.

## Extended Data scatters

Rows are the three primary relations. Columns are the four cohorts. Points are `source_data/units.tsv` rows with `included == True`. Axis limits are the included-unit min and max of that row, padded by 6% (x) and 8% (y), and are shared across the four cohorts in the row. Printed ρ, P and n match the cohort rows of panels a–c. No guide line is drawn.

## Extended Data sensitivity

| Panel | `score` | `relation` |
| --- | --- | --- |
| a | `junction_hallmark_no_cldn4` | `tacstd2_vs_junction` |
| b | `junction_hallmark_no_cldn4` | `junction_vs_tnk` |
| c | `junction_kegg` | `tacstd2_vs_junction` |
| d | `junction_kegg` | `junction_vs_tnk` |

TACSTD2 versus T/NK is not repeated: it does not use the junction gene list.

## Source files

| File | Role |
| --- | --- |
| `source_data/units.tsv` | One row per locked unit. Scores are blank when `included` is false. |
| `source_data/effects.tsv` | Cohort and pooled Spearman results for the primary score, the CLDN4-dropped score, and KEGG. |
| `source_data/genes_used.tsv` | Hallmark and KEGG membership and which genes entered each mean. |
| `source_data/calibration_tacstd2.tsv` | Max absolute difference versus the prior TACSTD2 log2(TMM-CPM+1) column. |
| `input/` | Malignant UMI sums, locked unit tables, frozen gene sets, and the prior TACSTD2 column used for calibration. |
| `KEY_STATS.json` | Score definition, inclusion, pooler, and the primary estimates. |

`units.tsv` columns used for the drawings: `patient`, `cohort`, `unit`, `frac_tnk`, `tacstd2`, `junction_hallmark`, `junction_hallmark_no_cldn4`, `junction_kegg`, `included`. `n_cells`, `n_malignant` and `n_tnk` are stored and not used as the sample size.

## Inclusion

Locked units, then keep those with a malignant pseudobulk column.

| Cohort | Unit | Malignant definition | Included |
| --- | --- | --- | ---: |
| GSE123902 | donor | marker, eligible primary or metastasis | 13 |
| GSE131907 | sample | author, tumour origin, ≥20 malignant cells | 21 |
| GSE205335 | patient | author | 21 |
| GSE189357 | patient | marker, eligible | 9 |

Excluded: P4001 (GSE205335), no malignant count column. That unit remains in the Part 2 n = 65 CLDN4 percent-positive table. Cohorts outside this lock (GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526) were not added.

## Score

Primary: `input/a8_sets.json` key `HALLMARK_APICAL_JUNCTION` (200 genes). TACSTD2 removed (not present). Mean of genes still in the matrix after count ≥ 10 in ≥ 3 units (194). Sensitivity drops CLDN4 as well (193). KEGG uses `KEGG_TIGHT_JUNCTION` with the same filter and the same mean (154 genes; 153 without CLDN4). Expression is malignant UMI-sum, TMM, log2(CPM+1), the same transform as the calibrated TACSTD2 values.

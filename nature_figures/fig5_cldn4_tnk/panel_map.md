# Figure 5 panel map

Concordant-4 malignant CLDN4 percent-positive versus T/NK fraction. Data visualization only. Every annotated estimate is read from the TSVs below. Sample size is tumour units (donor, sample, or patient, as defined per cohort), never cells or genes.

Regenerate with `python3 plot_fig5.py` from this directory (matplotlib, numpy). Text is set in Arimo, the Arial-metric face available in this environment; Arial itself is not installed. `pubstyle.py` supplies the palette and forest helpers and is otherwise unchanged.

## Source files

| File | Used by | Role |
| --- | --- | --- |
| `source_data/tnk_units.tsv` | a, b | One row per unit: `cldn4_pct`, `frac_tnk`, `cohort`, `quartile`, `unit`, `malig_def`, `in_count_matrix` |
| `source_data/tnk_pooled.tsv` | a, b, c | `score == pct`, `kind == concordant4` row only |
| `source_data/tnk_singles.tsv` | c | `score == pct`, `kind == single` rows |
| `source_data/family_de.tsv` | d | `split == q4q1` and `cohort == GSE123902+GSE131907+GSE205335+GSE189357` |

Columns present in those files but not drawn are listed at the end.

## Panel a — scatter

- x: `cldn4_pct` (percent scale, as stored). y: `frac_tnk`.
- All 65 rows. Cohort colours: GSE123902 `#0F4D92`, GSE131907 `#42949E`, GSE205335 `#9A4D8E`, GSE189357 `#8BCF8B`.
- Annotated ρ, P, n and I² come from the pooled `pct` row: ρ = −0.5311678045689989 displayed as −0.53; P = 1.6462232944573955×10⁻⁵ displayed as 1.6 × 10⁻⁵; n = 65; I² = 0.
- The pool is a DerSimonian–Laird combination of cohort Fisher-*z* Spearman coefficients (`tnk_pooled.tsv` note), not a Spearman fit to the superimposed cloud.
- Check, not annotated: Spearman on the 65 `tnk_units.tsv` rows is −0.549694055944056. Cohort rows in `tnk_singles.tsv` match Spearman values recomputed from the units.
- Grey line: ordinary least-squares fit of `frac_tnk` on `cldn4_pct` for those 65 points, drawn from the minimum to the maximum observed x. Slope and intercept are not annotated.

## Panel b — Q4 versus Q1

- Groups: `quartile == Q1` (n = 19) and `quartile == Q4` (n = 16), the same counts as `n_q1` and `n_q4` on the pooled `pct` row.
- Quartile rule, from that row: `within_cohort_rank_then_qcut`. Q1 is low CLDN4, Q4 is high CLDN4, within each cohort.
- y values are `frac_tnk`. Medians 0.4334841628959276 and 0.12380844872566725 are labeled 0.43 and 0.12.
- Δ median = −0.30967571417026035, labeled −0.31. Rank-biserial r = −0.7236842105263157 (`r_rb`), labeled −0.72. P = 2.8794861470847027×10⁻⁴ (`p_q4q1`), labeled 2.9 × 10⁻⁴.
- Boxes: median and interquartile range. Whiskers: 1.5 × IQR. Boxplot fliers are off because the strip draws every unit. Strip jitter is uniform on [−0.13, 0.13] with numpy Generator seed 21.
- Within-cohort Q1/Q4 counts: GSE123902 4/3, GSE131907 6/5, GSE205335 6/6, GSE189357 3/2. `tnk_singles.tsv` leaves Q4-versus-Q1 statistics blank for GSE189357 (`thin_q4`). No per-cohort Q4-versus-Q1 P is drawn.

## Panel c — forest

- Cohort ρ, P and n from `tnk_singles.tsv` (`score == pct`), displayed at two decimal places for ρ and with the shared P formatter.
- Pooled diamond: x = `rho`, horizontal extent = `ci95_lo` to `ci95_hi` (−0.6967711225150761 to −0.3118052568209761), labeled −0.70 to −0.31. This is the only 95% CI on the figure.
- `tnk_singles.tsv` has no confidence-interval columns. Cohort intervals are not calculated and not drawn.
- Point diameter is `1.48 × sqrt(n)`, so marker area scales with unit n. This is a display scale, not the DerSimonian–Laird weight.
- I² = 0 is taken from the pooled row.

## Panel d — malignant gene sets

- Four combined `q4q1` rows, in this order (not sorted by effect): IFN, chemokine, MHC-I/APM, TJ.
- x = `logFC`. Error bars = ± `se` from the same row (standard error, not a confidence interval).
- P labels are the `p` column. FDR is in the caption only (`fdr`).
- n = 34, n_q1 = 18, n_q4 = 16. The unit table has 19 Q1 units; `in_count_matrix` is False for one GSE205335 Q1 unit, which accounts for 18 versus 19. That unit remains in panels a and b.
- `n_genes` is the gene-set size (IFN 221, chemokine 25, MHC-I/APM 21, TJ 194). It is not a sample size.
- Reported `t` and `df` (df = 29) are not drawn.

## Present in the source files and not drawn

- `tnk_pooled.tsv` / `tnk_singles.tsv` rows with `score == mean` (pooled mean-score ρ = −0.403, P = 0.00186, I² = 0, with its own `ci95_lo`/`ci95_hi`). The figure uses percent-positive only.
- Pooled `stouffer_z` and `stouffer_p`.
- Per-cohort rows of `family_de.tsv`, and all `split == continuous` rows.
- `bulk_immune.tsv`, `family_barrier_inhibitory.tsv`, and `de_q4q1_combined_families.tsv` were in the source pack and are not used by this figure.

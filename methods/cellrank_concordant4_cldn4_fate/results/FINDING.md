# CellRank 2 fate toward barrier CLDN4-high (concordant-4)

ADDITIVE. **CLDN4-only.** This does not replace the locked concordant-4 patient result
(malignant CLDN4 %pos vs T/NK, n=65, ρ=−0.531, P=1.65×10⁻⁵, I²=0).
Cell counts below are not that n.

Question: in tumor epithelium from the same four public scRNA cohorts,
what is the CellRank 2 absorption probability of cells that are not already
in a terminal state, toward a pre-specified **barrier CLDN4-high** sink,
when the competing sink is cycle-high and CLDN4-low?

## Population

- GSE123902 (Laughney): tumor/metastasis donors; marker epithelium `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. NORMAL files dropped.
- GSE131907 (Kim): author `Cell_subtype==Malignant cells` in tumor-bearing samples with ≥20 malignant cells. nLung AT2/AT1/ciliated cells are not in the graph.
- GSE205335 (Ahn/Lee): author `lineage.sub==Malignant cells`, normal tissue dropped, patients with ≥20 malignant cells.
- GSE189357 (Zhu): TD1–TD9, same marker epithelium gate. No author cell-type column.
- Not included: GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, public mouse, private 8KL.
- Cap: ≤400 epithelial cells per unit after QC (n_genes≥200, n_UMI≥500, mito%<20), seed=1.
- Each dataset is its own graph. Harmony corrects the unit key inside the dataset. Datasets are not concatenated.

## State definitions (pre-specified)

- Barrier score **holds CLDN4 out**: CLDN3, CLDN7, CDH1, OCLN, TJP1, EPCAM, TACSTD2, KRT8, KRT19.
- `barrier_cldn4_high`: CLDN4 ≥ dataset 80th percentile (or CLDN4>0 if that percentile is 0) AND barrier score ≥ median AND airway score (FOXJ1, CAPS, PIFO, TPPP3, KRT5) < 75th percentile. Airway filter relaxes only if the terminal would otherwise have <25 cells.
- `cycle_cldn4_low`: cycle score (MKI67, TOP2A, STMN1, PCNA, HMGB2) ≥ 80th percentile AND CLDN4 ≤ median. Quantile relaxes to 70th once if n<25.
- `at2_like` root: AT2 score (SFTPC, SFTPB, SFTPA1, SFTPA2, NAPSA, SLC34A2) ≥ 75th percentile AND CLDN4 ≤ median, only if ≥30 cells. This root orients diffusion pseudotime. It is not a normal-lung atlas root.
- Primary kernel if that root exists **and** mean DPT of the barrier terminal exceeds mean DPT of AT2-like cells, with ≥95% finite DPT: **0.8 × PseudotimeKernel + 0.2 × ConnectivityKernel** (CellRank 2 default mixture). Otherwise primary kernel = **ConnectivityKernel** alone. Velocity is off (these GEO matrices have no spliced/unspliced counts).
- Fate of a unit = mean `P(barrier_cldn4_high)` over its non-terminal cells. Units with <15 non-terminal cells are listed and not tested. Terminal cells are left out of the mean.

## Honest n

| dataset | unit | cells QC / epithelial QC | cells used | units used |
|---|---|---:|---:|---:|
| GSE123902 | donor | 26056 / 4208 | 3090 | 13 |
| GSE131907 | sample | 24783 / 24783 | 6533 | 21 |
| GSE205335 | patient | 28512 / 28512 | 7552 | 22 |
| GSE189357 | patient | 119264 / 13789 | 3600 | 9 |

Units entering the fate test (non-terminal cells ≥15): **n=65** (GSE123902 13, GSE131907 21, GSE205335 22, GSE189357 9).
That unit split matches the locked concordant-4 table. The estimand here is fate, not the T/NK correlation. Do not quote cell counts as n.

## Primary result

DerSimonian–Laird meta-analysis of cohort mean unit fate minus 0.5 (k=4: GSE123902, GSE131907, GSE205335, GSE189357): **Δ = 0.045** (SE 0.048, descriptive z=0.94, p=0.3461, I²=79.3%).
Positive Δ means the unit-average fate sits above the two-sink coin-flip of 0.5.

| dataset | kernel | units | median fate | mean fate | fraction >0.5 | Wilcoxon p vs 0.5 | contact absorption (intermediate) | within-unit CellRank mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| GSE123902 | connectivity | 13 | 0.399 | 0.477 | 0.46 | 0.7869 | 0.492 | 0.427 |
| GSE131907 | connectivity | 21 | 0.690 | 0.643 | 0.86 | 0.0014 | 0.556 | 0.566 |
| GSE205335 | connectivity | 22 | 0.631 | 0.587 | 0.82 | 0.0190 | 0.523 | 0.540 |
| GSE189357 | connectivity | 9 | 0.452 | 0.451 | 0.22 | 0.3008 | 0.448 | 0.431 |

Wilcoxon p-values are descriptive (units inside one cohort). They are not a substitute for the meta-analytic Δ.

## Reading

The pre-specified 0.8/0.2 pseudotime mixture is the primary kernel in 0 of 4 estimable cohorts. In the other 4, AT2-like cells are present but the barrier terminal misses the downstream gate (barrier mean diffusion pseudotime must exceed the AT2-like mean by more than 0.02), so the primary kernel is ConnectivityKernel. Mean DPT of the cycle-high CLDN4-low sink exceeds mean DPT of the barrier sink in 4 of 4 cohorts.

The four cohort means disagree (I²=79.3%). Above 0.5: GSE131907 mean 0.643 (0.86 of units >0.5; within-unit mean 0.566); GSE205335 mean 0.587 (0.82 of units >0.5; within-unit mean 0.540). At or below 0.5: GSE123902 mean 0.477 (0.46 of units >0.5; within-unit mean 0.427); GSE189357 mean 0.451 (0.22 of units >0.5; within-unit mean 0.431). The meta-analytic Δ is compatible with an even split between the two sinks. A single concordant fate toward barrier CLDN4-high is not supported.

Within-unit ConnectivityKernel fates, recomputed inside each unit with no cross-unit edges, keep that split. Contact-matrix absorption of the intermediate state has the same sign relative to 0.5 and sits closer to 0.5.

GPCCA macrostates on the connectivity kernel were named by CLDN4 only after fitting. GSE123902 CLDN4 means 0.091–0.896 (gap 0.659), non-member fate toward the highest 0.633. GSE131907 CLDN4 means 1.293–1.690 (gap 0.151), non-member fate toward the highest 0.016. GSE205335 CLDN4 means 1.103–1.795 (gap 0.376), non-member fate toward the highest 0.178. GSE189357 CLDN4 means 0.362–1.590 (gap 0.565), non-member fate toward the highest 0.147. A small gap (GSE131907, GSE205335) means the macrostates are all CLDN4-high, so the highest one is not a distinct sink. A gap of at least 0.5 (GSE123902, GSE189357) separates the highest-CLDN4 macrostate from the rest.

The multinomial logistic puts a positive coefficient on the barrier score (CLDN4 held out) for the barrier-neighbor class in every cohort (`results/tables/multinomial_coefs.tsv`). That is co-localization on the kNN graph. CLDN4 was not a predictor, and the fit is cell-level, so it is not a unit-level p-value.

## What the terminals actually are

| dataset | n barrier / cycle / AT2 / intermediate | mean CLDN4 barrier vs cycle | mean barrier-score (CLDN4 held out) barrier vs rest | fraction of CLDN4-high barrier candidates removed as airway | DPT barrier vs AT2 | orientation ok |
|---|---|---|---|---:|---|---|
| GSE123902 | 344 / 397 / 169 / 2180 | 2.190 vs 0.011 | 0.579 vs 0.045 | 0.37 | 0.032 vs 0.046 | False |
| GSE131907 | 631 / 816 / 788 / 4298 | 2.660 vs 0.386 | 0.722 vs 0.278 | 0.36 | 0.024 vs 0.017 | False |
| GSE205335 | 779 / 792 / 630 / 5351 | 2.893 vs 0.270 | 1.073 vs 0.411 | 0.35 | 0.099 vs 0.114 | False |
| GSE189357 | 410 / 488 / 164 / 2538 | 2.484 vs 0.000 | 0.613 vs 0.024 | 0.34 | 0.089 vs 0.070 | False |

## Secondary descriptions (not the primary n)

- **Connectivity-only CellRank** is always computed. When orientation fails it is the primary kernel; otherwise it is a sensitivity that does not use pseudotime. Unit means are in `results/tables/unit_fate.tsv` (`mean_fate_connectivity`).
- **Within-unit CellRank** (ConnectivityKernel on that unit alone, states re-defined inside the unit, ≥80 epithelial cells) does not let Harmony stitch patients. Cohort means of those unit fates are in the table above. Missing units are too small or lack both sinks.
- **Contact matrix:** row-normalized kNN weights aggregated by the four states, then absorption of the intermediate (and AT2-like) state into the barrier sink. This is an undirected snapshot Markov model. It is the contact-tracing estimator. Table: `results/tables/contact_transition.tsv`.
- **Multinomial logistic:** one random graph neighbor per cell; class = neighbor state; predictors = AT2 score, barrier score with CLDN4 held out, and cycle score. CLDN4 itself is not a predictor. Coefficients are cell-level and dependent, so they are not a p-value claim. Table: `results/tables/multinomial_coefs.tsv`.
- **GPCCA macrostates** on the connectivity kernel are fit without using CLDN4 to choose them. The macrostate with the highest membership-weighted CLDN4 is named afterwards. Its fate among non-members is a sensitivity, not a test against 0.5 (there are 3–4 macrostates). Details are in `results/summary.json`.

## What this is not

- Not evidence that CLDN4 causes T/NK exclusion. The locked exclusion result is a different analysis.
- Not a timed trajectory and not RNA velocity. Public matrices here have no spliced/unspliced layer.
- Not a statement that AT2 becomes tumor. The root is an AT2-like program inside the tumor epithelial gate, and only when those cells exist and sit upstream.
- Not a cross-dataset Harmony object, and not a merge with mouse or with the private 8KL matrices.
- A fate near 0.5 means the two pre-specified sinks split the non-terminal mass. That is a result, not a failed audit of the n=65 correlation.

## Software

Python packages: cellrank 2.3.3, scanpy 1.12.4, anndata 0.13.4, numpy 2.4.4, pandas 3.0.6, scipy 1.16.3, scikit-learn 1.9.1, harmonypy 2.0.2.

## Reproduce

```bash
python3 methods/cellrank_concordant4_cldn4_fate/download.py
python3 methods/cellrank_concordant4_cldn4_fate/analyze.py
```

Figures: `results/figures/fig_cellrank_fate.png`.

# FINDING — max-effect barrier outgoing, concordant-four

ADDITIVE. The primary Q4/Q1 CellChat table in `FINDING.md` is unchanged.
CLDN4 only. No dual-high. No TACSTD2 gate. Concordant four only.
Honest unit = patient / locked sample.

Barrier/inhibitory outgoing is the JAM, NECTIN2–TIGIT, CDH1, and LGALS9 family.
IFN/recruit is a separate arm and is not added into the barrier score.

Search space for the probability-scale gate: 6912 barrier specs (split × receiver × mean × population.size × universe × aggregation × patient filter).
A spec clears the gate when all four cohort means are positive and the one-sided sign-flip p is ≤ 0.05 (B = 10000, seed 3979).
That p is the p of the reported spec. It is not multiplied by the number of specs.

Reference spec (14-pair universe, Q4/Q1, T/NK, truncatedMean trim 0.1, population.size TRUE, edge sum, all units): mean Δ = +0.00368135. This matches the primary table (max absolute patient-level difference below 1e-6).

## Largest probability-scale mean Δ that clears the gate

**Mean Δ = +0.190688** on 19 patients.
Spec: `barrier_pathway | top10 | CD8 | trimean | pop=FALSE | pathway_sum | filter=detected_ge3`.
Cliff's delta (zeros stay in the denominator) = +0.789474. Cliff's delta among nonzero pairs = +0.789474. Cohen's d (mean / sd) = +1.03628.
Median patient Δ = +0.207247. Sign-flip one-sided p = 2.00e-04. Sign-flip two-sided p = 3.00e-04. Wilcoxon p = 5.80e-04.
Cohort n (123902 / 131907 / 205335 / 189357) = 2 / 4 / 12 / 1.
Cohort mean Δ = +0.00170472 / +0.279353 / +0.200567 / +0.0954427.
A cohort mean is the mean of the patients who remain after the filter. One patient makes that cohort mean equal to that patient's delta.
Signs: 17 positive, 2 negative, 0 zero. Median CellChat expression max in this gene universe = +3.87566.

Largest eligible Δ that keeps the primary 14-pair gene universe (edge sum or mean of detected edges): mean Δ = +0.0914292, Cliff = +0.84, Cohen's d = +1.09212, n = 25, spec `edges14 | top10 | CD8 | trimean | pop=FALSE | edge_sum | filter=detected_ge3`.

Largest eligible Δ with at least 3 patients in every cohort: mean Δ = +0.150763, Cliff = +0.913043, Cohen's d = +1.29219, n = 46 (5/14/18/9), cohort means +0.0444458 / +0.178096 / +0.190435 / +0.0879635, spec `barrier_pathway | top10 | CD8 | thresh10 | pop=FALSE | pathway_sum | filter=detected_ge3`.

Largest eligible Δ that keeps every unit passing the cell-count floors (filter = all): mean Δ = +0.136588, Cliff = +0.90566, Cohen's d = +1.18388, n = 53 (8/15/21/9), spec `barrier_pathway | top10 | CD8 | thresh10 | pop=FALSE | pathway_sum | filter=all`.

### Next eligible probability specs

| mean Δ | Cliff | Cohen d | n | spec |
|---:|---:|---:|---:|---|
| +0.190688 | +0.789474 | +1.03628 | 19 | `barrier_pathway | top10 | CD8 | trimean | pop=FALSE | pathway_sum | filter=detected_ge3` |
| +0.176819 | +0.909091 | +1.13985 | 22 | `barrier_pathway | detected | CD8 | trimean | pop=FALSE | pathway_sum | filter=detected_ge3` |
| +0.159007 | +0.84 | +1.07919 | 25 | `barrier_pathway | top10 | CD8 | trimean | pop=FALSE | edge_sum | filter=detected_ge3` |
| +0.150763 | +0.913043 | +1.29219 | 46 | `barrier_pathway | top10 | CD8 | thresh10 | pop=FALSE | pathway_sum | filter=detected_ge3` |
| +0.145192 | +0.6 | +0.931699 | 20 | `barrier_pathway | top20 | CD8 | trimean | pop=FALSE | pathway_sum | filter=detected_ge3` |
| +0.142228 | +0.875 | +1.15153 | 48 | `barrier_pathway | detected | CD8 | thresh10 | pop=FALSE | pathway_sum | filter=detected_ge3` |
| +0.139214 | +0.923077 | +1.21175 | 52 | `barrier_pathway | top10 | CD8 | thresh10 | pop=FALSE | pathway_sum | filter=drop_zero_delta` |
| +0.138434 | +0.909091 | +1.13614 | 44 | `barrier_pathway | top10 | CD8 | thresh10 | pop=FALSE | pathway_sum | filter=both_arms_positive` |

## IFN/recruit on the same patients

Same patients as the barrier winner. The seven IFN/recruit edges are scored in the 14-pair universe, so an HLA-scale max-normalization is the one used for this arm even when the barrier winner uses the pathway universe.

IFN/recruit edge sum: mean Δ = +0.145431 (observed high>low; thesis expect low>high), n = 19, Cliff = +0.473684, Cohen's d = +0.634915.
Thesis-direction sign-flip (low>high) p = 0.992. Two-sided sign-flip p = 0.0141. Wilcoxon p = 0.023.
Cohort mean Δ (123902 / 131907 / 205335 / 189357) = +0.175325 / +0.250439 / +0.145751 / -0.338222.

### Pairs inside the winning probability

| pair | n | mean Δ | cohort means 123902 / 131907 / 205335 / 189357 | sign-flip greater p |
|---|---:|---:|---|---:|
| JAM1_ITGAL_ITGB2 | 19 | +0.0503172 | +0.000979615 / +0.0146611 / +0.0702822 / +0.0520362 | 1.00e-04 |
| NECTIN2_TIGIT | 19 | +0.0273971 | +0.0138756 / +0.0505051 / +0.0238226 / +0.00490191 | 1.00e-04 |
| CDH1_ITGAE_ITGB7 | 19 | +0.0138952 | +0.0115895 / +0.0133518 / +0.0156186 / +0 | 0.0011 |
| CDH1_KLRG1 | 19 | +0.0247294 | +0 / +0.0333025 / +0.0248453 / +0.0385046 | 2.00e-04 |
| LGALS9_HAVCR2 | 19 | +0.00206535 | -0.00273891 / +0.00819912 / +0.000993578 / +0 | 0.2 |
| LGALS9_CD44 | 19 | +0.0226562 | -0.00768638 / +0.0516132 / +0.019949 / +0 | 0.0194 |
| LGALS9_CD45 | 19 | +0.0366934 | -0.0126208 / +0.0772319 / +0.0344573 / +0 | 0.0163 |

Sum of these seven pair-mean deltas = +0.177754. The winning pathway sum is +0.190688. Other interactions in the JAM, NECTIN, CDH1, and GALECTIN pathways make up the difference.

## Largest fold / log-odds that clears the same gate

These are not probability differences. Log2 fold uses units with both arm sums positive. Log-odds uses edges positive on both arms, clipped to [1e-6, 1−1e-6].

**Mean = +3.29148** (log2_fold_pathway). Cliff = +0.756098. Cohen's d = +1.20344. Median = +2.80129. n = 41. Sign-flip p = 1.00e-04.
Spec: `barrier_pathway | detected | TNK | trunc10 | pop=TRUE | log2_fold_pathway | filter=detected_ge3`.

## Residual after total sender strength

On the winning split, receiver, mean, and population.size, one overexpressed full network (thresh.p = 0.05) supplies both the barrier edges and the total sender sum. They share one max-normalization. The seven barrier pairs are kept even when the overexpression filter would drop them. Total sender strength sums the overexpressed interactions only.

Within that network, on the winning patient set (n = 19): barrier-edge mean Δ = +0.0787093 (sign-flip p 2.00e-04, Cliff +0.789474, Cohen d +1.04591). Total sender mean Δ = +0.787071 (sign-flip p 1.00e-04, Cliff +0.684211, Cohen d +1.28438).
Mean barrier Δ / mean total Δ = +0.100003.
Proportional residual (barrier_high − barrier_low × total_high / total_low): mean = +0.074435, n = 19, Cliff = +0.789474, Cohen's d = +1.05518, sign-flip p = 2.00e-04.
Composition residual (barrier/total on high minus barrier/total on low): mean = +0.0315999, n = 19, sign-flip p = 1.00e-04, Cliff = +0.789474, Cohen's d = +1.22522.
OLS of within-network barrier Δ on total Δ: intercept = +0.0215182, slope = +0.0726632, R² = +0.350111. OLS R² of the headline sweep Δ on this total Δ = +0.241745.

## What this sweep is

Prepared patients with a matrix: 65.
Sender arms need at least 10 cells. Receivers need at least 20. Each group is capped at 200.
population.size TRUE multiplies by (n_sender / N) × (n_receiver / N) after the Hill probability.
The 14-pair universe is the normalization used in the primary table. The barrier-pathway universe renormalizes inside JAM, NECTIN, CDH1, and GALECTIN, so its absolute probabilities are a different scale.
Top 25% vs bottom 25% is the same rank rule as Q4/Q1.

## Reproduce

```bash
Rscript methods/cellchat_v2_network_roles_concordant4/scripts/run_max_effect_sweep.R --raw=/tmp/concordant4_raw
```


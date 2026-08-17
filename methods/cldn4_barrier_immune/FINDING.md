# CLDN4-high vs CLDN4-low marks a barrier / immune-low state

**CLDN4-first.** The gate is **CLDN4-high vs CLDN4-low** (median split, and Q4 vs Q1). TACSTD2 is a companion column only. Dual-high (TACSTD2-high AND CLDN4-high) is **not** required and is not tested here.

Thesis (given, not re-audited): CLDN4-high marks a tight-junction barrier / immune-low state. Patient is the unit. Numbers are computed Spearman ρ and Cliff δ on **existing** public patient-level scores. p-values are descriptive.

GSE207422 A3 TACSTD2 is taken as given (n=12, ρ=−0.490, p=0.106). CLDN4 on that same 12-patient table is one honest row, not an A3 re-cut. No slide was re-scored.

Sources: `scrna_meta_tnk`, `scrna_cldn4_combo` TLS extract, GSE207422 given table, GSE253013 9-patient tumor extract, GSE131907, GSE148071 TISCH + TLS, GSE241934 IIT/Real, GSE285029, GSE218989, CPTAC LUAD / LSCC.

## Full-pool continuous Spearman (one row, not the answer)

DerSimonian–Laird RE on Fisher-z of per-cohort Spearman ρ. This is the merge-everything row.

| endpoint | k | N | pooled ρ (p, I²) |
|---|---:|---:|---|
| T/NK | 5 | 92 | −0.086 (0.45, 0%) |
| CD8A | 7 | 878 | −0.062 (0.31, 55%) |
| ImmuneScore | 6 | 879 | +0.002 (0.96, 8%) |
| CXCL13 | 6 | 213 | −0.260 (0.030, 50%) |

The answer is the **per-cohort gated table** below, not this pool. CXCL13 is the only full-pool row with ρ<0 and p<0.05; T/NK / CD8A / ImmuneScore pools are near zero because they mix null bulk ImmuneScore with CD8-low tails.

## Per-cohort table — CLDN4-high vs CLDN4-low

Cliff δ = P(high>low) − P(high<low). **Negative δ = CLDN4-high is immune-low.** |δ|≈0.33 medium, ≥0.47 large. Continuous ρ is the same for both gates (CLDN4 vs endpoint, ungated).

### T/NK (scRNA fraction)

| cohort | n | ρ (p) | median δ (p) nH/nL | Q4 vs Q1 δ (p) nH/nL |
|---|---:|---|---|---|
| GSE207422 (A3 given) | 12 | −0.091 (0.78) | +0.056 (0.94) 6/6 | −0.333 (0.70) 3/3 |
| GSE253013 extract | 9 | −0.333 (0.38) | −0.300 (0.56) 5/4 | −0.556 (0.40) 3/3 |
| GSE131907 tumor-site epi | 36 | −0.210 (0.22) | **−0.444 (0.024)** 18/18 | −0.136 (0.66) 9/9 |
| GSE241934_IIT | 11 | −0.064 (0.85) | −0.267 (0.54) 6/5 | +0.111 (1.00) 3/3 |
| GSE241934_Real eligible | 24 | +0.178 (0.40) | +0.306 (0.21) 12/12 | +0.111 (0.82) 6/6 |

GSE131907 median split is the T/NK row that moves: CLDN4-high has lower T/NK (δ=−0.44, p=0.024, n=36). Q4 vs Q1 on the same 36 is weaker (tails are small). GSE207422 A3 CLDN4 vs T/NK stays near null (honest; the given TACSTD2 slide is the A3 result). GSE241934_Real is the one T/NK cohort with ρ>0.

### CD8A

| cohort | n | ρ (p) | median δ (p) nH/nL | Q4 vs Q1 δ (p) nH/nL |
|---|---:|---|---|---|
| GSE218989 | 355 | **−0.172 (0.0011)** | −0.124 (0.043) 178/177 | **−0.307 (0.00041)** 89/89 |
| CPTAC LUAD RNA | 110 | **−0.225 (0.018)** | **−0.256 (0.021)** 55/55 | −0.270 (0.084) 28/28 |
| GSE131907 | 36 | −0.094 (0.58) | −0.352 (0.074) 18/18 | −0.012 (1.00) 9/9 |
| GSE241934_IIT | 11 | −0.109 (0.75) | −0.467 (0.25) 6/5 | −0.111 (1.00) 3/3 |
| GSE241934_Real | 24 | +0.185 (0.39) | +0.264 (0.29) 12/12 | +0.056 (0.94) 6/6 |
| GSE285029 | 234 | +0.067 (0.30) | −0.011 (0.89) 117/117 | +0.176 (0.10) 59/59 |
| CPTAC LSCC RNA | 108 | +0.035 (0.72) | −0.110 (0.32) 54/54 | +0.021 (0.90) 27/27 |

GSE218989 Q4 vs Q1 is the large-n CD8 cut: n=89 vs 89, δ=−0.31, p=4.1×10⁻⁴. CPTAC LUAD RNA median split agrees (δ=−0.26, p=0.021). GSE285029 and CPTAC LSCC RNA do not.

### ImmuneScore

| cohort | n | score | ρ (p) | median δ (p) | Q4 vs Q1 δ (p) |
|---|---:|---|---|---|---|
| GSE218989 | 355 | ESTIMATE mean-z | +0.044 (0.41) | +0.091 (0.14) | +0.042 (0.63) |
| GSE285029 | 234 | immune8 | +0.080 (0.22) | −0.001 (0.99) | +0.198 (0.064) |
| CPTAC LUAD RNA | 110 | ESTIMATE | −0.079 (0.41) | −0.074 (0.51) | −0.102 (0.52) |
| CPTAC LSCC RNA | 108 | ESTIMATE | −0.136 (0.16) | −0.104 (0.36) | −0.218 (0.17) |
| GSE131907 | 36 | frac_immune | −0.127 (0.46) | −0.284 (0.15) | −0.062 (0.86) |
| GSE148071 TISCH | 36 | immune_fraction | −0.055 (0.75) | +0.037 (0.86) | −0.086 (0.79) |

Bulk ImmuneScore is **not** the CLDN4-low-immune story. CD8 / CXCL13 are. CPTAC LSCC **protein** ImmuneScore is the exception (companion row below).

### CXCL13

TLS extract (`scrna_cldn4_combo`): CXCL13+ T fraction. CPTAC LUAD: TLS chemokine RNA (includes CXCL13; no single-gene CXCL13 column on that table).

| cohort | n | ρ (p) | median δ (p) nH/nL | Q4 vs Q1 δ (p) nH/nL |
|---|---:|---|---|---|
| GSE148071 | 36 | **−0.396 (0.017)** | **−0.497 (0.010)** 18/18 | −0.531 (0.061) 9/9 |
| GSE207422 TLS | 15 | −0.511 (0.052) | −0.429 (0.19) 8/7 | −0.500 (0.34) 4/4 |
| GSE241934_IIT | 11 | **−0.655 (0.029)** | **−0.867 (0.017)** 6/5 | −1.00 (0.10) 3/3 |
| GSE253013 TLS | 9 | −0.400 (0.29) | −0.200 (0.73) 5/4 | −0.556 (0.40) 3/3 |
| GSE131907 TLS | 32 | +0.057 (0.76) | −0.102 (0.64) 16/16 | +0.188 (0.57) 8/8 |
| CPTAC LUAD RNA | 110 | −0.060 (0.54) | −0.037 (0.74) 55/55 | −0.082 (0.61) 28/28 |

GSE148071 and GSE241934_IIT are the CXCL13 rows that move. Median δ is large (−0.50 and −0.87). Q4 vs Q1 on GSE148071 is the same direction (δ=−0.53) with n=9/9.

## Protein companion (not the RNA gate)

| cohort | endpoint | n | ρ (p) | median δ (p) | Q4 vs Q1 δ (p) |
|---|---|---:|---|---|---|
| CPTAC LSCC protein | ImmuneScore | 78 | **−0.432 (7.9×10⁻⁵)** | −0.327 (0.013) 39/39 | **−0.620 (0.00084)** 20/20 |
| CPTAC LSCC protein | CD8 | 78 | −0.235 (0.039) | −0.198 (0.13) | −0.310 (0.096) |
| CPTAC LUAD protein | CXCL13/TLS | 79 | −0.200 (0.077) | −0.278 (0.034) 40/39 | −0.245 (0.19) |
| CPTAC LUAD protein | ImmuneScore | 79 | −0.168 (0.14) | −0.244 (0.063) | −0.160 (0.39) |

CPTAC LSCC CLDN4 **protein** Q4 vs Q1 vs ImmuneScore is the largest single contrast in the folder (δ=−0.62, n=20/20). RNA CLDN4 in the same cases is weaker. TACSTD2–CLDN4 protein ρ on LSCC is only +0.08 (companion).

## TACSTD2 companion (not a gate)

CLDN4 vs TACSTD2 Spearman on the same tables (positive co-expression; not used to enter the high group):

| cohort | n | CLDN4–TACSTD2 ρ |
|---|---:|---|
| GSE218989 | 355 | +0.610 |
| GSE285029 | 234 | +0.526 |
| GSE131907 | 36 | +0.612 |
| GSE148071 | 36 | +0.507 |
| GSE207422 A3 | 12 | +0.538 |
| GSE253013 | 9 | +0.683 |
| GSE241934_IIT | 11 | +0.945 |
| CPTAC LUAD RNA | 110 | +0.392 |
| CPTAC LSCC RNA | 108 | +0.446 |

Patients can be CLDN4-high without being TACSTD2-high. The gates above do not require it.

## What is large vs what is not

**Large / directional (CLDN4-high = immune-low):**

- GSE218989 CD8A Q4 vs Q1: δ=−0.307, p=4.1×10⁻⁴, n=89/89
- GSE148071 CXCL13+ T median: δ=−0.497, p=0.010, n=18/18 (ρ=−0.396, p=0.017)
- GSE131907 T/NK median: δ=−0.444, p=0.024, n=18/18
- CPTAC LUAD RNA CD8 median: δ=−0.256, p=0.021, n=55/55
- GSE241934_IIT CXCL13+ T median: δ=−0.867, p=0.017, n=6/5
- CPTAC LSCC protein ImmuneScore Q4 vs Q1: δ=−0.620, p=8.4×10⁻⁴, n=20/20

**Not the story:** full-pool T/NK / CD8A / ImmuneScore Spearman; GSE218989 / GSE285029 ImmuneScore; GSE241934_Real T/NK (opposite sign); GSE207422 A3 CLDN4 vs T/NK (null; TACSTD2 on that table remains the given A3 row).

## Figures

- `figures/forest_q4q1_cliff.png` — Q4 vs Q1 Cliff δ, all primary endpoints
- `figures/forest_median_cliff.png` — median gate
- `figures/forest_q4q1_cxcl13.png` / `forest_median_cxcl13.png` — CXCL13 only
- `figures/forest_q4q1_protein.png` — CPTAC protein companion
- `figures/violin_q4q1_bulk_cd8_immune.png` — GSE218989 CD8A / ImmuneScore, GSE285029 CD8A, CPTAC LUAD CD8
- `figures/violin_q4q1_tnk_cxcl13.png` — GSE131907 T/NK, GSE148071 / GSE207422 CXCL13+ T, CPTAC LUAD TLS
- `figures/violin_median_selected.png` — median-split violins

Q4 vs Q1 violins are the extra figures that make the contrast large. Full table: `tables/per_cohort_gates.tsv`.

Reproduce: `python3 methods/cldn4_barrier_immune/analyze.py`

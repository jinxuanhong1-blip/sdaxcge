# FINDING — scCODA-style composition, CLDN4 only

**Question.** Do **T/NK** or **B** fractions drop in patients with high
malignant CLDN4, on GSE207422 and already-extracted public lung ICI scRNA,
using a Dirichlet-multinomial / ALR composition model and honest n?

**Answer.** **None**, after BH on the pre-specified primary family
(12 patient-level permutation tests: T/NK and B × 6 ICI slices). Several
slices point in the hypothesized direction; none survive q < 0.10.
Bonferroni on the same 12 tests is 1 for every row.

TACSTD2 is **never a gate**. Prior TACSTD2 scCODA (PR #283 / #284) is taken
as given and is not re-audited.

## Primary contract

- Exposure: within-slice median split of malignant (or author-epithelial)
  CLDN4 mean. High = score ≥ slice median.
- Test: ALR slope of T/NK or B vs a non-lymphoid reference, **patient
  permutation p** (exact enumeration when `C(n, n_high) ≤ 3000`).
- Family: TNK + B × {GSE207422 post A3-malignant, GSE241934 IIT,
  GSE241934 RWC, GSE291670, GSE205335, GSE253013 tumor}.
- Recovery: ALR effect < 0 and BH q < 0.10.
- scCODA HMC: `not_available` (no anndata / tensorflow). Not faked.
- DM two-group LRT: diagnostic only. It scales with n_cells and is
  anti-conservative.

## Honest n

Unit = patient. Cells are library size only.

| Slice | Extracted | CLDN4-eligible | high vs low | CLDN4 definition | Notes |
|---|---:|---:|---|---|---|
| GSE207422 post A3-mal | 15 | **10** | 5 vs 5 | A3-malignant mean log1p(CP10k) | 3 pre excluded; P11/P14 empty malignant; P06/P13/P15 have 1–4 malignant cells |
| GSE207422 DRMref | 12 | 12 | 6 vs 6 | DRMref malignant mean | Sensitivity, not in the primary family |
| GSE241934 IIT | 11 | **11** | 6 vs 5 | author epithelial CLDN4 | EGFR-mut NEOTIDE; CopyKAT IDs not public |
| GSE241934 RWC | 34 | **29** | 15 vs 14 | author epithelial CLDN4 | drop n_epi < 10 (P33, P345, P348, P481, P579) |
| GSE291670 | 6 | **6** | 3 vs 3 | marker-malignant CLDN4 | exact Wilcoxon / 3 vs 3 floor p = 0.10 two-sided |
| GSE205335 | 22 | **22** | 11 vs 11 | author-malignant CLDN4 | B = B+plasma only |
| GSE253013 tumor | 9 | **9** | 5 vs 4 | marker-malignant-like CLDN4 | B from TLS extract; 9.3 GB GEO RDS not downloaded |
| ICI pool (cohort-centered) | — | **87** | 45 vs 42 | mixed annotations | Secondary; not independent of the singles |

Not included (already extracted, but not public lung **ICI** scRNA, or no
malignant CLDN4): GSE131907 (treatment-naive atlas); leftover GSE267108 /
GSE274595 (malignant CLDN4 empty).

## Composition table

Patient-level counts, fractions, CLDN4 score, and high/low call:

[`results/composition_table.tsv`](results/composition_table.tsv)

Long cell-type counts: [`results/composition_long.tsv`](results/composition_long.tsv).

### GSE207422 post, A3-malignant CLDN4 (primary n=10)

| Patient | Response | n_cells | n_mal | n_TNK | n_B | frac_TNK | frac_B | CLDN4 | high |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| P02 | NMPR | 5345 | 10 | 2580 | 1970 | 0.483 | 0.369 | 1.77 | 1 |
| P03 | MPR | 9259 | 746 | 4905 | 1050 | 0.530 | 0.113 | 1.02 | 0 |
| P04 | NMPR | 8398 | 51 | 5345 | 458 | 0.636 | 0.055 | 1.33 | 1 |
| P06 | MPR (pCR) | 4649 | 1 | 2809 | 94 | 0.604 | 0.020 | 0.00 | 0 |
| P07 | NMPR | 9131 | 3209 | 984 | 72 | 0.108 | 0.008 | 1.73 | 1 |
| P09 | NMPR | 4849 | 164 | 3133 | 141 | 0.646 | 0.029 | 1.86 | 1 |
| P10 | NMPR | 6643 | 189 | 3876 | 90 | 0.583 | 0.014 | 1.01 | 0 |
| P12 | NMPR | 6233 | 306 | 1272 | 1233 | 0.204 | 0.198 | 0.61 | 0 |
| P13 | NMPR | 6753 | 3 | 1286 | 117 | 0.190 | 0.017 | 2.54 | 1 |
| P15 | NMPR | 8086 | 4 | 1420 | 376 | 0.176 | 0.047 | 1.04 | 0 |

P11 and P14 (both MPR) have 0 A3-malignant cells and are not in this n=10.
The three pre-treatment biopsies are not in this n=10.

### Median fractions, CLDN4-high vs low

| Slice | n | T/NK high | T/NK low | B high | B low |
|---|---:|---:|---:|---:|---:|
| GSE207422 post A3-mal | 10 | 0.483 | 0.530 | 0.029 | 0.047 |
| GSE241934 IIT | 11 | 0.498 | 0.544 | 0.135 | 0.252 |
| GSE241934 RWC | 29 | 0.482 | 0.394 | 0.151 | 0.214 |
| GSE291670 | 6 | 0.031 | 0.098 | 0.003 | 0.005 |
| GSE205335 | 22 | 0.378 | 0.392 | 0.051† | 0.074† |
| GSE253013 tumor | 9 | 0.244 | 0.380 | 0.147 | 0.041 |
| ICI pool | 87 | 0.453 | 0.403 | 0.114 | 0.095 |

† GSE205335 B is B+plasma.

## Primary ALR tests (patient permutation)

| Slice | Compartment | n | high vs low | ALR effect | perm p | BH q | Recovery |
|---|---|---:|---|---:|---:|---:|---|
| GSE207422 post A3-mal | TNK | 10 | 5 vs 5 | +0.20 | 0.76 | 0.98 | no |
| GSE207422 post A3-mal | B | 10 | 5 vs 5 | +0.04 | 0.98 | 0.98 | no |
| GSE241934 IIT | TNK | 11 | 6 vs 5 | −0.74 | 0.40 | 0.98 | no |
| GSE241934 IIT | B | 11 | 6 vs 5 | −1.33 | 0.15 | 0.80 | no |
| GSE241934 RWC | TNK | 29 | 15 vs 14 | +0.10 | 0.87 | 0.98 | no |
| GSE241934 RWC | B | 29 | 15 vs 14 | −0.10 | 0.84 | 0.98 | no |
| GSE291670 | TNK | 6 | 3 vs 3 | −1.24 | 0.20 | 0.80 | no |
| GSE291670 | B | 6 | 3 vs 3 | −0.84 | 0.60 | 0.98 | no |
| GSE205335 | TNK | 22 | 11 vs 11 | −0.02 | 0.95 | 0.98 | no |
| GSE205335 | B | 22 | 11 vs 11 | +0.27 | 0.51 | 0.98 | no |
| GSE253013 tumor | TNK | 9 | 5 vs 4 | −0.10 | 0.92 | 0.98 | no |
| GSE253013 tumor | B | 9 | 5 vs 4 | +1.39 | 0.17 | 0.80 | no |

Full rows (OLS, fraction MWU, Spearman): [`results/tests_primary.tsv`](results/tests_primary.tsv).

## Closest unadjusted cells (not recoveries)

| Slice × compartment | ALR | perm p | BH q | n | Note |
|---|---:|---:|---:|---:|---|
| GSE241934 IIT × B | −1.33 | 0.145 | 0.80 | 11 | Naive fraction MWU p=0.0087 and Spearman ρ=−0.61 p=0.047; **not** the primary p |
| GSE291670 × TNK | −1.24 | 0.20 | 0.80 | 6 | 3 vs 3; exact ALR cannot go below 1/20 = 0.05 |
| GSE253013 tumor × B | +1.39 | 0.17 | 0.80 | 9 | Opposite of hypothesized B-down |

GSE207422 A3-malignant (the named series) is near null and slightly
**opposite** on ALR (T/NK +0.20, B +0.04). Continuous Spearman on the same
10 patients is ρ=−0.10 vs T/NK and vs B (p=0.78), matching the given
CLDN4-only GSE207422 slice (PR #318 ρ=−0.13 vs T/NK, n=10).

The ICI pool (n=87, mixed annotations, cohort-centered CLDN4) is null
(T/NK ALR −0.05 p_perm=0.86; B +0.06 p_perm=0.86).

DRMref sensitivity on GSE207422 (n=12, includes P11/P14) is also null
(T/NK ALR −0.08 p=0.93; B/TLS −0.60 p=0.64).

## What this does *not* say

- It does not refute a biological CLDN4–exclusion hypothesis; it says the
  already-extracted public lung ICI scRNA grid does not recover T/NK or B
  down after honest FDR / n.
- It does not use author CopyKAT IDs (not public on GSE207422).
- GSE241934 CLDN4 is author-epithelial, not a CopyKAT malignant call.
- GSE253013 B counts come from the TLS extract, not the 9 GB RDS.
- The ICI pool is not an independent replication of any single series.
- Cell-level p-values are not reported (pseudoreplication).
- scCODA credible effects are **not** claimed.

## Files

- `results/composition_table.tsv` — the composition table (done criterion)
- `results/composition_long.tsv`
- `results/honest_n.tsv` / `results/tests_primary.tsv` / `results/tests_full.tsv`
- `results/figures/fig1_fractions_high_vs_low.png`
- `results/figures/fig2_honest_n.png`
- `results/figures/fig3_alr_forest.png`
- `data/PROVENANCE.md`

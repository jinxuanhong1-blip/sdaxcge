# Methods — concordant-4 decoupler / DoRothEA

CLDN4 only. Patient / donor / sample is the unit. Quartile labels are taken from PR #503 (`data/sample_inventory.tsv`) and are not recomputed.

## Contrast

Malignant pseudobulk UMI sums from GSE123902 (marker-malignant donors), GSE131907 (author-malignant tumor-bearing samples), GSE205335 (author-malignant patients), and GSE189357 (marker-malignant patients). P4001 is absent from the count matrix.

Primary contrast: within-cohort CLDN4 %pos Q4 vs Q1, stacked, OLS with cohort indicators (reference GSE123902). Honest n = 18 low / 16 high. GSE189357 contributes 3/2 and is not tested alone.

Companion: the same activity vs continuous CLDN4 %pos (z-scored, n=64).

This is not the PR #454 within-patient cell-tail split.

## Expression

Genes are the intersection of the four matrices, then genes with count ≥10 in ≥3 units. Normalization is log2(TMM-CPM+1), with the same TMM implementation as PR #503, fit on all 64 units. The Hallmark IFN and MHC-I family scores on this matrix reproduce PR #503 (−0.584 and −0.779).

## Network

DoRothEA human A+B+C from OmniPath, as returned by `decoupler.op.dorothea` (academic license). Weights are mode of regulation divided by the confidence level: A = ±1, B = ±1/2, C = ±1/3. A TF is scored only if ≥5 targets are present (`tmin=5`).

Absent from A+B+C, and not imputed: NLRC5, CIITA, RFXANK, RFXAP, GRHL1, GRHL3, OVOL1, OVOL2. IRF7 has 2 targets in this matrix and is not scored. A ±1 weight sensitivity (signs only) is stored in `results/ulm_sign_only_weights.tsv`.

CollecTRI is a sensitivity network only. It still lacks NLRC5. It is not pooled with DoRothEA.

## Activity

decoupler 2.2.0, samples × genes = the log2(TMM-CPM+1) matrix.

| Method | What the score is |
|---|---|
| **ULM (primary)** | t-value of a univariate regression of expression on that TF’s weights |
| MLM | t-value from one multivariate regression on all TF weights |
| VIPER | aREA normalized enrichment, pleiotropy correction on |
| wmean | Σ(w·x) / Σ\|w\|, no permutation (`times=1`) |

ULM is invariant to a per-sample shift or scale of expression. The weighted mean is not, which is why a cell-level wmean can track depth (PR #454).

## Tests

For each TF and method, OLS of the per-unit activity on `CLDN4_Q4` plus cohort indicators. Positive coefficient = higher activity in CLDN4-high. BH FDR is computed across the 297 scored TFs within each method and is reported so the pre-specified p-values are not mistaken for a screen.

Within-cohort ULM refits TMM inside the cohort. Leave-one-cohort-out drops that cohort’s units from the stacked OLS and keeps the joint activity scores.

Gene-program scores are the mean log2(TMM-CPM+1) of the PR #503 Hallmark IFN union, the custom MHC-I/APM list, KEGG∪GO tight junction with CLDN4 removed, and Hallmark OXPHOS.

## Focused TFs

Pre-specified. IFN / MHC regulators are labeled against the KD expectation (lower in CLDN4-high). Junction TFs are reported with no expected KD sign.

| Program | TFs |
|---|---|
| IFN | STAT1, STAT2, IRF1, IRF2, IRF3, IRF7, IRF8, IRF9 |
| MHC | RFX5; NLRC5 and CIITA only as transcripts, or CollecTRI when present |
| junction | GRHL2, ELF3, KLF4, TFAP2A |

## Specification sweep

`scripts/sweep.py` repeats the activity estimation on the same log2(TMM-CPM+1) matrix. Networks are DoRothEA confidence tiers A, A+B, A+B+C, and A–D (D weight = ±1/4), CollecTRI, PROGENy top 500, and unsigned IFN / TJ / apical-junction gene sets. Methods are ULM, MLM, weighted mean, weighted sum, and permutation-normalized weighted mean (`times=100`) for A+B+C and CollecTRI. MLM was not fit on A–D. Cutoffs are within-cohort: the locked Q4/Q1 labels, equal-count quartiles, tertiles, halves, and the top and bottom 30% and 20% (a cohort enters only if each tail has at least 2 units; a test is reported only if each tail has at least 3). Each cutoff is tested as a stacked OLS with cohort indicators and as an inverse-variance meta of the cohorts that clear that floor. The meta treats cohort standard errors as known, so its p-values run smaller than the stacked model and are not a second independent cohort.

The sweep table is the record of that search. The highlighted joint panel is the stacked-or-meta row with both nominal p<0.05 and the smallest max(p). It is not a pre-specified primary test. All-epithelial matrices are not in this folder.

## Software

Python, decoupler 2.2.0, numpy, pandas, scipy, matplotlib. No R Bioconductor decoupleR. The Python port is the implementation that was run.

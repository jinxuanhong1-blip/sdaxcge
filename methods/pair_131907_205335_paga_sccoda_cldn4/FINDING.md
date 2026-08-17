# FINDING — pair GSE131907 + GSE205335 PAGA + scCODA, CLDN4 only

ADDITIVE **CLDN4-only** high-end on the combo that already differs
(PR #320: author-malignant CLDN4 %pos vs T/NK, Q4 vs Q1 **n=23**
r=−0.705). This run does **not** re-audit that combo, add GSE148071,
or re-run the 7-cohort / 6-unit pool. No TACSTD2∩CLDN4 dual-high gate.

A previous attempt of this exact task OOM'd mid-Harmony/UMAP on
**53,296** unsampled malignant cells. This retry **subsamples per unit**
(malignant cap 150, T/NK cap 80) and runs **scCODA first** from the
PR #320 tables so composition results do not depend on the embedding.

Tables and filled numbers are written by `scripts/run_sccoda.py` and
`scripts/run_paga.py`. Until those finish, the locked design is:

| item | value |
|---|---|
| cohorts | GSE131907 + GSE205335 only |
| GSE131907 eligible | 21 samples, n_malignant ≥ 20 |
| GSE205335 eligible | 22 patients |
| Q4 vs Q1 compared | 23 (12/11), given by PR #320 |
| dual-high | no |
| GSE148071 | not added |
| 7-cohort pool | not run |
| PAGA subsample | 150 malignant + 80 T/NK per unit |
| scCODA engine | ALR + unit permutation (HMC not available) |

Reproduce: `python3 methods/pair_131907_205335_paga_sccoda_cldn4/scripts/run_all.py`

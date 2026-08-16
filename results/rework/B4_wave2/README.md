# B4 extras: public lung ICI TJ / CLDN4 + TCGA TJ vs CD8/GEP

The GSE126044 TJ / NR result is the **index B4 cohort** and is taken as given.

This folder keeps that index computation plus **additive** paper figures from other public lung ICI matrices (not GSE126044) and TCGA-LUAD/LUSC.

**Paper extras:** [`extra/EXTRA_FOR_PAPER.md`](extra/EXTRA_FOR_PAPER.md)

## Index cohort (GSE126044), for context

Public GEO counts only. n=16 (5 responder / 11 non-responder). Score = log2(CPM+1) for single genes; modules = mean of per-gene z-scores.

| feature | n R / NR | median R | median NR | two-sided MW p | direction |
|---|---|---|---|---|---|
| TJ 7-gene (`CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN`) | 5 / 11 | −0.133 | 0.129 | 0.01923 | NR > R |
| CLDN4 | 5 / 11 | 2.54 | 3.77 | 0.115 | NR > R |
| CLDN1/4/7/F11R/PARD3 | 5 / 11 | −0.018 | 0.188 | 0.320 | NR > R |
| CD8A | 5 / 11 | — | — | 0.00092 | R > NR |

Full GSE126044 tables: `tables/`, `figures/`, `summary.json`.

## Extra public lung ICI + TCGA

```bash
python3 -m pip install -r scripts/rework/B4_wave2/requirements.txt
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/extra_cohorts.py
```

Outputs: `extra/tables/`, `extra/figures/`, `extra/EXTRA_FOR_PAPER.md`.

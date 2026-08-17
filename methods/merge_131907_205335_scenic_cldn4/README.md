# merge_131907_205335_scenic_cldn4

ADDITIVE **CLDN4-only** SCENIC/GRN proxy on the merged **GSE131907 + GSE205335**
author-malignant combo that already differs (PR #320).

Question: do CLDN4-high malignant cells show different IFN / MHC-I / TJ / keratin
regulons vs CLDN4-low, at the **patient** level?

A10 ELF3–CLDN4 is given and is not rediscovered as the headline.
No dual-high TACSTD2×CLDN4. No GSE207422.

Full pySCENIC cisTarget is **not** run (no motif DBs). The implemented method is
documented AUCell + public TF–target priors (TRRUST / DoRothEA / CollecTRI).

## Reproduce

```bash
python3 methods/merge_131907_205335_scenic_cldn4/scripts/download.py
python3 methods/merge_131907_205335_scenic_cldn4/scripts/analyze.py
```

Priors and gene sets are committed under `resources/`. GEO matrices stay in
`/tmp` (not git).

Primary table: `tables/regulons.tsv` (patient n). Writeup: `FINDING.md`.

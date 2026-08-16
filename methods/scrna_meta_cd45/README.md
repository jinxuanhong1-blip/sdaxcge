# Immune-compartment extra: TACSTD2/CLDN4 leak in CD45+ lung IO scRNA vs MPR

Additive only. User A3 (epithelial / malignant TACSTD2–T/NK) is taken as given.

This folder scores **immune libraries**, not malignant RNA. TACSTD2/CLDN4 detection
in CD45+ matrices is treated as residual/ambient/epithelial leak.

```bash
python3 run_all.py
```

Public inputs stay under `/workspace/data/` (gitignored). Results go to
`results/extra/scrna_meta_cd45/`.

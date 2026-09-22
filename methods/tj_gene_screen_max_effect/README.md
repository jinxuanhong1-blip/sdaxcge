# TJ gene screen max-effect (Part2)

Maximize |ρ| / win-margin of CLDN4 vs other TJ/claudin genes on
concordant-4 + TCGA + CosMx. Part2-ready; do not force into Part1.

## Quick start

```bash
python3 methods/tj_gene_screen_max_effect/scripts/sweep_concordant4.py
TCGA_CACHE=/tmp/tcga python3 methods/tj_gene_screen_max_effect/scripts/sweep_tcga.py
python3 methods/tj_gene_screen_max_effect/scripts/download_cosmx.py
COSMX_H5AD=methods/tj_gene_screen_max_effect/data/cosmx_nsclc/cosmx_human_nsclc_clustered.h5ad \
  python3 methods/tj_gene_screen_max_effect/scripts/sweep_cosmx.py
python3 methods/tj_gene_screen_max_effect/scripts/make_verdict.py
```

Read `FINDING.md`. CosMx h5ad is gitignored (~2.6 GB).

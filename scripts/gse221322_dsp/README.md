# GSE221322 GeoMx protein DSP (barrier vs immune)

Public NanoString GeoMx nCounter protein on an immunotherapy-treated NSCLC TMA (Monkman et al., *Immunology* 2023; GEO [GSE221322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221322)).

CLDN4 and TROP2 are not on the deposited 68-plex. The script inventories the panel and scores **EpCAM** (TACSTD1) and **PanCk** against CD8 / CD3 / PD-1 and neighborhood immune proteins.

```bash
bash scripts/gse221322_dsp/00_fetch.sh
python3 scripts/gse221322_dsp/analyze.py
```

Outputs land in `results/gse221322_dsp/`. Raw GEO files stay in `data/gse221322/` (gitignored).

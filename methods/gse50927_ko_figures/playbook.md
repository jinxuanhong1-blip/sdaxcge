# Methods: GSE50927 Cldn4 KO extra figures

**Additive only.** High-contrast redraw of locked GSE50927 numbers. No new *p*. Not lung cancer.

## In / out

| In | Out |
|---|---|
| Frozen tables copied from `cldn4_ko_gsea` + `fable_cldn4_kdko` | Re-running GSEA or edgeR |
| Author edgeR logFC / deposited FDR | Invented *p* |
| Naive whole-lung Cldn4 KO vs WT | VILI-only contrasts as primary |
| IFN/MHC 34/39 panel as already matched | New ortholog maps |

## Figures

`python3 scripts/gse50927_ko_figures/plot_extra.py` writes `figures/fig_extra1_…` through `fig_extra6_…` (PNG+PDF).

## Locked numbers

See `tables/locked_stats.json`. IFN-γ NES +1.508 is the already-published +1.51.

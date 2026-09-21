# CosMx bivariate Moran / Lee / Gi* — CLDN4 vs CD8A

Public He et al. 2022 CosMx NSCLC object on figshare [25976224](https://doi.org/10.6084/m9.figshare.25976224) (`cosmx_human_nsclc_clustered.h5ad`). Eight sections, five patients. No private 8-KL matrices.

This layer is neighbor-graph autocorrelation of the **CLDN4 and CD8A transcript fields**. It is a different estimand from CLDN4-high versus CLDN4-low tumor-cell cytotoxic **cell counts**.

## Estimands

Expression is `log1p(counts / n_counts × 10,000)` from `layers/counts`. The object `X` matrix is not used. Coordinates are `obsm['spatial']` in pixels, converted at 0.18 µm/px.

Within each FOV, and again on each whole section:

- Bivariate Moran's I: CLDN4 with the row-standardized lag of CD8A (esda `Moran_BV` definition).
- Lee's L (Lee 2001) on row-standardized weights.
- The same two statistics after within-unit OLS residualization of both genes on KRT8 and EPCAM.
- CLDN4 residual only (CD8A left raw).
- KRT8–CD8A and EPCAM–CD8A Moran on the same CD8A permutations (epithelial baselines).
- Reverse Moran: CD8A with the lag of CLDN4.
- Within-tumor Pearson r of CLDN4 (raw, and residualized inside the tumor set) with the CD8A lag. Tumor labels are the author classes `tumor 5/6/9/12/13`. The null shuffles CLDN4 among tumor cells.
- Tumor-cell CLDN4 top-minus-bottom quartile difference in mean neighbor CD8A.
- Getis-Ord Gi* z with binary star weights, for raw expression and for residuals, plus CLDN4-hot / CD8A-cold overlap.

Neighbor graphs, built inside the unit (isolates dropped for Moran and Lee): kNN k=6, 15, 30; radius 20, 50, 100 µm; Delaunay edges trimmed at the 99th percentile of edge length. Primary graph: **50 µm**. Co-primary: **kNN k=15**.

199 within-unit permutations, seed 20260921, independent of job order. Section Wilcoxon tests and the five-patient sign test are the confirmatory summaries. FOV inverse-variance meta-analysis assumes independent FOVs and is descriptive.

## Reproduce

```bash
bash methods/cosmx_bivar_autocorr/download.sh
python3 methods/cosmx_bivar_autocorr/analyze.py --self-test
python3 methods/cosmx_bivar_autocorr/analyze.py --h5ad data/cosmx_human_nsclc_clustered.h5ad
```

The h5ad is about 2.6 GB and is not committed. A gene/coordinate cache is written to `/tmp/cosmx/expr_cache_v1.npz`.

Numbers are in `RESULTS.md`. Tables are under `tables/`, figures under `figures/`.

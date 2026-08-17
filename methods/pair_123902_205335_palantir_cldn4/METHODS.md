# Methods — pair GSE123902+GSE205335 REAL Palantir (CLDN4-only)

## Scope

- Accessions: GSE123902 (Laughney 2020) + GSE205335 (Ahn/Lee 2024) only.
- Not GSE148071. Not the seven-cohort pool. Not a TACSTD2∩CLDN4 dual-high gate.
- Estimand: CLDN4. Barrier score holds CLDN4 out. IFN = Hallmark IFNα ∪ IFNγ
  (same family as the given pair IFN DE logFC −1.05).
- Inferential unit: GSE123902 **donor** / GSE205335 **patient**.

## Epithelium

- GSE123902: marker epithelium `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.
  One tumor/met library per donor (PRIMARY preferred). NORMAL libraries are
  kept only as the Palantir **root pool**.
- GSE205335: author `lineage.total == Epithelial cells`. Normal tissues dropped.
- Cap ≤250 cells / unit after a fixed seed.

## Palantir (real package)

1. log-normalize, 3000 HVG, 30 PCs, Harmony on `dataset`, Leiden 0.6.
2. `import palantir` is required. No DPT fallback.
3. Diffusion maps (10 components) + multiscale space on Harmony/PCA.
4. **Early cell is not CLDN4-high.** Prefer GSE123902 NORMAL cells in the
   bottom CLDN4 tertile of that pool, farthest from the malignant centroid.
5. Terminals are Palantir-auto. They are **not** the CLDN4-high quantile.
   Post-hoc labels (AT2 / malignant-like / IFN / barrier) are descriptive.
6. Waypoints = 500; knn = 30.

## Destinies

Each cell gets Palantir pseudotime, entropy, and a fate-probability vector.
Destiny = argmax fate. Unit means are the claim. Trends are unsmoothed
log-normalized score means in pseudotime bins on each destiny.

## Honest n

Report cells, tumor units, NORMAL-root cells, eligible units (≥10 cells),
paired tertile units (≥8 cells / arm), and Palantir version. Cell-level
p-values are not the claim.

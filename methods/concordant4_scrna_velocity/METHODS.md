# Methods

## Question

Do concordant-4 epithelial cells move toward a CLDN4-high barrier state and away from IFN / antigen-presentation programs?

## Object

Tumor epithelial cells from the locked 65 units, capped at 160 cells per unit (seed 4). Root-pool cells are GSE123902 NORMAL epithelium and GSE131907 nLung epithelium. They orient pseudotime and are excluded from the Spearman n.

Gates: GSE131907 and GSE205335, author epithelial label. GSE123902 and GSE189357, PTPRC = 0 and at least one of EPCAM, KRT8, KRT18, KRT19 detected. Counts are log1p of CP10k on the genes shared by all four matrices (11,805 genes).

Barrier score, CLDN4 excluded: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR. TJ score, CLDN4 excluded: OCLN, TJP1, TJP2, TJP3, CLDN3, CLDN7, CDH1, F11R, MARVELD2, CGN, CRB3. IFN: Hallmark IFNα ∪ IFNγ (211 genes present). APM: HLA-A/B/C, B2M, TAP1, TAP2, TAPBP, PSMB8/9/10, NLRC5.

Scores: gene-wise z-mean, UCell-like relative rank, AUCell-like recovery in the top 5% of genes, and Seurat-style AddModuleScore (24 expression bins).

## Clocks (expression only)

Neighbors on PCA of 1,000 or 2,000 HVGs, with or without CLDN4 in the HVG list, 20 or 40 PCs, 15 or 30 neighbors, Harmony on dataset or not.

- Palantir 1.4.5 (400 waypoints) on the Harmony multiscale diffusion map, early cell = the IFN-high root. One run, not the full grid.
- DPT from a root cell.
- PAGA shortest-path distance from the root cell's Leiden cluster (resolution 0.5). This is the principal-graph ordering. R / Monocle3 was not installed.
- MST on the cluster graph. One leaf is the highest-barrier cluster. Another leaf is the highest-CLDN4 cluster (that choice makes the CLDN4 endpoint partly circular). Units with no cells on the path are dropped.
- CytoTRACE-like: pseudotime = −(genes detected), smoothed on the kNN.
- Absorption probability of a random walk onto one terminal cluster.

Roots: (1) highest AT2 in the uninvolved-lung pool, (2) highest IFN among cells with CLDN4 at or below the tumor median, (3) lowest CLDN4 among cells with AT2 at or above the median.

The reported test is Spearman of unit means, n up to 65, not a cell-level p-value. A row matches the thesis when ρ(CLDN4, PT) > 0, ρ(barrier, PT) > 0, and ρ(IFN, PT) < 0. The headline is the matching row with all 65 units and the smallest Fisher combination of the three one-sided p-values, excluding leaves that were defined as CLDN4-high. BH q-values are computed inside the 1,920-row grid.

Within-unit contrasts split each unit at median, tertile, quartile, and the outer 20/25/30% of CLDN4 and test the paired median difference with a one-sided Wilcoxon.

## Splicing

`inventory.py` records that GEO supplements are single total-count matrices. scVelo needs spliced and unspliced layers. Those layers are not deposited. GSE131907 (PRJNA545296) and GSE205335 (PRJNA844398; authors cite EGA EGAD00001008703) have zero public SRA runs.

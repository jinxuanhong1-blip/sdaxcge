# Methods — GSE148071 PAGA scored by CLDN4

Additive. Does not rewrite `methods/scrna_paga/` or `methods/scrna_paga_cldn4/` (those are GSE131907).

**Question.** On public advanced-NSCLC epithelium (GSE148071), where does **CLDN4** sit on PAGA connectivities and Alveolar/AT2-rooted diffusion pseudotime, relative to AT2 / club / basal / barrier-keratin / malignant-like programs?

**Not a TACSTD2 redo.** TACSTD2 is a comparator Spearman only. The extra figure is gated on CLDN4. Barrier/keratin **excludes CLDN4**.

| Item | Choice |
| --- | --- |
| Cohort | GSE148071 (Wu et al. 2021); 42 advanced NSCLC biopsies |
| Counts | Public GEO raw UMI (`GSE148071_RAW.tar`) |
| Labels | TISCH2 major-lineage ∈ {Malignant, Alveolar, Basal, Epithelial} |
| Graph | log1p(CP10k) → Seurat HVG 3000 → PCA → k-NN (30/30) |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Cap | max 500 cells/patient if n>20k (seed 0) |
| Root | TISCH Alveolar (median AT2 score). Never CLDN4-high. |
| Score | CLDN4 continuous + tertiles |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | Patient (one sample each). Cell-level ρ is descriptive. |
| Extra figure | Within-sample CLDN4-high vs low program scores (min 8 cells/arm) |
| Unused | TISCH expression.h5 (1.1 GB; labels only). No ICI labels. No GEO histology. |

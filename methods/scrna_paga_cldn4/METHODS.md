# Methods — PAGA on epithelium scored by CLDN4

Additive. Does not rewrite `methods/scrna_paga/` (TACSTD2-primary) or `methods/trajectory/`.

**Question.** On public LUAD epithelium, where does **CLDN4** sit on PAGA connectivities and AT2-rooted diffusion pseudotime, relative to AT2 / club / basal / barrier-keratin / malignant-like programs?

**Not a TACSTD2 redo.** TACSTD2 is a comparator Spearman only. The extra figure is gated on CLDN4, not TACSTD2. Barrier/keratin **excludes CLDN4**.

| Item | Choice |
| --- | --- |
| Cohort | GSE131907 LUAD; `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung} |
| Counts | Author raw UMI text (public). Not the 2.9 GB log2TPM. |
| Graph | Seurat-v3 HVG 3000 → PCA → k-NN (30/30) recomputed on epithelium |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Root | nLung author AT2 (median AT2 score). Never CLDN4-high. |
| Score | CLDN4 continuous + tertiles on the epithelial object |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | Sample. Cell-level ρ is descriptive. |
| Extra figure | Within-sample CLDN4-high vs low program scores (min 8 cells/arm) |
| Unused | GSE207422 (NSCLC, not LUAD-only); GSE253013 (size cap) |

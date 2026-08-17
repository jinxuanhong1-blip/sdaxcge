# Methods — REAL Slingshot on GSE131907 epithelium/malignant, CLDN4 only

ADDITIVE. Does **not** rewrite `methods/scrna_paga_cldn4/` (PR #325, DPT/PAGA).
Does **not** rewrite `methods/winpair_131907_205335_slingshot_cldn4/` (PR #449,
DPT fallback). GSE207422 and GSE205335 are not added. No dual-high gate.

**Question.** On public GSE131907 LUAD **epithelium / malignant cells only**
(Kim et al., *Nat Commun* 2020, PMID 32385277), where do **CLDN4**, a
CLDN4-excluded barrier/keratin program, and a compact IFN/ISG program sit on
AT2-rooted Slingshot lineages?

## Dataset

| Item | Choice |
| --- | --- |
| Cohort | GSE131907 LUAD, public raw UMI + author annotation |
| Cells | `Cell_type` ∈ {Epithelial cells, Malignant cells} (or author tS1/tS2/tS3) |
| Origins | nLung, tLung, tL/B, mLN, mBrain. PE unlabeled epithelium dropped |
| Unused | Immune, stroma, nLN, GSE205335, GSE207422, 2.86 GB log2TPM, EGA FASTQ |

nLung AT2 is kept so the root is an external arrow, never a CLDN4-high cell.

## Trajectory clock

Slingshot (Street et al. 2018): MST on Leiden centroids in the first 5 PCs,
then Hastie–Stuetzle principal curves per lineage. If Bioconductor `slingshot`
is installed, that engine is used; otherwise the Python equivalent in
`scripts/slingshot_py.py` (same two steps). PAGA is companion geometry.
scanpy DPT is a **sensitivity** clock only and is not the primary claim.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → k-NN 30/30 |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Root | nLung author AT2 (median AT2 score). Never CLDN4-high |
| Score | CLDN4 continuous + tertiles |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| IFN | compact ISG panel in `scripts/gene_sets.py` (**no CLDN4**) |
| Inferential n | GSE131907 `Sample`. Cell-level ρ is descriptive |
| Cap | ≤400 cells / sample after protecting nLung AT2 |
| Extra figures | Within-sample CLDN4-high vs low; CLDN4/barrier/IFN along lineage bins |

## Tests

Primary: sample-level Spearman of mean CLDN4 vs mean Slingshot / AT2 /
barrier (CLDN4 excluded) / IFN (CLDN4 excluded). BH inside that list only.

Sensitivity (not BH): tLung-only, nLung-only, tumor/met (drop nLung),
nLung+tLung-only, and CLDN4 vs DPT.

The pooled Slingshot correlation mixes nLung, tLung, and metastatic sites.
It is **not** a within-tumor progression test.

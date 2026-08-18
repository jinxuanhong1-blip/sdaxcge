# Concordant-4 cell-level atlas (methods / identity)

ADDITIVE. **CLDN4-only.** This is the missing cell-level atlas figure
for how the four public scRNA sets were taken in, integrated, clustered,
and annotated. It is **methods/identity, not a new claim test.**

Datasets ONLY: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.
Not GSE148071 / GSE127465 / GSE154826 / GSE207422. No dual-high.

**PR #503 T/NK ρ and IFN/MHC numbers were not re-run and were not re-audited.**
The locked four-set T/NK result remains n=65, %pos ρ=−0.531. Do not replace
that patient/sample/donor number with a cell count from this atlas.

## Honest n

- **n_units = 65** (13 donors + 21 samples + 22 patients + 9 patients).
- **n_cells used (after QC + cap) = 22653**.
- **n_clusters (Leiden 0.6) = 26**.
- Plot used **all 22653** cells (below the 80k plot cap; no subsample).

| dataset | unit | n_units | n_cells in units (raw / author) | n_cells QC (pre-cap) | n_cells used | cap |
|---|---|---:|---:|---:|---:|---:|
| GSE123902 | donor | 13 | 30126 | 26056 | 4453 | 350 |
| GSE131907 | sample | 21 | 62612 | 62612* | 7350 | 350 |
| GSE205335 | patient | 22 | 80790 | 80790* | 7700 | 350 |
| GSE189357 | patient | 9 | 122373 | 119264 | 3150 | 350 |

\*GSE131907 / GSE205335 public matrices are already author-QC'd. The same min-gene / min-count / mito gates were applied after the 350/unit cap; every capped cell passed, so the honest pre-cap QC n is the author-in-units count.

Annotation (used cells): {'malignant': 7285, 'T': 6620, 'myeloid': 4187, 'B': 1663, 'NK': 1621, 'other': 1277}.

## QC / Harmony / Leiden

- QC: n_genes ≥ 200, n_counts ≥ 500, mitochondrial % < 20.0.
- Cap: ≤350 cells / unit when the QC-pass set is larger (memory).
- Concatenate (inner gene join) → HVG 2000 (Seurat flavor on log-normalized data, `batch_key=dataset`; Seurat v3 skipped — `skmisc` not installed) → PCA 30.
- Harmony: `batch=dataset`, theta=2.0, max_iter=20. Sample was not a second Harmony key.
- Neighbors k=15 on `X_pca_harmony`; UMAP; Leiden resolution **0.6**.
- Public processed matrices only. Author 36.5 GB GSE123902 H5 and GSE131907 log2TPM text were not used.

## Annotation rules

Cluster → label. Do not invent cell types. Mixed clusters are labelled mixed and mapped to **other** on the UMAP.

1. If ≥50% of cells in the cluster carry an **author malignant** label
   (GSE131907 `Cell_subtype==Malignant cells`; GSE205335 `lineage.sub==Malignant cells`) → malignant (author).
2. Else highest mean marker score on log-normalized expression:
   EPCAM/KRT8/KRT18/KRT19/CLDN4 epithelium-malignant; CD3D/CD3E T; NKG7/GNLY/KLRD1 NK;
   CD68/LYZ myeloid; MS4A1/CD79A B; PECAM1/VWF endothelial; COL1A1 fibroblast.
3. If the top score is <0.15 → other. If the top two scores differ by <0.05 and both ≥0.15 → mixed (other).
4. Endothelial and fibroblast clusters are **other** on the coarse UMAP; the fine label is in `annotation.tsv`.
5. GSE123902 and GSE189357 have no author cell-type column; they use markers only.
6. Cluster 15 top genes are plasma-like (MZB1, IGKC, JCHAIN) but the locked marker panel has no plasma set, so the cluster stays **mixed / other**. Cluster 19 is a low-score island (MBP/CRYAB; GSE131907 mBrain oligodendrocytes) and is left as **other**.

## Cluster table

| cluster | n | label | fine | top markers | rule |
|---|---:|---|---|---|---|
| 0 | 3436 | T | T | IL7R,RPS12,RPS15A,RPS27A,RPS27 | T score=1.501 |
| 1 | 2254 | T | T | CCL5,GZMA,CD3D,NKG7,IL32 | T score=1.932 |
| 10 | 862 | malignant | malignant (author) | EPCAM,MDK,ELF3,CLDN3,KRT8 | author_malignant_frac=0.88≥0.50 |
| 11 | 373 | other | fibroblast | DCN,IGFBP7,MGP,CALD1,COL1A2 | COL1A1 score=2.258 |
| 12 | 1663 | B | B | MS4A1,CD37,HLA-DRA,CD74,RPS23 | B score=2.007 |
| 13 | 210 | other | endothelial | RAMP2,GNG11,CLDN5,SPARCL1,EPAS1 | PECAM1/VWF score=1.911 |
| 14 | 795 | T | T | IL32,B2M,CD3D,TRAC,HLA-A | T score=1.888 |
| 15 | 518 | other | mixed endothelial/myeloid | MZB1,IGKC,SSR4,JCHAIN,DERL3 | mixed: endothelial=0.215 vs myeloid=0.175 (Δ<0.05) |
| 16 | 135 | T | T | HMGB2,HMGB1,TUBA1B,TUBB,STMN1 | T score=1.423 |
| 17 | 519 | malignant | malignant (author) | S100A6,NAPSA,MGST1,SLC9A3R2,S100A16 | author_malignant_frac=0.56≥0.50 |
| 18 | 121 | malignant | malignant (marker epithelium) | CAPS,TPPP3,RSPH1,MORN2,PIFO | EPCAM/KRT/CLDN4 score=1.687 |
| 19 | 176 | other | other (low marker scores) | MBP,CRYAB,PTGDS,S100B,GPM6B | top myeloid=0.097<0.15 |
| 2 | 1972 | malignant | malignant (author) | RPLP0,GAPDH,RPL8,KRT8,KRT18 | author_malignant_frac=0.74≥0.50 |
| 20 | 438 | malignant | malignant (author) | KRT19,SEC61G,ERO1A,CST6,S100A6 | author_malignant_frac=0.97≥0.50 |
| 21 | 289 | malignant | malignant (author) | DEFB1,AGR2,TMEM176A,HMGB3,CP | author_malignant_frac=0.98≥0.50 |
| 22 | 223 | malignant | malignant (author) | EFEMP1,SPINT2,SOX4,TACSTD2,IGF2BP2 | author_malignant_frac=0.87≥0.50 |
| 23 | 466 | malignant | malignant (marker epithelium) | MDK,PDCD6,EPHX1,KRT7,CYB5A | EPCAM/KRT/CLDN4 score=1.993 |
| 24 | 320 | malignant | malignant (author) | MAP1B,STMN1,CKB,H1FX,HMGN2 | author_malignant_frac=0.92≥0.50 |
| 25 | 329 | malignant | malignant (author) | RBP1,H3F3A,TUBA1A,STMN1,PSIP1 | author_malignant_frac=0.91≥0.50 |
| 3 | 218 | malignant | malignant (author) | PERP,SFN,KLF5,FXYD3,S100A2 | author_malignant_frac=0.97≥0.50 |
| 4 | 902 | malignant | malignant (marker epithelium) | SFTPB,SFTA2,NAPSA,MUC1,SLC34A2 | EPCAM/KRT/CLDN4 score=1.574 |
| 5 | 1621 | NK | NK | NKG7,CST7,CCL5,GNLY,KLRD1 | NK score=2.992 |
| 6 | 626 | malignant | malignant (author) | MUC1,S100P,SLPI,GPRC5A,TM4SF1 | author_malignant_frac=0.82≥0.50 |
| 7 | 376 | myeloid | myeloid | HPGDS,SRGN,LAPTM4A,CLU,VIM | myeloid score=0.336 |
| 8 | 3181 | myeloid | myeloid | TYROBP,FTL,AIF1,FCER1G,CD68 | myeloid score=2.261 |
| 9 | 630 | myeloid | myeloid | S100A9,S100A8,LST1,LYZ,TYROBP | myeloid score=2.321 |

## What this is not

- Not a re-audit of PR #503 T/NK ρ (n=65) or malignant IFN/MHC DE.
- Not a dual-high TACSTD2∩CLDN4 object.
- Not GSE148071 / GSE127465 / GSE154826 / GSE207422.
- Not evidence that CLDN4 *causes* T/NK exclusion. Identity figure only.
- Cell counts are not the inferential n. Inferential n for the thesis remains the PR #503 units.

## Files

- `results/tables/annotation.tsv` — cluster → label + top markers
- `results/tables/inventory.tsv` — honest n per dataset
- `results/figures/fig_atlas_umap.png` / `.svg` / `.pdf` — Figure 0 / Fig. atlas

Reproduce:

```bash
python3 methods/concordant4_atlas_umap_annotate/download.py
python3 methods/concordant4_atlas_umap_annotate/analyze.py
```

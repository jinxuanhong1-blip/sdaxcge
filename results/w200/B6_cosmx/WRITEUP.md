# B6 leftover: public CosMx NSCLC — CLDN4 vs immune niches

Honest bottom line: **CosMx leftover exists and was analyzed. MERFISH leftover has no processed cell-by-gene table under 2 GB without login (honest skip).**

Claim B6 (CLDN4-high epithelium sits in immune-poor niches) is **supported for a broad RNA-immune neighborhood at single-cell resolution**, **weaker when the neighborhood is CD45 protein**, and **not supported for T-like (CD3/CD8) neighbors** — those correlations are positive, the opposite of T-cell exclusion.

No ICI labels. These are the public NanoString CosMx SMI NSCLC demo samples (He et al. 2022, *Nat Biotechnol*).

## Why this is leftover

PR #67 (`fable_spatial`) and PR #110 (`B6_radius`) already tested B6 on **GeoMx** (GSE271689) and **Visium** (E-MTAB-13530). Both reported |median partial ρ| ≤ 0.06 at spot/hex-ring scale and did not rescue the claim. PR #8 skipped targeted CosMx-1K / Xenium-IO as instructed.

This slice takes the leftover **public CosMx** processed flat files and catalogues **public MERFISH/MERSCOPE lung**.

## MERFISH: honest skip

| Source | Processed cell-by-gene? | Decision |
|---|---|---|
| Vizgen MERSCOPE FFPE IO showcase (lung 1/2) | gated; GCS `vz-ffpe-showcase` returns HTTP 403 | skip — no anonymous download |
| Figshare 25983526 CellCharter MERSCOPE lung (2 patients) | yes, `*.h5ad` **3.67 GB** | skip — over 2 GB leftover budget |
| Zenodo 11198494 Chen 2024 (Pt38/43/67/73) | **no** — transcript tables only (0.26–3.9 GB) | skip — not processed cell-by-gene |

See `catalog.tsv`.

## CosMx data used

NanoString public S3 `nanostring-public-share/SMI-Compressed/` flat files (exprMat + metadata). Pixel size 0.18 µm. QC: ≥20 counts and ≥5 genes; drop `cell_ID==0`.

| Sample | QC cells | FOVs | RNA epithelial | RNA immune | RNA stromal |
|---|---:|---:|---:|---:|---:|
| Lung9_Rep1 | 87,617 | 20 | 41,595 | 12,836 | 32,161 |
| Lung12 | 71,440 | 28 | 25,599 | 9,247 | 35,970 |
| Lung13 | 81,248 | 20 | 26,713 | 19,700 | 34,585 |

Lung5/6/9_Rep2 exist on the same bucket and were not needed (3 samples, 68 FOVs). The 930 MB Giotto tarball was skipped as redundant.

The flat release has **no official cell types**. Compartments are RNA marker-score argmax (immune / epithelial / stromal), intersected with the 960-plex panel. Protein stains (Mean.PanCK, Mean.CD45) are used only for validation and a robustness split.

## A. CLDN4 / TACSTD2 are epithelial

| Sample | CLDN4 epi vs immune AUROC | CLDN4 frac+ epi / immune | TACSTD2 epi vs immune AUROC |
|---|---:|---|---:|
| Lung9_Rep1 | 0.747 | 0.608 / 0.100 | 0.595 |
| Lung12 | 0.609 | 0.321 / 0.097 | 0.609 |
| Lung13 | 0.809 | 0.689 / 0.093 | 0.645 |

RNA epithelial vs PanCK AUROC is usable (0.79 / 0.75 / 0.94). RNA immune vs CD45 AUROC is **weak** (0.57 / 0.59 / 0.62) — the RNA “immune” label is myeloid-heavy and only partly CD45-concordant. That is why the protein-stain split is required.

## B. Spatial: epithelial gene vs local immune fraction

Index = RNA-epithelial cells with ≥5 neighbors. Neighborhood = fraction of other cells inside a **50 µm** radius, computed **inside each FOV** (no cross-FOV stitching). Inference is FOV-level Spearman + Wilcoxon signed-rank on the 20–28 FOV ρ values (same style as PR #67).

### RNA-immune neighborhood (50 µm)

| Sample | n epi | pooled ρ | median FOV ρ | FOVs ρ<0 | Wilcoxon p |
|---|---:|---:|---:|---|---:|
| Lung9_Rep1 | 41,594 | −0.373 | **−0.250** | 19/20 | 3.8e-6 |
| Lung12 | 25,504 | −0.258 | **−0.204** | 28/28 | 7.5e-9 |
| Lung13 | 26,693 | −0.175 | **−0.138** | 20/20 | 1.9e-6 |

25 µm radius is the same sign and similar magnitude (median FOV ρ −0.270 / −0.196 / −0.166). Pre-specified; not tuned.

TACSTD2 vs the same RNA-immune neighborhood is also negative but smaller (median FOV ρ −0.097 / −0.222 / −0.039).

### Protein-stain robustness (PanCK-high ∩ CD45-low index vs CD45-high neighborhood)

| Sample | n index | pooled ρ | median FOV ρ | FOVs ρ<0 | Wilcoxon p |
|---|---:|---:|---:|---|---:|
| Lung9_Rep1 | 27,527 | −0.163 | **−0.170** | 19/20 | 3.8e-6 |
| Lung12 | 18,612 | −0.070 | **−0.100** | 28/28 | 7.5e-9 |
| Lung13 | 21,256 | −0.087 | **−0.074** | 16/20 | 0.0012 |

Still negative, but |ρ| drops toward the Visium band in Lung12/13.

### T-like neighborhood (any of CD3D/E/G, CD8A/B > 0) — opposite sign

| Sample | pooled ρ | median FOV ρ | FOVs ρ<0 | Wilcoxon p |
|---|---:|---:|---|---:|
| Lung9_Rep1 | **+0.318** | **+0.230** | 0/20 | 1.9e-6 |
| Lung12 | **+0.238** | **+0.202** | 1/28 | 7.5e-8 |
| Lung13 | **+0.115** | **+0.084** | 1/20 | 3.8e-6 |

CLDN4-high epithelium is **not** T-cell poor. If anything it sits in slightly T-like-richer 50 µm neighborhoods. This matches PR #67’s weakly **positive** CLDN4 vs T-effector hex rings.

## Honest verdict on claim B6

- **Does not skip.** Public CosMx processed leftover was real.
- **Supports** a moderate, FOV-replicated anti-correlation of epithelial CLDN4 with a **broad RNA-immune** (myeloid-heavy) 50 µm neighborhood. This is larger than the Visium |ρ| ≤ 0.06 band because CosMx is single-cell, not 55 µm spots.
- **Does not support T-cell exclusion.** T-like neighborhoods are positive in 3/3 samples.
- **Protein-stain split shrinks |ρ|** (median FOV −0.07 to −0.17). Do not treat the RNA-immune ρ ≈ −0.25 as a protein-level immune-exclusion effect size.
- **Cannot claim ICI biology.** No treatment, no response labels.
- **MERFISH remains a skip** until a processed cell-by-gene table is public, <2 GB, and anonymously downloadable.

## Reproduce

```bash
# data/ is gitignored; stream-extract flat files from NanoString S3
python3 scripts/b6_cosmx_cldn4_immune.py
```

Requires pandas, numpy, scipy, matplotlib. Raw CosMx tarballs are not committed.

## 中文摘要

CosMx 公开 NSCLC leftover **有处理后的 cell-by-gene，已分析**。MERFISH/MERSCOPE 肺：**无**可匿名下载且 <2 GB 的处理后细胞矩阵，诚实跳过。

上皮 CLDN4 与 50 µm RNA-immune 邻域在 3 个样本、68 个 FOV 上稳定负相关（中位 FOV ρ −0.14 至 −0.25）。换成 CD45 蛋白邻域后 |ρ| 变小。T 样（CD3/CD8）邻域为**正相关**，不支持 T 细胞排斥。无 ICI 标签。

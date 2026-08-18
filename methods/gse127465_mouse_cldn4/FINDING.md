# FINDING — GSE127465 mouse NSCLC / myeloid-rich lung (Cldn4-only)

ADDITIVE public **MOUSE**. **Cldn4 only.** Not a TACSTD2 gate. **Not merged with human GSE127465.**

Zilionis et al., *Immunity* 2019, PMID [30979687](https://pubmed.ncbi.nlm.nih.gov/30979687/). GEO [GSE127465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465). KP1.9 lung adenocarcinoma (Pfirschke et al. 2016) + matched healthy lung. Public processed object is **CD45+** inDrops: 15,939 cells after the authors' QC. Author major types are Neutrophils / B / MoMacDC / T / NK / pDC / Basophils. **Author epithelium = 0.**

Processed mouse matrix is public and <2 GB. Cldn4 is on the gene list. Epithelium is thin; mouse-level Cldn4 is still reported vs T fraction and myeloid. Human MTX / human metadata were not opened.

## Verdict

**GO (mouse-only; epithelium thin).** Public files exist (`GSE127465_mouse_counts_normalized_15939x28205.mtx.gz` 174,301,312 B, metadata 15,939 × 12, 28,205 genes). Cldn4 is present (column 3875 / 0-based index 3874). Author malignant/epithelial barcodes = **0**. Marker leftover (Epcam+ or Cdh1+∩Krt8+) = **270 / 15,939** cells — contamination/ambient in a CD45+ sort, not a tumor-epithelial compartment.

Cldn4 is **present but sparse** in CD45+ (**47 / 15,939** cells > 0). Overlap with marker leftover epithelium = **1 / 47** (leftover frac+ 0.0037 vs rest 0.0029) — not leftover-restricted. Inferential unit = **mouse**. Honest n = **4** (2 tumor-bearing, 2 healthy). Tumor-only n = **2** — Spearman not computed on n=2.

Mouse-level mean-log1p Cldn4 vs author T fraction: n=4, ρ = +0.400 [-0.911, +0.983], p = 0.600 (descriptive). vs myeloid fraction: n=4, ρ = -0.400 [-0.983, +0.911], p = 0.600. vs T+NK: n=4, ρ = +0.400 [-0.911, +0.983], p = 0.600. **n=4 is thin.** Do not write a precise effect. Do not treat 15,939 cells as n.

## Honest n

| Item | n | Note |
|---|---:|---|
| Biological mice deposited | **4** | Healthy 1, Healthy 2, Tumor-bearing 1, Tumor-bearing 2 |
| Tumor-bearing / healthy | **2 / 2** | KP1.9 lung vs tumor-free lung |
| Libraries (technical) | 14 | not the test n |
| Cells after author QC | 15,939 | CD45+; not the test n |
| Author major types | 7 | Neutrophils 8,022; B 2,813; MoMacDC 2,397; T 1,981; NK 630; pDC 62; Basophils 34 |
| Author epithelium / malignant | **0** | CD45+ design; no Type I/II / club / ciliated / KP labels |
| Marker leftover Epcam+ or (Cdh1+ and Krt8+) | **270** | Healthy 1 39 (0.013), Healthy 2 73 (0.019), Tumor-bearing 1 113 (0.020), Tumor-bearing 2 45 (0.013) |
| Cldn4 > 0 cells | **47** | sparse in immune cells |
| Mice with ≥20 T and ≥20 myeloid | **4** | all four pass a composition gate |
| Dual-high Tacstd2∩Cldn4 | **not defined** | Cldn4-only |
| Human GSE127465 cells used | **0** | not mega-merged |

Tumor-bearing mice (author cells): 9201. Healthy: 6738.

## Gate

| File | Bytes | Used |
|---|---:|---|
| `GSE127465_mouse_counts_normalized_15939x28205.mtx.gz` | 174,301,312 | **yes** — author-normalized; 15,939 × 28,205 |
| `GSE127465_mouse_cell_metadata_15939x12.tsv.gz` | 453,068 | **yes** |
| `GSE127465_gene_names_mouse_28205.tsv.gz` | 73,869 | **yes** — Cldn4 present |
| Human MTX / human metadata / human gene names | — | **no** |
| `GSE127465_RAW.tar` / SRA | — | **no** |

Stop-if-missing does **not** apply: processed mouse MTX is public and Cldn4 is present. Thin epithelium is not a no-go on this assignment.

## Per-mouse Cldn4 vs T and myeloid

Author-normalized MTX; Cldn4 score = mean log1p of that value in **all CD45+ cells of that mouse** (no epithelial subset exists). T fraction = author `T cells` / cells. T/NK = (`T cells` + `NK cells`) / cells. Myeloid = (`Neutrophils` + `MoMacDC` + `pDC` + `Basophils`) / cells.

| mouse | condition | n cells | n T | n T+NK | n myeloid | leftover epi | Cldn4+ | Cldn4 frac+ | Cldn4 mean log1p | frac T | frac T+NK | frac myeloid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Healthy 1 | healthy | 2950 | 239 | 396 | 2228 | 39 | 11 | 0.0037 | 0.0023 | 0.081 | 0.134 | 0.755 |
| Healthy 2 | healthy | 3788 | 251 | 324 | 2851 | 73 | 5 | 0.0013 | 0.0007 | 0.066 | 0.086 | 0.753 |
| Tumor-bearing 1 | tumor | 5741 | 933 | 1196 | 3712 | 113 | 14 | 0.0024 | 0.0015 | 0.163 | 0.208 | 0.647 |
| Tumor-bearing 2 | tumor | 3460 | 558 | 695 | 1724 | 45 | 17 | 0.0049 | 0.0030 | 0.161 | 0.201 | 0.498 |

Machine table: `results/tables/per_mouse.tsv`.

## Spearman (mouse is the unit)

**n=4 is thin.** Fisher-z CI on n=4 is almost the full [−1, 1] range. Tumor-only n=2 is reported as two points, not a ρ. Full sensitivity table (neutrophil / MoMacDC / tumor-only NA rows) is in `results/tables/spearman.tsv`.

| Contrast | n | ρ [95% CI] | p | Note |
|---|---:|---|---:|---|
| Cldn4_mean_log1p vs frac_T | 4 | +0.400 [-0.911, +0.983] | 0.600 | n=4 is thin; p descriptive |
| Cldn4_mean_log1p vs frac_TNK | 4 | +0.400 [-0.911, +0.983] | 0.600 | n=4 is thin; p descriptive |
| Cldn4_mean_log1p vs frac_myeloid | 4 | -0.400 [-0.983, +0.911] | 0.600 | n=4 is thin; p descriptive |
| Cldn4_frac_pos vs frac_T | 4 | +0.400 [-0.911, +0.983] | 0.600 | n=4 is thin; p descriptive |
| Cldn4_frac_pos vs frac_myeloid | 4 | -0.400 [-0.983, +0.911] | 0.600 | n=4 is thin; p descriptive |

Tumor-bearing only (n=2, no Spearman):
- Tumor-bearing 1: Cldn4 mean log1p = 0.0015, frac T = 0.163, frac T+NK = 0.208, frac myeloid = 0.647, Cldn4+ = 14/5741.
- Tumor-bearing 2: Cldn4 mean log1p = 0.0030, frac T = 0.161, frac T+NK = 0.201, frac myeloid = 0.498, Cldn4+ = 17/3460.

Healthy only (n=2, no Spearman):
- Healthy 1: Cldn4 mean log1p = 0.0023, frac T = 0.081, frac T+NK = 0.134, frac myeloid = 0.755, Cldn4+ = 11/2950.
- Healthy 2: Cldn4 mean log1p = 0.0007, frac T = 0.066, frac T+NK = 0.086, frac myeloid = 0.753, Cldn4+ = 5/3788.

## Cldn4 by author major type

No author epithelial row. Detection is in immune lineages (sparse).

| condition | major type | n cells | Cldn4+ | frac+ |
|---|---|---:|---:|---:|
| healthy | Neutrophils | 4429 | 9 | 0.0020 |
| healthy | B cells | 939 | 0 | 0.0000 |
| healthy | MoMacDC | 629 | 3 | 0.0048 |
| healthy | T cells | 490 | 3 | 0.0061 |
| healthy | NK cells | 230 | 1 | 0.0043 |
| healthy | Basophils | 11 | 0 | 0.0000 |
| healthy | pDC | 10 | 0 | 0.0000 |
| tumor | Neutrophils | 3593 | 9 | 0.0025 |
| tumor | B cells | 1874 | 6 | 0.0032 |
| tumor | MoMacDC | 1768 | 12 | 0.0068 |
| tumor | T cells | 1491 | 3 | 0.0020 |
| tumor | NK cells | 400 | 1 | 0.0025 |
| tumor | pDC | 52 | 0 | 0.0000 |
| tumor | Basophils | 23 | 0 | 0.0000 |

## What this can and cannot say

**Can say.** Public mouse processed MTX exists. Cldn4 is on the matrix. The deposit is CD45+ / myeloid-rich lung (author Neutrophils are the modal type). Mouse-level Cldn4 (all CD45+ cells) vs T and myeloid can be written with **n=4**. Epithelium is thin (author 0; marker leftover 270).

**Cannot say.** This is not an epithelial/malignant Cldn4 score. It is not a Type II / KP tumor-cell state. It is not n=15,939. It is not a precise ρ. It is not a human–mouse mega-merge. Tacstd2 was audited only and was not a gate. Healthy vs tumor is 2 vs 2.

## Methods

- Source: GEO series supplementary mouse MTX + metadata + gene names only.
- Matrix values are the author total-count-normalized inDrops counts (filename `*_counts_normalized_*`). Cldn4 score = mean log1p of that value.
- Lineages = author `Major cell type`. Myeloid = Neutrophils + MoMacDC + pDC + Basophils.
- Leftover epithelium = Epcam>0 or (Cdh1>0 and Krt8>0) on the same MTX. Not used as the Cldn4 denominator because it is not an author tumor compartment. Cldn4+ ∩ leftover = 1 of 47.
- Spearman on 4 mice. Fisher-z 95% CI. p-values descriptive.
- Human GSE127465 files were not downloaded for this folder.

## What this is not

- Not human GSE127465 (7 patients; 40,362 tumor cells). That analysis is a separate folder.
- Not a cell-level merge of mouse + human MTX.
- Not ICI / MPR / RECIST.
- Not dual-high Tacstd2×Cldn4.
- Not a Slingshot / PAGA / CellChat redo.
- Not a claim that Cldn4 is a myeloid marker. Sparse CD45+ detection is compatible with ambient epithelial RNA.

## Outputs

- `results/tables/per_mouse.tsv`
- `results/tables/by_majortype.tsv`
- `results/tables/spearman.tsv`
- `results/tables/geo_file_inventory.tsv`
- `results/figures/fig_mouse_cldn4_vs_T_myeloid.png`
- `results/figures/fig_honest_n.png`
- `results/figures/fig_cldn4_by_majortype.png`
- `results/summary.json`

## 结论

GSE127465 **小鼠** processed MTX 公开且 <2 GB，基因表有 **Cldn4**，因此是 GO，不是 no-go。对象是 CD45+ / 髓系为主的肺（作者上皮 **0**；Epcam/Cdh1∩Krt8 残留 **270**）。Cldn4 稀疏（**47/15,939**，与残留上皮重叠 **1**）。小鼠层面 Cldn4 vs T 分数 n=4 ρ=+0.40（p=0.60），vs 髓系 ρ=−0.40（p=0.60）。**n=4 太薄**，不能写成精确效应，也不能写成 n=15,939。肿瘤侧只有 2 只鼠。**没有**与人 GSE127465 合并。不是 dual-high。

## Reproduce

```bash
python3 methods/gse127465_mouse_cldn4/scripts/download.py --out /tmp/gse127465_mouse
python3 methods/gse127465_mouse_cldn4/scripts/analyze.py \
  --data /tmp/gse127465_mouse \
  --outdir methods/gse127465_mouse_cldn4/results \
  --finding methods/gse127465_mouse_cldn4/FINDING.md
```

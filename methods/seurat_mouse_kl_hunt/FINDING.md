# Leftover public mouse lung-tumor scRNA hunt (Seurat, Cldn4-only)

Additive catalog of **public** mouse lung-tumor **scRNA** with **Cldn4** plus immune or GEMM labels. Assigned accessions were cataloged, not re-scored as the only analysis. Private 8 KL mice were not used. FASTQ was not downloaded. Every scored table was produced with **R 4.3.3 + Seurat 5.0.1**.

## Assigned (catalog only; not re-scored here)

GSE179502, GSE154989, GSE267321, GSE154977, GSE165641, GSE180963, GSE127465, GSE274477.

## Catalog

`methods/seurat_mouse_kl_hunt/catalog.tsv` lists leftovers, assigned series, and no-gos (bulk-only, CD45-only, n=2, DEG-only, wrong tissue, 9.9 Gb Parse dump, private 8 KL).

## What was scored

Best unused processed leftovers (not the assigned eight):

| Accession | Why this leftover | Honest n | Folder |
|---|---|---|---|
| **GSE179501** | Lkb1-XTR **total-viable** sister of assigned GSE179502; mixed lineages so T/NK fraction is real | **4 mice** (2 Restored, 2 Non-Restored) | `GSE179501/mouse_level_scores.tsv` |
| **GSE201247** | Kras / ATTAC whole lung; largest unused easy n with processed h5 | **6 mice** (2 WT, 2 Kras, 2 ATTAC;Kras) | `GSE201247/mouse_level_scores.tsv` |
| **GSE266323** | KP LUAD TME with inflammatory labels; Cldn4 not at the floor | **4 mice** (2 Malat1 CRISPRa, 2 Tomato) | `GSE266323/mouse_level_scores.tsv` |

Unused but not scored: GSE149813 (honest n=2 mice after YFP collapse), GSE264739 (n=6 KP/KPP, 711 Mb RAW; kept for a later pass), GSE297023 (n=16 but 9.9 Gb h5ad), GSE322632 (Lkb1 hashed Seurat object 1.9 Gb / CD45-enriched), GSE194166 (KL vs KP immune-sorted so Cldn4 is the wrong compartment).

## Methods

**Cldn4 only** (no other claudins).

**Seurat:** `CreateSeuratObject` from GEO MTX/h5 → QC (200–10000 features, ≥500 UMI) → `NormalizeData` (log1p, scale 1e4) → `AddModuleScore`.

**T/NK genes (16/16 present in all three):** Cd3d, Cd3e, Cd3g, Cd2, Cd8a, Cd8b1, Cd4, Nkg7, Gzma, Gzmb, Prf1, Klrb1c, Ncr1, Klrd1, Klrc1, Ifng.

**IFN/MHC genes:** 30/30 in GSE179501 and GSE201247; 29/30 in GSE266323 (Ifnb1 absent). Stat1/2, Irf1/7/9, Isg15, Ifit1/2/3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1, B2m, H2-K1/D1/Q4/Q6/Q7, H2-Aa/Ab1/Eb1, Tap1/2, Psmb8/9, Nlrc5, Ciita.

**Cell calls (marker-positive, not a full annotation):** epithelial = Epcam>0 or (Krt8>0 and Ptprc==0); T/NK = Cd3e/Cd3d/Nkg7/Ncr1>0.

**Unit:** one biological mouse. Mouse-level Cldn4 is mean log-normalized Cldn4 in epithelial cells (require ≥10 epithelial cells). High vs low is a **median split** of that value. Tests: Spearman and two-sided Mann–Whitney. No multiple-testing theater.

Processed inputs only: GSE179501 combined MTX; GSE201247 filtered h5; GSE266323 GEX MTX (antibody capture dropped).

## Results (honest n)

Cldn4 was present in all three feature tables. No mouse-level test was significant. That is the result, not a power complaint after the fact.

### 1. GSE179501 — Lkb1-XTR total viable (n=4 mice)

20,351 cells deposited → 20,299 after QC. Epithelial Cldn4 is sparse (2.5–5.2% epithelial cells positive).

| Mouse | Cohort | n cells | n epi | frac T/NK | mean epi Cldn4 |
|---|---|---:|---:|---:|---:|
| CM2260 | Restored | 4279 | 388 | 0.355 | 0.071 |
| CM2319 | Restored | 4807 | 637 | 0.277 | 0.034 |
| CM2324 | Non-Restored | 6598 | 665 | 0.251 | 0.044 |
| CM2328 | Non-Restored | 4615 | 1208 | 0.192 | 0.051 |

| Endpoint | n | Spearman r | p | High-vs-low MWU p |
|---|---:|---:|---:|---:|
| T/NK fraction | 4 | +0.20 | 0.80 | 1.00 |
| IFN/MHC in epithelium | 4 | −0.80 | 0.20 | 0.25 |
| Restored vs Non-Restored epi Cldn4 | 2 vs 2 | — | — | 1.00 (0.053 vs 0.047) |

No Lkb1-restore Cldn4 difference at n=2 vs 2. IFN trend is the expected inverse direction and is not significant.

### 2. GSE201247 — Kras / ATTAC whole lung (n=6 mice)

70,480 cells deposited → 69,680 after QC. Epithelial Cldn4 is near the floor (0.5–1.5% epithelial cells positive; mean log-norm 0.002–0.009). This is not a usable high/low continuum.

| Endpoint | n | Spearman r | p | High-vs-low MWU p |
|---|---:|---:|---:|---:|
| T/NK fraction (all 6) | 6 | +0.49 | 0.33 | 0.19 |
| IFN/MHC in epithelium (all 6) | 6 | −0.20 | 0.70 | 0.66 |
| T/NK fraction (Kras-bearing only) | 4 | +0.80 | 0.20 | 0.25 |

Direction of Cldn4 vs T/NK is **positive**, opposite a “Cldn4-high / T-low” law, and not significant. Cldn4 abundance is too low to support a claim.

### 3. GSE266323 — KP LUAD TME (n=4 mice)

31,084 cells deposited → 29,922 after QC (GEX only). This is the only leftover with non-trivial epithelial Cldn4 (19–43% positive).

| Mouse | Group | n cells | n epi | frac T/NK | mean epi Cldn4 | % epi Cldn4+ |
|---|---|---:|---:|---:|---:|---:|
| d10_1 | Malat1 CRISPRa | 9905 | 960 | 0.231 | 0.297 | 18.9 |
| d10_2 | Malat1 CRISPRa | 4729 | 715 | 0.240 | 0.370 | 23.5 |
| dTom_1 | Tomato | 8379 | 1374 | 0.246 | 0.613 | 42.8 |
| dTom_2 | Tomato | 6909 | 868 | 0.472 | 0.299 | 21.9 |

| Endpoint | n | Spearman r | p | High-vs-low MWU p |
|---|---:|---:|---:|---:|
| T/NK fraction | 4 | +0.40 | 0.60 | 1.00 |
| IFN/MHC in epithelium | 4 | −0.40 | 0.60 | 1.00 |
| Malat1 vs Tomato epi Cldn4 | 2 vs 2 | — | — | 0.70 (0.334 vs 0.456) |

dTom_2 is a T/NK-high outlier with mid Cldn4. No Cldn4–immune law at n=4.

## Bottom line

Leftover **public mouse lung-tumor scRNA** with processed matrices, Cldn4, and immune or GEMM labels still exist after the assigned eight are set aside. The three scored here have honest mouse n of **4, 6, and 4**.

- **KL-related leftover:** GSE179501 (not GSE179502). Cldn4 is sparse; restore vs not is n=2 vs 2 and not different.
- **Largest unused easy n:** GSE201247 n=6, but Cldn4 is at the floor.
- **Best leftover Cldn4 dynamic range:** GSE266323 KP TME n=4. No significant Cldn4 vs T/NK or IFN association.

No result here is a KL Cldn4–T/NK law. The private 8 KL mice were not used.

## Re-run

```bash
# processed GEO files in /tmp/geo_dl (URLs in catalog.tsv). Do not fetch SRA FASTQ.
Rscript methods/seurat_mouse_kl_hunt/scripts/score_leftovers.R
```

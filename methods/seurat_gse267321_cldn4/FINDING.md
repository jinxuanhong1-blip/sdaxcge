# FINDING — Seurat GSE267321 LKR13 K / KK / KLK (Cldn4-only)

**ADDITIVE. Public mouse. Cldn4-only. No dual-high. Seurat primary.** Thesis is already correct and is not rewritten. This folder is the R + Seurat (`CreateSeuratObject`) run on one public GEO object: [GSE267321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267321) (Qian / Skoulidis / Heymach, MDACC). LKR13 syngeneic subcutaneous tumors in 129SV mice. Genotypes on GEO: **K** (KrasG12D), **KK** (KrasG12D + KEAP1 KO), **KLK** (KrasG12D + KEAP1/LKB1 KO). There is **no LKR13-KL** (STK11-only) library and **no private 8 KL** object is used. **KLK is the closest public stand-in for user KL** because it is the only STK11/LKB1-loss arm.

GEO title: *Gene expression profile at single cell level of **non-malignant cells** in KRAS syngeneic mouse tumor models harboring STK11 and/or KEAP1 co-mutation.* Processed file (corrected 12 Feb 2026): `GSE267321_Normalized_expression_matrix_02122026.csv.gz`. No author cluster / malignant column on GEO. No ICI arm (`treatment: no`).

**Engine.** R 4.6.1 + Seurat 5.5.1 / SeuratObject 5.4.0. `CreateSeuratObject(counts=..., project="GSE267321", min.cells=0, min.features=0)` then `NormalizeData` + `AddModuleScore` for IFN/MHC. Python is not the primary. Seurat **did run**.

**Verdict.** The public object is what the title says: a **non-malignant** digest. Marker leftover epithelium exists (**380** cells, 4.8% of the matrix) but is **not** an author-malignant compartment. Cldn4 is in the matrix and is **at the floor: 6 / 7956 cells > 0 (0.075%)**. Cldn4 cannot be tested as a genotype effect in epithelium. T/NK **can**: both KLK tumors are below both K tumors (sample-mean fraction **0.095 vs 0.371**, Δ = -0.276). Honest n = **2 vs 2 tumors**. IFN/MHC in leftover epithelium is mixed and KLK-2 has only 11 epi cells — not a claim. Host-lung epithelium (Sftpc/Scgb1a1/Ager ∩ epi markers) = **3** cells.

---

## Decision

| Question | Answer |
|---|---|
| Engine | **R + Seurat CreateSeuratObject** (not Python-only) |
| Public processed matrix | **yes** — 15.7 MB CSV.gz on GEO |
| Author malignant / epithelial labels | **no** — barcode matrix only; title says non-malignant |
| Marker tumor epithelium | **380** cells |
| Host-lung epithelium | **3** cells |
| Cldn4 row present | **yes** — **6 cells > 0** (floor) |
| Cldn4 usable as a genotype test | **no** — too few positive cells |
| T/NK by genotype | **yes** — KLK < K in both tumors |
| IFN/MHC in leftover epithelium | **scored in leftover marker epi; mixed / underpowered** |
| Closest to user KL | **KLK** (STK11/LKB1 loss). KK is KEAP1-only. No KL library. |
| Private 8 KL | **not used** |
| Unit | **tumor / replicate** (2 per genotype) |
| Dual-high TACSTD2 × Cldn4 | **not defined** |
| ICI / PD-1 | **no** |

Cldn4 was scored because the row exists. The score is a floor, not a biology test. T/NK fraction is the genotype readout this object can support.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO series | **1** | GSE267321 |
| Tumors (unit) | **6** | 2 K + 2 KK + 2 KLK |
| KLK vs K tumors | **2 vs 2** | MW p cannot beat 1/3 |
| Cells in matrix | **7956** | not the unit |
| Genes | **22497** | mm10 symbols |
| Author malignant labels | **0** | none deposited |
| Marker epithelial cells | **380** | not author-malignant |
| Host-epithelial cells | **3** | residual normal lung |
| T/NK cells | **1794** | Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c |
| Cldn4-positive cells (any) | **6** | 0.075% of matrix |
| Dual-high | **0** | not defined |
| ICI arms | **0** | untreated |
| Private KL libraries | **0** | public GEO only |

Do not write n = 7956. Do not write n = 3 genotypes as if they were biological replicates of KL. Do not import a private 8-KL object.

---

## Genotype table (unit = genotype, built from 2 tumors)

K = KrasG12D parental. KK = KEAP1-loss. KLK = STK11/LKB1 + KEAP1-loss (closest to user KL).

| genotype | STK11 | KEAP1 | n tumors | n cells | n epi | n host-epi | n T/NK | frac T/NK | Cldn4 all mean | Cldn4 all %pos | Cldn4 epi mean | Cldn4 epi %pos | IFN epi | MHC-I epi |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| K | WT | WT | 2 | 3090 | 207 | 3 | 1172 | 0.379 | 0.00453 | 0.0971 | 0.0628 | 0.966 | 0.288 | 1.115 |
| KK | WT | KO | 2 | 2853 | 102 | 0 | 411 | 0.144 | 0.000351 | 0.0351 | 0 | 0 | 0.394 | 1.179 |
| KLK | KO | KO | 2 | 2013 | 71 | 0 | 211 | 0.105 | 0.000994 | 0.0994 | 0 | 0 | 0.495 | 1.234 |

Source: `tables/genotype_table.tsv` (written by the Seurat script after `CreateSeuratObject`).

---

## Per-tumor table (the actual unit)

| sample | genotype | n cells | n epi | n host-epi | n T/NK | frac T/NK | Cldn4 all | Cldn4 %pos | Cldn4 epi | IFN epi | MHC-I epi |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LKR13-K-1 | K | 1767 | 127 | 0 | 753 | 0.426 | 0.00736 | 0.113 | 0.102 | 0.304 | 1.232 |
| LKR13-K-2 | K | 1323 | 80 | 3 | 419 | 0.317 | 0.000756 | 0.0756 | 0 | 0.264 | 0.929 |
| LKR13-KK-1 | KK | 2172 | 82 | 0 | 367 | 0.169 | 0.00046 | 0.046 | 0 | 0.423 | 1.275 |
| LKR13-KK-2 | KK | 681 | 20 | 0 | 44 | 0.0646 | 0 | 0 | 0 | 0.275 | 0.784 |
| LKR13-KLK-1 | KLK | 1232 | 60 | 0 | 171 | 0.139 | 0 | 0 | 0 | 0.542 | 1.355 |
| LKR13-KLK-2 | KLK | 781 | 11 | 0 | 40 | 0.0512 | 0.00256 | 0.256 | 0 | 0.241 | 0.571 |

---

## KLK vs K (primary)

Closest public STK11-loss arm versus parental K. Tumor-level means; n=2 vs 2. Wilcoxon / Mann-Whitney computed in R on the two tumor values.

| metric | n tumors | KLK | K | Δ (KLK−K) | MW p | KLK values | K values |
|---|---|---:|---:|---:|---:|---|---|
| frac_tnk | 2 vs 2 | 0.095 | 0.371 | -0.276 | 0.333 | 0.138799,0.0512164 | 0.426146,0.316704 |
| Cldn4_all_mean | 2 vs 2 | 0.00128 | 0.00406 | -0.00278 | 0.667 | 0,0.00256082 | 0.0073571,0.000755858 |
| Cldn4_all_pctpos | 2 vs 2 | 0.128 | 0.0944 | 0.0337 | 1.000 | 0,0.256082 | 0.113186,0.0755858 |
| Cldn4_epithelial_mean | 2 vs 2 | 0 | 0.0512 | -0.0512 | 1.000 | 0,0 | 0.102362,0 |
| Cldn4_epithelial_pctpos | 2 vs 2 | 0 | 0.787 | -0.787 | 1.000 | 0,0 | 1.5748,0 |
| ifn_isg_epithelial_mean | 2 vs 2 | 0.391 | 0.284 | 0.108 | 1.000 | 0.541537,0.241368 | 0.303881,0.263899 |
| mhc1_epithelial_mean | 2 vs 2 | 0.963 | 1.080 | -0.117 | 1.000 | 1.355,0.570821 | 1.23172,0.928504 |
| n_epithelial | 2 vs 2 | 35.500 | 103.5 | -68.000 | 0.333 | 60,11 | 127,80 |

**T/NK is the only directional genotype result.** Both KLK tumors sit below both K tumors (0.139 and 0.0512 vs 0.426 and 0.317). KK is also T/NK-low (0.169 and 0.0646) — KEAP1-loss alone already looks cold, and KLK is not warmer. Mann-Whitney at 2 vs 2 cannot beat p = 1/3; the sign is the result.

**Cldn4 is empty.** 6 positive cells in the whole matrix (`tables/cldn4_positive_cells.tsv`). That is dropout / floor, not “KLK down-regulates Cldn4.”

**IFN/MHC in leftover epithelium is not a claim.** KLK-1 leftover epi IFN is 0.542 (n=60); KLK-2 is 0.241 (n=11). Tumor-level MW p cannot support a genotype effect. Seurat `AddModuleScore` values are in `genotype_table.tsv` as `seurat_ifn_epithelial_mean` / `seurat_mhc1_epithelial_mean` and were not used to invent a leftover-epi IFN claim.

KK is the KEAP1-only extra arm. It is **not** user KL.

---

## Epithelium (honest)

GEO title and summary say **non-malignant cells**. The digest protocol is whole-tumor (collagenase / hyaluronidase / dispase; no CD45 sort on the GEO record). The public object has **no** author `Malignant` / `Epithelial` column.

This page therefore uses markers on the Seurat counts layer, not author calls:

- **Marker tumor epithelium** = Epcam+ **or** (Cdh1+ and Krt8+) **or** (Krt18+ and Krt19+), and **not** Sftpc / Scgb1a1 / Ager.
- **Host epithelium** = those same epi markers **and** Sftpc or Scgb1a1 or Ager.
- T/NK = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Epithelium wins if both fire.

Tumors are **subcutaneous** (GEO `source_name`: Subcutaneous tumor). There is no orthotopic lung, so host AT2 should be near zero. Sftpc and Sftpa1 are **missing from the matrix**; Sftpb is present and **0 / 7956**. Scgb1a1 = 5 cells. The 380 leftover Epcam/keratin cells are therefore residual LKR13 tumor epithelium the authors did not fully strip, not normal lung. KLK-2 has **11** such cells (below a ≥20 occupancy floor). Do not call this an author-malignant atlas.

Lineage genes present: Epcam, Cdh1, Krt8, Krt18, Krt19, Cldn4, Cd3d, Cd3e, Nkg7, Ptprc.
Missing: Sftpc.

---

## Methods (short)

1. Public only. Downloaded `GSE267321_Normalized_expression_matrix_02122026.csv.gz` and the series matrix from NCBI GEO FTP. No SRA / FASTQ. No private object. No private 8 KL.
2. Matrix is genes × cells, 10x-style barcodes `{K|KK|KLK}_{UMI}.1_LKR13.{geno}.{rep}`. Genotype and replicate parsed from the barcode. That is the only metadata join.
3. Ingest: `data.table::fread` → sparse `dgCMatrix` → **`CreateSeuratObject`**. Values: counts_log1p (frac_int=1.000); CreateSeuratObject treats matrix as counts. Positive = raw counts layer > 0. Simple module scores = mean log1p of present genes. Seurat `AddModuleScore` run after `NormalizeData` as a second IFN/MHC score.
4. Cldn4-only. Tacstd2 is inventory/audit, never a gate. No dual-high.
5. IFN ISG = mouse orthologs of the public type-I ISG core (Isg15, Ifit1/2/3, Mx1, Oas1a/Oas2, Stat1, Irf7, …). MHC-I APM = B2m, H2-K1/D1/Q4/Q7/T23, Tap1/2, Psmb8/9/10, Nlrc5, ….
6. Unit = tumor. Primary contrast = KLK vs K. KK is extra.
7. Thesis is not rewritten.

```bash
Rscript methods/seurat_gse267321_cldn4/scripts/analyze_gse267321_seurat.R
```

---

## How to read this

- **Additive public mouse**, not a human concordant-pool join. Additive to the existing public GSE267321 folder; this page is the Seurat engine.
- **KLK ≠ KL.** KLK is STK11-loss **plus** KEAP1-loss. It is the closest arm in *this* series, not a clean STK11-only replicate of user KL.
- **No private 8 KL.** Honest n is the six public tumors.
- **Title says non-malignant.** Leftover Epcam/keratin is residual tumor epithelium, not an author malignant call.
- **Cldn4 is present and empty (6 cells).** Score it; do not build a Cldn4-high story on this object.
- **T/NK fraction is the genotype table:** KLK (and KK) colder than K. n=2 vs 2.
- **Honest n is 2 vs 2 tumors.** Direction can be described. A p-value cannot.
- **No dual-high. No ICI endpoint. Thesis unchanged.**

## Files

- `tables/genotype_table.tsv` — required genotype roll-up (done criterion)
- `tables/per_tumor.tsv` — unit-level table
- `tables/klk_vs_k.tsv` — primary contrast
- `tables/cldn4_positive_cells.tsv` — every Cldn4>0 barcode
- `tables/compartment_by_genotype.tsv`
- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/summary.json`, `tables/sessionInfo.txt`
- `figures/fig1_tnk_fraction.png`, `fig2_cldn4_all.png`, `fig3_compartments.png`
- `scripts/analyze_gse267321_seurat.R`

## 结论

GSE267321 是公开的 LKR13 皮下同基因瘤 scRNA（K / KK / KLK；**无 KL 库**）。本页用 **R + Seurat CreateSeuratObject** 跑同一份 GEO normalized CSV，**不用私有 8 KL**。**KLK（STK11/LKB1+KEAP1 缺失）是本系列最接近用户 KL 的一臂。** GEO 标题为 **non-malignant cells**：无作者恶性标签；标记残留上皮 **380** 个细胞（皮下，不是正常肺）。Cldn4 行在，但全矩阵仅 **6** 个细胞 >0，**不能做 Cldn4 基因型检验**。T/NK 可以：KLK 两瘤都低于 K 两瘤（样本均数 0.095 vs 0.371）。IFN/MHC 在残留上皮上混合且 KLK-2 仅 11 个上皮细胞，不作结论。诚实 n = **2 vs 2 个瘤**。无 dual-high，无 ICI，不改 thesis。


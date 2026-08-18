# FINDING — Seurat/Harmony pair GSE154977 + GSE180963, Cldn4-only

**ADDITIVE. Public mouse pair. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. This folder integrates **two** public GEO objects only: [GSE154977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154977) (Marjanovic KP 30w 10x; [PMID 32707077](https://pubmed.ncbi.nlm.nih.gov/32707077/)) and [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) (Bai / Guo K vs KL whole-tumor 10x; [CAN-22-1740](https://doi.org/10.1158/0008-5472.can-22-1740)). **Not the triple.** GSE267321 was not opened. No human cohort. No private 8-KL matrix.

Primary engine is **R + Seurat 5.5.1 + Harmony 2.0.5** (`CreateSeuratObject` → `RunPCA` → `RunHarmony(group.by.vars = "dataset")`). No Python-only primary. If Seurat had failed to install, this folder would have stopped.

**Verdict.** The pair object exists. Honest n = **6 mice** (4 KP + 1 K + 1 KL), not 25,277 cells. **Cldn4 vs T/NK is a no-go** at the pair: only GSE180963 has an honest T/NK fraction (n = 2). GSE154977 is FACS `CD45−` tumor — that design is noted, not sold as T/NK, and those four mice are not given a T/NK fraction. **Cldn4 vs epithelial IFN/MHC** can be written at n = 6, but the unadjusted Spearman is a **dataset offset**, not a Cldn4 law. After a dataset residual / OLS covariate the sign is **not down** and is **not significant**. Leave-one-mouse-out keeps the same split: unadjusted ρ stays ≤ 0 while a GSE180963 mouse remains; dataset-adjusted ρ stays positive. Q4 vs Q1 is locked (needs n ≥ 8). Thesis unchanged.

---

## Decision

| Question | Answer |
|---|---|
| Pair (not triple) | **yes** — GSE154977 + GSE180963 only |
| Seurat / Harmony | **yes** — Seurat 5.5.1, Harmony 2.0.5, reduction `harmony`, 1 iteration to converge |
| Shared genes | **18,529** (154977 27,999 ∩ 180963 20,304) |
| Cldn4 row | **yes** |
| Honest unit | **mouse** (n = 6) |
| Cldn4 vs T/NK | **no-go** — honest T/NK mice = 2 (GSE180963 only) |
| Cldn4 vs epi IFN/MHC | **n = 6 written**; unadjusted not significant; dataset-adjusted not down |
| Dataset covariate | **yes** — residual Spearman + OLS `~ Cldn4 + dataset` |
| Leave-one-out | **yes** — mouse and dataset |
| Q4 vs Q1 | **locked** — n = 6 < 8 |
| FACS-epi-only as pair primary | **no** |
| GSE267321 / human / private 8 KL | **not used** |
| Dual-high TACSTD2 × Cldn4 | **not defined** |
| ICI / PD-1 | **no** |

Cldn4 was scored because the row exists in the shared universe. Six mice, two designs, two genotypes, and cisplatin on two KP libraries cannot test a law.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO series (pair, not triple) | **2** | GSE154977 + GSE180963 |
| Mice (unit) | **6** | 4 KP + 1 K + 1 KL |
| T/NK-honest mice | **2** | GSE180963 only |
| IFN/MHC mice (tight epi ≥10 cells) | **6** | all pass |
| Cells after pair QC | **25,277** | not the unit |
| Shared genes | **18,529** | intersected before merge |
| Tight epithelial cells | **10,485** | Epcam+ structural+ Ptprc−; Cldn4 not a caller |
| Marker T/NK cells | **7,453** | almost all from GSE180963 |
| Cldn4-positive cells | **5,636** | count > 0 |
| Q4 vs Q1 tails | **0** | locked |
| GSE267321 | **0** | excluded |
| Human cohorts | **0** | excluded |
| Private 8-KL mice | **0** | excluded |

Do not write n = 25,277. Do not write n = 2 datasets as if they were biological replicates.

---

## Per-mouse table (the actual unit)

Tight epithelium = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−. Cldn4 is not a caller. IFN / MHC = mean log-norm of the locked compact mouse sets (Ifnb1 absent from the shared universe; 14/15 IFN genes and 15/15 MHC genes present).

| mouse | GSM | dataset | genotype | treatment | n cells | n epi | n T/NK | frac T/NK | Cldn4 epi mean | Cldn4 epi %pos | IFN epi | MHC epi |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KP_30w_ND_m3 | GSM4685281 | GSE154977 | KP | ND | 3945 | 3668 | 0 | — | 0.323 | 37.6 | 0.039 | 0.449 |
| KP_30w_ND_m4 | GSM4685282 | GSE154977 | KP | ND | 4008 | 3778 | 0 | — | 0.642 | 61.2 | 0.053 | 0.455 |
| KP_30w_Cis72_m5 | GSM4685283 | GSE154977 | KP | Cis72 | 2065 | 1917 | 1 | — | 0.563 | 58.8 | 0.049 | 0.491 |
| KP_30w_Cis72_m6 | GSM4685284 | GSE154977 | KP | Cis72 | 999 | 920 | 0 | — | 0.865 | 66.8 | 0.079 | 0.540 |
| K | GSM5481386 | GSE180963 | KrasG12D/+ | ND | 6696 | 28 | 3634 | 0.543 | 0.027 | 3.57 | 0.124 | 0.827 |
| KL | GSM5481387 | GSE180963 | KrasG12D/+;Lkb1fl/fl | ND | 7564 | 174 | 3818 | 0.505 | 0.098 | 9.77 | 0.261 | 1.697 |

GSE154977 T/NK is **blank on purpose**. Those libraries are live `tdTomato+ / CD45− / CD11b− / TER119− / CD31−`. The one leftover T/NK call in m5 is leak, not a fraction. GSE180963 tight epithelium is **28** and **174** cells — the same thin epi as the single-dataset folder.

Pair IFN/MHC scores are **not** the Hallmark-mapped all-cell scores from `methods/seurat_gse154977_cldn4`. This folder locks the compact mouse IFN/MHC sets on **tight epithelium** so both datasets use one scale.

---

## Cldn4 vs T/NK (mouse)

| contrast | honest n | Spearman | note |
|---|---:|---|---|
| Cldn4 epi vs T/NK fraction | **2** | **locked** | GSE180963 only; n < 4 |

K: Cldn4 epi 0.027, T/NK 0.543. KL: Cldn4 epi 0.098, T/NK 0.505. Two points. Do not draw a slope. Do not borrow a T/NK fraction from FACS `CD45−` KP libraries. Pair T/NK is a **no-go**.

---

## Cldn4 vs epithelial IFN/MHC (mouse)

Spearman uses the **exact** test (n ≤ 10). Asymptotic p is recorded and is **not** quoted as the claim p. Q4 vs Q1 is locked.

### Unadjusted pair (n = 6)

| family | n | ρ | exact p | usable |
|---|---:|---:|---:|---|
| IFN | 6 | **−0.371** | 0.497 | yes |
| MHC | 6 | **−0.429** | 0.419 | yes |

The negative sign is **dataset geometry**. GSE180963 tight epi sits at low Cldn4 and high IFN/MHC; GSE154977 FACS tumor sits at high Cldn4 and low IFN/MHC on this locked set. That is two assays, not a Cldn4 gradient.

### Dataset residual (n = 6)

Residualize Cldn4 and the family score on `factor(dataset)`, then Spearman.

| family | n | ρ | exact p | usable |
|---|---:|---:|---:|---|
| IFN | 6 | **+0.771** | 0.103 | yes |
| MHC | 6 | **+0.657** | 0.175 | yes |

Sign **flips**. IFN/MHC are **not down** in Cldn4-high after the dataset offset is removed. Exact p is still not significant. Do not claim IFN-up.

### OLS covariate (descriptive)

`family ~ Cldn4_epi + dataset`. n = 6 is thin.

| outcome | β_Cldn4 | p_Cldn4 | β_GSE180963 | p_dataset | R² |
|---|---:|---:|---:|---:|---:|
| IFN epi | +0.103 | 0.508 | +0.193 | 0.114 | 0.76 |
| MHC epi | +0.352 | 0.719 | +0.967 | 0.186 | 0.69 |

The Cldn4 coefficient is positive and not significant. Most of the R² is the dataset intercept.

### Within dataset

| dataset | family | n | ρ | exact p | usable |
|---|---|---:|---:|---:|---|
| GSE154977 | IFN | 4 | +1.000 | 0.083 | yes |
| GSE154977 | MHC | 4 | +0.800 | 0.333 | yes |
| GSE180963 | IFN | 2 | — | — | **no** |
| GSE180963 | MHC | 2 | — | — | **no** |

GSE154977-alone matches the single-dataset folder: not down; n = 4 at the lock floor; exact p for ρ = 1 is 0.083, not 0.

### Leave-one-mouse-out

| left out | n | IFN ρ | MHC ρ | IFN ρ (dataset-adj) | MHC ρ (dataset-adj) |
|---|---:|---:|---:|---:|---:|
| KP_30w_Cis72_m5 | 5 | −0.50 | −0.50 | +0.80 | +0.80 |
| KP_30w_Cis72_m6 | 5 | −0.50 | −0.60 | +0.60 | +0.50 |
| KP_30w_ND_m3 | 5 | −0.50 | −0.60 | +0.60 | +0.50 |
| KP_30w_ND_m4 | 5 | −0.50 | −0.50 | +0.80 | +0.80 |
| K | 5 | 0.00 | −0.10 | +0.90 | +0.60 |
| KL | 5 | 0.00 | −0.10 | +0.90 | +0.60 |

Unadjusted IFN ρ range **−0.50 to 0**. Dataset-adjusted IFN ρ range **+0.60 to +0.90**. Leaving out either GSE180963 mouse collapses the unadjusted inverse to ~0. The negative unadjusted pair ρ is those two points.

### Leave-one-dataset-out

| left out | kept | n | IFN ρ (exact p) | MHC ρ (exact p) |
|---|---|---:|---|---|
| GSE154977 | GSE180963 | 2 | locked | locked |
| GSE180963 | GSE154977 | 4 | +1.000 (0.083) | +0.800 (0.333) |

A pair that disappears when one GEO series is dropped is not a Cldn4 law.

### ND-only sensitivity (drop Cis72; n = 4)

Unadjusted IFN/MHC ρ = −0.600 (exact p = 0.417). Dataset-residual ρ = +0.600 (exact p = 0.417). Same flip. Still not significant.

Cell-level Spearman inside epithelium is in `tables/cldn4_vs_ifn_mhc_cells_exploratory.tsv`. It is **pseudoreplication**. Do not cite it as n.

---

## Epithelium and T/NK (honest)

- **Tight epithelium (primary)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.
- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Tight epithelium wins if both fire.
- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.
- GSE154977 is FACS `CD45−` tumor. That is why T/NK is a design no-go there. The pair is **not** a FACS-epi-only analysis: GSE180963 whole-tumor cells stay in the Harmony object.
- GSE180963 Sftpc ambient is unchanged from the single-dataset folder. Tight epi is the IFN/MHC unit, not loose Epcam+.

Lineage genes present in the shared universe: Epcam, Cdh1, Krt8, Krt18, Krt19, Sftpc, Scgb1a1, Ager, Cd3d, Cd3e, Cd3g, Cd8a, Nkg7, Ncr1, Klrb1c, Cldn4, Ptprc, Stk11, Tacstd2, Cldn18, Nkx2-1.
Missing from the locked IFN set: **Ifnb1**.

---

## Methods (short)

1. Public only. GEO processed matrices. No SRA / FASTQ. No GSE267321. No human. No private 8-KL.
2. GSE154977: custom COO `GSE154977_mmLung10x_cis_dSp_rawCount.h5` (MATLAB rank-2 `i`/`j`/`v`) → `CreateSeuratObject`. Four 30w KP libraries (2 ND + 2 Cis72).
3. GSE180963: Cell Ranger v2 MTX in `GSE180963_RAW.tar` → `Read10X` → `CreateSeuratObject`. One K and one KL, mixed library, demultiplexed by label.
4. Intersect symbols (18,529). Merge. Light QC (`nFeature ≥ 200`, `percent.mt < 25`); all 25,277 depositor cells kept.
5. `NormalizeData` (LogNormalize, 1e4) → HVG 2000 → `ScaleData` → `RunPCA` (30) → **`RunHarmony(group.by.vars = "dataset")`** → neighbors / Leiden 0.4 / UMAP on Harmony dims 1–20. Harmony is for the embedding. Scores use joined RNA log-norm, not Harmony space.
6. Cldn4-only. Positive = raw count > 0. Module scores = mean log-norm of present locked IFN/MHC genes. `AddModuleScore` is audit only.
7. Unit = mouse. Spearman exact when n ≤ 10. Spearman lock n ≥ 4. Q4 vs Q1 lock n ≥ 8.
8. Dataset covariate = residual Spearman + OLS. LOO = leave one mouse; also leave one dataset.
9. Thesis is not rewritten.

```bash
bash methods/seurat_integrate_154977_180963_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_154977_180963_cldn4/scripts/download.sh /tmp/pair_154977_180963
Rscript methods/seurat_integrate_154977_180963_cldn4/scripts/analyze.R \
  --data /tmp/pair_154977_180963 \
  --out methods/seurat_integrate_154977_180963_cldn4
```

---

## How to read this

- **Additive public mouse pair**, not a human concordant-pool join, not a triple.
- **Honest n is 6 mice.** Do not quote 25,277 cells.
- **T/NK cannot be tested at the pair.** Two honest mice. Four FACS `CD45−` mice contribute none.
- **Unadjusted Cldn4–IFN/MHC ρ is negative because the two GEO objects sit in different corners.** That is a dataset offset.
- **After dataset residual / OLS, IFN/MHC are not down** and are not significant. LOO agrees.
- **GSE180963 tight epi is 28 and 174 cells.** Do not overweight those two IFN/MHC points.
- **No FACS-epi-only primary. No 267321. No human. No private 8 KL. No dual-high. No ICI. Thesis unchanged.**

## Files

- `tables/per_mouse.tsv` — unit-level table
- `tables/cldn4_vs_tnk.tsv` — no-go
- `tables/family_spearman.tsv` — unadjusted + dataset-residual
- `tables/dataset_covariate_lm.tsv`
- `tables/within_dataset_spearman.tsv`
- `tables/loo_mouse.tsv`, `tables/loo_dataset.tsv`
- `tables/honest_n.tsv`, `tables/gene_inventory.tsv`, `tables/summary.json`, `tables/sessionInfo.txt`
- `tables/cldn4_vs_ifn_mhc_cells_exploratory.tsv` — do not cite as n
- `figures/fig_honest_n.png`, `fig_cldn4_vs_ifn.png`, `fig_cldn4_vs_mhc.png`, `fig_cldn4_vs_tnk.png`, `fig_loo_*.png`, `fig_umap_*.png`
- `scripts/analyze.R` — R + Seurat + Harmony primary

## 结论

GSE154977 + GSE180963 是公开小鼠 pair（不是 triple）。用 **R + Seurat 5.5.1 + Harmony 2.0.5** 做了 Cldn4-only 整合。诚实 n = **6 只鼠**，不是 25277 个细胞。T/NK 只有 GSE180963 两只鼠能计，pair 上 **不能做检验**；GSE154977 是 FACS CD45−，不当成 T/NK。上皮 IFN/MHC 的未校正 Spearman 是负的，但是 **数据集偏移**（exact p > 0.4）；去掉 dataset 后符号翻成正、仍不显著。LOO 同一结论。未打开 GSE267321、人数据或私有 8 只 KL。无 dual-high，无 ICI，不改 thesis。

# FINDING — Seurat/Harmony integrate KP+K+KL 10x (Cldn4-only)

**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. Primary engine is **R + Seurat 5.5.1 + Harmony** (`CreateSeuratObject` on public processed matrices, then `RunHarmony` by **dataset**). No Python-only primary. Private 8-KL matrices were not opened. Human series were not opened. Locked-out GEO (GSE179502 FACS epi / no T/NK, GSE154989 plate epi, GSE267321 subq Cldn4 floor, GSE127465 CD45+, GSE50927) were not added.

**Requested trio.** GSE154977 (KP 30w 10x) + GSE180963 (K/KL) + GSE165641 (KL). Each series is kept **only** if a processed matrix exists. Dropped for no matrix: **none**. Kept: **GSE154977, GSE165641, GSE180963**.

**Honest unit = mouse**, not cell n. Dataset is a covariate. Within-genotype and leave-one-dataset-out are reported so a KP/K/KL mix is not sold as a Cldn4 effect.

---

## Decision

| Question | Answer |
|---|---|
| Public processed matrices | GSE154977 **h5 COO**; GSE180963 **10x MTX**; GSE165641 **Cell Ranger MTX** |
| CreateSeuratObject | **yes** — Seurat 5.5.1 |
| Harmony | **yes** — harmony::RunHarmony, `group.by.vars = dataset` |
| Integrated object | `objects/kpkl_10x_harmony_seurat.rds` |
| Mice (unit) | **8** |
| T/NK-usable mice | **4** (mixed digest only; GSE154977 is AT2-lineage FACS / tumor-state) |
| Epithelium-usable mice | **8** |
| Cldn4 vs T/NK (mouse) | n=4 mice, ρ=-1.000, p=0 |
| Cldn4 vs T/NK, dataset-adjusted | n=4 mice, ρ=-1.000, p=0 |
| Cldn4 vs epithelial IFN | n=8 mice, ρ=-0.048, p=0.911 |
| Cldn4 vs epithelial IFN, dataset-adjusted | n=8 mice, ρ=0.690, p=0.058 |
| Cldn4 vs epithelial MHC | n=8 mice, ρ=-0.619, p=0.102 |
| Cldn4 vs epithelial MHC, dataset-adjusted | n=8 mice, ρ=-0.095, p=0.823 |
| Dual-high TACSTD2 × Cldn4 | **not defined** |
| Private 8 KL | **not used** |
| Unit | **mouse** |

Do not write n = 31970 cells.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO series requested | **3** | GSE154977 + GSE180963 + GSE165641 |
| GEO series kept (processed matrix) | **3** | GSE154977,GSE165641,GSE180963 |
| GEO series dropped (no matrix) | **0** | none |
| mice (unit) | **8** | honest unit; do not write cell n |
| KP mice (GSE154977) | **4** | AT2-lineage / tumor-state FACS; T/NK design no-go |
| K mice (GSE180963) | **1** | mixed digest |
| KL mice (GSE180963 + GSE165641) | **3** | mixed digest |
| T/NK-usable mice | **4** | mixed digest AND T/NK floor |
| epithelium-usable mice | **8** | n_epi >= 10 |
| cells after QC (not the unit) | **31970** | nFeature>=200, nCount>=500, mt<25; before=32457 |
| genes (intersected symbols) | **18529** | symbol intersect across kept series |
| Cldn4-positive cells | **5933** | count > 0 |
| dual-high TACSTD2 x Cldn4 | **0** | not defined |
| private 8-KL mice used | **0** | public GEO only |
| human series | **0** | none |
| excluded GEO (locked) | **6** | GSE179502, GSE154989, GSE267321, GSE127465, GSE50927, human |

---

## Per-mouse table (the actual unit)

| mouse | dataset | GSM | genotype | treatment | digest | n cells | n epi | n T/NK | frac T/NK | Cldn4 epi | Cldn4 epi %pos | IFN epi | MHC epi | T/NK usable | epi usable |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| GSE154977_KP_30w_Cis72_m5 | GSE154977 | GSM4685283 | KP | Cis72 | AT2_lineage_FACS | 2065 | 2057 |    1 | 0.000 | 0.5415 | 56.10 | 0.048 | 0.483 | FALSE | TRUE |
| GSE154977_KP_30w_Cis72_m6 | GSE154977 | GSM4685284 | KP | Cis72 | AT2_lineage_FACS |  999 |  995 |    0 | 0.000 | 0.8315 | 63.42 | 0.078 | 0.537 | FALSE | TRUE |
| GSE154977_KP_30w_ND_m3 | GSE154977 | GSM4685281 | KP | ND | AT2_lineage_FACS | 3945 | 3897 |    0 | 0.000 | 0.3165 | 36.05 | 0.039 | 0.437 | FALSE | TRUE |
| GSE154977_KP_30w_ND_m4 | GSE154977 | GSM4685282 | KP | ND | AT2_lineage_FACS | 4008 | 3980 |    0 | 0.000 | 0.6322 | 59.80 | 0.053 | 0.447 | FALSE | TRUE |
| GSE165641_KL1 | GSE165641 | GSM5047302 | KL | untreated | mixed | 4199 |  207 |   97 | 0.023 | 1.0900 | 64.73 | 0.148 | 0.380 | TRUE | TRUE |
| GSE165641_KL2 | GSE165641 | GSM5047303 | KL | untreated | mixed | 2494 |  160 |  253 | 0.101 | 0.5683 | 51.25 | 0.082 | 0.616 | TRUE | TRUE |
| GSE180963_K | GSE180963 | GSM5481386 | K | untreated | mixed | 6696 |   28 | 3634 | 0.543 | 0.0272 | 3.57 | 0.124 | 0.827 | TRUE | TRUE |
| GSE180963_KL | GSE180963 | GSM5481387 | KL | untreated | mixed | 7564 |  174 | 3818 | 0.505 | 0.0978 | 9.77 | 0.261 | 1.697 | TRUE | TRUE |

T/NK fraction is a real fraction only in **mixed digest** (GSE180963, GSE165641). GSE154977 GEO source is *Lung tissue, AT2 cells* with a FACS step (Marjanovic / Cancer Cell 2020 cisplatin 10x). Those KP mice are epithelium for IFN/MHC and a **design no-go** for T/NK fraction. Cis72 (m5, m6) is a 72 h cisplatin arm — kept in the object, dropped in the `EPI_no_Cis72` sensitivity.

---

## Cldn4 vs T/NK (mouse unit)

Primary: epithelial Cldn4 mean vs T/NK fraction on **tnk_usable** mice. Spearman requires n≥4. Q4 vs Q1 requires n≥8.

**Pooled mixed-digest:** n=4 mice, ρ=-1.000, p=0.

**Dataset as covariate** (Spearman of residuals after `~ dataset`): n=4 mice, ρ=-1.000, p=0.

ρ=-1 at n=4 is a **rank identity**, not a Cldn4 law. The four mixed-digest points are two libraries: GSE180963 T/NK ≈ 0.50–0.54 with low epithelial Cldn4 (K 0.027, KL 0.098; K epithelium is only 28 cells) versus GSE165641 T/NK 0.023–0.101 with higher epithelial Cldn4 (KL2 0.568, KL1 1.090). Dataset-adjusted ρ stays −1 because both within-dataset pairs go the same way (higher Cldn4, lower T/NK), but that is still **1 K vs 3 KL** and GSE180963 is 1 vs 1. Within-KL T/NK is n=3 and is locked off. Asymptotic p=0 at n=4 is not a test of a law.

Pooled K+KL is a genotype mix. Do not sell it as a Cldn4 law. Within-genotype and leave-one-dataset-out:

| subset | x | y | n mice | Spearman | dataset-adjusted | Q4 vs Q1 | note |
|---|---|---|---:|---|---|---|---|
| TNK_all_usable | Cldn4_epi_mean | frac_tnk | 4 | ρ=-1.000 p=   0 | ρ=-1.000 p=   0 | locked off (need n>=8) | mixed-digest mice only; dataset as covariate |
| TNK_all_usable | Cldn4_epi_mean | IFN_epi_mean | 4 | ρ=0.000 p=   1 | ρ=0.600 p= 0.4 | locked off (need n>=8) | mixed-digest mice only; dataset as covariate |
| TNK_all_usable | Cldn4_epi_mean | MHC_epi_mean | 4 | ρ=-0.800 p= 0.2 | ρ=0.000 p=   1 | locked off (need n>=8) | mixed-digest mice only; dataset as covariate |
| TNK_within_K | Cldn4_epi_mean | frac_tnk | 1 | locked off (n=1) | locked off | locked off (need n>=8) | within-genotype K only |
| TNK_within_K | Cldn4_epi_mean | IFN_epi_mean | 1 | locked off (n=1) | locked off | locked off (need n>=8) | within-genotype K only |
| TNK_within_K | Cldn4_epi_mean | MHC_epi_mean | 1 | locked off (n=1) | locked off | locked off (need n>=8) | within-genotype K only |
| TNK_within_KL | Cldn4_epi_mean | frac_tnk | 3 | locked off (n=3) | locked off | locked off (need n>=8) | within-genotype KL only |
| TNK_within_KL | Cldn4_epi_mean | IFN_epi_mean | 3 | locked off (n=3) | locked off | locked off (need n>=8) | within-genotype KL only |
| TNK_within_KL | Cldn4_epi_mean | MHC_epi_mean | 3 | locked off (n=3) | locked off | locked off (need n>=8) | within-genotype KL only |
| TNK_within_KP | Cldn4_epi_mean | frac_tnk | 0 | locked off (n=0) | locked off | locked off (need n>=8) | within-genotype KP only empty |
| TNK_LODO_drop_GSE154977 | Cldn4_epi_mean | frac_tnk | 4 | ρ=-1.000 p=   0 | ρ=-1.000 p=   0 | locked off (need n>=8) | leave-one-dataset-out: drop GSE154977 |
| TNK_LODO_drop_GSE154977 | Cldn4_epi_mean | IFN_epi_mean | 4 | ρ=0.000 p=   1 | ρ=0.600 p= 0.4 | locked off (need n>=8) | leave-one-dataset-out: drop GSE154977 |
| TNK_LODO_drop_GSE154977 | Cldn4_epi_mean | MHC_epi_mean | 4 | ρ=-0.800 p= 0.2 | ρ=0.000 p=   1 | locked off (need n>=8) | leave-one-dataset-out: drop GSE154977 |
| TNK_LODO_drop_GSE165641 | Cldn4_epi_mean | frac_tnk | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE165641 |
| TNK_LODO_drop_GSE165641 | Cldn4_epi_mean | IFN_epi_mean | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE165641 |
| TNK_LODO_drop_GSE165641 | Cldn4_epi_mean | MHC_epi_mean | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE165641 |
| TNK_LODO_drop_GSE180963 | Cldn4_epi_mean | frac_tnk | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE180963 |
| TNK_LODO_drop_GSE180963 | Cldn4_epi_mean | IFN_epi_mean | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE180963 |
| TNK_LODO_drop_GSE180963 | Cldn4_epi_mean | MHC_epi_mean | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE180963 |

---

## Cldn4 vs epithelial IFN / MHC (mouse unit)

Primary: epithelial Cldn4 mean vs mean lognorm of the locked IFN or MHC sets, mice with n_epi ≥ 10.

**Pooled epithelium:** IFN n=8 mice, ρ=-0.048, p=0.911; MHC n=8 mice, ρ=-0.619, p=0.102.

**Dataset as covariate:** IFN n=8 mice, ρ=0.690, p=0.058; MHC n=8 mice, ρ=-0.095, p=0.823.

Unadjusted IFN is null; dataset-adjusted IFN flips toward positive (p=0.058). That flip is why dataset stays a covariate. Within-KP IFN ρ=1 (n=4) tracks the **Cis72 vs ND** rank (Cis72 mice are higher Cldn4 and higher IFN) — do not sell treatment as a Cldn4 effect. LODO drop GSE180963 IFN ρ=0.829 (n=6) is still KP+KL. Within-KL IFN is n=3 and locked off. Q4 vs Q1 on n=8 is written and is not significant.

Pooled KP+K+KL is a genotype mix and a dataset mix. Within-genotype and LODO:

| subset | x | y | n mice | Spearman | dataset-adjusted | Q4 vs Q1 | note |
|---|---|---|---:|---|---|---|---|
| EPI_all_usable | Cldn4_epi_mean | frac_tnk | 4 | ρ=-1.000 p=   0 | ρ=-1.000 p=   0 | locked off (need n>=8) | mice with epithelium; dataset as covariate |
| EPI_all_usable | Cldn4_epi_mean | IFN_epi_mean | 8 | ρ=-0.048 p=0.911 | ρ=0.690 p=0.058 | Δmed=-0.080 p=0.699 | mice with epithelium; dataset as covariate |
| EPI_all_usable | Cldn4_epi_mean | MHC_epi_mean | 8 | ρ=-0.619 p=0.102 | ρ=-0.095 p=0.823 | Δmed=-0.803 p=0.245 | mice with epithelium; dataset as covariate |
| EPI_no_Cis72 | Cldn4_epi_mean | frac_tnk | 4 | ρ=-1.000 p=   0 | ρ=-1.000 p=   0 | locked off (need n>=8) | drop GSE154977 Cis72; dataset as covariate |
| EPI_no_Cis72 | Cldn4_epi_mean | IFN_epi_mean | 6 | ρ=-0.143 p=0.787 | ρ=0.657 p=0.156 | locked off (need n>=8) | drop GSE154977 Cis72; dataset as covariate |
| EPI_no_Cis72 | Cldn4_epi_mean | MHC_epi_mean | 6 | ρ=-0.771 p=0.0724 | ρ=-0.200 p=0.704 | locked off (need n>=8) | drop GSE154977 Cis72; dataset as covariate |
| EPI_within_K | Cldn4_epi_mean | frac_tnk | 1 | locked off (n=1) | locked off | locked off (need n>=8) | within-genotype K only |
| EPI_within_K | Cldn4_epi_mean | IFN_epi_mean | 1 | locked off (n=1) | locked off | locked off (need n>=8) | within-genotype K only |
| EPI_within_K | Cldn4_epi_mean | MHC_epi_mean | 1 | locked off (n=1) | locked off | locked off (need n>=8) | within-genotype K only |
| EPI_within_KL | Cldn4_epi_mean | frac_tnk | 3 | locked off (n=3) | locked off | locked off (need n>=8) | within-genotype KL only |
| EPI_within_KL | Cldn4_epi_mean | IFN_epi_mean | 3 | locked off (n=3) | locked off | locked off (need n>=8) | within-genotype KL only |
| EPI_within_KL | Cldn4_epi_mean | MHC_epi_mean | 3 | locked off (n=3) | locked off | locked off (need n>=8) | within-genotype KL only |
| EPI_within_KP | Cldn4_epi_mean | IFN_epi_mean | 4 | ρ=1.000 p=   0 | ρ=1.000 p=   0 | locked off (need n>=8) | within-genotype KP only |
| EPI_within_KP | Cldn4_epi_mean | MHC_epi_mean | 4 | ρ=0.800 p= 0.2 | ρ=0.800 p= 0.2 | locked off (need n>=8) | within-genotype KP only |
| EPI_LODO_drop_GSE154977 | Cldn4_epi_mean | frac_tnk | 4 | ρ=-1.000 p=   0 | ρ=-1.000 p=   0 | locked off (need n>=8) | leave-one-dataset-out: drop GSE154977 |
| EPI_LODO_drop_GSE154977 | Cldn4_epi_mean | IFN_epi_mean | 4 | ρ=0.000 p=   1 | ρ=0.600 p= 0.4 | locked off (need n>=8) | leave-one-dataset-out: drop GSE154977 |
| EPI_LODO_drop_GSE154977 | Cldn4_epi_mean | MHC_epi_mean | 4 | ρ=-0.800 p= 0.2 | ρ=0.000 p=   1 | locked off (need n>=8) | leave-one-dataset-out: drop GSE154977 |
| EPI_LODO_drop_GSE165641 | Cldn4_epi_mean | frac_tnk | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE165641 |
| EPI_LODO_drop_GSE165641 | Cldn4_epi_mean | IFN_epi_mean | 6 | ρ=-0.371 p=0.468 | ρ=0.657 p=0.156 | locked off (need n>=8) | leave-one-dataset-out: drop GSE165641 |
| EPI_LODO_drop_GSE165641 | Cldn4_epi_mean | MHC_epi_mean | 6 | ρ=-0.429 p=0.397 | ρ=0.486 p=0.329 | locked off (need n>=8) | leave-one-dataset-out: drop GSE165641 |
| EPI_LODO_drop_GSE180963 | Cldn4_epi_mean | frac_tnk | 2 | locked off (n=2) | locked off | locked off (need n>=8) | leave-one-dataset-out: drop GSE180963 |
| EPI_LODO_drop_GSE180963 | Cldn4_epi_mean | IFN_epi_mean | 6 | ρ=0.829 p=0.0416 | ρ=0.943 p=0.0048 | locked off (need n>=8) | leave-one-dataset-out: drop GSE180963 |
| EPI_LODO_drop_GSE180963 | Cldn4_epi_mean | MHC_epi_mean | 6 | ρ=-0.143 p=0.787 | ρ=-0.257 p=0.623 | locked off (need n>=8) | leave-one-dataset-out: drop GSE180963 |

---

## Epithelium and T/NK (honest)

- **Tight epithelium (mixed digest)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.
- **GSE154977 epithelium** = Ptprc− and (Epcam+ or structural+); AT2-lineage / tumor-state FACS, not a mixed digest.
- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Epithelium wins if both fire.
- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.
- **T/NK usable** = mixed digest AND n_T/NK ≥ 20 AND frac ≥ 0.02.

Lineage genes present: Epcam, Cdh1, Krt8, Krt18, Krt19, Sftpc, Scgb1a1, Ager, Cd3d, Cd3e, Cd3g, Cd8a, Nkg7, Ncr1, Klrb1c, Cldn4, Ptprc, Stk11, Tacstd2, Cldn18, Nkx2-1.
Missing: none of the core set.

---

## Methods (short)

1. Public only. Processed matrices from NCBI GEO FTP. No SRA / FASTQ. No private 8-KL object.
2. GSE154977: author COO `rawCount.h5` + smp/gene tables. GSE180963: `Read10X` on K/ and KL/ MTX. GSE165641: `Read10X` on KL1/KL2 `filtered_feature_bc_matrix`. A series with no matrix is dropped.
3. Gene universe = **symbol intersect**. Cldn4 must be present.
4. `CreateSeuratObject` → light QC (nFeature≥200, nCount≥500, percent.mt<25) → LogNormalize 1e4 → VST 2000 → ScaleData on HVG → PCA 30.
5. **Harmony by dataset** (not by mouse). Neighbors / clusters / UMAP on Harmony dims 1:20.
6. Cldn4-only. Positive = raw count > 0. Module scores = mean lognorm of present locked genes.
7. Unit = mouse. Spearman n≥4. Q4 vs Q1 n≥8. Dataset-adjusted = Spearman of `lm(~ dataset)` residuals. Within-genotype and leave-one-dataset-out are written even when locked off.
8. Thesis is not rewritten.

```bash
bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/download.sh /tmp/kpkl_10x
Rscript methods/seurat_integrate_kpkl_10x_cldn4/scripts/analyze.R --data /tmp/kpkl_10x --out methods/seurat_integrate_kpkl_10x_cldn4
```

---

## How to read this

- **Additive public mouse**, not a human concordant-pool join.
- **GSE154977 is KP epithelium**, not a T/NK fraction series.
- **GSE180963 is 1 K + 1 KL.** **GSE165641 is 2 KL.** Together they give mixed-digest T/NK mice.
- **A pooled Spearman that mixes genotypes is not a Cldn4 effect.** Read the within-genotype and LODO rows.
- **Honest n is mice.** Cell counts are inventory.
- **No dual-high. No private 8 KL. No locked-out series. Thesis unchanged.**

## Files

- `objects/kpkl_10x_harmony_seurat.rds` — integrated Seurat (gitignored if huge; exists after the run)
- `tables/mouse_table.tsv` — unit-level table
- `tables/mouse_contrasts.tsv` — pooled / dataset-adjusted / within-genotype / LODO
- `tables/honest_n.tsv`, `tables/dataset_inventory.tsv`, `tables/gene_inventory.tsv`, `tables/summary.json`
- `figures/DimPlot_dataset.png` (and genotype / compartment / mouse)
- `figures/fig_cldn4_vs_tnk.png`, `fig_cldn4_vs_ifn.png`, `fig_cldn4_vs_mhc.png`
- `scripts/analyze.R` — R + Seurat / Harmony primary

## 结论

用 **R + Seurat / Harmony** 整合了公开小鼠肺 GEMM 10x：GSE154977（KP）、GSE180963（K/KL）、GSE165641（KL），缺矩阵的系列已丢。诚实单位是 **鼠**，不是细胞数。数据集作为协变量；同时报告基因型内和 leave-one-dataset-out，避免把基因型混合物写成 Cldn4 效应。Cldn4-only，无 dual-high，无私有 8 只 KL，不改 thesis。


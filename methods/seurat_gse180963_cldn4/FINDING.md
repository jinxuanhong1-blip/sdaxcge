# FINDING — GSE180963 K vs KL lung GEMM scRNA (Cldn4-only, Seurat)

**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. This folder scores one public GEO object: [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) (Bai / Guo / Zhang / Long / Dong, Southern Medical University; paper: [CAN-22-1740](https://doi.org/10.1158/0008-5472.can-22-1740)). KrasG12D/+ (**K**, Tomato lenti) vs KrasG12D/+;Lkb1-targeted (**KL**) FVB lung tumor nodules. Primary engine is **R + Seurat 5.5.1** (`Read10X` → `CreateSeuratObject`). No Python-only primary. Private 8-KL matrices were not opened.

**Verdict.** Processed 10x MTX **exists** on GEO (`GSE180963_RAW.tar` → `K/` and `KL/` `matrix.mtx` + `genes.tsv` + `barcodes.tsv`). This is **not** a no-go for matrix. It **is** a no-go for a mouse-level Cldn4 vs T/NK or Cldn4 vs epithelial IFN/MHC test: honest n = **2 mice (1 vs 1)**. Genotype = mouse = sample. The two samples were mixed in **one 10x library** and demultiplexed by label. Cell-level p-values are pseudoreplication. Cldn4 is in the matrix and is **epithelial-restricted and sparse**. Do not write n = 14260 cells.

---

## Decision

| Question | Answer |
|---|---|
| Public processed matrix | **yes** — 10x MTX in `GSE180963_RAW.tar` (87.9 MB) |
| CreateSeuratObject | **yes** — Seurat 5.5.1 |
| Author malignant / epithelial labels | **no** — barcode matrix only |
| Marker epithelium (tight) | **202** cells (K 28 / KL 174) |
| Loose Epcam+ (audit) | **2433** cells; Sftpc is ambient |
| Cldn4 row present | **yes** — **26** cells > 0 |
| Cldn4 vs T/NK at mouse unit | **no-go** — n = 1 vs 1 mice |
| Epithelial IFN/MHC at mouse unit | **no-go** — n = 1 vs 1 mice |
| T/NK fraction (descriptive) | K 0.543 vs KL 0.505 |
| Closest to user KL | **this KL arm** (Lkb1-targeted KrasG12D/+ lung GEMM) |
| Unit | **mouse** |
| Dual-high TACSTD2 × Cldn4 | **not defined** |
| Private 8 KL | **not used** |
| ICI / PD-1 | **no** |

Cldn4 was scored because the row exists. Two mice cannot test a law. Thesis unchanged.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO series | **1** | GSE180963 |
| Mice (unit) | **2** | 1 K + 1 KL |
| K vs KL mice | **1 vs 1** | MW / Spearman cannot be computed |
| 10x libraries | **1** | mixed, then demultiplexed by label |
| Cells in matrix | **14260** | not the unit |
| Genes | **20304** | mm10 symbols |
| Author malignant labels | **0** | none deposited |
| Marker epithelial cells (tight) | **202** | Epcam+ structural+ Ptprc- |
| Loose Epcam+ (audit) | **2433** | ambient-polluted |
| Sftpc+ cells (ambient flag) | **11001** | not host AT2 |
| T/NK cells | **7452** | Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c |
| Cldn4-positive cells (any) | **26** | count > 0 |
| Dual-high | **0** | not defined |
| Private 8-KL mice | **0** | public GEO only |
| ICI arms | **0** | untreated |

Do not write n = 14260. Do not write n = 2 genotypes as if they were biological replicates.

---

## Per-mouse table (the actual unit)

| mouse | GSM | genotype | n cells | n epi tight | n epi loose | n T/NK | frac T/NK | %Sftpc | Cldn4 all | Cldn4 %pos | Cldn4 epi | Cldn4 epi %pos | IFN epi | MHC epi |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| K | GSM5481386 | KrasG12D/+ | 6696 | 28 | 1313 | 3634 | 0.543 | 76.5 | 0.0001 | 0.01 | 0.0270 | 3.57 | 0.123 | 0.822 |
| KL | GSM5481387 | KrasG12D/+;Lkb1fl/fl | 7564 | 174 | 1120 | 3818 | 0.505 | 77.7 | 0.0031 | 0.33 | 0.0970 | 9.77 | 0.259 | 1.687 |

---

## KL vs K (descriptive; n = 1 vs 1)

| metric | n mice | KL | K | Δ (KL−K) | mouse-level test |
|---|---|---:|---:|---:|---|
| T/NK fraction | 1 vs 1 | 0.505 | 0.543 | -0.038 | **none** |
| Cldn4 mean, all cells | 1 vs 1 | 0.0031 | 0.0001 | 0.0030 | **none** |
| Cldn4 %pos, all cells | 1 vs 1 | 0.33 | 0.01 | 0.32 | **none** |
| Cldn4 mean, marker epithelium | 1 vs 1 | 0.0970 | 0.0270 | 0.0701 | **none** |
| Cldn4 %pos, marker epithelium | 1 vs 1 | 9.77 | 3.57 | 6.20 | **none** |
| IFN ISG mean, marker epithelium | 1 vs 1 | 0.259 | 0.123 | 0.136 | **none** |
| MHC/APM mean, marker epithelium | 1 vs 1 | 1.687 | 0.822 | 0.865 | **none** |
| n marker epithelial cells | 1 vs 1 | 174 | 28 | 146 | **none** |

**Cldn4 vs T/NK.** Two points. Spearman at the mouse unit is not defined. K: Cldn4 epi mean 0.0270, T/NK fraction 0.543. KL: Cldn4 epi mean 0.0970, T/NK fraction 0.505. Do not draw a slope.

**Cldn4 vs epithelial IFN/MHC.** Same two mice. Tight epithelium is **28** (K) and **174** (KL) cells. Cell-level Spearman inside KL tight epithelium (n=174 cells, ρ=0.160, p=0.0353 (pseudoreplication) for IFN; n=174 cells, ρ=-0.076, p=0.319 (pseudoreplication) for MHC) is exploratory and must not be cited as n.

The series summary frames LKB1 loss as an immune-desert TME. This object does **not** independently confirm that as a drop in T/NK fraction (K 0.543 vs KL 0.505). One library, one mouse per arm.

---

## Epithelium (honest)

GEO deposits **no** author `Malignant` / `Epithelial` column. Calling is markers:

- **Tight epithelium (primary)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.
- **Loose Epcam+** is audit only. Sftpc is an **ambient flag** in this digest, not host AT2.
- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Tight epithelium wins if both fire.
- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.

Do not treat K vs KL epithelial IFN/MHC as a powered contrast. n = 1 vs 1 mice.

Lineage genes present: Epcam, Cdh1, Krt8, Krt18, Krt19, Sftpc, Scgb1a1, Ager, Cd3d, Cd3e, Cd3g, Cd8a, Nkg7, Ncr1, Klrb1c, Cldn4, Ptprc, Stk11, Tacstd2, Cldn18, Nkx2-1.
Missing: none of the core set.

---

## Methods (short)

1. Public only. Downloaded `GSE180963_RAW.tar` and the series matrix from NCBI GEO FTP. No SRA / FASTQ. No private 8-KL object.
2. Each sample is a Cell Ranger v2-style MTX (`genes.tsv` + `barcodes.tsv` + `matrix.mtx`). Integer counts. Author QC already applied (Seurat 3.1.5: 500–6000 features, mito < 20%, genes in ≥3 cells).
3. Primary: `Seurat::Read10X` → `CreateSeuratObject` → `NormalizeData` (LogNormalize, 1e4). Light PCA / neighbors / Leiden / UMAP for figures only. **No integration** — genotype is the mouse.
4. Cldn4-only. Tacstd2 is inventory, never a gate. No dual-high.
5. Positive = raw count > 0. Module scores = mean lognorm of present genes in the locked IFN and MHC sets.
6. Unit = mouse. n = 1 vs 1. No mouse-level p-value.
7. Thesis is not rewritten.

```bash
bash methods/seurat_gse180963_cldn4/scripts/install_r.sh
bash methods/seurat_gse180963_cldn4/scripts/download.sh /tmp/gse180963
Rscript methods/seurat_gse180963_cldn4/scripts/analyze.R --data /tmp/gse180963 --out methods/seurat_gse180963_cldn4
```

---

## How to read this

- **Additive public mouse**, not a human concordant-pool join.
- **This is real KL** (Lkb1-targeted Kras lung GEMM), not KLK and not the private 8-KL cohort.
- **Matrix exists.** The no-go is the **mouse n**, not the file.
- **Cldn4 is present and sparse / epithelial.** Score it; do not build a Cldn4-high story on two mice.
- **Cldn4 vs T/NK and epithelial IFN/MHC cannot be tested at the unit.** Two points.
- **Honest n is 1 vs 1 mice.** Direction can be listed. A p-value cannot.
- **No dual-high. No ICI endpoint. Thesis unchanged.**

## Files

- `tables/per_mouse.tsv` — unit-level table
- `tables/kl_vs_k.tsv` — descriptive 1 vs 1
- `tables/cldn4_positive_cells.tsv`
- `tables/compartment_by_mouse.tsv`
- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/summary.json`
- `tables/cldn4_vs_ifn_mhc_cells_exploratory.tsv` — do not cite as n
- `figures/fig_honest_n.png`, `fig_tnk_fraction.png`, `fig_cldn4_epi.png`, `fig_epi_ifn_mhc.png`, `fig_cldn4_vs_tnk.png`, `fig_umap_*.png`
- `scripts/analyze.R` — R + Seurat primary

## 结论

GSE180963 是公开的 K vs KL 肺 GEMM scRNA，GEO **有** 处理后的 10x MTX，因此用 **R + Seurat `CreateSeuratObject`** 做了 Cldn4-only 计分。诚实 n = **2 只鼠（1 vs 1）**，不是 14260 个细胞。Cldn4 行在，主要在上皮且稀疏；Cldn4 对 T/NK、上皮 IFN/MHC 在鼠单位上 **不能做检验**。未打开私有 8 只 KL。无 dual-high，无 ICI，不改 thesis。


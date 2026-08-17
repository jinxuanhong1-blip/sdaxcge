# FINDING — GSE316655 CLDN / LILRB2 myeloid IO (public MTX/TSV)

**Species / model (do not drop):** *Homo sapiens* CD45+ FACS scRNA from **SK-MEL-5** (human melanoma cell line) tumors and peripheral blood of **NSG-SGM3** mice humanized with human cord-blood CD34+ cells. Dual-genome Cell Ranger `GRCh38_and_mm10-2020-A`. Liu et al. *Sci Immunol* 2026 (PMID 41931598, DOI 10.1126/sciimmunol.adt7832). GEO **GSE316655**.

**This is not a patient tumor. This is not a CLDN18.2 gastric xenograft.** The paper’s mechanistic CLDN18.2–LILRB2 work used other models; the public count matrices here are SK-MEL-5 immune cells.

---

## Verdict

The public deposit is **four 10x libraries, one labeled donor (`ND`), one library per tissue × treatment cell**. Myeloid cells are **rare** after a human-high QC (150 / 15,566). **CLDN18 is zero UMIs.** **CLDN4 is 6 cells total, each with 1 UMI** (ambient / leak, not a CLDN program). Cheap CLDN–LILRB ligand–receptor potential is **numerically ~0**.

What *is* present: classical **MHC-I is nearly ubiquitous** (HLA-A 85–99%, B2M ~99.5%). **LILRB2** is detectable in a minority of marker-called myeloid cells (tumor 10% anti-LILRB2 vs 23% isotype; n_myeloid = 50 vs 22). Tumor myeloid **NF-κB / STAT / suppressor** scores do **not** separate the two tumor libraries in a way that survives an honest n=1 vs 1 (cell-level p = 0.46 / 0.44 / 0.62). T/NK cytotoxicity is slightly higher and exhaustion slightly lower in the anti-LILRB2 tumor library, but that library is **T-heavy** and the isotype library is **B-heavy** — composition, not a replicated treatment effect.

**Do not use cell-level p-values as evidence.** Biological n = 1 mouse/library per arm.

---

## Honest n

| Item | n | Note |
| --- | --- | --- |
| Public 10x MTX/TSV libraries | **4** | GSM9457798–801; series matrix lists all four |
| Labeled humanized donors | **1** | GEO title prefix `ND` |
| Tumor anti-LILRB2 libraries | **1** | GSM9457798 `ND_anti_B2` |
| Tumor isotype libraries | **1** | GSM9457799 `ND_CTR` |
| PB anti-LILRB2 libraries | **1** | GSM9457800 |
| PB isotype libraries | **1** | GSM9457801 (FTP/MTX present; some GEO HTML views omit it) |
| Barcodes as deposited | 30,296 | Cell Ranger filtered matrices |
| QC cells (human UMI≥200, genes≥100, human fraction≥0.80, mito≤0.20) | **15,566** | 44.7% of barcodes have human UMI fraction <0.80 (dual-genome) |
| Tumor / PB QC cells | 11,201 / 4,365 | |
| Marker myeloid / T / NK / B | **150 / 8,873 / 498 / 5,074** | not author annotations |
| Tumor myeloid anti vs isotype | **50 vs 22** | too few for a myeloid-checkpoint claim |
| Biological replicates per arm | **1** | no sample-level test is valid |

Lineage is a compact marker argmax (`LYZ/CD14/…`, `CD3D/E/G`, `NKG7/GNLY`, `MLANA/PMEL/TYR`, …), not Seurat clusters and not the paper’s labels.

---

## 1. Claudins (including CLDN4)

FACS was **human CD45+**. Tumor claudin protein is not expected in this matrix except residual melanoma / ambient RNA.

| Gene | Tumor anti %UMI>0 | Tumor iso | PB anti | PB iso | QC cells UMI>0 |
| --- | ---: | ---: | ---: | ---: | ---: |
| **CLDN4** | 0.043 (3 cells) | 0.047 (2) | 0 | 0.034 (1) | **6**, all UMI=1 |
| **CLDN18** | 0 | 0 | 0 | 0 | **0** |
| CLDN1 | 0 | 0 | 0.069 | 0 | 1 |
| CLDN3 | 0.014 | 0.164 | 0 | 0.034 | few |
| CLDN7 | 0.058 | 0.117 | 0.276 | 0.034 | few |
| CLDN12 | 0.68 | 1.85 | 0.35 | 0.48 | highest CLDN (still <2%) |
| CLDN15 | — | — | — | — | max library 3.66% |

Residual marker-called **melanoma** cells (n=315; 153 anti / 157 iso tumor): CLDN4 in 2 cells (0.63%), CLDN18 in **0**. SK-MEL-5 in this deposit is not a CLDN4/CLDN18-high epithelial compartment.

---

## 2. LILRB family (myeloid)

LILR genes are human-only in this reference (no mouse `Lilrb` homologs under those symbols).

| Gene | Tumor anti myeloid n=50 | Tumor iso n=22 | PB anti n=69 | PB iso n=9 |
| --- | ---: | ---: | ---: | ---: |
| LILRB1 %pos | 10.0 | 22.7 | 7.2 | 0 |
| **LILRB2** %pos | **10.0** | **22.7** | 13.0 | 11.1 |
| LILRB3 %pos | 12.0 | 4.5 | 15.9 | 33.3 |
| LILRB4 %pos | 14.0 | 40.9 | 7.2 | 0 |
| **LILRB5** %pos | **0** | **0** | **0** | **0** |

Tumor-wide LILRB2+ cells: 8 (anti) vs 13 (isotype). Medians for myeloid LILRB2 log1p(CP10k) are **0 in both tumor arms** (a few positives pull the mean). Direction is *lower* LILRB2 detection after antibody, but **n_library = 1 vs 1** and n_myeloid = 50 vs 22.

---

## 3. Myeloid NF-κB / STAT

Compact scores = mean log1p(CP10k) of a priori gene lists (`gene_sets.py`). Tumor myeloid only:

| Score | Anti mean (n=50) | Iso mean (n=22) | Cell MWU p | How to read |
| --- | ---: | ---: | ---: | --- |
| NF-κB | 0.313 | 0.243 | 0.46 | no library-level test |
| STAT | 0.264 | 0.277 | 0.44 | no library-level test |
| suppressors (IDO1/IL10/CD274/…) | 0.229 | 0.243 | 0.62 | no library-level test |
| MHC-I cassette | 1.67 | 1.85 | 0.22 | both high |

The paper’s claim that CLDN–LILRB2 tunes myeloid NF-κB/STAT **cannot be confirmed or refuted** from these four libraries. Myeloid n is small, CLDN ligands are absent, and there is no replicate mouse.

---

## 4. T/NK and MHC-I

| Library | QC n | % T | % NK | % B | % myeloid | T/NK cytotox mean | T/NK exhaust mean | MHC-I (all) mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Tumor anti-LILRB2 | 6,940 | 56.4 | 3.6 | 31.3 | 0.72 | 0.293 | 0.161 | 1.82 |
| Tumor isotype | 4,261 | 27.2 | 1.6 | 65.0 | 0.52 | 0.252 | 0.223 | 1.62 |
| PB anti-LILRB2 | 1,449 | 79.6 | 3.4 | 3.3 | 4.76 | 0.254 | 0.200 | 2.02 |
| PB isotype | 2,916 | 90.9 | 4.3 | 2.9 | 0.31 | 0.227 | 0.159 | 2.09 |

Tumor anti vs iso T/NK cell-level MWU p-values are tiny (cytotox 5.6×10⁻⁸; exhaust 1.5×10⁻⁹) **because thousands of cells from one library are treated as replicates**. The libraries also differ in T vs B mix. Report the **direction and the composition confound**; do not report a treatment p-value.

MHC-I (HLA-A/B/C + B2M + TAP1/2 + TAPBP + NLRC5) is high in every library. This is the only ligand class with enough detection to pair with LILRB.

---

## 5. Cheap ligand–receptor (not CellChat)

Naive potential = (% ligand+ cells) × (% receptor+ cells) / 100. No permutation, no spatial data.

| Pair | Tumor anti product | Tumor iso product | Comment |
| --- | ---: | ---: | --- |
| CLDN4–LILRB2 | 5×10⁻⁵ | 1.4×10⁻⁴ | 5 and 2 CLDN4+ cells |
| CLDN4–LILRB5 | ~0 | 0 | LILRB5 ~0 |
| CLDN18–LILRB2 | 0 | 0 | CLDN18 absent |
| HLA-A–LILRB1 | 0.99 | 3.69 | HLA-A ~90%; LILRB1 mostly myeloid-rare |
| HLA-A–LILRB2 | 0.11 | 0.26 | limited by LILRB2 frequency |
| ANGPTL2–LILRB2 | 0 | 4×10⁻⁴ | ANGPTL2 almost absent |
| CD274–PDCD1 | 0.027 | 0.038 | both low |

**The only non-trivial cheap pair is classical MHC-I → LILRB1/2**, because MHC-I is on almost every human cell and LILRB is on a slice of myeloid cells. That is the known MHC-I–LILRB biology, not a CLDN axis.

---

## What this is not

- Not a slide / figure audit of Liu et al.
- Not SRA re-alignment, not Cell Ranger re-run, not author RDS/h5ad.
- Not CellChat / LIANA / NicheNet (pairs are a 26-row curated list).
- Not syngeneic LILRB2-transgenic mice and not MIA PaCa-2 (mentioned in the GEO design text; **no MIA PaCa-2 MTX is in this series**).
- Not evidence that anti-LILRB2 reprograms myeloid NF-κB/STAT in this deposit.
- Not evidence against the paper’s biochemistry or spatial human cohorts — those data are not in GSE316655.

---

## Methods (short)

Public GSM MTX + TSV only (`download.py`). Human vs mouse UMIs split on `GRCh38_` / `mm10___` feature prefixes. Scores = mean log1p(CP10k) of genes present in the matrix (missing genes dropped, not imputed; only `CD45` was absent from the panel, `PTPRC` used). Scripts: `analyze.py`, `gene_sets.py`.

---

## Files

`tables/n_table.tsv`, `library_summary.tsv`, `claudin_detection.tsv`, `lilr_by_lineage.tsv`, `ligand_receptor_naive.tsv`, `exploratory_cell_mwu.tsv`, `gene_coverage.tsv`, `summary.json`.

Figures: `figures/fig1_qc_species` … `fig8_score_overview` (png+pdf). Footer on every panel repeats the species/model string.

---

## 中文

**物种/模型：** 人源化 **NSG-SGM3**（人脐血 CD34+）皮下 **SK-MEL-5** 黑色素瘤，FACS **人 CD45+** 单细胞；双基因组 GRCh38+mm10。不是病人肿瘤，也不是论文里的 CLDN18.2 胃癌模型。

**n 必须写清楚：** 公开 MTX 只有 **4 个 10x 文库**，GEO 标题供者标记 **ND 1 个**，每个组织×处理 **1 个文库**。QC 后 15,566 细胞；髓系标记细胞只有 **150**（肿瘤抗 LILRB2 50 vs 同型 22）。

**结论：** **CLDN18 全程 0 UMI。CLDN4 一共 6 个细胞、各 1 UMI**，不能当屏障程序。CLDN–LILRB 廉价配体–受体乘积≈0。MHC-I 几乎全阳性，才是能和 LILRB1/2 配上的配体。肿瘤髓系 NF-κB/STAT 在 1 vs 1 文库上分不开。抗 LILRB2 肿瘤库偏 T、同型库偏 B，T/NK 分数差别不能当成处理效应。细胞水平 p 值不要当证据。

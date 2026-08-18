# FINDING — GSE137669 mouse Cldn4: NO-GO

**ADDITIVE. Cldn4-only.** Thesis already correct. Do not rewrite concordant-4
(GSE123902 + GSE131907 + GSE205335 + GSE189357). No dual-high Tacstd2×Cldn4.

Assigned accession: GEO [GSE137669](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137669).
Question: mouse-level epithelial **Cldn4 vs T/NK**, and epithelial IFN/MHC
high vs low, on public processed mouse lung-tumor scRNA. Honest n = mice.

**Verdict: NO-GO.** GSE137669 is not mouse lung-tumor scRNA. The series
matrix has **0** expression rows. The only processed count file is bulk
left-ventricle edgeR. **Cldn4 is present as a gene row and is all zeros.**
Epcam and Tacstd2 are also all zeros. There is no epithelial compartment
and no T/NK fraction to score. Stop.

---

## Decision

| Question | Answer |
|---|---|
| Is GSE137669 mouse lung-tumor scRNA? | **No.** Bulk polyA RNA-seq of **left ventricle** after TAC / sham |
| Species | Mus musculus (C57BL/6N male), n=8 libraries |
| Tissue | Left ventricle (heart failure diet study), PMID [36626303](https://pubmed.ncbi.nlm.nih.gov/36626303/) |
| Series-matrix expression rows | **0** |
| Public processed count file | `GSE137669_edgeR_normalized_counts.txt.gz` (bulk, 8 columns) |
| RAW.tar | 1005.7 MB of **bigWig** (RSeQC), not a cell-by-gene matrix |
| Cldn4 on the processed file? | Gene row **ENSMUSG00000047501 / Cldn4** exists |
| Cldn4 expression | **0 / 0 / 0 / 0 / 0 / 0 / 0 / 0** |
| Epcam / Tacstd2 | **all zeros** |
| Mouse-level Cldn4 vs T/NK | **not computed** — no epithelium, no scRNA T/NK fraction |
| Epithelial IFN/MHC Q4 vs Q1 | **not computed** — no epithelium; Q4 needs n_units≥8 anyway |
| Join concordant pool? | **No** — sign not tested |
| Closest substitute scored for T/NK? | GSE154977 peeked: Cldn4 yes, **T/NK = 0** (tumor-cell 10x) |

Stop rule from the assignment: *write a no-go if the series has no matrix/Cldn4.*
That is this page. The series matrix has no expression. Detectable Cldn4 is
absent (all-zero row). The edgeR file is not an scRNA matrix.

---

## What GSE137669 actually is

Ragni, Ruocco, et al., dietary BCAA / essential-amino-acid formula in
pressure-overload heart failure. Eight-week-old male C57BL/6N mice, SFA vs
SFA-EAA diet, sham vs transverse aortic constriction (TAC). RNA from **LV**
at 4 weeks. NextSeq 500, STAR → GRCm38.87, edgeR.

| Library | Diet | Surgery | Cldn4 | Epcam | Tacstd2 |
|---|---|---|---:|---:|---:|
| Sham_sfa.rep1 | SFA | sham | 0 | 0 | 0 |
| Sham_sfa.rep2 | SFA | sham | 0 | 0 | 0 |
| Sham_sfa.eaa.rep1 | SFA-EAA | sham | 0 | 0 | 0 |
| Sham_sfa.eaa.rep2 | SFA-EAA | sham | 0 | 0 | 0 |
| Tac_sfa.rep1 | SFA | TAC | 0 | 0 | 0 |
| Tac_sfa.rep2 | SFA | TAC | 0 | 0 | 0 |
| Tac_sfa.eaa.rep1 | SFA-EAA | TAC | 0 | 0 | 0 |
| Tac_sfa.eaa.rep2 | SFA-EAA | TAC | 0 | 0 | 0 |

Honest n if this were a bulk heart study = **8 mice** (2 per diet×surgery).
That n is **not** a lung-tumor epithelial unit. Do not quote n=8 as a Cldn4
vs T/NK n. Do not quote n=49,671 genes or the 1.0 GB RAW.tar as a cell n.

Low-level T-lineage genes are detectable in bulk LV (Ptprc max 48.5;
Cd8a max 4.07; Nkg7 max 6.11). That is heart-resident / circulating
immune RNA, not a tumor T/NK fraction.

---

## Assigned tests (empty)

| test | n | status | why |
|---|---:|---|---|
| Public processed scRNA matrix? | **0** | absent | series matrix 0 rows; RAW = bigWig |
| Detectable Cldn4? | 8 zeros | **absent** | gene row exists; expression = 0 |
| Mouse-level epithelial Cldn4 vs T/NK | **0** | not computed | no epithelium; not scRNA |
| Epithelial IFN/MHC Q4 vs Q1 | **0** | not run | no epithelium; n_units would be 8 heart mice, not tumor |
| Join concordant-4 if T/NK inverse | **0** | not joined | no sign to match |

---

## Closest public processed mouse lung-tumor scRNA (peek only)

The assignment allowed a closest substitute **if** GSE137669 failed and the
substitute had Cldn4 **and** immune. Public processed files were checked.
None of the closest peeks is a usable mouse-level Cldn4 vs T/NK test.
Already-scored leftover ICI libraries (GSE297632, GSE285606, GSE133604, …)
were **not** re-opened.

| Accession | What it is | Processed? | Cldn4 | Immune / T/NK | Why not scored |
|---|---|---|---|---|---|
| **GSE137669** (assigned) | Bulk LV, TAC diet | edgeR yes | row yes, **all 0** | bulk only | **NO-GO** |
| GSE137396 | Bulk KP vs KL nodules | yes (opus-gemm) | — | bulk | not scRNA; already used |
| **GSE154977** | Marjanovic KP 30w 10x ± cisplatin, PMID 32707077 | h5 + gene/cell tables | **yes** (5,610 / 11,017 cells >0) | **T/NK = 0** | tumor-cell 10x; no T/NK fraction |
| GSE154989 | Marjanovic SmartSeq2 timecourse | h5 | not required | tumor cells only (`typeID` = stage labels) | no immune compartment |
| GSE127465 mouse MTX | Zilionis *Immunity* 2019 | mouse MTX 166 MB | known present on human side | myeloid atlas | already used (human slingshot PR); not additive here |
| GSE218544 | SCLC GEMM + preSC allograft | RAW.tar MTX 434 MB | not peeked | title says tumor + immune | SCLC, n=4 libraries; skipped after assigned series failed the matrix/Cldn4 gate |

GSE154977 peek (public `GSE154977_mmLung10x_cis_dSp_rawCount.h5` +
`geneTable` / `smpTable`): 11,017 cells, 4 mice
(ND m3/m4, Cis72 m5/m6). Marker epithelium (Epcam or Cdh1+Krt8) = 10,617.
T/NK (Cd3d/e, Cd8a, Nkg7, Ncr1, Epcam wins) = **0**. Cldn4 is real in
tumor epithelium (%pos 0.37–0.66) but there is no same-mouse T/NK
denominator. IFN/MHC Q4 vs Q1 was not run (n_mice = 4 < 8).

---

## Honest n

| item | n | note |
|---|---:|---|
| Assigned GEO series | **1** | GSE137669 |
| Series-matrix expression rows | **0** | metadata only |
| Bulk LV libraries | **8** | 2×2 diet × surgery |
| Cldn4 non-zero libraries | **0** | |
| Epcam non-zero libraries | **0** | |
| Lung-tumor scRNA cells | **0** | |
| Mouse-level Cldn4 vs T/NK pairs | **0** | |
| Epithelial IFN/MHC Q4 vs Q1 units | **0** | |
| Concordant-pool join | **0** | |
| Dual-high Tacstd2 × Cldn4 | **0** | not defined |
| GSE154977 T/NK cells (peek) | **0** | 4 mice; epithelium present |

---

## Methods (locked)

- Cldn4 only. Tacstd2 is an audit gene, never a gate.
- Public processed GEO only. SRA / FASTQ / 1.0 GB bigWig RAW.tar not used.
- Files used: `GSE137669_series_matrix.txt.gz` (2.8 KB) and
  `GSE137669_edgeR_normalized_counts.txt.gz` (1.8 MB).
- Cldn4 = `ENSMUSG00000047501`. Zero in all 8 columns is “no Cldn4”.
- Spearman would require n≥4 finite **epithelial** pairs. Q4 vs Q1 would
  require n_units≥8 **epithelial** mice. Neither exists here.
- Concordant-4 vectors are the locked PR #459 tables. Not re-audited.
- Closest peek: GSE154977 public h5 + CSVs; marker gate as leftover mouse
  scRNA (Epcam or Cdh1+Krt8; T/NK = Cd3d/e, Cd8a, Nkg7, Ncr1; Epcam wins).
  Do **not** use Sftpc/Scgb1a1.
- Reproduce: `python3 methods/gse137669_mouse_cldn4/inventory.py`

---

## How to read this

- This is a **no-go**, not a negative Cldn4–T/NK correlation.
- Do not call GSE137669 mouse lung tumor. Do not call n=8 a tumor n.
- Do not treat the all-zero Cldn4 row as “Cldn4 present in heart.”
- Do not add GSE154977 to the concordant pool. It has no T/NK.
- No dual-high. No CellChat. Thesis unchanged.

---

## 中文摘要

GSE137669 **不是**小鼠肺肿瘤单细胞：是 C57BL/6N 左室 bulk RNA-seq（TAC/假手术 ×
SFA/SFA-EAA 饮食，PMID 36626303）。系列矩阵表达行 = **0**。公开 counts 为
edgeR 8 列。Cldn4（ENSMUSG00000047501）行存在但 **8/8 为 0**；Epcam、Tacstd2
同样全 0。无上皮室，无 T/NK 比例，未做小鼠水平 Spearman，未做 IFN/MHC
Q4–Q1，未加入 concordant-4。最近公开处理后选 GSE154977（KP 30 周 10x）有
Cldn4，但标记 T/NK = **0**（肿瘤细胞文库），不能做 Cldn4 vs T/NK。按任务：
**到此停止**。诚实 n = **0** 小鼠肺肿瘤单位。

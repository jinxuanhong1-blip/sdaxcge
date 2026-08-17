# FINDING — pairwise CLDN4-only: GSE154826 + GSE131907

**Verdict: stop. GSE154826 is unusable for patient-level malignant CLDN4 vs T/NK.** The pair was not built. ρ / Q4 were not computed. CellChat-style outgoing CLDN4-high → T/NK was not run. Honest pair n = **0**.

This page is **additive CLDN4-only**. Dual-high / TACSTD2 gates were not used. GSE131907 was **not** downloaded as a UMI matrix and was **not** scored, because the GSE154826 gate failed.

---

## Gate

Assigned test: **pairwise combo (not a pile)** of Leader CITE-seq **GSE154826 + GSE131907 if both have public processed matrices <2 GB and CLDN4+T/NK**, then patient-level malignant CLDN4 vs T/NK on this pair only. If ρ or Q4 differs, add CellChat-style outgoing CLDN4-high → T/NK. **If GSE154826 unusable, say so and stop.**

| Series | Public processed expression | Size | Usable for malignant CLDN4 vs T/NK |
|---|---|---:|---|
| [GSE154826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154826) | GEO series matrices are **metadata only** (8.6K + 6.3K). Supplementary is a **pile of 77 unlabeled 10x MTX tarballs** (3.567 GB / 3.32 GiB together). No labeled TME / malignant UMI on the series record. | each MTX <2 GB; **pile 3.567 GB** | **no** — CD45+ immune-enriched CITE-seq; epithelium gated as `epi_endo_fibro_doublet`; **0** barcodes labeled malignant |
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | `raw_UMI_matrix.txt.gz` / `.rds.gz` and `normalized_log2TPM_matrix.rds.gz` | **389.8 / 604.2 / 825.7 MiB** (all <2 GB). log2TPM txt is **2.86 GB** and was not needed. | would have been size-eligible **and** has epithelial / malignant + T/NK labels; **not used** |

Stop rule from the assignment: *If GSE154826 unusable, say so and stop.* That is this page.

---

## Why GSE154826 fails

Leader, Grout et al., *Cancer Cell* 2021 (PMID [34767762](https://pubmed.ncbi.nlm.nih.gov/34767762/)). GEO title: “CITEseq analysis of non-small-cell lung cancer lesions.” Series design (author text): tissues were **enriched for CD45+ cells** by bead positive selection or FACS; a minority of libraries used dead-cell depletion or CD2+ TCR capture. This is an **immune CITE-seq**, not a tumor-epithelium atlas.

1. **No public labeled processed matrix <2 GB on GEO.**
   - `GSE154826-GPL18573_series_matrix.txt.gz` = **8,816 B**; `GSE154826-GPL24676_series_matrix.txt.gz` = **6,500 B**. Both are sample metadata. Expression rows = **0**.
   - Supplementary expression is `GSE154826_amp_batch_ID_*.tar.gz` × **77** (Content-Length sum **3,567,029,508 B**). Each file is an unlabeled Cell Ranger feature-barcode MTX. Assignment said **not a pile**. Matching 77 MTX dumps to GitHub barcodes and re-clustering 362k cells is not a public processed TME matrix.
   - Author `lung_ldm.rd` (Dropbox, not GEO) HEAD = **2,123,085,302 B** (2.123 GB decimal / 1.977 GiB). Over a 2 GB decimal gate; not a GEO deposit; still the CD45+ Mount Sinai object.

2. **CLDN4+T/NK in the malignant sense is not there.**
   - Author `annots_list.csv` lineages on 361,929 barcodes: T 157,641; NK 11,456; `epi_endo_fibro_doublet` 33,759; **Malignant = 0**. T/NK exist. The failure is the malignant CLDN4 arm.
   - Author immune vs epithelial-gated DE (`immune_vs_ep_de.csv`): CLDN4 l2fc = **−7.42** (immune `fg_exprs` 3.21×10⁻⁶ vs epi `bg_exprs` 7.17×10⁻⁴). Residual immune CLDN4 is not tumor-cell CLDN4.
   - Real epithelium without CD45 depletion: **2 digest patients** (695 LUAD, 706 LUSC). That is not a patient-level malignant-CLDN4 vs T/NK cohort.
   - CITE-seq ADT panel is immune-focused. There is **no CLDN4 / Trop-2 antibody** in the public record used here.

Do not treat the leftover `epi_endo_fibro_doublet` gate as malignant CLDN4. Do not treat immune-cell CLDN4 UMI as tumor CLDN4.

---

## Honest n

Patient is the unit for the assigned pair. Do not write n = 361,929 or n = 35 as a malignant-CLDN4 n.

| Item | n | Note |
|---|---:|---|
| Public GSE154826 series | **1** | GEO 2020; PMID 34767762 |
| GSE154826 GEO labeled processed TME / malignant matrix | **0** | series matrices metadata-only; MTX pile unlabeled |
| GSE154826 GEO series-matrix expression rows | **0** | 8.6K + 6.3K |
| GSE154826 GEO MTX tarballs | 77 | unlabeled; pile 3.567 GB |
| GSE154826 author-annotated barcodes | 361,929 | `cell_metadata.csv` (GitHub, not a GEO matrix) |
| GSE154826 Mount Sinai patients (Table S1) | 35 | excludes Lambrechts / Zilionis re-use rows |
| GSE154826 digest patients (real epithelium) | **2** | 695, 706 |
| GSE154826 cells labeled malignant | **0** | no Malignant field |
| GSE154826 T + NK barcodes | 169,097 | present; not the failure |
| GSE154826 usable matrix for malignant CLDN4 | **0** | **stop** |
| GSE131907 size-eligible matrices <2 GB | 3 | raw UMI txt/rds + log2TPM rds; UMI **not opened** |
| GSE131907 annotation cells | 208,506 | annotation only |
| GSE131907 samples | 58 | paper: 44 patients |
| GSE131907 epithelial / malignant-subtype / T / NK | 36,467 / 24,784 / 79,676 / 11,551 | labels only |
| GSE131907 samples with malignant **and** T/NK | 21 | would have been pair-eligible if 154826 passed |
| Pair patients (154826 + 131907) | **0** | 154826 gate failed |
| Patient-level malignant CLDN4 vs T/NK | **0** | empty |
| CellChat-style CLDN4-high → T/NK | **0** | ρ / Q4 not computed |
| Dual-high TACSTD2 × CLDN4 | **0** | not defined |

---

## What was opened

**GSE154826 (compact public tables only):**

| File | Bytes | Role |
|---|---:|---|
| GEO series matrices (2 platforms) | 8,816 + 6,500 | metadata; no expression |
| `GSE154826_sample_annots.csv.gz` | 1,339 | prep = CD45+ beads / FACS / digest / CD2 |
| `cell_metadata.csv` + `annots_list.csv` + `table_s1_sample_table.csv` | GitHub | barcode → cluster → lineage; 35 Mount Sinai patients |
| `immune_vs_ep_de.csv` | GitHub | CLDN4 epithelial restriction |

GEO MTX tarballs and Dropbox `lung_ldm.rd` were **not** downloaded.

**GSE131907 (annotation + HEAD only):**

| File | Bytes | Opened? |
|---|---:|---|
| `GSE131907_Lung_Cancer_cell_annotation.txt.gz` | 1,886,187 | yes — labels |
| `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` | 408,736,818 | **no** |
| `GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz` | 633,500,069 | **no** |
| `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds.gz` | 865,847,230 | **no** |
| `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz` | 3,065,932,358 | **no** (>2 GB) |

Annotation `Cell_subtype = Malignant cells` is **24,784** barcodes in **mBrain (15,423) / tL/B (6,400) / mLN (2,961)** — not tLung (tLung malignant label n = **0**; tLung epithelial = 7,270). That is a Kim et al. label fact, not a CLDN4 score. 21 of 58 samples have both malignant and T/NK barcodes. None of that was used as a pair substitute.

---

## Assigned tests (empty)

| Test | n | Status | Why |
|---|---:|---|---|
| Pair 154826 + 131907 | **0** | stop | 154826 has no usable processed malignant-CLDN4 matrix |
| Patient malignant CLDN4 vs T/NK ρ | **0** | empty | no 154826 malignant CLDN4 vector |
| Q4 split | **0** | empty | same |
| CellChat-style CLDN4-high → T/NK | **0** | not run | requires ρ or Q4 to differ on a real pair |
| Dual-high | **0** | not defined | CLDN4-only |

No ligand table was invented. No GSE131907-only substitute was scored under this folder name.

---

## What this is not

- Not a claim that GSE154826 is private. The MTX pile and Dropbox object exist. They are the wrong object for **malignant** CLDN4 under a <2 GB labeled-matrix gate.
- Not a T/NK-intrinsic CLDN4 analysis. T/NK are abundant here; that is a different question.
- Not a GSE131907-only re-score. The assignment was a **pair** gated on both series.
- Not dual-high TACSTD2/CLDN4.
- Not CellChat R. The package was not called.

---

## 结论

GSE154826（Leader 2021 CITE-seq）GEO 上没有可用的 **恶性上皮 processed matrix**：系列 matrix 只有元数据（8.6K+6.3K）；补充文件是 **77 个未标注 MTX，合计 3.57 GB**（任务写明不要 pile）。作者对象在 Dropbox（2.12 GB），且实验是 **CD45+ 富集**，上皮被打成 `epi_endo_fibro_doublet`，**0** 个 malignant 标签。CLDN4 在作者 immune-vs-epi DE 里是上皮基因（l2fc −7.42），免疫里几乎测不到。按任务门控：**到此停止**。未与 GSE131907 做 pairwise ρ / Q4，未跑 CellChat-style。配对层面 malignant CLDN4 vs T/NK 的诚实 n = **0**。GSE131907 本身有 <2 GB 的 UMI/RDS 和上皮/恶性 + T/NK 标注（注释层 21/58 样本同时有恶性与 T/NK），但 **未打开表达矩阵、未单独顶替**。未编造表达值。

Files: `tables/honest_n.tsv`, `tables/geo_file_inventory.tsv`, `tables/summary.json`, `tables/gse154826_lineage.tsv`, `tables/gse131907_celltype.tsv`, `scripts/inventory.py`.

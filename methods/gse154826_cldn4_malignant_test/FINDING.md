# FINDING — GSE154826 Leader CITE-seq / LCAM: malignant CLDN4 test

**ADDITIVE. CLDN4-only. No dual-high.** This folder does not rewrite prior TACSTD2–LCAM hunts or the 154826+131907 pair stop. It answers one gate: do author annotations contain a **real malignant / epithelial compartment** with enough patients to score CLDN4 vs same-patient T/NK?

**Verdict: malignant absent. CD45+ / immune-only. Stop.**

Author malignant n_units = **0**. Real no-CD45 digest epithelium = **2** patients (695 LUAD, 706 LUSC), which is **&lt;8**. Score threshold (n_units ≥ 10) is not met. CLDN4 %pos vs T/NK was **not** computed. Concordant-pool join was **not** done. GSE123902 / GSE131907 were **not** opened and were **not** merged.

TACSTD2 and CLDN4 **are** in the public 10x features. That does not create a malignant compartment.

---

## Decision

| Question | Answer |
|---|---|
| Author lineage named Malignant / tumor cell / epithelial? | **No.** 60 clusters → T, MNP, B&plasma, NK, mast, pDC, and `epi_endo_fibro_doublet` |
| Author barcodes labeled malignant | **0** / 361,929 |
| Design | CD45+ bead enrichment or FACS on almost every Mount Sinai library; CITE-seq / LCAM immune atlas |
| Leftover `epi_endo_fibro_doublet` | 33,759 cells; `lig_rec_group = gate`. Mixed epi / endo / fibro / doublet. **Not** author-malignant |
| Tumor patients with leftover-gate ≥20 **and** T/NK ≥20 | 30 / 35 — temptation only; **not scored** |
| No-CD45 digest / dead-cell Tumor patients | **2** (695, 706) |
| Malignant n_units | **0** (&lt;8 → stop; ≥10 would have scored) |
| CLDN4-only %pos vs same-patient T/NK | **empty** |
| Malignant IFN/MHC Q4 vs Q1 | **not run** |
| Join concordant pool? | **No** — sign not tested; T/NK inverse not evaluated |
| Force-merge GSE123902 / GSE131907? | **No** |
| Dual-high TACSTD2 × CLDN4 | **not defined** |

Stop rule from the assignment: *If it is CD45+/immune-only or malignant n&lt;8, write that and STOP.* That is this page.

---

## What the public objects actually are

Leader, Grout et al., *Cancer Cell* 2021 ([PMID 34767762](https://pubmed.ncbi.nlm.nih.gov/34767762/)). GEO [GSE154826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154826): “CITEseq analysis of non-small-cell lung cancer lesions.”

| Public object | Size | Used here |
|---|---|---|
| `GSE154826_sample_annots.csv.gz` | 1,339 B | **yes** — prep / tissue / patient |
| Series matrices GPL18573 + GPL24676 | 8,816 + 6,500 B | **yes** — metadata; **0** expression rows |
| `GSE154826_amp_batch_ID_*.tar.gz` × **77** | **3,567,029,508 B** together | listing + HEAD of all; **one** smallest tar peeked for features only (`amp_batch_ID_26`, 6.27 MB). Pile **not** downloaded |
| Author `cell_metadata.csv` + `annots_list.csv` + Table S1 | GitHub [effiken/Leader_et_al](https://github.com/effiken/Leader_et_al) | **yes** — labels, not a UMI matrix |
| Author `immune_vs_ep_de.csv` | GitHub | **yes** — compact DE, not a matrix |
| Dropbox `lung_ldm.rd` | HEAD **2,123,085,302 B** | **no** — extra huge object |

GEO series matrices have no expression. Supplementary expression is 77 unlabeled Cell Ranger MTX tarballs. Matching that pile to GitHub barcodes and calling malignant cells is not an author annotation. It was not done.

---

## Author annotations (the compartment question)

`annots_list.csv` lineages on 361,929 barcodes:

| lineage | n cells | role |
|---|---:|---|
| T | 157,641 | T/NK present |
| MNP | 106,662 | immune |
| B&plasma | 41,422 | immune |
| **epi_endo_fibro_doublet** | **33,759** | leftover **gate**, not malignant |
| NK | 11,456 | T/NK present |
| mast | 9,082 | immune |
| pDC | 1,907 | immune |
| Malignant | **0** | **absent** |

Token search on every cluster row (malignant / tumor cell / cancer / epithelial / epi): the only hits are the 11 `epi_endo_fibro_doublet` clusters, all tagged `lig_rec_group=gate`. There is no Malignant field.

GEO `sample_annots` prep (human NSCLC rows):

| prep | Tumor sample rows | Tumor patients |
|---|---:|---:|
| CD45+ bead enrichment | 33 | 28 |
| CD45+ FACS | 6 | 5 |
| CD2+ bead enrichment | 2 | 2 |
| dead cell depletion kit (Table S1: digest/dead cell) | 11 | **2** |

This is an immune CITE-seq. Epithelium was gated out. T/NK abundance is not the failure.

### Leftover gate is not a malignant n=30

Thirty Tumor patients have ≥20 leftover-gate barcodes **and** ≥20 T/NK. That is CD45-enrichment leak plus mixed stroma/doublets. Authors did not split epithelium from endothelium / fibroblasts / doublets. Concordant-pool work already dropped GSE154826 gate rows as all-epithelial leftovers. Treating the gate as malignant CLDN4 would invent a compartment the authors did not call.

The two digest patients are the only libraries without CD45 depletion (real tissue digest):

| patient | histology | Tumor cells | leftover gate | T/NK | author malignant |
|---|---|---:|---:|---:|---:|
| 695 | LUAD | 15,757 | 2,590 | 5,725 | 0 |
| 706 | LUSC | 24,409 | 2,752 | 10,668 | 0 |

Two points. n=2 &lt; 8. Not a patient-level malignant-CLDN4 vs T/NK cohort.

---

## Features (peeked public batch tar)

Smallest public tarball: `GSE154826_amp_batch_ID_26.tar.gz` (6,266,201 B; patient 464 Normal, CD45+ FACS). Members are a standard 10x `features.tsv` / `barcodes.tsv` / `matrix.mtx`. Features = **33,694** Gene Expression rows.

Present: **TACSTD2, CLDN4, EPCAM, KRT19, PTPRC, CD3D, NKG7**.

The MTX was not scored. The other 76 tarballs were not downloaded.

Author immune-vs-epi DE (immune = foreground) on the leftover gate:

| gene | l2fc (immune vs epi) | immune expr | epi-gate expr |
|---|---:|---:|---:|
| CLDN4 | **−7.42** | 3.21×10⁻⁶ | 7.17×10⁻⁴ |
| TACSTD2 | −5.53 | 1.13×10⁻⁵ | 5.71×10⁻⁴ |
| EPCAM | −7.83 | 2.27×10⁻⁶ | 7.43×10⁻⁴ |
| KRT19 | −7.82 | 4.91×10⁻⁶ | 1.33×10⁻³ |
| PTPRC | +3.68 | 3.24×10⁻⁴ | 2.44×10⁻⁵ |

CLDN4 is an epithelial-restricted gene in the authors’ own DE. Residual immune CLDN4 is not tumor-cell CLDN4. That still does not mint an author-malignant cohort.

CITE ADT panel (from prior public record; not re-downloaded): EPCAM / CD8 / CD45 exist on hashed libraries. **No CLDN4 / Trop-2 antibody.**

---

## Honest n

Patient is the unit. Do not write n = 361,929 or n = 35 as a malignant-CLDN4 n. Do not write n = 30 leftover-gate patients as malignant.

| item | n | note |
|---|---:|---|
| Public GSE154826 series | **1** | GEO; PMID 34767762 |
| GEO series-matrix expression rows | **0** | 8.6K + 6.3K metadata |
| GEO amp-batch MTX tarballs | 77 | unlabeled; pile 3.567 GB; not downloaded as a pile |
| GEO sample_annots rows | 119 | includes PBMC / mouse extras |
| Author-annotated barcodes | 361,929 | GitHub labels |
| Author clusters | 60 | no Malignant lineage |
| Author cells labeled malignant | **0** | |
| Leftover epi_endo_fibro_doublet | 33,759 | gate, not malignant |
| T + NK barcodes | 169,097 | present; not the failure |
| Mount Sinai patients (Table S1) | 35 | excludes Lambrechts / Zilionis re-use |
| Tumor patients with author barcodes | 35 | |
| Digest / no-CD45 Tumor patients | **2** | 695, 706 |
| Author-malignant units | **0** | stop (n&lt;8); score bar was ≥10 |
| Author-malignant ∩ same-patient T/NK units | **0** | |
| Tumor patients leftover-gate≥20 and T/NK≥20 | 30 | **not** author-malignant; not scored |
| CLDN4 %pos vs T/NK rows | **0** | malignant absent |
| IFN/MHC Q4 vs Q1 | **0** | optional; not run |
| Concordant-pool join | **0** | sign not tested |
| Dual-high TACSTD2 × CLDN4 | **0** | not defined |
| GSE123902 / GSE131907 matrices opened | **0** | not merged |
| Dropbox lung_ldm.rd downloaded | **0** | 2.123 GB skipped |

---

## Assigned tests (empty)

| test | n | status | why |
|---|---:|---|---|
| Author-malignant compartment present? | **0** units | **absent** | no Malignant label; CD45+ design |
| CLDN4-only %pos vs same-patient T/NK | **0** | not computed | malignant n&lt;8 |
| Malignant IFN/MHC Q4 vs Q1 | **0** | not run | optional; requires malignant n≥10 |
| Join concordant pool if T/NK inverse | **0** | not joined | no sign to match |
| Merge GSE123902 / GSE131907 | **0** | not done | assignment: do not force |

No T/NK scoring table was written. The assignment said to write one **if** malignant is present.

---

## What this is not

- Not a claim that CLDN4 is missing from the features. It is present.
- Not a T/NK-intrinsic CLDN4 analysis. T/NK are abundant; that is a different question.
- Not a leftover-gate CLDN4 vs T/NK substitute. The gate is not author-malignant.
- Not a pair with GSE123902 or GSE131907.
- Not dual-high TACSTD2 × CLDN4.
- Not a download of the 3.57 GB MTX pile or the 2.12 GB Dropbox object.

---

## 结论

GSE154826（Leader 2021 CITE-seq / LCAM）作者注释里 **没有恶性上皮室**：60 个 cluster 无 Malignant；361,929 个 barcode 恶性标签 = **0**。实验是 **CD45+ 免疫富集**，残留上皮打成 `epi_endo_fibro_doublet`（`lig_rec_group=gate`），不是肿瘤上皮。无 CD45 的消化法只有 **2** 个病人（695 LUAD、706 LUSC），n&lt;8。公开 features 里有 TACSTD2/CLDN4（已用最小 batch tar 核对），作者 immune-vs-epi DE 里 CLDN4 l2fc=−7.42，但仍不能当成恶性室。按任务：**到此停止**。未算 CLDN4 %pos vs 同病人 T/NK，未做 IFN/MHC Q4–Q1，未因符号加入 concordant pool，未强行与 GSE123902/GSE131907 合并。无 dual-high。诚实恶性 n = **0**。

Files: `tables/honest_n.tsv`, `tables/summary.json`, `tables/per_patient_tumor_author_counts.tsv`, `tables/author_lineage.tsv`, `tables/geo_sample_prep.tsv`, `tables/features_peek.tsv`, `scripts/inventory.py`.

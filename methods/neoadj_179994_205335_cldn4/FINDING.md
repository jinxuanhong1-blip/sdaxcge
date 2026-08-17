# FINDING — neoadjuvant / ICI-adjacent merge GSE179994 + GSE205335 (CLDN4-only)

**Verdict: stop. GSE179994 has no usable processed matrix for patient-level malignant CLDN4 vs T/NK.** The merge was not built. CellChat-style and LIANA-style were not run. Honest merge n = **0**.

This page is **additive CLDN4-only**. Dual-high / TACSTD2 gates were not used. **GSE207422 was not added** (prior flat CLDN4–T/NK). GSE205335 was not downloaded or scored here because the GSE179994 gate failed.

---

## Gate

Assigned test: new merge of public neoadjuvant / ICI-adjacent scRNA **GSE179994 (LUAD pembro+chemo) + GSE205335 if both have usable processed matrices <2 GB**, then patient-level malignant CLDN4 vs T/NK, then CellChat-style + LIANA-style if n_patients ≥ 8.

| Series | Processed expression on GEO | Size | Usable for malignant CLDN4 vs T/NK |
|---|---|---:|---|
| [GSE179994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994) | `GSE179994_all.Tcell.rawCounts.rds.gz` only | **421.0 MiB** (441,470,188 B) | **no** — T-cell-restricted; no epithelial / malignant compartment |
| [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) | `GSE205335_Lung_IO_UMI_matrix.rds.gz` | **499.5 MiB** (523,720,227 B) | would have been size-eligible; **not used** |
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | — | — | **excluded by assignment** |

Both deposited count objects are <2 GB. Size is not the failure. **GSE179994 does not deposit a full TME / epithelial UMI.** Liu et al. *Nat Cancer* 2022 (PMID [35121991](https://pubmed.ncbi.nlm.nih.gov/35121991/)) released temporal scRNA + paired TCR on sorted / annotated **T cells** from 47 biopsies / 36 NSCLC patients on pembrolizumab + chemo. GEO states **raw data not provided**; processed data on the series record are T-cell counts, T-cell metadata, scTCR, PBMC bulk TCR, and a 44.1 MiB `RAW.tar` that is one blood T-cell Seurat (`GSM5444629` / P019), not a tumor epithelium matrix.

Stop rule from the assignment: *If GSE179994 has no usable processed matrix, say so and stop.* That is this page.

---

## Honest n

Patient is the unit for the assigned merge. Do not write n = 150,849 or n = 36 as a malignant-CLDN4 n.

| Item | n | Note |
|---|---:|---|
| Public GSE179994 series | **1** | GEO 2022; PMID 35121991 |
| GSE179994 processed TME / epithelial / malignant matrix | **0** | not on series supplementary; raw FASTQ not provided |
| GSE179994 T-cell barcodes (metadata) | 150,849 | `GSE179994_Tcell.metadata.tsv.gz` |
| GSE179994 patients in that T-cell table | 36 | P1–P30, P33–P38 |
| GSE179994 samples in that T-cell table | 47 | pre / post biopsy labels |
| GSE179994 cells labeled malignant or epithelial | **0** | `celltype` = CD8 59,656 / CD4 58,706 / NA 32,487 |
| GSE179994 usable matrix for malignant CLDN4 | **0** | stop |
| Merge patients (179994 + 205335) | **0** | 179994 gate failed; 205335 not opened |
| Patient-level malignant CLDN4 vs T/NK | **0** | empty |
| CellChat-style (n_patients ≥ 8) | **0** | not run |
| LIANA-style (n_patients ≥ 8) | **0** | not run |
| GSE207422 patients added | **0** | excluded |
| Dual-high cells / patients | **0** | not defined |

Do not treat the T-cell RDS as a stand-in malignant matrix. CLDN4 on CD4/CD8 barcodes is not tumor-cell CLDN4.

---

## What was opened (GSE179994)

GEO FTP `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/`:

| File | Bytes | Role |
|---|---:|---|
| `GSE179994_all.Tcell.rawCounts.rds.gz` | 441,470,188 | T-cell raw counts only |
| `GSE179994_Tcell.metadata.tsv.gz` | 935,750 | `cellid patient sample celltype cluster` |
| `GSE179994_all.scTCR.tsv.gz` | 6.5 MiB (listing) | TCR, not expression |
| `GSE179994_PBMC.bulkTCR.tsv.gz` | 2.6 MiB (listing) | bulk TCR |
| `GSE179994_RAW.tar` | 46,213,120 | `filelist.txt` lists only `GSM5444629_P019.blood.post.TCR.tsv.gz` and `GSM5444629_P019.blood.post.Tcell.Seurat.object.rds.gz` |

Metadata columns are only T-cell states (`Non-exhausted`, `Tex`, `CD4_C8-Treg`, naive / Tcm / Tem, `XCL1`, proliferating). There is no `Malignant`, `Epithelial`, `EPCAM+`, or lineage field that would let a malignant CLDN4 score be computed from the deposit.

GSE205335 HEAD: UMI RDS 523,720,227 B and `GSE205335_Lung_IO_CellIdentity.txt.gz` 735,956 B are present and <2 GB. They were **not** fetched. A prior public slice on that matrix (author `Malignant cells`, patient-level CLDN4 vs T/NK n = 22) is a different folder and is not re-run or merged here.

---

## Assigned tests (empty)

| Test | n | Status | Why |
|---|---:|---|---|
| Merge 179994 + 205335 | **0** | stop | 179994 has no usable processed TME matrix |
| Patient malignant CLDN4 vs T/NK fraction | **0** | empty | no malignant CLDN4 vector from 179994 |
| CellChat-style Mal-CLDN4-high → T/NK | **0** | not run | requires n_patients ≥ 8 on a real merge |
| LIANA-style LR | **0** | not run | same |

No ligand table was invented. No GSE205335-only substitute was scored under this folder name.

---

## What this is not

- Not a claim that GSE179994 is private or empty. The T-cell RDS is public and <2 GB. It is the wrong compartment for **malignant** CLDN4.
- Not a T-cell-intrinsic CLDN4 analysis. That would be a different question and is not this merge.
- Not a GSE205335-only re-score. The assignment was a **new merge** gated on both matrices.
- Not GSE207422. That series was left out on purpose.
- Not dual-high TACSTD2/CLDN4.
- Not CellChat R or LIANA Python. Those packages were not called.

---

## 结论

GSE179994 公开的 processed 表达是 **T 细胞 only**（421 MiB RDS；36 人 / 47 样本 / 150,849 barcode；标注 CD4/CD8）。**没有** 可用的恶性上皮 processed matrix，GEO 也未放 raw。按任务门控：**到此停止**。未与 GSE205335 合并，未跑 CellChat / LIANA。合并层面 malignant CLDN4 vs T/NK 的诚实 n = **0**。未编造表达值。

Files: `tables/honest_n.tsv`, `tables/geo_file_inventory.tsv`, `tables/summary.json`, `tables/gse179994_tcell_celltype.tsv`, `tables/gse179994_tcell_cluster.tsv`, `tables/gse179994_tcell_samples.tsv`.

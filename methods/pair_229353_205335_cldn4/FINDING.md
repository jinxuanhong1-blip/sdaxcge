# FINDING — pairwise combo GSE229353 + GSE205335 (CLDN4-only)

**Verdict: stop. GSE229353 has no malignant epithelium.** The assigned pair with GSE205335 was not built. Patient-level malignant CLDN4 vs T/NK was not scored. CellChat-style was not run. Honest pair **n = 0**.

This page is **additive CLDN4-only**. Dual-high / TACSTD2∩CLDN4 was not used. **Not a bigger pile** — GSE207422 and other series were not added. GSE205335 was not downloaded or re-scored here because the GSE229353 gate failed.

---

## Gate

Assigned test: pairwise combo of public **GSE229353 (neoadjuvant pembro+chemo CD45+ 10x) + GSE205335**. If both have malignant cells and T/NK, do patient-level CLDN4 vs T/NK, then CellChat-style if that contrast differs. Assignment stop: *CD45+ may lack malignant epithelium — if CLDN4/malignant absent, write that and stop.*

| Series | What GEO deposits | Size | Malignant epithelium? | Usable for malignant CLDN4 vs T/NK |
|---|---|---:|---|---|
| [GSE229353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229353) | `GSE229353_RAW.tar` — seven 10x MTX libraries (`*_C_{barcodes,features,matrix}`) | **155.5 MiB** (163,010,560 B) | **no** — every GSM is `cell type: CD45+ immune cells` | **no** |
| [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) | `GSE205335_Lung_IO_UMI_matrix.rds.gz` + `GSE205335_Lung_IO_CellIdentity.txt.gz` | **499.5 MiB** UMI (523,720,227 B) | yes (author `Malignant cells`) | would have been pair-eligible; **not used** |

Size is not the failure. Both deposits are <2 GB. **GSE229353 does not deposit an unsorted TME or epithelial UMI.** Hui et al. *npj Precis Oncol* 2023 (PMID [37231145](https://pubmed.ncbi.nlm.nih.gov/37231145/)) ran 10x on **CD45 Microbead-sorted immune cells** from seven resected NSCLC tumors (6 NAPC / 1 NAC). GEO sample titles are `Neoadjuvant-P0x_C`. There is no second `_E` / epithelial library, no cell-identity table, and no malignant barcode set.

CLDN4 the **gene** is on the 10x feature list (P01 `GSM7159183_P01_C_features.tsv.gz`: 33,538 features; `CLDN4` once, index 12425). That is not malignant CLDN4. Sparse CLDN4+ barcodes in a CD45+ library are compatible with ambient RNA / doublet / residual-tumor leak (already treated as immune-library leak in the CD45 extra; not re-opened here). Do not treat leak counts as a stand-in malignant score.

---

## Honest n

Patient is the unit for the assigned pair. Do not write n = 7 libraries or n = 31,203 barcodes as a malignant-CLDN4 n.

| Item | n | Note |
|---|---:|---|
| Public GSE229353 series | **1** | GEO 2023; PMID 37231145 |
| GSE229353 GSM libraries | 7 | P01–P07; all `cell type: CD45+ immune cells` |
| GSE229353 epithelial / malignant processed matrix | **0** | not on series supplementary |
| GSE229353 cell-identity / malignant labels | **0** | no identity file on GEO |
| GSE229353 patients with malignant + T/NK for CLDN4 | **0** | CD45+ sort removes epithelium by design |
| Pair patients (229353 + 205335) | **0** | 229353 gate failed; 205335 not opened |
| Patient-level malignant CLDN4 vs T/NK | **0** | empty |
| CellChat-style (only if the pair contrast exists and differs) | **0** | not run |
| Extra series added | **0** | not a bigger pile |
| Dual-high cells / patients | **0** | not defined |

GSE205335 HEAD confirms the author TME UMI (499.5 MiB) and cell-identity table (735,956 B). A prior public slice on that matrix (author `Malignant cells`, patient-level CLDN4 %pos vs T/NK n = 22, ρ = −0.435, p = 0.043; Q4 vs Q1 n = 6 vs 6) lives in another folder and is **not re-run or substituted** for this pair.

---

## What was opened (GSE229353)

GEO FTP `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229353/suppl/`:

| File | Bytes | Role |
|---|---:|---|
| `GSE229353_RAW.tar` | 163,010,560 | only series supplementary; 21 MTX/TSV members |

`filelist.txt` members (all `_C_` = the CD45+ Chromium library):

| GSM | Patient | Title | Treatment (GEO) | Pathology | Files |
|---|---|---|---|---|---|
| GSM7159183 | P01 | Neoadjuvant-P01_C | Chemo | SCC | barcodes / features / matrix |
| GSM7159184 | P02 | Neoadjuvant-P02_C | anti-PD1+Chemo | SCC | barcodes / features / matrix |
| GSM7159185 | P03 | Neoadjuvant-P03_C | anti-PD1+Chemo | SCC | barcodes / features / matrix |
| GSM7159186 | P04 | Neoadjuvant-P04_C | anti-PD1+Chemo | AD | barcodes / features / matrix |
| GSM7159187 | P05 | Neoadjuvant-P05_C | anti-PD1+Chemo | SCC | barcodes / features / matrix |
| GSM7159188 | P06 | Neoadjuvant-P06_C | anti-PD1+Chemo | SCC | barcodes / features / matrix |
| GSM7159189 | P07 | Neoadjuvant-P07_C | anti-PD1+Chemo | AD | barcodes / features / matrix |

Every GSM SOFT line is `cell type: CD45+ immune cells`. Extraction protocol (all seven): Human CD45 Microbeads Isolation Kit, then 10x Chromium 5' v1. Paper methods match: “Single-cell RNA sequencing was performed on CD45+ immune cells isolated from surgically resected fresh tumor tissues of seven NSCLC patients.”

No `Malignant`, `Epithelial`, `EPCAM+`, or unsorted-TME file is listed. Raw FASTQ is in SRA; this assignment uses public processed only.

GSE205335 HEAD: UMI RDS 523,720,227 B and identity 735,956 B are present and <2 GB. They were **not** fetched.

---

## Assigned tests (empty)

| Test | n | Status | Why |
|---|---:|---|---|
| Pair 229353 + 205335 | **0** | stop | 229353 has no malignant epithelium |
| Patient malignant CLDN4 vs T/NK | **0** | empty | no malignant CLDN4 vector from 229353 |
| CellChat-style Mal-CLDN4-high → T/NK | **0** | not run | pair contrast does not exist |
| Dual-high TACSTD2∩CLDN4 | **0** | not used | CLDN4-only |
| Add a third series | **0** | not used | not a bigger pile |

No ligand table was invented. No GSE205335-only substitute was scored under this folder name. No immune-leak CLDN4 vs T/NK was promoted to a malignant test.

---

## What this is not

- Not a claim that GSE229353 is private or empty. The 155.5 MiB TAR is public. It is the **wrong compartment** for malignant CLDN4.
- Not a claim that CLDN4 is missing from the 10x panel. The gene is present. **Malignant cells are absent.**
- Not a T-cell-intrinsic or CD45+ leak CLDN4 analysis. That is a different question (see the CD45 extra) and is not this pair.
- Not a GSE205335-only re-score. The assignment was a **pairwise combo** gated on both series having malignant + T/NK.
- Not GSE207422 or any other add-on.
- Not dual-high TACSTD2/CLDN4.
- Not CellChat R. That package was not called.

---

## 结论

GSE229353 公开的 processed 表达是 **CD45+ 免疫 10x only**（Hui 2023, PMID 37231145；7 个 GSM，全部 `cell type: CD45+ immune cells`；155.5 MiB TAR）。**没有** 恶性上皮 / unsorted TME processed matrix。CLDN4 基因在 feature 表里，但没有 malignant CLDN4 可评。按任务门控：**到此停止**。未与 GSE205335 合并，未跑 CellChat-style。配对层面 malignant CLDN4 vs T/NK 的诚实 n = **0**。未编造表达值，未把 CD45+ leak 当成上皮。

Files: `tables/honest_n.tsv`, `tables/geo_file_inventory.tsv`, `tables/gse229353_samples.tsv`, `tables/summary.json`.

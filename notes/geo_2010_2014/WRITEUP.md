# GEO 2010–2014 leftover human lung ICI / PD-1 / PD-L1 — TACSTD2 / CLDN4

> Parallel leftover slice. Outputs: `notes/geo_2010_2014/`,
> `scripts/geo_2010_2014/`, `results/w200/GEO_2010_2014/`.
> **No invented IDs**: every accession was returned by NCBI E-utilities and
> re-verified against its authoritative GEO SOFT record.

Bilingual write-up: **English first, 中文在后。**

---

## 1. Objective (EN)

Enumerate leftover **GEO Series (GSE) public 2010-01-01 … 2014-12-31** that
concern **lung** and **immune-checkpoint / immunotherapy** (PD-1/PDCD1,
PD-L1/CD274, CTLA-4, ICI drugs, or the broader word “immunotherapy”). This
window is leftover relative to already-mined GEO slices (2015–2018, 2019–2021,
2022–2023, 2026). For every real accession: verify it, download open processed
matrices **< 2 GB**, and measure **TACSTD2 (TROP2)** and **CLDN4** where the
platform actually carries them.

## 2. Pipeline (EN)

| step | script | output |
|---|---|---|
| 1 | `01_search.py` | `search_candidates.tsv`, `search_manifest.json` |
| 2 | `02_verify.py` | `verified_series.tsv`, `leftover_shortlist.tsv`, `suppl_files.tsv` |
| 3 | `03_download.py` | `download_manifest.tsv` + `data/` (git-ignored) |
| 4 | `04_build_clinical.py` | `clinical/*`, `clinical_label_index.tsv` |
| 5 | `05_target_gene_analysis.py` | `analysis/*` |
| 6 | `06_curate.py` | `master_summary.tsv`, `verdict.json` |

Public NCBI endpoints only. No API key. No hard-coded GSE/GPL/GSM IDs.

## 3. Search result (EN)

Live `gds` esearch (2026-08-16):

| query | UIDs |
|---|---|
| lung × ICI/immunotherapy × human × GSE × 2010–2014 | **4** |
| same without organism token | **9** |
| strict text `(lung) AND (PD-L1 OR PD-1 OR CD274 OR PDCD1)` × 2010–2014 | **0** |

Union of the ICI queries: **9 real GSE**. The strict PD-1/PD-L1 phrase probe is
empty — that terminology is largely post-2014. Accessions (NCBI only):

```
GSE27556 GSE35640 GSE41088 GSE44825 GSE46517
GSE54351 GSE54352 GSE54353 GSE54514
```

## 4. Verification / leftover classification (EN)

SOFT-brief + FTP listings for all 9. **0** human + lung + **checkpoint ICI**
(PD-1/PD-L1/CTLA-4 drug or blockade). That is historically expected: NSCLC
nivolumab approval is 2015.

| accession | leftover_reason | n | platform | note |
|---|---|---|---|---|
| **GSE27556** | leftover human lung immunotherapy, **not** checkpoint ICI | 10 | GPL96; GPL570 | HLA ligands / T-cell epitopes for lung-cancer immunotherapy. No response labels. |
| **GSE35640** | leftover human immunotherapy, **not** checkpoint ICI | 65 | GPL570 | recMAGE-A3 **cancer vaccine**. Title mentions NSCLC; **deposited samples are all melanoma** (PMID 23715562). Has responder / non-responder labels. |
| GSE54351 / 54352 / 54353 | mouse lung **Pdl1 expression**, not ICI treatment | 9 / 6 / 15 | GPL6246 | Lkb1/Pten SCC (Xu et al., PMID 24794706). SuperSeries + epithelial + stroma. |
| GSE41088 | non-human | 12 | GPL6887 | H1N1 mouse, not ICI |
| GSE44825 | non-human | 12 | GPL1261 | TB DNA vaccine, not lung cancer ICI |
| GSE46517 | human, not lung | 121 | GPL96 | melanoma progression |
| GSE54514 | human, not lung | 163 | GPL6947 | sepsis whole blood |

## 5. Downloads (EN)

Processed series matrices only (RAW.tar / CEL / filelist skipped as not
processed). **6 files downloaded, 0 failures, 0 ≥2 GB skips.**

## 6. TACSTD2 / CLDN4 (EN)

### 6a. GSE27556 — human lung tissue, no ICI labels (n=8 on GPL570)

All four genes measured. Values as deposited (Affymetrix). n=8 is too small
for a claim.

| gene | n | median | min–max |
|---|---|---|---|
| TACSTD2 | 8 | 5.455 | 4.576–6.748 |
| CLDN4 | 8 | 6.057 | 4.534–7.491 |
| CD274 | 8 | 6.283 | 4.862–9.475 |
| PDCD1 | 8 | 5.772 | 5.521–7.045 |

TACSTD2 vs CD274 Spearman **ρ = 0.24** (n=8). GPL96 split is 2 normal
tissues only; CD274 is absent on U133A.

### 6b. GSE35640 — melanoma MAGE-A3 vaccine, **not lung ICI** (n=65)

Deposited samples are melanoma biopsies before recMAGE-A3. **Not PD-1/PD-L1
ICI. Not lung samples.** Reported only because it is the only leftover with
both genome-wide TACSTD2/CLDN4 and response labels.

Mann-Whitney U, responder (n=22) vs non-responder (n=34); 9 not-evaluable
excluded:

| gene | median R | median NR | U | p |
|---|---|---|---|---|
| TACSTD2 | 2.960 | 2.565 | 417 | **0.47** |
| CLDN4 | 2.358 | 2.361 | 337 | **0.53** |
| CD274 | 5.368 | 4.175 | 523 | 0.012 |

TACSTD2/CLDN4 do not differ by MAGE-A3 response. CD274 is higher in
responders in this vaccine cohort (not an ICI claim). TACSTD2 vs CD274
Spearman **ρ = −0.08** (n=65). TACSTD2 vs CLDN4 **ρ = 0.47**.

### 6c. GSE54351 — mouse Lkb1/Pten lung SCC epithelium (n=9)

Not ICI treatment. Title reports elevated Pdl1. GPL6246 `gene_assignment`
maps Tacstd2 / Cldn4 / Cd274 / Pdcd1 (probes 10545168 / 10534395 / 10462390 /
10356866). Genotype medians (n=3 each):

| genotype | TACSTD2 | CLDN4 | CD274 |
|---|---|---|---|
| wild type | 8.98 | 9.40 | 8.21 |
| Kras G12D | 10.24 | 9.89 | 8.60 |
| Lkb1−/− Pten−/− | 12.58 | 13.02 | 11.02 |

Cd274 is higher in LP epithelium than WT/Kras, consistent with the paper’s
Pdl1 claim. TACSTD2 and CLDN4 track the same direction in this n=3-per-group
array. **Do not treat this as a human ICI result.**

## 7. Honest verdict (EN)

**There is no leftover GEO 2010–2014 human lung PD-1/PD-L1/CTLA-4 ICI
treatment-response cohort.** The window predates NSCLC ICI approval. The two
human “immunotherapy” leftovers are a 10-sample HLA-epitope lung series
without labels (GSE27556) and a melanoma MAGE-A3 vaccine series whose title
mentions NSCLC but deposits no lung samples (GSE35640). TACSTD2/CLDN4 are
measurable on those Affymetrix matrices and are null vs vaccine response.
Mouse GSE54351/2/3 is Pdl1-expression biology, not checkpoint blockade.

---

## 1. 目的（中文）

盘点 **2010-01-01 至 2014-12-31** 这一段、相对后续 GEO 切片（2015–2018 及以后）
而言的 **剩余（leftover）** 肺癌 × 免疫检查点 / 免疫治疗系列，并在平台真实
含有探针时测量 **TACSTD2 / CLDN4**。不编造任何 GSE/GPL/GSM。

## 2. 检索与核实（中文）

NCBI `gds` 现场检索：人源 4 条、不限物种 9 条；字面 PD-1/PD-L1/CD274/PDCD1
查询为 **0**。SOFT 核实后：**人 + 肺 + 检查点 ICI = 0**。符合史实（NSCLC
纳武利尤单抗获批于 2015）。

## 3. 剩余系列（中文）

- **GSE27556**：肺癌 HLA 配体 / T 细胞表位免疫治疗，n=10，无疗效标签。GPL570
  上 TACSTD2/CLDN4/CD274/PDCD1 均可测；n=8 不可外推。
- **GSE35640**：recMAGE-A3 **疫苗**（非 PD-1/PD-L1）。题名含 NSCLC，**入库样本
  全是黑色素瘤**。应答 vs 无应答：TACSTD2 p=0.47，CLDN4 p=0.53。
- **GSE54351/2/3**：小鼠 Lkb1/Pten 肺鳞癌 **Pdl1 表达**，不是 ICI 治疗。
  上皮中 LP 组 Cd274 高于野生型 / Kras，与原文一致。

## 4. 结论（中文）

**2010–2014 GEO 没有剩余的人肺癌 PD-1/PD-L1/CTLA-4 ICI 疗效队列。**
这是空窗，不是检索失败。疫苗 / 表位系列与小鼠 Pdl1 生物学已如实记录，
未把它们写成检查点治疗结果。

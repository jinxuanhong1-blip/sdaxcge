# Leftover GEO lung ICI series: TACSTD2 / CLDN4 vs response

Slice directory prefix: `notes/opus_geo_leftover/`, `scripts/opus_geo_leftover/`, `results/opus_geo_leftover/`.

This slice searched NCBI GEO for **remaining open human lung immune-checkpoint-inhibitor (ICI) bulk series** that are **not** the seven named exclusions (`GSE126044`, `GSE135222`, `GSE136961`, `GSE166449`, `GSE93157`, `GSE207422`, `GSE205335`), verified every downloaded accession against Entrez, skipped FASTQ / files > 2 GB, and tested **TACSTD2 (TROP2)** and **CLDN4** against real clinical endpoints where a usable matrix existed.

No accessions were invented. No p-values were fabricated. Tests that cannot be run (gene absent, no MPR cases in the deposited matrix, raw-only) are reported as such.

---

## 中文摘要

在排除 7 个已被其他切片覆盖的经典队列后，对 GEO 中剩余的**开放、人类、肺癌、ICI、bulk** 表达系列做了穷尽检索与登录号核验。真正同时满足「全转录组 + 肿瘤组织 + 可核验的 ICI 疗效标签 + 测到 TACSTD2/CLDN4」的剩余队列，只剩同一项新辅助 durvalumab ± SBRT 试验的两个系列：

| 系列 | 时点 | n | 可做的终点 |
|---|---|---|---|
| **GSE253564** | 治疗前 FFPE RNA-seq | 32 | MPR（11 vs 21）、病理缓解深度、复发/PFS（7 个事件） |
| **GSE248378** | 术后残留灶 RNA-seq | 29 | **不能做 MPR**（入库样本病理缓解全部 <90%）；可做复发（9 vs 20）与 PFS |

**治疗前 TACSTD2 在 MPR 中更高**（中位 log2(FPKM+1) 7.30 vs 5.89；Mann–Whitney *p* = 0.016；Cliff δ = 0.53；对两个基因的 BH-FDR *q* = 0.031）。方向与「TROP2 高 = ICI 冷/耐药」相反。该关联**不能独立于增殖**：对 MKI67 做 OLS 残差后 *p* = 0.25。阳性对照 MKI67 本身强烈升高（*p* = 2.2×10⁻⁴，δ = 0.81），与原文增殖签名一致，说明标签与表达对接是对的。

**CLDN4 治疗前 vs MPR 为阴性**（*p* = 0.34）。

**术后 TACSTD2 连续值与较差 PFS 相关**（Cox HR = 2.24，95% CI 1.16–4.30，*p* = 0.016；9 个事件）。复发二分 *p* = 0.056，未过 0.05。CLDN4 术后生存为阴性。CD8A 作为对照呈保护方向（HR = 0.53，*p* = 0.011），说明生存分析管道本身能检出已知方向。

其余「看起来像剩余 ICI 队列」的系列：HTG / NanoString / GeoMx CTA 面板**没有这两个基因**（GSE161537、GSE162520、GSE309652、GSE221733）；GSE182328 / GSE111414 等已在兄弟 PR 做完；GSE126045 是被排除系列 GSE126044 的 SuperSeries；GSE248249 的 GPL family SOFT 为 4.8 GB，按 >2 GB 规则跳过。

---

## English summary

After excluding the seven named series, an Entrez + GEO-FTP sweep left **one leftover whole-transcriptome tumour ICI trial** in which TACSTD2 and CLDN4 are actually measured and an outcome can be recovered from open-access papers (GEO itself stores only the randomisation arm):

- **GSE253564** (pre-treatment, n=32): TACSTD2 is **higher** in major pathologic response (MPR 11 vs no-MPR 21; two-sided Mann–Whitney *p* = 0.016; Cliff’s δ = 0.53; BH *q* = 0.031 across the two genes). The association **does not survive residualising on MKI67** (*p* = 0.25). CLDN4 is null (*p* = 0.34).
- **GSE248378** (post-treatment, n=29): the deposited matrix contains **zero MPR tumours** (deepest Table S1 reduction = −80%), so an MPR test is not estimable. Post-treatment TACSTD2 as a continuous Cox covariate associates with worse PFS (HR 2.24, 95% CI 1.16–4.30, *p* = 0.016; 9 events). Recurrence Mann–Whitney *p* = 0.056. CLDN4 is null.

This is **not** independent evidence that TROP2-high tumours are ICI-resistant. The only binary response test that can be run points the other way and is explained by the trial’s own proliferation effect.

---

## 1. Search, exclusions, and what “leftover” means

Scripts: `01_search_geo.py`, `01b_search_supplementary.py`, `02_triage_candidates.py`, `03_probe_metadata.py`, `04_select_datasets.py`, `06_inventory_response_cohorts.py`, `07_scan_supp_clinical.py`, `10_leftover_catalog.py`.

- Four complementary `db=gds` queries (lung ∩ ICI; +response words; +expression type; neoadjuvant) plus a second pass (trial names, biomarker phrasing, pan-cancer ICI).
- Pool after merge: **576 GSE** (427 + 149), of which **461** were expression series and were FTP-probed.
- Hard skips: FASTQ / SRA / BAM; any single supplementary file **> 2 GB** (Content-Length via HTTP HEAD).
- Named exclusions were dropped at search time and again at catalog time. `GSE126045` is the SuperSeries of excluded `GSE126044` and is not an independent cohort.

Sibling year-slices and blood/purity PRs had already finished several leftover-looking accessions (GSE161537/GSE162520 panel-null; GSE182328 Akkermansia surrogate; GSE111414 PBMC; GSE216297 platelets; GSE305086 blood). This slice does **not** re-analyse those. Disposition of every lung+ICI+outcome accession is in `results/opus_geo_leftover/leftover_disposition.tsv`.

PR 11 used GSE253564 / GSE248378 only for TACSTD2–CD8/NK correlation (“GEO carries no responder label”). The leftover contribution here is the **MPR / recurrence / PFS** join from the two open papers.

---

## 2. Accession verification

Every downloaded series was re-queried with `GSE########[Accession]` against Entrez `gds` (script `05_download.py`). All nine fetched accessions returned `verified=true` and the expected organism / sample count.

| Accession | Entrez n | Type | Taxon | PMIDs | GEO |
|---|---|---|---|---|---|
| GSE253564 | 32 | RNA-seq | *Homo sapiens* | 38401548 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564 |
| GSE248378 | 29 | RNA-seq | *Homo sapiens* | 38114518, 38401548 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378 |
| GSE161537 | 82 | targeted RNA (HTG) | *Homo sapiens* | 36111282, 36038492 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161537 |
| GSE162520 | 92 | targeted RNA (HTG) | *Homo sapiens* | 36111282, 36038492 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162520 |
| GSE309652 | 72 | NanoString metabolism | *Homo sapiens* | — | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309652 |
| GSE182328 | 44 | RNA-seq | *Homo sapiens* | 35115705 (paper; not in GEO pubmedids) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182328 |
| GSE248249 | 42 | Clariom D array | *Homo sapiens* | 38215748 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249 |
| GSE216297 | 286 | platelet RNA | *Homo sapiens* | 37270827 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216297 |
| GSE111414 | 20 | PBMC CD8 RNA-seq | *Homo sapiens* | 30765392 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111414 |

SHA-256 and byte sizes for every fetched file: `results/opus_geo_leftover/download_manifest.json`. Largest processed payload was the GSE248249 series matrix (20 MB). The GSE248249 GPL23126 family SOFT is **4,849,895,138 bytes** and was not downloaded.

---

## 3. Why most leftover “hits” cannot test TACSTD2/CLDN4 vs ICI response

| Accession | Why it is not a leftover TACSTD2/CLDN4 vs-response test |
|---|---|
| GSE161537, GSE162520 | HTG EdgeSeq ~2,560-gene immune panel. Grep of the public log2-CPM matrices: **0 TACSTD2, 0 CLDN4**. GSE162520 is surgically treated NSCLC, not ICI-treated (ICI-naive comparator in the same Data in Brief package). |
| GSE309652 | NanoString `NS_Hs_Metabolism_v1.0` RCC files. **0 TACSTD2, 0 CLDN4**. Has R/NR labels that therefore cannot be used. |
| GSE221733 | GeoMx CTA (~68 genes). Neither gene is on the panel. |
| GSE126045 | SuperSeries of excluded GSE126044. |
| GSE182328 | Whole transcriptome + genes present, but GEO has no RECIST/MPR field (Akkermansia only). Already analysed in PR 36. |
| GSE111414, GSE216297 | Blood / platelet. Already in PR 36 / 115. Epithelial genes are at/near zero or filtered. |
| GSE248249 | Acquired-resistance Clariom D; probe IDs only; annotation SOFT > 2 GB. No response field in GEO (pre/post timepoint only). |
| GSE266219 | Seurat RDS of sorted CD8/PBMC, not a bulk tumour matrix. |
| GSE206613 | Colorectal (inventory false positive from shared ICI words). |

---

## 4. Clinical join for the durvalumab ± SBRT trial

GEO sample characteristics for GSE253564 / GSE248378 contain only **tissue / histology / Arm1 vs Arm2**. Outcomes were taken from the two CC-licensed papers of the same trial:

- Cell Reports Medicine 2024, PMID 38401548, PMC10982989, Table S1 (`mmc2.xlsx`): arm, percent reduction in viable cancer cells (“Pathology Response”), pre/post RNA-seq flags.
- Nature Communications 2023, PMID 38114518, PMC10730562, figure source data: recurrence / PFS / DFS and a Major vs None pathologic-response call.

Join key = the integer embedded in the GEO sample title (`durva036`, `pod01`, `45-M-PO`, `POD20R` → `durva036` / `durva001` / `durva045` / `durva020`).

**Validation:** 32/32 GSE253564 and 29/29 GSE248378 titles join to Table S1. Randomisation arm agrees between GEO and Table S1 for **60/61** samples. The single conflict is GSE248378 `45-M-PO` / `durva045`: GEO says Arm2, Table S1 + Nat Commun + the same patient’s GSE253564 GSM + the title token `M` (monotherapy) all say Arm1. Canonical arm = Table S1. That sample is excluded from arm-stratified sensitivity tests. See `clinical_join_validation.json`.

**MPR definition:** Table S1 Pathology Response ≤ −90% (standard ≥90% reduction in viable cancer cells). This matches the Nat Commun “Major” call on every pre-treatment sample that has both fields.

**Structural limit of GSE248378:** every one of the 29 deposited post-treatment tumours has Pathology Response in [−10, −80]. There is **no MPR case** in this matrix. GEO’s own design line is “biomarkers of **recurrence**”, which is consistent with a residual-tumour series. An MPR contrast is therefore **not estimable** here; we do not invent one.

---

## 5. Statistics (real tests only)

Expression = submitter FPKM, analysed as **log2(FPKM + 1)**. Tests: two-sided Mann–Whitney U + Cliff’s δ; Spearman; log-rank (median split) and univariate Cox PH. BH-FDR is applied only to the two primary pre-treatment MPR U-tests (the only binary response tests that exist). No imputation.

### 5.1 GSE253564 pre-treatment (n = 32; MPR 11 / no-MPR 21; 7 PFS events)

| Gene | Test | Result |
|---|---|---|
| TACSTD2 | MWU MPR vs no-MPR | median 7.30 vs 5.89; U-test *p* = **0.016**; δ = 0.53; BH *q* = **0.031** |
| TACSTD2 | MWU Arm2 only (10 vs 6) | *p* = 0.022; δ = 0.70 |
| TACSTD2 | residual after MKI67 OLS | *p* = **0.25**; δ = 0.26 |
| TACSTD2 | residual after EPCAM OLS | *p* = 0.11; δ = 0.35 |
| TACSTD2 | Spearman vs response depth | ρ = 0.32; *p* = 0.072 |
| TACSTD2 | Cox PFS | HR 1.05 (0.64–1.70); *p* = 0.85; 7 events |
| CLDN4 | MWU MPR vs no-MPR | median 5.97 vs 5.66; *p* = 0.34; δ = 0.21; *q* = 0.34 |
| CLDN4 | all other pre-treatment tests | all *p* > 0.5 |
| MKI67 (control) | MWU MPR vs no-MPR | median 5.04 vs 3.13; *p* = **2.2×10⁻⁴**; δ = 0.81 |
| MKI67 | Spearman vs response depth | ρ = 0.58; *p* = 4.9×10⁻⁴ |
| TACSTD2 vs CLDN4 | Spearman | ρ = 0.69; *p* = 1.4×10⁻⁵ |

Interpretation: the TACSTD2–MPR gap is real in the raw bulk numbers and goes **toward responders**, but it is not separable from the trial’s published proliferation effect. It is not evidence for a TROP2-high ICI-resistant state.

### 5.2 GSE248378 post-treatment residual tumours (n = 29; 0 MPR; 9 recurrence / PFS events)

| Gene | Test | Result |
|---|---|---|
| TACSTD2 | MWU MPR vs no-MPR | **not estimable** (0 vs 29) |
| TACSTD2 | MWU recurrence vs no | median 6.34 vs 5.42; *p* = 0.056; δ = 0.46 |
| TACSTD2 | Cox PFS | HR **2.24** (1.16–4.30); *p* = **0.016**; 9 events |
| TACSTD2 | log-rank median split | *p* = 0.17 |
| CLDN4 | MWU recurrence | *p* = 0.31 |
| CLDN4 | Cox PFS | HR 1.41 (0.86–2.30); *p* = 0.17 |
| CD8A (control) | Cox PFS | HR 0.53 (0.32–0.87); *p* = 0.011 |
| TACSTD2 vs CLDN4 | Spearman | ρ = 0.77; *p* = 1.3×10⁻⁶ |

The Cox result is a **small-event** finding (9 events, no multiplicity correction across survival tests). Direction: higher residual TACSTD2, worse PFS. It should not be over-read as confirmatory of a pretreatment predictive biomarker.

Paired pre→post Δlog2(FPKM+1) exists for 20 patients, all no-MPR, so a paired MPR test is also not estimable. Spearman of Δ vs pathologic depth is ~0 (*p* > 0.97).

Figures: `results/opus_geo_leftover/GSE253564_MPR_boxplots.png`, `GSE248378_recurrence_boxplots.png`.

---

## 6. Reproduce

```bash
python3 scripts/opus_geo_leftover/01_search_geo.py
python3 scripts/opus_geo_leftover/01b_search_supplementary.py
python3 scripts/opus_geo_leftover/02_triage_candidates.py
python3 scripts/opus_geo_leftover/03_probe_metadata.py
python3 scripts/opus_geo_leftover/04_select_datasets.py
python3 scripts/opus_geo_leftover/05_download.py
python3 scripts/opus_geo_leftover/06_inventory_response_cohorts.py
# Table S1 + Nat Commun source data are fetched from Europe PMC OA packages
# (PMC10982989, PMC10730562) into /tmp/geo_dl/pmc/ as in the session log.
python3 scripts/opus_geo_leftover/08_build_clinical.py
python3 scripts/opus_geo_leftover/09_analyze_tacstd2_cldn4.py
python3 scripts/opus_geo_leftover/10_leftover_catalog.py
```

Dependencies: `pandas`, `numpy`, `scipy`, `lifelines`, `openpyxl`, `matplotlib`. No FASTQ, no >2 GB files.

---

## 7. Honest limits

- Leftover **whole-transcriptome + tumour + ICI + genes + response** public GEO data, after the seven named series and the sibling year-slices, is essentially **one trial**.
- n = 32 / 29 and 7 / 9 survival events. Only large effects are detectable; a null CLDN4 result does not rule out a modest true effect.
- MPR labels are from the papers, not from GEO. The sample-to-patient map is validated by arm concordance (60/61) but remains an inferred join.
- Bulk FPKM is purity-confounded. The MKI67 residual kills the TACSTD2–MPR gap; we do not claim a tumour-intrinsic pretreatment predictive effect.
- GSE248378 cannot be used as an MPR replication cohort. Anyone pooling “post-treatment RNA-seq vs MPR” from this accession would be fitting a contrast with zero MPR cases.

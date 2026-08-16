# TACSTD2 / CLDN4 vs ICI response in older GEO microarray cohorts

Slice paths: `notes/opus_microarray/`, `scripts/opus_microarray/`, `results/opus_microarray/`.

This write-up is bilingual. English first; 中文在后。

---

## English

### Question

Do **TACSTD2** (TROP2) and **CLDN4** microarray / targeted-panel probes associate with PD-1 / PD-L1 checkpoint-inhibitor **response** in older public GEO series, especially lung ICI cohorts and the named series **GSE93157**?

If a cohort cannot measure the genes, the absence is documented rather than guessed.

### Methods (what was actually done)

1. **Discovery.** NCBI E-utilities (`db=gds`) with four complementary queries (ICI drugs / PD-1 / PD-L1 × lung / NSCLC × human GSE; array-typed; response-worded; NanoString / nCounter / HTG). Seed accessions (GSE93157 and literature ICI series) were force-merged. Exact terms: `results/opus_microarray/geo_search_queries.json`. Hit table: `geo_search_hits.tsv` (520 GSE).
2. **Sample-level screen.** Brief GSM metadata for 243 plausible series. A series was shortlisted only if *per-sample* records mentioned an ICI drug **and** a response / RECIST / PFS field — titles alone are not trusted. Table: `series_screen.tsv`.
3. **Platform gate.** GEO platform tables were scanned for official symbols and aliases (TROP2 / M1S1 / EGP-1; CPE-R / CPETR1 / WBSCR8). Analysis used **official-symbol** matches only (alias-only hits such as MAL / p32 were dropped — they are not CD8A). Table: `platform_summary.tsv`, `platform_probes.tsv`.
4. **Extraction.** Series matrices were streamed; only the needed probe rows were kept. Cached downloads are checksummed in `download_manifest.tsv`.
5. **Association.** Probes were collapsed to the gene (median). Linear-looking matrices (`max > 30`) were `log2(x+1)`-transformed (rank tests are invariant to this monotone map). Two endpoints:
   - **OR** = CR/PR (or `R` / `response`) vs PD/NR (or `SD+PD` when only RECIST is available and OR is defined as CR/PR vs the rest).
   - **DCB** = CR/PR/SD vs PD, only when it differs from OR.
   Pre-specified subsets: GSE93157 **lung only**; GSE202417 **pre-treatment only**. Tests: two-sided Mann–Whitney U, AUC (= U / n1 n0), Cliff’s δ = 2·AUC − 1, 2000-shuffle permutation p, univariate logistic when both classes had n ≥ 3.
6. **Verification.** Eight independent checks in `verification.tsv` (all passed).

No RNA-seq series was analysed as a primary result (GSE135222 / GSE126044 / GSE136961 / GSE207422 are logged as out-of-slice).

### Discovery result

Keyword filters on titles produce many false hits (cell lines, pirfenidone, parasites). After sample-level screening, **array/panel series with ICI + response + lung material** that matter for this slice are:

| GSE | Platform | Assay | TACSTD2 / CLDN4 | Response in GEO? | Verdict |
|---|---|---|---|---|---|
| **GSE93157** | GPL19965 | nCounter PanCancer Immune 730 | **absent / absent** | RECIST `best.resp`, PFS | **Cannot test targets.** Named series. |
| GSE140901 | GPL19965 | same 730-gene panel (HCC) | **absent / absent** | RECIST + CBR | Same panel absence; not lung. |
| **GSE202417** | GPL23126 | Affymetrix Clariom D | present / present | phenotype R / NR | **Analysed.** Blood CD8+ T cells, nivo + bezafibrate, n=14 pre-tx. |
| GSE248249 | GPL23126 | Clariom D | present / present | **no** | Genes fire in NSCLC FFPE; labels not in GEO. |
| GSE141479 | GPL23126 | Clariom D | present / present | **no** | Blood CD8+; labels in the paper, not SOFT. |
| GSE305086 | GPL570 | HG-U133 Plus 2.0 | present / present | **no** (line / timepoint only) | Whole blood; TACSTD2 near floor. |
| **GSE67501** | GPL14951 | Illumina HT-12 WG-DASL | present / present | binary + RECIST | **Analysed.** RCC (not lung), n=11. |
| GSE99070 | GPL10558 | Illumina HT-12 v4 | present / present | treatment vs naive, no RECIST | Mesothelioma. |
| GSE180347 | GPL29738 | NanoString Cancer Immune | absent / absent | no ICI treatment | Excluded. |
| GSE261345 / 261348 | GPL30173 | NextSeq 2000 GeoMx | n/a (sequencing) | RECIST yes | Excluded: not microarray. |

Curated reasons: `catalog_decisions.tsv`.

### Platform verification (the gate)

**GSE93157 / GPL19965** — the series matrix lists **765** feature IDs (`GSE93157_panel_genes.txt`). Neither `TACSTD2` nor `CLDN4`, nor the aliases TROP2 / M1S1 / EGP-1 / CPE-R / CPETR1 / WBSCR8, appear. Immune controls **are** on the panel (`CD274`, `CD8A`, `PDCD1`, `GZMB`, `IFNG`, `CXCL9`, `STAT1`, `HLA-DRA`). This is expected: the nCounter PanCancer Immune Profiling panel is an immune-gene set, not an epithelium set.

The same ID-level absence holds for GSE140901’s supplementary processed matrix (785 IDs; `CD274`/`CD8A`/`PDCD1` present, targets absent). That series’ GEO matrix table is empty — expression lives in `GSE140901_processed_data.txt.gz`.

**Genome-wide arrays that do measure both genes**

| Platform | TACSTD2 probe | CLDN4 probe |
|---|---|---|
| GPL23126 Clariom D (transcript) | `TC0100014340.hg.1` (`NM_002353 // TACSTD2 // …`) | `TC0700007993.hg.1` (`NM_001305 // CLDN4 // …`) |
| GPL570 HG-U133 Plus 2.0 | `202285_s_at`, `202286_s_at`, `202287_s_at`, `227128_s_at` | `201428_at`, `1569421_at` |
| GPL14951 HT-12 WG-DASL | `ILMN_1739001` | `ILMN_2132458` |
| GPL10558 HT-12 v4 | `ILMN_1739001` | `ILMN_2132458` |

### Association results

Full table: `association_results.tsv`. Sample-level values used for the tests: `sample_level.tsv`. Figures: `figures/`.

#### Targets (the question)

| Series | Tissue | Endpoint | Gene | n_resp / n_non | AUC | MW p | perm p | Note |
|---|---|---|---|---|---|---|---|---|
| GSE202417 | NSCLC blood CD8+, pre-tx | OR (R vs NR) | TACSTD2 | 6 / 8 | 0.375 | 0.49 | 0.49 | Near detection floor after log2(x+1); epithelial gene in sorted T cells |
| GSE202417 | same | OR | CLDN4 | 6 / 8 | 0.729 | 0.18 | 0.18 | Same floor / compartment caveat |
| GSE67501 | RCC tumor | OR (response vs no_response) | TACSTD2 | 4 / 7 | 0.357 | 0.53 | 0.53 | Well measured (log2 mean ~9.4 vs 10.2) |
| GSE67501 | RCC tumor | OR | CLDN4 | 4 / 7 | 0.429 | 0.79 | 0.78 | Well measured (log2 mean ~6.63 vs 6.58) |

No target-gene test reaches p < 0.05. Effect sizes are small or sit on background. Summary plot: `figures/target_auc_summary.png`.

**Compartment matters.** TACSTD2 and CLDN4 are epithelial. In GSE202417 the Clariom matrix is linear intensity: GZMB / CD8A reach 10⁴, while TACSTD2 pre-treatment values are ~7–12 (background). Rank tests on that noise are not a test of *tumour* TROP2. GSE248249 (NSCLC FFPE, already-logged) shows TACSTD2 median 5.02 and CLDN4 median 5.14 — the probes *do* fire in lung tumour — but GEO has no responder label, so no test was run. GSE305086 whole-blood TACSTD2 mean 3.61 (U133 log2) is likewise near floor.

#### Named series GSE93157 — controls only (pipeline sanity)

Lung subset n=35 (22 non-squamous + 13 squamous). OR = CR/PR (9) vs SD/PD (26). DCB = CR/PR/SD (21) vs PD (14).

| Endpoint | Gene | AUC | MW p | perm p |
|---|---|---|---|---|
| OR | CD274 | 0.63 | 0.25 | 0.26 |
| OR | CD8A | 0.41 | 0.46 | 0.47 |
| DCB | CD274 | 0.63 | 0.21 | 0.21 |
| DCB | CD8A | 0.65 | 0.15 | 0.14 |
| DCB | PDCD1 | **0.73** | **0.025** | **0.024** |

PDCD1 is higher in DCB than PD in this lung subset (mean 5.77 vs 4.85). That is a weak, pre-specified *control* observation (immune receptor on an immune panel), not a target-gene finding. It shows the phenotype coding and the Mann–Whitney implementation agree with the expected direction for at least one immune gene. OR contrasts are underpowered (9 vs 26) and non-significant.

GSE140901 (HCC, same panel) was not re-tested after a supplementary-matrix column quirk; target absence is already proven at the ID level.

### Verification

`verification.tsv` — 8/8 passed:

| ID | Claim |
|---|---|
| V1 | GSE93157 matrix IDs contain no TACSTD2/CLDN4 alias |
| V2 | External NanoString 730 list was unavailable (404); recorded as skip, V1 is primary |
| V3 | Official-symbol probes present on GPL23126 / GPL570 / GPL14951 / GPL10558; absent on GPL19965 / GPL29738 |
| V4 | Extracted probe IDs match the platform map for GSE93157, GSE202417, GSE67501 |
| V5 | Mann–Whitney U and AUC recomputed from `sample_level.tsv` match `association_results.tsv` to 1e-6 |
| V6 | GSE93157 analysis rows are lung-only (n=35) |
| V7 | GSE202417 analysis rows are pre-treatment only (n=14) |
| V8 | GSE93157 extraction recovered CD274 and CD8A |

### Conclusions

1. **GSE93157 cannot answer the TACSTD2/CLDN4 question.** It is a 730-gene immune panel with excellent ICI + RECIST annotation (including 35 NSCLC tumours) and **no probe for either gene**. Documented by the series matrix ID list, not by a missing-download.
2. **Older public ICI *microarrays* that both measure these genes and carry a response label are scarce.** The two that do — GSE202417 (wrong cell compartment, n=14) and GSE67501 (RCC, n=11) — show **no significant association**.
3. The lung tumour array that *would* be informative (GSE248249, Clariom D FFPE NSCLC, probes expressed) **does not publish response in GEO**.
4. Most contemporary lung ICI expression studies on GEO are **RNA-seq** and were left to other slices.
5. Therefore this slice’s evidence is a **verified absence / underpowered null**, not a biomarker claim.

### Limitations

- Small n; no multiple-testing claim is made beyond the two pre-specified targets.
- GSE202417 treatment is nivolumab **plus bezafibrate**, not PD-1 monotherapy, and the analyte is CD8+ blood RNA.
- GSE67501 is RCC, archival interval to anti-PD-1 is months to years, n=11.
- Response coding follows GEO fields (RECIST / submitter R-NR), not a re-review of scans.
- Alias-based probe matching is noisy; analysis restricted to official-symbol columns.
- We did not impute paper-only labels for GSE141479.

### Reproduce

```bash
python3 scripts/opus_microarray/run_all.sh
# or, if discovery tables are already present:
python3 scripts/opus_microarray/03_platform_probes.py
python3 scripts/opus_microarray/04_extract_expression.py
python3 scripts/opus_microarray/05_analyze.py
python3 scripts/opus_microarray/06_verify.py
```

Requires outbound HTTPS to `eutils.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`, and `ftp.ncbi.nlm.nih.gov`. Bulky caches stay in `results/opus_microarray/cache/` (git-ignored).

---

## 中文

### 问题

在较早期的公开 GEO **芯片 / 靶向面板** 队列（含指定系列 **GSE93157** 及其他肺 ICI / PD-1 系列）中，**TACSTD2（TROP2）** 与 **CLDN4** 探针表达是否与免疫检查点抑制剂疗效相关？若平台测不到这两个基因，则如实记录缺失，而不是用邻近基因替代。

### 方法

1. 用 NCBI E-utilities 四组检索式收集人源 GSE（药物名 / PD-1 / PD-L1 × 肺 / NSCLC；限定 array；含疗效用词；NanoString / nCounter）。种子登录号（含 GSE93157）强制并入。检索式见 `geo_search_queries.json`，命中 520 个 GSE。
2. 对 243 个可疑系列下载 GSM 简要表型：必须在**样本级**同时出现 ICI 用药与疗效 / RECIST / PFS 字段才进入短名单（不信任标题关键词）。
3. 下载平台注释，按官方基因符号 + 别名（TROP2 / M1S1 / EGP-1；CPE-R / CPETR1 / WBSCR8）匹配探针；统计分析只保留官方符号列，去掉 MAL / p32 这类别名误配。
4. 流式读取 series matrix，只抽出目标 / 对照探针。下载带 sha256（`download_manifest.tsv`）。
5. 探针中位数合并为基因。若矩阵最大值 > 30 则 `log2(x+1)`（单调变换不改变秩和检验）。终点：OR = CR/PR（或 R）对 PD/NR；DCB = CR/PR/SD 对 PD。预指定子集：GSE93157 **仅肺**；GSE202417 **仅治疗前**。检验：双侧 Mann–Whitney U、AUC、Cliff’s δ、2000 次置换、n≥3 时的单变量 logistic。
6. 八项独立复核见 `verification.tsv`，全部通过。

RNA-seq 队列（GSE135222 等）只作排除记录，不作为本切片主结果。

### 发现摘要

标题检索噪声很大。样本级筛选后，与本切片真正相关的 array/面板系列见英文表。要点：

- **GSE93157（指定系列）**：nCounter 730 免疫基因面板，有 RECIST 与 PFS，**没有 TACSTD2 / CLDN4**。
- **能测到两基因且 GEO 带疗效标签**的旧芯片队列只有 **GSE202417**（NSCLC 外周血 CD8，纳武利尤单抗 + 苯扎贝特，治疗前 n=14）和 **GSE67501**（肾癌，n=11）。
- 肺肿瘤 Clariom D 系列 **GSE248249** 探针有表达，但 GEO **无疗效字段**。
- GSE261345 / GSE261348 虽有广泛期 SCLC 的 RECIST，但是 NextSeq / GeoMx **测序**，不是芯片。

### 平台核查

GSE93157 矩阵含 **765** 个特征 ID（`GSE93157_panel_genes.txt`），其中无 TACSTD2、CLDN4 及 TROP2 / M1S1 / EGP-1 / CPE-R 等别名；对照基因 CD274、CD8A、PDCD1 等在面板上。GSE140901 补充矩阵 785 个 ID，结论相同。这与 NanoString PanCancer Immune 面板的设计一致：免疫基因，不是上皮基因。

全基因组芯片上的官方探针：Clariom D 为 `TC0100014340.hg.1`（TACSTD2）与 `TC0700007993.hg.1`（CLDN4）；U133 Plus 2.0 为 `202285_s_at` 等与 `201428_at`；Illumina HT-12 为 `ILMN_1739001` 与 `ILMN_2132458`。

### 关联结果

**目标基因：无显著关联。**

- GSE202417：TACSTD2 AUC 0.375，p=0.49；CLDN4 AUC 0.729，p=0.18。矩阵为线性强度，GZMB/CD8A 达 10⁴，而 TACSTD2 治疗前仅约 7–12，属于检测下限。这是在**分选 T 细胞**里测上皮基因，不能代表肿瘤 TROP2。
- GSE67501：肿瘤组织中两基因表达充分（TACSTD2 log2 均值约 9.4 vs 10.2；CLDN4 约 6.63 vs 6.58），但 n=11，AUC 0.36 / 0.43，p>0.5。
- GSE248249 肺 FFPE 中 TACSTD2 中位数 5.02、CLDN4 5.14（已是 log2），说明探针在肺肿瘤里能工作，但无法做疗效检验。

**GSE93157 对照（证明流程而非目标结论）：** 肺亚组 n=35；DCB 下 PDCD1 AUC 0.73，p=0.025（置换 p=0.024），获益组更高。这与免疫面板的预期方向一致，说明表型编码与检验实现可用。OR（9 vs 26）下 CD274 / CD8A 均不显著。

### 核查

V1–V8 全部通过：矩阵 ID 无目标基因；平台有/无与目录一致；抽出的探针 ID 未漂移；从 `sample_level.tsv` 重算的 U 与 AUC 和结果表在 1e-6 内一致；GSE93157 仅肺、GSE202417 仅治疗前；对照基因抽取成功。

### 结论

1. **GSE93157 回答不了 TACSTD2/CLDN4 问题**——不是数据没下下来，而是 730 免疫面板上没有这两基因。
2. 旧公开 ICI **芯片**中，同时具备“能测到这两基因 + 有效疗标签”的队列极少；仅有的两套（血 CD8 n=14；肾癌 n=11）均为**阴性 / 效能不足**。
3. 真正合适的肺肿瘤芯片（GSE248249）未在 GEO 公布疗效。
4. 当代肺 ICI 表达研究多为 RNA-seq，不属于本切片。
5. 本切片的可核实结论是：**缺失已证实，可分析队列为小样本阴性**，不是生物标志物阳性声明。

### 限制

样本量小；GSE202417 为联合苯扎贝特且分析物是血 CD8 RNA；GSE67501 为肾癌、标本距用药间隔长；疗效完全按 GEO 字段，未重审影像；未从论文补全 GSE141479 的标签。

### 复现

见英文 “Reproduce” 一节。`run_all.sh` 可一键重跑。缓存目录 `results/opus_microarray/cache/` 不入库。

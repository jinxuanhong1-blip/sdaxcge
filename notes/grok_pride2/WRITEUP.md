# grok_pride2 WRITEUP — PRIDE / MassIVE / ProteomeXchange lung IO search

Parallel slice outputs only. Protein tables only; no raw MS.

Excluded as already known: **PXD042091**, **PXD059688**, **PXD019061**.

---

## 中文

### 任务

在 PRIDE、MassIVE、ProteomeXchange（含 iProX / jPOST / Panorama）中穷尽检索 **肺 + 免疫治疗 / PD-1 / TROP2 / CLDN4** 蛋白质组，排除上述三个已知号，只收集蛋白鉴定/定量表。

### 方法

1. PRIDE Archive v2 `search/projects?keyword=`：immunotherapy、pembrolizumab、nivolumab、durvalumab、PD-1/PD-L1、checkpoint、TROP2、TACSTD2、CLDN4 等 22 个词，共约 856 个去重 accession。
2. OmicsDI：`id:PXD*` / `id:MSV*` 与 NSCLC + immunotherapy/PD-1 交叉；iProX PROXI 补全未进 PRIDE 的项目文件。
3. ProteomeXchange GetDataset JSON 核对标题/摘要。
4. 文献种子：PSME4 系列、Lehtiö 肺蛋白基因组、iProX 血浆 ICI 队列等。
5. 自动关键词会把有丝分裂 checkpoint、多组织“lung”等算进来，因此人工审阅后只保留高置信列表（`results/grok_pride2/curated_lung_io.csv`）。
6. 下载策略：蛋白表 / mzTab / Spectronaut PG report / xlsx；跳过 `.raw/.wiff/.mzML` 以及 >80 MB 的 SEARCH 压缩包。

### 已排除（不作为新发现）

| Accession | 内容 |
| --- | --- |
| PXD042091 | 转移性 NSCLC 血浆 SWATH-MS，pembrolizumab 应答签名 |
| PXD059688 | 小鼠 LLC1，高剂量抗坏血酸 + anti-PD1 |
| PXD019061 | Panorama 靶向 MS 定量 NSCLC FFPE 中 PD-1/PD-L1 等检查点 |
| MSV000085049 | 上项的 MassIVE 全局蛋白组镜像，不是新研究 |

### A 级：ICI 治疗相关、肺来源、超出排除列表

| Accession | 宿主 | 设计 | 蛋白表 |
| --- | --- | --- | --- |
| **PXD062630** | PRIDE | 33 例晚期 NSCLC，anti-PD1；尿液 EV 宿主+细菌蛋白预测 PFS | 已下载 mzTab（11.6 MB） |
| **PXD031474** | iProX IPX0004062000 | 一线 anti-PD-1 NSCLC 血浆 DIA；S100A8/A9+SAA1/2 升高提示不响应 | 已下载 `Protein_identification.xlsx`；论文 Table 2 已转录 |
| **PXD039141** | iProX IPX0003162000 | 转移性肺鳞/腺癌 ICI 前后血清糖蛋白 | 已下载糖肽搜索 xlsx（非完整蛋白定量矩阵） |
| **PXD058967** | iProX IPX0010564000 | 不可切除 III/IV 期驱动基因阴性 NSCLC，免疫治疗前血浆 DIA | 已下载 Spectronaut protein pivot TSV |
| **PXD066804** | jPOST | 一线 anti-PD1/PD-L1 NSCLC，循环低密度中性粒细胞蛋白组 | PRIDE/iProX 无文件；仅元数据 |
| **PXD078247** | OmicsDI 有记录 | NSCLC ICB 单药血清蛋白组 | PRIDE/iProX API 尚无项目/文件（可能未完全公开） |
| **PXD065735** / MSV000098579 | PRIDE + MassIVE | SCLC，ERBB2 介导免疫逃逸与 ICI 耐药 | 已下载 Spectronaut `PG_Report.tsv` |
| **PXD040761** | PRIDE | 肺腺癌细胞系，durvalumab 增敏化疗 | 仅有 2.6 GB MaxQuant `txt.zip`，按体积上限未下 |

### B 级：肺免疫机制 / 免疫肽组（多数不是 ICI 治疗队列）

- **PXD044086**：人肺癌膜蛋白 TMT，明确写“寻找 PD-1/PD-L1 之外的免疫治疗新靶点”；仓库几乎只有 raw。
- **PXD034772 / PXD043057 / PXD022949 / PXD027766 / PXD058303**：肺癌或胸水 HLA/MHC-I 免疫肽组，与 ICB 失败或抗原呈递相关；文件以肽谱/mzid 为主。
- **PXD042991**：小鼠肺癌免疫应答重建，PISA 蛋白表已下载，含 **Cd274 (PD-L1)**。
- **PXD049438**：肺癌 TIL 细胞治疗时间序列蛋白组。
- **PXD081026**：NSCLC MUC1–MAPK–Elk-1–PD-L1 轴；已下载 `Proteins.tsv`，表内无 PD-L1/TROP2/CLDN4 行。
- **PXD019573 / PXD028364 / PXD037365**：PSME4/PA200 与抗原加工；临床关联 durvalumab。SEARCH tar 达 GB 级，只下了实验设计 xlsx。

### C 级：相邻肺蛋白组（非 ICI 治疗）

- **PXD020191 / PXD020548**（Lehtiö NSCLC 蛋白基因组 / DIA 验证）：与免疫逃逸亚型相关。已下载 PXD020548 组织 DIA 蛋白矩阵，见 TROP2/PD-L1。
- CPTAC LUAD/LSCC 由并行 slice 处理，本 slice 不重复下载。

### TROP2 / CLDN4

**没有**找到以肺 TROP2 或 CLDN4 为题目、且与 ICI 绑定的独立 PX 数据集。

PRIDE 关键词 TROP2=6、TACSTD2=0、CLDN4=0、claudin-4=0。命中的 TROP2 项目是前列腺癌（PXD039272、PXD065965）或 FRTACs 膜蛋白降解（PXD058514/523/539），不是肺 ICI。

在已下载蛋白表中的实际检出：

| 数据集 | TACSTD2 / TROP2 | CLDN4 | PD-1 / PD-L1 |
| --- | --- | --- | --- |
| PXD065735 SCLC ICI 耐药模型 | **有（鼠 Tacstd2）** | **有（鼠 Cldn4）** | 未见 Pdcd1/Cd274 行 |
| PXD020548 Lehtiö NSCLC 组织 DIA | **有（人 TACSTD2）** | 无 | **有 CD274**；无 PDCD1 |
| PXD042991 小鼠肺癌 | 无 | 无（仅 Cldn12 / Cldnd1） | **有 Cd274** |
| PXD031474 / PXD058967 血浆 DIA | 无 | 无 | 无 |
| PXD062630 尿 EV | 无基因名命中 | 无 | 无 |

结论：TROP2/CLDN4 的 MS 定量目前主要来自 **组织/模型蛋白组**（SCLC ICI 模型 + Lehtiö 肺组织），而不是 ICI 血浆/尿液液体活检。血浆 ICI 队列（PXD031474、PXD058967、已排除的 PXD042091）覆盖的是急性期/炎症蛋白，不是这两个膜靶点。

### 下载清单（仅蛋白表）

见 `results/grok_pride2/protein_tables/`：

- PXD062630 mzTab
- PXD031474 Protein_identification.xlsx + 论文 Table 2 CSV
- PXD039141 糖肽 xlsx
- PXD058967 Spectronaut protein pivot
- PXD065735 PG_Report.tsv
- PXD042991 mouse proteins xlsx
- PXD020548 DIA protein matrix（tsv）
- PXD081026 Proteins.tsv
- PXD019573 实验设计表（非定量矩阵）

### 缺口

1. PXD066804、PXD078247 尚无公开蛋白表。
2. PXD040761 / PXD019573 的 MaxQuant/PD 结果被打成 GB 级 SEARCH 包，本 slice 按“只要蛋白表、不要原始/大包”未解压。
3. 无肺 TROP2 ADC ± ICI 的 PX 蛋白组。
4. 无肺 CLDN4 专题 MS 数据集。

---

## English

### Task

Exhaustive PRIDE / MassIVE / ProteomeXchange (including iProX, jPOST, Panorama) search for **lung + immunotherapy / PD-1 / TROP2 / CLDN4**. Exclude PXD042091, PXD059688, PXD019061. Protein identification/quantification tables only; no raw MS.

### Methods

PRIDE Archive v2 keyword search (22 terms) plus OmicsDI `id:PXD*` / `id:MSV*` NSCLC+ICI queries, ProteomeXchange JSON, iProX PROXI file lists, and literature seeds. Automatic keyword scoring over-called mitotic “checkpoint” and multi-tissue “lung” hits (856 unique accessions → 88 auto-keep). Manual curation produced `curated_lung_io.csv`. Downloads were capped at 80 MB and restricted to protein-table-like files.

### Excluded known sets

PXD042091 (pembrolizumab plasma SWATH), PXD059688 (ascorbate + anti-PD1 mouse), PXD019061 (Panorama targeted checkpoints). MSV000085049 is the MassIVE global companion of PXD019061, not a new study.

### New high-confidence ICI-linked lung proteomes

**Patient / clinical ICI proteomes**

- **PXD062630** (PRIDE) — urine EV host+bacterial proteome, 33 advanced NSCLC on anti-PD1; mzTab downloaded. Human proteins MPP5, IGKV6-21, NT5E, KRT27 associate with long PFS; LMAN2, NUTF2, NID1, TNC, IGF1, BCR, GPHN, PPBP with short PFS (Koranyi et al., *Front Immunol* 2025).
- **PXD031474** (iProX IPX0004062000) — first-line anti-PD-1 NSCLC plasma DIA. Protein ID table downloaded. Published DEPs (S100A9/SAA2/S100A8/SAA1 up in non-responders) transcribed from Chao et al., *Clin Exp Immunol* 2022 Table 2.
- **PXD039141** (iProX) — serum glycoproteins before/during anti-PD-1/PD-L1 in metastatic LUAD/LUSC. Glycopeptide search table only.
- **PXD058967** (iProX IPX0010564000) — plasma DIA from unresectable driver-negative NSCLC before immunotherapy. Spectronaut protein pivot downloaded; no TACSTD2/CLDN4/PDCD1/CD274 rows.
- **PXD066804** (jPOST) — LDN myeloid proteome, first-line anti-PD1/PD-L1 resistance. Metadata only; no file API.
- **PXD078247** — serum proteome, NSCLC ICB monotherapy. Indexed by OmicsDI; PRIDE/iProX project/files not resolvable at search time.
- **PXD065735** / **MSV000098579** — SCLC ERBB2-driven immune evasion and ICI resistance. Spectronaut PG report downloaded; **mouse Tacstd2 and Cldn4 are quantified**.
- **PXD040761** — durvalumab + chemo in LUAD lines; only a 2.6 GB MaxQuant zip (not fetched).

**Mechanism / immunopeptidome (B)**

PXD044086 (membrane TMT for new IO targets; raw-only), PXD034772 / PXD043057 / PXD022949 / PXD027766 / PXD058303 (HLA/MHC-I peptidomes), PXD042991 (mouse PISA; **Cd274 present**), PXD049438 (TIL therapy), PXD081026 (MUC1–PD-L1 axis; protein table lacks the three targets), PXD019573 series (PSME4/PA200 + durvalumab association; GB-scale SEARCH archives skipped).

**Adjacent tissue atlases (C)**

PXD020191 / PXD020548 (Lehtiö NSCLC proteogenomics / DIA). PXD020548 tissue matrix contains **TACSTD2 and CD274**, not CLDN4, and is not an ICI-treated cohort. CPTAC LUAD/LSCC are out of scope for this slice.

### TROP2 / CLDN4

No PX dataset is titled as lung TROP2 or CLDN4 immunotherapy proteomics. PRIDE keyword hits for CLDN4/TACSTD2 are zero; TROP2 hits are prostate or FRTACs.

Protein-table detections:

1. **PXD065735** — murine **Tacstd2** and **Cldn4** in an SCLC ICI-resistance model (highest-value hit for this slice’s target genes).
2. **PXD020548** — human **TACSTD2** and **CD274** in untreated/early NSCLC tissue DIA.
3. Plasma/urine ICI datasets (PXD031474, PXD058967, PXD062630) do **not** report these membrane proteins — they report acute-phase and EV cargo proteins instead.

### Files written

| Path | Role |
| --- | --- |
| `scripts/grok_pride2/search_px.py` | repository sweep |
| `scripts/grok_pride2/curate_and_download.py` | curated download |
| `results/grok_pride2/px_catalog.csv` | 856-accession catalog |
| `results/grok_pride2/curated_lung_io.csv` | manual keep list |
| `results/grok_pride2/protein_tables/` | downloaded tables |
| `results/grok_pride2/target_protein_hits.csv` | TROP2/CLDN4/PD-L1 presence/absence |
| `results/grok_pride2/paper_tables/` | published DEP table for PXD031474 |

### Gaps

PXD066804 and PXD078247 lack public protein tables. Several PRIDE SEARCH archives are too large to treat as tables. There is still no public MS protein matrix for a **lung TROP2-ADC ± ICI** or **CLDN4-high NSCLC ICI** cohort.

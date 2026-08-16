# IMpower150/110/130 与 PACIFIC 公共组学审计 / Public-omics audit

检索截止 / Search cutoff: **2026-08-16 UTC**

## 中文

### 结论

严格按“样本明确来自原始临床试验”筛选后，四个队列中只有
**IMpower150** 找到试验专属的公共目录组学记录，而且全部是 EGA
受控访问；文件虽为处理后数据且均小于 2 GB，仍不能匿名公开下载。
IMpower110、IMpower130 和原始 PACIFIC 未找到试验专属的公开组学
accession。

另有一项跨 14 个 Roche 试验的 WGS 药物基因组学研究，明确包含
IMpower110/130/150。其两个汇总统计文件是真正公开、处理后且小于
2 GiB，故已纳入下载清单；但它们是合并分析，不能拆分为三个
IMpower 队列。PACIFIC 不在该研究中。

### Accession 清单

- IMpower150 ctDNA/临床/代码：**EGAS00001006703**
  - **EGAD00001009725**：ctDNA mutation calls、样本状态、患者级特征
  - **EGAD00001009726**：临床表
  - **EGAD00001009764**：R Markdown、PDF、RData
  - **EGAD50000000273**：311 基因列表
- IMpower150 基线肿瘤 RNA-seq TPM 与临床信息：
  **EGAS50000001272 / EGAD50000001814**
- 公开的跨试验 GWAS 汇总统计：LocusZoom **74850**（全队列）与
  **743668**（taxane 亚组）
- IMpower110：未找到直接 EGA accession；注册号 **NCT02409342**
- IMpower130：未找到直接 EGA accession；注册号 **NCT02367781**
- PACIFIC：未找到原始试验直接公共组学 accession；注册号
  **NCT02125461**

所有 EGAF 文件号、精确字节数、访问属性和链接见
`results/gpt_impower/accessions.tsv`。

### 下载判定

- 未下载任何 EGA 文件：EGA API 明确返回 `controlled`，小文件也不等于
  open。
- 下载器只允许显式标注为公开且处理后的两个 LocusZoom 文件，并在
  HEAD 阶段强制 `<2 GiB`，下载后核对长度与 SHA-256。
- `results/gpt_impower/open_processed/manifest.{tsv,json}` 记录下载结果。
  大型 `.gz` 本地缓存不进入 Git；可用脚本确定性重取。

### 边界

PRJNA1026052 是独立的小型 durvalumab-after-CRT TCR 队列，只引用
PACIFIC 作为治疗依据，不是 PACIFIC 原试验，因此排除。IMpower133
（EGAS00001004888）不是 IMpower130；IMvigor130 也不是 IMpower130。

## English

### Bottom line

Under a strict rule requiring explicit provenance from the original trial,
**IMpower150 is the only cohort with trial-specific public-catalog omics
records**. Every such record is controlled-access at EGA. Processed format and
sub-2-GB size do not override the EGA access policy, so none was anonymously
downloaded. No trial-specific public omics accession was found for IMpower110,
IMpower130, or the original PACIFIC trial.

A separate pharmacogenomic WGS analysis across 14 Roche trials explicitly
included IMpower110, IMpower130, and IMpower150. Its two pooled GWAS summary
files are genuinely public, processed, and below 2 GiB, so they are the only
eligible downloads. They do not expose trial-stratified effects or participant
labels and cannot be treated as three cohort datasets. PACIFIC was not part of
that analysis.

### Accession inventory

- IMpower150 ctDNA/clinical/code: **EGAS00001006703**
  - **EGAD00001009725**: ctDNA mutation calls, sample status, patient features
  - **EGAD00001009726**: clinical table
  - **EGAD00001009764**: R Markdown, PDF, and RData
  - **EGAD50000000273**: 311-gene list
- IMpower150 baseline tumor RNA-seq TPM plus clinical metadata:
  **EGAS50000001272 / EGAD50000001814**
- Open cross-trial GWAS summary statistics: LocusZoom **74850** (all patients)
  and **743668** (taxane subgroup)
- IMpower110: no direct EGA accession found; registry **NCT02409342**
- IMpower130: no direct EGA accession found; registry **NCT02367781**
- PACIFIC: no direct public-omics accession found for the original trial;
  registry **NCT02125461**

See `results/gpt_impower/accessions.tsv` for every EGAF file accession, exact
byte count, access classification, and source URL.

### Download decision

- No EGA payload was downloaded: the official API labels every candidate
  `controlled`.
- The downloader allowlists only the two explicitly public processed
  LocusZoom files, enforces `<2 GiB` from HTTP headers, and verifies byte count
  and SHA-256 after transfer.
- Download outcomes are recorded in
  `results/gpt_impower/open_processed/manifest.{tsv,json}`. Large `.gz` cache
  files are intentionally not committed and can be reproduced with the script.

### Scope guardrails

PRJNA1026052 is an independent small durvalumab-after-CRT TCR cohort that cites
PACIFIC; it is not the PACIFIC trial. IMpower133 (EGAS00001004888) is not
IMpower130, and IMvigor130 is a different urothelial-cancer trial.

## Reproduction

```bash
python3 scripts/gpt_impower/audit_ega.py
python3 scripts/gpt_impower/download_open_processed.py
```

Evidence and exclusions are documented in `notes/gpt_impower/SOURCES.md`.

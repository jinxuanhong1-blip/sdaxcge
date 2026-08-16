# Project-wide reproducibility playbook / 项目级可复现性手册

Version / 版本: 1.0  
Applies to / 适用范围: every data-analysis agent and every reported run / 所有数据分析代理及所有报告运行

This is a release gate, not a suggestion. A run is reportable only when its accessions are verified, inputs are cataloged, randomness and environments are captured, the FASTQ size policy is satisfied, and the bilingual write-up reports sample sizes and multiple-testing control.

本手册是发布门槛，并非建议。只有在登录号已核验、输入已编目、随机性与环境已记录、符合 FASTQ 大小规则，且双语报告说明样本量和多重检验控制时，运行结果才可报告。

## 1. Required project artifacts / 必需的项目文件

Paths are relative to the project root. Keep generated evidence as immutable run artifacts.

以下路径均相对于项目根目录。生成的证据应作为不可变运行产物保存。

```text
catalog.tsv
WRITEUP.md
repro/
├── seed.txt
├── conda-explicit.txt
├── pip-freeze.txt
└── accessions/
    └── <ACCESSION>.<retrieved-UTC>.<json|xml|tsv>
```

Rules / 规则:

1. `catalog.tsv` has one row per remote file, including files refused by policy. A biological sample may therefore occur on several rows.
2. `repro/seed.txt` contains exactly one non-negative base-10 integer and a trailing newline.
3. Both environment captures are required and must be generated from the final analysis environment; handwritten placeholders are not valid freezes.
4. `WRITEUP.md` follows the template in §8 with no unresolved `TBD`, `TODO`, `FIXME`, `待补`, or `待定`.
5. Raw verification responses are never edited after capture. If verification is repeated, add a new timestamped file.

1. `catalog.tsv` 每个远程文件一行，也必须包含因规则被拒绝的文件。因此同一个生物学样本可以出现多行。
2. `repro/seed.txt` 只能包含一个非负十进制整数及末尾换行。
3. 必须从最终分析环境生成两种环境记录；手写占位内容不能作为有效的环境冻结文件。
4. `WRITEUP.md` 必须遵循第 8 节模板，不能遗留 `TBD`、`TODO`、`FIXME`、`待补` 或 `待定`。
5. 原始核验响应一经保存不得修改。重新核验时应新增带时间戳的文件。

## 2. Accession verification; never invent identifiers / 登录号核验；严禁虚构标识符

An accession-looking string is not evidence that a record exists. Never guess, synthesize, increment, autocomplete, or “repair” an accession. Never substitute an accession found in a secondary paper without checking the authoritative registry.

形似登录号的字符串不能证明记录存在。严禁猜测、拼接、递增、自动补全或“修复”登录号。不得直接采用二手论文中的登录号，必须到权威数据库核验。

### Required procedure / 必须执行的流程

1. Copy the candidate accession verbatim from the primary source or authoritative search result.
2. Query the authoritative registry over HTTPS. Use exact matching, not a free-text near match.
3. Confirm that the returned record contains the exact accession, is public/retrievable, has the expected organism/assay/layout, and lists the exact remote file used.
4. Save the complete raw response under `repro/accessions/`. Record the request URL and UTC retrieval time in `catalog.tsv`.
5. Use the registry's file byte count before any transfer. Record the registry checksum when available. After an accepted download, compute SHA-256 locally.
6. If the service is unavailable, returns zero or ambiguous records, redirects to a different accession, or does not expose file size, mark the candidate unverified/refused. Do not download it and do not replace it with a plausible ID.

1. 从一手来源或权威搜索结果逐字复制候选登录号。
2. 通过 HTTPS 查询权威数据库，必须精确匹配，不能采用自由文本的相近结果。
3. 确认返回记录包含完全相同的登录号、公开且可获取、物种/实验类型/测序布局符合预期，并列出实际使用的远程文件。
4. 将完整原始响应保存到 `repro/accessions/`，并在 `catalog.tsv` 记录请求 URL 和 UTC 获取时间。
5. 下载前采用数据库给出的文件字节数；如数据库提供校验值也要记录。仅对获准下载的文件计算本地 SHA-256。
6. 若服务不可用、返回零条或含糊结果、跳转到其他登录号，或不提供文件大小，则标记为未核验/拒绝。不得下载，也不得用“看起来合理”的 ID 替代。

Example exact queries (encode the accession as a URL parameter) / 精确查询示例（将登录号编码为 URL 参数）:

```bash
# ENA run metadata; save the response exactly as received.
ACC=ERR000000
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
curl --fail --silent --show-error --location \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=run_accession,scientific_name,library_layout,fastq_ftp,fastq_bytes,fastq_md5&format=tsv" \
  --output "repro/accessions/${ACC}.${STAMP}.tsv"

# NCBI E-utilities exact lookup; inspect the returned count and IDs.
ACC=SRR000000
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
curl --fail --silent --show-error --location \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${ACC}%5BAccession%5D&retmode=json" \
  --output "repro/accessions/${ACC}.${STAMP}.json"
```

The values above are placeholders demonstrating command shape; they are not approved accessions and must never be copied into a real catalog. API fields and endpoints can change, so preserve the response and fail closed if required fields are absent.

以上值只是演示命令格式的占位符，并非已批准的登录号，严禁复制到真实目录中。API 字段和端点可能变化，因此必须保留响应；若缺少必需字段，应按拒绝处理。

## 3. `catalog.tsv` specification / `catalog.tsv` 编写规范

Use UTF-8, Unix newlines, a single header, literal tab separators, and no quoted multiline cells. Sort by `sample_id`, then `accession`, then `file_name`. Do not use spreadsheet date or number formatting. `bytes` is an integer count, not “MB/GB”.

使用 UTF-8、Unix 换行、单一表头、真实制表符分隔，不得使用带换行的引号单元格。按 `sample_id`、`accession`、`file_name` 排序。不得采用电子表格日期或数字格式。`bytes` 必须是整数字节数，不能写“MB/GB”。

Required columns, in this order / 必需列及顺序:

| Column | English definition | 中文定义 |
|---|---|---|
| `sample_id` | Stable project-local biological sample ID; repeats are allowed for multiple files | 项目内稳定的生物学样本 ID；多个文件时允许重复 |
| `accession` | Exact, verified public accession | 经精确核验的公开登录号 |
| `database` | Authoritative registry, e.g. `ENA`, `NCBI_SRA`, `GEO` | 权威数据库，例如 `ENA`、`NCBI_SRA`、`GEO` |
| `source_url` | Exact HTTPS verification/download URL | 精确的 HTTPS 核验/下载 URL |
| `verified_at_utc` | Retrieval time in `YYYY-MM-DDTHH:MM:SSZ` | `YYYY-MM-DDTHH:MM:SSZ` 格式的 UTC 核验时间 |
| `verification_evidence` | Project-relative path to saved raw response | 已保存原始响应的项目相对路径 |
| `file_name` | Registry file name, unchanged | 数据库中的原始文件名，不得修改 |
| `bytes` | Registry or HTTP `Content-Length` byte count | 数据库或 HTTP `Content-Length` 给出的字节数 |
| `decision` | Exactly `accepted` or `refused` | 只能是 `accepted` 或 `refused` |
| `refusal_reason` | `NA` for accepted rows; explicit reason for refused rows | 接受行写 `NA`；拒绝行写明原因 |
| `sha256` | 64 hex characters for accepted local file; `NA` for a file not downloaded | 已下载接受文件写 64 位十六进制 SHA-256；未下载写 `NA` |

Recommended columns after the required set / 建议追加列:

```text
organism	assay	layout	read_pair	registry_checksum	local_path	downloaded_at_utc	notes
```

Header and illustrative row shapes only / 仅演示表头与行结构:

```tsv
sample_id	accession	database	source_url	verified_at_utc	verification_evidence	file_name	bytes	decision	refusal_reason	sha256
sample-A	<VERIFIED_ACCESSION>	ENA	<EXACT_HTTPS_URL>	2026-01-01T00:00:00Z	repro/accessions/<RAW_RESPONSE>	<REMOTE_FILE.fastq.gz>	123456	accepted	NA	<64_HEX_SHA256>
sample-B	<VERIFIED_ACCESSION>	ENA	<EXACT_HTTPS_URL>	2026-01-01T00:00:00Z	repro/accessions/<RAW_RESPONSE>	<REMOTE_FILE.fastq.gz>	2500000001	refused	FASTQ exceeds 2000000000-byte limit	NA
```

Angle-bracket values are placeholders and must be replaced from verification evidence. Never present this example as data.

尖括号内容均为占位符，必须用核验证据中的真实值替换。严禁将此示例作为真实数据。

After download / 下载后:

```bash
sha256sum path/to/accepted.fastq.gz
```

Copy the full digest into `catalog.tsv`; do not copy only a prefix. If the local digest changes, quarantine the file and investigate rather than updating the catalog silently.

将完整摘要复制到 `catalog.tsv`，不得只写前缀。若本地摘要发生变化，应隔离文件并调查，不能静默修改目录。

## 4. FASTQ larger than 2 GB: mandatory refusal / 大于 2 GB 的 FASTQ：必须拒绝

The hard limit is **2,000,000,000 bytes per FASTQ object** (decimal 2 GB), applied to `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz`. The registry-reported object size or trustworthy HTTP `Content-Length` is evaluated before transfer.

硬限制为**每个 FASTQ 对象 2,000,000,000 字节**（十进制 2 GB），适用于 `.fastq`、`.fq`、`.fastq.gz` 和 `.fq.gz`。必须在传输前根据数据库报告的对象大小或可信的 HTTP `Content-Length` 判断。

Decision protocol / 决策流程:

1. Obtain byte size from the verified registry response. A `HEAD` request may corroborate it.
2. If size is greater than `2,000,000,000`, do not start or continue the transfer.
3. Add a `refused` catalog row with the observed byte count, evidence path, reason `FASTQ exceeds 2000000000-byte limit`, and `sha256=NA`.
4. If size is missing, inconsistent, or cannot be trusted, refuse with the actual reason until metadata is resolved.
5. If an oversized file is already present, stop processing, remove it from the analysis input set, retain metadata/evidence, and report the policy violation. Do not hash or inspect its biological contents merely to rescue the run.

1. 从已核验的数据库响应获得字节数，可用 `HEAD` 请求进行佐证。
2. 若大于 `2,000,000,000` 字节，不得开始或继续传输。
3. 在目录中加入 `refused` 行，记录实际字节数、证据路径、原因 `FASTQ exceeds 2000000000-byte limit`，并写 `sha256=NA`。
4. 若大小缺失、不一致或不可信，在元数据问题解决前按真实原因拒绝。
5. 若超限文件已存在，应停止处理、将其移出分析输入集合、保留元数据/证据并报告违规。不得仅为挽救该运行而对其生物学内容进行哈希或检查。

Do not bypass the rule by splitting, range-downloading, streaming into a pipeline, recompressing, renaming, or using Git LFS/object storage. Paired reads are assessed per object, and every accepted mate must independently satisfy the limit.

不得通过切分、分段下载、流式输入管道、重新压缩、重命名、Git LFS 或对象存储绕过规则。双端数据按每个对象分别判断，每个获准的 mate 都必须独立满足限制。

## 5. Random seed / 随机种子

Choose one project seed before analysis and write it once:

分析前选择一个项目种子并一次性写入：

```bash
mkdir -p repro
printf '%s\n' '20260816' > repro/seed.txt
```

`20260816` is an example, not a required seed. Read this value into every stochastic library and command. Record any tool that cannot accept a seed and describe the resulting nondeterminism.

`20260816` 只是示例，并非指定种子。所有含随机性的库和命令都必须读取该值。若工具不支持种子，应记录该工具并说明由此产生的非确定性。

Python example / Python 示例:

```python
from pathlib import Path
import random
import numpy as np

seed = int(Path("repro/seed.txt").read_text().strip())
random.seed(seed)
np.random.seed(seed)
```

Also set framework-, sampler-, cross-validation-, bootstrap-, permutation-, and parallel-worker seeds. A global seed does not automatically control independent generators or nondeterministic GPU kernels.

还必须设置框架、采样、交叉验证、自助法、置换检验和并行 worker 的种子。全局种子不会自动控制独立随机数生成器或非确定性 GPU 内核。

## 6. Environment capture: Conda and pip / 环境记录：Conda 与 pip

Capture after the final successful run, from the same activated environment:

在最终成功运行后，于同一已激活环境中记录：

```bash
mkdir -p repro
conda list --explicit > repro/conda-explicit.txt
python -m pip freeze --all > repro/pip-freeze.txt
```

Do not hand-edit either file. `conda list --explicit` captures exact platform builds; `pip freeze --all` captures packages pip sees, including packaging tools. In `WRITEUP.md`, also report OS/architecture, Python version, Conda environment name, key external tool versions, and accelerator/driver details when relevant.

不得手工编辑上述文件。`conda list --explicit` 记录精确的平台构建，`pip freeze --all` 记录 pip 可见的包（包括打包工具）。`WRITEUP.md` 还应报告操作系统/架构、Python 版本、Conda 环境名、关键外部工具版本，以及相关的加速器/驱动信息。

Useful commands / 常用命令:

```bash
uname -a
python --version
conda info --envs
git rev-parse HEAD
```

Environment capture supports reconstruction but does not prove deterministic results. The seed, immutable inputs, code commit, parameters, and output checksums remain necessary.

环境记录有助于重建，但不能证明结果具有确定性。随机种子、不可变输入、代码提交、参数和输出校验值仍然必需。

## 7. Reporting sample size (`n`) and FDR / 报告样本量（`n`）与 FDR

### Sample size / 样本量

Report `n` at the biological-unit level for every group and analysis stage. Distinguish:

每个分组和分析阶段都必须按生物学单位报告 `n`，并区分：

- enrolled/available / 纳入前可用数；
- passing accession and file policy / 通过登录号及文件规则的数量；
- passing quality control / 通过质量控制的数量；
- analyzed / 实际分析数量；
- independent biological replicates / 独立生物学重复；
- technical replicates, cells, reads, or repeated measures / 技术重复、细胞、reads 或重复测量。

Never call reads, cells, technical replicates, or repeated observations independent `n` when inference is at the donor/sample level. State exclusions with counts and predeclared reasons. For paired designs report complete pairs and unpaired exclusions. For models, state the exact number of observations and independent units used.

当推断单位是供体/样本时，不得把 reads、细胞、技术重复或重复观测当作独立 `n`。必须给出排除数量及预先规定的原因。配对设计应报告完整配对数和未配对排除数。模型应说明实际观测数与独立单位数。

Minimum phrasing / 最低报告格式:

> EN: Group A: `n=12` independent donors available, `n=10` analyzed; two were excluded before testing (one accession unverified, one QC failure). Technical duplicates were averaged within donor and were not counted as independent `n`.
>
> 中文：A 组有 `n=12` 名独立供体可用，实际分析 `n=10`；检验前排除 2 名（1 个登录号未核验，1 个未通过质控）。供体内技术重复先取平均，未计为独立 `n`。

### Multiple testing / 多重检验

For every family of more than one hypothesis, report:

对每个包含多个假设的检验族，必须报告：

1. the family definition and number of tests `m`;
2. raw p-value method and whether tests were one- or two-sided;
3. correction method (default: Benjamini–Hochberg unless justified otherwise);
4. prespecified FDR threshold `q` (for example `q < 0.05`);
5. number passing the threshold, with effect sizes and confidence intervals where appropriate;
6. software/function and version;
7. treatment of missing/non-finite p-values and ties.

1. 检验族定义及检验总数 `m`；
2. 原始 p 值方法，以及单侧或双侧；
3. 校正方法（除非有充分理由，默认 Benjamini–Hochberg）；
4. 预先设定的 FDR 阈值 `q`（例如 `q < 0.05`）；
5. 通过阈值的数量，并在适当时报告效应量和置信区间；
6. 软件/函数及版本；
7. 缺失/非有限 p 值与并列值的处理方式。

Do not write “FDR significant” without `m`, method, and threshold. Do not use raw `p < 0.05` across a multi-hypothesis family. Do not call an adjusted p-value “FDR”; label it `q-value` or `BH-adjusted p-value`. If no multiplicity correction is used, say so explicitly and label findings exploratory.

不得只写“FDR 显著”而不说明 `m`、方法和阈值。不得在多假设检验族中直接使用原始 `p < 0.05`。不得把校正 p 值本身称作“FDR”；应标为 `q 值` 或 `BH 校正 p 值`。若未做多重校正，必须明确说明，并将发现标为探索性。

Minimum phrasing / 最低报告格式:

> EN: We tested `m=18,742` genes as one family using two-sided Wald tests. P-values were adjusted by Benjamini–Hochberg (`statsmodels 0.x`, `multipletests`); `q<0.05` was prespecified. Non-finite p-values were excluded before adjustment and counted (`n=3`). `427/18,739` evaluable genes passed.
>
> 中文：我们将 `m=18,742` 个基因作为一个检验族，采用双侧 Wald 检验。使用 Benjamini–Hochberg 法（`statsmodels 0.x`，`multipletests`）校正 p 值，预设阈值为 `q<0.05`。校正前排除并计数非有限 p 值（`n=3`）。在 `18,739` 个可评估基因中，`427` 个通过阈值。

Version strings above are placeholders and must be replaced with the captured version.

以上版本号是占位符，必须替换为实际记录的版本。

## 8. Mandatory bilingual `WRITEUP.md` template / 强制双语 `WRITEUP.md` 模板

Other agents must copy this structure, replace every bracketed field with evidence-backed content, keep English and Chinese numerically consistent, and remove all instructional comments. Neither language may contain information absent from the other.

其他代理必须复制此结构，用有证据支持的内容替换所有方括号字段，确保中英文数字完全一致，并删除所有说明性注释。任一语言均不得包含另一语言中缺失的信息。

```markdown
# [Analysis title / 分析标题]

Analysis commit / 分析提交: `[full Git SHA]`  
Run UTC / 运行 UTC: `[YYYY-MM-DDTHH:MM:SSZ]`  
Catalog / 数据目录: `catalog.tsv`  
Seed / 随机种子: `[integer from repro/seed.txt]`

## English

### Question and design
[Primary question, outcome, exposure/groups, independent biological unit,
design, and prespecified contrasts.]

### Data provenance and accession verification
[Databases searched; exact verified accessions or catalog reference;
verification date; accepted/refused file counts; refusal reasons; statement
that no unverified identifiers entered analysis.]

### Sample size / 样本量
[For each group/stage: available n, policy-passing n, QC-passing n, analyzed
n, biological unit, technical/repeated measures, and exclusions by reason.]

### Methods and reproducibility
[Preprocessing, QC thresholds, model/test, sidedness, covariates, software
versions, seed propagation, nondeterministic operations, environment files,
and exact command/workflow entry point.]

### Multiple testing / 多重检验
[Hypothesis family, m, raw p-value method, correction method, prespecified q
threshold, handling of invalid values/ties, and number passing.]

### Results
[Effect sizes and uncertainty first; exact n and q/p values; distinguish
confirmatory from exploratory findings. Do not overstate causality.]

### Limitations
[Sampling, missingness, confounding, power, technical limitations,
nondeterminism, and generalizability.]

## 中文

### 研究问题与设计
[主要问题、结局、暴露/分组、独立生物学单位、设计和预设比较。]

### 数据来源与登录号核验
[检索的数据库；经精确核验的登录号或目录引用；核验日期；接受/拒绝
文件数；拒绝原因；声明分析未纳入未经核验的标识符。]

### Sample size / 样本量
[按各组/阶段报告：可用 n、通过规则 n、通过质控 n、实际分析 n、
生物学单位、技术/重复测量，以及按原因分类的排除数。]

### 方法与可复现性
[预处理、质控阈值、模型/检验、单/双侧、协变量、软件版本、种子传递、
非确定性操作、环境文件及精确命令/工作流入口。]

### Multiple testing / 多重检验
[假设检验族、m、原始 p 值方法、校正方法、预设 q 阈值、无效值/并列值
处理，以及通过阈值的数量。]

### 结果
[优先报告效应量和不确定性；给出精确 n 与 q/p 值；区分验证性与探索性
发现。不得夸大因果关系。]

### 局限性
[抽样、缺失、混杂、统计功效、技术限制、非确定性和可推广性。]
```

## 9. Checklist and release gate / 检查脚本与发布门槛

Run from the project root:

在项目根目录运行：

```bash
python methods/repro/check_repro.py --root .
```

The script uses only the Python standard library. It checks required catalog columns and values, accession syntax, saved verification evidence, accepted/refused state, the 2,000,000,000-byte FASTQ limit, accepted-file SHA-256 shape, local oversized FASTQs, seed format, Conda/pip captures, bilingual headings, `n`/multiple-testing sections, and unresolved placeholders.

脚本仅使用 Python 标准库。它检查目录必需列和值、登录号格式、已保存核验证据、接受/拒绝状态、2,000,000,000 字节 FASTQ 限制、接受文件 SHA-256 格式、本地超限 FASTQ、种子格式、Conda/pip 环境记录、双语标题、`n`/多重检验章节以及未解决占位符。

Exit status / 退出状态:

- `0`: automated checks pass / 自动检查通过；
- `1`: reproducibility violations found / 发现可复现性违规；
- `2`: invalid checker invocation / 检查脚本调用无效。

Passing the script is necessary but not sufficient. A reviewer must still compare each accession and file against its raw authoritative response, confirm the English and Chinese narratives are semantically equivalent, verify that reported `n` matches model inputs, and audit the multiple-testing family.

脚本通过是必要条件，但不是充分条件。审核者仍须将每个登录号和文件与权威原始响应逐一比对，确认中英文叙述语义一致，核实报告的 `n` 与模型输入一致，并审核多重检验族。

### Final human checklist / 最终人工检查表

- [ ] Every identifier was copied and verified against an authoritative registry; none was inferred or invented. / 每个标识符均从来源复制并经权威数据库核验，无推断或虚构。
- [ ] Every catalog evidence path opens to an immutable raw response containing the exact accession. / 每个目录证据路径均指向包含精确登录号的不可变原始响应。
- [ ] Every remote file has a byte count before transfer; oversized or unknown-size FASTQs were refused and recorded. / 每个远程文件传输前均有字节数；超限或大小未知的 FASTQ 已拒绝并记录。
- [ ] Every accepted local file has a matching full SHA-256 digest. / 每个获准本地文件都有匹配的完整 SHA-256。
- [ ] One seed is recorded and propagated to every stochastic step; exceptions are disclosed. / 已记录唯一种子并传递至所有随机步骤；例外已披露。
- [ ] Conda and pip captures come from the final successful environment and were not edited. / Conda 与 pip 记录来自最终成功环境且未经编辑。
- [ ] `n` is reported by group and stage at the correct independent biological level, with exclusions. / `n` 按组别和阶段、在正确独立生物学层级报告，并说明排除。
- [ ] Every multiple-testing family reports `m`, method, threshold, invalid-value handling, and discoveries. / 每个多重检验族均报告 `m`、方法、阈值、无效值处理和发现数。
- [ ] English and Chinese contain the same design, counts, thresholds, results, and limitations. / 中英文的设计、数量、阈值、结果和局限完全一致。
- [ ] The automated checklist exits `0`, and its output is retained with the run. / 自动检查退出码为 `0`，其输出已随运行保留。

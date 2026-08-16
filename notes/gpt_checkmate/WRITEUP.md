# CheckMate lung omics: open processed data and EGA catalog

Data cut-off / 数据截止：2026-08-16 UTC

## 中文

### 结论

在指定的 6 项肺癌试验（CheckMate-017、057、227、9LA、816、153）中，本次检索只确认
**CheckMate-153** 有可直接下载的开放处理后组学数据：论文的两个 XLSX 补充表。

- 补充表 1：差异表达分析汇总结果。
- 补充表 2：用于图形复现的完整新抗原免疫原性筛选表。

两者均为开放、处理后数据，远小于 2 GB，已下载到
`results/gpt_checkmate/open_processed/`；大小、SHA-256、来源 URL、工作表名称见
`results/gpt_checkmate/open_processed_manifest.tsv`。

其余五项试验未发现可核实的开放处理后组学文件。这里的“未发现”是有范围限定的：
EGA 全量公开元数据中的显式试验别名匹配和 GEO 精确试验名检索均为零；它不证明
未发布、未公开、命名不完整或今后提交的数据不存在。部分论文明确分析过 RNA-seq、
panel sequencing、TMB 或 ctDNA，但“论文报告做过组学分析”不等于“个体级处理后数据
已开放下载”。

### EGA 目录

EGA 中只匹配到 CheckMate-153：

| 试验 | EGA study | EGA dataset | 访问 |
|---|---|---|---|
| 017 | — | — | 未匹配 |
| 057 | — | — | 未匹配 |
| 227 | — | — | 未匹配 |
| 9LA | — | — | 未匹配 |
| 816 | — | — | 未匹配 |
| 153 | EGAS00001007508（WES）；EGAS00001007509（RNA-seq） | EGAD00001011302 | controlled |

EGAD00001011302 同时关联两个 study，EGA 标注 200 个样本、Illumina HiSeq 4000、
WES 与 RNA-seq。公开文件元数据显示 400 个 `fastq.gz`，合计
1,957,541,154,314 bytes。虽然其中 47 个文件单个小于十进制 2 GB，但它们仍是
**受控访问的原始 FASTQ**，同时违反“open”和“processed”两个下载条件，因此一个也
没有下载。完整机器可读目录在 `results/gpt_checkmate/ega_catalog.json`，简表在
`results/gpt_checkmate/ega_catalog.tsv`。

### 检索和复现

```bash
python3 scripts/gpt_checkmate/catalog_ega.py
python3 scripts/gpt_checkmate/search_geo.py
python3 scripts/gpt_checkmate/fetch_open_processed.py
```

下载脚本同时执行 2,000,000,000-byte 硬上限、流式计数和 XLSX/ZIP 完整性检查；
超过或达到上限即删除临时文件并失败。来源、检索词和局限见
`notes/gpt_checkmate/SOURCES.md`。

## English

### Result

Of the six requested lung trials (CheckMate-017, 057, 227, 9LA, 816, and 153),
only **CheckMate-153** yielded verified, directly downloadable open processed
omics data: two article supplementary XLSX workbooks.

- Supplementary Table 1: summarized differential-expression results.
- Supplementary Table 2: full neoantigen immunogenicity-screen data used for
  figure reproducibility.

Both are open processed data and well below 2 GB. They are stored under
`results/gpt_checkmate/open_processed/`. Their byte sizes, SHA-256 checksums,
source URLs, and worksheet names are recorded in
`results/gpt_checkmate/open_processed_manifest.tsv`.

No verifiable open processed omics file was found for the other five trials.
That is a scoped negative result: explicit trial aliases had no match in the
complete public EGA metadata scan, and exact trial-name searches had no match
in GEO DataSets. It does not prove that unreleased, private, incompletely named,
or future deposits do not exist. Several publications report RNA-seq, panel
sequencing, TMB, or ctDNA analyses, but reporting an omics analysis is not the
same as releasing downloadable participant-level processed data.

### EGA catalog

Only CheckMate-153 matched in EGA:

| Trial | EGA study | EGA dataset | Access |
|---|---|---|---|
| 017 | — | — | no match |
| 057 | — | — | no match |
| 227 | — | — | no match |
| 9LA | — | — | no match |
| 816 | — | — | no match |
| 153 | EGAS00001007508 (WES); EGAS00001007509 (RNA-seq) | EGAD00001011302 | controlled |

EGAD00001011302 is linked to both studies. EGA describes 200 samples, Illumina
HiSeq 4000, WES, and RNA-seq. Its public file metadata lists 400 `fastq.gz`
files totaling 1,957,541,154,314 bytes. Although 47 individual files are below
the decimal 2 GB threshold, they remain **controlled raw FASTQ**, failing both
the “open” and “processed” download criteria; none was downloaded. The complete
machine-readable record is `results/gpt_checkmate/ega_catalog.json`, with a
compact table in `results/gpt_checkmate/ega_catalog.tsv`.

### Search and reproduction

```bash
python3 scripts/gpt_checkmate/catalog_ega.py
python3 scripts/gpt_checkmate/search_geo.py
python3 scripts/gpt_checkmate/fetch_open_processed.py
```

The downloader enforces a 2,000,000,000-byte hard limit using declared and
streamed sizes and validates XLSX/ZIP integrity. A file at or above the limit
is rejected and its partial download removed. Sources, aliases, and limitations
are documented in `notes/gpt_checkmate/SOURCES.md`.

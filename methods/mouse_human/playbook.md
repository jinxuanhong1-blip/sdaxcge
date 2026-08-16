# Mouse ICI vs human ICI playbook: Tacstd2/Cldn4 ↔ TACSTD2/CLDN4

**Outputs live only under `methods/mouse_human/`.**  
**Companion files:** `templates/ortholog_focal_genes.tsv`, `templates/gse239485_style_metadata.tsv`, `templates/contrasts.yaml`, `templates/deseq2_gse239485_style.R`.

How to read: each numbered section is English first, then 中文. Do not mix contrast classes across species.

如何阅读：每一节先英文、后中文。跨物种时不要混用对照类型。

---

## 0. Scope / 范围

### English

This playbook is for **bulk (or pseudo-bulk) tumor RNA-seq** when the scientific question is:

> Does a mouse *Tacstd2* / *Cldn4* (TROP2 / claudin-4) program under immune-checkpoint inhibition (ICI) inform human *TACSTD2* / *CLDN4* biology in ICI-treated patients?

It covers five operational problems:

1. **Orthologs** — 1:1 maps, aliases, paralog traps (*Epcam* / *EPCAM*, other claudins).
2. **Batch** — within-species technical batch vs the species barrier (you cannot ComBat species away).
3. **Syngeneic vs GEMM** — which model class can support which claim.
4. **Treatment-arm DE vs human R/NR** — these are different questions; **you usually cannot compare them directly**.
5. **GSE239485-style templates** — nested combination arms (vehicle / ICI combo / ICI combo + add-on), no outcome labels.

**In scope:** LLC-like syngeneic multi-arm ICI RNA-seq; Tacstd2-KO / humanized-TROP2 add-on studies; human ICI cohorts with explicit pre- vs on-treatment and R vs NR labels.

**Out of scope:** writing an exploit of any system; claiming that mouse on-treatment logFC “validates” human pre-treatment response; treating *Cldn4* as the unique TROP2 binding partner (*Cldn7* / *CLDN7* is the better-documented interactor in TNBC immune-exclusion papers — report it beside *Cldn4*).

### 中文

本手册用于**肿瘤 bulk（或伪 bulk）RNA-seq**，核心问题是：

> 小鼠 ICI 条件下的 *Tacstd2* / *Cldn4*（TROP2 / claudin-4）程序，能否为人类 ICI 队列中的 *TACSTD2* / *CLDN4* 提供可辩护的证据？

覆盖五个操作问题：

1. **直系同源** — 1:1 映射、别名、旁系陷阱（*Epcam* / *EPCAM*、其他 claudin）。
2. **批次** — 种内技术批次 vs 物种屏障（不能用 ComBat 把物种“消掉”）。
3. **同源移植 vs GEMM** — 哪类模型能支撑哪类结论。
4. **处理臂差异表达 vs 人 R/NR** — 问题不同，**通常不能直接比较**。
5. **GSE239485 型模板** — 嵌套组合臂（溶媒 / ICI 联合 / ICI 联合+加药），无疗效标签。

**范围内：** LLC 类同源移植多臂 ICI RNA-seq；Tacstd2-KO / 人源化 TROP2 研究；带明确治疗前/治疗中与 R vs NR 标签的人类 ICI 队列。

**范围外：** 把小鼠治疗中 logFC 说成“验证”了人治疗前疗效；把 *Cldn4* 当成 TROP2 唯一结合伴侣（TNBC 免疫排斥文献中更常写 *Cldn7* / *CLDN7* — 与 *Cldn4* 并列报告，不要互换 ID）。

---

## 1. Orthologs / 直系同源

### English

#### 1.1 Focal 1:1 pairs

| Protein (common) | Mouse | Human | Homology | NCBI (Mm / Hs) | Ensembl (Mm / Hs) |
| --- | --- | --- | --- | --- | --- |
| TROP2 | *Tacstd2* | *TACSTD2* | 1:1 | 56753 / 4070 | ENSMUSG00000051397 / ENSG00000184292 |
| Claudin-4 | *Cldn4* | *CLDN4* | 1:1 | 12740 / 1364 | ENSMUSG00000047501 / ENSG00000189143 |

MGI / HGNC: *Tacstd2* MGI:1861606 ↔ *TACSTD2* HGNC:11530; *Cldn4* MGI:1313314 ↔ *CLDN4* HGNC:2044. UniProt: Q8BGV3 ↔ P09758; O35054 ↔ O14493.

Both *TACSTD2* and *CLDN4* are **intronless** in human; *Tacstd2* is intronless in mouse. That is a mapping convenience (one mature mRNA), not a reason to skip ID-based maps.

Always keep a **related interactor** row: *Cldn7* / *CLDN7* (NCBI 53624 / 1366). TROP2–claudin immune-exclusion work in TNBC ties TROP2 to **claudin-7** and a broader tight-junction program. *Cldn4* is a legitimate barrier-program gene and a frequent co-traveller; it is **not** a drop-in synonym for that interactor.

#### 1.2 Mapping rules (do this, in this order)

1. Map **NCBI Gene ID or Ensembl gene ID**, then attach the official symbol for the target genome.
2. Accept only **one-to-one** orthologs from NCBI Gene / Ensembl Compara / MGI / Alliance. Drop one-to-many and many-to-one for gene-level transfer. If a gene is needed in a set, keep it inside a **set-level** score, not as a 1:1 key.
3. Rebuild the map for the **annotation build** you counted against (GRCm39 / GRCh38, Ensembl release recorded in the methods).
4. After mapping, **assert** that *Tacstd2*→*TACSTD2* and *Cldn4*→*CLDN4* are present and unique. Fail the job if either is missing or collided.
5. Convert gene sets with the same table (`msigdbr`, biomaRt, or a frozen TSV). Do not case-fold symbols as a map (`tacstd2` is fine; `trop2` and `tacstd*` are not).

#### 1.3 Paralog and alias traps

| Trap | Why it fires | Rule |
| --- | --- | --- |
| *Epcam* / *EPCAM* (*Tacstd1* / *TACSTD1*, TROP1) | Same family, same “TACSTD*” / “TROP*” strings | Never collapse *Tacstd2* onto *Epcam*. Keep *Epcam* as a negative-control paralog. |
| *Cldn3*, *Cldn7*, *Cldn1* | Claudin family; TJ programs move together | Score a **tight-junction set**; quote *Cldn4* and *Cldn7* as named genes. |
| Protein aliases (GA733-1, EGP-1, CPE-R, WBSCR8) | GEO `gene_assignment` and vendor panels use aliases | Alias → NCBI ID → official symbol. |
| Humanized *TACSTD2* knock-in | Reads may be human *TACSTD2* in a mouse genome BAM | Count the KI allele as its own feature; do not average with mouse *Tacstd2*. |
| Antibody non-cross-reactivity | hRS7 / sacituzumab often need human TROP2 | Transcript of mouse *Tacstd2* ≠ target engagement of a human-specific ADC. |

#### 1.4 What to transfer across species

Transfer **ranks and programs**, not raw counts:

- 1:1-mapped ranked list (DESeq2 `stat`, or signed −log10 *p* × sign(log2FC)).
- Custom sets: TROP2/claudin TJ barrier; IFN-γ / CXCL9–11; cytotoxic T (*Cd8a*, *Gzmb*, *Prf1* and human counterparts); myeloid *C5ar1* if the mouse design used C5aR1 blockade.
- Single-gene claims for *TACSTD2* or *CLDN4* only after the ID assert in 1.2.

Do not transfer TPM matrices, voom objects, or “the top 50 mouse genes” as an untitled human signature.

### 中文

#### 1.1 焦点 1:1 配对

| 蛋白（俗称） | 小鼠 | 人 | 同源 | NCBI（Mm / Hs） | Ensembl（Mm / Hs） |
| --- | --- | --- | --- | --- | --- |
| TROP2 | *Tacstd2* | *TACSTD2* | 1:1 | 56753 / 4070 | ENSMUSG00000051397 / ENSG00000184292 |
| Claudin-4 | *Cldn4* | *CLDN4* | 1:1 | 12740 / 1364 | ENSMUSG00000047501 / ENSG00000189143 |

MGI / HGNC：*Tacstd2* MGI:1861606 ↔ *TACSTD2* HGNC:11530；*Cldn4* MGI:1313314 ↔ *CLDN4* HGNC:2044。UniProt：Q8BGV3 ↔ P09758；O35054 ↔ O14493。

人 *TACSTD2*、*CLDN4* 均为**无内含子**基因；小鼠 *Tacstd2* 亦无内含子。这只是注释上的便利，不能代替基于 ID 的映射。

始终保留**相关互作**行：*Cldn7* / *CLDN7*（NCBI 53624 / 1366）。TNBC 免疫排斥文献把 TROP2 与 **claudin-7** 及更宽的紧密连接程序联系起来。*Cldn4* 是合理的屏障程序基因、也常共变，但**不是**该互作蛋白的同义词。

#### 1.2 映射规则（按此顺序）

1. 先映射 **NCBI Gene ID 或 Ensembl gene ID**，再挂目标基因组的官方符号。
2. 只接受 NCBI Gene / Ensembl Compara / MGI / Alliance 的 **一对一**直系同源。一对多、多对一不得作为基因级转移键；如需保留，只进**基因集**评分。
3. 映射表必须与计数所用的 **注释版本**一致（GRCm39 / GRCh38，方法中写明 Ensembl release）。
4. 映射后**断言** *Tacstd2*→*TACSTD2*、*Cldn4*→*CLDN4* 存在且唯一。缺失或碰撞则任务失败。
5. 基因集用同一张表转换。禁止用大小写折叠当映射（`tacstd2` 可以；`trop2`、`tacstd*` 不行）。

#### 1.3 旁系与别名陷阱

| 陷阱 | 为何会中 | 规则 |
| --- | --- | --- |
| *Epcam* / *EPCAM*（*Tacstd1* / *TACSTD1*，TROP1） | 同一家族，“TACSTD*” / “TROP*” 字符串相同 | 禁止把 *Tacstd2* 塌缩到 *Epcam*。把 *Epcam* 当阴性对照旁系。 |
| *Cldn3*、*Cldn7*、*Cldn1* | claudin 家族；TJ 程序常一起动 | 评**紧密连接基因集**；*Cldn4* 与 *Cldn7* 作为具名基因引用。 |
| 蛋白别名（GA733-1、EGP-1、CPE-R、WBSCR8） | GEO `gene_assignment` 与试剂盒用别名 | 别名 → NCBI ID → 官方符号。 |
| 人源化 *TACSTD2* 敲入 | BAM 里可能是人 *TACSTD2* 读段 | KI 等位单独计数，不与小鼠 *Tacstd2* 平均。 |
| 抗体不交叉反应 | hRS7 / sacituzumab 常需人 TROP2 | 小鼠 *Tacstd2* 转录 ≠ 人源 ADC 的靶点占用。 |

#### 1.4 跨物种可以转移什么

转移**秩次与程序**，不转移原始计数：

- 1:1 映射后的排序列表（DESeq2 `stat`，或带符号的 −log10 *p* × sign(log2FC)）。
- 自定义基因集：TROP2/claudin TJ 屏障；IFN-γ / CXCL9–11；细胞毒 T；若小鼠设计含 C5aR1 阻断则加髓系 *C5ar1*。
- 对 *TACSTD2* 或 *CLDN4* 的单基因结论，必须先通过 1.2 的 ID 断言。

不要转移 TPM 矩阵、voom 对象，或把“小鼠 top 50 基因”当成未命名的人签名。

---

## 2. Batch / 批次

### English

Treat **species as a biological barrier**, not a batch factor you can estimate away.

#### 2.1 Never

- `cbind` / `rbind` mouse and human count matrices, then ComBat, limma `removeBatchEffect`, Harmony, or a shared `~ species + arm` DESeq2 design.
- Normalize mouse and human together (shared size factors, shared TMM, shared voom).
- Call a gene “replicated across species” because it sits on the same PCA after such a merge.
- Use symbol case-folding as the join key and then “correct batch”.

Species differences include gene content, 3′ bias vs intronless genes, immune-cell marker orthologs, and tumor-cell fraction. Those are not additive technical offsets.

#### 2.2 Within one GSE239485-style series

GSE239485 is 24 LLC tumors, one series, Illumina stranded mRNA, NextSeq 2000 (`GPL30172`). Still:

1. PCA / MDS on VST or log-CPM **before** DE. Color by `arm`, label outliers.
2. Check library size, mapped fraction, rRNA, and sex if the metadata allow.
3. If a hidden flow-cell or harvest-day split appears **and** is not collinear with `arm`, add `+ batch`. If it **is** collinear with `arm`, you cannot adjust; report the confound.
4. Do not add `response` — the series does not have it.

#### 2.3 Within species, multiple mouse series

When stacking LLC + 4T1 + MC38 + a GEMM:

| Factor | Treat as |
| --- | --- |
| GEO series / library kit / center | technical batch |
| Cell line or GEMM allele | **biology** (do not ComBat it away if the claim is about that model) |
| Host strain (C57BL/6 vs BALB/c) | biology + MHC |
| Implant site, harvest day, sex | covariate or batch, depending on balance |

Preferred: **analyze each series**, then meta-analyze ranks (Fisher / Stouffer on 1:1-mapped genes, or fixed gene-set scores). Second choice: one species-internal model `~ arm + series` only when the contrast is shared and series is not collinear with arm.

#### 2.4 Human ICI cohorts

R vs NR labels travel with **center, panel, RNA kit, and biopsy timing**. Model `~ response + cohort` only inside a pre-specified joint human analysis. Never use the mouse series as another “cohort” level.

#### 2.5 Allowed cross-species integration

Do these **after** separate QC and DE:

1. Ortholog-map the ranked lists → FGSEA of species A’s signature on species B’s ranks (matched contrast class; §4).
2. RRHO / rank-rank hypergeometric on the same class.
3. GSVA or ssGSEA of a frozen set in each species; compare **score–phenotype associations** (effect size, CI), not concatenated scores in one linear model without a species term that you then ignore.
4. Transfer learning on **program scores**, not on gene-level counts.

Record the map version, the contrast IDs, and the join key in every figure.

### 中文

把**物种当作生物学屏障**，而不是可以估计掉的批次因子。

#### 2.1 禁止

- 把小鼠与人的 count 矩阵 `cbind` / `rbind` 后做 ComBat、limma `removeBatchEffect`、Harmony，或共享 `~ species + arm` 的 DESeq2。
- 把小鼠与人放在同一套 size factor / TMM / voom 里归一化。
- 因为合并后 PCA 上“在一起”就称基因“跨物种重复”。
- 用符号大小写折叠当连接键再“校正批次”。

物种差异包括基因含量、3′ 偏好 vs 无内含子基因、免疫标记直系同源、肿瘤细胞比例。这些不是可加的技术偏移。

#### 2.2 单个 GSE239485 型系列内部

GSE239485：24 个 LLC 肿瘤、一个 series、Illumina stranded mRNA、NextSeq 2000（`GPL30172`）。仍然要：

1. 在 DE **之前**对 VST 或 log-CPM 做 PCA / MDS。按 `arm` 着色，标出离群点。
2. 查文库大小、比对率、rRNA；若元数据允许则查性别。
3. 若出现隐藏的 flow-cell / 取材日分裂，且**不与** `arm` 共线，可加 `+ batch`；若**与** `arm` 共线，不能校正，只能报告混杂。
4. 不要加 `response` — 该系列没有这个变量。

#### 2.3 同物种、多个小鼠系列

堆叠 LLC + 4T1 + MC38 + GEMM 时：

| 因子 | 当作 |
| --- | --- |
| GEO series / 建库试剂盒 / 中心 | 技术批次 |
| 细胞系或 GEMM 等位基因 | **生物学**（若结论针对该模型，禁止 ComBat 掉） |
| 宿主品系（C57BL/6 vs BALB/c） | 生物学 + MHC |
| 接种部位、取材日、性别 | 协变量或批次，取决于是否平衡 |

首选：**逐系列分析**，再对秩次做荟萃（1:1 基因上的 Fisher / Stouffer，或固定基因集分数）。次选：仅当对照共享且 series 不与 arm 共线时，用种内 `~ arm + series`。

#### 2.4 人类 ICI 队列

R vs NR 标签与**中心、芯片/panel、RNA 试剂盒、活检时点**绑定。只有在预先指定的人类联合分析里才建模 `~ response + cohort`。禁止把小鼠系列当成又一个 `cohort` 水平。

#### 2.5 允许的跨物种整合

在各自 QC 与 DE **之后**：

1. 直系同源映射排序列表 → 用 A 的签名对 B 的秩做 FGSEA（对照类型必须匹配；见 §4）。
2. 同一对照类型上的 RRHO / rank-rank。
3. 在各种属内对冻结基因集做 GSVA / ssGSEA；比较**分数–表型关联**（效应量、CI），而不是把分数拼进一个随后忽略物种项的线性模型。
4. 在**程序分数**上做迁移，不在基因级 count 上做。

每张图写明映射版本、对照 ID、连接键。

---

## 3. Syngeneic vs GEMM / 同源移植 vs 遗传工程小鼠

### English

Pick the model class to match the **claim**, not the GEO accession that was easiest to download.

#### 3.1 Syngeneic (GSE239485 is this class)

Transplantable lines in an immunocompetent, MHC-matched host (LLC on C57BL/6; 4T1 / EMT6 on BALB/c; MC38, CT26, B16).

| Property | Consequence for Tacstd2/Cldn4 × ICI |
| --- | --- |
| Fast, balanced arms, n ≈ 6–8 | Adequate for **treatment-arm DE** (the GSE239485 job). |
| Culture-adapted, clonal | Barrier / TJ programs can be plastic vs the human tumor they “represent”. |
| Ectopic (usually s.c.) | Architecture and immune exclusion are not autochthonous. |
| Known ICI “hot/cold” spectrum | LLC is typically **ICI-cold**; a PD signature here is not a typical human NSCLC response. |
| Host immune system is mouse | Fine for mechanism; not a human HLA / IO-experienced TME. |
| Human TROP2 ADC | Usually needs a **humanized TACSTD2** knock-in line; mouse *Tacstd2* mRNA will not tell you occupancy. |

Use syngeneic when the question is: *what does this regimen change in a controlled TME?* or *does Tacstd2 loss change ICI sensitivity in this line?* (e.g. 4T1 Trop2 WT vs KO).

#### 3.2 GEMM

Autochthonous, genetically driven tumors (KPC, KP lung, MMTV-PyMT, etc.).

| Property | Consequence |
| --- | --- |
| Stepwise evolution, native stroma | Better for **barrier-mediated exclusion** claims (TROP2/claudin TJ). |
| Variable burden and start time | Noisy treatment-arm DE; underpowered if you treat it like LLC n=8. |
| Often lower ICI response | Do not force an R/NR split with n=4. |
| Still a mouse immune system | Closer tissue, not a human trial. |

Use GEMM when the question is: *does the TROP2/claudin program exist in an autochthonous epithelium and track T-cell exclusion?*

#### 3.3 What each class can support vs human

| Claim | Syngeneic | GEMM | Human |
| --- | --- | --- | --- |
| Regimen X changes IFN / myeloid / TJ transcripts | Yes (powered) | Possible, smaller n | Only if on-tx vs pre / untreated exists |
| Tacstd2/Cldn4 track immune exclusion | Possible; confirm histology | Stronger face validity | Spatial / IHC / deconvolution |
| High TACSTD2 predicts ICI non-response | Not from treatment-arm DE | Not from treatment-arm DE | Pre-tx R vs NR or survival |
| ADC (hRS7/SG) + PD-1 synergy | Only with humanized TROP2 | Rare | Clinical / translational cohorts |

Write `model_class` on every figure: `syngeneic` | `GEMM` | `human_tumor`.

### 中文

按**结论类型**选模型，不要按最好下的 GEO 号选模型。

#### 3.1 同源移植（GSE239485 属于此类）

免疫健全、MHC 匹配宿主中的可移植细胞系（C57BL/6 上的 LLC；BALB/c 上的 4T1 / EMT6；以及 MC38、CT26、B16）。

| 性质 | 对 Tacstd2/Cldn4 × ICI 的含义 |
| --- | --- |
| 快、臂平衡、n ≈ 6–8 | 足以做**处理臂 DE**（GSE239485 的本职）。 |
| 培养适应、克隆性 | 屏障 / TJ 程序相对其“代表”的人类肿瘤可能已可塑性改变。 |
| 多为异位（皮下） | 结构与免疫排斥不是原位发生的。 |
| 已知 ICI“热/冷”谱 | LLC 通常 **ICI 冷**；此处的 PD 签名 ≠ 典型人 NSCLC 响应。 |
| 宿主免疫是小鼠 | 可做机制；不是人 HLA / 经 IO 治疗的 TME。 |
| 人 TROP2 ADC | 通常需要**人源化 TACSTD2** 敲入株；小鼠 *Tacstd2* mRNA 不能说明占用。 |

问题是“该方案在受控 TME 中改变了什么”或“该株中 Tacstd2 缺失是否改变 ICI 敏感性”时，用同源移植。

#### 3.2 GEMM

原位、基因驱动肿瘤（KPC、KP 肺、MMTV-PyMT 等）。

| 性质 | 含义 |
| --- | --- |
| 逐步演化、原位基质 | 更适合**屏障介导免疫排斥**（TROP2/claudin TJ）的主张。 |
| 负荷与起始时间不一 | 处理臂 DE 噪声大；不能按 LLC n=8 的方式当充足样本。 |
| ICI 响应常更低 | 禁止用 n=4 强行切 R/NR。 |
| 免疫仍是小鼠 | 组织更接近，仍不是人体试验。 |

问题是“TROP2/claudin 程序是否存在于原位上皮并与 T 细胞排斥共变”时，用 GEMM。

#### 3.3 各类模型相对人类能支撑什么

| 主张 | 同源移植 | GEMM | 人 |
| --- | --- | --- | --- |
| 方案 X 改变 IFN / 髓系 / TJ 转录 | 能（有功效） | 可能，n 更小 | 仅当存在治疗中 vs 治疗前/未治 |
| Tacstd2/Cldn4 与免疫排斥共变 | 可能；需组织学确认 | 表面效度更强 | 空间 / IHC / 反卷积 |
| 高 TACSTD2 预测 ICI 无响应 | 不能从处理臂 DE 得出 | 不能从处理臂 DE 得出 | 治疗前 R vs NR 或生存 |
| ADC（hRS7/SG）+ PD-1 协同 | 仅人源化 TROP2 | 少见 | 临床 / 转化队列 |

每张图写 `model_class`：`syngeneic` | `GEMM` | `human_tumor`。

---

## 4. Treatment-arm DE vs human R/NR — you usually cannot compare them directly

### English

This is the failure mode that produces false “cross-species validation”.

#### 4.1 The two contrasts are different scientific objects

| | Mouse **treatment-arm DE** (GSE239485-style) | Human **R vs NR** |
| --- | --- | --- |
| What is assigned | Treatment (vehicle / doublet / triple) | Outcome (RECIST, irRECIST, PFS cut, investigator call) |
| Typical question | *What did the drug change?* | *Who benefited?* |
| Typical time | Single on-treatment harvest | Often **pre-treatment** biopsy; sometimes on-tx |
| Label source | Cage card | Clinic + imaging + censoring |
| Composition | Treated tumors may shrink and become immune-rich | NR tumors often grow; R tumors may be residual |
| *TACSTD2* meaning | Pharmacodynamic or leftover epithelium | Baseline barrier / exclusion predictor (if pre-tx) |

GSE239485 arms are **C** vehicle (n=8), **D** poly I:C + anti-PD-1 (n=8), **T** poly I:C + anti-PD-1 + anti-C5aR1 (n=8). Those are treatment assignments. There is **no** responder field. Do not invent one from *Cd8a*, IFN, or tumor-size proxies unless the original study scored each mouse and you ingest that table.

#### 4.2 Why a direct gene-level compare is usually invalid

1. **Causal question mismatch.** logFC(treated / vehicle) is a **perturbation**. logFC(R / NR) at baseline is a **selection / state**. Overlap is not replication.
2. **Time-axis mismatch.** Mouse on-tx IFN-high is expected PD. Human pre-tx IFN-high is a common *predictor* of response. Same genes, opposite inferential role. Correlating the two logFC vectors will happily look “significant”.
3. **Direction can invert for compositional reasons.** A shrinking mouse tumor concentrates T-cell RNA; a progressing human NR concentrates tumor RNA. *Tacstd2* / *TACSTD2* (epithelial) will move with purity, not necessarily with pathway logic.
4. **Label mismatch.** RECIST R/NR ≠ “this cage received anti-PD-1”. Mouse volume cutoffs (e.g. 30% shrinkage) are not RECIST.
5. **Model mismatch.** LLC + poly I:C + PD-1 ± C5aR1 is not pembrolizumab monotherapy in human NSCLC or TNBC. Shared *TACSTD2* does not make the contrasts commensurate.
6. **Multiple testing after a bad join.** Case-folded symbols + top-N overlap statistics will report “replication” of housekeeping and interferon genes in almost any ICI pair.

**Rule:** two ranked lists may be compared only when you can fill this line with the same tokens:

```
species A: {pre|on-tx|KO} {treatment_arm|R_vs_NR|high_vs_low|KO_vs_WT}
species B: {pre|on-tx|KO} {treatment_arm|R_vs_NR|high_vs_low|KO_vs_WT}
```

If the second field differs, **do not** report gene-level concordance as validation. If the first field differs, say so in the caption and downgrade to “hypothesis-generating program overlap”.

#### 4.3 What to do instead (allowed hand-offs)

**A. Match the contrast class, then transfer a program.**

- Mouse `D_vs_C` or `T_vs_D` (on-tx treatment DE) → human **on-treatment vs pre** or **on-tx treated vs untreated**, same IO class if possible. Test FGSEA / RRHO on 1:1 orthologs.
- Mouse **vehicle-only** high vs low *Tacstd2* (or *Cldn4*) → human **pre-tx** high vs low *TACSTD2*, and **separately** human pre-tx R vs NR. The second test is the clinical claim; the first is the expression claim. Do not merge them.
- Mouse *Tacstd2* KO vs WT (perturbation; e.g. 4T1 Trop2 studies) → human *TACSTD2*-high vs low (state). This is still not R/NR, but it is a matched “TROP2 activity” axis.

**B. Score, don’t splice matrices.**

Compute a frozen TJ / TROP2-claudin score (include *CLDN7*) and a cytotoxic / IFN score in each species. Report: score vs T-cell estimates; score vs ICI outcome **in human only**; score vs treatment arm **in mouse only**. A three-line forest plot of associations is a valid cross-species figure. A Venn of DE genes from `T_vs_C` ∩ human R vs NR is not.

**C. If the mouse study actually scored responders.**

Some syngeneic ICI papers label mice by volume. That contrast is `mouse_R_vs_NR` (often on-tx). It is closer to human R vs NR but still differs in time, label, and species. Compare it to human **on-tx** R vs NR first; to pre-tx R vs NR only as a secondary, captioned as unmatched time.

**D. Language that is allowed vs forbidden**

| Forbidden | Allowed |
| --- | --- |
| “Mouse DE genes validate human ICI biomarkers.” | “On-treatment mouse programs were tested for enrichment in a time-matched human contrast.” |
| “Tacstd2 down after combo ⇒ TACSTD2-low patients respond.” | “Combo lowered *Tacstd2* in LLC; whether baseline *TACSTD2* predicts human R/NR is a separate pre-tx test.” |
| “Treated mice = responders.” | “Arm = assigned regimen; outcome was not scored.” |
| ComBat(mouse+human) “removes species.” | Separate DE + ortholog FGSEA. |

#### 4.4 Special note on *TACSTD2* / *CLDN4* in ICI

A TROP2/claudin **barrier** hypothesis predicts:

- High *TACSTD2* / tight-junction score ↔ lower T-cell infiltration ↔ **worse** ICI outcome in **human pre-tx** data.
- *Tacstd2* loss or TROP2-directed antibody ↔ disrupted TJ, more T cells, **better** PD-1 activity in **mouse intervention** data.

Those two sentences can both be true and still **must not** be proven by correlating mouse `treated_vs_vehicle` logFC with human `R_vs_NR` logFC. The mouse sentence is a perturbation; the human sentence is a baseline association. Test them as two pre-specified analyses.

### 中文

这是产生虚假“跨物种验证”的主要失败模式。

#### 4.1 两种对照是不同的科学对象

| | 小鼠**处理臂 DE**（GSE239485 型） | 人 **R vs NR** |
| --- | --- | --- |
| 被分配的是什么 | 处理（溶媒 / 双药 / 三药） | 结局（RECIST、irRECIST、PFS 切割、研究者判定） |
| 典型问题 | *药改变了什么？* | *谁获益了？* |
| 典型时点 | 单次治疗中取材 | 常为**治疗前**活检；有时为治疗中 |
| 标签来源 | 笼卡 | 临床 + 影像 + 删失 |
| 组成 | 治疗肿瘤可能缩小、免疫 RNA 变富 | NR 常增大；R 可能是残留灶 |
| *TACSTD2* 含义 | 药效学或残留上皮 | 基线屏障 / 排斥预测（若为治疗前） |

GSE239485 的臂是 **C** 溶媒（n=8）、**D** poly I:C + anti-PD-1（n=8）、**T** poly I:C + anti-PD-1 + anti-C5aR1（n=8）。这是处理分配，**没有**响应者字段。除非原研究给每只鼠打了分且你读入了该表，禁止用 *Cd8a*、IFN 或瘤体积代理去发明 R/NR。

#### 4.2 为何基因级直接比较通常无效

1. **因果问题不匹配。** logFC(处理/溶媒) 是**扰动**。基线 logFC(R/NR) 是**选择/状态**。重叠不是重复。
2. **时间轴不匹配。** 小鼠治疗中 IFN 高是预期 PD。人治疗前 IFN 高常是响应**预测**。同一套基因，推断角色相反。硬相关两条 logFC 向量会显得“显著”。
3. **组成效应可让方向反转。** 缩小的小鼠肿瘤浓缩 T 细胞 RNA；进展的人 NR 浓缩肿瘤 RNA。上皮基因 *Tacstd2* / *TACSTD2* 会随纯度移动，不一定随通路逻辑移动。
4. **标签不匹配。** RECIST R/NR ≠ “这笼打了 anti-PD-1”。小鼠体积阈值（如缩小 30%）不是 RECIST。
5. **模型不匹配。** LLC + poly I:C + PD-1 ± C5aR1 不是人 NSCLC / TNBC 的 PD-1 单药。共享 *TACSTD2* 不能让对照变得可通约。
6. **错误连接后的多重检验。** 大小写折叠 + top-N 重叠，几乎任意一对 ICI 列表都会“重复”到管家基因与干扰素基因。

**规则：** 只有当你能用**同一组记号**填完下面一行时，才可以比较两条排序列表：

```
物种A: {治疗前|治疗中|KO} {处理臂|R_vs_NR|高_vs_低|KO_vs_WT}
物种B: {治疗前|治疗中|KO} {处理臂|R_vs_NR|高_vs_低|KO_vs_WT}
```

若第二项不同，**不要**把基因级一致说成验证。若第一项不同，必须写进图注，并降级为“产生假说的程序重叠”。

#### 4.3 应该怎么做（允许的交接）

**A. 先匹配对照类型，再转移程序。**

- 小鼠 `D_vs_C` 或 `T_vs_D`（治疗中处理 DE）→ 人**治疗中 vs 治疗前**或**治疗中已治 vs 未治**，IO 类别尽量相同。在 1:1 直系同源上做 FGSEA / RRHO。
- 小鼠**仅溶媒**的 *Tacstd2*（或 *Cldn4*）高 vs 低 → 人**治疗前** *TACSTD2* 高 vs 低，并**另外**做人治疗前 R vs NR。第二个才是临床主张，第一个是表达主张。不要合并。
- 小鼠 *Tacstd2* KO vs WT（扰动；如 4T1 Trop2）→ 人 *TACSTD2* 高 vs 低（状态）。仍不是 R/NR，但是匹配的“TROP2 活性”轴。

**B. 评分，不要拼接矩阵。**

在各种属内计算冻结的 TJ / TROP2-claudin 分数（含 *CLDN7*）和细胞毒 / IFN 分数。报告：分数 vs T 细胞估计；分数 vs ICI 结局（**仅人**）；分数 vs 处理臂（**仅小鼠**）。三行关联的森林图是合法的跨物种图。`T_vs_C` ∩ 人 R vs NR 的 DE 基因韦恩图不合法。

**C. 若小鼠研究确实打了响应分。**

部分同源移植 ICI 按体积标小鼠。该对照是 `mouse_R_vs_NR`（常为治疗中）。更接近人 R vs NR，但时点、标签、物种仍不同。先与人**治疗中** R vs NR 比；与治疗前 R vs NR 只能作为次级、并注明时点不匹配。

**D. 允许与禁止的表述**

| 禁止 | 允许 |
| --- | --- |
| “小鼠 DE 基因验证了人 ICI 生物标志物。” | “小鼠治疗中程序在时点匹配的人对照中做了富集检验。” |
| “联合治疗后 Tacstd2 下降 ⇒ TACSTD2 低的患者会响应。” | “联合治疗降低了 LLC 中的 *Tacstd2*；基线 *TACSTD2* 是否预测人 R/NR 是另一项治疗前检验。” |
| “处理组小鼠 = 响应者。” | “臂 = 分配的方案；未评分结局。” |
| ComBat(小鼠+人)“去掉物种。” | 分开 DE + 直系同源 FGSEA。 |

#### 4.4 关于 ICI 中 *TACSTD2* / *CLDN4* 的特别说明

TROP2/claudin **屏障**假说预测：

- 高 *TACSTD2* / 紧密连接分数 ↔ 较低 T 细胞浸润 ↔ 在**人治疗前**数据中 ICI 结局**更差**。
- *Tacstd2* 缺失或 TROP2 抗体 ↔ TJ 破坏、更多 T 细胞、在**小鼠干预**中 PD-1 活性**更好**。

这两句可以同时成立，也**绝不能**靠把小鼠 `处理_vs_溶媒` logFC 与人 `R_vs_NR` logFC 相关来证明。小鼠句是扰动，人句是基线关联。应作为两项预先指定的分析分别检验。

---

## 5. Templates for GSE239485-style designs / GSE239485 型设计模板

### English

#### 5.1 What “GSE239485-style” means

Copy this checklist when a GEO series looks like GSE239485 even if the drugs differ.

| Item | GSE239485 fact | Template rule |
| --- | --- | --- |
| Organism | *Mus musculus* | Mouse genome; human only after §1 map |
| Model | LLC syngeneic, tumor bulk | `model_class=syngeneic`; not GEMM |
| Platform | NextSeq 2000, `GPL30172` | Record kit: Illumina Stranded mRNA Prep Ligation |
| Accessions | GSE239485 / SRP451881 / PRJNA999384 | GSM7666074–GSM7666097 (24 samples) |
| Arms | C vehicle; D poly I:C + aPD-1; T D + aC5aR1 | 3-level `arm`; **nested**, not 2×2×2 |
| n | 8 / 8 / 8 | Balanced; still check PCA |
| Outcome | None | `response_label=NA` |
| Science goal of the series | Effect of C5aR1 blockade **on top of** poly I:C + PD-1 | Primary contrast `T_vs_D` |
| Authors (GEO) | Ajona, Senent, Pio (Cima, University of Navarra) | Cite GEO + paper when used |

Filled tables: `templates/gse239485_style_metadata.tsv`, `templates/contrasts.yaml`.  
DE skeleton: `templates/deseq2_gse239485_style.R`.

#### 5.2 Metadata columns you must create

Minimum: `sample_id`, `geo_gsm`, `organism`, `model_class`, `model_name`, `host_strain`, `arm`, `timepoint`, `response_label`, `batch`.

For nested combos, also store **explicit 0/1 factors** (`factor_polyIC`, `factor_anti_PD1`, `factor_anti_C5aR1`) so nobody later fits a full factorial. In GSE239485, `factor_anti_C5aR1` is 1 only when both poly I:C and anti-PD-1 are 1. That empty-cell structure is the design.

`response_label` stays `NA` unless a **per-mouse outcome table** exists. Expression-derived “responders” are forbidden.

#### 5.3 Design formulae

**Use:** `~ arm` with levels `vehicle`, `polyIC_aPD1`, `polyIC_aPD1_aC5aR1`.

**Equivalent nested:** `~ combo + add_on` where `combo` is none vs polyIC+aPD1 and `add_on` is none vs aC5aR1 (identifiable only because add-on never appears alone).

**Do not use:** `~ polyIC * aPD1 * aC5aR1` (empty cells).  
**Do not use:** `~ response`.  
**Do not use:** a shared mouse+human design.

Contrasts to emit, in this order:

1. `D_vs_C` — combo vs vehicle (treatment-arm DE).
2. `T_vs_D` — add-on vs combo (primary scientific contrast of this series).
3. `T_vs_C` — full regimen vs vehicle (interpret only after 1 and 2).

Each result table gets `contrast_type=treatment_arm_DE` or `add_on_treatment_arm_DE`. Downstream human tests must read that field.

#### 5.4 QC and DE defaults

1. Count on GRCm39; keep NCBI/Ensembl IDs in the row names or a parallel column.
2. Filter low counts **within this series** (e.g. ≥10 counts in ≥8 samples — one full arm).
3. DESeq2 (or edgeR/QL) with `~ arm`; apeglm or ashr shrink for viz, **unshrunk or `stat`** for ranking / FGSEA.
4. Pre-specify focal genes: *Tacstd2*, *Cldn4*, *Cldn7*, *Epcam* (paralog control), plus IFN and cytotoxic anchors.
5. Pathways: MSigDB Hallmark / GO TJ and IFN; custom TROP2-claudin set. Convert to human only with the §1 table.
6. Deconvolution (CIBERSORTx mouse, mMCPcounter, or seqNMF) as a **composition check**, especially if *Tacstd2* moves opposite to *Cd8a*.

#### 5.5 Cross-species hand-off from this design (worked)

You have three mouse ranked lists, all **on-treatment, treatment-arm**. Allowed next steps:

| Mouse object | Human object you may test | Human object you may not call “validation” |
| --- | --- | --- |
| `T_vs_D` or `D_vs_C` ranks | On-tx vs pre (or treated vs untreated) ICI RNA-seq | Pre-tx R vs NR gene-level overlap |
| Vehicle-only *Tacstd2* or TJ score | Pre-tx *TACSTD2* or TJ score vs TIL / R vs NR | Using treated arms to define “high TROP2” |
| TJ score Δ (T − D) | On-tx change in a TJ or IFN score | Saying Δ equals a response biomarker |
| Nothing in GSE239485 | — | Any plot titled “mouse vs human ICI DE” that hides contrast class |

If you need a mouse **perturbation** signature for TROP2 itself, **do not** extract it from GSE239485 (no Tacstd2 allele). Use a Trop2 WT vs KO series (e.g. 4T1 RNA-seq such as GSE334497) and keep that contrast ID separate.

#### 5.6 Figure and caption contract

Every cross-species figure must state, in both languages if the paper is bilingual:

- species + `model_class`
- `contrast_type` and `timepoint` on **each** side
- gene ID space and ortholog map version
- that GSE239485-style arms are assigned treatments, not R/NR

Refuse captions of the form “conserved ICI DE genes (*Tacstd2*, *Cldn4*)” unless both sides are treatment-arm DE at a matched time.

#### 5.7 Operator checklist (print this)

- [ ] IDs: *Tacstd2* 56753 ↔ *TACSTD2* 4070; *Cldn4* 12740 ↔ *CLDN4* 1364; *Cldn7* listed; *Epcam* not merged.
- [ ] Mouse and human counted and DE’d **separately**; no cross-species ComBat.
- [ ] `model_class` recorded (syngeneic LLC ≠ GEMM ≠ human).
- [ ] `response_label` is NA for GSE239485-style series.
- [ ] Primary mouse contrasts are `D_vs_C`, `T_vs_D`, `T_vs_C` — all treatment-arm.
- [ ] Any human R/NR test is a **new**, pre-specified analysis (usually pre-tx), not a Venn with `T_vs_C`.
- [ ] Caption tokens match on contrast class or the figure is labeled unmatched.

### 中文

#### 5.1 “GSE239485 型”指什么

只要 GEO 系列看起来像 GSE239485（即使药物不同），就套用本清单。

| 项 | GSE239485 事实 | 模板规则 |
| --- | --- | --- |
| 物种 | *Mus musculus* | 小鼠基因组；人只在 §1 映射之后出现 |
| 模型 | LLC 同源移植，肿瘤 bulk | `model_class=syngeneic`；不是 GEMM |
| 平台 | NextSeq 2000，`GPL30172` | 记录试剂盒：Illumina Stranded mRNA Prep Ligation |
| 登录号 | GSE239485 / SRP451881 / PRJNA999384 | GSM7666074–GSM7666097（24 例） |
| 臂 | C 溶媒；D poly I:C + aPD-1；T = D + aC5aR1 | 三水平 `arm`；**嵌套**，不是 2×2×2 |
| n | 8 / 8 / 8 | 平衡；仍要看 PCA |
| 结局 | 无 | `response_label=NA` |
| 系列科学目标 | 在 poly I:C + PD-1 **之上**看 C5aR1 阻断 | 主对照 `T_vs_D` |
| GEO 作者 | Ajona, Senent, Pio（Navarra Cima） | 使用时同时引 GEO 与论文 |

填好的表：`templates/gse239485_style_metadata.tsv`、`templates/contrasts.yaml`。  
DE 骨架：`templates/deseq2_gse239485_style.R`。

#### 5.2 必须建立的元数据列

最低：`sample_id`、`geo_gsm`、`organism`、`model_class`、`model_name`、`host_strain`、`arm`、`timepoint`、`response_label`、`batch`。

嵌套组合还要存**显式 0/1 因子**（`factor_polyIC`、`factor_anti_PD1`、`factor_anti_C5aR1`），以免后人拟满因子设计。在 GSE239485 中，仅当 poly I:C 与 anti-PD-1 均为 1 时 `factor_anti_C5aR1` 才为 1。空单元格就是该设计。

除非存在**逐鼠结局表**，`response_label` 保持 `NA`。禁止用表达定义“响应者”。

#### 5.3 设计公式

**使用：** `~ arm`，水平为 `vehicle`、`polyIC_aPD1`、`polyIC_aPD1_aC5aR1`。

**等价嵌套：** `~ combo + add_on`。`combo` 为无 vs polyIC+aPD1，`add_on` 为无 vs aC5aR1（可识别是因为加药从不单独出现）。

**不要用：** `~ polyIC * aPD1 * aC5aR1`（空单元格）。  
**不要用：** `~ response`。  
**不要用：** 小鼠+人共享设计。

按此顺序输出对照：

1. `D_vs_C` — 联合 vs 溶媒（处理臂 DE）。
2. `T_vs_D` — 加药 vs 联合（本系列的科学主对照）。
3. `T_vs_C` — 全方案 vs 溶媒（必须在 1 和 2 之后解释）。

每张结果表带 `contrast_type=treatment_arm_DE` 或 `add_on_treatment_arm_DE`。下游人类检验必须读该字段。

#### 5.4 QC 与 DE 默认

1. 在 GRCm39 上计数；行名或并行列保留 NCBI/Ensembl ID。
2. **在本系列内**过滤低计数（例如 ≥10 counts 且至少 8 个样本——一整臂）。
3. DESeq2（或 edgeR/QL）用 `~ arm`；作图可用 apeglm / ashr 收缩，**排序 / FGSEA 用未收缩或 `stat`**。
4. 预指定焦点基因：*Tacstd2*、*Cldn4*、*Cldn7*、*Epcam*（旁系对照），外加 IFN 与细胞毒锚点。
5. 通路：MSigDB Hallmark / GO 的 TJ 与 IFN；自定义 TROP2-claudin 集。转人只用 §1 的表。
6. 反卷积（CIBERSORTx 小鼠、mMCPcounter 或 seqNMF）作为**组成检查**，尤其当 *Tacstd2* 与 *Cd8a* 反向移动时。

#### 5.5 从本设计出发的跨物种交接（实例）

你有三份小鼠排序列表，全部是**治疗中、处理臂**。允许的下一步：

| 小鼠对象 | 可以检验的人对象 | 不能称为“验证”的人对象 |
| --- | --- | --- |
| `T_vs_D` 或 `D_vs_C` 秩 | 治疗中 vs 治疗前（或已治 vs 未治）ICI RNA-seq | 治疗前 R vs NR 的基因级重叠 |
| 仅溶媒的 *Tacstd2* 或 TJ 分数 | 治疗前 *TACSTD2* 或 TJ 分数 vs TIL / R vs NR | 用处理臂来定义“高 TROP2” |
| TJ 分数差（T − D） | 治疗中 TJ 或 IFN 分数的变化 | 声称该差值等于疗效生物标志物 |
| GSE239485 中不存在的东西 | — | 任何隐瞒对照类型的“小鼠 vs 人 ICI DE”图 |

若需要 TROP2 本身的小鼠**扰动**签名，**不要**从 GSE239485 抽（无 Tacstd2 等位基因）。改用 Trop2 WT vs KO 系列（如 4T1 的 GSE334497），并分开保留对照 ID。

#### 5.6 图与图注契约

每张跨物种图必须写明（若论文双语则两种语言都写）：

- 物种 + `model_class`
- **两侧**的 `contrast_type` 与 `timepoint`
- 基因 ID 空间与直系同源映射版本
- GSE239485 型臂是分配的处理，不是 R/NR

除非两侧都是时点匹配的处理臂 DE，否则拒绝“保守的 ICI DE 基因（*Tacstd2*、*Cldn4*）”这类图注。

#### 5.7 操作清单（打印）

- [ ] ID：*Tacstd2* 56753 ↔ *TACSTD2* 4070；*Cldn4* 12740 ↔ *CLDN4* 1364；已列 *Cldn7*；未与 *Epcam* 合并。
- [ ] 小鼠与人**分开**计数与 DE；无跨物种 ComBat。
- [ ] 已记录 `model_class`（同源移植 LLC ≠ GEMM ≠ 人）。
- [ ] GSE239485 型系列的 `response_label` 为 NA。
- [ ] 小鼠主对照为 `D_vs_C`、`T_vs_D`、`T_vs_C` — 全是处理臂。
- [ ] 任何人 R/NR 检验都是**新的**、预先指定的分析（多为治疗前），不是与 `T_vs_C` 的韦恩图。
- [ ] 图注记号在对照类型上匹配，否则标明不匹配。

---

## 6. File map / 文件对照

| Path | Role |
| --- | --- |
| `methods/mouse_human/playbook.md` | This document (EN + 中文) |
| `methods/mouse_human/templates/ortholog_focal_genes.tsv` | Frozen ID table for TROP2 / CLDN4 / CLDN7 / EPCAM |
| `methods/mouse_human/templates/gse239485_style_metadata.tsv` | 24-row worked metadata (GSE239485) |
| `methods/mouse_human/templates/contrasts.yaml` | Allowed / forbidden contrasts and hand-off rules |
| `methods/mouse_human/templates/deseq2_gse239485_style.R` | `~ arm` DESeq2 skeleton; no R/NR |

| 路径 | 作用 |
| --- | --- |
| `methods/mouse_human/playbook.md` | 本文（中英） |
| `methods/mouse_human/templates/ortholog_focal_genes.tsv` | TROP2 / CLDN4 / CLDN7 / EPCAM 冻结 ID 表 |
| `methods/mouse_human/templates/gse239485_style_metadata.tsv` | GSE239485 的 24 行实例元数据 |
| `methods/mouse_human/templates/contrasts.yaml` | 允许 / 禁止的对照与交接规则 |
| `methods/mouse_human/templates/deseq2_gse239485_style.R` | `~ arm` DESeq2 骨架；无 R/NR |

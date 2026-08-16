# Tacstd2 (TROP2) 与 Cldn4 在小鼠肺癌 ICI 队列中的表达分析
# Tacstd2 (TROP2) and Cldn4 across mouse lung ICI cohorts

> 本文所有统计量均由 `notes/mouse/data/` 中的真实下载矩阵计算得到，脚本见 `scripts/mouse/`，
> 结果表见 `results/mouse/*.csv`。未编造任何数据或登录号。
> All statistics below are computed from the real processed matrices under `notes/mouse/data/`
> (scripts in `scripts/mouse/`, result tables in `results/mouse/*.csv`). No values or accessions were fabricated.

---

## 1. 摘要 / Executive summary

- **目标基因 / Target genes:** `Tacstd2` (TROP2, ENSMUSG00000051397)、`Cldn4` (ENSMUSG00000047501)。
- **范围 / Scope:** 仅小鼠 (Mus musculus) 肺肿瘤、免疫检查点抑制 (ICI) 相关数据集；仅下载处理后的矩阵；跳过 FASTQ / 原始质谱及 >2 GB 文件。
- **编目 / Catalog:** 原 11 个登录号 + 本轮补充的 GSE309199 / GSE330941（见 `notes/mouse/catalog.md`）。**GSE76628 已核实并排除**（胃/侧腹 Ad-VEGF 基质，非肺、非 ICI）。
- **对 ICB 应答（responder vs non-responder）的诚实结论 / Honest ICB-response finding:**
  **未找到任何带有逐鼠 ICB 应答/不应答标签的公开小鼠肺癌 RNA 处理矩阵。** 公开的小鼠 R vs NR RNA 存在于 CT26 结肠、SCC、MC38、黑色素瘤，**不是肺**。检索记录见 `notes/mouse/icb_response_search.md`。
  **No public mouse *lung* ICI RNA matrix with per-mouse ICB responder vs non-responder labels was found.** Public mouse R vs NR RNA exists for CT26 colon, SCC, MC38, and melanoma — not lung.
- **核心结论（治疗 ≠ 应答） / Key finding (treatment ≠ response):**
  - **Tacstd2/TROP2 在以抗 PD-1/PD-L1 为基础的治疗后一致性上调**，并在 bulk 肿瘤中与"T 细胞炎症/细胞毒"免疫评分**正相关**。这是**治疗/免疫表型**关联，**不能**外推为 ICB 临床应答。
    Tacstd2/TROP2 is **up-regulated after anti-PD-1/PD-L1–based therapy** and **positively correlates** with a cytotoxic/T-cell-inflamed immune score in bulk tumors. This is a **treatment / immune-phenotype** association, **not** ICB response.
  - **Cldn4 无一致上调**，在多个 ICI 联合臂中反而趋于下降，且与免疫评分无正相关。
    Cldn4 shows **no consistent induction**, trends downward in several ICI-combination arms, and is not positively associated with the immune score.

---

## 2. 方法 / Methods

- **表达标准化 / Normalization:** E-MTAB-13704 counts→log2 CPM；GSE239485 使用作者归一化 log2 值；GSE297630 使用 log2 RMA 每样本信号；GSE330658 / GSE197260 使用 log2(TPM+1)；scRNA 使用每细胞 log1p(CP10K) 与 pseudobulk log2 CPM。
- **分组比较 / Group comparisons:** 每个处理臂对对照计算 log2 fold change、Welch t 检验、Mann-Whitney U 检验；同一基因内多重比较用 BH-FDR 校正 (`treated_vs_control_stats.csv`)。
- **免疫结局轴 / Immune-outcome axis:** 以细胞毒/T 细胞炎症模块 (`Cd8a, Cd8b1, Gzmb, Gzmk, Prf1, Ifng, Nkg7, Cxcl9, Cxcl10, Cd3e, Pdcd1, Ptprc`) 的每样本 z 分数均值作为"免疫评分"，与目标基因做 Spearman/Pearson 相关 (`immune_correlation.csv`)。
- **基因 ID / Gene IDs:** 通过 Ensembl REST 校验（`notes/mouse/gene_map.json`）；大小写不敏感匹配以兼容 GSE297630 的大写符号。
- **单细胞 / Single cell:** 稀疏 mtx 流式解析，按每细胞 UMI 过滤（GSE129297 未过滤矩阵阈值 1000，其余 500），计算 pseudobulk、表达阳性细胞比例与每细胞 Spearman 相关 (`scrna_pseudobulk.csv`, `scrna_singlecell_immune_corr.csv`)。

---

## 3. 核心（可重复）队列结果 / Core replicated-cohort results

### 3.1 处理 vs 对照 — Tacstd2 / Treated vs control — Tacstd2

| Dataset | Comparison | n (trt/ctrl) | log2FC | Welch p | MWU p | BH-FDR |
|---|---|---|---|---|---|---|
| GSE239485 | PolyIC+anti-PD-1 vs vehicle | 8/8 | **+2.093** | 6.3e-4 | 6.2e-4 | **9.7e-4** |
| GSE239485 | PolyIC+anti-PD-1+anti-C5aR1 vs vehicle | 8/8 | **+2.576** | 9.7e-4 | 1.9e-3 | **9.7e-4** |
| GSE297630 | anti-PD-1 tolerant vs control | 3/3 | **+2.050** | **4.2e-5** | 0.100* | **4.2e-5** |
| E-MTAB-13704 | VEGFRi/aPD-L1 vs vehicle | 4/5 | **+0.414** | 7.6e-4 | 0.016 | **3.8e-3** |
| E-MTAB-13704 | aPD-L1 mono vs vehicle | 5/5 | -0.123 | 0.143 | 0.151 | 0.357 |
| E-MTAB-13704 | ATRi/aPD-L1 vs vehicle | 5/5 | +0.036 | 0.530 | 0.548 | 0.662 |
| E-MTAB-13704 | Cisplatin/aPD-L1/aCTLA4 vs vehicle | 4/5 | -0.009 | 0.918 | 0.556 | 0.918 |

\* n=3/组时 Mann-Whitney U 的最小可能双侧 p 即为 0.1（样本量限制，而非无效应）。/ With n=3 per group the smallest possible two-sided MWU p is 0.1 (a sample-size floor, not absence of effect).

**解读 / Interpretation:** 在 3 个独立、含重复的队列中，Tacstd2 在抗 PD-1/PD-L1 相关处理后显著升高（GSE239485、GSE297630 效应量大，log2FC≈2；E-MTAB-13704 的 VEGFRi/aPD-L1 联合臂 FDR=0.0038）。单用 aPD-L1（E-MTAB-13704）未见升高，提示上调更多出现在**联合治疗 / 存活（tolerant）肿瘤细胞**背景下。

### 3.2 处理 vs 对照 — Cldn4 / Treated vs control — Cldn4

| Dataset | Comparison | n (trt/ctrl) | log2FC | Welch p | BH-FDR |
|---|---|---|---|---|---|
| GSE297630 | anti-PD-1 tolerant vs control | 3/3 | -0.190 | 0.016 | 0.016 |
| E-MTAB-13704 | ATRi vs vehicle | 4/5 | -0.655 | 0.040 | 0.100 |
| E-MTAB-13704 | Cisplatin/aPD-L1/aCTLA4 vs vehicle | 4/5 | -0.504 | 0.038 | 0.100 |
| GSE239485 | PolyIC+anti-PD-1 vs vehicle | 8/8 | -0.771 | 0.238 | 0.357 |
| E-MTAB-13704 | aPD-L1 mono vs vehicle | 5/5 | -0.098 | 0.621 | 0.621 |

**解读 / Interpretation:** Cldn4 未出现与 Tacstd2 类似的一致性上调；在若干 ICI 联合臂（E-MTAB-13704 的 ATRi、Cisplatin/aPD-L1/aCTLA4，以及 GSE297630）中反而**下降**（多为未校正显著、校正后趋势）。因此 Tacstd2 与 Cldn4 在 ICI 背景下**行为相反/解耦**。

### 3.3 与免疫评分的相关性 / Correlation with the immune score

| Dataset | Gene | n | Spearman ρ | p | Pearson r | p |
|---|---|---|---|---|---|---|
| E-MTAB-13704 | Tacstd2 | 27 | **+0.518** | **0.0057** | +0.506 | 0.0070 |
| E-MTAB-13704 | Cldn4 | 27 | -0.279 | 0.159 | -0.264 | 0.183 |
| GSE239485 | Tacstd2 | 24 | **+0.482** | **0.0171** | +0.534 | 0.0072 |
| GSE239485 | Cldn4 | 24 | +0.242 | 0.255 | +0.171 | 0.424 |
| GSE297630 | Tacstd2 | 6 | -0.314 | 0.544 | +0.020 | 0.970 |
| GSE297630 | Cldn4 | 6 | -0.029 | 0.957 | -0.047 | 0.929 |

**解读 / Interpretation:** 在两个样本量较大的 bulk 队列（E-MTAB-13704、GSE239485）中，**Tacstd2 与细胞毒/T 细胞炎症评分显著正相关**（ρ≈0.5，p<0.02），即 TROP2 高的小鼠肺肿瘤更"热"（更富 CD8/细胞毒特征）。Cldn4 无此正相关（甚至在 E-MTAB-13704 中呈负向趋势）。GSE297630 因 n=6 且为经 anti-PD-1 筛选的存活细胞，功效不足，未见显著相关。

---

## 4. 单细胞与 pseudobulk 队列 / Single-cell & pseudobulk cohorts

### 4.1 Pseudobulk（每条件） / Pseudobulk per condition

| Dataset (model) | Condition | Tacstd2 log2CPM (%细胞+) | Cldn4 log2CPM (%细胞+) |
|---|---|---|---|
| GSE129297 (SCLC, RPM) | Ctrl | 1.62 (0.68%) | 3.77 (5.29%) |
| GSE129297 | anti-PD-1 | 1.97 (1.28%) | 3.71 (4.58%) |
| GSE129297 | CDK7i (YKL) | 4.15 (5.70%) | 3.78 (3.68%) |
| GSE129297 | anti-PD-1+CDK7i | 2.60 (2.07%) | 3.51 (3.37%) |
| GSE133604 (KP) | Ctrl | 3.64 (3.42%) | 0.49 (0.08%) |
| GSE133604 | Ctrl+anti-PD-1 | 3.63 (3.16%) | 0.12 (0.03%) |
| GSE133604 | Asf1a-KO | 4.11 (3.93%) | 0.30 (0.03%) |
| GSE133604 | Asf1a-KO+anti-PD-1 | 4.07 (3.72%) | 0.80 (0.07%) |
| GSE297632 (LLC, 全肿瘤) | Control | 0.61 (0.71%) | 0.46 (0.37%) |
| GSE297632 | anti-PD-1 | 0.36 (0.32%) | 0.90 (0.27%) |

**解读 / Interpretation:** scRNA 队列每臂仅 1 个样本，不做组间 p 值。GSE129297 中 Tacstd2 阳性细胞比例随 anti-PD-1（0.68%→1.28%）及 CDK7i 明显上升，方向与 bulk 一致。GSE297632（anti-PD-1 tolerant 的 scRNA 配套）为**全肿瘤**取样、上皮/肿瘤细胞比例低，故 pseudobulk 不同于其配套的纯肿瘤细胞芯片 GSE297630——凸显 bulk 与 scRNA 采样区室的差异。

### 4.2 单细胞层面的区室效应 / Single-cell compartment effect

- 在**全肿瘤 scRNA**（GSE129297、GSE133604、GSE297632）中，Tacstd2/Cldn4 在单细胞层面与免疫模块**负相关或近零**（例如 GSE129297 Cldn4 ρ≈-0.17~-0.19，p<1e-16）。这是**区室互斥**的表现：TROP2/Claudin-4 主要由上皮/肿瘤细胞表达，而免疫模块由免疫细胞表达，二者很少在同一细胞共表达。
  In whole-tumor scRNA, Tacstd2/Cldn4 anti-correlate (or are ~0) with the immune module at single-cell resolution — a **compartment mutual-exclusivity** effect (epithelial/tumor vs immune cells).
- 在 **FACS 分选的免疫细胞** GSE222158 中，Cldn4 几乎不表达（阳性细胞 ≈0.03–0.12%），证实其上皮/肿瘤限定性；Tacstd2 在部分髓系/DC 细胞有表达并与免疫模块弱正相关（如 CD3_Ctrl ρ=+0.135, p=3.6e-21）。
  In FACS-sorted immune cells (GSE222158), Cldn4 is essentially absent, confirming epithelial restriction; Tacstd2 is detectable in some myeloid/DC cells.
- 因此，第 3.3 节 bulk 层面的 **Tacstd2–免疫正相关** 反映的是"TROP2 高的肿瘤同时更 T 细胞炎症"这一**样本层面**关联，而非单细胞共表达。
  Hence the bulk-level positive Tacstd2–immune association is a **tumor-level** (sample-level) relationship, not single-cell co-expression.

---

## 5. 支持性 / 描述性证据 / Supporting & descriptive evidence

- **GSE330658 (Egfr-mut 肺癌, n=2/臂):** Tacstd2/Cldn4 各处理臂对对照差异均不显著（功效不足，见 `treated_vs_control_stats.csv`）。已沉积臂为 PTX / anti-VEGF（anti-PD-L1 臂未提供处理矩阵）。
- **GSE197260 (EGFR-TKI + aPD1/aVEGFR2, n=1/臂):** 仅描述性 log2FC；相对 gef_vehicle_d21，anti-PD-1 (4H2) 臂中 Tacstd2 与 Cldn4 均较低（分别 -5.6、-6.5 log2），但无重复无法检验。
- **GSE241978 (CMT167 LUAD, AhR 敲除):** 属基因层面的检查点通路（PD-L1/IDO1）扰动而非 ICI 药物治疗。DE 表中 Cldn4 fold change=-6.27（线性）、p=0.40（NS）；Tacstd2 未进入过滤后 DE 表。
- **PXD059688 (anti-PD-1 ± 高剂量维生素C, 蛋白质组):** 全部定量质谱文件为 >2 GB 原始数据（.raw/.mgf/.msf），按规则跳过；唯一 <2 GB 的处理文件为单一条件 (AA) 的 93 蛋白鉴定报告，**不含 Tacstd2/Cldn4**，无跨条件蛋白定量表，故**仅编目、无法分析**。
- **GSE76628:** **排除**。GEO 标题为 *Stromal-Based Signatures for the Classification of Gastric Cancer [part II]*；实验为侧腹 Ad-VEGF-A164 血管生成 + 抗 VEGFR，**不是肺肿瘤、不是 PD-1/PD-L1 ICB**。

### 5.1 本轮补充：最接近 ICB 应答的肺 RNA / Closest lung ICI RNA added for the ICB-response question

完整检索见 `notes/mouse/icb_response_search.md`。结果表：`results/mouse/icb_proxy_*.csv`。

| Dataset | Comparison | Gene | n | log2FC | Welch p | MWU p | What this is |
|---|---|---|---|---|---|---|---|
| GSE309199 | aPD-1 vs Ctrl | Tacstd2 | 3/3 | +0.491 | 0.474 | 0.184 | treatment, not R vs NR |
| GSE309199 | aPD-1 vs Ctrl | Cldn4 | 3/3 | +0.254 | 0.606 | 0.700 | treatment, not R vs NR |
| GSE309199 | aPD-1+entinostat vs Ctrl | Tacstd2 | 3/3 | +1.031 | 0.168 | 0.376 | treatment, not R vs NR |
| GSE309199 | aPD-1+entinostat vs Ctrl | Cldn4 | 3/3 | +1.467 | 0.034 | 0.100* | treatment; FDR 0.074 |
| GSE330941 | Ago2KO (ICI-sensitized) vs WT (ICI-refractory) | Tacstd2 | 4/4 | +1.480 | 0.082 | 0.114 | **no ICI on these RNA samples** |
| GSE330941 | Ago2KO vs WT | Cldn4 | 4/4 | +1.352 | **0.013** | 0.029 | model-level sensitivity proxy |

\* n=3 时 MWU 最小双侧 p=0.1。

**解读 / Interpretation:**
- GSE309199 中 Tacstd2 在 aPD-1 后方向与先前队列一致（升高），但 n=3 **不显著**。不能把“治疗后升高”写成“预测 ICB 应答”。
- GSE330941 中 Tacstd2 与 Cldn4 在 **ICI 敏感化（Ago2KO）** 模型中均高于难治 WT（Cldn4 Welch p=0.013）。RNA 取自接种后第 12 天、**未做 ICI**。这只说明 ICI 敏感基因型的基线转录，**不是**治疗后 R vs NR。
- 因此：**不能声称 Tacstd2 或 Cldn4 预测小鼠肺癌 ICB 应答**——公开数据里没有该标签。

---

## 6. 结论 / Conclusions

1. **完成标准 / Done-criteria:** 小鼠肺 ICI 目录已建立；≥2 个处理队列已用真实统计量分析。本轮补做了 **ICB 应答**专项检索。
2. **ICB 应答（诚实）:** **没有**公开的小鼠肺癌 RNA 处理矩阵带有逐鼠 ICB responder vs non-responder 标签。GSE76628 不是肺 ICI，已排除。公开 R/NR RNA 在 CT26 / SCC / MC38 / 黑色素瘤，不在肺。
   **No public mouse lung ICI RNA with per-mouse ICB R vs NR labels.** GSE76628 is gastric/flank VEGF stroma, excluded.
3. **Tacstd2/TROP2** 在抗 PD-1/PD-L1 **治疗**后上调，并标记更 T 细胞炎症的小鼠肺肿瘤。这是治疗/免疫表型，**不是 ICB 应答**。
4. **Cldn4** 在治疗臂中无一致上调；在 GSE330941 的 ICI 敏感化模型中反而更高。不能写成“Cldn4 与 ICB 应答负相关”——没有 R/NR 标签。

## 7. 局限 / Limitations

- 部分队列每臂样本量小（GSE297630 / GSE309199 n=3、GSE330658 n=2、GSE197260/scRNA n=1），限制统计功效。
- 不同平台/归一化方式差异，log2FC 绝对值不可跨数据集直接比较。
- 免疫评分为转录模块代理，非实测细胞比例或临床结局。
- scRNA 的负相关来自区室互斥，勿与 bulk 样本层面正相关混淆。
- **最大缺口 / Main gap:** 缺小鼠肺癌逐鼠 ICB 应答标签；不能用治疗 vs 对照代替应答。

## 8. 产物清单 / Output files

- `results/mouse/cohort_overview.csv` — 各队列样本数与分组 / cohort samples & groups
- `results/mouse/gene_group_summary.csv` — 每数据集每组 Tacstd2/Cldn4 均值±SD / per-group means
- `results/mouse/treated_vs_control_stats.csv` — 处理 vs 对照 log2FC + Welch/MWU + FDR
- `results/mouse/immune_correlation.csv` — 与免疫评分的 Spearman/Pearson
- `results/mouse/scrna_pseudobulk.csv` — scRNA pseudobulk 与阳性细胞比例
- `results/mouse/scrna_singlecell_immune_corr.csv` — 单细胞相关
- `results/mouse/icb_proxy_group_summary.csv` / `icb_proxy_stats.csv` / `icb_proxy_immune_corr.csv` — GSE309199 + GSE330941
- `notes/mouse/icb_response_search.md` — ICB 应答检索（含 GSE76628 排除理由）
- `notes/mouse/catalog.md` / `catalog.json` — 数据目录 / dataset catalog
- `notes/mouse/data/MANIFEST.tsv` — 下载来源/大小/校验 / download provenance
- `scripts/mouse/*.py` — 可重复流程 / reproducible pipeline

# GEO 2025 lung ICI leftovers: TACSTD2 / CLDN4

## English

### Scope

This is a strict **GEO-public-in-2025** follow-up of human lung-cancer ICI
series not in the previously analyzed reference set (GSE126044, GSE135222,
GSE136961, GSE166449, GSE93157, GSE207422, and GSE205335). A series was
eligible for an outcome test only if:

1. the deposited material contained tumor/epithelial cells;
2. both TACSTD2 and CLDN4 were measured; and
3. GEO supplied a patient-level ICI outcome without an inferred join.

The focused audit is in `series_audit.tsv`. Of 13 plausible series checked,
**GSE233203 was the only eligible series**. Exclusions are results, not missing
work: for example, GSE309652 has a clean 72-patient R/NR endpoint but its
768-gene panel measures neither target, while GSE292098 measures both targets
but supplies no per-ROI clinical outcome key.

### GSE233203 analysis

GSE233203 contains pleural-fluid scRNA-seq from seven stage-IV lung
adenocarcinoma patients (3 Response, 4 Non-response) receiving an
atezolizumab-containing combination. Raw 10x counts were summed across all
cells within each patient, converted to log2(CPM + 1), and tested at the
patient—not cell—level. P-values are exact two-sided permutation
Mann-Whitney values. Holm values adjust the two prespecified gene tests.

| Gene | Median Response | Median Non-response | Cliff's delta | Exact p | Holm p |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | 5.743 | 3.702 | 0.833 | 0.114 | 0.229 |
| CLDN4 | 5.543 | 3.184 | 0.500 | 0.400 | 0.400 |

### Honest conclusion

Both point estimates are higher in responders, but **neither gene has
statistically supported evidence of association with response**. The cohort is
too small to distinguish a reproducible biomarker signal from chance.

More importantly, these are whole-pleural-fluid pseudobulks, not
malignant-cell pseudobulks. Target detection varies sharply among patients:
TACSTD2-positive cells range from 2.3% to 48.7%, and CLDN4-positive cells from
2.3% to 48.3%. The largest responding sample drives much of the separation and
likely also differs in epithelial-cell abundance. Therefore, the observed
direction is hypothesis-generating and cannot be interpreted as a
tumor-cell-intrinsic ICI-response effect.

No pooled effect was calculated across excluded cohorts, no response labels
were reconstructed from paper figures, and no blood epithelial signal was
treated as tumor expression.

### Files and reproduction

- `series_audit.tsv`: reviewed series and explicit inclusion/exclusion reasons
- `GSE233203_pseudobulk.tsv`: patient-level expression and detection fractions
- `marker_outcome_tests.tsv`: exact tests and effect sizes
- `input_manifest.tsv`: source URLs, byte sizes, and SHA-256 checksums
- `GSE233203_*_response.svg`: individual patient values with medians

Run from repository root:

```bash
python3 scripts/w200/GEO_2025/analyze_gse233203.py
```

The script uses only the Python standard library. It downloads 545 MB into the
ignored `data/` directory and streams the sparse matrices without treating
cells as independent replicates.

## 中文

### 范围

本分析严格限定为 **2025 年在 GEO 公开**、且未进入既往参考集的人类肺癌 ICI
系列。只有同时满足以下条件的数据才做疗效检验：样本含肿瘤/上皮细胞、平台同时检测
TACSTD2 与 CLDN4、GEO 直接提供患者级 ICI 疗效标签且无需猜测式拼接。

`series_audit.tsv` 记录了 13 个重点候选及排除理由。最终只有 **GSE233203**
满足全部条件。排除本身也是结果：例如 GSE309652 有 72 例清晰 R/NR 标签，但
768 基因面板不含两个目标；GSE292098 检测两个目标，却没有公开 ROI 到临床结局的
对应键。

### 结果

GSE233203 含 7 例 IV 期肺腺癌胸水单细胞数据（缓解 3 例、未缓解 4 例）。
按患者汇总全部细胞原始计数，转为 log2(CPM + 1)，并在患者层面进行精确双侧置换
Mann-Whitney 检验；Holm P 值校正两个预设基因检验。

| 基因 | 缓解组中位数 | 未缓解组中位数 | Cliff's delta | 精确 P | Holm P |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | 5.743 | 3.702 | 0.833 | 0.114 | 0.229 |
| CLDN4 | 5.543 | 3.184 | 0.500 | 0.400 | 0.400 |

### 诚实结论

两个基因的点估计都在缓解组较高，但**均无统计学支持的疗效关联证据**。样本量只有
7 例，无法区分可重复的标志物信号与随机波动。

此外，这里是全胸水样本拟 bulk，而不是恶性细胞拟 bulk。不同患者 TACSTD2
阳性细胞比例为 2.3%–48.7%，CLDN4 为 2.3%–48.3%；一个高表达缓解样本对组间
差异影响很大，也很可能反映上皮细胞比例差异。因此当前方向仅能生成假设，不能解释
为肿瘤细胞内在的 ICI 疗效效应。

本分析未合并不合格队列、未从论文图片反推疗效标签，也未把血液中的上皮信号当作
肿瘤表达。

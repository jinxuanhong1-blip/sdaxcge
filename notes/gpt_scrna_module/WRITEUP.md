# TACSTD2/CLDN4 tumor-epithelial state versus immune co-occurrence in NSCLC scRNA-seq atlases

## English

### Question and scope

We tested whether a tumor-epithelial **TACSTD2/CLDN4** expression module covaries across specimens with **CD8**, **TLS-like**, or **LCAM** immune programs in the open processed NSCLC atlases GSE131907, GSE154826, and GSE148071. No ICI labels were required or used.

These are dissociated scRNA-seq data, so they contain no spatial coordinates. “Neighborhood” is therefore not spatially identifiable. The tested outcomes are explicitly **sample-level co-occurrence proxies**, not physical cell neighborhoods.

### Dataset eligibility and verification

We interpreted “<2 GB” as a strict decimal limit of 2,000,000,000 bytes.

| Atlas | Open processed input | Verified size | Decision |
|---|---:|---:|---|
| GSE131907 | Raw UMI matrix + author cell annotation | 408,736,818 + 1,886,187 bytes | Included |
| GSE148071 | GEO tar of 42 processed count matrices | 180,193,280 bytes | Included |
| GSE154826 | Complete author/HCA `lung_ldm.rds` | 2,123,085,302 bytes | Excluded: over limit |
| GSE154826 | 77 GEO processed sample archives, aggregate | 3,567,029,508 bytes | Excluded: no complete representation under limit |

SHA-256 values are in `results/gpt_scrna_module/source_manifest.csv`. The two included matrices contained 208,506 and 89,887 cells, respectively. GSE154826 metadata was inspected to verify its LCAM cluster definitions, but it did not enter any statistic.

### Predeclared analysis

1. Raw UMI counts were transformed cell-wise as `log1p(10,000 × count / total UMI)`.
2. GSE131907 used author annotations. Tumor epithelium was restricted to the author’s `Malignant cells` subtype. The cohort comprised all 21 malignant-cell-bearing specimens (`tL/B`, `mLN`, and `mBrain`), not the `tLung` TME-only samples.
3. GSE148071 has no GEO cell annotation. Cells were assigned by the maximum of predeclared epithelial, T, B, myeloid, endothelial, and fibroblast marker scores. Its epithelial cells are therefore **putative tumor epithelial cells**; malignancy was not established by CNV.
4. The epithelial exposure was each specimen’s median `TACSTD2`/`CLDN4` score among tumor/putative-tumor epithelial cells.
5. CD8 was scored in T cells (`CD8A`, `CD8B`, `CCL5`, `NKG7`, `GZMK`).
6. TLS-like co-occurrence averaged robust within-atlas z-scores for B (`MS4A1`, `CD79A`, `CD74`, `CD37`), Tfh-like (`CXCL13`, `CXCR5`, `ICOS`, `PDCD1`, `TOX`), and mregDC-like (`LAMP3`, `CCL19`, `CCL22`, `CD274`) programs in their respective lineages.
7. LCAM averaged robust z-scores for activated T (`PDCD1`, `CXCL13`, `TOX`, `HAVCR2`), IgG/plasma (`IGHG1`, `IGHG3`, `MZB1`, `JCHAIN`), and SPP1 macrophage (`SPP1`, `APOC1`, `APOE`, `LPL`) programs.
8. A specimen required at least 20 tumor epithelial, 20 T, 10 B, and 20 myeloid cells. This retained 18 GSE131907 and 25 GSE148071 specimens.
9. Primary tests were two-sided sample-level Spearman correlations with 100,000 label permutations, 20,000 sample bootstraps, and Benjamini–Hochberg correction across the three outcomes within each atlas. Sensitivity estimates residualized ranks for epithelial fraction, cell yield, median UMI depth, and—where applicable—tissue site.
10. A fixed-effect Fisher-z summary was reported as descriptive cross-atlas evidence; only two heterogeneous platforms were available.

All gene sets, thresholds, random seed (`20260816`), audit details, and per-sample values are machine-readable in `audit.json` and `sample_scores.csv`.

### Results

| Atlas | Outcome | n | Spearman ρ | 95% bootstrap CI | Permutation p | BH q | Adjusted rank ρ |
|---|---|---:|---:|---:|---:|---:|---:|
| GSE131907 | CD8 | 18 | 0.050 | −0.465, 0.546 | 0.844 | 0.892 | 0.373 |
| GSE131907 | TLS-like | 18 | −0.034 | −0.519, 0.480 | 0.892 | 0.892 | −0.048 |
| GSE131907 | LCAM | 18 | 0.093 | −0.429, 0.587 | 0.710 | 0.892 | 0.483 |
| GSE148071 | CD8 | 25 | −0.252 | −0.617, 0.187 | 0.221 | 0.332 | −0.314 |
| GSE148071 | TLS-like | 25 | −0.195 | −0.551, 0.229 | 0.350 | 0.350 | −0.154 |
| GSE148071 | LCAM | 25 | −0.462 | −0.707, −0.079 | 0.0211 | 0.0634 | −0.072 |

No association passed within-atlas FDR 0.05. GSE148071 showed a nominal inverse TACSTD2/CLDN4–LCAM association, but it missed FDR control and nearly vanished after adjustment, indicating sensitivity to specimen composition/depth rather than robust independent evidence.

Cross-atlas fixed-effect summaries were also null:

| Outcome | Total n | Summary ρ | 95% CI | p |
|---|---:|---:|---:|---:|
| CD8 | 43 | −0.132 | −0.426, 0.187 | 0.418 |
| TLS-like | 43 | −0.130 | −0.425, 0.189 | 0.425 |
| LCAM | 43 | −0.253 | −0.524, 0.063 | 0.115 |

### Conclusion

The eligible open processed atlases do **not** provide robust evidence that higher tumor-epithelial TACSTD2/CLDN4 is associated with CD8, TLS-like, or LCAM sample co-occurrence. The GSE148071 LCAM signal is exploratory only. Spatial data would be required to claim neighborhood exclusion or proximity, and independent malignant-cell annotation would strengthen GSE148071.

### Reproduction

```bash
python3 -m pip install -r scripts/gpt_scrna_module/requirements.txt
python3 scripts/gpt_scrna_module/download_data.py \
  --output results/gpt_scrna_module/source_data
python3 scripts/gpt_scrna_module/analyze_atlases.py \
  --source results/gpt_scrna_module/source_data \
  --output results/gpt_scrna_module
python3 scripts/gpt_scrna_module/verify_results.py \
  --results results/gpt_scrna_module
```

Primary sources:

- GSE131907: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907>
- GSE148071: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071>
- GSE154826: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154826>
- GSE154826 author metadata: <https://github.com/effiken/Leader_et_al>
- GSE154826 HCA project: <https://explore.data.humancellatlas.org/projects/b208466a-6fb0-4385-8cfb-8e03ff6b939e>

---

## 中文

### 问题与范围

本分析检验肿瘤上皮细胞的 **TACSTD2/CLDN4** 表达模块是否在样本层面与 **CD8**、**TLS 样**或 **LCAM** 免疫程序共同变化。分析对象为 GSE131907、GSE154826 和 GSE148071 的开放处理数据；不要求、也未使用 ICI 疗效标签。

这些数据来自解离后的 scRNA-seq，没有空间坐标。因此不能识别真正的空间“邻域”。本文检验的是明确标注的**样本层面共现代理分数**，不能解释为细胞间的物理邻近。

### 数据纳入与核验

“<2 GB”按严格十进制阈值 2,000,000,000 字节执行。

| 图谱 | 开放处理数据 | 核验大小 | 决策 |
|---|---:|---:|---|
| GSE131907 | 原始 UMI 矩阵及作者细胞注释 | 408,736,818 + 1,886,187 字节 | 纳入 |
| GSE148071 | GEO 中 42 个处理后计数矩阵的 tar 包 | 180,193,280 字节 | 纳入 |
| GSE154826 | 完整作者/HCA `lung_ldm.rds` | 2,123,085,302 字节 | 超过阈值，排除 |
| GSE154826 | 77 个 GEO 样本处理包合计 | 3,567,029,508 字节 | 无低于阈值的完整表示，排除 |

SHA-256 校验值见 `results/gpt_scrna_module/source_manifest.csv`。两个纳入矩阵分别包含 208,506 和 89,887 个细胞。GSE154826 的元数据用于核验 LCAM 簇定义，但未进入任何统计检验。

### 预先规定的分析方法

1. 每个细胞的原始 UMI 采用 `log1p(10,000 × count / total UMI)` 转换。
2. GSE131907 使用作者注释；肿瘤上皮严格限定为作者标注的 `Malignant cells`。纳入全部 21 个含恶性细胞的标本（`tL/B`、`mLN`、`mBrain`），不使用不含恶性细胞的 `tLung` TME 样本。
3. GSE148071 在 GEO 中没有细胞注释，因此按预先规定的上皮、T、B、髓系、内皮和成纤维标记模块最大值分类。其上皮细胞只能称为**推定肿瘤上皮细胞**，未通过 CNV 证明恶性。
4. 暴露变量为每个标本肿瘤/推定肿瘤上皮细胞中 `TACSTD2` 与 `CLDN4` 模块分数的中位数。
5. CD8 分数仅在 T 细胞中计算：`CD8A`、`CD8B`、`CCL5`、`NKG7`、`GZMK`。
6. TLS 样分数为三个谱系特异程序的图谱内稳健 z 分数均值：B 细胞程序、Tfh 样程序和 mregDC 样程序。
7. LCAM 分数为活化 T、IgG/浆细胞和 SPP1 巨噬细胞三个程序的稳健 z 分数均值。完整基因表见英文方法及 `audit.json`。
8. 每个标本至少需有 20 个肿瘤上皮、20 个 T、10 个 B 和 20 个髓系细胞；最终保留 GSE131907 的 18 个标本和 GSE148071 的 25 个标本。
9. 主检验为样本层面的双侧 Spearman 相关，使用 100,000 次标签置换、20,000 次样本 bootstrap，并在每个图谱内对三个结局做 Benjamini–Hochberg 校正。敏感性分析对上皮比例、细胞数、UMI 深度及适用时的组织部位进行秩残差化。
10. 两图谱 Fisher-z 固定效应汇总仅作描述，因为平台不同且图谱数只有两个。

所有基因集、阈值、随机种子（`20260816`）、审计信息和逐样本数值均保存在 `audit.json` 与 `sample_scores.csv`。

### 结果

| 图谱 | 结局 | n | Spearman ρ | 95% bootstrap CI | 置换 p | BH q | 校正后秩 ρ |
|---|---|---:|---:|---:|---:|---:|---:|
| GSE131907 | CD8 | 18 | 0.050 | −0.465, 0.546 | 0.844 | 0.892 | 0.373 |
| GSE131907 | TLS 样 | 18 | −0.034 | −0.519, 0.480 | 0.892 | 0.892 | −0.048 |
| GSE131907 | LCAM | 18 | 0.093 | −0.429, 0.587 | 0.710 | 0.892 | 0.483 |
| GSE148071 | CD8 | 25 | −0.252 | −0.617, 0.187 | 0.221 | 0.332 | −0.314 |
| GSE148071 | TLS 样 | 25 | −0.195 | −0.551, 0.229 | 0.350 | 0.350 | −0.154 |
| GSE148071 | LCAM | 25 | −0.462 | −0.707, −0.079 | 0.0211 | 0.0634 | −0.072 |

没有任何关联通过图谱内 FDR 0.05。GSE148071 中 TACSTD2/CLDN4 与 LCAM 呈名义上的负相关，但未通过 FDR，且校正后几乎消失，提示该信号可能受样本细胞组成或测序深度影响，不能视为稳健的独立证据。

跨图谱固定效应汇总同样不显著：CD8 的汇总 ρ = −0.132（95% CI −0.426 至 0.187，p = 0.418）；TLS 样 ρ = −0.130（−0.425 至 0.189，p = 0.425）；LCAM ρ = −0.253（−0.524 至 0.063，p = 0.115）。

### 结论

在满足大小限制的开放处理图谱中，**没有稳健证据**表明肿瘤上皮 TACSTD2/CLDN4 升高与 CD8、TLS 样或 LCAM 样本共现分数相关。GSE148071 的 LCAM 负相关只能作为探索性线索。若要提出免疫排斥或空间邻近结论，必须使用空间数据；GSE148071 还需要独立的恶性细胞注释验证。

复现命令与数据来源见英文部分。

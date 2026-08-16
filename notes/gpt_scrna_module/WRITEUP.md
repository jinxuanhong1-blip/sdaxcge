# TACSTD2/CLDN4 tumor-epithelial state versus immune co-occurrence in NSCLC scRNA-seq atlases

Machine-readable thesis table: `results/gpt_scrna_module/thesis_results.csv`.

## English

### Question and theses

Public processed lung scRNA-seq atlases GSE131907, GSE154826, and GSE148071 were tested. No ICI labels were required or used. Skip rule: EGA/dbGaP or a series missing **both** `TACSTD2` and `CLDN4`. All three series contain both genes.

These are dissociated scRNA-seq data, so they contain no spatial coordinates. “Neighborhood” is a **sample-level co-occurrence proxy**, not a physical neighborhood.

Predeclared theses:

1. **TROP2-high → immune-low after purity.** Higher tumor-epithelial `TACSTD2` associates with lower immune scores after residualizing immune ranks on epithelial fraction (purity proxy). Predicted: partial ρ < 0; OR of immune-low given TROP2-high > 1; log2FC of immune fraction (high/low TROP2) < 0.
2. **CLDN4 = tight-junction barrier.** Tumor-epithelial `CLDN4` tracks a CLDN4-independent TJ module (`CLDN3`, `CLDN7`, `OCLN`, `TJP1`, `F11R`) and the TROP2 epithelial state (predicted ρ > 0). A secondary barrier-immune test uses the same immune-low-after-purity predictions as thesis 1, with `CLDN4` as the exposure.

Match rule, applied to each row independently: **match** if the estimate is in the predicted direction and p < 0.05; **mismatch** if the opposite direction and p < 0.05; **inconclusive** if p ≥ 0.05 or the estimate is undefined.

### Dataset eligibility

| Atlas | Open processed input | Per-file size | Decision |
|---|---|---:|---|
| GSE131907 | Raw UMI matrix + author annotation | 408,736,818 + 1,886,187 | Included; both genes present |
| GSE148071 | GEO tar of 42 count matrices | 180,193,280 | Included; both genes present |
| GSE154826 | 29 clustering-used tumor MTX archives | each 11–153 MB | Included; both genes present |
| GSE154826 | Author/HCA `lung_ldm.rds` | 2,123,085,302 | Not used; per-sample MTX preferred |

SHA-256 values are in `source_manifest.csv`. Parsed cells: 208,506 (GSE131907), 89,887 (GSE148071), 164,000 matched GSE154826 barcodes. Eligible specimens after lineage minimums: 18, 25, and 19.

### Methods

1. Counts were transformed as `log1p(10,000 × count / total UMI)`.
2. GSE131907 used author labels. Tumor epithelium = `Malignant cells`. Cohort = 21 malignant-cell-bearing specimens (`tL/B`, `mLN`, `mBrain`).
3. GSE148071 has no GEO labels. Cells were assigned by predeclared marker-argmax. Epithelial cells are **putative**.
4. GSE154826 used author immune lineages. Gated `epi_endo_fibro_doublet` cells were classified by the same marker-argmax. Only clustering-used **Tumor** samples were scored.
5. Exposures are specimen medians of `TACSTD2`, `CLDN4`, their mean, and the TJ module among tumor/putative-tumor epithelial cells.
6. Immune outcomes: lineage-restricted CD8, TLS-like, and LCAM z-scores; their mean (`immune_score`); and `immune_fraction` = (T+B+myeloid)/n.
7. Specimen filters: ≥20 epithelial, ≥20 T, ≥10 B, ≥20 myeloid cells.
8. Combined-module tests remain two-sided Spearman with 100,000 permutations, 20,000 bootstraps, and within-atlas BH FDR.
9. Thesis tests use partial Spearman after purity, Fisher OR on median splits of purity-residual immune scores, and Mann–Whitney log2FC of immune fraction. Seed `20260816`.

### Thesis results

| Atlas | Hypothesis | n | Metric | Estimate | p | Match |
|---|---|---:|---|---:|---:|---|
| GSE131907 | TROP2-high → immune-low after purity | 18 | partial ρ | 0.280 | 0.266 | inconclusive |
| GSE131907 | TROP2-high → immune-low after purity | 18 | OR | 0.64 | 1.00 | inconclusive |
| GSE131907 | TROP2-high → immune-low after purity | 18 | log2FC | −0.116 | 0.724 | inconclusive |
| GSE148071 | TROP2-high → immune-low after purity | 25 | partial ρ | −0.286 | 0.172 | inconclusive |
| GSE148071 | TROP2-high → immune-low after purity | 25 | OR | 2.24 | 0.434 | inconclusive |
| GSE148071 | TROP2-high → immune-low after purity | 25 | log2FC | −0.459 | 0.109 | inconclusive |
| GSE154826 | TROP2-high → immune-low after purity | 19 | partial ρ | 0.152 | 0.720 | inconclusive |
| GSE154826 | TROP2-high → immune-low after purity | 19 | OR | 0.33 | 0.370 | inconclusive |
| GSE154826 | TROP2-high → immune-low after purity | 19 | log2FC | −0.103 | 0.055 | inconclusive |
| GSE131907 | CLDN4 barrier → immune-low after purity | 18 | partial ρ | −0.157 | 0.336 | inconclusive |
| GSE148071 | CLDN4 barrier → immune-low after purity | 25 | partial ρ | −0.493 | 0.0114 | **match** |
| GSE154826 | CLDN4 barrier → immune-low after purity | 19 | partial ρ | 0.188 | 0.954 | inconclusive |
| GSE131907 | CLDN4 tracks TJ module | 18 | ρ | 0.692 | 0.00147 | **match** |
| GSE148071 | CLDN4 tracks TJ module | 25 | ρ | 0.584 | 0.00218 | **match** |
| GSE154826 | CLDN4 tracks TJ module | 19 | ρ | 0.733 | 0.000358 | **match** |
| GSE131907 | CLDN4 co-occurs with TROP2 | 18 | ρ | 0.298 | 0.230 | inconclusive |
| GSE148071 | CLDN4 co-occurs with TROP2 | 25 | ρ | 0.374 | 0.0655 | inconclusive |
| GSE154826 | CLDN4 co-occurs with TROP2 | 19 | ρ | 0.673 | 0.00161 | **match** |

No TROP2-after-purity immune test reached p < 0.05. The only immune-low match is GSE148071 CLDN4 versus purity-adjusted `immune_score`. That same atlas’s CLDN4 OR and immune-fraction log2FC were inconclusive, so the immune-low claim is not internally consistent.

CLDN4 as a TJ-associated epithelial state is supported in all three atlases. TROP2–CLDN4 co-occurrence is significant only in GSE154826.

### Combined TACSTD2/CLDN4 module versus CD8/TLS/LCAM

| Atlas | Outcome | n | ρ | 95% CI | pperm | BH q |
|---|---|---:|---:|---:|---:|---:|
| GSE131907 | CD8 | 18 | 0.050 | −0.465, 0.546 | 0.844 | 0.892 |
| GSE131907 | TLS-like | 18 | −0.034 | −0.519, 0.480 | 0.892 | 0.892 |
| GSE131907 | LCAM | 18 | 0.093 | −0.429, 0.587 | 0.710 | 0.892 |
| GSE148071 | CD8 | 25 | −0.252 | −0.617, 0.187 | 0.221 | 0.332 |
| GSE148071 | TLS-like | 25 | −0.195 | −0.551, 0.229 | 0.350 | 0.350 |
| GSE148071 | LCAM | 25 | −0.462 | −0.707, −0.079 | 0.0211 | 0.0634 |
| GSE154826 | CD8 | 19 | 0.621 | 0.244, 0.844 | 0.00539 | 0.0162 |
| GSE154826 | TLS-like | 19 | 0.046 | −0.466, 0.597 | 0.852 | 0.852 |
| GSE154826 | LCAM | 19 | 0.381 | −0.115, 0.747 | 0.108 | 0.162 |

GSE154826 CD8 is a significant **positive** association (FDR 0.016). That is opposite an immune-low reading of the combined module. Three-atlas fixed-effect summaries are nonsignificant and heterogeneous for CD8 and LCAM.

### Conclusion

- **TROP2-high → immune-low after purity:** not supported. All nine TROP2 tests are inconclusive.
- **CLDN4 = TJ barrier:** supported as a transcriptional TJ-state marker in all three atlases. The secondary immune-low claim is a single-atlas match only (GSE148071 partial ρ) and does not replicate.
- Spatial exclusion cannot be claimed from dissociated data.

### Reproduction

```bash
python3 -m pip install -r scripts/gpt_scrna_module/requirements.txt
python3 scripts/gpt_scrna_module/download_data.py \
  --output results/gpt_scrna_module/source_data
python3 scripts/gpt_scrna_module/download_gse154826_batches.py \
  --source results/gpt_scrna_module/source_data \
  --output results/gpt_scrna_module/source_data/gse154826_batches
python3 scripts/gpt_scrna_module/analyze_atlases.py \
  --source results/gpt_scrna_module/source_data \
  --output results/gpt_scrna_module
python3 scripts/gpt_scrna_module/verify_results.py \
  --results results/gpt_scrna_module
```

---

## 中文

### 问题与假说

在 GSE131907、GSE154826、GSE148071 的开放处理肺 scRNA 图谱中检验肿瘤上皮 **TACSTD2/CLDN4**。不使用 ICI 标签。仅跳过 EGA/dbGaP，或同时缺少这两个基因的系列。三个系列均同时含有两基因。

数据为解离 scRNA-seq，没有空间坐标。“邻域”只是**样本层面共现代理**。

预先规定的假说：

1. **TROP2 高 → 校正纯度后免疫低。** 预测：偏 ρ < 0；TROP2 高组出现免疫低的 OR > 1；免疫细胞比例的 log2FC < 0。
2. **CLDN4 = 紧密连接屏障。** 预测 CLDN4 与不含 CLDN4 的 TJ 模块及 TROP2 正相关；次级检验为 CLDN4 高 → 校正纯度后免疫低。

判定规则：估计值方向符合且 p < 0.05 为 **match**；方向相反且 p < 0.05 为 **mismatch**；否则 **inconclusive**。

### 数据纳入

GSE131907 与 GSE148071 使用原先 GEO 处理矩阵。GSE154826 不使用超过 2,000,000,000 字节的完整 RDS，改为 29 个聚类所用肿瘤样本的 GEO MTX（单文件 11–153 MB）。过滤后标本数为 18、25、19。

### 论题结果

| 图谱 | 假说 | n | 指标 | 估计值 | p | 判定 |
|---|---|---:|---|---:|---:|---|
| GSE131907 | TROP2 高 → 校正纯度后免疫低 | 18 | 偏 ρ | 0.280 | 0.266 | inconclusive |
| GSE148071 | TROP2 高 → 校正纯度后免疫低 | 25 | 偏 ρ | −0.286 | 0.172 | inconclusive |
| GSE154826 | TROP2 高 → 校正纯度后免疫低 | 19 | 偏 ρ | 0.152 | 0.720 | inconclusive |
| GSE148071 | CLDN4 屏障 → 校正纯度后免疫低 | 25 | 偏 ρ | −0.493 | 0.0114 | **match** |
| GSE131907 | CLDN4 跟随 TJ 模块 | 18 | ρ | 0.692 | 0.00147 | **match** |
| GSE148071 | CLDN4 跟随 TJ 模块 | 25 | ρ | 0.584 | 0.00218 | **match** |
| GSE154826 | CLDN4 跟随 TJ 模块 | 19 | ρ | 0.733 | 0.000358 | **match** |
| GSE154826 | CLDN4 与 TROP2 共现 | 19 | ρ | 0.673 | 0.00161 | **match** |

完整 OR/log2FC 行见 `thesis_results.csv`。TROP2 的 9 项免疫检验全部 inconclusive。CLDN4 作为 TJ 转录状态在三个图谱中均成立。免疫低仅在 GSE148071 的偏相关中成立，OR 与 log2FC 并不支持，且另外两个图谱不重复。

联合模块方面，GSE154826 的 CD8 为显著正相关（ρ = 0.621，置换 p = 0.00539，q = 0.0162），与“免疫低”读法相反。三图谱固定效应汇总不显著。

### 结论

- **TROP2 高 → 校正纯度后免疫低：** 不被支持。
- **CLDN4 = TJ 屏障：** 作为 TJ 相关上皮状态被三个图谱支持；免疫低只是单图谱、单指标的 match。
- 解离数据不能证明空间排斥。

复现命令见英文部分。

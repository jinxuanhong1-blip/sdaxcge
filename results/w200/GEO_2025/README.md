# GEO 2025 lung ICI leftovers: TACSTD2 / CLDN4

## English

### Scope

This is a leftover slice: **GEO-public-in-2025** human lung ICI series that
were not in the previously analyzed reference set (GSE126044, GSE135222,
GSE136961, GSE166449, GSE93157, GSE207422, GSE205335). A leftover series was
eligible for an outcome test only if:

1. the deposited material contained tumor/epithelial cells;
2. both TACSTD2 and CLDN4 were measured; and
3. GEO supplied a patient-level ICI outcome without an inferred join.

The first pass searched only `anti-PD-1` phrasing and therefore missed
`GSE291670` (`PD-1 blockade` / camrelizumab). This continuation catalogs 61
leftover accessions in `leftover_catalog.tsv` and adds that series.

### Outcome tests

| Dataset | Setting | Gene | n | Median pos | Median neg | Cliff's δ | Exact p | Holm p |
|---|---|---|---:|---:|---:|---:|---:|---:|
| GSE233203 | pleural fluid, Response vs Non-response | TACSTD2 | 3 vs 4 | 5.743 | 3.702 | 0.833 | 0.114 | 0.229 |
| GSE233203 | pleural fluid, Response vs Non-response | CLDN4 | 3 vs 4 | 5.543 | 3.184 | 0.500 | 0.400 | 0.400 |
| GSE291670 | post-treatment lung, MPR vs Non-MPR | TACSTD2 | 3 vs 3 | 6.356 | 3.902 | 0.778 | 0.200 | 0.400 |
| GSE291670 | post-treatment lung, MPR vs Non-MPR | CLDN4 | 3 vs 3 | 5.563 | 3.864 | 0.556 | 0.400 | 0.400 |

### Honest conclusion

**No leftover 2025 series supports TACSTD2 or CLDN4 as an ICI-response
biomarker.** Both eligible leftover cohorts are n≤7. All four exact tests are
nonsignificant after a two-gene Holm correction.

GSE291670 is leftover-important and leftover-limited at the same time. It is
post-treatment surgical tissue after four cycles of camrelizumab plus
anlotinib. Epithelial genes in MPR versus Non-MPR are therefore sensitive to
residual tumor burden. The observed direction is higher, not lower, in MPR,
so a simple “MPR has fewer tumor cells” explanation does not fit these six
pseudobulks. That still does not make a biomarker: n=3 versus 3, exact
p=0.20/0.40, and detection fractions still vary from 2.5% to 10.6%.

GSE233203 remains composition-sensitive pleural-fluid leftover. One high
epithelial-fraction responder still drives most of the separation.

### Leftover series that were not outcome-tested

GSE309652 remains the cleanest leftover R/NR lung cohort (n=72), but GPL31904
measures neither gene. GSE206127 is leftover nivolumab response at n=213, but
it is serum miRNA. Blood leftover series (GSE202417, GSE295969, GSE285888,
GSE190905) were not treated as tumor expression.

GSE292098 leftover GeoMx WTA does measure both genes. GEO has no per-ROI
outcome key, so no response test was performed. Compartment medians are
honest leftover biology only:

| Compartment | TACSTD2 median raw | CLDN4 median raw | n segments |
|---|---:|---:|---:|
| CK (tumor) | 245 | 199 | 81 |
| CD45 (immune) | 10 | 10 | 72 |
| CD68 (macrophage) | 11 | 10 | 81 |

Both genes are about 20-fold higher in CK than in immune compartments. That
supports leftover epithelial restriction. It is not an ICI-outcome result.

### Files and reproduction

```bash
python3 scripts/w200/GEO_2025/analyze_scrna_leftovers.py
python3 scripts/w200/GEO_2025/describe_gse292098.py
```

Standard library only. Large raw archives stay in ignored `data/`.

## 中文

### 范围

本切片是 **2025 年 GEO 公开**、且未进入既往参考集的人类肺癌 ICI **剩余系列**。
只有同时满足肿瘤/上皮样本、平台同时检测 TACSTD2 与 CLDN4、以及 GEO
直接提供患者级疗效标签的剩余系列才做疗效检验。

第一轮检索只用了 `anti-PD-1` 措辞，因此漏掉了 `GSE291670`（`PD-1 blockade` /
卡瑞利珠单抗）。本轮补全 61 个剩余编号，并补做该系列。

### 诚实结论

**没有一个 2025 剩余系列能支持 TACSTD2 或 CLDN4 作为 ICI 疗效标志物。**
两个合格剩余队列都只有 ≤7 例；四个精确检验在两基因 Holm 校正后均不显著。

GSE291670 是治疗后手术标本。按常理，MPR 残存肿瘤细胞更少，上皮基因应更低；
这 6 个拟 bulk 的点估计却是 MPR 更高。这不能用“残瘤更少”一句话解释，但也
不能当成标志物：n=3 vs 3，精确 P=0.20/0.40。

GSE292098 剩余 GeoMx 能测到两个基因，但 GEO 没有 ROI 级疗效键，因此只报告
隔室中位数：CK 约比 CD45/CD68 高 20 倍。这只说明上皮限制，不是疗效结果。

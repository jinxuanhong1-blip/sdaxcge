# Public mouse KL vs KP, beyond GSE137244

KL means Kras plus Lkb1/Stk11 loss with Trp53 intact. KP means Kras plus Trp53 loss with Lkb1 intact. Nothing here is merged with the private 8KL single-cell matrices.

The locked GSE137244 cell-line result is reproduced, not replaced: log2(FPKM+1), KL n=5 (KL155mix, KL47-1, KLC, KLD, KLE) vs KP n=5 (B6AL10-1..5), two-sided exact Mann–Whitney. Tacstd2 Δ=+3.238, Cldn4 Δ=+5.570, p=0.00794 for both (complete separation; that p is the floor for n=5 vs 5). A seven-gene tight-junction mean (Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln) is Δ=+3.269 on the same matrix. The locked TJ figure of +3.03 used a different gene average; this seven-gene mean is reported under its own name and is not substituted for +3.03.

Usable n means both arms have at least 3 independent samples and a per-sample matrix that contains the gene. n=3 vs 3 cannot reach a two-sided exact Mann–Whitney p below 0.10.

## What was open

| Accession | Material | KL vs KP n | Scored |
|---|---|---|---|
| GSE137244 | GEMM cell lines | 5 vs 5 | Locked reference, reproduced |
| GSE137396 | GEMM lung nodules, same paper | 5 vs 5 | Rescored |
| GSE164758 | Primary bulk tumors (Eichner) | 9 vs 8 | New |
| GSE6135 | Ji 2007 microarray, mouse-level | 7 vs 5 | Tacstd2 new; Cldn4 matches the prior mouse-level result |
| GSE244452 | Syngeneic bulk tumors | 3 vs 3 | New; p floor 0.10 |
| GSE274352 | Empty-vector cell lines | 3 vs 3 | New; p floor 0.10 |
| GSE274351 | LCM early adenomas | 5 vs 5 | Rescored |
| GSE165640 | Bulk lungs, KL/KP plus Il1f knockouts | pure KP n=2 | Not scored; no KL-vs-KP matrix |
| GSE277929, GSE193895, GSE21581, GSE69552 | KL or K/KL only | no KP arm | Not a KL vs KP contrast |
| GSE180963, GSE165641, GSE154977 | scRNA | one genotype | Not merged with private 8KL |
| KP-Tracer (Yang, Cell 2022) | scRNA KP vs KPL | KPL is not KL | Not scored |
| Organoids | — | — | No public KL vs KP organoid RNA-seq or scRNA |

GEO queries were mouse lung expression with Stk11/Lkb1 and Trp53/p53, plus organoid and single-cell filters. No deposited mouse lung organoid series contains both KL and KP. Dost GSE150425 is Kras only. GSE227719 is KP organoids only.

## KL vs KP

Deltas are mean(KL) − mean(KP). p is two-sided exact Mann–Whitney.

| Accession | Scale | Tacstd2 | Cldn4 | TJ7 |
|---|---|---|---|---|
| GSE137244 cell lines | log2(FPKM+1) | **+3.238, p=0.00794** | **+5.570, p=0.00794** | **+3.269, p=0.00794** |
| GSE164758 primary tumors | log2(FPKM+1) | **+0.858, p=8.2×10⁻⁵** | −0.355, p=0.42 | +0.031, p=0.81 |
| GSE6135 all histologies, mouse | array value | **+1.795, p=0.00253** | +0.224, p=0.88 | +0.043, p=1 |
| GSE6135 adenocarcinoma only, mouse | array value | **+0.989, p=0.0159** | −0.590, p=0.19 | −0.175, p=0.41 |
| GSE137396 nodules | log2(abundance+1) | +1.144, p=0.095 | +0.557, p=0.84 | +0.311, p=0.55 |
| GSE244452 syngeneic tumors | log2(count+1) | +7.228, p=0.10 | +8.976, p=0.10 | +8.271, p=0.10 |
| GSE274352 empty-vector lines | log2(count+1) | −1.543, p=0.10 | −0.196, p=1 | −0.390, p=0.70 |
| GSE274351 LCM adenomas | log2(TPM+1) | −0.860, p=0.55 | −1.036, p=0.69 | −0.693, p=0.095 |

The two results with n large enough for p<0.05 both separate on Tacstd2. Raw Cldn4 overlaps in both. A later sweep of histology, probes, epithelial normalization, and tight-junction modules is in the next section. It does not add a second fair raw-Cldn4 result.

**GSE164758** is untreated primary tumors (KL n=9, KP n=8). Every KL tumor has higher Tacstd2 than every KP tumor (KL 6.21–6.99, KP 4.91–6.15 on log2(FPKM+1)), so p=8.2×10⁻⁵ is the minimum for 9 vs 8. Cldn4 overlaps (Δ=−0.355, p=0.42). These are bulk tumors, so stroma is in the RNA.

**GSE6135** (Ji et al. 2007) uses the same mouse cohort as the earlier Cldn4 contrast: Lkb1 L/L or L/− primaries, L/+ heterozygotes and the metastasis removed, multiple tumors from one mouse averaged. Cldn4 Δ=+0.224, p=0.88, which matches that earlier mouse-level number. Tacstd2 was not in that table. All 7 KL mice sit above all 5 KP mice (Δ=+1.795, p=0.00253, the minimum for 7 vs 5). Restricting to adenocarcinoma primaries (4 KL mice vs 5 KP) still separates Tacstd2 (Δ=+0.989, p=0.0159, again the minimum) and moves Cldn4 the other way (Δ=−0.590, p=0.19). Squamous and mixed KL tumors are the highest Tacstd2 samples, and they are not required for the separation.

**GSE137396** nodules from the GSE137244 paper are a trend for Tacstd2 (Δ=+1.144, p=0.095) and null for Cldn4 (Δ=+0.557, p=0.84). KL samples are labeled NA and KP samples are labeled vehicle.

**GSE244452** syngeneic tumors are completely separated for both genes (Tacstd2 Δ=+7.228, Cldn4 Δ=+8.976) but n=3 vs 3, so p=0.10 is the floor. Krt8 is high in both arms (KP mean 1.6×10⁴, KL mean 3.6×10⁴ deposited counts), so the KP samples are not empty of epithelium. Sftpc is detected only in KP. Cd8a is lower in KL (mean 13 vs 74). Cdh1 and F11r are absent from the 7,449-gene deposit, so the TJ7 row for this series is the mean of the five tight-junction genes that are present.

**GSE274352** empty-vector lines go the other way for Tacstd2: all three KP lines are above all three KL lines (Δ=−1.543, p=0.10). Cldn4 is flat (Δ=−0.196, p=1).

**GSE274351** early LCM adenomas do not show KL above KP (Tacstd2 Δ=−0.860, p=0.55; Cldn4 Δ=−1.036, p=0.69). Supplementary columns are K1–K4, KL1–KL5, KP1–KP5. GEO sample text instead lists five K adenomas and four KL adenomas. Stk11 TPM is near the floor in all five KL-labeled columns, which is why those labels are the ones used.

## Cldn4 limb sweep

`sweep.py` keeps the same two-sided exact Mann–Whitney and the same delta. Welch p is stored beside it. Benjamini–Hochberg FDR is computed inside each accession, cohort, and family. Pre-specified Tacstd2 and Cldn4 from the table above stay the primary calls; the gene-family FDR is the penalty for looking across the tight-junction list.

The sweep covers every scored KL vs KP matrix: histology filters, the single GPL8321 probe for each GSE6135 gene, Cldn4 minus Epcam, the residual of Cldn4 on Epcam and Krt8, and five tight-junction modules (mean of within-contrast gene z-scores). Module deltas are z-score means. They are not the raw seven-gene mean of +3.269. No public single-cell matrix contains both KL and KP with usable mouse n, so there is no epithelial-only single-cell contrast. Private 8KL matrices were not opened.

**Raw Cldn4.** The only fair specification with Mann–Whitney p<0.05 is the locked GSE137244 cell lines (Δ=+5.570, p=0.00794, FDR=0.0098). Everywhere else the fair raw call stays null:

| Specification | Cldn4 Δ | Mann–Whitney | Welch |
|---|---|---|---|
| GSE164758 all tumors | −0.355 | 0.42 | 0.28 |
| GSE164758, Cldn4 minus Epcam | −0.525 | 0.24 | 0.16 |
| GSE164758, residual on Epcam+Krt8 | −0.500 | 0.074 | 0.092 |
| GSE6135 all histologies, gene and the only probe `1418283_at` | +0.224 | 0.88 | 0.63 |
| GSE6135 adenocarcinoma only | −0.590 | 0.19 | 0.18 |
| GSE6135 adenocarcinoma plus mixed | −0.244 | 0.55 | 0.61 |
| GSE137396 | +0.557 | 0.84 | 0.48 |
| GSE274351 | −1.036 | 0.69 | 0.49 |
| GSE244452 n=3 vs 3 | +8.976 | 0.10 (floor) | 0.0070 |
| GSE244452, Cldn4 minus Epcam | −0.960 | 0.40 | 0.34 |
| GSE274352 n=3 vs 3 | −0.196 | 1 | 0.77 |

GSE244452 is completely separated on raw Cldn4, and n=3 vs 3 cannot go below p=0.10. Epcam is about a thousand-fold higher in the KL counts, so Cldn4 minus Epcam reverses the sign. The residual on Epcam and Krt8 is about zero (Δ=−0.019, p=0.70). GPL8321 carries one Cldn4 probe, so a probe swap cannot move GSE6135.

GSE164758 has no pathologist histology field. A marker proxy (squamous: Krt5, Krt14, Trp63; adeno: Nkx2-1, Sftpc, Sftpb) calls all 9 KL tumors adeno-like and 7 of 8 KP tumors squamous-like. The only adeno-like KP tumor is KP116, so a matched adeno-like contrast has KP n=1 and was not tested. There is no squamous-like KL subset. Lower Cldn4 in these KL tumors sits next to higher adeno-marker scores, not next to a squamous program.

One unfair GSE6135 contrast does reach p<0.05: squamous or mixed KL mice (n=3) versus adenocarcinoma KP mice (n=5), Cldn4 Δ=+1.310, p=0.036 (the floor for 3 vs 5), FDR=0.067 inside that cohort. The same rows are positive after the Epcam adjustments. The file marks `fair_genotype=False`. That contrast compares histologies.

**Tight-junction genes that do separate, fair cohorts.**

GSE164758 all tumors. These seven genes are completely separated (p=8.2×10⁻⁵, gene-family FDR=1.9×10⁻⁴):

| Direction | Gene | Δ |
|---|---|---|
| KL higher | Cldn3 | +3.843 |
| KL higher | Tjp3 | +1.687 |
| KL higher | Cldn7 | +1.211 |
| KL higher | Tacstd2 | +0.858 |
| KL lower | Cldn6 | −2.641 |
| KL lower | Ocln | −1.086 |
| KL lower | F11r | −0.533 |

Three more are KL lower with p<0.05 and are not complete separation: Marveld2 Δ=−0.321, p=0.00099, FDR=0.0020; Cldn1 Δ=−0.592, p=0.0037, FDR=0.0059; Tjp1 Δ=−0.400, p=0.0037, FDR=0.0059. Cldn4 overlaps (Δ=−0.355, p=0.42). Tjp2 is flat (Δ=−0.009, p=1). Cdh1, Cgn, Marveld3, and Cldn18 also have Mann–Whitney p>0.05.

The z-score modules follow the split: TJ7 Δ=−0.395, p=0.021, FDR=0.10 (the lower genes pull the seven-gene mean down). CLDN_STRAND (Cldn3, Cldn4, Cldn7, Ocln) Δ=+0.365, p=0.046, FDR=0.12. TJ_CORE Δ=−0.203, p=0.54. A whole-core mean is not a KL>KP result on this series. Cldn3, Cldn7, and Tjp3 are.

GSE6135 mouse-level, all histologies: Tjp1 Δ=+0.888, p=0.010, FDR=0.038; Marveld2 Δ=+0.844, p=0.010, FDR=0.038; Cdh1 Δ=+0.556, p=0.018, FDR=0.053. The plaque module (Tjp1, Tjp2, Tjp3, Cgn, Marveld2) is Δ=+0.971, p=0.030, FDR=0.15. Cldn3 and Cldn7 do not separate (all-histology Cldn7 Δ=+0.491, p=0.34; adenocarcinoma-only Cldn7 Δ=+0.008, p=0.90). Adenocarcinoma plus mixed still has Tjp1 Δ=+0.786 and Marveld2 Δ=+0.798, both p=0.032, FDR=0.12. Adenocarcinoma-only Tjp1 is Δ=+0.701, Mann–Whitney p=0.063, Welch p=0.035. Cldn18 has two probes. Gene-level adenocarcinoma-only Δ=+1.188, Mann–Whitney p=0.11, Welch p=0.045. One probe, `1449428_at`, is p=0.032 with FDR=0.28. That is one probe, not a gene call. Ocln is lower in KL (all-histology Welch p=0.049, Mann–Whitney p=0.073).

GSE137244 cell lines separate on the z-score modules as well (TJ7, TJ_CORE, TJ_CORE without Cldn4, CLDN_STRAND, PLAQUE; each Mann–Whitney p=0.00794). Those z-score deltas are not the raw +3.269. GSE137396 has nominal KL>KP Tjp1 and Cdh1 (both p=0.016, FDR=0.063) and KL-lower Cgn and Ocln (Ocln Δ=−0.462, p=0.016). Its modules are null. GSE244452 modules are completely separated at the n=3 floor (TJ_CORE Welch p=0.0012). GSE274351 modules trend KL lower (TJ7 Δ=−0.51, p=0.056).

Tacstd2 positives from the first table are unchanged by the sweep: GSE137244 +3.238 p=0.00794; GSE164758 +0.858 p=8.2×10⁻⁵; GSE6135 all histologies +1.795 p=0.00253; adenocarcinoma plus mixed +1.337 p=0.00794; adenocarcinoma only +0.989 p=0.0159; GSE244452 +7.228 with Welch p=1.9×10⁻⁴ at the Mann–Whitney floor of 0.10. GSE137396 remains a trend (+1.144, p=0.095). GSE274352 remains the opposite direction (−1.543, p=0.10).

## NHEJ, STING, and IFN

Pathway scores use the same z-score mean and the same two-sided exact Mann–Whitney, only on series where the rank-test floor can fall below 0.05. GSE244452 and GSE274352 (n=3) are omitted. Gene lists are the genes present in each matrix.

| Accession | NHEJ | STING | IFN |
|---|---|---|---|
| GSE164758 tumors, 9 vs 8 | **−1.074, p=8.2×10⁻⁵**, FDR=2.5×10⁻⁴ | +0.348, p=0.20 | +0.258, p=0.54 |
| GSE137244 cell lines, 5 vs 5 | **+0.493, p=0.032**, FDR=0.095 | +0.225, p=0.84 | −0.725, p=0.22 |
| GSE137396 nodules, 5 vs 5 | +0.108, p=0.84 | +0.450, p=0.55 | −0.003, p=0.84 |
| GSE274351 LCM, 5 vs 5 | −0.292, p=0.42 | −0.105, p=0.84 | −0.633, p=0.056 |
| GSE6135 all histologies, 7 vs 5 | −0.389, p=0.27 | +0.246, p=0.27 | +0.570, p=0.073 |
| GSE6135 adenocarcinoma, 4 vs 5 | −0.681, p=0.41 | +0.684, p=0.063 | +0.127, p=0.73 |

On GSE164758 the NHEJ module uses Xrcc4, Xrcc5, Xrcc6, Lig4, Prkdc, Nhej1, and Dclre1c. Xrcc4 (Δ=−0.447), Xrcc5 (−0.663), and Xrcc6 (−0.288) are each completely separated with KL lower (p=8.2×10⁻⁵, FDR=4.3×10⁻⁴). Prkdc Δ=−0.247, p=0.0025, FDR=0.0092. Dclre1c Δ=−0.207, p=0.015, FDR=0.036. Nhej1 is the exception inside the module: Δ=+0.221, p=0.011, FDR=0.029, Welch p=0.055. The STING module (Tmem173, Mb21d1, Tbk1, Irf3) is null, and two of its genes go up in KL: Irf3 Δ=+0.770, p=8.2×10⁻⁵; Tmem173 Δ=+0.621, p=0.0037, FDR=0.012. Tbk1 goes down (Δ=−0.160, p=0.00033, FDR=0.0014). The 15-gene IFN module is null. Inside it, Mx1 is higher in KL (Δ=+0.756, p=0.0055, FDR=0.016) and Ifih1 is lower (Δ=−0.473, complete separation). Stat2 (Δ=+0.611) and Isg15 (Δ=+0.924) have Mann–Whitney p=0.027 and FDR=0.055.

The cell-line NHEJ direction is the other way. GSE137244 Lig4 Δ=+0.694, p=0.016, FDR=0.061; Prkdc Δ=+0.637, p=0.032; Dclre1c Δ=+0.374, p=0.032. Paxx, present only in that matrix, is lower in KL (Δ=−0.953, p=0.032). Tbk1 is higher (Δ=+0.525, p=0.00794, FDR=0.061). STING and IFN modules on the cell lines are null.

GSE6135 all-histology IFN is a trend (p=0.073, Welch p=0.051). Genes that do separate there, KL higher: Tbk1 Δ=+0.890 and Irf1 Δ=+0.869, both p=0.00253, FDR=0.028; Stat1 Δ=+0.706, p=0.0051, FDR=0.037. Xrcc6 is lower (Δ=−0.580, p=0.030, FDR=0.13). Lig4, Cgas, and Mb21d1 are absent from the GPL8321 symbol match used here, so those modules are the genes that were present. Adenocarcinoma-only Tbk1 and Irf1 keep Mann–Whitney p=0.016 with FDR=0.17 inside the larger gene family.

## KPL is not KL

On the same GSE164758 tumors, KPL (Kras/p53/Lkb1, n=15) versus KP (n=8) is a different genotype. Tacstd2 Δ=+0.922, p=0.0019. Cldn4 Δ=+1.595, p=8.2×10⁻⁶. TJ7 Δ=+0.874, p=3.9×10⁻⁴. Bulk Cldn4 is higher when Lkb1 is lost on a p53-null background. It is not higher in the p53-intact KL tumors from that series.

## 中文

锁死的 GSE137244 细胞系（KL n=5 vs KP n=5）按 log2(FPKM+1) 复现：Tacstd2 Δ=+3.238，Cldn4 Δ=+5.570，精确 Mann–Whitney p=0.00794。七基因 TJ 均值是 +3.269，不等于原文的 TJ +3.03。

能做 KL vs KP、且 n 够用的新矩阵里，Tacstd2 在两套 n 足够的数据上分开，Cldn4 没有。GSE164758 原发瘤 KL n=9 vs KP n=8：Tacstd2 完全分开（Δ=+0.858，p=8.2×10⁻⁵），Cldn4 Δ=−0.355，p=0.42。GSE6135 小鼠水平（与先前 Cldn4 同一队列）Tacstd2 完全分开（7 vs 5，Δ=+1.795，p=0.00253）；只留腺癌仍分开（4 vs 5，Δ=+0.989，p=0.0159）。Cldn4 在该队列仍是 +0.224，p=0.88。GSE244452 方向相同但 n=3，p 下限 0.10。GSE274352 的 Tacstd2 方向相反（KP 更高，n=3）。早期 LCM（GSE274351）和结节（GSE137396）都不到 p<0.05。没有公开的 KL vs KP 类器官。单细胞里没有同时含 KL 和 KP、且能用的小鼠 n；不与私有 8KL 合并。GSE164758 里的 KPL vs KP 是三重突变，不是 KL。

Cldn4 扫描（组织学、探针、上皮归一化、TJ 模块）没有把公平对比里的原始 Cldn4 做成第二个显著结果。GSE164758 上完全分开、且 KL 更高的是 Cldn3（Δ=+3.843）、Tjp3（Δ=+1.687）、Cldn7（Δ=+1.211），同一 p=8.2×10⁻⁵；Cldn6、Ocln、F11r、Tjp1 则是 KL 更低。整段 TJ_CORE 均值不显著（Δ=−0.203，p=0.54），TJ7 均值是 KL 更低（Δ=−0.395，p=0.021）。GSE6135 公平对比里显著的是 Tjp1 与 Marveld2（全组织学 p=0.010），不是 Cldn4。鳞癌/混合 KL 对腺癌 KP 的 Cldn4 p=0.036 是组织学混杂，标为 unfair。GSE244452 的 Cldn4 完全分开但 n=3，减去 Epcam 后符号反过来。NHEJ 模块在 GSE164758 原发瘤是 KL 更低（Δ=−1.074，p=8.2×10⁻⁵），在 GSE137244 细胞系是 KL 更高（Δ=+0.493，p=0.032）。STING 与 IFN 模块在这些 n 足够的队列里都不到 p<0.05。GSE164758 上 Tmem173 与 Irf3 是 KL 更高。

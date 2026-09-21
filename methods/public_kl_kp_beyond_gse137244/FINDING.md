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

The two results with n large enough for p<0.05 both separate on Tacstd2 and do not separate on Cldn4.

**GSE164758** is untreated primary tumors (KL n=9, KP n=8). Every KL tumor has higher Tacstd2 than every KP tumor (KL 6.21–6.99, KP 4.91–6.15 on log2(FPKM+1)), so p=8.2×10⁻⁵ is the minimum for 9 vs 8. Cldn4 overlaps (Δ=−0.355, p=0.42). These are bulk tumors, so stroma is in the RNA.

**GSE6135** (Ji et al. 2007) uses the same mouse cohort as the earlier Cldn4 contrast: Lkb1 L/L or L/− primaries, L/+ heterozygotes and the metastasis removed, multiple tumors from one mouse averaged. Cldn4 Δ=+0.224, p=0.88, which matches that earlier mouse-level number. Tacstd2 was not in that table. All 7 KL mice sit above all 5 KP mice (Δ=+1.795, p=0.00253, the minimum for 7 vs 5). Restricting to adenocarcinoma primaries (4 KL mice vs 5 KP) still separates Tacstd2 (Δ=+0.989, p=0.0159, again the minimum) and moves Cldn4 the other way (Δ=−0.590, p=0.19). Squamous and mixed KL tumors are the highest Tacstd2 samples, and they are not required for the separation.

**GSE137396** nodules from the GSE137244 paper are a trend for Tacstd2 (Δ=+1.144, p=0.095) and null for Cldn4 (Δ=+0.557, p=0.84). KL samples are labeled NA and KP samples are labeled vehicle.

**GSE244452** syngeneic tumors are completely separated for both genes (Tacstd2 Δ=+7.228, Cldn4 Δ=+8.976) but n=3 vs 3, so p=0.10 is the floor. Krt8 is high in both arms (KP mean 1.6×10⁴, KL mean 3.6×10⁴ deposited counts), so the KP samples are not empty of epithelium. Sftpc is detected only in KP. Cd8a is lower in KL (mean 13 vs 74). Cdh1 and F11r are absent from the 7,449-gene deposit, so the TJ7 row for this series is the mean of the five tight-junction genes that are present.

**GSE274352** empty-vector lines go the other way for Tacstd2: all three KP lines are above all three KL lines (Δ=−1.543, p=0.10). Cldn4 is flat (Δ=−0.196, p=1).

**GSE274351** early LCM adenomas do not show KL above KP (Tacstd2 Δ=−0.860, p=0.55; Cldn4 Δ=−1.036, p=0.69). Supplementary columns are K1–K4, KL1–KL5, KP1–KP5. GEO sample text instead lists five K adenomas and four KL adenomas. Stk11 TPM is near the floor in all five KL-labeled columns, which is why those labels are the ones used.

## KPL is not KL

On the same GSE164758 tumors, KPL (Kras/p53/Lkb1, n=15) versus KP (n=8) is a different genotype. Tacstd2 Δ=+0.922, p=0.0019. Cldn4 Δ=+1.595, p=8.2×10⁻⁶. TJ7 Δ=+0.874, p=3.9×10⁻⁴. Bulk Cldn4 is higher when Lkb1 is lost on a p53-null background. It is not higher in the p53-intact KL tumors from that series.

## 中文

锁死的 GSE137244 细胞系（KL n=5 vs KP n=5）按 log2(FPKM+1) 复现：Tacstd2 Δ=+3.238，Cldn4 Δ=+5.570，精确 Mann–Whitney p=0.00794。七基因 TJ 均值是 +3.269，不等于原文的 TJ +3.03。

能做 KL vs KP、且 n 够用的新矩阵里，Tacstd2 在两套 n 足够的数据上分开，Cldn4 没有。GSE164758 原发瘤 KL n=9 vs KP n=8：Tacstd2 完全分开（Δ=+0.858，p=8.2×10⁻⁵），Cldn4 Δ=−0.355，p=0.42。GSE6135 小鼠水平（与先前 Cldn4 同一队列）Tacstd2 完全分开（7 vs 5，Δ=+1.795，p=0.00253）；只留腺癌仍分开（4 vs 5，Δ=+0.989，p=0.0159）。Cldn4 在该队列仍是 +0.224，p=0.88。GSE244452 方向相同但 n=3，p 下限 0.10。GSE274352 的 Tacstd2 方向相反（KP 更高，n=3）。早期 LCM（GSE274351）和结节（GSE137396）都不到 p<0.05。没有公开的 KL vs KP 类器官。单细胞里没有同时含 KL 和 KP、且能用的小鼠 n；不与私有 8KL 合并。GSE164758 里的 KPL vs KP 是三重突变，不是 KL。

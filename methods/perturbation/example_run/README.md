# Worked example — GSE207704 (CLDN4 knockout RNA-seq)

**中文摘要在下方 / Chinese summary below.**

This is one *runnable* companion to [`../playbook.md`](../playbook.md) (CLDN4 KO).
Sibling public TACSTD2 runs: [`../example_gse245459/`](../example_gse245459/)
(SKOV3 shTACSTD2) and [`../example_gse334497/`](../example_gse334497/) (4T1 Trop2 KO).
All accessions were verified live; none are invented.

| Field | Value |
|---|---|
| Accession | **GSE207704** (verified via NCBI E-utilities) |
| Title | *Claudin-4–adhesion signaling drives breast cancer metabolism and progression via liver X receptor β* |
| Perturbation | **CLDN4 CRISPR knockout** vs wild-type |
| Models | Human breast cancer lines **MCF7** and **T47D** |
| Assay | Expression profiling by high-throughput sequencing (RNA-seq) |
| Processed file | `GSE207704_CLDN4_RNAseq.txt.gz` (~1.0 MB, Cufflinks **FPKM**) |

Run it:

```bash
pip install pandas numpy gseapy
python run_gse207704.py
```

Outputs land in `results/` (already committed): `de_log2fc_table.tsv`,
`focus_genes.tsv`, `gsea_prerank_report.tsv`, `rank.rnk`, `run_summary.txt`.

## Important caveat (this is a teaching point, not a defect)

The **processed** matrix GEO distributes for this series is **FPKM with exactly
one column per condition per cell line** — i.e. **no within-group replicates**.
A negative-binomial DE test (DESeq2/edgeR) needs **raw integer counts** and
**≥2 replicates per group**, so we deliberately **do not** fabricate p-values
from it. Instead this example demonstrates what a replicate-free FPKM matrix
*legitimately* supports, and the DESeq2/edgeR/PyDESeq2 templates in
[`../templates/`](../templates/) show the full statistical pipeline for when you
have raw counts (which here would require re-quantifying the SRA FASTQs — out of
the <2 GB budget).

## What the real run shows

Ranking metric = mean of per-cell-line `log2((KO_FPKM+1)/(WT_FPKM+1))`.

1. **On-target check (sanity):** `CLDN4` mean log2FC = **−0.88** (down in the KO,
   concordant in both MCF7 and T47D). A knockout should drive its own transcript
   down; this passes. (It is only ~2-fold, not −∞, because a CRISPR frameshift
   KO can still transcribe an NMD-targeted mRNA and Cufflinks FPKM at the
   locus level dampens the estimate — a normal, expected observation.)

2. **"Opposite gene" hypothesis — CLDN4↓ → TACSTD2/TROP2↑?**
   **Not supported in this dataset.** `TACSTD2` mean log2FC = **−0.81** — TROP2
   goes **down together with CLDN4**, concordantly in both lines. This is the
   whole point of the check: a plausible compensation story must be *tested*,
   and here the data falsify it. (`EPCAM`, the TROP2 paralog, is mildly **up**,
   +0.22.) See [`../templates/opposite_gene_check.py`](../templates/opposite_gene_check.py)
   for the reusable, p-value-aware version.

3. **GSEA (pre-ranked, Hallmark + KEGG):** positive NES = up in KO.

   | Gene set | NES | NOM p | FDR q |
   |---|---:|---:|---:|
   | Hallmark Interferon Alpha Response | **−1.66** | 0.007 | **0.007** |
   | Hallmark Interferon Gamma Response | **−1.46** | 0.013 | **0.035** |
   | Hallmark Epithelial-Mesenchymal Transition | −1.24 | 0.10 | 0.14 |
   | KEGG Tight junction | −1.23 | 0.11 | 0.12 |

   The **interferon-α/γ** signatures are **significantly down** in CLDN4-KO
   (FDR < 0.05); EMT and tight-junction sets trend down but are not significant
   with this replicate-free ranking. Treat these as *hypothesis-generating*
   given the design limitation.

## 中文摘要

本目录是 [`../playbook.md`](../playbook.md) 的**可运行**示例，使用**真实且经核实**的
GEO 数据集 **GSE207704**（CLDN4 CRISPR 敲除 vs 野生型，人乳腺癌 MCF7/T47D，RNA-seq；
未编造任何 accession）。

**关键提醒（这是教学点，不是缺陷）：** GEO 提供的处理后矩阵是 **FPKM**，每种条件每个
细胞系**只有 1 列、没有组内重复**。负二项检验（DESeq2/edgeR）需要**原始整数 counts**
且**每组 ≥2 个重复**，因此本示例**不伪造 p 值**，而是演示无重复 FPKM 矩阵**合理支持**
的分析；完整统计流程见 [`../templates/`](../templates/)（需从 SRA 重新定量原始 counts，
超出 <2GB 预算）。

真实结果：(1) **靶基因验证**：`CLDN4` 平均 log2FC = **−0.88**（敲除后下调，两细胞系一致），
sanity check 通过。(2) **“反向基因”假设 CLDN4↓→TACSTD2↑ 在本数据中不成立**：`TACSTD2`
= **−0.81**，与 CLDN4 一起下调——说明这类“代偿”假设必须用数据检验，此处被证伪。
(3) **GSEA**：干扰素 α/γ 特征在敲除中**显著下调**（FDR<0.05）；EMT 与紧密连接集下调但不显著。
因设计限制，以上结论应作为**假设生成**看待。

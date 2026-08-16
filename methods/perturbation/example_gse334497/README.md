# Worked example — GSE334497 (4T1 Trop2 KO RNA-seq, mouse)

**中文摘要在下方.**

Verified public GEO series. No invented accession.

| Field | Value |
|---|---|
| Accession | **GSE334497** |
| Title | *TROP2/claudin program mediates immune exclusion to impede checkpoint blockade in breast cancer* |
| Perturbation | **Trop2 / Tacstd2 knockout** vs WT |
| Model | Mouse **4T1** tumors grown in vivo, **n=5 vs 5** |
| Assay | RNA-seq |
| Processed file | `GSE334497_normalized_counts.csv.gz` (~1.2 MB, **normalized** non-integer counts) |
| IDs | mouse Ensembl → MGI symbols via mygene.info (cached in `data/ensembl_to_mgi.json`) |

Sample mapping (from `!Sample_description` / `!Sample_title` in the series matrix):

| Group | Library names in the count table |
|---|---|
| KO | KO162, KO164, KO165, KO172, RESUB-KO163R |
| WT | control170, RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R |

```bash
pip install pandas numpy scipy gseapy
python run_gse334497.py
```

## Caveat

The supplement is **normalized counts**, not raw integers → **not DESeq2/edgeR**.
Welch t + BH is exploratory. GSEA uses `toupper(MGI)` against a **human** GMT
(playbook §9 limitation; acceptable for these conserved Hallmark/KEGG sets).

In-vivo tumors have high sample variance. `Tacstd2` mean drops 159 → 11
(log2FC −3.82) so the KO is real at the effect-size level, but Welch padj is
**not** significant (0.65) because WT tumors span 66–361. That is the n-small /
heterogeneity lesson: **rank and GSEA still work; a single-gene p-value may not.**

## Real results (KO vs WT)

| Check | Gene | log2FC | padj | Verdict |
|---|---|---:|---:|---|
| On-target | Tacstd2 | **−3.82** | 0.65 | mean drop is large; p is not |
| Opposite-gene (expect UP) | Cldn4 | **−0.82** | 0.94 | **DOWN in mean — not supported** |

GSEA (positive NES = up in KO):

| Set | NES | FDR |
|---|---:|---:|
| Hallmark IFN-γ | **+1.97** | ~0 |
| Hallmark IFN-α | **+1.67** | ~0 |
| KEGG Tight junction | **−1.90** | ~0 |
| Hallmark EMT | +0.95 | 0.62 |

IFN programs go **up** in this in-vivo mouse KO — the opposite direction from the
two human *cell-line* perturbations (GSE207704 CLDN4 KO, GSE245459 shTACSTD2).
Do not average across models; report each design.

## 中文摘要

**GSE334497** 为公开的小鼠 4T1 **Trop2/Tacstd2 敲除** RNA-seq（体内肿瘤，n=5）。
补充文件是**标准化 counts**（非整数），本运行**不是** DESeq2/edgeR。`Tacstd2` 均值
159→11（log2FC −3.82）说明敲除在效应量上成立，但组内方差大，Welch padj 不显著——
这正是小样本/异质性课：排序与 GSEA 仍可用，单基因 p 值未必。`Cldn4` 均值下调，
「代偿性上调」假设不成立。GSEA：干扰素 α/γ 在敲除中**上调**（与两个人细胞系扰动
方向相反），紧密连接下调。按设计分别报告，不要跨模型平均。

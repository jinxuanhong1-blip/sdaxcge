# Worked example — GSE245459 (SKOV3 shTACSTD2 RNA-seq)

**中文摘要在下方.**

Verified public GEO series (NCBI E-utilities + series matrix). No invented accession.

| Field | Value |
|---|---|
| Accession | **GSE245459** |
| Title | *Exploring the Mechanism of TACSTD2 Causing Cisplatin Resistance* |
| Perturbation | **shTACSTD2** (`genotype: sh`) vs **shNC** |
| Model | Human ovarian cancer **SKOV3**, n=3 per group |
| Extra arm | same KD under cisplatin (DDP): shDDP vs shNCDDP |
| Assay | RNA-seq (Illumina NovaSeq 6000) |
| Processed file | `GSE245459_fpkm.anno.txt.gz` (~14 MB **FPKM**) |

```bash
pip install pandas numpy scipy gseapy
python run_gse245459.py
```

## Caveat

GEO ships **FPKM**, not raw integer counts. This run is **not DESeq2/edgeR**.
Welch t-tests on `log2(FPKM+1)` + BH FDR are **exploratory**. Use the R/Python
count-model templates in `../templates/` if you re-quantify the SRA FASTQs.

## Real results (primary: sh vs shNC)

| Check | Gene | log2FC | padj | Verdict |
|---|---|---:|---:|---|
| On-target | TACSTD2 | **−2.63** | 0.004 | KD worked (5.40 → 0.03 FPKM) |
| Opposite-gene (expect UP) | CLDN4 | **−1.92** | 0.004 | **DOWN — hypothesis not supported** |

Other claudins (`CLDN1/3/7`) and `EPCAM` also fall. Hallmark **IFN-α / IFN-γ / EMT**
are down in the vehicle KD (GSEA FDR ≈ 0); KEGG tight junction is not significant
(FDR 0.20).

Under **cisplatin** the on-target KD is even stronger (TACSTD2 −4.03) but CLDN4
is already low in shNC+DDP, so the opposite-gene test is uninformative (log2FC
−0.14, padj 0.25). IFN signatures **flip up** in that arm — report the contrast,
do not pool.

## 中文摘要

**GSE245459** 为公开的 SKOV3 **shTACSTD2** RNA-seq（每组 n=3；另有顺铂臂）。GEO
补充文件是 **FPKM**，本运行**不是** DESeq2/edgeR。靶基因 `TACSTD2` 下调 −2.63
（padj 0.004），敲低成立；反向假设「TACSTD2↓ → CLDN4↑」**不成立**（CLDN4 −1.92）。
无药条件下 Hallmark 干扰素 α/γ 与 EMT 下调。

# Ku70/Ku80 loss and IFN / STING / antigen presentation

Public GEO and ENCODE RNA-seq only. Differential expression is DESeq2 (Wald statistic, KO/KD over the matched control). Enrichment is pre-ranked GSEA (1,000 gene permutations, weight 1, seed 123). Positive NES means the set sits among genes that rise after Ku loss.

Ku genes (`XRCC5`, `XRCC6`) and `PRKDC` were removed from every set before testing. Full tables are `gsea_prerank.tsv`, `contrast_inventory.tsv`, `focus_genes.tsv`, and `meta_stouffer.tsv`. Search notes are in `notes/ku70_ku80_ifn_gsea/dataset_inventory.md`.

## What was tested

Six Ku-loss contrasts enter the primary meta. Each is a separate biological system, and each passed an on-target check (mRNA down, except the AID degron, which removes protein):

| Contrast | Perturbation | On-target log2FC (FDR) |
|---|---|---|
| GSE294709 HCT116 Ku80-AID, CDKi day 4, 3 clones | Ku80 protein degron, proliferation-matched | XRCC5 −0.22 (FDR 0.24); XRCC6 −0.61 (FDR 7.6×10⁻⁵) |
| GSE294709 HCT116 Ku80 flox, Cre × 4OHT interaction | Ku80 excision above the 4OHT-only effect | XRCC5 −3.19 (FDR 8.4×10⁻¹⁸⁴) |
| GSE294709 HEK293 Ku70-KO, dox withdrawal, 3 clones | ectopic Ku70 removed | XRCC6 −2.74 (FDR 1.4×10⁻¹⁰⁹) |
| GSE180581 HEK293T Ku70+/− siKu70, n=3 | siRNA | XRCC6 −2.01 (FDR 2.1×10⁻⁵⁰) |
| GSE180581 HEK293T Ku80+/− siKu80, n=3 | siRNA | XRCC5 −1.99 (FDR 3.2×10⁻⁵⁰) |
| GSE247031 mouse Treg Xrcc6 cKO, LLC lung tumors, n=3 | genetic knockout | Xrcc6 −3.07 (FDR 1.9×10⁻⁵) |

Gene sets: MSigDB Hallmark interferon-alpha and interferon-gamma; Reactome STING-mediated induction (Ku/DNA-PKcs removed; 13 genes); a curated MHC-I antigen-presentation set (HLA-A/B/C/E/F/G, B2M, TAP1/2, TAPBP, TAPBPL, PSMB8/9/10, PSME1/2, NLRC5, ERAP1/2); KEGG antigen processing, kept as a broad set; Reactome DDX58/IFIH1 interferon induction, kept as the RNA-sensor comparator.

The same HCT116 AID clones were also tested at 24 h, day 2, and CDKi 48 h. Those arms repeat one system, so they are a time course, not extra meta units. ENCODE shRNA and CRISPR (K562 and HepG2, n=2 vs matched non-targeting controls, GRCh38 V29 expected counts) are a second meta. All seven ENCODE contrasts reduced the targeted mRNA.

## Result

**HCT116 Ku80 loss induces interferon and classical MHC-I genes. Other Ku-loss systems do not.**

In proliferation-matched HCT116 Ku80-AID cells (CDKi + IAA versus CDKi, day 4):

- Hallmark IFN-α NES 3.36, FDR 0. Hallmark IFN-γ NES 3.26, FDR 0.
- ISG15 log2FC +3.35 (FDR 3.4×10⁻⁶⁹), IFIT1 +4.22 (FDR 1.5×10⁻¹⁷⁵), STAT1 +2.07 (FDR 6.7×10⁻⁵⁸), DDX58 +2.40 (FDR 1.8×10⁻⁵⁷), IFIH1 +2.17 (FDR 5.6×10⁻⁴⁴).
- MHC-I set NES 2.51, FDR 0. Leading edge: HLA-A, HLA-C, HLA-B, TAP1, TAP2, HLA-E, HLA-F, B2M, NLRC5. HLA-A +1.57 (FDR 1.4×10⁻³⁵), TAP1 +1.36 (FDR 3.4×10⁻¹⁹), B2M +0.77 (FDR 6.1×10⁻¹⁰).
- STING1 log2FC −0.19 (FDR 0.56). TBK1 −0.05. The Reactome STING set NES is 1.45 (FDR 0.043) with leading edge TRIM21 and IRF3 only.

The same direction is already present at 24 h (IFN-α NES 1.90, FDR 0.002; MHC-I NES 1.59, FDR 0.022) and at day 2 and CDKi 48 h.

The floxed HCT116 system agrees after the 4OHT control is subtracted. 4OHT without Cre leaves XRCC5 unchanged (log2FC −0.04) and leaves ISG15 flat (−0.11), but the broad IFN-α hallmark still moves (NES 1.91, FDR 0.001; IFIT1 +0.77, FDR 0.093). The interaction contrast (4OHT effect in CreERT2 minus 4OHT effect without Cre) keeps the Ku-specific signal: IFN-α NES 2.63, FDR 0; MHC-I NES 2.13, FDR 0; ISG15 +2.18 (FDR 2.1×10⁻⁵⁰); TAP1 +1.28 (FDR 3.2×10⁻²¹); HLA-A +0.85 (FDR 4.1×10⁻¹⁴). STING1 stays flat (−0.16, FDR 0.56).

On the MAVS-knockout Ku80-AID background, the ISG spike shrinks: ISG15 +0.29 (FDR 0.25), IFIT1 −0.31 (FDR 0.46), IFIH1 +0.11 (FDR 0.84). HLA-A remains up (+1.18, FDR 1.8×10⁻¹⁴), and the MHC-I set NES is 2.32 (FDR 0). The hallmark IFN-α NES is still 1.75 (FDR 0.006), with a different leading edge. Canonical ISG15/IFIT1 induction is the MAVS-sensitive part.

**The other four primary systems do not reproduce that induction.**

- HEK293 Ku70 withdrawal: IFN-α NES −1.32 (FDR 0.10). MHC-I NES −1.69 (FDR 0.007). HLA-A −0.44 (FDR 0.0026), B2M −0.43 (FDR 0.0025). STING1 −1.25 (FDR 1.0×10⁻⁵). ISG15 −0.34 (FDR 0.36).
- HEK293T siKu70: IFN-α NES −0.95 (FDR 0.55). MHC-I NES −1.50 (FDR 0.12). ISG15 −0.70 (FDR 0.11).
- HEK293T siKu80: IFN-α NES −0.65 (FDR 0.96). MHC-I NES −0.81 (FDR 1). B2M −0.64 (FDR 0.0014).
- Treg Xrcc6 knockout: IFN-α NES 1.10 (FDR 0.29). MHC-I NES −0.71 (FDR 1). ISG15 +0.22 (FDR 0.99), Tap1 +0.03 (FDR 1.0).

KEGG antigen processing is positive in the siRNA and Treg contrasts (NES 1.55–1.76, FDR ≤ 0.04). The leading-edge genes there are heat-shock and chaperone genes (HSPA1B, HSPA5, CALR) or, in Tregs, MHC-II and lineage genes (H2-Eb1, Cd74, Ciita, Cd4, Cd8b1). That is why the curated MHC-I set is the APM readout used below.

DNA-PKcs siRNA in the same HEK293T experiment (not a Ku gene) lowers the IFN-α hallmark (NES −1.52, FDR 0.043) and does not raise MHC-I (NES −1.16, FDR 0.22).

ENCODE n=2 contrasts are mixed. IFN-α NES is positive in 4/7 and negative in 3/7. MHC-I NES is positive in 1/7. None of that meta is concordant (Stouffer one-sided p for IFN-α up = 0.19; for MHC-I up = 0.33).

## Meta across the six primary contrasts

Stouffer combination of the GSEA one-sided p-values for enrichment upward, weighted by sqrt(n). This p-value is a gene-permutation p conditional on each ranking. It is pulled by the two very large HCT116 effects. Concordance is the number that describes agreement.

| Set | Contrasts with NES > 0 | Contrasts with NES > 0 and FDR < 0.25 | Median NES | Stouffer p (up) |
|---|---|---|---|---|
| IFN-α | 3/6 | 2/6 | +0.22 | 1.7×10⁻³ |
| IFN-γ | 2/6 | 2/6 | −0.92 | 1.2×10⁻² |
| STING Reactome | 5/6 | 3/6 | +1.09 | 1.8×10⁻² |
| MHC-I APM | 2/6 | 2/6 | −0.76 | 0.11 |
| KEGG antigen processing | 5/6 | 5/6 | +1.65 | 2.8×10⁻⁵ |
| RIG-I/MDA5 | 3/6 | 3/6 | +0.34 | 5.1×10⁻³ |

The two IFN-α and MHC-I successes are the two HCT116 Ku80 contrasts. STING1 itself does not rise in those contrasts; the Reactome STING leading edges are TRIM21, IRF3, STAT6, or DDX41. A positive STING NES here is not evidence that cGAS–STING transcription turned on.

## Bottom line

Public Ku80 depletion in HCT116 induces an interferon program and classical MHC-I genes (HLA, TAP, B2M). The ISG15/IFIT1 spike shrinks on a MAVS-knockout background. STING1 transcript levels stay flat. The same Ku70 or Ku80 loss in HEK293, HEK293T, and mouse tumor Tregs does not induce that program; HEK293 Ku70 withdrawal lowers MHC-I. ENCODE K562/HepG2 knockdowns are inconsistent. The interferon/MHC-I result is real in HCT116 Ku80-loss RNA-seq and is not a general consequence of losing Ku.

## 结论（简）

HCT116 里敲掉 Ku80，干扰素和 MHC-I 抗原呈递基因上升（IFN-α NES 3.36；ISG15 log2FC +3.35；HLA-A +1.57）。MAVS 敲除背景下 ISG15/IFIT1 这支明显变弱。STING1 转录不变。HEK293/HEK293T 的 Ku70 或 Ku80 敲低，以及小鼠 Treg 的 Xrcc6 敲除，不重复这套上升；HEK293 撤掉 Ku70 后 MHC-I 反而下降。所以这是 HCT116 Ku80 缺失的公共转录组结果，不是所有 Ku 缺失的普遍现象。

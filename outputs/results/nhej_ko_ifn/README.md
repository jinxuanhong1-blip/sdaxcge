# NHEJ-gene knockout RNA-seq: interferon scores

Search date: **2026-09-21**. GEO DataSets only. Direction is **knockout minus matched control**. A positive score means higher interferon-stimulated gene expression in the knockout.

## What was scored

Prespecified primary signature: a cell-intrinsic type-I ISG core (`ISG_CORE` in `scripts/nhej_ko_ifn/analyze_nhej_ko_ifn.py`). Mouse stand-ins are Oas1a, Oasl2, Herc6, and Ifi27l2a. IFI44L and IFI6 are human-only.

Secondary signatures: MSigDB Hallmark interferon-alpha (97 genes) and interferon-gamma (200 genes), Enrichr `MSigDB_Hallmark_2020`, stored in `scripts/nhej_ko_ifn/hallmark_ifn.tsv`. Mouse hallmark symbols are title-case orthologs with the same stand-ins. Genes that are absent from the matrix are dropped.

For each contrast, a gene enters the score only if its maximum log2(abundance+1) in the contrast samples is at least 1 (CPM or FPKM at least 1). The sample score is the unweighted mean of those log2 values. The effect is the difference of arm means. The p-value is a two-sided Welch test on the sample scores. Genes are not treated as independent replicates. One library per arm has no p-value.

The knocked-out gene itself is a QC measurement, not part of the score. `rna_down` means its log2 delta is −0.5 or lower. `near_floor_both` means both arms sit below CPM/FPKM of 1. `rna_not_reduced` means the matrix does not show loss of that mRNA.

## Cancer-line result

There is no single interferon effect of “NHEJ knockout.” The four genes do not point the same way, so they are not pooled.

| Gene | Contrast | n | ISG-core Δ | Welch p | Genes up | KO mRNA | Call |
|---|---|---:|---:|---:|---:|---|---|
| XRCC4 | HeLa, no mirin, GSE135274 | 2 vs 2 | −0.95 | 0.11 | 2/27 | down (Δ −0.74) | lower, not significant |
| LIG4 | CH12F3 lymphoma, GSE154443 | 1 vs 1 | −0.34 | — | 5/16 | both near the floor | not interpretable |
| PRKDC | HCT116, normoxia, GSE285698 | 3 vs 3 | +0.56 | 8.9×10⁻⁴ | 17/19 | down (Δ −1.95) | higher |
| TP53BP1 | MCF-7, untreated, GSE84986 | 3 vs 6 | +0.71 | 0.0014 | 23/23 | down (Δ −0.94) | higher |

Hallmark interferon-alpha and interferon-gamma move in the same direction as the core in each of these four rows. Full precision is in `tables/primary_cancer_line_contrasts.tsv`.

**PRKDC / HCT116.** DNA-PK knockout raises the ISG core under vehicle and under CoCl2 hypoxia (hypoxia Δ +0.48, p = 2.2×10⁻⁴). The largest core increases are CMPK2, IRF7, IFIT1, IFIT2, and IFI6 (about +1.0 to +1.4 log2). ISG15 (−0.67) and IFI27 (−0.19) do not follow the set. This is a consistent moderate shift, not a uniform induction of every ISG.

**TP53BP1 / MCF-7.** Both null clones are above every wild-type replicate. Clone 1 alone is Δ +0.68 (p = 0.064) and clone 2 alone is Δ +0.75 (p = 0.012). The pooled p-value treats six knockout libraries as independent even though the clones share a parent, so 0.0014 is sharper than the clone-level tests. The 4-hour 5 Gy arm is similar (Δ +0.69, p = 8.9×10⁻⁴). Largest untreated increases: IFI6, IFIT1, and MX1 (about +1.6 to +1.9 log2).

**XRCC4 / HeLa.** Both experiments sit below both controls (sample scores about 4.0–4.1 versus 4.8–5.2). Twenty-five of 27 detected core genes are lower, median gene log2 fold-change −1.14. With two samples per arm the Welch p-value stays 0.11. XRCC4 mRNA is down. Mirin does not reverse the decrease. This is a large, consistent point estimate that is not a significant sample-level test.

**LIG4 / CH12F3.** The deposited file has one FPKM column per genotype. Lig4 itself is 0.13 FPKM in the parental column and 0.64 FPKM in the knockout column, so both are at the floor and the RNA-seq does not confirm the knockout. Isg15 is lower in the knockout column (109 versus 48 FPKM). The core is not coherently up. The GEO summary says Lig4 loss triggers cytosolic DNA-sensing and RIG-I-like pathways; that claim is not what this prespecified core shows in the deposited matrix.

## Scored, not in the cancer-line table

| Accession | Why it is not in the primary table | ISG-core Δ | KO mRNA |
|---|---|---:|---|
| GSE145148 MCF10A XRCC4 KO vs WT | MCF10A is not a cancer line | −1.59 (p = 0.026) | not reduced (Δ +0.08) |
| GSE145148 MCF10A XRCC4/TP53 DKO vs TP53 KO | same, non-malignant | −1.72 (p = 7.1×10⁻⁴) | not reduced (Δ +0.02) |
| GSE280049 hTERT-RPE1 sg53BP1 | not a cancer line | −0.27 (p = 0.88) | higher in the sg53BP1 columns (Δ +1.55) |
| GSE237615 ID8 53BP1-KO tumors | bulk allograft, immune cells included | +1.20 (p = 7.8×10⁻⁴) | Trp53bp1 down (Δ −1.21) |
| GSE237615 KPC 53BP1-KO tumors | bulk allograft | +0.37 (p = 0.26) | Trp53bp1 down (Δ −0.59) |

MCF10A interferon genes are lower in the samples labeled XRCC4 knockout, and TP53 mRNA does fall in the TP53-knockout samples, but XRCC4 mRNA does not. Those interferon numbers are not evidence that XRCC4 loss caused them.

RPE1 is not usable as a 53BP1-loss contrast. The columns named RPE-KO have more TP53BP1 counts (1743 and 1516) than the columns named RPE-NC (431 and 594). One control is also an interferon outlier (STAT1 counts 24192 versus 3681).

ID8 tumors are the strongest interferon increase in this screen, and Trp53bp1 mRNA is down, but the RNA is from implanted tumors. Infiltrate can produce that signature. KPC tumors move less: the core is not significant, while Hallmark interferon-gamma is Δ +0.53 (p = 0.042). Neither tumor row is a cell-intrinsic cancer-line measurement.

## Search

Queries, the union of 217 GEO series, and the series that were read and rejected are in `search_log.md` and `dataset_inventory.tsv`. Knockdown, inhibitor, and non-cancer series were not promoted into the primary table. No second replicated cancer-line RNA-seq knockout was found for any of the four genes, so there is no within-gene random-effects pool.

## Reproduce

```bash
python3 scripts/nhej_ko_ifn/analyze_nhej_ko_ifn.py
```

Dependencies are `pandas`, `numpy`, `scipy`, and `openpyxl`. The script downloads GEO supplementary files into `outputs/results/nhej_ko_ifn/_cache/` (git-ignored) and rewrites the tables.

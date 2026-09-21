# FINDING — GSE334497 Trop2 KO vs a CLDN4-KD direction

**Additive public evidence.** Wu *et al.*, *Journal for ImmunoTherapy of Cancer* 2026;14:e012265 ([JITC](https://jitc.bmj.com/content/14/4/e012265)). GEO [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497). CRISPR **Trop2 (Tacstd2)** knockout versus wild-type **4T1** tumors, 5 vs 5, grown 3 weeks in BALB/c, RNA from **frozen whole-tumor sections**. Breast, not lung. This is the depositing paper’s matrix, not an independent cohort, and it is **not a CLDN4 knockdown**.

## Answer

**No.** Trop2 KO does not match a CLDN4 knockdown. *Cldn4* is lower on average, and the 5-vs-5 test does not clear 0.05. The epithelial interferon genes stay flat. Bulk immune RNA does rise.

| Question | Result |
|---|---|
| Did Trop2 KO work? | Yes. *Tacstd2* log2FC **-3.821**, one-sided Welch *p* = **0.00054**. |
| Does *Cldn4* drop? | Direction only. log2FC **-0.817** (95% CI -2.37 to +0.73), one-sided Welch *p* = **0.126**, two-sided *p* = 0.253, MWU *p* = 0.421. Arm A does not pass. |
| Does epithelial IFN rise? | **No.** Epithelial ISG Δ **+0.075**, exact one-sided perm *p* = **0.452**. Arm B does not pass. |
| Does the full Hallmark IFN-γ score rise? | **Not at the sample level.** Δ **+0.287**, *d* = +0.83, exact perm *p* = **0.107**. |
| What part of IFN-γ does rise? | The infiltrate-leaning slice. Δ **+0.692**, perm *p* = **0.004**. The ISG slice inside the same hallmark does not (Δ +0.068, *p* = 0.448). |
| Does bulk immune RNA rise? | **Yes.** Δ **+0.934**, *d* = +1.86, exact perm *p* = **0.012**. |

Joint rule (both required): *Cldn4* down at one-sided *p* < 0.05, and epithelial ISG up at one-sided exact perm *p* < 0.05. Joint = **no**.

A prerank gene-set statistic can still put Hallmark IFN-γ at the KO end of the ranking, because *Cxcl9*, MHC-II, and other infiltrate genes sit at that end. That is a statement about gene order. It is not the same as the 10 tumors separating on the mean of all ~184 Hallmark genes. The pre-specified set test here is the exact sample permutation. On that test the full Hallmark score does not pass 0.05, and the epithelial ISG score is null. The score that passes is the bulk immune module, which is the paper’s T-cell infiltration result read out of the same frozen sections. OXPHOS, run with the same permutation, stays null.

## Design

| Item | Choice |
|---|---|
| Matrix | `GSE334497_normalized_counts.csv.gz` (author-normalized counts). No FASTQ. |
| Groups | KO: KO162, KO164, KO165, KO172, RESUB-KO163R. WT: RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R, control170. Library names are GEO `Sample_description`. |
| Transform | log2(normalized count + 1). Ensembl → NCBI symbol (`gene2ensembl` + `Mus_musculus.gene_info`). Duplicate symbols: keep the row with the higher mean count. |
| Gene test | Welch *t* on the 5 vs 5 log2 values. One-sided *p* is the pre-specified arm; two-sided and Mann–Whitney are reported with it. |
| Set test | Per-sample mean of gene-wise z-scores. Exact permutation of the 252 equal-sized label splits. One-sided up = KO higher. |
| Arm B genes | Isg15, Ifit1, Ifit2, Ifit3, Mx1, Mx2, Oas2, Oas3, Oasl1, Rsad2, Stat1, Irf7, Ifi27, Usp18, Bst2, Ifih1. |
| Bulk immune | Cd3e, Cd3d, Cd8a, Cd8b1, Cd4, Gzmb, Gzma, Prf1, Ifng, Nkg7, Klrd1, Cxcl9, Cxcl10. |
| Control | OXPHOS gene score, same permutation, expected null. |
| α | 0.05, pre-specified arms only. No genome-wide FDR claim. |

## Perturbation QC and the Cldn4 arm

| Gene | log2FC | 95% CI | two-sided *p* | one-sided down *p* | *d* |
|---|---:|---|---:|---:|---:|
| Tacstd2 | -3.821 | -5.52, -2.12 | 0.001 | 0.00054 | -3.35 |
| Cldn4 | -0.817 | -2.37, +0.73 | 0.253 | 0.126 | -0.79 |
| Cldn1 | -2.181 | -4.16, -0.20 | 0.035 | 0.017 | -1.61 |
| Cldn7 | -0.823 | -1.84, +0.19 | 0.091 | 0.045 | -1.34 |
| Ocln | -0.661 | -1.43, +0.11 | 0.078 | 0.039 | -1.40 |
| Tjp1 | -0.074 | -1.15, +1.00 | 0.875 | 0.437 | -0.10 |
| Epcam | -1.262 | -2.27, -0.25 | 0.021 | 0.011 | -1.89 |

*Cldn4* is abundant on both sides (mean normalized count KO 1674, WT 2605). The ten samples overlap, and the Welch interval crosses zero. The paper’s stated RNA mediators are **claudin 1, claudin 7, and occludin**, not claudin 4. *Cldn1* is the large drop (log2FC -2.181, two-sided *p* = 0.035, one-sided down *p* = 0.017, MWU *p* = 0.032). *Cldn7* has almost the same log2FC as *Cldn4* and separates more cleanly on the rank test (log2FC -0.823, two-sided *p* = 0.091, one-sided down *p* = 0.045, MWU *p* = 0.016). *Ocln* is smaller and one-sided only (log2FC -0.661, two-sided *p* = 0.078, one-sided down *p* = 0.039, MWU *p* = 0.095). The family table is below; the match arm is *Cldn4* only.

## IFN / immune scores

Positive Δ = higher in Trop2 KO. Permutation *p* is exact (252 splits).

| Set | n present | Δ | *d* | perm *p* up | perm *p* two | passes up |
|---|---:|---:|---:|---:|---:|---|
| Epithelial ISG (arm B) | 16/16 | +0.075 | +0.11 | 0.452 | 0.905 | no |
| Hallmark IFN-γ | 184/188 | +0.287 | +0.83 | 0.107 | 0.214 | no |
| Hallmark IFN-γ, infiltrate-leaning | 35/36 | +0.692 | +2.10 | 0.004 | 0.008 | yes |
| Hallmark IFN-γ ∩ epithelial ISG | 14/14 | +0.068 | +0.10 | 0.448 | 0.897 | no |
| Hallmark IFN-γ remainder | 135/138 | +0.204 | +0.58 | 0.187 | 0.373 | no |
| Bulk immune | 13/13 | +0.934 | +1.86 | 0.012 | 0.024 | yes |
| OXPHOS control | 16/16 | -0.235 | -0.28 | 0.671 | 0.667 | no |

Epithelial ISG genes, log2FC (KO − WT):

| Gene | log2FC | 95% CI | two-sided *p* | one-sided up *p* | *d* |
|---|---:|---|---:|---:|---:|
| Isg15 | -0.091 | -1.09, +0.91 | 0.839 | 0.581 | -0.13 |
| Ifit1 | -0.021 | -0.68, +0.64 | 0.943 | 0.529 | -0.05 |
| Mx1 | +0.160 | -0.70, +1.02 | 0.666 | 0.333 | +0.29 |
| Oas2 | -0.104 | -1.05, +0.84 | 0.805 | 0.597 | -0.16 |
| Stat1 | +0.460 | -0.54, +1.46 | 0.287 | 0.144 | +0.75 |
| Irf7 | -0.183 | -1.04, +0.67 | 0.617 | 0.692 | -0.33 |
| Ifi27 | +0.222 | -0.39, +0.84 | 0.416 | 0.208 | +0.55 |

Bulk immune genes:

| Gene | log2FC | 95% CI | two-sided *p* | one-sided up *p* | *d* |
|---|---:|---|---:|---:|---:|
| Cxcl9 | +1.082 | +0.41, +1.76 | 0.006 | 0.003 | +2.34 |
| Cxcl10 | +0.564 | -0.05, +1.18 | 0.067 | 0.034 | +1.36 |
| Cd8a | +0.674 | -0.19, +1.54 | 0.108 | 0.054 | +1.15 |
| Cd3e | +0.443 | -0.52, +1.41 | 0.322 | 0.161 | +0.67 |
| Gzmb | +0.763 | -0.76, +2.28 | 0.260 | 0.130 | +0.79 |
| Prf1 | +1.034 | +0.04, +2.03 | 0.043 | 0.021 | +1.52 |
| Ifng | +0.889 | -0.34, +2.12 | 0.132 | 0.066 | +1.08 |
| Nkg7 | +1.021 | -0.00, +2.04 | 0.050 | 0.025 | +1.48 |
| Cd274 | +0.772 | -0.22, +1.77 | 0.105 | 0.052 | +1.23 |

*Cxcl9* is the clearest single immune gene. *Cd274* (PD-L1) is in the table because it sits in Hallmark IFN-γ; it is not a pre-specified arm.

## Claudin family (every *Cldn* present)

| Gene | log2FC | 95% CI | two-sided *p* | one-sided down *p* | *d* |
|---|---:|---|---:|---:|---:|
| Cldn1 | -2.181 | -4.16, -0.20 | 0.035 | 0.017 | -1.61 |
| Cldn3 | -1.438 | -3.11, +0.23 | 0.082 | 0.041 | -1.27 |
| Cldn7 | -0.823 | -1.84, +0.19 | 0.091 | 0.045 | -1.34 |
| Cldn4 | -0.817 | -2.37, +0.73 | 0.253 | 0.126 | -0.79 |
| Cldn15 | -0.747 | -2.07, +0.57 | 0.220 | 0.110 | -0.86 |
| Cldn23 | -0.698 | -1.51, +0.12 | 0.082 | 0.041 | -1.31 |
| Cldn2 | -0.284 | -1.70, +1.13 | 0.652 | 0.326 | -0.30 |
| Cldnd1 | -0.246 | -0.63, +0.14 | 0.161 | 0.081 | -1.03 |
| Cldn12 | -0.245 | -0.48, -0.01 | 0.043 | 0.022 | -1.55 |
| Cldn20 | +0.057 | -0.95, +1.06 | 0.899 | 0.551 | +0.08 |
| Cldn5 | +0.317 | -0.76, +1.40 | 0.515 | 0.743 | +0.43 |
| Cldn10 | +0.568 | -0.31, +1.44 | 0.159 | 0.920 | +1.03 |

## What this is and is not

**Is**

- A pre-specified match test of “Cldn4 down and epithelial IFN up” on the public 4T1 Trop2 KO counts.
- A split of Hallmark IFN-γ. The infiltrate-leaning slice is higher in KO. The epithelial ISG slice is flat. The full 184-gene score falls between them and does not pass 0.05.
- Consistent with Wu *et al.*: Trop2 loss goes with less tight-junction RNA and more T-cell / inflammatory RNA in the bulk tumor.

**Is not**

- A CLDN4 knockdown or knockout.
- A lung model.
- Tumor-cell-intrinsic IFN (no sorted epithelium and no in-vitro 4T1 arm in this accession).
- An independent replication of the JITC paper. These are their tumors.
- A genome-wide discovery list. At n = 5 vs 5, gene-level tests are the pre-specified ones above.

Figure: `figures/fig_match.png`.

Reproduce: `python3 scripts/gse334497_cldn4_match/analyze.py`

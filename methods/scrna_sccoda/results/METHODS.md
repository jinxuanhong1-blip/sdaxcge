# Methods — public scRNA compositional analysis

## Series (public GEO only)

| Accession | Paper | What was used | Annotation |
|---|---|---|---|
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | Hu et al. 2023, neoadjuvant PD-1 + chemo, 15 NSCLC | Author UMI matrix (92,330 cells) + sample xlsx | Marker hierarchy. CopyKAT IDs are not on GEO. |
| [GSE241934](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241934) | Zhang et al. *Cell Rep Med* 2024, NEOTIDE/CTONG2104 + real-world | Author `major.cell.type` metadata + TACSTD2 streamed from MTX | Author labels (Epi/T/NK/B/Myeloid/Fibro/Endo/Mast) |
| [GSE291670](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291670) | Xia et al. *J Transl Med* 2025, anlotinib + camrelizumab, 6 tumors | 10x MTX in `GSE291670_RAW.tar` | Same marker hierarchy as GSE207422 |

No FASTQ. No EGA/dbGaP. IIT (11 EGFR-mutant, sintilimab) and REAL (34 EGFR-WT, mixed PD-1) in GSE241934 are kept as separate cohorts and also pooled.

## Lineages

Six parts: Epithelial, T/NK, B, Myeloid, Stromal (fibroblast + endothelial + mast), Other.

GSE241934 maps author `major.cell.type` onto those parts (T+NK → T/NK).

GSE207422 / GSE291670: exclusive marker presence, PTPRC+ for immune, PTPRC− for epithelial/stromal. Priority T → NK → B → Myeloid → Mast → Epithelial → Endo → Fibro. Marker lists are in `scripts/build_compositions.py`.

Malignant TACSTD2 = mean `log1p(UMI)` in that patient’s epithelial cells.

## Tests (patient = unit)

scCODA: **not run**. `pip install sccoda` failed (rpy2 needs system R).

Fallback, all on the 6-part count table with CLR/ILR pseudocount 0.5:

1. **H1** (TACSTD2-high → lower T/NK). Samples with ≥200 cells and ≥20 epithelial cells.
   - Spearman of malignant TACSTD2 vs fraction / CLR / ILR of T/NK.
   - Median-split high vs low: Wilcoxon on those T/NK coordinates (one-sided high < low, plus two-sided).
   - Exact enumeration when n ≤ 16; otherwise asymptotic MWU.
2. **H2** (NMPR → higher epithelial, lower T/NK). Samples with ≥200 cells and MPR/NMPR (pCR counted as MPR).
   - Wilcoxon NMPR vs MPR on epithelial and T/NK fraction/CLR/ILR.
3. **Dirichlet-multinomial** two-group LRT (shared α vs group-specific α) with permutation p (200–500 shuffles).
4. **Combinations:** each series, GSE241934 IIT/REAL/both, and every union. Multi-cohort models add `C(cohort)` or cohort-median-center the Spearman.

Primary GSE207422 slice is **post-treatment surgery** (n=12). Pre-treatment biopsies are only in the “all labeled” sensitivity.

## What “matches the user trend” means

Both **directions** hold: CLR ρ(TACSTD2, T/NK) < 0 **and** NMPR median epithelial fraction > MPR **and** NMPR median T/NK < MPR. This is not a p<0.05 rule. The write-up states the p-values next to the sign.

## Residual-tumor confound

MPR/pCR residuals are defined by low remaining viable tumor. Higher epithelial fraction in NMPR on post-drug resections is expected from that definition, especially in GSE241934.

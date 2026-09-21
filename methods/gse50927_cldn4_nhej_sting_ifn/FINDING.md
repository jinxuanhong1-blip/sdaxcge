# GSE50927 mouse VILI — Cldn4 NHEJ, STING, and IFN panels

**Additive public mouse. Non-cancer.** Kage / Borok *Am J Physiol Lung Cell Mol Physiol* 2014, PMID [25106430](https://pubmed.ncbi.nlm.nih.gov/25106430/); GEO [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927). Whole lung after ventilator-induced lung injury, Illumina HiSeq 2000 (GPL13112), author EdgeR on mm9. This folder scores three pre-specified panels on the four deposited contrasts. It does not replace the earlier IFN/MHC/TJ result in `methods/gse50927_mouse_cldn4`.

Cldn4-high = WT (Cldn4-intact). Cldn4-low = Cldn4 KO. **KOlow / KOhigh are BAL-protein injury strata, not Cldn4-expression strata.**

## Honest n

| Item | n | Note |
|---|---:|---|
| Author EdgeR tables | 4 | naive KO vs WT; WT VILI vs naive WT; KO VILIhigh vs WT VILI; KO VILIlow vs WT VILI |
| GSM per arm | **1 vs 1** | design text says duplicates; 2 SRA runs per GSM; no count matrix |
| Series-matrix expression rows | 0 | `Sample_data_row_count=0` |
| Tumour / LUAD / ICI | 0 | mixed-background whole lung |

Do not write n=5 or n=10. Prerank GSEA is gene-set permutation on the deposited logFC (1000 perms, seed 42, minimum set size 6, BH FDR within the 12 sets of that contrast). The Wilcoxon test asks whether the panel’s logFC values sit away from zero. Neither is a mouse-level test. Author FDR on an unreplicated design is reported gene-by-gene and is not used as a sample-size claim.

Author logFC is KO minus WT, or VILI minus naive. Compact means for genotype contrasts are flipped to **Cldn4-high minus Cldn4-low**. Positive NES means the set is enriched at the **second** group (KO, or VILI on the induction contrast).

## Panels

| panel | genes used | what it is |
|---|---:|---|
| NHEJ core | 9/9 | Ku70/80, DNA-PKcs, Artemis, Lig4, XRCC4, XLF, Pol λ, Pol μ |
| KEGG NHEJ (mmu03450) | 13/13 | core plus Mre11a, Rad50, Fen1, Dntt |
| STING core | 6/6 | Mb21d1 (Cgas), Tmem173 (Sting), Tbk1, Ikbke, Irf3, Ddx41 |
| Reactome STING minus DNA-PK | 11/13 | official R-HSA-1834941 without Prkdc / Xrcc5 / Xrcc6 |
| IFN compact | 24/26 on the naive table | same list as the earlier folder; logCPM < 0 dropped |
| Hallmark IFN-γ / IFN-α | 181/188, 89/94 | MSigDB 2023.2 mouse |

KEGG cytosolic DNA-sensing (60/61) is reported as a mixed companion. It contains Cxcl10, Ccl5, Il6 and Il1b together with the sensors, so it is not a STING-machinery score. Sources and the alias map are in `genesets/SOURCES.txt`.

## One-row table

Mean is Cldn4-high minus low, except the WT VILI row, which is VILI minus naive. NES is on the deposited rank (positive = up in KO, or up in VILI).

| contrast | n | Cldn4 logFC | NHEJ core mean (p) | STING core mean (p) | IFN compact mean (p) | KEGG NHEJ NES (FDR) | STING core NES (FDR) | IFN-γ NES (FDR) | IFN-α NES (FDR) |
|---|---|---:|---|---|---|---|---|---|---|
| naive KO vs WT (**primary**) | 1 vs 1 | **−6.061** | +0.029 (0.82) | −0.061 (0.84) | **−0.797 (5.1×10⁻⁶)** | −0.75 (0.50) | +0.52 (0.50) | **+1.488 (0.006)** | **+1.601 (0.008)** |
| WT VILI vs naive WT | 1 vs 1 | **+3.959** | −0.015 (0.91) | −0.168 (0.44) | −0.181 (0.038) | +0.64 (0.54) | −0.87 (0.53) | −0.897 (0.53) | **−1.582 (0.012)** |
| KO VILIhigh vs WT VILI | 1 vs 1 | **−10.16** | +0.062 (1.00) | −0.211 (0.56) | **−0.322 (0.018)** | −0.68 (0.46) | +0.95 (0.46) | **+1.670 (0.012)** | +1.229 (0.20) |
| KO VILIlow vs WT VILI | 1 vs 1 | **−12.69** | +0.058 (0.73) | −0.058 (1.00) | +0.046 (0.43) | −0.97 (0.50) | +0.75 (0.50) | +1.092 (0.50) | +0.613 (0.50) |

## Primary: naive Cldn4-intact vs Cldn4-null

`GSE50927_Cldn4lungWTvsKOgenes.csv.gz`. Cldn4 logFC = **−6.061** (logCPM 2.25; submitter FDR 4.1×10⁻²⁶).

**IFN.** Compact mean (high − low) = **−0.797** on 24 genes (Ifng and Cxcl11 have logCPM < 0 and are dropped). 21 of 24 genes are lower in WT. Wilcoxon p = 5.1×10⁻⁶. This is the same compact number as the earlier folder. Hallmark enrichment scores match that folder exactly (IFN-γ ES 0.506, IFN-α ES 0.591; 181 and 89 genes). NES in this run is +1.488 (FDR 0.006) and +1.601 (FDR 0.008). The earlier folder, which permuted Hallmark first from the same seed, reported NES +1.508 and +1.587. Same genes, same ES, same call: IFN is higher after Cldn4 loss. Lead genes include Cxcl9 (+2.88), Ccl5 (+2.71), Cxcl10 (+2.54), Isg15 (+1.08, FDR 0.0031).

**NHEJ.** The nine core enzymes do not move. Mean high − low = +0.029 (5 higher in WT, 4 higher in KO; Wilcoxon p = 0.82; background MWU p = 0.76). KEGG NHEJ NES = −0.75 (FDR 0.50). Deposited logFC (KO − WT): Xrcc6 +0.31, Xrcc5 −0.01, Prkdc +0.05, Lig4 +0.27, Xrcc4 −0.17, Nhej1 +0.21, Dclre1c −0.36, Poll −0.09, Polm −0.47. Every core gene has submitter FDR ≥ 0.21. GO NHEJ (38 genes) and Reactome NHEJ (32 genes; 19 histone symbols have no row) are also flat (NES +0.46 and +0.49, FDR 0.50).

**STING.** The six machinery genes do not move. Mean high − low = −0.061 (2 higher in WT, 4 higher in KO; Wilcoxon p = 0.84). STING-core NES = +0.52 (FDR 0.50). Reactome STING with DNA-PK removed: NES +0.76 (FDR 0.50). Deposited logFC: Mb21d1 +0.50 (FDR 0.67), Tmem173 +0.13, Tbk1 +0.12, Ikbke −0.27, Irf3 +0.05, Ddx41 −0.15. Trex1 and Ifnb1 are undetected (logCPM −2.07, logFC 0).

KEGG cytosolic DNA-sensing NES = +1.35 (nominal p = 0.028, FDR 0.084). That lean tracks the chemokine / cytokine members of the set. It is not a cGAS–STING transcript result.

Reading on this contrast: Cldn4-high lung is IFN-low, which agrees with the locked IFN result. The NHEJ enzyme panel and the STING machinery panel stay put. The ISG rise is present with Ifnb1 and Trex1 below detection.

## Companion: WT VILI induces Cldn4

`GSE50927_VILIwtGenes.csv.gz`. Cldn4 logFC = **+3.959** (FDR 7.2×10⁻⁹¹). This is injury versus the naive WT lung, not a genotype split.

NHEJ core mean (VILI − naive) = −0.015 (p = 0.91). STING core mean = −0.168 (p = 0.44); Irf3 −0.42 (FDR 0.071), Tbk1 +0.32 (FDR 0.25). Hallmark IFN-α ES = −0.548, NES = −1.582 (FDR 0.012), the same ES as the earlier folder. IFN compact mean = −0.181 on 23 genes (p = 0.038). Tnf (+3.27), Il1b (+1.64) and Egr1 (+0.73) rise with the injury program and are not used as the STING or IFN score. Ifnb1 and Trex1 remain undetected.

When WT lung induces Cldn4, IFN-α goes down and the NHEJ and STING machinery panels do not.

## Companion: VILI, Cldn4-null vs Cldn4-intact

Both KO arms are Cldn4-null (logFC −10.16 and −12.69).

**KO VILIhigh** (higher BAL leak). IFN compact mean high − low = −0.322 (p = 0.018), matching the earlier folder. Hallmark IFN-γ ES = 0.559, NES = +1.670 (FDR 0.012). NHEJ core mean = +0.062 (5 vs 4; p = 1.00). Lig4 is the largest core shift (KO − WT −1.06, FDR 0.11) and the other eight core genes stay within about ±0.5. STING core mean = −0.211 (p = 0.56). Tmem173 is +0.86 in the KO (FDR 0.22); Mb21d1 is −0.18. Paper injury genes are up (Tnf +2.01, FDR 3.8×10⁻⁸; Il1b +1.73, FDR 2.6×10⁻⁷; Egr1 +2.13, FDR 3.3×10⁻¹⁷).

**KO VILIlow** (WT-like leak). NHEJ, STING and IFN are all flat (IFN compact +0.046, p = 0.43; Hallmark FDR 0.50). Tnf, Il1b and Egr1 are not FDR-significant versus WT VILI.

The high-injury Cldn4-null lung is IFN-γ–high. The WT-like-injury Cldn4-null lung is not. NHEJ and STING machinery stay flat on both arms.

## What this adds

On the only public Cldn4-null mouse lung RNA-seq series, the IFN panel repeats the earlier result and the NHEJ and STING machinery panels do not move with genotype or with Cldn4 induction. That is additive public evidence from a non-cancer, whole-lung injury experiment. It cannot be written as a tumour NHEJ or STING law, and it cannot be written as n>1.

## What this does not test

- A tumour, LUAD, KL/KP, or ICI comparison.
- Sample-level Spearman or Cldn4 quartiles (no per-sample matrix).
- NHEJ repair activity, DNA-PK kinase activity, micronuclei, or cGAS/STING protein.
- Ifnb1 induction (the gene is below detection on all four tables).
- KOlow / KOhigh as Cldn4-expression classes.
- A failed audit of the IFN result. The naive IFN compact mean is −0.797, as before.

## Reproduce

```bash
python3 -m pip install -r methods/gse50927_cldn4_nhej_sting_ifn/requirements.txt
python3 methods/gse50927_cldn4_nhej_sting_ifn/analyze.py
```

Downloads (not committed) go to `$GSE50927_NHEJ_DATA` (default `/tmp/gse50927_nhej`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtGenes.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz`

## Files

- `analyze.py` — download, map panels, prerank GSEA, compact means, figures
- `gsea_core.py` — same weighted-KS prerank engine as the earlier GSE50927 folder
- `genesets/panels.json` — frozen NHEJ, STING, and IFN lists
- `genesets/SOURCES.txt`
- `tables/one_row.tsv`
- `tables/compact_panel_scores.tsv`
- `tables/gsea_panels.tsv`
- `tables/gene_level.tsv`
- `tables/gene_coverage.tsv`
- `tables/cldn4_logfc.tsv`
- `tables/summary.json`
- `figures/fig1_compact_panels.png`
- `figures/fig2_gsea_nes.png`
- `figures/fig3_focal_logfc.png`

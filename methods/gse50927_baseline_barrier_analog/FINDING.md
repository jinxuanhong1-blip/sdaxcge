# GSE50927 baseline Cldn4 KO lung: IFN / chemokine / MHC

**Barrier-loss analog only. Not cancer.** This is uninjured mouse whole lung after germline Cldn4 deletion (Kage / Flodby / Borok, *Am J Physiol Lung Cell Mol Physiol* 2014, PMID [25106430](https://pubmed.ncbi.nlm.nih.gov/25106430/); GEO [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927)). It is not a tumour, not LUAD, not ICI, and not evidence about CLDN4-high cancer cells.

The KD-match contrast is **naive Cldn4 KO vs WT (no VILI)**. Ventilator-induced lung injury (VILI) tables are used only to pull injury genes out of that contrast.

## Honest n

| Item | Value |
|---|---|
| Primary contrast | naive KO vs WT, author edgeR `GSE50927_Cldn4lungWTvsKOgenes.csv` |
| Sign | logFC = KO − WT |
| Honest n | **1 vs 1 GEO sample** |
| Design text | “all in duplicates”; those replicates are not two processed columns |
| Expression matrix | none (`Sample_data_row_count = 0`) |
| Cldn4 | logFC **−6.061**, author FDR **4.07×10⁻²⁶** (KO check) |
| Strain | mixed 129S6 / C57BL/6 / BALB/c |

Author FDR is an edgeR dispersion call on that unreplicated design. It confirms the deposited ranking. It is not a multi-mouse test. Gene-wise Wilcoxon and binomial p-values below treat genes as the units. Those genes are correlated, so the p-values are descriptive.

Hit rule, fixed before the class counts: **FDR < 0.05, logCPM > 0, logFC > 0**. logCPM > 0 drops the edgeR floor (about −2.07). A gene is **injury-shared** when it is also FDR-up in WT VILI vs naive WT. A gene is **injury-only** when that VILI contrast is the hit and the baseline contrast is not.

## Result

| Module | What is counted in the analog | What is separated as injury |
|---|---|---|
| IFN / ISG | **17** FDR-up genes with logCPM > 0, and all 17 are baseline-only (each logCPM ≥ 1). **0** are injury-shared. 8 of the 17 are FDR-**down** in WT VILI | Stat3 is the one panel gene that is a VILI hit and not a baseline hit |
| Chemokine | High-CPM baseline-only: **Ccl5, Cx3cl1**. Low-CPM baseline-only: Ccl8, Ccl20 | Injury-shared: **Cxcl1, Ccl4**, plus low-baseline-CPM **Cxcl2**. Injury-only: **Ccl2, Ccl3, Ccl7, Ccl9, Ccl11, Ccl17** |
| MHC-I classical | **B2m** and **Psmb9** only (2 / 10). H2-K1 is unchanged (logFC −0.071, FDR 1) | No classical MHC-I gene is injury-shared. B2m falls in WT VILI |
| MHC-II | Not up. H2-Oa and H2-Ob are FDR-down | — |
| Nonclassical H2 | **Not counted.** H2-K2 / H2-M2 / H2-Bl are large and track the KO line in both VILI arms | Haplotype / mixed-background risk, not a barrier-loss MHC call |
| Paper injury mediators | Not the baseline program | Tnf, Il1b, Il6 are injury-only. Egr1 is **down** at baseline and up with injury |

## IFN / ISG (chemokines held out)

Panel is the public IFN list already used for this accession, minus Cxcl9 / Cxcl10 / Cxcl11 (those sit in the chemokine module). 60 / 62 symbols are in the table (Ifi27 and Tgtp1 are absent). 52 genes have baseline logCPM > 0; 46 / 52 have logFC > 0 (median **+0.45**). Gene-wise Wilcoxon p = **1.9×10⁻⁸**. Dropping injury-class genes leaves p = **2.6×10⁻⁸**.

The 17 FDR hits are all baseline-only:

Isg15, Oas1g, Oas2, Oas3, Oasl2, Ifit1, Ifitm1, Ifitm3, Gbp3, Gbp4, Gbp5, Gbp6, Ifi44, Ifi47, Irgm2, Igtp, Tgtp2.

Eight of those seventeen are FDR-down when a WT lung is ventilated (Oasl2, Gbp3, Gbp4, Gbp5, Ifi47, Irgm2, Igtp, Tgtp2). The resting ISG cassette is not the ventilator program.

Not counted: Ubd (logFC +5.62, logCPM −0.75) and Irg1 (logFC +5.49, logCPM −0.82) are FDR < 0.05 only at the detection floor. Stat3 is injury-only (baseline FDR 0.67; VILI logFC +0.43, FDR 0.018).

## Chemokines

The ligand catalog has 32 symbols, all present; 22 have baseline logCPM > 0. A panel-wide shift does **not** survive the injury split: gene-wise Wilcoxon p = 0.010 on all 22 detected ligands, and p = **0.083** after injury-shared and injury-only genes are removed. The sign test on the 22 is p = 0.29. The analog is the individual baseline-only ligands, not a chemokine-set call.

| Gene | Class | Baseline logFC (FDR) | logCPM | WT VILI logFC (FDR) |
|---|---|---:|---:|---:|
| Ccl5 | baseline only | +2.710 (0.0011) | 4.96 | +0.32 (0.72) |
| Cx3cl1 | baseline only | +0.970 (8.5×10⁻⁷) | 6.70 | −0.46 (0.19) |
| Ccl8 | baseline only, low CPM | +1.962 (0.044) | 0.35 | +1.54 (0.30) |
| Ccl20 | baseline only, low CPM | +2.227 (0.027) | 0.35 | +2.50 (0.21) |
| Cxcl1 | injury-shared | +1.632 (0.021) | 1.72 | +3.21 (1.4×10⁻¹⁴) |
| Ccl4 | injury-shared, low CPM | +2.512 (0.047) | 0.10 | +2.92 (0.0017) |
| Cxcl2 | injury-shared, baseline logCPM ≤ 0 | +2.639 (0.037) | −0.12 | +5.10 (2.8×10⁻¹⁸) |
| Ccl2 | injury only | +0.31 (1) | −0.37 | +3.82 (6.3×10⁻¹⁶) |

Ccl3, Ccl7, Ccl9, Ccl11, and Ccl17 are the other injury-only ligands. Cxcl9 (logFC +2.88, FDR 0.062) and Cxcl10 (logFC +2.54, FDR 0.094) move up at baseline and are not FDR hits, so they are not counted.

## MHC

Classical MHC-I / APM (10 genes, all detected): **B2m** logFC +0.655, FDR 0.034, logCPM 9.73; **Psmb9** logFC +0.902, FDR 0.037, logCPM 5.01. **H2-K1** logFC −0.071, FDR 1. The other seven (H2-D1, Tap1, Tap2, Tapbp, Psmb8, Psmb10, Nlrc5) are not FDR < 0.05. Eight of ten logFC values are positive (median +0.36; Wilcoxon p = 0.0098; sign test p = 0.11). That is a small consistent tilt, not ten discoveries. B2m is FDR-down in WT VILI (logFC −0.59, FDR 3.1×10⁻⁴), so the B2m increase is not the injury program.

MHC-II is not opened. H2-Oa logFC −1.40 (FDR 0.001) and H2-Ob logFC −2.29 (FDR 1.1×10⁻²²).

Nonclassical H2 genes are excluded from the numerator. On this mixed 129 / B6 / BALB line the extreme calls are H2-K2 logFC +4.26 (FDR 1.4×10⁻⁶⁵), H2-M2 +5.18 (FDR 3.4×10⁻⁵⁴), and H2-Bl +3.90 (FDR 9.6×10⁻²⁵). The same genes stay up in both VILI KO arms (KOhigh and KOlow), which is the pattern of a line difference rather than ventilator injury. H2-Q4 is absent from the table.

## Injury confounds (not the analog)

WT VILI vs naive WT induces the paper’s mediators and leaves Cldn4 itself up (logFC **+3.959**, FDR 7.17×10⁻⁹¹). That contrast is induction by injury, not Cldn4 loss.

| Gene | Baseline KO vs WT | WT VILI vs naive WT | Class |
|---|---|---|---|
| Tnf | +0.24, FDR 1 | +3.27, FDR 3.3×10⁻¹⁰ | injury only |
| Il1b | +0.006, FDR 1 | +1.64, FDR 1.6×10⁻⁷ | injury only |
| Il6 | −0.96, FDR 1 | +4.65, FDR 3.8×10⁻²¹ | injury only |
| Egr1 | **−2.14, FDR 3.3×10⁻⁵** | +0.73, FDR 0.036 | down at rest, up with injury |

KOhigh vs WT VILI is the paper’s high-BAL-protein stratum. It is not a Cldn4-expression class. Those rows are in `tables/gene_level.tsv` and are not folded into the 17-gene IFN count.

## What this does not say

- Nothing about cancer, LUAD, ICI, or CLDN4-high tumour cells.
- No mouse-level p. n is 1 vs 1.
- No claim that the whole chemokine catalog opens after Cldn4 loss. After the injury split, it does not.
- No claim that classical MHC-I heavy chain (H2-K1) rises.
- No use of H2-K2 / H2-M2 / H2-Bl as antigen-presentation evidence.
- Ubd and Irg1 are not evidence. Their counts sit below logCPM 0.

## Reproduce

```bash
python3 -m pip install -r methods/gse50927_baseline_barrier_analog/requirements.txt
python3 methods/gse50927_baseline_barrier_analog/analyze.py
```

Downloads go to `$GSE50927_BASELINE_DATA` (default `/tmp/gse50927_baseline_barrier`):

- `GSE50927_Cldn4lungWTvsKOgenes.csv.gz` — baseline KO vs WT
- `GSE50927_VILIwtGenes.csv.gz` — WT VILI vs naive WT
- `GSE50927_VILIwtkohiGenes.csv.gz` — KO VILIhigh vs WT VILI
- `GSE50927_VILIwtkoloGenes.csv.gz` — KO VILIlow vs WT VILI

The script checks that baseline Cldn4 logFC is −6.061296179.

## Files

- `analyze.py` — download, classify, gene-wise tests, figure
- `tables/gene_level.tsv` — every panel gene across the four contrasts
- `tables/module_summary.tsv`
- `tables/cldn4_qc.tsv`
- `tables/label_inventory.tsv`
- `tables/summary.json`
- `figures/fig_baseline_vs_injury.png`

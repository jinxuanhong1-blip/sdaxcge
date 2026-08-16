# Methods: B5 leftover 2023–2026 open lung ICI (additive only)

## Scope

User claim **B5** (11-cohort CLDN4-high ICI **OR = 0.42**) is **taken as given**. This slice does **not** re-cut those cohorts and does **not** pool new rows into that OR.

The question here is only: among **new** open lung ICI **bulk** series from 2023–2026 that were **not** in the original 11, what are the honest extra rows (n, OR or rank-biserial, p) for CLDN4 and TACSTD2 versus DCB / ORR / PFS?

Public GEO / EuropePMC only. No EGA, no dbGaP, no FASTQ.

## Locked analysis (do not pick the cutoff nearest 0.42)

1. Score **CLDN4** and **TACSTD2** on the deposited processed matrix (STAR gene counts → `log2(CPM+1)` using the **full-transcriptome** library size; FPKM or log2(RPKM+1) when that is what GEO deposited).
2. Primary cutoff for a binary OR: **within-cohort median** (ties at the median assigned to low).
3. Sensitivity: **tertile T3 vs T1** (middle tertile dropped).
4. Continuous complement: two-sided Mann–Whitney U and **rank-biserial / Cliff’s δ** (same sign: positive if responders are higher).
5. OR uses the 2×2 of gene-high vs gene-low × responder vs non-responder. Zero cells get Haldane–Anscombe +0.5; p is Fisher’s exact. Woolf 95% CI on the log-OR.
6. DCB, when PFS is available: PFS ≥ 6 months. PFS time-to-event: median-split log-rank and univariable Cox per 1 SD. The Table S1 column `Response status available` is used as the event indicator (it tracks progression, not “label present”).
7. **NSCLC and SCLC are never mixed** into one OR.

## Must-try accessions

| Accession | What GEO actually is | What can be tested |
|---|---|---|
| **GSE328294** | H23 KRAS-mutant NSCLC cell line, DMSO vs IBI351, n=6 FPKM (public 2026) | Treatment contrast only. Abstract mentions KRAS G12C clinical ICI specimens; **those labels are not in the series**. |
| **GSE244944** | SCLC ChIP-seq | No expression-vs-ICB test. |
| **GSE244945** | Human/mouse SCLC cell-line bulk RNA-seq (NOTCH1 / LSD1 / REST) | Cell-line contrast only. IMpower133 RNA is **not** in GEO. |
| **GSE244946** | COR-L88 scRNA-seq, n=3 | Not patient ICB; not bulk tumor. |

## Other 2023–2026 series that were actually downloadable with response labels

| Accession | Why it is in / out |
|---|---|
| **GSE274975** (2024) | **In.** 61 FFPE lung tumors, STAR counts. RECIST + PFS from EuropePMC Table S1 (PMID 39723204). One GEO-annotated SCLC (`OB_pat_LuC_92`) is held out of the NSCLC OR. Four patients lack RECIST and are ORR-ineligible (kept for DCB/PFS). |
| **GSE283829** (2025) | **In.** 27 NSCLC tumors, raw counts. GEO `disease stage` is CR/SD/PD (no PR). |
| GSE202417 | Out of the tumor OR: Clariom D of **PBMC CD8 T cells**, not tumor bulk. |
| GSE317309, GSE291670 | Out: scRNA / PBMC, not bulk tumor. |
| GSE207422 | Out: already in the original 11; not recut. |
| GSE253564 / GSE248378 | Not recut here (separate leftover already exists). |
| IMpower133 trial RNA | Out: not on GEO; sponsor-controlled. |

## Reproduction

```bash
python3 -m pip install -r methods/B5_leftover_2024_26/requirements.txt
python3 methods/B5_leftover_2024_26/analyze.py
```

Outputs land in `results/B5_leftover_2024_26/`.

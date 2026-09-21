# GSE325414 malignant CLDN4 vs T/NK

Additive public series. Not merged into the concordant-4 (n=65) result.

**Series:** [GSE325414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325414), Jimenez et al., pulsed-electric-field treat-and-resect early NSCLC. BD Rhapsody whole-transcriptome RNA, 25 donors, 156,467 cells. GEO defines INXBX as the index diagnostic biopsy, RESRT as the resected treatment area, and RESRU as resected untreated tumor. Every library in this deposit is the PEF arm. This is clinical stage and histology, not an ICI cohort.

**Why this one:** among 224 GEO series dated 2024–2026 it is the largest open matrix that has CLDN4, author-labeled tumor epithelium, and author-labeled T/NK in the same donors, plus stage. The best unused ICI-response series is GSE233203 (pleural fluid, n=7), scored separately in `results/gse233203/FINDING.md`. GSE243013 (n=234, ICI) and GSE280232 (n=13, neoadjuvant nivolumab + ipilimumab) are immune-only or sorted T cells. Candidate table: `CANDIDATES.md`.

## LUAD stability sweep

Restricted to RESRU donors whose histology contains adenocarcinoma (mucinous included; strict LUAD is separate). Squamous donors stay out of this sweep. A result is a **stable inverse** only if n≥8, ρ<0, p<0.05, and every single-donor deletion still has ρ<0 and p<0.05. The grid is sensitivity, not a discovery p-value. Full grid: `luad_sweep_grid.tsv` (630 CLDN4 vs T/NK tests).

**CLDN4 mean is not stable.** Zero of those tests meet the rule. At the original gate (≥50 malignant, ≥30 T/NK, author tumor epithelium, n=12) the mean is ρ=−0.62, p=0.031, and dropping T13 moves it to ρ=−0.51, p=0.11. Raising the malignant floor does not rescue it: once T13 (88 malignant cells) falls out, ρ is weaker. The same fragility holds after dropping the top quartile of normal-lung score inside author tumor cells (without T13, p=0.20).

The only nominal mean result that does not contain T13 is author tumor cells with normal-lung score <0.25 (n=8, ρ=−0.71, p=0.047). That call already removed T13, and dropping T24 instead moves p to 0.18. It is one-donor fragile as well.

**The percent-positive call is the strongest stable inverse, and it is thin.** At the same pre-specified author gate, the fraction of malignant cells with CLDN4 UMI>0 is ρ=−0.73, p=0.0074, n=12. Every leave-one-out stays p<0.05. The weakest deletion is T24 (ρ=−0.65, p=0.032), not T13 (ρ=−0.71, p=0.026). Strict LUAD, dropping two mucinous donors, is the same story at n=10 (ρ=−0.77, p=0.0092; without T13, p=0.030; worst leave-one-out T24, p=0.042).

Searching the grid does not improve this. Of 168 nominal inverses, 63 pass the every-deletion rule, and all 63 are percent metrics, never the mean. The most negative stable cell (ρ=−0.76) only appears when the malignant floor is lowered to 1 cell. That is not a better result than the pre-specified percent test.

Marker epithelium (epithelial score ≥1, PTPRC log <0.3, normal-lung score <0.3) leaves n=6 at the original gate. Its percent correlation is ρ=−0.93, p=0.0077, and the worst leave-one-out is p=0.054. Below the stability rule. Not used.

Figure: `fig_luad_t13_sensitivity.png`.

## Definitions

- Malignant: author `sub.pop.level2 = EpithelialCells_TumorCells`
- T/NK: author `Tcell` or `NKcell` (CD8 is `sub.pop.level3 = Tcell_CD8`)
- CLDN4: mean log1p(CP10k) inside malignant cells; percent is UMI > 0
- Unit: donor. Gate: ≥50 malignant cells and ≥30 T/NK cells
- Primary compartment: RESRU. RESRT is expected to recruit lymphocytes after ablation, so it is not the primary test

Tumor epithelium expresses the epithelial program (CLDN4 mean log 1.69 vs 0.098 in other cells; 88% CLDN4-positive; EPCAM mean log 1.83 vs 0.001).

## Primary result (RESRU)

| Test | n | ρ | p |
|---|---:|---:|---:|
| CLDN4 mean vs T/NK, all histology | 17 | −0.26 | 0.31 |
| CLDN4 percent vs T/NK, all histology | 17 | −0.42 | 0.090 |
| CLDN4 mean vs CD8 fraction, all histology | 17 | −0.29 | 0.25 |
| TACSTD2 mean vs T/NK, all histology | 17 | −0.37 | 0.14 |
| CLDN4 mean vs T/NK, adenocarcinoma (includes mucinous) | 12 | −0.62 | 0.031 |
| CLDN4 percent vs T/NK, adenocarcinoma | 12 | −0.73 | 0.0074 |
| CLDN4 mean vs CD8 fraction, adenocarcinoma | 12 | −0.47 | 0.12 |
| CLDN4 mean vs T/NK, squamous | 5 | +0.60 | 0.28 |

Three RESRU donors fail the cell gate (T29, T30, T31) and are excluded.

The pre-specified all-histology test is negative and not significant. Squamous tumors (n=5) point the other way and dilute it. The adenocarcinoma slice is nominally negative. Leave-one-out: the mean correlation loses p<0.05 if T13 is removed (ρ=−0.51, p=0.11). The percent correlation stays ρ≤−0.65 and p≤0.032 on every single-donor deletion. n=12 is still a subset, not a second concordant cohort.

## Other compartments (not primary)

| Compartment | CLDN4 mean vs T/NK | n | ρ | p |
|---|---|---:|---:|---:|
| RESRT ablated area, all histology | mean | 13 | +0.10 | 0.75 |
| INXBX index biopsy, all histology | mean | 13 | +0.46 | 0.12 |

The negative RESRU adenocarcinoma association is not repeated in the ablated area or the index biopsy. Slices with n=4 are not interpreted.

## Read this as

The LUAD mean correlation is real only while T13 is kept. Document it as fragile to T13. The percent-positive correlation at the same gate survives every one-donor deletion, including T13, with the weakest deletion at p=0.032 (T24). That is the strongest stable inverse in the sweep. It is still n=12 and is not merged into concordant-4. No gate or malignant definition makes the mean stable. This series is not ICI evidence.

Figure: `fig_resru_cldn4_vs_tnk.png`. Tables: `per_donor_compartment.tsv`, `stats.tsv`, `leave_one_out_resru.tsv`.

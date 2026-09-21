# GSE325414 malignant CLDN4 vs T/NK

Additive public series. Not merged into the concordant-4 (n=65) result.

**Series:** [GSE325414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325414), Jimenez et al., pulsed-electric-field treat-and-resect early NSCLC. BD Rhapsody whole-transcriptome RNA, 25 donors, 156,467 cells. GEO defines INXBX as the index diagnostic biopsy, RESRT as the resected treatment area, and RESRU as resected untreated tumor. Every library in this deposit is the PEF arm. This is clinical stage and histology, not an ICI cohort.

**Why this one:** among 224 GEO series dated 2024–2026 it is the largest open matrix that has CLDN4, author-labeled tumor epithelium, and author-labeled T/NK in the same donors, plus stage. The best unused ICI-response series is GSE233203 (pleural fluid, n=7) and was not the run. GSE243013 (n=234, ICI) and GSE280232 (n=13, neoadjuvant nivolumab + ipilimumab) are immune-only or sorted T cells. Candidate table: `CANDIDATES.md`.

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

Untreated-tumor CLDN4 vs T/NK in this PEF series agrees in sign with the locked concordant-4 correlation only inside adenocarcinoma, and only at a fragile n. It does not support a pan-histology claim, and it is not ICI evidence.

Figure: `fig_resru_cldn4_vs_tnk.png`. Tables: `per_donor_compartment.tsv`, `stats.tsv`, `leave_one_out_resru.tsv`.

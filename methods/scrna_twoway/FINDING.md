# Two-way extra: public lung scRNA with timepoint × response

Public hunt for lung scRNA series that carry **both** a timepoint (pre/post) and a response label (MPR or RECIST), then malignant TACSTD2 / CLDN4 in that 2×2. Patient is the unit of n. No second series was stacked.

## Hunt

eUtils queries plus a curated seed of known lung ICI scRNA accessions (135 series scored). A series is two-way only if the **public** record has a timepoint axis, an MPR/RECIST (or equivalent) axis, and a processed matrix that includes epithelial / malignant cells.

| Accession | Timepoint | Response | Malignant matrix | Two-way? | Why not stacked |
|---|---|---|---|---|---|
| **GSE207422** | 3 pre + 12 post (unpaired) | MPR/NMPR/pCR/NE and RECIST on the GEO xlsx | Yes (UMI, 92,330 cells) | **Yes** | — |
| GSE291670 | Post only (n=6) | 3 MPR + 3 Non-MPR | Yes (MTX) | No | Paper TN samples are HRA001033 (GSA), not in this GEO series |
| GSE205335 | Response core is pre-ICI | RECIST on 14 core / 11 patients | Yes (RDS) | No | Post-ICI samples were excluded from the response core; cell-identity table has no timepoint/response columns |
| GSE243013 | Post only (n=234) | In the paper, not in GEO characteristics | No public UMI+label table | No | Raw held at NGDC; GEO sample traits are tissue / disease only |
| GSE146100 | Post only (3 nodules, 1 patient) | Radiographic R vs NR per nodule | Yes | No | No pre-treatment scRNA |
| GSE179994 | Paired pre/post | Good / poor | No (T cells only) | No | Cannot score malignant TACSTD2/CLDN4 |
| GSE229353 | scRNA is post-surgery CD45+ | Paper has MPR-style outcomes | Immune only | No | No malignant compartment |
| GSE253013 | Treatment-naive | None | Yes | No | No ICI labels |

**Stacked series: none.** Only GSE207422 is public two-way with malignant cells.

## GSE207422 2×2 occupancy

GEO `GSE207422_NSCLC_scRNAseq_metadata.xlsx`. pCR (P06) counted as MPR. The three pre biopsies are **different patients** from the twelve post resections (not paired deltas). P01 is pre / pathologic NE and is omitted from the MPR grid (kept on the RECIST grid as SD).

|  | MPR | NMPR |
|---|---|---|
| **pre** | **empty** (n=0) | P05, P08 (n=2) |
| **post** | P03, P06, P11, P14 (n=4) | P02, P04, P07, P09, P10, P12, P13, P15 (n=8) |

RECIST (PR vs SD; no CR/PD in the xlsx): pre-PR n=1 (P05), pre-SD n=2 (P01, P08), post-PR n=8, post-SD n=4.

**Missing:** pre-MPR; within-patient pre→post; a second public series to stack.

## Malignant TACSTD2 / CLDN4

Author CopyKAT barcodes are not on GEO. Malignant = epithelial lineage minus clear alveolar / club / ciliated programs, with a chromosome-dispersion support (stromal reference). All 15 patients have ≥10 malignant cells under this call (9,849 malignant cells). Metric = patient mean log1p(CP10k). Mann–Whitney two-sided.

### MPR grid (mean log1p CP10k)

| Gene | pre MPR | pre NMPR | post MPR | post NMPR | post NMPR vs MPR | NMPR post vs pre | MPR post vs pre |
|---|---|---|---|---|---|---|---|
| TACSTD2 | empty | 1.75 (n=2) | 1.33 (n=4) | 1.69 (n=8) | Δ=+0.36, p=0.21 | p=0.89 | empty pre |
| CLDN4 | empty | 1.23 (n=2) | 1.55 (n=4) | 1.52 (n=8) | Δ=−0.03, p=0.93 | p=0.40 | empty pre |

Percent-positive TACSTD2: post NMPR 84.3% vs post MPR 85.2%, p=0.81. CLDN4 percent-positive: 82.9% vs 90.2%, p=0.57.

### RECIST grid (mean log1p CP10k)

| Gene | pre PR | pre SD | post PR | post SD | post SD vs PR | pre SD vs PR |
|---|---|---|---|---|---|---|
| TACSTD2 | 1.86 (n=1) | 1.70 (n=2) | 1.56 (n=8) | 1.58 (n=4) | p=1.00 | n=2 vs 1 |
| CLDN4 | 0.89 (n=1) | 1.81 (n=2) | 1.55 (n=8) | 1.48 (n=4) | p=0.93 | n=2 vs 1 |

No testable 2×2 contrast has n≥2 in every cell. The only fully populated within-timepoint contrast is post MPR vs NMPR (4 vs 8). None of the patient-level tests above are p<0.05.

## Files

- Hunt: `results/scrna_twoway/hunt_series.tsv`
- Patient table: `results/scrna_twoway/per_patient.tsv`
- Grid / tests: `results/scrna_twoway/twoway_grid.tsv`, `twoway_tests.tsv`
- Figure: `results/scrna_twoway/fig_twoway_mpr.png`
- Scripts: `methods/scrna_twoway/`

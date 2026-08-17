# Extra — GSE207422 T/NK exhaustion / cytotoxicity vs malignant CLDN4

Additive extra on public [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) (Hu et al., *Genome Med* 2023, PMID 36869384). **Not a TACSTD2 redo.** The named question is T/NK **cytotoxicity** and **exhaustion** versus **malignant CLDN4**. Patient / one sample per patient is the unit. Cell-level p-values are not reported.

## Question

In the 12 post-resection NSCLC tumors (neoadjuvant PD-1 + chemo), does mean malignant CLDN4 track T/NK cytotoxicity (`GZMB`, `PRF1`, `GNLY`, `NKG7`) or exhaustion (`PDCD1`, `HAVCR2`, `LAG3`, `TIGIT`, `TOX`)?

## Honest n

- Matrix: **92,330** cells (paper post-QC count). QC here: UMI ≥ 200 → still **92,330**.
- Author CopyKAT / per-cell labels are **not on GEO**. Malignant is a marker proxy.
- **Primary malignant-like** = epithelial lineage **and** tumor-epi module > normal-lung module **and** normal-lung < 0.4 (same gate as the multi-cohort exhaustion meta).
- Primary filter: ≥10 malignant-like **and** ≥10 T/NK. Post-treatment: **7 / 12** patients enter Spearman.
- Dropped post (malignant-like < 10): **P02** (2, NMPR), **P06** (4, pCR/MPR), **P11** (1, MPR), **P13** (5, NMPR), **P14** (2, MPR).
- Under the primary gate, MPR with usable malignant CLDN4 is **n=1 (P03 only)**. Do not read an MPR contrast from n=1.
- Pre-treatment biopsies: 3 patients (P01 NE, P05 NMPR, P08 NMPR); unpaired with the 12 posts. Not stacked into the primary ρ.
- All 12 posts have ≥10 T/NK. All-epithelial sensitivity keeps every post (n=12).

### Per-patient occupancy

| Patient | Timing | MPR | Histology | n cells | n malig-module | n malig-zero | n epithelial | n T/NK | n CD8 | In primary? |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|
| P02 | post | NMPR | Adeno | 5345 | 2 | 14 | 118 | 2897 | 949 | no |
| P03 | post | MPR | Squamous | 9259 | 1210 | 978 | 1361 | 5060 | 2940 | yes |
| P04 | post | NMPR | Adeno | 8398 | 63 | 49 | 159 | 5625 | 2506 | yes |
| P06 | post | MPR (pCR) | Adeno | 4649 | 4 | 4 | 74 | 2945 | 1893 | no |
| P07 | post | NMPR | Squamous | 9131 | 4847 | 3225 | 5443 | 1078 | 399 | yes |
| P09 | post | NMPR | Squamous | 4849 | 290 | 185 | 299 | 3279 | 1890 | yes |
| P10 | post | NMPR | Squamous | 6643 | 327 | 320 | 389 | 4092 | 2664 | yes |
| P11 | post | MPR | Adeno | 3759 | 1 | 0 | 460 | 1687 | 673 | no |
| P12 | post | NMPR | Adeno | 6233 | 462 | 433 | 596 | 1480 | 674 | yes |
| P13 | post | NMPR | Squamous | 6753 | 5 | 14 | 48 | 1364 | 445 | no |
| P14 | post | MPR | Squamous | 5546 | 2 | 0 | 207 | 2850 | 912 | no |
| P15 | post | NMPR | Adeno | 8086 | 33 | 18 | 611 | 1532 | 513 | yes |
| P01 | pre | NE | Squamous | 4741 | 715 | 730 | 820 | 688 | 136 | no (pre) |
| P05 | pre | NMPR | Adeno | 3113 | 1505 | 1396 | 1513 | 808 | 204 | no (pre) |
| P08 | pre | NMPR | Adeno | 5825 | 9 | 21 | 199 | 3152 | 655 | no (pre) |

pCR (P06) is grouped as MPR. P01 is pathologic NE. That n is **not** 12 and **not** 92,330.

## Primary table (patient-level Spearman)

Score = mean log1p(CP10k) of the named genes in the named compartment.

| Contrast | Compartment | Subset | n | ρ | p | Patients |
|---|---|---|---:|---:|---:|---|
| CLDN4 vs T/NK cytotoxicity | malignant-like | post | **7** | **0.00** | **1** | P03, P04, P07, P09, P10, P12, P15 |
| CLDN4 vs T/NK exhaustion | malignant-like | post | **7** | **0.00** | **1** | P03, P04, P07, P09, P10, P12, P15 |

Primary n = 1 MPR (P03) + 6 NMPR.

## Sensitivity

| Contrast | Compartment | Subset | n | ρ | p |
|---|---|---|---:|---:|---:|
| CLDN4 vs T/NK cytotoxicity | zero-normal-UMI | post | 9 | −0.10 | 0.80 |
| CLDN4 vs T/NK exhaustion | zero-normal-UMI | post | 9 | +0.27 | 0.49 |
| CLDN4 vs T/NK cytotoxicity | all-epithelial | post | 12 | −0.31 | 0.32 |
| CLDN4 vs T/NK exhaustion | all-epithelial | post | 12 | −0.38 | 0.23 |
| CLDN4 vs CD8 cytotoxicity | malignant-like | post | 7 | −0.07 | 0.88 |
| CLDN4 vs CD8 exhaustion | malignant-like | post | 7 | −0.04 | 0.94 |
| CLDN4 vs CD8 cytotoxicity | all-epithelial | post | 12 | −0.48 | 0.11 |
| CLDN4 vs CD8 exhaustion | all-epithelial | post | 12 | −0.42 | 0.17 |
| CLDN4 %pos vs T/NK cytotoxicity | malignant-like | post | 7 | −0.04 | 0.94 |
| CLDN4 %pos vs T/NK exhaustion | malignant-like | post | 7 | +0.18 | 0.70 |
| CLDN4 vs T/NK cytotoxicity | malignant-like | post NMPR-only | 6 | −0.03 | 0.96 |
| CLDN4 vs T/NK exhaustion | malignant-like | post NMPR-only | 6 | +0.09 | 0.87 |
| CLDN4 vs T/NK cytotoxicity | all-epithelial | post adeno | 6 | −0.20 | 0.70 |
| CLDN4 vs T/NK cytotoxicity | all-epithelial | post squamous | 6 | −0.37 | 0.47 |
| CLDN4 vs T/NK cytotoxicity | malignant-like | pre+post | 9 | +0.40 | 0.29 |
| CLDN4 vs T/NK exhaustion | malignant-like | pre+post | 9 | −0.27 | 0.49 |
| CLDN4 vs T/NK cytotoxicity | all-epithelial | pre+post | 15 | −0.09 | 0.75 |
| CLDN4 vs T/NK exhaustion | all-epithelial | pre+post | 15 | −0.38 | 0.16 |

Secondary MPR vs NMPR (not the named question): malignant-like CLDN4 cannot be tested (MPR n=1). All-epithelial CLDN4, 8 vs 4, median 1.42 vs 1.46, MWU p=0.81. T/NK cytotoxicity p=0.93; exhaustion p=1.00.

## Extra figure

`results/fig_extra_gse207422_cldn4_exh.png` — four patient-level scatters (NMPR red, MPR blue; adeno circle, squamous square).

| Panel | Test | n | ρ | p |
|---|---|---:|---:|---:|
| A | malignant-like CLDN4 vs T/NK cytotoxicity | 7 | 0.00 | 1 |
| B | malignant-like CLDN4 vs T/NK exhaustion | 7 | 0.00 | 1 |
| C | all-epithelial CLDN4 vs T/NK cytotoxicity | 12 | −0.31 | 0.32 |
| D | all-epithelial CLDN4 vs T/NK exhaustion | 12 | −0.38 | 0.23 |

## Verdict

Primary post malignant-like (**n=7**, not 12, not 92,330): CLDN4 vs cytotoxicity **ρ=0.00, p=1**; vs exhaustion **ρ=0.00, p=1**. Five post tumors have fewer than 10 malignant-like cells (three of four MPR/pCR samples). The all-epithelial sensitivity uses every post sample and is still a weak negative, non-significant association (ρ≈−0.3 to −0.4). The strongest trend is all-epithelial CLDN4 vs CD8 cytotoxicity (n=12, ρ=−0.48, p=0.11) and is still not p<0.05. This does not support a GSE207422 claim that malignant CLDN4 tracks T/NK exhaustion or cytotoxicity.

What this is **not**: TACSTD2 vs T/NK fraction, NMPR>MPR TACSTD2, or a multi-cohort meta (those already exist elsewhere). TACSTD2 is extracted only so the panel is complete; it is not tested here.

## Methods (short)

- Public files: `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (92,330 cells) and the GEO xlsx. Raw HRA001033 was not used.
- Per cell: `log1p(CP10k) = log1p(UMI / library_size × 10⁴)`.
- Lineage = argmax of mean log1p marker modules (epithelial / T/NK / B / myeloid / mast / endo / fibro).
- T/NK = that lineage. CD8 sensitivity = T/NK with CD8A or CD8B UMI > 0.
- Malignant-like (primary) as above. Zero-normal-UMI = epithelial and zero UMI of `SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3` (A3-style). All-epithelial = no normal-lung gate.
- Spearman two-sided on patients. Mann–Whitney on MPR vs NMPR is recorded in `association_statistics.tsv` but is not a primary contrast when MPR n<2 under the malignant gate.

## Files

- `FINDING.md` — this note (table is the deliverable)
- `results/fig_extra_gse207422_cldn4_exh.png` / `.pdf` — extra figure
- `results/per_patient.tsv` — one row / sample, all scores and counts
- `results/honest_n.tsv` — occupancy
- `results/association_statistics.tsv` — every Spearman / MWU
- `results/summary.json`
- Scripts: `download.py`, `extract.py`, `analyze.py`

## Reproduce

```bash
python3 methods/gse207422_cldn4_exh/download.py
python3 methods/gse207422_cldn4_exh/extract.py
python3 methods/gse207422_cldn4_exh/analyze.py
```

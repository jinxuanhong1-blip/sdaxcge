# Claim A3 — GSE207422 (exact)

**Claim tested:** malignant-only TACSTD2 is higher in NMPR than MPR, and per-patient malignant TACSTD2 vs T/NK fraction Spearman ρ is about −0.40 to −0.50.

## Verdict (honest)

- **NMPR > MPR (malignant TACSTD2, 12 post-treatment patients): NOT SUPPORTED.** Median malignant mean-log1p-CPM NMPR=1.972 vs MPR=1.846 (n=8/4). Direction is weakly NMPR>MPR; one-sided Mann–Whitney p=0.533, two-sided p=1.000.
- **Malignant % TACSTD2+** is also n.s. (NMPR median 87.3% vs MPR 78.1%; p_greater=0.276).
- **Malignant pseudobulk CPM goes the other way** (NMPR median 919 vs MPR 1427; p_greater=0.923). Driven by tiny MPR malignant n (P06/P11/P14 have 4–7 cells).
- **All-epithelial sensitivity** is the closest to the claim (NMPR 1.660 vs MPR 1.285; p_greater=0.055) but still not p<0.05.
- **ρ ≈ −0.40 to −0.50: NOT SUPPORTED on the primary analysis.** All 12 post-treatment: Spearman ρ=-0.238 (p=0.457) for malignant mean log1p CPM vs T/NK fraction. %pos ρ=-0.263; pseudobulk ρ=6.99e-03; epithelial ρ=-0.161.
- One post-hoc slice (drop pCR P06) gives ρ=-0.509 (p=0.110, n=11) — inside the claimed band but n.s. and not pre-specified.
- The only MPR patient with a large malignant compartment (P03, 871 cells) has the **lowest** malignant TACSTD2 (0.83). Three other MPR/pCR samples have 4–7 marker-malignant cells because residual epithelium scores as normal lung — expected after MPR, and a hard limit without author CopyKAT labels.

## Data and annotations

- GEO **GSE207422** (Hu et al., *Genome Medicine* 2023, PMID 36869384): BD Rhapsody UMI matrix, 92,330 cells × 24,292 genes, 15 samples / 15 patients.
- **Author per-cell annotations are not public.** Searched: GEO suppl (UMI + sample xlsx only), paper Additional files 1/3/4 (clinical / module genes / steroids; no barcodes), TISCH2 NSCLC gallery (GSE207422 absent), CELLxGENE. Analysis uses marker-inferred labels.
- Lineage = argmax of mean log1p(CPM) across canonical panels (T, NK, B, plasma, myeloid, neutrophil, pDC, mast, stromal, epithelial).
- **Malignant** = assigned epithelial and *not* normal-lung (SFTPA2/SFTPC/AGER/SCGB1A1/TPPP3/FOXJ1), matching the authors' normal clusters (alveolar / club / ciliated). This is **not** CopyKAT.
- Primary test uses **post-treatment surgery** samples only (paper Fig. 1). pCR (P06) is grouped with MPR as in the paper (MPR n=4, NMPR n=8). Pre-treatment biopsies (P01 NE, P05/P08 NMPR-labeled) are excluded from the primary contrast.

## Post-treatment per-patient table

| patient | response | n_mal | n_epi | T/NK frac | mal TACSTD2 mean | mal %pos | mal pb CPM | epi TACSTD2 mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| P03 | MPR | 871 | 928 | 0.545 | 0.827 | 76.1 | 191 | 0.847 |
| P06 | pCR | 7 | 67 | 0.629 | 2.447 | 85.7 | 2075 | 1.427 |
| P11 | MPR | 5 | 447 | 0.449 | 1.871 | 80.0 | 1542 | 1.530 |
| P14 | MPR | 4 | 192 | 0.516 | 1.821 | 75.0 | 1311 | 1.143 |
| P02 | NMPR | 6 | 111 | 0.541 | 1.112 | 50.0 | 576 | 1.076 |
| P04 | NMPR | 80 | 150 | 0.669 | 2.155 | 86.2 | 1275 | 1.742 |
| P07 | NMPR | 4797 | 5318 | 0.124 | 2.302 | 99.1 | 1165 | 2.224 |
| P09 | NMPR | 283 | 283 | 0.676 | 1.789 | 88.3 | 875 | 1.789 |
| P10 | NMPR | 194 | 196 | 0.619 | 1.559 | 73.7 | 847 | 1.551 |
| P12 | NMPR | 409 | 481 | 0.224 | 0.872 | 60.4 | 292 | 0.952 |
| P13 | NMPR | 11 | 29 | 0.201 | 2.338 | 100.0 | 1016 | 1.940 |
| P15 | NMPR | 48 | 554 | 0.190 | 2.256 | 100.0 | 962 | 1.577 |

## QC (marker restriction)

- TACSTD2+ : malignant 90.1%, normal epithelium 83.7%, T 1.3%, NK 1.6%.
- EPCAM+ : malignant 91.8%, T 0.8%.
- PTPRC/CD45+ : malignant 4.8%, T 93.4%.
- Assigned counts: malignant 8724, normal epithelium 2209, T 30543, NK 7993.

## Files

- `per_patient_metrics.csv` — per-sample counts, fractions, TACSTD2 summaries (all 15 samples).
- `summary.json` / `tests.csv` — all tests (mean / %positive / pseudobulk; all-12 / min5 / min20 / drop-pCR / epithelial-wide).
- `claim_A3_summary.png` / `.pdf` — NMPR vs MPR, scatter vs T/NK, lineage restriction.
- `composition_post.png` — malignant and T/NK fractions by post-treatment patient.

## Reproduce

```bash
python3 scripts/claim_A3_gse207422.py
```

Downloads the GEO UMI matrix (~184 MB) and sample metadata if missing. Streams requested genes plus per-cell library size; caches `data/gse207422_marker_stream.npz` (gitignored).

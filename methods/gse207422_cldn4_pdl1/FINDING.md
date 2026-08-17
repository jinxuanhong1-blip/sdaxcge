# GSE207422 — malignant CLDN4 vs CD274/HLA (same cells) + T/NK

**CLDN4 only.** TACSTD2 is never a gate. Dual-high was not run. The given A3 TACSTD2 analysis on this public UMI is not re-argued.

**Verdict (honest n):** On A3-malignant **patient-pseudobulk**, CLDN4 vs CD274 directionally positive but NS (ρ=0.56, p=0.093, n=10) and CLDN4 vs MHC-I near-null (ρ=0.19, p=0.6, n=10). Pooled-cell ρ is descriptive only (n_cells=4,683 from 10 patients) and is near-zero / weakly negative. Patient-level CLDN4 vs T/NK is near-null (ρ=-0.13, p=0.73, n=10; epithelial complete-case ρ=-0.06, p=0.86, n=12).

## Data and n

- Public GEO UMI only: **92,330** cells × **24,292** genes. Raw GSA-Human HRA001033 was not used. Author CopyKAT / epithelium RDS barcodes are not on GEO.
- **Patient is the inference unit** for every test that reports a p-value (12 post-treatment patients: MPR n=4 including pCR P06; NMPR n=8). The three pre-treatment biopsies are excluded from tests.
- Marker lineages (Hu canonical argmax): epithelial 11,019; A3-malignant-like 6,627 (4,683 in the 12 post samples); T/NK 36,323.
- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the given A3 TACSTD2 slice; **not** CopyKAT). Post patients with 0 A3-malignant cells (dropped from malignant tests): **P11, P14**. Post patients with 1–9 A3-malignant cells (kept, noisy): **P06, P13, P15**. A3-malignant n=10/12. Epithelial complete-case n=12/12.
- MHC-I cassette used: **HLA-A/HLA-B/HLA-C/B2M**. Genes missing from the UMI: **none**.
- CD274 is sparse in A3-malignant post cells: 15.4% UMI≥1 (CLDN4 85.6% UMI≥1).

## Definitions

| Item | Rule |
|---|---|
| CLDN4 / CD274 / HLA | `log1p(CP10k)` from public UMI |
| MHC-I | mean of HLA-A, HLA-B, HLA-C, B2M `log1p(CP10k)` (genes present) |
| same-cell pooled | Spearman inside A3-malignant (or epithelial) cells; **p not used** |
| within-patient | Spearman per post patient with ≥20 A3-malignant cells; median ρ |
| patient-pseudobulk | Spearman of per-patient mean `log1p(CP10k)` (honest same-cell n) |
| T/NK | lineage T or NK / all cells, per patient |
| TACSTD2 | not used |

## Primary table — same-cell CLDN4 vs CD274/HLA (honest n)

Pooled-cell ρ treats thousands of cells as independent and is **not** the claim. Patient-pseudobulk and within-patient median are the honest rows.

| contrast | unit | n | result |
|---|---|---:|---|
| CLDN4 vs CD274 | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=-0.02; 4,683 cells (p not used) |
| CLDN4 vs MHC-I | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=-0.10; 4,683 cells (p not used) |
| CLDN4 vs HLA-A | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=-0.23; 4,683 cells (p not used) |
| CLDN4 vs HLA-B | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=-0.11; 4,683 cells (p not used) |
| CLDN4 vs HLA-C | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=-0.06; 4,683 cells (p not used) |
| CLDN4 vs B2M | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=0.15; 4,683 cells (p not used) |
| CLDN4 vs HLA-DRA | pooled A3-malignant cells (descriptive) | 10 patients / 4683 cells | ρ=-0.10; 4,683 cells (p not used) |
| CLDN4 vs CD274 | within-patient median (≥20 cells) | 6 | median ρ=0.13 (IQR 0.02–0.17); 2/6 negative; n=6, Wilcoxon p=0.16 |
| CLDN4 vs MHC-I | within-patient median (≥20 cells) | 6 | median ρ=0.12 (IQR 0.07–0.13); 1/6 negative; n=6, Wilcoxon p=0.062 |
| CLDN4 vs HLA-A | within-patient median (≥20 cells) | 6 | median ρ=0.10 (IQR 0.09–0.14); 1/6 negative; n=6, Wilcoxon p=0.062 |
| CLDN4 vs HLA-B | within-patient median (≥20 cells) | 6 | median ρ=0.08 (IQR 0.07–0.09); 1/6 negative; n=6, Wilcoxon p=0.062 |
| CLDN4 vs HLA-C | within-patient median (≥20 cells) | 6 | median ρ=0.17 (IQR 0.08–0.25); 1/6 negative; n=6, Wilcoxon p=0.062 |
| CLDN4 vs B2M | within-patient median (≥20 cells) | 6 | median ρ=0.02 (IQR -0.04–0.07); 3/6 negative; n=6, Wilcoxon p=0.84 |
| CLDN4 vs HLA-DRA | within-patient median (≥20 cells) | 6 | median ρ=0.12 (IQR 0.07–0.18); 1/6 negative; n=6, Wilcoxon p=0.062 |
| CLDN4 vs CD274 | patient-pseudobulk A3-malignant | 10 | ρ=0.56, p=0.093, n=10 |
| CLDN4 vs MHC-I | patient-pseudobulk A3-malignant | 10 | ρ=0.19, p=0.6, n=10 |
| CLDN4 vs HLA-A | patient-pseudobulk A3-malignant | 10 | ρ=0.58, p=0.082, n=10 |
| CLDN4 vs HLA-B | patient-pseudobulk A3-malignant | 10 | ρ=0.03, p=0.93, n=10 |
| CLDN4 vs HLA-C | patient-pseudobulk A3-malignant | 10 | ρ=0.04, p=0.91, n=10 |
| CLDN4 vs B2M | patient-pseudobulk A3-malignant | 10 | ρ=0.47, p=0.17, n=10 |
| CLDN4 vs HLA-DRA | patient-pseudobulk A3-malignant | 10 | ρ=0.02, p=0.96, n=10 |

Within-patient patients (≥20 A3-malignant cells): **P03, P04, P07, P09, P10, P12**.

## Primary table — patient-level CLDN4 vs T/NK

| CLDN4 score | vs T/NK | n |
|---|---|---:|
| A3-malignant mean log1p(CP10k) | ρ=-0.13, p=0.73, n=10 | 10 |
| A3-malignant %pos | ρ=-0.33, p=0.35, n=10 | 10 |
| epithelial mean log1p(CP10k) (complete-case) | ρ=-0.06, p=0.86, n=12 | 12 |
| epithelial %pos | ρ=-0.29, p=0.35, n=12 | 12 |

## Sensitivity — epithelial patient-pseudobulk (n=12 complete)

| contrast | n | result |
|---|---:|---|
| CLDN4 vs CD274 | 12 | ρ=0.03, p=0.91, n=12 |
| CLDN4 vs MHC-I | 12 | ρ=0.03, p=0.91, n=12 |
| CLDN4 vs HLA-DRA | 12 | ρ=0.38, p=0.23, n=12 |

## Honest limits

1. n=12 (4 vs 8) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. A3-malignant tests drop empty residual tumors (P11, P14); that n is not patched.
2. Pooled same-cell p-values are not reported as inference. Cells from one patient are not independent.
3. Within-patient Spearman needs ≥20 A3-malignant cells. Patients below that floor are excluded from the median-ρ row, not imputed.
4. CD274 (PD-L1) RNA is sparse in scRNA. A near-null same-cell ρ can be dropout, not biology.
5. This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak into A3-malignant-like.
6. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.
7. Expression is log1p(CP10k) from the public UMI, not author Seurat-normalized values, and not PD-L1 / MHC protein.

## Figures

- `figures/fig_same_cell.png` — pooled hexbin (descriptive) + patient-pseudobulk CLDN4 vs CD274
- `figures/fig_vs_tnk.png` — patient CLDN4 vs T/NK and vs MHC-I pseudobulk
- `figures/fig_within_patient.png` — within-patient median ρ strip

## Files

- `tables/per_patient.tsv` — 15 samples; tests use the 12 post rows
- `tables/same_cell_pooled.tsv` / `tables/within_patient.tsv` / `tables/within_patient_summary.tsv`
- `tables/patient_pseudobulk.tsv` / `tables/vs_tnk.tsv`
- `sanity.json` / `summary.json` / `file_manifest.tsv`
- Scripts: `scripts/download.py`, `scripts/extract.py`, `scripts/analyze.py`

## Reproduce

```bash
python3 methods/gse207422_cldn4_pdl1/scripts/download.py
python3 methods/gse207422_cldn4_pdl1/scripts/extract.py
python3 methods/gse207422_cldn4_pdl1/scripts/analyze.py
```

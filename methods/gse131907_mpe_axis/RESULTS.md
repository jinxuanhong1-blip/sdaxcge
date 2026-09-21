# CLDN4–TACSTD2–ELF3–EpCAM in public MPE scRNA

HRA006761 is controlled access, so the numbers below are from public GSE131907 pleural-effusion samples versus primary lung tumors (Kim et al., Nat Commun 2020).

In primary tumor-state cells, CLDN4 and TACSTD2 co-vary across samples (Spearman ρ = 0.770, n = 10, p = 0.0092). Inside MPE carcinoma-like cells, CLDN4 is co-detected with TACSTD2 and with ELF3. The MPE sample count is too small for a CLDN4-versus-T/NK correlation.

## HRA006761

GSA-Human accession HRA006761 (BioProject PRJCA023797; “Single-cell transcriptomics of malignant pleural effusion in patients with advanced non-small-cell lung cancer”) is marked controlled access / request data. The data-access committee is HDAC002197 (Caicun Zhou). No expression matrix was available to download, so no cell from that study was reanalyzed here.

The publication (Clin Transl Med 2024;14:e1649, doi:10.1002/ctm2.1649, PMID 38629624) states that CLDN4 was positively correlated with ELF3, EpCAM, and TACSTD2 in recurrent MPE. That sentence is the authors’ claim. It is not a coefficient recomputed in this repository.

## GSE131907 compartments

The matrix has 208,506 cells. Tests below use sample or cell counts inside a named compartment, not 208,506 as n.

Five pleural-effusion samples (EFFUSION_06, 11, 12, 13, 64; 20,304 cells) are the MPE subset. Eleven `tLung` resections are the primary subset. Kim et al. assigned `Cell_subtype` “Malignant cells” to 0 PE cells and 0 primary-tumor cells. Those labels sit in advanced biopsies and brain metastases and are the cells behind the already reported GSE131907 patient correlation. They are left untouched here.

PE epithelial cells (`Cell_type` Epithelial, subtype NA): **396**. They split on markers:

| Gate | Rule | Cells |
|---|---|---:|
| Carcinoma-like | EPCAM > 0, WT1 = 0, CALB2 = 0 | 259 |
| Mesothelial-like | EPCAM = 0 and (WT1 > 0 or CALB2 > 0) | 76 |
| Outside both gates | mostly EPCAM-negative without WT1/CALB2; 47 of these are in EFFUSION_06 | 61 |

Primary tumor states `tS1`/`tS2`/`tS3`: **6,352** cells. `LUNG_T09` has 5 of them and drops out of every sample test that requires at least 20 compartment cells (10 primary samples remain).

Mesothelial-like PE cells are CLDN4-detected in 2.6% of cells, TACSTD2 in 6.6%, EPCAM in 0%, CALB2 in 96%, and WT1 in 72%. Carcinoma-like PE cells are CLDN4-detected in 83.4%, TACSTD2 in 83.8%, and ELF3 in 80.7% (EPCAM is 100% by the gate). EFFUSION_13 has 46 mesothelial-like epithelial cells and 0 carcinoma-like cells.

Normal-lung AT2 cells also detect CLDN4 (72.0% of 2,020 cells). Mean log1p(CP10k) is 0.91 in AT2 and 1.62 in primary tS, so detection alone does not mark tumor cells in this atlas.

Normalization for every mean below is `log1p(UMI / full-library UMI × 10000)`. A cell is positive when the raw UMI count is greater than 0.

## Coexpression

Co-detection is a Fisher exact test on UMI > 0. Level correlation is Spearman on the log-normalized values.

### Primary tumor states (the quantitative CLDN4–TACSTD2 result)

Across the 10 primary samples with at least 20 tS cells, sample-mean CLDN4 versus sample-mean TACSTD2 is ρ = 0.770, p = 0.0092. The same sample correlation is ρ = 0.055 (p = 0.88) for ELF3 and ρ = 0.527 (p = 0.12) for EPCAM. Between patients, the partner that tracks CLDN4 is TACSTD2.

Inside the 6,352 tS cells:

| Partner | % detected in CLDN4+ | % detected in CLDN4− | Odds ratio | Fisher p | Spearman, all cells | Spearman, double-positive cells |
|---|---:|---:|---:|---:|---:|---:|
| TACSTD2 | 92.0 | 60.7 | 7.45 | 1.9×10⁻¹⁰² | 0.383 | 0.402 (n = 5,142) |
| ELF3 | 95.8 | 72.1 | 8.82 | 5.7×10⁻⁸⁶ | 0.328 | 0.336 (n = 5,354) |
| EPCAM | 90.6 | 59.1 | 6.69 | 4.9×10⁻⁹⁷ | 0.162 | 0.153 (n = 5,065) |

Cell-level ELF3 co-detection is strong. It does not reappear as a between-sample correlation.

### MPE carcinoma-like cells

Among 259 carcinoma-like PE cells:

| Partner | % detected in CLDN4+ | % detected in CLDN4− | Odds ratio | Fisher p | Spearman, all cells | Spearman, double-positive |
|---|---:|---:|---:|---:|---:|---:|
| TACSTD2 | 88.4 | 60.5 | 5.00 | 3.7×10⁻⁵ | 0.047 (p = 0.45) | −0.025 (n = 191, p = 0.74) |
| ELF3 | 90.7 | 30.2 | 22.6 | 3.2×10⁻¹⁶ | 0.266 (p = 1.4×10⁻⁵) | 0.143 (n = 196, p = 0.046) |

CLDN4 and TACSTD2 are co-detected in these MPE cells. Their log levels do not rise together once both genes are detected. ELF3 is the tighter MPE partner: co-detection and a positive level correlation.

EPCAM cannot be tested inside this gate, because the gate requires EPCAM detection. On the unfiltered 396 PE epithelial cells, which still include the mesothelial-like group, EPCAM is detected in 93.2% of CLDN4+ cells and 30.2% of CLDN4− cells (odds ratio 31.9, Fisher p = 1.1×10⁻⁴¹). That odds ratio mixes the carcinoma-like versus mesothelial-like split with within-carcinoma coexpression.

## MPE versus primary levels

Sample means, Mann–Whitney exact test. The symmetric contrast uses the same EPCAM+ / WT1− / CALB2− gate on both sites, and keeps samples with at least 20 gated cells: 3 MPE samples versus 10 primary samples.

| Gene | Median of sample means, MPE | Primary | Δ (MPE − primary) | p |
|---|---:|---:|---:|---:|
| CLDN4 | 1.71 | 1.60 | +0.12 | 0.69 |
| TACSTD2 | 1.94 | 1.73 | +0.21 | 0.69 |
| ELF3 | 1.52 | 1.93 | −0.41 | 0.16 |
| EPCAM | 2.07 | 1.58 | +0.48 | 0.049 |

CLDN4 and TACSTD2 sample means are not separated by site at this n. The EPCAM p-value is a level comparison among cells that already detect EPCAM, on 3 versus 10 samples.

An all-epithelial contrast that leaves the mesothelial-like cells inside the MPE mean makes ELF3 look lower in MPE (p = 0.014, 4 vs 10). That row is in `mpe_vs_primary.tsv` and is the mesothelial mix, not a carcinoma-program difference.

## CLDN4 and T/NK fraction

T/NK fraction is (T lymphocytes + NK cells) / all cells in the sample. CLDN4 is summarized inside the compartment named in the row. Spearman p-values for n ≤ 7 are exact two-sided permutations.

| Set | Samples | ρ, CLDN4+ % vs T/NK | p | ρ, CLDN4 mean vs T/NK | p |
|---|---:|---:|---:|---:|---:|
| All PE epithelial cells | 5 | −0.20 | 0.78 | +0.10 | 0.95 |
| PE epithelial, ≥20 cells | 4 | −0.80 | 0.33 | −0.80 | 0.33 |
| Carcinoma-like, at least 1 cell | 4 | −0.60 | 0.42 | +0.40 | 0.75 |
| Carcinoma-like, ≥20 cells | 3 | not computed |  | not computed |  |
| Primary tS, ≥20 cells | 10 | −0.38 | 0.28 | −0.19 | 0.60 |

The ≥20 carcinoma-like samples are the ones that can be read as MPE malignant-like:

| Sample | Patient | Carcinoma-like cells | CLDN4+ % | T/NK fraction |
|---|---|---:|---:|---:|
| EFFUSION_06 | P1006 | 158 | 76.6 | 0.674 |
| EFFUSION_11 | P1011 | 72 | 95.8 | 0.662 |
| EFFUSION_12 | P1012 | 20 | 95.0 | 0.494 |

EFFUSION_64 has 9 carcinoma-like cells (CLDN4+ 77.8%, T/NK 0.728). EFFUSION_13 has none. Its T/NK fraction is 0.713, and its epithelial cells are the mesothelial-like group (CLDN4+ 0%). Putting that sample in as a CLDN4-low malignant sample is what produces the epithelial-gate ρ of −0.80.

On the four samples with any carcinoma-like cells, the percent summary and the mean-log summary point in opposite directions (ρ = −0.60 versus ρ = +0.40). Three samples meet the 20-cell minimum, which is below the size this script uses for a Spearman test.

The primary tS row (ρ = −0.38, p = 0.28, n = 10) is a different compartment from the locked GSE131907 author-malignant patient result. It is reported because those primary samples are drawn on the same figure. It is not a re-estimate of that locked correlation.

Four patients have a paired advanced biopsy (EBUS or bronchoscopic lymph node) as well as a pleural effusion. Those biopsies are not primary resections. For P1006, P1011, and P1012 the effusion contains carcinoma-like cells and the paired biopsy expresses the same genes. P1013’s effusion epithelial cells are mesothelial-like; the paired EBUS_13 malignant cells are a separate population.

## What stays as previously reported

This file does not replace the locked public results: CosMx exclusion (not muzzling), the concordant-4 malignant CLDN4+ % versus T/NK correlation, the GSE131907 author-malignant patient row, GSE137244 KL versus KP, the TCGA keratin-adjusted correlations, or the TISMO Tacstd2 ICB count. No Visium same-spot correlation is treated as spatial exclusion. No private KL matrix is included.

## Files

- `results/tables/cell_codetection.tsv`, `cell_spearman.tsv`
- `results/tables/sample_pseudobulk_spearman.tsv`, `mpe_vs_primary.tsv`
- `results/tables/mpe_cldn4_vs_tnk.tsv`
- `results/tables/mpe_carcinoma_like_detail.tsv`, `mpe_mesothelial_like_detail.tsv`, `mpe_sample_detail.tsv`
- `results/tables/identity_cell_summary.tsv`, `sample_compartment.tsv`, `paired_mpe_vs_biopsy.tsv`
- `results/figures/fig_codetection.png`, `fig_cldn4_vs_tnk.png`, `fig_mpe_vs_primary_means.png`, `fig_cldn4_identity.png`
- `results/summary.json`

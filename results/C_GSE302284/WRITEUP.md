# C — GSE302284 public context (TROP2 ADC / EGFR DTP NSCLC)

**Slice:** `methods/gse302284/` + `results/C_GSE302284/` only.
**Thesis:** SKB264 / CLDN4 is taken as given. This slice adds a public GEO series; it does not audit the thesis.
**Accession:** [GSE302284](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284) (Baldacci / Brea / Jänne; public 2025-07-10). Cell Ranger v7.1 filtered H5, GRCh38.

## What is (and is not) on GEO

| Question | Answer |
|---|---|
| Sacituzumab / SKB264 / TROP2-ADC treated RNA? | **No.** Six libraries only. The series abstract discusses sacituzumab govitecan (SG) *in vivo* efficacy vs TROP2 CAR-T, but **no SG / SKB264 / TROP2-ADC treated sample is deposited.** |
| Public analog of SKB264? | **Not as a treatment RNA contrast.** SG is the same antibody lineage / TROP2-ADC class as SKB264; that analog cannot be scored here because those tumors are not on GEO. |
| What *can* be scored? | **Osimertinib vs vehicle** in DFCI282 and PC9 (n_library = 1 vs 1 each), plus descriptive residual tumor vs lymph node from one neoadjuvant EGFR-TKI patient (357). |

This is extra public context around TROP2-high EGFR DTP NSCLC, **not** a SKB264 on-treatment analog.

## Design (pre-specified)

- **Primary public contrast:** osimertinib − vehicle, separately in DFCI282 and PC9.
- **Unit:** library is the biological replicate. Each arm has **n = 1 library**. Cell-level MWU *p* is exploratory (cells are not independent biological replicates). Primary effect size = **pseudobulk log2FC** = log2((CPM_treat+1)/(CPM_ctrl+1)).
- **QC:** Cell Ranger filtered H5, then author-like filters (drop bottom 10% UMI, bottom 10% genes, mito > 25%). Scrublet was not re-run.
- **Patient 357:** epithelial cells only (EPCAM/KRT vs PTPRC). Cell-line models: all QC cells.
- **TJ signature:** mean log1p(CP10k) of CLDN1, CLDN4, CLDN7, F11R, PARD3, OCLN, TJP1.
- **IFN / MHC-I / APM:** PRIORITY6 (IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A), IFN_ISG, MHC1_APM (fixed lists in `methods/gse302284/gene_sets.py`).

## Coverage

| GSM | Library | Model | Treatment | n_raw | n_mito>25% | median mito | n_QC |
|---|---|---|---|---:|---:|---:|---:|
| GSM9101264 | N357 | patient_357 | neoadjuvant_EGFR_TKI_residual | 6469 | 226 | 0.068 | 5317 |
| GSM9101265 | T357 | patient_357 | neoadjuvant_EGFR_TKI_residual | 9274 | 2813 | 0.196 | 5455 |
| GSM9128162 | DFCI282_Ositreated | DFCI282 | osimertinib | 2608 | 1215 | 0.184 | 1175 |
| GSM9128163 | DFCI282_Vehicle | DFCI282 | vehicle | 2421 | 397 | 0.004 | 1793 |
| GSM9128164 | PC9_Osi | PC9 | osimertinib | 1688 | 385 | 0.065 | 1217 |
| GSM9128165 | PC9_Veh | PC9 | vehicle | 1631 | 67 | 0.013 | 1420 |

GEO note: `PC9_Veh` (GSM9128165) has `cell line: PC10` / source_name PC10; library name is PC9_Veh. It is used as the paired vehicle for PC9_Osi. DFCI282 osimertinib loses many cells to mito > 25% (DTP stress); that imbalance is real, not a coding error.

## Osimertinib vs vehicle — direction (honest)

Positive log2FC = higher in osimertinib. n_cells = QC cells used; n_library = **1 vs 1**.

| Model / feature | n_cells osi / veh | log2FC | direction | cell MWU *p* | rank-biserial |
|---|---:|---:|---|---:|---:|
| DFCI282 TACSTD2 | 1175 / 1793 | 2.147 | up | 4.62e-75 | 0.323 |
| DFCI282 CLDN4 | 1175 / 1793 | 1.688 | up | 3.77e-17 | 0.160 |
| DFCI282 TJ signature | 1175 / 1793 | 1.704 | up | 5.91e-82 | 0.411 |
| DFCI282 IFN/ISG | 1175 / 1793 | 0.109 | up | 0.0006 | -0.075 |
| DFCI282 MHC-I/APM | 1175 / 1793 | 0.096 | up | 0.0078 | -0.058 |
| DFCI282 PRIORITY6 | 1175 / 1793 | -0.675 | down | 6.17e-24 | -0.219 |
| PC9 TACSTD2 | 1217 / 1420 | 1.391 | up | 4.31e-76 | 0.405 |
| PC9 CLDN4 | 1217 / 1420 | 1.641 | up | 8.47e-106 | 0.444 |
| PC9 TJ signature | 1217 / 1420 | 0.805 | up | 1.06e-60 | 0.368 |
| PC9 IFN/ISG | 1217 / 1420 | 0.865 | up | 3.43e-65 | 0.382 |
| PC9 MHC-I/APM | 1217 / 1420 | 0.594 | up | 1.64e-07 | 0.118 |
| PC9 PRIORITY6 | 1217 / 1420 | 1.537 | up | 2.01e-23 | 0.221 |

Set-row log2FC is the **mean per-gene pseudobulk log2FC**. Cell MWU / rank-biserial use the per-cell mean log1p score. When those two summaries disagree, the analog table below says **discordant**.

### TJ genes (osi vs veh)

| Model / gene | n_cells osi / veh | log2FC | direction | cell MWU *p* | rank-biserial |
|---|---:|---:|---|---:|---:|
| DFCI282 CLDN1 | 1175 / 1793 | 3.472 | up | 1.98e-53 | 0.178 |
| PC9 CLDN1 | 1217 / 1420 | 1.562 | up | 9.53e-30 | 0.157 |
| DFCI282 CLDN4 | 1175 / 1793 | 1.688 | up | 3.77e-17 | 0.160 |
| PC9 CLDN4 | 1217 / 1420 | 1.641 | up | 8.47e-106 | 0.444 |
| DFCI282 CLDN7 | 1175 / 1793 | 1.228 | up | 4.49e-16 | 0.139 |
| PC9 CLDN7 | 1217 / 1420 | 0.828 | up | 7.34e-29 | 0.247 |
| DFCI282 F11R | 1175 / 1793 | 1.254 | up | 9.34e-11 | 0.092 |
| PC9 F11R | 1217 / 1420 | 0.521 | up | 1.38e-09 | 0.105 |
| DFCI282 PARD3 | 1175 / 1793 | 0.822 | up | 2.76e-13 | 0.109 |
| PC9 PARD3 | 1217 / 1420 | 0.187 | up | 0.0244 | 0.029 |
| DFCI282 OCLN | 1175 / 1793 | 1.897 | up | 6.33e-32 | 0.168 |
| PC9 OCLN | 1217 / 1420 | 0.644 | up | 9.07e-11 | 0.090 |
| DFCI282 TJP1 | 1175 / 1793 | 1.569 | up | 5.29e-68 | 0.310 |
| PC9 TJP1 | 1217 / 1420 | 0.250 | up | 3.05e-06 | 0.074 |

### PRIORITY6 genes (osi vs veh)

| Model / gene | n_cells osi / veh | log2FC | direction | cell MWU *p* | rank-biserial |
|---|---:|---:|---|---:|---:|
| DFCI282 IFI27 | 1175 / 1793 | -0.682 | down | 2.18e-31 | -0.252 |
| PC9 IFI27 | 1217 / 1420 | 2.768 | up | 6.68e-09 | 0.030 |
| DFCI282 OAS2 | 1175 / 1793 | -0.524 | down | 0.5454 | -0.005 |
| PC9 OAS2 | 1217 / 1420 | 1.355 | up | 4.28e-05 | 0.014 |
| DFCI282 IFIT1 | 1175 / 1793 | -0.837 | down | 0.0303 | -0.022 |
| PC9 IFIT1 | 1217 / 1420 | 0.904 | up | 0.0008 | 0.009 |
| DFCI282 MX1 | 1175 / 1793 | -0.085 | down | 0.2290 | 0.020 |
| PC9 MX1 | 1217 / 1420 | 1.983 | up | 1.90e-09 | 0.030 |
| DFCI282 ISG15 | 1175 / 1793 | -0.895 | down | 7.31e-10 | -0.125 |
| PC9 ISG15 | 1217 / 1420 | 1.286 | up | 1.33e-12 | 0.129 |
| DFCI282 HLA-A | 1175 / 1793 | -1.026 | down | 1.09e-17 | -0.163 |
| PC9 HLA-A | 1217 / 1420 | 0.929 | up | 5.45e-17 | 0.184 |

### Analog question (TJ down / IFN up) — not testable for SKB264 here

The SKB264 thesis is taken as given. That analog asks whether TROP2-ADC lowers TJ and raises IFN/MHC-I. **GSE302284 has no TROP2-ADC RNA**, so the analog is not scored. Extra context on **osimertinib vs vehicle** (n_lib = 1 vs 1):

| Model | TJ | IFN/ISG | MHC-I/APM | PRIORITY6 |
|---|---|---|---|---|
| DFCI282 | up (mean gene log2FC 1.704; cell-score Δ 0.254; rank-biserial 0.411); n_cells 1175/1793; cell *p*=5.91e-82 | discordant: mean gene log2FC up (0.109), cell-score Δ down (-0.017; rank-biserial -0.075); cell *p*=0.0006 | discordant: mean gene log2FC up (0.096), cell-score Δ down (-0.052; rank-biserial -0.058); cell *p*=0.0078 | down (mean gene log2FC -0.675; cell-score Δ -0.182; rank-biserial -0.219); cell *p*=6.17e-24 |
| PC9 | up (mean gene log2FC 0.805; cell-score Δ 0.249; rank-biserial 0.368); n_cells 1217/1420; cell *p*=1.06e-60 | up (mean gene log2FC 0.865; cell-score Δ 0.027; rank-biserial 0.382); cell *p*=3.43e-65 | up (mean gene log2FC 0.594; cell-score Δ 0.060; rank-biserial 0.118); cell *p*=1.64e-07 | up (mean gene log2FC 1.537; cell-score Δ 0.084; rank-biserial 0.221); cell *p*=2.01e-23 |

Do not read osi-vs-vehicle as an SKB264 / SG on-treatment result. On this contrast, TJ RNA goes **up** with osimertinib in both models (same direction as TACSTD2 / CLDN4), which is the deposited DTP context, not an ADC knockdown.

## Within-library Spearman (cell-level)

| Library | TACSTD2 vs CLDN4 | TACSTD2 vs TJ | CLDN4 vs IFN/ISG |
|---|---|---|---|
| DFCI282 osi | ρ=0.385, p=6.74e-43, n=1175 | ρ=0.434, p=2.72e-55, n=1175 | ρ=0.129, p=8.58e-06, n=1175 |
| DFCI282 veh | ρ=0.205, p=1.83e-18, n=1793 | ρ=0.270, p=3.24e-31, n=1793 | ρ=0.103, p=1.30e-05, n=1793 |
| PC9 osi | ρ=0.664, p=1.50e-155, n=1217 | ρ=0.765, p=7.18e-235, n=1217 | ρ=0.450, p=9.58e-62, n=1217 |
| PC9 veh | ρ=0.482, p=1.67e-83, n=1420 | ρ=0.588, p=5.37e-133, n=1420 | ρ=0.304, p=7.97e-32, n=1420 |
| pt357 tumor epi | ρ=0.288, p=3.10e-28, n=1401 | ρ=0.329, p=1.12e-36, n=1401 | ρ=0.197, p=1.10e-13, n=1401 |
| pt357 LN epi | undefined (constant gene), n=88 | ρ=0.397, p=0.0001, n=88 | undefined (constant gene), n=88 |

LN epithelial CLDN4 is undetectable (0/88 cells UMI>0), so TACSTD2–CLDN4 ρ is undefined there.

## Patient 357 residual tumor vs LN (descriptive only)

Epithelial cells after neoadjuvant EGFR-TKI. Not a treatment contrast.

| Feature | n_cells tumor / LN | log2FC | direction | cell MWU *p* | rank-biserial |
|---|---:|---:|---|---:|---:|
| TACSTD2 | 1401 / 88 | 3.744 | up | 0.0045 | 0.100 |
| CLDN4 | 1401 / 88 | 8.631 | up | 3.57e-05 | 0.166 |
| TJ signature | 1401 / 88 | 4.077 | up | 6.90e-05 | 0.204 |
| IFN/ISG | 1401 / 88 | -0.302 | down | 0.2106 | -0.079 |
| MHC-I/APM | 1401 / 88 | -0.791 | down | 0.5890 | -0.034 |

## Caveats

- **n_library = 1 vs 1.** No sample-level *p* exists. Cell-level *p* will be small whenever thousands of cells shift even slightly; use log2FC / rank-biserial.
- DFCI282 osimertinib has a large mito>25% drop; remaining cells are a stressed DTP subset.
- No doublet re-call (authors used scrublet).
- Patient 357 is one person, two sites, both residual on EGFR-TKI. LN epithelium n=88; CLDN4 is all zero there.
- mRNA ≠ TROP2 / claudin-4 protein.
- PC9 vehicle GEO annotation inconsistency (PC10 vs PC9).
- Public data only. SRA FASTQ were not re-aligned.

## Reproduce

```bash
python3 methods/gse302284/download.py --cache-dir /tmp/gse302284
python3 methods/gse302284/analyze.py --cache-dir /tmp/gse302284
```

Figures: `results/C_GSE302284/figures/`. Tables: `results/C_GSE302284/tables/`.

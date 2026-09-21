# GSE207422 — NHEJ, STING, and IFN scored separately

**FINAL discordant.** An exhaustive sweep (64 pre-specified contrasts) found no patient-level specification in which CLDN4-low cells have higher IFN, higher STING, and lower NHEJ together. The search is closed for this series. Details are in the sweep section below.

**Additive slice.** Public GEO only. TACSTD2 (TROP2) and CLDN4 are both rows in the UMI, so the split was run. CLDN4 is the primary splitter. TACSTD2 is a second splitter and is never a gate. The three modules share no genes. Dual-high was not run. GSE148071 was not merged. The locked T/NK-fraction result on this series is unchanged.

**Verdict (honest n=6 paired post patients, of 12):** in the same A3-malignant cells, CLDN4-high (within-tumor Q4) and CLDN4-low (Q1) separate the three programs.

| Module | Q4 vs Q1 | median Δ | high>low | exact p | BH q (3 tests) |
|---|---|---:|---:|---:|---:|
| IFN effectors | 0.462 vs 0.382 | +0.062 | 6/6 | 0.031 | 0.094 |
| STING | 0.187 vs 0.120 | +0.065 | 5/6 | 0.063 | 0.094 |
| NHEJ | 0.211 vs 0.198 | +0.029 | 4/6 | 0.84 | 0.84 |

IFN is the only module at the n=6 floor (exact two-sided Wilcoxon cannot go below 0.031 when every sign agrees). Across the three tests that floor is BH q=0.094. STING misses the floor by one patient. NHEJ is flat. Homologous recombination is also flat (3/6, p=1.0). Patient-mean Spearman, same 6 tumors: CLDN4 vs IFN ρ=+0.94 (p=0.0048), vs STING ρ=+0.60 (p=0.21), vs NHEJ ρ=+0.09 (p=0.87).

Removing the sensor/signaling genes from the old 40-gene ISG list does not change the IFN ranking (delta Spearman ρ=1). The earlier same-cell IFN lift is the effector program.

## Sweep — thesis direction was not found

Thesis hunted: CLDN4-low → IFN up, STING up, NHEJ down. Delta is CLDN4-high minus CLDN4-low. A hit required all three signs, with n≥4 patients and a strict majority on each module for paired tests, or n≥8 and |ρ|≥0.3 on each module for a single correlation. Bulk sample quartiles needed both arms n≥4.

**0 / 64 specifications are joint hits. 0 / 64 have the joint sign even before the majority rule.**

Families run on the public UMI, post-treatment A3-malignant cells unless noted:

| Family | What was varied |
|---|---|
| Quartile splits | Q4/Q1, median, tertile, top 20%, decile, detected vs undetected, and the same splits inside CLDN4-positive cells only. Epithelial as well as A3. |
| Continuous | Within-patient Spearman; patient-mean Spearman |
| Pre vs post | Pre biopsies alone (2 tumors pass the cell floor), post, and both timings together |
| Response | NMPR, MPR, squamous, adeno, RECIST PR, RECIST SD |
| UMI matching | Log-UMI calipers 0.10, 0.25, 0.50; UMI quintiles and deciles; residual on log UMI |
| Ribosomal depth | 99 cytosolic RPS/RPL genes (median fraction 0.073 in post A3 cells). Non-ribosomal library size, residual on ribosome fraction, residual on log UMI plus ribosome fraction, drop above-median ribosome fraction |
| Subtypes | Non-cycling, cycling, keratin-high, keratin-low, EPCAM-high, squamous-like (KRT5/KRT17/TP63), adeno-like (KRT7). KRT6A is absent |
| Double-low | TACSTD2-low ∩ CLDN4-low vs double-high, and CLDN4 split inside TACSTD2-low cells |
| Bulk | 24 pre-treatment biopsies: all, NMPR, MPR, squamous, adeno, EPCAM partial, sample Q4/Q1, double-low |

The original Q4 vs Q1 contrast stays IFN-high in CLDN4-high cells (median Δ +0.062, 6/6). In four of the six tumors the Q1 tail is mostly CLDN4-undetected cells (P04, P10, P12 are 100% zeros in Q1; P03 is 95%). P07 is the tumor whose Q1 is mostly detected (q25 = 1.19) and its IFN delta is +0.006.

Inside CLDN4-positive cells only, the median split flips IFN and STING (median Δ −0.017 and −0.003; 6/6 and 4/6) and also flips NHEJ (median Δ −0.048; 0/6 in the thesis direction). Q4 vs Q1 among positive cells is the same pattern (IFN 5/6 negative, STING 5/6 negative, NHEJ 0/6 thesis). CLDN4-low expressors are higher on all three modules. That is not IFN/STING up with NHEJ down.

Caliper UMI matching leaves IFN and STING higher in CLDN4-high cells (median IFN Δ +0.022 at caliper 0.10). Removing ribosomal depth does the same (non-ribosomal library, IFN median Δ +0.060, 6 tumors). NMPR-only and squamous-only repeats stay IFN-positive in the CLDN4-high tail. Pre-treatment scRNA has two eligible tumors, both with IFN higher in CLDN4-high cells. Bulk CLDN4 vs IFN is ρ=+0.17 (n=24); the squamous bulk IFN ρ is −0.16 and the STING ρ in that stratum is +0.15.

Per-specification signs are in `tables/sweep_specs.tsv` and `figures/fig_sweep_signs.png`. Patient deltas are in `tables/sweep_patient_deltas.tsv`.

## Depth

CLDN4-Q4 cells carry more UMI (14582 vs 6335, median Δ +8144, 5/6, p=0.063). P07 is the depth-matched tumor (3209 malignant cells; UMI 12938 vs 12966). Its IFN delta is +0.006.

A within-patient linear residual of each module on log1p(UMI), then the same Q4 vs Q1 contrast:

| Module | median residual Δ | residual >0 | p |
|---|---:|---:|---:|
| IFN effectors | +0.015 | 6/6 | 0.031 |
| STING | +0.012 | 4/6 | 0.22 |
| STING without IRF3 | +0.009 | 4/6 | 0.56 |
| NHEJ | −0.029 | 1/6 | 0.063 |
| OXPHOS | −0.064 | 1/6 | 0.31 |

The IFN direction remains after the residual, and the median delta shrinks from +0.062 to +0.015. STING does not. NHEJ leans the other way and does not clear the n=6 floor. P04 (51 malignant cells) carries the largest residual IFN delta (+0.116). The continuous within-patient partial Spearman of CLDN4 vs IFN, also controlling library size, is near zero (median ρ=+0.016, 3/6, p=0.44).

## TACSTD2 splitter (same 6 patients)

Unadjusted CP10k repeats the CLDN4 pattern: IFN 6/6 (median Δ +0.051, p=0.031), STING 5/6 (p=0.063), NHEJ 3/6 (p=1.0). Proliferation is lower in TACSTD2-Q4 (0/6 higher, p=0.031). UMI is higher in 6/6 (p=0.031).

After the UMI residual, TACSTD2-IFN is 4/6 (median Δ +0.008, p=0.44). TACSTD2-NHEJ residual is lower in 6/6 (median Δ −0.044, p=0.031), and the proliferation residual is lower in 6/6 as well.

## Immune compartments

Malignant CLDN4 vs the same modules inside other cells, post patients with ≥20 malignant and ≥20 cells in the compartment (n=6):

| Contrast | ρ | p |
|---|---:|---:|
| malignant CLDN4 vs myeloid STING | +0.09 | 0.87 |
| malignant CLDN4 vs myeloid IFN | +0.43 | 0.40 |
| malignant CLDN4 vs T/NK STING | +0.26 | 0.62 |
| malignant CLDN4 vs T/NK IFN | +0.43 | 0.40 |
| malignant CLDN4 vs T/NK NHEJ | +0.37 | 0.47 |

TMEM173 detection is 11% of post A3-malignant cells and 46% of post myeloid cells. CGAS is 26% and 30%. IRF3 is 47% of malignant cells, so the primary STING score is not a pure adaptor measurement. Dropping IRF3 leaves the unadjusted STING contrast at 5/6, p=0.094.

## Bulk baseline (secondary)

`GSE207422` bulk log2TPM is 24 pre-treatment biopsies. It is not the 12 post-surgery scRNA tumors. Response was not tested.

CLDN4 vs EPCAM ρ=+0.78 (p=5.7×10⁻⁶). TACSTD2 vs EPCAM ρ=+0.70 (p=1.4×10⁻⁴). CLDN4 vs TACSTD2 ρ=+0.86.

| Splitter | module | ρ | p | ρ partial EPCAM | p |
|---|---|---:|---:|---:|---:|
| CLDN4 | NHEJ | +0.25 | 0.23 | −0.21 | 0.31 |
| CLDN4 | STING | +0.17 | 0.43 | +0.07 | 0.76 |
| CLDN4 | IFN effectors | +0.17 | 0.43 | +0.06 | 0.79 |
| TACSTD2 | NHEJ | +0.35 | 0.090 | +0.04 | 0.86 |
| TACSTD2 | STING | +0.24 | 0.25 | +0.18 | 0.40 |
| TACSTD2 | IFN effectors | +0.24 | 0.26 | +0.17 | 0.43 |

The identical CLDN4–STING and CLDN4–IFN rank correlations (both ρ=+0.1696) are a tie in the ranks. The two module scores correlate with each other at ρ=+0.61, so they are not the same vector. Bulk CLDN4 vs OXPHOS is ρ=+0.50 (p=0.013) and falls to ρ=+0.08 (p=0.70) after EPCAM.

## Data and rules

- Public UMI: **92,330** cells × **24,292** genes. Author CopyKAT barcodes are not on GEO.
- Presence gate: TACSTD2 total UMI 114,638 across all cells (12.6% of cells; 88.9% of post A3-malignant). CLDN4 total UMI 80,868 (12.2% of cells; 85.6% of post A3-malignant).
- Symbols in the UMI: **CGAS** and **TMEM173**. Absent: NHEJ1, MB21D1, STING1, TREX1. NHEJ is 11 genes (no NHEJ1). STING negative-regulator module is ENPP1 only (0.6% of post malignant cells) and was not interpreted.
- Bulk symbols differ: MB21D1 is present, the CGAS symbol is absent, NHEJ1 is present, PAXX is absent.
- A3-malignant-like = epithelial marker-argmax AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3. Epithelial 11,019; A3-malignant 6,627. This is not Hu et al. CopyKAT.
- Unit of every primary test is the patient. Within each post patient, malignant cells are split on log1p(CP10k) Q4 vs Q1. Module score = mean log1p(CP10k). Paired Wilcoxon is two-sided exact.
- Floor: ≥20 A3-malignant cells and ≥8 cells in each tail, and the quartile must have spread. **Eligible n=6:** P03 (MPR), P04, P07, P09, P10, P12 (NMPR). Same six for TACSTD2.
- Dropped post (A3-malignant n): P02 (10), P06 (1), P11 (0), P13 (3), P14 (0), P15 (4). Three of four MPR are empty or one cell. Pre-treatment biopsies (P01, P05, P08) are not in the tests.
- Epithelial complete-case, CLDN4 tails, n=11: IFN 9/11 p=0.019; STING 8/11 p=0.067; NHEJ 7/11 p=0.46. Epithelial UMI is higher in 9/11 (p=0.0068). That sensitivity still includes residual normal epithelium.

## Modules

| Module | Genes |
|---|---|
| NHEJ | XRCC6, XRCC5, PRKDC, LIG4, XRCC4, DCLRE1C, PAXX, PNKP, APLF, POLL, POLM. NHEJ1 absent. No MRN. |
| STING | CGAS, TMEM173, TBK1, IKBKE, IRF3 |
| STING without IRF3 | CGAS, TMEM173, TBK1, IKBKE |
| IFN effectors | 34 type-I ISGs. DDX58, IFIH1, IRF7, IRF9, STAT1, STAT2 removed |
| IFN ISG-40 | the previous 40-gene core, kept for comparison |
| HR control | RAD51, BRCA1, BRCA2, PALB2, RAD51C, XRCC2, XRCC3, RPA1 |
| OXPHOS / proliferation | same control idea as the earlier IFN slice |

## Figures

- `figures/fig_paired_cldn4_split.png` — paired Q4 vs Q1 for the three modules
- `figures/fig_resid_umi.png` — same contrast after the UMI residual
- `figures/fig_delta_split.png` — patient deltas against each other
- `figures/fig_detection.png` — detection in malignant vs myeloid cells
- `figures/fig_honest_n.png` — A3-malignant occupancy; floor n=20
- `figures/fig_bulk_cldn4_modules.png` — baseline bulk, n=24
- `figures/fig_sweep_signs.png` — thesis sign for every sweep contrast

## Honest limits

1. Header n stays 12 attempted / 6 tested. MPR residual tumor is the occupancy hole.
2. Exact p=0.031 is the smallest two-sided p at n=6. P07 IFN is a near-tie, and P07 is the tumor whose library sizes match.
3. BH q=0.094 on the three primary paired tests.
4. Scores are CP10k. The residual check is a linear within-patient adjustment, not a claim that depth is fully removed.
5. STING in malignant cells is sparse at TMEM173. IRF3 is a general IRF.
6. Bulk is a different sample set (pre-treatment biopsies) and CLDN4 tracks EPCAM.

## Reproduce

```bash
python3 methods/gse207422_nhej_sting_ifn/scripts/download.py
python3 methods/gse207422_nhej_sting_ifn/scripts/extract.py
python3 methods/gse207422_nhej_sting_ifn/scripts/analyze.py
python3 methods/gse207422_nhej_sting_ifn/scripts/analyze_bulk.py
python3 methods/gse207422_nhej_sting_ifn/scripts/sweep_thesis.py
```

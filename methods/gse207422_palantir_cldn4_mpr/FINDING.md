# GSE207422 — Palantir destinies on A3-malignant-like epithelium, CLDN4-only, MPR labels

**Additive slice.** This is not dual-high. CLDN4 is the readout. TACSTD2 is a companion and is **never a gate**. Barrier / TJ scores **exclude CLDN4**. Terminals were **not** defined as CLDN4-high. The given A3 TACSTD2 / dual-high write-ups are not re-argued.

**Real Palantir, not DPT.** Diffusion maps + Palantir pseudotime + **branch / fate probabilities (destinies)** + entropy (Setty et al., 2019). PAGA is the Leiden graph. A parallel MPR-traj agent may emit DPT only; this run still produces destinies.

**Verdict (honest n).** Palantir auto-detected 4 destinies on 6627 QC A3-malignant-like cells (from 6627 extracted; matrix 92330 cells). Root is CLDN4-low, not CLDN4-high. Destinies largely recover **patient identity** (Cramér's V=0.959) — not a shared CLDN4 lineage. dest_0 n=98 (CLDN4 1.59, barrier 0.63, IFN 0.57), dest_1 n=1356 (CLDN4 0.86, barrier 0.90, IFN 0.51), dest_2 n=4874 (CLDN4 1.63, barrier 0.79, IFN 0.52), dest_3 n=299 (CLDN4 0.59, barrier 0.74, IFN 0.17). Patient-level CLDN4 vs Palantir PT: n=6, ρ=-0.714, p=0.1108. Barrier vs PT: ρ=-0.429, p=0.3965. IFN vs PT: ρ=-0.771, p=0.0724. MPR test n_MPR=1 after the occupancy floor — small n, stated. Eligible: P03(MPR),P04(NMPR),P07(NMPR),P09(NMPR),P10(NMPR),P12(NMPR).

## Data and n

- Public GEO UMI only: **92330** cells × **24292** genes. Author CopyKAT barcodes are not on GEO.
- A3-malignant-like extracted **n=6627**; after QC (n_genes≥200, UMI≥500) **n=6627**.
- Epithelial (pass 1) **n=11019**. A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (**not** CopyKAT).
- Unit of every primary test is the **patient**. Post attempted **n=12**. Eligible (≥20 A3-malignant after QC): **n=6** (5 NMPR, 1 MPR).
- Dropped / below floor post: P02 (NMPR, n=10), P06 (MPR, n=1), P11 (MPR, n=0), P13 (NMPR, n=3), P14 (MPR, n=0), P15 (NMPR, n=4).
- Eligible post: P03(MPR),P04(NMPR),P07(NMPR),P09(NMPR),P10(NMPR),P12(NMPR). dest_1 is TN-only (P05); post fate is identically 0 (Spearman NA).
- Genes absent from barrier module: ['KRT6A', 'KRT6B', 'KRT14']. IFN ISG absent: none. SFTPC / KRT6A / KRT6B / KRT14 are known holes on this public UMI.
- No dual-high (TACSTD2 AND CLDN4) gate. TACSTD2 is reported only as a companion.

## Definitions

| Item | Rule |
|---|---|
| Object | A3-malignant-like epithelium (marker, not CopyKAT) |
| CLDN4 | log1p(CP10k) from raw UMI; readout, not a terminal definition |
| Barrier | TJ (no CLDN4) + simple/basal keratins present on GEO |
| IFN ISG | 40-gene type-I ISG core |
| Root | CLDN4-low tertile, farthest from CLDN4-high PCA centroid; not CLDN4-high |
| Destinies | Palantir auto-detected terminals + fate probabilities |
| Primary n | post patients with ≥20 QC A3-malignant cells |
| MPR | Hu et al. pathologic response; pCR collapsed to MPR; TN has no MPR |

## Root (not CLDN4-high)

Root cell `BD_immune05_738899` in BD_immune05 (TN). CLDN4 log1p(CP10k)=0.735 (tertile **low**; q33=0.970, q66=1.890). Rule: CLDN4-low tertile, farthest from CLDN4-high PCA centroid. `use_early_cell_as_start=True`. Root is CLDN4-high: **False**.

## Palantir destinies

Auto-detected terminals: **4**. PAGA components at connectivity>0: **1** among 12 Leiden vertices.

| Destiny | Terminal cell | Sample | MPR | n_cells | mean CLDN4 | mean barrier | mean IFN |
|---|---|---|---|---:|---:|---:|---:|
| dest_0 | `BD_immune04_867797` | BD_immune04 | NMPR | 98 | 1.591 | 0.631 | 0.569 |
| dest_1 | `BD_immune05_371320` | BD_immune05 | TN | 1356 | 0.860 | 0.902 | 0.515 |
| dest_2 | `BD_immune07_678631` | BD_immune07 | NMPR | 4874 | 1.629 | 0.790 | 0.518 |
| dest_3 | `BD_immune12_381866` | BD_immune12 | NMPR | 299 | 0.594 | 0.739 | 0.172 |

Auto-terminals were **not** set as CLDN4-high. Palantir still placed three terminal cells in the CLDN4-high tertile (dest_0/2/3 terminals). Assigned dest_1 and dest_3 clouds are CLDN4-low (patient P05 TN and P12 NMPR). dest_2 is the majority multi-patient cloud. dest_0 is small and mostly P04.

## CLDN4 + barrier + IFN along destinies

Primary = patient-mean Spearman on post patients with ≥20 A3-malignant cells. Cell-level is exploratory.

| Contrast | n | ρ | p |
|---|---:|---:|---:|
| CLDN4_vs_palantir_pseudotime | 6 | -0.714 | 0.1108 |
| CLDN4_vs_fate_dest_0 | 6 | 0.143 | 0.7872 |
| CLDN4_vs_fate_dest_1 | 6 | NA | NA |
| CLDN4_vs_fate_dest_2 | 6 | 0.486 | 0.3287 |
| CLDN4_vs_fate_dest_3 | 6 | -0.655 | 0.1583 |
| barrier_no_cldn4_vs_palantir_pseudotime | 6 | -0.429 | 0.3965 |
| barrier_no_cldn4_vs_fate_dest_0 | 6 | -0.371 | 0.4685 |
| barrier_no_cldn4_vs_fate_dest_1 | 6 | NA | NA |
| barrier_no_cldn4_vs_fate_dest_2 | 6 | 0.143 | 0.7872 |
| barrier_no_cldn4_vs_fate_dest_3 | 6 | 0.131 | 0.8047 |
| ifn_isg_vs_palantir_pseudotime | 6 | -0.771 | 0.0724 |
| ifn_isg_vs_fate_dest_0 | 6 | -0.029 | 0.9572 |
| ifn_isg_vs_fate_dest_1 | 6 | NA | NA |
| ifn_isg_vs_fate_dest_2 | 6 | 0.600 | 0.2080 |
| ifn_isg_vs_fate_dest_3 | 6 | -0.655 | 0.1583 |

Cell-level (exploratory, post n_cells=4683): CLDN4_vs_palantir_pseudotime ρ=-0.401 p=4.23e-180; barrier_no_cldn4_vs_palantir_pseudotime ρ=-0.355 p=7.11e-139; ifn_isg_vs_palantir_pseudotime ρ=-0.390 p=2.17e-170. dest_1 fate is identically 0 on all post cells (TN-only destiny, P05), so those Spearman rows are NA. These p-values treat cells as independent and are not a claim.

## Destiny vs MPR

Post occupancy floor ≥20 leaves **n=6 patients (5 NMPR, 1 MPR)** of 12 post attempted. Destiny vs MPR is **descriptive**. Do not cite a MPR test as powered. Three of four MPR tumors are nearly empty of A3-malignant cells on GEO. Destiny vs patient Cramér's V=0.959 (p=0.00e+00).

| Contrast | n_NMPR | n_MPR | mean_NMPR | mean_MPR | p | note |
|---|---:|---:|---:|---:|---:|---|
| post_mean_fate_dest_0_NMPR_vs_MPR | 5 | 1 | 0.211 | 0.048 | NA | too_few_samples |
| post_mean_fate_dest_1_NMPR_vs_MPR | 5 | 1 | 0.000 | 0.000 | NA | too_few_samples |
| post_mean_fate_dest_2_NMPR_vs_MPR | 5 | 1 | 0.593 | 0.952 | NA | too_few_samples |
| post_mean_fate_dest_3_NMPR_vs_MPR | 5 | 1 | 0.195 | 0.000 | NA | too_few_samples |
| post_mean_CLDN4_NMPR_vs_MPR | 5 | 1 | 1.308 | 1.024 | NA | too_few_samples |
| post_mean_barrier_no_cldn4_NMPR_vs_MPR | 5 | 1 | 0.725 | 0.624 | NA | too_few_samples |
| post_mean_ifn_isg_NMPR_vs_MPR | 5 | 1 | 0.446 | 0.349 | NA | too_few_samples |
| post_mean_pseudotime_NMPR_vs_MPR | 5 | 1 | 0.516 | 0.474 | NA | too_few_samples |

## Honest limits

1. **Small n.** 12 post patients attempted; **6** meet the occupancy floor; **1** MPR among them. Do not cite n=12 as the tested n.
2. This is not Hu et al. CopyKAT. A3-malignant-like can leak unmarked epithelium.
3. There is no nLung AT2 root: the A3 rule zeros SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3. The external arrow is “not CLDN4-high,” not an AT2→tumor proof.
4. Palantir destinies on a multi-patient tumor subset can recover patient identity (Cramér's V=0.959). That is a limit, not a lineage.
5. Cell-level p-values are exploratory (pseudoreplication).
6. Dual-high (TACSTD2 AND CLDN4) was not run.
7. RNA velocity was not run (no spliced/unspliced on GEO).
8. A MPR-traj / DPT-only analysis is a different object. This FINDING is the Palantir-destiny one.

## Extra figures

- `results/figures/fig_palantir_paga.png — PAGA + UMAP destiny / PT / CLDN4 / barrier / IFN`
- `results/figures/fig_programs_along_destiny.png — CLDN4, barrier, IFN along Palantir PT and destinies`
- `results/figures/fig_destiny_vs_mpr.png — destiny composition and patient-mean fate vs MPR`
- `results/figures/fig_honest_n.png — A3-malignant occupancy; floor visible`
- `results/figures/fig_extra_root_not_cldn4high.png — root marked on CLDN4 / tertile UMAP`
- `results/figures/fig_extra_entropy.png — Palantir entropy vs barrier and IFN`

## Files

- `results/tables/destiny_vs_mpr.tsv — destinies × MPR (cells and patients)`
- `results/tables/patient_destiny.tsv — per-patient fate means and majority destiny`
- `results/tables/programs_along_destiny.tsv — CLDN4 / barrier / IFN vs destinies`
- `results/tables/honest_n.tsv / tests.tsv / destiny_by_sample.tsv`
- `results/tables/leiden_paga_vertices.tsv / paga_connectivities.tsv / root_info.json`
- `results/tables/cell_palantir_scores.tsv.gz`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/gse207422_palantir_cldn4_mpr/requirements.txt
python3 methods/gse207422_palantir_cldn4_mpr/scripts/download.py
python3 methods/gse207422_palantir_cldn4_mpr/scripts/extract.py
python3 methods/gse207422_palantir_cldn4_mpr/scripts/analyze.py
```


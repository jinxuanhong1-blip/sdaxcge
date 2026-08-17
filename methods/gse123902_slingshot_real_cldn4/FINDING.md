# Finding — GSE123902 Laughney tumor epithelium, CLDN4-only real Slingshot/PAGA

ADDITIVE. **CLDN4 only.** GSE123902 (Laughney et al., *Nat Med* 2020, PMID 32042191) human LUAD primary + metastasis epithelium. SuperSeries GSE123904 / mouse GSE123903 are not used. No TACSTD2∩CLDN4 dual-high gate. Author annotated H5 (36.5 GB) skipped; GEO dense UMI CSVs (90.4 MB) are the matrix.

Primary clock: **R_slingshot**. Inferential unit = **donor** (LX ID). Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. IFN is a locked Hallmark IFNα∩IFNγ core. Root is never CLDN4-high.

**What holds (n=13 tumor donors).** CLDN4 vs barrier/keratin (CLDN4 excluded) is n=13, ρ=0.670, p=0.0122. CLDN4 vs IFN is n=13, ρ=-0.280, p=0.354. Along Slingshot PT, barrier is n=13, ρ=0.566, p=0.0438 and IFN is n=13, ρ=-0.082, p=0.789. **What is underpowered / null.** CLDN4 vs Slingshot PT is n=13, ρ=0.385, p=0.194. CLDN4 vs AT2 is n=13, ρ=0.555, p=0.049. n=13 tumor donors is the honest ceiling.

## Verdict

Donor-level CLDN4 vs Slingshot PT: n=13, ρ=0.385, p=0.194. CLDN4 vs DPT: n=7, ρ=0.393, p=0.383. CLDN4 vs AT2 score: n=13, ρ=0.555, p=0.049. CLDN4 vs barrier/keratin (CLDN4 excluded): n=13, ρ=0.670, p=0.0122. CLDN4 vs IFN: n=13, ρ=-0.280, p=0.354. Barrier vs Slingshot PT: n=13, ρ=0.566, p=0.0438. IFN vs Slingshot PT: n=13, ρ=-0.082, p=0.789. Paired CLDN4-high vs low barrier/keratin: n=11, W=0.0, Δmed=0.468, p=0.000977. Paired CLDN4-high vs low IFN: n=11, W=1.0, Δmed=0.088, p=0.00195. Paired CLDN4-high vs low AT2: n=11, W=0.0, Δmed=0.247, p=0.000977. Slingshot engine=R_slingshot; 13 lineage(s) from Leiden 4. PAGA has 2 component(s) at connectivity>0 among 20 Leiden vertices. Root is not CLDN4-high. Donor is the unit. Not a TACSTD2 redo. No both-high gate. 36.5 GB H5 skipped.

## Honest n

- GEO catalog: **17 samples / 14 LX IDs** (8 primary + 5 metastasis + 4 matched-normal files; LX685 is normal-only). Paper QC atlas is 41,384 cells; we do not substitute that n.
- Marker epithelium before cap (EPCAM|KRT8|KRT18|KRT19>0 and PTPRC==0): **n_cells_uncapped = 5534** (tumor 4733, normal 801).
- Analysis cells after QC and ≤350/sample cap: **n_cells = 3610** (tumor 2852, normal 758).
- Donors on the object: **n_donors = 14**. Tumor donors: **13**. Donors with ≥10 tumor epithelial cells used for Spearman: **n = 13**.
- Donors with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 11**.
- CLDN4 tertile cells: {'low': 1661, 'high': 1203, 'mid': 746}.
- Tissue cells: {'PRIMARY': 1832, 'METASTASIS': 1020, 'NORMAL': 758}.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- Slingshot: available=True; engine=R_slingshot; Bioconductor slingshot 2.10.0.
- Lineages: **13**. Lineage table: `results/tables/lineage.tsv`.
- PAGA components at connectivity>0: **2** among 20 Leiden vertices.
- Marker epithelium ≠ author / inferCNV malignant. No ICI / MPR / RECIST labels.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- DPT / Slingshot root: matched-normal AT2-like, CLDN4 not high (median AT2) (root cell index 1045, donor LX675, tissue NORMAL, Leiden 4).
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN genes: STAT1/2, IRF1/7/9, ISG15/20, MX1/2, IFIT1/2/3, OAS1/2/L, IFI27/44/44L/6, RSAD2, CXCL9/10/11, B2M, TAP1, PSMB8/9.

## Lineage table (done criterion)

| lineage | start | end | n_clusters | n_cells | n_donors | mean CLDN4 | mean barrier | mean IFN | donor ρ CLDN4 vs PT |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Lineage1 | 4 | 5 | 6 | 1174 | 13 | 0.469 | 0.749 | 0.438 | n=12, ρ=-0.441 |
| Lineage2 | 4 | 1 | 5 | 929 | 13 | 0.580 | 0.826 | 0.424 | n=10, ρ=0.103 |
| Lineage3 | 4 | 11 | 5 | 866 | 13 | 1.084 | 1.077 | 0.426 | n=9, ρ=0.733 |
| Lineage4 | 4 | 2 | 4 | 843 | 13 | 0.802 | 0.920 | 0.380 | n=9, ρ=0.683 |
| Lineage5 | 4 | 10 | 4 | 781 | 13 | 0.888 | 1.023 | 0.420 | n=10, ρ=0.648 |
| Lineage6 | 4 | 7 | 4 | 726 | 13 | 0.897 | 1.109 | 0.457 | n=9, ρ=0.600 |
| Lineage7 | 4 | 16 | 4 | 844 | 13 | 1.048 | 1.287 | 0.410 | n=9, ρ=0.717 |
| Lineage8 | 4 | 17 | 4 | 897 | 13 | 0.870 | 1.154 | 0.456 | n=9, ρ=0.533 |
| Lineage9 | 4 | 18 | 4 | 680 | 13 | 0.868 | 1.034 | 0.414 | n=7, ρ=0.595 |
| Lineage10 | 4 | 19 | 4 | 962 | 13 | 0.944 | 1.255 | 0.447 | n=10, ρ=0.697 |
| Lineage11 | 4 | 3 | 3 | 393 | 13 | 0.341 | 0.879 | 0.532 | n=5, ρ=-0.051 |
| Lineage12 | 4 | 13 | 3 | 661 | 13 | 0.504 | 1.147 | 0.493 | n=6, ρ=-0.174 |
| Lineage13 | 4 | 14 | 3 | 378 | 13 | 0.349 | 0.858 | 0.502 | n=5, ρ=0.205 |

## Primary (donor-level Spearman on tumor epithelium, BH inside this list)

| Contrast | n_donors | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 13 | 0.385 | 0.194 | 0.389 |
| CLDN4 vs DPT | 7 | 0.393 | 0.383 | 0.548 |
| CLDN4 vs AT2 score | 13 | 0.555 | 0.049 | 0.152 |
| CLDN4 vs barrier/keratin (no CLDN4) | 13 | 0.670 | 0.0122 | 0.122 |
| CLDN4 vs IFN | 13 | -0.280 | 0.354 | 0.548 |
| CLDN4 vs TACSTD2 (comparator) | 13 | 0.533 | 0.0607 | 0.152 |
| barrier/keratin vs Slingshot PT | 13 | 0.566 | 0.0438 | 0.152 |
| IFN vs Slingshot PT | 13 | -0.082 | 0.789 | 0.877 |
| SFTPC vs Slingshot PT (control) | 13 | -0.028 | 0.928 | 0.928 |
| AT2 score vs Slingshot PT (control) | 13 | 0.104 | 0.734 | 0.877 |

## Sensitivity (not in the BH family)

| Contrast | n_donors | ρ | p |
| --- | ---: | ---: | ---: |
| primary-only CLDN4 vs Slingshot PT | 8 | 0.310 | 0.456 |
| met-only CLDN4 vs Slingshot PT | 5 | 0.300 | 0.624 |
| primary-only CLDN4 vs barrier/keratin | 8 | 0.595 | 0.12 |
| met-only CLDN4 vs barrier/keratin | 5 | 0.800 | 0.104 |
| primary-only CLDN4 vs IFN | 8 | -0.929 | 0.000863 |
| met-only CLDN4 vs IFN | 5 | 0.100 | 0.873 |
| tumor-donor CLDN4 vs DPT | 7 | 0.393 | 0.383 |
| all-tissues (include normal) CLDN4 vs Slingshot PT | 14 | 0.459 | 0.0985 |

## Extra figure — CLDN4-high vs CLDN4-low (donor-paired)

Emitted: **True**. Observed Spearman(CLDN4, barrier/keratin no CLDN4) n=13, ρ=0.67, p=0.0122.

| Paired contrast (high − low) | n_donors | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 11 | 0.468 | 0.000977 |
| IFN high vs low | 11 | 0.088 | 0.00195 |
| AT2 high vs low | 11 | 0.247 | 0.000977 |
| Slingshot PT high vs low | 11 | -0.656 | 0.898 |
| DPT high vs low | 5 | 0.011 | 0.188 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- Marker epithelium is **not** a CNV-malignant call. The 36.5 GB author H5 was skipped.
- n=13 tumor donors is small; CIs are wide. Do not write a multi-cohort pile.
- The graph mixes matched normal (root) and tumor/met. Pooled PT is not a within-tumor clock.
- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot/DPT is an ordering, not a developmental clock. No ICI language.

## Outputs

- `results/tables/lineage.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_extra_along_pt.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
bash methods/gse123902_slingshot_real_cldn4/scripts/install_tools.sh
python3 methods/gse123902_slingshot_real_cldn4/scripts/download.py --out /tmp/gse123902_slingshot
python3 methods/gse123902_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse123902_slingshot \
  --out /tmp/gse123902_slingshot/epithelium.h5ad
python3 methods/gse123902_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse123902_slingshot/epithelium.h5ad \
  --outdir methods/gse123902_slingshot_real_cldn4/results \
  --finding methods/gse123902_slingshot_real_cldn4/FINDING.md
```

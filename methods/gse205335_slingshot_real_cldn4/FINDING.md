# Finding — GSE205335 malignant-only REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. **CLDN4 only.** Ahn / Lee ICI biopsy/effusion cohort, GEO [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335) (eLife 98366). Author-malignant cells only (`lineage.sub == Malignant cells`) on non-normal tissues. This folder does **not** redo the winning-pair epithelium trajectory and does **not** substitute a RECIST table for a lineage. No TACSTD2∩CLDN4 dual-high gate. GSE131907 / GSE207422 are not added.

Primary clocks: **Slingshot + scanpy DPT**. Inferential unit = **patient**. Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. Root is lowest-CLDN4 Leiden 8 (median CLDN4 cell); never CLDN4-high cluster 7.

**What holds (n=22 patients).** CLDN4 vs barrier/keratin (CLDN4 excluded) n=22, ρ=0.293, p=0.186. CLDN4 vs IFN n=22, ρ=-0.099, p=0.662. Within-patient CLDN4-high vs low barrier n=21, W=0.0, Δmed=0.473, p=9.54e-07; IFN n=21, W=41.0, Δmed=0.053, p=0.00801. **Trajectory vs programs.** CLDN4 vs DPT n=22, ρ=0.435, p=0.0429; barrier vs DPT n=22, ρ=-0.312, p=0.157; IFN vs DPT n=22, ρ=-0.240, p=0.282. **RECIST.** DPT R vs NR n=6 vs 10, Δmed=0.035, p=0.181. n_R vs n_NR is thin; a null is inconclusive, not “CLDN4 is unrelated to ICI.”

## Verdict

Patient-level CLDN4 vs DPT: n=22, ρ=0.435, p=0.0429. CLDN4 vs Slingshot PT: n=22, ρ=0.634, p=0.00153. CLDN4 vs barrier/keratin (CLDN4 excluded): n=22, ρ=0.293, p=0.186. CLDN4 vs IFN: n=22, ρ=-0.099, p=0.662. barrier vs DPT: n=22, ρ=-0.312, p=0.157. IFN vs DPT: n=22, ρ=-0.240, p=0.282. DPT RECIST R vs NR: n=6 vs 10, Δmed=0.035, p=0.181. Slingshot PT RECIST R vs NR: n=6 vs 10, Δmed=2.473, p=0.263. Paired CLDN4-high vs low barrier/keratin: n=21, W=0.0, Δmed=0.473, p=9.54e-07. Paired CLDN4-high vs low IFN: n=21, W=41.0, Δmed=0.053, p=0.00801. PAGA has 1 component(s) at connectivity>0 among 17 Leiden vertices. Slingshot was run on Harmony-PCA with the non-CLDN4-high start cluster. Not a TACSTD2 redo. No both-high gate. Malignant-only. Patient unit.

## Honest n

- GEO patients / samples: **26 / 33**. Four normal-only patients (P2001, P2009, P2016, P3032) have 0 malignant cells and are out.
- Author-malignant cells on non-normal tissues (catalog): **28512**.
- Analysis cells after QC (capped ≤400/patient): **n_cells = 7552** in **n_patients = 22**.
- Patients with ≥10 malignant cells used for Spearman: **n = 22**.
- RECIST R (PR) / NR (SD+PD) / NE among analysis patients: **6 / 10 / 6**. MPR/NMPR is unlabeled (n=0).
- NSCLC ADC+SQ patients: **17**. SCLC+NUT remain in the primary n.
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 21**.
- Author subtypes (cells): {'Malignant cells': 7552}.
- Histology (cells): {'ADC': 4352, 'SCLC': 1600, 'SQ': 1200, 'NUT': 400}.
- CLDN4 tertile cells: low 2599, mid 2436, high 2517.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- Slingshot: ran=True; available=True; [1] ‘2.10.0’.
- PAGA components at connectivity>0: **1** among 17 Leiden vertices.
- Do not cite n_cells as the inferential n.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=patient_id.
- Root: lowest-CLDN4 Leiden 8 (median CLDN4 cell); never CLDN4-high cluster 7 (root cell index 16, patient P1006).
- Root cluster mean CLDN4 = 0.211; CLDN4-high cluster 7 mean = 2.189.
- Slingshot start.clus = 8.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN: compact Hallmark-like ISG panel (STAT1, IRF1, ISG15, MX1, OAS1, IFIT1, IFIT3, IFI6, BST2, GBP1, CXCL9, CXCL10, IFI27, RSAD2, IFITM1, OAS2, IRF7, IFI44L, MX2, ISG20).

## Lineage (required; not a RECIST-only table)

PAGA has **1** component among 17 Leiden vertices. Slingshot (v2.10.0) drew **9** lineages from start cluster 8 (lowest CLDN4; not cluster 7). These are graph paths, not nine biological fates.

- Lineage1: 8>2>11>5>4>9
- Lineage2: 8>2>11>5>13
- Lineage3: 8>2>1>15>16
- Lineage4: 8>2>6>7 (reaches the CLDN4-high cluster)
- Lineage5: 8>2>1>3
- Lineage6: 8>2>0
- Lineage7: 8>2>10
- Lineage8: 8>2>14
- Lineage9: 8>12

Figures: `results/figures/fig_lineage_paga_slingshot.png` (PAGA + Slingshot curves), `results/figures/fig_along_pseudotime.png` (CLDN4 + barrier + IFN vs DPT).

## Primary (patient-level Spearman, BH inside this list)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs DPT | 22 | 0.435 | 0.0429 | 0.171 |
| CLDN4 vs Slingshot PT | 22 | 0.634 | 0.00153 | 0.0122 |
| CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.293 | 0.186 | 0.297 |
| CLDN4 vs IFN | 22 | -0.099 | 0.662 | 0.662 |
| barrier/keratin vs DPT | 22 | -0.312 | 0.157 | 0.297 |
| IFN vs DPT | 22 | -0.240 | 0.282 | 0.322 |
| AT2 residual vs DPT (control) | 22 | -0.357 | 0.102 | 0.273 |
| CLDN4 vs TACSTD2 (comparator) | 22 | 0.265 | 0.234 | 0.312 |

## DPT / Slingshot vs RECIST (patient unit; not in the BH family)

| Contrast | n_R vs n_NR | median R | median NR | Δ (R−NR) | p |
| --- | ---: | ---: | ---: | ---: | ---: |
| DPT, RECIST R vs NR | 6 vs 10 | 0.317 | 0.282 | 0.035 | 0.181 |
| Slingshot PT, RECIST R vs NR | 6 vs 10 | 78.574 | 76.101 | 2.473 | 0.263 |
| CLDN4, RECIST R vs NR | 6 vs 10 | 1.429 | 0.947 | 0.482 | 0.492 |
| barrier/keratin, RECIST R vs NR | 6 vs 10 | 1.090 | 1.146 | -0.056 | 0.958 |
| IFN, RECIST R vs NR | 6 vs 10 | 0.096 | 0.452 | -0.356 | 0.368 |

RECIST R/NR/NE = **6 / 10 / 6**. NE is dropped from R vs NR. MPR is not labelled and is not substituted.

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| ADC+SQ CLDN4 vs DPT | 17 | 0.429 | 0.0858 |
| ADC+SQ CLDN4 vs barrier/keratin (no CLDN4) | 17 | 0.507 | 0.0376 |
| ADC+SQ CLDN4 vs IFN | 17 | -0.120 | 0.646 |
| ADC+SQ IFN vs DPT | 17 | 0.098 | 0.708 |
| drop n_cells<50 CLDN4 vs DPT | 21 | 0.373 | 0.0961 |
| CLDN4 %pos vs DPT | 22 | 0.444 | 0.0383 |
| ADC+SQ CLDN4 vs Slingshot PT | 17 | 0.586 | 0.0135 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)

Emitted: **True**. Observed Spearman(CLDN4, barrier_keratin_no_CLDN4) n=22, ρ=0.293, p=0.186.

| Paired contrast (high − low) | n_patients | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 21 | 0.473 | 9.54e-07 |
| IFN high vs low | 21 | 0.053 | 0.00801 |
| DPT high vs low | 21 | 0.029 | 1.34e-05 |
| AT2 residual high vs low | 21 | 0.011 | 0.00186 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a winning-pair (GSE131907+GSE205335) redo and not a GSE131907 PAGA redo.
- Author `Malignant cells` is **not CNV**. Residual AT2 score in malignant cells is not a normal AT2 root.
- DPT / Slingshot are orderings, not clocks. Harmony on patient removes a main-effect batch; it does not prove a within-tumor differentiation axis.
- RECIST 6 vs 10 is thin. MPR/NMPR is unlabeled (n=0). Do not write MPR on these figures.
- SCLC / NUT stay in the primary n. ADC+SQ is a sensitivity.
- No TACSTD2∩CLDN4 both-high gate. TACSTD2 is a comparator only.
- Do not write “malignant cells differentiate because PAGA is connected.”
- Do not write “CLDN4 marks the ICI-resistant terminal.”
- A sibling RECIST-traj folder may exist; this folder still emits a real lineage plot.

## Outputs

- `results/tables/lineage_vertices.tsv` — **lineage table (done criterion)**
- `results/tables/patient_means.tsv` — **patient table (done criterion)**
- `results/tables/patient_level_spearman.tsv`
- `results/tables/recist_vs_dpt.tsv`
- `results/tables/paga_connectivities.tsv`
- `results/tables/slingshot_lineages.tsv`
- `results/figures/fig_lineage_paga_slingshot.png`
- `results/figures/fig_along_pseudotime.png`
- `results/figures/fig_dpt_vs_recist.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/gse205335_slingshot_real_cldn4/requirements.txt
# R / Bioconductor slingshot (user library)
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER", "~/R/library"), recursive=TRUE);
  .libPaths(Sys.getenv("R_LIBS_USER", "~/R/library"));
  install.packages("BiocManager", repos="https://cloud.r-project.org");
  BiocManager::install(c("slingshot","SingleCellExperiment","DelayedMatrixStats"), ask=FALSE, update=FALSE)'
python3 methods/gse205335_slingshot_real_cldn4/scripts/download.py \
  --out /tmp/gse205335_slingshot_real
python3 methods/gse205335_slingshot_real_cldn4/scripts/extract_malignant.py \
  --data /tmp/gse205335_slingshot_real \
  --out /tmp/gse205335_slingshot_real/malignant.h5ad
python3 methods/gse205335_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse205335_slingshot_real/malignant.h5ad \
  --outdir methods/gse205335_slingshot_real_cldn4/results \
  --finding methods/gse205335_slingshot_real_cldn4/FINDING.md
```


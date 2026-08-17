# Finding — GSE207422 PAGA / DPT scored by CLDN4 vs MPR

Additive public slice. **CLDN4-only. Not a TACSTD2 redo. No dual-high.** Primary readout is **CLDN4** on marker-defined A3-malignant-like epithelium from GSE207422 (Hu et al., *Genome Med* 2023, PMID 36915183): neoadjuvant PD-1 + platinum, MPR vs NMPR. Counts are the public GEO processed UMI matrix. Author CopyKAT / epithelium barcodes are **not** on GEO. GSE207422-only T/NK was already flat and is **not** re-audited here.

PAGA + leftover/AT2-rooted diffusion pseudotime. Inferential unit = **patient**. Cell-level ρ is descriptive and is not the claim. Barrier/keratin score **excludes CLDN4**. Root is leftover (normal-lung-program) epithelium, never a CLDN4-high cell. R was not present; Slingshot was not run.

## Verdict

Patient-level A3-malignant CLDN4 vs leftover/AT2-rooted DPT: n=7, ρ=-0.214, p=0.645. CLDN4 vs AT2 score (malignant-like): n=7, ρ=0.786, p=0.0362. CLDN4 vs barrier/keratin (CLDN4 excluded): n=7, ρ=-0.071, p=0.879. A3-malignant DPT vs MPR: n=6 vs 1 (too few for exact MWU). Epithelial DPT vs MPR (complete-case companion, not the malignant-like n): mean 0.410 vs 0.545 (Δ=-0.134); exact p=0.214; n=8 vs 4. Paired CLDN4-high vs low barrier/keratin: n=11, W=5.0, Δmed=0.127, p=0.00977. Paired CLDN4-high vs low DPT: n=11, W=10.0, Δmed=0.030, p=0.042. PAGA has 1 component(s) at connectivity>0 among 17 Leiden vertices. Empty A3-malignant post patients: ['P11', 'P14']. Noisy 1–9 cell post patients: ['P06', 'P13', 'P15']. Honest n is the patient n after the ≥10 malignant-like gate. This set is small. Not a TACSTD2 redo. No both-high gate. T/NK was not re-tested.

## Honest n

- GEO samples: **15** (3 pre-treatment biopsies + 12 post-treatment resections). Catalog n is not the test n.
- Public UMI matrix: **92330** cells × **24292** genes.
- Marker epithelial cells: **11019**. A3-malignant-like: **6627**. Leftover epithelium: **4392**.
- After QC (min 200 genes, min 500 UMI): **n_cells_qc = 11019** in **15** samples.
- Graph cap: max 500/sample when n>8000 (seed 0). Analysis object: **n_cells = 4537** in **15** samples. Do not quote this as the inferential n.
- Lineage/compartment on the analysis object: {'A3_malignant': 2314, 'leftover_epi': 2223}.
- Inferential cohort: **12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). Pre-treatment biopsies stay on the graph for the manifold/root and are dropped from every test.
- Post patients with 0 A3-malignant cells on the analysis object (NaN, dropped from malignant-like tests): **P11, P14**.
- Post patients with 1–9 A3-malignant cells (kept only in the noisy sensitivity): **P06, P13, P15**.
- Post patients with ≥10 A3-malignant cells (primary malignant-like Spearman): **n = 7**.
- Post patients with ≥10 epithelial cells (complete-case epithelium companion): **n = 12**.
- Post patients with ≥8 cells in both CLDN4-high and CLDN4-low arms (paired extra): **n = 11**.
- CLDN4 tertile cells on the graph: low 1513, mid 1512, high 1512.
- Genes absent from locked sets: {'AT2': ['SFTPC'], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'single_genes': ['SFTPC']}.
- This set is small. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. Malignant-like MPR residual tumors are often empty. Do not inflate cell count.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- DPT root: leftover epithelium, not CLDN4-high, median AT2 score (root cell index 2121, patient P07, compartment leftover_epi, CLDN4=1.7845048904418945).
- PAGA components at connectivity>0: **1** among 17 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the given A3 slice; **not** CopyKAT).
- Slingshot: not run (R absent).

## Primary — patient-level CLDN4 vs DPT (A3-malignant-like, post, ≥10 cells)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| A3-malignant CLDN4 vs DPT | 7 | -0.214 | 0.645 | 0.806 |
| A3-malignant CLDN4 vs AT2 score | 7 | 0.786 | 0.0362 | 0.181 |
| A3-malignant CLDN4 vs barrier/keratin (no CLDN4) | 7 | -0.071 | 0.879 | 0.879 |
| A3-malignant CLDN4 vs TACSTD2 (companion) | 7 | 0.250 | 0.589 | 0.806 |
| A3-malignant AT2 vs DPT (control) | 7 | 0.357 | 0.432 | 0.806 |

## Primary — patient-level DPT vs MPR (exact two-sided MWU)

| Contrast | n_NMPR | n_MPR | mean NMPR | mean MPR | Δ (NMPR−MPR) | exact p | note |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| NMPR vs MPR: A3-malignant DPT, post ≥10 malignant-like | 6 | 1 | 0.379 | 0.326 | NA | NA | too_few_samples; primary malignant-like gate |
| NMPR vs MPR: A3-malignant DPT, post any malignant-like (≥1 cell) | 8 | 2 | 0.393 | 0.362 | +0.031 | 0.889 | noisy; includes 1–9 cell residuals |
| NMPR vs MPR: epithelial DPT, post complete-case | 8 | 4 | 0.410 | 0.545 | -0.134 | 0.214 | companion; leftover+malignant on the graph |
| NMPR vs MPR: leftover-epi DPT, post ≥10 leftover | 8 | 4 | 0.415 | 0.548 | -0.133 | 0.214 | leftover only |
| NMPR vs MPR: A3-malignant CLDN4, post ≥10 malignant-like | 6 | 1 | 1.375 | 1.014 | NA | NA | too_few_samples; CLDN4 vs MPR companion on the same n |
| NMPR vs MPR: epithelial CLDN4, post complete-case | 8 | 4 | 1.509 | 1.505 | +0.003 | 0.933 | CLDN4 vs MPR companion |

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| epithelium CLDN4 vs DPT (post complete-case) | 12 | -0.049 | 0.88 |
| leftover-epi CLDN4 vs DPT (post ≥10 leftover) | 12 | 0.028 | 0.931 |
| A3-malignant CLDN4 vs DPT including 1–9 cell patients | 10 | -0.152 | 0.676 |
| epithelium CLDN4 vs barrier/keratin (post) | 12 | 0.203 | 0.527 |
| epithelium CLDN4 vs TACSTD2 (companion, not a gate) | 12 | 0.392 | 0.208 |

## Extra figure — CLDN4-high vs CLDN4-low (sample-paired, epithelium on graph)

Emitted: **True**. Rule: always emit (requested extra figure). Paired tertile n=11.

| Contrast | n | W | Δmed (high−low) | p |
| --- | ---: | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 11 | 5.0 | 0.127 | 0.00977 |
| AT2 score high vs low | 11 | 23.0 | -0.021 | 0.413 |
| DPT high vs low | 11 | 10.0 | 0.030 | 0.042 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- A3-malignant-like is a marker rule, **not CNV**. Residual unmarked epithelium can leak in.
- Empty MPR residuals (normal-lung program) are reported, not patched by dual-high.
- GSE207422-only T/NK was flat in the given slice and is not re-tested here.
- No TACSTD2∩CLDN4 both-high gate. TACSTD2 is a companion Spearman only.
- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Do not write n=12 as the malignant-like test n when empty MPR were dropped.
- Mixed NSCLC (LUAD + LUSC). Not LUAD-only.

## Outputs

- `figures/fig_trajectory_cldn4.png` — PAGA + UMAP CLDN4 + patient CLDN4 vs DPT
- `figures/fig_dpt_vs_mpr.png` — patient DPT vs MPR (malignant-like and epithelium)
- `figures/fig_honest_n.png` — per-patient malignant vs leftover counts
- `figures/fig_extra_cldn4_tertile.png` — within-sample CLDN4-high vs low
- `tables/patient_means.tsv` — per-patient means (do not sum cells across patients)
- `tables/patient_cldn4_vs_dpt.tsv` — patient-level CLDN4 vs DPT
- `tables/patient_dpt_vs_mpr.tsv` — patient-level DPT vs MPR
- `tables/summary.json`

## Reproduce

```bash
pip install -r methods/gse207422_traj_cldn4_mpr/requirements.txt
python3 methods/gse207422_traj_cldn4_mpr/scripts/download.py
python3 methods/gse207422_traj_cldn4_mpr/scripts/extract_epithelium.py
python3 methods/gse207422_traj_cldn4_mpr/scripts/analyze_paga_dpt.py
```


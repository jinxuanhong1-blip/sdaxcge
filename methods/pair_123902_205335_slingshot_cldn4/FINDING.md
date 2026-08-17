# Finding — pair GSE123902+GSE205335 epithelium, CLDN4-only Palantir/PAGA

ADDITIVE. **CLDN4 only.** Pair that **differs** in PR #459 (GSE123902+GSE205335 %pos vs T/NK, n=35, ρ=−0.522, Q4 vs Q1 r=−0.802). That T/NK Spearman is **taken as given** and is not re-audited. PR #473 malignant IFN DE (IFN logFC −1.05, n=9/9) is a **different** contrast and is not re-run. No TACSTD2∩CLDN4 dual-high gate. GSE148071 is not added.

Primary clock: **Palantir (Setty et al. 2019); PAGA geometry; DPT companion**. Inferential unit = **GSE123902 donor + GSE205335 patient** (tumor units only). Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. IFN score = Hallmark IFNα ∪ IFNγ (same family as PR #473). Root is GSE123902 NORMAL AT2-like, **never CLDN4-high**.

**What holds (n=35 tumor units).** CLDN4 vs Palantir PT is n=35, ρ=-0.333, p=0.0504. IFN vs Palantir PT is n=35, ρ=0.255, p=0.139. barrier/keratin (no CLDN4) vs Palantir PT is n=35, ρ=-0.376, p=0.026. CLDN4 vs IFN is n=35, ρ=-0.139, p=0.426. CLDN4 vs barrier is n=35, ρ=0.375, p=0.0264. Within-unit CLDN4-high vs low IFN: n=32, W=64.0, Δmed=0.062, p=6.69e-05; barrier: n=32, W=0.0, Δmed=0.564, p=4.66e-10. **What this is not.** Given T/NK n=35 and IFN DE n=18 are different contrasts. Pooled PT mixes cohorts and is not a within-tumor progression test. Within-unit CLDN4-high cells are slightly *more* IFN (Δmed +0.062); that is **not** a copy of the PR #473 malignant pseudobulk IFN logFC −1.05.

## Verdict

Sample-level CLDN4 vs Palantir PT: n=35, ρ=-0.333, p=0.0504. IFN vs Palantir PT: n=35, ρ=0.255, p=0.139. barrier/keratin (CLDN4 excluded) vs Palantir PT: n=35, ρ=-0.376, p=0.026. CLDN4 vs IFN: n=35, ρ=-0.139, p=0.426. CLDN4 vs barrier/keratin (CLDN4 excluded): n=35, ρ=0.375, p=0.0264. GSE123902-only CLDN4 vs PT: n=13, ρ=-0.484, p=0.0941. GSE205335-only CLDN4 vs PT: n=22, ρ=-0.128, p=0.57. Paired CLDN4-high vs low barrier/keratin: n=32, W=0.0, Δmed=0.564, p=4.66e-10. Paired CLDN4-high vs low IFN: n=32, W=64.0, Δmed=0.062, p=6.69e-05. PAGA has 1 component(s) at connectivity>0 among 28 Leiden vertices. Root is GSE123902 NORMAL AT2-like and is not CLDN4-high. Given PR #459 T/NK and PR #473 IFN DE were not re-audited. No both-high gate. GSE148071 not added.

## Honest n

- Analysis cells after QC (capped ≤350/unit): **n_cells = 10334** (GSE123902 3403, GSE205335 6931).
- Tumor units (GSE123902 donor + GSE205335 patient): **n_tumor_units = 35** (GSE123902 13, GSE205335 22).
- Tumor units with ≥10 epithelial cells used for Spearman: **n = 35**.
- GSE123902 NORMAL cells in the object (root material, not primary units): **755** (4 libraries).
- Author / marker subtypes (cells): {'Malignant cells': 6228, 'marker_malignant': 2648, 'NORMAL_AT2like': 755, 'Non-malignant cells': 703}.
- CLDN4 tertile cells: {'low': 4121, 'high': 3445, 'mid': 2768}.
- Units with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 32**.
- Genes absent from locked sets: {'AT2': ['SFTPC', 'SFTPA1'], 'AT1': [], 'club': ['SCGB1A1', 'SCGB3A2', 'SCGB3A1'], 'basal': ['KRT5', 'TP63', 'NGFR'], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': ['CEACAM5', 'CEACAM6'], 'IFN': ['BATF2', 'CCL7', 'CXCL11', 'CXCL9', 'GBP6', 'IDO1', 'KLRK1', 'MARCHF1', 'RNF31', 'SELP', 'SLAMF7', 'TENT5A', 'WARS1'], 'single_genes': ['SFTPC', 'SCGB1A1', 'KRT5']}.
- Given PR #459 T/NK n=35 is **not** this n.
- Given PR #473 IFN DE n=18 (9/9) is **not** this n.
- GSE148071 not used. Dual-high not used. T/NK infiltrate not re-scored.
- Slingshot: available=False; Rscript not on PATH.
- Palantir: available=True; version=1.4.5. Palantir warned that some cells were unreachable at k=30; the clock is still Palantir PT on the reachable graph, not a fake ordering.
- SFTPC / SFTPA1 / club genes are absent from the shared GSE123902∩GSE205335 universe. AT2 score uses the remaining locked AT2 genes (SFTPB, NAPSA, LAMP3, ABCA3). The SFTPC-vs-PT control is therefore n=0.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- Root: GSE123902 NORMAL AT2-like (median AT2; CLDN4-high excluded) (root cell index 1178, unit GSE123902:LX675:NORMAL).
- PAGA components at connectivity>0: **1** among 28 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN genes: Hallmark IFNα ∪ IFNγ, n_present=211 / 224.

## Primary (sample-level Spearman, BH inside this list)

| Contrast | n_units | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Palantir PT | 35 | -0.333 | 0.0504 | 0.113 |
| IFN vs Palantir PT | 35 | 0.255 | 0.139 | 0.249 |
| barrier/keratin (no CLDN4) vs Palantir PT | 35 | -0.376 | 0.026 | 0.0792 |
| CLDN4 vs IFN | 35 | -0.139 | 0.426 | 0.639 |
| CLDN4 vs barrier/keratin (no CLDN4) | 35 | 0.375 | 0.0264 | 0.0792 |
| CLDN4 vs AT2 score | 35 | 0.054 | 0.757 | 0.851 |
| CLDN4 vs TACSTD2 (comparator) | 35 | 0.426 | 0.0108 | 0.0792 |
| SFTPC vs Palantir PT (control) | 0 | NA | NA | NA |
| AT2 score vs Palantir PT (control) | 35 | -0.092 | 0.6 | 0.771 |

## Sensitivity (not in the BH family)

| Contrast | n_units | ρ | p |
| --- | ---: | ---: | ---: |
| GSE123902-only CLDN4 vs Palantir PT | 13 | -0.484 | 0.0941 |
| GSE205335-only CLDN4 vs Palantir PT | 22 | -0.128 | 0.57 |
| GSE123902-only IFN vs Palantir PT | 13 | 0.489 | 0.0899 |
| GSE205335-only IFN vs Palantir PT | 22 | 0.137 | 0.543 |
| GSE123902-only barrier vs Palantir PT | 13 | 0.055 | 0.859 |
| GSE205335-only barrier vs Palantir PT | 22 | -0.525 | 0.0122 |
| GSE123902-only CLDN4 vs IFN | 13 | -0.143 | 0.642 |
| GSE205335-only CLDN4 vs IFN | 22 | -0.022 | 0.923 |
| include-NORMAL-units CLDN4 vs Palantir PT | 39 | -0.340 | 0.034 |
| GSE205335 ADC+SQ CLDN4 vs Palantir PT | 17 | -0.277 | 0.282 |
| tumor CLDN4 vs DPT (companion) | 35 | -0.115 | 0.512 |
| tumor IFN vs DPT (companion) | 35 | 0.082 | 0.638 |

## Extra figure — CLDN4 / IFN / barrier along Palantir PT

Emitted: **True**. Sample-level CLDN4 vs barrier ρ=0.375, p=0.0264.

| Paired contrast (high − low) | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 32 | 0.564 | 4.66e-10 |
| IFN high vs low | 32 | 0.062 | 6.69e-05 |
| AT2 high vs low | 32 | 0.165 | 5.99e-06 |
| Palantir PT high vs low | 32 | -0.025 | 0.000837 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a redo of PR #459 T/NK infiltrate or PR #473 malignant IFN DE.
- The pooled PT Spearman mixes two cohorts and is **not** a within-tumor progression test.
- Marker-malignant on GSE123902 is EPCAM|KRT>0 and PTPRC==0 — **not CNV**.
- GSE205335 is an ICI biopsy/effusion cohort, but this analysis is **not** an ICI / MPR / RECIST test.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “NORMAL AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Palantir/DPT is an ordering, not a clock of real time.
- GSE148071 was not added.

## Outputs

- `results/tables/lineage_leiden.tsv` — **done criterion (lineage)**
- `results/tables/sample_level_spearman.tsv` — **done criterion (sample-level)**
- `results/tables/sample_means.tsv`
- `results/tables/paga_connectivities.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_along_pt_programs.png`
- `results/figures/fig_extra_sample_along_pt.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/figures/fig_paga_ifn.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/pair_123902_205335_slingshot_cldn4/requirements.txt
python3 methods/pair_123902_205335_slingshot_cldn4/scripts/download.py \
  --outdir /tmp/pair_123902_205335_traj
python3 methods/pair_123902_205335_slingshot_cldn4/scripts/extract_epithelium.py \
  --data /tmp/pair_123902_205335_traj \
  --out /tmp/pair_123902_205335_traj/epithelium.h5ad
python3 methods/pair_123902_205335_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/pair_123902_205335_traj/epithelium.h5ad \
  --outdir methods/pair_123902_205335_slingshot_cldn4/results \
  --finding methods/pair_123902_205335_slingshot_cldn4/FINDING.md
```


# Finding — GSE189357 tumor epithelium, CLDN4-only REAL Slingshot/PAGA

ADDITIVE. **CLDN4 only.** Zhu et al., *Exp Mol Med* 2022 (DOI 10.1038/s12276-022-00896-9), GEO **GSE189357**: nine treatment-naïve resected LUAD lesions (TD1–TD9; AIS=3, MIA=3, IAC=3). This folder does **not** redo winning-pair Slingshot (GSE131907+GSE205335) or GSE131907-only PAGA. No TACSTD2∩CLDN4 dual-high gate. Spatial GSE189487 is not used.

Primary clock: **REAL Slingshot** (Street et al. 2018; Bioconductor `slingshot`). PAGA is geometry only. Inferential unit = **patient** (one 10x sample each). **n may be 9 — say so.** Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4**. IFN is Hallmark IFNα ∪ IFNγ. Root is the highest-AT2 Leiden cluster that is **not** the CLDN4-highest cluster; the root cell is never CLDN4-high.

**What holds (n=9 patients; n may be 9).** CLDN4 vs Slingshot PT n=9, ρ=-0.400, p=0.286. CLDN4 vs barrier/keratin (CLDN4 excluded) n=9, ρ=0.517, p=0.154. CLDN4 vs IFN n=9, ρ=0.517, p=0.154. Within-patient CLDN4-high vs low barrier n=9, W=0.0, Δmed=0.510, p=0.00391; IFN n=9, W=8.0, Δmed=0.034, p=0.0977. **Limits.** Catalog n=9 (3 AIS / 3 MIA / 3 IAC). Harmony on patient removes between-stage mean shifts. No nLung AT2. Marker-malignant ≠ CNV.

## Verdict

Patient-level CLDN4 vs Slingshot PT: n=9, ρ=-0.400, p=0.286 (n may be 9). CLDN4 vs barrier/keratin (CLDN4 excluded): n=9, ρ=0.517, p=0.154. CLDN4 vs IFN (Hallmark α∪γ): n=9, ρ=0.517, p=0.154. barrier vs Slingshot PT: n=9, ρ=-0.850, p=0.0037. IFN vs Slingshot PT: n=9, ρ=0.400, p=0.286. AT2 vs Slingshot PT (control): n=9, ρ=-0.433, p=0.244. Paired CLDN4-high vs low barrier/keratin: n=9, W=0.0, Δmed=0.510, p=0.00391. Paired CLDN4-high vs low IFN: n=9, W=8.0, Δmed=0.034, p=0.0977. Slingshot lineages=10; start=11; excluded CLDN4-high cluster=4. PAGA has 1 component(s) at connectivity>0 among 18 Leiden vertices. REAL Slingshot was run. Patient is the unit. n may be 9. Not a TACSTD2 redo. No both-high gate.

## Honest n

- GEO catalog: **n_patients = 9** (TD1–TD9). This is the catalog n and, after the ≥10-cell rule, the Spearman n may still be 9. **Say so.**
- Stage split: AIS=3, MIA=3, IAC=3 (n=3/stage — not a stage test).
- Marker-malignant cells after QC: **n_cells = 14168** (gate `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`).
- Patients with ≥10 epithelial cells used for Spearman: **n = 9**.
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 9**.
- CLDN4 tertile cells: low 6023, mid 3422, high 4723.
- Cells per patient: {'TD1': 2153, 'TD2': 1187, 'TD3': 663, 'TD4': 490, 'TD5': 1729, 'TD6': 1202, 'TD7': 849, 'TD8': 1018, 'TD9': 4877}.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': ['CCL7', 'KLRK1', 'MARCHF1', 'WARS1'], 'IFN_core': [], 'single_genes': ['PTPRC']}.
- Slingshot: available=True; version=2.10.0; n_lineages=10.
- Harmony: harmonypy on PCA, batch=patient. Each patient is one stage (AIS/MIA/IAC), so Harmony removes between-patient/stage mean shifts; the remaining axis is within-epithelium programs (n=9).
- No nLung / uninvolved AT2 exists in this accession. Root is AT2-high tumor epithelium, not CLDN4-high.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: Harmony on **patient** (n=9; each patient is one stage).
- Slingshot start cluster: **11** (excluded CLDN4-high cluster 4).
- Root cell: patient TD1 stage IAC tertile mid (never high).
- PAGA components at connectivity>0: **1** among 18 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN: Hallmark IFNα ∪ IFNγ (224 genes; CLDN4 not in the set).

## Slingshot lineages (done criterion)

| lineage | path | n_clusters | start | end | n_cells | mean PT | mean CLDN4 | mean barrier | mean IFN |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| Lineage1 | 11->6->0->12->4->3->1 | 7 | 11 | 1 | 7467 | 18.316 | 1.586 | 1.180 | 0.270 |
| Lineage2 | 11->6->0->12->4->3->8 | 7 | 11 | 8 | 7128 | 17.380 | 1.560 | 1.152 | 0.264 |
| Lineage3 | 11->6->0->12->10->16->15 | 7 | 11 | 15 | 6596 | 20.182 | 1.089 | 0.964 | 0.286 |
| Lineage4 | 11->6->0->12->10->5 | 6 | 11 | 5 | 7585 | 24.291 | 0.970 | 0.924 | 0.273 |
| Lineage5 | 11->6->0->12->10->7 | 6 | 11 | 7 | 6802 | 22.704 | 1.257 | 1.001 | 0.277 |
| Lineage6 | 11->6->0->12->10->9 | 6 | 11 | 9 | 7391 | 20.195 | 1.042 | 0.931 | 0.278 |
| Lineage7 | 11->6->0->12->10->13 | 6 | 11 | 13 | 6915 | 20.581 | 1.124 | 0.973 | 0.278 |
| Lineage8 | 11->6->0->12->10->14 | 6 | 11 | 14 | 6803 | 19.739 | 1.114 | 0.974 | 0.271 |
| Lineage9 | 11->6->0->12->10->17 | 6 | 11 | 17 | 6474 | 19.262 | 1.120 | 0.980 | 0.278 |
| Lineage10 | 11->6->0->2 | 4 | 11 | 2 | 3759 | 20.401 | 1.403 | 1.135 | 0.261 |

Machine table: `results/tables/slingshot_lineages.tsv`. CLDN4-ward lineage used for along-pseudotime plots: **Lineage1**.

## Primary (patient-level Spearman, BH inside this list)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 9 | -0.400 | 0.286 | 0.409 |
| CLDN4 vs Slingshot CLDN4-ward PT | 9 | 0.233 | 0.546 | 0.606 |
| CLDN4 vs AT2 score | 9 | 0.167 | 0.668 | 0.668 |
| CLDN4 vs barrier/keratin (no CLDN4) | 9 | 0.517 | 0.154 | 0.386 |
| CLDN4 vs IFN (Hallmark α∪γ) | 9 | 0.517 | 0.154 | 0.386 |
| barrier/keratin vs Slingshot PT | 9 | -0.850 | 0.0037 | 0.037 |
| IFN vs Slingshot PT | 9 | 0.400 | 0.286 | 0.409 |
| CLDN4 vs TACSTD2 (comparator) | 9 | 0.800 | 0.00963 | 0.0481 |
| SFTPC vs Slingshot PT (control) | 9 | 0.367 | 0.332 | 0.415 |
| AT2 score vs Slingshot PT (control) | 9 | -0.433 | 0.244 | 0.409 |

n may be 9. A significant p at n=9 is a small-n result, not a cohort.

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| drop IAC CLDN4 vs Slingshot PT | 6 | 0.029 | 0.957 |
| drop AIS CLDN4 vs Slingshot PT | 6 | -0.771 | 0.0724 |
| CLDN4 vs DPT (companion) | 9 | -0.533 | 0.139 |
| IFN core vs Slingshot PT | 9 | 0.800 | 0.00963 |
| CLDN4 vs IFN core | 9 | -0.100 | 0.798 |
| CLDN4 vs malignant-like | 9 | 0.517 | 0.154 |
| IFN vs Slingshot CLDN4-ward PT | 9 | -0.200 | 0.606 |
| barrier vs Slingshot CLDN4-ward PT | 9 | 0.700 | 0.0358 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired) and programs along Slingshot

Emitted: **True**. Rule: always emit extra figures on this accession. Observed Spearman(CLDN4, barrier) n=9, ρ=0.517, p=0.154.

| Paired contrast (high − low) | n_patients | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 9 | 0.510 | 0.00391 |
| IFN high vs low | 9 | 0.034 | 0.0977 |
| AT2 high vs low | 9 | 0.658 | 0.00391 |
| Slingshot PT high vs low | 9 | -5.500 | 0.0391 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- **n may be 9.** This is not a 65-unit winning-pair result and not a stage-powered test (n=3/stage).
- Harmony on patient removes between-stage mean shifts. Do not read Slingshot as AIS→IAC invasion time.
- No uninvolved nLung AT2 exists here. The root is AT2-high **tumor** epithelium.
- Malignant is a marker gate, **not CNV** (inferCNV was not run).
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into IAC because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot/DPT is an ordering, not a clock.

## Outputs

- `results/tables/slingshot_lineages.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_along_pseudotime.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
bash methods/gse189357_slingshot_real_cldn4/scripts/install_tools.sh
python3 methods/gse189357_slingshot_real_cldn4/scripts/download.py \
  --out /tmp/gse189357
python3 methods/gse189357_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --tar /tmp/gse189357/GSE189357_RAW.tar \
  --out /tmp/gse189357/epithelium.h5ad
python3 methods/gse189357_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse189357/epithelium.h5ad \
  --outdir methods/gse189357_slingshot_real_cldn4/results \
  --finding methods/gse189357_slingshot_real_cldn4/FINDING.md
```


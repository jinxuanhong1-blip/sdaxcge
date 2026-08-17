# Finding — ICI pair GSE205335 + GSE207422, CLDN4-only PAGA/DPT vs outcome

ADDITIVE. **CLDN4 only.** The two public ICI scRNA sets that carry outcome labels: **GSE205335** (Ahn/Lee eLife 2024; palliative ICI; RECIST) + **GSE207422** (Hu *Genome Med* 2023; neoadjuvant PD-1 + chemo; MPR). **Not a 7-pool. No GSE148071. No dual-high / TACSTD2 gate.** Patient is the unit. Honest n is not cell count.

PAGA + leftover-AT2-rooted diffusion pseudotime are run **per dataset**. Platforms differ (GSE205335 mixed 3′/5′ 10x biopsies/effusions; GSE207422 10x post-surgery). Harmony was **not** forced. Raw DPT is not commensurate across datasets; the stacked DPT tests use within-dataset rank / z.

Root is leftover (non-malignant / A3-normal-lung) epithelium in the giant PAGA component, high AT2, **not CLDN4-high.** Barrier/keratin excludes CLDN4.

## Verdict

GSE205335 patient-level CLDN4 vs leftover-AT2-rooted DPT: n=22, ρ=0.069, p=0.759. GSE205335 DPT RECIST PR vs PD: PR n=6 med=0.316 vs PD n=7 med=0.084, p=0.101. GSE207422 patient-level CLDN4 vs leftover-AT2-rooted DPT: n=12, ρ=-0.517, p=0.0849. GSE207422 DPT MPR vs NMPR: MPR n=4 med=0.491 vs NMPR n=8 med=0.524, p=0.808. Stacked CLDN4 vs within-dataset DPT rank: n=34, ρ=-0.213, p=0.225. Stacked DPT z, PR+MPR vs PD+NMPR: benefit n=10 med=0.249 vs no-benefit n=15 med=-0.299, p=0.233. DL Fisher-z of the two per-dataset CLDN4–DPT ρ: k=2, N=34, ρ=-0.203, p=0.517. Harmony was not run. Not a 7-pool. No GSE148071. No dual-high. Root is not CLDN4-high.

## Honest n

### GSE205335 (RECIST)

- GEO patients / samples: **26 / 33**. Catalog, not the test n.
- Tumor-tissue epithelium in the object: **n_cells = 7882** in **22** patients (normal LN/lung/brain dropped; 4 normal-only patients never enter).
- Patients with ≥10 epithelial cells (program Spearman): **n = 22**.
- Patients with finite leftover-rooted DPT: **n = 22**.
- Patients with ≥10 author malignant cells: **n = 22**.
- RECIST on eligible patients: {'PD': 7, 'NE': 6, 'PR': 6, 'SD': 3}.
- Primary outcome contrast is **PR vs PD** (SD and NE excluded). R vs NR (PR vs SD+PD) is sensitivity only.
- Lineage on the object: {'malignant': 6524, 'leftover': 1358}. Author subtypes: {'Malignant cells': 6524, 'Non-malignant cells': 1358}.
- Histology (patients): {'ADC': 14, 'NUT': 1, 'SCLC': 4, 'SQ': 3}. SCLC/NUT stay in the all-comer row.
- DPT root: GSE205335 leftover in giant PAGA component, CLDN4≤leftover tertile-2 (cut=0.000), max AT2 (80th percentile was 0) (patient P0031; CLDN4=0.000; AT2=4.463).
- PAGA components at connectivity>0: **1** among 26 Leiden vertices.
- MPR is not a GEO field here. RECIST is not substituted for MPR.

### GSE207422 (MPR)

- GEO samples: **15** (3 pre-treatment TN + 12 post-surgery). Catalog, not the test n.
- Post-treatment epithelium in the object: **n_cells = 2994** in **12** patients. The 3 pre-biopsies are out of the graph.
- Patients with ≥10 epithelial cells (program Spearman): **n = 12**.
- Patients with finite leftover-rooted DPT: **n = 12**. Root was placed in the giant PAGA component so DPT is defined for every eligible patient; PAGA still has 2 component(s).
- Patients with ≥10 A3-malignant-like cells: **n = 7**. A3 empties some MPR residuals (normal-lung program); that is reported, not patched.
- Pathologic response on eligible patients: {'NMPR': 8, 'MPR': 4} (pCR P06 = MPR).
- Lineage on the object: {'malignant': 1116, 'leftover': 1878}.
- DPT root: GSE207422 leftover in giant PAGA component, CLDN4≤leftover tertile-2 (cut=1.982), AT2 ~80th percentile (1.618) (patient P14; CLDN4=0.000; AT2=1.618).
- PAGA components at connectivity>0: **2** among 15 Leiden vertices.
- Author CopyKAT / DRMref barcodes are not on GEO. A3 marker malignant-like is used.

### Stacked (not a joint embedding)

- Stacked patient rows: **N = 34** (GSE205335 22 + GSE207422 12).
- Stacked CLDN4 vs within-dataset DPT rank: **n = 34**.
- Stacked benefit vs no-benefit (PR+MPR vs PD+NMPR) on DPT z: **n = 10 vs 15**. SD/NE/TN are out of the binary stack.
- Raw DPT is **not** pooled. Harmony was not run.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: none (per-dataset PCA/neighbors). Harmony not forced.
- DPT root: leftover epithelium in the giant PAGA component, not CLDN4-high, high AT2.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- GSE205335 malignant = author `Malignant cells`. GSE207422 malignant-like = epithelial AND zero UMI of SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3.

## Per-dataset — CLDN4 vs DPT (patient Spearman, BH inside each dataset list)

| Dataset | Contrast | n | ρ | p | q |
| --- | --- | ---: | ---: | ---: | ---: |
| GSE205335 | CLDN4 vs DPT | 22 | 0.069 | 0.759 | 0.759 |
| GSE205335 | CLDN4 vs AT2 score | 22 | 0.161 | 0.474 | 0.61 |
| GSE205335 | CLDN4 vs club score | 22 | 0.363 | 0.0968 | 0.145 |
| GSE205335 | CLDN4 vs basal score | 22 | 0.078 | 0.728 | 0.759 |
| GSE205335 | CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.390 | 0.0726 | 0.131 |
| GSE205335 | CLDN4 vs malignant-like score | 22 | 0.476 | 0.0251 | 0.113 |
| GSE205335 | CLDN4 vs TACSTD2 (comparator) | 22 | 0.395 | 0.0691 | 0.131 |
| GSE205335 | SFTPC vs DPT (control) | 22 | -0.390 | 0.0725 | 0.131 |
| GSE205335 | AT2 score vs DPT (control) | 22 | -0.727 | 0.000128 | 0.00115 |
| GSE207422 | CLDN4 vs DPT | 12 | -0.517 | 0.0849 | 0.255 |
| GSE207422 | CLDN4 vs AT2 score | 12 | 0.559 | 0.0586 | 0.255 |
| GSE207422 | CLDN4 vs club score | 12 | 0.273 | 0.391 | 0.678 |
| GSE207422 | CLDN4 vs basal score | 12 | -0.224 | 0.484 | 0.678 |
| GSE207422 | CLDN4 vs barrier/keratin (no CLDN4) | 12 | 0.203 | 0.527 | 0.678 |
| GSE207422 | CLDN4 vs malignant-like score | 12 | 0.028 | 0.931 | 1 |
| GSE207422 | CLDN4 vs TACSTD2 (comparator) | 12 | 0.399 | 0.199 | 0.449 |
| GSE207422 | SFTPC vs DPT (control) | 0 | NA | NA | NA |
| GSE207422 | AT2 score vs DPT (control) | 12 | -0.776 | 0.00299 | 0.0269 |

## Per-dataset — DPT vs ICI benefit (patient MWU)

| Dataset | Contrast | n_a | n_b | med_a | med_b | p |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| GSE205335 | DPT PR vs PD | 6 | 7 | 0.316 | 0.084 | 0.101 |
| GSE205335 | CLDN4 PR vs PD | 6 | 7 | 1.281 | 0.986 | 0.628 |
| GSE205335 | DPT RECIST R(PR) vs NR(SD+PD) [sensitivity] | 6 | 10 | 0.316 | 0.107 | 0.118 |
| GSE205335 | DPT PR vs PD ADC+SQ [sensitivity] | 4 | 6 | 0.143 | 0.065 | 0.257 |
| GSE207422 | DPT MPR vs NMPR | 4 | 8 | 0.491 | 0.524 | 0.808 |
| GSE207422 | CLDN4 MPR vs NMPR | 4 | 8 | 1.538 | 1.616 | 0.933 |

## Stacked patient table — CLDN4 vs DPT

DPT is rank- and z-transformed **within dataset** before stacking. Fisher-z DL is descriptive (k=2).

| Contrast | k | N | ρ | p | I² | note |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| GSE205335 CLDN4 vs DPT | 1 | 22 | 0.069 | 0.759 | — | per-dataset raw DPT |
| GSE207422 CLDN4 vs DPT | 1 | 12 | -0.517 | 0.0849 | — | per-dataset raw DPT |
| stacked CLDN4 vs DPT rank (within dataset) | 2 | 34 | -0.213 | 0.225 | — | pooled patients; DPT ranked within dataset |
| stacked CLDN4 vs DPT z (within dataset) | 2 | 34 | -0.092 | 0.607 | — | pooled patients; DPT z within dataset |
| DL Fisher-z of per-dataset CLDN4 vs DPT | 2 | 34 | -0.203 | 0.517 | 60% | DL Fisher-z |

## Stacked — DPT vs ICI benefit (PR+MPR vs PD+NMPR)

| Contrast | n_benefit | n_no_benefit | med_z_benefit | med_z_no | p |
| --- | ---: | ---: | ---: | ---: | ---: |
| stacked DPT z, benefit (PR+MPR) vs no-benefit (PD+NMPR) | 10 | 15 | 0.249 | -0.299 | 0.233 |
| stacked DPT rank, benefit vs no-benefit | 10 | 15 | 8.500 | 8.000 | 0.232 |
| stacked CLDN4, benefit (PR+MPR) vs no-benefit (PD+NMPR) | 10 | 15 | 1.497 | 1.515 | 1 |
| stacked DPT z, benefit vs no-benefit, drop GSE205335 SCLC/NUT | 8 | 14 | 0.192 | -0.453 | 0.474 |

## Sensitivity (not in the BH family)

| Dataset | Contrast | n | ρ | p |
| --- | --- | ---: | ---: | ---: |
| GSE205335 | malignant-only CLDN4 vs DPT | 22 | -0.037 | 0.871 |
| GSE205335 | malignant-only CLDN4 vs AT2 | 22 | 0.226 | 0.311 |
| GSE205335 | ADC+SQ CLDN4 vs DPT | 17 | -0.020 | 0.94 |
| GSE207422 | malignant-only CLDN4 vs DPT | 7 | -0.357 | 0.432 |
| GSE207422 | malignant-only CLDN4 vs AT2 | 7 | 0.571 | 0.18 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- This is not a 7-pool, not GSE148071, not a TACSTD2 redo, not dual-high.
- Per-dataset DPT clocks are not a shared latent time. Do not interpret stacked raw DPT.
- Do not write “AT2 differentiates into ICI-resistant NSCLC because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- GSE205335 RECIST PR vs PD is not MPR. GSE207422 MPR is not RECIST.
- A3 malignant-like is not CopyKAT. Author malignant on GSE205335 is not CNV re-called here.
- Harmony was not run. A null stacked test is not evidence that a joint embedding would be null.
- DPT is an ordering, not a clock. Slingshot was not required.

## Outputs (done criterion)

- `results/tables/gse205335_patient_means.tsv` — per-dataset patient table
- `results/tables/gse207422_patient_means.tsv` — per-dataset patient table
- `results/tables/per_dataset_cldn4_vs_dpt.tsv`
- `results/tables/per_dataset_dpt_vs_benefit.tsv`
- `results/tables/stacked_patient_table.tsv`
- `results/tables/stacked_cldn4_vs_dpt.tsv`
- `results/tables/stacked_dpt_vs_benefit.tsv`
- `results/tables/honest_n.tsv`

## Reproduce

```bash
pip install -r methods/ici_pair_205335_207422_traj_cldn4/requirements.txt
python3 methods/ici_pair_205335_207422_traj_cldn4/scripts/download.py \
  --out /tmp/ici_pair_205335_207422
python3 methods/ici_pair_205335_207422_traj_cldn4/scripts/extract.py \
  --data /tmp/ici_pair_205335_207422 \
  --out /tmp/ici_pair_205335_207422/extracted
python3 methods/ici_pair_205335_207422_traj_cldn4/scripts/analyze.py \
  --extracted /tmp/ici_pair_205335_207422/extracted \
  --outdir methods/ici_pair_205335_207422_traj_cldn4/results \
  --finding methods/ici_pair_205335_207422_traj_cldn4/FINDING.md
```


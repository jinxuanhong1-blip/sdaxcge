# Finding — GSE205335 epithelium trajectory, CLDN4 and RECIST

ADDITIVE. **CLDN4 only.** GSE205335 advanced NSCLC ICI (Ahn / Lee, *eLife* 2024), public processed UMI. **No GSE148071. No dual-high.** Does **not** re-audit PR #320 T/NK ρ. Patient is the unit.

Primary clock: **documented AT2-like / low-CLDN4 diffusion pseudotime (scanpy DPT; Slingshot R missing)**. Root is a biological AT2-like / low-CLDN4 state, **never CLDN4-high**. Barrier/keratin score **excludes CLDN4**. RECIST test is **PR vs PD/SD** (NE excluded). Cell-level ρ is descriptive.

**What holds (patient n=22; RECIST n=16, PR 6 vs PD/SD 10).** CLDN4 vs DPT: n=22, ρ=0.098, p=0.665. DPT vs RECIST (PR vs PD/SD): n=16 (PR 6 vs PD/SD 10), Δmed=+0.172, r=+0.467, p=0.147. CLDN4 vs barrier/keratin (CLDN4 excluded): n=22, ρ=0.233, p=0.296. **What does not hold.** Patient-level CLDN4 vs DPT is null, so this is **not** a lineage / differentiation claim. The RECIST split is still the deliverable. DPT does not significantly separate PR vs PD/SD at this n; that is reported, not hidden. Root used: author Non-malignant cells, CLDN4 ≤ object median, median AT2 among AT2>0 (n_low=684, n_AT2>0=170). Do not re-read this as a T/NK result (PR #320 is not re-audited).

## Verdict

Patient-level CLDN4 vs AT2-like/low-CLDN4 DPT: n=22, ρ=0.098, p=0.665. CLDN4 vs DPT is null at the patient — this is not lineage proof; RECIST is still reported. CLDN4 vs AT2 score: n=22, ρ=0.133, p=0.556. CLDN4 vs barrier/keratin (CLDN4 excluded): n=22, ρ=0.233, p=0.296. CLDN4 vs malignant-like: n=22, ρ=0.482, p=0.0232. DPT vs RECIST (PR vs PD/SD): n=16 (PR 6 vs PD/SD 10), Δmed=+0.172, r=+0.467, p=0.147. CLDN4 vs RECIST (PR vs PD/SD): n=16 (PR 6 vs PD/SD 10), Δmed=+0.447, r=+0.233, p=0.492. Paired CLDN4-high vs low barrier/keratin: n=22, W=0.0, Δmed=+0.559, p=4.77e-07. Paired CLDN4-high vs low DPT: n=22, W=123.0, Δmed=-0.002, p=0.924. PAGA has 1 component(s) at connectivity>0 among 25 Leiden vertices. Slingshot R was not run unless available; DPT is the documented clock. Not a TACSTD2 redo. No both-high gate. GSE148071 not used. PR #320 T/NK not re-tested.

## Honest n

- GEO catalog epithelium on non-normal tissues: **30448** cells / **22** patients.
- Analysis cells after QC (capped ≤500/patient): **n_cells = 9500**.
- Patients on the object: **n_patients = 22**.
- Patients with ≥10 epithelial cells (Spearman): **n = 22**.
- RECIST-evaluable (PR vs PD/SD, NE out): **n = 16** (PR 6 vs PD/SD 10; NE 6).
- Author malignant / non-malignant cells on the object: 8640 / 860.
- Patients with ≥10 author-malignant cells: **n = 22**.
- Author subtypes (cells): {'Malignant cells': 8640, 'Non-malignant cells': 860}.
- Histology (patients): {'ADC': 14, 'SCLC': 4, 'SQ': 3, 'NUT': 1}.
- RECIST (patients): {'PD': 7, 'NE': 6, 'PR': 6, 'SD': 3}.
- CLDN4 tertile cells: {'low': 3605, 'high': 3167, 'mid': 2728}.
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 22**.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': []}.
- Slingshot: available=False; Rscript not on PATH.
- GSE148071 not used. GSE131907 not used. Dual-high not used. PR #320 T/NK not re-tested.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- No Harmony (single cohort; Harmony on patient would erase the RECIST contrast).
- DPT root: author Non-malignant cells, CLDN4 ≤ object median, median AT2 among AT2>0 (n_low=684, n_AT2>0=170) (root cell index 2991, patient P1076, subtype Non-malignant cells, CLDN4=0.586, AT2=0.293).
- PAGA components at connectivity>0: **1** among 25 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- RECIST arms: PR vs PD+SD. NE is excluded from the RECIST table, kept in CLDN4–DPT.

## Primary — sample-level CLDN4 vs DPT (Spearman, BH inside this list)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs DPT | 22 | 0.098 | 0.665 | 0.749 |
| CLDN4 vs AT2 score | 22 | 0.133 | 0.556 | 0.715 |
| CLDN4 vs club score | 22 | 0.252 | 0.257 | 0.47 |
| CLDN4 vs basal score | 22 | 0.022 | 0.923 | 0.923 |
| CLDN4 vs barrier/keratin (no CLDN4) | 22 | 0.233 | 0.296 | 0.47 |
| CLDN4 vs malignant-like score | 22 | 0.482 | 0.0232 | 0.105 |
| CLDN4 vs TACSTD2 (comparator) | 22 | 0.289 | 0.193 | 0.47 |
| SFTPC vs DPT (control) | 22 | -0.225 | 0.313 | 0.47 |
| AT2 score vs DPT (control) | 22 | -0.523 | 0.0124 | 0.105 |

## Primary — sample-level DPT vs RECIST (PR vs PD/SD)

Mann–Whitney on patient means. Positive Δmed / r = PR higher than PD/SD. NE excluded. Honest n is patients, not cells.

| Contrast | n | n_PR | n_PD_SD | Δmed (PR−PD/SD) | r_rb | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DPT vs RECIST (PR vs PD/SD) | 16 | 6 | 10 | +0.172 | +0.467 | 0.147 |
| CLDN4 vs RECIST (PR vs PD/SD) | 16 | 6 | 10 | +0.447 | +0.233 | 0.492 |
| AT2 vs RECIST (PR vs PD/SD) | 16 | 6 | 10 | -0.080 | -0.133 | 0.713 |
| barrier/keratin vs RECIST (PR vs PD/SD) | 16 | 6 | 10 | +0.019 | -0.067 | 0.875 |

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| RECIST-evaluable only CLDN4 vs DPT | 16 | 0.153 | 0.572 |
| malignant-only CLDN4 vs DPT | 22 | 0.076 | 0.736 |
| ADC+SQ CLDN4 vs DPT | 17 | -0.093 | 0.722 |
| drop-SCLC CLDN4 vs DPT | 18 | -0.069 | 0.785 |

| Contrast | n | n_PR | n_PD_SD | Δmed (PR−PD/SD) | r_rb | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| malignant-only DPT vs RECIST | 16 | 6 | 10 | +0.176 | +0.467 | 0.147 |
| ADC+SQ DPT vs RECIST | 12 | 4 | 8 | +0.037 | +0.312 | 0.461 |
| drop-SCLC DPT vs RECIST | 13 | 4 | 9 | +0.038 | +0.278 | 0.503 |
| malignant-only CLDN4 vs RECIST | 16 | 6 | 10 | +0.413 | +0.200 | 0.562 |
| ADC+SQ CLDN4 vs RECIST | 12 | 4 | 8 | +0.248 | +0.188 | 0.683 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)

Emitted: **True**. Paired tertile n=22.

| Paired contrast (high − low) | n_patients | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 22 | +0.559 | 4.77e-07 |
| AT2 high vs low | 22 | +0.036 | 4.2e-05 |
| malignant-like high vs low | 22 | +0.242 | 9.87e-05 |
| DPT high vs low | 22 | -0.002 | 0.924 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- If CLDN4 vs DPT is null, that is **not** lineage proof. The RECIST split is still reported.
- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Author Malignant cells is the published label, **not CNV re-called here**.
- Q4 mixes SCLC with ADC/SQ; histology is reported, not hidden.
- This folder does **not** re-audit PR #320 T/NK ρ.
- No TACSTD2∩CLDN4 both-high gate. GSE148071 not used.
- Slingshot R was not the clock unless `available=True`. DPT is an ordering, not a developmental clock.

## Outputs

- `results/tables/sample_level_spearman.tsv` — CLDN4 vs DPT (done criterion)
- `results/tables/sample_level_recist.tsv` — DPT vs RECIST (done criterion)
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_trajectory_recist.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_extra_sample_cldn4_dpt_recist.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/gse205335_traj_cldn4_recist/requirements.txt
python3 methods/gse205335_traj_cldn4_recist/scripts/download.py \
  --out /tmp/gse205335_traj
python3 methods/gse205335_traj_cldn4_recist/scripts/extract_epithelium.py \
  --data /tmp/gse205335_traj \
  --out /tmp/gse205335_traj/epithelium.h5ad
python3 methods/gse205335_traj_cldn4_recist/scripts/analyze.py \
  --input /tmp/gse205335_traj/epithelium.h5ad \
  --outdir methods/gse205335_traj_cldn4_recist/results \
  --finding methods/gse205335_traj_cldn4_recist/FINDING.md
```

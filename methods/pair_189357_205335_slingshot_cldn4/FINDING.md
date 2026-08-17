# Finding — pair GSE189357+GSE205335, CLDN4-only REAL Slingshot/PAGA

ADDITIVE. **CLDN4 only.** Pair from PR #459 that already differs (malignant CLDN4 %pos vs T/NK, n=31, ρ=−0.478, Q4 r=−0.750). That T/NK cut is **not re-audited**. No dual-high TACSTD2∩CLDN4. GSE131907 is not added. Not CellChat.

Primary clock: **REAL Slingshot (Street 2018)**. Root / start cluster is **not** CLDN4-high (Leiden cluster with highest AT2 score among clusters not in the top tercile of mean CLDN4; prefer author AT2 if present. Never CLDN4-high.). Start Leiden = `4`; forbidden high-CLDN4 Leiden = `0`. Inferential unit = **patient**. Barrier/keratin and IFN scores **exclude CLDN4**. PAGA is the connectivity companion, not the clock.

**What holds (n=31 patients).** CLDN4 vs barrier/keratin (no CLDN4): ρ=0.319, p=0.0803. CLDN4 vs IFN: ρ=-0.077, p=0.68. Along Slingshot PT: CLDN4 ρ=-0.042, p=0.824; barrier ρ=-0.118, p=0.528; IFN ρ=0.116, p=0.535.

## Verdict

Patient-level CLDN4 vs Slingshot PT: n=31, ρ=-0.042, p=0.824. CLDN4 vs barrier/keratin (CLDN4 excluded): n=31, ρ=0.319, p=0.0803. CLDN4 vs IFN (Hallmark IFNα∩IFNγ): n=31, ρ=-0.077, p=0.68. barrier vs PT: n=31, ρ=-0.118, p=0.528. IFN vs PT: n=31, ρ=0.116, p=0.535. Slingshot engine=R_slingshot. Lineages=18. PAGA components at connectivity>0: 1 among 29 Leiden vertices. Not a TACSTD2 redo. No both-high gate. PR #459 not re-ranked.

## Honest n

- Analysis cells after QC (capped ≤350/patient): **n_cells = 10143** (GSE189357 3124, GSE205335 7019).
- Patients: **n_units = 31** (GSE189357 9, GSE205335 22). Do not write n=31 as one cohort.
- Patients with ≥10 epithelial cells used for Spearman: **n = 31**.
- GSE189357 gate: marker epithelium (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0. Author cell types are **absent**.
- GSE205335 gate: author `lineage.total == Epithelial cells` on non-normal tissue. Histology is mixed (ADC/SQ/SCLC/NUT).
- Author AT2 cells in the object: **0** (GSE205335 only).
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms: **n = 11**.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- Slingshot: engine=R_slingshot; 2.10.0.
- Dual-high: not used. GSE131907: not added. PR #459 T/NK Spearman: not re-audited.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: harmonypy on PCA, batch=dataset.
- Slingshot start cluster: `4` (Leiden cluster with highest AT2 score among clusters not in the top tercile of mean CLDN4; prefer author AT2 if present. Never CLDN4-high.).
- PAGA components at connectivity>0: **1** among 29 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN genes: Hallmark IFNα ∩ IFNγ (**CLDN4 out**).

## Lineage table (done criterion)

| lineage | path | n_cells | start | end | mean CLDN4 | mean barrier | mean IFN | ρ CLDN4~PT | ρ barrier~PT | ρ IFN~PT |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Lineage1 | `4->18->8->14->1->19->25->26->11` | 3968 | 4 | 11 | 1.276 | 0.352 | 0.091 | 0.016 | -0.302 | -0.359 |
| Lineage2 | `4->18->8->14->1->12->5->20` | 3829 | 4 | 20 | 1.209 | 0.507 | 0.102 | -0.051 | 0.167 | -0.247 |
| Lineage3 | `4->18->8->14->1->0->24` | 3343 | 4 | 24 | 1.382 | 0.633 | 0.159 | 0.160 | 0.495 | 0.477 |
| Lineage4 | `4->18->8->14->1->7` | 3204 | 4 | 7 | 1.042 | 0.444 | 0.131 | -0.400 | -0.036 | 0.167 |
| Lineage5 | `4->18->8->14->1->3` | 3134 | 4 | 3 | 1.102 | 0.492 | 0.128 | -0.254 | 0.032 | -0.025 |
| Lineage6 | `4->18->8->14->1->9` | 3343 | 4 | 9 | 1.168 | 0.559 | 0.130 | -0.127 | 0.365 | 0.237 |
| Lineage7 | `4->18->8->14->1->13` | 3181 | 4 | 13 | 1.094 | 0.547 | 0.145 | -0.265 | 0.328 | 0.344 |
| Lineage8 | `4->18->8->14->1->15` | 3199 | 4 | 15 | 1.049 | 0.453 | 0.140 | -0.360 | 0.032 | 0.316 |
| Lineage9 | `4->18->8->14->1->16` | 3498 | 4 | 16 | 1.096 | 0.527 | 0.151 | -0.239 | 0.266 | 0.350 |
| Lineage10 | `4->18->8->14->1->17` | 3396 | 4 | 17 | 1.067 | 0.443 | 0.119 | -0.283 | -0.044 | -0.092 |
| Lineage11 | `4->18->8->14->1->21` | 3374 | 4 | 21 | 1.046 | 0.486 | 0.139 | -0.326 | 0.141 | 0.273 |
| Lineage12 | `4->18->8->14->1->22` | 3233 | 4 | 22 | 1.119 | 0.508 | 0.141 | -0.208 | 0.275 | 0.337 |
| Lineage13 | `4->18->8->14->1->23` | 3351 | 4 | 23 | 1.090 | 0.577 | 0.145 | -0.259 | 0.356 | 0.317 |
| Lineage14 | `4->18->8->14->1->27` | 3338 | 4 | 27 | 1.095 | 0.447 | 0.107 | -0.267 | -0.087 | -0.233 |
| Lineage15 | `4->18->8->14->1->28` | 3163 | 4 | 28 | 1.125 | 0.516 | 0.129 | -0.209 | 0.301 | 0.232 |
| Lineage16 | `4->18->8->14->2` | 2034 | 4 | 2 | 1.094 | 0.308 | 0.115 | -0.413 | -0.213 | 0.145 |
| Lineage17 | `4->18->8->14->6` | 2532 | 4 | 6 | 0.931 | 0.133 | 0.165 | -0.534 | -0.512 | 0.588 |
| Lineage18 | `4->18->8->14->10` | 2099 | 4 | 10 | 1.418 | 0.401 | 0.134 | 0.111 | 0.109 | 0.393 |

Machine table: `results/tables/lineage_table.tsv` (18 lineages).

## Primary (patient-level Spearman, BH inside this list)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 31 | -0.042 | 0.824 | 0.824 |
| CLDN4 vs DPT (companion) | 31 | 0.135 | 0.47 | 0.654 |
| CLDN4 vs AT2 score | 31 | -0.202 | 0.275 | 0.504 |
| CLDN4 vs barrier/keratin (no CLDN4) | 31 | 0.319 | 0.0803 | 0.177 |
| CLDN4 vs IFN | 31 | -0.077 | 0.68 | 0.749 |
| CLDN4 vs malignant-like | 31 | 0.425 | 0.0172 | 0.0629 |
| CLDN4 vs TACSTD2 (comparator) | 31 | 0.400 | 0.0258 | 0.0709 |
| barrier/keratin vs Slingshot PT | 31 | -0.118 | 0.528 | 0.654 |
| IFN vs Slingshot PT | 31 | 0.116 | 0.535 | 0.654 |
| SFTPC vs Slingshot PT (control) | 31 | -0.616 | 2.25e-04 | 0.0023 |
| AT2 vs Slingshot PT (control) | 31 | -0.595 | 4.18e-04 | 0.0023 |

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| GSE189357-only CLDN4 vs Slingshot PT | 9 | -0.517 | 0.154 |
| GSE205335-only CLDN4 vs Slingshot PT | 22 | -0.116 | 0.608 |
| GSE189357-only CLDN4 vs IFN | 9 | 0.217 | 0.576 |
| GSE205335-only CLDN4 vs IFN | 22 | -0.153 | 0.497 |
| GSE189357-only CLDN4 vs barrier/keratin | 9 | 0.450 | 0.224 |
| GSE205335-only CLDN4 vs barrier/keratin | 22 | 0.243 | 0.275 |
| ADC/LUAD-only CLDN4 vs Slingshot PT | 23 | 0.006 | 0.979 |
| GSE205335-only IFN vs Slingshot PT | 22 | -0.086 | 0.702 |
| GSE189357-only IFN vs Slingshot PT | 9 | 0.383 | 0.308 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)

Emitted: **True**. Paired tertile n=11.

| Paired contrast (high − low) | n | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 11 | 0.320 | 0.00195 |
| IFN high vs low | 11 | 0.006 | 0.638 |
| AT2 high vs low | 11 | -0.023 | 0.175 |
| Slingshot PT high vs low | 11 | 0.202 | 0.365 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- PR #459 (combo enum vs T/NK) is not re-audited here.
- GSE189357 has no author AT2 labels; the start cluster is AT2-scored, not an author nLung AT2 root.
- GSE205335 is mixed histology (ADC/SQ/SCLC/NUT) and an ICI cohort; this is **not** an ICI / RECIST test.
- Marker epithelium ≠ CNV-called malignant. Author epithelial ≠ CNV.
- No TACSTD2∩CLDN4 both-high gate.
- Do not write “AT2 differentiates into LUAD because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot is an ordering on a Harmony-merged pair, not a developmental clock.

## Outputs

- `results/tables/lineage_table.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_along_pseudotime.png`
- `results/figures/fig_paga.png`
- `results/figures/fig_honest_n.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_lineages.png`
- `results/figures/fig_extra_umaps.png`
- `results/summary.json`

## Reproduce

```bash
pip install -r methods/pair_189357_205335_slingshot_cldn4/requirements.txt
python3 methods/pair_189357_205335_slingshot_cldn4/scripts/download.py \
  --out /tmp/geo_pair_189357_205335
python3 methods/pair_189357_205335_slingshot_cldn4/scripts/extract_epithelium.py \
  --data /tmp/geo_pair_189357_205335 \
  --out /tmp/geo_pair_189357_205335/epithelium.h5ad
python3 methods/pair_189357_205335_slingshot_cldn4/scripts/analyze.py \
  --input /tmp/geo_pair_189357_205335/epithelium.h5ad \
  --outdir methods/pair_189357_205335_slingshot_cldn4/results \
  --finding methods/pair_189357_205335_slingshot_cldn4/FINDING.md
```


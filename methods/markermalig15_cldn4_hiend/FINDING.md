# FINDING — CLDN4-only high-end on the marker-malignant combo (n=15)

**Additive. CLDN4 only. No dual-high.** The marker-malignant combo that already differs is taken as given from PR #290 / #279 and is **not re-audited**:

| combo | k | N | CLDN4 ρ (p, I²) |
|---|---:|---:|---|
| Marker-malignant %pos vs T/NK: GSE253013 + GSE291670 | 2 | **15** | **−0.714 (0.0217, 23%)** |

**n=15 is small. Do not inflate it.** Adjacent-normal GSE253013 libraries are not added. Extra GEO samples are not added. Cell-level p-values are not a 15-patient mixed model. The two cohorts are scored separately (CLDN4 scales differ).

Sender definition is **CLDN4-only** among marker-malignant cells. TACSTD2 is a companion column. Dual-high (TACSTD2-high ∩ CLDN4-high) is **not** the sender set (1,203 companion cells scored, unused).

- GSE253013: median `log1p(CP10k)` CLDN4 = 2.24 → 1,479 high / 1,478 low
- GSE291670: median = 0 (zero-inflated, matches the notes’ %pos scoring) → high = detected (500) / low = undetected (5,015)

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| Given combo (not re-audited) | **15** | 9 GSE253013 tumors + 6 GSE291670 |
| GSE253013 tumor patients | **9** | MRC001–004, MRC006–010 (no MRC005 in the notes) |
| GSE291670 tumors | **6** | MPR-1/2/3 vs Non-MPR-1/2/3 |
| Thin marker-malignant tumors | **2** | MRC004 n_mal=16, MRC007 n_mal=18 (stay in n=15) |
| Marker-malignant cells | **8,472** | CLDN4-high 1,979 / low 6,493 |
| T/NK cells | **77,239** | GSE253013 is CD45-sorted in 6/9 tumors — T/NK fraction is sampling-biased |
| Dual-high companion (not senders) | 1,203 | TACSTD2 scored only |
| CellChatDB v2 pairs in the panel | 2,226 | protein pairs |
| Significant differential directed tests | **129** | 63 GSE253013 + 66 GSE291670; 100 permutations; smallest p = 1/101 = 0.0099 |

Patient/sample table: `results/patient_table.tsv`. Combo membership: `data/combo/samples_n15.tsv`.

### Ligand–receptor table (CellChat-style, Mal CLDN4-high vs low → T/NK)

Hill / mass-action probability on 10% truncated means (CellChatDB v2). The CellChat R package is not run. Full significant table: `results/lr_table.tsv`. All scored pairs: `results/lr_contrasts_all.tsv`.

**GSE253013 outgoing (CLDN4-high higher, permutation p = 0.0099):**

| Pair | Pathway | Class | P high | P low | ΔP |
| --- | --- | --- | ---: | ---: | ---: |
| COL6A1–CD44 | COLLAGEN | other | 0.479 | 0.235 | +0.244 |
| CD55–ADGRE5 | ADGRE | other | 0.503 | 0.277 | +0.226 |
| LAMC2–CD44 | LAMININ | other | 0.441 | 0.256 | +0.185 |
| **LGALS9–PTPRC (CD45)** | GALECTIN | inhibitory | 0.474 | 0.314 | +0.159 |
| **LGALS9–CD44** | GALECTIN | inhibitory | 0.382 | 0.240 | +0.143 |
| APP–CD74 | APP | other | 0.768 | 0.627 | +0.140 |
| LAMB3–CD44 | LAMININ | other | 0.606 | 0.469 | +0.137 |
| CD99–CD99 | CD99 | other | 0.476 | 0.376 | +0.101 |

SFTPD–ADGRE5 and CLEC2B–KLRB1 are higher from CLDN4-low (ΔP −0.190 / −0.123; significant on the low arm).

**GSE291670 outgoing (zero-inflated split):**

| Pair | Pathway | Class | P high | P low | ΔP | p high |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| LAMB3–CD44 | LAMININ | other | 0.212 | 0.004 | +0.208 | 0.0099 |
| LAMC2–CD44 | LAMININ | other | 0.163 | 0.050 | +0.112 | 0.0099 |
| COL4A5–CD44 | COLLAGEN | other | 0.544 | 0.466 | +0.079 | 0.0099 |
| COL4A3–CD44 | COLLAGEN | other | 0.288 | 0.518 | −0.230 | 1 (sig on low) |
| COL4A4–CD44 | COLLAGEN | other | 0.136 | 0.354 | −0.219 | 1 (sig on low) |
| CADM1–CADM1 | CADM | barrier | 0.615 | 0.709 | −0.094 | 1 (sig on low) |

**Shared direction:** LAMB3–CD44 and LAMC2–CD44 are higher from CLDN4-high in **both** cohorts. LGALS9 inhibitory pairs are GSE253013-only (GSE291670 does not call them differential). CD274–PDCD1 is not detected (`expr_prop < 0.10`).

### NicheNet-style ligand activity (CLDN4-high malignant → T/NK IFN / cytotoxicity)

Unsigned NicheNet-v2 regulatory potential. Target set is the a-priori T/NK IFN + cytotoxicity list (18 genes in the T/NK-expressed background of 206). Ligands are those expressed in ≥10% of CLDN4-high marker-malignant cells (n=254 scored). Honest n remains 15.

| ligand | pearson | pearson_p | auroc | aupr |
| --- | ---: | ---: | ---: | ---: |
| LIF | 0.451 | 1.03e-11 | 0.676 | 0.445 |
| HLA-DRA | 0.372 | 3.57e-08 | 0.605 | 0.346 |
| VSIG10 | 0.360 | 1.05e-07 | 0.697 | 0.401 |
| HLA-A | 0.354 | 1.86e-07 | 0.678 | 0.382 |
| IL23A | 0.327 | 1.59e-06 | 0.702 | 0.308 |
| BST2 | 0.310 | 5.61e-06 | 0.682 | 0.400 |
| VTCN1 | 0.301 | 1.14e-05 | 0.600 | 0.277 |
| CLCF1 | 0.297 | 1.45e-05 | 0.705 | 0.374 |
| HLA-F | 0.296 | 1.52e-05 | 0.628 | 0.399 |
| HLA-E | 0.286 | 3.12e-05 | 0.647 | 0.357 |

This is prior activity, not a patient-level test that the ligand is higher in CLDN4-high tumors.

### What is not supported

- Re-auditing ρ=−0.714. That number is given.
- **n > 15.** ANT libraries, extra cohorts, or cell-as-N inflation are out of scope.
- Dual-high (TACSTD2-high and CLDN4-high) as the sender definition.
- A confirmatory ICI-response claim. GSE253013 has **no** public MPR/R labels. GSE291670 is 3 vs 3.
- Patient mixed-model communication. CellChat here is cell-pooled truncated means with a high/low label permutation, run **per cohort**.
- A single pooled Hill network across the two platforms.

### Extra figures

- `results/fig_given_combo_n15.png` — given combo scatter (notes, not re-audited)
- `results/fig_n_honesty.png` — per-sample malignant / T/NK counts (n=15)
- `results/fig_lr_delta.png` — CellChat-style ΔP
- `results/fig_nichenet_activity.png` — NicheNet ligand activity
- `results/fig_patient_cldn4_vs_tnk_ifn.png` — patient-level CLDN4 vs T/NK IFN

### Files

- `results/lr_table.tsv` — CellChat-style differential LR table (the PR done criterion)
- `results/nichenet_ligand_activity.tsv` — NicheNet-style ligand activity
- `results/patient_table.tsv` — honest n=15 sample counts
- `results/summary.json`

## 中文

在已经给出的 marker-malignant 组合（GSE253013 肿瘤 9 例 + GSE291670 6 例，**n=15**，CLDN4 %pos vs T/NK **ρ=−0.714**，不重审）上，做 **CLDN4-only** 高端：CellChat 风格的恶性 CLDN4-high vs low → T/NK，以及 NicheNet 风格的配体活性。**禁止 dual-high。** n=15 很小，不把 ANT 或其它队列加进来把 n 做大。完整 LR 表见 `results/lr_table.tsv`。

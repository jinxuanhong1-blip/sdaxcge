# Results — GSE207422 CellChat-style TACSTD2 vs T/NK

Public UMI only (Hu et al., *Genome Medicine* 2023, PMID 36869384). CellChat R and LIANA were not run. Pairs below are from the Hill probability + 100 high/low permutations in `scripts/analyze.py`. A pair is counted significant if it is detected (`expr_prop ≥ 0.10` on both sides), `P > 0`, and permutation `p < 0.05`. The smallest possible p with 100 permutations is 1/101 = 0.0099.

No edge is drawn unless it meets that rule. Full pair tables are in `results/lr_pairs_*.tsv` and `results/contrast_*.tsv`.

## Matrix and compartments

| Item | n |
| --- | --- |
| Genes × barcodes in the GEO UMI | 24,292 × 92,330 |
| CellChatDB v2 protein pairs with every subunit in the matrix | 1,537 (of 2,239 protein pairs) |
| Post-treatment samples used for splits | 12 (MPR+pCR = 4; NMPR = 8) |
| Post-treatment epithelial cells | 8,407 |
| Post-treatment T+NK | 33,760 |
| Pre-treatment biopsies | excluded from splits (3 samples) |

Epithelial cells are the malignant compartment (author CopyKAT IDs are not on GEO). T and NK are merged as TNK.

BD_immune07 (NMPR, squamous) is **5,121 / 8,407 (60.9%)** of post-treatment epithelial cells and **5,121 / 6,914 (74.1%)** of NMPR epithelial cells. Group means are cell-pooled, so this sample dominates the NMPR and pooled contrasts. Per-sample counts: `results/per_sample_post.tsv`.

## Splits compared

TACSTD2 = `log1p(CP10k)` on epithelial cells. Tertile high/low drops the middle third.

| Split | Mal high | Mal low | TNK | Detected LR | Sig LR (out / in) | Cold/barrier score | Recruit-up in high |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| tertile, all post | 2,803 (12 pts) | 2,803 | 33,760 | 220 | 117 (80 / 37) | 15 | 0 |
| median, all post | 4,203 (12 pts) | 4,204 | 33,760 | 218 | 105 (69 / 36) | 14 | 0 |
| **tertile, NMPR only (kept)** | **2,305 (8 pts)** | **2,305** | **21,223** | **204** | **93 (69 / 24)** | **15** | **0** |
| tertile, MPR only | 498 (4 pts) | 498 | 12,537 | 229 | 119 (86 / 33) | 12 | 3 |
| tertile × MPR (two arms) | 2,703+100 | 1,743+1,060 | 21,223+12,537 | 467 | 272 (200 / 72) | 22* | 2 |

\*Combo score adds NMPR and MPR arms and is not comparable to a single contrast. Global tertiles leave only **100** TACSTD2-high epithelial cells in MPR (MPR epithelium is shifted lower: tertile cuts 0.73 / 1.44 vs NMPR 1.75 / 2.41).

## Kept split

**NMPR-only TACSTD2 tertile.** MPR labels exist; the MPR-only tertile does not match the immune-cold / barrier pattern (CXCL16, CCL5, CCL4 outgoing and IFNG incoming are *higher* from/to TACSTD2-high). The NMPR-only tertile has recruit-up = 0 and CXCL16–CXCR6 lower in high. The pooled post tertile has the same score (15) on 12 patients and is the all-post view.

## Kept split — significant LR counts

Among 1,537 database pairs × 4 directed tests = 6,148 rows:

- Detected: 204
- Significant: 93
- Outgoing Mal → TNK: 69 significant (29 high-arm, 40 low-arm)
- Incoming TNK → Mal: 24 significant (13 high-arm, 11 low-arm)

## Outgoing Mal → T/NK (NMPR tertile)

Barrier / inhibitory pairs **higher** in TACSTD2-high (p = 0.0099 unless noted):

| Pair | Pathway | ΔP (high−low) | P high | P low |
| --- | --- | ---: | ---: | ---: |
| NECTIN2–TIGIT | NECTIN | +0.127 | 0.432 | 0.305 |
| CDH1–ITGAE/ITGB7 | CDH1 | +0.078 | 0.273 | 0.195 |
| JAM1–ITGAL/ITGB2 | JAM | +0.074 | 0.702 | 0.628 |
| NECTIN2–CD226 | NECTIN | +0.028 | 0.068 | 0.040 |
| HLA-E–CD8B | MHC-I | +0.022 | 0.429 | 0.407 |
| HLA-E–CD8A | MHC-I | +0.014 | 0.817 | 0.803 |
| LGALS9–CD45 | GALECTIN | +0.034 (p=0.020) | 0.424 | 0.390 |
| CDH1–KLRG1 | CDH1 | +0.005 | 0.015 | 0.010 |

Recruiting pair **lower** in TACSTD2-high:

| Pair | ΔP | P high | P low | p low |
| --- | ---: | ---: | ---: | ---: |
| CXCL16–CXCR6 | −0.012 | 0.143 | 0.155 | 0.030 |

Largest absolute outgoing deltas are **MHC-II → CD4**, all *lower* in TACSTD2-high (HLA-DRA/DRB1/DPA1/DPB1–CD4, ΔP −0.25 to −0.34). ICAM1–ITGAL/ITGB2 is also lower in high (ΔP −0.16). Those pairs are significant on the low arm; they are in `contrast_tertile_nmpr_outgoing.tsv`.

PVR–TIGIT and PVR–CD226 are significant on the low arm (higher in TACSTD2-low). CD274–PDCD1 is not significant on the NMPR-high arm.

## Incoming T/NK → Mal (NMPR tertile)

| Pair | ΔP | Note |
| --- | ---: | --- |
| FASLG–FAS | −0.0085 | lower into TACSTD2-high (p_low = 0.0099) |
| TGFB1–TGFBR complexes | −0.006 to −0.011 | detected only into TACSTD2-low |
| TNFSF10–TNFRSF10B | +0.005 | TRAIL slightly higher into high |
| CD8A–CEACAM5 | +0.318 | largest incoming |ΔP|; CEACAM5 higher on TACSTD2-high epithelium |

IFNG–IFNGR is **not** significant in the NMPR kept split. It **is** higher into TACSTD2-high in the MPR-only tertile (ΔP +0.077).

## MPR-only tertile (not kept)

498 vs 498 epithelial cells, 4 patients. Significant outgoing includes CXCL16–CXCR6 (+0.106), CCL5–CCR5 (+0.009), CCL4–CCR5 (+0.006) in TACSTD2-high, and IFNG incoming +0.077. That is the opposite recruit / attack direction from NMPR.

## What is not claimed

- These are permutation tests on cell-pooled truncated means, not a patient-level mixed model. n_patients = 8 (NMPR) or 4 (MPR).
- P07 supplies most NMPR epithelial UMIs.
- CellChat R visualizations were not generated. Figures show only pairs that pass p < 0.05.
- Author cell-type and CopyKAT objects were not used.

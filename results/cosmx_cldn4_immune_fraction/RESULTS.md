# CosMx NSCLC: CLDN4 detected vs absent immune fraction, smaller than 0.38×

ADDITIVE, **CLDN4-only**, He et al. 2022 CosMx 960-plex (figshare 25976224, `cosmx_human_nsclc_clustered.h5ad`; 8 sections / 5 patients). The estimand is the previous immune-cell fraction: neighbors inside a ball, index excluded, immune count divided by all other cells, and a malignant cell with no neighbor counted as fraction 0. The previous specification, recomputed here, is **0.377×** at 10 µm (count ≥ 1 vs 0, broad immune labels), 8/8 sections and 5/5 patients. This note does **not** replace the locked 50/100 µm cytotoxic exclusion result, and it does **not** say nearby effectors are muzzled. No private 8-KL. No ICI labels.

## Selection

The grid was fixed before ranking: radii 4, 5, 6, 7, 8, 9, 10, 12 µm; CLDN4 cuts ge1, ge2, ge3, ge5, ge10, pos_top50, pos_top25, pos_top10 (`geK` is raw count ≥ K versus count = 0; `pos_top50/25/10` keeps detected cells at or above the 50th/75th/90th percentile of that section's positive counts, versus count = 0); immune definitions broad, no_neutrophil, no_granulocyte, no_macrophage, no_macrophage_neutrophil, lymphoid, t_nk, t_nk_no_treg, cytotoxic, cd8, myeloid. 704 specs had both arms at n ≥ 30 in every section. A spec is eligible when the detected arm is strictly lower in **8/8 sections and 5/5 patients** and the absent-arm mean is ≥ 0.005 in every section (59 specs).

The smallest eligible ratio is **0.203** (9 µm, `ge5`, `broad`). It is not the headline. Its smallest detected-arm section mean is 0.0000. The contact fraction is not lower in every section. A section mean of zero, or a contact fraction that flips, shrinks the ratio of means without a compositional contrast in every section.

The headline is the smallest eligible ratio that also has a detected-arm mean above 0 in every section and a contact fraction lower in 8/8 sections and 5/5 patients (46 specs). Section p = 0.0078 and patient p = 0.0625 whenever every unit has the same sign, so those p-values do not choose the spec. FOV p-values are nominal (FOVs sit inside 5 patients; the grid was searched). Radii 4–8 µm produced no eligible spec.

## Winning contrast

- **Index:** patient-matched malignant cells (`tumor 5/6/9/12/13`).
- **Cutoff:** `pos_top25` within each section. Positive cells below the threshold are in neither arm.
- **Radius:** **9 µm** (global centroids × 0.18 µm/pixel).
- **Immune definition `broad`:** B-cell, NK, T CD4 memory, T CD4 naive, T CD8 memory, T CD8 naive, Treg, mDC, macrophage, mast, monocyte, neutrophil, pDC, plasmablast.
- **Ratio:** **0.229×** (detected 0.0112 / absent 0.0490), versus **0.377×** for 10 µm, count ≥ 1 vs 0, the same immune labels (detected 0.0242 / absent 0.0643).
- **Concordance:** sections 8/8, patients 5/5. Weakest section ratio 0.680.
- **Contact fraction** (degree ≥ 1 only): ratio 0.313, also 8/8 and 5/5.
- **Neighborhood size:** mean degree 0.583 (detected) and 0.964 (absent), degree ratio 0.605. Detected cells have fewer neighbors. The contact ratio above is the part that remains after those empty neighborhoods are dropped.
- **FOVs:** 190/204 lower, nominal Wilcoxon p = 5.96e-33. The reproduced 10 µm spec is 200/211, p = 4.13e-33.

| Section | Patient | Count ≥ | n detected | n absent | Detected | Absent | Ratio | Contact ratio | Degree ratio |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 4 | 3741 | 5226 | 0.0121 | 0.0244 | 0.496 | 0.648 | 0.617 |
| LUAD-5 R2 | Lung5 | 3 | 4917 | 5537 | 0.0111 | 0.0190 | 0.588 | 0.733 | 0.650 |
| LUAD-5 R3 | Lung5 | 4 | 2601 | 5908 | 0.0138 | 0.0203 | 0.680 | 0.912 | 0.629 |
| LUSC-6 | Lung6 | 2 | 4895 | 49873 | 0.0029 | 0.0053 | 0.539 | 0.786 | 0.647 |
| LUAD-9 R1 | Lung9 | 4 | 8450 | 11030 | 0.0047 | 0.0594 | 0.079 | 0.117 | 0.559 |
| LUAD-9 R2 | Lung9 | 3 | 20216 | 36954 | 0.0030 | 0.0495 | 0.060 | 0.079 | 0.669 |
| LUAD-12 | Lung12 | 3 | 1990 | 11007 | 0.0126 | 0.1205 | 0.105 | 0.204 | 0.443 |
| LUAD-13 | Lung13 | 4 | 6377 | 6853 | 0.0293 | 0.0933 | 0.314 | 0.461 | 0.562 |

Stable specs: 46 of 59 eligible. The ten smallest stable ratios:

| Radius | Cut | Immune definition | Ratio | Contact ratio | Max section ratio | Degree ratio | Min n detected |
|---:|---|---|---:|---:|---:|---:|---:|
| 9 | pos_top25 | broad | 0.229 | 0.313 | 0.680 | 0.605 | 1990 |
| 9 | ge3 | broad | 0.243 | 0.323 | 0.627 | 0.634 | 1719 |
| 9 | ge2 | broad | 0.273 | 0.345 | 0.563 | 0.706 | 3723 |
| 9 | pos_top50 | broad | 0.274 | 0.344 | 0.601 | 0.713 | 3723 |
| 10 | pos_top25 | broad | 0.276 | 0.348 | 0.770 | 0.653 | 1990 |
| 10 | ge3 | broad | 0.280 | 0.343 | 0.693 | 0.679 | 1719 |
| 10 | ge5 | no_neutrophil | 0.287 | 0.377 | 0.720 | 0.565 | 386 |
| 10 | ge3 | myeloid | 0.294 | 0.361 | 0.769 | 0.679 | 1719 |
| 10 | pos_top25 | no_neutrophil | 0.305 | 0.390 | 0.630 | 0.653 | 1990 |
| 10 | ge2 | broad | 0.307 | 0.365 | 0.592 | 0.745 | 3723 |

Holding the cut at count ≥ 1 vs 0 and the immune labels at the broad set, 9 µm is **0.354×** and 10 µm is **0.377×**. Shortening the radius by 1 µm does not produce the drop above. The drop comes from requiring a higher CLDN4 count in the detected arm. Narrowing the immune labels did not beat the broad set under the stability rule.

## Smaller ratios that missed the absent-arm floor

The grid contains smaller primary ratios that stay 8/8 and 5/5, including on the contact fraction, with a positive detected mean in every section. The smallest is **0.132** at 8 µm, cut `pos_top10`, immune definition `no_macrophage` (contact ratio 0.235). It is not eligible: the quietest absent-arm section mean is 0.0015, below the 0.005 floor. That floor is what keeps a near-empty absent arm from manufacturing a small ratio. The smallest contact ratio if a detected-arm section mean of zero is allowed is 0.215 at 8 µm, `ge5`, `no_macrophage` (detected-arm minimum 0.0000; primary ratio 0.127).

## What this does not claim

- Not a re-estimate of the locked 50/100 µm cytotoxic **count** ratio.
- Not muzzling of GZMB, PRF1, NKG7, or IFNG in the effector cells that are present.
- Not a smaller section-level p-value. Every 8/8 spec sits on the same Wilcoxon floor.
- FOV p-values are nominal. The sign that was required is 8/8 sections and 5/5 patients.
- The 0.203× grid minimum is a real computed ratio. It is not the headline, because one section's detected mean is zero and the contact fraction is not 8/8.
- No ICI labels. No private 8-KL.

```bash
python3 scripts/download_cosmx_nsclc_h5ad.py
python3 scripts/cosmx_cldn4_immune_fraction_grid.py
```


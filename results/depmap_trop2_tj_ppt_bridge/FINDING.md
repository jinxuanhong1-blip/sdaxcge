# DepMap/CCLE PPT funnel: TROP2–CLDN4/TJ protein + immune scores in lung lines

Numbers below are written by `analyze.py` from `tables/key_stats.json` and `tables/correlations.csv`. They are not typed by hand.

## Question

For a PPT resistance / tight-junction bridge slide, do **public DepMap/CCLE lung cell lines** support:

1. TROP2 coexpression with CLDN4 protein and locked TJ gene scores, and
2. TROP2 coexpression with cancer-cell IFN / MHC-I / immune-ligand scores,
3. with TJ scores explaining (attenuating) the TROP2–immune associations?

Cultured lines have no T-cell infiltrate and no IFN treatment. Partials are association contrasts, not mediation.

## Locked protein check (do not change)

Gygi TenPx `_LUNG_` columns, TACSTD2 and CLDN4 both quantified, replicates not collapsed: **Spearman ρ = 0.693** (n=45, p=1.31e-07). Rounds to 0.69 = True.

Model-mapped lung cell lines with both proteins: ρ = 0.704 (n=44, p=9.81e-08).

RPPA antibody info: mentions CLDN4 = False, TACSTD2/TROP2 = False, Claudin-7 = True. RPPA is not used for TROP2–CLDN4.

## Protein: TROP2 vs TJ genes / scores

Gygi lung cell lines mapped to DepMap models: n = 76. With TACSTD2 protein: n = 76. Protein TJ genes used in TJ_EPITHELIAL score: CLDN1, CLDN3, CLDN4, CLDN7, OCLN, MARVELD2, MARVELD3, TJP1, TJP2, TJP3, F11R, JAM2, JAM3, CGN, CGNL1, CRB3, LSR.

| cohort | partner | Spearman |
|---|---|---|
| lung | CLDN4 | +0.704 (n=44, p=9.81e-08); boot CI [+0.49, +0.83] |
| lung | TJ epithelial (18) | +0.648 (n=76, p=2.56e-10); boot CI [+0.48, +0.77] |
| lung | TJ TISMO (7) | +0.727 (n=76, p=1.03e-13) |
| lung | CLDN4 TJ edge (11) | +0.667 (n=76, p=4.83e-11) |
| lung | CLDN7 | +0.702 (n=65, p=7.48e-11) |
| lung | CLDN1 | +0.122 (n=69, p=0.3184) |
| lung | OCLN | +0.553 (n=76, p=2.26e-07) |
| lung | TJP1 | +0.340 (n=76, p=0.0027) |
| lung | F11R | +0.609 (n=76, p=5.28e-09) |
| lung | CDH1 | +0.702 (n=76, p=1.64e-12) |
| lung | EPCAM | +0.571 (n=76, p=7.13e-08) |
| NSCLC | CLDN4 | +0.727 (n=35, p=7.73e-07); boot CI [+0.52, +0.84] |
| NSCLC | TJ epithelial (18) | +0.710 (n=63, p=6.98e-11); boot CI [+0.54, +0.82] |
| NSCLC | TJ TISMO (7) | +0.784 (n=63, p=2.83e-14) |
| NSCLC | CLDN4 TJ edge (11) | +0.691 (n=63, p=3.71e-10) |
| NSCLC | CLDN7 | +0.787 (n=54, p=1.66e-12) |
| NSCLC | CLDN1 | +0.079 (n=57, p=0.5607) |
| NSCLC | OCLN | +0.542 (n=63, p=4.51e-06) |
| NSCLC | TJP1 | +0.242 (n=63, p=0.0564) |
| NSCLC | F11R | +0.664 (n=63, p=2.99e-09) |
| NSCLC | CDH1 | +0.750 (n=63, p=1.50e-12) |
| NSCLC | EPCAM | +0.659 (n=63, p=4.45e-09) |

## RNA: TROP2 vs TJ scores (DepMap 24Q4 lung cell lines)

Lung cell lines with TACSTD2 RNA: n = 214 (NSCLC 143, LUAD 80). TJ epithelial genes used: 18/18.

| cohort | TJ score | Spearman | after EPCAM | after keratins |
|---|---|---|---|---|
| lung | CLDN4 | +0.607 (n=214, p=5.64e-23) | +0.506 (n=214, p=3.23e-15); |partial|/|unadj|=0.83; retained | +0.087 (n=214, p=0.2107); |partial|/|unadj|=0.14; attenuated |
| lung | TJ epithelial (18) | +0.506 (n=214, p=2.67e-15) | +0.355 (n=214, p=1.05e-07); |partial|/|unadj|=0.70; retained | +0.069 (n=214, p=0.3213); |partial|/|unadj|=0.14; attenuated |
| lung | TJ TISMO (7) | +0.551 (n=214, p=2.33e-18) | +0.441 (n=214, p=1.46e-11); |partial|/|unadj|=0.80; retained | +0.105 (n=214, p=0.1277); |partial|/|unadj|=0.19; attenuated |
| lung | CLDN4 TJ edge (11) | +0.564 (n=214, p=2.10e-19) | +0.457 (n=214, p=2.06e-12); |partial|/|unadj|=0.81; retained | +0.122 (n=214, p=0.0782); |partial|/|unadj|=0.22; attenuated |
| NSCLC | CLDN4 | +0.667 (n=143, p=8.81e-20, q=2.06e-19); boot CI [+0.56, +0.75] | +0.474 (n=143, p=2.65e-09); boot CI [+0.32, +0.61]; |partial|/|unadj|=0.71; retained | +0.242 (n=143, p=0.0040); boot CI [+0.04, +0.43]; |partial|/|unadj|=0.36; reduced_but_still_associated |
| NSCLC | TJ epithelial (18) | +0.668 (n=143, p=7.45e-20, q=2.06e-19); boot CI [+0.55, +0.76] | +0.441 (n=143, p=4.01e-08); boot CI [+0.26, +0.59]; |partial|/|unadj|=0.66; retained | +0.304 (n=143, p=0.0003); boot CI [+0.08, +0.48]; |partial|/|unadj|=0.46; reduced_but_still_associated |
| NSCLC | TJ TISMO (7) | +0.663 (n=143, p=1.84e-19, q=3.22e-19); boot CI [+0.55, +0.76] | +0.429 (n=143, p=9.88e-08); boot CI [+0.27, +0.56]; |partial|/|unadj|=0.65; retained | +0.295 (n=143, p=0.0004); boot CI [+0.09, +0.48]; |partial|/|unadj|=0.44; reduced_but_still_associated |
| NSCLC | CLDN4 TJ edge (11) | +0.696 (n=143, p=4.92e-22, q=3.44e-21); boot CI [+0.59, +0.78] | +0.499 (n=143, p=2.54e-10); boot CI [+0.34, +0.62]; |partial|/|unadj|=0.72; retained | +0.341 (n=143, p=3.76e-05); boot CI [+0.15, +0.52]; |partial|/|unadj|=0.49; reduced_but_still_associated |
| LUAD | CLDN4 | +0.756 (n=80, p=5.11e-16) | +0.515 (n=80, p=1.19e-06); |partial|/|unadj|=0.68; retained | +0.398 (n=80, p=0.0003); |partial|/|unadj|=0.53; retained |
| LUAD | TJ epithelial (18) | +0.740 (n=80, p=4.63e-15) | +0.419 (n=80, p=0.0001); |partial|/|unadj|=0.57; retained | +0.391 (n=80, p=0.0004); |partial|/|unadj|=0.53; retained |
| LUAD | TJ TISMO (7) | +0.726 (n=80, p=2.70e-14) | +0.411 (n=80, p=0.0002); |partial|/|unadj|=0.57; retained | +0.367 (n=80, p=0.0010); |partial|/|unadj|=0.51; retained |
| LUAD | CLDN4 TJ edge (11) | +0.763 (n=80, p=1.82e-16) | +0.485 (n=80, p=5.90e-06); |partial|/|unadj|=0.64; retained | +0.432 (n=80, p=8.84e-05); |partial|/|unadj|=0.57; retained |

## RNA: TROP2 vs immune-related scores

Hallmark IFN-γ genes used: 198/200. Immune-ligand members used: CD274, PDCD1LG2, CXCL9, CXCL10, CXCL11, HLA-E, IDO1.

| cohort | score | TROP2 Spearman |
|---|---|---|
| lung | Hallmark IFN-γ | +0.482 (n=214, p=7.23e-14) |
| lung | MHC-I (HLA-A/B/C+B2M) | +0.314 (n=214, p=2.88e-06) |
| lung | Immune-ligand score | +0.550 (n=214, p=2.71e-18) |
| lung | Compact ISG | +0.449 (n=214, p=5.02e-12) |
| lung | CD274 | +0.537 (n=214, p=2.27e-17) |
| NSCLC | Hallmark IFN-γ | +0.339 (n=143, p=3.38e-05, q=3.95e-05); boot CI [+0.18, +0.48] |
| NSCLC | MHC-I (HLA-A/B/C+B2M) | +0.154 (n=143, p=0.0671, q=0.0671); boot CI [-0.03, +0.31] |
| NSCLC | Immune-ligand score | +0.405 (n=143, p=5.09e-07, q=7.13e-07); boot CI [+0.26, +0.53] |
| NSCLC | Compact ISG | +0.275 (n=143, p=0.0009) |
| NSCLC | CD274 | +0.400 (n=143, p=7.22e-07) |
| LUAD | Hallmark IFN-γ | +0.268 (n=80, p=0.0161) |
| LUAD | MHC-I (HLA-A/B/C+B2M) | +0.078 (n=80, p=0.4928) |
| LUAD | Immune-ligand score | +0.352 (n=80, p=0.0014) |
| LUAD | Compact ISG | +0.220 (n=80, p=0.0499) |
| LUAD | CD274 | +0.432 (n=80, p=6.38e-05) |

## Bridge: does TJ attenuate TROP2–immune (RNA NSCLC)?

Primary FDR family = the six NSCLC RNA partials (TROP2 vs IFNG/MHC1/IMMUNE given TJ_EPITHELIAL, and TROP2 vs TJ_EPITHELIAL/TJ_TISMO/CLDN4 given IFNG).

| immune score | TROP2 unadjusted | TROP2 \| TJ epithelial | TROP2 \| CLDN4 | TJ epithelial \| IFNG context |
|---|---|---|---|---|
| Hallmark IFN-γ | +0.339 (n=143, p=3.38e-05, q=3.95e-05); boot CI [+0.18, +0.48] | +0.277 (n=143, p=0.0009, q=0.0010); boot CI [+0.11, +0.43]; |partial|/|unadj|=0.82; retained | +0.251 (n=143, p=0.0026); boot CI [+0.08, +0.40]; |partial|/|unadj|=0.74; retained | +0.207 (n=143, p=0.0133); boot CI [+0.04, +0.37] |
| MHC-I (HLA-A/B/C+B2M) | +0.154 (n=143, p=0.0671, q=0.0671); boot CI [-0.03, +0.31] | +0.081 (n=143, p=0.3378, q=0.3378); boot CI [-0.08, +0.24]; |partial|/|unadj|=0.53; no_unadjusted_association | +0.058 (n=143, p=0.4915); boot CI [-0.10, +0.22]; |partial|/|unadj|=0.38; no_unadjusted_association | +0.140 (n=143, p=0.0943); boot CI [-0.03, +0.30] |
| Immune-ligand score | +0.405 (n=143, p=5.09e-07, q=7.13e-07); boot CI [+0.26, +0.53] | +0.405 (n=143, p=5.87e-07, q=8.22e-07); boot CI [+0.25, +0.54]; |partial|/|unadj|=1.00; retained | +0.347 (n=143, p=2.31e-05); boot CI [+0.20, +0.48]; |partial|/|unadj|=0.86; retained | +0.162 (n=143, p=0.0530); boot CI [-0.02, +0.32] |

TROP2 vs TJ scores given IFN-γ (does immune score explain the TJ coexpression?):

| TJ score | TROP2 unadjusted | TROP2 \| IFNG |
|---|---|---|
| CLDN4 | +0.667 (n=143, p=8.81e-20, q=2.06e-19); boot CI [+0.56, +0.75] | +0.642 (n=143, p=6.87e-18, q=1.20e-17); boot CI [+0.53, +0.74]; |partial|/|unadj|=0.96; retained |
| TJ epithelial (18) | +0.668 (n=143, p=7.45e-20, q=2.06e-19); boot CI [+0.55, +0.76] | +0.650 (n=143, p=2.06e-18, q=4.81e-18); boot CI [+0.52, +0.75]; |partial|/|unadj|=0.97; retained |
| TJ TISMO (7) | +0.663 (n=143, p=1.84e-19, q=3.22e-19); boot CI [+0.55, +0.76] | +0.657 (n=143, p=6.65e-19, q=2.33e-18); boot CI [+0.54, +0.75]; |partial|/|unadj|=0.99; retained |
| CLDN4 TJ edge (11) | +0.696 (n=143, p=4.92e-22, q=3.44e-21); boot CI [+0.59, +0.78] | +0.669 (n=143, p=9.90e-20, q=6.93e-19); boot CI [+0.56, +0.76]; |partial|/|unadj|=0.96; retained |

## Protein immune layer (honest n)

Protein immune-ligand members used: CD274, PDCD1LG2, HLA-E. CXCL9/10/11 are often absent from Gygi.

| cohort | score | TROP2 Spearman | TROP2 \| CLDN4 |
|---|---|---|---|
| lung | Hallmark IFN-γ | +0.416 (n=76, p=0.0002); boot CI [+0.20, +0.60] | +0.210 (n=44, p=0.1764); boot CI [-0.15, +0.53]; |partial|/|unadj|=0.51; ns_after_adjustment_point_estimate_still_half_or_more |
| lung | MHC-I (HLA-A/B/C+B2M) | +0.283 (n=76, p=0.0133) | +0.107 (n=44, p=0.4940); boot CI [-0.25, +0.43]; |partial|/|unadj|=0.38; attenuated |
| lung | Immune-ligand score | +0.402 (n=72, p=0.0005); boot CI [+0.17, +0.59] | +0.170 (n=41, p=0.2930); boot CI [-0.21, +0.52]; |partial|/|unadj|=0.42; attenuated |
| NSCLC | Hallmark IFN-γ | +0.333 (n=63, p=0.0077); boot CI [+0.10, +0.54] | +0.087 (n=35, p=0.6228); |partial|/|unadj|=0.26; attenuated |
| NSCLC | MHC-I (HLA-A/B/C+B2M) | +0.179 (n=63, p=0.1604) | -0.029 (n=35, p=0.8686); |partial|/|unadj|=0.16; no_unadjusted_association |
| NSCLC | Immune-ligand score | +0.338 (n=62, p=0.0073); boot CI [+0.08, +0.57] | +0.090 (n=34, p=0.6182); |partial|/|unadj|=0.27; attenuated |

## PPT claim language (honest)

- KEEP on PPT: Gygi lung TROP2–CLDN4 protein Spearman ρ=0.693 (n=45, p=1.31e-07).
- KEEP on PPT (RNA NSCLC): TROP2 coexpresses with CLDN4 ρ=+0.667 (n=143, p=8.81e-20); TJ epithelial (18) ρ=+0.668 (n=143, p=7.45e-20); TJ TISMO (7) ρ=+0.663 (n=143, p=1.84e-19); CLDN4 TJ edge (11) ρ=+0.696 (n=143, p=4.92e-22).
- KEEP on PPT (RNA NSCLC): TROP2 coexpresses with Hallmark IFN-γ ρ=+0.339 (n=143, p=3.38e-05); Immune-ligand score ρ=+0.405 (n=143, p=5.09e-07).
- Bridge (TROP2–immune | TJ epithelial, RNA NSCLC): Hallmark IFN-γ: unadj +0.339 → |TJ epithelial +0.277 (retained); Immune-ligand score: unadj +0.405 → |TJ epithelial +0.405 (retained).
- DO NOT write 'PPT resistance' as a DepMap result: these are untreated cultured lines without ICI labels.

## Not claimed

- ICI resistance in patients (no response labels in DepMap).
- Immune exclusion or spatial TJ barrier (cell lines only).
- RPPA n≈118 TROP2–CLDN4 (no TROP2/CLDN4 antibodies).
- Private KL scRNA, KD co-culture, or PDX protein.

## Reproduce

```bash
python3 scripts/depmap_trop2_tj_ppt_bridge/download.py
python3 scripts/depmap_trop2_tj_ppt_bridge/analyze.py
```


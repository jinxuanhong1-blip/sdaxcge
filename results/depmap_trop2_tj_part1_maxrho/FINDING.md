# Part1 bridge max ρ: DepMap/CCLE TROP2–CLDN4 protein + TROP2–TJ gene

Numbers below are written by `analyze.py` from `summary.json` / `correlations.csv`. They are not typed by hand.

## Question

For the **Part1 bridge** (TROP2 coexpresses with a tight-junction program), what is the largest pre-specified Spearman ρ for (1) Gygi TACSTD2–CLDN4 **protein** and (2) DepMap 24Q4 TACSTD2 versus locked **TJ gene** scores, with honest n?

## Rules (fixed before sorting)

- Minimum headline n = 15; minimum residual df for partials = 10.
- Eligible: histology slices and histology partials. Epithelial/keratin partials and peptide QC are reported and cannot win.
- Locked TJ gene sets only (TJ epithelial 18, TJ TISMO 7, CLDN4 TJ edge 11, plus CLDN4 alone). No gene dropping to raise ρ.
- No imputation. No removing individual cell lines to raise ρ.
- No n = 118 TROP2–CLDN4 protein row exists (RPPA lacks TROP2/CLDN4 antibodies).

## Protein maximum (TROP2–CLDN4)

- **Locked baseline (rounds to 0.69):** S1 Lung complete cases **ρ = 0.693**, n = 45, p = 1.31e-07, 95% CI 0.48–0.82.
- **Histology / partial maximum:** NSCLC Primary | OncotreeSubtype: **ρ = 0.846**, n = 17, p = 7.13e-05, 95% CI 0.54–0.96 (BH q = 9.87e-05).
- Peptide ≥ 2 sensitivity (not eligible for max): **ρ = 0.829**, n = 25, p = 3.01e-07, 95% CI 0.60–0.93; same-size null p = 2.30e-02.
- S1 Lung lines = 77; complete cases = 45. Any n = 118 in protein grid: **False**.

| Analysis | n | ρ | p | eligible |
|---|---:|---:|---:|:---:|
| NSCLC Primary | OncotreeSubtype | 17 | 0.846 | 7.13e-05 | yes |
| NSCLC Primary | 17 | 0.821 | 5.35e-05 | yes |
| Oncotree LUAD | PrimaryOrMetastasis | 22 | 0.782 | 2.79e-05 | yes |
| Oncotree LUAD | 22 | 0.782 | 1.71e-05 | yes |
| Oncotree NSCLC | OncotreeSubtype | 35 | 0.751 | 7.51e-07 | yes |
| S1 Lung Primary | 20 | 0.744 | 1.67e-04 | yes |
| Oncotree NSCLC | OncotreeSubtype + PrimaryOrMetastasis | 35 | 0.744 | 1.58e-06 | yes |
| S1 Lung Primary | OncotreePrimaryDisease | 20 | 0.728 | 4.13e-04 | yes |

## RNA / gene maximum (TROP2–TJ)

- Lung lines with TACSTD2 RNA: n = 214 (NSCLC 143, LUAD 80).
- TJ genes used: {"TJ_EPITHELIAL": 18, "TJ_TISMO": 7, "CLDN4_TJ_EDGE": 11}.
- Powered NSCLC baselines: CLDN4 **ρ = 0.667**, n = 143, p = 8.81e-20, 95% CI 0.56–0.75; CLDN4 TJ edge **ρ = 0.695**, n = 143, p = 5.98e-22, 95% CI 0.59–0.78.
- **TJ gene maximum (eligible):** LUAD Primary vs CLDN4: **ρ = 0.791**, n = 45, p = 1.05e-10, 95% CI 0.62–0.90 (BH q = 1.49e-10).

| Analysis | n | ρ | p | eligible |
|---|---:|---:|---:|:---:|
| LUAD Primary vs CLDN4 | 45 | 0.791 | 1.05e-10 | yes |
| LUAD Metastatic vs CLDN4 TJ edge (11) | 35 | 0.779 | 3.57e-08 | yes |
| LUAD Metastatic vs TJ epithelial (18) | 35 | 0.771 | 6.07e-08 | yes |
| LUAD vs CLDN4 TJ edge (11) | 80 | 0.765 | 1.46e-16 | yes |
| LUAD vs CLDN4 TJ edge (11) | PrimaryOrMetastasis | 80 | 0.763 | 2.91e-16 | yes |
| LUAD vs CLDN4 | PrimaryOrMetastasis | 80 | 0.758 | 6.29e-16 | yes |
| LUAD vs CLDN4 | 80 | 0.756 | 5.11e-16 | yes |
| LUAD vs TJ epithelial (18) | 80 | 0.753 | 8.15e-16 | yes |

Keratin partials on NSCLC (not eligible for max):

| Partner | unadjusted ρ | keratins partial ρ | n |
|---|---:|---:|---:|
| CLDN4 | 0.667 | 0.242 | 143 |
| TJ epithelial (18) | 0.676 | 0.315 | 143 |
| CLDN4 TJ edge (11) | 0.695 | 0.341 | 143 |

## Part1 bridge language (honest)

- KEEP: Gygi lung TROP2–CLDN4 protein ρ = 0.693 (n=45). Maximum eligible in the locked grid: NSCLC Primary | OncotreeSubtype ρ = 0.846 (n=17).
- KEEP: DepMap NSCLC RNA TROP2–CLDN4 ρ = 0.667 (n=143); TROP2–CLDN4 TJ edge ρ = 0.695 (n=143). Maximum eligible TJ gene coexpression: LUAD Primary vs CLDN4 ρ = 0.791 (n=45).
- DO NOT claim ICI resistance, immune exclusion, or RPPA n≈118 from DepMap/CCLE alone.
- Keratin/EPCAM partials shrink RNA TJ ρ; they are caveats, not the maximum.

## Reproduce

```bash
python3 -m pip install -r scripts/depmap_trop2_tj_part1_maxrho/requirements.txt
python3 scripts/depmap_trop2_tj_part1_maxrho/download.py
python3 scripts/depmap_trop2_tj_part1_maxrho/analyze.py
```


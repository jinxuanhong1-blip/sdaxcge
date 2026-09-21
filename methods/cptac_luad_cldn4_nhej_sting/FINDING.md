# Finding — CPTAC LUAD protein: CLDN4 vs DNA-PKcs / Ku / STING / TBK1 / IRF3 / HLA

**Protein vs protein.** Public CPTAC TMT freeze v1.2, **LUAD only**. Predictor is **CLDN4 protein**. **No TACSTD2 gate.** Honest pairwise-complete n.

Treatment-naive surgical LUAD (Gillette *Cell* 2020). **No ICI labels.** These correlations are not immunotherapy outcomes.

## Strongest match of the requested pattern

Requested direction: **CLDN4-low protein ↔ DNA-PKcs/Ku lower and STING/TBK1/IRF3/HLA higher.**

The pre-specified continuous 9-gene Spearman (next section) does not show that pattern at BH-FDR. A follow-up search changed the panel, the CLDN4 cut, and the correlation, without imputing CLDN4 and without filling in missing endpoint values. The strongest real match is a **lowest-quartile vs rest** contrast on measured CLDN4 protein.

**Selected cut.** Among 79 tumors with CLDN4 protein, Q1 (n=**20**) vs the other 59. Zero of those 20 are imputed. Mann–Whitney, two-sided. Panel values are the mean of within-cohort z-scores. Δ is median(Q1) − median(rest).

| Panel | Members | Δ median z | p | Direction vs CLDN4-low |
|---|---|---:|---:|---|
| DNA-PKcs + LIG4 | PRKDC, LIG4 | **−0.262** | **0.025** | lower |
| DNA-PKcs + Ku80 + LIG4 | PRKDC, XRCC5, LIG4 | **−0.519** | **0.033** | lower |
| STING + TBK1 + IRF3 | STING1, TBK1, IRF3 | **+0.418** | **0.0044** | higher |
| STING + TBK1 + IRF3 + IRF7 | those four | **+0.500** | **3.9×10⁻⁴** | higher |

DNA-PKcs+LIG4 is the NHEJ panel whose worse arm is smallest. Adding Ku80 keeps the same direction (p=0.033). **Ku70 does not:** XRCC6 Q1−rest Δ = +0.008, p=0.74. **HLA-A/B/C do not:** adding them to STING+TBK1+IRF3 wipes the immune contrast (Δ = −0.003, p=0.52). HLA-B and HLA-C alone are the wrong sign on this cut.

Single proteins, same Q1 vs rest, raw log2 abundance (not z):

| Protein | Δ median (Q1−rest) | p | Matches CLDN4-low pattern? |
|---|---:|---:|---|
| DNA-PKcs | −0.092 | 0.14 | yes, not significant alone |
| Ku80 | −0.033 | 0.21 | yes, not significant alone |
| Ku70 | +0.008 | 0.74 | no |
| LIG4 | −0.086 | 0.092 | yes |
| STING | +0.186 | 0.15 | yes |
| TBK1 | +0.069 | 0.052 | yes |
| IRF3 | +0.091 | 0.042 | yes |
| IRF7 | +0.215 | 0.012 | yes |
| HLA-A | +0.020 | 0.74 | flat |
| HLA-B | −0.131 | 0.65 | no |
| HLA-C | −0.031 | 0.46 | no |

**Search-wide check for this family.** All nonempty subsets of {DNA-PKcs, Ku80, Ku70, LIG4} × {STING, TBK1, IRF3, IRF7}, four CLDN4 cuts (Q1 vs rest, Q1 vs Q4, quintiles, median). 200 permutations of CLDN4 across the 110 tumors. The best worse-arm one-sided p in the real data is 0.013. Permutation p = **0.025** (null median 0.17). That permutation covers this family only, not every gene in the matrix.

Imputing the 31 missing CLDN4 values (minimum, minimum−1, or a Perseus downshift) did not beat this complete-case quartile. Q1 vs Q4 on DNA-PKcs+Ku80 vs STING+TBK1+IRF3 is the same direction and weaker (NHEJ p=0.31, immune p=0.041, 20 vs 20).

**Phospho layer (same tumors, gene-median of sites).** STING phosphosites are quantified in only 44/110 tumors, so STING is not in this phospho contrast. Ku80 phospho vs IRF3 phospho, Spearman n=79: ρ = +0.202 (p=0.074) and ρ = −0.234 (p=0.038). DNA-PKcs+Ku70+Ku80 phospho, Q1 vs Q4 (20 vs 20): Δ = −0.307, p=0.047; IRF3 phospho Δ = +0.252, p=0.091. Named functional sites do not carry it (TBK1 S172 n=61, ρ=+0.14, p=0.28; PRKDC S3205 n=79, ρ=−0.028, p=0.81; PRKDC T2609 n=12). A 200-permutation check of the fully observed phospho gene subsets gives p=0.14 for that layer’s best pair.

## Verdict on the pre-specified continuous test

The pre-specified **9-protein** family is **null at BH-FDR 0.05** (n=**79**). DNA-PKcs, Ku70, Ku80, and HLA-A/B/C protein do not track CLDN4 protein on a continuous Spearman. STING, TBK1, and IRF3 are the same sign (inverse) and each misses that FDR (closest: TBK1 ρ=−0.217, p=0.055, q=0.33; 95% CI includes 0).

A pre-specified mean of STING+TBK1+IRF3 is inverse before purity adjustment (ρ=−0.289, p=0.0098, q=0.030 among **3 composites**). Those three proteins are **not co-abundant** in this table (STING vs TBK1 ρ=0.07; STING vs IRF3 ρ=−0.02), so the mean is not a STING-pathway score. After a WES-purity residual the mean is still inverse (ρ=−0.255, p=0.025) but q=0.075 in the composite family. Do not write a single-protein STING, TBK1, or IRF3 result from this slice.

## Primary table (CLDN4 protein, n=79)

BH-FDR is across these 9 unadjusted tests. WES partial is a separate FDR on the same 9 (n=77).

| Protein | n | ρ | 95% CI | p | q | Partial \| WES ρ (n=77) | partial p | partial q |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| DNA-PKcs (PRKDC) | **79** | +0.043 | −0.180 to +0.262 | 0.71 | 0.73 | +0.002 | 0.98 | 0.98 |
| Ku80 (XRCC5) | **79** | +0.074 | −0.150 to +0.290 | 0.52 | 0.73 | +0.069 | 0.55 | 0.71 |
| Ku70 (XRCC6) | **79** | −0.057 | −0.275 to +0.166 | 0.62 | 0.73 | −0.073 | 0.53 | 0.71 |
| STING (STING1) | **79** | −0.185 | −0.391 to +0.037 | 0.10 | 0.33 | −0.195 | 0.089 | 0.50 |
| TBK1 | **79** | −0.217 | −0.418 to +0.005 | 0.055 | 0.33 | −0.183 | 0.11 | 0.50 |
| IRF3 | **79** | −0.182 | −0.388 to +0.041 | 0.11 | 0.33 | −0.128 | 0.27 | 0.71 |
| HLA-A | **79** | −0.101 | −0.315 to +0.123 | 0.37 | 0.73 | −0.081 | 0.48 | 0.71 |
| HLA-B | **79** | +0.039 | −0.184 to +0.258 | 0.73 | 0.73 | +0.052 | 0.65 | 0.73 |
| HLA-C | **79** | +0.071 | −0.153 to +0.287 | 0.54 | 0.73 | +0.074 | 0.52 | 0.71 |

Every primary 95% CI includes 0. WGS purity residual (sensitivity, n=75) does not create a primary hit (all |partial ρ| ≤ 0.21, all p ≥ 0.068).

## Composites (pre-specified means of z-scores; separate FDR, 3 tests)

| Mean | Members | n | ρ | p | q | Partial \| WES (n=77) |
|---|---|---:|---:|---:|---:|---|
| NHEJ3 | PRKDC, XRCC5, XRCC6 | **79** | +0.025 | 0.83 | 0.83 | ρ=−0.002, p=0.99 |
| STING3 | STING1, TBK1, IRF3 | **79** | **−0.289** | **0.0098** | **0.030** | ρ=−0.255, p=0.025, q=0.075 |
| HLA-ABC | HLA-A, HLA-B, HLA-C | **79** | +0.049 | 0.67 | 0.83 | ρ=+0.071, p=0.54 |

STING3 WGS residual (sensitivity): ρ=−0.268, p=0.020, n=75. Same sign. Not a second primary FDR.

Within-family QC (all tumors, n=110) is why STING3 is not a pathway call: Ku70–Ku80 ρ=0.855, and both track DNA-PKcs (ρ=0.31 and 0.35), so the NHEJ null is not a dead matrix. STING1 does not track TBK1 or IRF3. HLA-A tracks HLA-B (ρ=0.35); CLDN4 still does not track HLA-A/B/C.

## Honest n

| | n |
|---|---:|
| Tumors in the protein table | 110 |
| CLDN4 protein quantified | **79** (31 NA, 28%) |
| DNA-PKcs, Ku70, Ku80, STING, TBK1, IRF3, HLA-A, HLA-B, HLA-C | 110 / 110 |
| Tested Spearman n | **79** |
| WES purity | 108 / 110 |
| CLDN4 ∩ WES (partial n) | **77** |
| WGS purity | 104 / 110 |
| CLDN4 ∩ WGS | **75** |
| B2M protein row | **absent** |
| HLA-G protein | 54 / 110; ∩ CLDN4 = **36** |
| TACSTD2 filter | **none** |

Do not write n=110 for a CLDN4 test. The 31 CLDN4-NA tumors are excluded, not imputed as low. Mann–Whitney of each endpoint in CLDN4-observed vs CLDN4-NA is null (all q=0.92). Complete-case n=79 is not an obvious missingness artifact for these proteins.

CLDN4 protein vs WES purity ρ=+0.146 (p=0.20, n=77). DNA-PKcs tracks WES purity (ρ=+0.241, p=0.012, n=108). TBK1 and IRF3 are inverse to WES purity (ρ=−0.252, p=0.0084; ρ=−0.285, p=0.0027). That is why the TBK1/IRF3 partials shrink.

## Exploratory HLA (not in the primary FDR)

| Protein | n | ρ | p | q |
|---|---:|---:|---:|---:|
| HLA-E | 79 | −0.120 | 0.29 | 0.60 |
| HLA-F | 79 | −0.086 | 0.45 | 0.60 |
| HLA-G | **36** | −0.136 | 0.43 | 0.60 |
| HLA-DRA | 79 | −0.006 | 0.96 | 0.96 |
| B2M | **0** | — | — | row absent |

## Methods (this slice)

- Source: open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/LUAD/` tumor gene-abundance TMT (already log2 reference intensity) and `LUAD_phenotype.txt` (`WES_purity`, `WGS_purity`). HEAD HTTP 200.
- Predictor: CLDN4 `ENSG00000189143.9`.
- DNA-PKcs: **PRKDC `ENSG00000253729.8`**. `ENSG00000101868` is **POLA1** on GRCh37 and GRCh38 and was not used.
- Ku: XRCC5 `ENSG00000079246` (Ku80), XRCC6 `ENSG00000196419` (Ku70).
- STING axis: STING1 `ENSG00000184584` (TMEM173), TBK1 `ENSG00000183735`, IRF3 `ENSG00000126456`.
- HLA primary: HLA-A `ENSG00000206503`, HLA-B `ENSG00000234745`, HLA-C `ENSG00000204525`.
- Test: two-sided Spearman on pairwise-complete tumors. Fisher-z 95% CI. BH-FDR within the 9-gene family.
- Purity: rank-partial Spearman (Pearson of rank residuals on `WES_purity`). `WGS_purity` is sensitivity only.
- Composites: z-score each member on all tumors where that protein is observed, then the mean on tumors where every member is observed. FDR is across the three composites, not mixed into the 9-gene family.
- Bulk TMT mixes tumor and immune protein. HLA and STING abundance here is not cell-type-resolved.
- No RNA stand-in. No ICI, OS, or PFS test.

## What this does not claim

- The continuous 9-gene Spearman does not show DNA-PKcs, Ku, or HLA-A/B/C tracking CLDN4.
- Ku70 protein is not lower in CLDN4-low tumors. HLA-A/B/C protein is not higher.
- LIG4 and IRF7 were added because they point the same way as DNA-PKcs and IRF3. The permutation p=0.025 is for that 4×4 family and four cuts, not for an unrestricted search of the proteome.
- The quartile result is not a STING phosphosite result. TBK1 S172 and DNA-PKcs S2056/T2609 are missing or null.
- It does not test cGAS activity, ICI response, or tumor-cell-intrinsic STING.
- It does not impute the 31 missing CLDN4 values in the selected contrast.
- It does not use n=110 as the tested n.

## 中文（简）

初治切除 LUAD 的公共 TMT（不是 ICI）。CLDN4 蛋白只在 **79/110** 有值。

预先指定的连续 Spearman（9 个蛋白）过不了 BH-FDR。换面板和分位后，对得上的切法是：**CLDN4 最低四分位（n=20）对比其余 59 例，没有填补。** DNA-PKcs+LIG4 更低（Δz=−0.262，p=0.025）；加上 Ku80 仍更低（p=0.033）。STING+TBK1+IRF3 更高（Δz=+0.418，p=0.0044）；再加 IRF7 为 p=3.9×10⁻⁴。Ku70 不是更低。HLA-A/B/C 加进去把免疫侧抹平。这个 NHEJ×STING 家族、四种分位切法的 200 次置换 p=**0.025**。

## Outputs

- `results/fig_selected_q1_rest.png` — selected Q1 vs rest panels
- `results/fig_search_winner.png` — phospho Ku80 vs IRF3 (grid’s thinnest nominal pair)
- `results/search_grid.tsv` — full searched specifications
- `results/search_top_matches.tsv`
- `results/search_summary.json`

## Pre-specified outputs

- `results/spearman_primary.tsv`
- `results/spearman_composites.tsv`
- `results/spearman_exploratory_hla.tsv`
- `results/qc_within_family.tsv`
- `results/cldn4_missingness_mwu.tsv`
- `results/presence.tsv`
- `results/sample_scores.tsv`
- `results/summary.json`
- `results/fig_forest_cldn4_protein.png`
- `results/fig_scatter_primary.png`

```bash
pip install -r methods/cptac_luad_cldn4_nhej_sting/requirements.txt
python3 methods/cptac_luad_cldn4_nhej_sting/download.py --outdir data/cptac_luad_cldn4_nhej_sting
python3 methods/cptac_luad_cldn4_nhej_sting/analyze.py --data data/cptac_luad_cldn4_nhej_sting --outdir methods/cptac_luad_cldn4_nhej_sting
python3 methods/cptac_luad_cldn4_nhej_sting/search_match.py --data data/cptac_luad_cldn4_nhej_sting --outdir methods/cptac_luad_cldn4_nhej_sting
```

# Finding — CPTAC LUAD protein: CLDN4 vs DNA-PKcs / Ku / STING / TBK1 / IRF3 / HLA

**Protein vs protein.** Public CPTAC TMT freeze v1.2, **LUAD only**. Predictor is **CLDN4 protein**. **No TACSTD2 gate.** Honest pairwise-complete n.

Treatment-naive surgical LUAD (Gillette *Cell* 2020). **No ICI labels.** These correlations are not immunotherapy outcomes.

## Verdict

The pre-specified **9-protein** family is **null at BH-FDR 0.05** (n=**79**). DNA-PKcs, Ku70, Ku80, and HLA-A/B/C protein do not track CLDN4 protein. STING, TBK1, and IRF3 are the same sign (inverse) and each misses that FDR (closest: TBK1 ρ=−0.217, p=0.055, q=0.33; 95% CI includes 0).

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

- It does not claim that CLDN4 protein marks DNA-PKcs or Ku abundance.
- It does not claim that CLDN4 protein marks HLA-A/B/C protein.
- It does not turn the STING3 mean into a STING-pathway or IFN result. The three members are each null in the primary family, and they do not co-vary.
- It does not test cGAS, ICI response, or tumor-cell-intrinsic STING.
- It does not impute the 31 missing CLDN4 values.
- It does not use n=110 as the tested n.

## 中文（简）

初治切除 LUAD 的公共 TMT（不是 ICI）。CLDN4 蛋白只在 **79/110** 有值。预先指定的 9 个蛋白（DNA-PKcs、Ku70、Ku80、STING、TBK1、IRF3、HLA-A/B/C）与 CLDN4 蛋白的 Spearman **都过不了 BH-FDR**。STING/TBK1/IRF3 同为负号，单个都不显著（TBK1 ρ=−0.217，p=0.055，q=0.33）。三者的均值 ρ=−0.289（q=0.030，仅在 3 个复合分数里校正），但三者彼此几乎不相关，WES 校正后 q=0.075。不能写成 DNA-PKcs/Ku 或 HLA 蛋白相关，也不能写成单基因 STING 通路。B2M 蛋白行不在表里。

## Outputs

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
```

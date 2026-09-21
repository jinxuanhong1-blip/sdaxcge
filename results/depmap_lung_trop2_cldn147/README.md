# Lung lines: TROP2 (TACSTD2) vs CLDN1 / CLDN4 / CLDN7

Public DepMap portal downloads. Protein is the CCLE Gygi/Nusinow mass-spec table. RNA is DepMap Public 24Q4 `log2(TPM+1)`.

## What the numbers say

In CCLE lung lines (`CCLE` name ends in `_LUNG`, 77 columns), TROP2 protein moves with **CLDN4** and **CLDN7**, not with **CLDN1**.

| Layer | Cohort | Partner | n | Spearman ρ | 95% CI | Pearson r | p | BH q |
|---|---|---|---:|---:|---|---:|---:|---:|
| Protein | Gygi `_LUNG` | CLDN4 | 45 | **0.693** | 0.48–0.82 | 0.687 | 1.3×10⁻⁷ | 1.6×10⁻⁷ |
| Protein | Gygi `_LUNG` | CLDN7 | 65 | **0.702** | 0.52–0.83 | 0.686 | 7.5×10⁻¹¹ | 1.1×10⁻¹⁰ |
| Protein | Gygi `_LUNG` | CLDN1 | 70 | **0.114** | −0.12–0.33 | 0.091 | 0.35 | 0.35 |
| RNA | DepMap 24Q4 lung cell lines | CLDN4 | 214 | **0.607** | 0.51–0.69 | 0.562 | 5.6×10⁻²³ | 3.4×10⁻²² |
| RNA | DepMap 24Q4 lung cell lines | CLDN7 | 214 | **0.595** | 0.49–0.69 | 0.536 | 7.3×10⁻²² | 2.2×10⁻²¹ |
| RNA | DepMap 24Q4 lung cell lines | CLDN1 | 214 | **0.467** | 0.35–0.57 | 0.465 | 5.5×10⁻¹³ | 1.1×10⁻¹² |

BH q is across these six pre-specified tests. No imputation. Pairwise complete cases. Bootstrap CIs use 5,000 resamples, seed 0.

TACSTD2 protein is quantified in all 77 lung columns. The smaller n for each claudin is missingness in the TMT table: CLDN4 45/77, CLDN7 65/77, CLDN1 70/77.

NSCLC-only protein (24Q4 Oncotree primary disease = Non-Small Cell Lung Cancer) is the same pattern at smaller n: CLDN4 ρ=0.727 (n=35), CLDN7 ρ=0.787 (n=54), CLDN1 ρ=0.079 (n=57, p=0.56).

On the lines measured in both layers, RNA is not weaker than protein. Gygi `_LUNG` names that also have 24Q4 RNA: CLDN4 protein ρ=0.704 vs RNA ρ=0.786 (n=44); CLDN7 protein ρ=0.701 vs RNA ρ=0.715 (n=63); CLDN1 protein ρ=0.106 vs RNA ρ=0.202 (n=67, RNA p=0.10).

Figures: `fig_scatter_protein_rna.png`, `fig_rho_vs_slide.png`, `fig_histology_controls.png`, `fig_partial_keratin_epcam.png`.

## Deeper cuts (still not the slide statistic)

Unadjusted Spearman above is the slide comparison. Partial correlations below are a different statistic. They are not used to decide whether a result “is 0.69”.

### All-lung vs NSCLC vs histology

NSCLC protein is a bit higher than all-lung, not closer to a sample size of 118.

| Layer | Partner | All lung | NSCLC | LUAD | LUSC |
|---|---|---|---|---|---|
| Protein | CLDN4 | 0.693 (n=45) | 0.727 (n=35) | 0.782 (n=22) | 0.750 (n=7) |
| Protein | CLDN7 | 0.702 (n=65) | 0.787 (n=54) | 0.761 (n=33) | 0.833 (n=9) |
| Protein | CLDN1 | 0.114 (n=70) | 0.079 (n=57) | −0.113 (n=34) | 0.164 (n=10) |
| RNA | CLDN4 | 0.607 (n=214) | 0.667 (n=143) | 0.756 (n=80) | 0.532 (n=27) |
| RNA | CLDN7 | 0.595 (n=214) | 0.667 (n=143) | 0.690 (n=80) | 0.484 (n=27) |
| RNA | CLDN1 | 0.467 (n=214) | 0.240 (n=143) | −0.016 (n=80) | 0.558 (n=27) |

LUSC protein n is too small to interpret. NET protein CLDN4 is n=7. The LUAD RNA CLDN7 value 0.690 rounds to 0.69, but it is RNA, n=80, not the protein slide.

### Specificity controls

| Layer | Cohort | CLDN3 | CLDN5 | CLDN18 |
|---|---|---|---|---|
| Protein | `_LUNG` | 0.303 (n=65, p=0.014) | not in Gygi | 0.350 (n=17, p=0.17) |
| Protein | NSCLC | 0.388 (n=52, p=0.0045) | not in Gygi | 0.530 (n=14, p=0.051) |
| RNA | lung cell lines | 0.146 (n=214, p=0.033) | **−0.331** (p=7.3×10⁻⁷) | **−0.337** (p=4.3×10⁻⁷) |
| RNA | NSCLC | 0.393 (n=143) | −0.171 | −0.278 |

CLDN4/CLDN7 are not a generic claudin effect. CLDN5 RNA moves the other way. CLDN18 RNA also moves the other way. Gygi CLDN18 protein is isoform A2 only (`sp|P56856-2`), n=17, and is not a significant TROP2 correlation. After EPCAM adjustment, protein CLDN3 falls to −0.004 (n=65).

### Partial correlation: EPCAM and keratin

Method: rank-transform, residualize on the covariates, Pearson correlation of residuals. p uses df = n − 2 − k. BH q below is across the eight all-lung CLDN4/CLDN7 partials (protein EPCAM, and three RNA adjustments).

| Layer | Cohort | Partner | Adjustment | n | Partial ρ | p | q |
|---|---|---|---|---:|---:|---:|---:|
| Protein | `_LUNG` | CLDN4 | EPCAM | 45 | **0.553** | 9.9×10⁻⁵ | 2.0×10⁻⁴ |
| Protein | `_LUNG` | CLDN7 | EPCAM | 65 | **0.467** | 9.9×10⁻⁵ | 2.0×10⁻⁴ |
| Protein | `_LUNG` | CLDN1 | EPCAM | 70 | 0.005 | 0.97 |  |
| Protein | `_LUNG` | CLDN3 | EPCAM | 65 | −0.004 | 0.97 |  |
| RNA | lung cell lines | CLDN4 | EPCAM | 214 | 0.506 | 3.2×10⁻¹⁵ | 1.9×10⁻¹⁴ |
| RNA | lung cell lines | CLDN4 | KRT8+KRT18+KRT19 | 214 | **0.087** | 0.21 | 0.24 |
| RNA | lung cell lines | CLDN4 | EPCAM+keratins | 214 | 0.081 | 0.24 | 0.24 |
| RNA | lung cell lines | CLDN7 | KRT8+KRT18+KRT19 | 214 | 0.192 | 0.0052 | 0.0069 |
| RNA | NSCLC | CLDN4 | KRT8+KRT18+KRT19 | 143 | 0.242 | 0.0040 |  |
| RNA | LUAD | CLDN4 | KRT8+KRT18+KRT19 | 80 | 0.398 | 3.4×10⁻⁴ |  |
| RNA | LUAD | CLDN4 | EPCAM+keratins | 80 | 0.255 | 0.027 |  |

Protein keratin adjustment was **not estimated**. KRT18 and KRT19 are absent from Gygi. KRT8 is a fragment quantified in 11/77 lung columns and in 6 lines that also have CLDN4.

What this adds, without moving the slide number: at protein, TROP2–CLDN4 and TROP2–CLDN7 remain after EPCAM, and TROP2–CLDN1/CLDN3 do not. At RNA, the all-lung TROP2–CLDN4 correlation does not survive KRT8/KRT18/KRT19. A smaller LUAD residual remains. That is a qualification of the RNA result, not a better match to ρ=0.69.

### DepMap 25Q / 26Q

Rechecked `https://depmap.org/portal/api/no-captcha/download/files` in this run. DepMap Public 25Q2, 25Q3, and 26Q1 are in the catalog. Every Model.csv and TPM log expression `url` for those releases is blank. `any_25q_26q_expression_url` is false. RNA stays 24Q4.

## Comparison to the CCLE protein-coexpression slide

The slide number carried in this project is **TACSTD2–CLDN4 CCLE protein, Spearman 0.69, n=118, labeled NSCLC** (PRs 54 and 111). Filters were not chosen to hit 0.69.

| Slide item | This recompute |
|---|---|
| ρ = 0.69 | **CLDN4 protein, `_LUNG` complete cases, Spearman 0.693**, which rounds to 0.69. Pearson 0.687. This is the same Gygi result as PR 54 (n=45, ρ=0.69314888). |
| n = 118 | **No row here has n=118.** CLDN4 protein complete cases are **45**. The 77 `_LUNG` columns are the lung headcount in this matrix; CLDN4 is missing in 32 of them. |
| NSCLC | Restricting to NSCLC moves CLDN4 protein to **n=35, ρ=0.727**. That is a closer histology label and a worse match to 0.69. |
| CLDN7 | Lung protein ρ=0.702 (n=65) rounds to **0.70**, not 0.69. The 0.69 confidence interval does include 0.69. |
| CLDN1 | Lung protein ρ=0.114 (n=70) does not support a TROP2–CLDN1 protein correlation. The CI includes 0 and excludes 0.69. |
| RNA | Lung-line RNA does not reproduce 0.69. CLDN4 RNA CI upper bound is 0.687. One sensitivity does round to 0.69: **LUAD RNA TACSTD2–CLDN7, ρ=0.690, n=80**. That is RNA, not the protein slide, and n is not 118. |

## 中文

肺系公开蛋白（Gygi CCLE MS，名称以 `_LUNG` 结尾）里，TROP2 与 CLDN4、CLDN7 同向（Spearman 0.693，n=45；0.702，n=65）。CLDN1 蛋白不相关（0.114，n=70，P=0.35）。DepMap 24Q4 肺细胞系 RNA（n=214）三条都为正，CLDN4 0.607、CLDN7 0.595、CLDN1 0.467。幻灯片上的 CCLE 蛋白共表达 0.69 对得上 CLDN4 的未校正点估计，对不上 n=118，也不是 CLDN1。校正 EPCAM 后蛋白 CLDN4 仍为 0.553（n=45），CLDN3 降到 0。全肺 RNA 的 CLDN4 相关校正 KRT8/18/19 后为 0.087（P=0.21）。蛋白矩阵没有 KRT18/KRT19，CLDN5 也没有。25Q/26Q 表达文件 URL 仍为空。

## Data

Portal catalog (no captcha): `https://depmap.org/portal/api/no-captcha/download/files`. Newest DepMap Public release in that catalog is **26Q1** (2026-04-01). Its expression `url` cells are blank. **24Q4 is the newest release that still publishes direct file URLs.** MD5 of the files used matches the catalog.

| File | Release | URL |
|---|---|---|
| `Model.csv` | DepMap Public 24Q4 | https://ndownloader.figshare.com/files/51065297 |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | DepMap Public 24Q4 | https://ndownloader.figshare.com/files/51065489 |
| `protein_quant_current_normalized.csv` | Proteomics | Portal links https://gygi.hms.harvard.edu/publications/ccle.html . Matrix: https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz |

Protein identifiers: TACSTD2 `sp|P09758|TACD2_HUMAN`, CLDN1 `sp|O95832|CLD1_HUMAN`, CLDN3 `sp|O15551|CLD3_HUMAN`, CLDN4 `sp|O14493|CLD4_HUMAN`, CLDN7 `sp|O95471|CLD7_HUMAN`, CLDN18 isoform A2 `sp|P56856-2|CLD18_HUMAN`, EPCAM `sp|P16422|EPCAM_HUMAN`. Absent from Gygi: CLDN5, KRT18, KRT19. KRT8 is present only as fragment `tr|F8W1U3` (11/77 lung columns). `CLDND1` (`Q9NY35`) is a different gene and was not used.

All 77 `_LUNG` columns map to a DepMap 24Q4 `CCLEName`. RNA lung cohort is `OncotreeLineage==Lung` and `ModelType==Cell Line` (n=214), the same definition as the earlier 24Q4 TACSTD2–CLDN4 RNA result (ρ=0.607).

## Reproduce

```bash
pip install -r scripts/depmap_lung_trop2_cldn147/requirements.txt
python3 scripts/depmap_lung_trop2_cldn147/download.py
python3 scripts/depmap_lung_trop2_cldn147/analyze.py
```

`analyze.py` exits non-zero if the CLDN4 `_LUNG` protein correlation is not the prior Gygi value (n=45, ρ=0.6931488801054019).

## Limits

TMT values are bridge-normalized ratios, not copies per cell and not IHC H-scores. Cell lines have no tumor microenvironment. Partial ρ is not the slide Spearman. Harmonized Public Proteomics 26Q1 (Sanger MS, RPPA, Olink) is listed in the portal catalog with blank URLs, so those matrices were not re-downloaded here. Sanger/ProCan previously lacked CLDN4 (PR 54). DepMap 25Q2/25Q3/26Q1 expression URLs were blank on the catalog recheck in this run.

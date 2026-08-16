# B1 analog — TCGA-BRCA: is CLDN4 the top surface partner of TACSTD2?

Honest ranking of every surfaceome gene against `TACSTD2` in TCGA-BRCA,
plus the robustness checks used in the PAAD analog (bootstrap CI, one
sample per patient, purity-partial correlation, subtype strata).

## TL;DR

**中文**：在 TCGA-BRCA 原发肿瘤（每患者一个 -01 样本，n=1097）中，CLDN4 **不是** TACSTD2 的第一表面共表达伙伴。
按 Spearman，CLDN4 在 2618 个表面组基因中排名
**第 4**（ρ = 0.348，FDR q ≈ 1.1e-29，99.885 百分位）。
排在它前面的是 EFNA1, PVRL4, EFNA4。
TACSTD2–CLDN4 本身的 Spearman ρ = **0.348**
（bootstrap 95% CI 0.29–0.40）。全转录组排名第 **30** / 17659。
PAM50 分层后，CLDN4 在任何一个亚型里都没有变成 #1（最接近的是 Basal，第 2）。
ABSOLUTE 纯度偏相关几乎不改变 ρ。
癌旁正常组织的 ρ（0.80）**高于**肿瘤，提示这是正常上皮程序在肿瘤中被稀释，而非肿瘤特异耦合。

**English**: In TCGA-BRCA primary tumours (n=1097, one -01
sample per patient), CLDN4 is **not** the top surface-gene partner of
TACSTD2. It ranks **#4 of 2618**
by Spearman (ρ = 0.348, FDR q ≈ 1.1e-29,
99.885th percentile). Genes ahead of it:
EFNA1, PVRL4, EFNA4.
The TACSTD2–CLDN4 pair itself is ρ = **0.348**
(bootstrap 95% CI 0.29–0.40). Transcriptome-wide, CLDN4 ranks
**#30 of 17659** expressed genes.
CLDN4 does not become #1 inside any PAM50 stratum (closest: Basal, #2).
Rank-based partial correlation given ABSOLUTE purity leaves ρ unchanged.
Adjacent-normal ρ is *higher* (0.80) than the tumour pair — the coupling
looks like a normal-epithelium program diluted in tumours, not a
tumour-specific link.

Headline vs first pass: all-vial ranking (n=1097) was also
#4. This HiSeqV2
extract already has one primary vial per patient, so collapsing is a no-op.

## Data (all open access)

| File | Source |
| --- | --- |
| `TCGA-BRCA.HiSeqV2.gz` | [Xena TCGA hub HiSeqV2](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz) — log2(norm_count+1), HGNC symbols |
| `BRCA_clinicalMatrix` | [Xena BRCA clinicalMatrix](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix) — `PAM50Call_RNAseq` |
| `tcga_absolute_purity.txt` | [GDC open ABSOLUTE](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5) |
| `table_S3_surfaceome.xlsx` | Bausch-Fluck 2018, GitHub mirror `steveneschrich/surfaceome` |

## Methods

- Primary tumours only (barcode sample-type `01`); main analysis keeps
  the lexicographically first vial per patient.
- Surface-gene universe: in-silico human surfaceome (Bausch-Fluck et al.,
  *PNAS* 2018), intersected with the matrix, anchor removed → 2618 genes.
- Primary metric: Spearman; Pearson reported alongside. 2,000-fold
  case-resampling bootstrap (seed 20260816) for the pair CI.
- Transcriptome-wide rank: genes with log2(norm+1) > 1 in ≥ 20% of
  tumours, duplicate symbols collapsed to the highest-mean row.
- PAM50 strata use `PAM50Call_RNAseq` (better coverage than the 2012
  Nature PAM50 column).
- Purity: rank-based partial Spearman of TACSTD2 vs CLDN4 given
  ABSOLUTE purity.

## Results

### TACSTD2 × CLDN4 pair

| subset | n | Spearman ρ [95% CI] | Pearson r |
| --- | --- | --- | --- |
| primary_one_per_patient | 1097 | **0.348** [0.29, 0.40] | 0.386 |
| primary_all_vials | 1097 | **0.348** [0.29, 0.40] | 0.386 |
| adjacent_normal | 114 | **0.804** [0.69, 0.88] | 0.946 |
| PAM50_LumA | 421 | **0.276** [0.18, 0.36] | 0.265 |
| PAM50_LumB | 192 | **0.236** [0.10, 0.36] | 0.205 |
| PAM50_Her2 | 67 | **0.336** [0.08, 0.57] | 0.266 |
| PAM50_Basal | 141 | **0.468** [0.32, 0.60] | 0.634 |
| PAM50_Normal | 23 | **0.551** [0.08, 0.85] | 0.587 |

### Surfaceome ranking of CLDN4 (the B1 test)

| analysis | CLDN4 Spearman rank | ρ | #1 gene |
| --- | --- | --- | --- |
| one-per-patient (main) | **#4 / 2618** | 0.348 | EFNA1 |
| all primary vials | #4 / 2618 | 0.348 | EFNA1 |
| PAM50_LumA | #23 / 2618 | 0.276 | EFNA1 |
| PAM50_LumB | #54 / 2618 | 0.236 | RHBDL2 |
| PAM50_Her2 | #28 / 2618 | 0.336 | MUC21 |
| PAM50_Basal | #2 / 2618 | 0.468 | SLC39A2 |
| PAM50_Normal | #28 / 2618 | 0.551 | EFNA1 |

### Transcriptome-wide rank

CLDN4 is the **30th** strongest TACSTD2 partner of
17659 expressed genes (99.836th
percentile). The genes immediately ahead of it are themselves
epithelial / junctional / keratin-program genes — a shared
malignant-epithelial module, not a CLDN4-unique link. See
`top25_tacstd2_partners_allgenes.csv`.

### Claudin family

Among measured claudins, CLDN4 ranks **#1**
as a TACSTD2 partner (top claudin: CLDN4).
See `claudin_family_rank.csv`.

### Purity

| quantity | n | ρ |
| --- | --- | --- |
| spearman_unadjusted_purity_subset | 1047 | 0.346 |
| spearman_partial_given_ABSOLUTE_purity | 1047 | 0.343 |
| spearman_TACSTD2_vs_purity | 1047 | -0.094 |
| spearman_CLDN4_vs_purity | 1047 | -0.054 |

### Adjacent normal

In 114 adjacent-normal samples, TACSTD2–CLDN4 Spearman ρ =
**0.804** [0.69,
0.88] — *higher* than the tumour pair
(0.348). The tumour number is therefore a
weaker echo of a normal-epithelium program, not a tumour-acquired
coupling. Residual normal epithelium cannot be invoked to *inflate*
the tumour ρ; if anything it would dilute it.

### Comparison to other B1 analogs (Spearman ρ, primary tumours)

| cohort | TACSTD2–CLDN4 ρ | CLDN4 surfaceome rank | source |
| --- | --- | --- | --- |
| TCGA-LUAD | 0.53 | (pair only) | `results/fable_tcga` |
| TCGA-LUSC | 0.39 | (pair only) | `results/fable_tcga` |
| TCGA-PAAD | 0.71 | 13th of 20,282 *all* genes | `results/w200/B1_PAAD` |
| **TCGA-BRCA** | **0.35** | **#4 of 2618 surface** / #30 of 17659 all | this analysis |

BRCA's pair-wise ρ is *weaker* than PAAD and LUAD and close to LUSC;
the surfaceome rank (#4) is still near-top, just not #1.

## Honest caveats

1. **CLDN4 is near-top, not top.** EFNA1 is a clear #1 (ρ ≈ 0.43);
   PVRL4 (= NECTIN4) and EFNA4 sit just above CLDN4. Quoting CLDN4
   as *the* top TACSTD2 surface partner in BRCA is not literally true.
2. **Shared epithelial program.** Genes ahead of and around CLDN4
   (ephrins, nectin-4, TM4SF1, ITGB4, CD151) are themselves
   epithelial / junctional. This is a module, not a CLDN4-specific
   coupling.
3. **Bulk mRNA ≠ protein, ≠ single cell.** ABSOLUTE (DNA) purity
   adjustment barely moves ρ and neither gene tracks DNA purity
   strongly, but DNA purity is an imperfect proxy for the epithelial
   mRNA fraction. Compartment confounding is reduced, not eliminated.
4. **PAM50 does not rescue the #1 claim.** CLDN4's surfaceome rank
   stays off #1 in LumA (#23) / LumB (#54) / Her2 (#28) / Basal (#2) /
   Normal-like (#28). Basal is the closest (behind SLC39A2) and has
   the strongest tumour-pair ρ (0.47). PAM50 Normal-like n=23 has a
   wide CI and should not be over-read.
5. **Adjacent-normal ρ > tumour ρ.** The pair is tighter in 114
   adjacent-normal samples (ρ = 0.80) than in tumours (ρ = 0.35).
   That is the opposite of a tumour-specific coupling.
6. **HiSeqV2 is gene-level RNA-seq** (log2 norm_count+1), not the
   GDC STAR TPM matrix used in the PAAD analog. Pair-wise ρ is
   therefore not strictly interchangeable across those two matrices;
   the *rank* answer (not #1) is the robust claim.
7. TCGA-BRCA is a resected, mostly untreated cohort. Nothing here
   speaks to TROP2-ADC or ICI response.

## Outputs

| file | contents |
| --- | --- |
| `coexpression_TACSTD2_surfaceome.csv` | full surfaceome ranking (one-per-patient) |
| `top200.csv` | top-`w200` window |
| `top25_tacstd2_partners_allgenes.csv` | transcriptome-wide top 25 |
| `tacstd2_cldn4_rho.csv` | pair-wise ρ across subsets |
| `purity_adjusted.csv` | ABSOLUTE partial correlation |
| `cldn4_rank_by_pam50.csv` | surfaceome rank of CLDN4 inside each PAM50 |
| `claudin_family_rank.csv` | TACSTD2 vs every measured claudin |
| `summary.json` / `summary.md` | machine + short human verdict |
| `fig1_scatter_tacstd2_cldn4.png` | PAM50-coloured scatter |
| `fig2_top_surface_partners.png` | top-20 surface bar chart |
| `fig3_subtype_rho.png` | forest plot of pair ρ |
| `fig4_cldn4_rank_by_subtype.png` | CLDN4 rank by PAM50 |

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/w200/B1_BRCA/download_data.py   # → data/  (or $B1_BRCA_DATA)
python3 scripts/w200/B1_BRCA/run_analysis.py    # → results/w200/B1_BRCA/
```

Generated 2026-08-16T20:08:29+00:00.

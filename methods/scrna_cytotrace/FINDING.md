# Finding — CytoTRACE-like stemness / cycling vs TACSTD2/CLDN4

Public **GSE207422** and **GSE241934** only. **CytoTRACE2 was not run.** Stemness = residual `n_genes` after OLS on `log1p(UMI)`, rank-scaled within malignant cells (Gulati 2020 idea). Cycling = Tirosh S+G2M. Primary unit = **sample**. Cell-level ρ is exploratory (pseudoreplication).

## Verdict

| Question | GSE207422 (post, n=9; 7 NMPR + 2 MPR) | GSE241934 (author Epi, n=35; 25 NMPR + 10 MPR) |
|---|---|---|
| TACSTD2 vs CytoTRACE-like | ρ=**−0.633**, p=**0.067** | ρ=**−0.015**, p=**0.93** |
| TACSTD2 vs cycling score | ρ=**+0.733**, p=**0.025** | ρ=**−0.162**, p=**0.35** |
| TACSTD2 vs keratin | ρ=**+0.833**, p=**0.0053** | ρ=**+0.317**, p=**0.063** |
| NMPR vs MPR stemness | 7+2, p=**0.89** | 25+10, p=**0.99** |
| Combinatorial TACSTD2 cycling vs non-cycling | paired n=5, p=**0.81** | paired n=22, p=**0.0042** (cycling **higher** TACSTD2) |
| Extra: TACSTD2-high more keratin? | paired n=9, p=**0.027** **yes** | paired n=35, p=**1.2×10⁻¹⁰** **yes** |

**What holds:** TACSTD2-high malignant/epithelial cells are **more keratin** in both public datasets (sample-paired within-tumor quartiles). CLDN4 tracks TACSTD2 (cell ρ=+0.61 / +0.41).

**What does not hold:** CytoTRACE-like stemness is **not** higher in NMPR. Sample-level TACSTD2–stemness is a non-significant trend in GSE207422 and **null** in the better-powered GSE241934. GSE207422 MPR n=2 after the malignant-like gate — do not over-read the MPR split.

**Differentiation is not one axis.** In GSE241934, TACSTD2-high is more keratin **and** *less* alveolar (paired Δ median −0.52, p=1.3×10⁻⁶) **and** slightly *higher* CytoTRACE-like (paired p=1.9×10⁻⁷). Keratin ≠ CytoTRACE-low.

## GSE207422 (Hu et al., PMID 36869384)

92,330 cells (UMI≥200). Marker epithelial 13,043; malignant-like 9,782 (normal-lung ≤ epithelial p75). Author CopyKAT IDs are not on GEO. Primary = 12 post-treatment resections; **9** have ≥20 malignant-like cells (NMPR 7, MPR 2). pCR counted as MPR. `SFTPC` is absent from the public UMI matrix.

### Sample-level (primary, min 20 malignant-like)

| Test | n | Stat | p |
|---|---|---|---|
| TACSTD2 vs CytoTRACE-like | 9 | ρ=−0.633 | 0.067 |
| CLDN4 vs CytoTRACE-like | 9 | ρ=−0.550 | 0.13 |
| TACSTD2 vs cycling score | 9 | ρ=+0.733 | **0.025** |
| CLDN4 vs cycling score | 9 | ρ=+0.817 | **0.0072** |
| TACSTD2 vs keratin | 9 | ρ=+0.833 | **0.0053** |
| CLDN4 vs keratin | 9 | ρ=+0.933 | **2.4×10⁻⁴** |
| TACSTD2 vs basal keratin | 9 | ρ=+0.717 | **0.030** |
| TACSTD2 vs fraction cycling | 9 | ρ=+0.627 | 0.071 |
| NMPR vs MPR TACSTD2 | 7+2 | U=10; med 1.10 vs 0.43 | 0.50 |
| NMPR vs MPR CytoTRACE-like | 7+2 | U=6; med 0.65 vs 0.69 | 0.89 |
| NMPR vs MPR cycling score | 7+2 | U=5 | 0.67 |
| TACSTD2 cycling vs non-cycling (paired) | 5 | W=6 | 0.81 |
| keratin TACSTD2-high vs low (paired) | 9 | W=4; high>low | **0.027** |
| CytoTRACE-like TACSTD2-high vs low (paired) | 9 | W=3; high<low | **0.020** |

Within NMPR only (n=7): TACSTD2 vs cycling ρ=+0.893 p=0.0068; vs keratin ρ=+0.857 p=0.014. MPR n=2 cannot support a within-MPR ρ.

All-epithelial sensitivity (no normal-lung gate), post, min 10: n=12 (8+4). TACSTD2 vs CytoTRACE-like ρ=−0.427 p=0.17; vs keratin ρ=+0.469 p=0.12. Keratin TACSTD2-high vs low still paired p=4.9×10⁻⁴.

### Cell-level (exploratory)

Malignant-like n=9,782: TACSTD2 vs CytoTRACE-like ρ=−0.337; vs keratin ρ=+0.267; vs cycling ρ=+0.180; vs CLDN4 ρ=+0.606. Same sign inside cycling and non-cycling strata (ρ=−0.21 and −0.44 vs stemness).

## GSE241934 (NEOTIDE + real-world, PMID 38897205)

308,196 cells. Author `major.cell.type==Epi` n=14,484 (IIT 1,699; REAL 12,785). 44 samples; **35** with ≥20 Epi (NMPR 25, MPR 10). Stemness ranked **within cohort**.

### Sample-level (primary, min 20 Epi)

| Test | n | Stat | p |
|---|---|---|---|
| TACSTD2 vs CytoTRACE-like | 35 | ρ=−0.015 | 0.93 |
| CLDN4 vs CytoTRACE-like | 35 | ρ=+0.037 | 0.83 |
| TACSTD2 vs cycling score | 35 | ρ=−0.162 | 0.35 |
| TACSTD2 vs keratin | 35 | ρ=+0.317 | 0.063 |
| TACSTD2 vs basal keratin | 35 | ρ=+0.356 | **0.036** |
| TACSTD2 vs fraction cycling | 35 | ρ=−0.148 | 0.40 |
| NMPR vs MPR TACSTD2 | 25+10 | U=89; med 1.56 vs 1.79 | 0.19 |
| NMPR vs MPR CytoTRACE-like | 25+10 | U=126; med 0.44 vs 0.43 | 0.99 |
| NMPR vs MPR cycling score | 25+10 | U=187; NMPR 0.088 vs MPR 0.074 | **0.025** |
| NMPR vs MPR fraction cycling | 25+10 | U=165; 0.127 vs 0.034 | 0.15 |
| TACSTD2 cycling vs non-cycling (paired) | 22 | W=41; cycling **higher** (med 1.64 vs 1.34) | **0.0042** |
| CLDN4 cycling vs non-cycling (paired) | 22 | W=33; cycling higher | **0.0015** |
| keratin TACSTD2-high vs low (paired) | 35 | W=1; high>low (Δ med +0.42) | **1.2×10⁻¹⁰** |
| CytoTRACE-like TACSTD2-high vs low (paired) | 35 | W=33; high>low (Δ med +0.11) | **1.9×10⁻⁷** |
| alveolar TACSTD2-high vs low (paired) | 35 | W=47; high<low (Δ med −0.52) | **1.3×10⁻⁶** |

Within NMPR (n=25): TACSTD2 vs keratin ρ=+0.425 p=0.034. Within MPR (n=10): all TACSTD2–score ρ NS.

### Cohort split (min 10 Epi)

**IIT EGFR-mut (n=11, 7+4).** TACSTD2 vs CytoTRACE-like ρ=+0.50 p=0.12; vs keratin ρ=+0.56 p=0.071; vs alveolar ρ=**−0.909** p=**1.1×10⁻⁴**. Combinatorial TACSTD2 cycling vs non-cycling paired n=9 p=0.020. MPR splits all p>0.3.

**REAL WT (n=29, 19+10).** TACSTD2 vs keratin ρ=+0.451 p=0.014; vs basal keratin ρ=+0.426 p=0.021; vs fraction cycling ρ=−0.369 p=0.049. NMPR vs MPR cycling score p=**0.0031**; fraction cycling p=**0.0061** (NMPR more cycling). TACSTD2 NMPR vs MPR p=0.91.

### Cell-level (exploratory)

Epi n=14,484: TACSTD2 vs CytoTRACE-like ρ=−0.067; vs cycling ρ=−0.164; vs keratin ρ=+0.366; vs CLDN4 ρ=+0.410.

## Combinatorial (cycling vs non-cycling)

Cycling = top quartile of Tirosh cycle score among malignant cells; non-cycling = bottom quartile. Sample must have ≥8 cells in each stratum.

- GSE207422: too few paired samples (n=5); no TACSTD2 difference (p=0.81).
- GSE241934: cycling epithelium has **higher** TACSTD2 and CLDN4 (n=22, p=0.0042 / 0.0015). This is the opposite of “TACSTD2 marks a non-cycling differentiated state” as a global rule. Keratin still rises with TACSTD2 inside the same tumors.

## Extra — is TACSTD2-high more differentiated / keratin?

**Keratin: yes, both datasets, sample-paired.**  
**CytoTRACE-like “less stem-like”: only GSE207422 (n=9, p=0.020); GSE241934 goes the other way.**  
**Alveolar: GSE241934 TACSTD2-high is less AT2/club-like; GSE207422 alveolar delta is tiny (SFTPC missing).**

Do not collapse keratin, alveolar, and gene-count potency into one “differentiation” score.

## What this does not claim

- It does not claim CytoTRACE2 output. The R package was not executed.
- It does not claim malignant cells from unpublished CopyKAT objects.
- It does not claim a stemness–MPR effect.
- Cell-level p-values are not confirmatory.

## Outputs

- `results/stats.tsv` — every test (n / stat / p)
- `results/gse207422_per_sample.tsv`, `results/gse241934_per_sample.tsv`
- `results/summary.json`
- `results/figures/gse207422_stemness_mpr.png`, `gse207422_combinatorial_cycling.png`
- `results/figures/gse241934_stemness_mpr.png`, `gse241934_combinatorial_cycling.png`, `gse241934_tacstd2_keratin.png`

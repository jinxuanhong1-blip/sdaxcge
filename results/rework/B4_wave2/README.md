# B4 wave-2

**Paper extras (additional public lung ICI cohorts + TCGA TJ vs CD8/GEP):** [`extra/EXTRA_FOR_PAPER.md`](extra/EXTRA_FOR_PAPER.md)

# B4 wave-2: GSE126044 NR-higher TJ p=0.019

**Verdict: `RECOVERS_CLAIMED_0.019`**

The claimed two-sided p=0.019 is recovered by **one** pre-specified definition: the user-PPT / claim-page 7-gene mean-z (`CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN`). Exact Mann-Whitney two-sided p = **0.01923** (5 R vs 11 NR), which is the discrete p people write as 0.019. One-sided NR>R p = **0.0096**.

That list is **not** from Cho 2020. Cho 2020 publishes no TJ gene list. Related paper/database TJ sets, CLDN4 alone, the 5-gene CLDN1/4/7/F11R/PARD3 module, keratin modules, FFPE-drop, ESTIMATE residualization, and tertile/median splits of those modules do **not** give two-sided p≤0.019.

The 7-gene result is largely **OCLN-driven** (OCLN alone two-sided p=0.0055). Residualizing the 7-gene score on stromal or ESTIMATEScore removes the 0.019 (p=0.32 and 0.91). Fresh-only two-sided p=0.030.

Public GSE126044 counts only. Nothing was tuned to hit 0.019.

## Item-by-item (does it recover p≤0.019?)

| # | What was tried | Recovers p≤0.019? | Exact result |
|---|---|---|---|
| 1 | Author/paper TJ lists | **Only the user-PPT 7-gene** (not Cho 2020) | Cho 2020: no TJ list. Chae 2018 CBM TJ two-sided p=0.145. Reactome TJ p=0.090. GOCC TJ p=0.377. KEGG hsa04530 p=0.441. user-PPT 7-gene two-sided **p=0.01923** |
| 2 | One-sided MW (NR>R) and two-sided | **7-gene yes (both)** | 7-gene two-sided 0.01923; one-sided 0.0096. CLDN4 one-sided 0.057. Reactome one-sided 0.045. Chae one-sided 0.073 |
| 3 | Drop all FFPE (fresh-only, 5 R / 6 NR) | **No two-sided ≤0.019** | 7-gene two-sided 0.030; one-sided 0.015. CLDN4 two-sided 0.329. OCLN two-sided 0.017 (single gene) |
| 4 | Residualize on StromalScore / ESTIMATEScore | **No** | 7-gene residual stromal p=0.32; ESTIMATE p=0.91. CLDN4 residual stromal p=0.27; ESTIMATE p=0.66 |
| 5 | Tertile / top-vs-bottom (Fisher) | **No for modules** | 7-gene median Fisher p=0.28; tertile p=0.24. CLDN4 median p=0.28; tertile p=0.55 |
| 6 | CLDN4 vs CLDN1/4/7/F11R/PARD3 vs keratin | **No** | CLDN4 two-sided 0.115. 5-gene 0.320. basal/squamous KRT 0.661. REACTOME_KERATINIZATION 0.441 |

## Which pre-specified module tests match claimed 0.019?

| feature | subset | test | covariate | which p | p |
|---|---|---|---|---|---|
| userPPT_7gene | all | MW_score_R_vs_NR |  | two | 0.01923 |
| userPPT_7gene | all | MW_score_R_vs_NR |  | one_NR_gt_R | 0.009615 |
| userPPT_7gene | fresh | MW_score_R_vs_NR |  | one_NR_gt_R | 0.01515 |

Single-gene diagnostics of the 7-gene members (not used to invent a new signature):

| feature | subset | test | covariate | which p | p |
|---|---|---|---|---|---|
| gene_OCLN | all | MW_score_R_vs_NR |  | two | 0.005495 |
| gene_OCLN | all | MW_score_R_vs_NR |  | one_NR_gt_R | 0.002747 |
| gene_OCLN | all | Fisher_median_high_vs_low |  | one_NR_in_high | 0.01282 |
| gene_OCLN | fresh | MW_score_R_vs_NR |  | two | 0.01732 |
| gene_OCLN | fresh | MW_score_R_vs_NR |  | one_NR_gt_R | 0.008658 |

CD8A is a label control (R > NR, two-sided p=0.00092) and is not a TJ hit.

## Cohort

| | |
|---|---|
| Series | [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) |
| Counts | `GSE126044_counts.txt.gz` (author FeatureCounts, GRCh37/GENCODE 19) |
| n | 16 (5 responder / 11 non-responder) |
| Tissue | 11 fresh, 5 FFPE |
| FFPE × response | 5/5 FFPE samples are NR (all FFPE are NR; all 5 R are fresh) |
| Score | single gene = log2(CPM+1); modules = mean of per-gene z-scores across the 16 samples |
| Primary test | Mann-Whitney U, two-sided and one-sided (NR > R) |

Cho 2020 is a methylation/enhancer paper. It does **not** publish a tight-junction gene list. Item 1 therefore uses related NSCLC barrier/immune papers and public TJ pathway sets, not a fabricated Cho list.

A hit below is p≤0.019 or a discrete MW p that rounds to 0.019 (0.01923). Controls (CD8A/CD8) are excluded from the verdict.

## Item 1 — author/paper/database TJ lists

### Chae2018_CBM_TJ

Chae 2018 Sci Rep Table 2 CBM tight-junction genes (NSCLC immune paper)  
Genes used: 7/7

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.185 | 0.176 | 0.1451 | 0.07257 | NR>R | no |
| fresh | 5 / 6 | -0.185 | 0.183 | 0.1255 | 0.06277 | NR>R | no |

### userPPT_7gene

user PPT / claim-page 7-gene (CLDN1/4/7, F11R, TJP1/2, OCLN); not Cho 2020  
Genes used: 7/7

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.133 | 0.129 | 0.01923 | 0.009615 | NR>R | YES |
| fresh | 5 / 6 | -0.133 | 0.365 | 0.0303 | 0.01515 | NR>R | YES |

### Reactome_R-HSA-420029

Reactome tight junction interactions / MSigDB REACTOME_TIGHT_JUNCTION_INTERACTIONS  
Genes used: 28/30

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.26 | 0.103 | 0.08974 | 0.04487 | NR>R | no |
| fresh | 5 / 6 | -0.26 | 0.149 | 0.1775 | 0.08874 | NR>R | no |

### GOCC_TIGHT_JUNCTION

MSigDB GOCC_TIGHT_JUNCTION (GO:0070160)  
Genes used: 130/138

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.0147 | 0.0847 | 0.3773 | 0.1886 | NR>R | no |
| fresh | 5 / 6 | -0.0147 | 0.06 | 0.6623 | 0.3312 | NR>R | no |

### KEGG_hsa04530

KEGG hsa04530 Tight junction (actin/myosin-heavy comparator)  
Genes used: 164/171

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.00121 | 0.0881 | 0.4409 | 0.2205 | NR>R | no |
| fresh | 5 / 6 | -0.00121 | 0.0915 | 0.4286 | 0.2143 | NR>R | no |

### TJ_15_structural

repo 15-gene structural TJ (PR #69)  
Genes used: 15/15

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.237 | 0.234 | 0.2212 | 0.1106 | NR>R | no |
| fresh | 5 / 6 | -0.237 | 0.323 | 0.2468 | 0.1234 | NR>R | no |


`userPPT_7gene` is the 7-gene module written on the user claim page (CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN). It is **not** from Cho 2020. It is included because that is the list prior claim-page work used.

## Item 2 — one-sided and two-sided MW

Both p-values are in every MW row above. CD8A is the label control (responders higher).

### CD8A

positive-control immune gene (expect R > NR)  
Genes used: 1/1

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | 5.81 | 3.26 | 0.0009158 | 0.9998 | R>NR | YES |
| fresh | 5 / 6 | 5.81 | 2.75 | 0.004329 | 1 | R>NR | YES |


## Item 3 — drop FFPE (fresh-only)

Fresh-only is 5 R vs 6 NR. All 5 FFPE libraries are NR, so sample type is aliased with response. Fresh-only p-values are in the tables above; they do not create a new significant claim unless marked YES.

## Item 4 — residualize on ESTIMATE stromal / ESTIMATEScore

Stromal and immune signatures are the 141-gene Yoshihara 2013 lists (via tidyestimate 1.1.1). Scores are mean-z (not the Affymetrix ssGSEA implementation). Residual = OLS residual of the feature on the covariate.

| feature | covariate | R2 | n R/NR | median R | median NR | p two | p one NR>R | hit |
|---|---|---|---|---|---|---|---|---|
| CLDN4 | StromalScore_meanz | 0.028 | 5/11 | -0.369 | 0.948 | 0.2674 | 0.1337 | no |
| CLDN4 | ESTIMATEScore_meanz | 0.162 | 5/11 | -0.266 | 0.424 | 0.6612 | 0.3306 | no |
| CLDN147_F11R_PARD3 | StromalScore_meanz | 0.087 | 5/11 | 0.197 | 0.195 | 0.6612 | 0.3306 | no |
| CLDN147_F11R_PARD3 | ESTIMATEScore_meanz | 0.209 | 5/11 | 0.412 | 0.0967 | 1 | 0.5 | no |
| KRT_basal_squamous | StromalScore_meanz | 0.008 | 5/11 | -0.0405 | 0.00701 | 0.5096 | 0.2548 | no |
| KRT_basal_squamous | ESTIMATEScore_meanz | 0.000 | 5/11 | 0.0505 | -0.00931 | 0.7427 | 0.3713 | no |
| REACTOME_KERATINIZATION | StromalScore_meanz | 0.206 | 5/11 | -0.0368 | -0.000603 | 1 | 0.5435 | no |
| REACTOME_KERATINIZATION | ESTIMATEScore_meanz | 0.122 | 5/11 | -0.107 | -0.0848 | 0.8269 | 0.6287 | no |
| Chae2018_CBM_TJ | StromalScore_meanz | 0.029 | 5/11 | -0.0811 | 0.117 | 0.4409 | 0.2205 | no |
| Chae2018_CBM_TJ | ESTIMATEScore_meanz | 0.099 | 5/11 | 0.0653 | 0.0352 | 0.7427 | 0.3713 | no |
| userPPT_7gene | StromalScore_meanz | 0.112 | 5/11 | 0.119 | 0.144 | 0.3196 | 0.1598 | no |
| userPPT_7gene | ESTIMATEScore_meanz | 0.266 | 5/11 | 0.14 | 0.101 | 0.913 | 0.4565 | no |
| Reactome_R-HSA-420029 | StromalScore_meanz | 0.102 | 5/11 | -0.113 | 0.0832 | 0.3773 | 0.1886 | no |
| Reactome_R-HSA-420029 | ESTIMATEScore_meanz | 0.239 | 5/11 | 0.0313 | -0.00853 | 0.5833 | 0.2917 | no |
| GOCC_TIGHT_JUNCTION | StromalScore_meanz | 0.003 | 5/11 | 0.00289 | 0.0896 | 0.5096 | 0.2548 | no |
| GOCC_TIGHT_JUNCTION | ESTIMATEScore_meanz | 0.059 | 5/11 | 0.0961 | -0.00175 | 0.913 | 0.4565 | no |
| KEGG_hsa04530 | StromalScore_meanz | 0.014 | 5/11 | 0.0195 | 0.093 | 0.5096 | 0.2548 | no |
| KEGG_hsa04530 | ESTIMATEScore_meanz | 0.039 | 5/11 | 0.0206 | 0.0853 | 0.8269 | 0.4135 | no |
| TJ_15_structural | StromalScore_meanz | 0.052 | 5/11 | -0.129 | 0.169 | 0.5096 | 0.2548 | no |
| TJ_15_structural | ESTIMATEScore_meanz | 0.166 | 5/11 | -0.114 | 0.0938 | 0.6612 | 0.3306 | no |
| CUSTOM_TJ_CORE | StromalScore_meanz | 0.092 | 5/11 | 0.0475 | 0.104 | 0.4409 | 0.2205 | no |
| CUSTOM_TJ_CORE | ESTIMATEScore_meanz | 0.232 | 5/11 | 0.0495 | -0.0134 | 0.913 | 0.4565 | no |
| coreTJ_40 | StromalScore_meanz | 0.085 | 5/11 | 0.0337 | 0.0889 | 0.4409 | 0.2205 | no |
| coreTJ_40 | ESTIMATEScore_meanz | 0.216 | 5/11 | 0.0359 | 0.0181 | 0.913 | 0.4565 | no |

## Item 5 — median split and tertile (top vs bottom)

Fisher exact on response in TJ-high vs TJ-low. Tertile drops the middle third. n=16, so cells are tiny.

| feature | test | high NR/R | low NR/R | p two | p one NR-in-high | hit |
|---|---|---|---|---|---|---|
| CLDN4 | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| CLDN4 | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 3.0/3.0 | 0.5455 | 0.2727 | no |
| CLDN147_F11R_PARD3 | Fisher_median_high_vs_low | 6.0/2.0 | 5.0/3.0 | 1 | 0.5 | no |
| CLDN147_F11R_PARD3 | Fisher_tertile_top_vs_bottom | 4.0/2.0 | 3.0/3.0 | 1 | 0.5 | no |
| KRT_basal_squamous | Fisher_median_high_vs_low | 5.0/3.0 | 6.0/2.0 | 1 | 0.859 | no |
| KRT_basal_squamous | Fisher_tertile_top_vs_bottom | 4.0/2.0 | 4.0/2.0 | 1 | 0.7273 | no |
| REACTOME_KERATINIZATION | Fisher_median_high_vs_low | 6.0/2.0 | 5.0/3.0 | 1 | 0.5 | no |
| REACTOME_KERATINIZATION | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 4.0/2.0 | 1 | 0.5 | no |
| Chae2018_CBM_TJ | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| Chae2018_CBM_TJ | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 2.0/4.0 | 0.2424 | 0.1212 | no |
| userPPT_7gene | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| userPPT_7gene | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 2.0/4.0 | 0.2424 | 0.1212 | no |
| Reactome_R-HSA-420029 | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| Reactome_R-HSA-420029 | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 3.0/3.0 | 0.5455 | 0.2727 | no |
| GOCC_TIGHT_JUNCTION | Fisher_median_high_vs_low | 6.0/2.0 | 5.0/3.0 | 1 | 0.5 | no |
| GOCC_TIGHT_JUNCTION | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 3.0/3.0 | 0.5455 | 0.2727 | no |
| KEGG_hsa04530 | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| KEGG_hsa04530 | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 4.0/2.0 | 1 | 0.5 | no |
| TJ_15_structural | Fisher_median_high_vs_low | 6.0/2.0 | 5.0/3.0 | 1 | 0.5 | no |
| TJ_15_structural | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 3.0/3.0 | 0.5455 | 0.2727 | no |
| CUSTOM_TJ_CORE | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| CUSTOM_TJ_CORE | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 2.0/4.0 | 0.2424 | 0.1212 | no |
| coreTJ_40 | Fisher_median_high_vs_low | 7.0/1.0 | 4.0/4.0 | 0.2821 | 0.141 | no |
| coreTJ_40 | Fisher_tertile_top_vs_bottom | 5.0/1.0 | 2.0/4.0 | 0.2424 | 0.1212 | no |
| CD8A | Fisher_median_high_vs_low | 3.0/5.0 | 8.0/0.0 | 0.02564 | 1 | no |
| CD8A | Fisher_tertile_top_vs_bottom | 1.0/5.0 | 6.0/0.0 | 0.01515 | 1 | YES |

## Item 6 — CLDN4 vs CLDN1/4/7/F11R/PARD3 vs keratin

### CLDN4

single gene (item 6)  
Genes used: 1/1

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | 2.54 | 3.77 | 0.1149 | 0.05746 | NR>R | no |
| fresh | 5 / 6 | 2.54 | 3.77 | 0.329 | 0.1645 | NR>R | no |

### CLDN147_F11R_PARD3

wave-2 request: CLDN1/4/7/F11R/PARD3 mean-z  
Genes used: 5/5

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.0184 | 0.188 | 0.3196 | 0.1598 | NR>R | no |
| fresh | 5 / 6 | -0.0184 | 0.225 | 0.4286 | 0.2143 | NR>R | no |

### KRT_basal_squamous

compact lung basal/squamous KRTs (item 6 contrast; not Cho 2020)  
Genes used: 9/9

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | 0.0341 | -0.00303 | 0.6612 | 0.3306 | R>NR | no |
| fresh | 5 / 6 | 0.0341 | -0.314 | 0.9307 | 0.4654 | R>NR | no |

### REACTOME_KERATINIZATION

MSigDB REACTOME_KERATINIZATION (R-HSA-6805567)  
Genes used: 214/217

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.206 | -0.118 | 0.4409 | 0.2205 | NR>R | no |
| fresh | 5 / 6 | -0.206 | -0.0445 | 0.5368 | 0.2684 | NR>R | no |


## Prior-PR comparators (not expected to hit 0.019)

### CUSTOM_TJ_CORE

repo CUSTOM_TJ_CORE (prior PRs)  
Genes used: 24/26

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.122 | 0.17 | 0.06868 | 0.03434 | NR>R | no |
| fresh | 5 / 6 | -0.122 | 0.167 | 0.1255 | 0.06277 | NR>R | no |

### coreTJ_40

40-gene core TJ (prior PRs; two-sided MW ~0.115)  
Genes used: 40/42

| subset | n R / NR | median R | median NR | MW two-sided | MW one-sided NR>R | direction | hit 0.019 |
|---|---|---|---|---|---|---|---|
| all | 5 / 11 | -0.0762 | 0.157 | 0.1149 | 0.05746 | NR>R | no |
| fresh | 5 / 6 | -0.0762 | 0.196 | 0.1775 | 0.08874 | NR>R | no |


## Honest limits

- n=16. Only a large, well-aligned effect gives p=0.019.
- FFPE is completely nested in NR.
- Cho 2020 has no TJ list. Using the same cohort's DE genes as a "TJ list" would be circular and was not done.
- ESTIMATE here is mean-z of the published 141-gene lists, not the original Affymetrix ssGSEA + purity cosine.
- Multiple tests are reported. A single p≤0.019 among many definitions is not independent confirmation of the slide.

## Rerun

```bash
python3 -m pip install -r scripts/rework/B4_wave2/requirements.txt
python3 scripts/rework/B4_wave2/download.py
python3 scripts/rework/B4_wave2/analyze.py
```

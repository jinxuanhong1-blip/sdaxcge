# B5 wave 2: lung-first then all-open ICI, CLDN4 and TACSTD2

## User-reported number (not produced here)

CLDN4-high ICI **OR = 0.42 [0.18–0.95], k=11**.
Prior open tests: IMvigor210 OR=1.47 p=0.214 (opposite); open lung ICI bulk
all NS after purity residual (PR #114); LUSC-only ICI g=−0.29 p=0.46 (PR #205).
This wave recomputes **lung first**, then all independent open ICI, for
**CLDN4 and TACSTD2**, at median / tertile / continuous logistic, on raw
expression **and** ESTIMATEScore residuals, with DCB and RECIST ORR kept
as **separate** endpoints. It does not pick the cell closest to 0.42.

## Locked methods

- Public processed matrices only. No EGA. Liu/Braun raw are controlled;
  processed Zenodo ICB TSVs (CC-BY-4.0, 10.5281/zenodo.7058399) are used.
- ICB_Jung is the same patients as GSE135222 and is **not** an extra study.
- DCB = PFS ≥ 6 months (180 days). Early-censored before 6 months = NA.
- ORR = RECIST CR/PR vs SD+PD. Author/PredictIO R vs NR is a third label.
- ESTIMATE: Yoshihara 2013 rank-ssGSEA (R estimate v1.0.13 port). Residual
  = OLS of the gene on ESTIMATEScore. Residual skipped if a signature has
  <30 genes overlapping the matrix.
- Histology split only from pathology fields. GSE126044 / GSE135222 /
  GSE166449 have no public histology and are not inferred here.
- OR: Woolf + Haldane–Anscombe 0.5 if any 2×2 cell is 0; Fisher p.
- Continuous: logistic OR per 1 SD. Survival: Cox HR per 1 SD and median.
- Meta: DerSimonian–Laird RE + IVW FE on logOR / logHR.
- GSE207422 is neoadjuvant ICI+chemo (RECIST ORR; MPR is not ORR).
- GSE253564 has no public DCB/ORR in GEO; inventory only.
- Van_Allen / Nathanson are CTLA-4 and are not in the PD-1/PD-L1 metas.
- GSE190266 (France-4) is a reduced panel: CLDN4 present, TACSTD2 absent.
- GSE283829 RECIST labels in GEO are CR/SD/PD (no PR); ORR is CR vs SD+PD.
- GSE190265 histology is squamous vs non-squamous (not LUAD-named).
  Non-squamous is coded nonsquamous and enters the nonLUSC split only.

## Inventory

| Cohort | Cancer | CLDN4 | TACSTD2 | DCB n (R/NR) | ORR n (R/NR) | Author n | Histology | ESTIMATE overlap |
|---|---|---|---|---|---|---|---|---|
| GSE126044 | NSCLC | yes | yes | — | — | 16 (5/11) | NA:16 | stromal=139,immune=141 |
| GSE135222 | NSCLC | yes | yes | 27 (7/20) | — | — | NA:27 | stromal=140,immune=140 |
| GSE166449 | NSCLC | yes | yes | — | — | 22 (7/15) | NA:22 | stromal=140,immune=141 |
| GSE190265 | NSCLC | yes | yes | 43 (14/29) | — | — | NA:18,nonsquamous:14,LUSC:11 | stromal=141,immune=133 |
| GSE190266 | NSCLC | yes | no | 69 (17/52) | — | — | nonsquamous:54,LUSC:14,NA:2 | stromal=93,immune=97 |
| GSE207422 | NSCLC | yes | yes | — | 24 (17/7) | 24 (9/15) | LUSC:12,LUAD:8,NA:4 | stromal=141,immune=141 |
| GSE253564 | NSCLC | yes | yes | — | — | — | LUAD:20,LUSC:10,NA:2 | stromal=138,immune=138 |
| GSE283829 | NSCLC | yes | yes | — | 27 (7/20) | — | LUAD:16,LUSC:9,NA:2 | stromal=141,immune=141 |
| ICB_Jung | NSCLC | yes | yes | 27 (7/20) | 0 (0/0) | 27 (8/19) | NA:27 | stromal=141,immune=141 |
| Braun | ccRCC | yes | yes | 170 (65/105) | 172 (39/133) | 142 (39/103) | NA:181 | stromal=138,immune=140 |
| Gide | melanoma | yes | yes | 41 (22/19) | 41 (19/22) | 38 (19/19) | NA:41 | stromal=141,immune=141 |
| Hugo | melanoma | yes | yes | 0 (0/0) | 27 (14/13) | 27 (14/13) | NA:27 | stromal=141,immune=141 |
| IMvigor210 | urothelial | yes | yes | 0 (0/0) | 298 (68/230) | 235 (68/167) | NA:348 | stromal=138,immune=140 |
| Kim | gastric | yes | yes | 0 (0/0) | 45 (13/32) | 31 (13/18) | NA:45 | stromal=141,immune=141 |
| Liu | melanoma | yes | yes | 121 (57/64) | 121 (47/74) | 112 (52/60) | NA:121 | stromal=139,immune=140 |
| Miao1 | RCC | yes | yes | 14 (12/2) | 33 (8/25) | 28 (10/18) | NA:33 | stromal=140,immune=140 |
| Nathanson | melanoma | yes | yes | 0 (0/0) | 0 (0/0) | 24 (10/14) | NA:24 | stromal=141,immune=141 |
| Puch | melanoma | yes | yes | 0 (0/0) | 55 (14/41) | 49 (14/35) | NA:55 | stromal=138,immune=139 |
| Riaz | melanoma | yes | yes | 0 (0/0) | 44 (9/35) | 30 (9/21) | NA:46 | stromal=141,immune=141 |
| Shiuan | RCC | yes | yes | 0 (0/0) | 13 (6/7) | 13 (6/7) | NA:15 | stromal=137,immune=139 |
| Snyder | urothelial | yes | yes | 25 (9/16) | 21 (7/14) | 18 (8/10) | NA:25 | stromal=140,immune=140 |
| Van_Allen | melanoma | yes | yes | 42 (13/29) | 41 (7/34) | 39 (9/30) | NA:42 | stromal=140,immune=140 |

## Locked pooled results (residual, median)

| Pool | Gene | Endpoint | k | n | RE OR [95% CI] | p | I² | matches 0.42 | 0.42 in CI |
|---|---|---|---:|---:|---|---:|---:|---|---|
| lung_DCB | CLDN4 | DCB | 3 | 139 | 0.94 [0.44–1.99] | 0.873 | 0% | no | no |
| lung_ORR | CLDN4 | ORR | 2 | 51 | 0.64 [0.14–2.96] | 0.571 | 29% | no | yes |
| all_open_DCB | CLDN4 | DCB | 8 | 510 | 0.95 [0.66–1.36] | 0.766 | 0% | no | no |
| all_open_ORR | CLDN4 | ORR | 13 | 921 | 0.85 [0.63–1.14] | 0.274 | 0% | no | no |
| lung_DCB | TACSTD2 | DCB | 2 | 70 | 1.06 [0.38–2.98] | 0.913 | 0% | no | yes |
| lung_ORR | TACSTD2 | ORR | 2 | 51 | 0.64 [0.18–2.21] | 0.479 | 0% | no | yes |
| all_open_DCB | TACSTD2 | DCB | 7 | 441 | 0.92 [0.56–1.53] | 0.753 | 25% | no | no |
| all_open_ORR | TACSTD2 | ORR | 13 | 921 | 0.90 [0.65–1.25] | 0.536 | 8% | no | no |

## All pooled cells (do not cherry-pick)

Every gene × raw/residual × cutoff × endpoint × subset that was estimable.
`matches_user_0.42_at_2dp` is true only when the RE point estimate rounds to 0.42.

| Subset | Gene | Measure | Cutoff | Endpoint | k | RE [95% CI] | p | I² | 0.42@2dp | 0.42 in CI |
|---|---|---|---|---|---:|---|---:|---:|---|---|
| all_open_DCB | CLDN4 | raw | continuous | DCB | 8 | 1.08 [0.89–1.31] | 0.424 | 0% | no | no |
| all_open_DCB | CLDN4 | raw | median | DCB | 8 | 1.17 [0.81–1.70] | 0.407 | 1% | no | no |
| all_open_DCB | CLDN4 | raw | tertile | DCB | 8 | 1.02 [0.65–1.59] | 0.926 | 0% | no | no |
| all_open_DCB | CLDN4 | residual | continuous | DCB | 8 | 1.01 [0.84–1.22] | 0.915 | 0% | no | no |
| all_open_DCB | CLDN4 | residual | median | DCB | 8 | 0.95 [0.66–1.36] | 0.766 | 0% | no | no |
| all_open_DCB | CLDN4 | residual | tertile | DCB | 8 | 1.00 [0.64–1.56] | 0.995 | 0% | no | no |
| all_open_DCB | TACSTD2 | raw | continuous | DCB | 7 | 1.00 [0.78–1.29] | 0.99 | 23% | no | no |
| all_open_DCB | TACSTD2 | raw | median | DCB | 7 | 1.07 [0.63–1.82] | 0.805 | 31% | no | no |
| all_open_DCB | TACSTD2 | raw | tertile | DCB | 7 | 1.03 [0.51–2.10] | 0.935 | 36% | no | no |
| all_open_DCB | TACSTD2 | residual | continuous | DCB | 7 | 0.96 [0.78–1.18] | 0.693 | 5% | no | no |
| all_open_DCB | TACSTD2 | residual | median | DCB | 7 | 0.92 [0.56–1.53] | 0.753 | 25% | no | no |
| all_open_DCB | TACSTD2 | residual | tertile | DCB | 7 | 0.86 [0.44–1.69] | 0.666 | 32% | no | no |
| all_open_ORR | CLDN4 | raw | continuous | ORR | 13 | 0.91 [0.78–1.05] | 0.204 | 0% | no | no |
| all_open_ORR | CLDN4 | raw | median | ORR | 13 | 0.89 [0.66–1.20] | 0.437 | 0% | no | no |
| all_open_ORR | CLDN4 | raw | tertile | ORR | 13 | 0.82 [0.51–1.31] | 0.406 | 27% | no | no |
| all_open_ORR | CLDN4 | residual | continuous | ORR | 13 | 0.88 [0.76–1.03] | 0.108 | 0% | no | no |
| all_open_ORR | CLDN4 | residual | median | ORR | 13 | 0.85 [0.63–1.14] | 0.274 | 0% | no | no |
| all_open_ORR | CLDN4 | residual | tertile | ORR | 13 | 0.83 [0.54–1.29] | 0.415 | 19% | no | no |
| all_open_ORR | TACSTD2 | raw | continuous | ORR | 13 | 0.92 [0.80–1.07] | 0.298 | 0% | no | no |
| all_open_ORR | TACSTD2 | raw | median | ORR | 13 | 0.94 [0.70–1.26] | 0.686 | 0% | no | no |
| all_open_ORR | TACSTD2 | raw | tertile | ORR | 13 | 0.84 [0.58–1.20] | 0.33 | 0% | no | no |
| all_open_ORR | TACSTD2 | residual | continuous | ORR | 13 | 0.86 [0.74–1.00] | 0.0531 | 0% | no | no |
| all_open_ORR | TACSTD2 | residual | median | ORR | 13 | 0.90 [0.65–1.25] | 0.536 | 8% | no | no |
| all_open_ORR | TACSTD2 | residual | tertile | ORR | 13 | 0.58 [0.35–0.96] | 0.0332 | 33% | no | yes |
| all_open_PFS | CLDN4 | raw | continuous | PFS_Cox | 8 | 1.03 [0.93–1.13] | 0.587 | 0% | no | no |
| all_open_PFS | CLDN4 | raw | median | PFS_Cox | 8 | 0.94 [0.76–1.18] | 0.612 | 13% | no | no |
| all_open_PFS | CLDN4 | residual | continuous | PFS_Cox | 8 | 1.06 [0.96–1.17] | 0.228 | 0% | no | no |
| all_open_PFS | CLDN4 | residual | median | PFS_Cox | 8 | 1.07 [0.88–1.31] | 0.472 | 0% | no | no |
| all_open_PFS | TACSTD2 | raw | continuous | PFS_Cox | 7 | 1.07 [0.97–1.19] | 0.175 | 0% | no | no |
| all_open_PFS | TACSTD2 | raw | median | PFS_Cox | 7 | 1.07 [0.86–1.34] | 0.551 | 5% | no | no |
| all_open_PFS | TACSTD2 | residual | continuous | PFS_Cox | 7 | 1.10 [0.99–1.22] | 0.0656 | 0% | no | no |
| all_open_PFS | TACSTD2 | residual | median | PFS_Cox | 7 | 1.11 [0.90–1.37] | 0.329 | 0% | no | no |
| lung_DCB | CLDN4 | raw | continuous | DCB | 3 | 1.18 [0.76–1.83] | 0.452 | 9% | no | no |
| lung_DCB | CLDN4 | raw | median | DCB | 3 | 1.69 [0.79–3.63] | 0.177 | 0% | no | no |
| lung_DCB | CLDN4 | raw | tertile | DCB | 3 | 1.65 [0.67–4.02] | 0.274 | 0% | no | no |
| lung_DCB | CLDN4 | residual | continuous | DCB | 3 | 0.98 [0.67–1.44] | 0.921 | 0% | no | no |
| lung_DCB | CLDN4 | residual | median | DCB | 3 | 0.94 [0.44–1.99] | 0.873 | 0% | no | no |
| lung_DCB | CLDN4 | residual | tertile | DCB | 3 | 1.00 [0.41–2.43] | 0.997 | 0% | no | yes |
| lung_DCB | TACSTD2 | raw | continuous | DCB | 2 | 0.92 [0.55–1.54] | 0.75 | 0% | no | no |
| lung_DCB | TACSTD2 | raw | median | DCB | 2 | 1.06 [0.38–2.98] | 0.913 | 0% | no | yes |
| lung_DCB | TACSTD2 | raw | tertile | DCB | 2 | 1.00 [0.26–3.85] | 1 | 0% | no | yes |
| lung_DCB | TACSTD2 | residual | continuous | DCB | 2 | 0.97 [0.57–1.62] | 0.894 | 0% | no | no |
| lung_DCB | TACSTD2 | residual | median | DCB | 2 | 1.06 [0.38–2.98] | 0.913 | 0% | no | yes |
| lung_DCB | TACSTD2 | residual | tertile | DCB | 2 | 0.81 [0.22–3.01] | 0.754 | 0% | no | yes |
| lung_LUAD_ORR | CLDN4 | raw | continuous | ORR | 2 | 1.60 [0.57–4.45] | 0.37 | 0% | no | no |
| lung_LUAD_ORR | CLDN4 | raw | median | ORR | 2 | 0.59 [0.03–11.43] | 0.725 | 52% | no | yes |
| lung_LUAD_ORR | CLDN4 | raw | tertile | ORR | 2 | 0.99 [0.11–9.07] | 0.995 | 7% | no | yes |
| lung_LUAD_ORR | CLDN4 | residual | continuous | ORR | 2 | 1.71 [0.65–4.51] | 0.278 | 0% | no | no |
| lung_LUAD_ORR | CLDN4 | residual | median | ORR | 2 | 0.59 [0.03–11.43] | 0.725 | 52% | no | yes |
| lung_LUAD_ORR | CLDN4 | residual | tertile | ORR | 2 | 0.94 [0.07–12.66] | 0.963 | 15% | no | yes |
| lung_LUAD_ORR | TACSTD2 | raw | continuous | ORR | 2 | 1.57 [0.61–4.03] | 0.351 | 0% | no | no |
| lung_LUAD_ORR | TACSTD2 | raw | median | ORR | 2 | 0.25 [0.03–1.90] | 0.18 | 0% | no | yes |
| lung_LUAD_ORR | TACSTD2 | raw | tertile | ORR | 2 | 0.54 [0.06–5.19] | 0.594 | 0% | no | yes |
| lung_LUAD_ORR | TACSTD2 | residual | continuous | ORR | 2 | 1.70 [0.64–4.52] | 0.288 | 0% | no | no |
| lung_LUAD_ORR | TACSTD2 | residual | median | ORR | 2 | 4.03 [0.53–30.82] | 0.18 | 0% | no | no |
| lung_LUAD_ORR | TACSTD2 | residual | tertile | ORR | 2 | 1.85 [0.19–17.80] | 0.594 | 0% | no | yes |
| lung_LUSC_DCB | CLDN4 | raw | continuous | DCB | 2 | 0.57 [0.16–2.07] | 0.395 | 0% | no | yes |
| lung_LUSC_DCB | CLDN4 | raw | median | DCB | 2 | 1.89 [0.27–13.40] | 0.523 | 0% | no | yes |
| lung_LUSC_DCB | CLDN4 | raw | tertile | DCB | 2 | 0.49 [0.04–5.48] | 0.562 | 0% | no | yes |
| lung_LUSC_DCB | CLDN4 | residual | continuous | DCB | 2 | 0.33 [0.07–1.65] | 0.178 | 0% | no | yes |
| lung_LUSC_DCB | CLDN4 | residual | median | DCB | 2 | 0.82 [0.12–5.79] | 0.84 | 0% | no | yes |
| lung_LUSC_DCB | CLDN4 | residual | tertile | DCB | 2 | 0.31 [0.03–2.92] | 0.303 | 0% | no | yes |
| lung_LUSC_DCB | TACSTD2 | raw | continuous | DCB | 1 | 0.22 [0.03–1.79] | 0.156 | 0% | no | yes |
| lung_LUSC_DCB | TACSTD2 | raw | median | DCB | 1 | 0.33 [0.03–3.93] | 0.383 | 0% | no | yes |
| lung_LUSC_DCB | TACSTD2 | raw | tertile | DCB | 1 | 0.33 [0.02–6.65] | 0.472 | 0% | no | yes |
| lung_LUSC_DCB | TACSTD2 | residual | continuous | DCB | 1 | 0.26 [0.04–1.74] | 0.165 | 0% | no | yes |
| lung_LUSC_DCB | TACSTD2 | residual | median | DCB | 1 | 0.33 [0.03–3.93] | 0.383 | 0% | no | yes |
| lung_LUSC_DCB | TACSTD2 | residual | tertile | DCB | 1 | 0.33 [0.02–6.65] | 0.472 | 0% | no | yes |
| lung_LUSC_ORR | CLDN4 | raw | continuous | ORR | 2 | 2.24 [0.71–7.06] | 0.167 | 0% | no | no |
| lung_LUSC_ORR | CLDN4 | raw | median | ORR | 2 | 1.35 [0.09–20.09] | 0.828 | 37% | no | yes |
| lung_LUSC_ORR | CLDN4 | raw | tertile | ORR | 2 | 0.97 [0.08–11.23] | 0.979 | 11% | no | yes |
| lung_LUSC_ORR | CLDN4 | residual | continuous | ORR | 2 | 1.50 [0.52–4.36] | 0.456 | 0% | no | no |
| lung_LUSC_ORR | CLDN4 | residual | median | ORR | 2 | 1.35 [0.09–20.09] | 0.828 | 37% | no | yes |
| lung_LUSC_ORR | CLDN4 | residual | tertile | ORR | 2 | 0.97 [0.08–11.23] | 0.979 | 11% | no | yes |
| lung_LUSC_ORR | TACSTD2 | raw | continuous | ORR | 2 | 2.00 [0.66–6.01] | 0.219 | 0% | no | no |
| lung_LUSC_ORR | TACSTD2 | raw | median | ORR | 2 | 1.49 [0.19–11.76] | 0.705 | 0% | no | yes |
| lung_LUSC_ORR | TACSTD2 | raw | tertile | ORR | 2 | 2.22 [0.15–32.44] | 0.561 | 0% | no | yes |
| lung_LUSC_ORR | TACSTD2 | residual | continuous | ORR | 2 | 1.21 [0.43–3.36] | 0.719 | 0% | no | no |
| lung_LUSC_ORR | TACSTD2 | residual | median | ORR | 2 | 0.52 [0.07–4.14] | 0.54 | 0% | no | yes |
| lung_LUSC_ORR | TACSTD2 | residual | tertile | ORR | 2 | 0.29 [0.03–2.89] | 0.292 | 0% | no | yes |
| lung_ORR | CLDN4 | raw | continuous | ORR | 2 | 0.80 [0.43–1.47] | 0.467 | 0% | no | no |
| lung_ORR | CLDN4 | raw | median | ORR | 2 | 0.43 [0.12–1.55] | 0.197 | 0% | no | yes |
| lung_ORR | CLDN4 | raw | tertile | ORR | 2 | 0.77 [0.15–3.91] | 0.752 | 16% | no | yes |
| lung_ORR | CLDN4 | residual | continuous | ORR | 2 | 0.68 [0.36–1.29] | 0.242 | 0% | no | yes |
| lung_ORR | CLDN4 | residual | median | ORR | 2 | 0.64 [0.14–2.96] | 0.571 | 29% | no | yes |
| lung_ORR | CLDN4 | residual | tertile | ORR | 2 | 0.34 [0.08–1.53] | 0.161 | 0% | no | yes |
| lung_ORR | TACSTD2 | raw | continuous | ORR | 2 | 0.84 [0.36–1.92] | 0.672 | 43% | no | yes |
| lung_ORR | TACSTD2 | raw | median | ORR | 2 | 0.64 [0.18–2.21] | 0.479 | 0% | no | yes |
| lung_ORR | TACSTD2 | raw | tertile | ORR | 2 | 0.74 [0.16–3.44] | 0.701 | 0% | no | yes |
| lung_ORR | TACSTD2 | residual | continuous | ORR | 2 | 0.70 [0.34–1.46] | 0.339 | 22% | no | yes |
| lung_ORR | TACSTD2 | residual | median | ORR | 2 | 0.64 [0.18–2.21] | 0.479 | 0% | no | yes |
| lung_ORR | TACSTD2 | residual | tertile | ORR | 2 | 0.40 [0.08–2.02] | 0.265 | 0% | no | yes |
| lung_PFS | CLDN4 | raw | continuous | PFS_Cox | 3 | 0.96 [0.78–1.19] | 0.71 | 23% | no | no |
| lung_PFS | CLDN4 | raw | median | PFS_Cox | 3 | 0.83 [0.57–1.22] | 0.35 | 0% | no | no |
| lung_PFS | CLDN4 | residual | continuous | PFS_Cox | 3 | 1.05 [0.88–1.25] | 0.599 | 0% | no | no |
| lung_PFS | CLDN4 | residual | median | PFS_Cox | 3 | 1.10 [0.75–1.61] | 0.627 | 0% | no | no |
| lung_PFS | TACSTD2 | raw | continuous | PFS_Cox | 2 | 1.13 [0.85–1.50] | 0.398 | 0% | no | no |
| lung_PFS | TACSTD2 | raw | median | PFS_Cox | 2 | 1.10 [0.65–1.87] | 0.723 | 0% | no | no |
| lung_PFS | TACSTD2 | residual | continuous | PFS_Cox | 2 | 1.10 [0.83–1.44] | 0.507 | 0% | no | no |
| lung_PFS | TACSTD2 | residual | median | PFS_Cox | 2 | 1.10 [0.65–1.87] | 0.723 | 0% | no | no |
| lung_author | CLDN4 | raw | continuous | author_R_vs_NR | 2 | 0.72 [0.24–2.10] | 0.543 | 51% | no | yes |
| lung_author | CLDN4 | raw | median | author_R_vs_NR | 2 | 0.39 [0.09–1.69] | 0.209 | 0% | no | yes |
| lung_author | CLDN4 | raw | tertile | author_R_vs_NR | 2 | 0.51 [0.09–2.88] | 0.447 | 0% | no | yes |
| lung_author | CLDN4 | residual | continuous | author_R_vs_NR | 2 | 1.07 [0.53–2.18] | 0.852 | 0% | no | no |
| lung_author | CLDN4 | residual | median | author_R_vs_NR | 2 | 0.61 [0.15–2.44] | 0.487 | 0% | no | yes |
| lung_author | CLDN4 | residual | tertile | author_R_vs_NR | 2 | 1.00 [0.21–4.86] | 1 | 0% | no | yes |
| lung_author | TACSTD2 | raw | continuous | author_R_vs_NR | 2 | 0.81 [0.29–2.26] | 0.683 | 45% | no | yes |
| lung_author | TACSTD2 | raw | median | author_R_vs_NR | 2 | 1.53 [0.24–9.91] | 0.655 | 40% | no | yes |
| lung_author | TACSTD2 | raw | tertile | author_R_vs_NR | 2 | 1.47 [0.25–8.66] | 0.672 | 0% | no | yes |
| lung_author | TACSTD2 | residual | continuous | author_R_vs_NR | 2 | 0.94 [0.46–1.92] | 0.874 | 0% | no | no |
| lung_author | TACSTD2 | residual | median | author_R_vs_NR | 2 | 1.00 [0.25–3.98] | 0.997 | 0% | no | yes |
| lung_author | TACSTD2 | residual | tertile | author_R_vs_NR | 2 | 1.89 [0.39–9.17] | 0.43 | 0% | no | yes |
| lung_nonLUSC_DCB | CLDN4 | raw | continuous | DCB | 2 | 2.35 [0.69–8.00] | 0.172 | 37% | no | no |
| lung_nonLUSC_DCB | CLDN4 | raw | median | DCB | 2 | 4.05 [0.79–20.82] | 0.0943 | 24% | no | no |
| lung_nonLUSC_DCB | CLDN4 | raw | tertile | DCB | 2 | 4.43 [1.18–16.67] | 0.0276 | 0% | no | no |
| lung_nonLUSC_DCB | CLDN4 | residual | continuous | DCB | 2 | 2.05 [0.44–9.64] | 0.362 | 59% | no | no |
| lung_nonLUSC_DCB | CLDN4 | residual | median | DCB | 2 | 3.20 [0.24–42.99] | 0.38 | 61% | no | yes |
| lung_nonLUSC_DCB | CLDN4 | residual | tertile | DCB | 2 | 2.57 [0.20–32.86] | 0.468 | 54% | no | yes |
| lung_nonLUSC_DCB | TACSTD2 | raw | continuous | DCB | 1 | 1.40 [0.35–5.68] | 0.635 | 0% | no | yes |
| lung_nonLUSC_DCB | TACSTD2 | raw | median | DCB | 1 | 1.00 [0.10–10.17] | 1 | 0% | no | yes |
| lung_nonLUSC_DCB | TACSTD2 | raw | tertile | DCB | 1 | 2.67 [0.16–45.14] | 0.497 | 0% | no | yes |
| lung_nonLUSC_DCB | TACSTD2 | residual | continuous | DCB | 1 | 1.64 [0.39–6.85] | 0.5 | 0% | no | yes |
| lung_nonLUSC_DCB | TACSTD2 | residual | median | DCB | 1 | 1.00 [0.10–10.17] | 1 | 0% | no | yes |
| lung_nonLUSC_DCB | TACSTD2 | residual | tertile | DCB | 1 | 2.67 [0.16–45.14] | 0.497 | 0% | no | yes |
| lung_nonLUSC_ORR | CLDN4 | raw | continuous | ORR | 2 | 1.60 [0.57–4.45] | 0.37 | 0% | no | no |
| lung_nonLUSC_ORR | CLDN4 | raw | median | ORR | 2 | 0.59 [0.03–11.43] | 0.725 | 52% | no | yes |
| lung_nonLUSC_ORR | CLDN4 | raw | tertile | ORR | 2 | 0.99 [0.11–9.07] | 0.995 | 7% | no | yes |
| lung_nonLUSC_ORR | CLDN4 | residual | continuous | ORR | 2 | 1.71 [0.65–4.51] | 0.278 | 0% | no | no |
| lung_nonLUSC_ORR | CLDN4 | residual | median | ORR | 2 | 0.59 [0.03–11.43] | 0.725 | 52% | no | yes |
| lung_nonLUSC_ORR | CLDN4 | residual | tertile | ORR | 2 | 0.94 [0.07–12.66] | 0.963 | 15% | no | yes |
| lung_nonLUSC_ORR | TACSTD2 | raw | continuous | ORR | 2 | 1.57 [0.61–4.03] | 0.351 | 0% | no | no |
| lung_nonLUSC_ORR | TACSTD2 | raw | median | ORR | 2 | 0.25 [0.03–1.90] | 0.18 | 0% | no | yes |
| lung_nonLUSC_ORR | TACSTD2 | raw | tertile | ORR | 2 | 0.54 [0.06–5.19] | 0.594 | 0% | no | yes |
| lung_nonLUSC_ORR | TACSTD2 | residual | continuous | ORR | 2 | 1.70 [0.64–4.52] | 0.288 | 0% | no | no |
| lung_nonLUSC_ORR | TACSTD2 | residual | median | ORR | 2 | 4.03 [0.53–30.82] | 0.18 | 0% | no | no |
| lung_nonLUSC_ORR | TACSTD2 | residual | tertile | ORR | 2 | 1.85 [0.19–17.80] | 0.594 | 0% | no | yes |
| nonlung_DCB | CLDN4 | raw | continuous | DCB | 5 | 1.05 [0.85–1.31] | 0.627 | 0% | no | no |
| nonlung_DCB | CLDN4 | raw | median | DCB | 5 | 1.00 [0.57–1.77] | 0.987 | 29% | no | no |
| nonlung_DCB | CLDN4 | raw | tertile | DCB | 5 | 0.87 [0.52–1.46] | 0.6 | 0% | no | no |
| nonlung_DCB | CLDN4 | residual | continuous | DCB | 5 | 1.02 [0.82–1.26] | 0.859 | 0% | no | no |
| nonlung_DCB | CLDN4 | residual | median | DCB | 5 | 0.95 [0.62–1.44] | 0.802 | 0% | no | no |
| nonlung_DCB | CLDN4 | residual | tertile | DCB | 5 | 1.00 [0.60–1.67] | 0.992 | 0% | no | no |
| nonlung_DCB | TACSTD2 | raw | continuous | DCB | 5 | 1.00 [0.70–1.43] | 0.998 | 46% | no | no |
| nonlung_DCB | TACSTD2 | raw | median | DCB | 5 | 1.09 [0.53–2.23] | 0.824 | 50% | no | no |
| nonlung_DCB | TACSTD2 | raw | tertile | DCB | 5 | 1.00 [0.38–2.64] | 0.997 | 57% | no | yes |
| nonlung_DCB | TACSTD2 | residual | continuous | DCB | 5 | 0.93 [0.68–1.28] | 0.674 | 36% | no | no |
| nonlung_DCB | TACSTD2 | residual | median | DCB | 5 | 0.91 [0.46–1.78] | 0.773 | 45% | no | no |
| nonlung_DCB | TACSTD2 | residual | tertile | DCB | 5 | 0.86 [0.34–2.15] | 0.742 | 53% | no | yes |
| nonlung_ORR | CLDN4 | raw | continuous | ORR | 11 | 0.92 [0.78–1.07] | 0.26 | 0% | no | no |
| nonlung_ORR | CLDN4 | raw | median | ORR | 11 | 0.92 [0.68–1.26] | 0.607 | 2% | no | no |
| nonlung_ORR | CLDN4 | raw | tertile | ORR | 11 | 0.82 [0.49–1.38] | 0.463 | 34% | no | no |
| nonlung_ORR | CLDN4 | residual | continuous | ORR | 11 | 0.90 [0.77–1.05] | 0.171 | 0% | no | no |
| nonlung_ORR | CLDN4 | residual | median | ORR | 11 | 0.86 [0.63–1.17] | 0.331 | 0% | no | no |
| nonlung_ORR | CLDN4 | residual | tertile | ORR | 11 | 0.90 [0.57–1.42] | 0.644 | 21% | no | no |
| nonlung_ORR | TACSTD2 | raw | continuous | ORR | 11 | 0.93 [0.80–1.08] | 0.351 | 0% | no | no |
| nonlung_ORR | TACSTD2 | raw | median | ORR | 11 | 0.96 [0.71–1.30] | 0.808 | 0% | no | no |
| nonlung_ORR | TACSTD2 | raw | tertile | ORR | 11 | 0.79 [0.52–1.21] | 0.278 | 15% | no | no |
| nonlung_ORR | TACSTD2 | residual | continuous | ORR | 11 | 0.87 [0.75–1.02] | 0.083 | 0% | no | no |
| nonlung_ORR | TACSTD2 | residual | median | ORR | 11 | 0.88 [0.61–1.27] | 0.49 | 21% | no | no |
| nonlung_ORR | TACSTD2 | residual | tertile | ORR | 11 | 0.58 [0.33–1.02] | 0.0573 | 42% | no | yes |

## Per-cohort residual median (CLDN4 and TACSTD2)

| Cohort | Gene | Endpoint | n (R/NR) | High R/NR | Low R/NR | OR [95% CI] | p |
|---|---|---|---|---|---|---|---:|
| Braun | CLDN4 | DCB | 170 (65/105) | 32/53 | 33/52 | 0.95 [0.51–1.77] | 1 |
| GSE135222 | CLDN4 | DCB | 27 (7/20) | 4/10 | 3/10 | 1.33 [0.24–7.56] | 1 |
| GSE190265 | CLDN4 | DCB | 43 (14/29) | 7/15 | 7/14 | 0.93 [0.26–3.34] | 1 |
| GSE190266 | CLDN4 | DCB | 69 (17/52) | 8/27 | 9/25 | 0.82 [0.27–2.46] | 0.785 |
| Gide | CLDN4 | DCB | 41 (22/19) | 9/12 | 13/7 | 0.40 [0.11–1.43] | 0.215 |
| Liu | CLDN4 | DCB | 121 (57/64) | 29/32 | 28/32 | 1.04 [0.51–2.12] | 1 |
| Miao1 | CLDN4 | DCB | 14 (12/2) | 6/1 | 6/1 | 1.00 [0.05–19.96] | 1 |
| Snyder | CLDN4 | DCB | 25 (9/16) | 6/7 | 3/9 | 2.57 [0.47–14.10] | 0.411 |
| Van_Allen | CLDN4 | DCB | 42 (13/29) | 7/14 | 6/15 | 1.25 [0.34–4.64] | 1 |
| Braun | TACSTD2 | DCB | 170 (65/105) | 32/53 | 33/52 | 0.95 [0.51–1.77] | 1 |
| GSE135222 | TACSTD2 | DCB | 27 (7/20) | 3/11 | 4/9 | 0.61 [0.11–3.49] | 0.678 |
| GSE190265 | TACSTD2 | DCB | 43 (14/29) | 8/14 | 6/15 | 1.43 [0.40–5.16] | 0.747 |
| Gide | TACSTD2 | DCB | 41 (22/19) | 8/13 | 14/6 | 0.26 [0.07–0.97] | 0.0616 |
| Liu | TACSTD2 | DCB | 121 (57/64) | 28/33 | 29/31 | 0.91 [0.44–1.85] | 0.856 |
| Miao1 | TACSTD2 | DCB | 14 (12/2) | 6/1 | 6/1 | 1.00 [0.05–19.96] | 1 |
| Snyder | TACSTD2 | DCB | 25 (9/16) | 7/6 | 2/10 | 5.83 [0.90–37.82] | 0.0968 |
| Van_Allen | TACSTD2 | DCB | 42 (13/29) | 8/13 | 5/16 | 1.97 [0.52–7.49] | 0.505 |
| Braun | CLDN4 | ORR | 172 (39/133) | 16/70 | 23/63 | 0.63 [0.30–1.29] | 0.274 |
| GSE207422 | CLDN4 | ORR | 24 (17/7) | 7/5 | 10/2 | 0.28 [0.04–1.88] | 0.371 |
| GSE283829 | CLDN4 | ORR | 27 (7/20) | 4/10 | 3/10 | 1.33 [0.24–7.56] | 1 |
| Gide | CLDN4 | ORR | 41 (19/22) | 7/14 | 12/8 | 0.33 [0.09–1.19] | 0.121 |
| Hugo | CLDN4 | ORR | 27 (14/13) | 8/6 | 6/7 | 1.56 [0.34–7.11] | 0.706 |
| IMvigor210 | CLDN4 | ORR | 298 (68/230) | 38/111 | 30/119 | 1.36 [0.79–2.34] | 0.334 |
| Kim | CLDN4 | ORR | 45 (13/32) | 7/16 | 6/16 | 1.17 [0.32–4.25] | 1 |
| Liu | CLDN4 | ORR | 121 (47/74) | 20/41 | 27/33 | 0.60 [0.29–1.25] | 0.194 |
| Miao1 | CLDN4 | ORR | 33 (8/25) | 4/13 | 4/12 | 0.92 [0.19–4.54] | 1 |
| Puch | CLDN4 | ORR | 55 (14/41) | 7/21 | 7/20 | 0.95 [0.28–3.21] | 1 |
| Riaz | CLDN4 | ORR | 44 (9/35) | 4/18 | 5/17 | 0.76 [0.17–3.29] | 1 |
| Shiuan | CLDN4 | ORR | 13 (6/7) | 3/4 | 3/3 | 0.75 [0.08–6.71] | 1 |
| Snyder | CLDN4 | ORR | 21 (7/14) | 3/8 | 4/6 | 0.56 [0.09–3.52] | 0.659 |
| Van_Allen | CLDN4 | ORR | 41 (7/34) | 2/19 | 5/15 | 0.32 [0.05–1.86] | 0.238 |
| Braun | TACSTD2 | ORR | 172 (39/133) | 19/67 | 20/66 | 0.94 [0.46–1.91] | 1 |
| GSE207422 | TACSTD2 | ORR | 24 (17/7) | 8/4 | 9/3 | 0.67 [0.11–3.93] | 1 |
| GSE283829 | TACSTD2 | ORR | 27 (7/20) | 3/11 | 4/9 | 0.61 [0.11–3.49] | 0.678 |
| Gide | TACSTD2 | ORR | 41 (19/22) | 7/14 | 12/8 | 0.33 [0.09–1.19] | 0.121 |
| Hugo | TACSTD2 | ORR | 27 (14/13) | 6/8 | 8/5 | 0.47 [0.10–2.18] | 0.449 |
| IMvigor210 | TACSTD2 | ORR | 298 (68/230) | 42/107 | 26/123 | 1.86 [1.07–3.23] | 0.0379 |
| Kim | TACSTD2 | ORR | 45 (13/32) | 7/16 | 6/16 | 1.17 [0.32–4.25] | 1 |
| Liu | TACSTD2 | ORR | 121 (47/74) | 23/38 | 24/36 | 0.91 [0.44–1.89] | 0.853 |
| Miao1 | TACSTD2 | ORR | 33 (8/25) | 3/14 | 5/11 | 0.47 [0.09–2.42] | 0.438 |
| Puch | TACSTD2 | ORR | 55 (14/41) | 5/23 | 9/18 | 0.43 [0.12–1.53] | 0.227 |
| Riaz | TACSTD2 | ORR | 44 (9/35) | 3/19 | 6/16 | 0.42 [0.09–1.96] | 0.457 |
| Shiuan | TACSTD2 | ORR | 13 (6/7) | 3/4 | 3/3 | 0.75 [0.08–6.71] | 1 |
| Snyder | TACSTD2 | ORR | 21 (7/14) | 4/7 | 3/7 | 1.33 [0.21–8.29] | 1 |
| Van_Allen | TACSTD2 | ORR | 41 (7/34) | 5/16 | 2/18 | 2.81 [0.48–16.56] | 0.41 |
| Braun | CLDN4 | author_R_vs_NR | 142 (39/103) | 18/53 | 21/50 | 0.81 [0.39–1.69] | 0.707 |
| GSE126044 | CLDN4 | author_R_vs_NR | 16 (5/11) | 2/6 | 3/5 | 0.56 [0.06–4.76] | 1 |
| GSE166449 | CLDN4 | author_R_vs_NR | 22 (7/15) | 3/8 | 4/7 | 0.66 [0.11–4.00] | 1 |
| GSE207422 | CLDN4 | author_R_vs_NR | 24 (9/15) | 4/8 | 5/7 | 0.70 [0.13–3.68] | 1 |
| Gide | CLDN4 | author_R_vs_NR | 38 (19/19) | 7/12 | 12/7 | 0.34 [0.09–1.27] | 0.194 |
| Hugo | CLDN4 | author_R_vs_NR | 27 (14/13) | 8/6 | 6/7 | 1.56 [0.34–7.11] | 0.706 |
| IMvigor210 | CLDN4 | author_R_vs_NR | 235 (68/167) | 41/77 | 27/90 | 1.77 [1.00–3.15] | 0.0613 |
| Kim | CLDN4 | author_R_vs_NR | 31 (13/18) | 7/9 | 6/9 | 1.17 [0.28–4.87] | 1 |
| Liu | CLDN4 | author_R_vs_NR | 112 (52/60) | 22/34 | 30/26 | 0.56 [0.26–1.19] | 0.185 |
| Miao1 | CLDN4 | author_R_vs_NR | 28 (10/18) | 5/9 | 5/9 | 1.00 [0.21–4.69] | 1 |
| Nathanson | CLDN4 | author_R_vs_NR | 24 (10/14) | — | — | too_few | — |
| Puch | CLDN4 | author_R_vs_NR | 49 (14/35) | 7/18 | 7/17 | 0.94 [0.27–3.26] | 1 |
| Riaz | CLDN4 | author_R_vs_NR | 30 (9/21) | 4/11 | 5/10 | 0.73 [0.15–3.49] | 1 |
| Shiuan | CLDN4 | author_R_vs_NR | 13 (6/7) | 3/4 | 3/3 | 0.75 [0.08–6.71] | 1 |
| Snyder | CLDN4 | author_R_vs_NR | 18 (8/10) | 5/4 | 3/6 | 2.50 [0.37–16.89] | 0.637 |
| Van_Allen | CLDN4 | author_R_vs_NR | 39 (9/30) | 4/16 | 5/14 | 0.70 [0.16–3.13] | 0.716 |
| Braun | TACSTD2 | author_R_vs_NR | 142 (39/103) | 19/52 | 20/51 | 0.93 [0.45–1.95] | 1 |
| GSE126044 | TACSTD2 | author_R_vs_NR | 16 (5/11) | 3/5 | 2/6 | 1.80 [0.21–15.41] | 1 |
| GSE166449 | TACSTD2 | author_R_vs_NR | 22 (7/15) | 3/8 | 4/7 | 0.66 [0.11–4.00] | 1 |
| GSE207422 | TACSTD2 | author_R_vs_NR | 24 (9/15) | 5/7 | 4/8 | 1.43 [0.27–7.52] | 1 |
| Gide | TACSTD2 | author_R_vs_NR | 38 (19/19) | 7/12 | 12/7 | 0.34 [0.09–1.27] | 0.194 |
| Hugo | TACSTD2 | author_R_vs_NR | 27 (14/13) | 6/8 | 8/5 | 0.47 [0.10–2.18] | 0.449 |
| IMvigor210 | TACSTD2 | author_R_vs_NR | 235 (68/167) | 41/77 | 27/90 | 1.77 [1.00–3.15] | 0.0613 |
| Kim | TACSTD2 | author_R_vs_NR | 31 (13/18) | 7/9 | 6/9 | 1.17 [0.28–4.87] | 1 |
| Liu | TACSTD2 | author_R_vs_NR | 112 (52/60) | 25/31 | 27/29 | 0.87 [0.41–1.82] | 0.85 |
| Miao1 | TACSTD2 | author_R_vs_NR | 28 (10/18) | 4/10 | 6/8 | 0.53 [0.11–2.56] | 0.695 |
| Nathanson | TACSTD2 | author_R_vs_NR | 24 (10/14) | — | — | too_few | — |
| Puch | TACSTD2 | author_R_vs_NR | 49 (14/35) | 5/20 | 9/15 | 0.42 [0.12–1.50] | 0.217 |
| Riaz | TACSTD2 | author_R_vs_NR | 30 (9/21) | 3/12 | 6/9 | 0.38 [0.07–1.92] | 0.427 |
| Shiuan | TACSTD2 | author_R_vs_NR | 13 (6/7) | 3/4 | 3/3 | 0.75 [0.08–6.71] | 1 |
| Snyder | TACSTD2 | author_R_vs_NR | 18 (8/10) | 4/5 | 4/5 | 1.00 [0.16–6.42] | 1 |
| Van_Allen | TACSTD2 | author_R_vs_NR | 39 (9/30) | 5/15 | 4/15 | 1.25 [0.28–5.59] | 1 |

## Cox PFS (residual, continuous per 1 SD)

| Cohort | Gene | n | events | HR [95% CI] | p |
|---|---|---:|---:|---|---:|
| Braun | CLDN4 | 181 | 159 | 1.03 [0.88–1.20] | 0.725 |
| GSE135222 | CLDN4 | 27 | 21 | 1.07 [0.71–1.64] | 0.738 |
| GSE190265 | CLDN4 | 43 | 35 | 1.11 [0.81–1.53] | 0.523 |
| GSE190266 | CLDN4 | 70 | 52 | 1.01 [0.78–1.29] | 0.968 |
| Gide | CLDN4 | 41 | 29 | 1.28 [0.86–1.91] | 0.225 |
| Hugo | CLDN4 | 0 | 0 | not estimable | — |
| IMvigor210 | CLDN4 | 0 | 0 | not estimable | — |
| Kim | CLDN4 | 0 | 0 | not estimable | — |
| Liu | CLDN4 | 121 | 86 | 1.15 [0.93–1.42] | 0.19 |
| Miao1 | CLDN4 | 33 | 8 | 0.62 [0.16–2.41] | 0.494 |
| Nathanson | CLDN4 | 0 | 0 | not estimable | — |
| Puch | CLDN4 | 0 | 0 | not estimable | — |
| Riaz | CLDN4 | 0 | 0 | not estimable | — |
| Shiuan | CLDN4 | 0 | 0 | not estimable | — |
| Snyder | CLDN4 | 25 | 19 | 0.82 [0.49–1.36] | 0.442 |
| Van_Allen | CLDN4 | 42 | 36 | 0.99 [0.70–1.40] | 0.96 |
| Braun | TACSTD2 | 181 | 159 | 1.08 [0.93–1.24] | 0.332 |
| GSE135222 | TACSTD2 | 27 | 21 | 1.05 [0.67–1.64] | 0.833 |
| GSE190265 | TACSTD2 | 43 | 35 | 1.13 [0.80–1.60] | 0.499 |
| Gide | TACSTD2 | 41 | 29 | 1.40 [0.97–2.03] | 0.0712 |
| Hugo | TACSTD2 | 0 | 0 | not estimable | — |
| IMvigor210 | TACSTD2 | 0 | 0 | not estimable | — |
| Kim | TACSTD2 | 0 | 0 | not estimable | — |
| Liu | TACSTD2 | 121 | 86 | 1.17 [0.95–1.45] | 0.142 |
| Miao1 | TACSTD2 | 33 | 8 | 0.54 [0.21–1.38] | 0.195 |
| Nathanson | TACSTD2 | 0 | 0 | not estimable | — |
| Puch | TACSTD2 | 0 | 0 | not estimable | — |
| Riaz | TACSTD2 | 0 | 0 | not estimable | — |
| Shiuan | TACSTD2 | 0 | 0 | not estimable | — |
| Snyder | TACSTD2 | 25 | 19 | 0.88 [0.59–1.34] | 0.561 |
| Van_Allen | TACSTD2 | 42 | 36 | 0.89 [0.63–1.26] | 0.496 |

## Honest verdict

- All-open ORR, CLDN4 residual median: RE OR = **0.85 [0.63–1.14]**, k=13, p=0.274. Does **not** match user 0.42 (observed 0.85).
- Public independent ORR k=13, not the claimed k=11.
- Lung-first DCB, CLDN4 residual median: RE OR = **0.94 [0.44–1.99]**, k=3, p=0.873.
- We did not drop a cohort because it was null or opposite.
- We did not add extra series after seeing results in order to reach k=11.
- DCB and ORR were never pooled as if they were the same endpoint.
- Cross-cancer all-open pooling is not a lung-specific test.
- Small lung n makes tertile splits unstable; empty arms are not recoded.
- IMvigor210 CLDN4 raw median ORR is 1.47 p=0.214 (opposite of 0.42), matching
  the previously reported open urothelial test.
- One non-locked cell is nominally significant (all-open TACSTD2 residual
  tertile ORR OR=0.58 p=0.033). It is **not** adopted: tertile was not the
  locked cutoff, TACSTD2 is not the claimed gene, and median/continuous
  for the same gene/endpoint are NS.

## Reproduction

```bash
python3 -m pip install -r results/rework/B5_wave2/requirements.txt
python3 results/rework/B5_wave2/download.py
python3 results/rework/B5_wave2/analyze.py
```


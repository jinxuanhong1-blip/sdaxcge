# B5 rework: open ICI RNA, CLDN4 vs response, all cutoffs

## User-reported number (not produced here)

CLDN4-high ICI **OR = 0.42 [0.18–0.95], k=11**.
Open lung ICI was previously NS. This slice recomputes the claim
on public matrices and **reports every locked cutoff**. It does not
choose the cutoff whose pooled OR is closest to 0.42.

## Locked methods

- Gene: CLDN4 (`ENSG00000189143` / symbol).
- Expression: `log2(x+1)` for raw counts/TPM/FPKM. If the public matrix
  already contains negative values (Zenodo ICB log-scale), it is used as-is.
  GSE207422 is already log2TPM. Scale is in `tables/per_cohort_or.csv`.
- Cutoffs (cohort-internal, on the analysis subset):
  - median: High = ≥ median
  - tertile: High = T3 vs Low = T1 (middle tertile dropped)
  - quartile: High = Q4 vs Low = Q1 (Q2–Q3 dropped)
  - continuous: logistic OR per 1 SD of log expression
- Endpoints (all reported when the cohort has the labels):
  - `curated_R_vs_NR`: ICB/GEO native R vs NR (for IMvigor210/Gide/Liu/Riaz/Braun this is CR/PR vs PD; SD is already excluded in the curated field)
  - `CRPR_vs_PD`: RECIST CR/PR vs PD, SD dropped
  - `CRPR_vs_SDPD`: RECIST CR/PR vs SD+PD
- OR: Woolf interval; Haldane–Anscombe 0.5 if any 2×2 cell is 0; Fisher exact p.
- Meta: DerSimonian–Laird random effects and inverse-variance fixed effect on logOR.
- GSE135222_DCB180 is the same n as GSE135222_Jung with DCB = PFS ≥ 180 days. It is **not** an extra independent study.
- GSE207422 is neoadjuvant ICI + chemotherapy, pathologic MPR vs NMPR. RECIST has no PD.
- Liu raw data are dbGaP; Braun raw WES are EGA. Processed expression+response used here are the public Zenodo ICB TSVs (CC-BY-4.0).

## Inventory

| Cohort | Cancer | Matrix | CLDN4 | Binary n (R/NR) | Notes |
|---|---|---|---|---|---|
| Braun | ccRCC | public | yes | 142 (39/103) | Zenodo ICB_Braun; ccRCC nivolumab; patients with RNA |
| GSE126044 | NSCLC | public | yes | 16 (5/11) | GEO NSCLC anti-PD-1; GEO responder/non-responder |
| GSE135222_DCB180 | NSCLC | public | yes | 27 (7/20) | SAME patients as Jung; DCB=PFS≥180d; not independent |
| GSE135222_Jung | NSCLC | public | yes | 27 (8/19) | GEO GSE135222 / ICB_Jung; curated R/NR |
| GSE166449 | NSCLC | public | yes | 22 (7/15) | GEO NSCLC ICI; title responder/non-responder |
| GSE207422 | NSCLC | public | yes | 24 (9/15) | GEO neoadjuvant ICI+chemo; MPR vs NMPR |
| Gide | melanoma | public | yes | 38 (19/19) | Zenodo ICB_Gide; melanoma anti-PD-1; n=41 public TPM |
| IMvigor210_Mariathasan | urothelial | public | yes | 235 (68/167) | Zenodo ICB_Mariathasan; atezolizumab urothelial |
| Liu | melanoma | public | yes | 112 (52/60) | Zenodo ICB_Liu processed FPKM (raw dbGaP) |
| Riaz | melanoma | public | yes | 30 (9/21) | Zenodo ICB_Riaz; melanoma nivolumab; patients with RNA |
| GSE136961 | NSCLC | public panel | **no** | — | targeted immune panel; CLDN4 not measured |
| GSE93157 | NSCLC | public panel | **no** | — | targeted immune panel; CLDN4 not measured |

Independent studies with CLDN4 + binary response: **k = 9**, not 11.

## All pooled ORs (do not cherry-pick)

Every cutoff × endpoint × subset is listed. The user number is
shown only as a reference line. `matches_user_0.42_at_2dp` is
true only when the random-effects point estimate rounds to 0.42.

| Subset | Cutoff | Endpoint | k | RE OR [95% CI] | p | I² | matches 0.42 at 2dp | 0.42 in CI |
|---|---|---|---:|---|---:|---:|---|---|
| all_independent | continuous | CRPR_vs_PD | 5 | 0.92 [0.76–1.12] | 0.427 | 5% | no | no |
| all_independent | median | CRPR_vs_PD | 5 | 0.97 [0.57–1.65] | 0.908 | 43% | no | no |
| all_independent | quartile | CRPR_vs_PD | 5 | 0.74 [0.43–1.30] | 0.297 | 7% | no | no |
| all_independent | tertile | CRPR_vs_PD | 5 | 0.85 [0.49–1.48] | 0.571 | 25% | no | no |
| all_independent | continuous | CRPR_vs_SDPD | 6 | 0.86 [0.69–1.06] | 0.16 | 26% | no | no |
| all_independent | median | CRPR_vs_SDPD | 6 | 0.71 [0.41–1.22] | 0.211 | 49% | no | yes |
| all_independent | quartile | CRPR_vs_SDPD | 6 | 0.64 [0.35–1.18] | 0.154 | 25% | no | yes |
| all_independent | tertile | CRPR_vs_SDPD | 6 | 0.55 [0.27–1.14] | 0.108 | 55% | no | yes |
| all_independent | continuous | curated_R_vs_NR | 9 | 0.92 [0.78–1.09] | 0.354 | 0% | no | no |
| all_independent | median | curated_R_vs_NR | 9 | 0.70 [0.41–1.18] | 0.176 | 44% | no | yes |
| all_independent | quartile | curated_R_vs_NR | 9 | 0.80 [0.49–1.29] | 0.352 | 0% | no | no |
| all_independent | tertile | curated_R_vs_NR | 9 | 0.79 [0.51–1.22] | 0.282 | 6% | no | no |
| lung_only | continuous | CRPR_vs_SDPD | 1 | 0.97 [0.39–2.40] | 0.951 | 0% | no | yes |
| lung_only | median | CRPR_vs_SDPD | 1 | 0.28 [0.04–1.88] | 0.19 | 0% | no | yes |
| lung_only | quartile | CRPR_vs_SDPD | 1 | 0.50 [0.05–5.15] | 0.56 | 0% | no | yes |
| lung_only | tertile | CRPR_vs_SDPD | 1 | 0.33 [0.04–2.77] | 0.309 | 0% | no | yes |
| lung_only | continuous | curated_R_vs_NR | 4 | 0.91 [0.57–1.46] | 0.692 | 0% | no | no |
| lung_only | median | curated_R_vs_NR | 4 | 0.40 [0.15–1.04] | 0.0591 | 0% | no | yes |
| lung_only | quartile | curated_R_vs_NR | 4 | 1.00 [0.25–4.01] | 0.999 | 0% | no | yes |
| lung_only | tertile | curated_R_vs_NR | 4 | 0.55 [0.17–1.75] | 0.311 | 0% | no | yes |
| non_lung | continuous | CRPR_vs_PD | 5 | 0.92 [0.76–1.12] | 0.427 | 5% | no | no |
| non_lung | median | CRPR_vs_PD | 5 | 0.97 [0.57–1.65] | 0.908 | 43% | no | no |
| non_lung | quartile | CRPR_vs_PD | 5 | 0.74 [0.43–1.30] | 0.297 | 7% | no | no |
| non_lung | tertile | CRPR_vs_PD | 5 | 0.85 [0.49–1.48] | 0.571 | 25% | no | no |
| non_lung | continuous | CRPR_vs_SDPD | 5 | 0.84 [0.66–1.07] | 0.165 | 40% | no | no |
| non_lung | median | CRPR_vs_SDPD | 5 | 0.76 [0.43–1.32] | 0.329 | 53% | no | no |
| non_lung | quartile | CRPR_vs_SDPD | 5 | 0.62 [0.31–1.25] | 0.182 | 39% | no | yes |
| non_lung | tertile | CRPR_vs_SDPD | 5 | 0.57 [0.26–1.27] | 0.169 | 63% | no | yes |
| non_lung | continuous | curated_R_vs_NR | 5 | 0.91 [0.75–1.12] | 0.377 | 11% | no | no |
| non_lung | median | curated_R_vs_NR | 5 | 0.85 [0.48–1.50] | 0.571 | 52% | no | no |
| non_lung | quartile | curated_R_vs_NR | 5 | 0.73 [0.41–1.30] | 0.291 | 14% | no | yes |
| non_lung | tertile | curated_R_vs_NR | 5 | 0.77 [0.42–1.41] | 0.395 | 38% | no | no |

## Primary locked comparison (median, curated R vs NR, all independent)

Random-effects OR = **0.70 [0.41–1.18]**, k=9, p=0.176, I²=44%.

The primary point estimate does **not** match the user 0.42 (observed 0.70).
The user 0.42 does fall inside the primary 95% CI (wide CI, not confirmation).
k=9 public independent matrices, not the claimed k=11.

## Open lung ICI only (median, curated R vs NR)

Random-effects OR = **0.40 [0.15–1.04]**, k=4, p=0.0591, I²=0%.
This lung-only pool is **not significant** (NS), consistent with the prior open-lung note.

## Per-cohort median / curated R vs NR

| Cohort | n (R/NR) | High R/NR | Low R/NR | OR [95% CI] | Fisher p | Haldane |
|---|---|---|---|---|---:|---|
| Braun | 142.0 (39/103) | 18/53 | 21/50 | 0.81 [0.39–1.69] | 0.707 | no |
| GSE126044 | 16.0 (5/11) | 1/7 | 4/4 | 0.14 [0.01–1.76] | 0.282 | no |
| GSE135222_Jung | 27.0 (8/19) | 4/10 | 4/9 | 0.90 [0.17–4.70] | 1 | no |
| GSE166449 | 22.0 (7/15) | 3/8 | 4/7 | 0.66 [0.11–4.00] | 1 | no |
| GSE207422 | 24.0 (9/15) | 2/10 | 7/5 | 0.14 [0.02–0.96] | 0.0894 | no |
| Gide | 38.0 (19/19) | 6/13 | 13/6 | 0.21 [0.05–0.84] | 0.0502 | no |
| IMvigor210_Mariathasan | 235.0 (68/167) | 40/78 | 28/89 | 1.63 [0.92–2.88] | 0.114 | no |
| Liu | 112.0 (52/60) | 25/31 | 27/29 | 0.87 [0.41–1.82] | 0.85 | no |
| Riaz | 30.0 (9/21) | 4/11 | 5/10 | 0.73 [0.15–3.49] | 1 | no |

## Honest verdict

- The user number is a **user-reported** statistic. This slice does not treat it as ground truth.
- All cutoffs are in `tables/pooled_or.csv` and `tables/per_cohort_or.csv`.
- We did not drop a cohort because CLDN4 was null or opposite.
- We did not add extra melanoma/HNSCC series after seeing results in order to reach k=11.
- Small lung n (16–27) makes lung-only ORs unstable; tertile/quartile drop the middle and are even smaller.
- Cross-cancer pooling (bladder + melanoma + RCC + lung) is not a lung-specific test.
- Neoadjuvant GSE207422 is ICI+chemo pathologic response, not metastatic ORR.
- Median splits that land everyone on one side (empty High or Low arm) are reported as not estimable, not recoded after the fact to force an OR.

## Reproduction

```bash
python3 -m pip install -r results/rework/B5_OR/requirements.txt
python3 results/rework/B5_OR/download.py
python3 results/rework/B5_OR/analyze.py
```


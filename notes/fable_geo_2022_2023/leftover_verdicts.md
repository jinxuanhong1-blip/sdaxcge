# Leftover 2022–2023 GEO lung ICI series — finish log

This pass closed every relevant series that had an empty verdict after the
first run, plus four borderline series that the first-pass text filter missed
or rejected.

## Recovered / analysed (open processed data + gene + ICI outcome)

| GSE | Why leftover | Gene(s) | Outcome used | Test (computed, not invented) |
|---|---|---|---|---|
| **GSE221733** | First-pass `is_ici=False` because the regex required `immunotherap\b`, which does not match the word *immunotherapy*. Title/summary are NSCLC immunotherapy response groups. GeoMx DSP CTA, 93 AOIs / 41 patients. | TACSTD2 yes; **CLDN4 not on CTA panel** | GEO `response` (Responder / Non-responder); `followup` + `status` | Patient-level mean PanCK+ TACSTD2: R n=15 median 9.268 vs NR n=19 median 9.378; MWU p=**0.5324**, Cliff's δ=−0.13. PanCK− p=0.8067. Median-split OS log-rank p=**0.1545**. |
| **GSE248378** | Relevant but first pass stopped at “treatment = Arm1/Arm2”. Post-Rx FPKM has TACSTD2+CLDN4. | both | Published Nat Commun source data (MOESM6 Fig. 1b) joined on DurvaNNN = numeric GEO title. All 29 GEO samples are non-MPR (`Path Response_2Group=None`). Recurrence from paper Status. **45-M-PO**: GEO Arm2 vs paper Durvalumab monotherapy — kept for the n=29 paper set; dropped in a sensitivity test. | TACSTD2 no-recurrence n=20 median 41.98 FPKM vs recurrence n=9 median 80.01; MWU p=**0.0562**, δ=−0.456. CLDN4 p=**0.3108**, δ=−0.244. Sensitivity excluding 45-M-PO: TACSTD2 p=0.0766. Median-split DFS log-rank p=**0.1220**. |
| **GSE193049** | Relevant leftover; only RAW per-GSM CSVs. BALF, not tumour. | both | Sample titles `responder_*` / `non_responder_*` (4 vs 3) | TACSTD2 p=**0.6286**, δ=+0.333; CLDN4 p=**0.4000**, δ=+0.50. Opposite direction to tumour cohorts; n=7. |

## Borderline series inspected and rejected for this question

| GSE | n | Why it cannot test TACSTD2/CLDN4 vs ICI outcome |
|---|---|---|
| GSE221322 | 96 | Sister DSP **protein** panel; no TACSTD2/CLDN4 proteins. |
| GSE189045 | 23 | Serum **exosomal miRNA** with response labels; no mRNA. |
| GSE250262 | 51 | NanoString IO360 RCC of NSCLC tumours; panel has EPCAM, **not** TACSTD2/CLDN4; GEO has no response/survival field. |

## All other leftover relevant series (deep verdict: not usable)

Cell line / organoid / ATAC / infection / IPF / platelet-adjacent / sorted-immune / scRNA — no bulk tumour epithelium + per-patient ICI outcome + both genes:

GSE185204, GSE186446, GSE235048, GSE228419, GSE164146, GSE212622, GSE224099,
GSE193719, GSE189804, GSE217451, GSE224216, GSE150255, GSE194350, GSE195770,
GSE218402, GSE224246, GSE238006, GSE250254, GSE229353, GSE178521, GSE192591,
GSE192790, GSE197236, GSE198099, GSE213590, GSE213902, GSE223779, GSE193707.

One-line verdicts are in `results/fable_geo_2022_2023/geo_lung_ici_catalog.csv`.

## What was *not* done (anti-fabrication)

- No MPR contrast was computed for GSE248378: the deposited matrix is the paper’s non-MPR post-treatment set (29/29 `None`).
- No outcome labels were guessed from treatment arm.
- GSE248378 45-M-PO was not “corrected”; it is reported as arm-discordant and used only in the n=29 paper-faithful recurrence test.
- GSE221733 CLDN4 was not imputed.
- GSE193049 is reported as BALF, not tumour.

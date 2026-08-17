# FINDING — merged GSE131907+GSE205335 CLDN4-only SCENIC/GRN proxy

Additive **CLDN4-only** GRN on the **merged GSE131907 + GSE205335** author-malignant
combo that already differs vs T/NK (PR #320). **A10 ELF3–CLDN4 is taken as given**
and is not re-tested as a discovery. No dual-high TACSTD2×CLDN4. No GSE207422.

**Question:** do CLDN4-high malignant cells show different IFN / MHC-I / TJ / keratin
regulons vs CLDN4-low, **at the patient level**?

## Method (what was actually run)

Full **pySCENIC cisTarget was not run** (no motif ranking databases in this
environment; `pyscenic` importable=False).
No ChIP peaks were invented.

This is a documented **AUCell + public TF–target prior** proxy:

1. Public curated edges: TRRUST v2 + DoRothEA + CollecTRI (OmniPath). Not binding.
2. Program gene sets from the PR #267 A8 freeze (MSigDB Hallmark IFN-α/IFN-γ,
   GO tight-junction / keratinization, custom MHC-I antigen presentation,
   compact KRT panel). **CLDN4 is held out** of every set (it is the split gene).
3. Aibar-style AUCell on the **extracted gene universe** (program genes + TFs +
   prior targets), `auc_threshold=0.05`. Not R AUCell binary.
4. Pearson TF–program co-expression in CLDN4-high cells is descriptive only.

Patient is the unit. Cells are split **within patient** at the median of malignant
CLDN4 log1p(CP10k). If that median is 0 (zero-inflated), high = CLDN4>0 vs low = 0
so CLDN4-low patients are not dropped. Paired Wilcoxon on patient-mean AUCell
(high − low). Cohorts are scored separately, then patient deltas are stacked
(not Harmony). p-values are descriptive.

## Honest n

Primary paired test requires ≥20 CLDN4-high **and** ≥20 CLDN4-low
author-malignant cells after UMI≥200. That n is **patients**, not cells.

| item | n | note |
| --- | ---: | --- |
| GSE131907_author_malignant_cells | 31136 | Cell_subtype in {Malignant cells, tS1, tS2, tS3} |
| GSE131907_malignant_UMI>=200 | 31136 | cells |
| GSE131907_patients_ge20_malignant | 31 | patient_id; unit |
| GSE131907_patients_paired_ge20_high_and_low | 29 | primary paired AUCell unit |
| GSE205335_author_malignant_cells | 28512 | lineage.sub == Malignant cells |
| GSE205335_malignant_UMI>=200 | 28512 | cells |
| GSE205335_patients_ge20_malignant | 22 | SOFT patient; unit |
| GSE205335_patients_paired_ge20_high_and_low | 21 | primary paired AUCell unit |
| merged_patients_ge20_malignant | 53 | k=2 cohorts stacked; not cell-integrated |
| merged_patients_paired | 50 | honest n for high vs low regulon test |
| CLDN4_high_cells_in_paired_patients | 29286 | cells; not the claim |
| CLDN4_low_cells_in_paired_patients | 29875 | cells; not the claim |
| GSE131907_tLung_patients_paired | 10 | tS1/tS2/tS3 primary-tumor sensitivity |
| universe_genes_AUCell | 2020 | extracted programs + TFs + public prior targets |

**Paired unit = 50 patients** (29 GSE131907 + 21 GSE205335).
GSE131907 tLung author labels are tS1/tS2/tS3 (included). PR #320's T/NK extract
used `Malignant cells` only and therefore dropped tLung — that subset is a
sensitivity, not the GRN definition. GSE205335 Q4 in PR #320 is SCLC-heavy;
histology is reported, not hidden. Cell n is labeled as cells (pseudoreplication).

## Program AUCell (the question)

Patient-paired median Δ = CLDN4-high − CLDN4-low. Positive = higher in CLDN4-high.

| regulon | tf | program | kind | n_targets | elf3_given | n_patients_paired | n_GSE131907 | n_GSE205335 | delta_median_high_minus_low | p | fdr | delta_GSE131907 | p_GSE131907 | delta_GSE205335 | p_GSE205335 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| IFN |  | IFN | program_set | 222 | False | 50 | 29 | 21 | -0.000278 | 0.3326 | 0.3959 | -0.000177 | 0.7493 | -0.00098 | 0.3554 |
| IFN_IFNA |  | IFN_IFNA | program_set | 96 | False | 50 | 29 | 21 | -0.00095 | 0.0858 | 0.1129 | -0.000496 | 0.4420 | -0.00215 | 0.1032 |
| IFN_IFNG |  | IFN_IFNG | program_set | 198 | False | 50 | 29 | 21 | -0.000689 | 0.1139 | 0.1461 | -0.000545 | 0.3357 | -0.00108 | 0.2428 |
| MHC_I |  | MHC_I | program_set | 21 | False | 50 | 29 | 21 | 0.00623 | 0.0028 | 0.0042 | 0.00641 | 0.0075 | 0.00347 | 0.1111 |
| TJ |  | TJ | program_set | 80 | False | 50 | 29 | 21 | 0.00112 | 1.94e-04 | 3.35e-04 | 0.000631 | 0.0203 | 0.00133 | 0.0043 |
| keratin |  | keratin | program_set | 109 | False | 50 | 29 | 21 | 0.00199 | 1.91e-08 | 6.36e-08 | 0.00159 | 4.70e-05 | 0.00643 | 2.41e-04 |
| APICAL_JUNCTION |  | APICAL_JUNCTION | program_set | 199 | False | 50 | 29 | 21 | -0.00129 | 4.79e-10 | 2.40e-09 | -0.00108 | 1.10e-04 | -0.00289 | 1.91e-06 |

Stacked-patient program row: IFN Δ=-0.000278 p=0.3326 (n=50); MHC_I Δ=0.00623 p=0.0028 (n=50); TJ Δ=0.00112 p=1.94e-04 (n=50); keratin Δ=0.00199 p=1.91e-08 (n=50).

**Answer (descriptive, patient-paired):** IFN program AUCell does **not** differ
CLDN4-high vs low. Compact **TJ** and **keratin** are higher in CLDN4-high cells
in both cohorts. **MHC-I** is higher on the stacked n (same sign in both; GSE205335
alone is weaker). Hallmark **APICAL_JUNCTION** (200-gene mixed set) goes the other
way — it is not the compact TJ set. ELF3 is A10-given and is not the headline.
Between-patient CLDN4 %pos vs IFN (the PR #320 axis) is a different question;
GSE205335 leans IFN-low in CLDN4-high *patients*, with I².

Thin prior∩program intersections (n_targets < 8) are listed in `regulons.tsv`
but are not a binding claim. NLRC5 has no public prior edges in the snapshot.
Three patients fail the ≥20/20 tail rule (P1013 1 CLDN4+ cell; P3016 9 high;
P4001 13/14 of 27 cells) and are out of the paired n.

## TF prior AUCell (IFN / MHC-I / TJ / keratin TFs)

ELF3 rows are **A10-given** (`elf3_given=True`) and are not a discovery.

| regulon | tf | program | kind | n_targets | elf3_given | n_patients_paired | n_GSE131907 | n_GSE205335 | delta_median_high_minus_low | p | fdr | delta_GSE131907 | p_GSE131907 | delta_GSE205335 | p_GSE205335 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| STAT1_prior | STAT1 | IFN | public_prior | 856 | False | 50 | 29 | 21 | 1.73e-05 | 0.7090 | 0.7542 | 7.18e-05 | 0.0877 | -0.00012 | 0.2029 |
| STAT1_priorANDprog | STAT1 | IFN | prior_AND_program | 90 | False | 50 | 29 | 21 | -0.000614 | 0.2080 | 0.2599 | 8.51e-05 | 0.9152 | -0.00176 | 0.1193 |
| STAT2_prior | STAT2 | IFN | public_prior | 79 | False | 50 | 29 | 21 | 0.00247 | 9.75e-08 | 2.71e-07 | 0.00287 | 2.05e-07 | 0.00147 | 0.0175 |
| STAT2_priorANDprog | STAT2 | IFN | prior_AND_program | 36 | False | 50 | 29 | 21 | -0.00196 | 8.81e-04 | 0.0014 | -0.00143 | 0.0689 | -0.00405 | 0.0080 |
| IRF1_prior | IRF1 | IFN_MHC | public_prior | 169 | False | 50 | 29 | 21 | -0.000157 | 0.6947 | 0.7542 | -0.000179 | 1.0000 | -0.000134 | 0.6333 |
| IRF1_priorANDprog | IRF1 | IFN_MHC | prior_AND_program | 48 | False | 50 | 29 | 21 | -0.000149 | 0.9542 | 0.9542 | -0.000138 | 0.5362 | -0.000391 | 0.4733 |
| IRF7_prior | IRF7 | IFN | public_prior | 28 | False | 50 | 29 | 21 | -0.00107 | 8.81e-04 | 0.0014 | -0.00111 | 0.0190 | -0.000948 | 0.0263 |
| IRF7_priorANDprog | IRF7 | IFN | prior_AND_program | 5 | False | 50 | 29 | 21 | -0.00134 | 0.0663 | 0.0921 | -0.000785 | 0.2383 | -0.00211 | 0.1790 |
| IRF9_prior | IRF9 | IFN | public_prior | 35 | False | 50 | 29 | 21 | -0.000103 | 0.5985 | 0.6650 | 0.000506 | 0.2746 | -0.00191 | 0.0646 |
| IRF9_priorANDprog | IRF9 | IFN | prior_AND_program | 16 | False | 50 | 29 | 21 | -0.000551 | 0.3326 | 0.3959 | 0.000227 | 0.7013 | -0.00228 | 0.0958 |
| RFX5_prior | RFX5 | MHC_I | public_prior | 31 | False | 50 | 29 | 21 | -0.00395 | 0.0013 | 0.0020 | -0.0037 | 0.0258 | -0.00451 | 0.0239 |
| RFX5_priorANDprog | RFX5 | MHC_I | prior_AND_program | 4 | False | 50 | 29 | 21 | 5.51e-05 | 0.8109 | 0.8447 | -0.00106 | 0.9830 | 0.00255 | 0.8382 |
| GRHL2_prior | GRHL2 | TJ | public_prior | 10 | False | 50 | 29 | 21 | 0.0101 | 8.86e-10 | 3.69e-09 | 0.00701 | 4.70e-06 | 0.0188 | 1.05e-04 |
| GRHL2_priorANDprog | GRHL2 | TJ | prior_AND_program | 1 | False | 50 | 29 | 21 | 0.0824 | 1.14e-12 | 8.12e-12 | 0.0726 | 4.10e-07 | 0.168 | 1.91e-06 |
| OVOL1_prior | OVOL1 | TJ | public_prior | 9 | False | 50 | 29 | 21 | 0.00284 | 3.13e-05 | 5.79e-05 | 0.00333 | 6.07e-04 | 0.00211 | 0.0175 |
| OVOL2_prior | OVOL2 | TJ | public_prior | 1 | False | 50 | 29 | 21 | -5.09e-05 | 0.0056 | 0.0083 | 0 | 0.7946 | -0.00263 | 4.49e-04 |
| KLF5_prior | KLF5 | TJ | public_prior | 79 | False | 50 | 29 | 21 | 0.00116 | 1.49e-04 | 2.66e-04 | 0.00116 | 8.88e-04 | 0.00108 | 0.0319 |
| ELF3_prior | ELF3 | TJ_GIVEN | public_prior | 33 | True | 50 | 29 | 21 | 0.00387 | 1.29e-10 | 7.15e-10 | 0.00281 | 5.34e-05 | 0.00919 | 9.54e-07 |
| ELF3_priorANDprog | ELF3 | TJ_GIVEN | prior_AND_program | 1 | True | 50 | 29 | 21 | 0.0316 | 4.86e-09 | 1.87e-08 | 0.0259 | 2.43e-05 | 0.0513 | 1.05e-04 |
| TP63_prior | TP63 | keratin | public_prior | 207 | False | 50 | 29 | 21 | -0.00095 | 8.40e-06 | 1.68e-05 | -0.000965 | 1.56e-04 | -0.000934 | 0.0090 |
| TP63_priorANDprog | TP63 | keratin | prior_AND_program | 8 | False | 50 | 29 | 21 | 0.000137 | 0.3673 | 0.4271 | 0.000228 | 0.6545 | 4.63e-05 | 0.4948 |
| GRHL1_prior | GRHL1 | keratin | public_prior | 2 | False | 50 | 29 | 21 | 0 | 0.8652 | 0.8829 | 0 | 0.6874 | 0 | 0.4691 |
| KLF4_prior | KLF4 | keratin | public_prior | 131 | False | 50 | 29 | 21 | 0.00227 | 1.40e-07 | 3.51e-07 | 0.0013 | 6.07e-04 | 0.0045 | 1.05e-04 |
| KLF4_priorANDprog | KLF4 | keratin | prior_AND_program | 4 | False | 50 | 29 | 21 | 0.0243 | 4.44e-14 | 5.55e-13 | 0.0196 | 2.61e-08 | 0.0417 | 1.91e-06 |
| NKX2-1_prior | NKX2-1 | control | public_prior | 63 | False | 50 | 29 | 21 | -0.000797 | 0.4376 | 0.4972 | -0.00127 | 0.0920 | -0.000321 | 0.4524 |
| SOX2_prior | SOX2 | control | public_prior | 343 | False | 50 | 29 | 21 | -0.00101 | 1.37e-08 | 4.89e-08 | -0.000927 | 5.22e-08 | -0.00176 | 0.0049 |

Target lists are in [`tables/regulons.tsv`](tables/regulons.tsv)
and [`results/merge_131907_205335_scenic_cldn4/regulons.tsv`](../../results/merge_131907_205335_scenic_cldn4/regulons.tsv).

## Between-patient companion (not the cell-split)

Spearman of patient-level malignant CLDN4 %pos vs patient-mean program AUCell
(all cells in the patient; no high/low split). Fisher-z pool of the two cohorts.
This asks whether **CLDN4-high patients** (the PR #320 axis) also have higher
IFN/MHC-I/TJ/keratin programs — a different question from the within-patient split.

| regulon | n | k | rho | p | I2 | n_GSE131907 | rho_GSE131907 | p_GSE131907 | n_GSE205335 | rho_GSE205335 | p_GSE205335 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| IFN | 53 | 2 | -0.26 | 0.0685 | 61 | 31 | -0.0734 | 0.6948 | 22 | -0.5 | 0.0179 |
| IFN_IFNA | 53 | 2 | -0.147 | 0.3092 | 64 | 31 | 0.0516 | 0.7827 | 22 | -0.416 | 0.0541 |
| IFN_IFNG | 53 | 2 | -0.302 | 0.0328 | 37 | 31 | -0.159 | 0.3933 | 22 | -0.488 | 0.0211 |
| MHC_I | 53 | 2 | -0.0596 | 0.6826 | 73 | 31 | 0.171 | 0.3589 | 22 | -0.381 | 0.0801 |
| TJ | 53 | 2 | 0.152 | 0.2923 | 63 | 31 | -0.0427 | 0.8194 | 22 | 0.416 | 0.0541 |
| keratin | 53 | 2 | 0.0631 | 0.6648 | 23 | 31 | 0.197 | 0.2877 | 22 | -0.137 | 0.5426 |
| APICAL_JUNCTION | 53 | 2 | -0.21 | 0.1437 | 0 | 31 | -0.21 | 0.2557 | 22 | -0.209 | 0.3494 |

## What is / is not claimed

- **Claimed:** a descriptive patient-level regulon / program-AUCell table on
  author-malignant GSE131907+GSE205335 cells, with honest patient n.
- **Not claimed:** TF binding, pySCENIC cisTarget, or a new ELF3–CLDN4 discovery.
- **Not claimed:** dual-high TACSTD2×CLDN4, or any GSE207422-only result.
- **Not claimed:** cell-level p-values as the finding.
- **Not claimed:** ICI / MPR from GSE131907 (treatment-naive). GSE205335 has
  RECIST, not MPR; RECIST is not used as MPR.

Method flags: pyscenic_importable=False;
cisTarget_run=False;
chip_peaks_invented=False.

## Extra figures

- `figures/fig_program_paired_forest.png` — patient-paired program Δ
- `figures/fig_tf_prior_paired_forest.png` — TF-prior AUCell Δ (ELF3 marked given)
- `figures/fig_patient_delta_heatmap.png` — per-patient program Δ
- `figures/fig_paired_boxes.png` — IFN / MHC-I / TJ / keratin paired boxes
- `figures/fig_between_patient_scatter.png` — CLDN4 %pos vs program AUCell
- `figures/fig_n_patients.png` — honest n
- `figures/fig_high_specific_counts.png` — Pearson high-specific TF–program edges

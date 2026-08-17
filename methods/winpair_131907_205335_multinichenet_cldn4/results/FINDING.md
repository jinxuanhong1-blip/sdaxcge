# FINDING — MultiNicheNet / multi-sample NicheNet on winning pair GSE131907+GSE205335

**Verdict:** On the winning pair (eligible n=**51**: 29 GSE131907 + 22 GSE205335), same-patient T/NK **fraction** falls as malignant CLDN4 rises (mean ρ=-0.383 p=0.00551 (n=51); %pos ρ=-0.530 p=6.24e-05 (n=51)). T/NK IFN / cytotoxicity vs **mean** CLDN4 are the same sign but NS (IFN ρ=-0.251 p=0.0753 (n=51); cytotoxicity ρ=-0.241 p=0.0886 (n=51)). Vs CLDN4 **%pos** they are significant (IFN ρ=-0.365 p=0.00851 (n=51); cytotoxicity ρ=-0.358 p=0.0098 (n=51)). The ligand-activity table is a NicheNet-v2 prior ranking plus patient-paired CLDN4-high vs low ligand DE. Prior recovery of an IFN list is **not** evidence that CLDN4-high cells induce or repress IFN.

ADDITIVE. **CLDN4 only.** No dual-high. **GSE207422 was not used** (NicheNet #334 NS). Sender = CLDN4-high vs low malignant. Receiver = **same-patient T/NK**. **Patient is the unit.** `nichenetr` / `multinichenetr` were not installed; ligand activity is the published NicheNet-v2 ligand–target prior (Zenodo 7074291).

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| datasets | 2 | GSE131907 + GSE205335 only; GSE207422 excluded |
| GSE131907_cells | 208506 | Kim et al. author barcodes |
| GSE131907_malignant | 31136 | tumor-origin ∩ {Malignant cells, tS1, tS2, tS3} |
| GSE131907_cldn4_high_low | 15569 | high; low=15567; median=1.798 |
| GSE131907_TNK | 48012 | tumor-origin T lymphocytes + NK |
| GSE205335_cells | 96505 | Hu/Ahn/Lee processed UMI |
| GSE205335_malignant | 28512 | tumor sample ∩ lineage.sub == Malignant cells |
| GSE205335_cldn4_high_low | 14256 | high; low=14256; median=1.371 |
| GSE205335_TNK | 38265 | tumor sample ∩ lineage.total == T/NK cells |
| eligible_patients_GSE131907 | 29 | ≥10 high, ≥10 low, ≥20 T/NK; patient pooled |
| eligible_patients_GSE205335 | 22 | same gate; RECIST not used as MPR |
| eligible_patients_combined | 51 | unit of inference |
| potential_ligands_combined | 136 | CLDN4-high ligand + T/NK receptor + v2 prior |
| potential_ligands_GSE131907 | 114 | dataset-specific 10% rule |
| potential_ligands_GSE205335 | 161 | dataset-specific 10% rule |
| background_genes | 255 | T/NK-expressed ∩ prior targets (panel-restricted) |
| a_priori_ifn_in_prior | 16 |  |
| a_priori_cyto_in_prior | 14 |  |
| empirical_TNK_up_with_CLDN4 | 170 | patient Spearman p<0.15, ρ>0 |
| empirical_TNK_down_with_CLDN4 | 203 | patient Spearman p<0.15, ρ<0 |

GSE131907 is treatment-naive (Kim et al., *Nat Commun* 2020). It has **no ICI / MPR**. GSE205335 is palliative ICI biopsy/effusion (Hu / Ahn / Lee). RECIST is recorded and **not** used as MPR. Tumor samples are pooled per patient. A patient enters the paired ligand DE if it has ≥10 CLDN4-high, ≥10 CLDN4-low malignant cells and ≥20 T/NK.

### E1 / E2 — patient T/NK programs vs malignant CLDN4

Unit = eligible patient. Cells are not n. Q4 vs Q1 is the quartile tails of malignant CLDN4 mean (n_compared = n_Q1 + n_Q4).

| Dataset | Test | n | Result |
| --- | --- | ---: | --- |
| GSE131907 | spearman_mal_CLDN4_mean__tnk_ifn | 29 | ρ=-0.297 p=0.118 |
| GSE131907 | spearman_mal_CLDN4_mean__tnk_cyto | 29 | ρ=-0.309 p=0.102 |
| GSE131907 | spearman_mal_CLDN4_mean__tnk_ifn_minus_cyto | 29 | ρ=+0.270 p=0.156 |
| GSE131907 | spearman_mal_CLDN4_mean__tnk_exh | 29 | ρ=+0.022 p=0.909 |
| GSE131907 | spearman_mal_CLDN4_mean__tnk_frac | 29 | ρ=-0.576 p=0.00107 |
| GSE131907 | spearman_mal_CLDN4_pct_pos__tnk_ifn | 29 | ρ=-0.172 p=0.371 |
| GSE131907 | spearman_mal_CLDN4_pct_pos__tnk_cyto | 29 | ρ=-0.212 p=0.269 |
| GSE131907 | spearman_mal_CLDN4_pct_pos__tnk_frac | 29 | ρ=-0.552 p=0.0019 |
| GSE131907 | Q4_vs_Q1_tnk_ifn | 15 | Q1 mean +0.252 vs Q4 +0.209; p=0.281 |
| GSE131907 | Q4_vs_Q1_tnk_cyto | 15 | Q1 mean +0.897 vs Q4 +0.616; p=0.0939 |
| GSE131907 | Q4_vs_Q1_tnk_frac | 15 | Q1 mean +0.383 vs Q4 +0.124; p=0.00932 |
| GSE205335 | spearman_mal_CLDN4_mean__tnk_ifn | 22 | ρ=-0.155 p=0.49 |
| GSE205335 | spearman_mal_CLDN4_mean__tnk_cyto | 22 | ρ=-0.219 p=0.329 |
| GSE205335 | spearman_mal_CLDN4_mean__tnk_ifn_minus_cyto | 22 | ρ=+0.213 p=0.342 |
| GSE205335 | spearman_mal_CLDN4_mean__tnk_exh | 22 | ρ=-0.209 p=0.349 |
| GSE205335 | spearman_mal_CLDN4_mean__tnk_frac | 22 | ρ=-0.068 p=0.763 |
| GSE205335 | spearman_mal_CLDN4_pct_pos__tnk_ifn | 22 | ρ=-0.290 p=0.191 |
| GSE205335 | spearman_mal_CLDN4_pct_pos__tnk_cyto | 22 | ρ=-0.411 p=0.0577 |
| GSE205335 | spearman_mal_CLDN4_pct_pos__tnk_frac | 22 | ρ=-0.316 p=0.152 |
| GSE205335 | Q4_vs_Q1_tnk_ifn | 12 | Q1 mean +0.333 vs Q4 +0.245; p=0.394 |
| GSE205335 | Q4_vs_Q1_tnk_cyto | 12 | Q1 mean +0.899 vs Q4 +0.686; p=0.589 |
| GSE205335 | Q4_vs_Q1_tnk_frac | 12 | Q1 mean +0.368 vs Q4 +0.289; p=0.24 |
| combined | spearman_mal_CLDN4_mean__tnk_ifn | 51 | ρ=-0.251 p=0.0753 |
| combined | spearman_mal_CLDN4_mean__tnk_cyto | 51 | ρ=-0.241 p=0.0886 |
| combined | spearman_mal_CLDN4_mean__tnk_ifn_minus_cyto | 51 | ρ=+0.206 p=0.147 |
| combined | spearman_mal_CLDN4_mean__tnk_exh | 51 | ρ=-0.119 p=0.406 |
| combined | spearman_mal_CLDN4_mean__tnk_frac | 51 | ρ=-0.383 p=0.00551 |
| combined | spearman_mal_CLDN4_pct_pos__tnk_ifn | 51 | ρ=-0.365 p=0.00851 |
| combined | spearman_mal_CLDN4_pct_pos__tnk_cyto | 51 | ρ=-0.358 p=0.0098 |
| combined | spearman_mal_CLDN4_pct_pos__tnk_frac | 51 | ρ=-0.530 p=6.24e-05 |
| combined | Q4_vs_Q1_tnk_ifn | 26 | Q1 mean +0.264 vs Q4 +0.239; p=0.644 |
| combined | Q4_vs_Q1_tnk_cyto | 26 | Q1 mean +0.892 vs Q4 +0.704; p=0.166 |
| combined | Q4_vs_Q1_tnk_frac | 26 | Q1 mean +0.417 vs Q4 +0.235; p=0.021 |

### E3 — ligand activity (unsigned NicheNet-v2 prior)

Primary sender = CLDN4-high malignant (CLDN4 only); receiver = same-patient T/NK. Rank by Pearson of prior target scores vs gene-set membership on T/NK-expressed background genes (panel ∩ prior; **not** a full-transcriptome `nichenetr` run).

Full table: [`results/ligand_activity_table.tsv`](results/ligand_activity_table.tsv) (also `ligand_activity_all.tsv`).

**A priori IFN. Top 8 by Pearson (combined potential ligands):**

| rank_activity | ligand | pearson | auroc | pearson_p | n_patients_de |
|---|---|---|---|---|---|
| 1 | IFITM1 | 0.551 | 0.944 | 1.18e-21 | 51 |
| 2 | LIF | 0.472 | 0.905 | 1.61e-15 | 51 |
| 3 | VSIG10 | 0.451 | 0.898 | 3.43e-14 | 51 |
| 4 | HLA-F | 0.451 | 0.896 | 3.65e-14 | 51 |
| 5 | CRLF1 | 0.446 | 0.897 | 7.43e-14 | 51 |
| 6 | CLCF1 | 0.440 | 0.902 | 1.58e-13 | 51 |
| 7 | CD40 | 0.430 | 0.907 | 6.59e-13 | 51 |
| 8 | HLA-DRB5 | 0.423 | 0.891 | 1.68e-12 | 51 |

**MultiNicheNet-style prioritization (IFN set).** Combines scaled ligand activity, patient-paired CLDN4-high vs low ligand logFC, and T/NK receptor coverage. Top 12:

| rank_priority | ligand | pearson | auroc | mean_delta_high_minus_low | paired_p | frac_patients_receptor | prioritization |
|---|---|---|---|---|---|---|---|
| 1 | MDK | 0.354 | 0.851 | +0.242 | 8.13e-06 | 0.98 | 0.829 |
| 2 | HLA-B | 0.404 | 0.880 | +0.218 | 3.57e-05 | 0.94 | 0.819 |
| 3 | LGALS3 | 0.222 | 0.787 | +0.291 | 6.14e-07 | 1.00 | 0.804 |
| 4 | HLA-E | 0.376 | 0.865 | +0.193 | 4.58e-07 | 0.96 | 0.790 |
| 5 | HLA-A | 0.258 | 0.862 | +0.243 | 5.23e-06 | 1.00 | 0.784 |
| 6 | NECTIN2 | 0.217 | 0.778 | +0.275 | 9.54e-07 | 0.96 | 0.768 |
| 7 | CCN1 | 0.218 | 0.767 | +0.251 | 4.77e-06 | 1.00 | 0.766 |
| 8 | LGALS3BP | 0.155 | 0.793 | +0.290 | 1.05e-09 | 0.98 | 0.754 |
| 9 | CDH1 | 0.212 | 0.780 | +0.232 | 1.49e-09 | 1.00 | 0.746 |
| 10 | F11R | 0.215 | 0.754 | +0.229 | 7.35e-10 | 1.00 | 0.745 |
| 11 | CD55 | 0.222 | 0.794 | +0.238 | 3.56e-08 | 0.96 | 0.738 |
| 12 | LAMB3 | 0.229 | 0.770 | +0.187 | 8.65e-09 | 1.00 | 0.717 |

**A priori cytotoxicity. Top 12 by prioritization:**

| rank_priority | ligand | pearson | auroc | mean_delta_high_minus_low | paired_p | prioritization |
|---|---|---|---|---|---|---|
| 1 | HLA-A | 0.241 | 0.476 | +0.243 | 5.23e-06 | 0.897 |
| 2 | ICAM1 | 0.114 | 0.360 | +0.200 | 7.32e-09 | 0.766 |
| 3 | APP | 0.082 | 0.334 | +0.190 | 1.12e-07 | 0.734 |
| 4 | HLA-E | 0.082 | 0.443 | +0.193 | 4.58e-07 | 0.716 |
| 5 | LGALS3 | -0.067 | 0.336 | +0.291 | 6.14e-07 | 0.713 |
| 6 | CCN1 | -0.029 | 0.393 | +0.251 | 4.77e-06 | 0.705 |
| 7 | HLA-DRA | 0.297 | 0.400 | -0.026 | 0.311 | 0.701 |
| 8 | LGALS3BP | -0.073 | 0.330 | +0.290 | 1.05e-09 | 0.697 |
| 9 | MDK | -0.016 | 0.408 | +0.242 | 8.13e-06 | 0.697 |
| 10 | B2M | 0.110 | 0.459 | +0.133 | 0.000819 | 0.694 |
| 11 | HLA-B | 0.034 | 0.420 | +0.218 | 3.57e-05 | 0.694 |
| 12 | F11R | -0.028 | 0.358 | +0.229 | 7.35e-10 | 0.686 |

A high Pearson on the IFN list means the ligand sits next to IFN genes in the published prior (ISG / MHC ligands). It does **not** mean IFN is up or down in CLDN4-high patients. Direction is E1/E2 only.

### E4 — multi-sample (patient-paired) ligand DE

For each potential ligand, mean `log1p(CP10k)` in CLDN4-high vs CLDN4-low malignant cells is compared with a **paired Wilcoxon across patients** (combined eligible set, and each dataset alone). That is the MultiNicheNet sender-DE term. Receivers are the same patient's T/NK, so receiver gene-set membership does not split by sender state within a patient; receiver association is the between-patient Spearman of T/NK gene means vs malignant CLDN4 (empirical sets).

### What this is not

- Not GSE207422 and not a re-run of PR #334.
- Not TACSTD2∩CLDN4 dual-high senders.
- Not `nichenetr` / `multinichenetr` R (packages absent; documented v2 prior used).
- Not CellChat / LIANA communication probability.
- Not inferCNV/CopyKAT recomputed malignant IDs (author labels).
- Not ICI / MPR on GSE131907 (none labeled). RECIST on GSE205335 is not MPR.
- Cells are not the sample size.

### Files

`results/ligand_activity_table.tsv`, `ligand_activity_all.tsv`, `top_ligands_primary.tsv`, `ligand_de_paired.tsv`, `patient_level_tests.tsv`, `n_table.tsv`, `sample_metrics.tsv`, `fig1`–`fig7`.

---

## 中文

**结论：** On the winning pair (eligible n=**51**: 29 GSE131907 + 22 GSE205335), same-patient T/NK **fraction** falls as malignant CLDN4 rises (mean ρ=-0.383 p=0.00551 (n=51); %pos ρ=-0.530 p=6.24e-05 (n=51)). T/NK IFN / cytotoxicity vs **mean** CLDN4 are the same sign but NS (IFN ρ=-0.251 p=0.0753 (n=51); cytotoxicity ρ=-0.241 p=0.0886 (n=51)). Vs CLDN4 **%pos** they are significant (IFN ρ=-0.365 p=0.00851 (n=51); cytotoxicity ρ=-0.358 p=0.0098 (n=51)). The ligand-activity table is a NicheNet-v2 prior ranking plus patient-paired CLDN4-high vs low ligand DE. Prior recovery of an IFN list is **not** evidence that CLDN4-high cells induce or repress IFN.

只做胜出对 **GSE131907 + GSE205335**。**不用 GSE207422**（#334 NS）。发送端只用 CLDN4（恶性细胞全局中位数拆高/低），不用 TACSTD2 双高。接收端是**同一患者** T/NK。推断单位是患者：**合格 n=51（131907: 29；205335: 22）**。
配体活性是公开 NicheNet-v2 先验（Pearson / AUROC），加上患者配对的高 vs 低配体 DE （MultiNicheNet 式排序）。先验回收 IFN 基因集不是本队列里 CLDN4 诱导/抑制 IFN 的证据。GSE131907 无 ICI/MPR；GSE205335 的 RECIST 不代替 MPR。

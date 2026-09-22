# Public integrate cohorts: maximum |ρ| for mouse-level Tacstd2%

Public processed counts only: GSE154977 (KP, AT2-lineage FACS), GSE180963 (1 K + 1 KL, mixed), GSE165641 (2 KL, mixed). Private 8 KL matrices were not read and were not merged. Unit is the mouse. Tacstd2 is not an epithelial caller.

This file reports the maximum absolute Spearman correlation on a grid that was fixed in `scripts/analyze.py` before ranking. It is the extreme of that search, not a single prespecified test. Eligible tests need at least 4 mice. T-fraction tests that still contain a GSE154977 FACS library are a design no-go and are absent from the T argmax. Mixed-digest gates that call more than half of a mouse's cells epithelial drop that mouse.

Locked-gate replay against the integrate mouse table checks n_epi and T/NK fraction (same epithelial caller as the Cldn4 integrate). Tacstd2% is the new exposure; Cldn4% is recorded only for gate diagnostics. Epithelial IFN means may differ slightly from the integrate object because this run log-normalizes with each cell's full UMI total (median my/published IFN = 0.986).

Eligible tests searched: T_frac 1800, IFN/APM 10740, MHC (side, not in the IFN/APM argmax) 3580.

Epithelial gates that never produce an eligible T-fraction test, because GSE180963 then has fewer than 10 epithelial cells: epcam_krt_ptprc_neg_sftpc_neg, epcam_ptprc_neg_sftpc_neg, episcore_ge0.5_ptprc_neg_sftpc_neg, krt19_epcam_ptprc_neg_sftpc_neg. Those gates require Sftpc-negative epithelium. The T-fraction ceiling below is carried by gates that still keep Sftpc-positive cells.

## T fraction

Maximum |ρ| in T_frac: **n=4, ρ=-1.000, exact p=0.0833 (2/24)**.

Shown specification (study `drop_GSE154977`, genotype `all`, treatment `all`, epithelial gate `locked`, min epithelial cells `10`, Tacstd2% `pct_ge5`, endpoint `frac_tnk`). Dataset-residual Spearman ρ=-1.000 (asymptotic p=t approximation saturates at |ρ|=1). Genotype-residual Spearman ρ=-1.000 (asymptotic p=t approximation saturates at |ρ|=1).

Mice: GSE165641_KL1,GSE165641_KL2,GSE180963_K,GSE180963_KL.

Eligible tests in this family: 1800. Tests on the |ρ|=1 ceiling: 1104 (1104 negative, 0 positive). The row above is the ceiling tie that prefers the locked gate and Tacstd2 count>0. It is not a larger correlation than the other ceiling rows.

CD8 fraction does not reach the ceiling. Its maximum is n=4, ρ=-0.800, exact p=0.3333 (8/24) (study `drop_GSE154977`, genotype `all`, treatment `all`, epithelial gate `locked`, min epithelial cells `10`, Tacstd2% `pct_ge5`, endpoint `frac_cd8`).

Next rows under the |ρ| sort:

| |ρ| | ρ | exact p | n | studies | genotype filter | treatment | gate | min epi | Tacstd2% | endpoint | dataset residual |
|---:|---:|---:|---:|---:|---|---|---|---:|---|---|---:|
| 1.000 | -1.000 | 0.0833 | 4 | drop_GSE154977 | all | all | locked | 10 | pct_ge5 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | drop_GSE154977 | all | all | locked | 10 | pct_log1 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | drop_GSE154977 | all | no_Cis72 | locked | 10 | pct_ge5 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | drop_GSE154977 | all | no_Cis72 | locked | 10 | pct_log1 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | all3 | KL_K | all | locked | 10 | pct_ge5 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | all3 | KL_K | all | locked | 10 | pct_log1 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | drop_GSE154977 | KL_K | all | locked | 10 | pct_ge5 | frac_tnk | -1.000 |
| 1.000 | -1.000 | 0.0833 | 4 | drop_GSE154977 | KL_K | all | locked | 10 | pct_log1 | frac_tnk | -1.000 |


On the locked gate the four mixed mice are ordered the same way inside each study. Within-study Tacstd2% order versus T/NK is reported from the winner mice table. Each study alone has two mice, so a within-study Spearman is not an eligible test. The dataset residual stays −1 because those within-study orders agree. GSE154977 stays out of this endpoint: it is an AT2-lineage sort with essentially no T/NK cells.

## IFN/APM

IFN is the epithelial mean of the integrate ISG list, using only genes present in every library (Stat1, Stat2, Irf1, Irf7, Irf9, Isg15, Ifit1, Ifit2, Ifit3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58). Ifnb1 is absent from GSE180963, so it is not in the mean. APM is the epithelial mean of B2m, H2-K1, H2-D1, Tap1, Tap2, Tapbp, Psmb8, Psmb9, Psmb10, Nlrc5 (present in every library: B2m, H2-K1, H2-D1, Tap1, Tap2, Tapbp, Psmb8, Psmb9, Psmb10, Nlrc5). IFN_APM is the mean of those two scores. MHC-II genes stay in the side endpoint `MHC_published` and do not enter this argmax.

Maximum |ρ| in IFN_APM: **n=6, ρ=+1.000, exact p=0.0028 (2/720)**.

Shown specification (study `drop_GSE180963`, genotype `all`, treatment `all`, epithelial gate `locked`, min epithelial cells `10`, Tacstd2% `pct_log1`, endpoint `IFN`). Dataset-residual Spearman ρ=+1.000 (asymptotic p=t approximation saturates at |ρ|=1). Genotype-residual Spearman ρ=+1.000 (asymptotic p=t approximation saturates at |ρ|=1).

Mice: GSE154977_KP_30w_Cis72_m5,GSE154977_KP_30w_Cis72_m6,GSE154977_KP_30w_ND_m3,GSE154977_KP_30w_ND_m4,GSE165641_KL1,GSE165641_KL2.

Eligible tests in this family: 10740. Tests on the |ρ|=1 ceiling: 1953 (0 negative, 1953 positive). The row above is the ceiling tie that prefers the locked gate and Tacstd2 count>0. It is not a larger correlation than the other ceiling rows.

Largest mouse set in this family (8 mice): n=8, ρ=+0.833, exact p=0.0154 (620/40320). study `all3`, genotype `all`, treatment `all`, epithelial gate `nkx21_ptprc_neg`, min epithelial cells `10`, Tacstd2% `pct_log1`, endpoint `IFN`. Dataset-residual ρ=+0.857.

Next rows under the |ρ| sort:

| |ρ| | ρ | exact p | n | studies | genotype filter | treatment | gate | min epi | Tacstd2% | endpoint | dataset residual |
|---:|---:|---:|---:|---:|---|---|---|---:|---|---|---:|
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | all | all | locked | 10 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | KP_KL | all | locked | 10 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | all | all | locked | 30 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | all | all | locked | 50 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | KP_KL | all | locked | 30 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | KP_KL | all | locked | 50 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | all | all | epcam_ptprc_neg | 10 | pct_log1 | IFN | +1.000 |
| 1.000 | +1.000 | 0.0028 | 6 | drop_GSE180963 | all | all | epcam_sftpc_ptprc_neg | 10 | pct_log1 | IFN | +1.000 |


## MHC side endpoint (not the headline)

Maximum |ρ| in MHC_published: **n=4, ρ=+1.000, exact p=0.0833 (2/24)**.

Shown specification (study `only_GSE154977`, genotype `all`, treatment `all`, epithelial gate `locked`, min epithelial cells `10`, Tacstd2% `pct_ge5`, endpoint `MHC_published`). Dataset-residual Spearman ρ=+1.000 (asymptotic p=t approximation saturates at |ρ|=1). Genotype-residual Spearman ρ=+1.000 (asymptotic p=t approximation saturates at |ρ|=1).

Mice: GSE154977_KP_30w_Cis72_m5,GSE154977_KP_30w_Cis72_m6,GSE154977_KP_30w_ND_m3,GSE154977_KP_30w_ND_m4.

Eligible tests in this family: 3580. Tests on the |ρ|=1 ceiling: 350 (0 negative, 350 positive). The row above is the ceiling tie that prefers the locked gate and Tacstd2 count>0. It is not a larger correlation than the other ceiling rows.

Largest mouse set in this family (8 mice): n=8, ρ=-0.571, exact p=0.1511 (6094/40320). study `all3`, genotype `all`, treatment `all`, epithelial gate `broad_lung_ptprc_neg`, min epithelial cells `10`, Tacstd2% `pct_gt0`, endpoint `MHC_published`. Dataset-residual ρ=-0.238.

Smallest exact p in this family (not the |ρ| rule): n=5, ρ=+1.000, exact p=0.0167 (2/120). study `drop_GSE165641`, genotype `KP_KL`, treatment `all`, epithelial gate `locked`, min epithelial cells `10`, Tacstd2% `pct_ge5`, endpoint `MHC_published`. Dataset-residual ρ=+0.900.

Next rows under the |ρ| sort:

| |ρ| | ρ | exact p | n | studies | genotype filter | treatment | gate | min epi | Tacstd2% | endpoint | dataset residual |
|---:|---:|---:|---:|---:|---|---|---|---:|---|---|---:|
| 1.000 | +1.000 | 0.0833 | 4 | only_GSE154977 | all | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | all3 | KP | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | drop_GSE165641 | KP | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | drop_GSE180963 | KP | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | drop_GSE180963 | KP_K | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | only_GSE154977 | KP | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | only_GSE154977 | KP_KL | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |
| 1.000 | +1.000 | 0.0833 | 4 | only_GSE154977 | KP_K | all | locked | 10 | pct_ge5 | MHC_published | +1.000 |


## What the grid is allowed to change

Study filter, genotype filter, cisplatin drop, epithelial gate, epithelial-cell floor (10, 30, 50), and the Tacstd2 percent definition (count > 0, count ≥ 2, count ≥ 5, log1p CPM-10k ≥ 1). T endpoints are the T/NK call, CD3, and CD8 fractions of QC cells, with T/NK-marker cells that also pass the epithelial gate removed from the numerator.

Cells after the shared QC (nFeature ≥ 200, nCount ≥ 500, percent.mt < 25): 31970. Mice: 8. Private 8 KL mice used: 0.

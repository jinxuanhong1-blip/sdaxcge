# Kras-mutant mouse lung ICB: Cldn4 vs response and vs NHEJ/IFN

Public matrices only. TISMO Tacstd2 49/64 is not recomputed. No series here is a KL (Kras/Lkb1) ICB experiment. LLC and LL/2 are KrasG12C/Nras Lewis lung, not KL. CMT-167 (KrasG12V) has no ICB rows in TISMO. MLE12 is not Kras.

## Strongest thesis-aligned statistics

A headline row has the pre-set sign, no sample dropped, and a non-trivial effect (|ρ|≥0.30, percentile-score gap ≥0.02, other score gaps ≥0.20, resistant-minus-reference median gap ≥0.25, odds ratio ≥1.5). Response headlines also exclude the Setdb1-confounded TISMO label and need n≥3 per arm. Module headlines need n≥8. They may use any panel in that question, any score, Spearman, Pearson, or a median/tertile/quartile split. BH q is computed inside each question across every row of the sweep, aligned or not.

- **Cldn4 vs response.** No thesis-aligned statistic met the reporting rule.
- **Cldn4 vs IFN.** GSE262305_LLC (LLC), pearson / rank_within_sample / ifn_chemokine pct, filter all: effect -0.810 (rho), p 0.0082, BH q 0.4797, n 9. Cldn4 span on log2_value is 0.000 to 0.306 (width 0.306). That width is a near-floor gene, so the correlation is carried by a very small expression range. The module score itself spans 0.102.
- **Cldn4 vs NHEJ.** GSE246922_KP (KP), spearman / rank_within_sample / nhej_union pct, filter all: effect 0.607 (rho), p 0.0164, BH q 0.5394, n 15. Cldn4 span on author_scale is 8.636 to 9.957 (width 1.321). The module score itself spans 0.009.

Sweep size: cldn4_ifn 1372 tests (686 thesis-aligned), cldn4_nhej 1029 tests (395 thesis-aligned), cldn4_response 66 tests (3 thesis-aligned).

## Rows with BH q < 0.05

These are multiplicity-adjusted hits on the full sample set. The sign column says whether the hit matches the pre-set thesis direction.

| Question | Cohort | Test | Effect | p | q | n | Thesis sign |
|---|---|---|---:|---:|---:|---:|---|
| cldn4_nhej | GSE157880_HKP1 | nhej_extended pct pearson | -0.842 | 4.38e-05 | 0.0451 | 16 | opposite |
| cldn4_nhej | GSE157880_HKP1 | nhej_extended pct spearman | -0.824 | 8.84e-05 | 0.0455 | 16 | opposite |

The HKP1 BH rows are the extended NHEJ neighborhood (53BP1, ATM, PARP1, and the other genes outside the c-NHEJ core), scored as a within-sample percentile. c-NHEJ core on the same 16 samples is Spearman ρ=-0.288, p=0.2790. The extended-set correlation is in the IgG lungs (ρ=-0.881, n=8, p=0.0039) and in the anti-PD-1 lungs (ρ=-0.833, n=8, p=0.0102). 0 Gy alone is n=5, ρ=-0.600, p=0.2848. Radiated lungs are n=11, ρ=-0.827, p=0.0017.

## Pre-specified module correlations

Spearman of Cldn4 with the z-mean of c-NHEJ core or the IFN/chemokine panel, all samples, one expression scale per cohort. Thesis sign is positive for NHEJ and negative for IFN.

| Cohort | Module | n | ρ | p | q | Sign |
|---|---|---:|---:|---:|---:|---|
| GSE155972_LLC | nhej_core | 33 | -0.267 | 0.1328 | 0.5815 | opposite |
| GSE155972_LLC | ifn_chemokine | 33 | -0.083 | 0.6446 | 0.9766 | matches |
| GSE246922_KP | nhej_core | 15 | -0.029 | 0.9195 | 1.0000 | opposite |
| GSE246922_KP | ifn_chemokine | 15 | 0.046 | 0.8695 | 1.0000 | opposite |
| GSE246922_LLC1 | nhej_core | 9 | -0.226 | 0.5582 | 0.9004 | opposite |
| GSE246922_LLC1 | ifn_chemokine | 9 | 0.261 | 0.4974 | 0.9086 | opposite |
| GSE297630_LLC | nhej_core | 6 | 0.257 | 0.6228 | 0.9308 | matches |
| GSE297630_LLC | ifn_chemokine | 6 | 0.086 | 0.8717 | 1.0000 | opposite |
| GSE114601_KP | nhej_core | 8 | 0.190 | 0.6514 | 0.9441 | matches |
| GSE114601_KP | ifn_chemokine | 8 | 0.333 | 0.4198 | 0.8647 | opposite |
| GSE157880_HKP1 | nhej_core | 16 | -0.303 | 0.2541 | 0.6972 | opposite |
| GSE157880_HKP1 | ifn_chemokine | 16 | -0.176 | 0.5133 | 0.9161 | matches |
| GSE274960_LL2 | nhej_core | 13 | -0.514 | 0.0721 | 0.5394 | opposite |
| GSE274960_LL2 | ifn_chemokine | 13 | 0.602 | 0.0293 | 0.6701 | opposite |
| GSE262305_LLC | nhej_core | 9 | -0.310 | 0.4175 | 0.8094 | opposite |
| GSE262305_LLC | ifn_chemokine | 9 | -0.477 | 0.1942 | 0.7362 | matches |
| GSE239485_LLC | nhej_core | 24 | -0.002 | 0.9936 | 1.0000 | opposite |
| GSE239485_LLC | ifn_chemokine | 24 | 0.056 | 0.7962 | 1.0000 | opposite |
| GSE260596_KP | nhej_core | 7 | -0.214 | 0.6445 | 0.9367 | opposite |
| GSE260596_KP | ifn_chemokine | 7 | -0.571 | 0.1802 | 0.7362 | matches |

## Pre-specified response contrasts

These are Mann-Whitney tests on the pooled resistant-versus-reference contrast, all samples, one expression scale. They are reported whether or not they win the sweep.

| Cohort | Model | What the arms are | n resistant vs reference | median difference | p | BH q | sign |
|---|---|---|---:|---:|---:|---:|---|
| GSE155972_LLC | LLC | Setdb1 genotype, not an independent response call | 6 vs 7 | -0.111 | 0.4452 | 1.0000 | resistant lower or flat |
| GSE246922_KP | KP | acquired resistance | 9 vs 3 | -0.283 | 0.7273 | 1.0000 | resistant lower or flat |
| GSE246922_LLC1 | LLC1 | acquired resistance | 3 vs 3 | -0.195 | 0.1840 | 1.0000 | resistant lower or flat |
| GSE297630_LLC | LLC | tolerant proxy | 3 vs 3 | -0.180 | 0.1000 | 1.0000 | resistant lower or flat |

Pre-specified Mann-Whitney sign: Cldn4 is higher in the resistant arm in 0/4 and lower in 4/4.

GSE155972 is in that table and is not eligible for the response headline. TISMO calls Setdb1-KO ICB mice Responders and control-sgRNA ICB mice Non-responders.

## Cohorts

| Cohort | Model | Kras | Response label | Setting |
|---|---|---|---|---|
| GSE155972_LLC | LLC | KrasG12C and NrasQ61; not Stk11/Lkb1 | tismo_label; confounded | in vivo subcutaneous LLC, anti-PD-1 + anti-CTLA-4 |
| GSE246922_KP | KP | Kras/Trp53 lung cancer cells; not Lkb1 | acquired_resistance | CD45-negative cells from parental KP and from tumors that relapsed after immunotherapy |
| GSE246922_LLC1 | LLC1 | KrasG12C and Nras; not Stk11/Lkb1 | acquired_resistance | CD45-negative LLC1 cells, parental versus relapsed after immunotherapy |
| GSE297630_LLC | LLC | KrasG12C and Nras; not Stk11/Lkb1 | tolerant_proxy | subcutaneous LLC tumors; series title calls the anti-PD-1 arm tolerant cells |
| GSE114601_KP | KP | Kras/Trp53 GEMM lung tumor; not Lkb1 | none (modules and treatment context only) | KP lung nodules; anti-PD-1, JQ1, combination, vehicle. No per-mouse R/NR label. |
| GSE157880_HKP1 | HKP1 | KrasG12D/Trp53-null lung line; not Lkb1 | none (modules and treatment context only) | HKP1-bearing lungs, anti-PD-1 or IgG, with or without radiation. No per-mouse R/NR label. |
| GSE274960_LL2 | LL/2 | LL/2 is Lewis lung, KrasG12C; not Stk11/Lkb1 | none (modules and treatment context only) | LL/2 tumors, IgG, anti-PD-1, entrectinib, or the combination. shRNA arms are excluded. |
| GSE262305_LLC | LLC | KrasG12C and Nras; not Stk11/Lkb1 | none (modules and treatment context only) | subcutaneous LLC, isotype, anti-PD-L1, or bortezomib plus anti-PD-L1. No R/NR label. |
| GSE239485_LLC | LLC | KrasG12C and Nras; not Stk11/Lkb1 | none (modules and treatment context only) | LLC tumors. Treated arms are Poly I:C plus anti-PD-1, with or without anti-C5aR1. No monotherapy and no R/NR label. |
| GSE260596_KP | KP | refractory KP mouse lung tumors (series summary); not Lkb1 | none (modules and treatment context only) | All profiled tumors are on anti-PD-1, with control or anti-LAIR1. No untreated arm and no per-mouse R/NR label. |

## Screened and not scored

- TISMO CMT-167: KrasG12V lung carcinoma, no ICB expression rows.
- TISMO MLE12: lung adenocarcinoma driven by SV40 large T, not Kras, and no ICB rows.
- TISMO KPB25L: Kras/p53 mammary, not lung.
- GSE193895: Kras/Keap1/Lkb1 GEMM tumors, deposited matrix has no ICB arm.
- GSE277929: KL tumors treated with entinostat and trametinib, not ICB.
- GSE137244: locked KL-versus-KP cell-line result, not an ICB experiment, not rerun.

## Scales

Count matrices are tested as log2(count+1) and log2(CPM+1). FPKM and normalized counts use log2(value+1). Author VST, RMA, and the GSE239485 processed matrix are used as deposited. Within-sample percentile scores do not depend on a monotone per-sample rescaling, so they are emitted once per cohort.

Chronic IFNG arms in GSE246922 stay out of the resistant group. They are in the module correlations.


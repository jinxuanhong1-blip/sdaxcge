# RESULTS — pySCENIC on GSE131907 malignant cells, CLDN4-high vs low

Cohort: **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277), author-malignant cells.
Method that was actually run: **pySCENIC 0.12.1** (GRNBoost2 → cisTarget motif pruning → AUCell).
SCENIC+ was not run: this GEO deposit is scRNA-seq only (no matched scATAC).
Concordant-4 was not co-embedded. This machine has 15 GiB RAM; cisTarget rankings plus a four-cohort malignant matrix do not fit. GSE131907 is the largest single member of that concordant set.

This file separates three things that are easy to mix:

1. **Observational** regulon / program AUCell on tumor cells (this run).
2. **KD-like expectation** used by earlier AUCell work in this project: IFN and MHC activity **lower** in CLDN4-high cells (Δ high−low < 0), TJ **higher** in CLDN4-high. That is the direction in which CLDN4 loss would raise IFN. It is **not** estimated from these cells.
3. **Public CLDN4-loss transcriptomes**, which do not agree with each other. **GSE207704** (CLDN4 CRISPR in T47D/MCF7) has IFN trending **down** after loss, which predicts observational Δ > 0. **GSE50927** (whole-lung Cldn4 KO, n=1) has IFN/MHC **up** after loss (NES +1.51 in the earlier public summary), which predicts observational Δ < 0. GSE50927 is not a cancer-cell replicate.

A10 ELF3–CLDN4 is taken as given. ELF3 rows are labeled `elf3_given` and are not a new discovery.

## Observational IFN sign (read this before the rest)

The primary IFN board is STAT/IRF cisTarget regulons plus Hallmark IFN program AUCell (IFN-α, IFN-γ, and their union). Regulons that only overlap Hallmark IFN genes are **not** on this board and are **not** counted in the SAME / OPPOSITE totals below. They are printed in full in the next subsection so those signs are not dropped.

**PROGRAM::IFN_UNION** is higher in CLDN4-high (median Δ=0.003639, 26/29 patients Δ>0, p=4.70e-06; tLung Δ=0.003607, higher in CLDN4-high, p=0.02734); **OPPOSITE** versus the KD-like expectation, **SAME** versus GSE207704, **OPPOSITE** versus GSE50927. **PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE** is higher in CLDN4-high (median Δ=0.003594, 24/29 patients Δ>0, p=1.60e-05; tLung Δ=0.003379, higher in CLDN4-high, p=0.04883); **OPPOSITE** versus the KD-like expectation, **SAME** versus GSE207704, **OPPOSITE** versus GSE50927. **PROGRAM::HALLMARK_INTERFERON_ALPHA_RESPONSE** is higher in CLDN4-high (median Δ=0.004886, 23/29 patients Δ>0, p=1.56e-04; tLung Δ=0.004924, higher in CLDN4-high, p=0.08398); **OPPOSITE** versus the KD-like expectation, **SAME** versus GSE207704, **OPPOSITE** versus GSE50927. KD-like expectation for these rows: IFN **lower** in CLDN4-high (Δ < 0), the frame in which CLDN4 loss would raise IFN. That expectation was **not** fit on this matrix. GSE207704 (public CLDN4 CRISPR in T47D/MCF7) has IFN trending down after loss and predicts Δ > 0. GSE50927 (public whole-lung Cldn4 KO, n=1, not a cancer-cell test) has IFN/MHC up after loss and predicts Δ < 0. A row can be SAME versus one reference and OPPOSITE versus the other.

cisTarget STAT/IRF regulons at NES ≥ 3.0, CLDN4 held out (6 regulons). Observational sign **OPPOSITE** the KD-like expectation (higher in CLDN4-high): IRF1, STAT1, STAT2, IRF7, IRF3. Observational sign **SAME** as the KD-like expectation (higher in CLDN4-low): IRF2. Both lists are on the board. A non-significant p does not remove the sign.

**PROGRAM::MHC1_APM** is higher in CLDN4-high (median Δ=0.01588, 24/29 patients Δ>0, p=1.60e-05; tLung Δ=0.004258, higher in CLDN4-high, p=0.2754); **OPPOSITE** versus the KD-like expectation, **SAME** versus GSE207704, **OPPOSITE** versus GSE50927. KD-like expectation: MHC lower in CLDN4-high. cisTarget MHC TFs: RFXANK higher in CLDN4-low (Δ=-0.005909, p=7.71e-07, **SAME** KD-like); RFX5 higher in CLDN4-high (Δ=0.0011, p=0.09633, **OPPOSITE** KD-like). If the MHC program and a cisTarget MHC regulon disagree, both stay in this paragraph and in the MHC table.

**PROGRAM::TJ_STRUCT** is higher in CLDN4-high (median Δ=0.02886, 29/29 patients Δ>0, p=3.73e-09; tLung Δ=0.02567, higher in CLDN4-high, p=0.001953); **SAME** versus the KD-like expectation, **NA** versus GSE207704, **NA** versus GSE50927. **PROGRAM::KERATIN** is higher in CLDN4-high (median Δ=0.02998, 29/29 patients Δ>0, p=3.73e-09; tLung Δ=0.03048, higher in CLDN4-high, p=0.001953); **SAME** versus the KD-like expectation, **NA** versus GSE207704, **NA** versus GSE50927. KD-like expectation: TJ and keratin higher in CLDN4-high. cisTarget TJ-associated TFs: ELF3 (A10-given, not a discovery) higher in CLDN4-high (Δ=0.06082, p=3.73e-09, **SAME** KD-like); KLF4 higher in CLDN4-high (Δ=0.02801, p=3.73e-09, **SAME** KD-like); KLF5 higher in CLDN4-high (Δ=0.00766, p=1.56e-04, **SAME** KD-like); TFAP2C higher in CLDN4-low (Δ=-0.001727, p=0.02155, **OPPOSITE** KD-like); TFAP2A higher in CLDN4-low (Δ=-6.70e-04, p=0.3579, **OPPOSITE** KD-like; all-site direction and tLung direction disagree (tLung higher in CLDN4-high, p=0.6953)); OVOL2 higher in CLDN4-low (Δ=-1.24e-04, p=0.4946, **OPPOSITE** KD-like; all-site direction and tLung direction disagree (tLung higher in CLDN4-high, p=0.02734)). **PROGRAM::CTRL_RIBO** is higher in CLDN4-high (median Δ=0.00137, 17/29 patients Δ>0, p=0.7172; tLung Δ=-0.001861, higher in CLDN4-low, p=0.8457); **NA** versus the KD-like expectation, **NA** versus GSE207704, **NA** versus GSE50927 (null control; SAME/OPPOSITE is not applied).

Between-patient Spearman of malignant CLDN4+ fraction versus mean AUCell (a different question from the within-sample split): IFN_UNION ρ=0.08589, p=0.6459, n=31 patients; STAT1 ρ=0.02097, p=0.9109, n=31 patients; MHC1_APM ρ=0.1766, p=0.3419, n=31 patients; RFXANK ρ=-0.3964, p=0.02728, n=31 patients; TJ_STRUCT ρ=0.6323, p=1.36e-04, n=31 patients. A null between-patient IFN correlation does not replace the within-sample patient Δ, and the within-sample Δ does not replace the between-patient correlation.

Removing CLDN4 from IFN, MHC, and TJ cisTarget regulons did not change any observational direction relative to the as-returned regulon.

Primary IFN board, observational sign **OPPOSITE** the KD-like expectation: **8**. Primary IFN board, **SAME** sign as that expectation: **1**. These two counts are only the STAT/IRF and Hallmark/union rows below. They do not include Hallmark-overlap regulons.

- **IRF1(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.01605 (higher in CLDN4-high); 28/29 patients Δ>0; Wilcoxon p=1.12e-08; tLung median Δ = 0.02205 (higher in CLDN4-high, p=0.003906). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **STAT1(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.00481 (higher in CLDN4-high); 25/29 patients Δ>0; Wilcoxon p=7.58e-06; tLung median Δ = 0.003928 (higher in CLDN4-high, p=0.04883). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **STAT2(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.001502 (higher in CLDN4-high); 21/29 patients Δ>0; Wilcoxon p=0.00128; tLung median Δ = 7.91e-04 (higher in CLDN4-high, p=0.2324). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **IRF7(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.002458 (higher in CLDN4-high); 20/29 patients Δ>0; Wilcoxon p=0.003242; tLung median Δ = 0.003115 (higher in CLDN4-high, p=0.2324). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **IRF2(+)|CLDN4out** (cistarget_regulon): observational median Δ = -0.001529 (higher in CLDN4-low); 10/29 patients Δ>0; Wilcoxon p=0.02155; tLung median Δ = -0.001431 (higher in CLDN4-low, p=0.1055). Versus KD-like expectation (IFN lower in CLDN4-high): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **IRF3(+)|CLDN4out** (cistarget_regulon): observational median Δ = 2.00e-04 (higher in CLDN4-high); 15/29 patients Δ>0; Wilcoxon p=0.4946; tLung median Δ = 8.89e-05 (higher in CLDN4-high, p=0.7695). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **PROGRAM::IFN_UNION** (program_AUCell_not_cistarget): observational median Δ = 0.003639 (higher in CLDN4-high); 26/29 patients Δ>0; Wilcoxon p=4.70e-06; tLung median Δ = 0.003607 (higher in CLDN4-high, p=0.02734). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE** (program_AUCell_not_cistarget): observational median Δ = 0.003594 (higher in CLDN4-high); 24/29 patients Δ>0; Wilcoxon p=1.60e-05; tLung median Δ = 0.003379 (higher in CLDN4-high, p=0.04883). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **PROGRAM::HALLMARK_INTERFERON_ALPHA_RESPONSE** (program_AUCell_not_cistarget): observational median Δ = 0.004886 (higher in CLDN4-high); 23/29 patients Δ>0; Wilcoxon p=1.56e-04; tLung median Δ = 0.004924 (higher in CLDN4-high, p=0.08398). Versus KD-like expectation (IFN lower in CLDN4-high): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.

The KD-like column is the pre-specified “IFN lower in CLDN4-high” direction. The GSE207704 column is often the opposite reference, because that public CRISPR series did not open IFN. A row can be SAME versus one reference and OPPOSITE versus the other. That disagreement is the result, not a reason to drop a row.

Program AUCell rows are **not** cisTarget regulons. They use the same recovery curve (top 5% of the detected-gene ranking) on Hallmark IFN sets with CLDN4 removed. cisTarget rows on this board are STAT1, STAT2, or IRF1–IRF9 regulons that passed motif enrichment at NES ≥ 3.0, rescored after removing CLDN4 so the split gene cannot sit inside the signature.

### Primary IFN table

| score | tf | axis | kind | n_patients | delta_median_high_minus_low | direction | n_patients_delta_pos | n_patients_delta_neg | p_patient | fdr_patient | vs_kd_like | vs_GSE207704 | vs_GSE50927 | delta_median_tLung | direction_tLung | p_tLung |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| IRF1(+)/CLDN4out | IRF1 | IFN | cistarget_regulon | 29 | 0.01605 | higher in CLDN4-high | 28 | 1 | 1.12e-08 | 1.69e-07 | OPPOSITE | SAME | OPPOSITE | 0.02205 | higher in CLDN4-high | 0.003906 |
| STAT1(+)/CLDN4out | STAT1 | IFN | cistarget_regulon | 29 | 0.00481 | higher in CLDN4-high | 25 | 4 | 7.58e-06 | 4.42e-05 | OPPOSITE | SAME | OPPOSITE | 0.003928 | higher in CLDN4-high | 0.04883 |
| STAT2(+)/CLDN4out | STAT2 | IFN | cistarget_regulon | 29 | 0.001502 | higher in CLDN4-high | 21 | 8 | 0.00128 | 0.003861 | OPPOSITE | SAME | OPPOSITE | 7.91e-04 | higher in CLDN4-high | 0.2324 |
| IRF7(+)/CLDN4out | IRF7 | IFN | cistarget_regulon | 29 | 0.002458 | higher in CLDN4-high | 20 | 9 | 0.003242 | 0.00876 | OPPOSITE | SAME | OPPOSITE | 0.003115 | higher in CLDN4-high | 0.2324 |
| IRF2(+)/CLDN4out | IRF2 | IFN | cistarget_regulon | 29 | -0.001529 | higher in CLDN4-low | 10 | 19 | 0.02155 | 0.04383 | SAME | OPPOSITE | SAME | -0.001431 | higher in CLDN4-low | 0.1055 |
| IRF3(+)/CLDN4out | IRF3 | IFN | cistarget_regulon | 29 | 2.00e-04 | higher in CLDN4-high | 15 | 14 | 0.4946 | 0.5561 | OPPOSITE | SAME | OPPOSITE | 8.89e-05 | higher in CLDN4-high | 0.7695 |
| PROGRAM::IFN_UNION |  | IFN | program_AUCell_not_cistarget | 29 | 0.003639 | higher in CLDN4-high | 26 | 3 | 4.70e-06 | 1.10e-05 | OPPOSITE | SAME | OPPOSITE | 0.003607 | higher in CLDN4-high | 0.02734 |
| PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE |  | IFN | program_AUCell_not_cistarget | 29 | 0.003594 | higher in CLDN4-high | 24 | 5 | 1.60e-05 | 2.23e-05 | OPPOSITE | SAME | OPPOSITE | 0.003379 | higher in CLDN4-high | 0.04883 |
| PROGRAM::HALLMARK_INTERFERON_ALPHA_RESPONSE |  | IFN | program_AUCell_not_cistarget | 29 | 0.004886 | higher in CLDN4-high | 23 | 6 | 1.56e-04 | 1.82e-04 | OPPOSITE | SAME | OPPOSITE | 0.004924 | higher in CLDN4-high | 0.08398 |

### Hallmark-IFN target overlap (not IFN transcription factors)

These cisTarget regulons are **not** STAT or IRF. Each has ≥5 targets in the Hallmark IFN union and a hypergeometric FDR < 0.05. The SAME / OPPOSITE labels use the IFN KD-like direction so an opposite observational sign is still printed. That column does not rename FOS, JUNB, ATF3, or any other TF here as an interferon regulator.

Overlap rows if the IFN KD-like direction is used as the bookkeeping reference: **17** OPPOSITE, **7** SAME. These counts are **not** added to the primary IFN totals above.

- **FOS(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.03678 (higher in CLDN4-high); 29/29 patients Δ>0; Wilcoxon p=3.73e-09; tLung median Δ = 0.03574 (higher in CLDN4-high, p=0.001953). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **JUNB(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.04625 (higher in CLDN4-high); 29/29 patients Δ>0; Wilcoxon p=3.73e-09; tLung median Δ = 0.04566 (higher in CLDN4-high, p=0.001953). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **ATF3(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.02159 (higher in CLDN4-high); 28/29 patients Δ>0; Wilcoxon p=7.45e-09; tLung median Δ = 0.02179 (higher in CLDN4-high, p=0.001953). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **EGR1(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.01741 (higher in CLDN4-high); 27/29 patients Δ>0; Wilcoxon p=2.61e-08; tLung median Δ = 0.01706 (higher in CLDN4-high, p=0.003906). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **CEBPD(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.01069 (higher in CLDN4-high); 27/29 patients Δ>0; Wilcoxon p=7.08e-08; tLung median Δ = 0.01106 (higher in CLDN4-high, p=0.001953). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **KLF6(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.007133 (higher in CLDN4-high); 27/29 patients Δ>0; Wilcoxon p=7.08e-08; tLung median Δ = 0.007835 (higher in CLDN4-high, p=0.001953). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **ELK3(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.004996 (higher in CLDN4-high); 26/29 patients Δ>0; Wilcoxon p=3.28e-07; tLung median Δ = 0.01033 (higher in CLDN4-high, p=0.003906). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **EHF(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.01496 (higher in CLDN4-high); 25/29 patients Δ>0; Wilcoxon p=3.37e-06; tLung median Δ = 0.007928 (higher in CLDN4-high, p=0.005859). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **BHLHE40(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.004759 (higher in CLDN4-high); 25/29 patients Δ>0; Wilcoxon p=3.98e-06; tLung median Δ = 0.005026 (higher in CLDN4-high, p=0.009766). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **ELK1(+)|CLDN4out** (cistarget_regulon): observational median Δ = -0.001595 (higher in CLDN4-low); 3/29 patients Δ>0; Wilcoxon p=8.84e-06; tLung median Δ = -0.001308 (higher in CLDN4-low, p=0.005859). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **FOSL2(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.003977 (higher in CLDN4-high); 25/29 patients Δ>0; Wilcoxon p=1.19e-05; tLung median Δ = 0.00358 (higher in CLDN4-high, p=0.02734). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **E2F1(+)|CLDN4out** (cistarget_regulon): observational median Δ = -0.005407 (higher in CLDN4-low); 5/29 patients Δ>0; Wilcoxon p=2.11e-05; tLung median Δ = -0.006099 (higher in CLDN4-low, p=0.003906). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **CEBPB(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.0048 (higher in CLDN4-high); 24/29 patients Δ>0; Wilcoxon p=1.10e-04; tLung median Δ = 0.004063 (higher in CLDN4-high, p=0.009766). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **NFE2L3(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.001147 (higher in CLDN4-high); 19/29 patients Δ>0; Wilcoxon p=0.004436; tLung median Δ = 0.001339 (higher in CLDN4-high, p=0.04883). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **NFIC(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.002185 (higher in CLDN4-high); 21/29 patients Δ>0; Wilcoxon p=0.008008; tLung median Δ = 0.005695 (higher in CLDN4-high, p=0.02734). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **ERF(+)|CLDN4out** (cistarget_regulon): observational median Δ = -3.60e-04 (higher in CLDN4-low); 9/29 patients Δ>0; Wilcoxon p=0.02906; tLung median Δ = -3.88e-04 (higher in CLDN4-low, p=0.2754). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **SMAD1(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.001605 (higher in CLDN4-high); 19/29 patients Δ>0; Wilcoxon p=0.05623; tLung median Δ = 0.001265 (higher in CLDN4-high, p=0.7695). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **SPDEF(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.004261 (higher in CLDN4-high); 19/29 patients Δ>0; Wilcoxon p=0.06552; tLung median Δ = -2.40e-04 (higher in CLDN4-low, p=0.9219). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **NR3C1(+)|CLDN4out** (cistarget_regulon): observational median Δ = 0.001021 (higher in CLDN4-high); 19/29 patients Δ>0; Wilcoxon p=0.07236; tLung median Δ = 0.001502 (higher in CLDN4-high, p=0.009766). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **IKZF2(+)|CLDN4out** (cistarget_regulon): observational median Δ = -9.39e-05 (higher in CLDN4-low); 13/29 patients Δ>0; Wilcoxon p=0.5647; tLung median Δ = 2.00e-04 (higher in CLDN4-high, p=0.4922). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **ETV7(+)|CLDN4out** (cistarget_regulon): observational median Δ = -2.23e-04 (higher in CLDN4-low); 12/29 patients Δ>0; Wilcoxon p=0.6545; tLung median Δ = -0.001859 (higher in CLDN4-low, p=0.2754). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **ELF1(+)|CLDN4out** (cistarget_regulon): observational median Δ = 3.57e-05 (higher in CLDN4-high); 15/29 patients Δ>0; Wilcoxon p=0.8314; tLung median Δ = -9.46e-05 (higher in CLDN4-low, p=0.8457). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **OPPOSITE**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **SAME**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **OPPOSITE**.
- **CEBPG(+)|CLDN4out** (cistarget_regulon): observational median Δ = -8.19e-05 (higher in CLDN4-low); 13/29 patients Δ>0; Wilcoxon p=0.9152; tLung median Δ = -6.68e-05 (higher in CLDN4-low, p=0.5566). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.
- **PHF20(+)|CLDN4out** (cistarget_regulon): observational median Δ = -6.75e-04 (higher in CLDN4-low); 13/29 patients Δ>0; Wilcoxon p=0.9321; tLung median Δ = -6.80e-04 (higher in CLDN4-low, p=0.8457). Versus KD-like expectation (IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator): **SAME**. Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **OPPOSITE**. Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **SAME**.

| score | tf | axis | kind | n_patients | delta_median_high_minus_low | direction | n_patients_delta_pos | n_patients_delta_neg | p_patient | fdr_patient | vs_kd_like | vs_GSE207704 | vs_GSE50927 | delta_median_tLung | direction_tLung | p_tLung |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FOS(+)/CLDN4out | FOS | IFN_by_targets | cistarget_regulon | 29 | 0.03678 | higher in CLDN4-high | 29 | 0 | 3.73e-09 | 9.63e-08 | OPPOSITE | SAME | OPPOSITE | 0.03574 | higher in CLDN4-high | 0.001953 |
| JUNB(+)/CLDN4out | JUNB | IFN_by_targets | cistarget_regulon | 29 | 0.04625 | higher in CLDN4-high | 29 | 0 | 3.73e-09 | 9.63e-08 | OPPOSITE | SAME | OPPOSITE | 0.04566 | higher in CLDN4-high | 0.001953 |
| ATF3(+)/CLDN4out | ATF3 | IFN_by_targets | cistarget_regulon | 29 | 0.02159 | higher in CLDN4-high | 28 | 1 | 7.45e-09 | 1.69e-07 | OPPOSITE | SAME | OPPOSITE | 0.02179 | higher in CLDN4-high | 0.001953 |
| EGR1(+)/CLDN4out | EGR1 | IFN_by_targets | cistarget_regulon | 29 | 0.01741 | higher in CLDN4-high | 27 | 2 | 2.61e-08 | 3.37e-07 | OPPOSITE | SAME | OPPOSITE | 0.01706 | higher in CLDN4-high | 0.003906 |
| CEBPD(+)/CLDN4out | CEBPD | IFN_by_targets | cistarget_regulon | 29 | 0.01069 | higher in CLDN4-high | 27 | 2 | 7.08e-08 | 8.01e-07 | OPPOSITE | SAME | OPPOSITE | 0.01106 | higher in CLDN4-high | 0.001953 |
| KLF6(+)/CLDN4out | KLF6 | IFN_by_targets | cistarget_regulon | 29 | 0.007133 | higher in CLDN4-high | 27 | 2 | 7.08e-08 | 8.01e-07 | OPPOSITE | SAME | OPPOSITE | 0.007835 | higher in CLDN4-high | 0.001953 |
| ELK3(+)/CLDN4out | ELK3 | IFN_by_targets | cistarget_regulon | 29 | 0.004996 | higher in CLDN4-high | 26 | 3 | 3.28e-07 | 3.30e-06 | OPPOSITE | SAME | OPPOSITE | 0.01033 | higher in CLDN4-high | 0.003906 |
| EHF(+)/CLDN4out | EHF | IFN_by_targets | cistarget_regulon | 29 | 0.01496 | higher in CLDN4-high | 25 | 4 | 3.37e-06 | 2.10e-05 | OPPOSITE | SAME | OPPOSITE | 0.007928 | higher in CLDN4-high | 0.005859 |
| BHLHE40(+)/CLDN4out | BHLHE40 | IFN_by_targets | cistarget_regulon | 29 | 0.004759 | higher in CLDN4-high | 25 | 4 | 3.98e-06 | 2.40e-05 | OPPOSITE | SAME | OPPOSITE | 0.005026 | higher in CLDN4-high | 0.009766 |
| ELK1(+)/CLDN4out | ELK1 | IFN_by_targets | cistarget_regulon | 29 | -0.001595 | higher in CLDN4-low | 3 | 26 | 8.84e-06 | 5.00e-05 | SAME | OPPOSITE | SAME | -0.001308 | higher in CLDN4-low | 0.005859 |
| FOSL2(+)/CLDN4out | FOSL2 | IFN_by_targets | cistarget_regulon | 29 | 0.003977 | higher in CLDN4-high | 25 | 4 | 1.19e-05 | 6.54e-05 | OPPOSITE | SAME | OPPOSITE | 0.00358 | higher in CLDN4-high | 0.02734 |
| E2F1(+)/CLDN4out | E2F1 | IFN_by_targets | cistarget_regulon | 29 | -0.005407 | higher in CLDN4-low | 5 | 24 | 2.11e-05 | 1.09e-04 | SAME | OPPOSITE | SAME | -0.006099 | higher in CLDN4-low | 0.003906 |
| CEBPB(+)/CLDN4out | CEBPB | IFN_by_targets | cistarget_regulon | 29 | 0.0048 | higher in CLDN4-high | 24 | 5 | 1.10e-04 | 5.00e-04 | OPPOSITE | SAME | OPPOSITE | 0.004063 | higher in CLDN4-high | 0.009766 |
| NFE2L3(+)/CLDN4out | NFE2L3 | IFN_by_targets | cistarget_regulon | 29 | 0.001147 | higher in CLDN4-high | 19 | 10 | 0.004436 | 0.011 | OPPOSITE | SAME | OPPOSITE | 0.001339 | higher in CLDN4-high | 0.04883 |
| NFIC(+)/CLDN4out | NFIC | IFN_by_targets | cistarget_regulon | 29 | 0.002185 | higher in CLDN4-high | 21 | 8 | 0.008008 | 0.01907 | OPPOSITE | SAME | OPPOSITE | 0.005695 | higher in CLDN4-high | 0.02734 |
| ERF(+)/CLDN4out | ERF | IFN_by_targets | cistarget_regulon | 29 | -3.60e-04 | higher in CLDN4-low | 9 | 20 | 0.02906 | 0.05596 | SAME | OPPOSITE | SAME | -3.88e-04 | higher in CLDN4-low | 0.2754 |
| SMAD1(+)/CLDN4out | SMAD1 | IFN_by_targets | cistarget_regulon | 29 | 0.001605 | higher in CLDN4-high | 19 | 10 | 0.05623 | 0.09785 | OPPOSITE | SAME | OPPOSITE | 0.001265 | higher in CLDN4-high | 0.7695 |
| SPDEF(+)/CLDN4out | SPDEF | IFN_by_targets | cistarget_regulon | 29 | 0.004261 | higher in CLDN4-high | 19 | 10 | 0.06552 | 0.1098 | OPPOSITE | SAME | OPPOSITE | -2.40e-04 | higher in CLDN4-low | 0.9219 |
| NR3C1(+)/CLDN4out | NR3C1 | IFN_by_targets | cistarget_regulon | 29 | 0.001021 | higher in CLDN4-high | 19 | 10 | 0.07236 | 0.1202 | OPPOSITE | SAME | OPPOSITE | 0.001502 | higher in CLDN4-high | 0.009766 |
| IKZF2(+)/CLDN4out | IKZF2 | IFN_by_targets | cistarget_regulon | 29 | -9.39e-05 | higher in CLDN4-low | 13 | 16 | 0.5647 | 0.6233 | SAME | OPPOSITE | SAME | 2.00e-04 | higher in CLDN4-high | 0.4922 |
| ETV7(+)/CLDN4out | ETV7 | IFN_by_targets | cistarget_regulon | 29 | -2.23e-04 | higher in CLDN4-low | 12 | 17 | 0.6545 | 0.7051 | SAME | OPPOSITE | SAME | -0.001859 | higher in CLDN4-low | 0.2754 |
| ELF1(+)/CLDN4out | ELF1 | IFN_by_targets | cistarget_regulon | 29 | 3.57e-05 | higher in CLDN4-high | 15 | 14 | 0.8314 | 0.8648 | OPPOSITE | SAME | OPPOSITE | -9.46e-05 | higher in CLDN4-low | 0.8457 |
| CEBPG(+)/CLDN4out | CEBPG | IFN_by_targets | cistarget_regulon | 29 | -8.19e-05 | higher in CLDN4-low | 13 | 16 | 0.9152 | 0.9412 | SAME | OPPOSITE | SAME | -6.68e-05 | higher in CLDN4-low | 0.5566 |
| PHF20(+)/CLDN4out | PHF20 | IFN_by_targets | cistarget_regulon | 29 | -6.75e-04 | higher in CLDN4-low | 13 | 16 | 0.9321 | 0.9532 | SAME | OPPOSITE | SAME | -6.80e-04 | higher in CLDN4-low | 0.8457 |

### MHC, TJ, keratin, ribosome control

KD-like expectation: MHC lower in CLDN4-high; TJ and keratin higher in CLDN4-high; ribosome null. GSE207704 / GSE50927 columns are filled only for MHC (same public IFN/MHC notes as above). A cisTarget row that disagrees with its program is left in this table.

| score | tf | axis | kind | n_patients | delta_median_high_minus_low | direction | n_patients_delta_pos | n_patients_delta_neg | p_patient | fdr_patient | vs_kd_like | vs_GSE207704 | vs_GSE50927 | delta_median_tLung | direction_tLung | p_tLung |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PROGRAM::CTRL_RIBO |  | CTRL | program_AUCell_not_cistarget | 29 | 0.00137 | higher in CLDN4-high | 17 | 12 | 0.7172 | 0.7172 | NA | NA | NA | -0.001861 | higher in CLDN4-low | 0.8457 |
| PROGRAM::KERATIN |  | KERATIN | program_AUCell_not_cistarget | 29 | 0.02998 | higher in CLDN4-high | 29 | 0 | 3.73e-09 | 1.30e-08 | SAME | NA | NA | 0.03048 | higher in CLDN4-high | 0.001953 |
| RFXANK(+)/CLDN4out | RFXANK | MHC | cistarget_regulon | 29 | -0.005909 | higher in CLDN4-low | 2 | 27 | 7.71e-07 | 6.65e-06 | SAME | OPPOSITE | SAME | -0.006267 | higher in CLDN4-low | 0.01953 |
| RFX5(+)/CLDN4out | RFX5 | MHC | cistarget_regulon | 29 | 0.0011 | higher in CLDN4-high | 19 | 10 | 0.09633 | 0.1529 | OPPOSITE | SAME | OPPOSITE | 0.002121 | higher in CLDN4-high | 0.08398 |
| PROGRAM::MHC1_APM |  | MHC | program_AUCell_not_cistarget | 29 | 0.01588 | higher in CLDN4-high | 24 | 5 | 1.60e-05 | 2.23e-05 | OPPOSITE | SAME | OPPOSITE | 0.004258 | higher in CLDN4-high | 0.2754 |
| ELF3(+)/CLDN4out | ELF3 | TJ | cistarget_regulon | 29 | 0.06082 | higher in CLDN4-high | 29 | 0 | 3.73e-09 | 9.63e-08 | SAME | NA | NA | 0.05267 | higher in CLDN4-high | 0.001953 |
| KLF4(+)/CLDN4out | KLF4 | TJ | cistarget_regulon | 29 | 0.02801 | higher in CLDN4-high | 29 | 0 | 3.73e-09 | 9.63e-08 | SAME | NA | NA | 0.02944 | higher in CLDN4-high | 0.001953 |
| KLF5(+)/CLDN4out | KLF5 | TJ | cistarget_regulon | 29 | 0.00766 | higher in CLDN4-high | 23 | 6 | 1.56e-04 | 6.56e-04 | SAME | NA | NA | 0.007575 | higher in CLDN4-high | 0.005859 |
| TFAP2C(+)/CLDN4out | TFAP2C | TJ | cistarget_regulon | 29 | -0.001727 | higher in CLDN4-low | 8 | 21 | 0.02155 | 0.04383 | OPPOSITE | NA | NA | -0.001702 | higher in CLDN4-low | 0.1055 |
| TFAP2A(+)/CLDN4out | TFAP2A | TJ | cistarget_regulon | 29 | -6.70e-04 | higher in CLDN4-low | 14 | 15 | 0.3579 | 0.4261 | OPPOSITE | NA | NA | 0.001958 | higher in CLDN4-high | 0.6953 |
| OVOL2(+)/CLDN4out | OVOL2 | TJ | cistarget_regulon | 29 | -1.24e-04 | higher in CLDN4-low | 13 | 16 | 0.4946 | 0.5561 | OPPOSITE | NA | NA | 0.002218 | higher in CLDN4-high | 0.02734 |
| PROGRAM::TJ_STRUCT |  | TJ | program_AUCell_not_cistarget | 29 | 0.02886 | higher in CLDN4-high | 29 | 0 | 3.73e-09 | 1.30e-08 | SAME | NA | NA | 0.02567 | higher in CLDN4-high | 0.001953 |

## Split check

Within-sample CLDN4 log1p(CP10k), same patient aggregate as the regulon test: median Δ = 1.21, 29/29 patients positive, Wilcoxon p=3.73e-09. The high group is the upper tail of CLDN4 by construction. This row only checks that the labels were applied.

## Honest n

Primary unit = **patient**. Cells from a sample enter that patient only when the sample itself has ≥20 CLDN4-high and ≥20 CLDN4-low malignant cells. High/low is the within-sample median of CLDN4 log1p(CP10k). If that median is 0, high = CLDN4 > 0. Patient means are cell-weighted across that patient's paired samples. Cell counts are inventory, not the test.

| item | n | note |
| --- | --- | --- |
| cells_in_matrix | 208506 | GSE131907 UMI header |
| author_malignant | 31136 | Cell_subtype in {Malignant cells, tS1, tS2, tS3} |
| malignant_UMI_ge_200 | 31136 | QC cells scored by AUCell |
| genes_aucell_universe | 15075 | detected in ≥1% of QC malignant cells |
| grn_cells_subsample | 3000 | stratified by sample, seed 131907 |
| grn_genes | 3067 | HVG + TFs detected in ≥5% + program genes |
| tf_regulators_grnboost2 | 1060 | Lambert/Aerts TF list, detection ≥0.05 |
| cistarget_regulons_NES3 | 253 | activating modules kept by pySCENIC defaults |
| samples_ge20_malignant | 31 | before the high/low tail rule |
| samples_paired_ge20_each_tail | 29 | unit inside each patient |
| samples_split_by_gt0 | 4 | median CLDN4 log1p was 0 |
| patients_paired | 29 | primary Wilcoxon unit |
| paired_patients_single_sample | 29 | one sample supplied both tails |
| paired_patients_one_sample_gt_70pct_cells | 29 | cell-weighted mean can follow that sample |
| patients_tLung_paired | 10 | primary-tumor cells only; split recomputed |

Author-malignant cells by `Sample_Origin`: mBrain 15423, tL/B 6400, tLung 6352, mLN 2961.

mBrain can dominate a pooled cell test. The patient aggregate stops one metastasis from being hundreds of pseudo-replicates, but a patient whose cells are mostly one sample still follows that sample (29 paired patients have one sample with >70% of the paired cells).

Paired patients: P0006, P0008, P0018, P0019, P0020, P0025, P0028, P0030, P0031, P0034, P1006, P1010, P1011, P1012, P1015, P1019, P1028, P1049, P1051, P1058, P3002, P3003, P3004, P3006, P3007, P3012, P3013, P3017, P3019.

## cisTarget regulons (CLDN4 held out)

GRNBoost2 was fit on a stratified subsample of 3000 malignant cells and 3067 genes (log1p CP10k; library size from all genes; seed 131907). The subsample was not restricted to CLDN4-high cells. cisTarget used hg38 v10 cluster rankings for 500 bp upstream / 100 bp downstream and for ±10 kb, motif table `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`, rank threshold 1500, NES ≥ 3.0, pySCENIC default activating modules. AUCell was then run on **all** QC malignant cells, ranking genes detected in ≥1% of those cells. Within a cell, raw-count ranks match log1p(CP10k) ranks.

0 regulons contained CLDN4 before it was removed for the contrast. Full target lists: `tables/regulons.tsv`.

| score | tf | axis | nes | n_targets | elf3_given | n_patients | delta_median_high_minus_low | direction | n_patients_delta_pos | p_patient | fdr_patient | vs_kd_like | delta_median_samples | p_sample | delta_median_tLung | direction_tLung | p_tLung |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ELF3(+)/CLDN4out | ELF3 | TJ | 3.392 | 21 | True | 29 | 0.06082 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | SAME | 0.06082 | 3.73e-09 | 0.05267 | higher in CLDN4-high | 0.001953 |
| FOS(+)/CLDN4out | FOS | IFN_by_targets | 5.328 | 177 | False | 29 | 0.03678 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | OPPOSITE | 0.03678 | 3.73e-09 | 0.03574 | higher in CLDN4-high | 0.001953 |
| FOSB(+)/CLDN4out | FOSB | other | 4.728 | 65 | False | 29 | 0.02952 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | NA | 0.02952 | 3.73e-09 | 0.03248 | higher in CLDN4-high | 0.001953 |
| KLF4(+)/CLDN4out | KLF4 | TJ | 4.168 | 65 | False | 29 | 0.02801 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | SAME | 0.02801 | 3.73e-09 | 0.02944 | higher in CLDN4-high | 0.001953 |
| JUN(+)/CLDN4out | JUN | other | 4.401 | 161 | False | 29 | 0.04004 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | NA | 0.04004 | 3.73e-09 | 0.04015 | higher in CLDN4-high | 0.001953 |
| JUNB(+)/CLDN4out | JUNB | IFN_by_targets | 5.015 | 66 | False | 29 | 0.04625 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | OPPOSITE | 0.04625 | 3.73e-09 | 0.04566 | higher in CLDN4-high | 0.001953 |
| JUND(+)/CLDN4out | JUND | other | 3.202 | 36 | False | 29 | 0.02386 | higher in CLDN4-high | 29 | 3.73e-09 | 9.63e-08 | NA | 0.02386 | 3.73e-09 | 0.03111 | higher in CLDN4-high | 0.001953 |
| ATF3(+)/CLDN4out | ATF3 | IFN_by_targets | 3.639 | 91 | False | 29 | 0.02159 | higher in CLDN4-high | 28 | 7.45e-09 | 1.69e-07 | OPPOSITE | 0.02159 | 7.45e-09 | 0.02179 | higher in CLDN4-high | 0.001953 |
| GATA3(+)/CLDN4out | GATA3 | other | 3.596 | 70 | False | 29 | 0.004462 | higher in CLDN4-high | 28 | 1.12e-08 | 1.69e-07 | NA | 0.004462 | 1.12e-08 | 0.003411 | higher in CLDN4-high | 0.001953 |
| IRF1(+)/CLDN4out | IRF1 | IFN | 3.481 | 32 | False | 29 | 0.01605 | higher in CLDN4-high | 28 | 1.12e-08 | 1.69e-07 | OPPOSITE | 0.01605 | 1.12e-08 | 0.02205 | higher in CLDN4-high | 0.003906 |
| MAFF(+)/CLDN4out | MAFF | other | 4.224 | 41 | False | 29 | 0.0288 | higher in CLDN4-high | 28 | 1.12e-08 | 1.69e-07 | NA | 0.0288 | 1.12e-08 | 0.02647 | higher in CLDN4-high | 0.003906 |
| SOX4(+)/CLDN4out | SOX4 | other | 3.472 | 52 | False | 29 | 0.009978 | higher in CLDN4-high | 28 | 1.12e-08 | 1.69e-07 | NA | 0.009978 | 1.12e-08 | 0.01081 | higher in CLDN4-high | 0.001953 |
| EGR1(+)/CLDN4out | EGR1 | IFN_by_targets | 3.689 | 152 | False | 29 | 0.01741 | higher in CLDN4-high | 27 | 2.61e-08 | 3.37e-07 | OPPOSITE | 0.01741 | 2.61e-08 | 0.01706 | higher in CLDN4-high | 0.003906 |
| RARA(+)/CLDN4out | RARA | other | 4.55 | 17 | False | 29 | 0.01019 | higher in CLDN4-high | 28 | 2.61e-08 | 3.37e-07 | NA | 0.01019 | 2.61e-08 | 0.008206 | higher in CLDN4-high | 0.001953 |
| CEBPD(+)/CLDN4out | CEBPD | IFN_by_targets | 3.544 | 111 | False | 29 | 0.01069 | higher in CLDN4-high | 27 | 7.08e-08 | 8.01e-07 | OPPOSITE | 0.01069 | 7.08e-08 | 0.01106 | higher in CLDN4-high | 0.001953 |
| KLF6(+)/CLDN4out | KLF6 | IFN_by_targets | 4.506 | 212 | False | 29 | 0.007133 | higher in CLDN4-high | 27 | 7.08e-08 | 8.01e-07 | OPPOSITE | 0.007133 | 7.08e-08 | 0.007835 | higher in CLDN4-high | 0.001953 |
| SOX9(+)/CLDN4out | SOX9 | other | 4.446 | 36 | False | 29 | 0.008536 | higher in CLDN4-high | 26 | 2.05e-07 | 2.18e-06 | NA | 0.008536 | 2.05e-07 | 0.007023 | higher in CLDN4-high | 0.003906 |
| ELK3(+)/CLDN4out | ELK3 | IFN_by_targets | 4.039 | 62 | False | 29 | 0.004996 | higher in CLDN4-high | 26 | 3.28e-07 | 3.30e-06 | OPPOSITE | 0.004996 | 3.28e-07 | 0.01033 | higher in CLDN4-high | 0.003906 |
| NR4A1(+)/CLDN4out | NR4A1 | other | 3.163 | 11 | False | 29 | 0.0176 | higher in CLDN4-high | 25 | 4.10e-07 | 3.90e-06 | NA | 0.0176 | 4.10e-07 | 0.01011 | higher in CLDN4-high | 0.04883 |
| NFIA(+)/CLDN4out | NFIA | other | 3.611 | 35 | False | 29 | 0.01105 | higher in CLDN4-high | 27 | 6.30e-07 | 5.70e-06 | NA | 0.01105 | 6.30e-07 | 0.004565 | higher in CLDN4-high | 0.03711 |
| RFXANK(+)/CLDN4out | RFXANK | MHC | 4.251 | 38 | False | 29 | -0.005909 | higher in CLDN4-low | 2 | 7.71e-07 | 6.65e-06 | SAME | -0.005909 | 7.71e-07 | -0.006267 | higher in CLDN4-low | 0.01953 |
| PRDM1(+)/CLDN4out | PRDM1 | other | 4.073 | 10 | False | 29 | 0.006319 | higher in CLDN4-high | 25 | 1.14e-06 | 9.41e-06 | NA | 0.006319 | 1.14e-06 | 0.008412 | higher in CLDN4-high | 0.001953 |
| CLOCK(+)/CLDN4out | CLOCK | other | 3.241 | 21 | False | 29 | -0.002045 | higher in CLDN4-low | 3 | 1.38e-06 | 1.04e-05 | NA | -0.002045 | 1.38e-06 | -0.001778 | higher in CLDN4-low | 0.009766 |
| TEAD4(+)/CLDN4out | TEAD4 | other | 4.246 | 38 | False | 29 | -0.003022 | higher in CLDN4-low | 2 | 1.38e-06 | 1.04e-05 | NA | -0.003022 | 1.38e-06 | -0.001256 | higher in CLDN4-low | 0.06445 |
| MAF(+)/CLDN4out | MAF | other | 4.259 | 86 | False | 29 | 0.006539 | higher in CLDN4-high | 25 | 1.67e-06 | 1.16e-05 | NA | 0.006539 | 1.67e-06 | 0.005389 | higher in CLDN4-high | 0.01367 |
| CREB3L2(+)/CLDN4out | CREB3L2 | other | 3.216 | 11 | False | 29 | 0.005717 | higher in CLDN4-high | 26 | 1.67e-06 | 1.16e-05 | NA | 0.005717 | 1.67e-06 | 0.006021 | higher in CLDN4-high | 0.01953 |
| XBP1(+)/CLDN4out | XBP1 | other | 4.566 | 86 | False | 29 | 0.007141 | higher in CLDN4-high | 25 | 2.38e-06 | 1.60e-05 | NA | 0.007141 | 2.38e-06 | 0.001109 | higher in CLDN4-high | 0.04883 |
| MLX(+)/CLDN4out | MLX | other | 3.841 | 43 | False | 29 | -0.002763 | higher in CLDN4-low | 3 | 3.37e-06 | 2.10e-05 | NA | -0.002763 | 3.37e-06 | -0.003337 | higher in CLDN4-low | 0.001953 |
| EHF(+)/CLDN4out | EHF | IFN_by_targets | 3.945 | 46 | False | 29 | 0.01496 | higher in CLDN4-high | 25 | 3.37e-06 | 2.10e-05 | OPPOSITE | 0.01496 | 3.37e-06 | 0.007928 | higher in CLDN4-high | 0.005859 |
| BHLHE40(+)/CLDN4out | BHLHE40 | IFN_by_targets | 4.555 | 38 | False | 29 | 0.004759 | higher in CLDN4-high | 25 | 3.98e-06 | 2.40e-05 | OPPOSITE | 0.004759 | 3.98e-06 | 0.005026 | higher in CLDN4-high | 0.009766 |
| STAT1(+)/CLDN4out | STAT1 | IFN | 15.22 | 207 | False | 29 | 0.00481 | higher in CLDN4-high | 25 | 7.58e-06 | 4.42e-05 | OPPOSITE | 0.00481 | 7.58e-06 | 0.003928 | higher in CLDN4-high | 0.04883 |
| ELK1(+)/CLDN4out | ELK1 | IFN_by_targets | 5.887 | 41 | False | 29 | -0.001595 | higher in CLDN4-low | 3 | 8.84e-06 | 5.00e-05 | SAME | -0.001595 | 8.84e-06 | -0.001308 | higher in CLDN4-low | 0.005859 |
| FOSL2(+)/CLDN4out | FOSL2 | IFN_by_targets | 4.406 | 59 | False | 29 | 0.003977 | higher in CLDN4-high | 25 | 1.19e-05 | 6.54e-05 | OPPOSITE | 0.003977 | 1.19e-05 | 0.00358 | higher in CLDN4-high | 0.02734 |
| PBX1(+)/CLDN4out | PBX1 | other | 3.907 | 52 | False | 29 | 0.006334 | higher in CLDN4-high | 23 | 1.84e-05 | 9.79e-05 | NA | 0.006334 | 1.84e-05 | 0.006897 | higher in CLDN4-high | 0.01953 |
| E2F1(+)/CLDN4out | E2F1 | IFN_by_targets | 6.411 | 81 | False | 29 | -0.005407 | higher in CLDN4-low | 5 | 2.11e-05 | 1.09e-04 | SAME | -0.005407 | 2.11e-05 | -0.006099 | higher in CLDN4-low | 0.003906 |
| TAF1(+)/CLDN4out | TAF1 | other | 3.774 | 13 | False | 29 | 0.003148 | higher in CLDN4-high | 23 | 4.13e-05 | 2.08e-04 | NA | 0.003148 | 4.13e-05 | 0.003785 | higher in CLDN4-high | 0.003906 |
| STAT3(+)/CLDN4out | STAT3 | other | 4.056 | 16 | False | 29 | 0.003716 | higher in CLDN4-high | 25 | 7.73e-05 | 3.78e-04 | NA | 0.003716 | 7.73e-05 | 0.003219 | higher in CLDN4-high | 0.3223 |
| ZNF580(+)/CLDN4out | ZNF580 | other | 6.085 | 63 | False | 29 | -0.001694 | higher in CLDN4-low | 5 | 1.10e-04 | 5.00e-04 | NA | -0.001694 | 1.10e-04 | -9.81e-04 | higher in CLDN4-low | 0.1055 |
| E2F6(+)/CLDN4out | E2F6 | other | 4.533 | 32 | False | 29 | -0.002029 | higher in CLDN4-low | 6 | 1.10e-04 | 5.00e-04 | NA | -0.002029 | 1.10e-04 | -0.002565 | higher in CLDN4-low | 0.04883 |
| CEBPB(+)/CLDN4out | CEBPB | IFN_by_targets | 4.575 | 113 | False | 29 | 0.0048 | higher in CLDN4-high | 24 | 1.10e-04 | 5.00e-04 | OPPOSITE | 0.0048 | 1.10e-04 | 0.004063 | higher in CLDN4-high | 0.009766 |
| THAP1(+)/CLDN4out | THAP1 | other | 6.45 | 25 | False | 29 | -0.001055 | higher in CLDN4-low | 6 | 1.39e-04 | 6.14e-04 | NA | -0.001055 | 1.39e-04 | -6.51e-04 | higher in CLDN4-low | 0.375 |
| KLF5(+)/CLDN4out | KLF5 | TJ | 3.484 | 42 | False | 29 | 0.00766 | higher in CLDN4-high | 23 | 1.56e-04 | 6.56e-04 | SAME | 0.00766 | 1.56e-04 | 0.007575 | higher in CLDN4-high | 0.005859 |
| GMEB2(+)/CLDN4out | GMEB2 | other | 4.472 | 27 | False | 29 | -0.002089 | higher in CLDN4-low | 8 | 1.56e-04 | 6.56e-04 | NA | -0.002089 | 1.56e-04 | -0.001963 | higher in CLDN4-low | 0.01367 |
| PATZ1(+)/CLDN4out | PATZ1 | other | 4.333 | 39 | False | 29 | -6.80e-04 | higher in CLDN4-low | 7 | 1.74e-04 | 7.01e-04 | NA | -6.80e-04 | 1.74e-04 | -4.93e-04 | higher in CLDN4-low | 0.08398 |
| YY1(+)/CLDN4out | YY1 | other | 5.143 | 182 | False | 29 | 0.002997 | higher in CLDN4-high | 22 | 1.74e-04 | 7.01e-04 | NA | 0.002997 | 1.74e-04 | 0.003689 | higher in CLDN4-high | 0.001953 |
| ZNF331(+)/CLDN4out | ZNF331 | other | 3.564 | 15 | False | 29 | 0.002173 | higher in CLDN4-high | 25 | 1.95e-04 | 7.66e-04 | NA | 0.002173 | 1.95e-04 | 0.002366 | higher in CLDN4-high | 0.009766 |
| RXRB(+)/CLDN4out | RXRB | other | 3.757 | 23 | False | 29 | 0.002465 | higher in CLDN4-high | 23 | 2.17e-04 | 8.36e-04 | NA | 0.002465 | 2.17e-04 | 0.0026 | higher in CLDN4-high | 0.01953 |
| ETS2(+)/CLDN4out | ETS2 | other | 3.3 | 29 | False | 29 | 0.003187 | higher in CLDN4-high | 25 | 2.69e-04 | 9.95e-04 | NA | 0.003187 | 2.69e-04 | 6.92e-04 | higher in CLDN4-high | 0.6953 |
| MAFG(+)/CLDN4out | MAFG | other | 3.775 | 30 | False | 29 | 0.004634 | higher in CLDN4-high | 23 | 2.69e-04 | 9.95e-04 | NA | 0.004634 | 2.69e-04 | 0.006008 | higher in CLDN4-high | 0.003906 |
| HSF1(+)/CLDN4out | HSF1 | other | 4.05 | 31 | False | 29 | 0.004658 | higher in CLDN4-high | 23 | 3.68e-04 | 0.001307 | NA | 0.004658 | 3.68e-04 | 0.005425 | higher in CLDN4-high | 0.009766 |
| ZNF844(+)/CLDN4out | ZNF844 | other | 4.13 | 10 | False | 29 | -9.28e-04 | higher in CLDN4-low | 5 | 3.68e-04 | 0.001307 | NA | -9.28e-04 | 3.68e-04 | -5.93e-04 | higher in CLDN4-low | 0.04883 |
| PPARG(+)/CLDN4out | PPARG | other | 4.695 | 31 | False | 29 | 0.003076 | higher in CLDN4-high | 22 | 4.08e-04 | 0.00142 | NA | 0.003076 | 4.08e-04 | 0.005699 | higher in CLDN4-high | 0.003906 |
| ZBTB7A(+)/CLDN4out | ZBTB7A | other | 3.525 | 41 | False | 29 | 0.002466 | higher in CLDN4-high | 20 | 4.99e-04 | 0.001672 | NA | 0.002466 | 4.99e-04 | 0.002815 | higher in CLDN4-high | 0.009766 |
| ETV6(+)/CLDN4out | ETV6 | other | 3.429 | 18 | False | 29 | 0.002403 | higher in CLDN4-high | 21 | 4.99e-04 | 0.001672 | NA | 0.002403 | 4.99e-04 | 0.00354 | higher in CLDN4-high | 0.04883 |
| ETS1(+)/CLDN4out | ETS1 | other | 4.864 | 23 | False | 29 | 0.001361 | higher in CLDN4-high | 21 | 5.51e-04 | 0.001812 | NA | 0.001361 | 5.51e-04 | 0.002869 | higher in CLDN4-high | 0.01367 |
| ARID3A(+)/CLDN4out | ARID3A | other | 4.425 | 46 | False | 29 | 0.002016 | higher in CLDN4-high | 21 | 7.36e-04 | 0.002379 | NA | 0.002016 | 7.36e-04 | 0.002549 | higher in CLDN4-high | 0.02734 |
| TBL1XR1(+)/CLDN4out | TBL1XR1 | other | 3.792 | 38 | False | 29 | 0.004393 | higher in CLDN4-high | 22 | 8.09e-04 | 0.002569 | NA | 0.004393 | 8.09e-04 | 0.005144 | higher in CLDN4-high | 0.06445 |
| SRF(+)/CLDN4out | SRF | other | 3.616 | 14 | False | 29 | -8.76e-04 | higher in CLDN4-low | 5 | 0.00117 | 0.003651 | NA | -8.76e-04 | 0.00117 | -9.33e-04 | higher in CLDN4-low | 0.01953 |
| ZBTB5(+)/CLDN4out | ZBTB5 | other | 3.253 | 12 | False | 29 | 0.001436 | higher in CLDN4-high | 22 | 0.00128 | 0.003861 | NA | 0.001436 | 0.00128 | 0.002647 | higher in CLDN4-high | 0.1055 |
| STAT2(+)/CLDN4out | STAT2 | IFN | 7.356 | 79 | False | 29 | 0.001502 | higher in CLDN4-high | 21 | 0.00128 | 0.003861 | OPPOSITE | 0.001502 | 0.00128 | 7.91e-04 | higher in CLDN4-high | 0.2324 |
| DDIT3(+)/CLDN4out | DDIT3 | other | 7.336 | 27 | False | 29 | 0.006557 | higher in CLDN4-high | 21 | 0.001399 | 0.00415 | NA | 0.006557 | 0.001399 | 0.005078 | higher in CLDN4-high | 0.1055 |
| SP3(+)/CLDN4out | SP3 | other | 4.278 | 20 | False | 29 | -0.001879 | higher in CLDN4-low | 8 | 0.001527 | 0.004459 | NA | -0.001879 | 0.001527 | -4.42e-04 | higher in CLDN4-low | 0.5566 |
| ZNF528(+)/CLDN4out | ZNF528 | other | 3.474 | 10 | False | 29 | -0.00208 | higher in CLDN4-low | 9 | 0.001978 | 0.005683 | NA | -0.00208 | 0.001978 | 4.36e-04 | higher in CLDN4-high | 0.6953 |
| NFE2L1(+)/CLDN4out | NFE2L1 | other | 4.183 | 67 | False | 29 | 0.001585 | higher in CLDN4-high | 21 | 0.002152 | 0.006087 | NA | 0.001585 | 0.002152 | 0.002986 | higher in CLDN4-high | 0.03711 |
| SOX11(+)/CLDN4out | SOX11 | other | 3.823 | 21 | False | 29 | 0.006554 | higher in CLDN4-high | 20 | 0.002759 | 0.007683 | NA | 0.006554 | 0.002759 | 0.003446 | higher in CLDN4-high | 0.1055 |
| IRF7(+)/CLDN4out | IRF7 | IFN | 15.95 | 129 | False | 29 | 0.002458 | higher in CLDN4-high | 20 | 0.003242 | 0.00876 | OPPOSITE | 0.002458 | 0.003242 | 0.003115 | higher in CLDN4-high | 0.2324 |
| RXRA(+)/CLDN4out | RXRA | other | 3.654 | 23 | False | 29 | 7.92e-04 | higher in CLDN4-high | 23 | 0.003242 | 0.00876 | NA | 7.92e-04 | 0.003242 | 9.18e-04 | higher in CLDN4-high | 0.001953 |
| MAFB(+)/CLDN4out | MAFB | other | 5.192 | 37 | False | 29 | 0.002867 | higher in CLDN4-high | 20 | 0.003511 | 0.00921 | NA | 0.002867 | 0.003511 | 0.002295 | higher in CLDN4-high | 0.2754 |
| HOXA10(+)/CLDN4out | HOXA10 | other | 3.656 | 30 | False | 29 | 0.005223 | higher in CLDN4-high | 22 | 0.003511 | 0.00921 | NA | 0.005223 | 0.003511 | -1.64e-04 | higher in CLDN4-low | 0.5566 |
| THAP11(+)/CLDN4out | THAP11 | other | 5.369 | 23 | False | 29 | -0.003362 | higher in CLDN4-low | 10 | 0.003798 | 0.009822 | NA | -0.003362 | 0.003798 | -0.002292 | higher in CLDN4-low | 0.1309 |
| ZNF254(+)/CLDN4out | ZNF254 | other | 3.381 | 10 | False | 29 | 0.001489 | higher in CLDN4-high | 22 | 0.004106 | 0.01032 | NA | 0.001489 | 0.004106 | 2.48e-04 | higher in CLDN4-high | 0.4316 |
| HOXB7(+)/CLDN4out | HOXB7 | other | 3.871 | 22 | False | 29 | -0.004647 | higher in CLDN4-low | 6 | 0.004106 | 0.01032 | NA | -0.004647 | 0.004106 | -0.004572 | higher in CLDN4-low | 0.06445 |
| NFE2L3(+)/CLDN4out | NFE2L3 | IFN_by_targets | 3.604 | 41 | False | 29 | 0.001147 | higher in CLDN4-high | 19 | 0.004436 | 0.011 | OPPOSITE | 0.001147 | 0.004436 | 0.001339 | higher in CLDN4-high | 0.04883 |
| ETV5(+)/CLDN4out | ETV5 | other | 3.872 | 47 | False | 29 | -0.001985 | higher in CLDN4-low | 6 | 0.004788 | 0.01171 | NA | -0.001985 | 0.004788 | -0.003861 | higher in CLDN4-low | 0.001953 |
| NFATC3(+)/CLDN4out | NFATC3 | other | 3.834 | 44 | False | 29 | 0.00138 | higher in CLDN4-high | 22 | 0.005566 | 0.01343 | NA | 0.00138 | 0.005566 | -8.56e-05 | higher in CLDN4-low | 0.9219 |
| NFIC(+)/CLDN4out | NFIC | IFN_by_targets | 3.09 | 56 | False | 29 | 0.002185 | higher in CLDN4-high | 21 | 0.008008 | 0.01907 | OPPOSITE | 0.002185 | 0.008008 | 0.005695 | higher in CLDN4-high | 0.02734 |
| FOXK1(+)/CLDN4out | FOXK1 | other | 6.233 | 41 | False | 29 | 6.89e-04 | higher in CLDN4-high | 21 | 0.008594 | 0.0202 | NA | 6.89e-04 | 0.008594 | 8.49e-04 | higher in CLDN4-high | 0.04883 |
| GATA2(+)/CLDN4out | GATA2 | other | 4.389 | 14 | False | 29 | -0.001538 | higher in CLDN4-low | 9 | 0.009216 | 0.02139 | NA | -0.001538 | 0.009216 | -1.02e-04 | higher in CLDN4-low | 0.375 |
| REL(+)/CLDN4out | REL | other | 3.09 | 17 | False | 29 | 0.002253 | higher in CLDN4-high | 20 | 0.009877 | 0.02263 | NA | 0.002253 | 0.009877 | 0.003958 | higher in CLDN4-high | 0.06445 |
| FOSL1(+)/CLDN4out | FOSL1 | other | 3.94 | 38 | False | 29 | 0.002209 | higher in CLDN4-high | 21 | 0.01058 | 0.02393 | NA | 0.002209 | 0.01058 | 0.002625 | higher in CLDN4-high | 0.1602 |
| ELF5(+)/CLDN4out | ELF5 | other | 3.463 | 27 | False | 29 | 8.83e-04 | higher in CLDN4-high | 19 | 0.01211 | 0.02673 | NA | 8.83e-04 | 0.01211 | 4.53e-04 | higher in CLDN4-high | 0.6953 |
| ZNF354A(+)/CLDN4out | ZNF354A | other | 3.163 | 16 | False | 29 | -0.002492 | higher in CLDN4-low | 11 | 0.01211 | 0.02673 | NA | -0.002492 | 0.01211 | 3.68e-04 | higher in CLDN4-high | 1 |
| GABPB1(+)/CLDN4out | GABPB1 | other | 4.056 | 15 | False | 29 | -0.001202 | higher in CLDN4-low | 9 | 0.01383 | 0.03015 | NA | -0.001202 | 0.01383 | -8.54e-04 | higher in CLDN4-low | 0.2754 |
| CNOT3(+)/CLDN4out | CNOT3 | other | 4.247 | 12 | False | 29 | -6.39e-04 | higher in CLDN4-low | 8 | 0.01476 | 0.03143 | NA | -6.39e-04 | 0.01476 | -7.50e-04 | higher in CLDN4-low | 0.375 |
| ZNF200(+)/CLDN4out | ZNF200 | other | 4.017 | 12 | False | 29 | -4.75e-04 | higher in CLDN4-low | 10 | 0.01476 | 0.03143 | NA | -4.75e-04 | 0.01476 | -0.001369 | higher in CLDN4-low | 0.02734 |
| ZBTB25(+)/CLDN4out | ZBTB25 | other | 3.149 | 12 | False | 29 | -0.001961 | higher in CLDN4-low | 10 | 0.02027 | 0.04265 | NA | -0.001961 | 0.02027 | 3.34e-04 | higher in CLDN4-high | 0.5566 |
| TFAP2C(+)/CLDN4out | TFAP2C | TJ | 4.072 | 24 | False | 29 | -0.001727 | higher in CLDN4-low | 8 | 0.02155 | 0.04383 | OPPOSITE | -0.001727 | 0.02155 | -0.001702 | higher in CLDN4-low | 0.1055 |
| ILF2(+)/CLDN4out | ILF2 | other | 5.432 | 224 | False | 29 | -0.00101 | higher in CLDN4-low | 8 | 0.02155 | 0.04383 | NA | -0.00101 | 0.02155 | -0.001138 | higher in CLDN4-low | 0.2324 |
| IRF2(+)/CLDN4out | IRF2 | IFN | 4.411 | 47 | False | 29 | -0.001529 | higher in CLDN4-low | 10 | 0.02155 | 0.04383 | SAME | -0.001529 | 0.02155 | -0.001431 | higher in CLDN4-low | 0.1055 |
| KLF10(+)/CLDN4out | KLF10 | other | 3.297 | 19 | False | 29 | 0.001495 | higher in CLDN4-high | 18 | 0.02433 | 0.04894 | NA | 0.001495 | 0.02433 | 5.36e-04 | higher in CLDN4-high | 1 |
| MAX(+)/CLDN4out | MAX | other | 3.161 | 19 | False | 29 | 0.002482 | higher in CLDN4-high | 19 | 0.02583 | 0.05028 | NA | 0.002482 | 0.02583 | 0.002806 | higher in CLDN4-high | 0.1309 |
| TBP(+)/CLDN4out | TBP | other | 3.463 | 14 | False | 29 | -5.18e-04 | higher in CLDN4-low | 8 | 0.02583 | 0.05028 | NA | -5.18e-04 | 0.02583 | -7.37e-05 | higher in CLDN4-low | 0.6953 |
| CREB3(+)/CLDN4out | CREB3 | other | 4.45 | 48 | False | 29 | -0.001866 | higher in CLDN4-low | 9 | 0.02583 | 0.05028 | NA | -0.001866 | 0.02583 | -0.0024 | higher in CLDN4-low | 0.003906 |
| ERF(+)/CLDN4out | ERF | IFN_by_targets | 4.16 | 43 | False | 29 | -3.60e-04 | higher in CLDN4-low | 9 | 0.02906 | 0.05596 | SAME | -3.60e-04 | 0.02906 | -3.88e-04 | higher in CLDN4-low | 0.2754 |
| SIX5(+)/CLDN4out | SIX5 | other | 3.276 | 19 | False | 29 | -2.97e-04 | higher in CLDN4-low | 9 | 0.0308 | 0.05807 | NA | -2.97e-04 | 0.0308 | -4.36e-05 | higher in CLDN4-low | 0.9219 |
| ATF4(+)/CLDN4out | ATF4 | other | 4.801 | 210 | False | 29 | 0.003684 | higher in CLDN4-high | 20 | 0.0308 | 0.05807 | NA | 0.003684 | 0.0308 | 0.001843 | higher in CLDN4-high | 1 |
| ZBTB14(+)/CLDN4out | ZBTB14 | other | 4.16 | 31 | False | 29 | -3.62e-04 | higher in CLDN4-low | 10 | 0.03262 | 0.06087 | NA | -3.62e-04 | 0.03262 | -3.38e-04 | higher in CLDN4-low | 0.1934 |
| USF2(+)/CLDN4out | USF2 | other | 3.208 | 16 | False | 29 | 6.49e-04 | higher in CLDN4-high | 19 | 0.03453 | 0.06378 | NA | 6.49e-04 | 0.03453 | 0.001061 | higher in CLDN4-high | 0.1934 |
| HES6(+)/CLDN4out | HES6 | other | 3.777 | 43 | False | 29 | 0.002099 | higher in CLDN4-high | 20 | 0.03654 | 0.0668 | NA | 0.002099 | 0.03654 | 0.001021 | higher in CLDN4-high | 0.625 |
| PRRX2(+)/CLDN4out | PRRX2 | other | 3.231 | 12 | False | 29 | 0.001656 | higher in CLDN4-high | 19 | 0.04552 | 0.08239 | NA | 0.001656 | 0.04552 | 0.00143 | higher in CLDN4-high | 0.8457 |
| ZNF704(+)/CLDN4out | ZNF704 | other | 3.464 | 10 | False | 29 | 3.53e-04 | higher in CLDN4-high | 18 | 0.05338 | 0.09565 | NA | 3.53e-04 | 0.05338 | 2.81e-04 | higher in CLDN4-high | 0.375 |
| FOXP2(+)/CLDN4out | FOXP2 | other | 4.74 | 24 | False | 29 | 9.48e-04 | higher in CLDN4-high | 18 | 0.05623 | 0.09785 | NA | 9.48e-04 | 0.05623 | 0.001458 | higher in CLDN4-high | 0.1602 |
| SMAD1(+)/CLDN4out | SMAD1 | IFN_by_targets | 3.892 | 57 | False | 29 | 0.001605 | higher in CLDN4-high | 19 | 0.05623 | 0.09785 | OPPOSITE | 0.001605 | 0.05623 | 0.001265 | higher in CLDN4-high | 0.7695 |
| POU2F1(+)/CLDN4out | POU2F1 | other | 3.336 | 12 | False | 29 | 0.001337 | higher in CLDN4-high | 20 | 0.05623 | 0.09785 | NA | 0.001337 | 0.05623 | 9.72e-04 | higher in CLDN4-high | 0.5566 |
| ZBTB7B(+)/CLDN4out | ZBTB7B | other | 5.108 | 31 | False | 29 | 5.06e-04 | higher in CLDN4-high | 17 | 0.0592 | 0.102 | NA | 5.06e-04 | 0.0592 | -7.06e-05 | higher in CLDN4-low | 0.625 |
| TCF7(+)/CLDN4out | TCF7 | other | 3.206 | 11 | False | 29 | -4.52e-04 | higher in CLDN4-low | 9 | 0.06229 | 0.1054 | NA | -4.52e-04 | 0.06229 | -5.23e-04 | higher in CLDN4-low | 0.4316 |
| MZF1(+)/CLDN4out | MZF1 | other | 3.674 | 19 | False | 29 | -7.89e-04 | higher in CLDN4-low | 13 | 0.06229 | 0.1054 | NA | -7.89e-04 | 0.06229 | 9.65e-05 | higher in CLDN4-high | 0.8457 |
| SPDEF(+)/CLDN4out | SPDEF | IFN_by_targets | 5.348 | 75 | False | 29 | 0.004261 | higher in CLDN4-high | 19 | 0.06552 | 0.1098 | OPPOSITE | 0.004261 | 0.06552 | -2.40e-04 | higher in CLDN4-low | 0.9219 |
| NR3C1(+)/CLDN4out | NR3C1 | IFN_by_targets | 3.998 | 68 | False | 29 | 0.001021 | higher in CLDN4-high | 19 | 0.07236 | 0.1202 | OPPOSITE | 0.001021 | 0.07236 | 0.001502 | higher in CLDN4-high | 0.009766 |
| IRX5(+)/CLDN4out | IRX5 | other | 4.32 | 34 | False | 29 | 2.43e-04 | higher in CLDN4-high | 19 | 0.07976 | 0.1301 | NA | 2.43e-04 | 0.07976 | 8.18e-04 | higher in CLDN4-high | 0.06445 |
| ZBTB33(+)/CLDN4out | ZBTB33 | other | 3.036 | 10 | False | 29 | -6.55e-04 | higher in CLDN4-low | 7 | 0.07976 | 0.1301 | NA | -6.55e-04 | 0.07976 | -7.13e-04 | higher in CLDN4-low | 0.02734 |
| NFATC2(+)/CLDN4out | NFATC2 | other | 3.362 | 16 | False | 29 | -6.82e-04 | higher in CLDN4-low | 9 | 0.08774 | 0.1405 | NA | -6.82e-04 | 0.08774 | -5.65e-04 | higher in CLDN4-low | 0.4922 |
| ETV3(+)/CLDN4out | ETV3 | other | 4.654 | 21 | False | 29 | -6.26e-04 | higher in CLDN4-low | 10 | 0.08774 | 0.1405 | NA | -6.26e-04 | 0.08774 | 2.14e-04 | higher in CLDN4-high | 0.7695 |
| RFX5(+)/CLDN4out | RFX5 | MHC | 3.66 | 12 | False | 29 | 0.0011 | higher in CLDN4-high | 19 | 0.09633 | 0.1529 | OPPOSITE | 0.0011 | 0.09633 | 0.002121 | higher in CLDN4-high | 0.08398 |
| IRX3(+)/CLDN4out | IRX3 | other | 3.209 | 17 | False | 29 | 0.002416 | higher in CLDN4-high | 18 | 0.1207 | 0.1899 | NA | 0.002416 | 0.1207 | -0.002658 | higher in CLDN4-low | 0.5566 |
| PURA(+)/CLDN4out | PURA | other | 3.903 | 32 | False | 29 | -8.54e-04 | higher in CLDN4-low | 10 | 0.126 | 0.1967 | NA | -8.54e-04 | 0.126 | -5.53e-04 | higher in CLDN4-low | 0.5566 |
| THRB(+)/CLDN4out | THRB | other | 3.7 | 17 | False | 29 | -3.75e-04 | higher in CLDN4-low | 11 | 0.1433 | 0.2216 | NA | -3.75e-04 | 0.1433 | -2.91e-04 | higher in CLDN4-low | 0.4316 |
| ARID5B(+)/CLDN4out | ARID5B | other | 3.969 | 13 | False | 29 | -7.07e-04 | higher in CLDN4-low | 10 | 0.1494 | 0.2272 | NA | -7.07e-04 | 0.1494 | -2.91e-04 | higher in CLDN4-low | 0.8457 |
| BACH1(+)/CLDN4out | BACH1 | other | 4.274 | 15 | False | 29 | 0.001105 | higher in CLDN4-high | 18 | 0.1494 | 0.2272 | NA | 0.001105 | 0.1494 | 0.001368 | higher in CLDN4-high | 0.375 |
| GABPA(+)/CLDN4out | GABPA | other | 4.504 | 16 | False | 29 | -4.17e-04 | higher in CLDN4-low | 12 | 0.1557 | 0.2329 | NA | -4.17e-04 | 0.1557 | -1.07e-04 | higher in CLDN4-low | 0.7695 |
| CREM(+)/CLDN4out | CREM | other | 3.79 | 12 | False | 29 | 6.84e-04 | higher in CLDN4-high | 18 | 0.1557 | 0.2329 | NA | 6.84e-04 | 0.1557 | -5.66e-05 | higher in CLDN4-low | 0.625 |
| GATAD1(+)/CLDN4out | GATAD1 | other | 3.323 | 28 | False | 29 | 6.98e-04 | higher in CLDN4-high | 17 | 0.1622 | 0.2367 | NA | 6.98e-04 | 0.1622 | 0.002702 | higher in CLDN4-high | 0.2324 |
| ELK4(+)/CLDN4out | ELK4 | other | 3.981 | 39 | False | 29 | 3.03e-04 | higher in CLDN4-high | 17 | 0.1622 | 0.2367 | NA | 3.03e-04 | 0.1622 | 5.12e-04 | higher in CLDN4-high | 0.1934 |
| HMBOX1(+)/CLDN4out | HMBOX1 | other | 3.397 | 11 | False | 29 | 4.75e-04 | higher in CLDN4-high | 18 | 0.1622 | 0.2367 | NA | 4.75e-04 | 0.1622 | 4.47e-04 | higher in CLDN4-high | 0.4922 |
| ETV1(+)/CLDN4out | ETV1 | other | 4.118 | 33 | False | 29 | 4.60e-04 | higher in CLDN4-high | 17 | 0.1689 | 0.2407 | NA | 4.60e-04 | 0.1689 | -0.001348 | higher in CLDN4-low | 0.1934 |
| KLF9(+)/CLDN4out | KLF9 | other | 4.386 | 42 | False | 29 | -6.51e-04 | higher in CLDN4-low | 12 | 0.1689 | 0.2407 | NA | -6.51e-04 | 0.1689 | -2.89e-05 | higher in CLDN4-low | 1 |
| CHD1(+)/CLDN4out | CHD1 | other | 4.209 | 23 | False | 29 | 6.85e-04 | higher in CLDN4-high | 16 | 0.1689 | 0.2407 | NA | 6.85e-04 | 0.1689 | 0.001036 | higher in CLDN4-high | 0.1309 |
| HOXA3(+)/CLDN4out | HOXA3 | other | 3.16 | 17 | False | 29 | 4.13e-04 | higher in CLDN4-high | 18 | 0.1828 | 0.2585 | NA | 4.13e-04 | 0.1828 | 0.002502 | higher in CLDN4-high | 0.1309 |
| HDAC2(+)/CLDN4out | HDAC2 | other | 3.99 | 38 | False | 29 | 0.002377 | higher in CLDN4-high | 17 | 0.1901 | 0.2647 | NA | 0.002377 | 0.1901 | 9.51e-04 | higher in CLDN4-high | 0.8457 |
| HOXB9(+)/CLDN4out | HOXB9 | other | 4.247 | 66 | False | 29 | -0.001989 | higher in CLDN4-low | 10 | 0.1901 | 0.2647 | NA | -0.001989 | 0.1901 | -0.004656 | higher in CLDN4-low | 0.1309 |
| HOXD1(+)/CLDN4out | HOXD1 | other | 3.355 | 51 | False | 29 | 0.001855 | higher in CLDN4-high | 18 | 0.1976 | 0.2731 | NA | 0.001855 | 0.1976 | -0.005514 | higher in CLDN4-low | 0.1934 |
| MYBL2(+)/CLDN4out | MYBL2 | other | 3.615 | 28 | False | 29 | -0.005131 | higher in CLDN4-low | 11 | 0.2053 | 0.2816 | NA | -0.005131 | 0.2053 | -0.007563 | higher in CLDN4-low | 0.2754 |
| ZNF384(+)/CLDN4out | ZNF384 | other | 4.054 | 13 | False | 29 | -6.00e-04 | higher in CLDN4-low | 9 | 0.2297 | 0.3103 | NA | -6.00e-04 | 0.2297 | -7.13e-04 | higher in CLDN4-low | 0.375 |
| STAT6(+)/CLDN4out | STAT6 | other | 3.127 | 21 | False | 29 | -2.73e-04 | higher in CLDN4-low | 12 | 0.2297 | 0.3103 | NA | -2.73e-04 | 0.2297 | -0.001117 | higher in CLDN4-low | 0.3223 |
| KLF3(+)/CLDN4out | KLF3 | other | 4.358 | 39 | False | 29 | 2.11e-04 | higher in CLDN4-high | 15 | 0.2383 | 0.3194 | NA | 2.11e-04 | 0.2383 | 8.42e-04 | higher in CLDN4-high | 0.1602 |
| CEBPA(+)/CLDN4out | CEBPA | other | 4.908 | 44 | False | 29 | 3.40e-04 | higher in CLDN4-high | 15 | 0.256 | 0.3407 | NA | 3.40e-04 | 0.256 | -0.01029 | higher in CLDN4-low | 0.003906 |
| KLF16(+)/CLDN4out | KLF16 | other | 4.683 | 18 | False | 29 | 6.89e-04 | higher in CLDN4-high | 17 | 0.2652 | 0.3478 | NA | 6.89e-04 | 0.2652 | 8.18e-04 | higher in CLDN4-high | 0.08398 |
| ZNF148(+)/CLDN4out | ZNF148 | other | 3.393 | 55 | False | 29 | 6.35e-04 | higher in CLDN4-high | 16 | 0.2652 | 0.3478 | NA | 6.35e-04 | 0.2652 | 0.001188 | higher in CLDN4-high | 0.01367 |
| SP5(+)/CLDN4out | SP5 | other | 3.69 | 12 | False | 29 | -1.86e-04 | higher in CLDN4-low | 12 | 0.2843 | 0.3675 | NA | -1.86e-04 | 0.2843 | -5.00e-04 | higher in CLDN4-low | 0.1934 |
| HLTF(+)/CLDN4out | HLTF | other | 3.75 | 18 | False | 29 | -4.22e-04 | higher in CLDN4-low | 12 | 0.2843 | 0.3675 | NA | -4.22e-04 | 0.2843 | -7.21e-04 | higher in CLDN4-low | 0.6953 |
| TCF7L2(+)/CLDN4out | TCF7L2 | other | 4.321 | 22 | False | 29 | -5.37e-04 | higher in CLDN4-low | 12 | 0.2941 | 0.3723 | NA | -5.37e-04 | 0.2941 | 5.58e-05 | higher in CLDN4-high | 1 |
| ZNF282(+)/CLDN4out | ZNF282 | other | 3.083 | 12 | False | 29 | -3.96e-04 | higher in CLDN4-low | 10 | 0.2941 | 0.3723 | NA | -3.96e-04 | 0.2941 | -5.75e-04 | higher in CLDN4-low | 0.7695 |
| TCF7L1(+)/CLDN4out | TCF7L1 | other | 3.909 | 15 | False | 29 | -3.02e-04 | higher in CLDN4-low | 11 | 0.2941 | 0.3723 | NA | -3.02e-04 | 0.2941 | -5.56e-05 | higher in CLDN4-low | 0.6953 |
| NFIL3(+)/CLDN4out | NFIL3 | other | 4.824 | 10 | False | 29 | -7.59e-04 | higher in CLDN4-low | 10 | 0.3145 | 0.3899 | NA | -7.59e-04 | 0.3145 | 5.75e-04 | higher in CLDN4-high | 0.9219 |
| HOXB8(+)/CLDN4out | HOXB8 | other | 3.89 | 32 | False | 29 | -0.0015 | higher in CLDN4-low | 11 | 0.3145 | 0.3899 | NA | -0.0015 | 0.3145 | -0.003218 | higher in CLDN4-low | 0.03711 |
| BCLAF1(+)/CLDN4out | BCLAF1 | other | 3.304 | 10 | False | 29 | 5.41e-04 | higher in CLDN4-high | 15 | 0.3145 | 0.3899 | NA | 5.41e-04 | 0.3145 | -5.27e-04 | higher in CLDN4-low | 1 |
| ZFX(+)/CLDN4out | ZFX | other | 4.086 | 13 | False | 29 | 1.37e-05 | higher in CLDN4-high | 15 | 0.3357 | 0.4134 | NA | 1.37e-05 | 0.3357 | 0.001363 | higher in CLDN4-high | 0.01367 |
| HOXA5(+)/CLDN4out | HOXA5 | other | 3.455 | 19 | False | 29 | 6.48e-04 | higher in CLDN4-high | 15 | 0.3467 | 0.4211 | NA | 6.48e-04 | 0.3467 | -0.0019 | higher in CLDN4-low | 0.3223 |
| E2F3(+)/CLDN4out | E2F3 | other | 3.71 | 13 | False | 29 | -3.61e-04 | higher in CLDN4-low | 10 | 0.3467 | 0.4211 | NA | -3.61e-04 | 0.3467 | -3.03e-05 | higher in CLDN4-low | 1 |
| RELB(+)/CLDN4out | RELB | other | 5.668 | 29 | False | 29 | -6.03e-04 | higher in CLDN4-low | 12 | 0.3579 | 0.4261 | NA | -6.03e-04 | 0.3579 | 0.001395 | higher in CLDN4-high | 0.06445 |
| NR1D1(+)/CLDN4out | NR1D1 | other | 3.516 | 12 | False | 29 | 1.45e-04 | higher in CLDN4-high | 15 | 0.3579 | 0.4261 | NA | 1.45e-04 | 0.3579 | -3.23e-04 | higher in CLDN4-low | 1 |
| TFAP2A(+)/CLDN4out | TFAP2A | TJ | 3.837 | 11 | False | 29 | -6.70e-04 | higher in CLDN4-low | 14 | 0.3579 | 0.4261 | OPPOSITE | -6.70e-04 | 0.3579 | 0.001958 | higher in CLDN4-high | 0.6953 |
| ZBTB8A(+)/CLDN4out | ZBTB8A | other | 4.24 | 11 | False | 29 | -2.87e-04 | higher in CLDN4-low | 12 | 0.3692 | 0.4368 | NA | -2.87e-04 | 0.3692 | 7.21e-05 | higher in CLDN4-high | 0.8457 |
| EP300(+)/CLDN4out | EP300 | other | 4.094 | 35 | False | 29 | -3.94e-04 | higher in CLDN4-low | 13 | 0.3808 | 0.4476 | NA | -3.94e-04 | 0.3808 | 6.02e-04 | higher in CLDN4-high | 1 |
| NKX2-1(+)/CLDN4out | NKX2-1 | other | 3.67 | 65 | False | 29 | -0.001777 | higher in CLDN4-low | 11 | 0.3927 | 0.4585 | NA | -0.001777 | 0.3927 | -0.004694 | higher in CLDN4-low | 0.08398 |
| ZNF398(+)/CLDN4out | ZNF398 | other | 4.35 | 22 | False | 29 | 1.74e-04 | higher in CLDN4-high | 17 | 0.4169 | 0.4837 | NA | 1.74e-04 | 0.4169 | -6.48e-05 | higher in CLDN4-low | 0.8457 |
| BARX2(+)/CLDN4out | BARX2 | other | 3.381 | 18 | False | 29 | 1.93e-04 | higher in CLDN4-high | 16 | 0.442 | 0.5096 | NA | 1.93e-04 | 0.442 | 4.49e-04 | higher in CLDN4-high | 0.3223 |
| NFATC4(+)/CLDN4out | NFATC4 | other | 3.75 | 16 | False | 29 | 2.10e-04 | higher in CLDN4-high | 16 | 0.4679 | 0.5361 | NA | 2.10e-04 | 0.4679 | 7.52e-04 | higher in CLDN4-high | 0.06445 |
| IRF3(+)/CLDN4out | IRF3 | IFN | 4.495 | 38 | False | 29 | 2.00e-04 | higher in CLDN4-high | 15 | 0.4946 | 0.5561 | OPPOSITE | 2.00e-04 | 0.4946 | 8.89e-05 | higher in CLDN4-high | 0.7695 |
| OVOL2(+)/CLDN4out | OVOL2 | TJ | 3.855 | 24 | False | 29 | -1.24e-04 | higher in CLDN4-low | 13 | 0.4946 | 0.5561 | OPPOSITE | -1.24e-04 | 0.4946 | 0.002218 | higher in CLDN4-high | 0.02734 |
| CHURC1(+)/CLDN4out | CHURC1 | other | 3.67 | 49 | False | 29 | -0.001445 | higher in CLDN4-low | 12 | 0.4946 | 0.5561 | NA | -0.001445 | 0.4946 | -8.31e-04 | higher in CLDN4-low | 0.9219 |
| TFAP4(+)/CLDN4out | TFAP4 | other | 3.786 | 18 | False | 29 | -3.89e-04 | higher in CLDN4-low | 12 | 0.5083 | 0.5679 | NA | -3.89e-04 | 0.5083 | 3.71e-04 | higher in CLDN4-high | 0.9219 |
| NELFB(+)/CLDN4out | NELFB | other | 3.1 | 10 | False | 29 | -5.07e-04 | higher in CLDN4-low | 12 | 0.5647 | 0.6233 | NA | -5.07e-04 | 0.5647 | -2.88e-05 | higher in CLDN4-low | 0.9219 |
| IKZF2(+)/CLDN4out | IKZF2 | IFN_by_targets | 4.392 | 46 | False | 29 | -9.39e-05 | higher in CLDN4-low | 13 | 0.5647 | 0.6233 | SAME | -9.39e-05 | 0.5647 | 2.00e-04 | higher in CLDN4-high | 0.4922 |
| FOXJ3(+)/CLDN4out | FOXJ3 | other | 5.073 | 41 | False | 29 | -6.35e-05 | higher in CLDN4-low | 14 | 0.6391 | 0.6927 | NA | -6.35e-05 | 0.6391 | 1.82e-04 | higher in CLDN4-high | 0.9219 |
| HOXB2(+)/CLDN4out | HOXB2 | other | 3.439 | 66 | False | 29 | -0.0019 | higher in CLDN4-low | 12 | 0.6391 | 0.6927 | NA | -0.0019 | 0.6391 | -0.004271 | higher in CLDN4-low | 0.06445 |
| ZNF467(+)/CLDN4out | ZNF467 | other | 4.425 | 12 | False | 29 | -3.91e-05 | higher in CLDN4-low | 14 | 0.6391 | 0.6927 | NA | -3.91e-05 | 0.6391 | 1.17e-04 | higher in CLDN4-high | 0.625 |
| ETV7(+)/CLDN4out | ETV7 | IFN_by_targets | 10.78 | 63 | False | 29 | -2.23e-04 | higher in CLDN4-low | 12 | 0.6545 | 0.7051 | SAME | -2.23e-04 | 0.6545 | -0.001859 | higher in CLDN4-low | 0.2754 |
| NKX2-8(+)/CLDN4out | NKX2-8 | other | 3.437 | 16 | False | 29 | 0.001231 | higher in CLDN4-high | 16 | 0.6856 | 0.7343 | NA | 0.001231 | 0.6856 | -0.002811 | higher in CLDN4-low | 0.2754 |
| ELF2(+)/CLDN4out | ELF2 | other | 4.269 | 22 | False | 29 | -2.67e-04 | higher in CLDN4-low | 13 | 0.7332 | 0.7807 | NA | -2.67e-04 | 0.7332 | 2.84e-04 | higher in CLDN4-high | 0.7695 |
| HOXB6(+)/CLDN4out | HOXB6 | other | 4.203 | 11 | False | 29 | 4.07e-04 | higher in CLDN4-high | 15 | 0.7493 | 0.7931 | NA | 4.07e-04 | 0.7493 | -3.27e-04 | higher in CLDN4-low | 1 |
| FOXJ2(+)/CLDN4out | FOXJ2 | other | 3.12 | 19 | False | 29 | 3.40e-05 | higher in CLDN4-high | 15 | 0.8148 | 0.8525 | NA | 3.40e-05 | 0.8148 | 2.60e-04 | higher in CLDN4-high | 0.8457 |
| PITX1(+)/CLDN4out | PITX1 | other | 4.465 | 74 | False | 29 | -0.00122 | higher in CLDN4-low | 14 | 0.8148 | 0.8525 | NA | -0.00122 | 0.8148 | -0.001551 | higher in CLDN4-low | 1 |
| ELF1(+)/CLDN4out | ELF1 | IFN_by_targets | 5.096 | 44 | False | 29 | 3.57e-05 | higher in CLDN4-high | 15 | 0.8314 | 0.8648 | OPPOSITE | 3.57e-05 | 0.8314 | -9.46e-05 | higher in CLDN4-low | 0.8457 |
| ZNF189(+)/CLDN4out | ZNF189 | other | 3.035 | 10 | False | 29 | -2.72e-04 | higher in CLDN4-low | 13 | 0.848 | 0.8771 | NA | -2.72e-04 | 0.848 | -0.003082 | higher in CLDN4-low | 0.9219 |
| CEBPG(+)/CLDN4out | CEBPG | IFN_by_targets | 5.609 | 66 | False | 29 | -8.19e-05 | higher in CLDN4-low | 13 | 0.9152 | 0.9412 | SAME | -8.19e-05 | 0.9152 | -6.68e-05 | higher in CLDN4-low | 0.5566 |
| PHF20(+)/CLDN4out | PHF20 | IFN_by_targets | 5.607 | 117 | False | 29 | -6.75e-04 | higher in CLDN4-low | 13 | 0.9321 | 0.9532 | SAME | -6.75e-04 | 0.9321 | -6.80e-04 | higher in CLDN4-low | 0.8457 |
| ETV4(+)/CLDN4out | ETV4 | other | 4.723 | 39 | False | 29 | -6.35e-04 | higher in CLDN4-low | 14 | 0.966 | 0.9714 | NA | -6.35e-04 | 0.966 | -0.002863 | higher in CLDN4-low | 0.7695 |
| GLIS2(+)/CLDN4out | GLIS2 | other | 3.571 | 16 | False | 29 | 2.49e-04 | higher in CLDN4-high | 15 | 0.966 | 0.9714 | NA | 2.49e-04 | 0.966 | 0.001137 | higher in CLDN4-high | 0.375 |
| DPF2(+)/CLDN4out | DPF2 | other | 3.598 | 17 | False | 29 | -4.19e-04 | higher in CLDN4-low | 13 | 0.966 | 0.9714 | NA | -4.19e-04 | 0.966 | 5.48e-04 | higher in CLDN4-high | 0.4316 |
| DBP(+)/CLDN4out | DBP | other | 4.047 | 16 | False | 29 | -3.34e-04 | higher in CLDN4-low | 13 | 0.983 | 0.983 | NA | -3.34e-04 | 0.983 | 0.003224 | higher in CLDN4-high | 0.06445 |

FDR is Benjamini–Hochberg within the CLDN4-held-out cisTarget family. It is a descriptive multiplicity adjustment, not a binding claim.

### As-returned regulons (CLDN4 still inside, if it was a target)

Kept so a sign that appears only when CLDN4 is inside the signature is visible. Not the primary contrast.

| score | tf | n_patients | delta_median_high_minus_low | direction | p_patient | vs_kd_like |
| --- | --- | --- | --- | --- | --- | --- |
| ARID3A(+) | ARID3A | 29 | 0.002016 | higher in CLDN4-high | 7.36e-04 | NA |
| ARID5B(+) | ARID5B | 29 | -7.07e-04 | higher in CLDN4-low | 0.1494 | NA |
| ARNTL(+) | ARNTL | 29 | 1.36e-04 | higher in CLDN4-high | 0.8647 | NA |
| ARX(+) | ARX | 29 | -9.07e-04 | higher in CLDN4-low | 0.01679 | NA |
| ATF1(+) | ATF1 | 29 | -0.006061 | higher in CLDN4-low | 0.02155 | NA |
| ATF3(+) | ATF3 | 29 | 0.02159 | higher in CLDN4-high | 7.45e-09 | OPPOSITE |
| ATF4(+) | ATF4 | 29 | 0.003684 | higher in CLDN4-high | 0.0308 | NA |
| BACH1(+) | BACH1 | 29 | 0.001105 | higher in CLDN4-high | 0.1494 | NA |
| BARX2(+) | BARX2 | 29 | 1.93e-04 | higher in CLDN4-high | 0.442 | NA |
| BCL11A(+) | BCL11A | 29 | 1.99e-05 | higher in CLDN4-high | 0.2941 | NA |
| BCLAF1(+) | BCLAF1 | 29 | 5.41e-04 | higher in CLDN4-high | 0.3145 | NA |
| BHLHE40(+) | BHLHE40 | 29 | 0.004759 | higher in CLDN4-high | 3.98e-06 | OPPOSITE |
| CEBPA(+) | CEBPA | 29 | 3.40e-04 | higher in CLDN4-high | 0.256 | NA |
| CEBPB(+) | CEBPB | 29 | 0.0048 | higher in CLDN4-high | 1.10e-04 | OPPOSITE |
| CEBPD(+) | CEBPD | 29 | 0.01069 | higher in CLDN4-high | 7.08e-08 | OPPOSITE |
| CEBPG(+) | CEBPG | 29 | -8.19e-05 | higher in CLDN4-low | 0.9152 | SAME |
| CHD1(+) | CHD1 | 29 | 6.85e-04 | higher in CLDN4-high | 0.1689 | NA |
| CHURC1(+) | CHURC1 | 29 | -0.001445 | higher in CLDN4-low | 0.4946 | NA |
| CLOCK(+) | CLOCK | 29 | -0.002045 | higher in CLDN4-low | 1.38e-06 | NA |
| CNOT3(+) | CNOT3 | 29 | -6.39e-04 | higher in CLDN4-low | 0.01476 | NA |
| CREB3(+) | CREB3 | 29 | -0.001866 | higher in CLDN4-low | 0.02583 | NA |
| CREB3L2(+) | CREB3L2 | 29 | 0.005717 | higher in CLDN4-high | 1.67e-06 | NA |
| CREM(+) | CREM | 29 | 6.84e-04 | higher in CLDN4-high | 0.1557 | NA |
| DBP(+) | DBP | 29 | -3.34e-04 | higher in CLDN4-low | 0.983 | NA |
| DDIT3(+) | DDIT3 | 29 | 0.006557 | higher in CLDN4-high | 0.001399 | NA |
| DPF2(+) | DPF2 | 29 | -4.19e-04 | higher in CLDN4-low | 0.966 | NA |
| E2F1(+) | E2F1 | 29 | -0.005407 | higher in CLDN4-low | 2.11e-05 | SAME |
| E2F3(+) | E2F3 | 29 | -3.61e-04 | higher in CLDN4-low | 0.3467 | NA |
| E2F6(+) | E2F6 | 29 | -0.002029 | higher in CLDN4-low | 1.10e-04 | NA |
| E4F1(+) | E4F1 | 29 | -5.00e-04 | higher in CLDN4-low | 0.004436 | NA |
| EGR1(+) | EGR1 | 29 | 0.01741 | higher in CLDN4-high | 2.61e-08 | OPPOSITE |
| EHF(+) | EHF | 29 | 0.01496 | higher in CLDN4-high | 3.37e-06 | OPPOSITE |
| ELF1(+) | ELF1 | 29 | 3.57e-05 | higher in CLDN4-high | 0.8314 | OPPOSITE |
| ELF2(+) | ELF2 | 29 | -2.67e-04 | higher in CLDN4-low | 0.7332 | NA |
| ELF3(+) | ELF3 | 29 | 0.06082 | higher in CLDN4-high | 3.73e-09 | SAME |
| ELF5(+) | ELF5 | 29 | 8.83e-04 | higher in CLDN4-high | 0.01211 | NA |
| ELK1(+) | ELK1 | 29 | -0.001595 | higher in CLDN4-low | 8.84e-06 | SAME |
| ELK3(+) | ELK3 | 29 | 0.004996 | higher in CLDN4-high | 3.28e-07 | OPPOSITE |
| ELK4(+) | ELK4 | 29 | 3.03e-04 | higher in CLDN4-high | 0.1622 | NA |
| EP300(+) | EP300 | 29 | -3.94e-04 | higher in CLDN4-low | 0.3808 | NA |
| EPAS1(+) | EPAS1 | 29 | -0.004277 | higher in CLDN4-low | 0.004106 | NA |
| ERF(+) | ERF | 29 | -3.60e-04 | higher in CLDN4-low | 0.02906 | SAME |
| ETS1(+) | ETS1 | 29 | 0.001361 | higher in CLDN4-high | 5.51e-04 | NA |
| ETS2(+) | ETS2 | 29 | 0.003187 | higher in CLDN4-high | 2.69e-04 | NA |
| ETV1(+) | ETV1 | 29 | 4.60e-04 | higher in CLDN4-high | 0.1689 | NA |
| ETV3(+) | ETV3 | 29 | -6.26e-04 | higher in CLDN4-low | 0.08774 | NA |
| ETV4(+) | ETV4 | 29 | -6.35e-04 | higher in CLDN4-low | 0.966 | NA |
| ETV5(+) | ETV5 | 29 | -0.001985 | higher in CLDN4-low | 0.004788 | NA |
| ETV6(+) | ETV6 | 29 | 0.002403 | higher in CLDN4-high | 4.99e-04 | NA |
| ETV7(+) | ETV7 | 29 | -2.23e-04 | higher in CLDN4-low | 0.6545 | SAME |
| FOS(+) | FOS | 29 | 0.03678 | higher in CLDN4-high | 3.73e-09 | OPPOSITE |
| FOSB(+) | FOSB | 29 | 0.02952 | higher in CLDN4-high | 3.73e-09 | NA |
| FOSL1(+) | FOSL1 | 29 | 0.002209 | higher in CLDN4-high | 0.01058 | NA |
| FOSL2(+) | FOSL2 | 29 | 0.003977 | higher in CLDN4-high | 1.19e-05 | OPPOSITE |
| FOXJ2(+) | FOXJ2 | 29 | 3.40e-05 | higher in CLDN4-high | 0.8148 | NA |
| FOXJ3(+) | FOXJ3 | 29 | -6.35e-05 | higher in CLDN4-low | 0.6391 | NA |
| FOXK1(+) | FOXK1 | 29 | 6.89e-04 | higher in CLDN4-high | 0.008594 | NA |
| FOXK2(+) | FOXK2 | 29 | 0 | no difference | 0.648 | NA |
| FOXN2(+) | FOXN2 | 29 | -2.02e-04 | higher in CLDN4-low | 0.2746 | NA |
| FOXO3(+) | FOXO3 | 29 | 4.66e-04 | higher in CLDN4-high | 0.1207 | NA |
| FOXP1(+) | FOXP1 | 29 | -1.96e-04 | higher in CLDN4-low | 0.7983 | NA |
| FOXP2(+) | FOXP2 | 29 | 9.48e-04 | higher in CLDN4-high | 0.05623 | NA |
| GABPA(+) | GABPA | 29 | -4.17e-04 | higher in CLDN4-low | 0.1557 | NA |
| GABPB1(+) | GABPB1 | 29 | -0.001202 | higher in CLDN4-low | 0.01383 | NA |
| GATA2(+) | GATA2 | 29 | -0.001538 | higher in CLDN4-low | 0.009216 | NA |
| GATA3(+) | GATA3 | 29 | 0.004462 | higher in CLDN4-high | 1.12e-08 | NA |
| GATAD1(+) | GATAD1 | 29 | 6.98e-04 | higher in CLDN4-high | 0.1622 | NA |
| GLIS2(+) | GLIS2 | 29 | 2.49e-04 | higher in CLDN4-high | 0.966 | NA |
| GMEB2(+) | GMEB2 | 29 | -0.002089 | higher in CLDN4-low | 1.56e-04 | NA |
| HDAC2(+) | HDAC2 | 29 | 0.002377 | higher in CLDN4-high | 0.1901 | NA |
| HDAC8(+) | HDAC8 | 29 | -0.005466 | higher in CLDN4-low | 0.1316 | NA |
| HES6(+) | HES6 | 29 | 0.002099 | higher in CLDN4-high | 0.03654 | NA |
| HLF(+) | HLF | 29 | -0.002135 | higher in CLDN4-low | 0.04803 | NA |
| HLTF(+) | HLTF | 29 | -4.22e-04 | higher in CLDN4-low | 0.2843 | NA |
| HMBOX1(+) | HMBOX1 | 29 | 4.75e-04 | higher in CLDN4-high | 0.1622 | NA |
| HMG20B(+) | HMG20B | 29 | -0.002331 | higher in CLDN4-low | 0.003798 | NA |
| HMGA2(+) | HMGA2 | 29 | -5.63e-05 | higher in CLDN4-low | 0.8376 | NA |
| HNF1B(+) | HNF1B | 29 | 3.89e-04 | higher in CLDN4-high | 0.5504 | NA |
| HOXA10(+) | HOXA10 | 29 | 0.005223 | higher in CLDN4-high | 0.003511 | NA |
| HOXA3(+) | HOXA3 | 29 | 4.13e-04 | higher in CLDN4-high | 0.1828 | NA |
| HOXA5(+) | HOXA5 | 29 | 6.48e-04 | higher in CLDN4-high | 0.3467 | NA |
| HOXB2(+) | HOXB2 | 29 | -0.0019 | higher in CLDN4-low | 0.6391 | NA |
| HOXB3(+) | HOXB3 | 29 | 0 | no difference | 0.2213 | NA |
| HOXB5(+) | HOXB5 | 29 | -1.62e-05 | higher in CLDN4-low | 0.9491 | NA |
| HOXB6(+) | HOXB6 | 29 | 4.07e-04 | higher in CLDN4-high | 0.7493 | NA |
| HOXB7(+) | HOXB7 | 29 | -0.004647 | higher in CLDN4-low | 0.004106 | NA |
| HOXB8(+) | HOXB8 | 29 | -0.0015 | higher in CLDN4-low | 0.3145 | NA |
| HOXB9(+) | HOXB9 | 29 | -0.001989 | higher in CLDN4-low | 0.1901 | NA |
| HOXC9(+) | HOXC9 | 29 | -0.001206 | higher in CLDN4-low | 0.01383 | NA |
| HOXD1(+) | HOXD1 | 29 | 0.001855 | higher in CLDN4-high | 0.1976 | NA |
| HSF1(+) | HSF1 | 29 | 0.004658 | higher in CLDN4-high | 3.68e-04 | NA |
| IKZF2(+) | IKZF2 | 29 | -9.39e-05 | higher in CLDN4-low | 0.5647 | SAME |
| ILF2(+) | ILF2 | 29 | -0.00101 | higher in CLDN4-low | 0.02155 | NA |
| IRF1(+) | IRF1 | 29 | 0.01605 | higher in CLDN4-high | 1.12e-08 | OPPOSITE |
| IRF2(+) | IRF2 | 29 | -0.001529 | higher in CLDN4-low | 0.02155 | SAME |
| IRF3(+) | IRF3 | 29 | 2.00e-04 | higher in CLDN4-high | 0.4946 | OPPOSITE |
| IRF5(+) | IRF5 | 29 | -3.15e-04 | higher in CLDN4-low | 0.442 | SAME |
| IRF6(+) | IRF6 | 29 | -0.001936 | higher in CLDN4-low | 3.68e-04 | SAME |
| IRF7(+) | IRF7 | 29 | 0.002458 | higher in CLDN4-high | 0.003242 | OPPOSITE |
| IRF9(+) | IRF9 | 29 | 0.004201 | higher in CLDN4-high | 2.69e-04 | OPPOSITE |
| IRX3(+) | IRX3 | 29 | 0.002416 | higher in CLDN4-high | 0.1207 | NA |
| IRX5(+) | IRX5 | 29 | 2.43e-04 | higher in CLDN4-high | 0.07976 | NA |
| JUN(+) | JUN | 29 | 0.04004 | higher in CLDN4-high | 3.73e-09 | NA |
| JUNB(+) | JUNB | 29 | 0.04625 | higher in CLDN4-high | 3.73e-09 | OPPOSITE |
| JUND(+) | JUND | 29 | 0.02386 | higher in CLDN4-high | 3.73e-09 | NA |
| KLF10(+) | KLF10 | 29 | 0.001495 | higher in CLDN4-high | 0.02433 | NA |
| KLF11(+) | KLF11 | 29 | -4.62e-04 | higher in CLDN4-low | 0.3145 | NA |
| KLF16(+) | KLF16 | 29 | 6.89e-04 | higher in CLDN4-high | 0.2652 | NA |
| KLF3(+) | KLF3 | 29 | 2.11e-04 | higher in CLDN4-high | 0.2383 | NA |
| KLF4(+) | KLF4 | 29 | 0.02801 | higher in CLDN4-high | 3.73e-09 | SAME |
| KLF5(+) | KLF5 | 29 | 0.00766 | higher in CLDN4-high | 1.56e-04 | SAME |
| KLF6(+) | KLF6 | 29 | 0.007133 | higher in CLDN4-high | 7.08e-08 | OPPOSITE |
| KLF7(+) | KLF7 | 29 | 0.007342 | higher in CLDN4-high | 9.75e-04 | NA |
| KLF9(+) | KLF9 | 29 | -6.51e-04 | higher in CLDN4-low | 0.1689 | NA |
| MAF(+) | MAF | 29 | 0.006539 | higher in CLDN4-high | 1.67e-06 | NA |
| MAFB(+) | MAFB | 29 | 0.002867 | higher in CLDN4-high | 0.003511 | NA |
| MAFF(+) | MAFF | 29 | 0.0288 | higher in CLDN4-high | 1.12e-08 | NA |
| MAFG(+) | MAFG | 29 | 0.004634 | higher in CLDN4-high | 2.69e-04 | NA |
| MAFK(+) | MAFK | 29 | 0.006232 | higher in CLDN4-high | 2.69e-04 | NA |
| MAX(+) | MAX | 29 | 0.002482 | higher in CLDN4-high | 0.02583 | NA |
| MAZ(+) | MAZ | 29 | -9.39e-06 | higher in CLDN4-low | 0.4549 | NA |
| MECP2(+) | MECP2 | 29 | -2.62e-04 | higher in CLDN4-low | 0.594 | NA |
| MLX(+) | MLX | 29 | -0.002763 | higher in CLDN4-low | 3.37e-06 | NA |
| MXD4(+) | MXD4 | 29 | -0.003125 | higher in CLDN4-low | 0.05064 | NA |
| MYBL2(+) | MYBL2 | 29 | -0.005131 | higher in CLDN4-low | 0.2053 | NA |
| MZF1(+) | MZF1 | 29 | -7.89e-04 | higher in CLDN4-low | 0.06229 | NA |
| NCALD(+) | NCALD | 29 | -1.51e-04 | higher in CLDN4-low | 0.2214 | NA |
| NELFB(+) | NELFB | 29 | -5.07e-04 | higher in CLDN4-low | 0.5647 | NA |
| NFATC2(+) | NFATC2 | 29 | -6.82e-04 | higher in CLDN4-low | 0.08774 | NA |
| NFATC3(+) | NFATC3 | 29 | 0.00138 | higher in CLDN4-high | 0.005566 | NA |
| NFATC4(+) | NFATC4 | 29 | 2.10e-04 | higher in CLDN4-high | 0.4679 | NA |
| NFE2L1(+) | NFE2L1 | 29 | 0.001585 | higher in CLDN4-high | 0.002152 | NA |
| NFE2L3(+) | NFE2L3 | 29 | 0.001147 | higher in CLDN4-high | 0.004436 | OPPOSITE |
| NFIA(+) | NFIA | 29 | 0.01105 | higher in CLDN4-high | 6.30e-07 | NA |
| NFIB(+) | NFIB | 29 | 0.005945 | higher in CLDN4-high | 7.36e-04 | NA |
| NFIC(+) | NFIC | 29 | 0.002185 | higher in CLDN4-high | 0.008008 | OPPOSITE |
| NFIL3(+) | NFIL3 | 29 | -7.59e-04 | higher in CLDN4-low | 0.3145 | NA |
| NFKB2(+) | NFKB2 | 29 | 4.40e-04 | higher in CLDN4-high | 0.1494 | NA |
| NFYA(+) | NFYA | 29 | -5.79e-04 | higher in CLDN4-low | 0.03654 | NA |
| NKX2-1(+) | NKX2-1 | 29 | -0.001777 | higher in CLDN4-low | 0.3927 | NA |
| NKX2-8(+) | NKX2-8 | 29 | 0.001231 | higher in CLDN4-high | 0.6856 | NA |
| NR1D1(+) | NR1D1 | 29 | 1.45e-04 | higher in CLDN4-high | 0.3579 | NA |
| NR1H2(+) | NR1H2 | 29 | -0.004392 | higher in CLDN4-low | 0.0592 | NA |
| NR2C2(+) | NR2C2 | 29 | -6.75e-04 | higher in CLDN4-low | 0.02433 | NA |
| NR2F2(+) | NR2F2 | 29 | 0 | no difference | 0.3246 | NA |
| NR3C1(+) | NR3C1 | 29 | 0.001021 | higher in CLDN4-high | 0.07236 | OPPOSITE |
| NR4A1(+) | NR4A1 | 29 | 0.0176 | higher in CLDN4-high | 4.10e-07 | NA |
| OVOL2(+) | OVOL2 | 29 | -1.24e-04 | higher in CLDN4-low | 0.4946 | OPPOSITE |
| PATZ1(+) | PATZ1 | 29 | -6.80e-04 | higher in CLDN4-low | 1.74e-04 | NA |
| PBX1(+) | PBX1 | 29 | 0.006334 | higher in CLDN4-high | 1.84e-05 | NA |
| PHF20(+) | PHF20 | 29 | -6.75e-04 | higher in CLDN4-low | 0.9321 | SAME |
| PITX1(+) | PITX1 | 29 | -0.00122 | higher in CLDN4-low | 0.8148 | NA |
| PLAGL2(+) | PLAGL2 | 29 | 0.001272 | higher in CLDN4-high | 0.2652 | NA |
| POU2F1(+) | POU2F1 | 29 | 0.001337 | higher in CLDN4-high | 0.05623 | NA |
| PPARG(+) | PPARG | 29 | 0.003076 | higher in CLDN4-high | 4.08e-04 | NA |
| PRDM1(+) | PRDM1 | 29 | 0.006319 | higher in CLDN4-high | 1.14e-06 | NA |
| PRRX2(+) | PRRX2 | 29 | 0.001656 | higher in CLDN4-high | 0.04552 | NA |
| PURA(+) | PURA | 29 | -8.54e-04 | higher in CLDN4-low | 0.126 | NA |
| RARA(+) | RARA | 29 | 0.01019 | higher in CLDN4-high | 2.61e-08 | NA |
| REL(+) | REL | 29 | 0.002253 | higher in CLDN4-high | 0.009877 | NA |
| RELA(+) | RELA | 29 | -0.002373 | higher in CLDN4-low | 0.001666 | NA |
| RELB(+) | RELB | 29 | -6.03e-04 | higher in CLDN4-low | 0.3579 | NA |
| REST(+) | REST | 29 | -1.35e-04 | higher in CLDN4-low | 0.374 | NA |
| RFX2(+) | RFX2 | 29 | -3.30e-05 | higher in CLDN4-low | 0.6164 | NA |
| RFX5(+) | RFX5 | 29 | 0.0011 | higher in CLDN4-high | 0.09633 | OPPOSITE |
| RFXANK(+) | RFXANK | 29 | -0.005909 | higher in CLDN4-low | 7.71e-07 | SAME |
| RXRA(+) | RXRA | 29 | 7.92e-04 | higher in CLDN4-high | 0.003242 | NA |
| RXRB(+) | RXRB | 29 | 0.002465 | higher in CLDN4-high | 2.17e-04 | NA |
| SIRT6(+) | SIRT6 | 29 | -6.48e-04 | higher in CLDN4-low | 0.1161 | NA |
| SIX5(+) | SIX5 | 29 | -2.97e-04 | higher in CLDN4-low | 0.0308 | NA |
| SMAD1(+) | SMAD1 | 29 | 0.001605 | higher in CLDN4-high | 0.05623 | OPPOSITE |
| SOX11(+) | SOX11 | 29 | 0.006554 | higher in CLDN4-high | 0.002759 | NA |
| SOX4(+) | SOX4 | 29 | 0.009978 | higher in CLDN4-high | 1.12e-08 | NA |
| SOX9(+) | SOX9 | 29 | 0.008536 | higher in CLDN4-high | 2.05e-07 | NA |
| SP3(+) | SP3 | 29 | -0.001879 | higher in CLDN4-low | 0.001527 | NA |
| SP4(+) | SP4 | 29 | 6.10e-04 | higher in CLDN4-high | 0.3467 | NA |
| SP5(+) | SP5 | 29 | -1.86e-04 | higher in CLDN4-low | 0.2843 | NA |
| SPDEF(+) | SPDEF | 29 | 0.004261 | higher in CLDN4-high | 0.06552 | OPPOSITE |
| SREBF1(+) | SREBF1 | 29 | -0.00193 | higher in CLDN4-low | 1.60e-05 | NA |
| SRF(+) | SRF | 29 | -8.76e-04 | higher in CLDN4-low | 0.00117 | NA |
| STAT1(+) | STAT1 | 29 | 0.00481 | higher in CLDN4-high | 7.58e-06 | OPPOSITE |
| STAT2(+) | STAT2 | 29 | 0.001502 | higher in CLDN4-high | 0.00128 | OPPOSITE |
| STAT3(+) | STAT3 | 29 | 0.003716 | higher in CLDN4-high | 7.73e-05 | NA |
| STAT6(+) | STAT6 | 29 | -2.73e-04 | higher in CLDN4-low | 0.2297 | NA |
| TAF1(+) | TAF1 | 29 | 0.003148 | higher in CLDN4-high | 4.13e-05 | NA |
| TAF6(+) | TAF6 | 29 | -0.003963 | higher in CLDN4-low | 3.98e-06 | NA |
| TBL1XR1(+) | TBL1XR1 | 29 | 0.004393 | higher in CLDN4-high | 8.09e-04 | NA |
| TBP(+) | TBP | 29 | -5.18e-04 | higher in CLDN4-low | 0.02583 | NA |
| TCF7(+) | TCF7 | 29 | -4.52e-04 | higher in CLDN4-low | 0.06229 | NA |
| TCF7L1(+) | TCF7L1 | 29 | -3.02e-04 | higher in CLDN4-low | 0.2941 | NA |
| TCF7L2(+) | TCF7L2 | 29 | -5.37e-04 | higher in CLDN4-low | 0.2941 | NA |
| TEAD4(+) | TEAD4 | 29 | -0.003022 | higher in CLDN4-low | 1.38e-06 | NA |
| TEF(+) | TEF | 29 | -2.26e-05 | higher in CLDN4-low | 0.7493 | NA |
| TFAP2A(+) | TFAP2A | 29 | -6.70e-04 | higher in CLDN4-low | 0.3579 | OPPOSITE |
| TFAP2C(+) | TFAP2C | 29 | -0.001727 | higher in CLDN4-low | 0.02155 | OPPOSITE |
| TFAP4(+) | TFAP4 | 29 | -3.89e-04 | higher in CLDN4-low | 0.5083 | NA |
| TFDP1(+) | TFDP1 | 29 | -0.006162 | higher in CLDN4-low | 1.10e-04 | NA |
| THAP1(+) | THAP1 | 29 | -0.001055 | higher in CLDN4-low | 1.39e-04 | NA |
| THAP11(+) | THAP11 | 29 | -0.003362 | higher in CLDN4-low | 0.003798 | NA |
| THRA(+) | THRA | 29 | -7.58e-04 | higher in CLDN4-low | 0.01679 | NA |
| THRB(+) | THRB | 29 | -3.75e-04 | higher in CLDN4-low | 0.1433 | NA |
| TP53(+) | TP53 | 29 | 9.32e-04 | higher in CLDN4-high | 0.3357 | NA |
| USF2(+) | USF2 | 29 | 6.49e-04 | higher in CLDN4-high | 0.03453 | NA |
| VDR(+) | VDR | 29 | -2.15e-05 | higher in CLDN4-low | 0.4294 | NA |
| WRNIP1(+) | WRNIP1 | 29 | -4.46e-04 | higher in CLDN4-low | 0.1901 | NA |
| XBP1(+) | XBP1 | 29 | 0.007141 | higher in CLDN4-high | 2.38e-06 | NA |
| YY1(+) | YY1 | 29 | 0.002997 | higher in CLDN4-high | 1.74e-04 | NA |
| ZBED1(+) | ZBED1 | 29 | -3.21e-04 | higher in CLDN4-low | 0.1622 | NA |
| ZBTB14(+) | ZBTB14 | 29 | -3.62e-04 | higher in CLDN4-low | 0.03262 | NA |
| ZBTB17(+) | ZBTB17 | 29 | 0 | no difference | 0.6186 | NA |
| ZBTB25(+) | ZBTB25 | 29 | -0.001961 | higher in CLDN4-low | 0.02027 | NA |
| ZBTB33(+) | ZBTB33 | 29 | -6.55e-04 | higher in CLDN4-low | 0.07976 | NA |
| ZBTB4(+) | ZBTB4 | 29 | -0.001362 | higher in CLDN4-low | 0.3927 | NA |
| ZBTB5(+) | ZBTB5 | 29 | 0.001436 | higher in CLDN4-high | 0.00128 | NA |
| ZBTB7A(+) | ZBTB7A | 29 | 0.002466 | higher in CLDN4-high | 4.99e-04 | NA |
| ZBTB7B(+) | ZBTB7B | 29 | 5.06e-04 | higher in CLDN4-high | 0.0592 | NA |
| ZBTB8A(+) | ZBTB8A | 29 | -2.87e-04 | higher in CLDN4-low | 0.3692 | NA |
| ZFP62(+) | ZFP62 | 29 | 5.66e-04 | higher in CLDN4-high | 0.2652 | NA |
| ZFP91(+) | ZFP91 | 29 | -0.002105 | higher in CLDN4-low | 1.39e-04 | NA |
| ZFX(+) | ZFX | 29 | 1.37e-05 | higher in CLDN4-high | 0.3357 | NA |
| ZNF12(+) | ZNF12 | 29 | 3.94e-04 | higher in CLDN4-high | 0.1828 | NA |
| ZNF131(+) | ZNF131 | 29 | -3.64e-04 | higher in CLDN4-low | 0.145 | NA |
| ZNF133(+) | ZNF133 | 29 | 1.86e-04 | higher in CLDN4-high | 0.442 | NA |
| ZNF134(+) | ZNF134 | 29 | -2.68e-05 | higher in CLDN4-low | 0.7013 | NA |
| ZNF148(+) | ZNF148 | 29 | 6.35e-04 | higher in CLDN4-high | 0.2652 | NA |
| ZNF189(+) | ZNF189 | 29 | -2.72e-04 | higher in CLDN4-low | 0.848 | NA |
| ZNF200(+) | ZNF200 | 29 | -4.75e-04 | higher in CLDN4-low | 0.01476 | NA |
| ZNF217(+) | ZNF217 | 29 | -0.001619 | higher in CLDN4-low | 0.001816 | NA |
| ZNF224(+) | ZNF224 | 29 | -8.17e-04 | higher in CLDN4-low | 0.07976 | NA |
| ZNF254(+) | ZNF254 | 29 | 0.001489 | higher in CLDN4-high | 0.004106 | NA |
| ZNF282(+) | ZNF282 | 29 | -3.96e-04 | higher in CLDN4-low | 0.2941 | NA |
| ZNF331(+) | ZNF331 | 29 | 0.002173 | higher in CLDN4-high | 1.95e-04 | NA |
| ZNF354A(+) | ZNF354A | 29 | -0.002492 | higher in CLDN4-low | 0.01211 | NA |
| ZNF37A(+) | ZNF37A | 29 | -5.91e-04 | higher in CLDN4-low | 0.1976 | NA |
| ZNF384(+) | ZNF384 | 29 | -6.00e-04 | higher in CLDN4-low | 0.2297 | NA |
| ZNF398(+) | ZNF398 | 29 | 1.74e-04 | higher in CLDN4-high | 0.4169 | NA |
| ZNF426(+) | ZNF426 | 29 | -8.32e-04 | higher in CLDN4-low | 0.1316 | NA |
| ZNF467(+) | ZNF467 | 29 | -3.91e-05 | higher in CLDN4-low | 0.6391 | NA |
| ZNF518A(+) | ZNF518A | 29 | -2.13e-04 | higher in CLDN4-low | 0.5793 | NA |
| ZNF528(+) | ZNF528 | 29 | -0.00208 | higher in CLDN4-low | 0.001978 | NA |
| ZNF529(+) | ZNF529 | 29 | -8.78e-04 | higher in CLDN4-low | 0.1433 | NA |
| ZNF580(+) | ZNF580 | 29 | -0.001694 | higher in CLDN4-low | 1.10e-04 | NA |
| ZNF672(+) | ZNF672 | 29 | -0.003808 | higher in CLDN4-low | 5.34e-05 | NA |
| ZNF692(+) | ZNF692 | 29 | 8.99e-04 | higher in CLDN4-high | 0.7172 | NA |
| ZNF704(+) | ZNF704 | 29 | 3.53e-04 | higher in CLDN4-high | 0.05338 | NA |
| ZNF740(+) | ZNF740 | 29 | -4.04e-04 | higher in CLDN4-low | 0.4294 | NA |
| ZNF778(+) | ZNF778 | 29 | 0.001964 | higher in CLDN4-high | 0.02027 | NA |
| ZNF785(+) | ZNF785 | 29 | 6.17e-04 | higher in CLDN4-high | 0.02906 | NA |
| ZNF787(+) | ZNF787 | 29 | 2.54e-04 | higher in CLDN4-high | 0.2746 | NA |
| ZNF844(+) | ZNF844 | 29 | -9.28e-04 | higher in CLDN4-low | 3.68e-04 | NA |
| ZNF891(+) | ZNF891 | 29 | -2.14e-04 | higher in CLDN4-low | 0.1104 | NA |
| ZSCAN29(+) | ZSCAN29 | 29 | -9.97e-04 | higher in CLDN4-low | 0.06887 | NA |
| ZSCAN30(+) | ZSCAN30 | 29 | 0.004377 | higher in CLDN4-high | 0.003511 | NA |

## tLung sensitivity

Same rules, restricted to `Sample_Origin == tLung` before the median split. Columns `delta_median_tLung`, `direction_tLung`, and `p_tLung` are in the tables above. A sign that flips between the all-site patient test and tLung is a real disagreement and is not averaged away.

## Between-patient companion (not the high vs low split)

Spearman of the patient's malignant CLDN4+ fraction versus the mean AUCell of all that patient's malignant cells. This is the “CLDN4-high patients” question. It is not the within-sample high versus low contrast, and the two can disagree.

| score | n_patients | spearman_rho_frac_cldn4_pos | p | question |
| --- | --- | --- | --- | --- |
| ATF3(+)/CLDN4out | 31 | 0.2835 | 0.1223 | between_patient_not_the_within_sample_split |
| BHLHE40(+)/CLDN4out | 31 | 0.0496 | 0.791 | between_patient_not_the_within_sample_split |
| CEBPB(+)/CLDN4out | 31 | 0.1508 | 0.4181 | between_patient_not_the_within_sample_split |
| CEBPD(+)/CLDN4out | 31 | 0.3544 | 0.05042 | between_patient_not_the_within_sample_split |
| CEBPG(+)/CLDN4out | 31 | -0.2621 | 0.1543 | between_patient_not_the_within_sample_split |
| E2F1(+)/CLDN4out | 31 | -0.1153 | 0.5367 | between_patient_not_the_within_sample_split |
| EGR1(+)/CLDN4out | 31 | 0.1226 | 0.5112 | between_patient_not_the_within_sample_split |
| EHF(+)/CLDN4out | 31 | 0.3488 | 0.05447 | between_patient_not_the_within_sample_split |
| ELF1(+)/CLDN4out | 31 | -0.1093 | 0.5584 | between_patient_not_the_within_sample_split |
| ELF3(+)/CLDN4out | 31 | 0.6383 | 1.12e-04 | between_patient_not_the_within_sample_split |
| ELK1(+)/CLDN4out | 31 | -0.4335 | 0.01485 | between_patient_not_the_within_sample_split |
| ELK3(+)/CLDN4out | 31 | 0.2625 | 0.1537 | between_patient_not_the_within_sample_split |
| ERF(+)/CLDN4out | 31 | -0.5093 | 0.003433 | between_patient_not_the_within_sample_split |
| ETV7(+)/CLDN4out | 31 | 0.1133 | 0.5439 | between_patient_not_the_within_sample_split |
| FOS(+)/CLDN4out | 31 | 0.2532 | 0.1693 | between_patient_not_the_within_sample_split |
| FOSL2(+)/CLDN4out | 31 | 0.2214 | 0.2314 | between_patient_not_the_within_sample_split |
| IKZF2(+)/CLDN4out | 31 | -0.2706 | 0.141 | between_patient_not_the_within_sample_split |
| IRF1(+)/CLDN4out | 31 | 0.179 | 0.3352 | between_patient_not_the_within_sample_split |
| IRF2(+)/CLDN4out | 31 | 0.04839 | 0.796 | between_patient_not_the_within_sample_split |
| IRF3(+)/CLDN4out | 31 | -0.173 | 0.3521 | between_patient_not_the_within_sample_split |
| IRF7(+)/CLDN4out | 31 | 0.123 | 0.5098 | between_patient_not_the_within_sample_split |
| JUNB(+)/CLDN4out | 31 | 0.2202 | 0.234 | between_patient_not_the_within_sample_split |
| KLF4(+)/CLDN4out | 31 | 0.3246 | 0.07481 | between_patient_not_the_within_sample_split |
| KLF5(+)/CLDN4out | 31 | -0.1032 | 0.5805 | between_patient_not_the_within_sample_split |
| KLF6(+)/CLDN4out | 31 | 0.4101 | 0.02195 | between_patient_not_the_within_sample_split |
| NFE2L3(+)/CLDN4out | 31 | 0.1391 | 0.4555 | between_patient_not_the_within_sample_split |
| NFIC(+)/CLDN4out | 31 | 0.506 | 0.003679 | between_patient_not_the_within_sample_split |
| NR3C1(+)/CLDN4out | 31 | -0.1669 | 0.3694 | between_patient_not_the_within_sample_split |
| OVOL2(+)/CLDN4out | 31 | 0.3157 | 0.08359 | between_patient_not_the_within_sample_split |
| PHF20(+)/CLDN4out | 31 | -0.1823 | 0.3264 | between_patient_not_the_within_sample_split |
| PROGRAM::CTRL_RIBO | 31 | 0.08065 | 0.6663 | between_patient_not_the_within_sample_split |
| PROGRAM::HALLMARK_INTERFERON_ALPHA_RESPONSE | 31 | 0.1431 | 0.4424 | between_patient_not_the_within_sample_split |
| PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE | 31 | 0.04234 | 0.8211 | between_patient_not_the_within_sample_split |
| PROGRAM::IFN_UNION | 31 | 0.08589 | 0.6459 | between_patient_not_the_within_sample_split |
| PROGRAM::KERATIN | 31 | 0.2843 | 0.1212 | between_patient_not_the_within_sample_split |
| PROGRAM::MHC1_APM | 31 | 0.1766 | 0.3419 | between_patient_not_the_within_sample_split |
| PROGRAM::TJ_STRUCT | 31 | 0.6323 | 1.36e-04 | between_patient_not_the_within_sample_split |
| RFX5(+)/CLDN4out | 31 | 0.2073 | 0.2632 | between_patient_not_the_within_sample_split |
| RFXANK(+)/CLDN4out | 31 | -0.3964 | 0.02728 | between_patient_not_the_within_sample_split |
| SMAD1(+)/CLDN4out | 31 | 0.1419 | 0.4463 | between_patient_not_the_within_sample_split |
| SPDEF(+)/CLDN4out | 31 | 0.2718 | 0.1391 | between_patient_not_the_within_sample_split |
| STAT1(+)/CLDN4out | 31 | 0.02097 | 0.9109 | between_patient_not_the_within_sample_split |
| STAT2(+)/CLDN4out | 31 | -0.04919 | 0.7927 | between_patient_not_the_within_sample_split |
| TFAP2A(+)/CLDN4out | 31 | -0.3444 | 0.05783 | between_patient_not_the_within_sample_split |
| TFAP2C(+)/CLDN4out | 31 | -0.2302 | 0.2127 | between_patient_not_the_within_sample_split |

## What this is / is not

- **Is** a pySCENIC run (GRNBoost2 + cisTarget + AUCell) on GSE131907 author-malignant cells, with the CLDN4-high versus low contrast at patient level.
- **Is** an explicit report of observational IFN sign against the KD-like expectation and against GSE207704 and GSE50927, including opposite signs.
- **Is not** SCENIC+. There is no scATAC in this accession.
- **Is not** a concordant-4 integrated GRN. Memory stopped that merge. The other three cohorts were not scored here.
- **Is not** a knockdown, a coculture, or an estimate of the private CLDN4 KD. Those data are not in this public repository and were not re-fit.
- **Is not** ChIP, and cisTarget motif recovery is not proof that the TF binds the target in these tumors.
- **Is not** an ICI or MPR result. GSE131907 is treatment-naive.
- **Is not** a new ELF3–CLDN4 discovery.
- **Is not** a cell-level p-value. EBUS and brain-metastasis samples can dominate cell counts; the patient is the unit.
- numpy shim: `np.object = object` before importing pySCENIC 0.12.1, because NumPy 2 removed that alias. No GRN equation was changed.

## Figures

- `figures/fig_ifn_sign_board.png` — primary IFN board only (STAT/IRF cisTarget + Hallmark program AUCell), colored by SAME versus OPPOSITE the KD-like expectation
- `figures/fig_hallmark_overlap_not_ifn_tf.png` — regulons with Hallmark-IFN target overlap; not IFN transcription factors; signs kept
- `figures/fig_cistarget_forest.png` — cisTarget regulons, CLDN4 held out; red = STAT/IRF, gold = MHC TF, green = TJ-associated TF
- `figures/fig_patient_paired_ifn_tj.png` — per-patient high versus low for the IFN union and the TJ program
- `figures/fig_honest_n.png` — cell, sample, and patient counts

## Reproduce

```bash
python3 methods/gse131907_pyscenic_cldn4/scripts/01_prepare.py
python3 methods/gse131907_pyscenic_cldn4/scripts/02_pyscenic.py
python3 methods/gse131907_pyscenic_cldn4/scripts/03_contrast.py
```

Downloads (not committed): GEO UMI + annotation + series matrix; Aerts lab `allTFs_hg38.txt`, `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`, and the two hg38 v10 genes-vs-motifs ranking feathers; MSigDB Hallmark GMT `h.all.v2023.2.Hs.symbols.gmt`.

Run flags: cistarget_run=True; scenicplus_run=False; pyscenic=0.12.1.

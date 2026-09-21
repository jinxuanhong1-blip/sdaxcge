# CLDN4 virtual knockout in concordant-4 malignant cells

Additive public result. The private CLDN4 knockdown still pins CLDN4. This run does not replace it, and it does not use the private KL matrices.

Question asked of the public data: after a simulated CLDN4 knockout in malignant epithelial cells, do IFN, MHC-I, and chemokine programs go up, and do barrier / tight-junction ligands go down?

**Barrier ligands go down. IFN and MHC-I go down as well. The chemokine shift is about zero.**

The barrier-ligand decrease is larger than a housekeeping-gene knockout fit in the same model. The IFN and MHC-I decreases are not: a detection-matched housekeeping gene moves those two programs by a similar amount. Scores are mean log1p(CP10k) on the genes in the program. Full precision is in `results/tables/program_summary.tsv`.

## What ran

Stock CellOracle (0.20) was not executed. Its base GRN is motif links among transcription factors. CLDN4 is absent from the 1,339-TF list in `data/human_tfs.txt`, so a motif GRN gives CLDN4 no targets and a package knockout delta of zero. That cannot answer the question.

The analysis that ran is the same one-step shift CellOracle uses once a regulator has a coefficient. Inside each patient or sample, a ridge regression predicts program-gene expression from expressed TFs plus CLDN4. The knockout sets CLDN4 log1p(CP10k) to 0 and moves each target by `beta_CLDN4 * (0 - mean_CLDN4) / sd_CLDN4`. Ten bootstrap refits are averaged for the reported delta. A 40-permutation shuffle of the CLDN4 column, within the same cells, is the null for that coefficient. The smallest permutation p this design can return is 1/41 = 0.0244.

The GRN is fit separately in each unit, so the patient or sample is the biological replicate. Units with at least 80 QC-pass malignant cells are capped at 400 cells for the fit (seed fixed). CLDN4 has to be detected in at least 5% of those cells. Held-out R² of the whole ridge model is low (median 0.010). The penalty chosen by 3-fold CV is 1000 in 59/61 fits and 100 in 2/61. Sign consistency across units is the result. The absolute shifts are small.

## Honest n

Locked concordant-4 units: **65** (GSE123902 donors 13, GSE131907 tumor-bearing samples 21, GSE205335 patients 22, GSE189357 patients 9).

| dataset | malignant definition | QC malignant cells | units fit |
|---|---|---:|---:|
| GSE123902 | marker epithelium, PTPRC = 0; normal tissue dropped | 4,208 | 12 |
| GSE131907 | author subtype "Malignant cells" | 24,783 | 19 |
| GSE205335 | author "Malignant cells"; normal tissue dropped | 28,512 | 21 |
| GSE189357 | marker epithelium, PTPRC = 0 | 13,789 | 9 |
| fit units |  | 70,771 | **61** |

The GRN used 20,254 of those 70,771 cells after the per-unit cap of 400.

Not fit, and not imputed:

| unit | QC malignant cells | reason |
|---|---:|---|
| GSE123902 LX699 | 39 | below 80 |
| GSE131907 NS_16 | 79 | below 80 |
| GSE205335 P4001 | 27 | below 80 |
| GSE131907 EBUS_13 | 376 | CLDN4 detected in under 5% |

`unit_inventory.tsv` also lists 37 other GSE131907 samples (normal lung, normal lymph node, and tumor samples without the author malignant label). They are outside the locked 21 and were not modeled. GSE131907 primary `tLung` epithelial cells are not in the author "Malignant cells" label; this run keeps that locked definition. GSE205335 CLDN4 detection in malignant cells matches the locked patient table (for example P0031 59.2%, P1006 74.2%). QC removed no GSE205335 malignant cells; that matrix is already filtered.

One MHC-I score (GSE205335 P1056) and two chemokine scores (GSE131907 NS_12, NS_13) are left out of the patient tests because the bootstrap delta and the single-fit delta disagreed in sign. They remain in `patient_program_deltas.tsv` with `stable_sign = false`.

## Primary before / after

Positive delta means the program score rose when CLDN4 was set to 0. The KD-like direction asked for here is up for IFN, MHC-I, and chemokine, and down for barrier ligands. Wilcoxon p-values are paired across units and are descriptive of this simulation.

| program | n | median before | median after | median delta | units in the KD-like direction | Wilcoxon p | median housekeeping delta |
|---|---:|---:|---:|---:|---:|---:|---:|
| IFN | 61 | 0.2803 | 0.2786 | −0.00322 | 7/61 up | 1.96×10⁻⁹ | −0.00454 |
| MHC-I / APM | 60 | 0.8006 | 0.7890 | −0.00893 | 12/60 up | 5.32×10⁻⁸ | −0.0111 |
| chemokine | 59 | 0.0984 | 0.1007 | +0.00022 | 32/59 up | 0.264 | −0.00138 |
| barrier ligands | 61 | 0.3371 | 0.3276 | −0.0174 | 60/61 down | 1.17×10⁻¹¹ | −0.00509 |

Median of (delta / that unit's before score): IFN −0.012, MHC-I −0.012, chemokine +0.004, barrier ligands −0.056.

The barrier unit that did not go down is GSE205335 P1037 (delta +0.00081).

Cohort median barrier-ligand deltas, all down: GSE123902 −0.0112 (12/12), GSE131907 −0.0202 (19/19), GSE205335 −0.0263 (20/21), GSE189357 −0.0140 (9/9). Dropping GSE205335 histologies other than ADC leaves the same signs (`block = ALL_drop_GSE205335_nonADC` in the summary table).

## Housekeeping control

Each unit's model also contains one non-program gene chosen to match CLDN4's detection rate. Across the 61 fits that gene was UBB (24), GSTP1 (19), GAPDH (7), RPLP0 (5), TPT1 (2), FTL (2), or TMSB4X (2). It is knocked out with the same formula.

Paired test of CLDN4 delta minus control delta:

| program | what is specific to CLDN4 | Wilcoxon p |
|---|---|---:|
| barrier ligands | CLDN4 drop is larger (median −0.0174 vs −0.0051) | 3.31×10⁻⁶ |
| IFN | CLDN4 drop is not larger than the control | 0.027 |
| MHC-I / APM | no separation | 0.541 |
| chemokine | no separation | 0.526 |

GSE123902 is the cohort where the barrier result is not specific (CLDN4 median −0.0112, control −0.0139, p = 0.38). The separation comes from GSE131907, GSE205335, and GSE189357.

A broader tight-junction gene set (KEGG ∪ GOBP organization, plus CDH1/VIM/ZEB1, CLDN4 held out) also falls in all 61 units (median delta −0.00440), but the housekeeping knockout does the same (median −0.00414, paired p = 0.680). The CLDN4-linked decrease is the short barrier-ligand list, not the whole cytoskeletal junction set.

Partial coefficients (median across units, log1p per SD of CLDN4; positive means higher CLDN4, higher expression) line up with that split: EPCAM +0.042, CLDN7 +0.033, CDH1 +0.019, F11R +0.017, OCLN +0.015, and also HLA-A +0.020, STAT1 +0.0033. CXCL9 and CXCL10 are about zero (−0.00008). Setting CLDN4 to 0 therefore pulls the positive-coefficient genes down.

## Method notes

- Expression is log1p of counts per 10,000 using the full-transcriptome UMI total, not the gene-universe total.
- Regulators are the top variance TFs detected in ≥10% of the unit's fit cells (cap 15–40 depending on cell count), plus CLDN4, plus the control gene.
- Targets are the program genes with nonzero variance. CLDN4 is held out of every score.
- Programs: Hallmark IFNα ∪ IFNγ (224 genes); the locked MHC-I/APM list (21); the prior chemokine list (26); barrier ligands F11R, NECTIN2, CDH1, LGALS9, CLDN1, CLDN3, CLDN7, OCLN, TJP1, TJP2, TJP3, EPCAM, CGN, MARVELD2, JAM2, JAM3. A program is scored only when at least 5 of its genes enter the model.
- QC: ≥200 genes, ≥500 UMIs, mitochondrial fraction <20%.
- GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 were not added.

## What this is

A within-malignant-cell, within-unit partial association. It is not a knockdown, not a motif GRN, and not the between-patient CLDN4-high versus T/NK contrast already reported for these four cohorts. Low out-of-sample R² means the ridge model is a weak predictor of the transcriptome; the reproducible piece is the sign of the CLDN4 coefficient on barrier ligands.

Figures: `results/figures/before_after_programs.png`, `delta_by_cohort.png`, `cldn4_vs_control.png`.

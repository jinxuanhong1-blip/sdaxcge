# FINDING — GSE131907 SCENIC/GRN proxy, malignant CLDN4-high vs low

Additive public scRNA GRN on **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277).
**A10 ELF3–CLDN4 bulk RNA is taken as given** and is not re-tested as a discovery.

This is a **lightweight SCENIC proxy** (Pearson + public priors + AUCell). Full pySCENIC cisTarget was **not** run. No ChIP peaks are invented. Pearson neighborhoods are co-expression, not binding.

**CLDN4-high regulons only.** CLDN4-low edges are used only to mark high-specific targets. CLDN4 is held out of every regulon (it is the split gene).

Primary table: [`results/gse131907_scenic_cldn4/regulons.tsv`](../../results/gse131907_scenic_cldn4/regulons.tsv).

## Honest n

Primary unit = **sample**. Cell n is labeled as cells and is not the claim (pseudoreplication). EBUS_28 and NS_07 can dominate pooled cell counts.

| item | n | note |
| --- | ---: | --- |
| cells_in_matrix | 208506 | GSE131907 public UMI header |
| author_malignant | 31136 | Cell_subtype in {Malignant cells, tS1, tS2, tS3} |
| malignant_UMI>=200 | 31136 |  |
| samples_ge20_malignant | 31 | min 20 author-malignant after UMI QC; primary GRN samples |
| CLDN4_high_within_sample | 15984 | >= sample median CLDN4 log1p among author-malignant |
| CLDN4_low_within_sample | 15147 | complement in the same samples |
| samples_paired_ge20_high_and_low | 27 | unit of prior-AUCell high vs low |
| tLung_samples_ge20 | 10 | primary-tumor sensitivity |
| tLung_CLDN4_high | 3174 | cells |
| TFs_screened_in_high | 823 | human TF list ∩ det>=0.05 in CLDN4-high |

Author-labeled malignant cells = **31136**. Primary GRN uses samples with ≥20 such cells (**31 samples**): CLDN4-high **15984 cells**, CLDN4-low **15147 cells**. Paired prior-AUCell unit = **27 samples** with ≥20 high and low.

Unlabeled tumor epithelial cells (empty `Cell_subtype`) are **not** counted as malignant. PE / nLung / nLN contribute **0** author-malignant cells. LUNG_T09 has 5 malignant cells and is excluded from the primary sample set.

Four primary samples have **n_low = 0** after the ≥-median split (all malignant cells at/above that sample’s CLDN4 median): EBUS_13, EBUS_15, EBUS_49, NS_16. They stay in the GRN cell pool and are excluded from the paired AUCell unit.

Site mix of the 31,136 author-malignant cells: mBrain 15,423 / tL/B 6,400 / tLung 6,352 (of which 6,347 in the 10 tLung samples with ≥20) / mLN 2,961. EBUS_28 (4,640) and NS_07 (5,108) dominate pooled cell counts.

## CLDN4-high regulon table (focus TFs)

| regulon | tf | kind | n_targets | n_high_specific | mean_r | elf3_given |
| --- | --- | --- | --- | --- | --- | --- |
| ELF3_prior | ELF3 | public_prior | 33 |  |  | True |
| GRHL1_prior | GRHL1 | public_prior | 2 |  |  | False |
| GRHL2_prior | GRHL2 | public_prior | 10 |  |  | False |
| KLF4_prior | KLF4 | public_prior | 131 |  |  | False |
| KLF5_prior | KLF5 | public_prior | 79 |  |  | False |
| OVOL1_prior | OVOL1 | public_prior | 9 |  |  | False |
| OVOL2_prior | OVOL2 | public_prior | 1 |  |  | False |
| TFAP2A_prior | TFAP2A | public_prior | 317 |  |  | False |
| SP1_prior | SP1 | public_prior | 1316 |  |  | False |
| TP63_prior | TP63 | public_prior | 204 |  |  | False |
| ELF3_pearson_cldn4_high | ELF3 | pearson_cldn4_high | 50 | 11 | 0.293 | True |
| ELF3_high_specific | ELF3 | high_specific | 0 | 11 |  | True |
| ELF3_intersect | ELF3 | prior_AND_pearson_high | 1 |  |  | True |
| GRHL1_pearson_cldn4_high | GRHL1 | pearson_cldn4_high | 50 | 4 | 0.145 | False |
| GRHL1_high_specific | GRHL1 | high_specific | 4 | 4 |  | False |
| GRHL1_intersect | GRHL1 | prior_AND_pearson_high | 1 |  |  | False |
| GRHL2_pearson_cldn4_high | GRHL2 | pearson_cldn4_high | 50 | 0 | 0.076 | False |
| GRHL2_high_specific | GRHL2 | high_specific | 0 | 0 |  | False |
| GRHL2_intersect | GRHL2 | prior_AND_pearson_high | 0 |  |  | False |
| KLF4_pearson_cldn4_high | KLF4 | pearson_cldn4_high | 50 | 3 | 0.278 | False |
| KLF4_high_specific | KLF4 | high_specific | 0 | 3 |  | False |
| KLF4_intersect | KLF4 | prior_AND_pearson_high | 3 |  |  | False |
| KLF5_pearson_cldn4_high | KLF5 | pearson_cldn4_high | 50 | 0 | 0.185 | False |
| KLF5_high_specific | KLF5 | high_specific | 0 | 0 |  | False |
| KLF5_intersect | KLF5 | prior_AND_pearson_high | 3 |  |  | False |
| OVOL1_pearson_cldn4_high | OVOL1 | pearson_cldn4_high | 50 | 3 | 0.127 | False |
| OVOL1_high_specific | OVOL1 | high_specific | 3 | 3 |  | False |
| OVOL1_intersect | OVOL1 | prior_AND_pearson_high | 2 |  |  | False |
| OVOL2_pearson_cldn4_high | OVOL2 | pearson_cldn4_high | 27 | 0 | 0.119 | False |
| OVOL2_high_specific | OVOL2 | high_specific | 0 | 0 |  | False |
| OVOL2_intersect | OVOL2 | prior_AND_pearson_high | 0 |  |  | False |
| TFAP2A_pearson_cldn4_high | TFAP2A | pearson_cldn4_high | 50 | 2 | 0.143 | False |
| TFAP2A_high_specific | TFAP2A | high_specific | 2 | 2 |  | False |
| TFAP2A_intersect | TFAP2A | prior_AND_pearson_high | 5 |  |  | False |
| SP1_pearson_cldn4_high | SP1 | pearson_cldn4_high | 50 | 0 | 0.090 | False |
| SP1_high_specific | SP1 | high_specific | 0 | 0 |  | False |
| SP1_intersect | SP1 | prior_AND_pearson_high | 10 |  |  | False |
| TP63_pearson_cldn4_high | TP63 | pearson_cldn4_high | 0 | 0 |  | False |
| TP63_high_specific | TP63 | high_specific | 0 | 0 |  | False |
| TP63_intersect | TP63 | prior_AND_pearson_high | 0 |  |  | False |

Target lists (semicolon-separated) are in `regulons.tsv`. `elf3_given=True` means ELF3 is A10-given, not a discovery from this folder.

`n_high_specific` on a Pearson row is the count of high-specific edges for that TF (not necessarily inside the top-50). The `_high_specific` row’s `n_targets` is the overlap of those edges with the Pearson top list. TP63 Pearson is empty (0 targets above the r floor) — expected in this LUAD atlas; TP63 is not used as a LUSC claim.

## Sample-paired public-prior AUCell (CLDN4-high − low)

| regulon | kind | n_samples | delta_median_high_minus_low | p | fdr |
| --- | --- | --- | --- | --- | --- |
| ELF3_prior | public_prior_AUCell | 27 | 0.00375 | 1.04e-07 | 3.87e-07 |
| GRHL1_prior | public_prior_AUCell | 27 | 0 | 0.249 | 0.27 |
| GRHL2_prior | public_prior_AUCell | 27 | 0.004831 | 1.64e-06 | 4.26e-06 |
| KLF4_prior | public_prior_AUCell | 27 | 0.002964 | 1.49e-08 | 1.29e-07 |
| KLF5_prior | public_prior_AUCell | 27 | 0.002222 | 2.98e-08 | 1.55e-07 |
| OVOL1_prior | public_prior_AUCell | 27 | 0.001916 | 3.52e-05 | 6.54e-05 |
| OVOL2_prior | public_prior_AUCell | 27 | 0 | 0.0277 | 0.0343 |
| TFAP2A_prior | public_prior_AUCell | 27 | 0.0008087 | 7.99e-06 | 1.73e-05 |
| SP1_prior | public_prior_AUCell | 27 | 0.001208 | 1.49e-08 | 1.29e-07 |
| TP63_prior | public_prior_AUCell | 27 | 0.0006633 | 0.00355 | 0.00543 |
| NKX2-1_prior | public_prior_AUCell | 27 | 0.0009015 | 0.934 | 0.934 |
| SOX2_prior | public_prior_AUCell | 27 | 0.0001434 | 0.0772 | 0.0873 |
| FOXA1_prior | public_prior_AUCell | 27 | 0.0004197 | 3.77e-06 | 8.91e-06 |
| ELF3_RNA | tf_rna | 27 | 0.4834 | 1.49e-08 | 1.29e-07 |
| GRHL1_RNA | tf_rna | 27 | 0.04396 | 7.45e-08 | 3.23e-07 |
| GRHL2_RNA | tf_rna | 27 | 0.01645 | 2.59e-05 | 5.17e-05 |
| KLF4_RNA | tf_rna | 27 | 0.1265 | 2.98e-08 | 1.55e-07 |
| KLF5_RNA | tf_rna | 27 | 0.1344 | 2.83e-07 | 9.2e-07 |
| OVOL1_RNA | tf_rna | 27 | 0.02732 | 1.64e-06 | 4.26e-06 |
| OVOL2_RNA | tf_rna | 27 | 0.01142 | 0.0151 | 0.0196 |
| TFAP2A_RNA | tf_rna | 27 | 0.02333 | 0.00113 | 0.00184 |
| SP1_RNA | tf_rna | 27 | 0.01748 | 0.00011 | 0.00019 |
| TP63_RNA | tf_rna | 27 | 0 | 0.831 | 0.865 |
| NKX2-1_RNA | tf_rna | 27 | 0.04594 | 0.00502 | 0.00687 |
| SOX2_RNA | tf_rna | 27 | 0.005553 | 0.052 | 0.0615 |
| FOXA1_RNA | tf_rna | 27 | 0.0254 | 0.00388 | 0.0056 |

Prior-AUCell **deltas are small** (order 10⁻³) even where Wilcoxon p is low at n=27. Those p-values are not a binding claim and are not treated as a SUPPORT call. ELF3 RNA Δ ≈ 0.48 is **A10-given**, not a discovery from this folder.

## What is / is not claimed

- **Claimed:** a descriptive CLDN4-high regulon table on author-malignant GSE131907 cells, with honest sample n.
- **Not claimed:** TF binding, pySCENIC cisTarget, or a new ELF3–CLDN4 discovery (A10 is given).
- **Not claimed:** ICI / MPR association (this atlas is treatment-naive).
- **Not claimed:** cell-level p-values as the finding.
- **Not claimed:** that small prior-AUCell deltas at n=27 constitute a GRN mechanism.

Method flags: pyscenic importable=False; cisTarget_run=False; chip_peaks_invented=False.

# CLDN4 KD match wave: signed IFN/APM concordance

Thesis sign, taken from the public C4 test: **CLDN4 loss opens IFN and MHC-I/APM**. The signed effect is the median log2 fold change of the prespecified panel (KD/KO versus matched control). Positive matches the thesis.

This wave combines only the true CLDN4-loss transcriptomes. TACSTD2/TROP2 knockdowns are not included. GSE22493 stays in because GPL10555 probes map to gene symbols (IFN/ISG 27/39, APM 14/18; both above the 50% gate). GSE50927 contributes the uninjured baseline contrast only. VILI is not in the forest.

Numbers are the locked Claim C4 effects (PR #88), recomputed here only to attach gene-resample intervals and the cancer-versus-non-cancer labels. Medians, means, and up-counts are asserted against that table.

## Concordance

| Class | Organ | Contrast | CLDN4 log2FC | IFN/ISG | MHC-I/APM |
|---|---|---|---:|---|---|
| Non-cancer | Lung | GSE50927 naive KO vs WT | -6.06 | +0.44 (34/37; concordant) | +0.28 (14/16; concordant) |
| Cancer | Breast | GSE207704 T47D KO | -1.06 | -0.35 (7/21; discordant) | -0.22 (2/8; discordant) |
| Cancer | Breast | GSE207704 MCF-7 KO | -0.72 | -0.28 (6/17; unsigned) | -0.10 (2/8; unsigned) |
| Cancer | Ovary | GSE22493 SKOV-3 KD | -1.23 | -0.01 (12/27; unsigned) | -0.45 (4/14; discordant) |

A locked call is **concordant** only when Claim C4 called the panel `up` (median above 0 and a one-sided Wilcoxon or sign test P ≤ 0.05). **Discordant** means the locked call was `down`. **Unsigned** means the median did not clear that bar (`weak_down` here).

Non-cancer lung is concordant on both panels: IFN/ISG median +0.44 (34/37 up) and MHC-I/APM median +0.28 (14/16 up). Classical H2-K1, the mouse stand-in for HLA-A, is -0.07 in that same baseline table, so the APM concordance is the panel median, not MHC-I heavy-chain induction.

Every cancer median is negative, in breast and in ovary, on both panels. IFN/ISG locked calls: 0 concordant, 1 discordant, 2 unsigned. MHC-I/APM locked calls: 0 concordant, 2 discordant, 1 unsigned. No cancer contrast matches the thesis call. The median of the three cancer medians is -0.28 for IFN/ISG (range -0.35 to -0.01) and -0.22 for MHC-I/APM (range -0.45 to -0.10).

T47D clears a down call on both panels. Its IFN/ISG gene-resample interval still touches 0 (-0.62 to +0.02). That call is the locked Wilcoxon on the gene list (7 up, 14 down), which answers a different question than the median's resample interval. SKOV-3 clears a down call for APM and stays unsigned for IFN/ISG (median essentially zero, slightly negative). MCF-7 is unsigned on both, with a negative median and a thin gene list.

## Cross-organ gene signs

Genes counted below are observed in the lung baseline and in at least one cancer contrast. The cancer sign is the median of the cancer log2 fold changes available for that gene.

| Module | Genes in both | Same sign | Lung up, cancer down | Both up |
|---|---:|---:|---:|---:|
| IFN / ISG | 31 | 11 | 19 | 9 |
| MHC-I / APM | 15 | 4 | 11 | 2 |

Same-sign agreement across organ class is the minority. The common pattern among shared genes is lung up and cancer down.

## What the forest is

Thin bars are 2.5–97.5% gene-resample intervals of the locked panel median (10000 draws, seed 50927). They describe whether the genes in that panel move together. They are not confidence intervals from biological replicates. GSE207704 deposits pooled FPKM. GSE50927 baseline is one naive KO library versus one naive WT library. GSE22493 has three arrays, but the interval still resamples genes, so it matches the other rows.

The diamond is the median of the three cancer medians. Its width is the min-to-max of those medians. Cancer and non-cancer stay on separate rows.

A DerSimonian–Laird summary that weights each cancer median by its gene-resample SD is in `tables/stratum_summary.tsv`. IFN/ISG: -0.19 (Wald -0.43 to +0.04). MHC-I/APM: -0.23 (Wald -0.39 to -0.07). I² is 0 on both because those SDs are wide relative to the gap between medians. The weights are gene-resample SDs, so that summary stays in the table and off the figure.

## Files

- `tables/concordance_ifn_apm.tsv` — one row per contrast and panel
- `tables/gene_signed.tsv` — gene-level log2FC by contrast
- `tables/cross_organ_genes.tsv` — genes observed on both sides of the class split
- `tables/cross_organ_summary.tsv`
- `tables/stratum_summary.tsv`
- `tables/forest_effects.tsv`
- `figures/forest_ifn_apm.png`

Reproduce: `python3 scripts/cldn4_kd_match_wave/build_concordance.py`

Cancer-cell discordance is re-checked in `SWEEP.md`. Alternate collapses, pathway tests, outlier and probe rules, and a kallisto n=2 versus n=2 quantification of both GSE207704 lines (CLDN4 down in both; pseudoalignment 37.0–42.1%) do not call IFN/ISG or MHC-I/APM up on the combined cancer contrasts. The only non-optimistic up calls are one SKOV-3 array (GSM558701) for Hallmark IFN-α, Hallmark IFN-γ, NHEJ_EXT, and Reactome NHEJ. NHEJ_CORE stays weak in both classes. That sweep does not refit the lung IFN/APM baseline.

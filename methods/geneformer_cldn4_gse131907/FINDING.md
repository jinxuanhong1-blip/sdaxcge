# FINDING — Geneformer / scGPT embeddings of CLDN4-high vs low malignant cells

ADDITIVE. **CLDN4-only.** GSE131907 author-malignant cells only.
This does **not** re-estimate the concordant-4 patient Spearman
(n=65, ρ=−0.53). It does not use Visium, and it does not merge any
mouse matrix. Sample is the unit (Kim et al. sample id, not patient).
Cells are a **stratified** draw so both arms exist: up to 50
CLDN4 UMI>0 and 50 CLDN4 UMI=0 per sample after QC. They are
not a prevalence-weighted sample of the malignant compartment.

Primary embedding test: Geneformer **V1-10M** (the 10M model; V2-104M and
V2-316M are the heavier checkpoints), **CLDN4 token removed** before
ranking, leave-one-sample-out projection onto the other samples'
high-minus-low direction. Positive median Δ means the held-out sample's
CLDN4-high cells sit on the high side of a direction they did not help fit.

## Honest n

- Author-malignant samples in tumor-bearing origins with ≥20 malignant cells: **21**.
- Samples in the paired test (≥5 embedded cells in each arm): **20**.
- Embedded cells: **1874** (959 high / 915 low). This is not the test n.
- tLung epithelial states labeled tS1/tS2 are not `Cell_subtype == Malignant cells` and are not in this object.
- Full-sample CLDN4 %pos (all author-malignant cells, not the cap) is in `results/tables/sample_inventory.tsv`.

## Library size comes first

CLDN4 UMI>0 vs UMI=0 inside a sample is not a depth-balanced contrast. Unmatched embedded cells: median nUMI high 17189 vs low 2103; median genes detected 3902 vs 873 (20 samples).

| measure | n_samples | median_delta_high_minus_low | n_samples_delta_pos | p_wilcoxon |
|---|---|---|---|---|
| log1p_n_umi | 20 | 1.78 | 20 | 1.91e-06 |
| n_genes | 20 | 2.76e+03 | 20 | 1.91e-06 |

The unmatched program table below moves together with that gap: large gene sets rise in the deeper arm. That is not evidence for a specific program. The caliper match is a sensitivity added after this gap was measured, and it is the contrast used for program claims.

Caliper match |Δ log1p(nUMI)| ≤ 0.25, greedy within sample, ≥5 pairs: **297 pairs**, **19 samples**, median |Δ log1p(nUMI)| = 0.0401.

## Foundation-model separation (unmatched draw)

Primary result, Geneformer V1-10M holdout LOSO: median Δ = 0.335, 20/20 samples positive, Wilcoxon p = 1.91e-06.

| embedding | n_samples | median Δ | pos/n | p | q (holdout family) |
|---|---:|---:|---:|---:|---:|
| gf_v1_holdout_loso | 20 | 0.335 | 20/20 | 1.91e-06 | 1.91e-06 |
| gf_v1_full_loso | 20 | 0.336 | 20/20 | 1.91e-06 | NA |
| gf_v1_holdout_loso_depth_resid | 20 | 0.00715 | 18/20 | 2.67e-05 | NA |
| scgpt_holdout_loso | 20 | 0.33 | 20/20 | 1.91e-06 | 1.91e-06 |
| scgpt_full_loso | 20 | 0.331 | 20/20 | 1.91e-06 | NA |
| gf_v2_cancer_holdout_loso | NA | NA | NA | NA | NA |

`gf_v1_holdout_loso_depth_resid` residualizes the projection on log1p(nUMI) inside each sample. q is BH across the holdout tests that were fit (V1 and scGPT when both ran). Full-token and V2 rows are sensitivities and are not in that q.

CLDN4 itself (positive control) median paired Δ log1p(CP10k) = 1.65, p=1.91e-06, n=20.

Figure: `results/figures/fig_geneformer_v1_loso_paired.png`.

## Depth-matched contrast

Same embeddings, restricted to the caliper pairs. PCA and the diffusion map are refit on those cells with CLDN4 still held out of the features. Geneformer / scGPT vectors are not refit; only the leave-one-sample-out contrast is recomputed on the matched cells.

| embedding | n_samples | n_pairs | median_delta_high_minus_low | n_samples_delta_pos | p_wilcoxon | q_bh_holdouts |
|---|---|---|---|---|---|---|
| gf_v1_holdout | 19 | 297 | 0.0532 | 17 | 1.91e-05 | 3.81e-05 |
| gf_v1_full | 19 | 297 | 0.0528 | 17 | 1.91e-05 | NA |
| scgpt_holdout | 19 | 297 | 0.0951 | 17 | 5.34e-05 | 5.34e-05 |
| scgpt_full | 19 | 297 | 0.0946 | 17 | 5.34e-05 | NA |

This cell-level contrast is not the concordant-4 patient pseudobulk (IFN/MHC lower in CLDN4-high units). A positive cell-level Δ does not replace that result.

Depth-matched programs with BH q<0.10: IFN Δ=0.0242 q=0.00183; MHC-I/APM Δ=0.0752 q=0.00249; TJ Δ=0.0215 q=0.00183; EMT Δ=0.0254 q=0.00183; hypoxia Δ=0.0368 q=7.13e-04; G2M Δ=0.00406 q=0.0763; TNFA Δ=0.0633 q=3.24e-04; inflammatory Δ=0.0257 q=6.16e-04; apical_junction Δ=0.0174 q=0.00192; barrier_keratin Δ=0.0525 q=0.0763; ciliated Δ=0.0208 q=0.0729.

| program | n_samples | median_delta_high_minus_low | n_samples_delta_pos | p_wilcoxon | q_bh |
|---|---|---|---|---|---|
| IFN | 19 | 0.0242 | 15 | 6.45e-04 | 0.00183 |
| MHC-I/APM | 19 | 0.0752 | 15 | 0.00117 | 0.00249 |
| TJ | 19 | 0.0215 | 16 | 5.23e-04 | 0.00183 |
| chemokine | 19 | 0.00735 | 11 | 0.395 | 0.48 |
| EMT | 19 | 0.0254 | 15 | 6.45e-04 | 0.00183 |
| hypoxia | 19 | 0.0368 | 16 | 1.26e-04 | 7.13e-04 |
| G2M | 19 | 0.00406 | 14 | 0.0494 | 0.0763 |
| E2F | 19 | -0.000626 | 8 | 0.953 | 0.953 |
| TNFA | 19 | 0.0633 | 17 | 1.91e-05 | 3.24e-04 |
| inflammatory | 19 | 0.0257 | 17 | 7.25e-05 | 6.16e-04 |
| apical_junction | 19 | 0.0174 | 17 | 7.90e-04 | 0.00192 |
| barrier_keratin | 19 | 0.0525 | 13 | 0.0494 | 0.0763 |
| AT2 | 19 | 0.0234 | 11 | 0.441 | 0.5 |
| AT1 | 19 | 0 | 7 | 0.234 | 0.306 |
| club | 19 | 0 | 5 | 0.177 | 0.251 |
| basal | 19 | 0.00357 | 10 | 0.687 | 0.73 |
| ciliated | 19 | 0.0208 | 13 | 0.0386 | 0.0729 |

Depth-matched programs tracking the holdout axis at BH q<0.10: IFN median ρ=0.469 q=1.13e-04; MHC-I/APM median ρ=0.283 q=2.38e-04; TJ median ρ=0.453 q=1.13e-04; chemokine median ρ=0.281 q=0.0011; EMT median ρ=0.513 q=2.59e-05; hypoxia median ρ=0.459 q=7.57e-05; G2M median ρ=0.429 q=2.59e-05; E2F median ρ=0.287 q=0.00122; TNFA median ρ=0.665 q=2.16e-05; inflammatory median ρ=0.543 q=2.16e-05; apical_junction median ρ=0.544 q=2.16e-05; AT2 median ρ=0.188 q=0.0199; basal median ρ=0.109 q=0.0304.

Depth-matched genes with the highest mean within-sample Spearman vs the Geneformer holdout axis: NEAT1, MALAT1, WSB1, XIST, MUC1, RNF213, HOOK2, POLR2J3, N4BP2L2, MT-ND5, KLF6, MT-ND4.

| component | n_samples | median_delta_high_minus_low | n_samples_delta_pos | p_wilcoxon | q_bh |
|---|---|---|---|---|---|
| PC1 | 19 | 4.16 | 14 | 0.00117 | 0.0039 |
| PC2 | 19 | 3.91 | 17 | 2.67e-05 | 2.67e-04 |
| PC3 | 19 | 0.0283 | 10 | 0.312 | 0.347 |
| PC4 | 19 | 1.98 | 16 | 4.20e-04 | 0.0021 |
| PC5 | 19 | 0.181 | 11 | 0.86 | 0.86 |
| PC6 | 19 | 0.681 | 15 | 0.00823 | 0.0165 |
| PC7 | 19 | 0.895 | 15 | 0.0108 | 0.018 |
| PC8 | 19 | 0.52 | 13 | 0.0728 | 0.091 |
| PC9 | 19 | 0.55 | 13 | 0.0361 | 0.0515 |
| PC10 | 19 | 0.294 | 16 | 0.00392 | 0.00979 |
| DC1 | 19 | 0.000233 | 14 | 0.156 | 0.261 |
| DC2 | 19 | 0.000204 | 10 | 0.86 | 0.953 |
| DC3 | 19 | 0.000199 | 13 | 0.156 | 0.261 |
| DC4 | 19 | 0.000294 | 12 | 0.312 | 0.446 |
| DC5 | 19 | 0.000446 | 11 | 0.651 | 0.813 |
| DC6 | 19 | 0.000577 | 17 | 0.0024 | 0.024 |
| DC7 | 19 | 0.000676 | 14 | 0.0258 | 0.0861 |
| DC8 | 19 | 0.000544 | 15 | 0.0204 | 0.0861 |
| DC9 | 19 | 0.000396 | 12 | 0.0546 | 0.136 |
| DC10 | 19 | 1.4e-05 | 10 | 0.953 | 0.953 |

Figures: `fig_geneformer_v1_loso_paired_matched.png`, `fig_program_paired_matched.png`, `fig_nearest_genes_axis_matched.png`, `fig_pca_diffusion_paired_matched.png`.

## Nearest gene programs (unmatched draw)

Expression programs are mean log1p(CP10k). TJ / apical junction / every multi-gene set has **CLDN4 removed**. IFN is Hallmark IFNα ∪ IFNγ. MHC-I/APM is the repo custom set. Barrier/keratin matches the prior epithelial list and does not contain CLDN4. BH is across the 17 programs, not across the CLDN4 or TACSTD2 controls.

Programs with BH q<0.10: IFN Δ=0.0947 q=2.70e-06; MHC-I/APM Δ=0.227 q=2.70e-06; TJ Δ=0.0852 q=2.70e-06; chemokine Δ=0.0297 q=2.70e-06; EMT Δ=0.0728 q=2.70e-06; hypoxia Δ=0.0949 q=2.70e-06; G2M Δ=0.0757 q=2.70e-06; E2F Δ=0.0769 q=2.70e-06; TNFA Δ=0.149 q=2.70e-06; inflammatory Δ=0.0648 q=2.70e-06; apical_junction Δ=0.0693 q=2.70e-06; barrier_keratin Δ=0.304 q=2.70e-06; AT2 Δ=0.113 q=7.64e-05; AT1 Δ=0.0183 q=0.0167; basal Δ=0.0185 q=0.00724; ciliated Δ=0.0472 q=6.24e-05.

| program | n_genes | n_samples | median Δ | pos/n | p | q |
|---|---:|---:|---:|---:|---:|---:|
| IFN | 220 | 20 | 0.0947 | 20/20 | 1.91e-06 | 2.70e-06 |
| MHC-I/APM | 21 | 20 | 0.227 | 20/20 | 1.91e-06 | 2.70e-06 |
| TJ | 207 | 20 | 0.0852 | 20/20 | 1.91e-06 | 2.70e-06 |
| chemokine | 26 | 20 | 0.0297 | 20/20 | 1.91e-06 | 2.70e-06 |
| EMT | 198 | 20 | 0.0728 | 20/20 | 1.91e-06 | 2.70e-06 |
| hypoxia | 194 | 20 | 0.0949 | 20/20 | 1.91e-06 | 2.70e-06 |
| G2M | 190 | 20 | 0.0757 | 20/20 | 1.91e-06 | 2.70e-06 |
| E2F | 195 | 20 | 0.0769 | 20/20 | 1.91e-06 | 2.70e-06 |
| TNFA | 198 | 20 | 0.149 | 20/20 | 1.91e-06 | 2.70e-06 |
| inflammatory | 199 | 20 | 0.0648 | 20/20 | 1.91e-06 | 2.70e-06 |
| apical_junction | 192 | 20 | 0.0693 | 20/20 | 1.91e-06 | 2.70e-06 |
| barrier_keratin | 6 | 20 | 0.304 | 20/20 | 1.91e-06 | 2.70e-06 |
| AT2 | 6 | 20 | 0.113 | 17/20 | 6.29e-05 | 7.64e-05 |
| AT1 | 3 | 20 | 0.0183 | 16/20 | 0.0158 | 0.0167 |
| club | 3 | 20 | 0.00722 | 12/20 | 0.356 | 0.356 |
| basal | 4 | 20 | 0.0185 | 16/20 | 0.00639 | 0.00724 |
| ciliated | 3 | 20 | 0.0472 | 19/20 | 4.77e-05 | 6.24e-05 |

Along the Geneformer holdout axis (within-sample Spearman of the program score vs the LOSO projection, Wilcoxon across samples): Programs tracking the holdout axis at BH q<0.10: IFN median ρ=0.78 q=2.49e-06; MHC-I/APM median ρ=0.54 q=2.49e-06; TJ median ρ=0.825 q=2.49e-06; chemokine median ρ=0.448 q=2.49e-06; EMT median ρ=0.763 q=2.49e-06; hypoxia median ρ=0.722 q=2.49e-06; G2M median ρ=0.716 q=2.49e-06; E2F median ρ=0.686 q=2.49e-06; TNFA median ρ=0.748 q=2.49e-06; inflammatory median ρ=0.792 q=2.49e-06; apical_junction median ρ=0.814 q=2.49e-06; barrier_keratin median ρ=0.354 q=2.49e-06; AT2 median ρ=0.311 q=4.32e-06; AT1 median ρ=0.363 q=1.22e-05; club median ρ=0.218 q=2.14e-04; basal median ρ=0.345 q=4.32e-06; ciliated median ρ=0.369 q=2.49e-06.

Genes with the highest mean within-sample Spearman vs that axis (detection ≥5%; CLDN4's own ρ = 0.566): TAF10, MTDH, CLDN4, VWA1, TSR3, DSG2, ACBD3, ATP2A2, GNAI2, RBM47, ECI1, CTSZ.

Lowest: RPL31, RPS24, RPL21, RPL37, RPL35, RPL26, RPS3A, RPL39.

Token-space probe (same layer as the cell vector; cosine of the mean contextual token state to the high−low vector; tokens seen in ≥30 cells): `results/tables/nearest_genes_token_cosine.tsv`. The holdout run cannot list CLDN4 because that token was removed. CLDN4 is not among the top 25 full-token cosine genes, so that probe is not used as evidence that the marker token was recovered.

## PCA / diffusion map (CLDN4 held out of the features)

2000 HVGs by variance of log1p(CP10k), CLDN4 excluded, z-scored, 30 PCs (full SVD). Diffusion map: k=15 adaptive Gaussian on those PCs, symmetric normalization, 10 nontrivial components. Components are sign-oriented so the median paired Δ is ≥0; p-values are two-sided and unchanged by that flip. BH is within the 10 PCs and, separately, within the 10 DCs.

PCs with BH q<0.10: PC1 Δ=12.7 q=1.91e-05; PC2 Δ=5.81 q=2.38e-04; PC3 Δ=1.5 q=0.00358; PC4 Δ=2 q=4.45e-04; PC5 Δ=0.377 q=0.0666; PC9 Δ=0.75 q=0.0655. None of the tested diffusion components survive BH q<0.10.

| component | n | median Δ | pos/n | p | q |
|---|---:|---:|---:|---:|---:|
| PC1 | 20 | 12.7 | 20/20 | 1.91e-06 | 1.91e-05 |
| PC2 | 20 | 5.81 | 19/20 | 4.77e-05 | 2.38e-04 |
| PC3 | 20 | 1.5 | 16/20 | 0.00143 | 0.00358 |
| PC4 | 20 | 2 | 16/20 | 1.34e-04 | 4.45e-04 |
| PC5 | 20 | 0.377 | 14/20 | 0.04 | 0.0666 |
| PC6 | 20 | 0.858 | 14/20 | 0.154 | 0.171 |
| PC7 | 20 | 0.0697 | 11/20 | 0.812 | 0.812 |
| PC8 | 20 | 0.448 | 12/20 | 0.154 | 0.171 |
| PC9 | 20 | 0.75 | 15/20 | 0.0328 | 0.0655 |
| PC10 | 20 | 0.514 | 12/20 | 0.0897 | 0.128 |
| DC1 | 20 | 0.000131 | 16/20 | 0.0136 | 0.133 |
| DC2 | 20 | 9.17e-05 | 10/20 | 0.622 | 0.777 |
| DC3 | 20 | 0.000711 | 16/20 | 0.0266 | 0.133 |
| DC4 | 20 | 0.000444 | 13/20 | 0.294 | 0.736 |
| DC5 | 20 | 0.000285 | 12/20 | 0.452 | 0.777 |
| DC6 | 20 | 0.000332 | 15/20 | 0.0696 | 0.232 |
| DC7 | 20 | 0.000271 | 13/20 | 0.596 | 0.777 |
| DC8 | 20 | 0.000324 | 11/20 | 0.784 | 0.871 |
| DC9 | 20 | 0.000114 | 11/20 | 0.956 | 0.956 |
| DC10 | 20 | 0.000239 | 12/20 | 0.571 | 0.777 |

Trivial diffusion eigenvalue (component dropped): 1.

## What was installed and what was skipped

| model | decision |
|---|---|
| Geneformer-V1-10M | Ran. 6 layers, hidden 256, max 2048, no CLS. This is the lite checkpoint relative to V2. |
| Geneformer-V2-104M and V2-104M_CLcancer | 418MB, 12×768, max 4096, CLS+EOS. Cancer weights were benchmarked; the full cohort was run only if the probe was ≤1.5 s/cell. |
| Geneformer-V2-316M | Not loaded (1.26GB). Too heavy for this 16GB CPU host once V2-104M is the larger model under test. |
| scGPT whole-human (`wanglab/scGPT-human`) | Weights loaded. Official flash-attn kernel is not installed; the Wqkv projection runs as eager scaled-dot-product attention, post-norm, ReLU FFN, CLS pooling, L2-normalized. Binning is the official 51-bin within-cell quantile. Sequences longer than 1199 genes are a seeded subset, not a new sample of cells. |

| model | status | sec/cell | note |
|---|---|---:|---|
| Geneformer-V1-10M | self_check_ok | NA | BertForMaskedLM 6x256, 41MB safetensors, CPU |
| scGPT-human | self_check_ok | NA | eager attention in place of flash-attn; weights 196MB |
| Geneformer-V1-10M-holdout | ran | 0.148 | {"n_genes_in_vocab": 17413, "n_empty_cells": 0, "holdout": "ENSG00000189143", "special_tokens": false, "max_len": 2048} |
| Geneformer-V1-10M-full | ran | 0.158 | {"n_genes_in_vocab": 17414, "n_empty_cells": 0, "holdout": null, "special_tokens": false, "max_len": 2048} |
| scGPT-human-holdout | ran | 0.478 | CLS L2-normalized; eager attention; CLDN4 held out |
| scGPT-human-full | ran | 0.585 | CLS L2-normalized; eager attention; CLDN4 token kept |
| Geneformer-V2-104M_CLcancer | skipped_time | 2.7 | benchmark only; full cohort not run because probe > 1.5 s/cell on CPU |
| Geneformer-V2-316M | skipped_weight | NA | 1.26GB root checkpoint not loaded. Heavier than V2-104M on a 16GB CPU host; not run. |

## Reproduce

```bash
python3 -m venv /tmp/venv
/tmp/venv/bin/pip install -r methods/geneformer_cldn4_gse131907/requirements.txt
# torch CPU wheel if the default torch build is unwanted
GF_CACHE=/tmp/gf_cache /tmp/venv/bin/python methods/geneformer_cldn4_gse131907/analyze.py
```

GEO inputs (not committed): `GSE131907_Lung_Cancer_cell_annotation.txt.gz`, `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`.
Model weights (not committed): `ctheodoris/Geneformer` `Geneformer-V1-10M` plus the gc30M dictionaries; optional V2-104M_CLcancer and gc104M dictionaries; `wanglab/scGPT-human`.

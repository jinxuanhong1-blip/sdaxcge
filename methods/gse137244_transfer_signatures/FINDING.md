# FINDING — GSE137244 transferable KL vs KP signatures

Public cell-line RNA-seq only. Deng et al., *Nature Cancer* 2021 ([PMID 34142094](https://pubmed.ncbi.nlm.nih.gov/34142094/)), [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244). Matrices re-downloaded from the GEO supplementary files. Normal lung is held out. No private 8-KL matrices. No TISMO LLC.

Primary scale is the mean of log2(FPKM+1), KL minus KP. That is the scale of the locked Tacstd2 and Cldn4 deltas. The rank test is the two-sided exact Mann–Whitney test. For 5 vs 5 with complete separation the p-value is 2/252 = 0.00794, and it cannot go lower.

The deposited FPKM and normalized files store March- and Sept-family symbols as Excel dates (`1-Mar`, `1-Sep`). Row order matches the raw-count file, which still has `March1` and `Sept1`. Those 26 symbols were restored from the raw-count row order before any test. None of them is in the pre-specified TJ, NHEJ, STING, IFN, or APM lists.

## Samples

KP (KrasG12D; Trp53, GEO description) is B6AL10-1 through B6AL10-5. KL (KrasG12D; Lkb1) is KL155mix-control-2, KL47-1-untreated-1, KLC, KLD, and KLE. Titles that say control or untreated are the baseline libraries in this series. Honest n for the primary contrast is 5 libraries vs 5 libraries. Several libraries inside an arm are near-replicates; the maximum within-arm Pearson r of log2(FPKM+1) on genes with mean FPKM ≥ 1 is **0.990**.

Genotype check on the same scale: Stk11 Δ = **-2.523** (exact p 0.00794). Trp53 Δ = **+2.433** (exact p 0.00794). Every KL library is below every KP library for Stk11, and every KL library is above every KP library for Trp53.

## Locked genes

| gene | Δ log2(FPKM+1) | exact MW p | separation |
|---|---:|---:|---|
| Tacstd2 | **+3.238** | 0.00794 | KL>KP |
| Cldn4 | **+5.570** | 0.00794 | KL>KP |

Handoff rounding is Tacstd2 +3.24 and Cldn4 +5.57. Both p-values are the 5-vs-5 floor.

## Signatures to score later

Lists are in `tables/signatures.gmt` (mouse symbols as they appear in this mm10 matrix). The gene-level KL-versus-KP table is `tables/de_kl_vs_kp.tsv.gz`. On another mouse scRNA matrix, score each set as the mean of log-normalized expression of the genes that set and the matrix share. Do not re-cut the gene list on the query.

| signature | genes used | Δ mean log2(FPKM+1) | exact MW p | separation |
|---|---:|---:|---:|---|
| TJ_EPITHELIAL | 18 | +1.627 | 0.00794 | KL>KP |
| TJ_TISMO | 7 | +3.269 | 0.00794 | KL>KP |
| CLDN4_TJ_EDGE | 11 | +2.449 | 0.00794 | KL>KP |
| CLDN4_NEIGHBORHOOD | 107 | +0.165 | 0.4206 | overlap |
| NHEJ_CORE | 11 | +0.165 | 0.0317 | overlap |
| NHEJ_KEGG | 13 | +0.302 | 0.0159 | overlap |
| NHEJ_EXTENDED | 14 | +0.485 | 0.00794 | KL>KP |
| STING_CORE | 5 | +0.475 | 0.00794 | KL>KP |
| IFN_COMPACT | 15 | -0.871 | 0.0317 | overlap |
| IFN_HALLMARK_ALPHA | 90 | -0.631 | 0.0317 | overlap |
| IFN_HALLMARK_GAMMA | 186 | -0.598 | 0.0317 | overlap |
| IFN_OAS | 7 | -0.705 | 0.3095 | overlap |
| APM_MHCI | 16 | -0.500 | 0.3095 | overlap |
| APM_KEGG | 85 | -0.823 | 0.0159 | overlap |

Pre-specified TJ means on this scale are +1.627 (18 epithelial genes) and +3.269 (7-gene TISMO list). The handoff’s TJ figure of +3.03 is a different average.

NHEJ core Δ = +0.165. KEGG NHEJ Δ = +0.302. STING core Δ = +0.475. Compact IFN Δ = -0.871. Hallmark IFN-α Δ = -0.631. Hallmark IFN-γ Δ = -0.598. MHC-I APM Δ = -0.500.

## Cldn4 neighborhood

The co-expression list is built inside each genotype, on log2(FPKM+1), so the KL-versus-KP axis does not pick the genes. A gene qualifies when mean FPKM is at least 1 and the Spearman with Cldn4 is positive in both arms, with the mean of the two arm correlations at or above the threshold. Cldn4 itself is excluded. The operating threshold is the highest grid point with at least 8 genes: ρ ≥ 0.9, 107 genes. The mean score is Δ +0.165, exact p 0.4206, overlap. Epithelial TJ genes inside this list: F11r, Ocln, Tjp3. Epithelial TJ genes outside it: Cgn, Cgnl1, Cldn1, Cldn3, Cldn7, Crb3, Ildr1, Jam2, Jam3, Lsr, Marveld2, Marveld3, Tjp1, Tjp2. Within-arm correlations at n = 5 are unstable, so this list is the co-expression export, while the vendored 11-gene edge and the 18-gene epithelial core are the lists whose mean scores separate KL from KP.

## Sweep: KL>KP on Cldn4 and Tacstd2, with NHEJ lower and IFN higher

The grid has 192 cells: 4 scales (log2 FPKM+1, deposited normalized counts, CPM, median-of-ratios) × mean or median × 2 cohorts (all 10 libraries, or KP without B6AL10-3) × 3 NHEJ lists × 4 IFN lists. **0 cells** have Cldn4 higher, Tacstd2 higher, NHEJ lower, and IFN higher at the same time. NHEJ Δ is positive in 160/192 cells. IFN Δ is positive in 0/192 cells. The pre-specified primary cell (FPKM, mean, all 10, NHEJ core, compact IFN) has Cldn4 +5.570, Tacstd2 +3.238, NHEJ +0.165, IFN -0.871 (2 of 4 limbs in the requested direction). NHEJ Δ is negative in 32 cells, all of them median summaries, and IFN Δ is negative in every one of those cells. The highest Cliff strength on the grid is 1.840 (cpm, median, all_10, NHEJ_CORE, IFN_OAS; NHEJ Δ -0.363, IFN Δ -2.354, 3 of 4 limbs). That cell is the top of a ranking. It is not a selected signature.

## What a later mouse scRNA score can use

- `tables/signatures.gmt` — frozen lists.
- `tables/signature_genes.tsv` — one row per gene per list, including symbols that were absent from this matrix.
- `tables/de_kl_vs_kp.tsv.gz` — every gene, with Δ log2(FPKM+1), Δ log2(normalized+1), Δ log2(CPM+1), Δ log2(median-of-ratios+1), exact MW p, Cliff’s δ, and Welch p. Benjamini–Hochberg columns are on the full gene table. The MW p-values pile up at 0.00794, so that FDR is not a discovery list for Tacstd2 or Cldn4.
- `tables/cldn4_spearman.tsv.gz` — within-arm and 10-library Spearman with Cldn4, so the neighborhood threshold can be moved without re-downloading GEO.
- `tables/sweep_pattern.tsv` — every scale × mean/median × cohort × NHEJ list × IFN list.

Library-level p = 0.00794 is the floor for a 5-vs-5 rank test. It is not an independent-line p-value when libraries inside an arm are near-replicates. Cell-line IFN and APM are transcripts in cultured lines. They are not an immune-exclusion measurement.

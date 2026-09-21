# NHEJ–STING logic wave — TCGA-OV and TCGA-LUAD

Cross-sectional test of one chain, on public tumors only:

CLDN4-low → lower c-NHEJ transcription → higher mutation count → higher STING / IFN / antigen presentation → higher CD8.

The mutation-count piece is the claim in Yamamoto et al., *Molecular Cancer Therapeutics* 2022 (Fig. 4F): in TCGA ovarian tumors, lower CLDN4 tracked higher mutational burden. This note reproduces that comparison and asks the same question in LUAD. LUAD is also reported after a keratin partial correlation, because bulk CLDN4 moves with KRT19.

These cohorts were accrued before immune-checkpoint inhibitors. Nothing here is an ICI-response result.

## What was tested

Two layers, decided before looking at the p-values.

**Layer A. MCT mutation count.** cBioPortal `MUTATION_COUNT` joined to CLDN4 RNA-seq V2 RSEM.

- The published figure splits mutation count at 70 (low 2–69, n=95; high 70–1899, n=107) and also shows CPTAC protein (n=82).
- On the current public files, Pan-Cancer Atlas `MUTATION_COUNT` still has that long tail (OV maximum 1,893). Today's Firehose Legacy `MUTATION_COUNT` tops out at 137 in the OV RNA overlap, so it cannot be the x-axis of Fig. 4F. Both are reported. The figure-scale count is PanCan.
- The 70 cut is the median of the OV PanCan counts in this extract (median = 69). In LUAD the median is 193, so the same numeric cut is an OV cut applied to lung, not a lung median split. Lung is summarized by the continuous Spearman and by its own median split.

**Layer B. Logic-wave scores.** Xena GDC STAR log2(TPM+1), primary tumor (sample code 01), one aliquot per patient (lexicographically first). GSVA follows Hänzelmann et al. 2013 as implemented in Bioconductor GSVA with a Gaussian kernel: bandwidth = gene sd/4, row values = logit of the kernel CDF, column ranks with `ties.method="last"`, tau = 1, `maxDiff=TRUE`, `absRanking=FALSE`. Genes detected in fewer than 10% of tumors were dropped from the background (OV 59,427 → 47,005; LUAD 59,427 → 44,069). Curated set members were kept if present under a GENCODE v36 symbol.

Pre-specified primary endpoints, BH-FDR within cohort (eight tests):

| Endpoint | Genes | Sign if the chain is right |
|---|---|---|
| log10 mutation count (PanCan, count ≥ 1) | — | negative |
| c-NHEJ core GSVA | XRCC6, XRCC5, PRKDC, LIG4, XRCC4, NHEJ1, PAXX, DCLRE1C | positive |
| STING proximal GSVA | CGAS, STING1, TBK1, IKBKE, IRF3, IRF7 | negative |
| IFN-α GSVA | Hallmark Interferon Alpha Response (97/97) | negative |
| IFN-γ GSVA | Hallmark Interferon Gamma Response (200/200) | negative |
| APM GSVA | HLA-A/B/C, B2M, TAP1/2, TAPBP, PSMB8/9/10, PSME1/2, ERAP1/2, NLRC5, CALR, CANX, PDIA3 | negative |
| CD8A | single gene | negative |
| MCP-counter CD8 T | the author file lists **CD8B only** | negative |

GENCODE v36 calls cGAS `CGAS` and STING `STING1`. Both were used. KEGG non-homologous end-joining (includes DNTT) and the single genes TP53BP1 and XRCC1 are sensitivities, outside the FDR family. XRCC1 is BER, not c-NHEJ.

LUAD primary adjustment: partial Spearman of the ranks, residualizing CLDN4 and the endpoint on KRT18 and KRT19 (the locked keratin pair). Sensitivity: KRT8+KRT18+KRT19, ABSOLUTE purity, and log10 mutation count. A second FDR is applied to the KRT18/KRT19 p-values across the same eight endpoints.

## Layer A — mutation burden

### OV mRNA, figure-scale count

PanCan `MUTATION_COUNT` × Firehose RNA-seq V2 RSEM. 290 primary tumors join; 86 have a recorded count of 0 and drop out of log10 and out of the 2–1,899 window. Among the remaining **204**, the median count is **69**.

| Bin | n | Median CLDN4 RSEM | Median log2(RSEM+1) |
|---|---|---|---|
| 2–69 | 103 | 10,840 | 13.40 |
| 70–1,893 | 101 | 9,322 | 13.19 |

The high-count half has the lower CLDN4. Mann–Whitney p = 0.0076. Welch t-test on log2(RSEM+1) p = 0.019. Continuous Spearman rho = **−0.132**, p = 0.060 (n = 204).

The published bins were 95 and 107, with an upper count of 1,899. This extract is 103 and 101, upper count 1,893. The cut, the direction, and the group test agree with Fig. 4F. The sample list is not nucleotide-identical to the December 2020 portal pull, and the continuous correlation is milder than the median split.

The same count joined to PanCan batch-normalized RSEM (n = 205 in the window) gives Mann–Whitney p = 0.0062 and rho = −0.134, p = 0.055.

On Xena STAR, where the logic-wave scores live, the continuous association is a bit stronger: n = 284 tumors with count ≥ 1, rho = **−0.190** (95% CI −0.300 to −0.075), p = 0.0013, q = 0.010. CLDN4 Q1 versus Q4 median log10 counts are 1.92 versus 1.79 (Mann–Whitney p = 0.0089). Purity adjustment leaves rho = −0.190 (n = 277, q = 0.012).

### OV mRNA, today's Firehose count

Firehose Legacy `MUTATION_COUNT` × Firehose RSEM: n = 185, counts 8–137, median 40. Spearman rho = −0.069, p = 0.35. Mann–Whitney on the 70 cut p = 0.34 (156 vs 29, because 70 is no longer the median). This attribute does not carry the Fig. 4F result.

### OV protein

CPTAC CLDN4 from `ov_tcga_protein_quantification`.

| Mutation count | n | Spearman rho | p | Notes |
|---|---|---|---|---|
| Firehose | 84 | −0.183 | 0.095 | nearest to the published n = 82; median count in this overlap is 41, so the 70 cut is 70 vs 14 and is not a median split (Mann–Whitney p = 0.67) |
| PanCan | 112 with count ≥ 1 | −0.040 | 0.67 | 54 vs 58 on the 70 cut, Mann–Whitney p = 0.57 |

The mRNA group difference reproduces. The protein sentence does not clear p < 0.05 on the current CPTAC matrix.

### LUAD extension

PanCan count × Firehose RSEM, n = 505 with count ≥ 1 (5 zeros dropped; maximum count 1,582). Median count = **193**.

- Spearman rho = **+0.080**, p = 0.073. The sign is the opposite of the OV claim.
- OV cut (2–69 vs 70–1,899): 110 vs 395, Mann–Whitney p = 0.72.
- LUAD's own median split: 253 at or below 193 versus 252 above. Median log2(RSEM+1) is 13.12 versus 13.18. Mann–Whitney p = 0.19. The above-median half has the higher CLDN4.

STAR gives the same picture: n = 506, rho = +0.069, p = 0.12, q = 0.14. After KRT18/KRT19, rho = +0.072, p = 0.11. After the keratin triad, rho = +0.075, p = 0.091. Firehose counts in LUAD (n = 230, rho = +0.117, p = 0.077) agree in sign.

The inverse CLDN4–mutation relationship from the ovarian figure does not extend to lung adenocarcinoma.

## Layer B — c-NHEJ, STING/IFN/APM, CD8

Spearman of STAR CLDN4 versus each score. Positive rho means the endpoint is higher when CLDN4 is higher. "Chain" is the pre-specified sign above. q is BH across the eight primary endpoints.

### TCGA-OV (421 primary tumors)

| Endpoint | n | rho | p | q | Chain |
|---|---|---|---|---|---|
| log10 mutation count | 284 | −0.190 | 0.0013 | 0.010 | yes |
| c-NHEJ core GSVA | 421 | +0.017 | 0.73 | 0.73 | flat |
| STING proximal GSVA | 421 | +0.116 | 0.017 | 0.046 | opposite |
| IFN-α GSVA | 421 | +0.083 | 0.089 | 0.14 | opposite, not FDR |
| IFN-γ GSVA | 421 | +0.066 | 0.18 | 0.24 | flat |
| APM GSVA | 421 | +0.032 | 0.52 | 0.59 | flat |
| CD8A | 421 | −0.093 | 0.057 | 0.11 | direction only |
| MCP CD8 (CD8B) | 421 | −0.140 | 0.0040 | 0.016 | yes |

c-NHEJ mRNA does not track CLDN4 in OV. TP53BP1 itself is weakly positive (rho = +0.096, p = 0.049, not in the FDR family). XRCC1 rho = +0.090, p = 0.064. A transcriptomic NHEJ score is not the 53BP1-foci assay in the MCT paper; the flat GSVA means that assay is not visible as a bulk mRNA program here.

STING proximal is higher, not lower, when CLDN4 is higher (q = 0.046). IFN and APM are flat. Partialling log10 mutation count does not turn STING negative (partial rho = +0.152, n = 284, p = 0.011).

CD8B is higher when CLDN4 is lower. That survives purity adjustment (rho = −0.111, n = 409, q = 0.040) and is almost unchanged after adjusting for mutation count (rho = −0.138, n = 284, p = 0.020; q across the seven non-mutation partials = 0.071). CD8A is the same direction and does not meet FDR. The CD8B association is not accounted for by the mutation-count difference.

### TCGA-LUAD (516 primary tumors), with keratin adjustment

CLDN4 is collinear with keratins: Spearman with KRT8 +0.299, KRT18 +0.212, KRT19 **+0.428** (n = 516). That is why the partials are required.

| Endpoint | raw rho (q) | KRT18+KRT19 rho (q) | Chain after keratin adjustment |
|---|---|---|---|
| log10 mutation count | +0.069 (0.14) | +0.072 (0.14) | no inverse association |
| c-NHEJ core GSVA | **+0.174 (2.0×10⁻⁴)** | **+0.133 (0.010)** | yes |
| STING proximal GSVA | **+0.336 (3.6×10⁻¹⁴)** | **+0.220 (3.9×10⁻⁶)** | opposite, and it stays |
| IFN-α GSVA | +0.190 (5.6×10⁻⁵) | +0.071 (0.14) | raw opposite; gone after keratin |
| IFN-γ GSVA | +0.127 (0.0076) | +0.031 (0.48) | raw opposite; gone after keratin |
| APM GSVA | +0.012 (0.79) | −0.063 (0.17) | flat |
| CD8A | −0.098 (0.039) | −0.082 (0.14) | raw yes; not after keratin or purity |
| MCP CD8 (CD8B) | −0.096 (0.039) | −0.074 (0.14) | raw yes; not after keratin or purity |

c-NHEJ also survives the keratin triad (rho = +0.127, q = 0.016), ABSOLUTE purity (rho = +0.155, n = 503, q = 9.9×10⁻⁴), and mutation count (rho = +0.163, n = 506, q = 5.6×10⁻⁴). KEGG NHEJ, which adds DNTT, FEN1, POLL, POLM, MRE11 and RAD50, is the same direction (raw rho = +0.158, p = 3.1×10⁻⁴). CLDN4 Q1 versus Q4: c-NHEJ median −0.233 versus −0.005, Mann–Whitney p = 2.2×10⁻⁴. In LUAD, CLDN4-low tumors have lower c-NHEJ GSVA, and keratin content does not explain it.

STING stays positive after KRT18/KRT19 (rho = +0.220) and after the triad (rho = +0.219). CLDN4-low LUAD has lower, not higher, proximal STING GSVA (Q1−Q4 median rank gap is large and negative; Mann–Whitney p = 1.1×10⁻¹² on the raw score). IFN-α and IFN-γ move with CLDN4 only until keratins are held constant; those two raw associations are epithelial-program confounding. APM is flat either way.

CD8A and CD8B are weakly negative raw (q = 0.039) and lose the FDR after KRT18/KRT19. Purity adjustment removes them (CD8A partial rho = −0.052, p = 0.24). The triad partial for CD8A is rho = −0.087, p = 0.049, q = 0.13.

## How the chain reads

In OV the mutation-count link from Fig. 4F is there for mRNA: the high-count half has lower CLDN4, and on STAR the continuous correlation survives purity adjustment. c-NHEJ mRNA does not accompany it. STING does not rise in CLDN4-low tumors. CD8B does rise, and mutation count does not account for that CD8B difference. The steps are not one chain.

In LUAD the c-NHEJ mRNA link is there and keratin-stable: CLDN4-low tumors score lower for classical NHEJ genes. The mutation-count link is not there (point estimate positive). STING is higher when CLDN4 is higher, with or without keratin adjustment. IFN-α/γ look higher with CLDN4 only in the unadjusted bulk correlation. CD8 does not survive keratin or purity adjustment. Lung does not extend the ovarian mutation-burden claim.

GSVA is a rank score of mRNA. It is not NHEJ repair activity, not 53BP1 foci, and not phosphorylated STING.

## Limits

1. No ICI treatment in either cohort.
2. Bulk RNA. Keratin adjustment is a linear rank residual, not a pure malignant-cell measurement. LUAD CLDN4–KRT19 rho = 0.43, so residual CLDN4 is a thinner variable than raw CLDN4.
3. PanCan mutation counts reproduced the Fig. 4F scale; the current Firehose Legacy count did not. Protein n = 84 is close to 82 and the correlation is not significant.
4. MCP-counter "CD8 T cells" in the author file used here is CD8B alone. CD8A is reported separately.
5. Hallmark sets were taken from Enrichr MSigDB Hallmark 2020; KEGG from Enrichr KEGG 2021 Human. Every requested gene in those sets was present in the filtered STAR matrix.
6. Multiple testing is BH within the eight primary endpoints, per cohort and per adjustment family. Sensitivities are unadjusted p-values.

## Reproduce

```bash
pip install pandas numpy scipy matplotlib
python3 scripts/nhej_sting_wave/download_data.py   # /tmp/nhej_sting_data, or NHEJ_STING_DATA
python3 scripts/nhej_sting_wave/run_analysis.py    # results/nhej_sting_wave/
```

Tables: `mct_mutation_burden.csv`, `mct_ov_protein.csv`, `endpoint_correlations.csv`, `quartile_contrasts.csv`, `geneset_coverage.csv`, `sample_scores_OV.csv`, `sample_scores_LUAD.csv`. Figures: `fig1_mct_mutation_bins.png`, `fig2_spearman_forest.png`, `fig3_cldn4low_minus_high.png`.

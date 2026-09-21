# Open CLDN1 / CLDN3 / CLDN7 knockdown RNA profiles versus CLDN4

Question: in public data, is NHEJ-down plus IFN-up a CLDN4-only knockdown signature?

Answer: no open CLDN1, CLDN3, or CLDN7 series with a confirmed mRNA loss shows that pair. The open CLDN4 knockout profiles do not show that pair either. One confirmed CLDN1 knockdown raises type I IFN while NHEJ stays flat.

Private CLDN1/4/7 co-culture matrices are not in this repository and were not used.

## How a contrast was called

Effect size is the difference in mean log2 expression, perturbed minus control. For counts that is log2(CPM+1). For the Illumina and Agilent arrays it is log2(signal+1). For TPM it is log2(TPM+1). For group-mean FPKM it is log2(FPKM+0.25). For GSE50927 it is the author edgeR logFC.

The replicate-level test is a Welch t-test on the per-sample module mean. A module counts as down when the mean log2 difference is ≤ −0.10 and the median gene is ≤ −0.05. It counts as up when the mean is ≥ 0.10 and the median gene is ≥ 0.05. The joint pattern needs NHEJ core down and Hallmark interferon-alpha up, and the on-target claudin log2 difference must be < −0.5. A nominal joint call would also need ≥3 samples per arm and both Welch p-values < 0.05. No contrast met that bar.

NHEJ core is the canonical c-NHEJ set: PRKDC, XRCC5, XRCC6, XRCC4, LIG4, NHEJ1, DCLRE1C, PAXX, APLF, PNKP, APTX, POLL, POLM, CYREN. IFN-alpha and IFN-gamma are MSigDB hallmark sets (v2023.2). Hallmark DNA repair is reported beside NHEJ so a broad repair drop is visible. BH q-values are across the five modules inside each contrast.

## Confirmed mRNA loss

| Contrast | n | On-target log2 | NHEJ core | IFN-alpha | Joint pattern |
|---|---|---|---|---|---|
| GSE312713 CCA organoid shCLDN1 | 3 vs 3 | −1.83 (p=0.0085) | +0.033 (p=0.59) | **+0.250 (p=3.5×10⁻⁴, q=0.0018)** | absent |
| GSE296175 GIST siCLDN1 | 3 vs 3 | −1.83 (p=3.3×10⁻⁵) | +0.051 (p=0.0036, q=0.018) | −0.019 (p=0.25) | absent |
| GSE234513 DLD-1 shCLDN1 hairpin 1, liver tumors | 2 vs 2 | −0.80 (p=0.19) | −0.57 (p=0.50) | −0.40 (p=0.29) | absent |
| GSE234513 LoVo shCLDN1 hairpin 1, liver tumors | 2 vs 2 | −0.85 (p=0.064) | −1.26 (p=0.23) | −0.66 (p=0.19) | absent |
| GSE159914 Cldn3−/− vs WT adult liver | 3 vs 3 | −3.50 (p=0.0026) | −0.046 (p=0.53) | −0.005 (p=0.89) | absent |
| GSE159914 Cldn3−/− 48 h after hepatectomy | 3 vs 3 | −3.63 (p=9.4×10⁻⁵) | −0.046 (p=0.25) | −0.029 (p=0.61) | absent |
| GSE159914 Cldn3−/− aged liver | 3 vs 3 | −3.31 (p=0.0033) | +0.070 (p=0.71) | +0.074 (p=0.89) | absent |
| GSE273512 MMTV-Neu organoid shCldn7, both hairpins | 8 vs 4 | −2.12 (p=3.7×10⁻⁵) | −0.012 (p=0.81) | +0.100 (p=0.56) | absent |
| GSE273512 sh1 alone | 4 vs 4 | −2.02 (p=4.4×10⁻⁵) | −0.020 (p=0.70) | +0.016 (p=0.94) | absent |
| GSE273512 sh2 alone | 4 vs 4 | −2.21 (p=2.9×10⁻⁵) | −0.004 (p=0.94) | +0.184 (p=0.36) | absent |
| GSE26055 OVCA420 siCLDN7 | 2 vs 2 | −2.68 (p=0.024) | +0.058 (p=0.81) | +0.049 (p=0.87) | absent |
| GSE26055 OVCAR-2 siCLDN7 | 2 vs 2 | −3.16 (p=0.012) | −0.24 (p=0.25) | +0.14 (median +0.030, p=0.33) | absent |
| GSE256329 Cldn7−/− P3 small intestine | 2 vs 2 | −2.85 (p=0.061) | −0.44 (p=0.12) | −0.13 (p=0.72) | absent |
| GSE256329 Cldn7−/− P3 large intestine | 2 vs 2 | −2.29 (p=0.032) | −0.050 (p=0.54) | +0.17 (p=0.34) | absent |

IFN-gamma agrees with IFN-alpha in the organoid CLDN1 knockdown (GSE312713 +0.15, p=0.0040, q=0.010). It does not turn any other confirmed contrast into the joint pattern.

## What each paralog actually does

**CLDN1.** The clean replicated knockdown is the CCA organoid shRNA (GSE312713). Type I IFN goes up and NHEJ does not go down. The GIST siRNA (GSE296175) lowers CLDN1 and moves NHEJ slightly up (+0.05), with IFN flat. In the liver-tumor TPM matrix (GSE234513), columns follow GEO sample order. The first hairpin lowers CLDN1 in DLD-1 and LoVo; the second hairpin does not (log2 differences +0.08 and +0.15). On the working hairpin the point estimates for NHEJ and IFN-alpha are both negative, and Hallmark DNA repair falls with them (−0.42 and −0.62). Those shifts are n=2 and the Welch p-values are 0.19–0.50, so they are not a significant NHEJ-down / IFN-up result.

Two CLDN1 knockout series from the same cholangiocarcinoma study do not lose CLDN1 mRNA. In vitro (GSE312708) control counts are 58,519–71,118 and knockout counts are 43,024–67,825 (log2 difference −0.30). In the xenografts (GSE312711) knockout counts remain 50,574–197,758 against 102,432–264,345 in controls (log2 difference −0.28, p=0.20). Those matrices are not CLDN1-loss transcriptomes. The in vitro matrix does show IFN-alpha up (+0.56, p=0.0085, q=0.021) and a small NHEJ drop (−0.10, p=0.043, q=0.072) while CLDN1 mRNA is retained.

**CLDN3.** GEO has no cancer CLDN3 knockdown RNA-seq. The germline liver knockout (GSE159914) removes Cldn3 (log2 difference about −3.3 to −3.6) in adult liver, regenerating liver, and aged liver. NHEJ and IFN-alpha stay within ±0.08. This is liver, not a tumor knockdown.

**CLDN7.** The replicated cancer knockdown is MMTV-Neu organoids (GSE273512). Both hairpins lower Cldn7 about 4-fold. Pooled NHEJ is −0.012 (p=0.81). Pooled IFN-alpha is +0.10 (p=0.56). Hairpin 2 alone has a larger IFN-alpha point estimate (+0.18, p=0.36) that does not clear the replicate test. Ovarian siRNA arrays (GSE26055; the GEO series matrix is empty, so non-normalized signal was used) lower CLDN7 in both lines. OVCA420 modules are flat. OVCAR-2 NHEJ is −0.24 (p=0.25) and the IFN-alpha median is +0.030, with Hallmark DNA repair also down (−0.37), at n=2. Neonatal intestine arrays (GSE256329) lower Cldn7. Small intestine NHEJ is −0.44 (p=0.12) with IFN-alpha down; large intestine NHEJ is flat and IFN-alpha is +0.17 (p=0.34). The two tissues do not agree, and neither is a cancer knockdown.

## Open CLDN4, scored the same way

These are the public CLDN4-loss expression sets already catalogued. IFN behavior of GSE207704 and GSE50927 was reported before; NHEJ core was not.

| Contrast | CLDN4 | NHEJ core | IFN-alpha |
|---|---|---|---|
| GSE207704 T47D CRISPR, group-mean FPKM, 60 IFN-alpha genes measured | −1.05 | +0.055 (median +0.071) | **−0.37 (median −0.22)** |
| GSE207704 MCF7 CRISPR, group-mean FPKM | −0.75 | +0.053 (median +0.11) | +0.015 (median +0.047) |
| GSE50927 naive lung germline KO, author edgeR | −6.06 | −0.026 (median −0.015); 0 of 12 genes FDR<0.05 | **+0.30 (median +0.19)** |
| GSE22493 SKOV-3 two-color array, knockdown versus CLDN4 overexpression | −1.23 | mean −0.080, median **+0.36** (p=0.75) | −0.054 (p=0.83) |

T47D interferon-alpha falls. MCF7 interferon-alpha on this 60-gene hallmark overlap is flat; an earlier 22-gene IFN-I panel on the same FPKM file was down in both lines. Neither breast line lowers NHEJ mRNA. The lung knockout raises IFN, which matches the earlier mouse Cldn4 audit, and the NHEJ core logFC values sit on zero (POLM −0.47 is the largest drop, author FDR 0.21). The ovarian array compares knockdown with overexpression, and the NHEJ mean and median disagree in sign.

## Held out of the knockdown claim

- GSE10309 is CLDN1 overexpression in CL1-5 lung-adenocarcinoma cells, not a knockdown.
- CLDN1 antibody series (GSE262166 and its treatment subseries, including GSE298420) are not genetic knockdowns. Inside GSE296175 the PDS-0330 arm raises CLDN1 mRNA (+0.23, p=3.2×10⁻⁴) and raises IFN-alpha (+0.66, p=2.7×10⁻⁶). That arm is an antibody exposure.
- GSE274940 is an EpH4 line described as claudin-null. In the count file Cldn3 and Cldn7 fall, and Cldn4 does not (WT 12,380–18,757; null 12,982–21,063). It is not a CLDN4-loss transcriptome. NHEJ is flat (+0.019, p=0.84).
- GSE67164 is a dendritic-cell CRISPR screen, not a CLDN7 knockdown RNA-seq. GSE244629 and GSE180418 mention CLDN7 and are not CLDN7 knockdowns.
- No human lung CLDN4, CLDN1, CLDN3, or CLDN7 knockdown RNA-seq was found.

## Files

- `contrast_summary.tsv` — one row per contrast
- `module_stats.tsv` — NHEJ core, 53BP1/Shieldin, IFN-alpha, IFN-gamma, DNA repair
- `nhej_and_target_genes.tsv` — per-gene log2 differences
- `GSE50927_NHEJ_CORE_author_edgeR.tsv` — author logFC and FDR
- `fig_nhej_ifn_deltas.png` — mean log2 differences
- Script and gene lists: `scripts/cldn_paralog_nhej_ifn/`

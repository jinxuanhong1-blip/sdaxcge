# STING agonist and STING-knockout RNA-seq: CLDN4 in lung cancer lines

Public GEO screen for cGAS/STING agonist or STING-knockout RNA-seq in lung cancer lines, and the CLDN4 result in every open contrast. No private 8-KL matrices. Locked CosMx, concordant-4, GSE137244 KL-versus-KP, TCGA, and TISMO numbers are left as they are.

Searched 21 Sep 2026 (NCBI GEO): STING / STING1 / TMEM173, diABZI, cGAMP, ADU-S100, MSA-2, DMXAA, and “STING agonist” / “STING knockout”, crossed with lung, NSCLC, SCLC, or named lines (A549, H596, H82, Calu-3, LLC).

## Answer

In the one open **direct STING-agonist** RNA-seq of a lung cancer line, Calu-3 plus diABZI, **CLDN4 stays flat** while the interferon program turns on.

| Contrast | Gene | log2FC | p | padj | n |
|---|---|---:|---:|---:|---:|
| Calu-3 diABZI vs DMSO, 6 h | CLDN4 | +0.346 | 0.105 | 0.843 | 2 vs 2 |
| Calu-3 diABZI vs DMSO, 12 h | CLDN4 | −0.057 | 0.716 | 0.956 | 2 vs 2 |
| Calu-3 diABZI vs DMSO, 12 h | CXCL10 | +5.712 | 3.6×10⁻⁸⁶ | 1.1×10⁻⁸² | 2 vs 2 |
| Calu-3 diABZI vs DMSO, 12 h | IFIT1 | +3.273 | 2.4×10⁻⁸⁷ | 8.6×10⁻⁸⁴ | 2 vs 2 |
| Calu-3 diABZI vs DMSO, 12 h | TACSTD2 | +0.245 | 0.154 | 0.667 | 2 vs 2 |
| Calu-3 diABZI vs DMSO, 12 h | KRT8 | −0.0005 | 0.998 | 1.000 | 2 vs 2 |

Source: GSE166209 author DESeq2 tables (`Compare_Treatment_diABZI_*_vs_DMSO`). Calu-3 is a lung adenocarcinoma line. The same baseMean at 6 h and 12 h is the joint model. KRT8 does not move with the agonist, so the flat CLDN4 result is not a keratin-normalization artifact of this matrix.

## STING knockout

**GSE271679, NCI-H596** (lung adenosquamous), unstimulated WT versus CRISPR STING KO, n=4. PyDESeq2 on the deposited featureCounts. This is loss of STING, not an agonist time course.

| Gene | log2FC (KO / WT) | padj | What the counts do |
|---|---:|---:|---|
| CLDN4 | +0.744 | 1.0×10⁻³ | normalized counts 68–116 in WT, 141–167 in KO; the eight samples do not overlap |
| CLDN7 | +1.853 | 2.6×10⁻⁷ | up in KO |
| CLDN1 | −0.519 | 3.7×10⁻⁴ | down in KO |
| TACSTD2 | +1.728 | 1.6×10⁻¹³ | up in KO |
| KRT8 | +0.099 | 0.89 | flat |
| STING1 | −1.340 | 8.7×10⁻⁸ | mRNA down, not abolished |
| CXCL10 | −7.855 | 5.8×10⁻⁵ | WT ~32–44, KO 0 |
| IFIT1 | −2.032 | 8.8×10⁻⁶⁶ | down in KO |

3,678 of 20,120 tested genes have padj < 0.05. Basal ISGs collapse, so the knockout is functional. CLDN4 is higher without STING in this one unstimulated line. Acute HT-DNA in H596 (below) also moves CLDN4 up, so the knockout and the ligand experiment are not a single “STING represses CLDN4” curve.

**GSE244945, SCLC** (NCI-H82, CORL88), author log2 RPKM, n=3, Welch tests. TAS1440 and doxycycline-NOTCH1 are NOTCH activators in that paper, not STING agonists. STING1/TMEM173 mRNA is not lower in the knockout (CORL88 DMSO, Δ log2 −0.09, p=0.44). The functional check is ISG15: NOTCH-on versus NOTCH-off is +4.03 in H82 STING-WT (p=3.5×10⁻⁵) and +2.51 in the STING knockout (p=0.0014). The STING-dependent piece of that dox effect is +1.52 (OLS interaction p=6.1×10⁻⁴).

CLDN4 on the same matrix:

- CORL88 is already high (mean log2 RPKM 6.93–7.65). STING KO versus WT at DMSO: +0.21, p=0.17. TAS1440 raises CLDN4 in WT (+0.72, p=0.013) and in the knockout (+0.40, p=0.030). The STING-dependent piece of TAS1440 is +0.32, p=0.13.
- H82 CLDN4 sits near the floor (means 0.44–1.20, and the matrix contains exact zeros). A NOTCH-by-STING interaction is numerically present (+1.07, p=0.014) on that floor. It is not used as evidence of a CLDN4 program.

CXCL10 is essentially off in both SCLC lines, so it is not the pathway control here. ISG15 is.

## cGAS ligand (HT-DNA), six NSCLC lines

**GSE288796.** Herring-testis DNA versus untreated, n=2. This is a cGAS ligand, not diABZI, cGAMP, ADU-S100, or MSA-2. log2FC is the difference of mean log2(TPM+1). Column numbers follow the GEO sample-title prefixes. Every line induces CXCL10.

| Line | CLDN4 log2FC | CLDN4 Welch p | CXCL10 log2FC | CLDN4 TPM, untreated → DNA |
|---|---:|---:|---:|---|
| HCC827 | −0.164 | 0.24 | +5.29 | 212, 232 → 199, 196 |
| H2228 | +1.295 | 0.0093 | +3.20 | 28, 28 → 70, 72 |
| H1650 | +0.713 | 0.0069 | +3.67 | 143, 150 → 245, 237 |
| H596 | +0.773 | 0.066 | +7.77 | 17, 19 → 31, 31 |
| H1975 | −0.328 | 0.024 | +3.60 | 130, 125 → 102, 101 |
| H358 | +0.002 | 0.98 | +5.51 | 535, 549 → 567, 519 |

The two replicates inside each arm agree. The sign of CLDN4 does not. H358 is already high and stays high. n=2, so these p-values are descriptive.

## Adjacent series

**GSE252340, NCI-H1944**, BAY1217389 (MPS1 inhibitor, 48 h) versus DMSO, n=2. The authors use this compound to make micronuclei and thereby turn on cGAS-STING. It is not a STING agonist and not a STING knockout. PyDESeq2 log2FC for CLDN4 is +0.64 (lfcSE 0.087). Raw CLDN4 counts are 1,403 and 1,291 (DMSO) versus 2,357 and 2,259 (BAY). IFIT1 log2FC is +3.99 and CXCL10 log2FC is +8.42. PyDESeq2 warned that with residual df < 3 the dispersion prior is poorly estimated, so the very small padj values are not the claim. H1944 is a KRAS/LKB1 line; this contrast is still an MPS1-inhibitor treatment.

**GSE134129.** In vivo LLC tumors, NanoString PanCancer panel (750 genes), cGAMP-treated tumors in WT mice versus tumors from STING-knockout mice. **Cldn4 is not on the panel**, so there is no CLDN4 test. Cxcl10, Ifnb1, Ifit1, and Tmem173 are present. Mean log2(count+1) for Cxcl10 is higher in the STING columns than the control columns by 1.60, and lower in the knockout columns. That only shows the panel sees the pathway.

**GSE166209 mouse lung** (K18-hACE2 tissue, diABZI versus PBS, n=5) is the same series as Calu-3 and is not a cancer line. Cldn4 at 12 h: log2FC +1.03, p=0.023, padj=0.32, baseMean 29.5. At 6 h: log2FC +0.21, padj=0.86. Low expression, and the 12 h nominal p does not survive the author FDR.

## Not tested

No open RNA-seq of diABZI, cGAMP, ADU-S100, or MSA-2 was found in A549, H1299, H1944, or a KL/KP pair. Left out of the CLDN4 test on purpose:

- GSE254174, U2OS and BJ fibroblasts, cGAMP/DMXAA
- GSE305239, H4 glioblastoma, diABZI
- GSE308610, THP-1 and RPMI-8226, ADU-S100
- GSE147085, FaDu HNSCC STING knockout
- GSE269551, H2122 SMARCA4 knockout (STING is the paper’s mechanism, not the contrast)
- GSE137244, mouse KL versus KP baseline, already reported
- GSE271679 non-lung lines (SCC25, OE21, Detroit 562) and Calu-3 WT-only samples in that series

## Files

- `results/tables/inventory.tsv`
- `results/tables/contrasts.tsv` (every gene and contrast)
- `results/tables/h596_stingko_normalized_counts.tsv`
- `results/tables/highlights.json`
- `results/figures/cldn4_sting_lung_log2fc.png` (and `.pdf`)

Figure: CLDN4 and CXCL10 for Calu-3 diABZI, H596 STING knockout, and the six HT-DNA lines. Error bars are the author or PyDESeq2 lfcSE, or the SE of the log2(TPM+1) difference for HT-DNA. SCLC and the MPS1-inhibitor arm are in the tables only.

Reproduce:

```bash
python3 scripts/sting_lung_cldn4.py
```

Requires pandas, scipy, matplotlib, openpyxl, and pydeseq2. The script downloads GEO supplementary files into `data/geo_cache/` (gitignored).

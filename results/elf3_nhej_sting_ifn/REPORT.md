# ELF3 loss vs NHEJ, STING, and IFN

Public ELF3 ChIP-seq / CUT&Tag and ELF3-loss transcriptomes. CLDN4 is the positive-control node of the ELF3 epithelial network. Private KL matrices are not used.

## Verdict

classical NHEJ: no QC-pass primary dataset meets the UP or DOWN rule (4 contrasts: A549_shELF3 -0.010 (NULL), HBDEC2_KO 0.158 (WEAK), HBDEC2_KD -0.107 (NULL), JEG3_KO 0.007 (NULL)). cGAS–STING: no QC-pass primary dataset meets the UP or DOWN rule (4 contrasts: A549_shELF3 0.119 (NULL), HBDEC2_KO -0.008 (NULL), HBDEC2_KD -0.051 (NULL), JEG3_KO 0.101 (INSUFFICIENT)). interferon-alpha: no QC-pass primary dataset meets the UP or DOWN rule (4 contrasts: A549_shELF3 0.091 (WEAK), HBDEC2_KO -0.075 (NULL), HBDEC2_KD -0.011 (NULL), JEG3_KO -0.186 (WEAK)). interferon-gamma: no QC-pass primary dataset meets the UP or DOWN rule (4 contrasts: A549_shELF3 0.074 (WEAK), HBDEC2_KO -0.116 (NULL), HBDEC2_KD -0.022 (WEAK), JEG3_KO 0.004 (NULL)). Apical-junction control on the same contrasts: A549_shELF3 0.080 (WEAK), HBDEC2_KO -0.152 (NULL), HBDEC2_KD 0.001 (WEAK), JEG3_KO -0.030 (NULL). CLDN4 log2FC in those contrasts: A549_shELF3 -0.312; HBDEC2_KO -2.224; HBDEC2_KD 0.091; JEG3_KO -0.096. The one cytokine contrast that passes the log2FC gate (synovial siELF3, 16 h TNF+IL-17A, not lung) is NHEJ_CORE 0.464 (NULL), STING_CORE -0.004 (NULL), HALLMARK_INTERFERON_ALPHA_RESPONSE -0.266 (DOWN), HALLMARK_INTERFERON_GAMMA_RESPONSE -0.150 (WEAK). Interferon-alpha moves down there, not up. ELF3 raw p in that contrast is 0.21, so the knockdown itself is not significant and the interferon call is not secure. The 16 h TNF-only arm moves the same way but fails the ELF3 log2FC gate (−0.22) and is not counted. Mouse tubular Elf3 deletion on TNFα+IFNγ (ELF3 log2FC −2.45, raw p 0.002) does not call NHEJ, STING, or interferon UP or DOWN.

Binding is a separate result from those expression calls. ELF3 CUT&Tag TSS±2 kb overlap for classical NHEJ is ECC4 5/13 (p=0.025); A99 5/13 (p=0.009); DMS53 7/13 (p=0.003). STING and interferon-alpha sets are not enriched at that window. The 232-gene direct signature (near a peak and down in at least two lines) is NHEJ: DCLRE1C; STING: none; interferon-alpha: none; CLDN4: CLDN4. Peaks on a few NHEJ genes do not become a downregulated NHEJ program after ELF3 loss.

Calls use detected genes only. UP/DOWN requires median |log2FC| ≥ 0.25, same-direction fraction ≥ 0.60, and two-sided Mann–Whitney p < 0.05 versus the rest of the detected transcriptome. WEAK means p < 0.05 with a smaller shift. NULL means the set does not move relative to the transcriptome.

## ELF3 knockdown QC

| Dataset | Role | n KD / n ctrl | ELF3 log2FC | ELF3 raw p | QC | CLDN4 log2FC |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| A549_shELF3 | primary | 3 / 3 | -2.287 | 2.94e-04 | PASS | -0.312 |
| HBDEC2_KO | primary | 3 / 3 | -3.126 | 4.91e-04 | PASS | -2.224 |
| HBDEC2_KD | primary | 3 / 3 | -1.156 | 0.045 | PASS | 0.091 |
| JEG3_KO | primary | 3 / 3 | -1.340 | 0.012 | PASS | -0.096 |
| CRPC_shELF3 | supporting_n2 | 2 / 2 | -1.982 | 6.44e-04 | PASS | -1.347 |
| synov_siELF3_basal | primary | 4 / 4 | 0.360 | 0.045 | FAIL | 0.014 |
| synov_siELF3_TNF16 | stimulated | 4 / 4 | -0.222 | 0.523 | FAIL | 0.171 |
| synov_siELF3_TNF_IL17_16 | stimulated | 3 / 3 | -0.860 | 0.211 | PASS | 0.498 |
| mTEC_KO_TNFa_IFNg | stimulated | 3 / 3 | -2.452 | 0.002 | PASS | -0.308 |

QC PASS requires ELF3 log2FC ≤ −0.5. That is a point-estimate gate, not an FDR gate: A549 ELF3 raw p is small, and the genome-wide BH FDR can still sit above 0.05 because thousands of genes are tested at n=3. Contrasts that fail the gate are shown and are excluded from the basal cross-dataset sentence.

GSE303076 has two replicates per arm. A Mann–Whitney p-value there uses genes as the sample, so tight replicates produce very small p-values. Read that row as a direction.

## Program shifts (detected genes)

| Dataset | QC | Set | n | median log2FC | frac up | MW p | Wilcoxon p | Call |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| A549_shELF3 | PASS | NHEJ_CORE | 12 | -0.010 | 0.417 | 0.916 | 0.970 | NULL |
| A549_shELF3 | PASS | STING_CORE | 11 | 0.119 | 0.727 | 0.082 | 0.083 | NULL |
| A549_shELF3 | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 93 | 0.091 | 0.613 | 0.024 | 0.026 | WEAK |
| A549_shELF3 | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 194 | 0.074 | 0.613 | 0.006 | 0.005 | WEAK |
| A549_shELF3 | PASS | HALLMARK_DNA_REPAIR | 141 | -0.037 | 0.440 | 0.041 | 0.153 | WEAK |
| A549_shELF3 | PASS | HALLMARK_APICAL_JUNCTION | 191 | 0.080 | 0.618 | 5.60e-04 | 3.21e-04 | WEAK |
| A549_shELF3 | PASS | HALLMARK_E2F_TARGETS | 188 | -0.061 | 0.378 | 4.77e-05 | 8.03e-04 | WEAK |
| HBDEC2_KO | PASS | NHEJ_CORE | 13 | 0.158 | 0.615 | 0.002 | 0.068 | WEAK |
| HBDEC2_KO | PASS | STING_CORE | 11 | -0.008 | 0.455 | 0.284 | 0.966 | NULL |
| HBDEC2_KO | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 97 | -0.075 | 0.454 | 0.065 | 0.151 | NULL |
| HBDEC2_KO | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 198 | -0.116 | 0.379 | 0.254 | 0.002 | NULL |
| HBDEC2_KO | PASS | HALLMARK_DNA_REPAIR | 148 | -0.003 | 0.493 | 3.97e-08 | 0.860 | WEAK |
| HBDEC2_KO | PASS | HALLMARK_APICAL_JUNCTION | 200 | -0.152 | 0.325 | 0.836 | 2.01e-04 | NULL |
| HBDEC2_KO | PASS | HALLMARK_E2F_TARGETS | 199 | 0.186 | 0.749 | 3.42e-38 | 4.14e-15 | WEAK |
| HBDEC2_KD | PASS | NHEJ_CORE | 13 | -0.107 | 0.154 | 0.420 | 0.008 | NULL |
| HBDEC2_KD | PASS | STING_CORE | 11 | -0.051 | 0.182 | 0.481 | 0.067 | NULL |
| HBDEC2_KD | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 97 | -0.011 | 0.485 | 0.091 | 0.475 | NULL |
| HBDEC2_KD | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 198 | -0.022 | 0.455 | 0.007 | 0.372 | WEAK |
| HBDEC2_KD | PASS | HALLMARK_DNA_REPAIR | 148 | -0.063 | 0.365 | 0.105 | 0.002 | NULL |
| HBDEC2_KD | PASS | HALLMARK_APICAL_JUNCTION | 200 | 0.001 | 0.505 | 5.87e-05 | 0.735 | WEAK |
| HBDEC2_KD | PASS | HALLMARK_E2F_TARGETS | 199 | -0.080 | 0.271 | 0.322 | 6.60e-12 | NULL |
| JEG3_KO | PASS | NHEJ_CORE | 12 | 0.007 | 0.500 | 0.712 | 0.970 | NULL |
| JEG3_KO | PASS | STING_CORE | 4 | 0.101 | 0.500 | NA | NA | INSUFFICIENT |
| JEG3_KO | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 63 | -0.186 | 0.333 | 0.005 | 0.002 | WEAK |
| JEG3_KO | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 128 | 0.004 | 0.500 | 0.409 | 0.512 | NULL |
| JEG3_KO | PASS | HALLMARK_DNA_REPAIR | 140 | -0.028 | 0.421 | 0.820 | 0.073 | NULL |
| JEG3_KO | PASS | HALLMARK_APICAL_JUNCTION | 140 | -0.030 | 0.479 | 0.260 | 0.749 | NULL |
| JEG3_KO | PASS | HALLMARK_E2F_TARGETS | 197 | -0.045 | 0.431 | 0.396 | 0.004 | NULL |
| CRPC_shELF3 | PASS | NHEJ_CORE | 11 | 0.152 | 0.727 | 0.259 | 0.240 | NULL |
| CRPC_shELF3 | PASS | STING_CORE | 7 | -0.222 | 0.429 | 0.435 | NA | NULL |
| CRPC_shELF3 | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 64 | -0.264 | 0.359 | 1.92e-04 | 0.004 | DOWN |
| CRPC_shELF3 | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 120 | -0.222 | 0.383 | 1.91e-07 | 9.03e-05 | WEAK |
| CRPC_shELF3 | PASS | HALLMARK_DNA_REPAIR | 141 | 0.102 | 0.603 | 0.054 | 0.006 | NULL |
| CRPC_shELF3 | PASS | HALLMARK_APICAL_JUNCTION | 117 | -0.038 | 0.479 | 0.015 | 0.135 | WEAK |
| CRPC_shELF3 | PASS | HALLMARK_E2F_TARGETS | 198 | 0.606 | 0.848 | 7.10e-43 | 1.36e-26 | UP |
| synov_siELF3_basal | FAIL | NHEJ_CORE | 13 | 0.156 | 0.615 | 0.792 | 0.376 | NULL |
| synov_siELF3_basal | FAIL | STING_CORE | 10 | -0.145 | 0.200 | 0.170 | 0.375 | NULL |
| synov_siELF3_basal | FAIL | HALLMARK_INTERFERON_ALPHA_RESPONSE | 90 | -0.054 | 0.444 | 0.165 | 0.979 | NULL |
| synov_siELF3_basal | FAIL | HALLMARK_INTERFERON_GAMMA_RESPONSE | 166 | -0.013 | 0.482 | 0.204 | 0.561 | NULL |
| synov_siELF3_basal | FAIL | HALLMARK_DNA_REPAIR | 145 | -0.085 | 0.400 | 1.31e-04 | 0.031 | WEAK |
| synov_siELF3_basal | FAIL | HALLMARK_APICAL_JUNCTION | 145 | 0.153 | 0.634 | 0.019 | 6.23e-05 | WEAK |
| synov_siELF3_basal | FAIL | HALLMARK_E2F_TARGETS | 196 | 0.158 | 0.679 | 6.20e-05 | 1.19e-08 | WEAK |
| synov_siELF3_TNF16 | FAIL | NHEJ_CORE | 13 | 0.041 | 0.615 | 0.750 | 0.635 | NULL |
| synov_siELF3_TNF16 | FAIL | STING_CORE | 10 | -0.506 | 0.200 | 0.003 | 0.037 | DOWN |
| synov_siELF3_TNF16 | FAIL | HALLMARK_INTERFERON_ALPHA_RESPONSE | 95 | -0.693 | 0.074 | 2.25e-39 | 4.33e-16 | DOWN |
| synov_siELF3_TNF16 | FAIL | HALLMARK_INTERFERON_GAMMA_RESPONSE | 173 | -0.416 | 0.179 | 2.09e-36 | 2.34e-17 | DOWN |
| synov_siELF3_TNF16 | FAIL | HALLMARK_DNA_REPAIR | 146 | -0.013 | 0.479 | 0.167 | 0.741 | NULL |
| synov_siELF3_TNF16 | FAIL | HALLMARK_APICAL_JUNCTION | 146 | 0.029 | 0.527 | 0.690 | 0.172 | NULL |
| synov_siELF3_TNF16 | FAIL | HALLMARK_E2F_TARGETS | 196 | 0.281 | 0.770 | 2.19e-16 | 3.19e-19 | UP |
| synov_siELF3_TNF_IL17_16 | PASS | NHEJ_CORE | 13 | 0.464 | 0.846 | 0.160 | 0.003 | NULL |
| synov_siELF3_TNF_IL17_16 | PASS | STING_CORE | 11 | -0.004 | 0.455 | 0.298 | 0.831 | NULL |
| synov_siELF3_TNF_IL17_16 | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 95 | -0.266 | 0.347 | 2.98e-17 | 4.10e-06 | DOWN |
| synov_siELF3_TNF_IL17_16 | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 175 | -0.150 | 0.423 | 4.82e-14 | 0.005 | WEAK |
| synov_siELF3_TNF_IL17_16 | PASS | HALLMARK_DNA_REPAIR | 148 | 0.230 | 0.655 | 0.451 | 6.60e-04 | NULL |
| synov_siELF3_TNF_IL17_16 | PASS | HALLMARK_APICAL_JUNCTION | 146 | 0.224 | 0.692 | 0.944 | 2.00e-07 | NULL |
| synov_siELF3_TNF_IL17_16 | PASS | HALLMARK_E2F_TARGETS | 198 | 0.567 | 0.843 | 2.02e-17 | 3.94e-25 | UP |
| mTEC_KO_TNFa_IFNg | PASS | NHEJ_CORE | 13 | -0.144 | 0.385 | 0.350 | 0.216 | NULL |
| mTEC_KO_TNFa_IFNg | PASS | STING_CORE | 9 | -0.022 | 0.444 | 0.486 | 0.734 | NULL |
| mTEC_KO_TNFa_IFNg | PASS | HALLMARK_INTERFERON_ALPHA_RESPONSE | 76 | -0.136 | 0.342 | 0.002 | 7.09e-04 | WEAK |
| mTEC_KO_TNFa_IFNg | PASS | HALLMARK_INTERFERON_GAMMA_RESPONSE | 138 | -0.102 | 0.420 | 0.017 | 0.004 | WEAK |
| mTEC_KO_TNFa_IFNg | PASS | HALLMARK_DNA_REPAIR | 140 | 0.011 | 0.507 | 0.507 | 0.964 | NULL |
| mTEC_KO_TNFa_IFNg | PASS | HALLMARK_APICAL_JUNCTION | 155 | 0.090 | 0.619 | 1.91e-05 | 7.73e-04 | WEAK |
| mTEC_KO_TNFa_IFNg | PASS | HALLMARK_E2F_TARGETS | 192 | -0.130 | 0.344 | 6.25e-08 | 3.75e-09 | WEAK |

Apical junction and E2F/G2M are controls. The NEC paper’s ELF3 signature is G2M/E2F, not interferon. Hallmark DNA repair is the broad proliferation-linked repair set, not classical NHEJ.

## Focus genes

log2FC (KD − control). FDR is Benjamini–Hochberg within detected genes of that contrast. Blank FDR means the gene was below the detection floor.

| Gene | A549_shELF3 | HBDEC2_KO | HBDEC2_KD | JEG3_KO | CRPC_shELF3 | synov_siELF3_basal | synov_siELF3_TNF16 | synov_siELF3_TNF_IL17_16 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ELF3 | -2.287 | -3.126* | -1.156 | -1.340 | -1.982 | 0.360 | -0.222 | -0.860 |
| CLDN4 | -0.312 | -2.224* | 0.091 | -0.096 | -1.347 | 0.014 | 0.171 | 0.498 |
| TACSTD2 | 0.292 | -1.471* | -0.208 | 0.380 | -0.845 | -0.616 | 0.020 | -0.128 |
| XRCC5 | -0.187 | -0.010 | -0.013 | -0.036 | 0.744 | -0.175 | 0.010 | 0.468 |
| XRCC6 | -0.019 | 0.158 | -0.012 | -0.058 | 0.934 | -0.348 | -0.289 | 0.748 |
| PRKDC | -0.141 | 0.370* | -0.422* | 0.107 | 1.242 | 0.083 | 0.390 | 0.164 |
| XRCC4 | -0.639 | 0.463 | -0.165 | 0.581 | 0.133 | 0.265 | 0.160 | 0.563 |
| LIG4 | 0.352 | -0.206 | 0.150 | 0.131 | 0.178 | 0.354 | 0.293 | 0.854 |
| NHEJ1 | 0.548 | -0.085 | -0.508 | -0.064 | 0.094 | 0.156 | -0.434 | 0.331 |
| CGAS | -0.122 | -1.112* | -0.358 | 0.190 | 0.032 | 0.387 | 0.322 | 0.413 |
| STING1 | 0.087 | 0.742 | -0.016 | 0.148 | -0.363 | -0.618 | -0.829 | -0.004 |
| TBK1 | -0.091 | -0.008 | 0.134 | -0.069 | 0.231 | -0.043 | -0.170 | 0.996 |
| IRF3 | 0.522 | 0.128 | -0.048 | -0.120 | -0.222 | -0.026 | -0.107 | 0.262 |
| IFNB1 | 0.034 | -0.037 | -0.057 | 0.000 | NA | -0.071 | -0.065 | 0.000 |
| ISG15 | 0.226 | -0.789 | 0.132 | 0.088 | -1.216 | -0.745 | -2.083 | -2.235 |
| MX1 | 0.063 | -3.028 | 0.484 | -0.195 | -0.167 | -0.313 | -1.774 | -0.612 |
| OAS1 | -1.174 | -0.867 | -1.185 | 0.044 | -1.535 | -0.208 | -1.462 | -0.612 |
| IFIT1 | -0.113 | -0.913 | 0.418 | -1.299 | -0.708 | -0.460 | -1.950 | -0.614 |
| STAT1 | 0.927 | 0.156 | 0.009 | 0.298 | 0.348 | -0.160 | -0.923 | 0.068 |
| CXCL10 | -0.038 | 0.275 | -0.194 | 0.000 | -0.332 | 0.012 | -1.838 | -1.696 |
| CD274 | 0.326 | -0.504 | 0.053 | 0.149 | -0.117 | 1.038 | -0.174 | 0.075 |

\* FDR < 0.05 among detected genes. Mouse TEC (stimulated) is in `gene_log2fc.csv`, not this human table.

## Binding

ChIP-Atlas hg38 ELF3, 5 kb gene window, score > 0. Individual libraries jump from 0 to a MACS score above 50, so the >0 call is a peak call, not a tiny-score artifact. Every row in that table has ELF3|Average > 0: the file is the set of genes with a peak in at least one catalogued experiment (about 14,000 genes), not the whole genome. Fisher tests below are inside that table. There is no lung / NSCLC / A549 ELF3 ChIP-seq in it (CFPAC-1, ESO-26, HBDEC2, HepG2 only). CFPAC-1 and HepG2 peaks cover most of the table, so a high bound fraction there is broad occupancy, not a selective NHEJ or IFN program.

| Experiment | Set | bound / tested | fraction | background fraction | Fisher p |
| --- | --- | ---: | ---: | ---: | ---: |
| CFPAC1_any | NHEJ_CORE | 11/12 | 0.917 | 0.826 | 0.704 |
| CFPAC1_any | STING_CORE | 9/9 | 1.000 | 0.826 | 0.375 |
| CFPAC1_any | HALLMARK_APICAL_JUNCTION | 126/147 | 0.857 | 0.826 | 0.381 |
| CFPAC1_any | HALLMARK_INTERFERON_ALPHA_RESPONSE | 88/92 | 0.957 | 0.825 | 2.67e-04 |
| CFPAC1_any | HALLMARK_INTERFERON_GAMMA_RESPONSE | 157/167 | 0.940 | 0.825 | 1.98e-05 |
| CFPAC1_any | CLDN4 | 1/1 | 1.000 | 0.826 | 1.000 |
| CFPAC1_any | TACSTD2 | 1/1 | 1.000 | 0.826 | 1.000 |
| ESO26 | NHEJ_CORE | 2/12 | 0.167 | 0.050 | 0.121 |
| ESO26 | STING_CORE | 0/9 | 0.000 | 0.051 | 1.000 |
| ESO26 | HALLMARK_APICAL_JUNCTION | 8/147 | 0.054 | 0.051 | 0.849 |
| ESO26 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 4/92 | 0.043 | 0.051 | 1.000 |
| ESO26 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 9/167 | 0.054 | 0.051 | 0.858 |
| ESO26 | CLDN4 | 0/1 | 0.000 | 0.051 | 1.000 |
| ESO26 | TACSTD2 | 0/1 | 0.000 | 0.051 | 1.000 |
| HBDEC2 | NHEJ_CORE | 1/12 | 0.083 | 0.019 | 0.209 |
| HBDEC2 | STING_CORE | 0/9 | 0.000 | 0.019 | 1.000 |
| HBDEC2 | HALLMARK_APICAL_JUNCTION | 6/147 | 0.041 | 0.019 | 0.066 |
| HBDEC2 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 1/92 | 0.011 | 0.019 | 1.000 |
| HBDEC2 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 4/167 | 0.024 | 0.019 | 0.568 |
| HBDEC2 | CLDN4 | 0/1 | 0.000 | 0.019 | 1.000 |
| HBDEC2 | TACSTD2 | 0/1 | 0.000 | 0.019 | 1.000 |
| HepG2_any | NHEJ_CORE | 12/12 | 1.000 | 0.811 | 0.140 |
| HepG2_any | STING_CORE | 5/9 | 0.556 | 0.811 | 0.072 |
| HepG2_any | HALLMARK_APICAL_JUNCTION | 97/147 | 0.660 | 0.812 | 1.67e-05 |
| HepG2_any | HALLMARK_INTERFERON_ALPHA_RESPONSE | 66/92 | 0.717 | 0.811 | 0.031 |
| HepG2_any | HALLMARK_INTERFERON_GAMMA_RESPONSE | 119/167 | 0.713 | 0.812 | 0.002 |
| HepG2_any | CLDN4 | 0/1 | 0.000 | 0.811 | 0.189 |
| HepG2_any | TACSTD2 | 1/1 | 1.000 | 0.811 | 1.000 |

Genome-wide test: called peaks overlapped with any RefSeq TSS ±2 kb (hg38). DMS53 is SCLC; ECC4 and A99 are gastrointestinal NEC. HepG2 is ENCODE IDR (ENCFF080FAU); that file has ~40,000 peaks and a ~40% background rate, so its enrichment is not selective. The NEC/SCLC CUT&Tag background is 11–17%.

TSS ±2 kb peaks in at least two of ECC4, A99, and DMS53: PRKDC (3/3), XRCC4 (3/3), DCLRE1C (3/3), PNKP (3/3), POLM (3/3), STING1 (2/3), TREX1 (2/3), CLDN4 (2/3), TACSTD2 (2/3). A99 CLDN4 is outside ±2 kb and inside ±10 kb.

| Experiment | Peaks | Set | TSS±2 kb bound / tested | fraction | background | Fisher p |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| ECC4_CUTTag | 5067 | NHEJ_CORE | 5/13 | 0.385 | 0.139 | 0.025 |
| ECC4_CUTTag | 5067 | STING_CORE | 1/11 | 0.091 | 0.139 | 1.000 |
| ECC4_CUTTag | 5067 | HALLMARK_APICAL_JUNCTION | 31/198 | 0.157 | 0.139 | 0.470 |
| ECC4_CUTTag | 5067 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 11/97 | 0.113 | 0.139 | 0.557 |
| ECC4_CUTTag | 5067 | CLDN4 | 1/1 | 1.000 | 0.139 | 0.139 |
| ECC4_CUTTag | 5067 | TACSTD2 | 0/1 | 0.000 | 0.139 | 1.000 |
| A99_CUTTag | 3837 | NHEJ_CORE | 5/13 | 0.385 | 0.108 | 0.009 |
| A99_CUTTag | 3837 | STING_CORE | 2/11 | 0.182 | 0.108 | 0.335 |
| A99_CUTTag | 3837 | HALLMARK_APICAL_JUNCTION | 24/198 | 0.121 | 0.108 | 0.492 |
| A99_CUTTag | 3837 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 9/97 | 0.093 | 0.108 | 0.744 |
| A99_CUTTag | 3837 | CLDN4 | 0/1 | 0.000 | 0.108 | 1.000 |
| A99_CUTTag | 3837 | TACSTD2 | 1/1 | 1.000 | 0.108 | 0.108 |
| DMS53_CUTTag | 6415 | NHEJ_CORE | 7/13 | 0.538 | 0.170 | 0.003 |
| DMS53_CUTTag | 6415 | STING_CORE | 3/11 | 0.273 | 0.170 | 0.413 |
| DMS53_CUTTag | 6415 | HALLMARK_APICAL_JUNCTION | 45/198 | 0.227 | 0.170 | 0.037 |
| DMS53_CUTTag | 6415 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 14/97 | 0.144 | 0.171 | 0.589 |
| DMS53_CUTTag | 6415 | CLDN4 | 1/1 | 1.000 | 0.170 | 0.170 |
| DMS53_CUTTag | 6415 | TACSTD2 | 1/1 | 1.000 | 0.170 | 0.170 |
| HBDEC2_ChIP | 13960 | NHEJ_CORE | 1/13 | 0.077 | 0.033 | 0.355 |
| HBDEC2_ChIP | 13960 | STING_CORE | 0/11 | 0.000 | 0.033 | 1.000 |
| HBDEC2_ChIP | 13960 | HALLMARK_APICAL_JUNCTION | 9/198 | 0.045 | 0.033 | 0.315 |
| HBDEC2_ChIP | 13960 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 4/97 | 0.041 | 0.033 | 0.566 |
| HBDEC2_ChIP | 13960 | CLDN4 | 0/1 | 0.000 | 0.033 | 1.000 |
| HBDEC2_ChIP | 13960 | TACSTD2 | 0/1 | 0.000 | 0.033 | 1.000 |
| HepG2_ENCODE_IDR | 39978 | NHEJ_CORE | 12/13 | 0.923 | 0.403 | 1.50e-04 |
| HepG2_ENCODE_IDR | 39978 | STING_CORE | 5/11 | 0.455 | 0.403 | 0.765 |
| HepG2_ENCODE_IDR | 39978 | HALLMARK_APICAL_JUNCTION | 91/198 | 0.460 | 0.403 | 0.110 |
| HepG2_ENCODE_IDR | 39978 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 58/97 | 0.598 | 0.402 | 1.12e-04 |
| HepG2_ENCODE_IDR | 39978 | CLDN4 | 0/1 | 0.000 | 0.403 | 1.000 |
| HepG2_ENCODE_IDR | 39978 | TACSTD2 | 0/1 | 0.000 | 0.403 | 1.000 |

Author gene lists from the same NEC study (supplement tables S7–S9) are one-sided for expression: S8 is fold change < 0.75 only. Full siELF3 count matrices were not deposited in GSE190618 (that accession has baseline RNA-seq plus CUT&Tag peaks).

| List | n genes | Set | overlap | genes |
| --- | ---: | --- | ---: | --- |
| down_ECC4 | 2455 | NHEJ_CORE | 4/13 | APTX,DCLRE1C,XRCC5,XRCC6 |
| down_ECC4 | 2455 | STING_CORE | 2/11 | ENPP1,IKBKE |
| down_ECC4 | 2455 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 12/97 | ADAR,CMTR1,EIF2AK2,ELF1,HLA-C,NCOA7,PNPT1,SAMD9,TDRD7,TMEM140,TRIM25,USP18 |
| down_ECC4 | 2455 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 24/200 | ADAR,BPGM,CASP3,CASP4,CMTR1,EIF2AK2,GCH1,HIF1A,HLA-A,HLA-B,HLA-G,IFNAR2,MX2,NFKB1,PDE4B,PNPT1,PTPN1,RAPGEF6,SPPL2A,TDRD7,TNFSF10,TOR1B,TRIM25,USP18 |
| down_ECC4 | 2455 | HALLMARK_APICAL_JUNCTION | 25/200 | AKT2,AKT3,AMIGO2,CDH8,CLDN4,CNTN1,EXOC4,FLNC,INPPL1,INSIG1,LIMA1,MAP3K20,MPZL1,MYH10,NFASC,PIK3CB,PIK3R3,RAC2,RASA1,SHROOM2,SKAP2,SYK,TMEM8B,VCL,VWF |
| down_ECC4 | 2455 | HALLMARK_E2F_TARGETS | 54/200 | ASF1B,ATAD2,AURKA,AURKB,BIRC5,CBX5,CCNE1,CCP110,CDC25A,CDC25B,CDCA3,CDKN2C,CENPE,CENPM,CNOT9,DEK,DEPDC1,DIAPH3,DSCC1,DUT,E2F8,EED,ESPL1,GINS3,KIF18B,KIF2C,KIF4A,MAD2L1,MCM2,MCM5... |
| down_ECC4 | 2455 | CLDN4 | 1/1 | CLDN4 |
| down_ECC4 | 2455 | TACSTD2 | 0/1 |  |
| down_A99 | 199 | NHEJ_CORE | 0/13 |  |
| down_A99 | 199 | STING_CORE | 0/11 |  |
| down_A99 | 199 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 2/97 | IFITM2,TRIM14 |
| down_A99 | 199 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 2/200 | IFITM2,TRIM14 |
| down_A99 | 199 | HALLMARK_APICAL_JUNCTION | 3/200 | AKT3,CDH11,CLDN11 |
| down_A99 | 199 | HALLMARK_E2F_TARGETS | 2/200 | POLA2,RRM2 |
| down_A99 | 199 | CLDN4 | 0/1 |  |
| down_A99 | 199 | TACSTD2 | 0/1 |  |
| down_DMS53 | 2701 | NHEJ_CORE | 5/13 | APTX,DCLRE1C,NHEJ1,PRKDC,XRCC6 |
| down_DMS53 | 2701 | STING_CORE | 5/11 | ENPP1,IKBKE,IRF3,STING1,TREX1 |
| down_DMS53 | 2701 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 19/97 | CASP8,CMTR1,GBP2,HERC6,HLA-C,IFI35,IFIT3,IFITM2,IFITM3,OAS1,OGFR,PARP12,PLSCR1,PROCR,PSMB9,SAMD9,TAP1,TRIM25,TRIM26 |
| down_DMS53 | 2701 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 42/200 | CASP3,CASP4,CASP8,CD274,CMTR1,HERC6,HLA-A,HLA-B,HLA-DRB1,IFI35,IFIT3,IFITM2,IFITM3,IFNAR2,IL18BP,IL6,ITGB7,METTL7B,MYD88,NLRC5,NOD1,NUP93,OAS3,OGFR,PARP12,PIM1,PLSCR1,PNP,PSMA2,... |
| down_DMS53 | 2701 | HALLMARK_APICAL_JUNCTION | 33/200 | ACTB,ACTN1,ADAM9,ATP1A3,CD274,CDH8,CERCAM,CLDN4,CLDN7,CNN2,CRB3,DHX16,EPB41L2,EXOC4,ICAM5,ITGB4,JUP,LAMB3,MAP3K20,MPZL1,MYH9,MYL9,NECTIN4,NFASC,PCDH1,PFN1,PIK3CB,PTEN,RRAS,SKAP2... |
| down_DMS53 | 2701 | HALLMARK_E2F_TARGETS | 131/200 | ANP32E,ASF1B,ATAD2,AURKA,AURKB,BARD1,BIRC5,BRCA1,BUB1B,CBX5,CCNB2,CCNE1,CDC20,CDC25A,CDC25B,CDCA3,CDCA8,CDK1,CDKN2A,CDKN2C,CDKN3,CENPE,CENPM,CHEK1,CHEK2,CKS2,CSE1L,CTPS1,DCTPP1,... |
| down_DMS53 | 2701 | CLDN4 | 1/1 | CLDN4 |
| down_DMS53 | 2701 | TACSTD2 | 1/1 | TACSTD2 |
| direct_232 | 232 | NHEJ_CORE | 1/13 | DCLRE1C |
| direct_232 | 232 | STING_CORE | 0/11 |  |
| direct_232 | 232 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 0/97 |  |
| direct_232 | 232 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1/200 | PTPN1 |
| direct_232 | 232 | HALLMARK_APICAL_JUNCTION | 2/200 | CLDN4,SKAP2 |
| direct_232 | 232 | HALLMARK_E2F_TARGETS | 11/200 | ASF1B,ATAD2,AURKA,CBX5,CDC25B,CDCA3,KIF18B,POLD1,RACGAP1,RPA3,SYNCRIP |
| direct_232 | 232 | CLDN4 | 1/1 | CLDN4 |
| direct_232 | 232 | TACSTD2 | 0/1 |  |
| near_ge2 | 3288 | NHEJ_CORE | 5/13 | DCLRE1C,PNKP,POLM,PRKDC,XRCC4 |
| near_ge2 | 3288 | STING_CORE | 2/11 | STING1,TREX1 |
| near_ge2 | 3288 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 11/97 | CASP8,IRF1,LGALS3BP,MVB12A,PSMB8,PSMB9,PSME2,RNF31,TAP1,TENT5A,TRIM14 |
| near_ge2 | 3288 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 27/200 | BPGM,CASP8,CDKN1A,HLA-DMA,IRF1,LGALS3BP,MTHFD2,NAMPT,NFKBIA,PELI1,PSMB10,PSMB8,PSMB9,PSME2,PTPN1,PTPN6,RAPGEF6,RBCK1,RIPK1,RNF31,SOD2,SPPL2A,STAT3,TAP1,TRIM14,UPP1,ZNFX1 |
| near_ge2 | 3288 | HALLMARK_APICAL_JUNCTION | 28/200 | ACTB,ACTG1,ADAM15,AKT2,ARPC2,CAP1,CLDN4,CLDN7,CRB3,GNAI2,GTF2F1,JUP,KCNH2,LIMA1,MDK,MYH9,MYL12B,PFN1,PIK3R3,PTK2,RASA1,RRAS,SHC1,SKAP2,SYMPK,THBS3,TIAL1,VASP |
| near_ge2 | 3288 | HALLMARK_E2F_TARGETS | 62/200 | ASF1B,ATAD2,AURKA,BRCA1,BRCA2,CBX5,CDC25B,CDCA3,CDCA8,CDKN1A,CKS1B,CKS2,CTCF,CTPS1,DCLRE1B,DCTPP1,EIF2S1,GSPT1,HNRNPD,ING3,KIF18B,LIG1,LMNB1,MCM4,MCM7,MLH1,MRE11,MTHFD2,NASP,NBN... |
| near_ge2 | 3288 | CLDN4 | 1/1 | CLDN4 |
| near_ge2 | 3288 | TACSTD2 | 1/1 | TACSTD2 |
| down_ge2 | 1121 | NHEJ_CORE | 3/13 | APTX,DCLRE1C,XRCC6 |
| down_ge2 | 1121 | STING_CORE | 2/11 | ENPP1,IKBKE |
| down_ge2 | 1121 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 5/97 | CMTR1,HLA-C,IFITM2,SAMD9,TRIM25 |
| down_ge2 | 1121 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 10/200 | CASP3,CASP4,CMTR1,HLA-A,HLA-B,IFITM2,IFNAR2,PTPN1,TOR1B,TRIM25 |
| down_ge2 | 1121 | HALLMARK_APICAL_JUNCTION | 11/200 | AKT3,CDH8,CLDN4,EXOC4,MAP3K20,MPZL1,NFASC,PIK3CB,SKAP2,SYK,VCL |
| down_ge2 | 1121 | HALLMARK_E2F_TARGETS | 49/200 | ASF1B,ATAD2,AURKA,AURKB,BIRC5,CBX5,CCNE1,CDC25A,CDC25B,CDCA3,CDKN2C,CENPE,CENPM,DEK,DEPDC1,DIAPH3,DSCC1,DUT,E2F8,EED,ESPL1,GINS3,KIF18B,KIF2C,KIF4A,MAD2L1,MCM2,MCM5,MCM6,MELK,MK... |
| down_ge2 | 1121 | CLDN4 | 1/1 | CLDN4 |
| down_ge2 | 1121 | TACSTD2 | 0/1 |  |

## What this does not say

- No public lung ELF3 ChIP, so a LUAD binding claim at NHEJ/STING/IFN promoters is not testable. DMS53 CUT&Tag is SCLC.
- NEC siELF3 RNA-seq is only the downregulated-gene lists, so those tables cannot show interferon induction.
- Mouse TEC and the synovial TNF contrasts are cytokine-on. They are not basal ELF3-loss tests.
- Co-expression of ELF3 with CLDN4 (earlier public catalogs) is not binding and is not re-tested here.
- This is not an ICI or TROP2-ADC outcome analysis.

## Datasets

- **A549_shELF3** (PASS, primary): GSE137479 A549 LUAD stable shELF3 clones vs empty-vector clones, n=3. HuGene 2.0 ST RMA. The only public lung ELF3-loss transcriptome in this set.
- **HBDEC2_KO** (PASS, primary): GSE148105 HBDEC2 biliary epithelial ELF3 CRISPR KO vs WT, n=3. Agilent 8x60K. Not lung. Paired with GSE156165 ELF3 ChIP in the same line. Platform symbols TMEM173/MB21D1/C9orf142 are renamed to STING1/CGAS/PAXX.
- **HBDEC2_KD** (PASS, primary): GSE148106 HBDEC2 ELF3 miRNA KD vs control, n=3. Same platform and line as the KO, independent reagent.
- **JEG3_KO** (PASS, primary): GSE241792 JEG-3 trophoblast CRISPR sgELF3 vs WT, n=3. Salmon counts, log2(CPM+1). Not lung.
- **CRPC_shELF3** (PASS, supporting_n2): GSE303076 CRPC shELF3 vs PLKO, n=2. Author log2FoldChange is stored for ELF3/CLDN4 checks. Not in the 5-dataset verdict.
- **synov_siELF3_basal** (FAIL, primary): GSE129487 RA synovial fibroblasts, siELF3 vs siCtrl, unstimulated (t=0). Four donors, technical duplicates averaged. Not lung. Donors used: RA5,RA6,RA7,RA8.
- **synov_siELF3_TNF16** (FAIL, stimulated): GSE129487 siELF3 vs siCtrl at 16 h TNF. Same donors. Donors used: RA5,RA6,RA7,RA8.
- **synov_siELF3_TNF_IL17_16** (PASS, stimulated): GSE129487 siELF3 vs siCtrl at 16 h TNF+IL-17A. Same donors. Donors used: RA6,RA7,RA8.
- **mTEC_KO_TNFa_IFNg** (PASS, stimulated): GSE319933 primary mouse renal tubular cells, AdvCre Elf3 deletion vs control, all samples already on TNFα+IFNγ. Tests whether Elf3 loss changes programs during cytokine stimulation, not basal loss.
- **GSE190618** ELF3 CUT&Tag peaks: ECC4, A99, DMS53. Binding only.
- **GSE156165** HBDEC2 ELF3 ChIP peaks.
- **ENCODE ENCFF080FAU** HepG2 ELF3 IDR peaks.
- **ChIP-Atlas** hg38 ELF3 5 kb target matrix.
- **GSE319848** mouse TEC ELF3 CUT&RUN under TNFα+IFNγ (stimulated binding, not basal).

## Methods

Microarray: A549 HuGene 2.0 ST is GEO RMA (already log2), probes collapsed to HGNC symbols with Bioconductor `hugene20sttranscriptcluster.db` (unambiguous probes only) and the highest-mean probe kept. HBDEC2 Agilent intensities are log2(signal+1), same collapse using the platform gene-symbol column. Retired symbols (TMEM173, MB21D1, C9orf142, and other HGNC previous/alias symbols that are not themselves current symbols) are renamed to the current symbol before the tests. RNA-seq counts are log2(CPM+1). Mouse TPM is log2(TPM+1). Synovial technical duplicates are averaged inside donor, then a paired t-test is taken across donors. Detection floor for RNA-seq / TPM is mean abundance ≥ 1 (CPM or TPM); synovial detection is mean log2(CPM+1) ≥ 1 on either side. Set tests are two-sided Mann–Whitney of set log2FC versus other detected genes, plus a two-sided Wilcoxon signed-rank versus 0. Gene FDR is BH inside detected genes. Binding: ChIP-Atlas score > 0, or any RefSeq TSS ±2 kb overlapping a called peak. Fisher exact tests are two-sided against all genes in that universe.

Reproduce: `pip install -r scripts/elf3_nhej_sting_ifn/requirements.txt` then `python scripts/elf3_nhej_sting_ifn/analyze.py`.


# GSE22493 — symbol remap, then NHEJ / STING / IFN vs CLDN4-overexpression control

GSE22493 only. SKOV-3-IP-Luc ovarian cells, three two-color Operon arrays (GSM558700–702). Channel 1 is Cy3 **CLDN4 overexpression (control)**. Channel 2 is Cy5 **CLDN4 lentiviral siRNA**. The deposited series-matrix `VALUE` is log2(knockdown / control). This is not lung, and it is not knockdown versus scramble or wild type.

> **Bottom line.** After remapping GPL10555 tokens onto current HGNC symbols, the three panels do not line up as NHEJ-down, STING-up, IFN-up. NHEJ (10 genes scored) is mixed: median of per-gene mean log2 **+0.15**, Mann-Whitney p = **0.98** versus background. STING (7/9) sits on zero (**+0.04**, p = **0.62**); **STING1 and TBK1 are not on the array**, and the sensor **CGAS is down**. IFN (63/73) is not up: **28 up / 35 down**, median of means **−0.16**, p = **0.47**. No NHEJ or STING gene has a within-panel q < 0.05. CLDN4 itself is not cleanly knocked down on the array.
>
> **一句话。** 把 GPL10555 的旧符号映到现行 HGNC 之后，这三张卵巢芯片上看不到「NHEJ 下降、STING 上升、IFN 上升」的配套变化。NHEJ 10 个可评分基因中位均值 log2 为 **+0.15**（对背景 MW p = **0.98**）。STING 为 **+0.04**（p = **0.62**），而且 **STING1 和 TBK1 不在芯片上**，传感器 **CGAS 向下**。IFN 63/73 为 **28 上 / 35 下**，中位均值 **−0.16**（p = **0.47**）。CLDN4 敲低在阵列上不能干净确认。

---

## What was remapped

GPL10555 (BWH Human Release 3.0, Operon 60-mer) stores 2010-era ORF tokens. A current-symbol match misses genes that the array still calls by a retired name. Tokens were resolved against NCBI `Homo_sapiens.gene_info` (protein-coding, nomenclature status O; md5 `de642eac2f7a24354e080c249a44f609`):

| Rule | Probes |
|---|---:|
| Token is a current symbol: keep | 15,529 |
| Token is an unambiguous synonym of one current gene: remap | 4,974 |
| Token matches several genes and is not itself current: drop | 407 |
| Token not in the NCBI table: drop | 4,192 |
| No ORF and no `SYMBOL--` prefix: drop | 11,186 |

Resolved genes: 13,556. Background for the set test: 12,814 genes with a value on at least two arrays (median of per-gene mean log2 = **−0.027**).

A current symbol that is also a retired synonym of a different gene stays with the current symbol. On this panel that affects TAP1 (also a synonym of SEC14L2), TAP2 (SEC14L3), and ZBP1 (IGF2BP1). The GPL descriptions say ABC transporter and Z-DNA binding protein, so those probes stay TAP1, TAP2, and ZBP1 (`symbol_collisions_kept.tsv`).

Panel genes that exist only under a retired token:

| Current | Platform token | Where it matters |
|---|---|---|
| XRCC6 (Ku70) | G22P1 | NHEJ; scored, median −1.56 |
| APLF | C2orf13 | NHEJ; one array only, not scored |
| CGAS | C6orf150 | STING; scored, median −0.49 |
| H2AX | H2AFX | break marker; scored, median +0.40 |
| RIGI | DDX58 | IFN; scored, median −1.45 |
| ISG15 | G1P2 | IFN |
| IFI6 | G1P3 | IFN |
| IFI44L | C1orf29 | IFN |
| XAF1 | HSXIAPAF1 | IFN |
| EIF2AK2 | PRKR | IFN |
| NLRC5 | NOD27 | IFN |
| IL32 | NK4 | IFN |
| RSAD2 | cig5 | annotated, all three values missing |

The earlier GSE22493 IFN slice remapped only G1P2 and G1P3. It was a different question (IFN / MHC-I / APM). This run does not replace that table.

Free-text probes with an empty ORF were not folded into gene medians. Probe 8697 (“ARTEMIS PROTEIN”) is negative on all three arrays (−1.79 / −3.84 / −1.25), the same direction as the DCLRE1C ORF median. Probe 25010 (Ku80 / Ku86 text) is −4.32 / +1.18 / −2.84 and was left out.

## Panels (locked before the scores)

- **NHEJ (13):** XRCC6, XRCC5, PRKDC, XRCC4, LIG4, NHEJ1, PAXX, DCLRE1C, APLF, PNKP, APTX, POLL, POLM. Classical end-joining. Not MRN, not alt-EJ (PARP1, LIG3, XRCC1), not 53BP1-Shieldin.
- **STING (9):** CGAS, STING1, TBK1, IKBKE, IRF3, DDX41, IFI16, TREX1, ENPP1. Signaling axis, including two negative regulators. The IFN program is separate, so a STING score is not an ISG score.
- **IFN (73):** the same IFN_IMMUNE list as the earlier slice, with DDX58 entered as the current symbol RIGI. IFNB1 and CXCL10 stay here as STING transcriptional outputs.
- **Break markers (adjunct, not in the three-panel score):** H2AX, TP53BP1.

The pattern this wave is built to see is NHEJ down, STING up, IFN up on the deposited log2 ratio. The test against background is two-sided Mann-Whitney on the per-gene mean log2. One-sample t-tests versus 0 at n = 2–3 are descriptive. BH-FDR is within panel. Missing genes stay blank.

## CLDN4 check

Probe 17169. Deposited log2 is missing, **−1.74**, **−0.71**. ScanArray background-subtracted Cy5 is negative on GSM558700, and GSM558701 raw log2(Cy5/Cy3) is **+1.22**, opposite the deposited value. The control is overexpression. Knockdown is not confirmed on the array. Everything below sits on that design.

## NHEJ

Absent: **NHEJ1, PAXX**. APLF (C2orf13) has a value on GSM558702 only (+0.87) and is not scored.

| Gene | Map | Three arrays | Median | Direction | p |
|---|---|---|---:|---|---:|
| XRCC6 | G22P1 → XRCC6 | −1.74 / +0.07 / −1.56 | −1.56 | DOWN | 0.20 |
| XRCC5 | current | −1.18 / +0.42 / −0.84 | −0.84 | DOWN | 0.39 |
| DCLRE1C | current | −0.77 / −1.98 / −4.19 | −1.98 | DOWN | 0.15 |
| POLM | current | −0.99 / +1.10 / −1.56 | −0.99 | DOWN | 0.61 |
| LIG4 | current | −0.64 / +2.46 / 0 | 0 | ZERO | 0.59 |
| PNKP | current | −0.62 / +0.14 / +2.24 | +0.14 | UP | 0.56 |
| POLL | current | +1.66 / −0.01 / +0.15 | +0.15 | UP | 0.38 |
| XRCC4 | current | −2.94 / +0.26 / +2.49 | +0.26 | UP | 0.97 |
| PRKDC | current | +0.43 / −0.13 / +0.77 | +0.43 | UP | 0.31 |
| APTX | current | +1.44 / +0.73 / −0.60 | +0.73 | UP | 0.47 |

5 up / 4 down / 1 zero. Median of the per-gene means **+0.15** versus background **−0.027**, MW p = **0.98**. No within-panel q is below 0.68. Ku70 is down only because the remap recovered G22P1, and one of its three arrays is positive. Artemis is the most negative gene and still p = 0.15. This is not a coordinated NHEJ drop.

ScanArray Cy5/Cy3 (sensitivity, not the primary ratio): the four deposited downs stay down, and the five deposited ups (PRKDC, XRCC4, PNKP, APTX, POLL) flip to down. LIG4, deposited zero, is down on all three raw arrays. Raw intensities do not rescue an NHEJ-up reading. They also do not turn this into a clean NHEJ-down result: the arrays disagree with each other, and the CLDN4 channel is not a confirmed knockdown.

## STING

**STING1 and TBK1 are absent** (no TMEM173, MPYS, MITA, NAK, or T2K token either). The adaptor and the main kinase cannot be scored. cGAS is present only as C6orf150.

| Gene | Role | Arrays | Median | Direction | p |
|---|---|---|---:|---|---:|
| CGAS | sensor (C6orf150) | −0.64 / NA / −0.34 | −0.49 | DOWN | 0.19 |
| ENPP1 | negative regulator | −3.18 / +0.74 / −1.79 | −1.79 | DOWN | 0.34 |
| IFI16 | sensor | −0.04 / +1.12 / −1.03 | −0.04 | DOWN | 0.98 |
| DDX41 | sensor | +0.25 / +0.18 / −0.32 | +0.18 | UP | 0.86 |
| TREX1 | negative regulator | −0.51 / NA / +1.71 | +0.60 | UP | 0.69 |
| IKBKE | kinase | −0.14 / +1.09 / +0.77 | +0.77 | UP | 0.26 |
| IRF3 | transcription factor | +1.88 / −1.18 / +2.50 | +1.88 | UP | 0.45 |

4 up / 3 down. Median of means **+0.035**, MW p = **0.62**. The positive sign is a few hundredths of a log2 and is not an up-shift. CGAS, the sensor that is actually measured, is down on both arrays that have a value (ScanArray agrees, median −1.72). IRF3 is the most consistently positive gene, including on raw Cy5/Cy3 (all three arrays up), and its deposited p is 0.45. IKBKE and DDX41 flip to down on ScanArray. TREX1 up and ENPP1 down are opposite directions for two negative regulators. Without STING1 and TBK1 there is no STING-axis call on this platform.

## IFN

66/73 genes are on the platform after remap; 63 have at least two arrays. Median direction **28 up / 35 down**. Median of per-gene means **−0.157** versus background **−0.027**, MW p = **0.47**. Not an IFN up-shift.

The only nominal p < 0.05 in the whole IFN list is JAK2 (median +1.80, p = 0.033, within-panel q = 0.95). That one gene is not a panel result. RIGI, recovered from DDX58, is down on both available arrays (median −1.45, p = 0.32). IFNB1 median is +0.72 with the two arrays at +4.49 and −3.06 (p = 0.88). CXCL10 is +6.50 and −0.45 (p = 0.54). The proximal STING outputs are not a consistent rise.

Still absent after remap: MX2, IFITM1, IFITM3, DDX60, CMPK2, IRF9, IFNL1. RSAD2 (cig5) is annotated and has no deposited values. CXCL9 and CCL5 have one array each and are not scored.

## Break markers (adjunct)

H2AX (H2AFX) is positive on all three deposited arrays (median +0.40, p = 0.28). TP53BP1 median is −0.42 (p = 0.52), and the ScanArray median flips positive. A single histone gene at p = 0.28 is not a DNA-break signature, and it is not part of the NHEJ enzyme score.

## What this accession can support

On these three arrays, after the symbol remap, NHEJ is mixed, the measurable STING genes sit on zero with the hub missing and CGAS down, and IFN is not up. That is not the coordinated NHEJ-down / STING-up / IFN-up pattern. It is also not a clean demonstration that CLDN4 siRNA closes those programs: the control is CLDN4 overexpression, the CLDN4 probe contradicts itself, n = 3, and the platform is a 2010 custom array.

Do not export this to lung, to ICI, or to a private knockdown series.

## Reproduce

```bash
pip install numpy scipy statsmodels matplotlib
python3 scripts/nhej_sting_gse22493/download_data.py
python3 scripts/nhej_sting_gse22493/run_analysis.py
```

Outputs: `results/nhej_sting_gse22493/`. Series-matrix md5 `8f8a07cb8f3f3a90396c3feb1eb3ea0f`. Family SOFT md5 `5d88cef4c873370970e6ba6be82d7b48`.

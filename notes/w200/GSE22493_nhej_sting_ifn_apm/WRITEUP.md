# GSE22493 ovarian CLDN4：NHEJ / STING / IFN / APM

# GSE22493 ovarian CLDN4: NHEJ / STING / IFN / APM

**STATUS: FINAL for GSE22493.** One method sweep is in `method_sweep.tsv`. The primary deposited analysis and the three-array replicate test stay null. Alternate summaries that cross p < 0.05 disagree with that replicate test, and the STING set median flips sign between the deposited ratio and ScanArray. This accession stops here.

Additive public slice of **GSE22493 only** (SKOV-3-IP-Luc; lentiviral CLDN4 siRNA vs CLDN4-overexpression control; three two-color Operon arrays). IFN and APM use the same gene lists as the earlier C4 slice of this accession. NHEJ and STING are new pre-specified panels. This file does not revise CosMx, concordant-4, GSE137244, TCGA, or TISMO.

> **一句话 / FINAL**
>
> GSE22493 到此为止。主分析（作者沉积 log2，探针中位数再对阵列取均值）上 NHEJ / STING / IFN / APM 都没有 p < 0.05 的集合检验，也没有面板内 q < 0.05 的基因。以三张芯片为重复单位的 t 检验同样是空的（四个主面板 p = 0.28–0.75）。换汇总方式后出现的名义 p < 0.05 与这个重复单位检验不一致；STING 集合中位数在沉积值上是 +0.58，在 ScanArray 上是 −0.79。唯一 q < 0.05 的基因是 ScanArray 上的 HLA-B（q = 0.033），而沉积值在 GSM558701 上是 +2.33，与原始强度 −2.02 反号。
>
> **Bottom line / FINAL**
>
> GSE22493 stops here. On the primary analysis (deposited log2, probe median, then the mean across arrays) the NHEJ, STING, IFN, and APM set tests are all p ≥ 0.16, and no panel gene has q < 0.05. The three-array replicate test is also null (primary-panel p = 0.28–0.75). Nominal p < 0.05 values from other summaries disagree with that replicate test. The STING set median is +0.58 on the deposited ratios and −0.79 on ScanArray. The only gene with q < 0.05 is HLA-B on ScanArray (q = 0.033); its deposited value on GSM558701 is +2.33 against a ScanArray value of −2.02.

---

## Design

| Item | Fact |
|---|---|
| Accession | [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) (no companion paper; submitter Zhijian Gao, BWH) |
| Platform | [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555) BWH Human Release 3.0, Operon 60-mer, ~36k probes |
| Cells | SKOV-3-IP-Luc, human ovarian cancer. This is an ovarian line. |
| Samples | GSM558700, GSM558701, GSM558702 |
| Channels | Cy3 = CLDN4 overexpression (control); Cy5 = CLDN4 lentiviral siRNA |
| Deposited value | Series-matrix `VALUE`. The numbers are symmetric about 0 (background median −0.0265), so this slice reads them as log2(knockdown/control), the same reading as the C4 slice. GEO text says “normalized sample to control ratios” and does not print the word log2. |
| What was not done | No scramble or WT arm. GEO text mentions dye-swap; all three GSMs carry the same Cy3/Cy5 assignment. |

The contrast is overexpression versus knockdown. A change in a panel gene is a difference between those two arms.

## Panels (fixed before looking at log-ratios)

Primary:

- **NHEJ (8):** XRCC6 (Ku70), XRCC5 (Ku80), PRKDC (DNA-PKcs), LIG4, XRCC4, NHEJ1 (XLF), DCLRE1C (Artemis), PAXX
- **STING axis (5):** CGAS, STING1, TBK1, IRF3, IKBKE
- **IFN (73):** the C4 IFN_IMMUNE list (type-I ISGs, IFN signaling, a short chemokine set)
- **APM (16):** HLA-A/B/C, B2M, NLRC5, TAP1, TAP2, TAPBP, PSMB8, PSMB9, PSMB10, ERAP1, ERAP2, CALR, CANX, PDIA3

Secondary, kept out of the primary set medians:

- **NHEJ accessory (7):** PNKP, APLF, APTX, POLL, POLM, DNTT, WRN
- **STING regulators (8):** TREX1, ENPP1, IFI16, DDX41, SAMHD1, RNASEH2A, RNASEH2B, RNASEH2C

Mapping rules, same conservatism as the C4 slice:

- A probe is used when the ORF field is filled, or when DESCRIPTION starts with `SYMBOL--`.
- Explicit aliases: `G22P1 → XRCC6`, `C6orf150 → CGAS`, `G1P2 → ISG15`, `G1P3 → IFI6`.
- Near-misses left unmapped: TTBK1 (tau-tubulin kinase, not TBK1), KUB3 (not Ku70), DNTTIP1 (not DNTT), WRNIP1 (not WRN).
- Empty-ORF free text was recorded and excluded: probe 25010 (Ku80/Ku86 description) and probe 8697 (Artemis description). XRCC5 and DCLRE1C already have ORF probes.

## Method

1. Gene value on an array = median of mapped probes. Gene median and gene mean are then taken across arrays that have a value.
2. Set test: gene mean log2 for genes with a value on ≥2 arrays, versus all other mapped genes with ≥2 arrays. Two-sided Mann–Whitney. A Wilcoxon signed-rank of those gene means against 0 is reported when at least six genes are in the set and the signs are not all the same. One-sample t tests and within-panel BH-FDR are descriptive at n = 2–3.
3. A set whose median sits within 0.05 log2 of the background median is labeled FLAT.
4. Sensitivity: ScanArray log2 of (Cy5 Median−B) / (Cy3 Median−B), both channels required to be > 0. This is not the author’s deposited ratio.

## CLDN4 on the array

Probe 17169, from `cldn4_diagnostic.tsv`:

| Array | Deposited log2 | ScanArray Cy3 / Cy5 (Median−B) | ScanArray log2 |
|---|---:|---|---:|
| GSM558700 | missing | 89 / −39 | unusable (Cy5 background-subtracted median is negative) |
| GSM558701 | −1.74 | 251 / 584 | +1.22 |
| GSM558702 | −0.71 | 1031 / 330 | −1.64 |

Deposited median on the two arrays with a value is −1.23. GSM558701 deposited and ScanArray disagree in sign. The control arm is CLDN4 overexpression. The array does not confirm a clean knockdown.

## NHEJ

`nhej_genes.tsv`. Six of eight core genes are measured. NHEJ1 and PAXX are absent from GPL10555.

| Gene | Deposited log2 (700 / 701 / 702) | Median | Mean | p | ScanArray median | Signs agree |
|---|---|---:|---:|---:|---:|---|
| XRCC6 (G22P1) | −1.74 / +0.07 / −1.56 | −1.56 | −1.07 | 0.20 | −2.16 | yes, down |
| XRCC5 | −1.18 / +0.42 / −0.84 | −0.84 | −0.53 | 0.39 | −1.28 | yes, down |
| PRKDC | +0.43 / −0.13 / +0.77 | +0.43 | +0.36 | 0.31 | −0.48 | no |
| LIG4 | −0.64 / +2.46 / 0.00 | 0.00 | +0.60 | 0.59 | −1.63 | deposited median is zero |
| XRCC4 | −2.94 / +0.26 / +2.49 | +0.26 | −0.06 | 0.97 | −0.38 | no |
| DCLRE1C | −0.77 / −1.98 / −4.19 | −1.98 | −2.31 | 0.15 | −1.75 | yes, down |

Set test (`geneset_stats.tsv`): 6 genes, gene-mean direction 2 up / 4 down, median of gene means −0.298, background −0.027, delta −0.271, Mann–Whitney p = 0.36, Wilcoxon vs 0 p = 0.44. Within-panel q values are all > 0.5.

Ku70, Ku80, and Artemis are the three core genes whose deposited and ScanArray medians are both negative. Their nominal p values are 0.20, 0.39, and 0.15. PRKDC’s two ORF probes disagree in sign (probe 20613 median +0.83; probe 32709 median −0.56), and the gene-level ScanArray median is negative. XRCC4 swings from −2.94 to +2.49 across the three arrays. LIG4’s ScanArray values are negative on all three arrays (−1.63 / −2.96 / −1.13) while the deposited GSM558701 value is +2.46, so the deposited median of 0 is an unstable summary.

Accessory factors are a separate row. APLF is absent. DNTT has a value on only one array and is out of the set test. The five genes in the set have median of means +0.584, 4 up / 1 down, Mann–Whitney p = 0.14. ScanArray reverses PNKP, APTX, POLL, and WRN. The only accessory gene with agreeing signs is POLM, and it is down (deposited median −0.99, ScanArray −1.90).

## STING

`sting_genes.tsv`. STING1 (TMEM173 / MPYS / MITA / ERIS) and TBK1 are not annotated on GPL10555. TTBK1 is a different kinase and was not substituted. The axis that can be scored is CGAS, IRF3, and IKBKE.

| Gene | Deposited log2 | Median | Mean | p | ScanArray median | Signs agree |
|---|---|---:|---:|---:|---:|---|
| CGAS (C6orf150) | −0.64 / missing / −0.34 | −0.49 | −0.49 | 0.19 | −1.72 | yes, down |
| IRF3 | +1.88 / −1.18 / +2.50 | +1.88 | +1.06 | 0.45 | +2.11 | yes, up |
| IKBKE | −0.14 / +1.09 / +0.77 | +0.77 | +0.58 | 0.26 | −0.59 | no |

Set test: 3/5 genes, 2 up / 1 down, median of gene means +0.575, delta vs background +0.602, Mann–Whitney p = 0.36. No within-panel q < 0.05. The positive set median is the middle of three discordant genes, and the kinase that is present (IKBKE) flips on raw intensities. CGAS, the sensor that is present, is down on both readouts (two arrays).

Regulators (`sting_reg_genes.tsv`): 6/8 measured on ≥2 arrays (RNASEH2B and RNASEH2C absent). Median of gene means −0.025 versus background −0.027 (delta +0.0015), Mann–Whitney p = 0.50, Wilcoxon p = 0.56. Labeled FLAT. Gene means are 3 up / 3 down.

## IFN and APM

Same lists and the same collapse as the C4 slice. Numbers match that slice.

| Set | ≥2 arrays | Gene-mean up / down | Median of means | Mann–Whitney p | Shift |
|---|---:|---|---:|---:|---|
| IFN | 58/73 | 25 / 33 | −0.232 | 0.29 | DOWN |
| APM | 12/16 | 3 / 9 | −0.454 | 0.16 | DOWN |

Background median −0.027, 15,853 genes. CXCL9 and CCL5 have a value on only one array and stay out of the set test, as in the C4 slice. NLRC5, ERAP1, ERAP2, and PDIA3 are absent. APM medians are negative for HLA-A/B/C, B2M, TAP1, TAP2, TAPBP, PSMB10, and CANX. Positive APM medians are PSMB8, PSMB9, and CALR. PSMB8’s ScanArray median is negative, so 11/12 measured APM genes agree in sign across deposited and ScanArray, and nine of those agreements are down.

The C4 priority ISGs that look up on the deposited median (OAS2, MX1, ISG15) are down on ScanArray. IFIT1 is down on both. That sensitivity result is unchanged.

## What this accession can carry

On GSE22493, the pre-specified NHEJ, STING, IFN, and APM panels are null at the gene-FDR and set-test thresholds used here. Point estimates: NHEJ, IFN, and APM sit below the background median; the incomplete STING axis sits above it; STING regulators sit on the background median. The Ku70 / Ku80 / Artemis deposited medians are negative on both the deposited ratio and the ScanArray ratio, with nominal p from 0.15 to 0.39. Those three genes are a description of this array, and they do not survive multiplicity or the design limits below.

Limits that stay attached to every number:

1. Ovarian SKOV-3-IP-Luc. This slice is not a lung result and is not an ICI result.
2. The control is CLDN4 overexpression. There is no scramble or WT arm.
3. CLDN4 itself is missing on GSM558700, and GSM558701 deposited versus ScanArray disagree in sign.
4. n = 3 arrays from 2010. One-sample t tests are descriptive.
5. STING1 and TBK1 are absent, so a cGAS–STING axis test is not available on this platform.
6. NHEJ1 and PAXX are absent. PRKDC probes disagree with each other.
7. IFN and APM here reproduce the earlier C4 table. They do not add a second cohort.

## Method sweep (then stop)

`method_sweep.py` reruns the same gene lists under seven summaries. The gene lists do not change. `PRE_VALUE` in the SOFT sample tables is the linear ratio: |VALUE − log2(PRE_VALUE)| has maximum 5.0×10⁻⁵ on every array (`value_equals_log2_pre.tsv`), which is the rounding of VALUE. The log2 reading is the deposited scale.

Primary deposited result, reproduced by the sweep:

| Panel | Median of gene means | Mann–Whitney p | Min gene q |
|---|---:|---:|---:|
| NHEJ | −0.298 | 0.36 | 0.58 |
| STING | +0.575 | 0.36 | 0.39 |
| IFN | −0.232 | 0.29 | 0.90 |
| APM | −0.454 | 0.16 | 0.79 |

Treating each array as the replicate (mean of the panel on that array, one-sample t, n = 3), from `method_sweep_array_level.tsv`:

| Panel | Array means (700 / 701 / 702) | p | p after subtracting the array-wide median |
|---|---|---:|---:|
| NHEJ | −1.14 / +0.18 / −0.55 | 0.32 | 0.31 |
| STING | +0.37 / −0.05 / +0.98 | 0.28 | 0.31 |
| IFN | −0.47 / +0.58 / −0.50 | 0.75 | 0.78 |
| APM | −0.90 / +0.91 / −0.79 | 0.70 | 0.71 |

GSM558701 is the array that sits on the other side of zero for NHEJ, IFN, and APM. Transcriptome-wide medians on the three arrays are −0.059, +0.070, and −0.118.

Summaries that cross 0.05 somewhere, all from `method_sweep.tsv`:

| Summary | What crosses 0.05 | What stays null next to it |
|---|---|---|
| Deposited, gene median instead of gene mean | APM set p = 0.028 (median −0.91, still down) | APM min gene q = 0.79. Array-level APM p = 0.70 |
| ScanArray log2(Cy5/Cy3) | IFN set p = 0.036 (median −0.83). HLA-B q = 0.033 | IFN min gene q = 0.053. APM set p = 0.54. Array-level IFN p = 0.75 |
| Vendor Ch2 Log Ratio, \|value\| > 15 dropped | NHEJ set p = 0.043. IFN set p = 0.0098. HLA-B q = 0.034 | No NHEJ or IFN gene has q < 0.05 (IFN min q = 0.085). APM set p = 0.39. Array-level tests unchanged |

ScanArray and the capped vendor log ratio are almost the same column (median absolute difference about 0.003 on the bulk of spots). The vendor file keeps a few spots the background-subtracted log2 drops, which is why the IFN and NHEJ set p-values move. Seven methods times four primary panels is 28 set tests. The smallest of those p-values is 0.0098.

HLA-B is the only gene with within-panel q < 0.05 in the whole sweep, and only on the raw channel ratios. Deposited HLA-B is −1.56 / +2.33 / −2.12. ScanArray HLA-B is −1.68 / −2.02 / −1.90. The q = 0.033 call uses the raw ratios and conflicts with the deposited ratio on GSM558701.

STING set median by summary: deposited primary +0.58 (p = 0.36), ScanArray −0.79 (p = 0.76). Same genes, opposite sign. STING1 and TBK1 are absent in every method.

Adding the two empty-ORF spots (probe 25010 onto XRCC5, probe 8697 onto DCLRE1C) moves the NHEJ median from −0.30 to −0.57 and the set p from 0.36 to 0.24. Still null.

Secondary STING-regulator array means are −0.49 / −0.34 / −0.27, one-sample p = 0.029. After subtracting each array’s transcriptome-wide median the p is 0.065. The regulator set versus background remains flat (median −0.025, Mann–Whitney p = 0.50). That row is not a cGAS–STING result, and it is not a reason to keep scoring this accession.

`key_stats.json` field `accession_status` is `FINAL_DISCORDANT`.

## Rerun

```bash
python3 scripts/w200/GSE22493_nhej_sting_ifn_apm/download_data.py
python3 scripts/w200/GSE22493_nhej_sting_ifn_apm/run_analysis.py
python3 scripts/w200/GSE22493_nhej_sting_ifn_apm/method_sweep.py
```

Cache directory: `W200_GSE22493_PANEL_DATA` (default `/tmp/gse22493`).

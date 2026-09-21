# GSE22493: CLDN4 siRNA vs CLDN4 overexpression after HGNC symbol repair

## 中文

GSE22493（SKOV-3-IP-Luc，CLDN4 siRNA / CLDN4 过表达，3 张双色芯片）先前把 GPL10555 的 2006 年 ORF 字符串当成现行基因符号，只手工改了 G1P2→ISG15、G1P3→IFI6。NOD27、GRP58、ARTS-1、LRAP、cig5、PRKR、C1orf29、FLJ20035 等因此进不了 IFN/APM 检验。本分析用 HGNC complete set（批准符号，否则唯一的曾用名，否则唯一的别名）重映射后再测用户 KD 方向：CLDN4 降低后 IFN 与 APM 应升高，其它紧密连接基因应降低。

符号补上之后，结论没有翻向。以中位 log2(siRNA/OE) 为准，优先级 6 基因 3 升 3 降（中位 −0.744），APM 5/15 升（中位 −0.560），IFN 集 29/64 升（中位 −0.101）。相对背景基因的单侧 Mann-Whitney 都不支持“升高”（p = 0.881、0.973、0.888；四组 BH q = 0.973）。CLDN4 本身中位 −1.225（2 张芯片），方向与标记一致。去掉 CLDN4 后的紧密连接核心 8/14 为负（中位 −0.247），但和背景无差别（p = 0.516），基因均值符号相反。这是卵巢 siRNA 对过表达，不是肺，也不是亲本对照；三张芯片彼此 Spearman 只有 −0.24 到 +0.24。它不支持“CLDN4 敲低打开 IFN/APM”，也不拿来否定别的公共或私有 KD。

最后一轮又查了原始 ScanArray、染料交换符号、三张芯片留一、探针级、效应量分位数、ssGSEA 和 n=3 的 GSVA 式评分。没有任何一个方向正确的实验级处理能让 IFN 和 APM 同时升高且单侧 p<0.05。停在 FINAL discordant。下文表格是这一轮的全部 IFN/APM 结果。

## Contrast

[GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) on [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555) (BWH Human Release 3.0, 60-mer). Three two-colour arrays, GSM558700–GSM558702. GEO labels Cy3 / channel 1 as CLDN4 overexpression (control) and Cy5 / channel 2 as CLDN4 lentiviral siRNA. The series matrix VALUE is the deposited normalized sample-to-control ratio. CLDN4 itself is negative, so the ratio is oriented as siRNA relative to overexpression.

The processing text says a dye swap was done. The deposited channel labels are the same on all three samples. The primary analysis uses those deposited values as published. The raw ScanArray files show that GSM558701 was already sign-reversed in the deposit; that check is in the final sweep below.

## What the old symbol map missed

The earlier C4 slice (PR #104) and the prerank GSEA (PR #287) kept the ORF string, or the text before `--` when ORF was empty, and only rewrote `G1P2` and `G1P3`. Of 24,986 probes with a token, **9,334** tokens are not a current HGNC symbol.

This run resolves each token against the HGNC complete set downloaded 2026-09-21 (45,083 rows):

| Class | Probes |
|---|---|
| Current approved symbol | 15,650 |
| Unique previous symbol | 3,937 |
| Unique alias | 3,930 |
| Ambiguous (dropped) | 339 |
| Token with no HGNC hit | 1,130 |
| No symbol token | 11,302 |

Genes with values on at least two arrays: 14,738 after collapsing aliases, 15,768 when obsolete tokens were left as separate names.

Panel genes that the old map did not score under their current name, and that this map does:

| Current symbol | Array token | Median log2 | Direction |
|---|---|---|---|
| NLRC5 | NOD27 | +0.671 | up |
| ERAP2 | LRAP | +0.057 | up, near zero |
| PDIA3 | GRP58 | −0.737 | down |
| IFI44L | C1orf29 | −1.029 | down |
| DDX60 | FLJ20035 | +0.972 | up |
| XAF1 | HSXIAPAF1 | +0.444 | up |
| EIF2AK2 | PRKR | −0.179 | down |
| IL32 | NK4 | +0.355 | up |
| CGNL1 | JACOP | +0.251 | up |

Also annotated, not entered in the set test: ERAP1 is `ARTS-1` and has a value on GSM558702 only (−0.667). RSAD2 is `cig5` and has no deposited value. DDX58 is a previous symbol of **RIGI** in this HGNC file; the probe is kept and scored as DDX58 in the IFN panel (median −1.448, down). TACSTD2 is not on the platform.

## Thesis test

User KD thesis used as the alternative, not as a result: IFN and antigen presentation go **up** after CLDN4 loss; other tight-junction genes go **down**, because CLDN4-high is the barrier. CLDN4 is the perturbation check and is not inside the TJ set.

Gene statistic: median across arrays of the per-array probe median. A gene needs two arrays. Primary test: one-sided Mann-Whitney against all other mapped genes. BH is across the four panels below. The background median is −0.044, so the array is nearly centered and a set has to move past that.

| Panel | Thesis | Measured | Up / down | Set median | One-sided MW | Two-sided MW | BH q | Call |
|---|---|---|---|---|---|---|---|---|
| Priority (IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A) | up | 6/6 | 3 / 3 | −0.744 | 0.881 | 0.239 | 0.973 | opposite |
| APM (16 genes) | up | 15/16 | 5 / 10 | −0.560 | 0.973 | 0.054 | 0.973 | opposite |
| IFN / ISG (73 genes) | up | 64/73 | 29 / 35 | −0.101 | 0.888 | 0.225 | 0.973 | opposite |
| TJ core without CLDN4 | down | 14/15 | 6 / 8 | −0.247 | 0.516 | 0.968 | 0.973 | same sign, not significant |

Using the mean across arrays instead of the median does not change any call. The old map’s set medians were −0.744 (priority), −0.906 (APM, 12 genes), −0.326 (IFN, 58 genes). Putting the missing genes back moved APM and IFN toward zero. It did not turn them up.

Named genes, median log2(siRNA/OE):

| Gene | Median | On the two same-sign arrays | vs thesis (up) |
|---|---|---|---|
| IFI27 | −2.322 | −2.483 | opposite |
| OAS2 | +0.111 | −0.088 | small up; flips if GSM558701 is dropped |
| IFIT1 | −2.120 | −2.478 | opposite on all three arrays |
| MX1 | +0.358 | +0.358 | up (missing on GSM558701) |
| ISG15 (probe G1P2) | +0.714 | −0.173 | up; flips if GSM558701 is dropped |
| HLA-A (11 probes) | −1.600 | −1.668 | opposite; 7/11 probe medians negative |
| TAP1 | −1.472 | −1.687 | opposite |
| TAP2 | −0.556 | −0.362 | opposite |

The three priority genes that are up are the smaller effects. IFI27, IFIT1, and HLA-A are largely down. Classical MHC-I (HLA-A/B/C, B2M, TAP1, TAP2, TAPBP, CANX) is down. The clearest APM increase is PSMB9 (+1.709). Recovered NLRC5 is up (+0.671).

## Tight junction

CLDN4 probe 17169 is missing on GSM558700 and negative on the other two (−1.737, −0.713; median −1.225). That matches siRNA lower than overexpression. It does not, by itself, say the barrier program moved.

Of the other core genes, CLDN1, CLDN3, CLDN7, CLDN18, OCLN, TJP1, MARVELD2, and MARVELD3 are down. TJP2, TJP3, F11R, JAM3, CGN, and CGNL1 are up. The set median is slightly negative and the mean of those medians is slightly positive (+0.073). Against the background, p = 0.516. JAM2 is present on one array only and is not in the test.

## Array agreement

Spearman across 12,219 genes present on all three arrays: GSM558700 vs GSM558702 ρ = +0.244; GSM558700 vs GSM558701 ρ = −0.238; GSM558701 vs GSM558702 ρ = −0.197. GSM558701 is the array that often carries the opposite sign on IFN and MHC-I genes. It was not inverted.

Sensitivity that keeps only GSM558700 and GSM558702 (both values required): priority 1 up / 5 down, median −0.920; APM 4 / 11, median −0.959; IFN 21 / 37, median −0.351. All three stay opposite the “up” thesis (one-sided p = 0.961, 0.966, 0.972). TJ core on that pair is 2 up / 11 down, median −0.359, one-sided p = 0.092 (q = 0.369). The sign test on that pair is 0.011 and the Wilcoxon test is 0.029, so the genes lean down, but they are not separated from the rest of the array at p < 0.05.

## Limits

Ovarian SKOV-3, not lung. The control is CLDN4 overexpression, not parental or scramble cells, so a low TJ ratio can be overexpression-high rather than siRNA-low. n = 3, and the arrays disagree. HLA-A is an 11-probe median, not a single oligo. This accession does not show IFN or APM opening after the symbol repair. It is not a substitute for a lung knockout or for a private knockdown, and it is not used here to retract those.

## Reproduce

```bash
python3 scripts/gse22493_cldn4_remap/analyze.py
```

Downloads the GEO series matrix, GPL10555 SOFT, and the HGNC complete set into `/tmp/gse22493` (override with `GSE22493_CACHE`). Needs numpy, pandas, scipy, matplotlib.

The last sweep is below. It is the stopping point.

## FINAL: still discordant

No orientation-valid, experiment-level processing puts both IFN and APM up at p < 0.05. Stop.

The rule was fixed before reading these extra results. A processing counts only if, on that same processing, the IFN set and the APM set each have a positive effect, a one-sided up-test p < 0.05, and CLDN4 is negative on every array where it is measured. One array by itself does not count. Negating GSM558701 back to the raw Ch2/Ch1 sign does not count: that undoes the dye-swap correction already in the GEO matrix.

### What the raw files are

GEO has one supplement, `GSE22493_RAW.tar` (file list dated 2013-01-17): `GSM558700.txt.gz`, `GSM558701.txt.gz`, `GSM558702.txt.gz`. ScanArray Express exports. No other processed table and no author gene list.

Deposited VALUE is rank-identical to the ScanArray normalized Ch2/Ch1 log ratio on GSM558700 and GSM558702 (Spearman +1). On GSM558701 it is rank-identical to the negation of that column (Spearman +1). CLDN4 on that array is raw log2(Ch2/Ch1) = +1.22 and deposited = −1.737. The deposit is the dye-swap-corrected ratio. Flag 3 spots are the quantified spots (median background-subtracted intensity in the hundreds). Flag 1 spots are dim (median about 20–30).

### Sweep table

Effect is the set median log2(siRNA/OE) unless the row says otherwise. `p_up` is one-sided Mann-Whitney against the other genes or probes, except ssGSEA and the n=3 GSVA score, where it is a 999-permutation one-sided p (seed 42). Full rows, including TJ, STING, and NHEJ, are in `tables/sweep_final.tsv`.

| Processing | IFN n (up/down) | IFN effect | IFN p_up | APM n (up/down) | APM effect | APM p_up | CLDN4 negative on measured arrays |
|---|---|---|---|---|---|---|---|
| Deposited ratios | 64 (29/35) | −0.101 | 0.888 | 15 (5/10) | −0.560 | 0.973 | yes |
| Deposited, array median-centered | 64 (31/33) | −0.079 | 0.883 | 15 (4/11) | −0.482 | 0.976 | yes |
| Deposited, quantile-normalized | 64 (29/35) | −0.125 | 0.888 | 15 (4/11) | −0.551 | 0.975 | yes |
| Raw Ch2/Ch1, no flip | 62 (13/49) | −0.892 | 0.966 | 15 (3/12) | −0.796 | 0.726 | no (GSM558701 CLDN4 +1.22) |
| Raw flag 3, no flip | 55 (10/45) | −0.895 | 0.967 | 15 (3/12) | −0.796 | 0.763 | no |
| Raw, GSM558701 flipped to the deposited sign | 62 (20/42) | −0.589 | 0.724 | 15 (4/11) | −0.559 | 0.588 | yes |
| Raw flag 3, GSM558701 flipped | 55 (16/39) | −0.412 | 0.761 | 15 (4/11) | −0.559 | 0.586 | yes |
| Raw flag 3, flipped, median-centered | 55 (29/26) | +0.066 | 0.780 | 15 (8/7) | +0.070 | 0.507 | yes |
| Leave out GSM558700 | 59 (33/26) | +0.084 | 0.347 | 15 (8/7) | +0.062 | 0.490 | yes |
| Leave out GSM558701 | 58 (21/37) | −0.351 | 0.972 | 15 (4/11) | −0.959 | 0.966 | yes |
| Leave out GSM558702 | 57 (31/26) | +0.013 | 0.533 | 15 (7/8) | −0.089 | 0.342 | yes |
| Probe-level, no gene collapse | 98 (47/51) | −0.037 | 0.901 | 52 (14/38) | −0.850 | 1.000 | yes |
| \|effect\| ≥ genome q50 (0.600) | 43 (16/27) | −0.713 | 0.902 | 9 (2/7) | −1.252 | 0.972 | yes |
| \|effect\| ≥ genome q75 (1.104) | 21 (6/15) | −1.434 | 0.818 | 7 (1/6) | −1.472 | 0.869 | yes |
| \|effect\| ≥ genome q90 (1.711) | 8 (4/4) | −0.159 | 0.589 | 1 (0/1) | −2.556 | not tested | yes |
| \|effect\| ≥ \|CLDN4\| (1.225) | 20 (6/14) | −1.441 | 0.681 | 7 (1/6) | −1.472 | 0.756 | yes |
| Within-array percentile (up would be > 0.5) | 64 (31/33) | 0.472 | 0.916 | 15 (4/11) | 0.346 | 0.979 | yes |
| ssGSEA, α=0.25, median of 3 array ES | 64 | +77.3 | 0.773 | 15 | −936 | 0.955 | yes |
| GSVA-style KS on n=3 across-sample ranks | 55 | −452 | 0.971 | 15 | −937 | 0.981 | yes |

Two processings have both medians slightly positive. Neither is close to p < 0.05.

- Dropping GSM558700 leaves IFN at +0.084 (p = 0.347) and APM at +0.062 (p = 0.490).
- Median-centering the flag-3 raw ratios, after putting GSM558701 back on the deposited sign, leaves IFN at +0.066 (p = 0.780) and APM at +0.070 (p = 0.507).

Genes with a large absolute ratio are the ones that go down. Above the CLDN4 bar (\|log2\| ≥ 1.225), IFN is 6 up / 14 down (median −1.441) and APM is 1 up / 6 down (median −1.472).

ssGSEA on the deposited ranking does not rescue this. The IFN enrichment score is positive on GSM558700 (+77) and GSM558701 (+1759) and negative on GSM558702 (−130). The median score is positive and the permutation p on the gene-median ranking is 0.773. APM enrichment is negative on two arrays and the permutation p is 0.955. The n=3 GSVA-style score is negative for both sets.

### The one array that goes up

GSM558701 alone, using the deposited (dye-swap-corrected) ratios, does look like the thesis: IFN median +0.640 (37/24 up, p = 0.021), APM +0.872 (12/3 up, p = 0.011), priority 4/5 up (p = 0.035), CLDN4 −1.737. GSM558702 is the opposite (IFN −0.537, p_up = 0.977; APM −0.644). GSM558700 has no CLDN4 value and its IFN median is −0.360. GSM558701 is the array that anti-correlates with the other two (Spearman about −0.24 and −0.20). Keeping any second array wipes the up-shift out. That replicate is reported. It is not the result of the experiment.

### STING and NHEJ

These were scored. They are not the IFN/APM call.

On deposited gene medians, STING (12/14 genes; STING1 and TBK1 absent) median +0.508, up-test p = 0.127. ssGSEA is positive on all three arrays, permutation p = 0.118. Probe-level STING is the only up-test in this sweep with p < 0.05 (17 probes, 12 up, median +0.526, p = 0.042). That is one probe test among the sweep, and the gene-level test does not agree at p < 0.05. NHEJ (11/12 genes; PAXX absent) median +0.138, p = 0.690.

### Stop

Symbol repair, dye-swap-aware raw ratios, flag filtering, median centering, quantile normalization, within-array ranks, probe-level collapse, leave-one-out of each array, absolute-effect quantiles, ssGSEA, and an n=3 GSVA-style score were all run on the only raw files GEO holds. IFN/APM do not come out up. This accession stays discordant with that thesis. No further cut of these three arrays is planned.

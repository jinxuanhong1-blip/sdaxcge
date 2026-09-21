# GSE22493: CLDN4 siRNA vs CLDN4 overexpression after HGNC symbol repair

## 中文

GSE22493（SKOV-3-IP-Luc，CLDN4 siRNA / CLDN4 过表达，3 张双色芯片）先前把 GPL10555 的 2006 年 ORF 字符串当成现行基因符号，只手工改了 G1P2→ISG15、G1P3→IFI6。NOD27、GRP58、ARTS-1、LRAP、cig5、PRKR、C1orf29、FLJ20035 等因此进不了 IFN/APM 检验。本分析用 HGNC complete set（批准符号，否则唯一的曾用名，否则唯一的别名）重映射后再测用户 KD 方向：CLDN4 降低后 IFN 与 APM 应升高，其它紧密连接基因应降低。

符号补上之后，结论没有翻向。以中位 log2(siRNA/OE) 为准，优先级 6 基因 3 升 3 降（中位 −0.744），APM 5/15 升（中位 −0.560），IFN 集 29/64 升（中位 −0.101）。相对背景基因的单侧 Mann-Whitney 都不支持“升高”（p = 0.881、0.973、0.888；四组 BH q = 0.973）。CLDN4 本身中位 −1.225（2 张芯片），方向与标记一致。去掉 CLDN4 后的紧密连接核心 8/14 为负（中位 −0.247），但和背景无差别（p = 0.516），基因均值符号相反。这是卵巢 siRNA 对过表达，不是肺，也不是亲本对照；三张芯片彼此 Spearman 只有 −0.24 到 +0.24。它不支持“CLDN4 敲低打开 IFN/APM”，也不拿来否定别的公共或私有 KD。

## Contrast

[GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) on [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555) (BWH Human Release 3.0, 60-mer). Three two-colour arrays, GSM558700–GSM558702. GEO labels Cy3 / channel 1 as CLDN4 overexpression (control) and Cy5 / channel 2 as CLDN4 lentiviral siRNA. The series matrix VALUE is the deposited normalized sample-to-control ratio. CLDN4 itself is negative, so the ratio is oriented as siRNA relative to overexpression.

The processing text says a dye swap was done. The deposited channel labels are the same on all three samples. No array was sign-flipped here.

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

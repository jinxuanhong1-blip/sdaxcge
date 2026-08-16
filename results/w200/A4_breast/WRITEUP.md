# A4 · TISMO mammary (breast) models — Tacstd2 after ICB

Honest test of whether **Tacstd2** (Trop2) changes after immune-checkpoint blockade in the TISMO syngeneic mammary models. Output directory: `results/w200/A4_breast/`.

> **Claim under test:** Tacstd2 is induced, or otherwise systematically shifted, by ICB in TISMO breast / mammary models.
>
> **Verdict: not supported.** There is no consistent after-ICB shift. The TISMO website stars two contrasts; neither is a clean breast-wide ICB effect.

---

## 中文一句话

TISMO 七个有 ICB 对照的乳腺模型里，**Tacstd2 在 ICB 后没有一致变化**。29 个队列内 Mann-Whitney 全部 FDR 不显著；唯一未校正 p<0.05 的是 p53-2336R 终点非应答（2/5 只肿瘤把均值拉高）。TISMO 给 EMT6 两个臂标了同一颗星（p=0.022 被复用），单抗 PDL1 臂几乎不动。高表达模型（4T1 / KPB25L / p53-2225L）治疗后仍高，低表达模型（T11 / EMT6 / E0771）仍低。

## English TL;DR

In the seven TISMO mammary models that actually have a within-study ICB contrast, **Tacstd2 does not move consistently after ICB**.

| Check | Result |
|---|---|
| Independent MWU, 29 TISMO boxes, BH-FDR | **0 / 29** q < 0.05 |
| Unadjusted MWU p < 0.05 | **1 / 29**: p53-2336R end-stage NR vs baseline (p = 0.045) |
| That one signal | Mean Δ = +2.02, but **2 of 5** treated tumors sit at ~5 log TPM; the other three stay near baseline. Leave-one-out shift = 0.71 |
| TISMO's own stars | EMT6 GSE107801 anti-PDL1 **and** M7824 share **one** p = 0.022. Displayed-value MWU: anti-PDL1 Δ = +0.045, p = 0.29; combo Δ = +0.65, p = 0.054, one tumor at 4.55 |
| Latest-timepoint paired means (n = 15) | 8 up / 7 down; mean Δ = +0.023; Wilcoxon p = 0.98; paired t p = 0.90 |
| Model-clustered bootstrap of that mean Δ | 95% CI −0.31 to +0.62; p = 0.90 |
| Responder vs non-responder in the same box | **None.** TISMO splits R and NR into separate boxes, each vs its own baseline |
| High vs low expressors | 4T1 / KPB25L / p53-2225L stay high (~2–4). T11 / EMT6 / E0771 stay low (~0–0.6) |

This is **not** evidence that ICB induces Trop2 in breast models, and it is **not** evidence that Trop2-high tumors fail ICB. It is a within-study, bulk RNA, treated-vs-control comparison of one gene.

---

## What was asked, and what was actually testable

TISMO (Zeng et al., *NAR* 2022; https://tismo.pku-genomics.org) uniformly processed syngeneic RNA-seq. The Gene module's in-vivo ICB view documents seven **Mammary cancer** models: **4T1, E0771, EMT6, T11, KPB25L, p53-2225L, p53-2336R**.

The annotation table has more mammary lines (4T07, 67NR, 6DT1, Met1, Mvt1, Py230, …). Those have **no ICB arm** (`ICB=0`) and were not forced into a fake contrast. The seven Gene-module models are exactly the seven mammary lines with `ICB=1`.

TISMO returns 29 boxes (cell-line × study × condition × regimen). Many are timepoints of the same GSE124821 tumors (day3 / day7 / end) or GSE130472 arms that **share one isotype baseline**. They are not 29 independent experiments.

---

## Methods

- **Source:** public TISMO Gene-module table (`/rtismo/gene/downVivoExprn`) for Tacstd2, all six ICB regimen labels, the seven mammary models. Official metadata: `TISMO_vivosample_annotations.csv`. Official TISMO boxplot: `tismo_official_plot.jpg`.
- **Units:** TISMO quantile-normalised, ComBat-corrected log-scale TPM, as served. Not raw counts. TISMO's attached `pvalue` is the site's precomputed within-cohort statistic (paper: DESeq2 Wald on counts). Those p-values are recorded and then **audited**, not trusted.
- **Primary test:** two-sided Mann-Whitney U on the displayed values, baseline (`Baseline=1`) vs ICB (`Baseline=0`), one TISMO box at a time. BH-FDR across the 29 boxes. Cliff's delta as the effect size.
- **Honesty filters:** (1) same numeric TISMO p on two different boxes = reused, not a per-arm test; (2) leave-one-out shift of the treated mean ≥ 0.5 with fewer than two treated tumors > 2 log TPM = outlier-driven; (3) no pooling of day3/day7/end as independent; latest timepoint kept as a secondary table; (4) no R-vs-NR test where TISMO did not put both in one box.
- **Secondary:** paired Wilcoxon / t-test and a cell-line clustered bootstrap on the 15 latest-timepoint mean deltas. These are reported so a reader who wants a single number can see that it is null. They are **not** a breast-wide causal estimate.

Reproduce:

```bash
python3 scripts/w200/A4_breast/run_analysis.py
```

Requires network access to `tismo.pku-genomics.org` and the TISMO Aliyun share. Writes `results/w200/A4_breast/`.

---

## Results, model by model

### 4T1 — high baseline, no ICB shift

Seven latest contrasts (GSE130472 young/old × anti-CTLA4/anti-PDL1; GSE132529 anti-PD1; GSE137818 Brca1/Brca2 KO anti-PD1). Baseline means ~1.8–4.6. Deltas from −0.60 to +0.19. All MWU p ≥ 0.054. GSE130472 young anti-CTLA4 and young anti-PDL1 **share TISMO p = 0.482**. Old anti-CTLA4 and old anti-PDL1 share the same 7 isotype-old baselines (mean 2.51); they are not independent.

### E0771 — low/mid, null

One contrast, GSE174053 high-fat-diet anti-PD1, all non-responders. Baseline 0.60 vs treated 0.61. MWU p = 0.16.

### EMT6 — TISMO star is not a clean anti-PDL1 effect

GSE107801, both arms labelled Responders, both given TISMO p = **0.0223871957549608** and a green star.

| Arm | n | mean baseline | mean ICB | Δ | MWU p | Honest read |
|---|---|---|---|---|---|---|
| anti-PDL1 vs isotype | 10+10 | 0.148 | 0.193 | +0.045 | 0.29 | null |
| M7824 (PDL1×TGFβ trap) vs mutant-trap control | 10+10 | 0.097 | 0.743 | +0.647 | 0.054 | one tumor (SRX3453192) = 4.55; the other nine treated values are 0.01–0.93 |

Do not quote the TISMO star as “EMT6 Tacstd2 goes up after anti-PDL1.” The anti-PDL1-only arm does not move. The combo arm is a different drug and is outlier-sensitive.

### T11 — essentially off, stays off

Nine GSE124821 boxes (parental / Apobec / UV × day3 / day7 / end). All means < 0.45. Latest parental end NR Δ = −0.007; Apobec end R Δ = +0.25 (MWU p = 0.11). TISMO p-values here are ~0.06–1.0 and do not agree with a stable induction.

### KPB25L — high, noisy, no consistent direction

Six GSE124821 timepoint boxes, all responders to anti-CTLA4+anti-PD1. Latest UV-end Δ = −1.34 (MWU p = 0.17); other timepoints flip sign. Not a replicable ICB effect.

### p53-2225L — high, flat

Day3 and end, non-responders. Means stay ~3.3–3.9. Δ ≈ −0.1. Null.

### p53-2336R — the only unadjusted p < 0.05, still not a clean call

| Timepoint | n base / ICB | mean base | mean ICB | MWU p | TISMO p |
|---|---|---|---|---|---|
| day3 | 3 / 5 | 0.00 | 0.18 | 0.33 | NA |
| end | 8 / 5 | 0.15 | 2.17 | **0.045** | **2.3×10⁻⁴** (***) |

End-stage treated values: **4.999, 4.625, 0.785, 0.414, 0.025**. Two tumors account for the mean. BH-FDR q = 1.0 among 29 tests. Day3 in the same model is null. Treat as a **single heterogeneous end-stage NR box**, not as “ICB induces Tacstd2 in p53 mammary models.”

---

## What this does **not** show

- It does not test whether baseline Tacstd2 **predicts** ICB response. TISMO never puts responders and non-responders in the same box, so R vs NR was not estimated.
- It does not test protein / Trop2 IHC, or ADC activity.
- It does not transfer to human breast ICB. These are mouse syngeneic bulk tumors, several with extra factors (age, high-fat diet, Brca KO, UV, Apobec, TGFβ trap).
- It does not license a 49/64-style “most models go up” headline. Restricted to mammary models, the latest-timepoint score is 8 up / 7 down and the paired test is null.

---

## Files

| File | What |
|---|---|
| `summary.json` | Machine-readable headline and denominators |
| `cohort_stats.csv` | All 29 TISMO boxes |
| `cohort_stats_latest_timepoint.csv` | One row per model×study×regimen×status, latest timepoint |
| `model_summary.csv` | Per-model collapse of the latest table |
| `sample_level_tacstd2.csv` | Per-sample values used for MWU |
| `mammary_model_census.csv` | Every mammary line in the annotation table, ICB yes/no |
| `tismo_tacstd2_breast_raw.csv` | Raw Gene-module download |
| `tismo_official_plot.jpg` | TISMO's own boxplot (stars = their p, including the reused EMT6 p) |
| `fig1_delta_forest_latest.png` | Forest of latest-timepoint mean Δ |
| `fig2_model_paired_means.png` | Baseline vs ICB cohort means by model |
| `fig3_flagged_samples.png` | Individual tumors for p53-2336R end and both EMT6 arms |
| `TISMO_vivosample_annotations.csv` | Official in-vivo metadata |
| `NOTES.md` | Short machine-oriented note |

---

## Limitations

1. Displayed values are batch-corrected TPMs; TISMO's stars come from a different (count-level) test and, in two studies, are copied onto more than one arm.
2. Small n per box (often 3–8). MWU cannot be very small; DESeq2 can, which is why p53-2336R looks more dramatic on the TISMO plot than on ranks of the served values.
3. GSE124821 dominates the mammary ICB set (T11, KPB25L, p53-2225L, p53-2336R). One paper, one dual CTLA4+PD1 regimen.
4. Shared baselines in GSE130472 inflate the number of 4T1 boxes.
5. 14 samples appear in more than one TISMO box (317 rows, 303 unique sample IDs).
6. No multiplicity correction can rescue a true effect that is not in the data; conversely, calling p53-2336R “significant” after 29 looks is a false-discovery risk.

---

## Bottom line

**Tacstd2 after ICB in TISMO breast models is a null / model-idiosyncratic result, not a directional ICB effect.** Quote the per-model table. Do not quote the TISMO EMT6 stars. Do not average the 29 boxes.

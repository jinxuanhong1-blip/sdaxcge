# Open NSCLC ICI bulk: CLDN4-high versus response

Public GEO and ArrayExpress tumor bulk only. One locked binary endpoint per cohort. Keratin residual is the within-cohort OLS residual of log2 CLDN4 on KRT8, KRT18, and KRT19. Odds ratios below 1 mean CLDN4-high tumors had lower odds of response. Hazard ratios above 1 mean CLDN4-high tumors had worse survival.

## Does this support the PDF 11-cohort figure?

The figure is CLDN4-high **OR = 0.42 (95% CI 0.18–0.95), k = 11**.

**Supported: No.**

Open NSCLC bulk primary median-split OR is 0.96 [0.69–1.33] at k=9. Point estimate rounds to 0.42: False. k equals 11: False. 0.42 lies inside the 95% CI: False.

This NSCLC meta is not that figure. Nine independent open tumor-bulk series have CLDN4 and a clear response label. That is not 11, and the random-effects interval does not contain 0.42. No cohort was added or dropped to move the point estimate.

Leaving out GSE218989 (the only large series) gives k=8, n=283, RE 1.05 [0.61–1.80], p=0.852, I²=0%. That interval is wide enough to contain 0.42, and the point estimate is not 0.42. The small-series meta is underpowered. It does not recover the PDF figure either.

## Locked primary results

| Analysis | Result |
|---|---|
| Response OR, raw CLDN4, median split | k=9, n=638, RE 0.96 [0.69–1.33], p=0.794, I²=0% |
| Response OR, keratin residual, median split | k=9, n=638, RE 1.00 [0.66–1.51], p=0.99, I²=11% |
| PFS HR, raw CLDN4, median split | k=3, n=129, RE 1.12 [0.74–1.68], p=0.592, I²=0% |
| PFS HR, keratin residual, median split | k=3, n=129, RE 0.94 [0.63–1.42], p=0.783, I²=0% |

Forests: `figures/forest_or_raw_median.png`, `figures/forest_or_keratin_residual_median.png`, `figures/forest_hr_raw_median.png`, `figures/forest_hr_keratin_residual_median.png`.

### Primary response endpoint (median split)

One endpoint per cohort. Extra DCB and MPR rows are below and are not in the pooled OR.

| Cohort | Measure | Endpoint | n | 2×2 or events | Effect [95% CI] | Fisher or Cox p |
|---|---|---|---:|---|---|---:|
| GSE126044 | keratin_residual | author_R_vs_NR | 16 | 1/7/4/4 (hiR/hiNR/loR/loNR) | 0.14 [0.01–1.76] | 0.282 |
| GSE126044 | raw | author_R_vs_NR | 16 | 1/7/4/4 (hiR/hiNR/loR/loNR) | 0.14 [0.01–1.76] | 0.282 |
| GSE135222 | keratin_residual | DCB_PFS_ge_180d | 27 | 4/10/3/10 (hiR/hiNR/loR/loNR) | 1.33 [0.24–7.56] | 1 |
| GSE135222 | raw | DCB_PFS_ge_180d | 27 | 4/10/3/10 (hiR/hiNR/loR/loNR) | 1.33 [0.24–7.56] | 1 |
| GSE166449 | keratin_residual | author_R_vs_NR | 22 | 4/7/3/8 (hiR/hiNR/loR/loNR) | 1.52 [0.25–9.30] | 1 |
| GSE166449 | raw | author_R_vs_NR | 22 | 3/8/4/7 (hiR/hiNR/loR/loNR) | 0.66 [0.11–4.00] | 1 |
| GSE190265 | keratin_residual | DCB_PFS_ge_6mo | 43 | 7/15/7/14 (hiR/hiNR/loR/loNR) | 0.93 [0.26–3.34] | 1 |
| GSE190265 | raw | DCB_PFS_ge_6mo | 43 | 8/14/6/15 (hiR/hiNR/loR/loNR) | 1.43 [0.40–5.16] | 0.747 |
| GSE190266 | keratin_residual | DCB_PFS_ge_6mo | 69 | 11/24/6/28 (hiR/hiNR/loR/loNR) | 2.14 [0.69–6.65] | 0.265 |
| GSE190266 | raw | DCB_PFS_ge_6mo | 69 | 11/24/6/28 (hiR/hiNR/loR/loNR) | 2.14 [0.69–6.65] | 0.265 |
| GSE207422 | keratin_residual | RECIST_CRPR_vs_SDPD | 24 | 6/6/11/1 (hiR/hiNR/loR/loNR) | 0.09 [0.01–0.94] | 0.0686 |
| GSE207422 | raw | RECIST_CRPR_vs_SDPD | 24 | 7/5/10/2 (hiR/hiNR/loR/loNR) | 0.28 [0.04–1.88] | 0.371 |
| GSE218989 | keratin_residual | author_R_vs_NR | 355 | 84/94/84/93 (hiR/hiNR/loR/loNR) | 0.99 [0.65–1.50] | 1 |
| GSE218989 | raw | author_R_vs_NR | 355 | 82/96/86/91 (hiR/hiNR/loR/loNR) | 0.90 [0.60–1.37] | 0.671 |
| GSE274975 | keratin_residual | RECIST_CRPR_vs_SDPD | 55 | 8/20/6/21 (hiR/hiNR/loR/loNR) | 1.40 [0.41–4.76] | 0.758 |
| GSE274975 | raw | RECIST_CRPR_vs_SDPD | 55 | 8/20/6/21 (hiR/hiNR/loR/loNR) | 1.40 [0.41–4.76] | 0.758 |
| GSE283829 | keratin_residual | RECIST_CR_vs_SDPD | 27 | 3/11/4/9 (hiR/hiNR/loR/loNR) | 0.61 [0.11–3.49] | 0.678 |
| GSE283829 | raw | RECIST_CR_vs_SDPD | 27 | 3/11/4/9 (hiR/hiNR/loR/loNR) | 0.61 [0.11–3.49] | 0.678 |

### Other response labels (not pooled)

| Cohort | Measure | Endpoint | n | 2×2 or events | Effect [95% CI] | Fisher or Cox p |
|---|---|---|---:|---|---|---:|
| GSE207422 | keratin_residual | MPR_vs_NMPR | 24 | 4/8/5/7 (hiR/hiNR/loR/loNR) | 0.70 [0.13–3.68] | 1 |
| GSE207422 | raw | MPR_vs_NMPR | 24 | 2/10/7/5 (hiR/hiNR/loR/loNR) | 0.14 [0.02–0.96] | 0.0894 |
| GSE274975 | keratin_residual | DCB | 53 | 17/10/16/10 (hiR/hiNR/loR/loNR) | 1.06 [0.35–3.23] | 1 |
| GSE274975 | raw | DCB | 53 | 17/10/16/10 (hiR/hiNR/loR/loNR) | 1.06 [0.35–3.23] | 1 |

GSE207422 major pathologic response has a raw median OR whose Woolf interval excludes 1 (0.14 [0.02–0.96]) and a Fisher p of 0.089. That disagreement is the small-sample case noted above. It is not called a significant result, and it is not the primary endpoint (RECIST). Keratin residual MPR is null (OR 0.70, p=1).

### PFS or OS (median split)

PFS Cox is limited to cohorts with a deposited time and an event indicator, and with follow-up that is not capped at 6 months. GSE190266 is DCB-only for that reason. GSE218989 deposits PFS days without an event indicator, so it contributes an OS HR (Death) and not a PFS HR. GSE274975 PFS drops the two non-NSCLC samples (small-cell and neuroendocrine).

| Cohort | Measure | Endpoint | n | 2×2 or events | Effect [95% CI] | Fisher or Cox p |
|---|---|---|---:|---|---|---:|
| GSE135222 | keratin_residual | PFS | 27 | 21 events | 0.76 [0.32–1.81] | 0.541 |
| GSE135222 | raw | PFS | 27 | 21 events | 1.18 [0.50–2.79] | 0.705 |
| GSE190265 | keratin_residual | PFS | 43 | 35 events | 1.11 [0.57–2.17] | 0.759 |
| GSE190265 | raw | PFS | 43 | 35 events | 0.86 [0.44–1.70] | 0.671 |
| GSE218989 | keratin_residual | OS | 355 | 265 events | 0.98 [0.77–1.24] | 0.853 |
| GSE218989 | raw | OS | 355 | 265 events | 1.08 [0.85–1.37] | 0.539 |
| GSE274975 | keratin_residual | PFS | 59 | 40 events | 0.92 [0.49–1.72] | 0.787 |
| GSE274975 | raw | PFS | 59 | 40 events | 1.36 [0.72–2.58] | 0.34 |

## What is underpowered, and what is a null

Primary raw median-split tests with Fisher p < 0.05: 0. Directions are mixed: some cohorts have OR above 1 and some below 1. Keratin regression R² for CLDN4 ranges from 0.31 to 0.83, so the residual is a real adjustment, and the residual meta is still null.

Small series cannot rule a moderate association in or out. A confidence interval that includes 1 is a null at α = 0.05. In the smallest 2×2 tables the Woolf interval and the Fisher p can disagree; the p used here is Fisher exact, and the meta-analysis uses the Woolf log-OR.

- GSE126044 (author_R_vs_NR): n=16, OR 0.14 [0.01–1.76], Fisher p=0.282 (n<40; fewer than 8 patients in one response class; CI includes 1)
- GSE135222 (DCB_PFS_ge_180d): n=27, OR 1.33 [0.24–7.56], Fisher p=1 (n<40; fewer than 8 patients in one response class; CI includes 1)
- GSE166449 (author_R_vs_NR): n=22, OR 0.66 [0.11–4.00], Fisher p=1 (n<40; fewer than 8 patients in one response class; CI includes 1)
- GSE190265 (DCB_PFS_ge_6mo): n=43, OR 1.43 [0.40–5.16], Fisher p=0.747 (CI includes 1)
- GSE190266 (DCB_PFS_ge_6mo): n=69, OR 2.14 [0.69–6.65], Fisher p=0.265 (CI includes 1)
- GSE207422 (RECIST_CRPR_vs_SDPD): n=24, OR 0.28 [0.04–1.88], Fisher p=0.371 (n<40; fewer than 8 patients in one response class; CI includes 1)
- GSE283829 (RECIST_CR_vs_SDPD): n=27, OR 0.61 [0.11–3.49], Fisher p=0.678 (n<40; fewer than 8 patients in one response class; CI includes 1)
- GSE274975 (RECIST_CRPR_vs_SDPD): n=55, OR 1.40 [0.41–4.76], Fisher p=0.758 (CI includes 1)
- GSE218989 (author_R_vs_NR): n=355, OR 0.90 [0.60–1.37], Fisher p=0.671 (CI includes 1)

GSE218989 (Samsung–KAIST) is the only open NSCLC bulk series large enough for a moderate OR to have a relatively narrow interval. Its raw median-split result is OR 0.90 [0.60–1.37] (Fisher p=0.671, 168 responders / 187 non-responders). That interval does not contain 0.42. The series dominates the random-effects pool. The leave-GSE218989 row in `tables/pooled.csv` is the underpowered small-cohort meta.

Neoadjuvant PD-1 plus chemotherapy (GSE207422) uses pathologic response only as a sensitivity (`MPR_vs_NMPR`). The primary binary for that series is RECIST. Dropping it does not create an 11-cohort OR of 0.42 (`primary_OR_raw_median_no_neoadjuvant`).

## Duplicate samples

Spearman correlation on the shared gene panel (CLDN4, keratins, and the other extracted genes). A smaller cohort is non-independent when at least half of its samples match a larger cohort at ρ ≥ 0.99. Non-independent series would stay in the per-cohort table and leave the primary pool.

36 cohort pairs compared on the shared gene panel. Pairs declared non-independent (at least half the smaller cohort at Spearman ≥ 0.99): 0. Highest best-match Spearman was 0.97 (GSE283829 vs GSE274975; median best-match in that pair 0.92). No sample pair reached 0.99, so all nine series stay in the primary pool. A shared patient sequenced on a different pipeline could be missed by this fingerprint. Full pairs: `tables/sample_overlap.csv`.

## Cohorts opened and not put in the forest

| Accession | Why it is not in the OR/HR forest |
|---|---|
| GSE136961 | DCB labels, CLDN4 not on the panel (confirmed in this run) |
| GSE161537 | RECIST and PFS, CLDN4 not on the HTG panel (confirmed in this run; keratins are present) |
| GSE93157 | NSCLC RECIST/PFS; NanoString 730 without CLDN4 |
| GSE309652 | R/NR on GEO; CLDN4 absent from the targeted panel |
| GSE253564 | Pre-treatment bulk FPKM, durvalumab ± SBRT; no MPR/RECIST/PFS on GEO |
| GSE248378 | Post-treatment resected tumors, not baseline R/NR |
| GSE248249 | Pre/post immunotherapy FFPE; no RECIST on GEO |
| GSE182328 | Tumor RNA on PD-1, label is Akkermansia, not response |
| GSE208858 | 4-1BB agonist, not PD-1/PD-L1 |
| GSE202417, GSE260770 | Blood or exosome, not tumor bulk |
| ArrayExpress | No leftover open NSCLC ICI tumor matrix with R/NR or PFS |

GSE135222 is not double-counted as ICB_Jung. GSE207422 single-cell libraries are not in this bulk meta.

## Methods

- Expression: counts become log2(CPM+1). TPM and FPKM become log2(x+1). GSE207422 is the author log2TPM and is not logged again.
- Median: high = at or above the cohort median on the analysis subset. Tertile (T3 vs T1) and quartile (Q4 vs Q1) are in `tables/effects.csv` and are not used to choose a pooled number.
- Binary OR: Woolf interval; Haldane–Anscombe 0.5 if a 2×2 cell is 0; Fisher exact p. The interval and the Fisher p answer different questions and can disagree when n is small.
- Continuous OR: logistic regression per 1 SD. Keratin-adjusted OR: same model plus z-scored KRT8, KRT18, and KRT19.
- Keratin residual: within cohort, CLDN4 ~ KRT8 + KRT18 + KRT19. The residual R² is in `tables/effects.csv` (`keratin_r2`).
- DCB: PFS ≥ 180 days or ≥ 6 months. Patients censored before that landmark are excluded from the 2×2, not called non-responders.
- Cox models use `statsmodels` PHReg. GSE274975 uses the Table S1 column “Response status available” as the event indicator because no other censor column is deposited.
- Meta-analysis: DerSimonian–Laird random effects and the I² from Cochran’s Q.
- Sample-level genes and labels: `tables/patient_level.tsv`.

## Reproduce

```bash
python3 -m pip install -r results/ici_nsclc_cldn4_bulk/requirements.txt
python3 results/ici_nsclc_cldn4_bulk/analyze.py
```

The script expects the public files already fetched under `/tmp/ici_raw` (GEO FTP, EuropePMC PMC11669362 Table S1, Springer Supplementary Data 8 for GSE218989). It does not re-download if those files are present.

## 中文

开放 NSCLC ICI 肿瘤 bulk（GEO / ArrayExpress）里，能同时拿到 CLDN4 和明确 R/NR 或 PFS 的独立队列是 9 个，不是 11 个。锁定的中位数切分、随机效应合并：应答 OR 0.96 [0.69–1.33]，p=0.79，I²=0%，n=638。角蛋白残差（CLDN4 ~ KRT8+KRT18+KRT19）OR 1.00 [0.66–1.51]，p=0.99。PFS HR 只有 3 个有事件时间且未被 6 个月截尾的队列，HR 1.12 [0.74–1.68]，p=0.59。各队列主终点 Fisher 检验无一 p<0.05。GSE218989（n=355）的区间不含 0.42；去掉它之后的小队列合并区间很宽，点估计也不是 0.42。因此这组公开数字不支持 PDF 上的 11 队列 OR=0.42。小队列本身检验效能不足，合并结果是阴性而不是把 0.42 证伪到每一项灵敏度上。

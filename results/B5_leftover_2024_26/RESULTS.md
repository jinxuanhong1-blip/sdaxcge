# Results: extra 2023–2026 open lung ICI rows (B5 not recut)

User B5 (11-cohort CLDN4-high ICI **OR = 0.42**) is **taken as given**. These rows are **additive only**. They are **not** pooled into that OR.

## Extra NSCLC rows (do not mix SCLC)

Locked primary = within-cohort **median** split. Continuous = rank-biserial (responders minus non-responders). Honest n after dropping missing RECIST / missing gene.

| Series | Gene | Endpoint | n | Effect | Value | 95% CI | p |
|---|---|---|---:|---|---:|---|---:|
| GSE274975 NSCLC | CLDN4 | ORR (CR/PR vs SD/PD), median | 56 | OR high vs low | **1.00** | 0.30–3.35 | 1.00 |
| GSE274975 NSCLC | TACSTD2 | ORR (CR/PR vs SD/PD), median | 56 | OR high vs low | **1.00** | 0.30–3.35 | 1.00 |
| GSE274975 NSCLC | CLDN4 | ORR, continuous | 56 (14 R / 42 NR) | rank-biserial | −0.034 | — | 0.86 |
| GSE274975 NSCLC | TACSTD2 | ORR, continuous | 56 (14 R / 42 NR) | rank-biserial | +0.150 | — | 0.41 |
| GSE274975 NSCLC | CLDN4 | DCB (PFS ≥ 6 mo), median | 60 | OR high vs low | **1.00** | 0.36–2.78 | 1.00 |
| GSE274975 NSCLC | TACSTD2 | DCB (PFS ≥ 6 mo), median | 60 | OR high vs low | 0.76 | 0.27–2.12 | 0.79 |
| GSE274975 NSCLC | CLDN4 | PFS, Cox per 1 SD | 60 | HR | 1.02 | — | 0.89 |
| GSE274975 NSCLC | TACSTD2 | PFS, Cox per 1 SD | 60 | HR | 0.98 | — | 0.90 |
| GSE283829 NSCLC | CLDN4 | CR vs PD, median | 17 | OR high vs low | **0.75** | 0.11–5.24 | 1.00 |
| GSE283829 NSCLC | TACSTD2 | CR vs PD, median | 17 | OR high vs low | **0.75** | 0.11–5.24 | 1.00 |
| GSE283829 NSCLC | CLDN4 | CR vs PD, continuous | 17 (7 CR / 10 PD) | rank-biserial | −0.029 | — | 0.96 |
| GSE283829 NSCLC | TACSTD2 | CR vs PD, continuous | 17 (7 CR / 10 PD) | rank-biserial | −0.086 | — | 0.81 |

Tertile T3-vs-T1 and CR/PR-vs-PD sensitivities are in `tables/all_tests.tsv`. All remain null. No cell is 0.42.

### What these series are

- **GSE274975** (Poddubskaya et al., *Front Immunol* 2024, PMID 39723204). Pretreatment FFPE bulk RNA-seq, STAR counts. RECIST 14 responders (3 CR + 11 PR) vs 42 SD/PD among NSCLC with both labels and genes. Regimens are mixed PD-1/PD-L1 ± chemo ± CTLA-4 ± PARPi. One GEO-annotated SCLC (`OB_pat_LuC_92`; Table S1 histotype “mixed”) is **excluded** from this NSCLC table. Four patients have PFS but no RECIST and are ORR-ineligible.
- **GSE283829** (Lindberg et al., *J Thorac Oncol* 2025). GEO `disease stage` is RECIST-like best response (CR 7 / SD 10 / PD 10). No PR labels, no PFS.

## Extra SCLC table (separate; not in the NSCLC OR)

Must-try GSE244944 / GSE244945 / GSE244946 **do not contain patient ICB response labels**. IMpower133 RNA (the clinical claim in that paper) is not on GEO.

| Series | Gene | Contrast | n | Effect | Value | p | Use |
|---|---|---|---:|---|---:|---:|---|
| GSE244945 H82 cell line | CLDN4 | NOTCH1-dox vs parental | 3 vs 3 | rank-biserial | +1.00 | 0.10 | mechanistic only |
| GSE244945 H82 cell line | TACSTD2 | NOTCH1-dox vs parental | 3 vs 3 | rank-biserial | +1.00 | 0.064 | mechanistic only |
| GSE244945 H69 cell line | CLDN4 | NOTCH1-dox vs parental | 3 vs 3 | rank-biserial | −1.00 | 0.10 | mechanistic only |
| GSE244945 H69 cell line | TACSTD2 | NOTCH1-dox vs parental | 3 vs 3 | rank-biserial | +1.00 | 0.077 | mechanistic only |
| GSE244944 | — | ChIP-seq | 18 | — | — | — | no expression-vs-ICB test |
| GSE244946 | — | COR-L88 scRNA | 3 | — | — | — | not patient ICB |

No SCLC **patient** OR/PFS row can be computed from these public files.

## Must-try GSE328294 (not an ICI OR)

H23 KRAS-mutant NSCLC **cell line**, 48 h DMSO vs IBI351, n=6 FPKM. No patient ICI labels in GEO.

| Gene | Median log2(FPKM+1) DMSO vs IBI351 | rank-biserial (IBI351 − DMSO) | exact MWU p |
|---|---|---:|---:|
| CLDN4 | 4.68 vs 3.79 (down) | −1.00 | 0.10 |
| TACSTD2 | 5.30 vs 8.47 (up) | +1.00 | 0.10 |

Direction is a drug-on-cell-line effect, not a predictor of ICI response. n=3 vs 3 cannot reject a null at α=0.05 even when ranks do not overlap.

## Bottom line

New open lung ICI **tumor bulk** from 2023–2026 that we could actually download and label (GSE274975 n=56 ORR; GSE283829 n=17 CR vs PD) are **null** for both CLDN4 and TACSTD2. Median ORs are 1.00 and 0.75. They do **not** recover user B5 OR=0.42, and they are **not** added to that 11-cohort pool.

## 中文摘要

B5 的 11 队列 CLDN4-high ICI OR=0.42 **按用户给定，不重切、不并入**。2023–2026 年可下载的新 NSCLC 肿瘤 bulk：GSE274975（ORR n=56）CLDN4 中位数 OR=1.00，p=1.00；GSE283829（CR vs PD n=17）OR=0.75，p=1.00。GSE328294 是 H23 细胞系 ± IBI351，不是患者 ICI。GSE244944/5/6 是 SCLC 细胞系/ChIP/scRNA，**没有**患者 ICB 标签，不进入 NSCLC OR。

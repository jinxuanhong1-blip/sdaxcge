# Additive GSE11969 LUAD — CLDN4 vs CD8A / CD274

**Additive public cohort.** Takeuchi / Tomida **GSE11969** (PMID 16549822; Tomida follow-up PMID 21465578): resected Japanese lung series on **GPL7015** Agilent Homo sapiens 21.6K custom array. Series text says 149 NSCLC including **90 adenocarcinomas (AD)**, plus 9 SCLC and 5 normal-lung mixtures. This slice uses only **Histology = AD**.

Primary question: **CLDN4 vs CD8A** and **CLDN4 vs CD274**. TACSTD2 is a same-run companion, not an audit of prior TACSTD2 claims.

**CD274 mapping.** GPL7015 has no `CD274` symbol. The deposited symbol is **PDCD1LG1** (probe `A_23_P256487`, EST `BU621474`, title “Programmed cell death 1 ligand 1”). That is the CD274 row. Do not write “CD274 absent.”

## One-row table

| dataset | histology | platform | n series | n NSCLC (text) | n LUAD | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | verdict CD8A / CD274 |
|---|---|---|---:|---:|---:|---|---|---|---|---|---|---|---|
| GSE11969 Takeuchi | LUAD (AD) | GPL7015 Agilent 21.6K log10(R/G) | 163 | 149 | **90** | `A_23_P19944` | `A_23_P68110` | `A_23_P256487` (PDCD1LG1) | **-0.403 (8.36e-05)** | -0.392 (1.44e-04) | **+0.025 (0.813)** | -0.048 (0.655) | **HOLDS / NO_EVIDENCE** |

Do not write n=163 or n=149 for the correlations. The computable LUAD n is **90** arrays / **90** patients (AD001–AD090, 1:1).

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **163** | 21,619 probes × 163 GSM; 0 missing VALUE cells |
| NSCLC (AD+SQ+LA+AS+LCNEC) | yes | **149** | matches series text “149 patients with NSCLC” |
| **LUAD (Histology = AD)** | yes | **90** | titles all “Patient … with adenocarcinoma”; unique Annotation AD001–AD090 |
| SQ / LA / AS / LCNEC | yes | 35 / 18 / 4 / 2 | not in the LUAD slice |
| SCLC | yes | 9 | excluded |
| Normal lung mixtures | yes | 5 | Cy5 channel; excluded |
| Unique LUAD patients | yes | **90** | one array per AD annotation |
| ICI / PD-1 / PD-L1 treatment | no | **0** | 8/90 LUAD marked Gefitinib=Y; not an ICI series |
| Tumor % / ESTIMATE on GEO | no | **0** | only public proxy is an RNA epithelial score |
| CLDN4 finite (`A_23_P19944`) | yes | **90** | named GPL7015 CLDN4 / NM_001305 |
| CD8A finite (`A_23_P68110`) | yes | **90** | named GPL7015 CD8A / NM_001768 |
| CD274 finite (`A_23_P256487`) | yes | **90** | GPL7015 symbol **PDCD1LG1** / EST BU621474 |
| **Primary pairwise n (CLDN4 + CD8A + CD274)** | yes | **90** | this is the n used below |
| CLDN4 Q4 vs Q1 | yes | **23 vs 23** | not 90; Q2+Q3 dropped |

The author patient-info supplement has extra blank rows if read naively (183 lines). Histology counts above come from the series matrix, not that file.

## Verdict

| Test | n | ρ | 95% CI | p | ρ_adj (epithelial 5/6) | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| **CLDN4 vs CD8A** | **90** | **-0.403** | -0.568 to -0.211 | **8.36e-05** | -0.392 | 1.44e-04 | **HOLDS** |
| **CLDN4 vs CD274 (PDCD1LG1)** | **90** | **+0.025** | -0.181 to +0.228 | **0.813** | -0.048 | 0.655 | **NO_EVIDENCE** |
| CLDN4 vs TACSTD2 (companion) | 90 | +0.232 | +0.014 to +0.423 | 0.028 | +0.159 | 0.137 | companion |
| CLDN4 vs epithelial mean-z | 90 | +0.403 | — | 8.30e-05 | — | — | companion |

**What holds.** Higher **CLDN4** tracks **lower CD8A** on LUAD arrays (n=90, ρ=-0.403, p=8.36e-05). The inverse remains after residualising on the 5-gene epithelial mean-z (EPCAM absent; used KRT8, KRT18, KRT19, CDH1, KRT7; ρ_adj=-0.392, p_adj=1.44e-04). HOLDS rule: n≥40, ρ<0, p<0.05.

**What does not hold.** **CLDN4 vs CD274/PDCD1LG1** is null (n=90, ρ=+0.025, p=0.813; adj ρ=-0.048, p=0.655). Do not upgrade the CD8A anti-correlation into a PD-L1-RNA claim. The CD274 probe is an EST, not NM_014143; that is a platform limit, not a reason to impute a different n.

## Q4 vs Q1 (not n=90)

`pd.qcut(rank(method="first"), 4)` on the 90 LUAD CLDN4 values. Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | r | p |
|---|---:|---:|---:|---:|---:|---:|
| **CD8A** | **23 / 23** | 0.069 | 0.333 | -0.264 | **-0.67** | **9.21e-05** |
| CD274 (PDCD1LG1) | 23 / 23 | 0.039 | 0.040 | -0.002 | +0.10 | 0.568 |

The CD8A quartile cut is the same immune-low direction as the continuous test. It is not a second cohort. Do not write n=90 for 23 vs 23.

## Methods

- **Matrix:** GEO `GSE11969_series_matrix.txt.gz`. Deposited VALUE = LOWESS-normalized, background-subtracted **log10(processed Red / processed Green)**. Cy5 = sample; Cy3 = 20-cell-line lung RNA reference. Spearman is rank-based, so the log10 ratio scale does not change ρ.
- **LUAD filter:** `Histology = AD` on the series matrix (n=90). Titles match. SQ/LA/AS/LCNEC/SCLC/normal are counted above and then dropped.
- **Probes:** official GPL7015 table (`Gene symbol`, `Title`, `GB_LIST`). Primary genes use the **named single probes** in `tables/probe_confirm.tsv`. No max-mean collapse for CLDN4 / CD8A / CD274 because each has one probe.
- **CD274:** GPL7015 symbol **PDCD1LG1**. Alias search (CD274, PD-L1, PDL1, B7-H1) found no second probe.
- **Epithelial residual:** unweighted mean of gene-wise z for KRT8, KRT18, KRT19, CDH1, KRT7 (5/6; missing: EPCAM). Max-mean collapse only for this score. Partial Spearman = Pearson of rank residuals; df = n − 3.
- **HOLDS:** n≥40, ρ (or ρ_adj) < 0, p < 0.05. Used for the immune-low claim only.
- **Not tested:** OS/relapse (labels exist on GEO but are outside this slice), ICI response (none), pathologist CD8/PD-L1 IHC, purity/ESTIMATE (not deposited).

## What this does not claim

- It does not treat n=163 or n=149 as the LUAD correlation n.
- It does not claim CLDN4 vs PD-L1 protein. The CD274 row is one EST probe.
- It does not claim an ICI-response effect. 8 LUAD arrays have Gefitinib=Y; that is EGFR TKI, not PD-1/PD-L1 blockade.
- TACSTD2 numbers are a same-run companion only.

## Reproduce

```bash
python3 -m pip install -r methods/gse11969_cldn4/requirements.txt
python3 methods/gse11969_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE11969_CLDN4_DATA` (default `/tmp/gse11969_cldn4`).

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/n_table.tsv` — honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / PDCD1LG1 / TACSTD2 probes
- `tables/spearman_primary.tsv` — CLDN4 vs CD8A / CD274 / TACSTD2
- `tables/highlow_q4q1.tsv` — CLDN4 Q4 vs Q1
- `tables/sample_annotation.tsv` — 90 LUAD arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`

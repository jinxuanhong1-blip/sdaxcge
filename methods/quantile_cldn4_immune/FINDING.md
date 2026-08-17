# FINDING — extreme-quartile recut (CLDN4 / TACSTD2 Q4 vs Q1)

**Additive only. No new matrices. No invented signatures.**

Continuous Spearman clouds on these same tables are weak-to-moderate.
This folder **only recuts** existing public sample-level scores as **Q4 vs Q1**
and reports **δ, p, n per arm**.

**Verdict.** CLDN4-high (Q4) is **immune-low** vs Q1 for CD8A / ImmuneScore /
GEP18 on GSE218989, GSE68465, GSE31210, and CPTAC-LSCC protein. TCGA-LUAD is
the same sign, smaller. TCGA-LUSC CLDN4 is **null**. GSE285029 is the
**exception**: CLDN4 Q4 is IFN / MHC-I / GEP18 **high** (replicates that
table’s own Q4 file). TACSTD2 Q4 is CD8-low on GSE218989, GSE31210, OncoSG,
TCGA-LUSC, and GSE68465; null on GSE285029 and both CPTAC protein tables.
**OS is not significant** wherever survival was already on the table.

---

## What was recut (nothing invented)

Copies of existing sample-level tables; sources in `harvested/PROVENANCE.md`.

| Cohort | n on table | CLDN4 | TACSTD2 | CD8A | ImmuneScore | GEP18 | IFN | MHC-I | OS |
|---|---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| GSE218989 | 355 | yes | yes | yes | ESTIMATE ssGSEA | **no** | Ayers6 | HLA-A/B/C mean-z | yes |
| GSE285029 | 234 | yes | yes | yes | A1 8-gene (`immune8`) | yes | IFN-compact | MHC1 | **no** |
| TCGA-LUAD (Xena STAR) | 502 | yes | yes | yes | ESTIMATE | yes | Wolf IFNγ (n=444) | **no** | **no** |
| TCGA-LUSC (Xena) | 502 | yes | yes | CD8 score | ESTIMATE | yes | **no** | **no** | **no** |
| CPTAC LUAD **protein** | 110 (CLDN4 prot. 79) | protein | protein | CD8 T-cell RNA | ESTIMATE RNA | GEP RNA | IFNG-6 RNA | **no** | yes |
| CPTAC LSCC **protein** | 108 (CLDN4 prot. 78) | protein | protein | CD8A RNA | ESTIMATE RNA | GEP RNA | Hallmark IFNγ | **no** | yes |
| GSE68465 | 443 | yes | yes | yes | ESTIMATE | 15/18 | **no** | **no** | **no** |
| GSE31210 | 226 | yes | yes | yes | ESTIMATE | yes | **no** | **no** | **no** |
| OncoSG (public z) | 169 | **no** | yes | CD8A z | A1 8-gene | 17/18 | IMSIG IFN | **no** | **no** |

OncoSG CLDN4 is **not on the public cBioPortal z-score matrix** used in PR
OncoSG A1. It was not imputed.

---

## Methods (fixed before looking at p)

- Anchor quartiles: `rank(method='first')` then `qcut` into 4 equal-count bins.
  Q1 = lowest 25%, Q4 = highest 25%. Complete-case per (anchor, endpoint).
- Continuous endpoints: two-sided Mann–Whitney U (Q4 vs Q1).
  **δ_median** = median(Q4) − median(Q1) on the stored scale.
  **δ_rb** = rank-biserial \(2U/(n_4 n_1)-1\) (positive = higher in Q4).
  95% CI on δ_rb = percentile bootstrap, 2000 resamples, seed 20260817.
- OS (only if `time` + `event` already exist): log-rank; Mantel–Haenszel HR
  (Q4 vs Q1). Median OS is KM; “NA” = not reached.
- No purity residual, no new gene-set, no download.

GSE285029 CLDN4 Q4 vs Q1 medians match that branch’s `q4_vs_q1.csv`
(CD8A +0.354, GEP18 +0.418, IFN-compact +0.796, MHC1 +0.572; n=59/59).

---

## CLDN4 Q4 vs Q1

δ = median(Q4) − median(Q1). n = n_Q1 / n_Q4.

### CD8A

| Cohort | n | δ_median | δ_rb | p |
|---|---:|---:|---:|---:|
| GSE218989 | 89 / 89 | −0.856 | −0.307 | **4.1×10⁻⁴** |
| GSE285029 | 59 / 59 | +0.354 | +0.176 | 0.10 |
| TCGA-LUAD | 126 / 126 | −0.156 | −0.153 | **0.036** |
| TCGA-LUSC | 126 / 126 | −0.256 | −0.082 | 0.26 |
| CPTAC-LUAD-protein | 20 / 20 | −0.770 | −0.290 | 0.12 |
| CPTAC-LSCC-protein | 20 / 20 | −1.15 | −0.580 | **0.0018** |
| GSE68465 | 111 / 111 | −99.8 | −0.304 | **9.3×10⁻⁵** |
| GSE31210 | 57 / 57 | −268 | −0.483 | **8.9×10⁻⁶** |

### ImmuneScore

| Cohort | n | δ_median | δ_rb | p |
|---|---:|---:|---:|---:|
| GSE218989 | 89 / 89 | −0.134 | −0.372 | **1.8×10⁻⁵** |
| GSE285029 | 59 / 59 | +0.392 | +0.198 | 0.064 |
| TCGA-LUAD | 126 / 126 | −310 | −0.195 | **0.0074** |
| TCGA-LUSC | 126 / 125 | −27.0 | −0.049 | 0.50 |
| CPTAC-LUAD-protein | 20 / 20 | −606 | −0.160 | 0.39 |
| CPTAC-LSCC-protein | 20 / 20 | −2189 | −0.620 | **8.4×10⁻⁴** |
| GSE68465 | 111 / 111 | −627 | −0.354 | **5.3×10⁻⁶** |
| GSE31210 | 57 / 57 | −799 | −0.528 | **1.2×10⁻⁶** |

### GEP18 (absent on GSE218989)

| Cohort | n | δ_median | δ_rb | p |
|---|---:|---:|---:|---:|
| GSE285029 | 59 / 59 | +0.418 | +0.334 | **0.0018** |
| TCGA-LUAD | 126 / 126 | −0.059 | −0.045 | 0.54 |
| TCGA-LUSC | 126 / 126 | −0.172 | −0.096 | 0.19 |
| CPTAC-LUAD-protein | 20 / 20 | −0.290 | −0.285 | 0.13 |
| CPTAC-LSCC-protein | 20 / 20 | −0.785 | −0.670 | **3.1×10⁻⁴** |
| GSE68465 | 111 / 111 | −0.517 | −0.331 | **2.1×10⁻⁵** |
| GSE31210 | 57 / 57 | −0.420 | −0.415 | **1.4×10⁻⁴** |

### IFN / MHC-I (only if the column already existed)

| Cohort | endpoint | n | δ_median | δ_rb | p |
|---|---|---:|---:|---:|---:|
| GSE218989 | IFN Ayers6 | 89 / 89 | +0.013 | +0.019 | 0.83 |
| GSE218989 | MHC-I | 89 / 89 | +0.035 | +0.092 | 0.29 |
| GSE285029 | IFN-compact | 59 / 59 | +0.796 | +0.356 | **8.6×10⁻⁴** |
| GSE285029 | MHC-I | 59 / 59 | +0.572 | +0.411 | **1.2×10⁻⁴** |
| TCGA-LUAD | Wolf IFNγ | 111 / 111 | −0.002 | +0.003 | 0.97 |
| CPTAC-LUAD-protein | IFNG-6 RNA | 20 / 20 | −0.222 | −0.245 | 0.19 |
| CPTAC-LSCC-protein | Hallmark IFNγ | 20 / 20 | −3.54 | −0.570 | **0.0021** |

### OS (only GSE218989 + CPTAC)

| Cohort | n_Q1 / n_Q4 | events | HR (Q4 vs Q1) | 95% CI | log-rank p |
|---|---:|---:|---:|---|---:|
| GSE218989 | 89 / 89 | 64 / 69 | 1.32 | 0.93–1.85 | 0.12 |
| CPTAC-LUAD-protein | 19 / 19 | 4 / 4 | 1.85 | 0.42–8.09 | 0.41 |
| CPTAC-LSCC-protein | 17 / 17 | 5 / 4 | 0.76 | 0.21–2.81 | 0.68 |

KM median OS on GSE218989: Q1 471 d vs Q4 379 d. CPTAC medians mostly not reached. Underpowered protein OS.

---

## TACSTD2 Q4 vs Q1

### CD8A

| Cohort | n | δ_median | δ_rb | p |
|---|---:|---:|---:|---:|
| GSE218989 | 89 / 89 | −0.900 | −0.365 | **2.6×10⁻⁵** |
| GSE285029 | 59 / 59 | −0.273 | −0.117 | 0.28 |
| TCGA-LUAD | 126 / 126 | −0.287 | −0.140 | 0.054 |
| TCGA-LUSC | 126 / 126 | −0.766 | −0.294 | **5.5×10⁻⁵** |
| CPTAC-LUAD-protein | 28 / 28 | −0.621 | −0.255 | 0.10 |
| CPTAC-LSCC-protein | 27 / 27 | −0.420 | −0.130 | 0.42 |
| GSE68465 | 111 / 111 | −81.3 | −0.244 | **0.0017** |
| GSE31210 | 57 / 57 | −272 | −0.435 | **6.3×10⁻⁵** |
| OncoSG | 43 / 42 | −1.02 | −0.478 | **1.5×10⁻⁴** |

### ImmuneScore / GEP18 (short)

| Cohort | ImmuneScore δ_rb (p) | GEP18 δ_rb (p) |
|---|---|---|
| GSE218989 | −0.384 (**1.0×10⁻⁵**) | — |
| GSE285029 | −0.073 (0.50) | −0.007 (0.95) |
| TCGA-LUAD | −0.019 (0.79) | +0.017 (0.82) |
| TCGA-LUSC | −0.071 (0.33) | −0.220 (**0.0026**) |
| CPTAC-LUAD-protein | −0.230 (0.14) | −0.253 (0.11) |
| CPTAC-LSCC-protein | −0.010 (0.96) | +0.004 (0.99) |
| GSE68465 | −0.029 (0.71) | −0.150 (0.054) |
| GSE31210 | −0.420 (**1.1×10⁻⁴**) | −0.361 (**9.0×10⁻⁴**) |
| OncoSG | −0.503 (**6.7×10⁻⁵**) | −0.553 (**1.2×10⁻⁵**) |

TACSTD2 IFN / MHC-I: all p ≥ 0.16. OS: GSE218989 HR 1.28 (0.92–1.79) p=0.15;
CPTAC protein OS null (few events).

---

## How to read this vs the continuous cloud

The prior tables already showed negative Spearmans that look like a smear.
Q4 vs Q1 is the contrast that matches a “high vs low” claim:

- **Surgical LUAD microarray** (GSE31210, GSE68465): CLDN4 Q4 is clearly
  CD8 / ImmuneScore / GEP18 low (δ_rb −0.30 to −0.53).
- **ICI bulk GSE218989 (n=355):** CLDN4 and TACSTD2 Q4 are CD8-low and
  ESTIMATE-immune-low; IFN and MHC-I stay **null** (same as the continuous
  result on that table).
- **ICI bulk GSE285029 (n=234):** TACSTD2 Q4 is null. CLDN4 Q4 sits on
  **high** IFN / MHC-I / GEP18 — do not pool this cohort with the
  immune-low surgical sets.
- **TCGA-LUSC:** TACSTD2 Q4 is CD8 / GEP18 low; CLDN4 Q4 is not.
- **CPTAC protein:** LSCC CLDN4 protein Q4 is immune-low (n=20/20, large
  δ_rb). LUAD protein CLDN4 is the same sign, n=20/20, not significant.
- **OncoSG:** TACSTD2-only; Q4 is CD8 / GEP / A1-immune low (public z-scores).

---

## Figures (Q4 vs Q1 only; no scatter cloud)

| File | Content |
|---|---|
| `figures/fig1_forest_cldn4.png` | CLDN4 rank-biserial forest, all endpoints except OS |
| `figures/fig2_forest_tacstd2.png` | TACSTD2 forest |
| `figures/fig3_forest_os.png` | OS HR forest (3 cohorts with survival columns) |
| `figures/fig4_violins_cldn4_cd8a.png` | CLDN4 Q1 vs Q4 violins — CD8A |
| `figures/fig5_violins_cldn4_immunescore.png` | CLDN4 — ImmuneScore |
| `figures/fig6_violins_cldn4_gep18.png` | CLDN4 — GEP18 |
| `figures/fig7_violins_tacstd2_cd8a.png` | TACSTD2 — CD8A |
| `figures/fig8_violins_tacstd2_immunescore.png` | TACSTD2 — ImmuneScore |
| `figures/fig9_violins_cldn4_ifn.png` | CLDN4 — IFN |
| `figures/fig10_violins_cldn4_mhci.png` | CLDN4 — MHC-I (2 cohorts) |

Full numeric dump: `tables/q4_vs_q1.tsv`.

---

## Reproduce

```
pip install -r methods/quantile_cldn4_immune/requirements.txt
python methods/quantile_cldn4_immune/analyze.py
```

Reads only `harvested/`. Does not touch GEO/Xena/CPTAC.

---

## What this is not

- Not a new download of GSE218989 / GSE285029 / TCGA / CPTAC / OncoSG.
- Not a recompute of GEP18 / ESTIMATE / IFN from raw counts.
- Not a purity-residual or histology split (GSE218989 has no public histology).
- Not an ICI-response OR (GSE218989 response was already tested on that table
  and is null for CLDN4 / TACSTD2).
- Not a claim that CLDN4-high is immune-low **everywhere** — GSE285029
  contradicts that for IFN / MHC-I / GEP18.

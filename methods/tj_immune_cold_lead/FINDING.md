# PAPER FUNNEL: TJ pathway → then which gene tracks immune-cold?

**Never fabricate.** Every ρ / ratio / rank below is written by
`analyze_concordant4.py`, `analyze_cosmx.py`, or `analyze_tcga.py`, or copied
from a locked prior PR with the source named.

Pre-specified panel after Tacstd2 DEG / TJ up: **CLDN4 vs CLDN3, CLDN7, OCLN, F11R, CDH1**.

## Slide line (honest)

> **TJ first, then CLDN4 as the concordant-4 lead.** Malignant CLDN4 %pos is the most negative keratin-partial Spearman vs patient T/NK (ρ = −0.478, p = 2.7×10⁻⁴, I² = 0%). CosMx can only score CLDN4 vs CDH1 (CLDN3/7/OCLN/F11R off the 960 panel); on the same Q4/Q1 cut CDH1 is colder (8/8) and CLDN4 is not. TCGA keratin-partial CD8 ranks CLDN4 **5th of 6**. Public multimodal data do **not** uniquely pin CLDN4 over every TJ gene outside concordant-4.

## 0. TJ context (locked prior; not recomputed)

GSE137244 KL vs KP n=5 vs 5 (PR #694 / handoff): Tacstd2 Δ = +3.24, Cldn4 Δ = +5.57, TJ7 mean Δ = +3.27, MW p = 0.00794. That is the “TJ after Tacstd2 DEG” step. This PR does not re-run GSE137244.

## 1. Concordant-4 (n = 65) — **CLDN4 leads**

Malignant % positive vs same-patient T/NK fraction. Keratin = mean log1p of KRT8/18/19 on malignant cells, partialled on both sides. Calibration vs locked PR #644 units: max |Δ CLDN4 %| ≈ 0; unadjusted CLDN4 meta ρ = −0.531 (locked −0.531).

| rank | gene | KRT-partial ρ | p | I² |
|---:|---|---:|---:|---:|
| **1** | **CLDN4** | **−0.478** | **2.7×10⁻⁴** | **0%** |
| 2 | CLDN7 | −0.352 | 0.10 | 56% |
| 3 | CLDN3 | −0.346 | 0.012 | 0% |
| 4 | CDH1 | −0.277 | 0.073 | 16% |
| 5 | OCLN | −0.157 | 0.27 | 0% |
| 6 | F11R | −0.070 | 0.63 | 0% |

Label-swap head-to-head (CLDN4 − other, KRT, %pos): separable from **F11R** (Δ = −0.412, perm p = 0.040). Not separable from CLDN7 (p = 0.35), CLDN3 (p = 0.22), OCLN (p = 0.16), or CDH1 (p = 0.26).

## 2. CosMx He2022 — panel-limited; **CDH1 colder than CLDN4 on same cut**

| gene | on 960 / in h5ad | scored |
|---|---|---|
| CLDN4 | yes | yes |
| CDH1 | yes | yes |
| CLDN3, CLDN7, OCLN, F11R | **no** | **OFF — not fabricated** |

Primary cut = Q4 vs Q1 among author-tumor cells; CD8+NK neighbors at 50 µm (0.18 µm/px; median NN ≈ 7.8 µm).

| gene | ratio @ 50 µm | sections Δ&lt;0 | donors Δ&lt;0 | 8/8 & 5/5? |
|---|---:|---|---|---|
| **CDH1** | **0.786** | **8/8** | **5/5** | **yes** |
| CLDN4 | 0.800 | 4/8 | 3/5 | no |

Locked CLDN4 cytotoxic ratios **0.36 / 0.52** (50/100 µm, 8/8, 5/5, sign P = 0.031) are **not** recomputed or replaced here. That lock uses a different contrast than this same-cut Q4/Q1 head-to-head.

## 3. TCGA lung + keratin funnel — **CLDN4 does not lead**

Keratin-partial Spearman vs CD8 score (mean CD8A/CD8B); covariates KRT8+KRT18+KRT19. Cohorts: LUAD, LUSC, BRCA, CESC, KIRC, STAD, BLCA, PAAD.

| rank | gene | meta ρ | n_neg / 8 | p | I² |
|---:|---|---:|---|---:|---:|
| 1 | F11R | −0.131 | 7/8 | 0.0010 | 83% |
| 2 | CDH1 | −0.126 | 7/8 | 2.2×10⁻⁶ | 61% |
| 3 | CLDN7 | −0.109 | 7/8 | 0.00094 | 75% |
| 4 | OCLN | −0.101 | 6/8 | 0.00015 | 61% |
| **5** | **CLDN4** | **−0.081** | **7/8** | **0.0062** | **68%** |
| 6 | CLDN3 | −0.060 | 5/8 | 0.28 | 91% |

Matches the direction of PR #593 for CLDN4/CLDN7; extends the panel. Bulk surface ranking still does not pin CLDN4.

## PPT paste block

```
TJ after Tacstd2 DEG (GSE137244): TJ7 Δ=+3.27; Cldn4 Δ=+5.57 (locked)

Then immune-cold gene ranking (CLDN4 vs CLDN3/7/OCLN/F11R/CDH1):
  Concordant-4  → CLDN4 lead  ρ=−0.478 (I²=0%); not > CLDN7 by perm test
  CosMx         → only CLDN4 & CDH1 on panel; CDH1 Q4/Q1 colder (8/8); CLDN3/7/OCLN/F11R OFF
  TCGA          → CLDN4 rank 5/6 (F11R/CDH1/CLDN7/OCLN colder)

Lead for the slide: CLDN4 from concordant-4 patient T/NK.
Do not say CosMx or TCGA uniquely pin CLDN4 over the full TJ panel.
```

## Reproduce

```bash
# GEO matrices → /tmp/geo_c4 (see download notes in README)
GEO_DIR=/tmp/geo_c4 python3 methods/tj_immune_cold_lead/analyze_concordant4.py
TCGA_CACHE=/tmp/tcga python3 methods/tj_immune_cold_lead/analyze_tcga.py
# CosMx h5ad via download_cosmx_nsclc_h5ad.py → methods/data/cosmx_nsclc/
python3 methods/tj_immune_cold_lead/analyze_cosmx.py
```

Private 8-KL / KD co-culture are out of scope (handoff). They are what still functionally pin CLDN4 when public ranking is mixed.

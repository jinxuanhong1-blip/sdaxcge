# DepMap lung lines — CLDN4 vs Hallmark IFN-γ / MHC-I / CD274

**Additive only.** TACSTD2–CLDN4 co-expression on DepMap/CCLE lung RNA is treated as **already known** and is **not** restated as the finding. This page adds the **IFN / MHC-I / CD274** axis on the same public matrix.

**Data.** DepMap Public 24Q4 (Figshare+ [10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1)): `OmicsExpressionProteinCodingGenesTPMLogp1.csv` = log2(TPM+1); `Model.csv` for lineage. Hallmark IFN-γ gene set: MSigDB 2024.1 Hs `HALLMARK_INTERFERON_GAMMA_RESPONSE` (**198/200** genes present; `RIGI`, `TMT1B` absent).

**Honest n.** `OncotreeLineage == Lung` and `ModelType == Cell Line` with complete RNA for CLDN4 and the scored genes: **n = 214** (of 260 lung models in `Model.csv`; 46 have no 24Q4 RNA and are excluded). Split: LUAD 80 / LUSC 27 / other NSCLC 36 / SCLC+NET 60 / other lung 11. NSCLC (Oncotree primary disease) **n = 143**. These are cultured lines. There is **no T-cell infiltrate** and **no IFN-γ treatment**.

**Scatter.** `figures/fig_cldn4_vs_ifng.png` (primary) and `figures/fig_cldn4_vs_ifn_mhc_cd274.png` (three additive axes).

## Finding (IFN / MHC axis)

In cancer-cell **basal** RNA, CLDN4 is **positively** associated with Hallmark IFN-γ, MHC-I, and CD274. The sign is the opposite of “CLDN4-high = IFN-cold” at the cell-intrinsic level. MHC-I is the weakest of the three; in NSCLC the MHC-I 95% CI includes 0.

| cohort | n | CLDN4 vs Hallmark IFN-γ ρ (p, q) | CLDN4 vs MHC-I ρ (p, q) | CLDN4 vs CD274 ρ (p, q) |
|---|---:|---|---|---|
| **all lung (primary)** | **214** | **+0.284** (2.5×10⁻⁵, q=3.7×10⁻⁵); CI [0.16, 0.41] | **+0.217** (0.0014, q=0.0014); CI [0.08, 0.34] | **+0.330** (7.8×10⁻⁷, q=2.3×10⁻⁶); CI [0.20, 0.45] |
| NSCLC | 143 | +0.237 (0.0044, q=0.0066); CI [0.07, 0.39] | +0.166 (0.047, q=0.047); CI [−0.01, 0.33] | +0.275 (0.00087, q=0.0026); CI [0.10, 0.43] |
| LUAD | 80 | +0.249 (0.026) | +0.109 (0.34) | +0.350 (0.0015) |
| LUSC | 27 | +0.183 (0.36) | +0.123 (0.54) | +0.003 (0.99) |
| SCLC+NET | 60 | +0.356 (0.0052) | +0.293 (0.023) | +0.257 (0.048) |

Spearman, two-sided, complete cases. Hallmark IFN-γ and MHC-I (HLA-A, HLA-B, HLA-C, B2M; 4/4) are the mean of **within-cohort** gene-wise z-scores. BH *q* is within the three new axes (IFN-γ, MHC-I, CD274) on that cohort. Bootstrap 95% CI: 5,000 resamples, seed 0, all-lung and NSCLC only. LUSC n=27 is underpowered; do not read the LUSC nulls as a histology-specific rule.

All-lung ρ is slightly higher than NSCLC-only because SCLC/NET sit toward the low-CLDN4 / low-IFN / low-MHC-I corner. That mix is why both n=214 and n=143 are reported. Do not cite the all-lung number as “NSCLC.”

## What this is not

- **Not** the TACSTD2–CLDN4 correlation. That pair was scored only as a companion (same n=214 complete cases). It is not the finding.
- **Not** a hot/cold tumor, ICI-response, or IFN-stimulation result. No immune cells are in the dish.
- **Not** surface MHC or PD-L1 protein. CD274 here is RNA.
- **Not** a claim that CLDN4 *induces* IFN. Association only.
- Hallmark IFN-γ includes genes a pure cancer line may not use (e.g. *GZMA*). The MHC-I cassette is cleaner and weaker.

## How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib
python3 methods/depmap_cldn4_ifn/download.py --outdir data/depmap_cldn4_ifn
python3 methods/depmap_cldn4_ifn/analyze.py --data data/depmap_cldn4_ifn --outdir methods/depmap_cldn4_ifn
```

Tables: `tables/correlations.tsv`, `tables/lung_cell_lines.tsv`, `tables/cohort_counts.tsv`, `tables/key_stats.json`.

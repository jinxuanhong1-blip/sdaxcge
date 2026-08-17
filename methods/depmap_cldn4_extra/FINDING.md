# Extra cut: DepMap lung CLDN4 Chronos vs IFN / MHC-I / CD274

**Additive only. New cut only.** The all-lung CLDN4-**RNA** Spearman (n=214, Hallmark IFN-γ ρ=+0.28, CD274 ρ=+0.33) is already reported in `methods/depmap_cldn4_ifn` and is **not** restated as the finding. This page is the CRISPR∩RNA cut: CLDN4 **dependency** (Chronos) versus the same three axes, plus CLDN4 RNA on that smaller honest n.

**Data.** DepMap Public 24Q4 (Figshare+ [10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1)). RNA = `OmicsExpressionProteinCodingGenesTPMLogp1.csv` log2(TPM+1). Dependency = `CRISPRGeneEffect.csv` Chronos (more negative = stronger dependency). Hallmark IFN-γ: MSigDB 2024.1 Hs, **198/200** genes (`RIGI`, `TMT1B` missing). MHC-I = mean-z of HLA-A/B/C + B2M (4/4). Signatures are the mean of **within-cut** gene-wise z-scores.

**Honest n.** `OncotreeLineage == Lung`, `ModelType == Cell Line`: **260** models in `Model.csv` → **214** with complete RNA for CLDN4 / CD274 / MHC-I → **91** of those have no finite CLDN4 Chronos → **n = 123** CRISPR∩RNA (the new cut). NSCLC **n = 95** (LUAD 51 / LUSC 20 / other NSCLC 24). SCLC+NET n=26; other lung n=2. Cultured lines: no infiltrate, no IFN treatment.

**CLDN4 is not a lung essential.** On n=123, Chronos median **+0.050** (IQR −0.036 to +0.107; range −0.29 to +0.39). **0 / 123** lines have Chronos < −0.5. There is almost no CLDN4-dependency variance to correlate.

## Extra table (this cut)

Master file: `tables/extra_table.tsv`. Spearman, two-sided, complete cases. BH *q* is within the three axes on that cut × predictor. Bootstrap 95% CI: 5,000 resamples, seed 0.

| cut | predictor | target | n | ρ | p | q | 95% CI |
|---|---|---|---:|---:|---:|---:|---|
| CRISPR∩RNA lung | CLDN4 Chronos | Hallmark IFN-γ | **123** | **+0.077** | 0.397 | 0.595 | [−0.10, +0.25] |
| CRISPR∩RNA lung | CLDN4 Chronos | MHC-I | 123 | +0.098 | 0.279 | 0.595 | [−0.08, +0.27] |
| CRISPR∩RNA lung | CLDN4 Chronos | CD274 | 123 | −0.038 | 0.679 | 0.679 | [−0.22, +0.14] |
| CRISPR∩RNA lung | CLDN4 RNA | Hallmark IFN-γ | 123 | +0.330 | 1.9×10⁻⁴ | 2.9×10⁻⁴ | [+0.16, +0.48] |
| CRISPR∩RNA lung | CLDN4 RNA | MHC-I | 123 | +0.257 | 0.0041 | 0.0041 | [+0.08, +0.42] |
| CRISPR∩RNA lung | CLDN4 RNA | CD274 | 123 | +0.338 | 1.3×10⁻⁴ | 2.9×10⁻⁴ | [+0.16, +0.50] |
| CRISPR∩RNA NSCLC | CLDN4 Chronos | Hallmark IFN-γ | **95** | +0.106 | 0.305 | 0.725 | [−0.10, +0.31] |
| CRISPR∩RNA NSCLC | CLDN4 Chronos | MHC-I | 95 | +0.037 | 0.725 | 0.725 | [−0.16, +0.24] |
| CRISPR∩RNA NSCLC | CLDN4 Chronos | CD274 | 95 | +0.046 | 0.656 | 0.725 | [−0.16, +0.25] |
| CRISPR∩RNA NSCLC | CLDN4 RNA | Hallmark IFN-γ | 95 | +0.323 | 0.0014 | 0.0042 | [+0.12, +0.50] |
| CRISPR∩RNA NSCLC | CLDN4 RNA | MHC-I | 95 | +0.236 | 0.021 | 0.021 | [+0.03, +0.43] |
| CRISPR∩RNA NSCLC | CLDN4 RNA | CD274 | 95 | +0.301 | 0.0030 | 0.0045 | [+0.08, +0.50] |

## Finding (dependency, not the n=214 RNA ρ)

CLDN4 **Chronos is null** versus basal Hallmark IFN-γ, MHC-I, and CD274 on the CRISPR∩RNA lung set (n=123) and on NSCLC (n=95). Every Chronos 95% CI includes 0. That is the extra result: the already-reported RNA association is **not** a CLDN4-dependency association.

CLDN4 **RNA on this new n=123** still tracks IFN-γ (ρ=+0.330) and CD274 (ρ=+0.338). Sign matches the reported n=214 ρ=+0.28/+0.33; this is the screened subset, not a second discovery of the same all-lung Spearman. MHC-I is again the weakest of the three RNA axes.

Q4 vs Q1 of CLDN4 RNA on the same n=123 (31 vs 31; thresholds log2(TPM+1) 7.92 / 2.37): IFN-γ median-z +0.14 vs −0.19 (Cliff’s δ=+0.51, p=5.9×10⁻⁴, q=7.7×10⁻⁴); MHC-I +0.41 vs −0.28 (δ=+0.50, p=7.7×10⁻⁴); CD274 3.72 vs 0.67 (δ=+0.56, p=1.8×10⁻⁴). Same direction as the overlap Spearman. Table: `tables/extra_q4q1.tsv`.

Histology splits on Chronos are underpowered (LUAD n=51, LUSC n=20, SCLC+NET n=26) and stay in `tables/extra_correlations.tsv`. Do not read those nulls as a histology rule.

Scatter: `figures/fig_cldn4_chronos_vs_ifn_mhc_cd274.png`.

## What this is not

- **Not** the n=214 all-lung CLDN4-RNA ρ=+0.28/+0.33. That cut is already reported.
- **Not** a claim that CLDN4 knockout changes IFN/MHC/CD274. Chronos here is basal fitness, and CLDN4 is not essential in these lines.
- **Not** a hot/cold tumor or ICI-response result. No immune cells are in the dish.
- **Not** surface MHC or PD-L1 protein. CD274 is RNA.
- A near-zero Chronos distribution has little dynamic range. The Chronos null is compatible with “no CLDN4 dependency to correlate,” not with a precise zero biological effect.

## How to rerun

```bash
python3 -m pip install pandas numpy scipy matplotlib
python3 methods/depmap_cldn4_extra/download.py --outdir data/depmap_cldn4_extra
python3 methods/depmap_cldn4_extra/analyze.py --data data/depmap_cldn4_extra --outdir methods/depmap_cldn4_extra
```

Tables: `tables/extra_table.tsv` (PR extra table), `tables/extra_correlations.tsv`, `tables/extra_q4q1.tsv`, `tables/crispr_rna_lung_lines.tsv`, `tables/cohort_counts.tsv`, `tables/key_stats.json`.

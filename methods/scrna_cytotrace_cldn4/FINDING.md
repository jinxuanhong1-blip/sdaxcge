# Finding — CytoTRACE-like potency vs CLDN4 in GSE131907 epithelium

Additive slice under `methods/scrna_cytotrace_cldn4/`. **CLDN4 only.** Public **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277). **CytoTRACE2 was not run.** Potency = residual `n_genes` after OLS on `log1p(UMI)`, rank-scaled to `[0, 1]` within the epithelial universe of that contrast (Gulati 2020 gene-count idea). Primary unit = **sample**. Cell-level ρ is exploratory (pseudoreplication).

GSE207422 was not re-run here. This atlas is treatment-naive LUAD; it has **no ICI / MPR labels**.

## Verdict

| Contrast | Honest n | Stat | p |
|---|---|---|---|
| **Primary:** tLung sample-mean CLDN4 vs CytoTRACE-like | **11** (all tLung samples; min 42 epi cells) | ρ=**+0.155** | **0.65** |
| Sensitivity: tumor-site sample-mean CLDN4 vs CytoTRACE-like | **36** (≥20 epi; 1 PE sample dropped) | ρ=**−0.113** | **0.51** |
| nLung sample-mean CLDN4 vs CytoTRACE-like | **11** | ρ=**+0.464** | **0.15** |
| All-epithelial sample-mean (mixes nLung + tumor) | **47** | ρ=**+0.140** | **0.35** |
| Extra: paired CytoTRACE-like, CLDN4 Q4 vs Q1, tLung | **11** | W=6; median Δ=+0.122 | **0.014** |
| Extra: paired CytoTRACE-like, CLDN4 Q4 vs Q1, tumor sites | **35** | W=95; median Δ=+0.041 | **5.4×10⁻⁴** |

**What holds:** Between-sample CLDN4 does **not** track CytoTRACE-like potency. Primary tLung n=11 is underpowered for a moderate ρ and is a **null**. The better-powered tumor-site test (n=36) is also **null** and opposite in sign.

**What is only within-tumor:** CLDN4-high epithelial cells are slightly **more** CytoTRACE-like than CLDN4-low cells from the same sample (tLung Δ median +0.12; tumor-site Δ +0.041). That is a small paired shift, not a sample-level stemness association. Do not promote it to “CLDN4 marks potency.”

**What does not hold:** A sample-level CLDN4–potency correlation in GSE131907 epithelium. Cell-level p-values (tLung n_cells=7,270, ρ=+0.139) are not confirmatory.

## Honest n

Public UMI matrix: **208,506** barcodes × **29,634** genes. Author QC already has UMI≥1,000 and n_genes≥200; the UMI≥200 filter drops **0** cells.

| Universe | Cells | Samples with ≥20 epi | Notes |
|---|---:|---:|---|
| All cells | 208,506 | 58 | 44 patients |
| Author epithelial (`Cell_type == Epithelial cells`) | 36,467 | 47 | includes nLung |
| **tLung epithelial (primary)** | **7,270** | **11 / 11** | min 42, median 316, max 2,486 |
| Tumor-site epithelial (tLung + tL/B + mLN + mBrain + PE) | 32,764 | **36** | PE 4, mBrain 10, mLN 7, tL/B 4, tLung 11 |
| nLung epithelial | 3,703 | 11 | normal-lung contrast, not mixed into primary |
| Excluded PE | 15 | EFFUSION_64 | <20 epi; not in n=36 |
| Paired tumor-site quartiles | — | **35** | EFFUSION_12 has 25 epi (Q1/Q4 n=7 < 8) |

Origin-stratified sample-level ρ (potency ranked in all tumor epithelium; all NS): tLung n=11 ρ=+0.127 p=0.71; mBrain n=10 ρ=+0.333 p=0.35; mLN n=7 ρ=+0.321 p=0.48; tL/B n=4 ρ=−0.400 p=0.60; PE n=4 ρ=−0.400 p=0.60. n=4 is not a test.

Two tumor samples have essentially no CLDN4 (EBUS_13 0.26% positive; EFFUSION_13 0%). They stay in the n=36 table.

## Methods

1. Stream the public GEO UMI TSV. Keep **CLDN4** (EPCAM for QC only). Accumulate per-cell `total_UMI` and `n_genes` (UMI>0). Full dense matrix is not written.
2. Join author cell annotation. Epithelial = `Cell_type == Epithelial cells`. No marker re-annotation. No CopyKAT.
3. CLDN4 = `log1p(CP10k)` from raw UMI and the streamed totals.
4. Potency, per contrast universe: OLS `n_genes ~ log1p(total_UMI)`; residual rank-scaled to `[0, 1]`. Higher = more genes than expected for depth = more stem-like. **Not `cytotrace2()`.**
5. Sample mean of those cell scores. Spearman two-sided. Primary = tLung. Sensitivity = tumor sites. nLung is a contrast, not a discovery set.
6. Extra paired Wilcoxon: within-sample CLDN4 top vs bottom quartile, min 8 cells per stratum.

No TACSTD2 test. No cycling / keratin / MPR split in this slice.

## What this does not claim

- It does not claim CytoTRACE2 output. The R package was not executed.
- It does not claim a sample-level CLDN4–stemness correlation.
- It does not claim ICI / MPR biology. GSE131907 is treatment-naive.
- It does not use GSE207422 numbers from other modules as if they were computed here.
- Cell-level p-values are not the claim. n_cells is large by construction.
- Within-sample Q4>Q1 on potency is a small paired delta, not a between-tumor effect.

## Outputs

- `results/stats.tsv` — every test (n / stat / p)
- `results/summary.json`
- `results/per_sample_tLung_epithelial.tsv`
- `results/per_sample_tumor_epithelial.tsv`
- `results/per_sample_nLung_epithelial.tsv`
- `results/per_sample_all_epithelial.tsv`
- `results/figures/gse131907_cldn4_vs_potency.png`
- `results/figures/gse131907_cldn4_quartile_potency.png`

## Reproduce

```bash
pip install -r methods/scrna_cytotrace_cldn4/requirements.txt
bash methods/scrna_cytotrace_cldn4/scripts/run_all.sh
```

Matrices stay in `/tmp/gse131907/` and are not committed.

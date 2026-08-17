# FINDING — winning-pair DoRothEA TF activity (CLDN4-only)

ADDITIVE **CLDN4 only** on the PR #320 winning pair **GSE131907 + GSE205335**.
No TACSTD2∩CLDN4 dual-high. **GSE207422 and GSE148071 are not merged in.**
Patient is the unit. p-values are descriptive.

decoupleR / dorothea R were **not installed**. Activity is documented
DoRothEA **wmean** (Badia-i-Mompel et al. 2022):
wmean = sum(w_g * x_g) / sum(|w_g|) on `log1p(CP10k)` targets,
w = +1 / -1 from OmniPath A+B+C. No VIPER NES.

## Honest n

| Item | n | Note |
|---|---:|---|
| GSE131907 patients with any tumor-origin malignant cells | 32 | author `Malignant cells` + tS1/tS2/tS3; PE/nLung/nLN out |
| GSE205335 patients with any author-malignant cells | 22 | `lineage.sub == Malignant cells` |
| **Paired Q4 vs Q1 (floor ≥20 malignant, ≥8/tail)** | **52** | **31 + 21** |
| GSE131907 eligible | 31 / 32 | dropped if a tail <8 |
| GSE205335 eligible | 21 / 22 | dropped if a tail <8 |
| Cells as n | 0 | not used |

Do **not** cite the attempted header n as the test n. Dropped: GSE131907
**P0009** (5 malignant cells) and GSE205335 **P4001** (27 cells; `n//4=6<8`).
Equal-count tails are used because zero-inflated CLDN4 collapses `pd.qcut`
bins. GSE205335 eligible histology is mixed: 13 ADC + 4 SCLC + 3 SQ + 1 NUT.

GSE131907 matrix: 208506 barcodes streamed,
31136 malignant tumor cells kept.
GSE205335 matrix: 96505 barcodes,
28512 author-malignant kept.
Shared DoRothEA A+B+C targets in both matrices: 8780.
TFs scored (≥5 shared targets): 297.

Absent from DoRothEA A+B+C (not imputed): CIITA, NLRC5, GRHL1, GRHL3, OVOL1, OVOL2, RFXANK, RFXAP.

## Verdict (pooled Q4 vs Q1)

- **IFN / STAT1:** n=52 paired patients, Δ=+0.046, high>low 50/52, p=1.18e-09.
- **IFN∩MHC / IRF1:** n=52 paired patients, Δ=+0.062, high>low 50/52, p=2.63e-08.
- **MHC / RFX5:** n=52 paired patients, Δ=+0.049, high>low 48/52, p=3.24e-08.
- **TJ / GRHL2:** n=52 paired patients, Δ=+0.122, high>low 52/52, p=3.50e-10.
- **keratin / TP63:** n=52 paired patients, Δ=+0.063, high>low 52/52, p=3.50e-10.

This is TF activity in the **same malignant cells** that define the CLDN4
split, not a T/NK test. The PR #320 T/NK association on this pair is
given and is not re-ranked (Q4 vs Q1 n=23 r=−0.705 p=0.0003; continuous
n=43 ρ=−0.479).

## Primary TF table — pooled patients, Q4 vs Q1

| program | TF | n_targets | n | high vs low | Δ mean | high>low | W | p |
|---|---|---:|---:|---|---:|---:|---:|---|
| IFN | STAT1 | 1091 | 52 | 0.175 vs 0.129 | +0.046 | 50/52 | 21 | **1.18e-09** |
| IFN | STAT2 | 104 | 52 | 0.210 vs 0.154 | +0.056 | 50/52 | 90 | **4.90e-08** |
| IFN | IRF1 | 79 | 52 | 0.217 vs 0.155 | +0.062 | 50/52 | 78 | **2.63e-08** |
| IFN | IRF2 | 74 | 52 | 0.275 vs 0.213 | +0.062 | 50/52 | 50 | **5.91e-09** |
| IFN | IRF3 | 70 | 52 | 0.206 vs 0.172 | +0.034 | 48/52 | 109 | **1.28e-07** |
| IFN | IRF7 | 0 | 0 | NA | NA | NA | NA | not scored |
| IFN | IRF8 | 8 | 52 | 0.094 vs 0.059 | +0.035 | 48/52 | 13 | **7.45e-10** |
| IFN | IRF9 | 46 | 52 | 0.241 vs 0.175 | +0.066 | 50/52 | 24 | **1.39e-09** |
| IFN | STAT3 | 173 | 52 | 0.302 vs 0.208 | +0.093 | 52/52 | 0 | **3.50e-10** |
| MHC | RFX5 | 87 | 52 | 0.316 vs 0.267 | +0.049 | 48/52 | 82 | **3.24e-08** |
| TJ | GRHL2 | 74 | 52 | 0.295 vs 0.173 | +0.122 | 52/52 | 0 | **3.50e-10** |
| TJ | KLF4 | 48 | 52 | 0.295 vs 0.225 | +0.070 | 52/52 | 0 | **3.50e-10** |
| TJ | ELF3 | 52 | 52 | 0.233 vs 0.147 | +0.085 | 52/52 | 0 | **3.50e-10** |
| TJ | TFAP2A | 126 | 52 | 0.196 vs 0.135 | +0.061 | 52/52 | 0 | **3.50e-10** |
| keratin | TP63 | 86 | 52 | 0.147 vs 0.083 | +0.063 | 52/52 | 0 | **3.50e-10** |
| keratin | KLF5 | 78 | 52 | 0.418 vs 0.248 | +0.170 | 52/52 | 0 | **3.50e-10** |
| keratin | SOX2 | 391 | 52 | 0.106 vs 0.073 | +0.033 | 52/52 | 0 | **3.50e-10** |

High = within-patient malignant CLDN4 Q4; low = Q1. Mean is the mean of
patient-bin wmean values. **Bold p** is p<0.05 (descriptive).

### Same table by cohort

**GSE131907 (n=31)**

| program | TF | n_targets | n | high vs low | Δ mean | high>low | W | p |
|---|---|---:|---:|---|---:|---:|---:|---|
| IFN | STAT1 | 1091 | 31 | 0.173 vs 0.131 | +0.042 | 31/31 | 0 | **9.31e-10** |
| IFN | STAT2 | 104 | 31 | 0.197 vs 0.142 | +0.056 | 31/31 | 0 | **9.31e-10** |
| IFN | IRF1 | 79 | 31 | 0.213 vs 0.153 | +0.060 | 31/31 | 0 | **9.31e-10** |
| IFN | IRF2 | 74 | 31 | 0.268 vs 0.210 | +0.058 | 31/31 | 0 | **9.31e-10** |
| IFN | IRF3 | 70 | 31 | 0.198 vs 0.165 | +0.033 | 30/31 | 6 | **1.30e-08** |
| IFN | IRF7 | 0 | 0 | NA | NA | NA | NA | not scored |
| IFN | IRF8 | 8 | 31 | 0.097 vs 0.063 | +0.034 | 28/31 | 8 | **2.33e-08** |
| IFN | IRF9 | 46 | 31 | 0.239 vs 0.178 | +0.061 | 31/31 | 0 | **9.31e-10** |
| IFN | STAT3 | 173 | 31 | 0.307 vs 0.223 | +0.084 | 31/31 | 0 | **9.31e-10** |
| MHC | RFX5 | 87 | 31 | 0.330 vs 0.278 | +0.052 | 29/31 | 7 | **1.77e-08** |
| TJ | GRHL2 | 74 | 31 | 0.287 vs 0.188 | +0.100 | 31/31 | 0 | **9.31e-10** |
| TJ | KLF4 | 48 | 31 | 0.280 vs 0.218 | +0.062 | 31/31 | 0 | **9.31e-10** |
| TJ | ELF3 | 52 | 31 | 0.227 vs 0.163 | +0.065 | 31/31 | 0 | **9.31e-10** |
| TJ | TFAP2A | 126 | 31 | 0.190 vs 0.140 | +0.050 | 31/31 | 0 | **9.31e-10** |
| keratin | TP63 | 86 | 31 | 0.134 vs 0.085 | +0.049 | 31/31 | 0 | **9.31e-10** |
| keratin | KLF5 | 78 | 31 | 0.405 vs 0.264 | +0.141 | 31/31 | 0 | **9.31e-10** |
| keratin | SOX2 | 391 | 31 | 0.098 vs 0.072 | +0.026 | 31/31 | 0 | **9.31e-10** |

**GSE205335 (n=21)**

| program | TF | n_targets | n | high vs low | Δ mean | high>low | W | p |
|---|---|---:|---:|---|---:|---:|---:|---|
| IFN | STAT1 | 1091 | 21 | 0.178 vs 0.127 | +0.051 | 19/21 | 10 | **4.10e-05** |
| IFN | STAT2 | 104 | 21 | 0.228 vs 0.171 | +0.056 | 19/21 | 33 | **0.00286** |
| IFN | IRF1 | 79 | 21 | 0.224 vs 0.160 | +0.064 | 19/21 | 27 | **0.00118** |
| IFN | IRF2 | 74 | 21 | 0.285 vs 0.217 | +0.068 | 19/21 | 16 | **1.61e-04** |
| IFN | IRF3 | 70 | 21 | 0.218 vs 0.183 | +0.035 | 18/21 | 38 | **0.00554** |
| IFN | IRF7 | 0 | 0 | NA | NA | NA | NA | not scored |
| IFN | IRF8 | 8 | 21 | 0.089 vs 0.053 | +0.035 | 20/21 | 1 | **1.91e-06** |
| IFN | IRF9 | 46 | 21 | 0.244 vs 0.170 | +0.074 | 19/21 | 11 | **5.25e-05** |
| IFN | STAT3 | 173 | 21 | 0.294 vs 0.186 | +0.108 | 21/21 | 0 | **9.54e-07** |
| MHC | RFX5 | 87 | 21 | 0.294 vs 0.250 | +0.045 | 19/21 | 29 | **0.0016** |
| TJ | GRHL2 | 74 | 21 | 0.308 vs 0.152 | +0.155 | 21/21 | 0 | **9.54e-07** |
| TJ | KLF4 | 48 | 21 | 0.318 vs 0.236 | +0.082 | 21/21 | 0 | **9.54e-07** |
| TJ | ELF3 | 52 | 21 | 0.240 vs 0.124 | +0.116 | 21/21 | 0 | **9.54e-07** |
| TJ | TFAP2A | 126 | 21 | 0.205 vs 0.128 | +0.077 | 21/21 | 0 | **9.54e-07** |
| keratin | TP63 | 86 | 21 | 0.165 vs 0.082 | +0.084 | 21/21 | 0 | **9.54e-07** |
| keratin | KLF5 | 78 | 21 | 0.438 vs 0.224 | +0.214 | 21/21 | 0 | **9.54e-07** |
| keratin | SOX2 | 391 | 21 | 0.118 vs 0.075 | +0.042 | 21/21 | 0 | **9.54e-07** |

## Companion gene-set modules (not the TF table)

Mean `log1p(CP10k)` of the listed genes. CLDN4 is excluded from TJ.

| module | n | high vs low | Δ mean | high>low | p |
|---|---:|---|---:|---:|---|
| IFN ISG core | 52 | 0.333 vs 0.272 | +0.061 | 46/52 | 4.98e-07 |
| MHC-I APM | 52 | 0.907 vs 0.769 | +0.138 | 48/52 | 2.09e-07 |
| TJ (no CLDN4) | 52 | 0.421 vs 0.209 | +0.212 | 52/52 | 3.50e-10 |
| keratin | 52 | 1.005 vs 0.773 | +0.233 | 51/52 | 4.97e-10 |
| OXPHOS control | 52 | 0.878 vs 0.841 | +0.037 | 33/52 | 0.0147 |

## Median-split companion

Same patients are not required. Eligible median-split n = 51 (29 GSE131907 + 22 GSE205335). See `results/tf_tests.tsv` (`split=median`).

## Honest limits

1. **Winning pair only.** GSE207422 and GSE148071 were not downloaded.
2. **No dual-high.** TACSTD2 is not a gate.
3. **wmean, not decoupleR NES.** No permutation / VIPER normalization.
4. **Author labels.** CopyKAT was not re-run. GSE131907 tLung uses tS1/tS2/tS3.
5. **GSE205335 histology mix** (13 ADC + 4 SCLC + 3 SQ + 1 NUT) is part of the honest n.
6. IRF7 is not scored (4 ABC targets, below the ≥5 floor). IRF8 has 8 targets.
7. The OXPHOS control also rises (33/52, Δ=+0.037, p=0.015). Depth is not fully ruled out; TJ/keratin/IFN module deltas are larger.
8. p-values are descriptive.

## Extra figures

- `figures/fig_extra_honest_n.png`
- `figures/fig_extra_paired_tf.png`
- `figures/fig_extra_forest.png`
- `figures/fig_extra_heatmap_delta.png`
- `figures/fig_extra_modules.png`

## Files

- `results/tf_table.tsv` — compact primary TF table (pooled Q4 vs Q1)
- `results/tf_tests.tsv` — all cohort × split tests
- `results/per_patient.tsv` — occupancy + bin activities
- `results/tf_meta.tsv` — targets used per TF
- `results/module_tests.tsv` — companion gene-set tests
- `resources/dorothea_hs_ABC.tsv` — OmniPath A+B+C weights
- `METHODS.md` — wmean formula, floors, labels

## Reproduce

```bash
python3 methods/winpair_131907_205335_dorothea_cldn4/scripts/download.py
python3 methods/winpair_131907_205335_dorothea_cldn4/scripts/analyze.py
```

# TISMO in vivo ICB pairing: Tacstd2, Cldn4, and TJ score

Public recompute from TISMO (Zeng et al., *Nucleic Acids Research* 2022, PMID 34534350) plus four GEO KP/KL/KPM anti-PD-1 series. Values are not fabricated: every count and *p* below is in `results/tismo/tables/` or `results/gemm/tables/`.

## Paper sentence

After ICB, **Tacstd2 rose in 49/64 TISMO models** (Wilcoxon signed-rank *p* = 5.84×10⁻⁵; two-sided binomial *p* = 2.44×10⁻⁵); **Cldn4 rose in 34/64** (28 down, 2 ties; Wilcoxon *p* = 0.077; binomial *p* = 0.53); **TJ score rose in 37/64** (27 down; Wilcoxon *p* = 0.087; binomial *p* = 0.26). **Lung subset** (TISMO lung-carcinoma ICB = LLC only, 2/64 pairs; LLC is not KL): Tacstd2 2/2 up, Cldn4 2/2 up, TJ 0/2 (both down). TISMO has no KL or KP ICB pairs. In four public KP/KL/KPM anti-PD-1 GEO contrasts, Cldn4 mean-delta was down in 4/4; Tacstd2 was mixed (2/4 up).

## Pairing (TISMO 64)

TISMO’s Gene module (`POST https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn`, `type=3`, all ICB treatments × all tumors) returns one named slice per `cell_line` label (cell line × study × condition × ICB regimen). That Tacstd2 export has **exactly 64 slices** — the July PPT pairing.

- **Naive / baseline:** `Baseline == 1` (isotype / vehicle / no treatment / matched control).
- **ICB:** `Baseline == 0` (anti-PD1, anti-PD-L1, anti-CTLA4, or combo; responders + non-responders pooled).
- **Delta:** ICB − naive. Primary location = **mean** (this is the rule that yields 49/64 and *p* = 5.84×10⁻⁵). Median deltas are reported alongside.
- **TJ score:** per-sample mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln on the Tacstd2 sample lock (mean 6.99 / 7 genes present), then the same pair rule.
- Cldn4 and TJ are locked to the same 64 Tacstd2 stems. Cldn4’s native export has one extra stem (`MOC22_RU31562203_antiPD1`), not used in the 64.

The 64 slices are **17 cell lines / 22 studies**. GSE124821 alone is 19/64. Cancer mix: mammary 29, melanoma 14, colorectal 10, gastric 5, liver 2, **lung 2**, sarcoma 2.

TISMO `cellLineMeta` / `vivoMeta` contain **no STK11, LKB1, KL, or KP** ICB models. KPB25L is mammary, not lung. **LLC is Lewis lung carcinoma, not KL.**

## All models (n = 64)

| Marker | Location | Up | Down | Tie | Wilcoxon *p* (two-sided) | Binomial *p* (two-sided) | Mean Δ |
|---|---|---:|---:|---:|---:|---:|---:|
| Tacstd2 | mean | **49** | 15 | 0 | **5.84×10⁻⁵** | 2.44×10⁻⁵ | +0.265 |
| Tacstd2 | median | 40 | 18 | 6 | 0.0071 | 0.0054 | +0.187 |
| Cldn4 | mean | **34** | 28 | 2 | **0.077** | 0.53 | +0.129 |
| Cldn4 | median | 35 | 21 | 8 | 0.084 | 0.081 | +0.140 |
| TJ | mean | **37** | 27 | 0 | **0.087** | 0.26 | +0.048 |
| TJ | median | 35 | 29 | 0 | 0.24 | 0.53 | +0.029 |

Wilcoxon is signed-rank on the 64 (or non-zero) paired deltas. Binomial is a sign test vs 0.5 after dropping ties.

Tacstd2 and Cldn4 are not interchangeable on this pairing: Tacstd2’s up-count is 49/64; Cldn4’s is 34/64 with a Wilcoxon *p* of 0.077. The TJ mean (seven junction genes) sits between them (37/64, *p* = 0.087).

Figures: `results/tismo/figures/waterfall_Tacstd2.png`, `waterfall_Cldn4.png`, `waterfall_TJ.png`; paired-mean slopes `paired_Tacstd2.png`, `paired_Cldn4.png`, `paired_TJ.png`. Lung LLC slices are outlined in black.

## Lung / KL / KP / LLC

### TISMO lung and KL/KP/LLC strata

| Stratum | What is in TISMO ICB | Tacstd2 | Cldn4 | TJ |
|---|---|---|---|---|
| Lung | LLC GSE155972 only (2 slices) | 2/2 up (Wilcoxon *p* = 0.50) | 2/2 up (*p* = 0.50) | 0/2 (*p* = 0.50) |
| KL / KP / LLC | same as LLC (no KL, no KP) | 2/2 up | 2/2 up | 0/2 |

LLC GSE155972 (anti-PD1 + anti-CTLA4; Griffin et al., *Nature* 2021):

| Slice | n naive / ICB | Tacstd2 Δmean | Cldn4 Δmean | TJ Δmean |
|---|---|---:|---:|---:|
| LLC WT | 10 / 6 | +0.133 | +0.081 | −0.148 |
| LLC Setdb1-KO | 10 / 7 | +0.725 | +0.334 | −0.142 |

Both single-gene means go up; the seven-gene TJ mean goes down (other junction genes offset Cldn4). Expression of Tacstd2 and Cldn4 in LLC is low (naive means ~0.13–0.28 on the TISMO log scale). With n = 2, Wilcoxon/binomial *p* = 0.50 is the exact two-point value, not a powered test.

Paired-mean plots: `results/tismo/figures/paired_lung_Tacstd2.png`, `paired_lung_Cldn4.png`.

### Public GEMM / GEMM-derived ICB (because TISMO lung pairing is LLC-only)

Processed GEO matrices, gene-subset tables in `results/gemm/tables/`. Primary contrast per series = ICB vs its matched control. Platforms differ (normalized counts, FPKM, depositor log2); deltas are **within-series** only.

| Series | Model | Contrast | n | Tacstd2 Δmean | Cldn4 Δmean | TJ Δmean |
|---|---|---|---|---:|---:|---:|
| GSE182228 | LKB1-deficient LUAD s.c. (**KL**) | aPD-1 vs vehicle | 3 / 3 | −0.041 | −0.024 | +0.106 |
| GSE114601 | Kras;Trp53 nodules (**KP**) | aPD-1 vs vehicle | 2 / 2 | −1.155 | −1.323 | −0.381 |
| GSE157880 | HKP1 orthotopic lung (**KP**) | aPD-1 vs IgG, 0 Gy | 3 / 2 | +0.348 | −0.324 | +0.199 |
| GSE169194 | Kras;p53;Msh2 total viable (**KPM**) | A2V+aPD-1 vs A2V (no PD-1 mono) | 3 / 3 | +0.225 | −0.811 | −0.591 |

Sign count across these four pairs: **Tacstd2 2 up / 2 down** (binomial *p* = 1.0); **Cldn4 0 up / 4 down** (*p* = 0.125); **TJ 2 up / 2 down** (*p* = 1.0). Wilcoxon on four deltas is *p* = 1.0 / 0.125 / 0.625. GSE182228 Tacstd2/Cldn4 sit near the FPKM floor (vehicle means log2(FPKM+1) = 0.24 and 0.04). GSE114601 is the only KP nodule series where both genes are clearly expressed, and both means fall after aPD-1 (n = 2).

Secondary GSE182228 contrast (aPD-1+palbociclib vs palbociclib): Tacstd2 Δ = −0.65, Cldn4 Δ = −0.30, TJ Δ = −0.16 (n = 3 / 3).

Figures: `results/gemm/figures/paired_gemm_kp_kl.png`, `waterfall_gemm.png`.

## IFNγ in vitro (optional)

TISMO `downVitroExprn` for IFNg × IFNg-annotated lines: **13 pairs**.

| Gene | Up | Down | Tie | Wilcoxon *p* | Notes |
|---|---:|---:|---:|---:|---|
| Tacstd2 | 4 | 4 | 5 | 0.84 | 5 pairs are 0 vs 0 (CT26, EMT6, KPC, LLC, MC38 in RTM28723893 / XW33589424) |
| Cldn4 | 5 | 6 | 2 | 0.76 | LLC 0 vs 0; KPC (pancreatic, not lung KP) Δmean = +1.07 |

Neither gene shows a consistent IFNγ-up direction across lines. LLC in vitro Tacstd2 and Cldn4 are at the floor.

## Sources

- TISMO: http://tismo.cistrome.org → https://tismo.pku-genomics.org/ (metadata `/tismo/metaData/{vivoMeta,vitroMeta,cellLineMeta}`; expression `/rtismo/gene/downVivoExprn` and `downVitroExprn`). Downloaded 2026-08-19.
- GSE182228 (PMID 36871040), GSE114601 (PMID 30087114), GSE157880 (Ban et al. club-cell RT + PD-1), GSE169194 (PMID 34380768).

## Reproduce

```bash
pip install -r scripts/requirements.txt
python3 scripts/download_tismo.py
python3 scripts/analyze_tismo.py
python3 scripts/analyze_gemm.py
```

## Files

| Path | Content |
|---|---|
| `data/tismo/vivo/*.csv` | Official per-gene ICB tables |
| `data/tismo/vitro/*.csv` | Official IFNg tables |
| `data/tismo/*Meta.json` | TISMO metadata |
| `results/tismo/tables/pairs_merged.tsv` | 64-slice Tacstd2 / Cldn4 / TJ deltas |
| `results/tismo/tables/stats_by_stratum.tsv` | All / lung / KL-KP-LLC stats |
| `results/tismo/figures/` | Waterfalls + paired-mean plots |
| `results/gemm/tables/gemm_contrasts.tsv` | GEO ICB contrasts |
| `results/tismo/summary.json`, `results/gemm/summary.json` | Machine-readable numbers |

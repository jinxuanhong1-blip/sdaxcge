# Finding — leftover mouse ICI sc/spatial Cldn4 (not Tacstd2-only)

Additive only. **TISMO 49/64** and **PR #285** (no n≥2 epithelial ICB vs a true no-ICB control in that scRNA pool) are taken as given and are not re-scored here.

Public processed GEO only: **GSE297632**, **GSE261890**, **GSE285606**. **GSE303943 skipped** (CMT167R subcutaneous PKCi vs solvent — not ICB).

Question: epithelial **Cldn4** versus T-cell / exclusion scores. Honest library n. Extra figures use a **Cldn4-detected + T-low** cut (Cldn4>0 and T ≤ median T among Cldn4+). Median-median cuts are not used as the claim: both genes are sparse, so both medians are 0.

Primary table: `results/mouse_ici_cldn4_leftover/per_series_table.tsv`.

## Per-series table

| Series | What the public files are | Contrast | n | ICB vs true control? | Epi Cldn4 | T / exclusion | p |
|---|---|---|---|---|---|---|---|
| **GSE297632** | LLC subcutaneous scRNA; residual 7d after anti-PD-1 vs untreated | aPD-1 residual vs control | **1 vs 1** library (epi cells 455 vs 133) | yes, but n=1; **subcutaneous** | mean 0.027 → **0.137**; %pos 3.1% → **8.3%** | T/NK fraction 0.226 → **0.101** | NA (n<2) |
| **GSE261890** | BMKMANU S1000 L7 bins; titles Sen.Res.Veh vs Sen.Res.ICI | Veh vs ICI (not a 4-group sen/res split) | **1 vs 1** section (10,793 vs 10,797 bins) | Veh vs ICI; overall-design lists 4 groups, **2 GSM deposited** | bin mean 0.058 → 0.088; %pos 7.5% → 10.9% | T score 0.010 → 0.024; same-bin ρ(Cldn4,T)= **+0.093 / +0.105** | NA (n<2); bin ρ p≪0.05 but **positive** |
| **GSE285606** | Kras/p53 344SQ parental vs PD1R1 (140P) acquired aPD-1-resistant | PD1R1 vs WT; **both IgG** | **2 vs 2** (wk4 + wk6) | **no** | epi mean WT 0.211 vs PD1R1 **0.021**; %pos 13.0% vs **1.3%** | T/NK fraction WT 0.239 vs PD1R1 **0.360** | Welch Cldn4 **0.25**; T/NK frac **0.27** |
| **GSE303943** | CMT167R subcutaneous PKCi vs solvent | — | skip | **no — not ICB** | — | — | — |

Same-library epithelial vs T/NK Cldn4 (6 leftover libraries with ≥20 cells each): **6/6** epi > T/NK; Wilcoxon signed-rank p=**0.031**. That is a compartment restriction, not an ICB contrast.

## Verdict

No leftover series here is a powered epithelial **ICB vs true control** test of Cldn4 (the only n=2 vs 2 is GSE285606, and both arms are IgG).

- **GSE297632** (n=1 vs 1, subcutaneous LLC): residual/tolerant library is Cldn4-higher and T/NK-lower than control. Direction only.
- **GSE261890** (n=1 vs 1): public files are Veh vs ICI, not a separable sensitive-vs-resistant 2×2. Same-bin and neighbor ρ(Cldn4, T) are **weakly positive**. Cldn4+ bins are *less* often T=0 than the section background (Veh 79% vs 89%; ICI 63% vs 76%). That is not a Cldn4-high / T-excluded neighborhood.
- **GSE285606** (n=2 vs 2, both IgG): the acquired-resistant PD1R1 line is Cldn4-**low** and T/NK-**higher** than parental 344SQ. Opposite a Cldn4-high / immune-low pattern. Under-powered (Welch p=0.25 / 0.27).

## GSE297632 — LLC ± anti-PD-1 (tolerant residual)

Subcutaneous LLC. Control vs cells remaining 7 days after anti-PD-1. GEM-X Flex 10x MTX. Epithelial = Epcam+ or (Cdh1+ and Krt8+); T/NK = Cd3d/e, Cd8a, Nkg7, Ncr1 (Epcam wins over ambient T UMIs). min UMI 200.

| Library | n cells | n epi | n T/NK | epi Cldn4 mean | epi Cldn4 %pos | T/NK fraction | epi minus T/NK Cldn4 |
|---|---|---|---|---|---|---|---|
| LLC_Control | 8,441 | 455 | 1,909 | 0.027 | 3.1% | 0.226 | +0.024 |
| LLC_aPD1_residual | 9,423 | 133 | 949 | 0.137 | 8.3% | 0.101 | +0.134 |

Detected cut inside epithelium: Control 7/14 Cldn4+ cells are T-low (50%; background T=0 in epi = 36%). Residual 9/11 Cldn4+ are T-low (82%; background 59%). Cell n for Cldn4+ epithelium is **14 and 11** — do not treat the cut as a sample-level test.

Epithelial cell ρ(Cldn4, Tscore): Control +0.008 (p=0.87); residual −0.14 (p=0.11).

## GSE261890 — NSCLC ICI spatial (sensitive vs resistant as deposited)

BMKMANU S1000. Overall-design text lists immunotherapy-sensitive+veh, sensitive+anti-PD1, resistant+veh, resistant+anti-PD1. **Only two GSM are public:** Sen.Res.Veh and Sen.Res.ICI. They are scored as Veh vs ICI, n=1 vs 1. L7 is the primary bin (not L1 subcellular).

| Section | n L7 bins | Cldn4 mean | Cldn4 %pos | T score | same-bin ρ | neighbor ρ (k=6) | Cldn4+ & T-low |
|---|---|---|---|---|---|---|---|
| Sen.Res.Veh | 10,793 | 0.058 | 7.5% | 0.010 | **+0.093** (p=2.8×10⁻²²) | **+0.152** (p=4.0×10⁻⁵⁷) | 638 / 811 Cldn4+ |
| Sen.Res.ICI | 10,797 | 0.088 | 10.9% | 0.024 | **+0.105** (p=1.0×10⁻²⁷) | **+0.175** (p=7.9×10⁻⁷⁵) | 747 / 1,182 Cldn4+ |

Both Cldn4 and T are higher on the ICI section (descriptive). The extra cut figure highlights Cldn4+ / T=0 bins; that fraction is **below** the section-wide P(T=0). Do not read the red points as exclusion independent of sparsity.

No public spot-level sensitive vs resistant label is in the processed MTX. The 730 MB RDS was not used (no R/Seurat in this run).

## GSE285606 — Kras/p53 ICI-resistant scRNA (344SQ vs PD1R1)

344SQ parental (WT) vs 344SQ_PD1R1 / 140P, derived after in vivo anti-PD-1 resistance. **Both libraries are IgG**, week 4 (T1) and week 6 (T2). This is a genotype contrast on IgG, not ICB vs control.

| Library | n cells | n epi | n T/NK | epi Cldn4 mean | epi Cldn4 %pos | T/NK fraction |
|---|---|---|---|---|---|---|
| WT_IgG_T1 | 7,947 | 2,723 | 1,417 | 0.133 | 9.5% | 0.178 |
| WT_IgG_T2 | 11,380 | 2,422 | 3,410 | 0.289 | 16.4% | 0.300 |
| PD1R1_IgG_T1 | 5,506 | 1,027 | 2,260 | 0.015 | 1.1% | 0.410 |
| PD1R1_IgG_T2 | 9,219 | 2,050 | 2,849 | 0.026 | 1.5% | 0.309 |

PD1R1 − WT: epi Cldn4 Δ = −0.190 (Welch p=0.25, MWU p=0.33); T/NK fraction Δ = +0.121 (Welch p=0.27). Cldn4+ epithelial cells are few in PD1R1 (11 and 31).

## Extra figures (Cldn4-high + immune-low cut)

`results/mouse_ici_cldn4_leftover/figures/`

- `scrna_cldn4_epi_vs_tnk.png` — library epithelial vs T/NK Cldn4
- `scrna_cldn4_vs_tnk_frac.png` — exclusion-style scatter (library n)
- `cut_<sample>_epithelial.png` — epithelial Cldn4+ & T-low
- `spatial_cut_Sen.Res.*.png` / `spatial_map_cut_Sen.Res.*.png` — L7 detected cut and XY map

## Methods (short)

- Cldn4 = ENSMUSG00000047501 (symbol match). Tacstd2 scored but not the claim.
- Unsorted 10x: epithelial/tumor = Epcam+ or (Cdh1+ and Krt8+); T/NK = Cd3d/e, Cd8a, Nkg7, Ncr1. Do **not** use Sftpc/Scgb1a1. Epcam wins over ambient T UMIs.
- Expression = log1p(UMI). T score = mean log1p of Cd3d/e/g, Cd2, Cd8a, Cd8b1, Nkg7, Gzmb, Prf1, Ifng (present genes).
- Welch / MWU only if n≥2 libraries/arm. n=1 vs 1 reports direction only.
- Spatial primary bin = L7_heAuto. Neighbor T = mean of 6 nearest L7 bins.
- Public processed MTX/TSV/tar only. Large MTX/TAR stay gitignored.

## 中文摘要

TISMO 49/64 与 PR #285（该池无 n≥2 上皮 ICB vs 真对照）视为已知，不重算。本切片只补 GSE297632 / GSE261890 / GSE285606。GSE303943 为皮下 PKCi，已跳过。

**GSE297632** 皮下 LLC，n=1 vs 1：残存抗 PD-1 库上皮 Cldn4 更高（0.137 vs 0.027）、T/NK 比例更低（0.101 vs 0.226），只报方向。**GSE261890** 公开只有 Veh vs ICI 两张切片，不是可拆的敏感/耐药 2×2；同 bin 与邻域 ρ(Cldn4,T) 为弱阳性，不是排斥。**GSE285606** 是 IgG 上的亲本 vs 获得性耐药株（n=2 vs 2）：耐药株上皮 Cldn4 更低、T/NK 更高（Welch p=0.25 / 0.27），与 Cldn4 高 / 免疫低相反，也不是 ICB vs 对照。

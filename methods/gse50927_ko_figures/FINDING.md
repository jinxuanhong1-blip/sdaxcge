# FINDING — GSE50927 Cldn4 KO extra figure pack

**Additive only.** Redraws a high-contrast extra pack from **locked** numbers. Does **not** recompute GSEA, does **not** invent new *p*, and does **not** audit or retract any slide.

GSE50927 is **naive mouse whole lung** (mixed-cell bulk, n=1 vs 1). It is **not lung cancer** and it is **not** a tumour KO.

Locked sources:

- NES / FDR: `methods/cldn4_ko_gsea/tables/gsea_headline.tsv` (contrast `GSE50927_naive_lung`)
- 34/39 IFN+MHC panel: `results/fable_cldn4_kdko/panel_ifn_apm_results.json` (EXTENDED)
- Cldn4 QC: author edgeR table `GSE50927_Cldn4lungWTvsKOgenes.csv` (same file as the GSEA extra)

---

## 一句话 / TL;DR

On GSE50927 Cldn4 KO vs WT (whole lung, **not** lung cancer):

| Locked item | Value | Source |
|---|---|---|
| Cldn4 logFC | **−6.061** (author FDR 4.07×10⁻²⁶) | deposited edgeR |
| Hallmark IFN-γ NES | **+1.51** (FDR 0.005) | `cldn4_ko_gsea` |
| Hallmark IFN-α NES | **+1.59** (FDR 0.010) | `cldn4_ko_gsea` |
| MHC-I / APM NES | **+1.42** (FDR 0.050) | `cldn4_ko_gsea` |
| Extended IFN+MHC panel | **34 / 39 genes UP** (median logFC +0.43) | `fable_cldn4_kdko` |
| Core 8 | **6 / 8 UP** | same panel |

No new *p* is reported on the figures. Author edgeR FDR and the already-published GSEA FDR are copied as deposited / locked.

---

## What was redrawn

High-contrast extra pack under `methods/gse50927_ko_figures/figures/`:

1. `fig_extra1_ifn_mhc_heatmap` — IFN/MHC heatmap of the **39 detected** panel genes (**34 UP**). ISG vs MHC-I/APM arms. Color = author edgeR logFC.
2. `fig_extra2_nes_bar` — locked NES bars (IFN-γ **+1.51**, IFN-α +1.59, MHC-I +1.42; TJ / keratinization negative, FDR≥0.05).
3. `fig_extra3_cldn4_qc` — Cldn4 QC vs other expressed claudins / Tacstd2 / Epcam. Only Cldn4 collapses.
4. `fig_extra4_panel_lollipop` — same 39 genes as a ranked lollipop (34 up / 5 down).
5. `fig_extra5_core8` — user core 8 (IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A, TAP1, TAP2) → mouse matches; 6/8 up.
6. `fig_extra6_cldn4_qc_cpm` — second QC: author logCPM vs logFC (Cldn4 is detected, not a low-count drop).

PNG + PDF for each.

---

## Design (locked; not re-run)

| Item | Choice |
|---|---|
| Series | [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927) |
| Tissue | Naive **whole lung**, mixed 129S6/C57BL/6/BALB/c. **Not lung cancer.** |
| n | **1 vs 1** (author edgeR; no count matrix on GEO) |
| Rank / logFC | Author edgeR KO minus WT |
| GSEA | Already run: weighted KS *p*=1, 1000 gene-set permutations, seed=42. **Not recomputed.** |
| Panel | EXTENDED ISG + MHC-I/APM from `scripts/fable_cldn4_kdko/panel_ifn_apm.py` (43 human symbols; 39 detected on the mouse table) |
| New statistics | **None.** Figures plot frozen tables only. |

Mouse matching is the same as the locked 34/39 count (including shared orthologs: HLA-A and HLA-B both → H2-K1; GBP1 and GBP2 both → Gbp2). Four human symbols are not on the table: IFI44L, IFI6, HERC5, HLA-F.

Five detected genes are down: IFIT2, IFIT3, HLA-A (H2-K1), HLA-B (H2-K1), TAP2.

---

## Locked GSEA (copy; not recomputed)

| Set | NES | FDR | nom p | n in rank |
|---|---:|---:|---:|---:|
| Hallmark IFN-γ | +1.508 | 0.005 | 0.001 | 181 |
| Hallmark IFN-α | +1.587 | 0.010 | 0.004 | 89 |
| MHC-I / APM | +1.423 | 0.050 | 0.030 | 21 |
| KEGG tight junction | −0.963 | 0.283 | 0.283 | 154 |
| GO keratinization | −1.060 | 0.202 | 0.162 | 69 |

Positive NES = up after Cldn4 loss.

---

## Cldn4 QC (author edgeR, deposited)

| Gene | logFC | logCPM | author FDR |
|---|---:|---:|---:|
| **Cldn4** | **−6.061** | 2.251 | **4.07×10⁻²⁶** |
| Cldn3 | +0.091 | 7.126 | 1 |
| Cldn7 | +0.113 | 5.421 | 1 |
| Tacstd2 | −0.145 | 5.894 | 1 |
| Epcam | −0.088 | 5.990 | 1 |
| Actb | −0.023 | 10.687 | 1 |

Cldn4 is abundant enough (logCPM 2.25) and uniquely collapsed. Neighbor claudins and Tacstd2 are flat on this unreplicated table.

---

## Previously reported directional *p* (citation only)

`fable_cldn4_kdko` already reported one-sided tests on the same 34/39 panel (MW vs background *p* = 7.2×10⁻¹¹; Wilcoxon *p* = 3.1×10⁻⁷; sign test *p* = 1.2×10⁻⁶). Those values are **not recomputed** and are **not printed on the extra figures**.

---

## What this is not

- Not lung cancer. Not a LUAD/LUSC/NSCLC KO.
- Not SKB264 / TROP2-ADC.
- Not a new GSEA run and not a new DE *p*.
- Not sample-permutation GSEA (n=1 vs 1 cannot support that).
- Not a re-cut of GSE207704 or GSE22493 (those stay in `methods/cldn4_ko_gsea`).

---

## Reproduce

```bash
pip install -r methods/gse50927_ko_figures/requirements.txt
python3 scripts/gse50927_ko_figures/plot_extra.py
```

The script only reads `methods/gse50927_ko_figures/tables/` and writes `figures/fig_extra*.png` (+ PDF). It asserts 34/39 and refuses to finish with fewer than 4 extra PNGs.

---

## 中文摘要

只补 GSE50927 Cldn4 KO 的高对比额外图，不重算、不新造 *p*、不审不撤已有页。

这是小鼠**全肺**（naive，1 vs 1），**不是肺癌**。Cldn4 logFC −6.06（作者 FDR 4.07×10⁻²⁶）。已锁定的 prerank GSEA：IFN-γ NES **+1.51**（FDR 0.005），IFN-α +1.59，MHC-I +1.42。IFN/MHC 扩展面板 **34/39** 上调。图：热图、NES 条形、Cldn4 QC、lollipop、core 8、logCPM QC。

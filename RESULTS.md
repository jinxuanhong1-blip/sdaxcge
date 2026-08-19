# CLDN4-only leftover Visium lung adenocarcinoma / NSCLC (GEO)

Additive CLDN4-only screen. No private 8-gene keratin-like (8-KL) signature. No claim-failed language.

## Search

NCBI GEO (`gds`, entry type GSE), 2026-08-19:

```
visium AND lung AND (adenocarcinoma OR NSCLC OR LUAD)
```

25 series. Broader `"spatial transcriptomics" AND lung AND (adenocarcinoma OR NSCLC OR LUAD)` (61 series) was used only to catch Visium series that omit the word “Visium” in the series record.

**Already assigned (excluded from run):** GSE307534, GSE277206, GSE221733 (GeoMx), GSE299786, GSE300007, GSE299886, GSE311609.

## Metrics (CLDN4 only)

For every leftover slide with a processed matrix and CLDN4+CD8A+KRT8:

1. **Spot Spearman** — Spearman ρ(CLDN4, CD8A) on in-tissue spots with ≥100 UMI, after log1p(10⁴ × CPM).
2. **Q4 vs Q1 nearest CD8 distance** — CD8⁺ = raw CD8A > 0. Q4/Q1 = top/bottom CLDN4 quartile. Median Euclidean pixel distance to nearest CD8⁺ spot; Mann–Whitney Q4 vs Q1. Positive `q4_minus_q1` means CLDN4-high spots are farther from CD8⁺ spots.
3. **KRT8 residual** — OLS residual of log-normalized CLDN4 ~ KRT8, then Spearman ρ(residual, CD8A).

## Series opened

### Already assigned (opened, not re-run)

| GSE | Why skipped |
| --- | --- |
| GSE307534 | Already assigned. Public Visium LUAD (normal / AAH / AIS / invasive). |
| GSE277206 | Already assigned. Visium early LUAD in never-smokers. |
| GSE221733 | Already assigned. GeoMx DSP NSCLC, not Visium. |
| GSE299786 | Already assigned. |
| GSE300007 | Already assigned. |
| GSE299886 | Already assigned. |
| GSE311609 | Already assigned. Xenium + Visium comparison (lung and breast). |

### Leftover — ran

Numbers filled after `scripts/analyze_cldn4.py`. Per-slide table: `results/cldn4_visium_leftover.csv`.

| GSE | Decision |
| --- | --- |
| GSE189487 | **RAN.** Open 10x Visium LUAD (AIS/MIA/IAC), 6 slides, processed MTX + coordinates. Not in the literal `visium` keyword query (titled `[ST]`) but confirmed Visium-format files. |
| GSE273378 | **RAN.** Open Visium FFPE stage I LUAD, 16 slides, processed MTX + coordinates. |
| GSE322553 | **RAN (NSCLC only).** Open Visium; 2/7 slides are NSCLC (`GSM9554209`, `GSM9554210`). CRC/HCC slides not analyzed. |

### Leftover — opened and skipped

| GSE | Why skipped |
| --- | --- |
| GSE322553 CRC/HCC | Same series as above; not lung. |
| GSE303162 | Mouse CAF / Treg spatial, not human. |
| GSE292299 | Open Visium NSCLC (16) + DLBCL (2), processed H5 + spatial tars (~3.1 GB). Leftover and eligible; not downloaded in this pass (size). Prefer GSE189487 + GSE273378. |
| GSE325196 | Brain-metastasis Visium/Xenium/snRNA (mixed primaries). Not primary lung adeno; lung-origin subset not isolated here. |
| GSE308260 | Mouse thoracic adhesions, not tumor. |
| GSE325706 | Visium HD ASPS / GPNMB CAR, not lung. |
| GSE282057 | Same GPNMB CAR study, not lung adeno. |
| GSE301973 | Open Visium HD EGFR-mutant NSCLC, processed tar (3.1 GB, bin-level). Leftover and eligible; skipped this pass (HD size). |
| GSE288758 | Open Visium HD EGFR NSCLC pre/post osimertinib (2.3 GB). Same skip as HD size. |
| GSE300676 | Open standard Visium LUAD micropapillary, 8 slides, processed H5 (1.5 GB). Leftover and eligible; not downloaded in this pass (size; two other LUAD series already running). |
| GSE288479 | Mostly scRNA/WES/WTS LUAD. One Visium sample (`GSM8768597`) exists; series is not Visium-first. Not run. |
| GSE270431 | Spatial samples are mouse / xenograft LUAD models, not open human tumor Visium. |
| GSE268426 | snPATHO-seq / snRNA-seq, not Visium. |
| GSE246011 | Visium of bladder and gastric adenocarcinoma, not lung. |
| GSE222901 | Mouse scRNA-seq tobacco LUAD, not Visium. |
| GSE205354 | PDAC, not lung. |
| GSE215361 / GSE202159 | Mouse Treg Visium, not human. |
| GSE215858 | Mouse AT1 → LUAD Visium. |
| GSE193460 | Mouse Perturb-map LUAD. |
| GSE179572 | Visium brain metastases; primary site not restricted to LUAD/NSCLC in the series record. |
| GSE300827 | Open Visium-style H5 of pre-invasive LUAD (20 slides) but **no tissue coordinates** on GEO, so nearest-CD8 distance cannot be run. Not in the literal `visium` keyword query. |
| GSE172416 | Normal human lung Visium (CellDART), not adenocarcinoma. |

## Run numbers

_Pending download + analysis._

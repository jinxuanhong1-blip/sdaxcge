# CLDN4-only leftover Visium lung adenocarcinoma / NSCLC (GEO)

Additive CLDN4-only screen. No private 8-gene keratin-like (8-KL) signature.

Search date: 2026-08-19. NCBI GEO `gds`, entry type GSE.

Primary query:

```
visium AND lung AND (adenocarcinoma OR NSCLC OR LUAD)
```

25 series. A broader `"spatial transcriptomics" AND lung AND (adenocarcinoma OR NSCLC OR LUAD)` query (61 series) was used only to catch Visium series that omit the word “Visium” in the series record (GSE189487, GSE300827).

**Already assigned (not re-run):** GSE307534, GSE277206, GSE221733 (GeoMx), GSE299786, GSE300007, GSE299886, GSE311609.

## Metrics (CLDN4 only)

In-tissue spots with ≥100 UMI. Expression = log1p(10⁴ × CPM). CLDN4, CD8A, and KRT8 were present on every run slide.

1. **Spot Spearman** — ρ(CLDN4, CD8A).
2. **Q4 vs Q1 nearest CD8 distance** — CD8⁺ = raw CD8A count > 0. Q4/Q1 = top/bottom CLDN4 quartile. Median Euclidean full-resolution pixel distance to the nearest CD8⁺ spot; two-sided Mann–Whitney Q4 vs Q1. Positive `Δ = med(Q4) − med(Q1)` means CLDN4-high spots are farther from CD8⁺ spots. ~200 px is one Visium center-to-center spacing on these slides.
3. **KRT8 residual** — OLS residual of log-normalized CLDN4 ~ KRT8, then ρ(residual, CD8A).

Scripts: `scripts/download_geo.py`, `scripts/analyze_cldn4.py`. Per-slide CSV: `results/cldn4_visium_leftover.csv`.

## Leftover series that were run (n = 3)

CLDN4+CD8A+KRT8 present on all 24 slides. 20/24 slides have negative spot Spearman. Median ρ = −0.029. KRT8 residual median ρ = −0.007. Median Q4−Q1 distance = −0.15 px (CLDN4-high not farther from CD8 at the median; several GSE273378 slides are ~1 spot closer).

### GSE189487 — RAN (6/6 Visium LUAD)

Open 10x Visium of treatment-naïve LUAD (AIS / MIA / IAC). Titled `[ST]` so it missed the literal `visium` keyword query; files are Visium MTX + `tissue_positions_list`. Processed matrix present. Not controlled.

| sample | n_spots | CD8⁺ | Spearman ρ (p) | KRT8-resid ρ (p) | Q1 med px | Q4 med px | Δ Q4−Q1 | MW p |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | --- |
| TD1 | 4094 | 541 | −0.071 (5.6e-6) | −0.028 (0.075) | 201.2 | 200.7 | −0.50 | 0.20 |
| TD2 | 4414 | 327 | −0.050 (8.0e-4) | −0.064 (1.8e-5) | 346.8 | 348.0 | +1.19 | 0.0042 |
| TD3 | 1139 | 109 | −0.074 (0.012) | −0.023 (0.43) | 346.8 | 346.8 | 0.00 | 0.94 |
| TD5 | 1700 | 314 | −0.161 (2.8e-11) | −0.077 (0.0016) | 200.2 | 201.1 | +0.87 | 2.3e-6 |
| TD6 | 3911 | 196 | −0.052 (0.0010) | −0.015 (0.35) | 401.4 | 400.0 | −1.38 | 0.063 |
| TD8 | 1760 | 257 | −0.030 (0.22) | −0.039 (0.10) | 201.2 | 201.1 | −0.13 | 0.074 |

Series: median Spearman ρ = **−0.062** (6/6 negative). Median KRT8-resid ρ = **−0.033** (6/6 negative). Median Δ = **−0.06 px** (2/6 Q4 farther). Stouffer p for signed Δ = 0.28.

### GSE273378 — RAN (16/16 Visium FFPE stage I LUAD)

Open Visium FFPE, 16 tumors from the VI-signature validation cohort (stage IA/IB except one upstaged IIA). Processed MTX + coordinates. Not controlled.

| sample | n_spots | CD8⁺ | Spearman ρ (p) | KRT8-resid ρ (p) | Q1 med px | Q4 med px | Δ Q4−Q1 | MW p |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | --- |
| LM_SD_1216_1 | 2401 | 322 | −0.009 (0.66) | +0.034 (0.097) | 397.0 | 199.7 | −197.3 | 9.4e-21 |
| LM_SD_16 | 2755 | 54 | −0.030 (0.12) | −0.018 (0.34) | 690.0 | 718.3 | +28.3 | 0.0051 |
| LM_SD_11 | 2724 | 143 | −0.031 (0.10) | +0.058 (0.0025) | 397.0 | 398.5 | +1.52 | 1.7e-4 |
| LM_SD_2 | 1728 | 232 | −0.021 (0.37) | −0.003 (0.90) | 343.9 | 199.7 | −144.2 | 0.11 |
| LM_SD_3 | 3902 | 853 | −0.010 (0.53) | +0.020 (0.22) | 199.0 | 199.0 | −0.05 | 1.9e-4 |
| LM_SD_4 | 3524 | 1409 | −0.106 (2.4e-10) | −0.043 (0.011) | 198.8 | 198.8 | 0.00 | 0.47 |
| LM_SD_5 | 3230 | 154 | −0.038 (0.029) | −0.030 (0.089) | 525.4 | 525.7 | +0.32 | 0.30 |
| LM_SD_6 | 3194 | 335 | +0.064 (2.8e-4) | +0.049 (0.0058) | 344.7 | 344.7 | 0.00 | 0.065 |
| LM_SD_7 | 2822 | 758 | +0.028 (0.14) | +0.087 (3.3e-6) | 271.5 | 198.8 | −72.7 | 7.6e-31 |
| LM_SD_1216_8 | 2836 | 377 | −0.033 (0.076) | +0.009 (0.64) | 344.4 | 199.7 | −144.7 | 1.2e-7 |
| LM_SD_9 | 2801 | 290 | −0.006 (0.75) | +0.017 (0.38) | 597.4 | 344.9 | −252.5 | 3.1e-23 |
| LM_SD_10 | 2620 | 379 | +0.002 (0.94) | +0.044 (0.024) | 199.7 | 199.0 | −0.69 | 0.025 |
| LM_SD_1216_12 | 3604 | 94 | −0.019 (0.26) | −0.011 (0.50) | 528.1 | 526.9 | −1.22 | 0.021 |
| LM_SD_13 | 2912 | 577 | +0.019 (0.31) | +0.176 (9.9e-22) | 199.0 | 198.8 | −0.17 | 3.8e-11 |
| LM_SD_1216_14 | 1915 | 357 | −0.120 (1.4e-7) | −0.085 (2.0e-4) | 199.0 | 199.0 | +0.05 | 0.0077 |
| LM_SD_15 | 1744 | 134 | −0.016 (0.51) | −0.003 (0.89) | 398.2 | 346.0 | −52.1 | 0.029 |

Series: median Spearman ρ = **−0.017** (12/16 negative). Median KRT8-resid ρ = **+0.013** (7/16 negative). Median Δ = **−0.43 px**; mean Δ = −52 px because several slides have Q4 ~1 Visium spacing closer to CD8⁺. Stouffer p for signed Δ = 9.3e-29. LM_SD_6 has heavy CLDN4 ties (quartile collapse). LM_SD_16 has only 54 CD8⁺ spots.

### GSE322553 — RAN NSCLC only (2/7 slides)

Open Visium of CRC (n=4) + HCC (n=1) + NSCLC (n=2). CRC/HCC not analyzed. Processed H5 + coordinates. Not controlled.

| sample | n_spots | CD8⁺ | Spearman ρ (p) | KRT8-resid ρ (p) | Q1 med px | Q4 med px | Δ Q4−Q1 | MW p |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | --- |
| NSCLC_936_1 | 2894 | 775 | −0.029 (0.11) | +0.046 (0.014) | 199.0 | 199.0 | 0.00 | 2.0e-4 |
| NSCLC_41485_2 | 3945 | 39 | −0.024 (0.13) | −0.012 (0.46) | 527.0 | 456.4 | −70.6 | 3.5e-8 |

Series: both NSCLC slides negative Spearman (median ρ = **−0.027**). NSCLC_41485_2 has only 39 CD8⁺ spots; distance is unstable.

## Every GSE opened

### Already assigned

| GSE | Why skipped |
| --- | --- |
| GSE307534 | Already assigned. Public Visium LUAD (normal / AAH / AIS / invasive). |
| GSE277206 | Already assigned. Visium early LUAD in never-smokers. |
| GSE221733 | Already assigned. GeoMx DSP NSCLC, not Visium. |
| GSE299786 | Already assigned. |
| GSE300007 | Already assigned. |
| GSE299886 | Already assigned. |
| GSE311609 | Already assigned. Xenium + Visium comparison (lung and breast). |

### Leftover — opened and skipped

| GSE | Why skipped |
| --- | --- |
| GSE322553 CRC/HCC | Same series as the NSCLC run; those 5 slides are not lung. |
| GSE303162 | Mouse CAF / Treg spatial, not human. |
| GSE292299 | Open Visium NSCLC (16) + DLBCL (2), processed H5 + spatial (~3.1 GB). Eligible leftover; not downloaded this pass after 24 slides from three other series were already run. |
| GSE325196 | Brain-metastasis Visium/Xenium/snRNA, mixed primaries. Not primary lung adeno. |
| GSE308260 | Mouse thoracic adhesions, not tumor. |
| GSE325706 | Visium HD ASPS / GPNMB CAR, not lung. |
| GSE282057 | Same GPNMB CAR study, not lung adeno. |
| GSE301973 | Open Visium HD EGFR-mutant NSCLC, 3.1 GB bin-level tars. Eligible leftover; skipped this pass (HD size). |
| GSE288758 | Open Visium HD EGFR NSCLC pre/post osimertinib, 2.3 GB. Same HD-size skip. |
| GSE300676 | Open standard Visium LUAD micropapillary, 8 slides, 1.5 GB H5. Eligible leftover; not downloaded this pass (two other LUAD series already run). |
| GSE288479 | Mostly scRNA/WES/WTS LUAD. One Visium sample (`GSM8768597`); series is not Visium-first. |
| GSE270431 | Spatial samples are mouse / xenograft LUAD models, not open human tumor Visium. |
| GSE268426 | snPATHO-seq / snRNA-seq, not Visium. |
| GSE246011 | Visium of bladder and gastric adenocarcinoma, not lung. |
| GSE222901 | Mouse scRNA-seq tobacco LUAD, not Visium. |
| GSE205354 | PDAC, not lung. |
| GSE215361 / GSE202159 | Mouse Treg Visium, not human. |
| GSE215858 | Mouse AT1 → LUAD Visium. |
| GSE193460 | Mouse Perturb-map LUAD. |
| GSE179572 | Visium brain metastases; series record does not restrict primaries to LUAD/NSCLC. |
| GSE300827 | Open Visium-style H5 of pre-invasive LUAD (20 slides) but **no tissue coordinates** on GEO, so nearest-CD8 distance cannot be run. |
| GSE172416 | Normal human lung Visium (CellDART), not adenocarcinoma. |

No leftover series in the primary query was raw-FASTQ-only. None of the leftover human Visium series above were controlled-access.

## Combined leftover run

| series | slides | median Spearman ρ | n ρ<0 | median KRT8-resid ρ | median Δ Q4−Q1 (px) | n Q4 farther |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GSE189487 | 6 | −0.062 | 6/6 | −0.033 | −0.06 | 2/6 |
| GSE273378 | 16 | −0.017 | 12/16 | +0.013 | −0.43 | 4/16 |
| GSE322553 NSCLC | 2 | −0.027 | 2/2 | +0.017 | −35.3 | 0/2 |
| **all run** | **24** | **−0.029** | **20/24** | **−0.007** | **−0.15** | **6/24** |

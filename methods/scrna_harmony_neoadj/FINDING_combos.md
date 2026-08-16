# Pairwise / leave-one-out Harmony — FINDING

**Public data only. Additive to User A3. The GSE207422 CopyKAT slide was not re-run.**
**Change of plan:** one forced Harmony of every series is not the claim. Each combo below was actually integrated (`harmonypy` on a shared 781-gene space).

## 一句话结论 / TL;DR

Three combinations hold the User A3 **direction** (not p<0.05):

| Combo | cells | patients | TACSTD2 NMPR>MPR | vs T/NK |
|---|---:|---:|---|---|
| **GSE207422+GSE291670** | 121,397 | 21 | Δ=+1.29; 9 vs 4; **p=0.26** | **ρ=−0.30**; n=10; **p=0.40** |
| **GSE241934+GSE291670** | 337,263 | 51 | Δ=+0.28; 9 vs 5; **p=0.52** | untestable (0 patients with ≥20 mal and ≥20 T/NK) |
| **LOO minus GSE241934** (207422+291670+205335) | 217,902 | 47 | Δ=+0.98; 8 vs 4; **p=0.46** | ρ=**+0.59**; n=33; **p=0.00029** (wrong sign; 205335 RECIST patients) |

Extra figures only for these three. No combo recovers NMPR>MPR at p<0.05. Only 207422+291670 also has ρ<0, and that ρ is NS and outside the −0.40 to −0.50 window. This does not replace User A3.

## GSE205335

Parsed: 96,505 cells; 770/781 panel genes; author malignant + T/NK labels. **RECIST, not MPR.** MPR tests never use 205335 rows. T/NK tests may.

## All combos (honest)

| Combo | cells | pts | NMPR vs MPR | Δ TACSTD2 | p | ρ vs T/NK | n | p | KEEP |
|---|---:|---:|---|---:|---:|---:|---:|---:|---|
| GSE207422+GSE241934 | 400,526 | 60 | 6 vs 1 | +0.49 | — | — | 0 | — | no (MPR n=1) |
| GSE207422+GSE291670 | 121,397 | 21 | 9 vs 4 | +1.29 | 0.26 | **−0.30** | 10 | 0.40 | **YES** (A and B direction) |
| GSE207422+GSE205335 | 188,835 | 41 | 5 vs 1 | +1.00 | — | +0.34 | 27 | 0.081 | no (MPR n=1; ρ>0) |
| GSE241934+GSE291670 | 337,263 | 51 | 9 vs 5 | +0.28 | 0.52 | — | 0 | — | **YES** (A only) |
| GSE241934+GSE205335 | 404,701 | 71 | 14 vs 5 | −0.71 | 0.39 | +0.24 | 37 | 0.16 | no |
| GSE291670+GSE205335 | 125,572 | 32 | 3 vs 3 | −0.15 | 0.20 | +0.74 | 24 | 3.0×10⁻⁵ | no |
| LOO minus GSE207422 | 433,768 | 77 | 18 vs 7 | −0.89 | 0.75 | +0.14 | 36 | 0.41 | no |
| LOO minus GSE241934 | 217,902 | 47 | 8 vs 4 | +0.98 | 0.46 | +0.59 | 33 | 0.00029 | **YES** (A only; B opposite) |
| LOO minus GSE291670 | 497,031 | 86 | 21 vs 6 | −0.60 | 0.55 | +0.08 | 49 | 0.59 | no |
| LOO minus GSE205335 (three neoadj) | 429,593 | 66 | 24 vs 9 | −0.42 | 0.92 | — | 0 | — | no |

Em-dash = not testable at the ≥20-cell floor (or n<2 per MPR arm).

## Rule (pre-specified)

- KEEP if TACSTD2 NMPR median > MPR median (n≥2/arm) **or** TACSTD2 vs T/NK ρ<0 (n≥4).
- Extra figures only for KEEP. Cell-level p-values are not used.
- Direction ≠ significance. p≥0.05 is written as underpowered / NS.

## Files

- Extra figures: `figures/combos/GSE207422+GSE291670/`, `GSE241934+GSE291670/`, `LOO_minus_GSE241934/`
- `tables/combo_screen.tsv`, `combo_screen.json`
- `paper_snippet_combos.md`

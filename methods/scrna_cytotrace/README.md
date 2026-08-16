# Additive scRNA: CytoTRACE-like stemness / cycling vs TACSTD2/CLDN4

Public **GSE207422** and **GSE241934** only. Malignant (or author epithelial) cells.

- Stemness = residual gene-count (CytoTRACE-**like**; CytoTRACE2 R package not run)
- Cycling = Tirosh S + G2M
- Correlate with TACSTD2 / CLDN4; split by MPR; combinatorial cycling vs non-cycling
- Extra: is TACSTD2-high more keratin / differentiated?

Primary unit is the **sample**. Cell-level ρ is exploratory.

**Headline (honest):** TACSTD2-high is more keratin in both datasets (GSE207422 n=9 p=0.027; GSE241934 n=35 p=1.2×10⁻¹⁰). CytoTRACE-like vs TACSTD2 is a non-significant trend in GSE207422 (ρ=−0.63, p=0.067, n=9) and **null** in GSE241934 (ρ=−0.015, p=0.93, n=35). Stemness is not MPR-associated. See `FINDING.md`.

See `METHODS.md`. Numbers: `results/stats.tsv`.

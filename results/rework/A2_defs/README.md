# Rework A2: TACSTD2 vs immune definitions in public durvalumab NSCLC RNA

**Self-contained recompute. Every requested definition is reported. Nothing is hidden.**

## Claim and prior

| Source | n | Immune def | Raw ρ | Raw p | Purity-adj ρ | Adj p |
|---|---|---|---|---|---|---|
| User | ? | unspecified | **-0.65** | 2e-4 (which test unspecified) | **-0.46** | paired with 2e-4 in the slide |
| Prior A2 (this repo) | 29 post | ESTIMATE ImmuneScore | -0.665 | 8.3e-5 | -0.474 | **0.011** |
| Prior A2 (this repo) | 32 pre | ESTIMATE ImmuneScore | -0.716 | 4.0e-6 | -0.530 | 0.0022 |

Public durvalumab RNA used here (and in the prior A2):

- [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) pre-treatment, n=32
- [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378) post-treatment, n=29
- Same randomized phase II neoadjuvant durvalumab ± SBRT trial (Altorki / NCT02904954).

## Results — all five definitions

Partial Spearman always controls **ESTIMATEScore ranks** (n−3 df).

| Dataset | n | Feature | Raw ρ | Raw p | Adj ρ | Adj p |
|---|---|---|---|---|---|---|
| GSE253564_pretreatment | 32 | **ESTIMATE_ImmuneScore** | -0.7163 | 4.03e-06 | -0.5295 | 0.00219 |
| GSE253564_pretreatment | 32 | **CYT** | -0.4611 | 0.0079 | -0.1652 | 0.375 |
| GSE253564_pretreatment | 32 | **CD8A** | -0.5363 | 0.00156 | -0.2425 | 0.189 |
| GSE253564_pretreatment | 32 | **GEP18** | -0.4699 | 0.00665 | -0.0827 | 0.658 |
| GSE253564_pretreatment | 32 | **xCell_CD8** | -0.6525 | 5.19e-05 | -0.5707 | 8.02e-04 |
| GSE248378_posttreatment | 29 | **ESTIMATE_ImmuneScore** | -0.665 | 8.30e-05 | -0.4739 | 0.0109 |
| GSE248378_posttreatment | 29 | **CYT** | -0.8059 | 1.32e-07 | -0.7103 | 2.29e-05 |
| GSE248378_posttreatment | 29 | **CD8A** | -0.7059 | 1.89e-05 | -0.5773 | 0.0013 |
| GSE248378_posttreatment | 29 | **GEP18** | -0.6296 | 2.53e-04 | -0.404 | 0.033 |
| GSE248378_posttreatment | 29 | **xCell_CD8** | -0.8177 | 6.10e-08 | -0.7283 | 1.12e-05 |

### Supplement (not the five named scores; still reported)

| Dataset | Feature | Raw ρ | Raw p | Adj ρ | Adj p |
|---|---|---|---|---|---|
| GSE253564_pretreatment | xCell_CD8_raw_ssgsea | -0.8024 | 3.33e-08 | -0.6929 | 1.56e-05 |
| GSE253564_pretreatment | xCell_CD8_Tcm | -0.5117 | 0.00276 | -0.2228 | 0.228 |
| GSE253564_pretreatment | xCell_CD8_Tem | -0.4956 | 0.00392 | -0.2907 | 0.113 |
| GSE253564_pretreatment | xCell_CD8_naive_T-cells | -0.1404 | 0.443 | -0.2395 | 0.194 |
| GSE248378_posttreatment | xCell_CD8_raw_ssgsea | -0.8236 | 4.06e-08 | -0.757 | 3.13e-06 |
| GSE248378_posttreatment | xCell_CD8_Tcm | -0.6901 | 3.43e-05 | -0.5126 | 0.00528 |
| GSE248378_posttreatment | xCell_CD8_Tem | -0.6478 | 1.45e-04 | -0.4889 | 0.00829 |
| GSE248378_posttreatment | xCell_CD8_naive_T-cells | -0.5374 | 0.00264 | -0.4913 | 0.00794 |

## Verdict (honest — all five definitions)

**Closest match to the user numbers is post-treatment ESTIMATE ImmuneScore (coefficients) plus a raw p-value (not the adjusted one).**

1. **ESTIMATE ImmuneScore** — reproduces prior A2 exactly. Post n=29: raw ρ=-0.665 (p=8.30e-05), adj ρ=-0.4739 (p=0.0109 = prior 0.011). Pre n=32: raw ρ=-0.7163, adj ρ=-0.5295 (p=0.00219). User raw −0.65 ≈ post raw −0.665; user adj −0.46 ≈ post adj −0.474. User p=2e-4 is the **raw** post p-order, not the adjusted p.
2. **CYT** — post remains strong after purity (raw -0.8059, adj -0.7103, p=2.29e-05). Pre collapses (raw -0.4611, adj -0.1652, p=0.375).
3. **CD8A** — same pattern. Post adj ρ=-0.5773 (p=0.0013); pre adj ρ=-0.2425 (p=0.189, NS).
4. **GEP18** (all 18 genes present in both matrices) — post raw ρ=-0.6296, raw p=2.53e-04 is the single closest match to the user's **p=2e-4**, but adj ρ=-0.404 (p=0.033) is weaker than −0.46. Pre adj ρ=-0.0827 (p=0.658) is null — almost all of the pre GEP18 signal is purity.
5. **xCell CD8** (spillover-adjusted CD8+ T-cells) — negative in both cohorts after purity: pre adj ρ=-0.5707 (p=8.02e-04); post adj ρ=-0.7283 (p=1.12e-05). Stronger than the user −0.46; does not reproduce that coefficient. Raw ssGSEA (no spillover) is even stronger (see supplement).

- **p=2e-4 attached to adj ρ=−0.46 does not reproduce** at n=29 or n=32 (that pairing needs n≈55+). The prior A2 adj p=0.011 is the correct p for the ESTIMATE post adjusted coefficient.
- **Do not cite a single ρ.** Pre-treatment CYT / CD8A / GEP18 lose significance after ESTIMATEScore adjustment; ESTIMATE ImmuneScore and xCell CD8 do not.

## Methods

- **Expression:** GEO supplementary FPKM; duplicate symbols collapsed by max; TACSTD2 / CD8A / CYT on log2(FPKM+1).
- **ESTIMATE:** Python port of v1.0.13 `estimateScore` (official `SI_geneset.gmt`).
- **CYT:** Rooney *Cell* 2015, (log2(GZMA+1)+log2(PRF1+1))/2.
- **GEP18:** Ayers *JCI* 2017 18 genes listed in `scripts/rework_A2_defs.py`. Score = mean of per-gene z-scores of log2(FPKM+1). Unweighted; NanoString TIS regression weights are platform-specific and were not applied.
- **xCell CD8:** Aran *Genome Biol* 2017. ssGSEA (tau=0.25, no GSVA normalization) on the 15 official CD8+ T-cells signatures, mean-aggregated, RNA-seq power/calibration transform, spillover α=0.5 via non-negative least squares. Gene universe = all genes in each FPKM matrix (pre 20,187; post 15,165), not the locked ~10.8k xCell training list. Raw aggregated ssGSEA (rank-identical to the transform, before spillover) is in the supplement.
- **Partial Spearman:** Pearson on ranks of (TACSTD2, feature, ESTIMATEScore); p from Student t with n−3 df. Same estimator as prior A2.
- **No multiple-testing correction** is applied across the five definitions — this is a pre-specified sensitivity panel, and every p-value is shown raw.

## Caveats

1. ESTIMATE TumorPurity cosine is out of bounds for many FPKM samples (pre 23/32, post 6/29); only score *ranks* are used.
2. Both cohorts are small and from one early-stage neoadjuvant trial. The user slide did not name the accession or the immune definition.
3. xCell on n=29–32 is below the size where spillover is stable; that is why raw ssGSEA CD8 is also reported.
4. GEP18 here is an unweighted z-mean, not the commercial TIS.
5. Bessede et al. 2024 (TROP2 vs ICI resistance) used atezolizumab OAK/POPLAR under EGA access — not these public durvalumab series.

## Reproduce

```bash
pip install pandas numpy scipy matplotlib
python3 scripts/rework_A2_defs.py
```

## Files

- `results.json` — machine-readable metrics + methods notes
- `correlations.tsv` — one row per dataset × feature
- `scores_GSE253564_pretreatment.csv`, `scores_GSE248378_posttreatment.csv`
- `figures/scatter_all_defs.png`, `figures/bar_partial_rho.png`
- `scripts/rework_A2_defs.py`


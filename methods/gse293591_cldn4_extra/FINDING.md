# FINDING — GSE293591 leftover extra: CLDN4 vs CD274 and IFN beyond MHC-I

**Additive only. CLDN4 only.** Public BostonGene solid-tumor RNA-seq (Kushnarev et al., *Sci Rep* 2025, PMID [40753302](https://pubmed.ncbi.nlm.nih.gov/40753302/); GEO [GSE293591](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293591)). The paper scores nine IHC biomarkers (ESR1, PGR, AR, MKI67, ERBB2, **CD274**, CDX2, KRT7, KRT20). It does **not** score CLDN4. MHC-I / APM on this TPM matrix is treated as **already known (partial)** and is a reference row only. This leftover is **CD274** and **IFN after MHC-I genes are removed**. Extra cuts are median / tertile / Q4 vs Q1.

Unit is the **GEO sample / TPM column**. No FASTQ. No IHC PD-L1 on GEO. No ICI arm.

---

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| TPM columns / unique GSM | yes | **365** | `GSE293591_TPM_all_samples.tsv.gz`; 20,062 genes |
| Unique `Sample_*` IDs | yes | **365** | titles run `Sample_001`…`Sample_377` with gaps; **do not write n=377** |
| FFPE / total RNA / target enrichment | yes | **313** | Procrustes batch-corrected (author) |
| FF / polyA | yes | **52** | Cureline fresh-frozen |
| Breast cancer | yes | **201** | GEO `diagnosis` |
| Lung Adenocarcinoma | yes | **20** | GEO `diagnosis` |
| Squamous Cell Carcinoma of Lung | yes | **14** | GEO `diagnosis` |
| **Lung (LUAD+LUSC)** | yes | **34** | does **not** include non-lung SCC (n=13) |
| GI (CRC + gastric + eso + panc + liver/CCA) | yes | **71** | locked bag; not a single histology |
| CLDN4 / CD274 finite | yes | **365 / 365** | both present on every column |
| PD-L1 IHC (CPS / TC) | **no** | **0** | paper has IHC; not deposited |
| ICI / response | **no** | **0** | biomarker series, not a trial |
| FASTQ / SRA | **no** | **0** | author withheld identifiable sequence |

Primary leftover n is **365** (all deposited tumors). FFPE-only is **313**. Breast is **201**. Lung is **34** (exploratory; n<40). Do not write n=365 for a lung-only test.

---

## One-row table (all tumors, n=365)

log2(TPM+1). Signature = mean of gene-wise z. Spearman vs continuous CLDN4. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6). BH *q* is inside the **six leftover** endpoints on this cohort (CD274, Hallmark IFN-γ no MHC-I, same minus CD274, Ayers 6, ISG no MHC-I, IFNG). MHC-I rows are **not** in that BH.

| pair | leftover? | n | ρ | 95% CI | p | q (BH leftover) | ρ_adj epi | p_adj | extra-cut Q4 vs Q1 | verdict |
|---|---|---:|---:|---|---:|---:|---:|---:|---|---|
| CLDN4 vs **CD274** | yes | **365** | **−0.112** | −0.216 to −0.005 | 0.032 | 0.096 | −0.032 | 0.55 | 92 vs 92, δ=−0.16, MWU 0.067 | **NO_EVIDENCE** after BH / residual |
| CLDN4 vs Hallmark IFN-γ **minus MHC-I** | yes | 365 | +0.052 | −0.064 to +0.159 | 0.33 | 0.37 | +0.048 | 0.36 | δ=+0.11, p=0.21 | **NO_EVIDENCE** |
| CLDN4 vs IFN-γ minus MHC-I **and CD274** | yes | 365 | +0.053 | −0.062 to +0.162 | 0.31 | 0.37 | +0.048 | 0.36 | δ=+0.11, p=0.20 | **NO_EVIDENCE** |
| CLDN4 vs Ayers 6-gene IFN-γ | yes | 365 | −0.047 | −0.159 to +0.055 | 0.37 | 0.37 | −0.016 | 0.77 | δ=−0.08, p=0.38 | **NO_EVIDENCE** |
| CLDN4 vs compact ISG (no MHC-I) | yes | 365 | +0.048 | −0.068 to +0.150 | 0.36 | 0.37 | +0.036 | 0.49 | δ=+0.08, p=0.35 | **NO_EVIDENCE** |
| CLDN4 vs **IFNG** | yes | 365 | **−0.123** | −0.228 to −0.017 | 0.019 | 0.096 | −0.066 | 0.21 | δ=−0.16, MWU 0.062 | **NO_EVIDENCE** after BH / residual |
| CLDN4 vs MHC-I / APM (21) | **already known** | 365 | **+0.238** | +0.130 to +0.342 | 4.1×10⁻⁶ | — | +0.178 | 6.4×10⁻⁴ | δ=+0.33, p=1.4×10⁻⁴ | **POSITIVE** (reference) |
| CLDN4 vs MHC-I core (HLA-A/B/C+B2M+TAP1/2+TAPBP) | already known | 365 | **+0.289** | +0.183 to +0.388 | 2.0×10⁻⁸ | — | +0.249 | 1.4×10⁻⁶ | δ=+0.39, p=6.0×10⁻⁶ | **POSITIVE** (reference) |

**Leftover call:** MHC-I stays positive. CD274 and IFN **beyond MHC-I** do not. The two crude leftover *p*<0.05 (CD274, IFNG) fail leftover BH (*q*=0.096) and die after epithelial residual. Hallmark / Ayers / ISG signatures with MHC-I genes removed are null on n=365.

Full numbers: `tables/spearman_cldn4_vs_endpoints.tsv`, `tables/extra_cuts_cldn4.tsv`.

---

## Extra cuts (all tumors)

CLDN4 log2(TPM+1) cuts on n=365: median 6.624 (183 vs 182); tertile ≤6.117 vs ≥7.204 (122 vs 122); Q1 ≤5.808 vs Q4 ≥7.549 (92 vs 92). Cliff δ > 0 means CLDN4-high is endpoint-high.

| endpoint | median δ (p) | tertile δ (p) | Q4 vs Q1 δ (p) |
|---|---|---|---|
| **CD274** | −0.12 (0.054) | **−0.15 (0.046)** | −0.16 (0.067) |
| Hallmark IFN-γ no MHC-I | +0.04 (0.49) | +0.06 (0.39) | +0.11 (0.21) |
| Ayers 6 | −0.04 (0.50) | −0.08 (0.27) | −0.08 (0.38) |
| ISG no MHC-I | +0.05 (0.44) | +0.04 (0.58) | +0.08 (0.35) |
| **IFNG** | **−0.12 (0.043)** | **−0.17 (0.020)** | −0.16 (0.062) |
| MHC-I / APM (known) | **+0.25 (4.1×10⁻⁵)** | **+0.29 (9.1×10⁻⁵)** | **+0.33 (1.4×10⁻⁴)** |

CD274 / IFNG extra cuts are the same weak negative as the Spearman. They do not become a Q4 law. MHC-I extra cuts stay positive on every cut.

---

## Subsets (honest n; leftover vs known)

| cohort | n | CLDN4–CD274 ρ (p) | CLDN4–IFN-γ no MHC-I ρ (p) | CLDN4–MHC-I/APM ρ (p) |
|---|---:|---|---|---|
| all | **365** | −0.112 (0.032); q=0.096; adj p=0.55 | +0.052 (0.33) | **+0.238 (4.1×10⁻⁶)** |
| FFPE only | **313** | −0.015 (0.80) | +0.124 (0.029); q=0.086 | **+0.213 (1.5×10⁻⁴)** |
| breast | **201** | −0.142 (0.045); q=0.13; adj p=0.97 | −0.033 (0.64) | +0.143 (0.043) |
| GI | **71** | −0.168 (0.16) | −0.117 (0.33) | +0.256 (0.031) |
| lung LUAD+LUSC | **34** | −0.354 (0.040) **EXPLORATORY** | −0.245 (0.16) | +0.096 (0.59) |

FFPE-only Hallmark IFN-γ no MHC-I is a weak positive (ρ=+0.12, *q*=0.086). It is **not** the leftover headline and is not called HOLDS. Lung n=34 is below the n≥40 leftover rule; the CD274 ρ=−0.35 is written as exploratory, not a lung law. Non-lung SCC (n=13) is not in the lung bag.

---

## What “IFN beyond MHC-I” is

Hallmark IFN-γ is 200 genes (198/200 present; `MARCHF1`, `WARS1` absent). MHC-I / APM is the locked 21-gene panel from `methods/gse285029_cldn4_gsea` (HLA-A/B/C/E/F/G, B2M, TAP1/2, TAPBP/TAPBPL, NLRC5, PSMB8/9/10, ERAP1/2, CALR, CANX, PDIA3, IRF1). **Beyond MHC-I** = that Hallmark list with those 21 genes removed (187/189). A second leftover score also drops CD274 so the IFN test is not the PD-L1 test. Ayers 6 = IDO1, CXCL9, CXCL10, HLA-DRA, STAT1, IFNG (6/6; HLA-DRA is MHC-II, kept). Compact ISG = STAT1, IRF1, CXCL9, CXCL10, CXCL11, IDO1, IFNG, GBP1, ISG15, MX1, IFIT1, OAS2 (12/12; no classical MHC-I).

Single leftover genes on n=365 are mixed and are **not** an IFN program: ISG15 ρ=+0.305 (*q*=3.5×10⁻⁸) is the one strong leftover positive; IRF1 / MX1 are weak positives; STAT1 / CD274 / IFNG are weak negatives. The signature mean eats that mix and is null.

---

## Already known (partial): MHC-I

CLDN4-high tracks MHC-I / APM on this matrix (ρ=+0.24, residual +0.18). Core classical genes are the same direction (HLA-A +0.22, HLA-B +0.21, HLA-C +0.23, TAP1 +0.31, TAPBP +0.31). That is the partial already on the table. It is **not** evidence that CLDN4-high is IFN-high or PD-L1-high. Full Hallmark IFN-γ (which still contains MHC-I + CD274) is null on n=365 (ρ=+0.070, p=0.18).

Companion only: CLDN4 vs TACSTD2 ρ=+0.308; vs epithelial mean-z ρ=+0.562. CLDN4 sits on the epithelial / TACSTD2 side.

---

## What this does not test

- PD-L1 IHC CPS / TC (not on GEO; paper IHC is not re-scored here).
- ICI response, PFS, or OS.
- A lung-only law (n=34).
- Cell-intrinsic IFN. This is mixed-histology bulk FFPE/FF TPM.
- A claim that CLDN4 loss opens IFN / PD-L1. This is a cross-sectional correlation, not a KD/KO.

---

## Reproduce

```bash
python3 -m pip install -r methods/gse293591_cldn4_extra/requirements.txt
python3 methods/gse293591_cldn4_extra/analyze.py
```

Downloads (not committed) go to `$GSE293591_CLDN4_DATA` (default `/tmp/gse293591_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293591/suppl/GSE293591_TPM_all_samples.tsv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293591/matrix/GSE293591_series_matrix.txt.gz`

---

## Files

- `analyze.py` — download, scores, Spearman / partial Spearman, extra cuts, figures
- `tables/label_inventory.tsv` — honest n
- `tables/diagnosis_counts.tsv`
- `tables/gene_coverage.tsv`
- `tables/spearman_cldn4_vs_endpoints.tsv` — leftover + MHC-I reference, all cohorts
- `tables/spearman_cldn4_vs_genes.tsv` — CD274 / IFN genes / MHC-I core
- `tables/extra_cuts_cldn4.tsv` — median / tertile / Q4 vs Q1
- `tables/one_row.tsv` — all-tumor leftover + MHC-I
- `tables/per_sample.tsv`
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd274.png`
- `figures/fig2_spearman_forest.png`
- `figures/fig3_q4q1_violins.png`
- `figures/fig4_subset_forest.png`

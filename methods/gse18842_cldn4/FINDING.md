# GSE18842 NSCLC array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** No dual-high. Public Sanchez-Palencia resected NSCLC Affymetrix U133 Plus 2.0 series (Sanchez-Palencia et al., *Int J Cancer* 2011, PMID 20878980; GEO [GSE18842](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18842)). Unit is the **array**. No slide was re-scored. No ICI arm.

Primary pairs are **CLDN4 vs CD8A**, **CLDN4 vs CD274**, and **CLDN4 vs ESTIMATE ImmuneScore** on **tumor** arrays. Adjacent-normal / control arrays are inventoried and are **not** mixed into those tests.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **91** | 54,675 probes × 91 GSM; RMA; 0 missing values |
| Unique GSM / unique titles | yes | **91** | `Tumor N` / `Control N` / `NP`; all unique |
| Tumor arrays (`sample type: tumor`) | yes | **46** | **PRIMARY n**; `source_name` = Human Lung Tumor |
| Control / adjacent-normal arrays | yes | **45** | `sample type: control`; **dropped** |
| Title-paired tumor+control | yes | **44** | numeric pair IDs shared by Tumor and Control titles |
| Unpaired tumors (`Tumor NP_1`, `Tumor NP_2`) | yes | **2** | GEO `overall_design`: “All samples are paired except three” |
| Unpaired controls (`Control NP`) | yes | **1** | the third unpaired array |
| LUAD / LUSC / histology | **no** | 0 | not a GEO characteristic; **no split invented** |
| Stage | no | 0 | paper discusses stage; not deposited per array |
| OS / DFS | no | 0 | not deposited |
| ICI / treatment | no | 0 | resected atlas, not an ICI series |
| Tumor % / ABSOLUTE purity | no | 0 | public proxies: epithelial mean-z and ESTIMATE |
| CLDN4 finite on tumors (`201428_at`) | yes | **46** | named Plus2 probe; official GPL570 = CLDN4 |
| CD8A finite on tumors (`205758_at`) | yes | **46** | named Plus2 probe; official GPL570 = CD8A |
| CD274 finite on tumors (`227458_at`) | yes | **46** | named Plus2 probe (higher-mean of two clean CD274 probes) |
| ESTIMATE ImmuneScore on tumors | computed | **46** | Yoshihara Immune141 ssGSEA (137/141) |
| **Primary pairwise n** | yes | **46** | complete-case CLDN4 + CD8A + CD274 + ImmuneScore |

Do not write n=91 for the correlations. Mixing tumor with control would inflate n and is not an NSCLC-tumor test. The computable public tumor n is **46 arrays**.

Paper / secondary text often cites 14 adenocarcinoma + 32 squamous among the 46 tumors. Those counts are **not** per-array GEO characteristics. This slice does not assign LUAD/LUSC from titles, lab IDs, or literature. No histology split.

Deposited matrix: RMA. Spearman is rank-based, so a further monotone transform would not change ρ.

## One-row table

| dataset | histology | platform | n tumor | n control dropped | CLDN4 | CD8A | CD274 | CLDN4–CD8A ρ (p) | adj ρ (p) | CLDN4–CD274 ρ (p) | adj ρ (p) | CLDN4–ImmuneScore ρ (p) | adj ρ (p) | verdict | ICI |
|---|---|---|---:|---:|---|---|---|---|---|---|---|---|---|---|---|
| GSE18842 Sanchez-Palencia | NSCLC mixed; LUAD/LUSC **not deposited** | GPL570 U133 Plus 2.0 RMA | **46** | 45 | `201428_at` | `205758_at` | `227458_at` | **+0.260 (0.081)** | +0.217 (0.15) | **−0.000 (0.998)** | +0.055 (0.72) | **+0.190 (0.21)** | +0.279 (0.064) | **NO_EVIDENCE** | not deposited |

Full numbers: `tables/one_row.tsv`, `tables/spearman_cldn4_vs_cd8a_cd274.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary)

Named Plus2 probes. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6** present). HOLDS rule (same as PR 229 / GSE4573): n≥40, ρ_adj<0, p_adj<0.05.

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **46** | **+0.260** | −0.005 to +0.501 | 0.081 | +0.217 | 0.15 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | **46** | **−0.000** | −0.316 to +0.310 | 0.998 | +0.055 | 0.72 | **NO_EVIDENCE** |
| CLDN4 vs **ImmuneScore** | **46** | **+0.190** | −0.106 to +0.462 | 0.21 | +0.279 | 0.064 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 12 vs 12 | — | — | MWU 0.069 | — | — | null (rank-biserial **+0.44**) |
| CLDN4 Q4 vs Q1 on CD274 | 12 vs 12 | — | — | MWU 0.93 | — | — | null (rank-biserial −0.03) |
| CLDN4 Q4 vs Q1 on ImmuneScore | 12 vs 12 | — | — | MWU 0.67 | — | — | null (rank-biserial +0.11) |

CLDN4-high is not CD8A-low, not CD274-low, and not ImmuneScore-low on this NSCLC array. The CD8A crude interval is a weak positive that includes zero; after the epithelial residual it stays positive and non-significant. Q4 vs Q1 on CD8A leans the same way (CLDN4-high has *higher* CD8A, rank-biserial +0.44) and is still not a claim.

Alternate CD274 probe `223834_at` (also GPL570 = CD274; not the headline): ρ = −0.060, p = 0.69; ρ_adj = +0.016, p = 0.92. Same verdict.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z | 46 | **+0.562** | 0.268 to 0.787 | 4.8×10⁻⁵ |

CLDN4 sits on the epithelial side of this array. That is compatible with a tight-junction / tumour-cell program. It is not evidence that CLDN4-high NSCLC is immune-low or PD-L1-low.

## Positive controls (not the claim)

| pair | n | ρ | p | ρ_adj | p_adj |
|---|---:|---:|---:|---:|---:|
| CD8A vs ImmuneScore | 46 | +0.724 | 1.3×10⁻⁸ | +0.745 | 4.4×10⁻⁹ |
| CD274 vs CD8A | 46 | +0.458 | 0.0014 | +0.476 | 9.4×10⁻⁴ |

CD8A tracks ImmuneScore. The ImmuneScore axis is not inert on this matrix. CLDN4 still does not.

ESTIMATE TumorPurity (Yoshihara cosine) falls outside [0, 1] for **35 / 46** tumors. Those values are not used as a 0–1 purity fraction. The independent residual reported above is **epithelial mean-z**, not the cosine.

## What this does not test

- A LUAD vs LUSC split (histology is not deposited per array; not invented from the paper’s 14 / 32 sentence).
- ICI response or PFS (not on GEO).
- OS as a CLDN4 endpoint (not deposited).
- Pathologist CD8 / PD-L1 IHC.
- n=91 (tumor+control mix).
- Dual-high TACSTD2×CLDN4 (this slice is CLDN4-only).

## Reproduce

```bash
python3 -m pip install -r methods/gse18842_cldn4/requirements.txt
python3 methods/gse18842_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE18842_CLDN4_DATA` (default `/tmp/gse18842_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18842/matrix/GSE18842_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz`

## Files

- `analyze.py` — GEO download, tumor-only filter, ESTIMATE ImmuneScore, Spearman / epithelial residual, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv` — headline row
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/probe_confirm.tsv` — named CLDN4 / CD8A / CD274 probes
- `tables/spearman_cldn4_vs_cd8a_cd274.tsv` — primary pairs, positive controls, probe sensitivity
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 (12 vs 12)
- `tables/gene_coverage.tsv`
- `tables/sample_annotation.tsv` — 91 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274_immunescore.png`
- `figures/fig2_cldn4_forest.png`

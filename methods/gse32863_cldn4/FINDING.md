# GSE32863 LUAD array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** Public Selamat / Laird-Offringa resected LUAD Illumina HumanWG-6 v3.0 series (Selamat et al., *Genome Res* 2012, PMID 22613842; GEO [GSE32863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32863); GPL6884). Unit is the **tumor array**. Matched adjacent non-tumor arrays are dropped. No slide was re-scored. No ICI arm. No dual-high. TACSTD2 is not a claim.

Primary question: **CLDN4 vs CD8A**, **CLDN4 vs CD274**, and **CLDN4 vs ESTIMATE ImmuneScore**, with an **epithelial residual**.

## Honest n

Do **not** write n=60. The GEO series *summary* says “60 lung adenocarcinoma tumors”; the *overall design* and the deposited series matrix are **58 LUAD + 58 adjacent non-tumor**. Do **not** write n=59 (that is the methylation companion [GSE32867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32867) / paper Infinium 27k set). Do **not** write n=116 (that pool includes matched normal).

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **116** | 48,803 probes × 116 GSM; log2 + lumi RSN as deposited |
| Unique GSM / unique titles | yes | **116** | all unique |
| Adjacent non-tumor lung | yes | **58** | `source_name_ch1` / `tissue: Normal lung`; **dropped** |
| **Tumor LUAD arrays (primary n)** | yes | **58** | `source_name_ch1` = Lung adenocarcinoma; `tissue: Lung tumor` |
| Unique tumor patients (title id) | yes | **58** | one tumor array per id (`05L*` / `30xx` / `L*`) |
| Matched pairs (shared title id) | yes | **57** | unpaired tumor `3023`; unpaired normal `3035`. Pairing is inventory only |
| Stage | yes | 58 | deposited; **not** the claim |
| Smoking / pack-years | yes | 58 | deposited; **not** the claim |
| KRAS / EGFR / LKB1 | yes (tumor slot) | 58 | deposited; **not** the claim |
| Recurrence yes/no | yes | 58 | deposited; **not** a survival claim (no time) |
| OS / DSS time or event | no | 0 | not a GEO characteristic |
| ICI / treatment | no | 0 | 2012 methylation/expression atlas, not an ICI series |
| Numeric tumor % / ABSOLUTE purity | no | 0 | protocol is macrodissection; values not deposited |
| CLDN4 finite (`ILMN_2132458`) | yes | **58** | single GPL6884 probe; max-mean collapse |
| CD8A finite (max-mean of 3 probes) | yes | **58** | collapse chose `ILMN_2353732` |
| CD274 finite (`ILMN_1701914`) | yes | **58** | single GPL6884 probe |
| Epithelial mean-z | yes | **58** | 6/6 of EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |
| ESTIMATE ImmuneScore | computed | **58** | Yoshihara Immune141 ssGSEA (141/141 genes) |
| **Primary pairwise n** | yes | **58** | complete-case tumor CLDN4 + CD8A + CD274 + ImmuneScore |

The computable public tumor n is **58 LUAD arrays**. That is the n in the one-row table. Title ids pair **57** of the 58 tumors to an adjacent array; tumor `3023` has no same-id normal and normal `3035` has no same-id tumor. Tests do not use pairing.

Tumor-only clinical inventory (not tested): smoking Never 29, Current 29; stage Stage IB 18, Stage IA 16, Stage IIIA 12, Stage IIA 9, Stage IIB 2, Stage IV 1; gender Female 45, Male 13; ethnicity White 36, Asian 22; KRAS WT 36, MUT 22; EGFR WT 41, MUT 17; tissue source nci-edrn 44, ttr-data 14; expression batch Batch 2 35, Batch 1 23.

## One-row table

| dataset | histology | platform | n tumor | CLDN4–CD8A ρ (p) | ρ\|epi (p) | CLDN4–CD274 ρ (p) | ρ\|epi (p) | CLDN4–ImmuneScore ρ (p) | ρ\|epi (p) | verdict |
|---|---|---|---:|---|---|---|---|---|---|---|
| GSE32863 Selamat | LUAD (normals dropped) | GPL6884 HumanWG-6 v3.0 log2-RSN | **58** | **+0.084 (0.531)** | +0.096 (0.477) | **+0.130 (0.332)** | +0.164 (0.223) | **+0.047 (0.728)** | +0.047 (0.727) | **NO_EVIDENCE / NO_EVIDENCE / NO_EVIDENCE** |

Full numbers: `tables/one_row.tsv`, `tables/spearman.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary)

Deposited log2 + Robust Spline Normalization (lumi). Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM/KRT8/KRT18/KRT19/CDH1/KRT7; **6/6** present). HOLDS rule (same as PR 229 / GSE4573 / GSE10245): n≥40, ρ_adj<0, p_adj<0.05. No dual-high cut.

| pair | subset | n | ρ | 95% CI | p | ρ_adj epi (p) | verdict |
|---|---|---:|---:|---|---:|---|---|
| CLDN4 vs **CD8A** | tumor LUAD | **58** | **+0.084** | -0.169 to +0.335 | 0.531 | +0.096 (0.477) | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | tumor LUAD | **58** | **+0.130** | -0.139 to +0.383 | 0.332 | +0.164 (0.223) | **NO_EVIDENCE** |
| CLDN4 vs **ImmuneScore** | tumor LUAD | **58** | **+0.047** | -0.200 to +0.294 | 0.728 | +0.047 (0.727) | **NO_EVIDENCE** |
| CLDN4 vs epithelial mean-z | tumor LUAD | 58 | -0.110 | -0.337 to +0.120 | 0.412 | — | context |
| CD8A vs ImmuneScore | tumor LUAD | 58 | +0.862 | +0.774 to +0.912 | 3.83e-18 | — | positive control |

ImmuneScore is Yoshihara 2013 Immune141 ssGSEA (Barbie/GSVA, τ=0.25; 141/141 genes present after first-symbol max-mean collapse). Stromal141 coverage is 137/141 (not used in the residual). There is no deposited numeric purity; the residual is the **epithelial mean-z**, not ESTIMATE TumorPurity. On this LUAD matrix CLDN4 does **not** track the epithelial mean-z (see context row), so the residual barely moves ρ. CD8A vs ImmuneScore is the positive control.

Named-probe sensitivity (CLDN4 `ILMN_2132458` vs CD8A `ILMN_1768482` / CD274 `ILMN_1701914`) is in `tables/spearman.tsv`.

### Q4 vs Q1 (descriptive)

See `tables/highlow_cldn4.tsv`. Quartiles are descriptive only and are not a dual-high claim.

## Matrix and labels (nothing invented)

| field | public? | n | what is there |
|---|---|---:|---|
| tissue / source | yes | 116 | 58 tumor / 58 adjacent; tests use tumor only |
| histology | yes | 58 | LUAD only on the tumor arrays |
| ICI response | **no** | 0 | surgical methylation/expression series |
| OS time | **no** | 0 | recurrence yes/no is deposited; unused |
| ESTIMATE published scores | **no** | 0 | ImmuneScore computed here from Yoshihara lists |

Max-mean unique-symbol GPL6884 collapse for the named genes:

| gene | probe used |
|---|---|
| CLDN4 | `ILMN_2132458` |
| CD8A | `ILMN_2353732` |
| CD274 | `ILMN_1701914` |

## What is not done

- No dual-high (no TACSTD2∩CLDN4 high intersection).
- No pooling of matched normal into the immune Spearman.
- No re-use of the series-summary “60 tumors” or the methylation-paper n=59.
- No ICI ORR / PFS model (labels are not deposited).
- No OS model (no time-to-event).

## Files

- `analyze.py` — GEO download, tumor filter, epithelial residual, ImmuneScore, Spearman, figures
- `estimate_yoshihara_2013.tsv` — Yoshihara 2013 Supp Data 1 (141+141)
- `tables/one_row.tsv`, `spearman.tsv`, `label_inventory.tsv`, `per_sample.tsv`, `probe_confirm.tsv`, `gene_coverage.tsv`, `highlow_cldn4.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_cd274.png`
- `figures/fig3_cldn4_vs_immunescore.png`
- `figures/fig4_spearman_forest.png`

```bash
python3 methods/gse32863_cldn4/analyze.py
```

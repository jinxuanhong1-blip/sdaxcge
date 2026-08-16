# REWORK OncoSG A1 — TACSTD2 vs CD8 / GEP / immune after purity

**Self-contained. Public cBioPortal `luad_oncosg_2020` only. Written to be read without the rest of the repo.**

**Verdict: YES — TACSTD2 is negatively associated with CD8, GEP18, and the A1 8-gene immune score after published purity (n=169).**

This is **East-Asian surgical LUAD**, not an ICI-response cohort. A negative
bulk correlation is not evidence that TROP2-high tumors fail checkpoint blockade.

## Why this rework exists

User Claim A1 includes OncoSG (with TCGA). A prior OncoSG-only slice (PR 82)
reported TACSTD2 vs **IMSIG neutrophils** partial Spearman **ρ = −0.456, n = 169**.
That is a myeloid IMSIG score, not CD8, not Ayers GEP, and not the A1 8-gene
immune signature.

This folder recomputes the three endpoints the A1 rework language actually
uses — **CD8, GEP18, immune** — on the same public matrix, after the same
published `PURITY` field. The neutrophil number is recomputed as an
**anchor / sanity check**, not as the primary claim. Filters were not tuned
to recover −0.456.

## Analysis set

| Item | n | Note |
|---|---:|---|
| cBioPortal RNA sample list `luad_oncosg_2020_rna_seq_v2_mrna` | 181 | Portal description says “181 samples” |
| Public z-score matrix columns | **169** | This is the open expression table |
| RNA-list IDs absent from the public matrix | 12 | A008, A114, A122, A136, A139, A147, A184, A302, A435, A484, A489, A507 |
| Clinical samples | 305 | `data_clinical_sample.txt` |
| Clinical PURITY non-NA | 302 | includes many RNA-absent tumors |
| Clinical IMSIG T-cells non-NA | 172 | 169 matrix + 3 IMSIG-only (A184, A484, A489; those three lack PURITY) |
| **Complete-case analysis (TACSTD2 + feature + PURITY)** | **169** | All 169 matrix samples have PURITY and all IMSIG immune scores |

The 12 portal-listed RNA IDs are **not in the public z-score file**. Nine of
them have PURITY but no IMSIG; three have IMSIG but no PURITY. They cannot
enter a TACSTD2 correlation. **n = 169 is the honest public n**, not 181.

## Pre-specified features

| Name | Definition | Honest limitation |
|---|---|---|
| **CD8** | `CD8A` all-sample z-score | Single gene. Not a deconvolution fraction. Sensitivity: mean(`CD8A`,`CD8B`). |
| **GEP18** | unweighted mean of **17/18** Ayers 2017 T-cell-inflamed GEP gene z-scores | **`CCL5` is absent** from the public matrix and from the cBioPortal REST API (0 values). Merck NanoString TIS weights are not public. This is not the clinical assay. |
| **Immune** | original A1 8-gene T-cell effector: `CD8A GZMA GZMB IFNG EOMES CXCL9 CXCL10 TBX21` (all 8 present) | A gene-mean, not ESTIMATE ImmuneScore. |

**Covariate:** published sample-level `PURITY` (OncoSG clinical attribute).
This is **not** TCGA ABSOLUTE and **not** ESTIMATE cosine purity.

**Statistic:** first-order partial Spearman = Pearson correlation of
average-ranks after algebraic adjustment for ranked purity; two-sided t,
df = n−3; Fisher-z 95% CI with SE = 1/√(n−4). Unadjusted Spearman is
reported beside it. BH-FDR is within the 3 primary tests only.

### What cannot be computed (honest skip)

cBioPortal exposes **only z-score** mRNA profiles for this study
(`luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores` and the
diploid-referenced twin). Datahub has no `data_mrna_seq_v2_rsem.txt`.

| Score | Why skipped |
|---|---|
| ESTIMATE ImmuneScore / TumorPurity | Official ssGSEA ranks genes **within a sample** on raw/log expression. Ranking z-scores (standardized **across** samples) is not the same transform. The Affymetrix cosine purity formula is also invalid on z-scores. |
| xCell CD8 / immune | Needs counts or TPM plus RNA-seq spillover calibration. n=169 z-scores are not that input. |
| MCP-counter | Published as mean log2(TPM+1) of marker genes. Z-scores change within-sample ranks. |

Do not pretend those scores exist here.

## Direct answer

After published purity, TACSTD2 is **negatively** associated with all three
primary readouts in this public OncoSG LUAD matrix:

| Feature | Unadj ρ | Unadj p | Partial ρ | 95% CI | Partial p | BH q (3 tests) | n |
|---|---:|---:|---:|---|---:|---:|---:|
| CD8 (`CD8A`) | -0.380 | 3.38e-07 | -0.309 | -0.440 to -0.165 | 4.69e-05 | 4.69e-05 | 169 |
| GEP18 (17/18) | -0.414 | 2.27e-08 | -0.349 | -0.475 to -0.208 | 3.62e-06 | 1.09e-05 | 169 |
| Immune (A1 8-gene) | -0.387 | 1.98e-07 | -0.318 | -0.448 to -0.175 | 2.70e-05 | 4.05e-05 | 169 |

- CD8: ρ = **-0.309** (95% CI -0.440 to -0.165; p = 4.69e-05; n = 169)
- GEP18: ρ = **-0.349** (95% CI -0.475 to -0.208; p = 3.62e-06; n = 169)
- Immune: ρ = **-0.318** (95% CI -0.448 to -0.175; p = 2.70e-05; n = 169)

Purity adjustment **does not create** the sign: unadjusted Spearman is already
negative. Partialling purity shrinks |ρ| (TACSTD2 vs PURITY Spearman
ρ = 0.238, p = 0.0018) but leaves all three primary
tests negative.

**Do not quote the neutrophil −0.456 as if it were CD8 or GEP.** Those are
different endpoints. CD8/GEP/immune |ρ| here is smaller than the IMSIG
neutrophil |ρ|.

## Secondary / original A1 signatures

| Feature | Unadj ρ | Partial ρ | Partial p | n |
|---|---:|---:|---:|---:|
| CYT (GZMA+PRF1) | -0.421 | -0.359 | 1.72e-06 | 169 |
| CD8A+CD8B | -0.372 | -0.299 | 8.44e-05 | 169 |
| A1 cytotoxic 8-gene | -0.425 | -0.363 | 1.32e-06 | 169 |
| A1 exhaustion 12-gene | -0.444 | -0.389 | 1.92e-07 | 169 |

These **exactly match** the original OncoSG A1 per-cohort rows in PR 89
(immune / cytotoxic / exhaustion partial ρ = −0.3178 / −0.3631 / −0.3887,
n=169). CD8 and GEP18 were not primary endpoints in that run.

## Anchor: prior IMSIG neutrophils ρ = −0.456

Recomputed on the **same** Datahub revision and checksums as PR 82:

| | Prior (PR 82) | This run |
|---|---:|---:|
| n | 169 | 169 |
| Neutrophils raw ρ | −0.503 | -0.503 |
| Neutrophils partial ρ | **−0.456** | **-0.456** |
| Partial p | 5.06e-10 | 5.06e-10 |
| BH q (8 IMSIG) | 4.05e-09 | 4.05e-09 |

`matches_3dp` = true.
The neutrophil result **reproduces**. It is **not** a CD8 or GEP result.

Published IMSIG (not recomputed from genes):

| Signature | Unadj ρ | Unadj p | Partial ρ | 95% CI | Partial p | BH q |
|---|---:|---:|---:|---|---:|---:|
| IMSIG B cells | -0.322 | 2.00e-05 | -0.240 | -0.378 to -0.092 | 0.0017 | 1.96e-03 |
| IMSIG interferon | -0.130 | 0.0917 | -0.063 | -0.213 to 0.089 | 0.4161 | 4.16e-01 |
| IMSIG macrophages | -0.462 | 2.49e-10 | -0.408 | -0.527 to -0.274 | 3.98e-08 | 1.06e-07 |
| IMSIG monocytes | -0.440 | 2.13e-09 | -0.381 | -0.504 to -0.244 | 3.40e-07 | 5.43e-07 |
| IMSIG neutrophils | -0.503 | 3.18e-12 | -0.456 | -0.568 to -0.328 | 5.06e-10 | 4.05e-09 |
| IMSIG NK cells | -0.471 | 1.03e-10 | -0.421 | -0.538 to -0.288 | 1.29e-08 | 5.16e-08 |
| IMSIG plasma cells | -0.329 | 1.22e-05 | -0.285 | -0.419 to -0.140 | 0.0002 | 2.36e-04 |
| IMSIG T cells | -0.455 | 4.94e-10 | -0.404 | -0.524 to -0.269 | 5.55e-08 | 1.11e-07 |

IMSIG T cells (closest published “T-cell immune” column): ρ = **-0.404** (95% CI -0.524 to -0.269; p = 5.55e-08; n = 169).
IMSIG interferon is the weak / null member of that set, as in PR 82.

## Data

- Study: [luad_oncosg_2020](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020) — Chen et al., *Nat Genet* 2020
- Expression: Datahub `165bd77077b03038f9c2ee104959eb474770b2a9` `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` (sha256 `44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f`, 23,493,696 bytes)
- Clinical: same revision `data_clinical_sample.txt` (sha256 `55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368`)
- Profile id: `luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores`
- TACSTD2 Entrez 4070; API vs file check: max|API-file|=0.00e+00 on 8 samples

Raw downloads stay in the cache directory and are not committed.

## Caveats

1. **Not ICI.** OncoSG is a genomic LUAD series. No PD-1/PD-L1 treatment labels are used here. Do not call these tumors hot/cold on IO.
2. **Bulk composition.** TACSTD2 is epithelial. A negative correlation with CD8/GEP partly reflects tumor-cell content. That is why the purity-adjusted number is the one to quote — and why it is smaller than the raw ρ.
3. **Purity is a published scalar**, method not re-derived (not ABSOLUTE, not ESTIMATE). Residualizing on it is not cell-intrinsic causality.
4. **GEP18 is 17/18.** CCL5 is missing. We did not impute it and did not swap in TLR2.
5. **No ESTIMATE / xCell / MCP** on this public deposit. If a slide quotes those names for OncoSG, the public files cannot support that sentence.
6. **Z-scores.** Signature means are means of already z-scored genes. Spearman is rank-based, so a monotone per-gene transform would not change ρ; a within-sample ssGSEA (ESTIMATE) would.
7. **n = 169 ≠ 181.** The portal RNA list is larger than the open matrix. We did not drop samples to hit a target n.
8. **No multiple-testing story beyond the stated BH sets.** Secondary rows are descriptive.
9. **East-Asian LUAD.** Do not pool this ρ with TCGA-LUSC and call it “NSCLC.”

## Reproduce

```bash
pip install -r scripts/rework/OncoSG_A1/requirements.txt
python3 scripts/rework/OncoSG_A1/analyze.py \
  --cache-dir /tmp/oncosg_a1 \
  --out-dir results/rework/OncoSG_A1
```

## Artifacts

| File | What |
|---|---|
| `correlations.tsv` | Every feature, raw + partial ρ / p / CI / FDR |
| `sample_table.tsv` | 169 rows: TACSTD2, PURITY, scores, light clinical |
| `gene_coverage.tsv` | Present/missing genes per signature |
| `results.json` | Machine-readable primary numbers + methods flags |
| `provenance.json` | URLs, sha256, bytes |
| `figures/forest_partial_rho.png` | Primary + IMSIG T / neutrophils |
| `figures/scatter_partial_residuals.png` | Rank residuals for CD8 / GEP18 / immune |

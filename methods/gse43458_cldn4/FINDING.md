# GSE43458 never-smoker LUAD array — CLDN4 vs CD8A / CD274 / ImmuneScore

**Additive CLDN4-only.** This slice is **CLDN4** on the public Kadara / Kabbout never-smoker lung adenocarcinoma Affymetrix Human Gene 1.0 ST series (Kabbout et al., *Clin Cancer Res* 2013, PMID 23659968; GEO [GSE43458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43458)). Unit is the **array**. Tumor only. No dual-high class. No slide was re-scored. No ICI arm.

Primary set is GEO `tissue: Lung tumor` + `histology: Lung adenocarcinoma` + `smoking status: Never-smoker`. The 30 paired normals and the 40 smoker tumors are **not** that n.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| Arrays in the series matrix | yes | **110** | 33,252 transcript clusters × 110 GSM; RMA log2 as deposited (BRB-Array Tools) |
| Unique GSM | yes | **110** | all unique |
| **never-smoker LUAD tumor (PRIMARY)** | yes | **40** | GEO Never-smoker + Lung tumor; **this is the n used below** |
| smoker LUAD tumor (extra) | yes | 40 | **not** folded into never-smoker n |
| all LUAD tumors mixed | yes | 80 | 40 never + 40 smoker; composition, not the claim |
| paired normal lung | yes | 30 | paired to 30 never-smoker cases; **excluded** (tumor-only) |
| ICI / treatment | no | 0 | resected atlas, not an ICI series |
| Tumor % / ABSOLUTE purity | no | 0 | only public proxy is an RNA epithelial score |
| OS / stage on GEO characteristics | no | 0 | not deposited as sample characteristics |
| CLDN4 finite (`8133360`) | yes | **40** | single HuGene 1.0 ST transcript cluster |
| CD8A finite (`8053584`) | yes | **40** | single cluster |
| CD274 finite (`8154233`) | yes | **40** | single cluster |
| **Primary pairwise n (never-smoker LUAD tumor)** | yes | **40** | this is the n used below |
| ESTIMATE ImmuneSignature genes present | yes | **132 / 141** | missing HLA-E, GBP1, HLA-B, CD74, KLRK1, HLA-DRA, HLA-F, HLA-G, CD302 |
| dual-high TACSTD2×CLDN4 class | no | 0 | not scored |

Do not write n=110 for the tests. Do not write n=80 by adding smoker tumors. Do not write n=70 by adding the 30 paired normals.

## One-row table

| dataset | histology | platform | n arrays | CLDN4 | CD8A | CD274 | ImmuneScore | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–ImmuneScore ρ (p) | adj ImmuneScore | verdict | ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE43458 Kadara/Kabbout | never-smoker LUAD tumor | GPL6244 HuGene 1.0 ST RMA log2 | **40** | `8133360` | `8053584` | `8154233` | ESTIMATE ImmuneSignature mean-z, 132/141 | **+0.034 (0.83)** | **+0.162 (0.32)** | **−0.032 (0.85)** | −0.004 (0.98) | **NO_EVIDENCE** | none |

Full numbers: `tables/spearman_cldn4_vs_partners.tsv`, `tables/label_inventory.tsv`, `tables/n_table.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore (primary, never-smoker LUAD tumor n=40)

RMA log2 as deposited. Spearman is rank-based. Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; 6/6 present). HOLDS rule (same as PR 229 / PR 313): n≥40, ρ_adj<0, p_adj<0.05. n=40 is the floor of that rule; CIs are wide.

ImmuneScore here is **not** the official R `estimate::estimateScore` ssGSEA. It is the mean of gene-wise z-scores (z within the analysis set) of the Yoshihara 2013 ImmuneSignature (141 genes from `estimate` 1.0.13 `SI_geneset.gmt`).

| pair | n | ρ | 95% CI | p | ρ_adj | p_adj | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| CLDN4 vs **CD8A** | **40** | **+0.034** | −0.312 to +0.393 | 0.83 | +0.055 | 0.74 | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | **40** | **+0.162** | −0.190 to +0.508 | 0.32 | +0.208 | 0.20 | **NO_EVIDENCE** |
| CLDN4 vs **ImmuneScore** | **40** | **−0.032** | −0.382 to +0.309 | 0.85 | −0.004 | 0.98 | **NO_EVIDENCE** |
| CLDN4 Q4 vs Q1 on CD8A | 10 vs 10 | — | — | MWU 0.97 | — | — | null (rank-biserial −0.02) |
| CLDN4 Q4 vs Q1 on CD274 | 10 vs 10 | — | — | MWU 0.68 | — | — | null (rank-biserial +0.12) |
| CLDN4 Q4 vs Q1 on ImmuneScore | 10 vs 10 | — | — | MWU 0.52 | — | — | null (rank-biserial −0.18) |

CLDN4-high is not CD8A-low, CD274-low, or ImmuneScore-low on this never-smoker LUAD array. Under the locked HOLDS rule, nothing holds. Q4 vs Q1 is **10 vs 10**, not 40.

## Extra smoking strata (not the never-smoker claim)

Do not quote the mixed row as a never-smoker result. Smoker n=40 is the same size as the primary set; it is a different biology and is extra only.

| cohort | n | CLDN4–CD8A ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–ImmuneScore ρ (p) | ImmuneScore ρ_adj (p) |
|---|---:|---|---|---|---|
| never-smoker LUAD tumor | **40** | +0.034 (0.83) | +0.162 (0.32) | −0.032 (0.85) | −0.004 (0.98) |
| smoker LUAD tumor extra | 40 | −0.310 (0.052) | −0.293 (0.067) | −0.455 (0.003) | −0.236 (0.15) |
| all LUAD tumor mixed | 80 | −0.225 (0.045) | −0.140 (0.21) | −0.270 (0.015) | −0.160 (0.16) |

Mixed CD8A / ImmuneScore crude negatives are smoker-driven composition. Both vanish after the epithelial residual. That is not a never-smoker confirmation. The smoker ImmuneScore crude HOLDS and the smoker CD274 residual is a borderline HOLDS (ρ_adj=−0.322, p=0.046); neither is this claim.

## What CLDN4 *does* track here

| pair | n | ρ | 95% CI | p |
|---|---:|---:|---|---:|
| CLDN4 vs epithelial mean-z (never-smoker) | 40 | +0.226 | −0.090 to +0.503 | 0.16 |
| CLDN4 vs TACSTD2 (never-smoker) | 40 | +0.264 | −0.058 to +0.552 | 0.10 |
| CLDN4 vs epithelial mean-z (smoker extra) | 40 | +0.642 | 0.384 to 0.823 | 8.1×10⁻⁶ |
| CLDN4 vs TACSTD2 (smoker extra) | 40 | +0.755 | 0.543 to 0.884 | 1.8×10⁻⁸ |

On never-smoker tumors, CLDN4 does **not** significantly track the epithelial mean-z or TACSTD2. On smoker tumors it does. That is a companion observation, not an immune claim, and not a dual-high class.

## What this does not test

- ICI response (not an ICI series).
- Official ESTIMATE ssGSEA ImmuneScore from the R package (this is ImmuneSignature mean-z).
- Pathologist CD8 / PD-L1 IHC.
- A LUAD-wide CLDN4–immune law. This is one public HuGene 1.0 ST series, n=40 never-smoker tumors.
- Dual-high TACSTD2×CLDN4 (not scored).
- OS / stage (not on GEO characteristics).

## Reproduce

```bash
python3 -m pip install -r methods/gse43458_cldn4/requirements.txt
python3 methods/gse43458_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE43458_CLDN4_DATA` (default `/tmp/gse43458_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43458/matrix/GSE43458_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz`
- `https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz` (ImmuneSignature GMT only)

## Files

- `analyze.py` — download, probe map, Spearman / partial Spearman, figures
- `tables/label_inventory.tsv` — public vs missing fields; honest n
- `tables/n_table.tsv` — one-row table
- `tables/histology_counts.tsv` — GEO tissue × smoking
- `tables/probe_confirm.tsv` — HuGene 1.0 ST transcript clusters
- `tables/probe_all_mapped.tsv` — mapped CLDN4 / CD8A / CD274 / TACSTD2 / epithelial probes
- `tables/immunescore_gene_coverage.tsv` — 132/141 ImmuneSignature genes
- `tables/spearman_cldn4_vs_partners.tsv` — never-smoker + smoker extra + mixed tumors
- `tables/highlow_cldn4.tsv` — CLDN4 Q4 vs Q1 (10 vs 10 on the primary set)
- `tables/sample_annotation.tsv` — 110 arrays
- `tables/summary.json`
- `figures/fig1_cldn4_vs_partners.png`
- `figures/fig2_cldn4_forest.png`

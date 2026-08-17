# GSE40419 Korean LUAD RNA-seq (Seo) — CLDN4-only vs CD8A / CD274 / ImmuneScore / IFN

**Additive only. CLDN4-only. No dual-high.** Public Seo Korean lung adenocarcinoma RNA-seq (Seo et al., *Genome Res* 2012, PMID 22975805; GEO [GSE40419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419)). Unit is the **tumor RNA-seq library**. No slide was re-scored. No ICI arm. TACSTD2 is a same-run companion only and is **never a gate**.

The GEO series matrix has **0 expression rows** (`!Sample_data_row_count = 0`). Expression is the author processed RPKM file `GSE40419_LC-87_RPKM_expression.txt.gz` (~10 MB, under 2 GB). No FASTQ was used.

Primary question: **CLDN4 vs CD8A / CD274 / ImmuneScore / IFN** on **tumors only**, with an **epithelial residual**.

## Honest n

Do not write n=200 (paper clinical cohort before RNA-seq subset). Do not write n=164 (87 tumors + 77 adjacent normal). Do not write an expression n from the empty series matrix.

| Item | Public? | n | Note |
|---|---|---:|---|
| GSM in the series matrix | yes | **164** | 87 tumor + 77 adjacent normal titles |
| Expression rows in the series matrix | yes | **0** | empty table; `data_row_count=0` |
| Genes × libraries in author RPKM | yes | **36742 × 164** | RefSeq rows; 87 tumor + 77 `*_nor` |
| Unique HUGO after max-mean | yes | **22427** | isoform collapse on the tumor matrix |
| Tumors (`Lung cancer cells`) | yes | **87** | titles `LC_C*` / `LC_S*`; all unique |
| Adjacent normal (`*_nor`) | yes | **77** | **dropped** (tumor-only slice) |
| Tumors with a paired `*_nor` | yes | **77** | pairing not used |
| Unique tumor titles / patients | yes | **87** | 1 title = 1 library; no duplicate-patient sentence |
| Paper text “200 cancer patients” | paper only | 200 | broader clinical set; **not** the RNA-seq n |
| Stage on tumors | yes | **85** | GEO `Stage`; 2 tumor NA |
| Smoking on tumors | yes | **83** | never / smoker / current / NA |
| OS / ICI / response | no | **0** | surgical atlas, not an ICI series |
| Numeric tumor % / ABSOLUTE purity | no | **0** | not deposited |
| CLDN4 finite | yes | **87** | single RefSeq row `NM_001305` |
| CD8A finite | yes | **87** | protein-coding `NM_001768` (max-mean among NM_); 4 isoforms on the locus |
| CD274 finite | yes | **87** | single RefSeq row `NM_014143` |
| ImmuneScore (Yoshihara Immune141) | computed | **87** | ssGSEA; 141/141 genes present |
| IFN Ayers-6 | computed | **87** | 6/6 genes present |
| Epithelial mean-z (6/6) | computed | **87** | EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7 |
| **Primary pairwise n (tumor LUAD)** | yes | **87** | complete-case CLDN4 + CD8A + CD274 + scores |

The computable public n is **87 tumor libraries**. Use 164 only when counting GSM metadata. Use 77 only for the unused adjacent-normal columns.

## One-row table

| dataset | histology | platform | n tumors | CLDN4–CD8A ρ (p) | epi adj ρ (p) | CLDN4–CD274 ρ (p) | epi adj ρ (p) | CLDN4–ImmuneScore ρ (p) | epi adj ρ (p) | CLDN4–IFN ρ (p) | epi adj ρ (p) | verdict | OS / ICI |
|---|---|---|---:|---|---|---|---|---|---|---|---|---|---|
| GSE40419 Seo | LUAD | HiSeq 2000 author RPKM | **87** | -0.309 (0.00364) | -0.186 (0.0867) | -0.196 (0.0683) | -0.077 (0.479) | -0.381 (0.000276) | -0.143 (0.19) | -0.112 (0.303) | -0.035 (0.752) | **NO_EVIDENCE** | not deposited |

Full numbers: `tables/spearman.tsv`, `tables/label_inventory.tsv`.

## CLDN4 vs CD8A / CD274 / ImmuneScore / IFN (primary)

Author RPKM as deposited. Spearman is rank-based (log2 does not change ρ). Bootstrap 95% CI, 2,000 resamples, seed `20260817`. Partial Spearman residualises on the pan-epithelial mean-z (EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7; **6/6** present). HOLDS rule (same as PR 229 / GSE4573 / GSE10245): n≥40, ρ_adj<0, p_adj<0.05.

| pair | subset | n | ρ | 95% CI | p | ρ_adj epi (p) | verdict |
|---|---|---:|---:|---|---:|---|---|
| CLDN4 vs **CD8A** | LUAD tumor | **87** | **-0.309** | -0.493 to -0.102 | 0.00364 | -0.186 (0.0867) | **NO_EVIDENCE** |
| CLDN4 vs **CD274** | LUAD tumor | **87** | **-0.196** | -0.397 to +0.009 | 0.0683 | -0.077 (0.479) | **NO_EVIDENCE** |
| CLDN4 vs **ImmuneScore** | LUAD tumor | **87** | **-0.381** | -0.576 to -0.167 | 0.000276 | -0.143 (0.19) | **NO_EVIDENCE** |
| CLDN4 vs **IFN (Ayers-6)** | LUAD tumor | **87** | **-0.112** | -0.307 to +0.095 | 0.303 | -0.035 (0.752) | **NO_EVIDENCE** |
| CLDN4 vs IFN type-I ISG (sensitivity) | LUAD tumor | 87 | +0.108 | -0.114 to +0.326 | 0.319 | -0.016 (0.884) | sensitivity |
| CLDN4 vs CD8A `NR_027353` (sensitivity) | LUAD tumor | 87 | -0.306 | -0.492 to -0.098 | 0.0039 | -0.182 (0.093) | sensitivity |

Positive-control (same tumors): CD8A vs ImmuneScore ρ=+0.817 (p=4.71e-22); IFN Ayers-6 vs CD8A ρ=+0.761 (p=1.22e-17).

CLDN4 vs epithelial mean-z: ρ=+0.508 (p=5.24e-07, n=87). The residual asks whether any immune association is more than epithelial content.

**What is inverse unadjusted.** CLDN4 vs CD8A and vs ImmuneScore are negative at n=87 before the residual (CIs exclude 0). Q4 vs Q1 (22 vs 22) matches that direction for CD8A and ImmuneScore. That is a coarsened unadjusted test, not a residual claim.

**What does not HOLDS.** After the specified epithelial residual, every primary pair is NS (CD8A adj p=0.0867; ImmuneScore adj p=0.19; CD274 and IFN already NS unadjusted). HOLDS requires ρ_adj<0 and p_adj<0.05. Do not upgrade the unadjusted inverses into an epithelial-independent immune-low claim.

### Q4 vs Q1 (descriptive)

Quartiles are `pd.qcut(rank(method="first"), 4)` on the 87 tumors with finite CLDN4. Honest n is the quartile arms, not 87.

| endpoint | n Q4 vs Q1 | MWU p | rank-biserial (Q4>Q1) |
|---|---|---:|---:|
| CD8A | 22 vs 22 | 0.0109 | -0.450 |
| CD274 | 22 vs 22 | 0.156 | -0.252 |
| ImmuneScore | 22 vs 22 | 0.00503 | -0.496 |
| IFN Ayers-6 | 22 vs 22 | 0.296 | -0.186 |

## Context (not the claim)

This is a surgical Korean LUAD transcriptome atlas (Seo 2012). It is **not** an ICI-response series. Adjacent-normal columns exist in the RPKM file and were not used. Smoking and stage are deposited for tumors and were not residualised (the specified residual is epithelial). No dual-high (CLDN4-high AND TACSTD2-high) split was run.

## TACSTD2 companion (not the claim)

| pair | subset | n | ρ (p) | ρ_adj epi (p) | verdict |
|---|---|---:|---|---|---|
| TACSTD2 vs CD8A | LUAD tumor | 87 | -0.061 (0.574) | +0.005 (0.965) | companion |
| TACSTD2 vs CD274 | LUAD tumor | 87 | +0.059 (0.587) | +0.121 (0.266) | companion |
| CLDN4 vs TACSTD2 | LUAD tumor | 87 | +0.383 (0.000248) | +0.327 (0.00213) | companion |

Do not treat CLDN4 as interchangeable with TACSTD2. Dual-high was not computed.

## Methods

- **Matrix:** GEO supplementary `GSE40419_LC-87_RPKM_expression.txt.gz` (author RPKM, NCBI build 37.1 / GSNAP). The series matrix expression table is empty and was used only for sample characteristics.
- **Tumor rule:** GEO `source_name = Lung cancer cells` (titles without `_nor`). Adjacent normal dropped.
- **Gene collapse:** max-mean of RefSeq rows that share a HUGO symbol, computed on the tumor matrix. Named genes (CLDN4 / CD8A / CD274 / TACSTD2) take the highest-mean **NM_** accession when one exists.
- **CD8** = protein-coding `CD8A` `NM_001768` (highest-mean NM_ on the locus). The non-coding `NR_027353` isoform is a sensitivity row only. **CD274** = PD-L1 RNA (not protein).
- **ImmuneScore:** Yoshihara 2013 Immune141 ssGSEA (Barbie/GSVA, τ=0.25); ranks scaled to 1…10000. Stromal141 is computed only to form ESTIMATEScore for the coverage audit; the specified residual is epithelial, not ESTIMATE purity.
- **IFN (primary):** Ayers 2017 IFNG 6-gene mean of gene-wise z (`IFNG`, `STAT1`, `CXCL9`, `CXCL10`, `IDO1`, `HLA-DRA`).
- **IFN (sensitivity):** 10-gene type-I ISG mean-z (`ISG15`, `MX1`, `OAS1`, `OAS2`, `IFIT1`, `IFIT3`, `STAT1`, `IRF7`, `IFI27`, `IFI44L`).
- **Epithelial residual:** unweighted mean of gene-wise z for EPCAM / KRT8 / KRT18 / KRT19 / CDH1 / KRT7.
- **Partial Spearman:** Pearson of rank residuals; df = n − 2 − k.
- **Not done:** no FASTQ / recount. No dual-high. No ICI / OS model. No gene-set fishing beyond the pre-specified IFN lists.

## Files

- `tables/spearman.tsv` — all n / ρ / p / CI / epithelial residuals
- `tables/label_inventory.tsv` — honest n
- `tables/sample_annotation.tsv` — per-tumor genes + scores + GEO labels
- `tables/gene_coverage.tsv` / `highlow_cldn4.tsv` / `isoform_used.tsv` (named + set genes only) / `summary.json`
- `figures/fig1_cldn4_vs_cd8a_cd274.png`
- `figures/fig2_cldn4_vs_immunescore_ifn.png`
- `figures/fig3_forest.png`

## Reproduce

```bash
python3 -m pip install -r methods/gse40419_cldn4/requirements.txt
export GSE40419_CLDN4_DATA=/tmp/gse40419_cldn4
python3 methods/gse40419_cldn4/analyze.py
```

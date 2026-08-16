# TACSTD2 surface-gene co-expression in TCGA-READ (B1 analog)

Honest ranking of how well a candidate surface marker co-expresses with an
anchor gene, across the entire cell-surface proteome, in TCGA rectal
adenocarcinoma (READ). Same question and method as the B1_BRCA analog;
only the cohort changes.

**Driving question:** In TCGA-READ, is `CLDN4` the *top* co-expression partner
of `TACSTD2` (TROP2) among surface genes?

**Answer: No — and it is not even close.** Across **94 primary-tumour samples**
and the **2,618 surfaceome genes** present in the matrix, `CLDN4` ranks
**#1579** (Spearman ρ = −0.007, FDR q = 0.99; 39.7th percentile). Pearson is
the same null (ρ = 0.060, rank #951). This is a mid-pack / uncorrelated result,
not a near-miss of #1.

The actual top surface partner of TACSTD2 in READ is `PVRL4` (= NECTIN4,
ρ = 0.604). CLDN4 is constitutively high in this cohort (mean log2 14.04,
range 12.36–15.20) while TACSTD2 varies (mean 8.15, range 3.58–13.72), so
there is almost no CLDN4 dynamic range left to correlate.

| rank | gene | Spearman ρ | note |
| --- | --- | --- | --- |
| 1 | PVRL4 | 0.604 | = **NECTIN4** |
| 2 | LYPD3 | 0.465 | |
| 3 | AMIGO2 | 0.464 | |
| … | … | … | |
| **1579** | **CLDN4** | **−0.007** | focus gene; FDR q = 0.99 |

This is the opposite of the B1_BRCA analog, where CLDN4 was a genuine top-4
partner (ρ = 0.348). The BRCA finding does **not** transfer to READ.

Full machine-readable results live in
[`results/w200/B1_READ/`](results/w200/B1_READ/).

## What "honest ranking" means here

The analysis scores **every** surface gene against the anchor and reports where
the focus gene actually lands. It does not restrict to a hand-picked gene
neighbourhood, and it reports both Spearman (primary) and Pearson so a single
favourable metric cannot be cherry-picked. The `#1` gene is always named
explicitly in the summary next to the focus gene's rank.

## Data sources

- **Expression:** UCSC Xena TCGA hub, dataset `TCGA.READ.sampleMap/HiSeqV2`
  — RNA-seq, log2(norm_count + 1), 20,530 genes × 105 samples. Primary tumours
  (barcode sample-type code `01`, n = 94) are used by default. Same hub/pipeline
  as B1_BRCA. The remaining columns are 10 solid-tissue normals (`11`) and 1
  recurrent tumour (`02`).
- **Surface-gene universe:** the "in silico human surfaceome" (Bausch-Fluck
  et al., *PNAS* 2018; table S3, "in silico surfaceome only" sheet), 2,886
  predicted surface proteins. The workbook is pulled from the
  `steveneschrich/surfaceome` GitHub mirror because the original ETH host now
  serves a single-page app in place of the file.

Gene symbols follow the surfaceome table (legacy HGNC symbols), so a few genes
appear under older names — most notably `PVRL4` for `NECTIN4`.

## Sensitivity (not the primary ranking)

On the newer GDC-harmonized Xena matrix `TCGA-READ.star_tpm` (log2(TPM+0.001),
166 primary tumours), TACSTD2 vs CLDN4 is still null (Spearman ρ = 0.031,
p = 0.69). PVRL4 remains a strong partner (ρ = 0.54). The HiSeqV2 null is
not an artifact of the older freeze or the n = 94 sample size.

## Reproduce

```bash
pip install -r requirements.txt

# 1) fetch inputs into data/ (git-ignored)
python scripts/download_data.py

# 2) run the ranking (defaults reproduce the tables above)
python scripts/coexpression.py
```

Useful flags:

```bash
python scripts/coexpression.py \
  --anchor TACSTD2 --focus CLDN4 \
  --sample-types 01 \
  --window 200 \
  --outdir results/w200/B1_READ
```

`w200` is the reporting window (size of `top200.csv`). It does not change the
full-universe ranking.

## Outputs (`results/w200/B1_READ/`)

| file | contents |
| --- | --- |
| `coexpression_TACSTD2_surfaceome.csv` | full ranked table: Spearman/Pearson ρ, p, BH-FDR q, ranks, mean expression, % expressed for all 2,618 surface genes |
| `top200.csv` | the top-`window` partners (the "w200" window) |
| `summary.json` | machine-readable run parameters + focus-gene result + top-10 |
| `summary.md` | the human-readable verdict and top-15 table |
| `scatter_TACSTD2_vs_CLDN4.png` | anchor-vs-focus expression scatter |
| `top_partners_TACSTD2.png` | top-20 bar chart, focus gene highlighted |

## Notes / caveats

- `HiSeqV2` is gene-level RNA-seq; co-expression here is bulk-tumour mRNA and
  mixes tumour, stromal, and immune compartments — it is not protein-level or
  single-cell co-expression.
- Spearman is the primary metric (robust to the log-count scale and outliers);
  Pearson is reported alongside. Both say CLDN4 is uncorrelated with TACSTD2
  in READ.
- 169 surfaceome genes had zero variance in this small cohort and therefore
  undefined Spearman; they sort to the bottom and do not affect CLDN4's rank.
- The ranking excludes the anchor gene itself.
- READ is a small TCGA cohort relative to BRCA; the CLDN4 null is still
  well-powered enough to reject a BRCA-like ρ ≈ 0.35 (that would be obvious
  at n = 94).

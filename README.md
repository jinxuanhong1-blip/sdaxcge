# TACSTD2 surface-gene co-expression in TCGA-STAD (B1 analog)

Honest ranking of how well a candidate surface marker co-expresses with an
anchor gene, across the entire cell-surface proteome, in TCGA stomach
adenocarcinoma.

**Driving question:** In TCGA-STAD, is `CLDN4` the *top* co-expression partner
of `TACSTD2` (TROP2) among surface genes?

**Answer: No.** Using Spearman correlation across 415 primary-tumour samples
and 2,617 rankable surfaceome genes, `CLDN4` ranks **#33** (ρ = 0.310,
FDR q ≈ 8.7e-9; 98.8th percentile). That is a real, FDR-significant
association, but it is **not** rank 1. Thirty-two surface genes correlate
more strongly with TACSTD2. The actual #1 partner is `PVRL4` (= NECTIN4;
ρ = 0.578), well ahead of CLDN4.

By Pearson correlation `CLDN4` is #26 (r = 0.319). The two metrics agree
that CLDN4 sits in the top ~1–1.3% of the surfaceome, not at the top.

| rank | gene | Spearman ρ | note |
| --- | --- | --- | --- |
| 1 | PVRL4 | 0.578 | = **NECTIN4** |
| 2 | TM4SF1 | 0.438 | |
| 3 | TGFA | 0.416 | |
| 4 | GPR110 | 0.405 | = ADGRF1 (legacy SURFY symbol) |
| 5 | LYPD3 | 0.396 | |
| … | … | … | |
| 32 | CLDN18 | 0.310 | another claudin, just above the focus gene |
| **33** | **CLDN4** | **0.310** | focus gene |

Full machine-readable results live in
[`results/w200/B1_STAD/`](results/w200/B1_STAD/).

## What "honest ranking" means here

The analysis scores **every** surface gene against the anchor and reports
where the focus gene actually lands. It does not restrict to a hand-picked
gene neighbourhood, and it reports both Spearman (primary) and Pearson so a
single favourable metric cannot be cherry-picked. The `#1` gene is always
named explicitly in the summary next to the focus gene's rank.

One surfaceome gene (`OR2T29`) had zero variance in this matrix (all zeros)
and therefore has an undefined correlation; it is left in the full table
with missing ranks and is excluded from the 2,617-gene rank denominator.

## Data sources

- **Expression:** UCSC Xena TCGA hub, dataset `TCGA.STAD.sampleMap/HiSeqV2`
  — RNA-seq, log2(norm_count + 1), rows are HGNC gene symbols, 20,530 genes ×
  450 samples. Primary tumours (barcode sample-type code `01`; n = 415) are
  used by default.
- **Surface-gene universe:** the "in silico human surfaceome" (Bausch-Fluck
  et al., *PNAS* 2018; table S3, "in silico surfaceome only" sheet), 2,886
  predicted surface proteins. The workbook is pulled from the
  `steveneschrich/surfaceome` GitHub mirror because the original ETH host now
  serves a single-page app in place of the file.

Gene symbols follow the surfaceome table (legacy HGNC symbols), so a few
genes appear under older names — most notably `PVRL4` for `NECTIN4` and
`GPR110` for `ADGRF1`.

## Reproduce

```bash
pip install -r requirements.txt

# 1) fetch inputs into data/ (git-ignored; ~25 MB expression matrix)
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
  --outdir results/w200/B1_STAD
```

## Outputs (`results/w200/B1_STAD/`)

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
  Pearson is reported alongside. Both reject the "CLDN4 is #1" claim.
- The ranking excludes the anchor gene itself.
- Compared with the BRCA analog (CLDN4 Spearman rank #4, ρ = 0.348), the
  STAD association is weaker and farther from the top. Do not treat the
  BRCA rank as transferable to stomach adenocarcinoma.

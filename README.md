# TACSTD2 surface-gene co-expression in TCGA-BRCA (B1 analog)

Honest ranking of how well a candidate surface marker co-expresses with an
anchor gene, across the entire cell-surface proteome, in TCGA breast cancer.

**Driving question:** In TCGA-BRCA, is `CLDN4` the *top* co-expression partner
of `TACSTD2` (TROP2) among surface genes?

**Answer: No — but it is very near the top.** Using Spearman correlation across
1,097 primary-tumour samples and the 2,618 surfaceome genes present in the
matrix, `CLDN4` ranks **#4** (ρ = 0.348, FDR q ≈ 1e-29; 99.9th percentile).
The genes ahead of it are:

| rank | gene | Spearman ρ | note |
| --- | --- | --- | --- |
| 1 | EFNA1 | 0.429 | ephrin-A1 |
| 2 | PVRL4 | 0.351 | = **NECTIN4** (enfortumab vedotin target) |
| 3 | EFNA4 | 0.350 | ephrin-A4 |
| 4 | **CLDN4** | 0.348 | claudin-4 (focus gene) |
| 5 | TM4SF1 | 0.341 | |

By Pearson correlation `CLDN4` is #3. So the headline "CLDN4 is the top
TACSTD2-correlated surface gene in BRCA" is **not** literally true; the honest
statement is that CLDN4 is a top-4 (top ~0.15%) surface-gene partner of TACSTD2,
just behind EFNA1 and NECTIN4.

Full machine-readable results live in
[`results/w200/B1_BRCA/`](results/w200/B1_BRCA/).

## What "honest ranking" means here

The analysis scores **every** surface gene against the anchor and reports where
the focus gene actually lands. It does not restrict to a hand-picked gene
neighbourhood, and it reports both Spearman (primary) and Pearson so a single
favourable metric cannot be cherry-picked. The `#1` gene is always named
explicitly in the summary next to the focus gene's rank.

## Data sources

- **Expression:** UCSC Xena TCGA hub, dataset `TCGA.BRCA.sampleMap/HiSeqV2`
  — RNA-seq, log2(norm_count + 1), rows are HGNC gene symbols, 20,530 genes ×
  1,218 samples. Primary tumours (barcode sample-type code `01`) are used by
  default.
- **Surface-gene universe:** the "in silico human surfaceome" (Bausch-Fluck
  et al., *PNAS* 2018; table S3, "in silico surfaceome only" sheet), 2,886
  predicted surface proteins. The workbook is pulled from the
  `steveneschrich/surfaceome` GitHub mirror because the original ETH host now
  serves a single-page app in place of the file.

Gene symbols follow the surfaceome table (legacy HGNC symbols), so a few genes
appear under older names — most notably `PVRL4` for `NECTIN4`.

## Reproduce

```bash
pip install -r requirements.txt

# 1) fetch inputs into data/ (git-ignored; ~64 MB expression matrix)
python scripts/download_data.py

# 2) run the ranking (defaults reproduce the tables above)
python scripts/coexpression.py
```

Useful flags:

```bash
python scripts/coexpression.py \
  --anchor TACSTD2 --focus CLDN4 \
  --sample-types 01 \        # TCGA sample-type codes to keep
  --window 200 \             # size of the top-N companion CSV ("w200")
  --outdir results/w200/B1_BRCA
```

## Outputs (`results/w200/B1_BRCA/`)

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
  Pearson is reported alongside. The two agree that CLDN4 sits at rank 3–4.
- The ranking excludes the anchor gene itself.

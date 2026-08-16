# TACSTD2 surface-gene co-expression in TCGA-CESC (B1 analog)

Honest ranking of how well a candidate surface marker co-expresses with an
anchor gene, across the entire cell-surface proteome, in TCGA cervical
squamous cell carcinoma and endocervical adenocarcinoma (CESC).

This is the CESC analog of the BRCA B1 surface-rank analysis
(`results/w200/B1_BRCA/`): same question, same surfaceome universe, same
HiSeqV2 / Spearman protocol, different cohort.

**Driving question:** In TCGA-CESC, is `CLDN4` the *top* co-expression partner
of `TACSTD2` (TROP2) among surface genes?

**Answer:** see [`results/w200/B1_CESC/summary.md`](results/w200/B1_CESC/summary.md)
after running the pipeline. The ranking is not pre-committed; it is computed
from the public matrix.

## What "honest ranking" means here

The analysis scores **every** surface gene against the anchor and reports where
the focus gene actually lands. It does not restrict to a hand-picked gene
neighbourhood, and it reports both Spearman (primary) and Pearson so a single
favourable metric cannot be cherry-picked. The `#1` gene is always named
explicitly in the summary next to the focus gene's rank.

## Data sources

- **Expression:** UCSC Xena TCGA hub, dataset `TCGA.CESC.sampleMap/HiSeqV2`
  — RNA-seq, log2(norm_count + 1), rows are HGNC gene symbols. Primary tumours
  (barcode sample-type code `01`) are used by default.
- **Surface-gene universe:** the "in silico human surfaceome" (Bausch-Fluck
  et al., *PNAS* 2018; table S3, "in silico surfaceome only" sheet), 2,886
  predicted surface proteins. The workbook is pulled from the
  `steveneschrich/surfaceome` GitHub mirror because the original ETH host now
  serves a Git-LFS pointer / SPA in place of the file. The downloaded xlsx
  sha256 matches that official LFS pointer.

Gene symbols follow the surfaceome table (legacy HGNC symbols).

## Reproduce

```bash
pip install -r requirements.txt

# 1) fetch inputs into data/ (git-ignored)
python scripts/download_data.py

# 2) run the ranking (defaults write results/w200/B1_CESC/)
python scripts/coexpression.py
```

Useful flags:

```bash
python scripts/coexpression.py \
  --anchor TACSTD2 --focus CLDN4 \
  --sample-types 01 \
  --window 200 \
  --outdir results/w200/B1_CESC
```

`w200` is the reporting window (size of `top200.csv`). It does not change the
full-universe ranking.

## Outputs (`results/w200/B1_CESC/`)

| file | contents |
| --- | --- |
| `coexpression_TACSTD2_surfaceome.csv` | full ranked table: Spearman/Pearson ρ, p, BH-FDR q, ranks, mean expression, % expressed |
| `top200.csv` | the top-`window` partners |
| `summary.json` | machine-readable run parameters + focus-gene result + top-10 |
| `summary.md` | the human-readable verdict and top-15 table |
| `scatter_TACSTD2_vs_CLDN4.png` | anchor-vs-focus expression scatter |
| `top_partners_TACSTD2.png` | top-20 bar chart, focus gene highlighted |

## Notes / caveats

- `HiSeqV2` is gene-level RNA-seq; co-expression here is bulk-tumour mRNA and
  mixes tumour, stromal, and immune compartments — it is not protein-level or
  single-cell co-expression.
- Spearman is the primary metric (robust to the log-count scale and outliers);
  Pearson is reported alongside.
- The ranking excludes the anchor gene itself.
- CESC is a mixed squamous / adenocarcinoma cohort; this first pass does not
  split histology.

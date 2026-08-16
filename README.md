# TACSTD2 surface-gene co-expression in TCGA-CESC (B1 analog)

Honest ranking of how well a candidate surface marker co-expresses with an
anchor gene, across the entire cell-surface proteome, in TCGA cervical
squamous cell carcinoma and endocervical adenocarcinoma (CESC).

This is the CESC analog of the BRCA B1 surface-rank analysis
(`results/w200/B1_BRCA/`): same question, same surfaceome universe, same
HiSeqV2 / Spearman protocol, different cohort.

**Driving question:** In TCGA-CESC, is `CLDN4` the *top* co-expression partner
of `TACSTD2` (TROP2) among surface genes?

**Answer: No.** In the pooled primary-tumour cohort (the BRCA-analog protocol)
`CLDN4` ranks **#193 of 2,575** (ρ = 0.279, FDR q = 5.7e-6). That is
significant but mid-pack — not a top partner. The actual #1 is `PVRL4`
(= **NECTIN4**, ρ = 0.690).

The pooled rank is **histology-diluted**. CESC is ~82% squamous. Restricting
to squamous tumours moves CLDN4 to **#12** (ρ = 0.491). In adenocarcinoma
(n = 52) the link is weak and not FDR-significant (ρ = 0.203, q = 0.40,
rank #317). Even in squamous tumours CLDN4 is not #1 (`MUC21` leads).

| setting | n | CLDN4 Spearman rank | ρ | FDR q | actual #1 |
| --- | ---: | ---: | ---: | ---: | --- |
| **primary (protocol)** | 303 | **#193 / 2575** | 0.279 | 5.7e-6 | PVRL4 (0.690) |
| squamous only | 250 | #12 / 2560 | 0.491 | 2.9e-14 | MUC21 (0.601) |
| adeno only | 52 | #317 / 2464 | 0.203 | 0.40 | PVRL4 (0.725) |

Compare BRCA B1 analog: CLDN4 was #4 (ρ = 0.348). CESC does **not** reproduce
that near-top rank unless the cohort is restricted to squamous histology, and
even then it is #12, not #1.

Full machine-readable results live in
[`results/w200/B1_CESC/`](results/w200/B1_CESC/).

## What "honest ranking" means here

The analysis scores **every** surface gene against the anchor and reports where
the focus gene actually lands. It does not restrict to a hand-picked gene
neighbourhood, and it reports both Spearman (primary) and Pearson so a single
favourable metric cannot be cherry-picked. The `#1` gene is always named
explicitly in the summary next to the focus gene's rank.

The primary ranking pools all primary tumours, matching the BRCA analog. The
histology split is reported as robustness, not as a replacement, because
hiding the squamous-only #12 or the pooled #193 would both be selective.

## Data sources

- **Expression:** UCSC Xena TCGA hub, dataset `TCGA.CESC.sampleMap/HiSeqV2`
  — RNA-seq, log2(norm_count + 1), 20,530 genes × 308 samples. Primary tumours
  (barcode sample-type code `01`; n = 303) are used by default. Two metastases
  (`06`) and three normals (`11`) are excluded.
- **Surface-gene universe:** the "in silico human surfaceome" (Bausch-Fluck
  et al., *PNAS* 2018; table S3, "in silico surfaceome only" sheet), 2,886
  predicted surface proteins. 2,618 symbols are present in HiSeqV2; 43 with
  undefined Spearman (constant / zero-variance) are excluded from ranking,
  leaving 2,575. The workbook is pulled from the `steveneschrich/surfaceome`
  GitHub mirror because the original ETH host now serves a Git-LFS pointer;
  the downloaded xlsx sha256 matches that official LFS pointer.
- **Histology:** Xena `TCGA.CESC.sampleMap/CESC_clinicalMatrix`, field
  `histological_type`. Squamous = "Cervical Squamous Cell Carcinoma"
  (n = 250 among primary tumours in the matrix). Adeno = any
  `histological_type` containing "Adeno" (n = 52). One primary sample has
  missing histology.

Gene symbols follow the surfaceome table (legacy HGNC symbols), so a few genes
appear under older names — most notably `PVRL4` for `NECTIN4`.

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
| `coexpression_TACSTD2_surfaceome.csv` | full ranked table (pooled primary tumours) |
| `coexpression_TACSTD2_surfaceome_squamous.csv` | same ranking, squamous only |
| `coexpression_TACSTD2_surfaceome_adeno.csv` | same ranking, adenocarcinoma only |
| `top200.csv` | the top-`window` partners of the pooled ranking |
| `summary.json` | machine-readable run parameters + focus-gene result + histology split |
| `summary.md` | the human-readable verdict, top-15, and histology note |
| `histology_split.json` | squamous / adeno focus ranks |
| `scatter_TACSTD2_vs_CLDN4.png` | pooled anchor-vs-focus scatter |
| `top_partners_TACSTD2.png` | pooled top-20 bar chart, focus gene highlighted |

## Notes / caveats

- `HiSeqV2` is gene-level RNA-seq; co-expression here is bulk-tumour mRNA and
  mixes tumour, stromal, and immune compartments — it is not protein-level or
  single-cell co-expression.
- Spearman is the primary metric (robust to the log-count scale and outliers);
  Pearson agrees that pooled CLDN4 is ~#200 (Pearson rank #204, ρ = 0.258).
- The ranking excludes the anchor gene itself and 43 constant surface genes.
- Among claudins in the pooled set, `CLDN1` outranks `CLDN4` (#24, ρ = 0.522);
  `CLDN3` is inversely correlated (ρ = −0.294).
- Top pooled partners (PVRL4/NECTIN4, DUOXA1, CLCA4, TMPRSS11D, SDC1, LYPD3,
  DSG3, DSC3) are themselves squamous/epithelial-surface genes. The TACSTD2
  neighbourhood in CESC is a squamous-epithelial program, not a CLDN4-specific
  link.
- Adenocarcinoma n = 52 is small; the null there is a real lack of evidence,
  not proof of no association.
- No ICI / ADC treatment annotation is used. This is expression co-ranking
  only.

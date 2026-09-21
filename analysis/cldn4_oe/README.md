# Public CLDN4-high versus CLDN4-low transcriptomes

## Search

GEO and SRA were searched for a CLDN4 (or Cldn4) cDNA overexpression RNA-seq series. The queries and the series that were set aside are in `tables/search_log.tsv`.

There is no public RNA-seq series in which CLDN4 cDNA was overexpressed and the transcriptome was deposited. SRA queries for `CLDN4 overexpression` and `Cldn4 overexpression` return no runs. The only GEO series whose title says CLDN4 overexpression is **GSE22493**, a 2010 two-color microarray. Its design is parental SKOV3 cells, called “CLDN4 overexpression (control)” because the line already expresses CLDN4, versus CLDN4 siRNA.

The open RNA-seq that changes CLDN4 itself is **GSE207704** (PRJNA856719): parental MCF7 and T47D versus CLDN4 knockout, two replicates each (Murakami, Chiba and colleagues, Breast Cancer Research 2023). The paper’s rescue lines were profiled by RT-qPCR and are not in the RNA-seq series.

## Contrast

Both datasets are read in the CLDN4-high direction, parental (or the array control channel) over CLDN4-low. That is the direction an overexpression contrast would have: NHEJ up and interferon down.

Pre-specified sets, frozen from Enrichr (MSigDB Hallmark 2020, KEGG 2021 Human, Reactome 2022):

- KEGG non-homologous end-joining
- Hallmark interferon alpha response
- Hallmark interferon gamma response

Sensitivity sets: Reactome NHEJ with histone genes removed, and Reactome interferon alpha/beta signaling with the IFN ligand genes removed. Symbols are matched through the frozen NCBI synonym table in `genesets/hgnc_aliases.tsv`.

## GSE207704 deposited FPKM

`GSE207704_CLDN4_RNAseq.txt.gz` is a Cufflinks FPKM table with one column per genotype. Replicates are already collapsed, and many genes have no row (MX1, B2M, IFIT3, EGFR among them). Log2 fold change is log2((WT + 0.1) / (KO + 0.1)) for the dominant locus of each symbol with FPKM at least 1 in one genotype. The rank p-value is a two-sided Mann-Whitney test of set genes against the other detected genes. It is a gene-ranking statistic on a collapsed table, not a biological-replicate test. Benjamini-Hochberg q values are computed across the three primary sets inside each cell line.

CLDN4 mRNA is lower in the knockout but still present: MCF7 85.9 versus 50.8 FPKM, T47D 42.5 versus 20.4 FPKM.

| Cell line | Set | Genes detected | Mean log2(WT/KO) | Median | Fraction > 0 | Rank p | q |
|---|---|---:|---:|---:|---:|---:|---:|
| MCF7 | KEGG NHEJ | 12 / 13 | −0.029 | −0.058 | 0.417 | 0.87 | 0.87 |
| MCF7 | Hallmark IFN alpha | 58 / 97 | −0.065 | −0.044 | 0.431 | 0.85 | 0.87 |
| MCF7 | Hallmark IFN gamma | 116 / 200 | +0.004 | −0.023 | 0.474 | 0.87 | 0.87 |
| T47D | KEGG NHEJ | 12 / 13 | −0.067 | −0.027 | 0.333 | 0.72 | 0.72 |
| T47D | Hallmark IFN alpha | 63 / 97 | +0.370 | +0.211 | 0.683 | 5.2×10⁻⁴ | 0.0016 |
| T47D | Hallmark IFN gamma | 119 / 200 | +0.218 | +0.086 | 0.605 | 0.019 | 0.029 |

NHEJ sits on zero in both lines. Interferon alpha and gamma are higher in T47D parental cells than in the CLDN4 knockout. MCF7 does not show that shift. The global median log2 fold change is about zero in both lines (−0.027 in MCF7, −0.006 in T47D), so the T47D interferon shift is not a genome-wide offset. About a third of each hallmark interferon list has no row in this FPKM table.

## GSE22493 microarray

Channel 1 (Cy3) is the parental control and channel 2 (Cy5) is CLDN4 siRNA. The deposited VALUE tracks CLDN4: probe 17169 is missing, −1.74, and −0.71 across the three arrays, so VALUE is treated as log2(siRNA / parental). The CLDN4-high contrast is the sign flip of that ratio. Probes are averaged within a gene symbol. The replicate test is a one-sample t-test of the three array-level set means.

| Set | Array means, parental / siRNA | Mean | t p | q |
|---|---|---:|---:|---:|
| KEGG NHEJ | +0.719, −0.535, +0.290 | +0.158 | 0.71 | 0.99 |
| Hallmark IFN alpha | +0.263, −0.409, +0.157 | +0.004 | 0.99 | 0.99 |
| Hallmark IFN gamma | +0.436, −0.570, +0.410 | +0.092 | 0.81 | 0.99 |

The second array flips sign for every primary set. The pathway means are compatible with no shift. CLDN4 itself moves in the knockdown direction on the two arrays where the probe is present.

## What this adds

The open CLDN4-high versus CLDN4-low profiles do not show NHEJ up. The one interferon shift that survives the primary-set q filter is in the other direction, and only in T47D: parental cells, which retain more CLDN4 mRNA, score higher on interferon genes than the knockout. The SKOV3 siRNA array does not move either set. These are loss-of-function profiles, not a cDNA overexpression RNA-seq, and the deposited FPKM table collapses the two replicates.

## Reproduce

```bash
python3 analysis/cldn4_oe/analyze_public_cldn4.py --geo /tmp/geo
```

The script expects the GEO supplementary files under `--geo` (`GSE207704_CLDN4_RNAseq.txt.gz`, `GSE22493_series_matrix.txt.gz`, `GPL10555_family.soft.gz`). If Salmon `quant.sf` files are present in `--geo/quant/{MCF7,T47D}_{WT,KO}_rep{1,2}/` and `--geo/ref/pc.fa` is the Ensembl cDNA used to build the index, the script also writes replicate-level Salmon scores. Dependencies: numpy, scipy.

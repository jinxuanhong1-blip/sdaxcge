# Public CLDN4-high versus CLDN4-low transcriptomes

## Search

GEO and SRA were searched for a CLDN4 (or Cldn4) cDNA overexpression RNA-seq series. Queries and the series that were set aside are in `tables/search_log.tsv`.

There is no public RNA-seq series in which CLDN4 cDNA was overexpressed and the transcriptome was deposited. SRA queries for `CLDN4 overexpression` and `Cldn4 overexpression` return no runs. The only GEO series whose title says CLDN4 overexpression is **GSE22493**, a 2010 two-color microarray. Its design is parental SKOV3 cells, called “CLDN4 overexpression (control)” because the line already expresses CLDN4, versus CLDN4 siRNA.

The open RNA-seq that changes CLDN4 itself is **GSE207704** (PRJNA856719): parental MCF7 and T47D versus CLDN4 knockout, two replicates each (Murakami, Chiba and colleagues, Breast Cancer Research 2023). The paper’s rescue lines were profiled by RT-qPCR and are not in the RNA-seq series.

## Contrast

Both datasets are read in the CLDN4-high direction, parental (or the array control channel) over CLDN4-low. That is the direction an overexpression contrast would have: NHEJ up and interferon down.

Pre-specified sets, frozen from Enrichr (MSigDB Hallmark 2020, KEGG 2021 Human, Reactome 2022):

- KEGG non-homologous end-joining
- Hallmark interferon alpha response
- Hallmark interferon gamma response

Sensitivity sets: Reactome NHEJ with histone genes removed, and Reactome interferon alpha/beta signaling with the IFN ligand genes removed. Symbols are matched through `genesets/hgnc_aliases.tsv`. When a query symbol is absent, synonyms are tried in alphabetical order so the match does not depend on hash order.

## GSE207704 Salmon counts

The GEO FPKM table collapses the two replicates and has no row for many genes (MX1, B2M, IFIT3, EGFR). The eight SRA runs were therefore quantified separately: 51 bp single-end reads, Salmon 1.12.1, Ensembl 113 protein-coding cDNA, no genome decoys. Counts were summed to gene symbols and converted to CPM. Log2 fold change is log2((mean CPM of the two parental replicates + 0.1) / (mean CPM of the two knockout replicates + 0.1)) for genes with mean CPM at least 1 in one genotype. The rank p-value is a two-sided Mann-Whitney test of those gene log2 fold changes against the other detected genes. Benjamini-Hochberg q values are across the three primary sets inside each cell line. That rank test treats genes as the sampling unit. The replicate check is the mean log2(CPM + 0.1) of the set in each library.

Mapping to the protein-coding transcriptome was 36.3–41.2% (`tables/gse207704_salmon_mapping.tsv`). CLDN4 mRNA is lower in the knockout in every library and is not abolished. CPM is about 142 versus 84 in MCF7 and about 62 versus 31 in T47D (`tables/gse207704_cldn4_salmon_cpm.tsv`).

| Cell line | Set | Genes | Mean log2(WT/KO) | Median | Fraction > 0 | Rank p | q |
|---|---|---:|---:|---:|---:|---:|---:|
| MCF7 | KEGG NHEJ | 12 / 13 | −0.062 | −0.023 | 0.417 | 0.53 | 0.57 |
| MCF7 | Hallmark IFN alpha | 70 / 97 | +0.015 | +0.008 | 0.500 | 0.57 | 0.57 |
| MCF7 | Hallmark IFN gamma | 141 / 200 | +0.039 | +0.054 | 0.525 | 0.21 | 0.57 |
| T47D | KEGG NHEJ | 12 / 13 | +0.032 | −0.009 | 0.417 | 0.69 | 0.69 |
| T47D | Hallmark IFN alpha | 74 / 97 | +0.285 | +0.137 | 0.662 | 0.0010 | 0.0029 |
| T47D | Hallmark IFN gamma | 143 / 200 | +0.170 | +0.094 | 0.608 | 0.0019 | 0.0029 |

NHEJ is centered on zero in both lines. The replicate means overlap between parental and knockout (T47D delta of mean log2 CPM +0.031; MCF7 −0.060).

Interferon is flat in MCF7. In T47D both parental libraries sit above both knockout libraries: interferon-alpha mean log2(CPM+0.1) is 3.94 and 3.89 in parental versus 3.66 and 3.63 in knockout (delta +0.271). Interferon-gamma is 4.09 and 4.06 versus 3.89 and 3.90 (delta +0.177). The genome-wide median log2 fold change is −0.017 in MCF7 and −0.022 in T47D, so the T47D interferon shift is not a global offset.

The deposited FPKM table, scored the same way, agrees in direction: T47D interferon alpha mean log2 fold change +0.362 (63 genes present), NHEJ about −0.07 and −0.03. It is the secondary table because replicates are pre-collapsed and gene coverage is lower.

## GSE22493 microarray

Channel 1 (Cy3) is the parental control and channel 2 (Cy5) is CLDN4 siRNA. The deposited VALUE tracks CLDN4: probe 17169 is missing, −1.74, and −0.71 across the three arrays, so VALUE is treated as log2(siRNA / parental). The CLDN4-high contrast is the sign flip of that ratio. Probes are averaged within a gene symbol. The replicate test is a one-sample t-test of the three array-level set means.

| Set | Array means, parental / siRNA | Mean | t p | q |
|---|---|---:|---:|---:|
| KEGG NHEJ | +0.719, −0.535, +0.290 | +0.158 | 0.71 | 0.99 |
| Hallmark IFN alpha | +0.263, −0.409, +0.157 | +0.004 | 0.99 | 0.99 |
| Hallmark IFN gamma | +0.426, −0.565, +0.399 | +0.087 | 0.81 | 0.99 |

The second array flips sign for every primary set. The pathway means are compatible with no shift. CLDN4 itself moves in the knockdown direction on the two arrays where the probe is present.

## What this adds

The open CLDN4-high versus CLDN4-low profiles do not show NHEJ up. Interferon is not down. In T47D, parental cells, which keep more CLDN4 mRNA, score higher on interferon-alpha and interferon-gamma genes than the knockout, and the two replicates agree. MCF7 and the SKOV3 siRNA array do not move either set. These are loss-of-function profiles, not a cDNA overexpression RNA-seq.

## Reproduce

```bash
python3 analysis/cldn4_oe/analyze_public_cldn4.py --geo /tmp/geo
```

The script expects the GEO supplementary files under `--geo` (`GSE207704_CLDN4_RNAseq.txt.gz`, `GSE22493_series_matrix.txt.gz`, `GPL10555_family.soft.gz`). Salmon output is read from `--geo/quant/{MCF7,T47D}_{WT,KO}_rep{1,2}/quant.sf` with transcript-to-gene names taken from `--geo/ref/pc.fa` (Ensembl cDNA headers containing `gene_symbol:`). Dependencies: numpy, scipy.

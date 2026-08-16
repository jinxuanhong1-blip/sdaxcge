# DepMap lung TACSTD2 / IFN

Reproducible analysis of TACSTD2 expression and dependency against Hallmark
interferon signatures and IFN-gene CRISPR dependencies in DepMap lung cancer
cell lines. The generated interpretation, including negative results and
limitations, is in [`REPORT.md`](REPORT.md).

## Data

This analysis uses the stable DepMap Public 24Q4 Figshare archive
([DOI](https://doi.org/10.25452/figshare.plus.27993248.v1)). Raw matrices are
about 936 MB and are gitignored.

```bash
mkdir -p results/w200/DepMap_IFN/data
wget https://ndownloader.figshare.com/files/51065297 \
  -O results/w200/DepMap_IFN/data/Model.csv
wget https://ndownloader.figshare.com/files/51065489 \
  -O results/w200/DepMap_IFN/data/OmicsExpressionProteinCodingGenesTPMLogp1.csv
wget https://ndownloader.figshare.com/files/51064667 \
  -O results/w200/DepMap_IFN/data/CRISPRGeneEffect.csv
```

The script verifies the official MD5 checksums before reading the files.

## Run

```bash
python3 -m pip install -r results/w200/DepMap_IFN/requirements.txt
python3 results/w200/DepMap_IFN/analyze.py
```

The primary cohort is fixed as malignant DepMap models satisfying:

- `OncotreeLineage == "Lung"`
- `ModelType == "Cell Line"`
- `OncotreePrimaryDisease != "Non-Cancerous"`

The two primary tests correlate TACSTD2 log2(TPM+1) with the MSigDB Hallmark
IFN-alpha and IFN-gamma scores. A score is the mean of gene-wise z-scores
within the lung cohort. Spearman correlations, bootstrap confidence intervals,
and BH correction across the two signatures are reported.

Secondary analyses test TACSTD2 Chronos gene effect against both signatures,
adjust expression associations for model type, estimate major-subtype
associations, and screen the fixed union of Hallmark IFN genes for CRISPR
gene-effect correlation with TACSTD2. Candidate-screen FDR is controlled
separately for each TACSTD2 predictor.

Hallmark definitions were retrieved from MSigDB:

- [HALLMARK_INTERFERON_ALPHA_RESPONSE](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE)
- [HALLMARK_INTERFERON_GAMMA_RESPONSE](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_GAMMA_RESPONSE)

The analysis is cross-sectional and does not support causal, synthetic-lethal,
patient-response, or ADC-efficacy claims.

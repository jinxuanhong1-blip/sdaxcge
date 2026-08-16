#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pydeseq2_template.py  -  DESeq2 in pure Python (no R) for KD/KO RNA-seq.

Use this when you only have Python. It reproduces the DESeq2 negative-binomial
workflow via PyDESeq2. Requires RAW integer counts + >= 2 replicates per group.

    pip install pydeseq2 pandas numpy

INPUT
    counts.tsv   genes (rows) x samples (cols), RAW integer counts.
    coldata.tsv  index = sample, columns include `condition` (and `batch`).

USAGE
    python pydeseq2_template.py counts.tsv coldata.tsv out_dir
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def main(counts_file, coldata_file, out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(counts_file, sep="\t", index_col=0)
    counts = counts.round().astype(int).T          # PyDESeq2 wants samples x genes
    meta = pd.read_csv(coldata_file, sep="\t", index_col=0)
    meta = meta.loc[counts.index]                   # align order (critical)

    # Reference level = control. Design blocks on batch if present & informative.
    factors = ["condition"]
    if "batch" in meta.columns and meta["batch"].nunique() > 1:
        factors = ["batch", "condition"]

    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design_factors=factors,
        ref_level=["condition", "control"],
        refit_cooks=True,
    )
    # light pre-filter: drop genes with <10 total counts
    dds = dds[:, counts.sum(axis=0) >= 10]
    dds.deseq2()

    stat = DeseqStats(dds, contrast=["condition", "knockdown", "control"])
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "gene"
    res = res.sort_values("padj")
    res.to_csv(out / "pydeseq2_results.tsv", sep="\t")

    # ranked list for GSEA (gseapy prerank)
    rnk = res["log2FoldChange"].replace([np.inf, -np.inf], np.nan).dropna()
    rnk.sort_values(ascending=False).to_csv(
        out / "ranked_for_gsea.rnk", sep="\t", header=False
    )

    focus = ["CLDN4", "TACSTD2", "EPCAM", "CLDN1", "CLDN3", "CLDN7"]
    print(res.loc[res.index.intersection(focus), ["log2FoldChange", "padj"]])
    print("Done ->", out)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else "counts.tsv",
         a[1] if len(a) > 1 else "coldata.tsv",
         a[2] if len(a) > 2 else "pydeseq2_out")

#!/usr/bin/env python3
"""GSE289287 - human T-47D breast, Trop-2 KO vs WT xenografts (author DESeq2).

Uses the author-provided DESeq2 tables directly (log2FC, pvalue, padj).
Contrast Trop2KO_tumors_vs_WT: log2FC>0 => up in Trop-2 KO.
Also analyses DSG2 KO (desmosome partner) tumours & cells for comparison.
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis_common as AC

DATA = AC.DATA
jobs = [
    ("GSE289287", "Human T-47D Trop-2 KO vs WT xenografts (author DESeq2)",
     "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz",
     "4 Trop-2 KO vs 3 WT xenografts. log2FC>0 = up in Trop-2 KO."),
    ("GSE289287_DSG2tumor", "Human T-47D DSG2 KO vs WT xenografts (desmosome partner, author DESeq2)",
     "GSE289287_DESeq2-DSG2KO_tumors_vs_WT.tsv.gz",
     "DSG2 KO vs WT xenografts (comparison: desmosomal partner)."),
    ("GSE289287_DSG2cell", "Human T-47D DSG2 KO vs WT cells (desmosome partner, author DESeq2)",
     "GSE289287_DESeq2-DSG2KO_cells_vs_DSG2WT.tsv.gz",
     "DSG2 KO vs WT cells (comparison: desmosomal partner)."),
]

for dsid, desc, fname, notes in jobs:
    df = pd.read_csv(DATA / fname, sep="\t")
    print("="*80)
    print(dsid, df.shape, "cols:", list(df.columns)[:10])
    summary, pergene, setdf = AC.run_deseq2_table(
        dsid, desc, df,
        symbol_col="Feature_name", l2fc_col="log2FoldChange",
        p_col="pvalue", padj_col="padj", stat_col="stat", basemean_col="baseMean",
        notes=notes,
    )
    print("QC:", summary["perturbation_qc"])
    print("CLDN4:", summary["cldn4"])
    print(setdf[["gene_set", "n_detected", "comp_direction", "comp_rank_biserial",
                 "comp_p", "sc_mean_log2FC", "sc_n_up", "sc_n_down", "sc_wilcoxon_p"]].to_string(index=False))

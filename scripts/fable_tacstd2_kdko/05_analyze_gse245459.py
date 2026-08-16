#!/usr/bin/env python3
"""GSE245459 - human ovarian cancer, shTACSTD2 vs shNC (RNA-seq FPKM).

Primary contrast (no drug): sh1/2/3 (shTACSTD2) vs shNC1/2/3.
Stats: Welch t-test + BH FDR + Mann-Whitney on log2(FPKM + 1). n=3/3 -> exploratory.
log2FC>0 => up in shTACSTD2.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis_common as AC

DATA = AC.DATA
df = pd.read_csv(DATA / "GSE245459_fpkm.anno.txt.gz", sep="\t", low_memory=False)
print("raw shape:", df.shape)

case_cols = ["sh1", "sh2", "sh3"]
ctrl_cols = ["shNC1", "shNC2", "shNC3"]
symbols = df["GeneName"].astype(str).values
expr = df[case_cols + ctrl_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values
mat = np.log2(expr + 1.0)

# drop all-zero genes, collapse to gene level
nz = mat.sum(axis=1) > 0
symbols, mat = symbols[nz], mat[nz]
symbols, mat, _ = AC.collapse_to_gene(symbols, mat)
print("gene-level matrix:", mat.shape)

case_idx = [0, 1, 2]
ctrl_idx = [3, 4, 5]
summary, pergene, setdf = AC.run_expression_dataset(
    "GSE245459",
    "Human ovarian cancer cells, shTACSTD2 vs shNC (FPKM)",
    symbols, mat, case_idx, ctrl_idx,
    value_kind="log2(FPKM+1)",
    notes="3 shTACSTD2 vs 3 shNC (no drug). n=3/3 exploratory. log2FC>0 = up in shTACSTD2.",
)
print("QC perturbation:", summary["perturbation_qc"])
print("CLDN4:", summary["cldn4"])
print(setdf[["gene_set", "n_detected", "comp_direction", "comp_rank_biserial",
             "comp_p", "sc_mean_log2FC", "sc_n_up", "sc_n_down", "sc_wilcoxon_p"]].to_string(index=False))

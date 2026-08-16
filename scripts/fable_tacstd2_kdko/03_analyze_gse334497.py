#!/usr/bin/env python3
"""GSE334497 - mouse 4T1 TNBC, Trop2 KO vs WT tumours (RNA-seq, normalized counts).

Contrast: KO (5) vs WT (5). log2FC > 0 => up in Trop2 KO.
Stats: Welch t-test + BH FDR + Mann-Whitney on log2(normalized counts + 0.5).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis_common as AC

DATA = AC.DATA
df = pd.read_csv(DATA / "GSE334497_normalized_counts.csv.gz", index_col=0)
print("matrix:", df.shape, "cols:", list(df.columns))

# group by column name: contains 'KO' -> knockout, else control
case_idx = [i for i, c in enumerate(df.columns) if "KO" in c]
ctrl_idx = [i for i, c in enumerate(df.columns) if "KO" not in c]
print("KO cols  :", [df.columns[i] for i in case_idx])
print("WT cols  :", [df.columns[i] for i in ctrl_idx])
assert len(case_idx) == 5 and len(ctrl_idx) == 5

# map Ensembl mouse -> symbol
id2sym = AC.map_ensembl_to_symbol(list(df.index), "mouse", "GSE334497_id2symbol.json")
symbols = np.array([id2sym.get(str(e).split(".")[0], "") for e in df.index])
mapped = np.sum(symbols != "")
print(f"mapped {mapped}/{len(symbols)} Ensembl IDs to symbols")

mat = np.log2(df.values.astype(float) + 0.5)
# keep only mapped rows, then collapse to gene level
keep = symbols != ""
symbols, mat = symbols[keep], mat[keep]
symbols, mat, _ = AC.collapse_to_gene(symbols, mat)
print("gene-level matrix:", mat.shape)

summary, pergene, setdf = AC.run_expression_dataset(
    "GSE334497",
    "Mouse 4T1 TNBC tumours, Trop2 CRISPR KO vs WT (normalized counts)",
    symbols, mat, case_idx, ctrl_idx,
    value_kind="log2(normalized counts+0.5)",
    notes="5 KO vs 5 WT. log2FC>0 = up in Trop2 KO.",
)
print("QC perturbation:", summary["perturbation_qc"])
print("CLDN4:", summary["cldn4"])
print(setdf.to_string(index=False))

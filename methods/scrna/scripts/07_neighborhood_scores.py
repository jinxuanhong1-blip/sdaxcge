"""Immune neighborhood / niche scores: CD8, TLS, CXCL13 at the sample level.

NOTE: true neighborhood/TLS structure needs SPATIAL data. In dissociated scRNA we
APPROXIMATE via co-occurrence of B cells, CXCL13+ T cells, and CCL19/21 stroma per
sample, and via signature scores. State the approximation clearly.
注意：真正的邻域/TLS结构需空间数据；解离scRNA只能用共现与信号分数近似，须明确说明。

Usage: python 07_neighborhood_scores.py annotated.h5ad out_prefix
"""
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from gene_sets import CD8_CYTOTOXIC, CD8_EXHAUSTION, TLS_CORE, TLS_CHEMOKINE, CXCL13

def main(in_h5ad, out_prefix):
    adata = sc.read_h5ad(in_h5ad)
    use_raw = adata.raw is not None
    vnames = adata.raw.var_names if use_raw else adata.var_names

    for name, gs in {"CD8_cytotoxic": CD8_CYTOTOXIC, "CD8_exhaustion": CD8_EXHAUSTION,
                     "TLS_core": TLS_CORE, "TLS_chemokine": TLS_CHEMOKINE,
                     "CXCL13": CXCL13}.items():
        present = [g for g in gs if g in vnames]
        if present:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=use_raw)

    # sample-level summary. CD8 scores restricted to T/NK cells if lineage present.
    obs = adata.obs
    rows = []
    for s, d in obs.groupby("sample"):
        row = {"sample": s, "n_cells": len(d)}
        for m in ["CD8_cytotoxic", "CD8_exhaustion", "TLS_core", "TLS_chemokine", "CXCL13"]:
            if m in d:
                row[f"{m}_mean"] = float(d[m].mean())
        if "lineage" in d:
            row["frac_Bcell"] = float((d["lineage"] == "B/Plasma").mean())
            row["frac_TNK"] = float((d["lineage"] == "T/NK").mean())
        rows.append(row)
    summ = pd.DataFrame(rows).set_index("sample")
    summ.to_csv(f"{out_prefix}_neighborhood_by_sample.csv")
    print(summ)
    adata.write(f"{out_prefix}.h5ad")

if __name__ == "__main__":
    main(*sys.argv[1:])

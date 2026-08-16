#!/usr/bin/env python3
"""Peek at the structure of the two scRNA h5ad files without loading full X."""
import sys
from pathlib import Path
import anndata as ad
import pandas as pd

RAW = Path(__file__).resolve().parents[2] / "results" / "fable_irae" / "data" / "raw"
GENES = ["TACSTD2", "CLDN4"]


def peek(path):
    print("=" * 90)
    print(path.name)
    a = ad.read_h5ad(path, backed="r")
    print("shape (obs x var):", a.shape)
    print("X dtype:", a.X.dtype if hasattr(a.X, "dtype") else type(a.X))
    print("\nobs columns:")
    for c in a.obs.columns:
        vals = a.obs[c]
        nun = vals.nunique(dropna=True)
        sample = ", ".join(map(str, list(pd.unique(vals.dropna()))[:8]))
        print(f"  - {c} (nunique={nun}): {sample[:160]}")
    print("\nvar (n=%d) index name: %s" % (a.n_vars, a.var.index.name))
    print("var columns:", list(a.var.columns))
    print("first var names:", list(a.var_names[:8]))
    # Are our genes present? var_names may be symbols or ids
    vnames = set(a.var_names.astype(str))
    for g in GENES:
        present = g in vnames
        print(f"  gene {g} in var_names: {present}")
    # check any gene-symbol column
    for col in a.var.columns:
        colvals = set(a.var[col].astype(str))
        hits = [g for g in GENES if g in colvals]
        if hits:
            print(f"  genes found in var['{col}']: {hits}")
    # layers / raw
    print("layers:", list(a.layers.keys()) if a.layers else None)
    print("raw present:", a.raw is not None)
    if a.raw is not None:
        rn = set(map(str, a.raw.var_names))
        print("  raw n_vars:", a.raw.n_vars, "TACSTD2 in raw:", "TACSTD2" in rn, "CLDN4 in raw:", "CLDN4" in rn)
    print("obsm keys:", list(a.obsm.keys()))
    a.file.close()


def main():
    for name in ["GSE206300_ircolitis-tissue-epithelial.h5ad", "GSE277136_BLFplusNM.h5ad"]:
        try:
            peek(RAW / name)
        except Exception as e:
            print(f"ERROR peeking {name}: {e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()

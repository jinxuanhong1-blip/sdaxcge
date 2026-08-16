"""Claim A10 / P5+P6: cell-type-resolved public lung scRNA (CELLxGENE Census).

Two questions, answered separately:

1. Is the ELF3 module / TACSTD2 / CLDN4 program a *cell-type identity* (high in
   airway epithelium, low in AT2 where NKX2-1 lives)?  If yes, bulk correlations
   are composition.
2. *Within* a single epithelial type, do the genes still co-vary in the claimed
   directions?  That is the only scRNA evidence that would support a within-cell
   regulatory module rather than a lineage marker set.

We query only the claim genes + lineage markers. No genome-wide DE.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

HUMAN_GENES = C.CLAIM_GENES_H + list(C.MARKERS_H)
MOUSE_GENES = C.CLAIM_GENES_M + list(C.MARKERS_M)

HUMAN_EPITHELIAL = [
    "pulmonary alveolar type 2 cell",
    "pulmonary alveolar type 1 cell",
    "pulmonary alveolar epithelial cell",
    "club cell",
    "ciliated cell",
    "basal cell",
    "lung secretory cell",
    "epithelial cell of lung",
    "epithelial cell of lower respiratory tract",
    "epithelial cell",
    "multiciliated columnar cell of tracheobronchial tree",
    "goblet cell",
    "respiratory basal cell",
    "bronchial goblet cell",
    "ionocyte",
    "tuft cell",
    "malignant cell",
    "type II pneumocyte",
    "type I pneumocyte",
    "respiratory goblet cell",
    "nasal mucosa goblet cell",
    "tracheobronchial serous cell",
    "lung ciliated cell",
    "secretory cell",
]

# Max cells per (cell_type, disease) for within-type correlations.
MAX_PER_STRATUM = 8000
MIN_CELLS = 80
RNG = np.random.default_rng(C.RNG_SEED)


def _quote_list(xs):
    # tiledbsoma QueryCondition requires `attr in ['a', 'b']` (square brackets).
    return "[" + ", ".join(repr(x) for x in xs) + "]"


def fetch(organism: str, genes: list[str], cell_types: list[str] | None,
          diseases: list[str] | None):
    import cellxgene_census

    obs_filter = "tissue_general == 'lung' and is_primary_data == True"
    if cell_types:
        obs_filter += f" and cell_type in {_quote_list(cell_types)}"
    if diseases:
        obs_filter += f" and disease in {_quote_list(diseases)}"
    var_filter = f"feature_name in {_quote_list(genes)}"
    print(f"[{organism}] query: {obs_filter}", flush=True)
    with cellxgene_census.open_soma(census_version="stable") as census:
        adata = cellxgene_census.get_anndata(
            census,
            organism=organism,
            measurement_name="RNA",
            obs_value_filter=obs_filter,
            var_value_filter=var_filter,
            obs_column_names=["cell_type", "disease", "assay", "dataset_id",
                              "tissue", "suspension_type", "raw_sum"],
        )
    print(f"[{organism}] got {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
    return adata


def to_log1p_cpm(adata):
    """Library-size normalise from raw counts using Census `raw_sum`.

    `raw_sum` is the full-transcriptome UMI/count total. Using the 20-gene
    subset sum as a library size would inflate CPM in cells that express only
    the claim genes and produce NaNs in cells that express none of them.
    """
    X = adata.X
    if "raw_sum" in adata.obs.columns:
        lib = pd.to_numeric(adata.obs["raw_sum"], errors="coerce").to_numpy(dtype=float)
    else:
        lib = np.asarray(X.sum(axis=1)).ravel() if hasattr(X, "toarray") else X.sum(axis=1)
    lib = np.where(lib > 0, lib, np.nan)
    return X, lib, list(adata.var["feature_name"].astype(str))


def gene_col(X, lib, var_names, gene):
    if gene not in var_names:
        return None
    j = var_names.index(gene)
    col = X[:, j]
    if hasattr(col, "toarray"):
        col = np.asarray(col.toarray()).ravel()
    else:
        col = np.asarray(col).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.log1p(col / lib * 1e4)
    out = np.where(np.isfinite(out), out, 0.0)
    return out


def celltype_means(adata, genes, organism):
    X, lib, var_names = to_log1p_cpm(adata)
    obs = adata.obs.reset_index(drop=True)
    rows = []
    for (ct, dis), idx in obs.groupby(["cell_type", "disease"], observed=True).groups.items():
        idx = np.asarray(list(idx))
        if len(idx) < 20:
            continue
        rec = {"organism": organism, "cell_type": ct, "disease": dis, "n_cells": int(len(idx))}
        for g in genes:
            v = gene_col(X, lib, var_names, g)
            rec[g] = float(np.nanmean(v[idx])) if v is not None else np.nan
            rec[f"{g}_frac_pos"] = float((v[idx] > 0).mean()) if v is not None else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def within_type_corrs(adata, genes, organism, diseases_keep):
    X, lib, var_names = to_log1p_cpm(adata)
    obs = adata.obs.reset_index(drop=True)
    have = [g for g in genes if g in var_names]
    rows = []
    for (ct, dis), idx in obs.groupby(["cell_type", "disease"], observed=True).groups.items():
        if dis not in diseases_keep:
            continue
        idx = np.asarray(list(idx))
        if len(idx) < MIN_CELLS:
            continue
        if len(idx) > MAX_PER_STRATUM:
            idx = RNG.choice(idx, size=MAX_PER_STRATUM, replace=False)
        mat = {g: gene_col(X, lib, var_names, g)[idx] for g in have}
        # skip genes that are almost never detected in this type
        use = [g for g in have if np.nanmean(mat[g] > 0) >= 0.02]
        for i, ga in enumerate(use):
            for gb in use[i + 1:]:
                x, y = mat[ga], mat[gb]
                if np.nanstd(x) == 0 or np.nanstd(y) == 0:
                    continue
                rho, p = stats.spearmanr(x, y, nan_policy="omit")
                rows.append({"organism": organism, "cell_type": ct, "disease": dis,
                             "n_cells": int(len(idx)), "gene_a": ga, "gene_b": gb,
                             "spearman_rho": float(rho), "p_value": float(p),
                             "frac_pos_a": float((x > 0).mean()),
                             "frac_pos_b": float((y > 0).mean())})
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = C.bh_fdr(df["p_value"].to_numpy())
    return df


def main():
    # Human: normal + LUAD epithelium (the two contexts the claim is about)
    human_types = HUMAN_EPITHELIAL
    adata_h = fetch("homo_sapiens", HUMAN_GENES, human_types,
                    ["normal", "lung adenocarcinoma"])
    means_h = celltype_means(adata_h, HUMAN_GENES, "human")
    within_h = within_type_corrs(adata_h, C.CLAIM_GENES_H, "human",
                                 {"normal", "lung adenocarcinoma"})
    C.write_table(means_h, "census_human_celltype_means.csv")
    C.write_table(within_h, "census_human_within_type_spearman.csv")

    # Mouse: all lung primary cells (the atlas is small; airway types may be sparse)
    adata_m = fetch("mus_musculus", MOUSE_GENES, None, None)
    means_m = celltype_means(adata_m, MOUSE_GENES, "mouse")
    within_m = within_type_corrs(adata_m, C.CLAIM_GENES_M, "mouse",
                                 set(adata_m.obs["disease"].astype(str).unique()))
    C.write_table(means_m, "census_mouse_celltype_means.csv")
    C.write_table(within_m, "census_mouse_within_type_spearman.csv")
    print("done", flush=True)


if __name__ == "__main__":
    main()

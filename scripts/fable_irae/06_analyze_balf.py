#!/usr/bin/env python3
"""GSE277136 - bronchoalveolar lavage fluid (BALF) scRNA-seq, ICI-related pneumonitis.

Lung compartment. 4 ICI-pneumonitis patients ("AE", source='our data') integrated with
3 healthy controls ("HC", source='Nature Medicine' reference). The processed .h5ad stores
scaled HVGs in X, but .raw holds log1p-normalized expression for 22,126 genes including
TACSTD2, CLDN4 and EPCAM -> both markers are measurable here.

We (1) localize TACSTD2/CLDN4/EPCAM across BAL cell types (are epithelial cells present
and do they carry the markers?) and (2) compare AE vs HC by per-sample pseudobulk.

Caveat: AE vs HC is fully confounded with data source/batch (all AE from one lab, all HC
from a published reference), so the AE-vs-HC contrast is descriptive only.

Outputs:
  results/fable_irae/tables/balf_GSE277136_by_celltype.tsv
  results/fable_irae/tables/balf_GSE277136_pseudobulk_condition.tsv
  results/fable_irae/figures/balf_GSE277136_celltype.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import anndata as ad
from scipy import sparse, stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "results" / "fable_irae"
RAW = ROOT / "data" / "raw" / "GSE277136_BLFplusNM.h5ad"
TAB = ROOT / "tables"; TAB.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"; FIG.mkdir(parents=True, exist_ok=True)

GENES = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC"]


def main():
    a = ad.read_h5ad(RAW, backed="r")
    obs = a.obs[["sample", "source", "cell_label", "predicted_labels",
                 "majority_voting", "lung_altas_predict"]].copy()
    a.file.close()

    with h5py.File(RAW, "r") as f:
        var_idx = np.array([x.decode() if isinstance(x, bytes) else str(x)
                            for x in f["raw/var/_index"][:]])
        gpos = {g: int(np.where(var_idx == g)[0][0]) for g in GENES if g in set(var_idx)}
        shape = tuple(f["raw/X"].attrs["shape"])
        X = sparse.csr_matrix((f["raw/X/data"][:], f["raw/X/indices"][:],
                               f["raw/X/indptr"][:]), shape=shape)
    print("raw X loaded:", X.shape, "genes found:", gpos)
    Xc = X[:, [gpos[g] for g in gpos]].toarray()
    expr = pd.DataFrame(Xc, columns=list(gpos.keys()), index=obs.index)
    df = expr.join(obs)
    df["condition"] = np.where(df["source"].astype(str).str.contains("our data"), "AE_pneumonitis", "HC")

    # (1) per cell type: mean log1p expr + detection rate (use predicted_labels: has epithelial)
    ct_col = "predicted_labels"
    agg = {}
    for g in gpos:
        agg[f"{g}_mean"] = (g, "mean")
        agg[f"{g}_det"] = (g, lambda s: float((s > 0).mean()))
    byct = df.groupby(ct_col, observed=True).agg(n_cells=(list(gpos)[0], "size"), **agg)
    byct = byct.sort_values("n_cells", ascending=False)
    byct.to_csv(TAB / "balf_GSE277136_by_celltype.tsv", sep="\t")
    print("[write] BALF by cell type")
    epi_like = [c for c in byct.index if any(k in c for k in
                ["Club", "Multiciliated", "Suprabasal", "Basal", "Goblet", "AT1", "AT2",
                 "Ionocyte", "Deuterosomal", "Secretory", "epithel", "Ciliated"])]
    print("epithelial-like cell types detected:", epi_like)
    cols = ["n_cells"] + [c for c in byct.columns if c.startswith(("TACSTD2", "CLDN4", "EPCAM"))]
    print(byct[cols].head(20).to_string())

    # (2) per-sample pseudobulk, AE vs HC
    rows = []
    for samp, g in df.groupby("sample", observed=True):
        rec = {"sample": samp, "condition": g["condition"].iloc[0], "n_cells": len(g)}
        for gene in gpos:
            rec[f"{gene}_mean_log1p"] = g[gene].mean()
            rec[f"{gene}_detrate"] = float((g[gene] > 0).mean())
        rows.append(rec)
    pb = pd.DataFrame(rows)
    pb.to_csv(TAB / "balf_GSE277136_pseudobulk_condition.tsv", sep="\t", index=False)
    print("[write] BALF pseudobulk by sample")
    print(pb.to_string(index=False))

    de_rows = []
    for gene in gpos:
        ae = pb.loc[pb.condition == "AE_pneumonitis", f"{gene}_mean_log1p"].values
        hc = pb.loc[pb.condition == "HC", f"{gene}_mean_log1p"].values
        if len(ae) and len(hc):
            u, p = stats.mannwhitneyu(ae, hc, alternative="two-sided")
            de_rows.append(dict(gene=gene, mean_AE=ae.mean(), mean_HC=hc.mean(),
                                n_AE=len(ae), n_HC=len(hc), mwu_U=u, p_value=p))
    de = pd.DataFrame(de_rows)
    de.to_csv(TAB / "balf_GSE277136_de_condition.tsv", sep="\t", index=False)
    print("[write] BALF AE-vs-HC pseudobulk DE (batch-confounded, descriptive)")
    print(de.to_string(index=False))

    # figure: TACSTD2/CLDN4/EPCAM mean by cell type (top 15)
    top = byct.head(15)
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(top))
    w = 0.25
    for i, g in enumerate(["TACSTD2", "CLDN4", "EPCAM"]):
        if f"{g}_mean" in top:
            ax.bar(x + (i - 1) * w, top[f"{g}_mean"], w, label=g)
    ax.set_xticks(x)
    ax.set_xticklabels(top.index, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("mean log1p expression")
    ax.set_title("GSE277136 BALF (ICI pneumonitis) — epithelial markers by cell type")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "balf_GSE277136_celltype.png", dpi=140)
    print("[write] BALF figure")


if __name__ == "__main__":
    main()

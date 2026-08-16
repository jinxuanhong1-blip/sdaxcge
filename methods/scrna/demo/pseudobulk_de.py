#!/usr/bin/env python3
"""Demo pseudobulk DE on GSE207422 epithelial cells: MPR-like vs NMPR.

EXPLORATORY / DESCRIPTIVE ONLY. Aggregates raw counts per sample within
epithelial cells (from run_demo.py) and fits DESeq2 (pydeseq2) with
sample as the unit of replication. n is small; results are
hypothesis-generating, NOT confirmatory. No stats fabricated.

Grouping: pCR+MPR -> "MPR_like" vs NMPR. Pre-treatment / NE excluded
when they lack a usable pathologic-response call.

Run after run_demo.py:
    python pseudobulk_de.py
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> int:
    counts = pd.read_csv(
        os.path.join(HERE, "pseudobulk_epithelial_counts.csv"), index_col=0
    )
    design = pd.read_csv(
        os.path.join(HERE, "pseudobulk_design.csv"), index_col=0
    )

    resp_map = {"pCR": "MPR_like", "MPR": "MPR_like", "NMPR": "NMPR"}
    design["grp"] = design["Pathologic Response"].map(resp_map)
    post = design["Resource"].astype(str).str.contains("Post", na=False)
    keep = design.index[
        design["grp"].isin(["MPR_like", "NMPR"])
        & post
        & (design["n_cells_epithelial"] >= 20)
        & design.index.isin(counts.index)
    ]
    counts = counts.loc[keep]
    design = design.loc[keep]
    print(
        "Samples in DE:\n",
        design[
            ["Pathologic Response", "grp", "n_cells_epithelial", "Pathology"]
        ].to_string(),
    )
    print("\nGroup sizes:", design["grp"].value_counts().to_dict())

    if design["grp"].nunique() < 2 or design["grp"].value_counts().min() < 2:
        print("Not enough samples per group for a DE contrast; skipping.")
        return 0

    counts = counts.loc[:, counts.sum(axis=0) >= 10].astype(int)

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    meta = design[["grp"]].copy()
    meta["grp"] = pd.Categorical(meta["grp"], categories=["NMPR", "MPR_like"])
    dds = DeseqDataSet(
        counts=counts, metadata=meta, design="~grp", refit_cooks=True
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["grp", "MPR_like", "NMPR"])
    stat.summary()
    res = stat.results_df.sort_values("padj")
    res.to_csv(os.path.join(HERE, "de_epithelial_MPRlike_vs_NMPR.csv"))

    goi = [
        "TACSTD2", "CLDN4", "CLDN3", "CLDN7", "CLDN18", "CDH1", "EPCAM",
        "TJP1", "OCLN", "F11R", "ELF3", "CXCL13",
    ]
    sub = res[res.index.isin(goi)]
    print("\nGenes of interest (epithelial MPR_like vs NMPR):")
    print(sub.to_string())
    sub.to_csv(os.path.join(HERE, "de_epithelial_genes_of_interest.csv"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    top = [g for g in ["TACSTD2", "CLDN4"] if g in cpm.columns]
    figdir = os.path.join(HERE, "figures")
    os.makedirs(figdir, exist_ok=True)
    if top:
        rng = np.random.default_rng(0)
        fig, axes = plt.subplots(
            1, len(top), figsize=(4.5 * len(top), 4.2), squeeze=False
        )
        for ax, g in zip(axes[0], top):
            for i, grp in enumerate(["MPR_like", "NMPR"]):
                s = design.index[design["grp"] == grp]
                y = np.log1p(cpm.loc[s, g].values)
                ax.scatter(
                    rng.normal(i, 0.06, len(y)), y, s=40,
                    color="#c00000" if grp == "MPR_like" else "#4472c4",
                )
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["MPR_like", "NMPR"])
            ax.set_title(g)
            ax.set_ylabel("log1p CPM (pseudobulk)")
        fig.suptitle("Epithelial pseudobulk, per-sample points (exploratory)")
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "06_de_goi_points.png"), dpi=130)
        plt.close(fig)

    print("\nDE tables written. EXPLORATORY: small n, no external replication here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

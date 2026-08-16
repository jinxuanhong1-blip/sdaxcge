#!/usr/bin/env python3
"""Box + strip of patient-level detection fraction. Immune libraries only."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RESULTS


def _panel(ax, df, gene, title):
    mpr = df.loc[df["mpr_any"] == "yes", f"{gene}_frac_pos"].to_numpy()
    non = df.loc[df["mpr_any"] == "no", f"{gene}_frac_pos"].to_numpy()
    data = [non, mpr]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.45, patch_artist=True,
                    showfliers=False, medianprops={"color": "black"})
    colors = ["#C44E52", "#4C72B0"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data):
        x = rng.normal(i, 0.06, size=len(vals))
        ax.scatter(x, vals, s=14 if len(vals) > 20 else 36, c=colors[i],
                   edgecolors="white", linewidths=0.4, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"non-MPR\nn={len(non)}", f"MPR-any\nn={len(mpr)}"])
    ax.set_ylabel(f"{gene} fraction positive")
    ax.set_title(title, loc="left", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    g243 = pd.read_csv(RESULTS / "gse243013_patient_level.tsv", sep="\t")
    g229 = pd.read_csv(RESULTS / "gse229353_patient_level.tsv", sep="\t")
    g243 = g243[g243["mpr_any"].isin(["yes", "no"])]
    if "mtx_complete" in g229.columns:
        g229 = g229[g229["mtx_complete"].astype(str).isin(["True", "true", "1"])]

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.2))
    _panel(axes[0, 0], g243, "TACSTD2", "GSE243013 CD45+ immune atlas")
    _panel(axes[0, 1], g243, "CLDN4", "GSE243013 CD45+ immune atlas")
    _panel(axes[1, 0], g229, "TACSTD2", "GSE229353 CD45+ bead-sort (P03 MTX dropped)")
    _panel(axes[1, 1], g229, "CLDN4", "GSE229353 CD45+ bead-sort (P03 MTX dropped)")
    fig.suptitle(
        "Immune-library TACSTD2/CLDN4 detection vs MPR\n"
        "Not malignant RNA. Residual/ambient/epithelial leak in CD45+ libraries.",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_extra_immune_leak_vs_mpr.png", dpi=180)
    fig.savefig(RESULTS / "fig_extra_immune_leak_vs_mpr.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()

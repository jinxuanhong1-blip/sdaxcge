"""Figures for the replication / interaction follow-up."""

from __future__ import annotations

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from opus_tls_lib import RESULTS_DIR, fmt_p

FIG = os.path.join(RESULTS_DIR, "figures")
TAB = os.path.join(RESULTS_DIR, "tables")


def fig_interaction() -> None:
    d = pd.read_csv(os.path.join(TAB, "tcga_tacstd2_histology_interaction.csv"))
    d = d.set_index("signature").loc[
        ["TLS_Cabrita", "TLS_12chemokine", "TLS_imprint", "B_cell",
         "Plasma_cell", "Tfh", "T_cell_CD8", "IFNg_Ayers"]
    ]
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    y = np.arange(len(d))
    ax.scatter(d["beta_tac_luad"], y - 0.12, s=28, color="#888888", label="TACSTD2 slope in LUAD")
    ax.scatter(d["beta_tac_luad"] + d["beta_interaction"], y + 0.12, s=28,
               color="#b2182b", label="implied slope in LUSC")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(d.index)
    ax.set_xlabel("OLS β  (inverse-normal signature ~ TACSTD2 + LUSC + TACSTD2×LUSC + ABSOLUTE)")
    ax.set_title("TCGA n=1008  TACSTD2 × histology interaction\n"
                 "negative interaction = steeper anti-immune slope in LUSC")
    for i, sig in enumerate(d.index):
        ax.text(0.02, i, f"p_int={fmt_p(d.loc[sig, 'p_interaction'])}",
                transform=ax.get_yaxis_transform(), va="center", fontsize=7, color="#b2182b")
    ax.legend(frameon=False, loc="lower left")
    fig.savefig(os.path.join(FIG, "fig7_tcga_interaction.png"), bbox_inches="tight",
                facecolor="white")
    plt.close()


def fig_replication_grid() -> None:
    d = pd.read_csv(os.path.join(TAB, "replication_lusc_histology.csv"))
    d = d[d.gene == "TACSTD2"]
    sigs = ["TLS_Cabrita", "B_cell", "Plasma_cell", "T_cell_CD8"]
    cohs = ["GSE4573", "GSE17710", "GSE103584_SCC", "GSE50081_SCC",
            "GSE19188_SCC", "GSE103584_ADC", "GSE50081_ADC", "GSE19188_ADC"]
    color = {"HOLDS": "#2166ac", "NO_EVIDENCE": "#f0f0f0", "OPPOSITE": "#b2182b",
             "UNDERPOWERED": "#dddddd"}
    fig, ax = plt.subplots(figsize=(8.8, 3.2))
    ax.set_xlim(-0.5, len(cohs) - 0.5)
    ax.set_ylim(-0.5, len(sigs) - 0.5)
    ax.set_xticks(range(len(cohs)))
    ax.set_xticklabels(cohs, rotation=30, ha="right")
    ax.set_yticks(range(len(sigs)))
    ax.set_yticklabels(sigs)
    for i, sig in enumerate(sigs):
        for j, coh in enumerate(cohs):
            row = d[(d.signature == sig) & (d.cohort == coh)]
            if row.empty:
                continue
            v = row.iloc[0]["verdict"]
            ax.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                       facecolor=color.get(v, "#eee"), edgecolor="white"))
            ax.text(j, i, f"{row.iloc[0]['rho_adj_epithelial']:.2f}",
                    ha="center", va="center", fontsize=8,
                    color="white" if v in ("HOLDS", "OPPOSITE") else "black")
    ax.set_title("Independent public arrays/RNA  TACSTD2 after epithelial correction\n"
                 "blue=HOLDS  red=OPPOSITE  grey=no evidence / n<40")
    ax.invert_yaxis()
    fig.savefig(os.path.join(FIG, "fig8_replication_grid.png"), bbox_inches="tight",
                facecolor="white")
    plt.close()


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig_interaction()
    fig_replication_grid()
    print("wrote fig7 fig8")

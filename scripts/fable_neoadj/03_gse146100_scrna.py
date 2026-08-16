"""GSE146100 -- descriptive (single multiprimary-LUAD patient, 3 nodules).

scRNA-seq (normalized) of three synchronous lung adenocarcinoma nodules from
ONE patient treated with induction pembrolizumab (Zhang et al., PMID 33820821):
  W1 = non-responded (EGFR L858R)
  W2 = responded     (KRAS G12C)
  W3 = non-responded (EGFR L858R/R77H)

n=1 patient, so this is a purely descriptive per-nodule comparison of
epithelial (EPCAM+) TACSTD2 / CLDN4 -- not a statistical test.
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import GENES, DATA, FIG, TAB

FILE = os.path.join(DATA, "GSE146100_NormData.txt.gz")
NODULE_RESPONSE = {"W1": "non-responded", "W2": "responded", "W3": "non-responded"}
TARGETS = GENES + ["EPCAM"]


def main():
    with gzip.open(FILE, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        nodule = np.array([c.split("_")[0] for c in cells])
        want = set(TARGETS)
        mat = {}
        for line in fh:
            g = line.split("\t", 1)[0].strip('"')
            if g in want:
                mat[g] = np.array(line.rstrip("\n").split("\t")[1:], float)
                if len(mat) == len(want):
                    break

    epcam = mat["EPCAM"]
    rows = []
    for w in ["W1", "W2", "W3"]:
        cell_mask = nodule == w
        epi_mask = cell_mask & (epcam > 0)
        n_epi = int(epi_mask.sum())
        row = {"nodule": w, "response": NODULE_RESPONSE[w],
               "n_cells": int(cell_mask.sum()), "n_epithelial": n_epi}
        for g in GENES:
            v = mat[g][epi_mask]
            row[f"{g}_mean_epi"] = float(np.mean(v)) if n_epi else np.nan
            row[f"{g}_pct_pos_epi"] = float(100.0 * np.mean(v > 0)) if n_epi else np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB, "GSE146100_per_nodule.csv"), index=False)
    print(df.to_string())

    _bars(df)

    summary = {
        "dataset": "GSE146100_scRNA",
        "design": "1 patient, 3 synchronous LUAD nodules, induction pembrolizumab",
        "note": "descriptive only (n=1 patient); EPCAM+ epithelial cells",
    }
    with open(os.path.join(TAB, "GSE146100_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def _bars(df):
    fig, axes = plt.subplots(1, len(GENES), figsize=(4.2 * len(GENES), 4.0))
    colors = {"responded": "#d1495b", "non-responded": "#8aa1b1"}
    for ax, g in zip(axes, GENES):
        vals = df[f"{g}_mean_epi"].values
        cols = [colors[r] for r in df["response"]]
        ax.bar(df["nodule"], vals, color=cols)
        for i, (v, r) in enumerate(zip(vals, df["response"])):
            ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
        ax.set_ylabel(f"{g} mean (EPCAM+ cells)")
        ax.set_title(g)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors.values()]
    fig.legend(handles, list(colors.keys()), loc="upper right")
    fig.suptitle("GSE146100 epithelial ADC targets per nodule (1 patient)", y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "GSE146100_per_nodule.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

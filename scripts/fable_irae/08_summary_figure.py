#!/usr/bin/env python3
"""Cross-dataset summary figure: measurability (detection rate) of TACSTD2/CLDN4
across the three ICI compartments."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "results" / "fable_irae"
TAB = ROOT / "tables"; FIG = ROOT / "figures"

m = pd.read_csv(TAB / "master_summary.tsv", sep="\t")
labels = {
    "GSE206300": "Colon epithelium\n(colitis)",
    "GSE277136": "BALF / lung\n(pneumonitis)",
    "GSE319496": "Whole blood\n(irAE Y/N)",
}
order = ["GSE206300", "GSE277136", "GSE319496"]
genes = ["TACSTD2", "CLDN4"]

fig, ax = plt.subplots(figsize=(7.5, 4.5))
x = np.arange(len(order)); w = 0.38
for i, g in enumerate(genes):
    vals = []
    for ds in order:
        r = m[(m.dataset == ds) & (m.gene == g)]
        v = float(r["detection_rate"].values[0]) if len(r) else 0.0
        vals.append(v * 100)
    bars = ax.bar(x + (i - 0.5) * w, vals, w, label=g,
                  color=["#2c7fb8", "#c0392b"][i])
    for b, ds in zip(bars, order):
        r = m[(m.dataset == ds) & (m.gene == g)]
        meas = bool(r["measurable"].values[0]) if len(r) else False
        if not meas:
            ax.text(b.get_x() + b.get_width() / 2, 0.5, "absent\nfrom matrix",
                    ha="center", va="bottom", fontsize=7, rotation=0, color="grey")
ax.set_xticks(x); ax.set_xticklabels([labels[d] for d in order])
ax.set_ylabel("% cells / samples with detected expression")
ax.set_title("Measurability of TACSTD2 (TROP2) & CLDN4 across ICI irAE compartments")
ax.legend(title="gene")
fig.tight_layout()
fig.savefig(FIG / "summary_measurability.png", dpi=140)
print("[write] summary_measurability.png")

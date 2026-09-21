#!/usr/bin/env python3
"""Plot pair-level log2 fold-changes from the evidence table. No new statistics."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

src = Path("/workspace/results/tacstd2_cldn4_partial/evidence_table.tsv")
out = Path("/workspace/results/tacstd2_cldn4_partial/fig_pair_log2fc.png")
table = pd.read_csv(src, sep="\t", keep_default_na=False)
table = table[table["call"].isin(["supports", "null", "opposite", "perturbation_failed", "unreplicated"])].copy()
for col in ("tacstd2_log2fc", "cldn4_log2fc", "ifn_score_delta"):
    table[col] = pd.to_numeric(table[col], errors="coerce")
table = table.dropna(subset=["tacstd2_log2fc", "cldn4_log2fc", "ifn_score_delta"])
table["label"] = table["accession"] + "  " + table["contrast"]
table = table.iloc[::-1]

fig, ax = plt.subplots(figsize=(10.5, 8.2))
y = range(len(table))
ax.axvline(0, color="#888888", lw=0.8)
ax.scatter(table["tacstd2_log2fc"], y, s=36, c="#1f4e79", label="TACSTD2 / Tacstd2", zorder=3)
ax.scatter(table["cldn4_log2fc"], y, s=36, c="#c45911", label="CLDN4 / Cldn4", zorder=3)
ax.scatter(table["ifn_score_delta"], y, s=36, marker="D", c="#548235", label="IFN-score delta", zorder=3)
colors = {
    "supports": "#1f4e79",
    "null": "#666666",
    "opposite": "#c00000",
    "perturbation_failed": "#833c0c",
    "unreplicated": "#8064a2",
}
for i, call in enumerate(table["call"]):
    ax.plot(
        [table["tacstd2_log2fc"].iloc[i], table["cldn4_log2fc"].iloc[i], table["ifn_score_delta"].iloc[i]],
        [i, i, i],
        color=colors.get(call, "#999999"),
        lw=0.6,
        zorder=1,
    )
ax.set_yticks(list(y))
ax.set_yticklabels(table["label"], fontsize=7)
ax.set_xlabel("log2 difference (treated or perturbed minus control)")
ax.set_title("TACSTD2, CLDN4, and the IFN score in the same contrast")
ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(out, dpi=160)
fig.savefig("/opt/cursor/artifacts/tacstd2_cldn4_partial_log2fc.png", dpi=160)
print("wrote", out)

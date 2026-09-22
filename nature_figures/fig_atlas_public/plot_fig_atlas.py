#!/usr/bin/env python3
"""Public concordant-4 atlas: Harmony UMAP colored by cluster cell class.

Coordinates and labels are read only from data/umap_coordinates.tsv.gz,
which was extracted from the committed Seurat embeddings object
(methods/seurat_concordant4_cldn4/results/objects/seurat_harmony_embeddings.rds
on branch cursor/seurat-concordant4-cldn4-802a). Per-cell CLDN4 and TACSTD2
are not in that object, so this script does not draw feature panels.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
TABLE = HERE / "data" / "umap_coordinates.tsv.gz"
STEM = HERE / "fig_atlas_cellclass"

CLASS_ORDER = ["malignant", "T", "NK", "myeloid", "B", "other"]
DISPLAY = {
    "malignant": "Malignant",
    "T": "T",
    "NK": "NK",
    "myeloid": "Myeloid",
    "B": "B",
    "other": "Other",
}
# Distinct on white, readable in print. Other is light so it sits behind.
COLORS = {
    "malignant": "#C23B22",
    "T": "#2C5F9E",
    "NK": "#1B7F6E",
    "myeloid": "#C6860A",
    "B": "#6A4C9A",
    "other": "#C2C2C2",
}
# Draw abundant clouds first so smaller classes stay visible.
DRAW_ORDER = ["other", "myeloid", "malignant", "T", "B", "NK"]


def _use_arimo() -> None:
    """Arimo is the Arial-metric face installed here; Arial itself is not."""
    regular = None
    for path in fm.findSystemFonts():
        if path.endswith("Arimo-Regular.ttf") or path.endswith("Arimo[wght].ttf"):
            regular = path
            break
        if "Arimo" in path and "Italic" not in path and "Bold" not in path:
            regular = path
    if regular is None:
        mpl.rcParams["font.family"] = "Liberation Sans"
        return
    fm.fontManager.addfont(regular)
    name = fm.FontProperties(fname=regular).get_name()
    mpl.rcParams["font.family"] = name


def _density_peak(xy: np.ndarray, bins: int = 48) -> tuple[float, float]:
    hist, xe, ye = np.histogram2d(xy[:, 0], xy[:, 1], bins=bins)
    # Light box blur so a single crowded bin does not pin the label.
    kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], dtype=float)
    kernel /= kernel.sum()
    padded = np.pad(hist, 1, mode="constant")
    smooth = np.zeros_like(hist)
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            smooth[i, j] = np.sum(padded[i : i + 3, j : j + 3] * kernel)
    i, j = np.unravel_index(int(np.argmax(smooth)), smooth.shape)
    return float(0.5 * (xe[i] + xe[i + 1])), float(0.5 * (ye[j] + ye[j + 1]))


def main() -> None:
    df = pd.read_csv(TABLE, sep="\t")
    need = {"umap_1", "umap_2", "cell_class", "dataset", "unit_id"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"coordinate table missing columns: {sorted(missing)}")
    unknown = sorted(set(df["cell_class"]) - set(CLASS_ORDER))
    if unknown:
        raise SystemExit(f"unexpected cell_class values: {unknown}")

    _use_arimo()
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "font.size": 7,
        }
    )

    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(5.15, 4.85))
    # Point clouds rasterize; axis text, class names, and the legend stay vectors.
    ax.set_rasterization_zorder(1)

    for lab in DRAW_ORDER:
        sub = df.loc[df["cell_class"] == lab, ["umap_1", "umap_2"]].to_numpy(float)
        sub = sub[rng.permutation(len(sub))]
        ax.scatter(
            sub[:, 0],
            sub[:, 1],
            s=2.6,
            c=COLORS[lab],
            alpha=0.78,
            linewidths=0,
            zorder=0,
        )

    # Labels at the density peak of each class.
    for lab in CLASS_ORDER:
        xy = df.loc[df["cell_class"] == lab, ["umap_1", "umap_2"]].to_numpy(float)
        x, y = _density_peak(xy)
        color = "#333333" if lab == "other" else COLORS[lab]
        ax.text(
            x,
            y,
            DISPLAY[lab],
            color=color,
            fontsize=8,
            ha="center",
            va="center",
            zorder=6,
            bbox={
                "boxstyle": "round,pad=0.12",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.88,
            },
        )

    ax.set_xlabel("UMAP 1", fontsize=8, labelpad=2)
    ax.set_ylabel("UMAP 2", fontsize=8, labelpad=2)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.04)

    n_cells = len(df)
    n_units = int(df["unit_id"].nunique())
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COLORS[lab],
            markeredgecolor="none",
            markersize=5.5,
            label=f"{DISPLAY[lab]}  {int((df['cell_class'] == lab).sum()):,}",
        )
        for lab in CLASS_ORDER
    ]
    leg = ax.legend(
        handles=handles,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        fontsize=7,
        handletextpad=0.35,
        borderaxespad=0.0,
        labelspacing=0.45,
        title="Cell class",
    )
    leg.set_zorder(6)
    leg.get_title().set_fontsize(7)

    ax.text(
        0.0,
        1.02,
        "a",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
        zorder=6,
    )
    fig.text(
        0.01,
        0.012,
        f"{n_cells:,} cells plotted   ·   {n_units} units",
        fontsize=6.5,
        color="#333333",
        ha="left",
        va="bottom",
    )

    fig.subplots_adjust(left=0.06, right=0.78, bottom=0.08, top=0.94)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(STEM.with_suffix(f".{ext}"), dpi=400)
    plt.close(fig)

    # Echo the numbers the figure states, so a mismatch fails loudly.
    print(f"wrote {STEM}.*")
    print(f"cells={n_cells} units={n_units}")
    print(df.groupby("dataset")["unit_id"].nunique().to_string())
    print(df["cell_class"].value_counts().reindex(CLASS_ORDER).to_string())


if __name__ == "__main__":
    main()

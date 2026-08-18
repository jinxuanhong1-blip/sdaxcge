#!/usr/bin/env python3
"""Split the concordant-4 Harmony UMAP by unit CLDN4 quartile.

ADDITIVE visualization of the already-reported inverse association
(higher malignant CLDN4 ↔ fewer T/NK). Not a new test.

Datasets ONLY: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Same QC / cap≤350 / Harmony batch=dataset / annotations as the
atlas pipeline (PR #507). Quartiles = PR #503 locked unit table
(malignant CLDN4 %pos, within-cohort rank then Q1/Q4).

Does NOT re-audit PR #503 T/NK ρ (n=65, ρ=−0.531) or Q4 vs Q1 r.
CLDN4-only. No dual-high. No GSE148071.
"""
from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

import anndata as ad
import scanpy as sc

HERE = Path(__file__).resolve().parent
ATLAS_PY = HERE.parent / "concordant4_atlas_umap_annotate" / "analyze.py"
RESULTS = HERE / "results"
TABLES = RESULTS / "tables"
FIGS = RESULTS / "figures"
DEFAULT_DATA = Path("/tmp/geo_atlas")
LOCKED_UNITS = HERE / "data" / "tnk_units_pr503.tsv"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
TEAL = "#0f766e"
TEAL_FILL = "#14b8a6"
RED_HI = "#b91c1c"
GREY = "#c8c8c8"
MALIG_CMAP = LinearSegmentedColormap.from_list(
    "cldn4_red",
    ["#fde8e8", "#f5b7b1", "#e74c3c", "#922b21", "#641e16"],
)
COMP_COLORS = {
    "T": "#1d4e89",
    "NK": TEAL,
    "malignant": "#c0392b",
    "other": "#8a8a8a",
}
MM = 1.0 / 25.4

# Locked PR #503 numbers. Quoted, not re-tested.
PR503_N = 65
PR503_RHO = -0.531
PR503_Q1 = 19
PR503_Q4 = 16
PR503_R_RB = -0.724


def load_atlas_module():
    spec = importlib.util.spec_from_file_location("concordant4_atlas_analyze", ATLAS_PY)
    if spec is None or spec.loader is None:
        raise SystemExit(f"missing atlas pipeline: {ATLAS_PY}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def assign_quartiles(values: pd.Series) -> pd.Series:
    """PR #503 rule: rank (average ties) then qcut into Q1–Q4."""
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def nature_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def despine(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax, letter: str, x: float = -0.10, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
    )


def umap_limits(xy: np.ndarray, pad: float = 0.06) -> tuple[tuple[float, float], tuple[float, float]]:
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    dx = max(xmax - xmin, 1e-6)
    dy = max(ymax - ymin, 1e-6)
    return (xmin - pad * dx, xmax + pad * dx), (ymin - pad * dy, ymax + pad * dy)


def apply_umap_axes(ax, xlim, ylim) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    despine(ax)
    ax.set_aspect("equal", adjustable="box")


def load_locked_units(path: Path) -> pd.DataFrame:
    u = pd.read_csv(path, sep="\t")
    need = {"patient", "cohort", "quartile", "cldn4_pct", "frac_tnk"}
    missing = need - set(u.columns)
    if missing:
        raise SystemExit(f"{path} missing columns: {sorted(missing)}")
    u = u.rename(columns={"patient": "unit_id", "cohort": "dataset"})
    u["unit_id"] = u["unit_id"].astype(str)
    u["dataset"] = u["dataset"].astype(str)
    u["quartile"] = u["quartile"].astype(str)
    n_q1 = int((u["quartile"] == "Q1").sum())
    n_q4 = int((u["quartile"] == "Q4").sum())
    if n_q1 != PR503_Q1 or n_q4 != PR503_Q4:
        raise SystemExit(
            f"locked unit table Q1/Q4 = {n_q1}/{n_q4}, expected {PR503_Q1}/{PR503_Q4}"
        )
    if len(u) != PR503_N:
        raise SystemExit(f"locked unit table n={len(u)}, expected {PR503_N}")
    return u


def attach_quartiles(adata: ad.AnnData, locked: pd.DataFrame) -> pd.DataFrame:
    """Attach PR #503 quartiles to every cell of that unit. Also score atlas cells."""
    key = adata.obs["dataset"].astype(str) + "\t" + adata.obs["unit_id"].astype(str)
    lock_key = locked["dataset"].astype(str) + "\t" + locked["unit_id"].astype(str)
    qmap = locked.set_index(lock_key)["quartile"]
    pctmap = locked.set_index(lock_key)["cldn4_pct"]
    tnkmap = locked.set_index(lock_key)["frac_tnk"]
    adata.obs["cldn4_quartile"] = key.map(qmap)
    adata.obs["pr503_cldn4_pct"] = key.map(pctmap)
    adata.obs["pr503_frac_tnk"] = key.map(tnkmap)
    n_miss = int(adata.obs["cldn4_quartile"].isna().sum())
    if n_miss:
        miss = sorted(key[adata.obs["cldn4_quartile"].isna()].unique().tolist())
        raise SystemExit(f"{n_miss} cells have no PR #503 quartile; units={miss}")

    rows = []
    for (ds, uid), g in adata.obs.groupby(["dataset", "unit_id"], observed=True):
        lab = g["annotation"].astype(str)
        cldn4 = g["cldn4"].to_numpy(float)
        mal = lab.eq("malignant")
        n_mal = int(mal.sum())
        atlas_pct = float(100.0 * np.mean(cldn4[mal.to_numpy()] > 0)) if n_mal else np.nan
        n = int(len(g))
        rows.append(
            {
                "dataset": str(ds),
                "unit_id": str(uid),
                "pr503_quartile": str(g["cldn4_quartile"].iloc[0]),
                "pr503_cldn4_pct": float(g["pr503_cldn4_pct"].iloc[0]),
                "pr503_frac_tnk": float(g["pr503_frac_tnk"].iloc[0]),
                "n_cells_atlas": n,
                "n_malignant": n_mal,
                "n_T": int(lab.eq("T").sum()),
                "n_NK": int(lab.eq("NK").sum()),
                "n_other": int((~lab.isin(["T", "NK", "malignant"])).sum()),
                "atlas_cldn4_pct_pos": atlas_pct,
                "frac_T": float(lab.eq("T").mean()),
                "frac_NK": float(lab.eq("NK").mean()),
                "frac_malignant": float(lab.eq("malignant").mean()),
                "frac_other": float((~lab.isin(["T", "NK", "malignant"])).mean()),
                "frac_TNK": float(lab.isin(["T", "NK"]).mean()),
            }
        )
    units = pd.DataFrame(rows)
    units["atlas_quartile"] = ""
    for ds, idx in units.groupby("dataset").groups.items():
        sub = units.loc[idx]
        ok = sub["atlas_cldn4_pct_pos"].notna()
        if ok.sum() >= 4:
            qs = assign_quartiles(sub.loc[ok, "atlas_cldn4_pct_pos"])
            units.loc[qs.index, "atlas_quartile"] = qs.to_numpy()
    units["quartile_match"] = units["atlas_quartile"].eq(units["pr503_quartile"])
    return units


def composition_table(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for q in ("Q1", "Q4"):
        sub = units.loc[units["pr503_quartile"] == q]
        for col, name in (
            ("frac_T", "T"),
            ("frac_NK", "NK"),
            ("frac_malignant", "malignant"),
            ("frac_other", "other"),
        ):
            v = sub[col].to_numpy(float)
            rows.append(
                {
                    "quartile": q,
                    "compartment": name,
                    "n_units": int(len(sub)),
                    "mean_unit_frac": float(np.mean(v)),
                    "sd_unit_frac": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
                    "sem_unit_frac": float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
                    "median_unit_frac": float(np.median(v)),
                }
            )
    return pd.DataFrame(rows)


def _count_kde(xy: np.ndarray, grid_x, grid_y) -> np.ndarray:
    if len(xy) < 8:
        z = np.zeros(grid_x.shape, dtype=float)
        return z
    kde = gaussian_kde(xy.T, bw_method="scott")
    pos = np.vstack([grid_x.ravel(), grid_y.ravel()])
    return len(xy) * kde(pos).reshape(grid_x.shape)


def scatter_quartile_umap(ax, xy, anno, cldn4, vmax: float, s_other=1.0, s_mal=2.2, s_tnk=2.4):
    other = ~np.isin(anno, ["T", "NK", "malignant"])
    mal = anno == "malignant"
    tnk = np.isin(anno, ["T", "NK"])
    if other.any():
        ax.scatter(
            xy[other, 0],
            xy[other, 1],
            s=s_other,
            c=GREY,
            linewidths=0,
            alpha=0.45,
            rasterized=True,
            zorder=1,
        )
    sca = None
    if mal.any():
        sca = ax.scatter(
            xy[mal, 0],
            xy[mal, 1],
            s=s_mal,
            c=cldn4[mal],
            cmap=MALIG_CMAP,
            vmin=0,
            vmax=vmax,
            linewidths=0,
            alpha=0.92,
            rasterized=True,
            zorder=2,
        )
    if tnk.any():
        ax.scatter(
            xy[tnk, 0],
            xy[tnk, 1],
            s=s_tnk,
            c=TEAL,
            linewidths=0,
            alpha=0.80,
            rasterized=True,
            zorder=3,
        )
    return sca


def make_figure(
    adata: ad.AnnData,
    units: pd.DataFrame,
    path: Path,
) -> dict:
    nature_style()
    xy_all = np.asarray(adata.obsm["X_umap"])
    xlim, ylim = umap_limits(xy_all)
    obs = adata.obs
    anno = obs["annotation"].astype(str).to_numpy()
    quart = obs["cldn4_quartile"].astype(str).to_numpy()
    cldn4 = obs["cldn4"].to_numpy(float)
    mal = anno == "malignant"
    q14_mal = mal & np.isin(quart, ["Q1", "Q4"])
    pos = cldn4[q14_mal] if q14_mal.any() else cldn4[mal]
    pos = pos[pos > 0] if (pos > 0).any() else np.array([1.0])
    vmax = float(max(np.quantile(pos, 0.95), 0.5))

    n_q1_u = int((units["pr503_quartile"] == "Q1").sum())
    n_q4_u = int((units["pr503_quartile"] == "Q4").sum())
    n_q1_c = int((quart == "Q1").sum())
    n_q4_c = int((quart == "Q4").sum())
    n_q1_tnk = int(((quart == "Q1") & np.isin(anno, ["T", "NK"])).sum())
    n_q4_tnk = int(((quart == "Q4") & np.isin(anno, ["T", "NK"])).sum())

    fig = plt.figure(figsize=(180 * MM, 148 * MM), facecolor="white")
    gs = fig.add_gridspec(
        2,
        2,
        left=0.055,
        right=0.98,
        top=0.90,
        bottom=0.16,
        wspace=0.22,
        hspace=0.38,
        height_ratios=[1.12, 1.0],
    )

    # a / b — Q1 vs Q4 UMAP, shared xy
    handles_info = []
    for i, (letter, q, n_u, n_c) in enumerate(
        (
            ("a", "Q1", n_q1_u, n_q1_c),
            ("b", "Q4", n_q4_u, n_q4_c),
        )
    ):
        ax = fig.add_subplot(gs[0, i])
        panel_label(ax, letter, x=-0.08, y=1.10)
        m = quart == q
        sca = scatter_quartile_umap(ax, xy_all[m], anno[m], cldn4[m], vmax)
        apply_umap_axes(ax, xlim, ylim)
        ax.set_title(f"CLDN4-{q} units", pad=3)
        ax.text(
            0.0,
            1.015,
            f"n_units={n_u}  ·  n_cells={n_c:,}",
            transform=ax.transAxes,
            fontsize=6,
            va="bottom",
            ha="left",
            color="0.25",
        )
        if i == 1 and sca is not None:
            cb = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02, shrink=0.78)
            cb.set_label("malignant CLDN4", fontsize=5.5)
            cb.ax.tick_params(labelsize=5)
        handles_info.append(sca)

    # shared legend under a
    ax_a = fig.axes[0]
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=TEAL, markersize=4.5, label="T + NK"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#e74c3c",
            markersize=4.5,
            label="malignant (CLDN4)",
        ),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY, markersize=4.5, label="other"),
    ]
    ax_a.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        fontsize=5.5,
        handletextpad=0.25,
        borderpad=0.1,
        labelspacing=0.12,
        markerscale=1.0,
    )

    # c — T/NK count-density on the same UMAP (shared scale)
    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "c", x=-0.08, y=1.10)
    gx = np.linspace(xlim[0], xlim[1], 120)
    gy = np.linspace(ylim[0], ylim[1], 120)
    xx, yy = np.meshgrid(gx, gy)
    tnk = np.isin(anno, ["T", "NK"])
    z1 = _count_kde(xy_all[(quart == "Q1") & tnk], xx, yy)
    z4 = _count_kde(xy_all[(quart == "Q4") & tnk], xx, yy)
    zmax = float(max(z1.max() if z1.size else 0.0, z4.max() if z4.size else 0.0, 1e-9))
    levels = np.linspace(0.08 * zmax, zmax, 6)
    ax.scatter(
        xy_all[:, 0],
        xy_all[:, 1],
        s=0.4,
        c="#ececec",
        linewidths=0,
        alpha=0.35,
        rasterized=True,
        zorder=0,
    )
    ax.contourf(xx, yy, z1, levels=levels, cmap="YlGnBu", alpha=0.45, extend="max", zorder=1)
    ax.contour(xx, yy, z1, levels=levels, colors=TEAL, linewidths=0.55, zorder=2)
    ax.contourf(
        xx,
        yy,
        z4,
        levels=levels,
        colors=["#fde0dc", "#f4a582", "#d6604d", "#b2182b"],
        alpha=0.28,
        extend="max",
        zorder=3,
    )
    ax.contour(xx, yy, z4, levels=levels, colors=RED_HI, linewidths=0.7, linestyles="-", zorder=4)
    apply_umap_axes(ax, xlim, ylim)
    ax.set_title("T/NK density  ·  same UMAP, count scale", pad=3)
    ax.text(
        0.0,
        1.015,
        f"Q1 T/NK cells={n_q1_tnk:,}   Q4 T/NK cells={n_q4_tnk:,}",
        transform=ax.transAxes,
        fontsize=6,
        va="bottom",
        color="0.25",
    )
    ax.legend(
        handles=[
            mpatches.Patch(facecolor=TEAL_FILL, edgecolor=TEAL, alpha=0.7, label="Q1 T/NK"),
            mpatches.Patch(facecolor="#d6604d", edgecolor=RED_HI, alpha=0.7, label="Q4 T/NK"),
        ],
        loc="upper right",
        frameon=False,
        fontsize=5.5,
        handletextpad=0.35,
    )

    # d — unit-averaged composition
    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "d", x=-0.10, y=1.10)
    comps = ["T", "NK", "malignant", "other"]
    cols = ["frac_T", "frac_NK", "frac_malignant", "frac_other"]
    q1 = units.loc[units["pr503_quartile"] == "Q1"]
    q4 = units.loc[units["pr503_quartile"] == "Q4"]
    x = np.arange(len(comps))
    w = 0.36
    means1 = [float(q1[c].mean()) for c in cols]
    means4 = [float(q4[c].mean()) for c in cols]
    sem1 = [float(q1[c].sem()) for c in cols]
    sem4 = [float(q4[c].sem()) for c in cols]
    ax.bar(
        x - w / 2,
        means1,
        width=w,
        yerr=sem1,
        color="#5b9aa0",
        edgecolor="none",
        error_kw={"elinewidth": 0.6, "ecolor": "0.25", "capsize": 1.4},
        label=f"Q1  n_units={n_q1_u}",
        zorder=2,
    )
    ax.bar(
        x + w / 2,
        means4,
        width=w,
        yerr=sem4,
        color="#c0392b",
        edgecolor="none",
        error_kw={"elinewidth": 0.6, "ecolor": "0.25", "capsize": 1.4},
        label=f"Q4  n_units={n_q4_u}",
        zorder=2,
    )
    rng = np.random.default_rng(1)
    for i, c in enumerate(cols):
        for sub, shift, col in ((q1, -w / 2, "#1f4e5f"), (q4, w / 2, "#7b241c")):
            yy = sub[c].to_numpy(float)
            xx = np.full(len(yy), x[i] + shift) + rng.uniform(-0.05, 0.05, size=len(yy))
            ax.scatter(xx, yy, s=7, c=col, alpha=0.55, linewidths=0, zorder=3, rasterized=True)
    ax.set_xticks(x)
    ax.set_xticklabels(comps)
    ax.set_ylabel("unit fraction")
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0, decimals=0))
    despine(ax)
    ax.set_title("composition  ·  unit-averaged", pad=3)
    ax.legend(frameon=False, loc="upper right", fontsize=5.5)
    ax.text(
        0.0,
        1.015,
        "points = units; bars = mean ± SEM across units",
        transform=ax.transAxes,
        fontsize=5.5,
        va="bottom",
        color="0.35",
    )

    # e — schematic, no invented numbers
    e_ax = fig.add_axes([0.20, 0.035, 0.60, 0.07])
    e_ax.set_xlim(0, 1)
    e_ax.set_ylim(0, 1)
    e_ax.axis("off")
    e_ax.text(0.02, 0.55, "e", fontsize=8, fontweight="bold", va="center")
    e_ax.annotate(
        "",
        xy=(0.78, 0.50),
        xytext=(0.22, 0.50),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.1),
    )
    e_ax.text(0.22, 0.82, "CLDN4-high units", fontsize=6.5, ha="left", va="center", color=RED_HI)
    e_ax.text(0.78, 0.82, "T/NK down", fontsize=6.5, ha="right", va="center", color=TEAL)
    e_ax.text(
        0.50,
        0.12,
        "schematic of the already-reported inverse association  ·  not a new test",
        fontsize=5.5,
        ha="center",
        va="center",
        color="0.4",
    )

    fig.savefig(path.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), dpi=600, facecolor="white")
    fig.savefig(path.with_suffix(".svg"), dpi=600, facecolor="white")
    plt.close(fig)
    print(f"wrote {path}.png/pdf/svg", flush=True)
    return {
        "n_units_q1": n_q1_u,
        "n_units_q4": n_q4_u,
        "n_cells_q1": n_q1_c,
        "n_cells_q4": n_q4_c,
        "n_cells_tnk_q1": n_q1_tnk,
        "n_cells_tnk_q4": n_q4_tnk,
        "cldn4_vmax": vmax,
    }


def write_finding(units: pd.DataFrame, comp: pd.DataFrame, fig_stats: dict, path: Path) -> None:
    n_q1 = fig_stats["n_units_q1"]
    n_q4 = fig_stats["n_units_q4"]
    n_c1 = fig_stats["n_cells_q1"]
    n_c4 = fig_stats["n_cells_q4"]
    n_match = int(units["quartile_match"].sum())
    n_atlas_q = int((units["atlas_quartile"] != "").sum())
    by_ds = (
        units.groupby(["dataset", "pr503_quartile"])
        .size()
        .unstack(fill_value=0)
        .reindex(COHORTS)
        .reindex(columns=["Q1", "Q2", "Q3", "Q4"], fill_value=0)
    )

    def _m(q, name):
        r = comp[(comp["quartile"] == q) & (comp["compartment"] == name)]
        return float(r["mean_unit_frac"].iloc[0]) if len(r) else float("nan")

    lines = [
        "# Concordant-4 UMAP split by unit CLDN4 quartile",
        "",
        "ADDITIVE visualization only. **CLDN4-only.**",
        "",
        "This figure shows the already-reported inverse association",
        "(higher malignant CLDN4 ↔ fewer T/NK) on the cell-level Harmony UMAP.",
        "A single mixed UMAP cannot show it: T/NK and malignant occupy different",
        "clusters. The same embedding is therefore split by the **unit's** CLDN4",
        "quartile.",
        "",
        "**This is not a new test.** PR #503 T/NK numbers were not re-audited:",
        f"n={PR503_N}, %pos ρ={PR503_RHO}, stacked Q4 vs Q1 n={PR503_Q1}/{PR503_Q4}",
        f"r={PR503_R_RB}. No new Spearman is quoted.",
        "",
        "Datasets ONLY: **GSE123902 + GSE131907 + GSE205335 + GSE189357**.",
        "Not GSE148071 / GSE127465 / GSE154826 / GSE207422. No dual-high.",
        "",
        "## Honest n",
        "",
        f"- **n_units Q1 = {n_q1}**, **n_units Q4 = {n_q4}**",
        "  (PR #503 locked labels: within-cohort malignant CLDN4 %pos rank then qcut).",
        f"- **n_cells shown Q1 = {n_c1}**, **n_cells shown Q4 = {n_c4}**",
        "  (atlas cells after QC + cap ≤350/unit, from those units).",
        f"- T/NK cells in the density panel: Q1={fig_stats['n_cells_tnk_q1']},",
        f"  Q4={fig_stats['n_cells_tnk_q4']} (descriptive cell counts, not the inferential n).",
        f"- Inferential n remains the PR #503 units (n={PR503_N}; tails {PR503_Q1}/{PR503_Q4}).",
        "",
        "| dataset | Q1 | Q2 | Q3 | Q4 |",
        "|---|---:|---:|---:|---:|",
    ]
    for ds in COHORTS:
        r = by_ds.loc[ds]
        lines.append(f"| {ds} | {int(r['Q1'])} | {int(r['Q2'])} | {int(r['Q3'])} | {int(r['Q4'])} |")
    lines += [
        "",
        "## Quartile rule",
        "",
        "Unit-level malignant CLDN4 **%pos**, ranked **within each cohort**, then",
        "`qcut` on average-tie ranks → Q1/Q2/Q3/Q4. Same function as PR #503",
        "`assign_quartiles`. The figure uses the **locked PR #503 unit table**",
        "(`data/tnk_units_pr503.tsv`) so Q1/Q4 membership is the reported 19/16,",
        "not a re-cut on the 350/unit atlas subsample.",
        "",
        f"Atlas-capped malignant CLDN4 %pos, run through the same within-cohort",
        f"rank-then-qcut, matched the locked label for **{n_match}/{n_atlas_q}** units",
        "that had ≥1 malignant cell. Mismatches are expected: the atlas %pos is a",
        "capped subsample. The figure follows the locked labels.",
        "",
        "Every cell of a unit inherits that unit's quartile.",
        "",
        "## What the panels show",
        "",
        "- **a** UMAP of cells from CLDN4-Q1 units only. T+NK teal; malignant by",
        "  cell-level CLDN4 (red scale); other grey. Shared xy limits with b.",
        "- **b** Same for CLDN4-Q4 units. Same colour scale. Q4 is the hotter",
        "  malignant CLDN4 / thinner T/NK cloud.",
        "- **c** T/NK count-density (KDE × n) of Q1 vs Q4 on the same UMAP",
        "  coordinates and the same density levels — the “fewer T/NK” picture.",
        "- **d** Unit-averaged compartment fractions (T, NK, malignant, other) in",
        f"  Q1 vs Q4. Points are units. n_units = {n_q1} / {n_q4}. Not cell-pooled.",
        "- **e** Schematic only: CLDN4-high units → T/NK down. No new number.",
        "",
        "Unit-averaged fractions (descriptive; not a new test):",
        "",
        "| compartment | Q1 mean | Q4 mean |",
        "|---|---:|---:|",
        f"| T | {_m('Q1','T'):.3f} | {_m('Q4','T'):.3f} |",
        f"| NK | {_m('Q1','NK'):.3f} | {_m('Q4','NK'):.3f} |",
        f"| malignant | {_m('Q1','malignant'):.3f} | {_m('Q4','malignant'):.3f} |",
        f"| other | {_m('Q1','other'):.3f} | {_m('Q4','other'):.3f} |",
        "",
        "## Pipeline (reused, not re-audited)",
        "",
        "Same as `methods/concordant4_atlas_umap_annotate` (PR #507):",
        "",
        "- QC: n_genes ≥ 200, n_counts ≥ 500, mitochondrial % < 20.",
        "- Cap ≤350 cells / unit.",
        "- Inner-join genes → HVG 2000 → PCA 30 → Harmony `batch=dataset`",
        "  (theta=2.0). Sample is not a second Harmony key.",
        "- Neighbors k=15 on `X_pca_harmony`; one UMAP; Leiden 0.6.",
        "- Annotation: malignant / T / NK / myeloid / B / other",
        "  (author malignant where present, else marker scores).",
        "",
        "## What this is not",
        "",
        "- Not a re-audit of PR #503 T/NK ρ or IFN/MHC DE.",
        "- Not a new Spearman, MWU, or Q4 vs Q1 test.",
        "- Not a dual-high TACSTD2∩CLDN4 object.",
        "- Not GSE148071 / GSE127465 / GSE154826 / GSE207422.",
        "- Not evidence that CLDN4 *causes* T/NK exclusion.",
        "- Cell counts are not the inferential n.",
        "",
        "## Files",
        "",
        "- `results/figures/fig_umap_cldn4_tnk_split.png` / `.svg` / `.pdf`",
        "- `results/tables/unit_quartiles.tsv` — locked Q + atlas composition",
        "- `results/tables/composition_unit_avg.tsv` — panel d means",
        "- `data/tnk_units_pr503.tsv` — locked PR #503 unit table",
        "",
        "Reproduce:",
        "",
        "```bash",
        "python3 methods/concordant4_atlas_umap_annotate/download.py",
        "python3 methods/concordant4_umap_cldn4_tnk_split/analyze.py",
        "```",
        "",
    ]
    path.write_text("\n".join(lines))
    print(f"wrote {path}", flush=True)


def run_atlas(atlas, data: Path, cap: int) -> tuple[ad.AnnData, list[dict]]:
    rng = np.random.default_rng(atlas.SEED)
    sc.settings.verbosity = 2
    loaders = {
        "GSE123902": atlas.load_gse123902,
        "GSE131907": atlas.load_gse131907,
        "GSE205335": atlas.load_gse205335,
        "GSE189357": atlas.load_gse189357,
    }
    cache = data / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    adatas = []
    infos = []
    for name, fn in loaders.items():
        h5 = cache / f"{name}.h5ad"
        js = cache / f"{name}.json"
        if h5.exists() and js.exists():
            print(f"==== load {name} (cache) ====", flush=True)
            a = ad.read_h5ad(h5)
            info = json.loads(js.read_text())
        else:
            print(f"==== load {name} ====", flush=True)
            a, info = fn(data, cap, rng)
            a.write_h5ad(h5)
            js.write_text(json.dumps(info, indent=2, default=str))
        adatas.append(a)
        infos.append(info)
        gc.collect()
    print("==== concat / Harmony / Leiden ====", flush=True)
    adata = atlas.align_and_concat(adatas)
    del adatas
    gc.collect()
    adata = atlas.integrate(adata)
    print("==== annotate ====", flush=True)
    ann = atlas.annotate(adata)
    TABLES.joinpath("atlas_annotation.tsv").write_text(ann.to_csv(sep="\t", index=False))
    return adata, infos


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=DEFAULT_DATA)
    p.add_argument("--cap", type=int, default=350)
    p.add_argument("--h5ad", type=Path, default=None, help="skip rebuild if this integrated h5ad exists")
    args = p.parse_args()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    locked = load_locked_units(LOCKED_UNITS)
    atlas = load_atlas_module()
    integ = args.h5ad or (args.data / "cache" / "concordant4_integrated.h5ad")

    if integ.exists():
        print(f"==== load integrated {integ} ====", flush=True)
        adata = ad.read_h5ad(integ)
        infos = []
    else:
        adata, infos = run_atlas(atlas, args.data, args.cap)
        integ.parent.mkdir(parents=True, exist_ok=True)
        print(f"write {integ}", flush=True)
        adata.write_h5ad(integ)

    missing_obs = {"annotation", "cldn4", "dataset", "unit_id"} - set(adata.obs.columns)
    if missing_obs:
        raise SystemExit(f"integrated object missing obs: {sorted(missing_obs)}")
    if "X_umap" not in adata.obsm:
        raise SystemExit("integrated object has no X_umap")

    print("==== attach PR #503 quartiles ====", flush=True)
    units = attach_quartiles(adata, locked)
    units.to_csv(TABLES / "unit_quartiles.tsv", sep="\t", index=False)
    comp = composition_table(units)
    comp.to_csv(TABLES / "composition_unit_avg.tsv", sep="\t", index=False)

    print("==== figure ====", flush=True)
    fig_stats = make_figure(adata, units, FIGS / "fig_umap_cldn4_tnk_split")

    summary = {
        "additive_visualization_only": True,
        "new_spearman": False,
        "pr503_quoted": {
            "n": PR503_N,
            "rho_pct": PR503_RHO,
            "n_q1": PR503_Q1,
            "n_q4": PR503_Q4,
            "r_rb_q4q1": PR503_R_RB,
        },
        "n_units_q1": fig_stats["n_units_q1"],
        "n_units_q4": fig_stats["n_units_q4"],
        "n_cells_q1": fig_stats["n_cells_q1"],
        "n_cells_q4": fig_stats["n_cells_q4"],
        "n_cells_tnk_q1": fig_stats["n_cells_tnk_q1"],
        "n_cells_tnk_q4": fig_stats["n_cells_tnk_q4"],
        "n_units_total": int(len(units)),
        "n_cells_atlas": int(adata.n_obs),
        "quartile_source": "PR503_locked_tnk_units",
        "atlas_quartile_n_match": int(units["quartile_match"].sum()),
        "dual_high": False,
        "datasets": COHORTS,
        "gse148071": False,
        "composition_unit_avg": comp.to_dict(orient="records"),
        "infos": infos,
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_finding(units, comp, fig_stats, HERE / "FINDING.md")
    print(
        json.dumps(
            {k: summary[k] for k in ("n_units_q1", "n_units_q4", "n_cells_q1", "n_cells_q4")},
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

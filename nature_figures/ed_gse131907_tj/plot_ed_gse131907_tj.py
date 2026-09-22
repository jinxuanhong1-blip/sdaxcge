#!/usr/bin/env python3
"""Extended Data figure: GSE131907 patient-level TJ increase vs null TACSTD2–T.

Numbers are read from source tables copied from PR #773. The script recomputes
patient-level means, Wilcoxon tests and Spearman correlations and stops if they
disagree with the published values.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from scipy import stats

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source"

# Arimo is the metric-compatible Arial face shipped with the environment.
for _face in (
    "/usr/share/fonts/truetype/croscore/Arimo-Regular.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Bold.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-Italic.ttf",
    "/usr/share/fonts/truetype/croscore/Arimo-BoldItalic.ttf",
):
    font_manager.fontManager.addfont(_face)

NAVY = "#0072B2"
TEAL = "#009E73"
VERM = "#D55E00"
GREY = "#4D4D4D"
INK = "#1A1A1A"
ORANGE = "#E69F00"
HAIR = "#B5B5B5"

MINUS = "\u2212"
TIMES = "\u00d7"
SUP = str.maketrans("-0123456789", "\u207b⁰¹²³⁴⁵⁶⁷⁸⁹")

plt.rcParams.update(
    {
        "font.family": "Arimo",
        "font.size": 7,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.edgecolor": INK,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 400,
    }
)


def nfmt(x: float, nd: int = 3) -> str:
    return f"{x:.{nd}f}".replace("-", MINUS)


def pfmt(p: float) -> str:
    if p < 1e-3:
        exp = int(np.floor(np.log10(p)))
        mant = p / 10**exp
        return f"{mant:.1f} {TIMES} 10{str(exp).translate(SUP)}"
    if p < 0.01:
        return f"{p:.5f}"
    if p < 0.1:
        return f"{p:.3f}"
    return f"{p:.2f}"


def load_constants() -> dict[str, float]:
    tab = pd.read_csv(SRC / "published_constants.tsv", sep="\t")
    return {str(r.key): float(r.value) for r in tab.itertuples(index=False)}


def partial_spearman(x, y, z) -> tuple[float, float, int]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])

    def resid(a, b):
        design = np.column_stack([np.ones(len(b)), b])
        coef, _, _, _ = np.linalg.lstsq(design, a, rcond=None)
        return a - design @ coef

    xr, yr = resid(rx, rz), resid(ry, rz)
    res = stats.pearsonr(xr, yr)
    return float(res.statistic), float(res.pvalue), n


def mean_ci(x: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(x, float)
    n = x.size
    mean = float(x.mean())
    half = float(stats.t.ppf(0.975, n - 1) * stats.sem(x))
    return mean, mean - half, mean + half


def spearman_boot(x, y, rng: np.random.Generator, n_boot: int = 5000):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    res = stats.spearmanr(x, y)
    rho = float(res.statistic)
    p = float(res.pvalue)
    n = x.size
    draw = rng.integers(0, n, size=(n_boot, n))
    boots = np.empty(n_boot)
    for i in range(n_boot):
        boots[i] = stats.spearmanr(x[draw[i]], y[draw[i]]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return rho, p, n, float(lo), float(hi)


def style_data_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=6.5)


def style_forest(ax):
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(INK)
    ax.tick_params(axis="y", length=0, labelsize=6.5)
    ax.tick_params(axis="x", labelsize=6.5)


def prepare() -> dict:
    const = load_constants()
    contrasts = pd.read_csv(SRC / "tj_contrasts_for_figure.tsv", sep="\t")
    deltas = pd.read_csv(SRC / "tj_patient_deltas.tsv", sep="\t")
    author = pd.read_csv(SRC / "t_prespec_author_pct_T_patients.tsv", sep="\t")
    ts = pd.read_csv(SRC / "t_prespec_ts_pct_T_patients.tsv", sep="\t")
    searched = pd.read_csv(SRC / "t_winner_patients.tsv", sep="\t")
    prespec = pd.read_csv(SRC / "t_prespecified.tsv", sep="\t")

    winner = deltas.loc[deltas.tag == "winner"].copy()
    if not (winner.n_high == winner.n_low).all():
        raise SystemExit("Q4 and Q1 cell counts differ; the arm label would be wrong")
    w_delta = winner.delta.to_numpy(float)
    w_resid = winner.delta_resid_epi.to_numpy(float)
    w_mean, w_lo, w_hi = mean_ci(w_delta)
    r_mean, r_lo, r_hi = mean_ci(w_resid)
    w_p = float(stats.wilcoxon(w_delta, alternative="two-sided", method="exact").pvalue)
    r_p = float(stats.wilcoxon(w_resid, alternative="two-sided", method="exact").pvalue)

    if abs(w_mean - 0.5419839859008789) > 1e-12:
        raise SystemExit(f"claudin mean {w_mean} != published 0.542")
    if len(winner) != 10 or not np.all(w_delta > 0) or not np.all(w_resid > 0):
        raise SystemExit("winner panel is not 10/10 positive")
    if abs(w_p - 0.001953125) > 1e-15 or abs(r_p - 0.001953125) > 1e-15:
        raise SystemExit(f"Wilcoxon mismatch raw {w_p} residual {r_p}")
    if f"{w_mean:.3f}" != "0.542":
        raise SystemExit("display rounding for the claudin mean drifted")

    prespec_core = deltas.loc[deltas.tag == "prespec_author_median_pre_core"]
    core_delta = prespec_core.delta.to_numpy(float)
    c_mean, c_lo, c_hi = mean_ci(core_delta)
    c_p = float(stats.wilcoxon(core_delta, alternative="two-sided", method="exact").pvalue)
    if abs(c_mean - 0.13130335062742232) > 1e-12:
        raise SystemExit("pre-specified core mean drifted")
    if abs(c_p - 9.5367431640625e-06) > 1e-18:
        raise SystemExit(f"pre-specified Wilcoxon {c_p}")
    if int((core_delta > 0).sum()) != 18 or len(prespec_core) != 20:
        raise SystemExit("expected 18/20 positive pre-specified differences")

    # Published module rows that have no patient vector must match the file exactly.
    by_id = {r.id: r for r in contrasts.itertuples(index=False)}
    checks = {
        "cldn4_gene": 0.851188097681318,
        "claudins": 0.5419839859008789,
        "ctrl_epi": 0.5135294675827027,
        "ctrl_krt": 0.48922340869903563,
        "claudins_no4": 0.4170176848769188,
        "residual": 0.39787852401852264,
        "pre_core_q4": 0.2922197014093399,
        "pre_core_author": 0.13130335062742232,
    }
    for key, expect in checks.items():
        got = float(by_id[key].mean_delta)
        if abs(got - expect) > 1e-12:
            raise SystemExit(f"{key} file mean {got} != {expect}")
    if int(by_id["cldn4_gene"].n_patients) != 21:
        raise SystemExit("CLDN4 n is not 21")
    if f"{float(by_id['cldn4_gene'].mean_delta):.3f}" != "0.851":
        raise SystemExit("CLDN4 display rounding drifted")
    if abs(r_mean - float(by_id["residual"].mean_delta)) > 1e-12:
        raise SystemExit("residual recomputation does not match the published mean")

    rng = np.random.default_rng(131907)

    def row_from(frame, xcol, ycol, label, kind):
        rho, p, n, lo, hi = spearman_boot(frame[xcol], frame[ycol], rng)
        return {
            "label": label,
            "kind": kind,
            "xcol": xcol,
            "ycol": ycol,
            "rho": rho,
            "p": p,
            "n": n,
            "lo": lo,
            "hi": hi,
            "frame": frame,
        }

    assoc = [
        row_from(author, "tac_pct", "T_frac", "Mal. TACSTD2 % vs T", "primary"),
        row_from(author, "tac_mean", "T_frac", "Mal. TACSTD2 mean vs T", "prespec"),
        row_from(author, "tac_pct", "CD8_frac", "Mal. TACSTD2 % vs CD8", "prespec"),
        row_from(ts, "tac_pct", "T_frac", "Primary tS % vs T", "prespec"),
        row_from(ts, "tac_mean", "T_frac", "Primary tS mean vs T", "prespec"),
        row_from(ts, "tac_pct", "CD8_frac", "Primary tS % vs CD8", "prespec"),
        row_from(
            searched,
            "tac_pct2",
            "cyto_frac",
            "Searched, not confirmatory",
            "searched",
        ),
    ]

    for item in assoc:
        if item["kind"] == "searched":
            if abs(item["rho"] - const["searched_rho"]) > 1e-9:
                raise SystemExit(f"searched rho {item['rho']} != {const['searched_rho']}")
            if abs(item["p"] - const["searched_p"]) > 1e-9:
                raise SystemExit(f"searched p {item['p']}")
            continue
        gate = "author" if item["frame"] is author else "ts"
        hit = prespec[
            (prespec.gate == gate)
            & (prespec.predictor == item["xcol"])
            & (prespec.outcome == item["ycol"])
            & (prespec.denom == "gate_samples")
        ]
        if len(hit) != 1:
            raise SystemExit(f"no published row for {item['label']}")
        if abs(item["rho"] - float(hit.rho.iloc[0])) > 1e-9:
            raise SystemExit(f"rho mismatch {item['label']}: {item['rho']} vs {float(hit.rho.iloc[0])}")
        if abs(item["p"] - float(hit.p_spearman.iloc[0])) > 1e-9:
            raise SystemExit(f"p mismatch {item['label']}")
        if int(item["n"]) != int(hit.n_patients.iloc[0]):
            raise SystemExit(f"n mismatch {item['label']}")

    if abs(assoc[0]["rho"] - const["prespec_author_rho"]) > 1e-12:
        raise SystemExit("primary Tac% vs T rho drifted")
    if assoc[0]["n"] != 21:
        raise SystemExit("primary association n is not 21")
    if f"{assoc[0]['rho']:.2f}" != "-0.16":
        raise SystemExit("primary rho does not round to -0.16")

    pr_author, pp_author, _ = partial_spearman(author.tac_pct, author.T_frac, author.mal_frac)
    pr_ts, pp_ts, _ = partial_spearman(ts.tac_pct, ts.T_frac, ts.mal_frac)
    pr_sea, pp_sea, _ = partial_spearman(searched.tac_pct2, searched.cyto_frac, searched.mal_frac)
    for got, exp, name in (
        (pr_author, const["partial_rho_author_T"], "partial author rho"),
        (pp_author, const["partial_p_author_T"], "partial author p"),
        (pr_ts, const["partial_rho_ts_T"], "partial ts rho"),
        (pp_ts, const["partial_p_ts_T"], "partial ts p"),
        (pr_sea, const["partial_rho_searched"], "partial searched rho"),
        (pp_sea, const["partial_p_searched"], "partial searched p"),
    ):
        if abs(got - exp) > 1e-8:
            raise SystemExit(f"{name}: recomputed {got} != published {exp}")

    # Leave-one-out on the searched panel, and the cytotoxic-cell total.
    sea_x = searched.tac_pct2.to_numpy(float)
    sea_y = searched.cyto_frac.to_numpy(float)
    loo = []
    for i in range(len(searched)):
        m = np.ones(len(searched), dtype=bool)
        m[i] = False
        loo.append(float(stats.spearmanr(sea_x[m], sea_y[m]).pvalue))
    cyto_cells = float((searched.cyto_frac * searched.n_cells).sum())
    n_cyto_zero = int((searched.cyto_frac == 0).sum())

    return {
        "const": const,
        "contrasts": contrasts,
        "winner": winner.sort_values("delta", ascending=False),
        "w_mean": w_mean,
        "w_lo": w_lo,
        "w_hi": w_hi,
        "w_median": float(np.median(w_delta)),
        "w_p": w_p,
        "r_mean": r_mean,
        "r_lo": r_lo,
        "r_hi": r_hi,
        "r_median": float(np.median(w_resid)),
        "r_p": r_p,
        "c_mean": c_mean,
        "c_lo": c_lo,
        "c_hi": c_hi,
        "c_median": float(np.median(core_delta)),
        "c_p": c_p,
        "c_pos": int((core_delta > 0).sum()),
        "author": author,
        "assoc": assoc,
        "pr_author": pr_author,
        "pp_author": pp_author,
        "pr_ts": pr_ts,
        "pp_ts": pp_ts,
        "pr_sea": pr_sea,
        "pp_sea": pp_sea,
        "loo_min": float(min(loo)),
        "loo_max": float(max(loo)),
        "cyto_cells": cyto_cells,
        "n_cyto_zero": n_cyto_zero,
    }


def draw(ctx: dict) -> plt.Figure:
    # Nature double-column width. The caption is set apart from the art.
    w_mm, h_mm = 183.0, 222.0
    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4))

    def add(x, y, w, h):
        return fig.add_axes([x / w_mm, y / h_mm, w / w_mm, h / h_mm])

    ax_a = add(30, 160, 138, 48)
    ax_b = add(52, 80, 100, 60)
    ax_c = add(14, 16, 66, 44)
    ax_d = add(134, 16, 42, 44)
    cax = add(82.5, 22, 2.5, 32)

    _panel_a(ax_a, ctx)
    _panel_b(ax_b, ctx)
    _panel_c(ax_c, cax, ctx)
    _panel_d(ax_d, ctx)

    def letter(x, y, s):
        fig.text(x / w_mm, y / h_mm, s, fontsize=11, fontweight="bold", va="top", ha="left")

    def title(x, y, s):
        fig.text(x / w_mm, y / h_mm, s, fontsize=8, va="top", ha="left")

    letter(3.0, 218.5, "a")
    title(12, 218.5, "Within-patient claudin change, TACSTD2 Q4 minus Q1")
    letter(3.0, 146.5, "b")
    title(12, 146.5, "Mean paired change by gene set")
    letter(3.0, 66.5, "c")
    title(12, 66.5, "TACSTD2 percent versus T-cell fraction")
    letter(96, 66.5, "d")
    title(104, 66.5, "Spearman ρ")
    return fig


def _panel_a(ax, ctx):
    style_data_ax(ax)
    winner = ctx["winner"].reset_index(drop=True)
    y = np.arange(len(winner))[::-1]
    deltas = winner.delta.to_numpy(float)
    resids = winner.delta_resid_epi.to_numpy(float)
    for yi, d, r in zip(y, deltas, resids):
        ax.plot([0, d], [yi, yi], color=NAVY, lw=1.05, zorder=2, solid_capstyle="butt")
        ax.plot(
            r,
            yi,
            marker="o",
            ms=4.6,
            mfc="white",
            mec=TEAL,
            mew=1.05,
            zorder=3,
        )
        ax.plot(d, yi, marker="o", ms=4.6, mfc=NAVY, mec=NAVY, zorder=4)
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    ax.axvline(ctx["w_mean"], color=NAVY, lw=0.7, ls=(0, (3, 1.6)), zorder=1)
    ax.axvline(ctx["r_mean"], color=TEAL, lw=0.7, ls=(0, (1.2, 1.4)), zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(winner.patient.tolist(), fontsize=6.5)
    ax.set_ylim(-0.55, len(winner) - 0.35)
    ax.set_xlim(-0.04, 1.02)
    ax.set_xlabel("Paired Δ, log1p(CP10k)", fontsize=7, labelpad=2)
    ax.tick_params(axis="x", labelsize=6.5)

    # Cells per arm sit outside the axes so they are not read as a second effect.
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    for yi, n_arm in zip(y, winner.n_high.to_numpy(int)):
        ax.text(
            1.03,
            yi,
            f"{n_arm:,}",
            transform=trans,
            va="center",
            ha="left",
            fontsize=6,
            color=GREY,
            clip_on=False,
        )
    ax.text(
        1.03,
        1.04,
        "Cells/arm",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6,
        color=GREY,
        clip_on=False,
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=NAVY,
            mfc=NAVY,
            ms=4.2,
            lw=1.05,
            label="CLDN1/3/4/7, mean 0.542",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=TEAL,
            mfc="white",
            mew=1.05,
            ms=4.2,
            lw=0.8,
            ls=(0, (1.2, 1.4)),
            label="Residual, mean 0.398",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.01),
        ncol=2,
        frameon=False,
        fontsize=6.5,
        borderaxespad=0.0,
        handlelength=1.5,
        columnspacing=1.0,
    )


def _panel_b(ax, ctx):
    style_forest(ax)
    contrasts = {r.id: r for r in ctx["contrasts"].itertuples(index=False)}
    # Top to bottom. Section labels occupy their own rows so a different gate
    # cannot be read as part of the n = 10 contrast.
    layout = [
        ("section", "Single gene, not a module · epithelial cells, brain excluded"),
        ("row", "cldn4_gene"),
        ("section", "Same contrast · epithelial cells, tL/B + mLN + effusion · n = 10"),
        ("row", "claudins"),
        ("row", "ctrl_epi"),
        ("row", "ctrl_krt"),
        ("row", "claudins_no4"),
        ("row", "residual"),
        ("row", "pre_core_q4"),
        ("section", "Pre-specified · author-malignant cells · median split · n = 20"),
        ("row", "pre_core_author"),
    ]
    # reversed() puts the last layout entry at the bottom and the first at the top.
    y = 0.0
    placed = []
    for kind, key in reversed(layout):
        placed.append((y, kind, key))
        y += 0.72 if kind == "section" else 1.0

    ci_of = {
        "claudins": (ctx["w_lo"], ctx["w_hi"]),
        "residual": (ctx["r_lo"], ctx["r_hi"]),
        "pre_core_author": (ctx["c_lo"], ctx["c_hi"]),
    }
    style = {
        "cldn4_gene": dict(marker="o", mfc=VERM, mec=VERM),
        "claudins": dict(marker="o", mfc=NAVY, mec=NAVY),
        "ctrl_epi": dict(marker="s", mfc="white", mec=GREY),
        "ctrl_krt": dict(marker="s", mfc="white", mec=GREY),
        "claudins_no4": dict(marker="o", mfc=NAVY, mec=NAVY),
        "residual": dict(marker="o", mfc=TEAL, mec=TEAL),
        "pre_core_q4": dict(marker="o", mfc=NAVY, mec=NAVY),
        "pre_core_author": dict(marker="o", mfc=NAVY, mec=NAVY),
    }
    labels = {
        "cldn4_gene": "CLDN4 alone",
        "claudins": "CLDN1/3/4/7",
        "ctrl_epi": "EPCAM/KRT control",
        "ctrl_krt": "KRT8/18/19 control",
        "claudins_no4": "CLDN1/3/7",
        "residual": "Residual on EPCAM/KRT",
        "pre_core_q4": "12-gene core",
        "pre_core_author": "12-gene core",
    }

    yticks, ylabels = [], []
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    data_ys = []
    for yi, kind, key in placed:
        if kind == "section":
            ax.text(
                0.0,
                yi,
                key,
                transform=trans,
                ha="left",
                va="center",
                fontsize=6,
                color="#3F3F3F",
                fontstyle="italic",
                clip_on=False,
            )
            continue
        row = contrasts[key]
        mean = float(row.mean_delta)
        med = float(row.median_delta)
        data_ys.append(yi)
        yticks.append(yi)
        ylabels.append(labels[key])
        if key in ci_of:
            lo, hi = ci_of[key]
            color = TEAL if key == "residual" else NAVY
            ax.plot([lo, hi], [yi, yi], color=color, lw=1.0, zorder=3, solid_capstyle="butt")
            ax.plot([lo, lo], [yi - 0.16, yi + 0.16], color=color, lw=0.8, zorder=3)
            ax.plot([hi, hi], [yi - 0.16, yi + 0.16], color=color, lw=0.8, zorder=3)
        ax.plot([med, med], [yi - 0.18, yi + 0.18], color=INK, lw=0.65, zorder=4)
        st = style[key]
        ax.plot(
            mean,
            yi,
            marker=st["marker"],
            ms=4.8,
            mfc=st["mfc"],
            mec=st["mec"],
            mew=0.9,
            zorder=5,
        )
        n = int(row.n_patients)
        n_pos = int(round(float(row.frac_pos) * n))
        num_color = {"cldn4_gene": VERM, "residual": TEAL, "ctrl_epi": GREY, "ctrl_krt": GREY}.get(key, NAVY)
        ax.text(
            1.02,
            yi,
            f"{mean:.3f}   {n_pos}/{n}",
            transform=trans,
            va="center",
            ha="left",
            fontsize=6,
            color=num_color,
            clip_on=False,
        )

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(min(v[0] for v in placed) - 0.55, max(v[0] for v in placed) + 0.45)
    ax.set_xlim(-0.02, 1.0)
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    ax.set_xlabel("Mean paired Δ, log1p(CP10k)", fontsize=7, labelpad=2)


def _panel_c(ax, cax, ctx):
    style_data_ax(ax)
    author = ctx["author"]
    x = author.tac_pct.to_numpy(float) * 100
    y = author.T_frac.to_numpy(float)
    mal = author.mal_frac.to_numpy(float)
    sc = ax.scatter(
        x,
        y,
        c=mal,
        cmap="cividis",
        vmin=0,
        vmax=1,
        s=28,
        zorder=3,
        edgecolors="white",
        linewidths=0.4,
    )
    # Descriptive least-squares line. The test is the Spearman correlation.
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(0, 100, 50)
    ax.plot(xs, intercept + slope * xs, color=GREY, lw=0.7, ls=(0, (3, 1.8)), zorder=2)
    ax.set_xlim(-3, 105)
    ax.set_ylim(-0.02, 0.68)
    ax.set_ylabel("T-cell fraction", fontsize=7, labelpad=2)
    rho = ctx["assoc"][0]["rho"]
    p = ctx["assoc"][0]["p"]
    # Below axes-y 0.30 and left of axes-x 0.55 there are no patients.
    ax.text(
        0.03,
        0.22,
        f"ρ = {nfmt(rho, 3)}\nP = {pfmt(p)}\nn = 21\n"
        f"Partial ρ = {nfmt(ctx['pr_author'], 3)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        linespacing=1.25,
        zorder=4,
    )
    cb = plt.colorbar(sc, cax=cax)
    cb.outline.set_linewidth(0.4)
    cb.set_ticks([0, 0.5, 1])
    cb.ax.tick_params(labelsize=6, width=0.4, length=2)
    cb.set_label("Malignant-cell\nfraction", fontsize=6, labelpad=1)
    ax.set_xlabel("TACSTD2-positive malignant cells (%)\nGrey line, descriptive OLS. Test is Spearman.", fontsize=6.5, labelpad=2)


def _panel_d(ax, ctx):
    style_forest(ax)
    assoc = ctx["assoc"]
    # Primary test at the top. Extra gap before the searched row at the bottom.
    y_pos = []
    cursor = (len(assoc) - 1) + 0.7
    for item in assoc:
        if item["kind"] == "searched":
            cursor -= 0.7
        y_pos.append(cursor)
        cursor -= 1.0

    partial_at = {
        0: ctx["pr_author"],
        3: ctx["pr_ts"],
        6: ctx["pr_sea"],
    }
    for i, (item, yi) in enumerate(zip(assoc, y_pos)):
        if item["kind"] == "primary":
            ax.axhspan(yi - 0.42, yi + 0.42, color="#E6F1F8", zorder=0, lw=0)
            color = NAVY
        elif item["kind"] == "searched":
            ax.axhspan(yi - 0.42, yi + 0.42, color="#FBF4E6", zorder=0, lw=0)
            color = ORANGE
        else:
            color = NAVY
        ax.plot([item["lo"], item["hi"]], [yi, yi], color=color, lw=1.0, zorder=3, solid_capstyle="butt")
        ax.plot([item["lo"], item["lo"]], [yi - 0.14, yi + 0.14], color=color, lw=0.8, zorder=3)
        ax.plot([item["hi"], item["hi"]], [yi - 0.14, yi + 0.14], color=color, lw=0.8, zorder=3)
        ax.plot(item["rho"], yi, marker="o", ms=4.4, mfc=color, mec=color, zorder=4)
        if i in partial_at:
            ax.plot(
                partial_at[i],
                yi,
                marker="D",
                ms=3.6,
                mfc="white",
                mec=color,
                mew=0.8,
                zorder=5,
            )
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    ax.set_yticks(y_pos)
    tick_labels = ax.set_yticklabels([a["label"] for a in assoc])
    tick_labels[0].set_fontweight("bold")
    tick_labels[-1].set_color("#8A5A00")
    span = max(abs(a["lo"]) for a in assoc + [{"lo": -1}]) 
    span = max(span, max(abs(a["hi"]) for a in assoc), 1.0)
    ax.set_xlim(-span - 0.08, span + 0.08)
    ax.set_ylim(min(y_pos) - 0.7, max(y_pos) + 0.7)
    ax.set_xlabel("Spearman ρ\nBar, 95% bootstrap CI. Diamond, partial ρ.", fontsize=6.5, labelpad=2)
    if span <= 1.05:
        ax.set_xlim(-1.08, 1.08)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])


def write_caption(ctx: dict) -> None:
    cldn4 = ctx["contrasts"].loc[ctx["contrasts"].id == "cldn4_gene"].iloc[0]
    krt = ctx["contrasts"].loc[ctx["contrasts"].id == "ctrl_krt"].iloc[0]
    epi = ctx["contrasts"].loc[ctx["contrasts"].id == "ctrl_epi"].iloc[0]
    no4 = ctx["contrasts"].loc[ctx["contrasts"].id == "claudins_no4"].iloc[0]
    core = ctx["contrasts"].loc[ctx["contrasts"].id == "pre_core_q4"].iloc[0]
    primary = ctx["assoc"][0]
    ts = next(a for a in ctx["assoc"] if a["xcol"] == "tac_pct" and a["ycol"] == "T_frac" and a["n"] == 10)
    cd8 = next(a for a in ctx["assoc"] if a["xcol"] == "tac_pct" and a["ycol"] == "CD8_frac" and a["n"] == 10)
    sea = ctx["assoc"][-1]
    const = ctx["const"]
    text = f"""# Extended Data figure caption

**Tight-junction scores increase between TACSTD2 Q1 and Q4 in GSE131907, and TACSTD2 abundance does not track T-cell fraction.**

**a**, Within-patient change in the claudin score (mean of *CLDN1*, *CLDN3*, *CLDN4* and *CLDN7*; log1p of counts per 10,000 UMIs) for epithelial cells in metastatic lung (tL/B), lymph node (mLN) and pleural effusion. The contrast is TACSTD2 quartile 4 minus quartile 1 inside each patient. Filled circles are that paired difference. Open circles are the same score after a cell-level linear residual on *EPCAM*, *KRT8*, *KRT18* and *KRT19*. Dashed lines mark the unweighted means, {nfmt(ctx['w_mean'])} for the claudin score and {nfmt(ctx['r_mean'])} for the residual. All 10 of 10 patients are positive on both scores (exact two-sided Wilcoxon signed-rank *P* = {pfmt(ctx['w_p'])}). The median claudin difference is {nfmt(ctx['w_median'])}. Grey numbers are epithelial cells in each arm; quartile arms are the same size, and the mean is not weighted by cell count. *n* = 10 patients.

**b**, Mean paired difference. The *CLDN4*-only row is a different contrast: epithelial cells with brain metastasis excluded, TACSTD2 Q4 versus Q1, *n* = 21, mean {nfmt(float(cldn4.mean_delta))}, median {nfmt(float(cldn4.median_delta))}, 21/21 positive, exact *P* = {pfmt(float(cldn4.wilcoxon_p))}. That row is one gene, not a module, and the released table does not contain the 21 patient-level differences, so no interval is drawn. The middle block is one shared contrast (*n* = 10): *CLDN1/3/4/7* mean {nfmt(ctx['w_mean'])}; epithelial control *EPCAM/KRT8/18/19* mean {nfmt(float(epi.mean_delta))}; *KRT8/18/19* alone mean {nfmt(float(krt.mean_delta))}; *CLDN1/3/7* mean {nfmt(float(no4.mean_delta))}; claudin score residualized on the epithelial control, mean {nfmt(ctx['r_mean'])}; pre-specified 12-gene core mean {nfmt(float(core.mean_delta))}. The 12-gene core is *CLDN1*, *CLDN3*, *CLDN4*, *CLDN7*, *OCLN*, *TJP1*, *TJP2*, *TJP3*, *F11R*, *CGN*, *MARVELD2* and *CRB3*. The bottom row is the pre-specified core on author-malignant cells with a median split (*n* = 20, mean {nfmt(ctx['c_mean'])}, 18/20 positive, exact *P* = {pfmt(ctx['c_p'])}). Bars are 95% *t* intervals of the patient-level mean and are drawn only for rows whose patient-level differences are in the source table. Every *n* = 10 contrast in which all 10 differences are positive has the same exact signed-rank *P* = {pfmt(ctx['w_p'])}; the panel compares effect size. *TACSTD2* is not a member of any score. Epithelial controls are open squares.

**c**, Pre-specified test on author-malignant cells: percentage of cells with *TACSTD2* UMI > 0 versus the T-cell fraction of the sample. One point is one patient (*n* = 21). Colour is the malignant-cell fraction. Spearman ρ = {nfmt(primary['rho'], 3)}, *P* = {pfmt(primary['p'])}. The grey line is ordinary least squares and is not the test. Partial Spearman ρ given the malignant-cell fraction is {nfmt(ctx['pr_author'], 3)} (*P* = {pfmt(ctx['pp_author'])}). Author-annotated malignant cells in this matrix are metastatic (mBrain, tL/B, mLN), not primary lung.

**d**, Spearman ρ for pre-specified patient-level tests, with 95% percentile bootstrap intervals (5,000 resamples of patients). The highlighted row is panel **c**. Primary tumour epithelium (tS) versus T-cell fraction is ρ = {nfmt(ts['rho'], 3)}, *P* = {pfmt(ts['p'])}, *n* = 10. The same predictor versus CD8 fraction is ρ = {nfmt(cd8['rho'], 3)}, *P* = {pfmt(cd8['p'])} (the point estimate is positive). The orange row is the smallest nominal depletion *P* in the search: author-malignant brain metastases, *TACSTD2* fraction with UMI ≥ 2, versus the cytotoxic CD8 fraction (ρ = {nfmt(sea['rho'], 3)}, nominal *P* = {pfmt(sea['p'])}, *n* = 10). It is not a confirmatory result. Partial ρ on the malignant-cell fraction is {nfmt(ctx['pr_sea'], 3)} (*P* = {pfmt(ctx['pp_sea'])}). Across these 10 patients the cytotoxic outcome sums to {ctx['cyto_cells']:.0f} cells, and {ctx['n_cyto_zero']} of 10 patients are zero. Leave-one-out Spearman *P* ranges from {ctx['loo_min']:.3f} to {ctx['loo_max']:.2f}. A permutation *P* for this single panel is {pfmt(const['permutation_p_searched'])} ({int(const['permutation_n']):,} shuffles). The Bonferroni 0.05 threshold for the Spearman family ({int(const['spearman_family_n']):,} tests) is {pfmt(const['bonferroni_0.05'])}. Open diamonds are partial ρ given the malignant-cell fraction.

Inferential *n* is the patient. The matrix contains {int(const['n_cells_not_inferential']):,} cells; that count is not a sample size. Dissociated counts are not a spatial exclusion test.
"""
    (ROOT / "CAPTION.md").write_text(text)
    # The three numbers this figure is not allowed to drift away from.
    required = ("0.542", "0.851", "−0.157", "n* = 10", "n* = 21")
    # The caption uses markdown italics around n, so check the plain phrases.
    plain_checks = ("0.542", "0.851", "−0.157", "*n* = 10", "*n* = 21")
    missing = [s for s in plain_checks if s not in text]
    if missing:
        raise SystemExit(f"caption missing {missing}")
    _ = required


def write_audit(ctx: dict) -> None:
    audit = {
        "claudin_mean_delta": ctx["w_mean"],
        "claudin_median_delta": ctx["w_median"],
        "claudin_mean_ci95": [ctx["w_lo"], ctx["w_hi"]],
        "claudin_n": 10,
        "claudin_wilcoxon_p": ctx["w_p"],
        "residual_mean_delta": ctx["r_mean"],
        "residual_mean_ci95": [ctx["r_lo"], ctx["r_hi"]],
        "cldn4_mean_delta": 0.851188097681318,
        "cldn4_n": 21,
        "tac_vs_T_rho": ctx["assoc"][0]["rho"],
        "tac_vs_T_p": ctx["assoc"][0]["p"],
        "tac_vs_T_n": 21,
        "tac_vs_T_bootstrap_ci95": [ctx["assoc"][0]["lo"], ctx["assoc"][0]["hi"]],
        "partial_rho_author_T": ctx["pr_author"],
        "partial_p_author_T": ctx["pp_author"],
        "searched_rho": ctx["assoc"][-1]["rho"],
        "searched_p": ctx["assoc"][-1]["p"],
        "searched_loo_p_range": [ctx["loo_min"], ctx["loo_max"]],
        "searched_cyto_cells": ctx["cyto_cells"],
        "searched_cyto_zero_patients": ctx["n_cyto_zero"],
        "prespec_core_mean": ctx["c_mean"],
        "prespec_core_positive": f"{ctx['c_pos']}/20",
        "associations": [
            {
                "label": a["label"],
                "rho": a["rho"],
                "p": a["p"],
                "n": a["n"],
                "ci95": [a["lo"], a["hi"]],
            }
            for a in ctx["assoc"]
        ],
    }
    (SRC / "figure_audit.json").write_text(json.dumps(audit, indent=2) + "\n")


def main():
    ctx = prepare()
    fig = draw(ctx)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(ROOT / f"ed_gse131907_tj.{ext}")
    plt.close(fig)
    write_caption(ctx)
    write_audit(ctx)
    print("wrote", ROOT / "ed_gse131907_tj.pdf")
    print(
        f"claudin mean {ctx['w_mean']:.3f} n=10; "
        f"CLDN4 0.851 n=21; "
        f"Tac vs T rho {ctx['assoc'][0]['rho']:.3f} n=21"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""TCGA-LUSC CLDN4 Q4 vs Q1 vs ImmuneScore / CD8A / CD274.

Additive LUSC-only slice. Predictor is CLDN4 RNA (Xena HiSeqV2).
No TACSTD2 gate. Quartiles among primary tumors with non-NA CLDN4.
Endpoints: official ESTIMATE ImmuneScore, CD8A, CD274 (PD-L1).
Honest pairwise-complete n.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GENES = ("CLDN4", "CD8A", "CD274")
ENDPOINTS = ("ImmuneScore", "CD8A", "CD274")
RNG = np.random.default_rng(20260817)
N_BOOT = 2000


def is_primary_01(barcode: str) -> bool:
    parts = str(barcode).split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def sample15(barcode: str) -> str:
    return "-".join(str(barcode).split("-")[:4])[:15]


def load_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t", index_col=0)
    missing = [g for g in GENES if g not in expr.index]
    if missing:
        raise SystemExit(f"missing genes in Xena LUSC matrix: {missing}")
    keep = [c for c in expr.columns if is_primary_01(c)]
    slim = expr.loc[list(GENES), keep]
    slim.columns = [sample15(c) for c in slim.columns]
    slim = slim.T.groupby(level=0).mean().T
    return slim


def load_estimate(path: Path) -> pd.DataFrame:
    est = pd.read_csv(path, sep="\t")
    est.columns = [c.strip() for c in est.columns]
    idcol = est.columns[0]
    rename = {}
    for c in est.columns:
        cl = c.lower().replace(" ", "_")
        if cl in {"stromal_score", "stromalscore"}:
            rename[c] = "ESTIMATE_StromalScore"
        elif cl in {"immune_score", "immunescore"}:
            rename[c] = "ESTIMATE_ImmuneScore"
        elif cl in {"estimate_score", "estimatescore"}:
            rename[c] = "ESTIMATE_Score"
    est = est.rename(columns={idcol: "sample_raw", **rename})
    if "ESTIMATE_ImmuneScore" not in est.columns:
        raise SystemExit(f"ImmuneScore column not found: {list(est.columns)}")
    est = est[est["sample_raw"].astype(str).map(is_primary_01)].copy()
    est["sample"] = est["sample_raw"].map(sample15)
    out = est.groupby("sample")[["ESTIMATE_ImmuneScore"]].mean()
    if "ESTIMATE_StromalScore" in est.columns:
        out["ESTIMATE_StromalScore"] = est.groupby("sample")["ESTIMATE_StromalScore"].mean()
    if "ESTIMATE_Score" in est.columns:
        out["ESTIMATE_Score"] = est.groupby("sample")["ESTIMATE_Score"].mean()
    return out


def quartiles(s: pd.Series) -> pd.Series:
    """Equal-count Q1–Q4 on non-NA values. rank-first so ties still fill 4 bins."""
    out = pd.Series(np.nan, index=s.index, dtype=object)
    ok = s.dropna()
    if len(ok) < 8:
        return out
    labels = pd.qcut(ok.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    out.loc[ok.index] = labels.astype(str)
    return out


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    if n < 4 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": int(n), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def bootstrap_rb(a: np.ndarray, b: np.ndarray, n_boot: int = N_BOOT) -> tuple[float, float]:
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    vals = []
    for _ in range(n_boot):
        aa = RNG.choice(a, size=len(a), replace=True)
        bb = RNG.choice(b, size=len(b), replace=True)
        U, _ = stats.mannwhitneyu(aa, bb, alternative="two-sided")
        vals.append((2.0 * U) / (len(aa) * len(bb)) - 1.0)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def mannwhitney(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    rec = {
        "n_q4": int(len(a)),
        "n_q1": int(len(b)),
        "median_q4": float(a.median()) if len(a) else np.nan,
        "median_q1": float(b.median()) if len(b) else np.nan,
        "delta_median_q4_minus_q1": np.nan,
        "U": np.nan,
        "p": np.nan,
        "rank_biserial_q4_gt_q1": np.nan,
        "rank_biserial_lo95": np.nan,
        "rank_biserial_hi95": np.nan,
    }
    if len(a) and len(b):
        rec["delta_median_q4_minus_q1"] = float(a.median() - b.median())
    if len(a) < 3 or len(b) < 3:
        return rec
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = (2.0 * U) / (len(a) * len(b)) - 1.0
    lo, hi = bootstrap_rb(a.to_numpy(dtype=float), b.to_numpy(dtype=float))
    rec.update(
        {
            "U": float(U),
            "p": float(p),
            "rank_biserial_q4_gt_q1": float(r),
            "rank_biserial_lo95": lo,
            "rank_biserial_hi95": hi,
        }
    )
    return rec


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_q4_vs_q1(core: pd.DataFrame, mw: pd.DataFrame, path: Path) -> None:
    titles = {
        "ImmuneScore": "ESTIMATE ImmuneScore",
        "CD8A": "CD8A (CD8)",
        "CD274": "CD274 (PD-L1)",
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.8))
    sub = core[core["quartile"].isin(["Q1", "Q4"])].copy()
    for ax, ep in zip(axes, ENDPOINTS):
        a = sub.loc[sub["quartile"] == "Q1", ep].dropna()
        b = sub.loc[sub["quartile"] == "Q4", ep].dropna()
        parts = ax.violinplot(
            [a.to_numpy(dtype=float), b.to_numpy(dtype=float)],
            positions=[1, 2],
            showextrema=False,
            widths=0.7,
        )
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor("#d9e8f5" if i == 0 else "#f5d6d0")
            body.set_edgecolor("#2c3e50")
            body.set_alpha(0.85)
        ax.boxplot(
            [a, b],
            widths=0.18,
            patch_artist=True,
            boxprops=dict(facecolor="white", edgecolor="#2c3e50"),
            medianprops=dict(color="#c0392b", linewidth=1.6),
            whiskerprops=dict(color="#2c3e50"),
            capprops=dict(color="#2c3e50"),
            flierprops=dict(marker="o", markersize=2.5, markerfacecolor="#7f8c8d", alpha=0.6),
        )
        row = mw[mw.endpoint == ep].iloc[0]
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"Q1\nn={len(a)}", f"Q4\nn={len(b)}"])
        ax.set_title(f"{titles[ep]}\nMWU p={row['p']:.3g}", fontsize=9)
        ax.set_ylabel(ep, fontsize=8)
        style(ax)
    fig.suptitle(
        "TCGA-LUSC  ·  CLDN4 RNA Q4 vs Q1  ·  no TACSTD2 gate",
        fontsize=11,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_forest_extra(mw: pd.DataFrame, path: Path) -> None:
    """Extra figure: rank-biserial forest with honest n per arm."""
    labels = {
        "ImmuneScore": "ImmuneScore",
        "CD8A": "CD8A (CD8)",
        "CD274": "CD274 (PD-L1)",
    }
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    y = np.arange(len(ENDPOINTS))[::-1]
    for i, ep in enumerate(ENDPOINTS):
        row = mw[mw.endpoint == ep].iloc[0]
        r = row["rank_biserial_q4_gt_q1"]
        lo = row["rank_biserial_lo95"]
        hi = row["rank_biserial_hi95"]
        color = "#c0392b" if row["p"] < 0.05 else "#2c3e50"
        ax.plot([lo, hi], [y[i], y[i]], color=color, lw=1.8)
        ax.plot(r, y[i], "o", color=color, ms=7)
        ax.text(
            1.02,
            y[i],
            f"δ_rb={r:+.3f}  p={row['p']:.3g}  n={int(row['n_q1'])}/{int(row['n_q4'])}",
            va="center",
            ha="left",
            fontsize=8,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color="#7f8c8d", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([labels[ep] for ep in ENDPOINTS], fontsize=9)
    ax.set_xlabel("Rank-biserial  (positive = CLDN4 Q4 higher)", fontsize=9)
    ax.set_xlim(-0.45, 0.45)
    ax.set_title(
        "Extra: TCGA-LUSC CLDN4 Q4 vs Q1  ·  rank-biserial with 95% bootstrap CI",
        fontsize=10,
    )
    style(ax)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/tcga_lusc_cldn4_q4")
    ap.add_argument("--outdir", default="methods/tcga_lusc_cldn4_q4")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    res = out / "results"
    figdir = out / "figures"
    res.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    expr_path = data / "TCGA.LUSC.HiSeqV2.gz"
    est_path = data / "ESTIMATE_LUSC_RNAseqV2.txt"
    expr = load_expression(expr_path)
    est = load_estimate(est_path)

    n_expr_primary = int(expr.shape[1])
    n_estimate_primary = int(est.shape[0])

    core = expr.T.copy()
    core.index.name = "sample"
    core = core.join(est, how="left")
    core = core.rename(columns={"ESTIMATE_ImmuneScore": "ImmuneScore"})
    core["CLDN4"] = pd.to_numeric(core["CLDN4"], errors="coerce")
    core["CD8A"] = pd.to_numeric(core["CD8A"], errors="coerce")
    core["CD274"] = pd.to_numeric(core["CD274"], errors="coerce")
    core["ImmuneScore"] = pd.to_numeric(core["ImmuneScore"], errors="coerce")

    n_intersect_expr_est = int(core["ImmuneScore"].notna().sum())
    # Quartiles among CLDN4-complete primaries (not among ImmuneScore-complete only).
    core["quartile"] = quartiles(core["CLDN4"])
    core["tacstd2_gate"] = "none"

    n_cldn4 = int(core["CLDN4"].notna().sum())
    n_cldn4_na = int(core["CLDN4"].isna().sum())
    q_counts = core["quartile"].value_counts(dropna=True).to_dict()
    cuts = core.loc[core["CLDN4"].notna(), "CLDN4"]
    q_cuts = cuts.quantile([0.25, 0.5, 0.75]).to_dict() if n_cldn4 else {}

    endpoint_n = {}
    for ep in ENDPOINTS:
        endpoint_n[ep] = {
            "n_non_na_among_cldn4": int(core.loc[core["CLDN4"].notna(), ep].notna().sum()),
            "n_na_among_cldn4": int(core.loc[core["CLDN4"].notna(), ep].isna().sum()),
        }

    mw_rows = []
    sp_rows = []
    for ep in ENDPOINTS:
        q4 = core.loc[core["quartile"] == "Q4", ep]
        q1 = core.loc[core["quartile"] == "Q1", ep]
        r = mannwhitney(q4, q1)
        r.update(
            {
                "cohort": "TCGA-LUSC",
                "predictor": "CLDN4_RNA",
                "endpoint": ep,
                "test": "MWU_Q4_vs_Q1",
                "tacstd2_gate": "none",
            }
        )
        mw_rows.append(r)
        s = spearman(core["CLDN4"], core[ep])
        s.update(
            {
                "cohort": "TCGA-LUSC",
                "predictor": "CLDN4_RNA",
                "endpoint": ep,
                "test": "spearman_continuous",
                "tacstd2_gate": "none",
            }
        )
        sp_rows.append(s)

    mw = pd.DataFrame(mw_rows)
    sp = pd.DataFrame(sp_rows)

    n_row = {
        "cohort": "TCGA-LUSC",
        "n_xena_primary_01": n_expr_primary,
        "n_estimate_primary_01": n_estimate_primary,
        "n_cldn4": n_cldn4,
        "n_cldn4_na": n_cldn4_na,
        "n_ImmuneScore_join": n_intersect_expr_est,
        "n_Q1": int(q_counts.get("Q1", 0)),
        "n_Q2": int(q_counts.get("Q2", 0)),
        "n_Q3": int(q_counts.get("Q3", 0)),
        "n_Q4": int(q_counts.get("Q4", 0)),
        "n_ImmuneScore_among_cldn4": endpoint_n["ImmuneScore"]["n_non_na_among_cldn4"],
        "n_CD8A_among_cldn4": endpoint_n["CD8A"]["n_non_na_among_cldn4"],
        "n_CD274_among_cldn4": endpoint_n["CD274"]["n_non_na_among_cldn4"],
        "n_ImmuneScore_NA_among_cldn4": endpoint_n["ImmuneScore"]["n_na_among_cldn4"],
        "n_CD8A_NA_among_cldn4": endpoint_n["CD8A"]["n_na_among_cldn4"],
        "n_CD274_NA_among_cldn4": endpoint_n["CD274"]["n_na_among_cldn4"],
        "tacstd2_gate": "none",
        "quartile_rule": "pd.qcut(rank(method='first'), 4) among CLDN4-complete primaries",
    }
    ntab = pd.DataFrame([n_row])

    core.to_csv(res / "sample_scores.tsv", sep="\t")
    mw.to_csv(res / "q4_vs_q1.tsv", sep="\t", index=False)
    sp.to_csv(res / "spearman.tsv", sep="\t", index=False)
    ntab.to_csv(res / "n_table.tsv", sep="\t", index=False)

    plot_q4_vs_q1(core, mw, figdir / "fig1_q4_vs_q1.png")
    plot_forest_extra(mw, figdir / "fig2_forest_extra.png")

    info = {
        "cohort": "TCGA-LUSC",
        "predictor": "CLDN4 RNA (Xena HiSeqV2 log2(RSEM+1))",
        "tacstd2_gate": "none",
        "n_xena_primary_01": n_expr_primary,
        "n_estimate_primary_01": n_estimate_primary,
        "n_cldn4": n_cldn4,
        "n_cldn4_na": n_cldn4_na,
        "quartile_counts": {k: int(v) for k, v in q_counts.items()},
        "cldn4_q25": float(q_cuts.get(0.25, np.nan)),
        "cldn4_q50": float(q_cuts.get(0.50, np.nan)),
        "cldn4_q75": float(q_cuts.get(0.75, np.nan)),
        "endpoint_n": endpoint_n,
        "genes_present": {g: bool(g in expr.index) for g in GENES},
    }
    summary = {
        "question": (
            "TCGA-LUSC CLDN4 RNA Q4 vs Q1 versus ImmuneScore, CD8A, CD274. "
            "No TACSTD2 gate. Honest pairwise-complete n."
        ),
        "predictor": "CLDN4 RNA",
        "tacstd2_gate": "none",
        "quartile_rule": n_row["quartile_rule"],
        "endpoints": {
            "ImmuneScore": "Official MD Anderson ESTIMATE RNAseqV2 Immune_score",
            "CD8A": "Xena HiSeqV2 CD8A (CD8 gene; not a deconvolution fraction)",
            "CD274": "Xena HiSeqV2 CD274 (PD-L1 RNA; not IHC CPS)",
        },
        "cohort": info,
        "q4_vs_q1": mw.to_dict(orient="records"),
        "spearman": sp.to_dict(orient="records"),
        "n_table": n_row,
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"n": n_row, "q4_vs_q1": mw.to_dict(orient="records"), "spearman": sp.to_dict(orient="records")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

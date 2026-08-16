#!/usr/bin/env python3
"""Figures for the TF-network hunt. No extra statistics — just the tables."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TFS = ["ELF3", "GRHL1", "KLF4", "TFAP2A", "NKX2-1"]
TARGETS = ["TACSTD2", "CLDN4"]


def heatmap(pair: pd.DataFrame, scores: pd.DataFrame, out: Path) -> None:
    raw = pair[(pair["scale"] == "raw") & (pair["tf"].isin(TFS)) & (pair["target"].isin(TARGETS))]
    if raw.empty:
        return
    raw = raw.copy()
    raw["pair"] = raw["tf"] + " vs " + raw["target"]
    # keep cohorts that are not UNINFORMATIVE
    keep = scores.loc[scores["call_raw"] != "UNINFORMATIVE", "cohort_id"]
    raw = raw[raw["cohort_id"].isin(keep)]
    if raw.empty:
        return
    order = scores.set_index("cohort_id").loc[raw["cohort_id"].unique()].sort_values("score", ascending=False).index
    pair_order = [f"{tf} vs {tgt}" for tf in TFS for tgt in TARGETS]
    mat = (
        raw.pivot_table(index="cohort_id", columns="pair", values="rho", aggfunc="mean")
        .reindex(index=order, columns=pair_order)
    )
    fig_h = max(4.0, 0.28 * len(mat) + 2.2)
    fig, ax = plt.subplots(figsize=(10.5, fig_h))
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(mat.index)))
    labels = []
    for cid in mat.index:
        r = scores.set_index("cohort_id").loc[cid]
        labels.append(f"{cid}  {r['organism'][:1]}/{r['material'][:3]} n={int(r['n_used'])} {r['call_raw']}")
    ax.set_yticklabels(labels, fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Spearman ρ")
    ax.set_title("Raw Spearman: TFs vs TACSTD2 / CLDN4")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def score_forest(scores: pd.DataFrame, out: Path) -> None:
    df = scores[scores["call_raw"] != "UNINFORMATIVE"].copy()
    if df.empty:
        return
    df = df.sort_values("score")
    colors = {"SUPPORTED": "#1b9e77", "PARTIAL": "#d95f02", "NOT_SUPPORTED": "#666666"}
    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.28 * len(df) + 1.5)))
    y = np.arange(len(df))
    ax.axvline(0, color="0.7", lw=1)
    ax.scatter(df["score"], y, c=[colors.get(c, "0.4") for c in df["call_raw"]], s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.cohort_id} ({r.organism[0]}/{r.material[:3]}, n={int(r.n_used)})" for r in df.itertuples()],
        fontsize=7,
    )
    ax.set_xlabel("score = mean ρ(TFs, targets) − mean ρ(NKX2-1, targets)")
    ax.set_title("Pre-specified network score (raw Spearman)")
    for call, col in colors.items():
        ax.scatter([], [], c=col, label=call)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def raw_vs_partial(pair: pd.DataFrame, out: Path) -> None:
    raw = pair[(pair["scale"] == "raw") & (pair["tf"].isin(TFS[:4]))]
    ep = pair[(pair["scale"] == "epcam") & (pair["tf"].isin(TFS[:4]))]
    m = raw.merge(ep, on=["cohort_id", "tf", "target"], suffixes=("_raw", "_epcam"))
    if m.empty:
        return
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.axhline(0, color="0.8", lw=1)
    ax.axvline(0, color="0.8", lw=1)
    ax.plot([-1, 1], [-1, 1], color="0.7", lw=1)
    ax.scatter(m["rho_raw"], m["rho_epcam"], s=12, alpha=0.55, c="#345995")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("raw Spearman ρ")
    ax.set_ylabel("EPCAM-residual Spearman ρ")
    ax.set_title("Positive TF–target pairs: raw vs EPCAM-partial")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def composition_bars(scores: pd.DataFrame, out: Path) -> None:
    cols = [
        "comp_TACSTD2__SCGB1A1",
        "comp_TACSTD2__SFTPC",
        "comp_CLDN4__SCGB1A1",
        "comp_CLDN4__SFTPC",
        "comp_NKX2-1__SFTPC",
        "comp_NKX2-1__SCGB1A1",
    ]
    have = [c for c in cols if c in scores.columns]
    if not have:
        return
    df = scores[scores["call_raw"] != "UNINFORMATIVE"]
    if df.empty:
        return
    means = df[have].mean(numeric_only=True)
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axhline(0, color="0.7", lw=1)
    ax.bar(range(len(means)), means.to_numpy(), color="#4c78a8")
    ax.set_xticks(range(len(means)))
    ax.set_xticklabels([c.replace("comp_", "") for c in means.index], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("mean Spearman ρ across informative cohorts")
    ax.set_title("Composition diagnostic: targets / NKX2-1 vs club vs AT2 markers")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="results/hunt_tf/tables/cohort_scores.tsv")
    ap.add_argument("--pairs", default="results/hunt_tf/tables/pair_correlations.tsv")
    ap.add_argument("--fig-dir", default="results/hunt_tf/figures")
    args = ap.parse_args()
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(args.scores, sep="\t")
    pairs = pd.read_csv(args.pairs, sep="\t")
    heatmap(pairs, scores, fig_dir / "heatmap_tf_target_rho.png")
    score_forest(scores, fig_dir / "score_forest.png")
    raw_vs_partial(pairs, fig_dir / "raw_vs_epcam_partial.png")
    composition_bars(scores, fig_dir / "composition_diagnostic.png")
    print("wrote", fig_dir)


if __name__ == "__main__":
    main()

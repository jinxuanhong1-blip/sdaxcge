#!/usr/bin/env python3
"""Figures for the combinatorial composition grid."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig_recovery_heatmap(rec: pd.DataFrame, out: Path) -> None:
    rec = rec[rec["method"].isin(["alr_perm", "dm_glm"])].copy()
    rec = rec[rec["method"] == "alr_perm"] if "alr_perm" in set(rec["method"]) else rec
    rec = rec[rec["compartment"].isin(["TNK", "TLS", "T", "NK", "B", "plasma"])]
    if rec.empty:
        return
    rec["key"] = rec["cohort"] + " × " + rec["annotation"]
    rec["lab"] = rec["compartment"].where(rec["compartment"].isin(["TNK", "TLS"]), rec["celltype"])
    # one row per key × lab: prefer collapsed names
    rec = rec.sort_values(["key", "lab", "celltype"]).drop_duplicates(["key", "lab"])
    keys = sorted(rec["key"].unique())
    labs = [c for c in ["TNK", "TLS", "T", "NK", "B", "plasma"] if c in set(rec["lab"])]
    mat = np.full((len(keys), len(labs)), np.nan)
    qmat = np.full_like(mat, np.nan)
    for i, k in enumerate(keys):
        for j, lab in enumerate(labs):
            hit = rec[(rec["key"] == k) & (rec["lab"] == lab)]
            if hit.empty:
                continue
            mat[i, j] = hit["effect"].iloc[0]
            qmat[i, j] = hit["q_bh"].iloc[0]
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1
    vmax = max(float(vmax), 0.2)
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.38 * len(keys) + 1.5)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(labs, rotation=0)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(keys, fontsize=8)
    for i in range(len(keys)):
        for j in range(len(labs)):
            if not np.isfinite(mat[i, j]):
                ax.text(j, i, "·", ha="center", va="center", color="0.5", fontsize=8)
                continue
            star = ""
            q = qmat[i, j]
            if np.isfinite(q) and q < 0.05:
                star = "**"
            elif np.isfinite(q) and q < 0.10:
                star = "*"
            ax.text(j, i, f"{mat[i, j]:+.2f}{star}", ha="center", va="center", fontsize=7)
    ax.set_title("DM-GLM ALR slope: TACSTD2-high vs low\n(negative = down in TACSTD2-high; * q<0.10, ** q<0.05, BH in family)")
    fig.colorbar(im, ax=ax, fraction=0.03, label="ALR slope (vs reference)")
    _save(fig, out / "fig1_recovery_heatmap.png")


def fig_n_bars(cov: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    rows = []
    for cohort, g in cov.groupby("cohort"):
        rows.append(
            {
                "cohort": cohort,
                "n_total": len(g),
                "n_tacstd2": int(g["eligible_tacstd2"].fillna(0).astype(bool).sum()),
                "n_mpr": int(g["eligible_mpr"].fillna(0).astype(bool).sum()),
            }
        )
    d = pd.DataFrame(rows).sort_values("cohort")
    x = np.arange(len(d))
    ax.bar(x - 0.2, d["n_total"], 0.2, label="all samples", color="#4C72B0")
    ax.bar(x, d["n_tacstd2"], 0.2, label="TACSTD2-eligible", color="#55A868")
    ax.bar(x + 0.2, d["n_mpr"], 0.2, label="MPR-eligible", color="#C44E52")
    ax.set_xticks(x)
    ax.set_xticklabels(d["cohort"], rotation=25, ha="right")
    ax.set_ylabel("n patients")
    ax.set_title("Honest n (patients, not cells)")
    ax.legend(frameon=False)
    _save(fig, out / "fig2_honest_n.png")


def fig_scatter(cov: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=False)
    cohorts = [c for c in ["GSE207422", "GSE241934_IIT", "GSE241934_RWC", "GSE253013"] if c in set(cov["cohort"])]
    for ax, cohort in zip(axes, cohorts[:3] if len(cohorts) >= 3 else cohorts + [None] * 3):
        if cohort is None:
            ax.axis("off")
            continue
        g = cov[cov["cohort"] == cohort].copy()
        g = g[g["eligible_tacstd2"].fillna(0).astype(bool)]
        x = g["tacstd2_malig_mean_log1p"]
        y = g["n_tnk"] / g["n_cells"]
        ax.scatter(x, y, c=np.where(g["mpr"] == "MPR", "#C44E52", "#4C72B0"), s=36, edgecolor="k", linewidth=0.3)
        ax.set_title(f"{cohort}\nn={len(g)}")
        ax.set_xlabel("malignant TACSTD2 mean log1p")
        ax.set_ylabel("T/NK fraction")
    fig.suptitle("Descriptive only: T/NK fraction vs malignant TACSTD2 (not the compositional test)")
    _save(fig, out / "fig3_tnk_vs_tacstd2_descriptive.png")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args()
    res = args.root / "results"
    rec = pd.read_csv(res / "grid_effects.tsv", sep="\t")
    cov = pd.read_csv(res / "sample_covariates.tsv", sep="\t")
    fig_recovery_heatmap(rec, res)
    fig_n_bars(cov, res)
    fig_scatter(cov, res)
    print("wrote figures in", res)


if __name__ == "__main__":
    main()

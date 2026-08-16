#!/usr/bin/env python3
"""Pool per-dataset Spearman ρ (epithelial TACSTD2/CLDN4 vs T/NK fraction).

Primary pool: Fisher z inverse-variance meta-analysis of dataset-level Spearman.
Secondary: rank-within-dataset then one pooled Spearman (batch-removed).
Neither is an ICI-response test.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from config import MIN_N_FOR_FISHER_Z, MIN_N_FOR_SPEARMAN


def fisher_z_meta(rows: pd.DataFrame) -> dict:
    use = rows.loc[rows["ok"] & (rows["n"] >= MIN_N_FOR_FISHER_Z)].copy()
    rec = {
        "n_datasets": int(len(use)),
        "n_units_sum": int(use["n"].sum()) if len(use) else 0,
        "rho": np.nan,
        "p": np.nan,
        "z": np.nan,
        "se": np.nan,
        "ci95_lo": np.nan,
        "ci95_hi": np.nan,
        "Q": np.nan,
        "I2": np.nan,
        "note": "",
    }
    if len(use) < 2:
        rec["note"] = f"need >=2 datasets with n>={MIN_N_FOR_FISHER_Z}"
        return rec
    z = np.arctanh(use["rho"].to_numpy(dtype=float))
    w = (use["n"].to_numpy(dtype=float) - 3.0)
    zbar = float(np.sum(w * z) / np.sum(w))
    se = float(1.0 / math.sqrt(np.sum(w)))
    zstat = zbar / se
    p = float(2 * stats.norm.sf(abs(zstat)))
    rec["z"] = zbar
    rec["se"] = se
    rec["rho"] = float(math.tanh(zbar))
    rec["p"] = p
    rec["ci95_lo"] = float(math.tanh(zbar - 1.96 * se))
    rec["ci95_hi"] = float(math.tanh(zbar + 1.96 * se))
    q = float(np.sum(w * (z - zbar) ** 2))
    rec["Q"] = q
    df = len(use) - 1
    rec["I2"] = float(max(0.0, (q - df) / q * 100)) if q > 0 else 0.0
    rec["note"] = "Fisher z IVW of dataset Spearman ρ"
    return rec


def rank_within_pool(units: pd.DataFrame, xcol: str) -> dict:
    parts = []
    for ds, sub in units.groupby("dataset"):
        if len(sub) < MIN_N_FOR_SPEARMAN:
            continue
        s = sub.copy()
        s["_x"] = s[xcol].rank()
        s["_y"] = s["frac_tnk"].rank()
        parts.append(s[["dataset", "_x", "_y", "unit_id"]])
    rec = {"n": 0, "n_datasets": 0, "rho": np.nan, "p": np.nan, "note": ""}
    if not parts:
        rec["note"] = "no dataset with n>=5 eligible units"
        return rec
    cat = pd.concat(parts, ignore_index=True)
    rec["n"] = int(len(cat))
    rec["n_datasets"] = int(cat["dataset"].nunique())
    if np.nanstd(cat["_x"]) == 0 or np.nanstd(cat["_y"]) == 0:
        rec["note"] = "zero variance after ranking"
        return rec
    rho, p = stats.spearmanr(cat["_x"], cat["_y"])
    rec["rho"] = float(rho)
    rec["p"] = float(p)
    rec["note"] = "Spearman after within-dataset ranks (batch removed)"
    return rec


def forest_plot(df: pd.DataFrame, title: str, out: Path, pooled: dict | None = None) -> None:
    use = df.sort_values("dataset").reset_index(drop=True)
    fig_h = max(2.8, 0.45 * (len(use) + (1 if pooled else 0)) + 1.4)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    y = np.arange(len(use))
    for i, r in use.iterrows():
        if not r.get("ok", False) or pd.isna(r["rho"]):
            ax.plot(0, i, marker="x", color="0.5")
            continue
        # approximate 95% CI via Fisher z
        n = float(r["n"])
        if n > 3:
            se = 1.0 / math.sqrt(n - 3)
            lo, hi = math.tanh(math.atanh(r["rho"]) - 1.96 * se), math.tanh(
                math.atanh(r["rho"]) + 1.96 * se
            )
            ax.plot([lo, hi], [i, i], color="0.25", lw=1.4)
        ax.plot(r["rho"], i, "o", color="#1f4e79", ms=6)
        ax.text(1.02, i, f"n={int(r['n'])}  ρ={r['rho']:+.2f}  p={r['p']:.3g}", va="center", fontsize=8)
    labels = list(use["dataset"].str.replace("NSCLC_", "", regex=False))
    if pooled and pooled.get("n_datasets", 0) >= 2 and not pd.isna(pooled.get("rho", np.nan)):
        y0 = len(use)
        ax.plot([pooled["ci95_lo"], pooled["ci95_hi"]], [y0, y0], color="#8b1e3f", lw=2)
        ax.plot(pooled["rho"], y0, "D", color="#8b1e3f", ms=7)
        ax.text(
            1.02,
            y0,
            f"pool n_ds={pooled['n_datasets']}  ρ={pooled['rho']:+.2f}  p={pooled['p']:.3g}",
            va="center",
            fontsize=8,
        )
        labels.append("Fisher-z pool")
        y = np.arange(len(labels))
    ax.axvline(0, color="0.6", lw=0.8)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Spearman ρ (epithelial score vs T/NK fraction)")
    ax.set_title(title, fontsize=10)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def scatter_plot(units: pd.DataFrame, xcol: str, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    datasets = sorted(units["dataset"].unique())
    cmap = plt.get_cmap("tab10")
    for i, ds in enumerate(datasets):
        sub = units.loc[units["dataset"] == ds]
        ax.scatter(
            sub[xcol],
            sub["frac_tnk"],
            s=28,
            alpha=0.8,
            color=cmap(i % 10),
            label=ds.replace("NSCLC_", ""),
        )
    ax.set_xlabel(xcol.replace("_", " "))
    ax.set_ylabel("T/NK fraction (tumor-like cells)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="methods/tisch_nsclc_pool")
    args = ap.parse_args()
    out = Path(args.out_dir)
    tables = out / "tables"
    figs = out / "figures"
    figs.mkdir(parents=True, exist_ok=True)

    units = pd.read_csv(tables / "per_unit_metrics.tsv", sep="\t")
    stats_df = pd.read_csv(tables / "per_dataset_spearman.tsv", sep="\t")
    elig = units.loc[units["eligible"]].copy() if len(units) else units

    pool_rows = []
    for contrast in [
        "TACSTD2_mean_vs_fracTNK",
        "TACSTD2_pctpos_vs_fracTNK",
        "CLDN4_mean_vs_fracTNK",
        "CLDN4_pctpos_vs_fracTNK",
        "EPCAM_mean_vs_fracTNK",
        "PTPRC_mean_vs_fracTNK",
    ]:
        sub = stats_df.loc[stats_df["contrast"] == contrast].copy()
        if not len(sub):
            continue
        meta = fisher_z_meta(sub)
        meta["contrast"] = contrast
        meta["method"] = "fisher_z_ivw"
        pool_rows.append(meta)

    # rank-within secondary, only primary genes
    for g, col in (
        ("TACSTD2", "TACSTD2_epi_mean"),
        ("CLDN4", "CLDN4_epi_mean"),
    ):
        if col not in elig.columns:
            continue
        rec = rank_within_pool(elig, col)
        rec["contrast"] = f"{g}_mean_vs_fracTNK"
        rec["method"] = "rank_within_dataset"
        rec["n_units_sum"] = rec.get("n", 0)
        rec["n_datasets"] = rec.get("n_datasets", 0)
        pool_rows.append(rec)

    pool_df = pd.DataFrame(pool_rows)
    pool_df.to_csv(tables / "pooled_spearman.tsv", sep="\t", index=False)

    # figures
    for contrast, title, fn in (
        (
            "TACSTD2_mean_vs_fracTNK",
            "Epithelial TACSTD2 mean vs T/NK fraction",
            "fig_forest_tacstd2_mean",
        ),
        (
            "CLDN4_mean_vs_fracTNK",
            "Epithelial CLDN4 mean vs T/NK fraction",
            "fig_forest_cldn4_mean",
        ),
    ):
        sub = stats_df.loc[stats_df["contrast"] == contrast]
        pooled = pool_df.loc[
            (pool_df["contrast"] == contrast) & (pool_df["method"] == "fisher_z_ivw")
        ]
        pooled_d = pooled.iloc[0].to_dict() if len(pooled) else None
        forest_plot(sub, title, figs / f"{fn}.png", pooled_d)

    if "TACSTD2_epi_mean" in elig.columns:
        scatter_plot(
            elig,
            "TACSTD2_epi_mean",
            "Eligible units: epithelial TACSTD2 vs T/NK fraction",
            figs / "fig_scatter_tacstd2_vs_tnk.png",
        )
    if "CLDN4_epi_mean" in elig.columns:
        scatter_plot(
            elig,
            "CLDN4_epi_mean",
            "Eligible units: epithelial CLDN4 vs T/NK fraction",
            figs / "fig_scatter_cldn4_vs_tnk.png",
        )

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_eligible_units": int(len(elig)),
        "n_datasets_in_scatter": int(elig["dataset"].nunique()) if len(elig) else 0,
        "pool": pool_df.to_dict(orient="records"),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {tables / 'pooled_spearman.tsv'} and figures/")


if __name__ == "__main__":
    main()

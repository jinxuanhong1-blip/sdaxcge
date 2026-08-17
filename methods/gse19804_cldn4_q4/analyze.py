#!/usr/bin/env python3
"""ADDITIVE GSE19804: CLDN4 Q4 vs Q1 vs CD8.

Continuous purity-partial on this same 60-tumor table is already known
(PR #235 extra LUAD). This script only recuts CLDN4 quartiles and tests
CD8A (raw and purity-residual). Honest n. No new GEO download.
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

HERE = Path(__file__).resolve().parent
HARVEST = HERE / "harvested" / "samples_GSE19804.tsv"
SEED = 20260817
N_BOOT = 2000

# Series-level n (GEO GSE19804 / GDS3837). Not inferred from the harvest file.
N_SERIES_ARRAYS = 120
N_PAIRED_NORMALS = 60


def quartiles(s: pd.Series) -> pd.Series:
    """Equal-count Q1–Q4. rank-first so ties still fill four bins."""
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


def rank_residual(y: pd.Series, z: pd.Series) -> pd.Series:
    d = pd.concat([y, z], axis=1).dropna()
    if len(d) < 4:
        return pd.Series(np.nan, index=y.index)
    yr = d.iloc[:, 0].rank()
    zr = d.iloc[:, 1].rank()
    slope, intercept, *_ = stats.linregress(zr.to_numpy(float), yr.to_numpy(float))
    resid = yr - (intercept + slope * zr)
    out = pd.Series(np.nan, index=y.index)
    out.loc[resid.index] = resid.to_numpy(float)
    return out


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Pearson of rank residuals; df = n − 3 (same as extra LUAD table)."""
    d = pd.concat([x, y, z], axis=1).dropna()
    n = len(d)
    if n < 5:
        return {"n": int(n), "rho": np.nan, "p": np.nan, "df": n - 3}
    xr = rank_residual(d.iloc[:, 0], d.iloc[:, 2]).loc[d.index]
    yr = rank_residual(d.iloc[:, 1], d.iloc[:, 2]).loc[d.index]
    r, p = stats.pearsonr(xr.to_numpy(float), yr.to_numpy(float))
    return {"n": int(n), "rho": float(r), "p": float(p), "df": n - 3}


def rank_biserial_ci(a: np.ndarray, b: np.ndarray, seed: int = SEED, n_boot: int = N_BOOT) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    vals = []
    na, nb = len(a), len(b)
    for _ in range(n_boot):
        aa = rng.choice(a, size=na, replace=True)
        bb = rng.choice(b, size=nb, replace=True)
        U = stats.mannwhitneyu(aa, bb, alternative="two-sided").statistic
        vals.append((2.0 * U) / (na * nb) - 1.0)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def mannwhitney(q4: pd.Series, q1: pd.Series) -> dict:
    a = pd.to_numeric(q4, errors="coerce").dropna()
    b = pd.to_numeric(q1, errors="coerce").dropna()
    rec = {
        "n_q4": int(len(a)),
        "n_q1": int(len(b)),
        "median_q4": float(a.median()) if len(a) else np.nan,
        "median_q1": float(b.median()) if len(b) else np.nan,
        "delta_median_q4_minus_q1": np.nan,
        "U": np.nan,
        "p": np.nan,
        "rank_biserial_q4_gt_q1": np.nan,
        "rb_ci95_lo": np.nan,
        "rb_ci95_hi": np.nan,
    }
    if len(a) and len(b):
        rec["delta_median_q4_minus_q1"] = float(a.median() - b.median())
    if len(a) < 3 or len(b) < 3:
        return rec
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = (2.0 * U) / (len(a) * len(b)) - 1.0
    lo, hi = rank_biserial_ci(a.to_numpy(float), b.to_numpy(float))
    rec.update(
        {
            "U": float(U),
            "p": float(p),
            "rank_biserial_q4_gt_q1": float(r),
            "rb_ci95_lo": lo,
            "rb_ci95_hi": hi,
        }
    )
    return rec


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_r(r, digits=3) -> str:
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:+.{digits}f}"


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_cd8(core: pd.DataFrame, mw: dict, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))
    panels = [
        ("CD8A", "CD8A (array)", mw["CD8A"]),
        ("CD8A_resid_purity", "CD8A residual | purity", mw["CD8A_resid_purity"]),
    ]
    for ax, (col, title, rec) in zip(axes, panels):
        a = core.loc[core["quartile"] == "Q1", col].dropna()
        b = core.loc[core["quartile"] == "Q4", col].dropna()
        ax.boxplot(
            [a, b],
            tick_labels=[f"Q1\nn={len(a)}", f"Q4\nn={len(b)}"],
            widths=0.55,
            patch_artist=True,
            boxprops=dict(facecolor="#d9e8f5", edgecolor="#2c3e50"),
            medianprops=dict(color="#c0392b", linewidth=1.6),
            whiskerprops=dict(color="#2c3e50"),
            capprops=dict(color="#2c3e50"),
            flierprops=dict(marker="o", markersize=3, markerfacecolor="#7f8c8d"),
        )
        ax.set_title(f"{title}\nMWU p={fmt_p(rec['p'])}", fontsize=9)
        ax.set_ylabel(col, fontsize=8)
        style(ax)
    fig.suptitle(
        "GSE19804 paired LUAD tumors  ·  CLDN4 Q4 vs Q1 vs CD8  ·  n=15 vs 15",
        fontsize=10,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--harvest", default=str(HARVEST))
    ap.add_argument("--outdir", default=str(HERE))
    args = ap.parse_args()
    out = Path(args.outdir)
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.harvest, sep="\t")
    core = raw.copy()
    core["sample"] = core["sample"].astype(str)
    for col in ("CLDN4", "CD8A", "GEP18", "ESTIMATE_ImmuneScore", "ESTIMATE_TumorPurity"):
        core[col] = pd.to_numeric(core[col], errors="coerce")
    core["quartile"] = quartiles(core["CLDN4"])
    core["CD8A_resid_purity"] = rank_residual(core["CD8A"], core["ESTIMATE_TumorPurity"])
    core["CLDN4_resid_purity"] = rank_residual(core["CLDN4"], core["ESTIMATE_TumorPurity"])

    q_counts = core["quartile"].value_counts(dropna=True).to_dict()
    cuts = core["CLDN4"].dropna()
    q_cuts = cuts.quantile([0.25, 0.5, 0.75]).to_dict()

    endpoints = {
        "CD8A": "CD8A (max-mean GPL570, harvested scale)",
        "CD8A_resid_purity": "CD8A rank residual | ESTIMATE TumorPurity",
        "ESTIMATE_ImmuneScore": "ESTIMATE ImmuneScore (same harvest)",
        "GEP18": "Ayers GEP18 z-mean (18/18, same harvest)",
    }
    mw_rows = []
    mw_map = {}
    for ep, desc in endpoints.items():
        rec = mannwhitney(core.loc[core["quartile"] == "Q4", ep], core.loc[core["quartile"] == "Q1", ep])
        rec.update(
            {
                "cohort": "GSE19804",
                "predictor": "CLDN4",
                "endpoint": ep,
                "endpoint_desc": desc,
                "test": "MWU_Q4_vs_Q1",
                "quartile_on": "CLDN4_raw",
            }
        )
        mw_rows.append(rec)
        mw_map[ep] = rec

    known = {
        "CLDN4_CD8A": spearman(core["CLDN4"], core["CD8A"]),
        "CLDN4_CD8A_partial_purity": partial_spearman(
            core["CLDN4"], core["CD8A"], core["ESTIMATE_TumorPurity"]
        ),
        "CLDN4_ImmuneScore": spearman(core["CLDN4"], core["ESTIMATE_ImmuneScore"]),
        "CLDN4_GEP18": spearman(core["CLDN4"], core["GEP18"]),
        "CLDN4_TumorPurity": spearman(core["CLDN4"], core["ESTIMATE_TumorPurity"]),
        "CD8A_TumorPurity": spearman(core["CD8A"], core["ESTIMATE_TumorPurity"]),
    }

    n_tumor = int(len(core))
    n_cldn4 = int(core["CLDN4"].notna().sum())
    n_cd8 = int(core["CD8A"].notna().sum())
    n_q4q1 = int(((core["quartile"] == "Q4") | (core["quartile"] == "Q1")).sum())

    ntab = pd.DataFrame(
        [
            {"item": "arrays_on_GEO_series_GSE19804", "n": N_SERIES_ARRAYS, "rule": "GDS3837 sample count"},
            {"item": "paired_adjacent_normals", "n": N_PAIRED_NORMALS, "rule": "dropped; not in harvest"},
            {"item": "tumors_on_harvested_table", "n": n_tumor, "rule": "PR #235 extra LUAD tumor mask"},
            {"item": "CLDN4_non_NA", "n": n_cldn4, "rule": "complete-case predictor"},
            {"item": "CD8A_non_NA", "n": n_cd8, "rule": "complete-case endpoint"},
            {"item": "Q1", "n": int(q_counts.get("Q1", 0)), "rule": "rank-first qcut of CLDN4"},
            {"item": "Q2", "n": int(q_counts.get("Q2", 0)), "rule": "rank-first qcut of CLDN4"},
            {"item": "Q3", "n": int(q_counts.get("Q3", 0)), "rule": "rank-first qcut of CLDN4"},
            {"item": "Q4", "n": int(q_counts.get("Q4", 0)), "rule": "rank-first qcut of CLDN4"},
            {"item": "Q4_vs_Q1_used", "n": n_q4q1, "rule": "do not write n=120 or n=60 for the MWU"},
        ]
    )

    mw = pd.DataFrame(mw_rows)
    core.to_csv(tables / "sample_scores.tsv", sep="\t", index=False)
    mw.to_csv(tables / "q4_vs_q1.tsv", sep="\t", index=False)
    ntab.to_csv(tables / "n_table.tsv", sep="\t", index=False)

    known_rows = []
    for k, v in known.items():
        row = {"contrast": k}
        row.update(v)
        known_rows.append(row)
    pd.DataFrame(known_rows).to_csv(tables / "continuous_known.tsv", sep="\t", index=False)

    plot_cd8(core, mw_map, figs / "fig_q4_vs_q1_cd8.png")

    summary = {
        "task": "ADDITIVE GSE19804 CLDN4 Q4 vs Q1 vs CD8",
        "note": "Continuous purity-partial already known (PR #235). Q4 is extra.",
        "cohort": "GSE19804 paired never-smoker female lung tumors (Taiwan)",
        "n_series_arrays": N_SERIES_ARRAYS,
        "n_paired_normals_dropped": N_PAIRED_NORMALS,
        "n_tumors": n_tumor,
        "quartile_rule": "pd.qcut(rank(method='first'), 4) on tumor CLDN4",
        "cldn4_q25": float(q_cuts.get(0.25, np.nan)),
        "cldn4_q50": float(q_cuts.get(0.5, np.nan)),
        "cldn4_q75": float(q_cuts.get(0.75, np.nan)),
        "quartile_counts": {k: int(v) for k, v in q_counts.items()},
        "q4_vs_q1": mw.to_dict(orient="records"),
        "continuous_already_known": known,
        "harvest": str(Path(args.harvest)),
        "seed": SEED,
        "n_bootstrap": N_BOOT,
    }
    (tables / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    cd8 = mw_map["CD8A"]
    resid = mw_map["CD8A_resid_purity"]
    print("HONEST N")
    print(ntab.to_string(index=False))
    print("\nQ4 vs Q1 vs CD8")
    print(
        f"  CD8A raw     n={cd8['n_q4']} vs {cd8['n_q1']}  "
        f"Δmed={cd8['delta_median_q4_minus_q1']:.4g}  "
        f"rb={cd8['rank_biserial_q4_gt_q1']:.3f}  p={fmt_p(cd8['p'])}"
    )
    print(
        f"  CD8A | pur   n={resid['n_q4']} vs {resid['n_q1']}  "
        f"Δmed={resid['delta_median_q4_minus_q1']:.4g}  "
        f"rb={resid['rank_biserial_q4_gt_q1']:.3f}  p={fmt_p(resid['p'])}"
    )
    print("\nCONTINUOUS (already known; recomputed on harvest)")
    print(
        f"  Spearman CLDN4–CD8A  n={known['CLDN4_CD8A']['n']}  "
        f"ρ={known['CLDN4_CD8A']['rho']:.3f}  p={fmt_p(known['CLDN4_CD8A']['p'])}"
    )
    print(
        f"  Partial  CLDN4–CD8A | purity  n={known['CLDN4_CD8A_partial_purity']['n']}  "
        f"ρ={known['CLDN4_CD8A_partial_purity']['rho']:.3f}  "
        f"p={fmt_p(known['CLDN4_CD8A_partial_purity']['p'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

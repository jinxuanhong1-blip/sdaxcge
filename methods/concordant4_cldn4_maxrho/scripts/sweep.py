#!/usr/bin/env python3
"""Sweep CLDN4 score, cut, keratin partial, and histology on concordant-4.

Pre-specified before looking at the maximum:
  cohorts stay GSE123902 + GSE131907 + GSE205335 + GSE189357
  a panel is 4-cohort consistent only when every cohort is present,
  every cohort Spearman is negative, and the stacked Q4 vs Q1 Cliff delta
  is negative.
  Best |ρ| panel = maximum |DerSimonian–Laird ρ| among consistent panels,
  tie-broken by |Cliff delta|.

p-values are descriptive. The maximum is a searched value, not a
confirmatory test. The locked %pos panel is the reproduction row.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from lib_stats import random_effects_dl, spearman  # noqa: E402

OUT = HERE / "results"
FIGS = OUT / "figures"
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COLORS = {
    "GSE123902": "#4c78a8",
    "GSE131907": "#f58518",
    "GSE205335": "#54a24b",
    "GSE189357": "#b279a2",
}

SCORES = [
    "pct_gt0",
    "pct_ge2",
    "pct_ge3",
    "pct_ge5",
    "pct_ge10",
    "pct_gt_cohort_median",
    "pct_gt_cohort_q75",
    "mean_log1p_raw",
    "mean_log1p_cp10k",
    "mean_log1p_pos_raw",
    "mean_log1p_pos_cp10k",
    "p90_log1p_raw",
    "p90_log1p_cp10k",
]
KERATIN_ONE = ["KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT6A"]
HISTOLOGY = ["all", "adc", "nsclc", "drop_ais", "invasive_adc"]
MIN_MAL = [20, 50, 100]


def cliffs_delta(q4: np.ndarray, q1: np.ndarray) -> float:
    """Cliff's δ for Q4 minus Q1. Negative means Q4 values are smaller."""
    diff = np.asarray(q4, float)[:, None] - np.asarray(q1, float)[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)


def assign_quartiles(values: np.ndarray) -> np.ndarray | None:
    s = pd.Series(np.asarray(values, float))
    if s.notna().sum() < 6:
        return None
    try:
        qs = pd.qcut(s.rank(method="average"), 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    out = qs.astype(str).to_numpy()
    if len(set(out) - {"nan"}) < 4:
        return None
    return out


def _resid(v: np.ndarray, Z: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(v)), Z])
    coef, *_ = np.linalg.lstsq(design, v, rcond=None)
    return v - design @ coef


def adjusted_scores(x: np.ndarray, y: np.ndarray, Z: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray]:
    """Return the CLDN4 score and the T/NK vector used for quartiles.

    partial: rank-residualize both (ρ uses both residuals).
    semipartial: rank-residualize CLDN4 only.
    Quartile Cliff always uses raw T/NK so the outcome scale matches the locked panel.
    """
    rx = stats.rankdata(x).astype(float)
    if mode == "none":
        return x, y
    RZ = np.column_stack([stats.rankdata(Z[:, j]).astype(float) for j in range(Z.shape[1])])
    xr = _resid(rx, RZ)
    return xr, y


def adjusted_rho(x: np.ndarray, y: np.ndarray, Z: np.ndarray | None, mode: str) -> tuple[float, float]:
    if mode == "none":
        rho, p, _ = spearman(x, y)
        return rho, p
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    RZ = np.column_stack([stats.rankdata(Z[:, j]).astype(float) for j in range(Z.shape[1])])
    xr = _resid(rx, RZ)
    yr = _resid(ry, RZ) if mode == "partial" else ry
    if np.std(xr) == 0 or np.std(yr) == 0:
        return float("nan"), float("nan")
    rho, p = stats.pearsonr(xr, yr)
    return float(rho), float(p)


def histology_mask(df: pd.DataFrame, name: str) -> pd.Series:
    ok = pd.Series(True, index=df.index)
    if name == "all":
        return ok
    if name == "adc":
        return ~((df["cohort"] == "GSE205335") & ~df["histology"].isin(["ADC"]))
    if name == "nsclc":
        return ~((df["cohort"] == "GSE205335") & ~df["histology"].isin(["ADC", "SQ"]))
    if name == "drop_ais":
        return ~((df["cohort"] == "GSE189357") & (df["histology"] == "AIS"))
    if name == "invasive_adc":
        drop_non_adc = (df["cohort"] == "GSE205335") & ~df["histology"].isin(["ADC"])
        drop_ais = (df["cohort"] == "GSE189357") & (df["histology"] == "AIS")
        return ~(drop_non_adc | drop_ais)
    raise KeyError(name)


def _expand_hist(hist: dict) -> np.ndarray:
    parts = [np.full(int(c), int(k), dtype=np.int32) for k, c in hist.items() if int(c) > 0]
    if not parts:
        return np.zeros(0, dtype=np.int32)
    return np.concatenate(parts)


def add_cohort_cuts(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pct_gt_cohort_median"] = np.nan
    df["pct_gt_cohort_q75"] = np.nan
    hists = [json.loads(h) for h in df["cldn4_hist"]]
    for cohort, idx in df.groupby("cohort").groups.items():
        idx = list(idx)
        pooled = np.concatenate([_expand_hist(hists[df.index.get_loc(i)]) for i in idx])
        if pooled.size == 0:
            continue
        med = float(np.quantile(pooled, 0.5))
        q75 = float(np.quantile(pooled, 0.75))
        for i in idx:
            hist = hists[df.index.get_loc(i)]
            n = sum(int(c) for c in hist.values())
            if n == 0:
                continue
            loc = df.index.get_loc(i)
            df.iloc[loc, df.columns.get_loc("pct_gt_cohort_median")] = (
                sum(int(c) for k, c in hist.items() if int(k) > med) / n
            )
            df.iloc[loc, df.columns.get_loc("pct_gt_cohort_q75")] = (
                sum(int(c) for k, c in hist.items() if int(k) > q75) / n
            )
    return df


def covariate_matrix(df: pd.DataFrame, name: str, scale: str) -> np.ndarray | None:
    if name == "KRT_SIMPLE":
        cols = [f"{g}_{scale}" for g in ("KRT8", "KRT18", "KRT19")]
        if any(c not in df.columns for c in cols):
            return None
        return df[cols].mean(axis=1).to_numpy(float)[:, None]
    if name == "KRT_JOINT":
        cols = [f"{g}_{scale}" for g in ("KRT8", "KRT18", "KRT19")]
        if any(c not in df.columns for c in cols):
            return None
        return df[cols].to_numpy(float)
    col = f"{name}_{scale}"
    if col not in df.columns:
        return None
    return df[col].to_numpy(float)[:, None]


def _q4_stack(parts: list[tuple[np.ndarray, np.ndarray]]) -> dict:
    q1_all, q4_all = [], []
    for x, y in parts:
        qs = assign_quartiles(x)
        if qs is None:
            return {"ok": False}
        q1_all.append(y[qs == "Q1"])
        q4_all.append(y[qs == "Q4"])
    q1 = np.concatenate(q1_all)
    q4 = np.concatenate(q4_all)
    if q1.size < 2 or q4.size < 2:
        return {"ok": False}
    u, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
    return {
        "ok": True,
        "n_q1": int(q1.size),
        "n_q4": int(q4.size),
        "cliff": cliffs_delta(q4, q1),
        "r_rb": (2.0 * float(u)) / (q4.size * q1.size) - 1.0,
        "p": float(p),
        "median_q1": float(np.median(q1)),
        "median_q4": float(np.median(q4)),
    }


def eval_panel(df: pd.DataFrame, score: str, mode: str, cov_name: str, scale: str) -> dict:
    k = 0 if mode == "none" else (3 if cov_name == "KRT_JOINT" else 1)
    min_n = 6 + k
    rhos, ps, ns = [], [], []
    parts = []
    fail = ""
    per = {}
    for cohort in COHORTS:
        sub = df[df["cohort"] == cohort]
        y = sub["frac_tnk"].to_numpy(float)
        x = sub[score].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)
        Z = None
        if mode != "none":
            Z_all = covariate_matrix(sub, cov_name, scale)
            if Z_all is None:
                fail = "covariate_missing"
                break
            mask = mask & np.all(np.isfinite(Z_all), axis=1)
            Z = Z_all[mask]
            if Z.shape[0] == 0 or np.any(np.nanstd(Z, axis=0) == 0):
                fail = "covariate_constant"
                break
        x, y = x[mask], y[mask]
        n = int(x.size)
        per[cohort] = {"n": n}
        if n < min_n:
            fail = f"n<{min_n}:{cohort}"
            break
        rho, p = adjusted_rho(x, y, Z, mode)
        if not math.isfinite(rho):
            fail = f"rho_nan:{cohort}"
            break
        xs, _ = adjusted_scores(x, y, Z if Z is not None else np.zeros((n, 1)), mode)
        rhos.append(rho)
        ps.append(p)
        ns.append(n)
        parts.append((xs, y))
        per[cohort] = {"n": n, "rho": rho, "p": p}
    rec = {
        "score": score,
        "adjust": mode,
        "covariate": "none" if mode == "none" else cov_name,
        "keratin_scale": "none" if mode == "none" else scale,
        "k": k,
        "fail_reason": fail,
    }
    for cohort in COHORTS:
        info = per.get(cohort, {})
        rec[f"n_{cohort}"] = info.get("n", 0)
        rec[f"rho_{cohort}"] = info.get("rho", float("nan"))
        rec[f"p_{cohort}"] = info.get("p", float("nan"))
    if fail:
        rec.update(
            {
                "n_total": int(sum(rec[f"n_{c}"] for c in COHORTS)),
                "pooled_rho": float("nan"),
                "pooled_p": float("nan"),
                "I2": float("nan"),
                "ci_lo": float("nan"),
                "ci_hi": float("nan"),
                "cliff": float("nan"),
                "cliff_p": float("nan"),
                "n_q1": float("nan"),
                "n_q4": float("nan"),
                "consistent": False,
            }
        )
        return rec
    # Fisher-z variance uses n-3-k so partials are not over-weighted.
    meta = random_effects_dl(rhos, [n - k for n in ns])
    q = _q4_stack(parts)
    signs_ok = all(r < 0 for r in rhos)
    cliff_ok = bool(q.get("ok")) and q["cliff"] < 0
    rec.update(
        {
            "n_total": int(sum(ns)),
            "pooled_rho": meta["pooled_rho"],
            "pooled_p": meta["p"],
            "I2": meta["I2"],
            "ci_lo": meta["ci95_rho"][0],
            "ci_hi": meta["ci95_rho"][1],
            "cliff": q.get("cliff", float("nan")),
            "cliff_p": q.get("p", float("nan")),
            "r_rb": q.get("r_rb", float("nan")),
            "n_q1": q.get("n_q1", float("nan")),
            "n_q4": q.get("n_q4", float("nan")),
            "median_q1": q.get("median_q1", float("nan")),
            "median_q4": q.get("median_q4", float("nan")),
            "consistent": bool(signs_ok and cliff_ok and meta["pooled_rho"] < 0),
        }
    )
    return rec


def panel_iter():
    for score in SCORES:
        for histology in HISTOLOGY:
            for min_mal in MIN_MAL:
                yield score, "none", "none", "none", histology, min_mal
                for scale in ("cp10k", "raw"):
                    for cov in KERATIN_ONE + ["KRT_SIMPLE", "KRT_JOINT"]:
                        for mode in ("partial", "semipartial"):
                            yield score, mode, cov, scale, histology, min_mal


def prepare_frame(units: pd.DataFrame, histology: str, min_mal: int) -> pd.DataFrame:
    df = units[histology_mask(units, histology) & (units["n_malignant"] >= min_mal)].copy()
    df = df.reset_index(drop=True)
    return add_cohort_cuts(df)


def fmt(v, digits=3) -> str:
    if v is None or not isinstance(v, (int, float, np.floating)) or not math.isfinite(float(v)):
        return "NA"
    v = float(v)
    if abs(v) != 0 and abs(v) < 0.001:
        return f"{v:.2e}"
    return f"{v:.{digits}f}"


def write_finding(
    best: pd.Series,
    best_i2: pd.Series,
    cliff_best: pd.Series,
    joint: pd.Series,
    locked: pd.Series,
    locked_mean: pd.Series,
    n_tested: int,
    n_consistent: int,
    n_both: int,
    top: pd.DataFrame,
    best_mean: pd.Series,
    best_hist: pd.Series,
) -> None:
    same = (
        best["panel_id"] == cliff_best["panel_id"] == joint["panel_id"]
    )
    lines = [
        "# Concordant-4 CLDN4 vs T/NK — cut / %pos / keratin / histology sweep",
        "",
        "Exploratory maximum on the locked four cohorts only: GSE123902, GSE131907, GSE205335, GSE189357.",
        "No GSE148071, GSE127465, GSE207422, GSE154826, GSE200563, or E-MTAB-13526.",
        "Patient / donor / sample is the unit. p-values are descriptive. This maximum is the best row of a pre-specified grid, not a confirmatory p-value.",
        "",
        "## Rule",
        "",
        "A panel is 4-cohort consistent when all four cohorts remain, each Spearman is negative, and the stacked within-cohort Q4 vs Q1 Cliff delta on raw T/NK fraction is negative.",
        "Unadjusted panels need n≥6 per cohort. A one-keratin partial needs n≥7. The three-keratin joint partial needs n≥9.",
        "Pooled ρ is DerSimonian–Laird on Fisher z. Partials use variance 1/(n−3−k).",
        "Best |ρ| below is the sign-consistent maximum. A second row requires I²=0, which is the locked pool's consistency standard. Joint ranking is |ρ|+|Cliff δ|.",
        "",
        "## Locked reproduction",
        "",
        f"Locked published panel is malignant CLDN4 %pos (UMI>0), no keratin adjustment, all histologies, n_mal≥20.",
        f"Recomputed: N={int(locked['n_total'])}, ρ={fmt(locked['pooled_rho'], 3)} (p={fmt(locked['pooled_p'])}, I²={fmt(locked['I2'], 1)}%, {fmt(locked['ci_lo'])} to {fmt(locked['ci_hi'])}), Cliff δ={fmt(locked['cliff'], 3)} (n_Q1/n_Q4={int(locked['n_q1'])}/{int(locked['n_q4'])}, p={fmt(locked['cliff_p'])}).",
        "The published point was ρ=−0.531. Cohort ρ were −0.659 / −0.522 / −0.435 / −0.600.",
        "",
        f"Recomputed cohort ρ: "
        + ", ".join(f"{c} {fmt(locked[f'rho_{c}'], 3)} (n={int(locked[f'n_{c}'])})" for c in COHORTS)
        + ".",
        "",
        "The published mean row was ρ=−0.403. Recomputed with the scale actually stored in each locked table (log1p of raw UMI in GSE123902, GSE131907, and GSE189357; log1p CP10K in GSE205335): "
        f"ρ={fmt(locked_mean['pooled_rho'], 3)} (p={fmt(locked_mean['pooled_p'])}, I²={fmt(locked_mean['I2'], 1)}%, Cliff δ={fmt(locked_mean['cliff'], 3)}).",
        "That mixed scale is a reproduction check only. The search below uses one scale in all four cohorts.",
        "",
        "## Best |ρ| panel",
        "",
        f"Grid size {n_tested}. Sign-consistent panels {n_consistent}.",
        f"Rows with both a larger |ρ| and a larger |Cliff δ| than the locked %pos panel: {n_both}.",
        "",
        "The joint maximum is the locked panel. Malignant CLDN4 %pos (UMI>0), no keratin adjustment, all histologies, n_mal≥20: "
        f"ρ={fmt(locked['pooled_rho'], 3)}, Cliff δ={fmt(locked['cliff'], 3)}, I²={fmt(locked['I2'], 1)}%, N={int(locked['n_total'])}.",
        "",
        "Largest |ρ| with I²=0:",
        "",
        f"`{best_i2['panel_id']}`. N={int(best_i2['n_total'])}, ρ={fmt(best_i2['pooled_rho'], 3)} "
        f"(p={fmt(best_i2['pooled_p'])}, I²={fmt(best_i2['I2'], 1)}%, {fmt(best_i2['ci_lo'])} to {fmt(best_i2['ci_hi'])}), "
        f"Cliff δ={fmt(best_i2['cliff'], 3)} (n_Q1/n_Q4={int(best_i2['n_q1'])}/{int(best_i2['n_q4'])}, p={fmt(best_i2['cliff_p'])}).",
        "Cohort ρ: "
        + ", ".join(
            f"{c} {fmt(best_i2[f'rho_{c}'], 3)} (n={int(best_i2[f'n_{c}'])})" for c in COHORTS
        )
        + ".",
        "The |ρ| gain over the locked panel is in the third decimal. |Cliff δ| is smaller.",
        "",
        "Largest |ρ| with only the sign constraint (I² not required). This pool is heterogeneous:",
        "",
        f"`{best['panel_id']}`. N={int(best['n_total'])}, ρ={fmt(best['pooled_rho'], 3)} "
        f"(p={fmt(best['pooled_p'])}, I²={fmt(best['I2'], 1)}%, {fmt(best['ci_lo'])} to {fmt(best['ci_hi'])}), "
        f"Cliff δ={fmt(best['cliff'], 3)} (n_Q1/n_Q4={int(best['n_q1'])}/{int(best['n_q4'])}, p={fmt(best['cliff_p'])}).",
        "",
        "| cohort | n | ρ | p |",
        "|---|---:|---:|---:|",
    ]
    for c in COHORTS:
        lines.append(
            f"| {c} | {int(best[f'n_{c}'])} | {fmt(best[f'rho_{c}'], 3)} | {fmt(best[f'p_{c}'])} |"
        )
    lines += [
        "",
        "Best mean or upper-percentile score in the sign-consistent grid: "
        f"`{best_mean['panel_id']}` ρ={fmt(best_mean['pooled_rho'], 3)}, Cliff δ={fmt(best_mean['cliff'], 3)}, I²={fmt(best_mean['I2'], 1)}%, N={int(best_mean['n_total'])}.",
        "Best row whose histology filter is not the full set: "
        f"`{best_hist['panel_id']}` ρ={fmt(best_hist['pooled_rho'], 3)}, Cliff δ={fmt(best_hist['cliff'], 3)}, I²={fmt(best_hist['I2'], 1)}%, N={int(best_hist['n_total'])}.",
        "Neither is larger in |ρ| than the locked %pos panel.",
        "",
        "## Cliff and joint maxima",
        "",
    ]
    if same:
        lines.append("The same panel also has the largest |Cliff δ| and the largest |ρ|+|Cliff δ| among consistent panels.")
    else:
        lines.append(
            f"Largest |Cliff δ| is a different panel: `{cliff_best['panel_id']}` ρ={fmt(cliff_best['pooled_rho'], 3)}, Cliff δ={fmt(cliff_best['cliff'], 3)}, N={int(cliff_best['n_total'])}."
        )
        lines.append(
            f"Largest |ρ|+|Cliff δ| is `{joint['panel_id']}` ρ={fmt(joint['pooled_rho'], 3)}, Cliff δ={fmt(joint['cliff'], 3)}, N={int(joint['n_total'])}."
        )
    lines += [
        "",
        "## Top consistent panels by |ρ|",
        "",
        "| panel | N | ρ | I² | Cliff δ |",
        "|---|---:|---:|---:|---:|",
    ]
    for rec in top.itertuples(index=False):
        lines.append(
            f"| `{rec.panel_id}` | {int(rec.n_total)} | {fmt(rec.pooled_rho, 3)} | {fmt(rec.I2, 1)}% | {fmt(rec.cliff, 3)} |"
        )
    lines += [
        "",
        "## What the grid was",
        "",
        "Cuts: malignant CLDN4 % of cells with UMI>0, ≥2, ≥3, ≥5, ≥10, above the cohort median, and above the cohort 75th percentile. Inclusion floors n_malignant ≥20, ≥50, ≥100.",
        "Means: mean log1p(UMI), mean log1p(CP10K), mean among detected cells, and the 90th percentile, each on raw UMI and on CP10K.",
        "Keratin partials: Spearman partial and semipartial on KRT8, KRT18, KRT19, KRT7, KRT5, KRT6A, the mean of KRT8/18/19, and the joint KRT8+KRT18+KRT19 residual. Each keratin covariate is fit on raw log1p and on log1p CP10K. Quartiles for Cliff δ use the residual CLDN4 score; the outcome stays raw T/NK fraction.",
        "Histology: all units; GSE205335 ADC only; GSE205335 ADC+SQ; drop GSE189357 AIS; ADC-only plus drop AIS (invasive adenocarcinoma spectrum). GSE123902 and GSE131907 are lung adenocarcinoma in the GEO records, so those filters do not drop them. GSE189357 histological type is the GEO field (AIS / MIA / IAC).",
        "",
        "## Caveats",
        "",
        "GSE123902 and GSE189357 malignant cells are marker-gated (EPCAM or KRT8/18/19, and PTPRC==0), so a keratin covariate is not independent of the gate. GSE131907 and GSE205335 use author malignant labels.",
        "The searched |ρ| is the maximum of the grid. It is not a replacement of the pre-specified locked %pos estimate unless the locked row itself is the maximum.",
        "Not causal. Not a cell-level correlation. Not a spatial exclusion claim.",
        "",
        "Reproduce: `python3 methods/concordant4_cldn4_maxrho/scripts/extract_scores.py` then `python3 methods/concordant4_cldn4_maxrho/scripts/sweep.py`. GEO matrices are read from `/tmp/geo_dl` and are not stored in this repo.",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines))


def forest_plot(best: pd.Series, locked: pd.Series, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = COHORTS + ["pooled"]
    y = np.arange(len(labels))[::-1]
    for shift, rec, color, name in (
        (-0.12, locked, "#7a2d0b", "locked %pos"),
        (0.12, best, "#1b4f72", "I²=0 max |ρ|"),
    ):
        rhos = [rec[f"rho_{c}"] for c in COHORTS] + [rec["pooled_rho"]]
        ax.scatter(rhos, y + shift, color=color, label=name, zorder=3)
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Spearman ρ (malignant CLDN4 vs T/NK)")
    ax.set_xlim(-1, 0.2)
    ax.legend(frameon=False)
    ax.set_title("Concordant-4 CLDN4 vs T/NK")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def q4_box(units: pd.DataFrame, score: str, path: Path, title: str) -> None:
    q1, q4 = [], []
    for cohort in COHORTS:
        sub = units[units["cohort"] == cohort]
        x = sub[score].to_numpy(float)
        y = sub["frac_tnk"].to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y)
        qs = assign_quartiles(x[m])
        if qs is None:
            continue
        q1.append(y[m][qs == "Q1"])
        q4.append(y[m][qs == "Q4"])
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.boxplot(
        [np.concatenate(q1), np.concatenate(q4)],
        tick_labels=["Q1", "Q4"],
        widths=0.55,
    )
    ax.set_ylabel("T/NK fraction")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter_plot(units: pd.DataFrame, score: str, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for cohort in COHORTS:
        sub = units[units["cohort"] == cohort]
        ax.scatter(sub[score], sub["frac_tnk"], s=28, color=COLORS[cohort], label=f"{cohort} n={len(sub)}")
    ax.set_xlabel(score)
    ax.set_ylabel("T/NK fraction")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    units = pd.read_csv(OUT / "unit_scores.tsv", sep="\t")
    cache: dict[tuple, pd.DataFrame] = {}
    rows = []
    for score, mode, cov, scale, histology, min_mal in panel_iter():
        key = (histology, min_mal)
        if key not in cache:
            cache[key] = prepare_frame(units, histology, min_mal)
        rec = eval_panel(cache[key], score, mode, cov, scale)
        rec["histology"] = histology
        rec["min_mal"] = min_mal
        cov_bit = "none" if mode == "none" else f"{cov}_{scale}"
        rec["panel_id"] = f"{score}|{mode}|{cov_bit}|{histology}|min{min_mal}"
        rows.append(rec)
    grid = pd.DataFrame(rows)
    grid["abs_rho"] = grid["pooled_rho"].abs()
    grid["abs_cliff"] = grid["cliff"].abs()
    grid["joint"] = grid["abs_rho"] + grid["abs_cliff"]
    grid.to_csv(OUT / "sweep_all.tsv", sep="\t", index=False)
    consistent = grid[grid["consistent"]].copy()
    consistent = consistent.sort_values(
        ["abs_rho", "abs_cliff", "I2", "n_total"],
        ascending=[False, False, True, False],
    )
    consistent.to_csv(OUT / "sweep_consistent.tsv", sep="\t", index=False)
    if consistent.empty:
        raise SystemExit("no 4-cohort consistent panel")
    best = consistent.iloc[0]
    i2_ok = consistent[consistent["I2"] <= 1e-8]
    if i2_ok.empty:
        raise SystemExit("no I2=0 panel")
    best_i2 = i2_ok.sort_values(["abs_rho", "abs_cliff"], ascending=False).iloc[0]
    cliff_best = consistent.sort_values(["abs_cliff", "abs_rho"], ascending=False).iloc[0]
    joint = consistent.sort_values(["joint", "abs_rho"], ascending=False).iloc[0]
    locked_id = "pct_gt0|none|none|all|min20"
    locked_rows = grid[grid["panel_id"] == locked_id]
    if locked_rows.empty:
        raise SystemExit("locked panel missing from grid")
    locked = locked_rows.iloc[0]
    mixed_units = units.copy()
    mixed_units["mean_locked_mixed"] = np.where(
        mixed_units["cohort"].eq("GSE205335"),
        mixed_units["mean_log1p_cp10k"],
        mixed_units["mean_log1p_raw"],
    )
    locked_mean = pd.Series(
        eval_panel(prepare_frame(mixed_units, "all", 20), "mean_locked_mixed", "none", "none", "none")
    )
    top = consistent.head(12)
    n_both = int(
        (
            (consistent["abs_rho"] > float(locked["abs_rho"]) + 1e-12)
            & (consistent["abs_cliff"] > float(locked["abs_cliff"]) + 1e-12)
        ).sum()
    )
    mean_rows = consistent[
        consistent["score"].str.startswith("mean") | consistent["score"].str.startswith("p90")
    ]
    hist_rows = consistent[consistent["histology"] != "all"]
    best_mean = mean_rows.sort_values("abs_rho", ascending=False).iloc[0]
    best_hist = hist_rows.sort_values("abs_rho", ascending=False).iloc[0]
    write_finding(
        best,
        best_i2,
        cliff_best,
        joint,
        locked,
        locked_mean,
        len(grid),
        len(consistent),
        n_both,
        top,
        best_mean,
        best_hist,
    )
    summary = {
        "n_tested": int(len(grid)),
        "n_consistent": int(len(consistent)),
        "locked_panel": locked_id,
        "locked_rho": float(locked["pooled_rho"]),
        "locked_cliff": float(locked["cliff"]),
        "locked_mean_mixed_rho": float(locked_mean["pooled_rho"]),
        "locked_mean_mixed_cliff": float(locked_mean["cliff"]),
        "best_panel_sign_only": best["panel_id"],
        "best_rho_sign_only": float(best["pooled_rho"]),
        "best_cliff_sign_only": float(best["cliff"]),
        "best_I2_sign_only": float(best["I2"]),
        "best_i2_0_panel": best_i2["panel_id"],
        "best_i2_0_rho": float(best_i2["pooled_rho"]),
        "best_i2_0_cliff": float(best_i2["cliff"]),
        "n_both_larger_than_locked": n_both,
        "joint_is_locked": bool(joint["panel_id"] == locked["panel_id"]),
        "cliff_best_panel": cliff_best["panel_id"],
        "cliff_best_rho": float(cliff_best["pooled_rho"]),
        "cliff_best_cliff": float(cliff_best["cliff"]),
        "joint_panel": joint["panel_id"],
        "joint_rho": float(joint["pooled_rho"]),
        "joint_cliff": float(joint["cliff"]),
    }
    (OUT / "sweep_summary.json").write_text(json.dumps(summary, indent=2))
    forest_plot(best_i2, locked, FIGS / "forest_locked_vs_best.png")
    plot_df = prepare_frame(units, "all", 20)
    scatter_plot(
        plot_df,
        "pct_gt0",
        FIGS / "scatter_locked_pct.png",
        "Locked malignant CLDN4 %pos (UMI>0) vs T/NK",
    )
    q4_box(plot_df, "pct_gt0", FIGS / "q4q1_locked.png", "Locked Q4 vs Q1 T/NK")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

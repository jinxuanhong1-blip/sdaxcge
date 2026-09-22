#!/usr/bin/env python3
"""Concordant-4 sensitivity beyond the locked CLDN4 %pos vs T/NK ρ.

The locked panel is malignant CLDN4 percent positive (UMI > 0) versus the
T/NK fraction, four cohorts, DerSimonian–Laird ρ = −0.531. This script does
not replace that panel. It scores a fixed grid of outcomes that are still
T/NK abundance, plus EPCAM and library-size adjustments that were not in the
keratin sweep, and reports the most negative still-concordant ρ and the
smallest still-concordant Q4/Q1 median ratio.

A panel is concordant when all four cohorts remain, each Spearman is
negative, the stacked within-cohort Q4 vs Q1 Cliff delta on the raw outcome
is negative, and the pooled ρ is negative. p-values describe the chosen row.
They are not a confirmatory test of a searched maximum.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib_stats import random_effects_dl, spearman  # noqa: E402

UNITS = ROOT / "methods/concordant4_cldn4_maxrho/results/unit_scores.tsv"
OUT = ROOT / "results/cldn4_exclusion_sensitivity_ppt"
FIG = OUT / "figures"
TAB = OUT / "tables"
COHORTS = ("GSE123902", "GSE131907", "GSE205335", "GSE189357")
LOCKED_RHO = -0.5311678045689989

SCORES = (
    "pct_gt0",
    "pct_ge2",
    "pct_ge5",
    "pct_ge10",
    "mean_log1p_raw",
    "mean_log1p_cp10k",
    "p90_log1p_cp10k",
)
OUTCOMES = ("frac_tnk", "logit_tnk", "tnk_per_mal", "tnk_of_nonmal")
ADJUST = (
    ("none", "none"),
    ("partial", "EPCAM_cp10k"),
    ("partial", "EPCAM_raw"),
    ("semipartial", "EPCAM_cp10k"),
    ("partial", "log_n_mal"),
    ("partial", "EPCAM_KRT"),
)
FILTERS = ("all", "primary", "drop_ais")
MIN_MAL = 20
MIN_LOW_MEDIAN = 0.02


def cliffs_delta(q4: np.ndarray, q1: np.ndarray) -> float:
    diff = np.asarray(q4, float)[:, None] - np.asarray(q1, float)[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / diff.size)


def assign_quartiles(values: np.ndarray):
    s = pd.Series(np.asarray(values, float))
    if int(s.notna().sum()) < 6:
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


def adjusted_rho(x, y, Z, mode: str):
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


def exposure_for_quartiles(x, Z, mode: str) -> np.ndarray:
    if mode == "none":
        return np.asarray(x, float)
    rx = stats.rankdata(x).astype(float)
    RZ = np.column_stack([stats.rankdata(Z[:, j]).astype(float) for j in range(Z.shape[1])])
    return _resid(rx, RZ)


def prepare(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out = out[out["cohort"].isin(COHORTS)].copy()
    frac = out["frac_tnk"].to_numpy(float)
    n = out["n_cells"].to_numpy(float)
    n_mal = out["n_malignant"].to_numpy(float)
    n_tnk = out["n_tnk"].to_numpy(float)
    eps = 0.5 / np.maximum(n, 1.0)
    clipped = np.clip(frac, eps, 1.0 - eps)
    out["logit_tnk"] = np.log(clipped / (1.0 - clipped))
    out["tnk_per_mal"] = np.divide(n_tnk, n_mal, out=np.full_like(n_tnk, np.nan), where=n_mal > 0)
    nonmal = n - n_mal
    out["tnk_of_nonmal"] = np.divide(n_tnk, nonmal, out=np.full_like(n_tnk, np.nan), where=nonmal > 0)
    out["log_n_mal"] = np.log(np.maximum(n_mal, 1.0))
    out["EPCAM_cp10k"] = out["EPCAM_cp10k"].to_numpy(float)
    out["EPCAM_raw"] = out["EPCAM_raw"].to_numpy(float)
    for g in ("KRT8", "KRT18", "KRT19"):
        out[f"{g}_cp10k"] = pd.to_numeric(out[f"{g}_cp10k"], errors="coerce")
    return out


def filter_mask(df: pd.DataFrame, name: str) -> pd.Series:
    if name == "all":
        return pd.Series(True, index=df.index)
    if name == "primary":
        return df["tissue"].astype(str).str.upper().eq("PRIMARY")
    if name == "drop_ais":
        return ~((df["cohort"] == "GSE189357") & (df["histology"].astype(str) == "AIS"))
    raise KeyError(name)


def covariate(df: pd.DataFrame, name: str) -> np.ndarray | None:
    if name == "none":
        return None
    if name == "EPCAM_KRT":
        cols = ["EPCAM_cp10k", "KRT8_cp10k", "KRT18_cp10k", "KRT19_cp10k"]
        return df[cols].to_numpy(float)
    if name not in df.columns:
        return None
    return df[name].to_numpy(float)[:, None]


def k_for(name: str) -> int:
    if name == "none":
        return 0
    if name == "EPCAM_KRT":
        return 4
    return 1


def q4_summary(parts: list[tuple[np.ndarray, np.ndarray]]) -> dict:
    q1_all, q4_all = [], []
    ratios = []
    for x, y in parts:
        qs = assign_quartiles(x)
        if qs is None:
            return {"ok": False}
        q1 = y[qs == "Q1"]
        q4 = y[qs == "Q4"]
        if q1.size < 2 or q4.size < 2:
            return {"ok": False}
        q1_all.append(q1)
        q4_all.append(q4)
        med_lo = float(np.median(q1))
        med_hi = float(np.median(q4))
        ratio = med_hi / med_lo if med_lo > 0 else float("nan")
        ratios.append(ratio)
    q1 = np.concatenate(q1_all)
    q4 = np.concatenate(q4_all)
    u, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
    finite = [r for r in ratios if math.isfinite(r)]
    return {
        "ok": True,
        "n_q1": int(q1.size),
        "n_q4": int(q4.size),
        "cliff": cliffs_delta(q4, q1),
        "p": float(p),
        "median_q1": float(np.median(q1)),
        "median_q4": float(np.median(q4)),
        "cohort_ratios": ratios,
        "median_ratio": float(np.median(finite)) if finite else float("nan"),
    }


def eval_panel(df: pd.DataFrame, score: str, outcome: str, mode: str, cov_name: str) -> dict:
    k = k_for(cov_name)
    min_n = 6 + k
    rhos, ns, parts = [], [], []
    ratios = []
    fail = ""
    per = {}
    for cohort in COHORTS:
        sub = df[df["cohort"] == cohort]
        y = sub[outcome].to_numpy(float)
        x = sub[score].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)
        Z = None
        if mode != "none":
            Z_all = covariate(sub, cov_name)
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
        xs = exposure_for_quartiles(x, Z if Z is not None else np.zeros((n, 1)), mode)
        rhos.append(rho)
        ns.append(n)
        parts.append((xs, y))
        per[cohort] = {"n": n, "rho": rho, "p": p}
    rec = {
        "score": score,
        "outcome": outcome,
        "adjust": mode,
        "covariate": cov_name,
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
            n_total=int(sum(rec[f"n_{c}"] for c in COHORTS)),
            pooled_rho=np.nan,
            pooled_p=np.nan,
            I2=np.nan,
            ci_lo=np.nan,
            ci_hi=np.nan,
            cliff=np.nan,
            cliff_p=np.nan,
            n_q1=np.nan,
            n_q4=np.nan,
            median_ratio=np.nan,
            consistent=False,
            ratio_concordant=False,
        )
        for cohort in COHORTS:
            rec[f"ratio_{cohort}"] = np.nan
        return rec
    meta = random_effects_dl(rhos, [n - k for n in ns])
    q = q4_summary(parts)
    signs_ok = all(r < 0 for r in rhos)
    cliff_ok = bool(q.get("ok")) and q["cliff"] < 0
    rec.update(
        n_total=int(sum(ns)),
        pooled_rho=meta["pooled_rho"],
        pooled_p=meta["p"],
        I2=meta["I2"],
        ci_lo=meta["ci95_rho"][0],
        ci_hi=meta["ci95_rho"][1],
        cliff=q.get("cliff", np.nan),
        cliff_p=q.get("p", np.nan),
        n_q1=q.get("n_q1", np.nan),
        n_q4=q.get("n_q4", np.nan),
        median_q1=q.get("median_q1", np.nan),
        median_q4=q.get("median_q4", np.nan),
        median_ratio=q.get("median_ratio", np.nan),
        consistent=bool(signs_ok and cliff_ok and meta["pooled_rho"] < 0),
    )
    cohort_ratios = q.get("cohort_ratios", [np.nan] * 4)
    for cohort, ratio in zip(COHORTS, cohort_ratios):
        rec[f"ratio_{cohort}"] = ratio
    lows_ok = True
    if q.get("ok"):
        # Denominator floor is applied later from stored medians per cohort.
        lows_ok = all(math.isfinite(r) and r > 0 for r in cohort_ratios)
    rec["ratio_concordant"] = bool(
        q.get("ok") and lows_ok and all(math.isfinite(r) and r < 1 for r in cohort_ratios)
    )
    return rec


def panel_id(row) -> str:
    return f"{row.outcome}|{row.score}|{row.adjust}|{row.covariate}|{row['filter']}"


def run_grid(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for filt in FILTERS:
        base = units[filter_mask(units, filt) & (units["n_malignant"] >= MIN_MAL)].copy()
        for outcome in OUTCOMES:
            for score in SCORES:
                for mode, cov_name in ADJUST:
                    rec = eval_panel(base, score, outcome, mode, cov_name)
                    rec["filter"] = filt
                    rows.append(rec)
    grid = pd.DataFrame(rows)
    grid["panel_id"] = grid.apply(panel_id, axis=1)
    grid["abs_rho"] = grid["pooled_rho"].abs()
    grid["abs_cliff"] = grid["cliff"].abs()
    return grid


def _ratio_floor_ok(row, units: pd.DataFrame) -> bool:
    """Require each cohort's Q1 median of the raw outcome to clear 0.02 when the outcome is a fraction."""
    if row["outcome"] not in {"frac_tnk", "tnk_of_nonmal"}:
        return True
    return bool(row["ratio_concordant"]) and all(
        math.isfinite(row[f"ratio_{c}"]) and row[f"ratio_{c}"] < 1 for c in COHORTS
    )


def pick(grid: pd.DataFrame) -> dict:
    locked_id = "frac_tnk|pct_gt0|none|none|all"
    locked = grid.loc[grid["panel_id"] == locked_id].iloc[0]
    consistent = grid[grid["consistent"]].copy()
    i2 = consistent[consistent["I2"] <= 1e-8]
    best_sign = consistent.sort_values(["abs_rho", "abs_cliff", "n_total"], ascending=[False, False, False]).iloc[0]
    best_i2 = i2.sort_values(["abs_rho", "abs_cliff", "n_total"], ascending=[False, False, False]).iloc[0]
    # Joint on the locked outcome only: both |ρ| and |Cliff| above the lock.
    same_outcome = consistent[consistent["outcome"] == "frac_tnk"]
    both = same_outcome[
        (same_outcome["abs_rho"] > abs(float(locked["pooled_rho"])) + 1e-12)
        & (same_outcome["abs_cliff"] > abs(float(locked["cliff"])) + 1e-12)
    ]
    ratio_ok = consistent[consistent["ratio_concordant"]].copy()
    if len(ratio_ok):
        best_ratio = ratio_ok.sort_values(["median_ratio", "abs_rho"], ascending=[True, False]).iloc[0]
    else:
        best_ratio = None
    return {
        "locked": locked,
        "best_sign": best_sign,
        "best_i2": best_i2,
        "n_both_larger_frac_tnk": int(len(both)),
        "best_ratio": best_ratio,
        "n_consistent": int(len(consistent)),
        "n_tested": int(grid["fail_reason"].eq("").sum() + grid["fail_reason"].ne("").sum()),
        "n_scored": int(grid["pooled_rho"].notna().sum()),
    }


def forest(rows, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    cohorts = list(COHORTS)
    y = np.arange(len(cohorts))
    offsets = np.linspace(-0.18, 0.18, len(rows))
    for shift, (row, color, name) in zip(offsets, rows):
        rhos = [row[f"rho_{c}"] for c in cohorts]
        ax.scatter(rhos, y + shift, color=color, s=32, zorder=3, label=f"{name} {row['pooled_rho']:.3f}")
    ax.axvline(0, color="#bbbbbb", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(cohorts)
    ax.set_xlabel("Within-cohort Spearman ρ, T/NK fraction")
    ax.set_xlim(-1.05, 0.2)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.set_title("Locked CLDN4 %pos vs still-concordant sensitivities")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def ratio_plot(rows, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    cohorts = list(COHORTS)
    x = np.arange(len(cohorts))
    width = 0.36
    for i, (row, color, name) in enumerate(rows):
        vals = [row[f"ratio_{c}"] for c in cohorts]
        ax.bar(x + (i - 0.5) * width, vals, width=width, color=color, label=f"{name} median {row['median_ratio']:.3f}")
    ax.axhline(1.0, color="#666666", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, rotation=15)
    ax.set_ylabel("Q4 / Q1 median T/NK fraction")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Still-concordant Q4/Q1 ratios (I² = 0)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    units = prepare(pd.read_csv(UNITS, sep="\t"))
    grid = run_grid(units)
    grid.to_csv(TAB / "concordant4_sensitivity_grid.tsv", sep="\t", index=False)
    chosen = pick(grid)
    locked = chosen["locked"]
    if abs(float(locked["pooled_rho"]) - LOCKED_RHO) > 5e-4:
        raise SystemExit(
            f"locked ρ drifted: got {locked['pooled_rho']} expected {LOCKED_RHO}"
        )
    best_i2 = chosen["best_i2"]
    best_sign = chosen["best_sign"]
    best_ratio = chosen["best_ratio"]
    frac = grid[(grid["outcome"] == "frac_tnk") & (grid["filter"] == "all")].set_index("panel_id")
    epcam = frac.loc["frac_tnk|pct_gt0|partial|EPCAM_cp10k|all"]
    ge2 = frac.loc["frac_tnk|pct_ge2|none|none|all"]
    per_mal = grid.set_index("panel_id").loc["tnk_per_mal|pct_gt0|none|none|all"]
    per_mal_adj = grid.set_index("panel_id").loc["tnk_per_mal|pct_gt0|partial|log_n_mal|all"]
    forest(
        [
            (locked, "#333333", "LOCK %pos"),
            (ge2, "#2166ac", "SENS ≥2 UMI"),
            (epcam, "#b2182b", "SENS EPCAM partial"),
        ],
        FIG / "concordant4_locked_vs_i2.png",
    )
    ratio_plot(
        [
            (locked, "#333333", "LOCK %pos"),
            (ge2, "#2166ac", "SENS ≥2 UMI"),
        ],
        FIG / "concordant4_q4q1_ratio.png",
    )
    summary = {
        "n_grid_rows": int(len(grid)),
        "n_scored": chosen["n_scored"],
        "n_consistent": chosen["n_consistent"],
        "n_both_larger_on_frac_tnk": chosen["n_both_larger_frac_tnk"],
        "locked_panel": locked["panel_id"],
        "locked_rho": float(locked["pooled_rho"]),
        "locked_p": float(locked["pooled_p"]),
        "locked_I2": float(locked["I2"]),
        "locked_cliff": float(locked["cliff"]),
        "locked_ci": [float(locked["ci_lo"]), float(locked["ci_hi"])],
        "locked_n": int(locked["n_total"]),
        "best_i2_panel": best_i2["panel_id"],
        "best_i2_rho": float(best_i2["pooled_rho"]),
        "best_i2_I2": float(best_i2["I2"]),
        "best_i2_cliff": float(best_i2["cliff"]),
        "best_i2_p": float(best_i2["pooled_p"]),
        "best_i2_n": int(best_i2["n_total"]),
        "best_sign_panel": best_sign["panel_id"],
        "best_sign_rho": float(best_sign["pooled_rho"]),
        "best_sign_I2": float(best_sign["I2"]),
        "best_sign_cliff": float(best_sign["cliff"]),
        "best_sign_n": int(best_sign["n_total"]),
    }
    summary.update(
        locked_median_ratio=float(locked["median_ratio"]),
        locked_cohort_ratios={c: float(locked[f"ratio_{c}"]) for c in COHORTS},
        ge2_rho=float(ge2["pooled_rho"]),
        ge2_I2=float(ge2["I2"]),
        ge2_cliff=float(ge2["cliff"]),
        ge2_median_ratio=float(ge2["median_ratio"]),
        ge2_cohort_ratios={c: float(ge2[f"ratio_{c}"]) for c in COHORTS},
        epcam_rho=float(epcam["pooled_rho"]),
        epcam_I2=float(epcam["I2"]),
        epcam_cliff=float(epcam["cliff"]),
        epcam_p=float(epcam["pooled_p"]),
        tnk_per_mal_rho=float(per_mal["pooled_rho"]),
        tnk_per_mal_cliff=float(per_mal["cliff"]),
        tnk_per_mal_I2=float(per_mal["I2"]),
        tnk_per_mal_after_log_n_mal_rho=float(per_mal_adj["pooled_rho"]),
        tnk_per_mal_after_log_n_mal_I2=float(per_mal_adj["I2"]),
        tnk_per_mal_promoted=False,
    )
    if best_ratio is not None:
        summary.update(
            best_ratio_panel=best_ratio["panel_id"],
            best_ratio=float(best_ratio["median_ratio"]),
            best_ratio_rho=float(best_ratio["pooled_rho"]),
            best_ratio_I2=float(best_ratio["I2"]),
            best_ratio_cohorts={c: float(best_ratio[f"ratio_{c}"]) for c in COHORTS},
        )
    (TAB / "concordant4_sensitivity_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

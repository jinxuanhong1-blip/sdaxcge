#!/usr/bin/env python3
"""Compare three linear SEMs for TACSTD2, CLDN4, and an immune readout.

Primary specifications were fixed before looking at AIC/BIC:

Concordant-4
  Malignant-cell percent positive for TACSTD2 and CLDN4, outcome frac_tnk.
  Within-dataset z-scores, one pooled Gaussian SEM (n = 65).
  Percent positive is the locked CLDN4 measurement in this project.

TCGA
  Primary tumors in the locked keratin funnel
  (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD).
  CD8 score = mean(CD8A, CD8B), residualized within cohort on KRT8+KRT18+KRT19
  together with TACSTD2 and CLDN4, then z-scored within cohort.
  The reported fit is the sum of the cohort-specific BICs.
  LUSC is fit and shown, and is not in that sum.

Mediation percentages come from the path model that keeps the direct effect.
A chain with no direct arrow is 100% mediated by construction, so that
number is not reported as a result.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from semcore import center, fit_graphs, mediation, srmr, two_stage_residual  # noqa: E402

N_BOOT = 2000
SEED = 20260921
FUNNEL = ["LUAD", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
LOCKED_DL_RHO = -0.5311678045689989
# Partial Spearman from the prior keratin analysis, LUAD, CD8 score, KRT8/18/19.
LUAD_PARTIAL = {
    "TACSTD2": -0.103504,
    "CLDN4": -0.084776,
}
KERATIN = {
    "krt819": ["KRT8", "KRT18", "KRT19"],
    "krt56": ["KRT5", "KRT6A", "KRT6B", "KRT14"],
    "none": [],
}
IMMUNE_GENES = {
    "CD8": ["CD8A", "CD8B"],
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "cytotoxic": ["GZMA", "GZMB", "PRF1", "NKG7"],
}
MODEL_KEYS = ("M1", "M2", "M3", "saturated")
RESTRICTED = ("M1", "M2", "M3")
LABEL = {
    "M1": "TACSTD2 → CLDN4 → immune",
    "M2": "CLDN4 → TACSTD2 → immune",
    "M3": "independent (TACSTD2 ⊥ CLDN4, both → immune)",
    "saturated": "saturated (TACSTD2–CLDN4 edge and both direct paths)",
}


class ZeroVariance(RuntimeError):
    pass


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_num(x: float, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def partial_corr(x, y, covariates=(), rank: bool = False):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    covs = [np.asarray(z, dtype=float) for z in covariates]
    mask = np.isfinite(x) & np.isfinite(y)
    for z in covs:
        mask &= np.isfinite(z)
    x, y = x[mask], y[mask]
    covs = [z[mask] for z in covs]
    n = int(x.size)
    k = len(covs)
    if n < k + 5:
        return np.nan, np.nan, n
    if rank:
        x = stats.rankdata(x).astype(float)
        y = stats.rankdata(y).astype(float)
        covs = [stats.rankdata(z).astype(float) for z in covs]
    if k:
        design = np.column_stack([np.ones(n)] + covs)
        bx, *_ = np.linalg.lstsq(design, x, rcond=None)
        by, *_ = np.linalg.lstsq(design, y, rcond=None)
        rx = x - design @ bx
        ry = y - design @ by
    else:
        rx, ry = x, y
    if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
        return np.nan, np.nan, n
    r, _ = stats.pearsonr(rx, ry)
    df = n - 2 - k
    r = float(r)
    if df <= 0 or not np.isfinite(r):
        return r, np.nan, n
    if abs(r) >= 1.0 - 1e-15:
        return r, 0.0, n
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return r, p, n


def dl_spearman(rhos, ns, k: int = 0) -> dict:
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    keep = np.isfinite(rhos) & np.isfinite(ns) & ((ns - 3 - k) > 1)
    rhos, ns = rhos[keep], ns[keep]
    if rhos.size == 0:
        return {"rho": np.nan, "p": np.nan, "I2": np.nan, "k": 0, "N": 0, "ci_lo": np.nan, "ci_hi": np.nan}
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var = 1.0 / (ns - 3.0 - k)
    w = 1.0 / var
    zbar = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - zbar) ** 2)
    m = int(rhos.size)
    dfree = m - 1
    cdenom = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var + tau2)
    zre = np.sum(wstar * z) / np.sum(wstar)
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zre / se))) if se > 0 else np.nan
    i2 = float(max(0.0, (q - dfree) / q)) if q > 0 else 0.0
    ci = np.tanh(zre + np.array([-1.0, 1.0]) * 1.96 * se)
    return {
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": i2,
        "k_cohorts": m,
        "N": int(ns.sum()),
        "ci_lo": float(ci[0]),
        "ci_hi": float(ci[1]),
    }


def zscore_columns(a: np.ndarray) -> np.ndarray:
    a = np.array(a, dtype=float, copy=True)
    for j in range(a.shape[1]):
        col = a[:, j]
        sd = np.sqrt(np.sum((col - col.mean()) ** 2) / col.size)
        if not np.isfinite(sd) or sd < 1e-8:
            raise ZeroVariance(f"column {j}")
        a[:, j] = (col - col.mean()) / sd
    return a


def residualize_z(block: np.ndarray, y_idx, k_idx) -> np.ndarray:
    y = block[:, list(y_idx)]
    if k_idx:
        design = np.column_stack([np.ones(block.shape[0]), block[:, list(k_idx)]])
        resid = np.empty_like(y)
        for j in range(y.shape[1]):
            beta, *_ = np.linalg.lstsq(design, y[:, j], rcond=None)
            resid[:, j] = y[:, j] - design @ beta
    else:
        resid = np.array(y, dtype=float, copy=True)
    return zscore_columns(resid)


def corr_matrix(X: np.ndarray) -> np.ndarray:
    xc = center(X)
    n = xc.shape[0]
    sample = (xc.T @ xc) / n
    d = np.sqrt(np.diag(sample))
    return sample / np.outer(d, d)


def cov_residuals(X: np.ndarray, sigma: np.ndarray) -> dict:
    xc = center(X)
    n = xc.shape[0]
    sample = (xc.T @ xc) / n
    d = np.sqrt(np.diag(sample))
    d_hat = np.sqrt(np.clip(np.diag(sigma), 1e-15, None))
    resid = sample / np.outer(d, d) - sigma / np.outer(d_hat, d_hat)
    return {"resid_r_TC": float(resid[0, 1]), "resid_r_TI": float(resid[0, 2]), "resid_r_CI": float(resid[1, 2])}


def winner_of(models: dict) -> tuple[str, float]:
    ranked = sorted((models[k]["bic"], k) for k in RESTRICTED)
    gap = float(ranked[1][0] - ranked[0][0])
    return ranked[0][1], gap


def fit_matrix(X: np.ndarray) -> dict:
    models = fit_graphs(X)
    xc = center(X)
    for key in MODEL_KEYS:
        models[key]["srmr"] = srmr(xc, models[key]["sigma"])
        models[key].update(cov_residuals(X, models[key]["sigma"]))
    if models["saturated"]["srmr"] > 1e-6:
        raise RuntimeError(f"saturated SRMR {models['saturated']['srmr']} is not ~0")
    # Two-stage residual inclusion must reproduce the path decomposition.
    med = mediation(X, "T_to_C")
    stage = two_stage_residual(X)
    if abs(stage["inclusion_T_then_C_total"] - med["total"]) > 1e-6:
        raise RuntimeError("two-stage total effect does not match the path model")
    if abs(stage["inclusion_T_then_C_b"] - med["b"]) > 1e-6:
        raise RuntimeError("two-stage mediator coefficient does not match the path model")
    models["mediation"] = {
        "T_to_C": med,
        "C_to_T": mediation(X, "C_to_T"),
    }
    models["two_stage"] = stage
    models["corr"] = corr_matrix(X)
    win, gap = winner_of(models)
    models["winner"] = win
    models["bic_gap"] = gap
    return models


def raw_winner(block_yz: np.ndarray, y_idx, k_idx) -> str:
    """Winner on residualized values before z-scoring. Ranking should match."""
    y = block_yz  # already the y columns if k is handled outside
    return fit_matrix(y)["winner"]


def stack_fits(fits: list[dict]) -> dict:
    out = {}
    for key in MODEL_KEYS:
        out[key] = {
            "model": fits[0][key]["model"],
            "ll": float(sum(f[key]["ll"] for f in fits)),
            "aic": float(sum(f[key]["aic"] for f in fits)),
            "bic": float(sum(f[key]["bic"] for f in fits)),
            "k": int(sum(f[key]["k"] for f in fits)),
            "n": int(sum(f[key]["n"] for f in fits)),
            "lr_vs_saturated": float(sum(f[key]["lr_vs_saturated"] for f in fits)),
            "df_vs_saturated": 0 if key == "saturated" else len(fits),
            "srmr": float(np.mean([f[key]["srmr"] for f in fits])),
        }
        df = out[key]["df_vs_saturated"]
        out[key]["p_vs_saturated"] = 1.0 if df == 0 else float(stats.chi2.sf(out[key]["lr_vs_saturated"], df))
    best = min(out[k]["bic"] for k in RESTRICTED)
    for key in MODEL_KEYS:
        out[key]["delta_bic_vs_best_restricted"] = float(out[key]["bic"] - best)
        out[key]["aicc"] = np.nan
        out[key]["resid_r_TC"] = np.nan
        out[key]["resid_r_TI"] = np.nan
        out[key]["resid_r_CI"] = np.nan
    win, gap = winner_of(out)
    out["winner"] = win
    out["bic_gap"] = gap
    return out


def bootstrap_percents(groups: list[np.ndarray], prepare, n_boot: int, seed: int) -> dict:
    point_X = prepare(groups)
    point = {
        "T_to_C": mediation(point_X, "T_to_C"),
        "C_to_T": mediation(point_X, "C_to_T"),
    }
    rng = np.random.default_rng(seed)
    draws = {k: [] for k in point}
    n_skip = 0
    for _ in range(n_boot):
        resampled = [g[rng.choice(g.shape[0], size=g.shape[0], replace=True)] for g in groups]
        try:
            xb = prepare(resampled)
        except ZeroVariance:
            n_skip += 1
            continue
        for key in draws:
            draws[key].append(mediation(xb, key))
    summary = {"n_boot": n_boot, "n_skip": n_skip, "orderings": {}}

    def _ci(series: np.ndarray) -> tuple[float, float, int]:
        finite = series[np.isfinite(series)]
        if finite.size == 0:
            return np.nan, np.nan, 0
        lo, hi = np.quantile(finite, [0.025, 0.975])
        return float(lo), float(hi), int(finite.size)

    for key, med in point.items():
        packed = {name: np.array([d[name] for d in draws[key]], dtype=float) for name in ("percent", "indirect", "direct", "total", "a", "b")}
        total = packed["total"]
        ci = {name: _ci(packed[name]) for name in packed}
        summary["orderings"][key] = {
            **med,
            "ci_lo": ci["percent"][0],
            "ci_hi": ci["percent"][1],
            "indirect_ci_lo": ci["indirect"][0],
            "indirect_ci_hi": ci["indirect"][1],
            "direct_ci_lo": ci["direct"][0],
            "direct_ci_hi": ci["direct"][1],
            "total_ci_lo": ci["total"][0],
            "total_ci_hi": ci["total"][1],
            "n_finite": ci["percent"][2],
            "frac_small_total": float(np.mean(np.abs(total) < 0.02)) if total.size else np.nan,
            "unstable": bool(abs(med["total"]) < 0.02),
        }
    summary["point_X_n"] = int(point_X.shape[0])
    return summary


def fit_rows(spec: str, estimator: str, cohort: str, models: dict) -> list[dict]:
    rows = []
    for key in MODEL_KEYS:
        m = models[key]
        rows.append(
            {
                "spec": spec,
                "estimator": estimator,
                "cohort": cohort,
                "model": key,
                "model_name": m["model"],
                "n": m["n"],
                "k": m["k"],
                "ll": m["ll"],
                "aic": m["aic"],
                "bic": m["bic"],
                "aicc": m.get("aicc", np.nan),
                "delta_bic": m["delta_bic_vs_best_restricted"],
                "lr_vs_saturated": m["lr_vs_saturated"],
                "df": m["df_vs_saturated"],
                "p_vs_saturated": m["p_vs_saturated"],
                "srmr": m.get("srmr", np.nan),
                "resid_r_TC": m.get("resid_r_TC", np.nan),
                "resid_r_TI": m.get("resid_r_TI", np.nan),
                "resid_r_CI": m.get("resid_r_CI", np.nan),
                "is_winner": key == models["winner"],
            }
        )
    return rows


def prepare_z(groups: list[np.ndarray]) -> np.ndarray:
    return np.vstack([zscore_columns(g) for g in groups])


def make_tcga_prepare(k_idx: list[int]):
    def prepare(groups: list[np.ndarray]) -> np.ndarray:
        blocks = [residualize_z(g, (0, 1, 2), k_idx) for g in groups]
        return np.vstack(blocks)

    return prepare


def stratum_arrays(df: pd.DataFrame, columns: list[str], by: str, order: list[str]) -> list[np.ndarray]:
    groups = []
    for key in order:
        sub = df[df[by] == key]
        if sub.empty:
            raise RuntimeError(f"no rows for {key}")
        groups.append(sub[columns].to_numpy(dtype=float))
    return groups


def calibrate_concordant4(df: pd.DataFrame) -> dict:
    rows = []
    for col in ("pct_CLDN4", "mean_CLDN4", "pct_TACSTD2", "mean_TACSTD2"):
        rhos, ns = [], []
        per = []
        for ds, g in df.groupby("dataset", sort=True):
            r, p, n = partial_corr(g[col], g["frac_tnk"], rank=True)
            rhos.append(r)
            ns.append(n)
            per.append({"dataset": ds, "score": col, "rho": r, "p": p, "n": n})
        meta = dl_spearman(rhos, ns, k=0)
        meta["score"] = col
        rows.append(meta)
        rows[-1]["per"] = per
    locked = next(r for r in rows if r["score"] == "pct_CLDN4")
    if abs(locked["rho"] - LOCKED_DL_RHO) > 1e-9 or locked["N"] != 65:
        raise SystemExit(f"concordant-4 calibration failed: {locked}")
    return {"marginal": rows}


def calibrate_tcga(df: pd.DataFrame) -> None:
    luad = df[df.cohort == "LUAD"]
    if len(luad) != 516:
        raise SystemExit(f"LUAD n={len(luad)}")
    cd8 = luad[["CD8A", "CD8B"]].mean(axis=1).to_numpy()
    cov = [luad[c].to_numpy() for c in KERATIN["krt819"]]
    for gene, expected in LUAD_PARTIAL.items():
        rho, p, n = partial_corr(luad[gene].to_numpy(), cd8, cov, rank=True)
        if abs(rho - expected) > 1e-4 or n != 516:
            raise SystemExit(f"LUAD partial Spearman {gene}={rho} expected {expected} n={n}")
        print(f"[cal] LUAD {gene}–CD8 | KRT8/18/19 partial Spearman {rho:.6f} (published {expected})", flush=True)


def load_inputs():
    c4 = pd.read_csv(ROOT / "data" / "concordant4_patient_scores.tsv", sep="\t")
    tcga = pd.read_csv(ROOT / "data" / "tcga_primary01_genes.tsv.gz", sep="\t")
    if len(c4) != 65:
        raise SystemExit(f"concordant-4 n={len(c4)}")
    if int(tcga.isna().sum().sum()) != 0:
        raise SystemExit("TCGA table contains NA")
    funnel_n = int(tcga[tcga.cohort.isin(FUNNEL)].shape[0])
    if funnel_n != 3444:
        raise SystemExit(f"funnel n={funnel_n}")
    return c4, tcga


def c4_groups(df: pd.DataFrame, t_col: str, c_col: str, outcome: np.ndarray) -> tuple[list[str], list[np.ndarray]]:
    work = df[["dataset"]].copy()
    work["T"] = df[t_col].to_numpy(dtype=float)
    work["C"] = df[c_col].to_numpy(dtype=float)
    work["I"] = outcome
    order = sorted(work["dataset"].unique())
    return order, stratum_arrays(work, ["T", "C", "I"], "dataset", order)


def summarize_boot(spec: str, boot: dict) -> list[dict]:
    rows = []
    for key, med in boot["orderings"].items():
        direct = med["direct"]
        indirect = med["indirect"]
        inconsistent = bool(np.isfinite(direct) and np.isfinite(indirect) and direct * indirect < 0 and abs(direct) > 1e-8 and abs(indirect) > 1e-8)
        rows.append(
            {
                "spec": spec,
                "ordering": med["ordering"],
                "a": med["a"],
                "b": med["b"],
                "direct": direct,
                "indirect": indirect,
                "total": med["total"],
                "percent": med["percent"],
                "ci_lo": med["ci_lo"],
                "ci_hi": med["ci_hi"],
                "indirect_ci_lo": med["indirect_ci_lo"],
                "indirect_ci_hi": med["indirect_ci_hi"],
                "direct_ci_lo": med["direct_ci_lo"],
                "direct_ci_hi": med["direct_ci_hi"],
                "total_ci_lo": med["total_ci_lo"],
                "total_ci_hi": med["total_ci_hi"],
                "unstable": med["unstable"],
                "inconsistent": inconsistent,
                "frac_small_total": med["frac_small_total"],
                "n_boot": boot["n_boot"],
                "n_skip": boot["n_skip"],
                "n": boot["point_X_n"],
            }
        )
    return rows


def gap_word(gap: float) -> str:
    if gap < 2:
        return "small"
    if gap < 6:
        return "moderate"
    if gap < 10:
        return "large"
    return "larger than 10"


def percent_phrase(row: dict) -> str:
    if not np.isfinite(row["percent"]):
        return "not defined because the total effect is zero"
    base = (
        f"{row['percent']:.1f}% "
        f"(bootstrap 95% CI {row['ci_lo']:.1f} to {row['ci_hi']:.1f}; "
        f"standardized indirect {row['indirect']:.3f}, direct {row['direct']:.3f}, total {row['total']:.3f})"
    )
    notes = []
    if row["unstable"]:
        notes.append("total effect is under 0.02 in standardized units, so the percentage is unstable")
    if row["inconsistent"]:
        notes.append("direct and indirect paths have opposite signs (inconsistent mediation), so the percentage is not a share of a one-direction effect")
    if row.get("frac_small_total", 0) and row["frac_small_total"] >= 0.02:
        notes.append(
            f"{100 * row['frac_small_total']:.1f}% of bootstrap draws have a standardized total effect below 0.02 in absolute value, which stretches the ratio"
        )
    if notes:
        return base + "; " + "; ".join(notes)
    return base


def write_finding(path: Path, c4_primary, tcga_primary, sens_lines, cal_lines, two_stage_lines, reading_lines) -> None:
    c4_win = c4_primary["winner"]
    tc_win = tcga_primary["winner"]
    lines = []
    lines.append("# Which graph fits: TACSTD2, CLDN4, and immune")
    lines.append("")
    lines.append(
        "Three linear Gaussian graphs were fit to the same three observed variables. "
        "Each restricted graph has 5 free covariance parameters, so AIC and BIC rank them by likelihood. "
        "The saturated graph adds the direct path (6 parameters) and is the reference for a likelihood-ratio test."
    )
    lines.append("")
    lines.append(
        "M1 is TACSTD2 → CLDN4 → immune, with no direct TACSTD2 → immune arrow. "
        "M2 is CLDN4 → TACSTD2 → immune, with no direct CLDN4 → immune arrow. "
        "M3 is independent: TACSTD2 and CLDN4 are uncorrelated, and each has its own arrow into the immune readout."
    )
    lines.append("")
    lines.append(
        "A lower BIC is a better description of this covariance. "
        "It is not a knockdown, not an instrument, and not by itself a causal order. "
        "Concordant-4 and TCGA are not entered into one likelihood."
    )
    lines.append("")
    lines.append("## Concordant-4")
    lines.append("")
    lines.append(
        "Locked units only: GSE123902, GSE131907, GSE205335, and GSE189357 (n = 65). "
        "Primary variables are the malignant-cell detection percentages `pct_TACSTD2` and `pct_CLDN4`, "
        "and the unit T/NK fraction `frac_tnk`. "
        "Each variable is z-scored inside its dataset (maximum-likelihood standard deviation), then the four datasets are stacked. "
        "Percent positive is the locked CLDN4 measurement: the same table reproduces the locked DerSimonian–Laird Spearman."
    )
    lines.append("")
    lines.append(c4_primary["prose"])
    lines.append("")
    lines.append("## TCGA")
    lines.append("")
    lines.append(
        "UCSC Xena GDC STAR log2(TPM+1), primary solid tumor only, replicate aliquots averaged. "
        "Primary cohorts are the locked keratin funnel: LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD (n = 3,444). "
        "LUSC is estimated and is not in the sum. "
        "The immune readout is the CD8 score, the mean of CD8A and CD8B. "
        "Inside each cohort, TACSTD2, CLDN4, and the CD8 score are residualized on KRT8, KRT18, and KRT19, then z-scored. "
        "The fit that decides the TCGA comparison is the sum of those cohort-specific likelihoods. "
        "A pooled z-score fit is reported beside it."
    )
    lines.append("")
    lines.append(tcga_primary["prose"])
    lines.append("")
    lines.append("## Mediation percentage")
    lines.append("")
    lines.append(
        "The percentage is the product-of-coefficients indirect path divided by the total effect, "
        "with the direct path retained. "
        "On z-scored variables this is a standardized path decomposition. "
        "The bootstrap resamples patients inside each dataset or cohort, repeats the z-scoring "
        "(and, for TCGA, the keratin residualization), and refits. "
        f"There are {N_BOOT} replicates and the seed is {SEED}. "
        "The interval is the 2.5 and 97.5 percentiles."
    )
    lines.append("")
    lines.append(c4_primary["mediation_prose"])
    lines.append("")
    lines.append(tcga_primary["mediation_prose"])
    lines.append("")
    lines.extend(reading_lines)
    lines.append("")
    lines.append("## Two-stage residual")
    lines.append("")
    lines.append(
        "Conditional-independence check: residualize the immune variable on the hypothesized mediator, "
        "then correlate that residual with the other gene. "
        "M1 says the TACSTD2 residual correlation should be the one near zero. "
        "M2 says the CLDN4 residual correlation should be the one near zero. "
        "M3 says the TACSTD2–CLDN4 correlation itself should be near zero."
    )
    lines.append("")
    lines.append(
        "Effect decomposition: residualize the mediator on the exposure, then regress the immune variable "
        "on the exposure and that residual. "
        "The exposure coefficient equals the total effect, and the residual coefficient equals the mediator coefficient. "
        "On these matrices the identity holds to numerical error (the script checks it before writing results)."
    )
    lines.append("")
    for line in two_stage_lines:
        lines.append(line)
        lines.append("")
    lines.append("## Sensitivities")
    lines.append("")
    lines.append(
        "The primary specifications above were fixed first. "
        "The rows below use the same three graphs. "
        "A sensitivity that picks a different graph is listed as a disagreement, not averaged away."
    )
    lines.append("")
    lines.extend(sens_lines)
    lines.append("")
    lines.append("## Calibration")
    lines.append("")
    lines.extend(cal_lines)
    lines.append("")
    lines.append("## What this does not identify")
    lines.append("")
    lines.append(
        "Once both the TACSTD2–CLDN4 arrow and the direct arrow into immune are free, the saturated Gaussian model "
        "fits the 3×3 covariance exactly. The two directions of the gene–gene arrow are then the same likelihood, "
        "and AIC cannot choose between them. "
        "Only the restricted graphs, which drop one association, are separated by AIC and BIC."
    )
    lines.append("")
    lines.append(
        "Concordant-4 percent-positive and T/NK fraction are bounded. "
        "The linear SEM is a covariance approximation, which is why the mean log1p and logit-fraction sensitivities are in the table. "
        "TCGA CD8 is a bulk transcript score after a linear keratin residual, not a spatial immune count and not a purity-adjusted fraction."
    )
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def model_sentence(tag: str, models: dict, estimator_phrase: str) -> str:
    win = models["winner"]
    gap = models["bic_gap"]
    second = sorted((models[k]["bic"], k) for k in RESTRICTED)[1][1]
    sat_delta = models["saturated"]["delta_bic_vs_best_restricted"]
    bits = [
        f"On {tag}, {estimator_phrase}, the lowest BIC is {LABEL[win]} "
        f"(ΔBIC versus {LABEL[second]} = {gap:.2f}, a {gap_word(gap)} gap). "
        f"AIC picks the same graph." if _aic_same(models, win) else
        f"On {tag}, {estimator_phrase}, the lowest BIC is {LABEL[win]} "
        f"(ΔBIC versus {LABEL[second]} = {gap:.2f}, a {gap_word(gap)} gap). "
        f"AIC picks {LABEL[_aic_winner(models)]} instead."
    ]
    # The ternary above is a bit awkward if I embed it wrong. Let me build cleanly below.
    return bits[0] if False else _model_sentence(tag, models, estimator_phrase)


def _aic_winner(models: dict) -> str:
    return min(RESTRICTED, key=lambda k: (models[k]["aic"], k))


def _aic_same(models: dict, bic_win: str) -> bool:
    return _aic_winner(models) == bic_win


def _model_sentence(tag: str, models: dict, estimator_phrase: str) -> str:
    win = models["winner"]
    gap = models["bic_gap"]
    ranked = sorted((models[k]["bic"], k) for k in RESTRICTED)
    second = ranked[1][1]
    third = ranked[2][1]
    sat = models["saturated"]
    sat_delta = sat["delta_bic_vs_best_restricted"]
    aic_win = _aic_winner(models)
    aic_clause = "AIC selects the same graph." if aic_win == win else f"AIC selects {LABEL[aic_win]} instead."
    gap_clause = f"ΔBIC {gap:.2f} versus {LABEL[second]}"
    if gap >= 10:
        gap_clause += " (greater than 10)"
    else:
        gap_clause += f" (a {gap_word(gap)} gap)"
    text = (
        f"On {tag}, {estimator_phrase}, the lowest BIC is **{LABEL[win]}** "
        f"(BIC {models[win]['bic']:.1f}; {gap_clause}; "
        f"ΔBIC {models[third]['delta_bic_vs_best_restricted']:.2f} versus {LABEL[third]}). "
        f"{aic_clause} "
    )
    aic_gap = float(sat["aic"] - models[win]["aic"])
    if sat_delta < 0:
        text += (
            f"The saturated model has a still lower BIC (ΔBIC {sat_delta:.2f} relative to the best restricted graph; "
            f"likelihood-ratio versus the winner χ² = {models[win]['lr_vs_saturated']:.2f} on {models[win]['df_vs_saturated']} df, "
            f"p = {fmt_p(models[win]['p_vs_saturated'])}). "
            "The one-constraint graphs leave residual association. "
        )
    else:
        text += (
            f"The saturated model does not improve BIC (ΔBIC {sat_delta:.2f}; "
            f"likelihood-ratio versus the winner χ² = {models[win]['lr_vs_saturated']:.2f} on {models[win]['df_vs_saturated']} df, "
            f"p = {fmt_p(models[win]['p_vs_saturated'])}). "
        )
    if aic_gap < -2 and sat_delta > 0:
        text += (
            f"AIC still prefers the saturated model (AIC {sat['aic']:.1f} versus {models[win]['aic']:.1f} for the BIC winner). "
            "The extra direct-path parameters improve the likelihood enough for AIC and not enough for BIC."
        )
    elif aic_gap > 2 and sat_delta < 0:
        text += f"AIC agrees with BIC that the saturated model is worse (ΔAIC {aic_gap:.1f})."
    return text.strip()


def mediation_sentence(title: str, rows: list[dict]) -> str:
    by = {r["ordering"]: r for r in rows}
    t = by["TACSTD2_to_CLDN4_to_immune"]
    c = by["CLDN4_to_TACSTD2_to_immune"]
    return (
        f"{title}: TACSTD2 → CLDN4 → immune is {percent_phrase(t)}. "
        f"CLDN4 → TACSTD2 → immune is {percent_phrase(c)}."
    )


def two_stage_sentence(title: str, stage: dict, rank_note: str) -> str:
    return (
        f"{title}: Pearson residual correlation of TACSTD2 with immune after CLDN4 is "
        f"{stage['partial_r_TACSTD2_immune_given_CLDN4']:.3f} "
        f"(p = {fmt_p(stage['partial_p_TACSTD2_immune_given_CLDN4'])}). "
        f"Pearson residual correlation of CLDN4 with immune after TACSTD2 is "
        f"{stage['partial_r_CLDN4_immune_given_TACSTD2']:.3f} "
        f"(p = {fmt_p(stage['partial_p_CLDN4_immune_given_TACSTD2'])}). "
        f"Pearson correlation of TACSTD2 with CLDN4 is {stage['pearson_TACSTD2_CLDN4']:.3f} "
        f"(p = {fmt_p(stage['pearson_p_TACSTD2_CLDN4'])}). "
        f"{rank_note}"
    )


def plot_bic(path: Path, panels: list[tuple[str, dict]]) -> None:
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.4), squeeze=False)
    colors = {"M1": "#0072B2", "M2": "#E69F00", "M3": "#009E73", "saturated": "#666666"}
    short = {"M1": "M1\nT→C→immune", "M2": "M2\nC→T→immune", "M3": "M3\nindependent", "saturated": "saturated"}
    for ax, (title, models) in zip(axes[0], panels):
        keys = list(MODEL_KEYS)
        delta = [models[k]["delta_bic_vs_best_restricted"] for k in keys]
        cap = 40.0
        shown = [min(cap, v) if v > 0 else v for v in delta]
        ax.bar([short[k] for k in keys], shown, color=[colors[k] for k in keys], width=0.78)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel("ΔBIC versus best restricted graph")
        ax.set_title(title, fontsize=10)
        top = max(shown + [1])
        ax.set_ylim(min(-1.0, min(shown) - 2), top * 1.28 + 1)
        for i, val in enumerate(delta):
            label = f"{val:.1f}" if val <= cap else f"{val:.0f}\n(capped)"
            y = shown[i]
            ax.text(i, y + top * 0.03, label, ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_mediation(path: Path, series: list[tuple[str, list[dict]]]) -> None:
    """Standardized indirect and direct paths. The percentage is their ratio and is plotted separately only when the interval is finite and narrow."""
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    labels = []
    y = 0.0
    positions = []
    for title, rows in series:
        for row in rows:
            short = "T→C→immune" if row["ordering"].startswith("TACSTD2") else "C→T→immune"
            for key, color, marker in (("indirect", "#0072B2", "o"), ("direct", "#D55E00", "s")):
                labels.append(f"{title} {short}\n{key}")
                positions.append(y)
                lo, hi = row[f"{key}_ci_lo"], row[f"{key}_ci_hi"]
                est = row[key]
                if np.isfinite(est) and np.isfinite(lo) and np.isfinite(hi):
                    ax.errorbar(est, y, xerr=[[est - lo], [hi - est]], fmt=marker, color=color, capsize=3, markersize=5)
                y += 1.0
            y += 0.35
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Standardized coefficient")
    ax.set_title("Direct path retained")
    ax = axes[1]
    y = 0.0
    positions = []
    labels = []
    for title, rows in series:
        for row in rows:
            short = "T→C→immune" if row["ordering"].startswith("TACSTD2") else "C→T→immune"
            labels.append(f"{title}\n{short}")
            positions.append(y)
            width = row["ci_hi"] - row["ci_lo"] if np.isfinite(row["ci_lo"]) and np.isfinite(row["ci_hi"]) else np.inf
            if np.isfinite(row["percent"]) and width < 200:
                ax.errorbar(
                    row["percent"],
                    y,
                    xerr=[[row["percent"] - row["ci_lo"]], [row["ci_hi"] - row["percent"]]],
                    fmt="o",
                    color="#0072B2",
                    capsize=3,
                )
            elif np.isfinite(row["percent"]):
                labels[-1] = labels[-1] + "\n(ratio CI too wide to draw)"
            y += 1.0
    ax.axvline(0, color="black", lw=0.7)
    ax.axvline(100, color="#999999", lw=0.6, ls="--")
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Mediation % (indirect / total × 100)")
    ax.set_title("Ratio, with the direct path kept")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_partial(path: Path, rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    y = np.arange(len(rows))
    ax.barh(y - 0.15, [abs(r["r_T_given_C"]) for r in rows], height=0.3, color="#0072B2", label="|r| TACSTD2–immune | CLDN4")
    ax.barh(y + 0.15, [abs(r["r_C_given_T"]) for r in rows], height=0.3, color="#E69F00", label="|r| CLDN4–immune | TACSTD2")
    xmax = max(max(abs(r["r_T_given_C"]), abs(r["r_C_given_T"])) for r in rows)
    ax.set_xlim(0, xmax * 1.35 + 0.02)
    for i, r in enumerate(rows):
        ax.text(abs(r["r_T_given_C"]) + 0.01, i - 0.15, f"{r['r_T_given_C']:+.3f}", va="center", fontsize=7, color="#0072B2")
        ax.text(abs(r["r_C_given_T"]) + 0.01, i + 0.15, f"{r['r_C_given_T']:+.3f}", va="center", fontsize=7, color="#8a5a00")
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlabel("Absolute residual Pearson correlation (labels are signed)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Two-stage residual: which gene still tracks immune")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def md_table(headers: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return out


def rank_meta_sentence(df_groups: list[tuple[str, np.ndarray]]) -> str:
    """DL partial Spearman across strata. Each array is T, C, I."""
    summaries = []
    for name, x_i, y_i, c_i in (
        ("TACSTD2–immune | CLDN4", 0, 2, 1),
        ("CLDN4–immune | TACSTD2", 1, 2, 0),
        ("TACSTD2–CLDN4", 0, 1, None),
    ):
        rhos, ns = [], []
        for _, block in df_groups:
            if c_i is None:
                r, p, n = partial_corr(block[:, x_i], block[:, y_i], rank=True)
            else:
                r, p, n = partial_corr(block[:, x_i], block[:, y_i], [block[:, c_i]], rank=True)
            rhos.append(r)
            ns.append(n)
        meta = dl_spearman(rhos, ns, k=0 if c_i is None else 1)
        summaries.append(f"{name} DL ρ = {meta['rho']:.3f} (p = {fmt_p(meta['p'])}, I² = {100 * meta['I2']:.0f}%)")
    return "Within-stratum partial Spearman, DerSimonian–Laird: " + "; ".join(summaries) + "."


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    fig_dir = ROOT / "results" / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    c4, tcga = load_inputs()
    cal = calibrate_concordant4(c4)
    calibrate_tcga(tcga)
    print("[cal] concordant-4 pct_CLDN4 DL rho matches locked value", flush=True)

    fit_table = []
    med_table = []
    stage_table = []
    corr_table = []
    partial_plot_rows = []

    # ----- Concordant-4 specs -----
    c4_specs = []
    outcome_frac = c4["frac_tnk"].to_numpy(dtype=float)
    logit_frac = np.log(outcome_frac / (1.0 - outcome_frac))
    c4_specs.append(("c4_pct", "pct_TACSTD2", "pct_CLDN4", outcome_frac, c4.index == c4.index, True))
    c4_specs.append(("c4_mean", "mean_TACSTD2", "mean_CLDN4", outcome_frac, c4.index == c4.index, False))
    c4_specs.append(("c4_pct_logit", "pct_TACSTD2", "pct_CLDN4", logit_frac, c4.index == c4.index, False))
    min30 = c4["n_malignant"] >= 30
    c4_specs.append(("c4_pct_min30_malignant", "pct_TACSTD2", "pct_CLDN4", outcome_frac, min30, False))

    c4_done = {}
    for name, t_col, c_col, outcome, mask, is_primary in c4_specs:
        sub = c4.loc[mask].copy()
        sub["_outcome"] = outcome[mask.to_numpy()] if hasattr(mask, "to_numpy") else outcome
        # outcome was built on c4 row order. Align by position.
        aligned = c4[["dataset", t_col, c_col, "n_malignant"]].copy()
        aligned["I"] = outcome
        aligned = aligned.loc[mask].copy()
        order = sorted(aligned["dataset"].unique())
        groups = stratum_arrays(aligned.rename(columns={t_col: "T", c_col: "C"}), ["T", "C", "I"], "dataset", order)
        if any(g.shape[0] < 8 for g in groups):
            print(f"[skip] {name}: a dataset has fewer than 8 units", flush=True)
            continue
        pooled = fit_matrix(prepare_z(groups))
        fit_table.extend(fit_rows(name, "pooled_within_z", "ALL", pooled))
        # Per-dataset fits and the sum.
        per = []
        for ds, block in zip(order, groups):
            z = zscore_columns(block)
            fitted = fit_matrix(z)
            raw = fit_matrix(center(block))
            if fitted["winner"] != raw["winner"] and min(
                abs(fitted[a]["bic"] - fitted[b]["bic"]) for a in RESTRICTED for b in RESTRICTED if a < b
            ) > 0.05:
                # Ranking is invariant to scaling when the gap is not numerical dust.
                gap_here = fitted["bic_gap"]
                if gap_here > 0.05 and raw["bic_gap"] > 0.05 and fitted["winner"] != raw["winner"]:
                    raise RuntimeError(f"{name} {ds}: z-score winner {fitted['winner']} != raw {raw['winner']}")
            per.append(fitted)
            fit_table.extend(fit_rows(name, "within_dataset", ds, fitted))
        summed = stack_fits(per)
        fit_table.extend(fit_rows(name, "dataset_bic_sum", "ALL", summed))
        boot = bootstrap_percents(groups, prepare_z, N_BOOT, SEED)
        med_rows = summarize_boot(name, boot)
        med_table.extend(med_rows)
        stage = pooled["two_stage"]
        stage_table.append(
            {
                "spec": name,
                "estimator": "pooled_within_z",
                **{k: stage[k] for k in stage},
            }
        )
        r = pooled["corr"]
        corr_table.append(
            {
                "spec": name,
                "estimator": "pooled_within_z",
                "r_TACSTD2_CLDN4": r[0, 1],
                "r_TACSTD2_immune": r[0, 2],
                "r_CLDN4_immune": r[1, 2],
                "n": pooled["M1"]["n"],
            }
        )
        rank_note = rank_meta_sentence(list(zip(order, groups)))
        c4_done[name] = {
            "pooled": pooled,
            "summed": summed,
            "boot_rows": med_rows,
            "rank_note": rank_note,
            "stage": stage,
            "n_units": int(sum(g.shape[0] for g in groups)),
            "per_winners": {ds: fitted["winner"] for ds, fitted in zip(order, per)},
        }
        print(f"[c4] {name} BIC winner {pooled['winner']} gap {pooled['bic_gap']:.2f}; sum winner {summed['winner']}", flush=True)

    # ----- TCGA specs -----
    tcga_specs = [
        ("tcga_cd8_krt819", "CD8", "krt819", True, FUNNEL),
        ("tcga_cd3_krt819", "CD3", "krt819", False, FUNNEL),
        ("tcga_cytotoxic_krt819", "cytotoxic", "krt819", False, FUNNEL),
        ("tcga_cd8_none", "CD8", "none", False, FUNNEL),
        ("tcga_cd8_krt56", "CD8", "krt56", False, FUNNEL),
        ("tcga_cd3_krt56", "CD3", "krt56", False, FUNNEL),
        ("tcga_cd8_krt819_plus_lusc", "CD8", "krt819", False, FUNNEL + ["LUSC"]),
    ]
    tcga_done = {}
    for name, immune, keratin, is_primary, cohorts in tcga_specs:
        genes = IMMUNE_GENES[immune]
        kgenes = KERATIN[keratin]
        work = tcga[tcga.cohort.isin(cohorts)].copy()
        work["T"] = work["TACSTD2"]
        work["C"] = work["CLDN4"]
        work["I"] = work[genes].mean(axis=1)
        for i, gname in enumerate(kgenes):
            work[f"K{i}"] = work[gname]
        kcols = [f"K{i}" for i in range(len(kgenes))]
        cols = ["T", "C", "I"] + kcols
        groups = stratum_arrays(work, cols, "cohort", cohorts)
        k_idx = list(range(3, 3 + len(kgenes)))
        prepare = make_tcga_prepare(k_idx)
        pooled = fit_matrix(prepare(groups))
        fit_table.extend(fit_rows(name, "pooled_within_z", "ALL", pooled))
        per = []
        per_names = []
        for cohort, block in zip(cohorts, groups):
            z = residualize_z(block, (0, 1, 2), k_idx)
            fitted = fit_matrix(z)
            # Unscaled residuals, same adjustment.
            if k_idx:
                design = np.column_stack([np.ones(len(block)), block[:, k_idx]])
                raw = np.column_stack(
                    [
                        block[:, j] - design @ np.linalg.lstsq(design, block[:, j], rcond=None)[0]
                        for j in (0, 1, 2)
                    ]
                )
            else:
                raw = block[:, :3]
            raw_fit = fit_matrix(raw)
            if fitted["winner"] != raw_fit["winner"] and fitted["bic_gap"] > 0.05 and raw_fit["bic_gap"] > 0.05:
                raise RuntimeError(f"{name} {cohort}: scaling changed the BIC winner")
            per.append(fitted)
            per_names.append(cohort)
            est = "within_cohort_primary" if cohort in FUNNEL else "within_cohort_extra"
            fit_table.extend(fit_rows(name, est, cohort, fitted))
        summed = stack_fits(per)
        fit_table.extend(fit_rows(name, "cohort_bic_sum", "ALL", summed))
        boot = bootstrap_percents(groups, prepare, N_BOOT, SEED)
        med_rows = summarize_boot(name, boot)
        med_table.extend(med_rows)
        stage = pooled["two_stage"]
        stage_table.append({"spec": name, "estimator": "pooled_within_z", **{k: stage[k] for k in stage}})
        r = pooled["corr"]
        corr_table.append(
            {
                "spec": name,
                "estimator": "pooled_within_z",
                "r_TACSTD2_CLDN4": r[0, 1],
                "r_TACSTD2_immune": r[0, 2],
                "r_CLDN4_immune": r[1, 2],
                "n": pooled["M1"]["n"],
            }
        )
        # Per-cohort mediation point estimates (no extra bootstrap).
        for cohort, block in zip(cohorts, groups):
            z = residualize_z(block, (0, 1, 2), k_idx)
            for ordering in ("T_to_C", "C_to_T"):
                med = mediation(z, ordering)
                med_table.append(
                    {
                        "spec": name,
                        "ordering": med["ordering"],
                        "cohort": cohort,
                        "a": med["a"],
                        "b": med["b"],
                        "direct": med["direct"],
                        "indirect": med["indirect"],
                        "total": med["total"],
                        "percent": med["percent"],
                        "ci_lo": np.nan,
                        "ci_hi": np.nan,
                        "indirect_ci_lo": np.nan,
                        "indirect_ci_hi": np.nan,
                        "direct_ci_lo": np.nan,
                        "direct_ci_hi": np.nan,
                        "total_ci_lo": np.nan,
                        "total_ci_hi": np.nan,
                        "unstable": abs(med["total"]) < 0.02,
                        "inconsistent": bool(med["direct"] * med["indirect"] < 0),
                        "frac_small_total": np.nan,
                        "n_boot": 0,
                        "n_skip": 0,
                        "n": int(z.shape[0]),
                    }
                )
        rank_note = rank_meta_sentence(list(zip(cohorts, [g[:, :3] for g in groups])))
        # The rank meta above is on raw scores, not keratin residuals.
        # Recompute on residualized (not z-scored) values so the rank check matches the model.
        resid_blocks = []
        for cohort, block in zip(cohorts, groups):
            if k_idx:
                design = np.column_stack([np.ones(len(block)), block[:, k_idx]])
                raw = np.column_stack(
                    [
                        block[:, j] - design @ np.linalg.lstsq(design, block[:, j], rcond=None)[0]
                        for j in (0, 1, 2)
                    ]
                )
            else:
                raw = block[:, :3]
            resid_blocks.append((cohort, raw))
        rank_note = rank_meta_sentence(resid_blocks)
        tcga_done[name] = {
            "pooled": pooled,
            "summed": summed,
            "boot_rows": [r for r in med_rows],
            "rank_note": rank_note,
            "stage": stage,
            "per_winners": {ds: fitted["winner"] for ds, fitted in zip(per_names, per)},
            "n": int(sum(g.shape[0] for g in groups)),
        }
        votes = pd.Series(tcga_done[name]["per_winners"]).value_counts().to_dict()
        print(
            f"[tcga] {name} sum winner {summed['winner']} gap {summed['bic_gap']:.2f}; "
            f"pooled winner {pooled['winner']}; votes {votes}",
            flush=True,
        )

    fit_df = pd.DataFrame(fit_table)
    med_df = pd.DataFrame(med_table)
    stage_df = pd.DataFrame(stage_table)
    corr_df = pd.DataFrame(corr_table)
    fit_df.to_csv(table_dir / "model_fit.tsv", sep="\t", index=False)
    med_df.to_csv(table_dir / "mediation.tsv", sep="\t", index=False)
    stage_df.to_csv(table_dir / "two_stage_residual.tsv", sep="\t", index=False)
    corr_df.to_csv(table_dir / "analysis_correlations.tsv", sep="\t", index=False)

    # Calibration table
    cal_rows = []
    for meta in cal["marginal"]:
        cal_rows.append(
            {
                "layer": "concordant4_dl_spearman",
                "score": meta["score"],
                "rho": meta["rho"],
                "p": meta["p"],
                "I2": meta["I2"],
                "N": meta["N"],
            }
        )
    pd.DataFrame(cal_rows).to_csv(table_dir / "calibration.tsv", sep="\t", index=False)

    c4p = c4_done["c4_pct"]
    tcp = tcga_done["tcga_cd8_krt819"]
    c4_primary = {
        "winner": c4p["pooled"]["winner"],
        "prose": _model_sentence(
            "the concordant-4 primary specification (percent positive, n = 65)",
            c4p["pooled"],
            "using one SEM on the within-dataset z-scores",
        )
        + " "
        + _model_sentence(
            "the sum of the four dataset-specific BICs",
            c4p["summed"],
            "with a separate 5-parameter model in each dataset",
        ),
        "mediation_prose": mediation_sentence("Concordant-4 primary", c4p["boot_rows"]),
    }
    tc_primary = {
        "winner": tcp["summed"]["winner"],
        "prose": _model_sentence(
            "the TCGA primary specification (CD8 score, KRT8/18/19 residuals, seven cohorts, n = 3,444)",
            tcp["summed"],
            "summing cohort-specific BICs",
        )
        + " "
        + _model_sentence(
            "the same residuals stacked after within-cohort z-scoring",
            tcp["pooled"],
            "as one pooled SEM",
        ),
        "mediation_prose": mediation_sentence("TCGA primary, pooled within-cohort z-scores", tcp["boot_rows"]),
    }

    # Sensitivity markdown
    sens_rows = []
    for name, blob, estimator in (
        *((n, c4_done[n], "pooled_within_z") for n in c4_done),
        *((n, tcga_done[n], "cohort_bic_sum") for n in tcga_done),
    ):
        models = blob["pooled"] if estimator == "pooled_within_z" else blob["summed"]
        meds = {r["ordering"]: r for r in blob["boot_rows"]}
        tmed = meds["TACSTD2_to_CLDN4_to_immune"]
        cmed = meds["CLDN4_to_TACSTD2_to_immune"]
        role = "primary" if name in ("c4_pct", "tcga_cd8_krt819") else "sensitivity"
        sens_rows.append(
            [
                role,
                name,
                estimator,
                str(models["M1"]["n"]),
                models["winner"],
                f"{models['bic_gap']:.2f}",
                f"{models['saturated']['delta_bic_vs_best_restricted']:.2f}",
                fmt_num(tmed["percent"], 1),
                fmt_num(cmed["percent"], 1),
            ]
        )
    sens_lines = md_table(
        ["role", "spec", "estimator", "n", "BIC winner", "ΔBIC 2nd", "saturated ΔBIC", "M1 mediation %", "M2 mediation %"],
        sens_rows,
    )
    # Vote line for TCGA primary cohorts
    vote = pd.Series(tcp["per_winners"])
    vote = vote[vote.index.isin(FUNNEL)]
    sens_lines.append("")
    sens_lines.append(
        "TCGA primary, cohort-by-cohort BIC winners: "
        + ", ".join(f"{k} {int(v)}" for k, v in vote.value_counts().items())
        + f" of {len(FUNNEL)} funnel cohorts. "
        + "Per-cohort labels: "
        + ", ".join(f"{k} {v}" for k, v in tcp["per_winners"].items())
        + "."
    )
    sens_lines.append("")
    sens_lines.append(
        "Concordant-4 primary, dataset-by-dataset BIC winners: "
        + ", ".join(f"{k} {v}" for k, v in c4p["per_winners"].items())
        + "."
    )

    cal_lines = [
        f"Concordant-4 DerSimonian–Laird Spearman of `pct_CLDN4` versus `frac_tnk` is {cal['marginal'][0]['rho']:.6f} "
        if False
        else "",
    ]
    by_score = {m["score"]: m for m in cal["marginal"]}
    cal_lines = [
        (
            f"Concordant-4 DerSimonian–Laird Spearman of malignant CLDN4 percent positive versus T/NK fraction "
            f"is {by_score['pct_CLDN4']['rho']:.6f} (p = {fmt_p(by_score['pct_CLDN4']['p'])}, I² = {100 * by_score['pct_CLDN4']['I2']:.0f}%, N = 65). "
            "That matches the locked value −0.5311678045689989 from the concordant-4 patient table."
        ),
        (
            f"The same meta-analysis on mean log1p is {by_score['mean_CLDN4']['rho']:.3f} for CLDN4 and "
            f"{by_score['mean_TACSTD2']['rho']:.3f} for TACSTD2; percent-positive TACSTD2 is {by_score['pct_TACSTD2']['rho']:.3f}."
        ),
        (
            "LUAD partial Spearman of TACSTD2 versus the CD8 score after KRT8/18/19 is "
            f"{LUAD_PARTIAL['TACSTD2']}, and CLDN4 versus the CD8 score is {LUAD_PARTIAL['CLDN4']}. "
            "This extract reproduces both published values within 0.0001, and the funnel patient count is 3,444."
        ),
    ]
    # Use computed LUAD rhos rather than only the published constants. Recompute for the sentence.
    luad = tcga[tcga.cohort == "LUAD"]
    cd8 = luad[["CD8A", "CD8B"]].mean(axis=1).to_numpy()
    cov = [luad[c].to_numpy() for c in KERATIN["krt819"]]
    rho_t, _, _ = partial_corr(luad["TACSTD2"].to_numpy(), cd8, cov, rank=True)
    rho_c, _, _ = partial_corr(luad["CLDN4"].to_numpy(), cd8, cov, rank=True)
    cal_lines[2] = (
        f"LUAD partial Spearman of TACSTD2 versus the CD8 score after KRT8/18/19 is {rho_t:.6f} "
        f"(published −0.103504), and CLDN4 versus the CD8 score is {rho_c:.6f} (published −0.084776). "
        "Patient counts match the earlier Xena extract (funnel n = 3,444)."
    )

    two_stage_lines = [
        two_stage_sentence("Concordant-4 primary", c4p["stage"], c4p["rank_note"]),
        two_stage_sentence("TCGA primary pooled residuals", tcp["stage"], tcp["rank_note"]),
    ]

    def _med_pair(spec: str, cohort: str | None = None) -> dict[str, pd.Series]:
        sub = med_df[med_df.spec == spec]
        if cohort is None:
            sub = sub[sub.n_boot > 0]
        else:
            sub = sub[sub.cohort == cohort]
        return {str(r.ordering): r for _, r in sub.iterrows()}

    def _coef_phrase(row: pd.Series, key: str) -> str:
        return (
            f"{row[key]:.3f} (bootstrap 95% CI {row[key + '_ci_lo']:.3f} to {row[key + '_ci_hi']:.3f})"
        )

    c4_t = _med_pair("c4_pct")["TACSTD2_to_CLDN4_to_immune"]
    c4_c = _med_pair("c4_pct")["CLDN4_to_TACSTD2_to_immune"]
    tc_t = _med_pair("tcga_cd8_krt819")["TACSTD2_to_CLDN4_to_immune"]
    tc_c = _med_pair("tcga_cd8_krt819")["CLDN4_to_TACSTD2_to_immune"]
    reading_lines = [
        "## Reading the percentage and the cohort split",
        "",
        (
            "On concordant-4, TACSTD2 → CLDN4 → T/NK splits into a standardized indirect path of "
            f"{_coef_phrase(c4_t, 'indirect')} and a direct path of {_coef_phrase(c4_t, 'direct')}. "
            f"The total effect is {_coef_phrase(c4_t, 'total')}. "
            "The indirect path is negative and the direct path is positive, so they oppose each other. "
            f"The ratio of the indirect path to that net total is {c4_t['percent']:.1f}%, "
            f"and the bootstrap interval on the ratio is {c4_t['ci_lo']:.0f} to {c4_t['ci_hi']:.0f}. "
            "The total-effect interval includes zero. The indirect-path interval does not. "
            "The ratio is not an estimate of a mediation share."
        ),
        "",
        (
            "The reverse order, CLDN4 → TACSTD2 → T/NK, has indirect "
            f"{_coef_phrase(c4_c, 'indirect')}, direct {_coef_phrase(c4_c, 'direct')}, "
            f"and total {_coef_phrase(c4_c, 'total')}. "
            f"The ratio is {c4_c['percent']:.1f}% (bootstrap 95% CI {c4_c['ci_lo']:.1f} to {c4_c['ci_hi']:.1f}). "
            "The CLDN4 association with T/NK stays on the direct path. It is not carried by TACSTD2."
        ),
        "",
        (
            "On the pooled TCGA CD8 residuals both paths in TACSTD2 → CLDN4 → CD8 are negative: indirect "
            f"{_coef_phrase(tc_t, 'indirect')}, direct {_coef_phrase(tc_t, 'direct')}, "
            f"total {_coef_phrase(tc_t, 'total')}. "
            f"The ratio is {tc_t['percent']:.1f}% (bootstrap 95% CI {tc_t['ci_lo']:.1f} to {tc_t['ci_hi']:.1f}). "
            "The interval extends above 100, so the split between the indirect path and the direct path is not pinned down. "
            f"The reverse order is {tc_c['percent']:.1f}% (CI {tc_c['ci_lo']:.1f} to {tc_c['ci_hi']:.1f}): "
            "most of the CLDN4 association with the CD8 score is direct, and that percentage's interval includes zero."
        ),
        "",
        "The TCGA sum is not a 7-cohort vote. Cohort-specific BIC gaps on the primary CD8 specification:",
        "",
    ]
    cohort_rows = []
    within = fit_df[
        (fit_df.spec == "tcga_cd8_krt819")
        & (fit_df.estimator.str.startswith("within"))
        & (fit_df.model.isin(RESTRICTED))
    ]
    for cohort_name in FUNNEL:
        ranked = within[within.cohort == cohort_name].sort_values(["delta_bic", "model"])
        winrow = ranked.iloc[0]
        gap = float(ranked.iloc[1].delta_bic)
        pair = _med_pair("tcga_cd8_krt819", cohort_name)

        def _pct(row: pd.Series) -> str:
            text = fmt_num(row["percent"], 1)
            notes = []
            if bool(row["unstable"]):
                notes.append("unstable total")
            if bool(row["inconsistent"]):
                notes.append("opposite signs")
            if notes:
                text += " (" + "; ".join(notes) + ")"
            return text

        cohort_rows.append(
            [
                cohort_name,
                str(int(winrow.n)),
                str(winrow.model),
                f"{gap:.2f}",
                _pct(pair["TACSTD2_to_CLDN4_to_immune"]),
                _pct(pair["CLDN4_to_TACSTD2_to_immune"]),
            ]
        )
    reading_lines.extend(
        md_table(
            ["cohort", "n", "BIC winner", "ΔBIC to 2nd", "T→C→immune %", "C→T→immune %"],
            cohort_rows,
        )
    )
    lusc_ranked = fit_df[
        (fit_df.spec == "tcga_cd8_krt819_plus_lusc") & (fit_df.cohort == "LUSC") & (fit_df.model.isin(RESTRICTED))
    ].sort_values(["delta_bic", "model"])
    lusc_win = lusc_ranked.iloc[0]
    lusc_gap = float(lusc_ranked.iloc[1].delta_bic)
    krt = tcga_done["tcga_cd8_krt56"]
    krt_stage = krt["stage"]
    gse_bits_list = []
    gse_all = fit_df[
        (fit_df.spec == "c4_pct") & (fit_df.estimator == "within_dataset") & (fit_df.model.isin(list(RESTRICTED)))
    ]
    for dataset_name, sub in gse_all.groupby("cohort", sort=True):
        ranked = sub.sort_values(["delta_bic", "model"])
        gse_bits_list.append(
            f"{dataset_name} {ranked.iloc[0].model} (ΔBIC to 2nd {float(ranked.iloc[1].delta_bic):.2f}, n={int(ranked.iloc[0].n)})"
        )
    gse_bits = ", ".join(gse_bits_list)
    reading_lines.extend(
        [
            "",
            (
                f"LUSC alone, same CD8 score and KRT8/18/19 residual, prefers {lusc_win.model} by ΔBIC {lusc_gap:.2f} "
                f"(n = {int(lusc_win.n)}). Adding LUSC to the seven-cohort sum selects "
                f"{tcga_done['tcga_cd8_krt819_plus_lusc']['summed']['winner']} "
                f"(ΔBIC {tcga_done['tcga_cd8_krt819_plus_lusc']['summed']['bic_gap']:.2f}). "
                "LUSC was not in the pre-specified funnel."
            ),
            "",
            (
                "Replacing KRT8/18/19 with KRT5+KRT6A+KRT6B+KRT14 on the same seven cohorts leaves the cohort-sum gap at "
                f"ΔBIC {krt['summed']['bic_gap']:.2f} for {krt['summed']['winner']}, while the pooled SEM prefers "
                f"{krt['pooled']['winner']} (ΔBIC {krt['pooled']['bic_gap']:.2f}). "
                f"The pooled residual correlations are {krt_stage['partial_r_TACSTD2_immune_given_CLDN4']:.3f} for TACSTD2 after CLDN4 and "
                f"{krt_stage['partial_r_CLDN4_immune_given_TACSTD2']:.3f} for CLDN4 after TACSTD2. "
                "Under that keratin specification the two directions are not separated."
            ),
            "",
            (
                "Concordant-4 datasets, fit separately: "
                + gse_bits
                + ". GSE205335 is the one dataset whose BIC winner is not M1, and that gap is small. "
                "The pooled n=65 result is where the ΔBIC exceeds 10."
            ),
        ]
    )

    write_finding(
        ROOT / "FINDING.md",
        c4_primary,
        tc_primary,
        sens_lines,
        cal_lines,
        two_stage_lines,
        reading_lines,
    )

    plot_bic(
        fig_dir / "fig_delta_bic.png",
        [
            ("Concordant-4 % positive", c4p["pooled"]),
            ("TCGA CD8, sum of cohort BICs", tcp["summed"]),
            ("TCGA CD8, pooled z-scores", tcp["pooled"]),
        ],
    )
    plot_mediation(
        fig_dir / "fig_mediation.png",
        [
            ("Concordant-4", c4p["boot_rows"]),
            ("TCGA CD8", tcp["boot_rows"]),
        ],
    )
    plot_partial(
        fig_dir / "fig_two_stage.png",
        [
            {
                "label": "Concordant-4 %pos",
                "r_T_given_C": c4p["stage"]["partial_r_TACSTD2_immune_given_CLDN4"],
                "r_C_given_T": c4p["stage"]["partial_r_CLDN4_immune_given_TACSTD2"],
            },
            {
                "label": "Concordant-4 mean log1p",
                "r_T_given_C": c4_done["c4_mean"]["stage"]["partial_r_TACSTD2_immune_given_CLDN4"],
                "r_C_given_T": c4_done["c4_mean"]["stage"]["partial_r_CLDN4_immune_given_TACSTD2"],
            },
            {
                "label": "TCGA CD8 | KRT8/18/19",
                "r_T_given_C": tcp["stage"]["partial_r_TACSTD2_immune_given_CLDN4"],
                "r_C_given_T": tcp["stage"]["partial_r_CLDN4_immune_given_TACSTD2"],
            },
            {
                "label": "TCGA CD3 | KRT8/18/19",
                "r_T_given_C": tcga_done["tcga_cd3_krt819"]["stage"]["partial_r_TACSTD2_immune_given_CLDN4"],
                "r_C_given_T": tcga_done["tcga_cd3_krt819"]["stage"]["partial_r_CLDN4_immune_given_TACSTD2"],
            },
        ],
    )

    summary = {
        "concordant4_primary_winner": c4p["pooled"]["winner"],
        "concordant4_primary_bic_gap": c4p["pooled"]["bic_gap"],
        "concordant4_sum_winner": c4p["summed"]["winner"],
        "tcga_sum_winner": tcp["summed"]["winner"],
        "tcga_sum_bic_gap": tcp["summed"]["bic_gap"],
        "tcga_pooled_winner": tcp["pooled"]["winner"],
        "tcga_votes": tcp["per_winners"],
        "n_boot": N_BOOT,
        "seed": SEED,
    }
    (table_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

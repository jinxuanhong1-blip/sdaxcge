#!/usr/bin/env python3
"""Study-level REML meta-analysis of TISMO ICB paired pre/post expression.

The locked reference is the 64 TISMO Gene-module slices (cell line x study x
condition x regimen): Tacstd2 rose in 49/64 (Wilcoxon p = 5.84e-5). This
script recomputes that lock and does not replace it. The new analysis pools
those paired mean differences at the study (GEO/ENA accession) level with a
random-effects model whose heterogeneity variance is estimated by REML,
then draws forest plots and leave-one-study-out fits.

TISMO has no KL (Kras/Stk11) or KP (Kras/Trp53) lung ICB model. LLC stays LLC.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "tismo"
VIVO = DATA / "vivo"
OUT = ROOT / "results" / "tismo_icb_study_meta"
FIG = OUT / "figures"
TABLES = OUT / "tables"

TJ_GENES = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
MARKERS = ["Tacstd2", "Cldn4", "TJ"]
STEM_RE = re.compile(r"\(n=\d+\)$")
Z95 = float(stats.norm.ppf(0.975))

# TISMO cellLineMeta cancerType, collapsed. KPB25L is mammary, not lung KP.
CANCER = {
    "4T1": "Mammary",
    "E0771": "Mammary",
    "EMT6": "Mammary",
    "KPB25L": "Mammary",
    "T11": "Mammary",
    "p53-2225L": "Mammary",
    "p53-2336R": "Mammary",
    "B16": "Melanoma",
    "YUMM1.7": "Melanoma",
    "D3UV2": "Melanoma",
    "D4M.3A.3": "Melanoma",
    "CT26": "Colorectal",
    "MC38": "Colorectal",
    "LLC": "Lung",
    "402230": "Sarcoma",
    "BNL-MEA": "Liver",
    "YTN16": "Gastric",
}
CANCER_ORDER = ["Mammary", "Melanoma", "Colorectal", "Gastric", "Liver", "Lung", "Sarcoma", "Other"]
CANCER_COLOR = {
    "Mammary": "#4C78A8",
    "Melanoma": "#F58518",
    "Colorectal": "#54A24B",
    "Gastric": "#E45756",
    "Liver": "#B279A2",
    "Lung": "#111111",
    "Sarcoma": "#72B7B2",
    "Other": "#9D755D",
}

# Locked slice-level reference (PR #542). Tests fail if a recompute drifts.
LOCKED = {
    "Tacstd2": {"n_up": 49, "n_down": 15, "n_tie": 0, "wilcoxon_stat": 439.0, "wilcoxon_p": 5.8398482923347564e-05},
    "Cldn4": {"n_up": 34, "n_down": 28, "n_tie": 2, "wilcoxon_stat": 724.0, "wilcoxon_p": 0.07667785136224303},
    "TJ": {"n_up": 37, "n_down": 27, "n_tie": 0, "wilcoxon_stat": 784.0, "wilcoxon_p": 0.08689644662558298},
}


def read_expression_csv(path: Path, gene: str) -> pd.DataFrame:
    raw = path.read_text(encoding="utf-8", errors="replace")
    reader = csv.reader(io.StringIO(raw))
    header = next(reader)
    good = [row for row in reader if len(row) == len(header)]
    df = pd.DataFrame(good, columns=header)
    df = df[df["geneID"] == gene]
    for col in ("value", "Baseline"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["value", "Baseline"]).reset_index(drop=True)
    df["stem"] = df["cell_line"].map(lambda s: STEM_RE.sub("", str(s)).rstrip())
    df["line"] = df["stem"].map(lambda s: s.split("_")[0])
    return df


def file_md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_num(x, digits: int = 3) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    return f"{x:.{digits}f}"


def sign_test(values: np.ndarray) -> dict:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n_up = int(np.sum(v > 0))
    n_down = int(np.sum(v < 0))
    n_tie = int(np.sum(v == 0))
    n = n_up + n_down
    if n == 0:
        p_two = p_greater = math.nan
    else:
        p_two = float(stats.binomtest(n_up, n, 0.5, alternative="two-sided").pvalue)
        p_greater = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    return {
        "n": int(len(v)),
        "n_up": n_up,
        "n_down": n_down,
        "n_tie": n_tie,
        "n_tested": n,
        "binomial_p_two_sided": p_two,
        "binomial_p_greater": p_greater,
    }


def wilcoxon_signed(values: np.ndarray) -> dict:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    nz = v[v != 0]
    if len(nz) == 0:
        return {"n": 0, "statistic": math.nan, "p_two_sided": math.nan, "p_greater": math.nan}
    two = stats.wilcoxon(nz, alternative="two-sided", zero_method="wilcox")
    greater = stats.wilcoxon(nz, alternative="greater", zero_method="wilcox")
    return {
        "n": int(len(nz)),
        "statistic": float(two.statistic),
        "p_two_sided": float(two.pvalue),
        "p_greater": float(greater.pvalue),
    }


def reml_tau2(y: np.ndarray, v: np.ndarray) -> float:
    """REML estimate of between-study variance for the normal-normal model.

    Maximizes the restricted log-likelihood in Viechtbauer (2005) / metafor.
    Returns 0 when the maximum is on the boundary.
    """
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    if len(y) < 2:
        return 0.0
    if np.any(v <= 0) or not np.all(np.isfinite(v)) or not np.all(np.isfinite(y)):
        raise ValueError("REML inputs must be finite with positive variances")

    def nll(tau2: float) -> float:
        w = 1.0 / (v + tau2)
        sw = np.sum(w)
        mu = np.sum(w * y) / sw
        return 0.5 * (
            float(np.sum(np.log(v + tau2)))
            + math.log(float(sw))
            + float(np.sum(w * (y - mu) ** 2))
        )

    span = float(np.max(y) - np.min(y))
    upper = max(1.0, span ** 2 * 10.0, float(np.max(v)) * 100.0)
    grid = np.unique(
        np.concatenate(
            [
                [0.0],
                np.linspace(0.0, upper, 500),
                np.geomspace(max(upper * 1e-8, 1e-12), upper, 200),
            ]
        )
    )
    vals = np.array([nll(float(t)) for t in grid])
    tau_grid = float(grid[int(np.argmin(vals))])
    # Local polish around the grid minimum.
    from scipy.optimize import minimize_scalar

    lo = 0.0
    hi = upper
    bracket_left = max(0.0, tau_grid - upper / 50.0)
    res = minimize_scalar(nll, bounds=(lo, hi), method="bounded", options={"xatol": 1e-14})
    candidates = [0.0, tau_grid, float(res.x)]
    # Keep the search inside a neighborhood of the grid hit as well.
    if bracket_left < hi:
        res2 = minimize_scalar(
            nll, bounds=(bracket_left, min(hi, tau_grid + upper / 50.0)), method="bounded", options={"xatol": 1e-14}
        )
        candidates.append(float(res2.x))
    best = min(candidates, key=nll)
    if nll(0.0) <= nll(best) + 1e-10:
        return 0.0
    return float(best)


def der_simonian_laird(y: np.ndarray, v: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    k = len(y)
    if k < 2:
        return 0.0
    w = 1.0 / v
    mu = np.sum(w * y) / np.sum(w)
    q = float(np.sum(w * (y - mu) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    if c <= 0:
        return 0.0
    return max(0.0, (q - (k - 1)) / c)


def fit_meta(y: np.ndarray, v: np.ndarray, method: str = "REML") -> dict:
    """Random-effects meta-analysis. method is REML or DL for tau^2."""
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    k = int(len(y))
    if k == 0:
        raise ValueError("no studies")
    w_fe = 1.0 / v
    mu_fe = float(np.sum(w_fe * y) / np.sum(w_fe))
    se_fe = float(1.0 / math.sqrt(float(np.sum(w_fe))))
    q = float(np.sum(w_fe * (y - mu_fe) ** 2))
    df = k - 1
    p_q = float(stats.chi2.sf(q, df)) if df > 0 else math.nan
    i2 = float(max(0.0, (q - df) / q)) if q > 0 and df > 0 else 0.0
    h2 = float(q / df) if df > 0 else math.nan
    if method == "REML":
        tau2 = reml_tau2(y, v) if k >= 2 else 0.0
    elif method == "DL":
        tau2 = der_simonian_laird(y, v)
    else:
        raise ValueError(method)
    w = 1.0 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    se = float(1.0 / math.sqrt(float(np.sum(w))))
    z = mu / se if se > 0 else math.nan
    p_z = float(2 * stats.norm.sf(abs(z))) if se > 0 else math.nan
    ci_low = mu - Z95 * se
    ci_high = mu + Z95 * se
    # Knapp-Hartung, plus a modified version that never undercuts the Wald SE.
    # The slice-replication variance is a conservative plug-in. Raw KH can then
    # shrink the interval. Modified KH (scale floored at 1, t critical value)
    # does not turn that conservatism into a tighter claim.
    q_w = float(np.sum(w * (y - mu) ** 2))
    if k >= 2:
        tcrit = float(stats.t.ppf(0.975, k - 1))
        scale = math.sqrt(q_w / (k - 1)) if q_w > 0 else 0.0
        se_kh = se * scale
        if se_kh > 0:
            kh_low = mu - tcrit * se_kh
            kh_high = mu + tcrit * se_kh
            p_kh = float(2 * stats.t.sf(abs(mu / se_kh), k - 1))
        else:
            kh_low = kh_high = mu
            p_kh = 0.0 if mu != 0 else 1.0
        se_mkh = se * max(1.0, scale)
        mkh_low = mu - tcrit * se_mkh
        mkh_high = mu + tcrit * se_mkh
        p_mkh = float(2 * stats.t.sf(abs(mu / se_mkh), k - 1)) if se_mkh > 0 else math.nan
    else:
        se_kh = kh_low = kh_high = p_kh = math.nan
        se_mkh = mkh_low = mkh_high = p_mkh = math.nan
    if k >= 2:
        tcrit = float(stats.t.ppf(0.975, k - 1))
        se_pred = math.sqrt(tau2 + se ** 2)
        pi_low = mu - tcrit * se_pred
        pi_high = mu + tcrit * se_pred
    else:
        pi_low = pi_high = math.nan
    weights = w / np.sum(w)
    return {
        "k": k,
        "method": method,
        "mu": mu,
        "se": se,
        "z": float(z),
        "p_wald": p_z,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "tau2": float(tau2),
        "i2": i2,
        "h2": h2,
        "Q": q,
        "Q_df": df,
        "Q_p": p_q,
        "mu_fe": mu_fe,
        "se_fe": se_fe,
        "fe_low": mu_fe - Z95 * se_fe,
        "fe_high": mu_fe + Z95 * se_fe,
        "se_kh": float(se_kh) if se_kh == se_kh else math.nan,
        "kh_low": float(kh_low) if kh_low == kh_low else math.nan,
        "kh_high": float(kh_high) if kh_high == kh_high else math.nan,
        "p_kh": float(p_kh) if p_kh == p_kh else math.nan,
        "se_mkh": float(se_mkh) if se_mkh == se_mkh else math.nan,
        "mkh_low": float(mkh_low) if mkh_low == mkh_low else math.nan,
        "mkh_high": float(mkh_high) if mkh_high == mkh_high else math.nan,
        "p_mkh": float(p_mkh) if p_mkh == p_mkh else math.nan,
        "pi_low": float(pi_low) if pi_low == pi_low else math.nan,
        "pi_high": float(pi_high) if pi_high == pi_high else math.nan,
        "weights": weights,
    }


def leave_one_out(y: np.ndarray, v: np.ndarray, labels: list[str]) -> list[dict]:
    rows = []
    for i, label in enumerate(labels):
        mask = np.ones(len(y), dtype=bool)
        mask[i] = False
        fit = fit_meta(y[mask], v[mask], method="REML")
        rows.append(
            {
                "left_out": label,
                "k": fit["k"],
                "mu": fit["mu"],
                "se": fit["se"],
                "ci_low": fit["ci_low"],
                "ci_high": fit["ci_high"],
                "p_wald": fit["p_wald"],
                "tau2": fit["tau2"],
                "i2": fit["i2"],
                "kh_low": fit["kh_low"],
                "kh_high": fit["kh_high"],
                "p_kh": fit["p_kh"],
                "Q_p": fit["Q_p"],
            }
        )
    return rows


def _ids(series: pd.Series) -> frozenset:
    return frozenset(series.astype(str).tolist())


def cov_factor(a: dict, b: dict) -> float:
    """Coefficient of sigma^2 in Cov(delta_a, delta_b) under independent mice."""

    def cm(ids1: frozenset, ids2: frozenset) -> float:
        if not ids1 or not ids2:
            return 0.0
        return len(ids1 & ids2) / (len(ids1) * len(ids2))

    return (
        cm(a["icb_ids"], b["icb_ids"])
        + cm(a["base_ids"], b["base_ids"])
        - cm(a["icb_ids"], b["base_ids"])
        - cm(a["base_ids"], b["icb_ids"])
    )


def slice_table(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    rows = []
    for stem, g in df.groupby("stem", sort=True):
        base = g.loc[g["Baseline"] == 1]
        icb = g.loc[g["Baseline"] == 0]
        if len(base) == 0 or len(icb) == 0:
            continue
        bval = base[value_col].astype(float)
        tval = icb[value_col].astype(float)
        line = str(g["line"].iloc[0])
        rows.append(
            {
                "stem": stem,
                "line": line,
                "gse_id": str(g["GSE_ID"].iloc[0]),
                "cancer_group": CANCER.get(line, "Other"),
                "is_llc": line.upper() == "LLC",
                "is_lung": CANCER.get(line, "Other") == "Lung",
                "n_baseline": int(len(bval)),
                "n_icb": int(len(tval)),
                "delta_mean": float(tval.mean() - bval.mean()),
                "mean_baseline": float(bval.mean()),
                "mean_icb": float(tval.mean()),
                "sd_baseline": float(bval.std(ddof=1)) if len(bval) >= 2 else math.nan,
                "sd_icb": float(tval.std(ddof=1)) if len(tval) >= 2 else math.nan,
                "base_ids": _ids(base["Samples"]),
                "icb_ids": _ids(icb["Samples"]),
                "base_values": bval.to_numpy(dtype=float),
                "icb_values": tval.to_numpy(dtype=float),
            }
        )
    return pd.DataFrame(rows)


def borrowed_arm_variance(slices: pd.DataFrame) -> float:
    vars_ = []
    for _, row in slices.iterrows():
        for sd, n in (
            (row["sd_baseline"], row["n_baseline"]),
            (row["sd_icb"], row["n_icb"]),
        ):
            if n >= 2 and np.isfinite(sd) and sd > 0:
                vars_.append(float(sd) ** 2)
    if not vars_:
        raise RuntimeError("no positive within-arm variance to borrow")
    return float(np.median(vars_))


def pooled_slice_variance(slices: pd.DataFrame) -> float:
    """Pooled variance of slice deltas inside studies that have at least two slices.

    This is the replication variance of a locked contrast, not the mouse-level
    technical variance. Dividing by k_i is the variance of that study's mean.
    """
    num = 0.0
    den = 0
    for _, g in slices.groupby("gse_id"):
        k = len(g)
        if k < 2:
            continue
        num += (k - 1) * float(np.var(g["delta_mean"].to_numpy(dtype=float), ddof=1))
        den += k - 1
    if den == 0:
        raise RuntimeError("no multi-slice study to estimate slice-replication variance")
    return num / den


def pooled_sigma2(group: pd.DataFrame) -> tuple[float, int]:
    num = 0.0
    den = 0
    for _, row in group.iterrows():
        for vals in (row["base_values"], row["icb_values"]):
            if len(vals) >= 2:
                num += (len(vals) - 1) * float(np.var(vals, ddof=1))
                den += len(vals) - 1
    if den == 0:
        return 0.0, 0
    return num / den, den


def study_estimates(slices: pd.DataFrame, borrowed: float) -> pd.DataFrame:
    """One row per study.

    Primary variance uses the shared-mouse covariance of the slice deltas and
    the unweighted mean of those deltas (each locked slice counts once).
    """
    rows = []
    for gse, g in slices.groupby("gse_id", sort=True):
        g = g.reset_index(drop=True)
        k = len(g)
        deltas = g["delta_mean"].to_numpy(dtype=float)
        sigma2, df_sigma = pooled_sigma2(g)
        borrowed_used = False
        if not np.isfinite(sigma2) or sigma2 <= 0:
            sigma2 = borrowed
            borrowed_used = True
        records = g.to_dict(orient="records")
        sig = np.zeros((k, k), dtype=float)
        for j in range(k):
            for m in range(k):
                sig[j, m] = sigma2 * cov_factor(records[j], records[m])
        # Numerical floor so a singular shared-sample matrix cannot zero a weight.
        eig = np.linalg.eigvalsh(sig)
        if eig.min() < 1e-15:
            sig = sig + np.eye(k) * (1e-12 - min(0.0, float(eig.min())))
        ones = np.ones(k)
        y = float(np.mean(deltas))
        v = float(ones @ sig @ ones) / (k ** 2)
        if v <= 0 or not np.isfinite(v):
            v = borrowed / max(k, 1)
            borrowed_used = True
        # GLS (precision-weighted) companion inside the study.
        try:
            sig_inv = np.linalg.inv(sig)
            prec = float(ones @ sig_inv @ ones)
            y_gls = float(ones @ sig_inv @ deltas) / prec
            v_gls = 1.0 / prec
        except np.linalg.LinAlgError:
            y_gls, v_gls = y, v
        # Independence (off-diagonal zeroed) for the same unweighted mean.
        sig_ind = np.diag(np.diag(sig))
        v_ind = float(ones @ sig_ind @ ones) / (k ** 2)
        if k >= 2:
            v_disp = float(np.var(deltas, ddof=1)) / k
        else:
            v_disp = v
        v_disp = max(v_disp, v)
        lines = sorted(set(g["line"]))
        groups = sorted(set(g["cancer_group"]))
        cancer = groups[0] if len(groups) == 1 else "mixed"
        sample_ids = set()
        for rec in records:
            sample_ids |= set(rec["base_ids"]) | set(rec["icb_ids"])
        n_shared_base_pairs = 0
        for j in range(k):
            for m in range(j + 1, k):
                if records[j]["base_ids"] & records[m]["base_ids"]:
                    n_shared_base_pairs += 1
        rows.append(
            {
                "gse_id": gse,
                "cancer_group": cancer,
                "lines": ",".join(lines),
                "is_llc": bool(g["is_llc"].any()),
                "is_lung": bool(g["is_lung"].any()),
                "n_slices": k,
                "n_mice": len(sample_ids),
                "n_baseline_mice": len(set().union(*[set(r["base_ids"]) for r in records])),
                "n_icb_mice": len(set().union(*[set(r["icb_ids"]) for r in records])),
                "sigma2": sigma2,
                "sigma2_df": df_sigma,
                "borrowed_sigma2": borrowed_used,
                "n_shared_baseline_pairs": n_shared_base_pairs,
                "y": y,
                "v": v,
                "y_gls": float(y_gls),
                "v_gls": float(v_gls),
                "v_independence": float(v_ind),
                "v_dispersion": float(v_disp),
                "slice_delta_min": float(np.min(deltas)),
                "slice_delta_max": float(np.max(deltas)),
            }
        )
    out = pd.DataFrame(rows)
    out["cancer_rank"] = out["cancer_group"].map({c: i for i, c in enumerate(CANCER_ORDER)}).fillna(99)
    return out.sort_values(["cancer_rank", "gse_id"]).reset_index(drop=True)


def load_locked_slices() -> dict[str, pd.DataFrame]:
    tac = read_expression_csv(VIVO / "Tacstd2.csv", "Tacstd2")
    lock_samples = pd.Index(tac["Samples"].unique())
    frames = []
    for gene in ["Tacstd2", *TJ_GENES]:
        df = read_expression_csv(VIVO / f"{gene}.csv", gene)
        part = (
            df[["Samples", "value"]]
            .drop_duplicates("Samples")
            .rename(columns={"value": gene})
            .set_index("Samples")
        )
        frames.append(part)
    mat = frames[0]
    for frame in frames[1:]:
        mat = mat.join(frame, how="outer")
    mat = mat.reindex(lock_samples)
    tac = tac.copy()
    tac["Cldn4_value"] = tac["Samples"].map(mat["Cldn4"])
    tac["TJ"] = tac["Samples"].map(mat[TJ_GENES].mean(axis=1, skipna=True))
    tac["n_tj_genes"] = tac["Samples"].map(mat[TJ_GENES].notna().sum(axis=1))
    lock_stems = set(tac["stem"])
    slices = {
        "Tacstd2": slice_table(tac, "value"),
        "Cldn4": slice_table(tac.dropna(subset=["Cldn4_value"]).assign(value=lambda d: d["Cldn4_value"]), "value"),
        "TJ": slice_table(tac.dropna(subset=["TJ"]).assign(value=lambda d: d["TJ"]), "value"),
    }
    for marker, table in slices.items():
        slices[marker] = table[table["stem"].isin(lock_stems)].copy()
    meta = {
        "n_tacstd2_samples": int(tac["Samples"].nunique()),
        "mean_tj_genes_per_sample": float(tac["n_tj_genes"].mean()),
        "min_tj_genes_per_sample": int(tac["n_tj_genes"].min()),
    }
    return {"slices": slices, "sample_meta": meta}


def locked_reference_row(marker: str, slices: pd.DataFrame) -> dict:
    signs = sign_test(slices["delta_mean"].to_numpy())
    wil = wilcoxon_signed(slices["delta_mean"].to_numpy())
    expected = LOCKED[marker]
    return {
        "marker": marker,
        "unit": "slice",
        "n_slices": signs["n"],
        "n_up": signs["n_up"],
        "n_down": signs["n_down"],
        "n_tie": signs["n_tie"],
        "binomial_p_two_sided": signs["binomial_p_two_sided"],
        "wilcoxon_n": wil["n"],
        "wilcoxon_stat": wil["statistic"],
        "wilcoxon_p_two_sided": wil["p_two_sided"],
        "mean_delta": float(slices["delta_mean"].mean()),
        "median_delta": float(slices["delta_mean"].median()),
        "matches_lock": (
            signs["n_up"] == expected["n_up"]
            and signs["n_down"] == expected["n_down"]
            and signs["n_tie"] == expected["n_tie"]
            and abs(wil["statistic"] - expected["wilcoxon_stat"]) < 1e-8
            and abs(wil["p_two_sided"] - expected["wilcoxon_p"]) < 1e-12
        ),
    }


def kl_audit() -> dict:
    cell = json.loads((DATA / "cellLineMeta.json").read_text())
    vivo = json.loads((DATA / "vivoMeta.json").read_text())
    cell_blob = json.dumps(cell).lower()
    vivo_blob = json.dumps(vivo).lower()
    lung_lines = [
        {"cellLine": r.get("cellLine"), "cancerType": r.get("cancerType"), "background": r.get("background")}
        for r in cell
        if "lung" in str(r.get("cancerType", "")).lower()
    ]
    return {
        "n_cell_lines": len(cell),
        "n_vivo_annotation_rows": len(vivo),
        "cellLineMeta_mentions_stk11": "stk11" in cell_blob,
        "cellLineMeta_mentions_lkb1": "lkb1" in cell_blob,
        "vivoMeta_mentions_stk11": "stk11" in vivo_blob,
        "vivoMeta_mentions_lkb1": "lkb1" in vivo_blob,
        "lung_lines_in_catalog": lung_lines,
        "kpb25l_cancer": next((r.get("cancerType") for r in cell if r.get("cellLine") == "KPB25L"), None),
        "kpc_cancer": next((r.get("cancerType") for r in cell if r.get("cellLine") == "KPC"), None),
        "llc_cancer": next((r.get("cancerType") for r in cell if r.get("cellLine") == "LLC"), None),
        "note": (
            "No STK11/LKB1 annotation in TISMO cellLineMeta or vivoMeta. "
            "LLC is lung carcinoma, not KL. KPB25L is mammary. KPC is pancreatic."
        ),
    }


def attach_weights(studies: pd.DataFrame, fit: dict, prefix: str = "") -> pd.DataFrame:
    out = studies.copy()
    out[f"{prefix}weight"] = fit["weights"]
    out[f"{prefix}weight_pct"] = 100.0 * fit["weights"]
    return out


def analyze_marker(marker: str, slices: pd.DataFrame) -> dict:
    ref = locked_reference_row(marker, slices)
    if not ref["matches_lock"]:
        raise SystemExit(f"{marker} drifted from the locked 64-slice reference: {ref}")
    borrowed = borrowed_arm_variance(slices)
    sigma_slice = pooled_slice_variance(slices)
    studies = study_estimates(slices, borrowed)
    # Technical mouse-level variance lives in `v`. Primary variance is the
    # larger of that and the pooled slice-replication variance / n_slices.
    # The max() stops an unexpressed gene (technical variance near 0) from
    # taking essentially all of the weight.
    studies["v_technical"] = studies["v"]
    studies["v_slice"] = sigma_slice / studies["n_slices"].to_numpy(dtype=float)
    studies["v_primary"] = np.maximum(studies["v_technical"].to_numpy(dtype=float), studies["v_slice"].to_numpy(dtype=float))
    y = studies["y"].to_numpy(dtype=float)
    v = studies["v_primary"].to_numpy(dtype=float)
    fit = fit_meta(y, v, method="REML")
    fit_dl = fit_meta(y, v, method="DL")
    fit_technical = fit_meta(y, studies["v_technical"].to_numpy(dtype=float), method="REML")
    fit_slice_only = fit_meta(y, studies["v_slice"].to_numpy(dtype=float), method="REML")
    fit_gls = fit_meta(studies["y_gls"].to_numpy(dtype=float), studies["v_gls"].to_numpy(dtype=float), method="REML")
    fit_ind = fit_meta(y, studies["v_independence"].to_numpy(dtype=float), method="REML")
    studies = attach_weights(studies, fit)
    studies["weight_technical"] = fit_technical["weights"]
    studies["weight_technical_pct"] = 100.0 * fit_technical["weights"]
    signs = sign_test(y)
    wil = wilcoxon_signed(y)
    # One-study-one-vote: Student interval on the 22 study means.
    mean_vote = float(np.mean(y))
    median_vote = float(np.median(y))
    sd_vote = float(np.std(y, ddof=1))
    se_vote = sd_vote / math.sqrt(len(y))
    tcrit = float(stats.t.ppf(0.975, len(y) - 1))
    tstat = mean_vote / se_vote
    p_vote = float(2 * stats.t.sf(abs(tstat), len(y) - 1))
    loso = leave_one_out(y, v, studies["gse_id"].tolist())
    full_mu = fit["mu"]
    for row in loso:
        row["delta_mu"] = row["mu"] - full_mu
        row["wald_excludes_zero"] = bool(row["ci_low"] > 0 or row["ci_high"] < 0)
        row["marker"] = marker
    subgroups = []
    for cancer, sub in studies.groupby("cancer_group"):
        if len(sub) < 3:
            continue
        sub_fit = fit_meta(sub["y"].to_numpy(dtype=float), sub["v_primary"].to_numpy(dtype=float), method="REML")
        sub_signs = sign_test(sub["y"].to_numpy(dtype=float))
        subgroups.append(
            {
                "marker": marker,
                "cancer_group": cancer,
                "k": sub_fit["k"],
                "n_up": sub_signs["n_up"],
                "n_down": sub_signs["n_down"],
                "n_tie": sub_signs["n_tie"],
                "binomial_p_two_sided": sub_signs["binomial_p_two_sided"],
                "mu": sub_fit["mu"],
                "se": sub_fit["se"],
                "ci_low": sub_fit["ci_low"],
                "ci_high": sub_fit["ci_high"],
                "p_wald": sub_fit["p_wald"],
                "tau2": sub_fit["tau2"],
                "i2": sub_fit["i2"],
                "mkh_low": sub_fit["mkh_low"],
                "mkh_high": sub_fit["mkh_high"],
                "p_mkh": sub_fit["p_mkh"],
            }
        )
    lung = studies[studies["is_lung"]]
    tech_idx = studies["weight_technical"].idxmax()
    return {
        "marker": marker,
        "reference": ref,
        "borrowed_arm_variance": borrowed,
        "sigma_slice": sigma_slice,
        "n_studies_borrowed_sigma2": int(studies["borrowed_sigma2"].sum()),
        "studies": studies,
        "fit": fit,
        "fit_dl": fit_dl,
        "fit_technical": fit_technical,
        "fit_slice_only": fit_slice_only,
        "fit_gls": fit_gls,
        "fit_independence": fit_ind,
        "study_signs": signs,
        "study_wilcoxon": wil,
        "one_vote": {
            "mean": mean_vote,
            "median": median_vote,
            "se": se_vote,
            "ci_low": mean_vote - tcrit * se_vote,
            "ci_high": mean_vote + tcrit * se_vote,
            "p_two_sided": p_vote,
        },
        "loso": loso,
        "subgroups": subgroups,
        "lung_studies": lung.to_dict(orient="records"),
        "max_weight_study": str(studies.loc[studies["weight"].idxmax(), "gse_id"]),
        "max_weight_pct": float(studies["weight"].max() * 100.0),
        "gse124821_weight_pct": float(studies.loc[studies["gse_id"] == "GSE124821", "weight"].iloc[0] * 100.0)
        if (studies["gse_id"] == "GSE124821").any()
        else math.nan,
        "technical_max_weight_study": str(studies.loc[tech_idx, "gse_id"]),
        "technical_max_weight_pct": float(studies.loc[tech_idx, "weight_technical"] * 100.0),
    }


def _ci_text(fit: dict) -> str:
    return (
        f"{fmt_num(fit['mu'])} (Wald 95% CI {fmt_num(fit['ci_low'])} to {fmt_num(fit['ci_high'])}; "
        f"modified Knapp–Hartung {fmt_num(fit['mkh_low'])} to {fmt_num(fit['mkh_high'])})"
    )


def write_forest(marker: str, block: dict, path: Path) -> None:
    studies = block["studies"]
    fit = block["fit"]
    k = len(studies)
    # Display top-to-bottom in cancer order (frame is already sorted).
    y_pos = np.arange(k, 0, -1)
    fig_h = max(6.8, 0.38 * k + 2.4)
    fig = plt.figure(figsize=(11.6, fig_h))
    ax = fig.add_axes([0.34, 0.10, 0.38, 0.78])
    ax.axvline(0, color="0.35", lw=0.8, zorder=0)
    for pos, (_, row) in zip(y_pos, studies.iterrows()):
        color = CANCER_COLOR.get(row["cancer_group"], "#333333")
        se = math.sqrt(row["v_primary"])
        ax.plot([row["y"] - Z95 * se, row["y"] + Z95 * se], [pos, pos], color=color, lw=1.3, solid_capstyle="round")
        size = 18 + 220 * float(row["weight"])
        edge = "#111111" if row["is_llc"] else color
        lw = 1.4 if row["is_llc"] else 0.4
        ax.scatter([row["y"]], [pos], s=size, color=color, edgecolors=edge, linewidths=lw, zorder=3, marker="s")
    # Pooled diamond and prediction interval.
    diamond_y = -0.15
    pi_y = -1.05
    ax.plot([fit["pi_low"], fit["pi_high"]], [pi_y, pi_y], color="0.45", lw=4, solid_capstyle="butt", alpha=0.35, zorder=2)
    _diamond(ax, fit["mu"], fit["mkh_low"], fit["mkh_high"], diamond_y, "#1F4E79")
    labels = []
    for _, row in studies.iterrows():
        llc = " LLC, not KL" if row["is_llc"] else ""
        labels.append(f"{row['gse_id']}  {row['cancer_group']}{llc}\n{row['lines']}  ·  {int(row['n_slices'])} slices, {int(row['n_mice'])} mice")
    ax.set_yticks(list(y_pos) + [diamond_y, pi_y])
    ax.set_yticklabels(labels + ["REML pooled (mKH CI)", "95% prediction interval"], fontsize=7.5)
    ax.set_ylim(-1.8, k + 1.15)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Mean difference after ICB (TISMO log2(TPM+1))")
    ref = block["reference"]
    ax.set_title(
        f"{marker}  ·  study-level REML\n"
        f"locked reference {ref['n_up']}/{ref['n_slices']} slices up, Wilcoxon p={fmt_p(ref['wilcoxon_p_two_sided'])}",
        loc="left",
        fontsize=11,
    )
    # Right-hand estimates.
    x_right = ax.get_xlim()[1]
    # Fix xlim from data first, then place text in axes fraction via a second pass.
    xs = list(studies["y"]) + [fit["pi_low"], fit["pi_high"], fit["mkh_low"], fit["mkh_high"]]
    pad = 0.08 * (max(xs) - min(xs) + 1e-6)
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    fig.canvas.draw()
    for pos, (_, row) in zip(y_pos, studies.iterrows()):
        se = math.sqrt(row["v_primary"])
        txt = f"{row['y']:+.2f}  [{row['y'] - Z95 * se:+.2f}, {row['y'] + Z95 * se:+.2f}]   {row['weight_pct']:.1f}%"
        ax.annotate(
            txt,
            xy=(1.02, pos),
            xycoords=("axes fraction", "data"),
            va="center",
            ha="left",
            fontsize=7,
            annotation_clip=False,
            color="0.15",
        )
    ax.annotate(
        f"{fit['mu']:+.2f}  [{fit['mkh_low']:+.2f}, {fit['mkh_high']:+.2f}]   I²={100 * fit['i2']:.0f}%",
        xy=(1.02, diamond_y),
        xycoords=("axes fraction", "data"),
        va="center",
        ha="left",
        fontsize=7.5,
        fontweight="regular",
        annotation_clip=False,
    )
    ax.annotate(
        f"[{fit['pi_low']:+.2f}, {fit['pi_high']:+.2f}]",
        xy=(1.02, pi_y),
        xycoords=("axes fraction", "data"),
        va="center",
        ha="left",
        fontsize=7.5,
        annotation_clip=False,
        color="0.35",
    )
    ax.annotate(
        "estimate [95% CI]    weight",
        xy=(1.02, k + 0.72),
        xycoords=("axes fraction", "data"),
        va="center",
        ha="left",
        fontsize=7.5,
        annotation_clip=False,
        color="0.25",
    )
    fig.text(
        0.34,
        0.015,
        "Square area is the random-effects weight. Black edge marks LLC (GSE155972), the only lung ICB study. LLC is not KL.",
        fontsize=7.5,
        color="0.25",
    )
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _diamond(ax, mu, lo, hi, y, color) -> None:
    height = 0.38
    poly = plt.Polygon(
        [[mu, y + height], [hi, y], [mu, y - height], [lo, y]],
        closed=True,
        facecolor=color,
        edgecolor=color,
        zorder=4,
    )
    ax.add_patch(poly)


def write_panel(blocks: dict[str, dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 7.6), sharey=True)
    # Common study order from Tacstd2 (same 22 studies).
    order = blocks["Tacstd2"]["studies"]["gse_id"].tolist()
    y_pos = np.arange(len(order), 0, -1)
    for ax, marker in zip(axes, MARKERS):
        block = blocks[marker]
        studies = block["studies"].set_index("gse_id").loc[order].reset_index()
        fit = block["fit"]
        ax.axvline(0, color="0.4", lw=0.7)
        for pos, row in zip(y_pos, studies.itertuples(index=False)):
            color = CANCER_COLOR.get(row.cancer_group, "#333")
            se = math.sqrt(row.v_primary)
            ax.plot([row.y - Z95 * se, row.y + Z95 * se], [pos, pos], color=color, lw=1.15)
            ax.scatter([row.y], [pos], s=12 + 140 * row.weight, marker="s", color=color, zorder=3, linewidths=0)
        _diamond(ax, fit["mu"], fit["mkh_low"], fit["mkh_high"], -0.2, "#1F4E79")
        ax.set_title(
            f"{marker}\n{fmt_num(fit['mu'])} [{fmt_num(fit['mkh_low'])}, {fmt_num(fit['mkh_high'])}]",
            fontsize=10,
        )
        ax.set_xlabel("Δ after ICB")
        xs = []
        for row in studies.itertuples(index=False):
            se = math.sqrt(row.v_primary)
            xs.extend([row.y - Z95 * se, row.y + Z95 * se])
        xs.extend([fit["mkh_low"], fit["mkh_high"]])
        pad = 0.06 * (max(xs) - min(xs))
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(-1.3, len(order) + 0.8)
    labels = []
    tac = blocks["Tacstd2"]["studies"].set_index("gse_id")
    for gse in order:
        row = tac.loc[gse]
        tag = " LLC" if row["is_llc"] else ""
        labels.append(f"{gse}{tag}")
    axes[0].set_yticks(list(y_pos) + [-0.2])
    axes[0].set_yticklabels(labels + ["REML"], fontsize=7)
    fig.suptitle("TISMO ICB paired pre/post, one estimate per study (REML weights)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_loso(blocks: dict[str, dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 7.4), sharey=True)
    order = [row["left_out"] for row in blocks["Tacstd2"]["loso"]]
    y_pos = np.arange(len(order), 0, -1)
    for ax, marker in zip(axes, MARKERS):
        by = {row["left_out"]: row for row in blocks[marker]["loso"]}
        fit = blocks[marker]["fit"]
        ax.axvline(fit["mu"], color="#1F4E79", lw=1.0, ls="--")
        ax.axvline(0, color="0.45", lw=0.7)
        for pos, gse in zip(y_pos, order):
            row = by[gse]
            ax.plot([row["ci_low"], row["ci_high"]], [pos, pos], color="0.35", lw=1.1)
            ax.scatter([row["mu"]], [pos], s=18, color="#1F4E79", zorder=3)
        ax.set_title(marker, fontsize=10)
        ax.set_xlabel("REML μ with study left out")
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels([f"without {gse}" for gse in order], fontsize=7)
    fig.suptitle("Leave-one-study-out REML (dashed line = full-data estimate)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_pooled_comparison(blocks: dict[str, dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.axvline(0, color="0.4", lw=0.8)
    ypos = [3, 2, 1]
    for y, marker in zip(ypos, MARKERS):
        fit = blocks[marker]["fit"]
        ax.plot([fit["pi_low"], fit["pi_high"]], [y - 0.18, y - 0.18], color="0.65", lw=6, solid_capstyle="butt", alpha=0.5)
        ax.plot([fit["mkh_low"], fit["mkh_high"]], [y + 0.18, y + 0.18], color="#F2C14E", lw=2.0)
        _diamond(ax, fit["mu"], fit["ci_low"], fit["ci_high"], y, "#1F4E79")
        ax.text(
            fit["ci_high"],
            y + 0.42,
            f"  {fit['mu']:+.3f}   Wald p={fmt_p(fit['p_wald'])}   mKH p={fmt_p(fit['p_mkh'])}   I²={100 * fit['i2']:.0f}%",
            va="bottom",
            fontsize=8,
            color="0.2",
        )
    ax.set_yticks(ypos)
    ax.set_yticklabels(MARKERS)
    ax.set_xlabel("Pooled mean difference after ICB (log2(TPM+1) units)")
    ax.set_ylim(0.3, 3.9)
    ax.set_title("REML summary: Wald CI (diamond), modified Knapp–Hartung CI (gold), prediction interval (gray)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def flat_fit(marker: str, name: str, fit: dict) -> dict:
    return {
        "marker": marker,
        "model": name,
        "k": fit["k"],
        "mu": fit["mu"],
        "se": fit["se"],
        "ci_low": fit["ci_low"],
        "ci_high": fit["ci_high"],
        "p_wald": fit["p_wald"],
        "tau2": fit["tau2"],
        "i2": fit["i2"],
        "h2": fit["h2"],
        "Q": fit["Q"],
        "Q_df": fit["Q_df"],
        "Q_p": fit["Q_p"],
        "mu_fe": fit["mu_fe"],
        "se_fe": fit["se_fe"],
        "fe_low": fit["fe_low"],
        "fe_high": fit["fe_high"],
        "kh_low": fit["kh_low"],
        "kh_high": fit["kh_high"],
        "p_kh": fit["p_kh"],
        "mkh_low": fit["mkh_low"],
        "mkh_high": fit["mkh_high"],
        "p_mkh": fit["p_mkh"],
        "pi_low": fit["pi_low"],
        "pi_high": fit["pi_high"],
    }


def render_results(blocks: dict[str, dict], audit: dict, sample_meta: dict) -> str:
    tac, cld, tj = blocks["Tacstd2"], blocks["Cldn4"], blocks["TJ"]
    tr, cr, jr = tac["reference"], cld["reference"], tj["reference"]
    tf, cf, jf = tac["fit"], cld["fit"], tj["fit"]
    ts, cs, js = tac["study_signs"], cld["study_signs"], tj["study_signs"]

    def loso_span(block: dict) -> str:
        mus = [row["mu"] for row in block["loso"]]
        return f"{min(mus):+.3f} to {max(mus):+.3f}"

    def crossed(block: dict) -> list[str]:
        """Studies whose removal lets the Wald CI cover 0, if the full CI does not, or vice versa."""
        out = []
        full_pos = block["fit"]["ci_low"] > 0
        for row in block["loso"]:
            pos = row["ci_low"] > 0
            if pos != full_pos:
                out.append(row["left_out"])
        return out

    lines = []
    a = lines.append
    a("# TISMO ICB paired pre/post: study-level REML meta-analysis")
    a("")
    a("Public TISMO Gene-module tables (Zeng et al., *Nucleic Acids Research* 2022, PMID 34534350). Expression is the served log2(TPM+1) scale. Every count below is recomputed in `results/tismo_icb_study_meta/tables/`.")
    a("")
    a("## Locked reference (not replaced)")
    a("")
    a("The pairing unit stays the TISMO slice: one cell-line × study × condition × ICB-regimen label. Baseline is `Baseline == 1`. ICB is `Baseline == 0` (responders and non-responders pooled). Delta is mean(ICB) − mean(baseline). Cldn4 and the TJ score are locked to the same 64 Tacstd2 slices. The TJ score is the per-sample mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, and Ocln "
      f"(mean {sample_meta['mean_tj_genes_per_sample']:.2f} of 7 genes present; minimum {sample_meta['min_tj_genes_per_sample']}).")
    a("")
    a("| Marker | Up | Down | Tie | Wilcoxon *p* | Mean Δ |")
    a("|---|---:|---:|---:|---:|---:|")
    for block in (tac, cld, tj):
        r = block["reference"]
        a(
            f"| {r['marker']} | **{r['n_up']}/{r['n_slices']}** | {r['n_down']} | {r['n_tie']} | "
            f"{fmt_p(r['wilcoxon_p_two_sided'])} | {r['mean_delta']:+.3f} |"
        )
    a("")
    a(
        f"**Tacstd2 49/64 up, Wilcoxon *p* = {fmt_p(tr['wilcoxon_p_two_sided'])}**, is the reference sentence. "
        f"Cldn4 is weaker on that same lock: **{cr['n_up']}/{cr['n_slices']} up, {cr['n_down']} down, {cr['n_tie']} ties, *p* = {fmt_p(cr['wilcoxon_p_two_sided'])}**. "
        f"TJ is **{jr['n_up']}/{jr['n_slices']}**, *p* = {fmt_p(jr['wilcoxon_p_two_sided'])}. "
        "This meta-analysis does not revise those three rows."
    )
    a("")
    a("## Why pool by study")
    a("")
    n_studies = int(tac["fit"]["k"])
    gse_slices = int(tac["studies"].loc[tac["studies"]["gse_id"] == "GSE124821", "n_slices"].iloc[0])
    a(
        f"The 64 slices sit in **{n_studies} studies**. GSE124821 alone is {gse_slices}/64 slices, so a slice-level Wilcoxon treats many contrasts from one mammary experiment as independent. "
        "Fourteen sample IDs in GSE130472 are shared baselines across anti-PD-L1 and anti-CTLA4 arms (old and young 4T1). No other study reuses a sample ID across slices."
    )
    a("")
    a(
        "Primary study effect: unweighted mean of the locked slice deltas in that accession, so each slice still counts once inside its study. "
        "The mouse-level sampling variance treats mice as independent, uses a within-arm pooled residual variance, and adds a covariance when two slices share mice (GSE130472 only). "
        f"That technical variance is not the primary weight. Pooled slice-replication variances are Tacstd2 {tac['sigma_slice']:.3f}, Cldn4 {cld['sigma_slice']:.3f}, TJ {tj['sigma_slice']:.3f} "
        "(studies with at least two slices, denominator Σ(k−1)). "
        "Primary vᵢ is the larger of the technical variance and (slice-replication variance / kᵢ). "
        "An unexpressed gene has a technical variance near zero and a delta near zero; without the floor it would dominate the meta-analysis. "
        f"That is what happens to Cldn4 under technical weights alone: {cld['technical_max_weight_study']} carries {cld['technical_max_weight_pct']:.1f}% of the weight. "
        "That fit is reported as a sensitivity and is not the estimate used below."
    )
    a("")
    a("Across studies the model is normal-normal, yᵢ ~ N(μ, vᵢ + τ²), with τ² fit by REML (Viechtbauer 2005). The reported interval is a modified Knapp–Hartung interval: t critical value on k−1 degrees of freedom, and the Hartung scale is not allowed to shrink below the Wald standard error. Raw Knapp–Hartung, DerSimonian–Laird, technical-SE weights, slice-variance-only weights, within-study GLS, and ignoring shared mice are sensitivities. Leave-one-study-out refits the primary REML model 22 times.")
    a("")
    a("## Study-level sign count")
    a("")
    a("Each study contributes the sign of its mean slice delta. This is the 49/64 question asked once per study.")
    a("")
    a("| Marker | Studies up | Down | Tie | Binomial *p* | Wilcoxon *p* | Median Δ | Mean Δ (Student 95% CI, *p*) |")
    a("|---|---:|---:|---:|---:|---:|---:|---|")
    for block in (tac, cld, tj):
        s = block["study_signs"]
        w = block["study_wilcoxon"]
        vote = block["one_vote"]
        a(
            f"| {block['marker']} | **{s['n_up']}/{s['n']}** | {s['n_down']} | {s['n_tie']} | "
            f"{fmt_p(s['binomial_p_two_sided'])} | {fmt_p(w['p_two_sided'])} | {vote['median']:+.3f} | "
            f"{vote['mean']:+.3f} ({vote['ci_low']:+.3f} to {vote['ci_high']:+.3f}, *p* = {fmt_p(vote['p_two_sided'])}) |"
        )
    a("")
    a(
        f"Tacstd2 is up in **{ts['n_up']}/{ts['n']} studies** (binomial *p* = {fmt_p(ts['binomial_p_two_sided'])}; Wilcoxon *p* = {fmt_p(tac['study_wilcoxon']['p_two_sided'])}). "
        f"Cldn4 is up in **{cs['n_up']}/{cs['n']}** (binomial *p* = {fmt_p(cs['binomial_p_two_sided'])}; Wilcoxon *p* = {fmt_p(cld['study_wilcoxon']['p_two_sided'])}). "
        f"TJ is up in **{js['n_up']}/{js['n']}** (*p* = {fmt_p(js['binomial_p_two_sided'])}). "
        f"Cldn4’s equal-weight mean ({cld['one_vote']['mean']:+.3f}, Student *p* = {fmt_p(cld['one_vote']['p_two_sided'])}) sits above its median ({cld['one_vote']['median']:+.3f}) because a few studies have large positive deltas. "
        "The rank test and the sign count do not call that a consistent increase. Tacstd2 is positive on the sign count, the rank test, and the mean."
    )
    a("")
    a("## REML mean difference")
    a("")
    a("| Marker | μ | Wald 95% CI | Wald *p* | mKH 95% CI | mKH *p* | τ² | I² | Q *p* | 95% PI | max weight |")
    a("|---|---:|---|---:|---|---:|---:|---:|---:|---|---|")
    for block in (tac, cld, tj):
        f = block["fit"]
        a(
            f"| {block['marker']} | {f['mu']:+.3f} | {f['ci_low']:+.3f} to {f['ci_high']:+.3f} | {fmt_p(f['p_wald'])} | "
            f"{f['mkh_low']:+.3f} to {f['mkh_high']:+.3f} | {fmt_p(f['p_mkh'])} | {f['tau2']:.3f} | {100 * f['i2']:.1f}% | "
            f"{fmt_p(f['Q_p'])} | {f['pi_low']:+.3f} to {f['pi_high']:+.3f} | "
            f"{block['max_weight_study']} {block['max_weight_pct']:.1f}% |"
        )
    a("")
    a(
        f"Tacstd2 REML μ = {_ci_text(tf)}. "
        f"Cldn4 REML μ = {_ci_text(cf)}. "
        f"TJ REML μ = {_ci_text(jf)}."
    )
    a("")
    a(
        f"Tacstd2’s modified Knapp–Hartung interval lies above 0 (*p* = {fmt_p(tf['p_mkh'])}). "
        f"I² is {100 * tf['i2']:.0f}%. The prediction interval ({tf['pi_low']:+.3f} to {tf['pi_high']:+.3f}) includes decreases, so a new study need not be positive. "
        f"GSE124821 carries {tac['gse124821_weight_pct']:.1f}% of the Tacstd2 random-effects weight."
    )
    a("")
    ytn = tac["studies"].loc[tac["studies"]["gse_id"] == "GSE146027"].iloc[0]
    b16 = tac["studies"].loc[tac["studies"]["gse_id"] == "GSE149825"].iloc[0]
    a(
        f"The Tacstd2 mean is pulled by large study effects. GSE146027 (YTN16) is {ytn['y']:+.3f} "
        f"(slice deltas {ytn['slice_delta_min']:+.3f} to {ytn['slice_delta_max']:+.3f}). "
        f"GSE149825 (B16) is {b16['y']:+.3f}, the average of slice deltas {b16['slice_delta_min']:+.3f} and {b16['slice_delta_max']:+.3f} "
        "(plain dual checkpoint blockade versus birinapant plus dual blockade). Both slices stay inside the locked 64."
    )
    a("")
    a(
        f"Cldn4’s pooled mean is smaller ({cf['mu']:+.3f} versus Tacstd2 {tf['mu']:+.3f}) and its Wald and modified Knapp–Hartung intervals both include 0 "
        f"(*p* = {fmt_p(cf['p_wald'])} and {fmt_p(cf['p_mkh'])}). "
        f"The locked slice test remains {cr['n_up']}/{cr['n_slices']}, *p* = {fmt_p(cr['wilcoxon_p_two_sided'])}, and the study sign count remains {cs['n_up']}/{cs['n']}. "
        "Cldn4 is weaker than Tacstd2 on the locked slice test, the 22-study sign count, the study-level Wilcoxon test, and the primary REML interval. "
        "The equal-weight Student interval for the Cldn4 mean just excludes 0; that is the summary the right skew can move, and it is not the locked claim."
    )
    a("")
    a(
        f"TJ follows Cldn4 rather than Tacstd2: μ = {jf['mu']:+.3f}, interval includes 0, study signs {js['n_up']}/{js['n']}, slice signs {jr['n_up']}/{jr['n_slices']}."
    )
    a("")
    a("Forest plots: `results/tismo_icb_study_meta/figures/forest_Tacstd2.png`, `forest_Cldn4.png`, `forest_TJ.png`, and `forest_panel.png`. Pooled comparison: `pooled_comparison.png`.")
    a("")
    a("## Leave-one-study-out")
    a("")
    a(
        f"Refitting REML after dropping each study moves Tacstd2 μ across {loso_span(tac)}, Cldn4 across {loso_span(cld)}, and TJ across {loso_span(tj)}. "
        "Full table: `results/tismo_icb_study_meta/tables/loso.tsv`. Figure: `figures/loso_panel.png`."
    )
    a("")
    for marker, block in blocks.items():
        changed = crossed(block)
        if changed:
            a(f"{marker}: leaving out {', '.join(changed)} changes whether the Wald interval lies entirely above 0.")
        else:
            a(
                f"{marker}: no single study flips whether the Wald interval lies entirely above 0 "
                f"(full interval {'does' if block['fit']['ci_low'] > 0 else 'does not'})."
            )
    a("")
    # Name the largest |delta_mu| study.
    for marker, block in blocks.items():
        top = max(block["loso"], key=lambda r: abs(r["delta_mu"]))
        a(
            f"{marker} moves most when {top['left_out']} is removed (μ {block['fit']['mu']:+.3f} → {top['mu']:+.3f})."
        )
    a("")
    a("## Sensitivities")
    a("")
    a("Same 22 studies. μ is the REML pooled mean except the one-vote row, which is the unweighted mean of study effects with a Student interval.")
    a("")
    a("| Marker | Primary | Slice variance only | Technical SE | DL on primary v | One-vote mean | Technical max weight |")
    a("|---|---:|---:|---:|---:|---:|---|")
    for block in (tac, cld, tj):
        a(
            f"| {block['marker']} | {block['fit']['mu']:+.3f} | {block['fit_slice_only']['mu']:+.3f} | "
            f"{block['fit_technical']['mu']:+.3f} | {block['fit_dl']['mu']:+.3f} | "
            f"{block['one_vote']['mean']:+.3f} | {block['technical_max_weight_study']} {block['technical_max_weight_pct']:.1f}% |"
        )
    a("")
    a(
        "DerSimonian–Laird on the primary variances agrees with REML at the reported precision. "
        "Slice-variance-only weights are close to the primary fit. "
        f"Technical-SE weights are not usable for Cldn4: {cld['technical_max_weight_study']} (B16, Cldn4 at the expression floor, study Δ {cld['studies'].loc[cld['studies']['gse_id']==cld['technical_max_weight_study'], 'y'].iloc[0]:+.4f}) "
        f"takes {cld['technical_max_weight_pct']:.1f}% of the weight and pulls μ to {cld['fit_technical']['mu']:+.4f}. "
        "That number is a measurement-floor artifact. It is not evidence that Cldn4 is unchanged in the studies where it is expressed, and it is not used as the result."
    )
    a("")
    a("### Cancer groups with at least 3 studies")
    a("")
    a("Lung is one study (LLC) and is not pooled as its own meta-analysis. There is no KL subgroup.")
    a("")
    a("| Marker | Group | Studies up | μ | Wald 95% CI | Wald *p* | I² |")
    a("|---|---|---:|---:|---|---:|---:|")
    any_sub = False
    for block in (tac, cld, tj):
        for row in block["subgroups"]:
            any_sub = True
            a(
                f"| {row['marker']} | {row['cancer_group']} | {row['n_up']}/{row['k']} | {row['mu']:+.3f} | "
                f"{row['ci_low']:+.3f} to {row['ci_high']:+.3f} | {fmt_p(row['p_wald'])} | {100 * row['i2']:.0f}% |"
            )
    if not any_sub:
        a("| — | — | — | — | — | — | — |")
    a("")
    n_clear = sum(1 for block in (tac, cld, tj) for row in block["subgroups"] if row["ci_low"] > 0)
    mel = next(row for row in tac["subgroups"] if row["cancer_group"] == "Melanoma")
    a(
        f"Subgroup Wald intervals entirely above 0: {n_clear}. "
        f"Melanoma Tacstd2 is up in {mel['n_up']}/{mel['k']} studies "
        f"(μ {mel['mu']:+.3f}, Wald {mel['ci_low']:+.3f} to {mel['ci_high']:+.3f}). "
        "The overall Tacstd2 mean is not a within-histology certainty."
    )
    a("")
    a("## Lung is LLC, not KL")
    a("")
    lung_rows = tac["lung_studies"]
    if len(lung_rows) != 1 or lung_rows[0]["gse_id"] != "GSE155972":
        raise SystemExit(f"unexpected lung studies: {lung_rows}")
    lung = lung_rows[0]
    cld_lung = cld["lung_studies"][0]
    tj_lung = tj["lung_studies"][0]
    a(
        f"TISMO ICB in this 64-slice lock contains one lung study: **{lung['gse_id']}, line LLC** "
        f"({int(lung['n_slices'])} slices, {int(lung['n_mice'])} mice). "
        f"Study mean Δ: Tacstd2 {lung['y']:+.3f}, Cldn4 {cld_lung['y']:+.3f}, TJ {tj_lung['y']:+.3f}. "
        "LLC is Lewis lung carcinoma. It is not a Kras/Stk11 (KL) or Kras/Trp53 (KP) model."
    )
    a("")
    a(
        f"Catalog check: cellLineMeta has {audit['n_cell_lines']} lines; vivoMeta has {audit['n_vivo_annotation_rows']} rows. "
        f"STK11 mentioned in cellLineMeta: {audit['cellLineMeta_mentions_stk11']}. "
        f"LKB1 mentioned: {audit['cellLineMeta_mentions_lkb1']}. "
        f"Same two strings in vivoMeta: {audit['vivoMeta_mentions_stk11']}, {audit['vivoMeta_mentions_lkb1']}. "
        f"Lung lines on the catalog: "
        + ", ".join(f"{r['cellLine']} ({r['cancerType']})" for r in audit["lung_lines_in_catalog"])
        + f". KPB25L is {audit['kpb25l_cancer']} and stays in the mammary group. KPC is {audit['kpc_cancer']} and is not in this ICB lock."
    )
    a("")
    a("No KL or KP row was added. The older label `kl_kp_llc` is not used here, because that stratum was LLC twice.")
    a("")
    a("## What this does not say")
    a("")
    a("- It does not replace Tacstd2 49/64.")
    a("- It does not say Cldn4 rises after ICB in the same way Tacstd2 does. The locked slice test (34/64, *p* = 0.077), the 22-study sign count (12/22), and the primary REML interval (includes 0) all leave Cldn4 weaker.")
    a("- It does not treat a technical standard error of an unexpressed gene as infinite information. That sensitivity is shown and set aside.")
    a("- It does not say the next study will be positive. Prediction intervals cover negative deltas.")
    a("- It does not split responders from non-responders. The lock pools them.")
    a("- It does not use private 8-KL matrices, and it does not call LLC a KL line.")
    a("")
    a("## Reproduce")
    a("")
    a("```bash")
    a("pip install -r methods/tismo_icb_study_meta/requirements.txt")
    a("python3 methods/tismo_icb_study_meta/analyze.py")
    a("python3 -m unittest tests/test_tismo_icb_study_meta.py")
    a("```")
    a("")
    a("Inputs are the vendored TISMO Gene-module CSVs in `data/tismo/vivo/` (same export as the 49/64 lock) plus `cellLineMeta.json` and `vivoMeta.json` for the KL audit.")
    a("")
    a("## Files")
    a("")
    a("| Path | Content |")
    a("|---|---|")
    a("| `results/tismo_icb_study_meta/tables/locked_reference.tsv` | 49/64, 34/64, 37/64 recomputed |")
    a("| `results/tismo_icb_study_meta/tables/slice_effects.tsv` | 64 paired deltas per marker |")
    a("| `results/tismo_icb_study_meta/tables/study_effects.tsv` | 22 study effects, variances, weights |")
    a("| `results/tismo_icb_study_meta/tables/reml_summary.tsv` | REML and sensitivity fits |")
    a("| `results/tismo_icb_study_meta/tables/loso.tsv` | Leave-one-study-out |")
    a("| `results/tismo_icb_study_meta/tables/subgroup_reml.tsv` | Cancer groups with ≥3 studies |")
    a("| `results/tismo_icb_study_meta/figures/` | Forests, LOSO, pooled comparison |")
    a("")
    return "\n".join(lines) + "\n"


def drop_id_columns(slices: pd.DataFrame) -> pd.DataFrame:
    keep = [c for c in slices.columns if c not in {"base_ids", "icb_ids", "base_values", "icb_values"}]
    return slices[keep].copy()


def run() -> dict:
    loaded = load_locked_slices()
    audit = kl_audit()
    if audit["cellLineMeta_mentions_stk11"] or audit["cellLineMeta_mentions_lkb1"]:
        raise SystemExit("unexpected STK11/LKB1 hit in cellLineMeta; do not invent a KL call without reading it")
    if audit["vivoMeta_mentions_stk11"] or audit["vivoMeta_mentions_lkb1"]:
        raise SystemExit("unexpected STK11/LKB1 hit in vivoMeta")
    blocks = {}
    for marker, slices in loaded["slices"].items():
        blocks[marker] = analyze_marker(marker, slices)
    return {"blocks": blocks, "audit": audit, "sample_meta": loaded["sample_meta"]}


def write_outputs(payload: dict) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    blocks = payload["blocks"]
    ref_rows = [blocks[m]["reference"] for m in MARKERS]
    pd.DataFrame(ref_rows).to_csv(TABLES / "locked_reference.tsv", sep="\t", index=False)

    loaded = load_locked_slices()
    slice_frames = []
    for marker, slices in loaded["slices"].items():
        part = drop_id_columns(slices)
        part.insert(0, "marker", marker)
        slice_frames.append(part)
    pd.concat(slice_frames, ignore_index=True).to_csv(TABLES / "slice_effects.tsv", sep="\t", index=False)

    study_frames = []
    fit_rows = []
    loso_rows = []
    sub_rows = []
    for marker in MARKERS:
        block = blocks[marker]
        part = block["studies"].drop(columns=["cancer_rank"]).copy()
        part.insert(0, "marker", marker)
        study_frames.append(part)
        for name, fit in (
            ("REML_primary", block["fit"]),
            ("DL_primary", block["fit_dl"]),
            ("REML_technical_SE", block["fit_technical"]),
            ("REML_slice_variance_only", block["fit_slice_only"]),
            ("REML_ignore_shared_mice", block["fit_independence"]),
            ("REML_within_study_GLS", block["fit_gls"]),
        ):
            fit_rows.append(flat_fit(marker, name, fit))
        vote = block["one_vote"]
        fit_rows.append(
            {
                "marker": marker,
                "model": "one_study_one_vote",
                "k": block["fit"]["k"],
                "mu": vote["mean"],
                "se": vote["se"],
                "ci_low": vote["ci_low"],
                "ci_high": vote["ci_high"],
                "p_wald": vote["p_two_sided"],
                "tau2": math.nan,
                "i2": math.nan,
                "h2": math.nan,
                "Q": math.nan,
                "Q_df": block["fit"]["k"] - 1,
                "Q_p": math.nan,
                "mu_fe": math.nan,
                "se_fe": math.nan,
                "fe_low": math.nan,
                "fe_high": math.nan,
                "kh_low": math.nan,
                "kh_high": math.nan,
                "p_kh": math.nan,
                "mkh_low": math.nan,
                "mkh_high": math.nan,
                "p_mkh": math.nan,
                "pi_low": math.nan,
                "pi_high": math.nan,
            }
        )
        loso_rows.extend(block["loso"])
        sub_rows.extend(block["subgroups"])
    pd.concat(study_frames, ignore_index=True).to_csv(TABLES / "study_effects.tsv", sep="\t", index=False)
    pd.DataFrame(fit_rows).to_csv(TABLES / "reml_summary.tsv", sep="\t", index=False)
    pd.DataFrame(loso_rows).to_csv(TABLES / "loso.tsv", sep="\t", index=False)
    pd.DataFrame(sub_rows).to_csv(TABLES / "subgroup_reml.tsv", sep="\t", index=False)

    for marker in MARKERS:
        write_forest(marker, blocks[marker], FIG / f"forest_{marker}.png")
    write_panel(blocks, FIG / "forest_panel.png")
    write_loso(blocks, FIG / "loso_panel.png")
    write_pooled_comparison(blocks, FIG / "pooled_comparison.png")

    text = render_results(blocks, payload["audit"], payload["sample_meta"])
    (ROOT / "RESULTS.md").write_text(text)
    (OUT / "RESULTS.md").write_text(text)

    summary = {
        "locked_reference": ref_rows,
        "input_md5": {p.name: file_md5(VIVO / p.name) for p in sorted(VIVO.glob("*.csv"))},
        "kl_audit": payload["audit"],
        "sample_meta": payload["sample_meta"],
        "reml": [
            {
                "marker": m,
                **{k: blocks[m]["fit"][k] for k in (
                    "mu", "se", "ci_low", "ci_high", "p_wald", "tau2", "i2", "Q", "Q_p",
                    "kh_low", "kh_high", "p_kh", "mkh_low", "mkh_high", "p_mkh", "pi_low", "pi_high", "k",
                )},
                "sigma_slice": blocks[m]["sigma_slice"],
                "technical_max_weight_study": blocks[m]["technical_max_weight_study"],
                "technical_max_weight_pct": blocks[m]["technical_max_weight_pct"],
                "study_n_up": blocks[m]["study_signs"]["n_up"],
                "study_n": blocks[m]["study_signs"]["n"],
                "study_binomial_p": blocks[m]["study_signs"]["binomial_p_two_sided"],
                "max_weight_study": blocks[m]["max_weight_study"],
                "max_weight_pct": blocks[m]["max_weight_pct"],
            }
            for m in MARKERS
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (TABLES / "kl_audit.json").write_text(json.dumps(payload["audit"], indent=2) + "\n")


def main() -> int:
    payload = run()
    write_outputs(payload)
    for marker in MARKERS:
        ref = payload["blocks"][marker]["reference"]
        fit = payload["blocks"][marker]["fit"]
        signs = payload["blocks"][marker]["study_signs"]
        print(
            f"{marker}: locked {ref['n_up']}/{ref['n_slices']} p={ref['wilcoxon_p_two_sided']:.3e}; "
            f"studies {signs['n_up']}/{signs['n']}; "
            f"REML μ={fit['mu']:+.3f} Wald[{fit['ci_low']:+.3f},{fit['ci_high']:+.3f}] "
            f"p={fit['p_wald']:.3e} mKH p={fit['p_mkh']:.3e} I2={100*fit['i2']:.1f}%"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Does the TISMO Tacstd2–infiltrate link attenuate after holding Cldn4 fixed?

Pre-specified scores, not a search over the immune export:

- granulocyte / neutrophil scores, including Granulocytes_mMCPcounter
  (the study-level Tacstd2 association in the prior TISMO infiltrate table)
  and Neutrophils_mMCPcounter (the mouse-level neutrophil score)
- the six CD8 scores reported with that table

Units:

- study: one GEO or ArrayExpress series, unweighted median of its baseline mice
- mouse: one SRX run
- line: one syngeneic model, median of its baseline mice
- within-line: mouse values after subtracting the line median

The partial Spearman residualizes the ranks of both variables on the rank of
the control gene. The rank-regression coefficient is the Tacstd2 slope in
immune ~ Tacstd2 + Cldn4 on z-scored ranks, which is the attenuation estimand.
Neither estimand is a causal mediation test.

ICI responder labels are an arm property. The study-level response rate is the
unit used for the partial. Arm-level odds ratios are reported beside it.
Mouse-level response tests are not used: 62 of 64 arms repeat one label.

Nothing is imputed. The locked Tacstd2 49/64 ICB-versus-control count is not
recomputed. LLC is not labeled KL.
"""

from __future__ import annotations

import json
import math
import re
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "tismo_max_effect"
OUT = ROOT / "results" / "tismo_partial_cldn4"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

N_BOOT = 2000
SEED = 0

GRANULOCYTE = (
    "Granulocytes_mMCPcounter",
    "Neutrophils_mMCPcounter",
    "Neutrophil_TIMER",
    "Neutrophils_CIBERSORT_abs",
    "Neutrophil_quanTIseq",
    "Neutrophil_xCell",
)
CD8 = (
    "T CD8_TIMER",
    "T CD8_CIBERSORT_abs",
    "T CD8_EPIC",
    "T CD8_quanTIseq",
    "T CD8_xCell",
    "CD8 T_mMCPcounter",
)
FEATURES = GRANULOCYTE + CD8
PRIMARY_GRAN = "Granulocytes_mMCPcounter"
PRIMARY_NEUT = "Neutrophils_mMCPcounter"
PRIMARY_CD8 = "T CD8_CIBERSORT_abs"

SHORT = {
    "Granulocytes_mMCPcounter": "Granulocytes mMCP",
    "Neutrophils_mMCPcounter": "Neutrophils mMCP",
    "Neutrophil_TIMER": "Neutrophil TIMER",
    "Neutrophils_CIBERSORT_abs": "Neutrophils CIBERSORT",
    "Neutrophil_quanTIseq": "Neutrophil quanTIseq",
    "Neutrophil_xCell": "Neutrophil xCell",
    "T CD8_TIMER": "CD8 TIMER",
    "T CD8_CIBERSORT_abs": "CD8 CIBERSORT",
    "T CD8_EPIC": "CD8 EPIC",
    "T CD8_quanTIseq": "CD8 quanTIseq",
    "T CD8_xCell": "CD8 xCell",
    "CD8 T_mMCPcounter": "CD8 mMCP",
}

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
    "MOC22": "HeadNeck",
}
CANCER_COLOR = {
    "Mammary": "#b45309",
    "Melanoma": "#1d4ed8",
    "Colorectal": "#15803d",
    "Lung": "#b91c1c",
    "Sarcoma": "#7c3aed",
    "Liver": "#0f766e",
    "Gastric": "#be185d",
    "HeadNeck": "#44403c",
    "mixed": "#737373",
}


def stem_label(value: str) -> str:
    return re.sub(r"\(n=\d+\)$", "", str(value)).rstrip()


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, (np.floating, float)):
        value = float(value)
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def fmt(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_ci(lo, hi, digits=3) -> str:
    return f"{fmt(lo, digits)} to {fmt(hi, digits)}"


def safe_corr(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def n_unique(a) -> int:
    a = np.asarray(a, dtype=float)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return 0
    return int(np.unique(np.round(a, 10)).size)


def estimate_arrays(x, y, z) -> dict:
    """Spearman, partial Spearman, and the adjusted rank-regression slope."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    x, y, z = x[mask], y[mask], z[mask]
    n = int(len(x))
    empty = {
        "n": n,
        "rho": float("nan"),
        "partial": float("nan"),
        "beta": float("nan"),
        "semipartial": float("nan"),
        "rho_xz": float("nan"),
        "rho_yz": float("nan"),
        "analytic_p": float("nan"),
    }
    if n < 5 or min(n_unique(x), n_unique(y), n_unique(z)) < 2:
        return empty
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    rho = safe_corr(rx, ry)
    rho_xz = safe_corr(rx, rz)
    rho_yz = safe_corr(ry, rz)
    denom_x = 1.0 - rho_xz**2
    denom_y = 1.0 - rho_yz**2
    if rho != rho or not (denom_x > 1e-10 and denom_y > 1e-10):
        empty["n"] = n
        empty["rho"] = float(rho) if rho == rho else float("nan")
        empty["rho_xz"] = float(rho_xz) if rho_xz == rho_xz else float("nan")
        empty["rho_yz"] = float(rho_yz) if rho_yz == rho_yz else float("nan")
        return empty
    numer = rho - rho_xz * rho_yz
    partial = numer / math.sqrt(denom_x * denom_y)
    beta = numer / denom_x
    # Semipartial: Cldn4-residualized exposure rank vs the raw outcome rank.
    design = np.column_stack([np.ones(n), rz])
    coef, *_ = np.linalg.lstsq(design, rx, rcond=None)
    semipartial = safe_corr(rx - design @ coef, ry)
    if abs(partial) < 1:
        tstat = partial * math.sqrt((n - 3) / (1.0 - partial**2))
        analytic_p = float(2 * stats.t.sf(abs(tstat), n - 3))
    else:
        analytic_p = float("nan")
    return {
        "n": n,
        "rho": float(rho),
        "partial": float(partial),
        "beta": float(beta),
        "semipartial": float(semipartial),
        "rho_xz": float(rho_xz),
        "rho_yz": float(rho_yz),
        "analytic_p": analytic_p,
    }


def _draw_indices(rng, n: int, clusters: np.ndarray | None) -> np.ndarray:
    if clusters is None:
        return rng.integers(0, n, n)
    codes, inverse = np.unique(clusters, return_inverse=True)
    drawn = rng.integers(0, len(codes), len(codes))
    buckets = [np.flatnonzero(inverse == i) for i in range(len(codes))]
    return np.concatenate([buckets[i] for i in drawn])


def bootstrap_unit(frame: pd.DataFrame, xcol: str, zcol: str, features: tuple[str, ...], cluster_col: str | None, seed_key: str) -> dict[str, dict]:
    """Bootstrap every feature on the same resamples. Returns feature -> draws."""
    seed = zlib.adler32(f"{SEED}:{seed_key}".encode("utf-8")) & 0xFFFFFFFF
    rng = np.random.default_rng(seed)
    x = frame[xcol].to_numpy(dtype=float)
    z = frame[zcol].to_numpy(dtype=float)
    clusters = frame[cluster_col].to_numpy() if cluster_col else None
    ys = {feat: frame[feat].to_numpy(dtype=float) for feat in features}
    n = len(frame)
    store = {
        feat: {"partial": [], "rho": [], "beta": [], "delta": []}
        for feat in features
    }
    for _ in range(N_BOOT):
        idx = _draw_indices(rng, n, clusters)
        xb, zb = x[idx], z[idx]
        if n_unique(xb) < 2 or n_unique(zb) < 2:
            continue
        for feat, y in ys.items():
            est = estimate_arrays(xb, y[idx], zb)
            if est["rho"] != est["rho"]:
                continue
            store[feat]["rho"].append(est["rho"])
            if est["partial"] != est["partial"]:
                continue
            store[feat]["partial"].append(est["partial"])
            store[feat]["beta"].append(est["beta"])
            store[feat]["delta"].append(est["rho"] - est["partial"])
    out = {}
    for feat, draws in store.items():
        rec = {"n_boot": len(draws["partial"])}
        for key in ("partial", "rho", "beta", "delta"):
            arr = np.asarray(draws[key], dtype=float)
            if arr.size < 0.8 * N_BOOT:
                rec[f"{key}_lo"] = float("nan")
                rec[f"{key}_hi"] = float("nan")
            else:
                lo, hi = np.percentile(arr, [2.5, 97.5])
                rec[f"{key}_lo"] = float(lo)
                rec[f"{key}_hi"] = float(hi)
        out[feat] = rec
    return out


def leave_one_cluster_out(frame: pd.DataFrame, xcol: str, zcol: str, features: tuple[str, ...], cluster_col: str) -> dict[str, dict]:
    out = {
        feat: {
            "loo_partial_min": float("nan"),
            "loo_partial_max": float("nan"),
            "loo_influence_id": None,
            "loo_partial_at_influence": float("nan"),
            "loo_rho_min": float("nan"),
            "loo_rho_max": float("nan"),
        }
        for feat in features
    }
    clusters = frame[cluster_col].drop_duplicates().tolist()
    if len(clusters) < 8:
        return out
    collected = {feat: [] for feat in features}
    for cid in clusters:
        sub = frame[frame[cluster_col] != cid]
        for feat in features:
            est = estimate_arrays(sub[xcol], sub[feat], sub[zcol])
            if est["partial"] == est["partial"]:
                collected[feat].append((cid, est["partial"], est["rho"]))
    for feat, rows in collected.items():
        if len(rows) < 8:
            continue
        partials = [r[1] for r in rows]
        rhos = [r[2] for r in rows]
        out[feat]["loo_partial_min"] = float(min(partials))
        out[feat]["loo_partial_max"] = float(max(partials))
        out[feat]["loo_rho_min"] = float(min(rhos))
        out[feat]["loo_rho_max"] = float(max(rhos))
    return out


def attach_influence(frame, xcol, zcol, feature, cluster_col, full_partial) -> tuple[str | None, float]:
    if full_partial != full_partial:
        return None, float("nan")
    worst_id, worst_partial, worst_gap = None, float("nan"), -1.0
    for cid in frame[cluster_col].drop_duplicates().tolist():
        sub = frame[frame[cluster_col] != cid]
        est = estimate_arrays(sub[xcol], sub[feature], sub[zcol])
        if est["partial"] != est["partial"]:
            continue
        gap = abs(est["partial"] - full_partial)
        if gap > worst_gap:
            worst_id, worst_partial, worst_gap = str(cid), float(est["partial"]), gap
    return worst_id, worst_partial


def load_gene(path: Path, gene: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    if int(df.groupby("Samples")["value"].nunique(dropna=True).max()) > 1:
        raise RuntimeError(f"{gene}: one sample has conflicting expression values")
    # Shared isotype controls are copied onto two arms. The value matches; keep one row.
    df = df.drop_duplicates("Samples").dropna(subset=["value"]).copy()
    df["stem"] = df["cell_line"].map(stem_label)
    df["line"] = df["stem"].str.split("_").str[0]
    df["cancer"] = df["line"].map(CANCER).fillna("Other")
    df = df.rename(columns={"value": gene})
    return df


def load_immune(features: tuple[str, ...]) -> pd.DataFrame:
    raw = pd.read_csv(DATA / "infiltrate_icb.csv", usecols=["Samples", "geneID", "value"])
    raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
    raw = raw[raw["geneID"].isin(features)].copy()
    nunique = raw.groupby(["Samples", "geneID"])["value"].nunique(dropna=True)
    if int((nunique > 1).sum()) > 0:
        raise RuntimeError("infiltrate table has conflicting values")
    wide = raw.pivot_table(index="Samples", columns="geneID", values="value", aggfunc="mean")
    missing = [feat for feat in features if feat not in wide.columns]
    if missing:
        raise RuntimeError(f"missing infiltrate scores: {missing}")
    return wide[list(features)]


def series_row(g: pd.DataFrame, features: tuple[str, ...]) -> dict:
    lines = sorted(g["line"].unique())
    rec = {
        "unit_id": g["GSE_ID"].iloc[0],
        "n_mice": int(len(g)),
        "n_lines": int(len(lines)),
        "line": lines[0] if len(lines) == 1 else "mixed",
        "lines": ",".join(lines),
        "cancer": g["cancer"].iloc[0] if g["cancer"].nunique() == 1 else "mixed",
        "Tacstd2": float(g["Tacstd2"].median()),
        "Cldn4": float(g["Cldn4"].median()),
    }
    for feat in features:
        rec[feat] = float(g[feat].median()) if g[feat].notna().any() else float("nan")
    return rec


def build_frames(paired: pd.DataFrame, features: tuple[str, ...]) -> dict[str, pd.DataFrame]:
    base = paired[paired["Baseline"] == 1].copy()
    study = pd.DataFrame([series_row(g, features) for _, g in base.groupby("GSE_ID", sort=True)])
    line_rows = []
    for line, g in base.groupby("line", sort=True):
        rec = {
            "unit_id": line,
            "n_mice": int(len(g)),
            "n_lines": 1,
            "line": line,
            "lines": line,
            "cancer": g["cancer"].iloc[0],
            "Tacstd2": float(g["Tacstd2"].median()),
            "Cldn4": float(g["Cldn4"].median()),
            "cluster": line,
        }
        for feat in features:
            rec[feat] = float(g[feat].median()) if g[feat].notna().any() else float("nan")
        line_rows.append(rec)
    line = pd.DataFrame(line_rows)
    within = base.copy()
    counts = within["line"].value_counts()
    within = within[within["line"].isin(counts[counts >= 3].index)].copy()
    for col in ("Tacstd2", "Cldn4", *features):
        within[col] = within[col] - within.groupby("line")[col].transform("median")
    within["unit_id"] = within["Samples"]
    within["cluster"] = within["line"]
    base["unit_id"] = base["Samples"]
    base["cluster"] = base["GSE_ID"]
    study["cluster"] = study["unit_id"]
    ge8 = study[study["n_mice"] >= 8].copy()
    drop_big = base[base["GSE_ID"] != "GSE124821"].copy()
    return {
        "study": study,
        "study_n_ge_8": ge8,
        "mouse": base,
        "mouse_without_GSE124821": drop_big,
        "line": line,
        "within_line": within,
    }


def score_direction(frame: pd.DataFrame, exposure: str, control: str, features: tuple[str, ...], cluster_for_boot: str | None, cluster_for_loo: str | None, seed_key: str) -> list[dict]:
    boots = bootstrap_unit(frame, exposure, control, features, cluster_for_boot, seed_key)
    loo = {}
    if cluster_for_loo is not None:
        loo = leave_one_cluster_out(frame, exposure, control, features, cluster_for_loo)
    rows = []
    for feat in features:
        est = estimate_arrays(frame[exposure], frame[feat], frame[control])
        boot = boots[feat]
        influence_id, influence_partial = (None, float("nan"))
        if cluster_for_loo is not None:
            influence_id, influence_partial = attach_influence(
                frame, exposure, control, feat, cluster_for_loo, est["partial"]
            )
        y = frame[feat].to_numpy(dtype=float)
        mask = ~(np.isnan(frame[exposure].to_numpy(dtype=float)) | np.isnan(y) | np.isnan(frame[control].to_numpy(dtype=float)))
        rec = {
            "exposure": exposure,
            "control": control,
            "feature": feat,
            "family": "granulocyte" if feat in GRANULOCYTE else "CD8",
            "n": est["n"],
            "n_unique_x": n_unique(frame[exposure].to_numpy(dtype=float)[mask]),
            "n_unique_y": n_unique(y[mask]),
            "n_unique_z": n_unique(frame[control].to_numpy(dtype=float)[mask]),
            "rho": est["rho"],
            "partial": est["partial"],
            "beta": est["beta"],
            "semipartial": est["semipartial"],
            "rho_xz": est["rho_xz"],
            "rho_yz": est["rho_yz"],
            "delta": est["rho"] - est["partial"] if est["rho"] == est["rho"] and est["partial"] == est["partial"] else float("nan"),
            "analytic_p": est["analytic_p"],
            "rho_lo": boot["rho_lo"],
            "rho_hi": boot["rho_hi"],
            "partial_lo": boot["partial_lo"],
            "partial_hi": boot["partial_hi"],
            "beta_lo": boot["beta_lo"],
            "beta_hi": boot["beta_hi"],
            "delta_lo": boot["delta_lo"],
            "delta_hi": boot["delta_hi"],
            "n_boot": boot["n_boot"],
        }
        rec.update(loo.get(feat, {}))
        rec["loo_influence_id"] = influence_id
        rec["loo_partial_at_influence"] = influence_partial
        rows.append(rec)
    return rows


def cliffs_delta(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2 * u / (len(x) * len(y)) - 1)


def median_split_or(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 8:
        return {"n": n, "eligible": False, "a": 0, "b": 0, "c": 0, "d": 0, "OR": float("nan"), "OR_low": float("nan"), "OR_high": float("nan"), "fisher_p": float("nan")}
    hx = x >= np.median(x)
    hy = y >= np.median(y)
    a = int((hx & hy).sum())
    b = int((hx & ~hy).sum())
    c = int((~hx & hy).sum())
    d = int((~hx & ~hy).sum())
    margins_ok = min(int(hx.sum()), int((~hx).sum()), int(hy.sum()), int((~hy).sum())) >= 3
    cells_ok = min(a, b, c, d) >= 1
    out = {"n": n, "a": a, "b": b, "c": c, "d": d, "eligible": bool(margins_ok and cells_ok)}
    if not out["eligible"]:
        out.update(OR=float("nan"), OR_low=float("nan"), OR_high=float("nan"), fisher_p=float("nan"))
        return out
    odd = (a / b) / (c / d)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    log_or = math.log(odd)
    _, fisher_p = stats.fisher_exact([[a, b], [c, d]])
    out.update(
        OR=float(odd),
        OR_low=float(math.exp(log_or - 1.96 * se)),
        OR_high=float(math.exp(log_or + 1.96 * se)),
        fisher_p=float(fisher_p),
    )
    return out


def mantel_haenszel(tables: list[tuple[int, int, int, int]]) -> dict:
    """OR of response given high exposure, stratified by the control-gene split.

    Each table is (a, b, c, d) = high-R, high-NR, low-R, low-NR.
    Variance is the Robins-Breslow-Greenland variance. A stratum with a zero
    margin is dropped rather than given a continuity correction.
    """
    used = []
    for a, b, c, d in tables:
        if min(a + b, c + d, a + c, b + d) < 1:
            continue
        used.append((a, b, c, d))
    if len(used) < 2:
        return {"OR": float("nan"), "OR_low": float("nan"), "OR_high": float("nan"), "n_strata": len(used), "tables": tables}
    num = 0.0
    den = 0.0
    var = 0.0
    for a, b, c, d in used:
        n = a + b + c + d
        num += a * d / n
        den += b * c / n
    if num <= 0 or den <= 0:
        return {"OR": float("nan"), "OR_low": float("nan"), "OR_high": float("nan"), "n_strata": len(used), "tables": tables}
    # Robins, Breslow, Greenland 1986.
    term_p = term_q = term_r = 0.0
    for a, b, c, d in used:
        n = a + b + c + d
        term_p += (a + d) * a * d / (n * n)
        term_q += ((a + d) * b * c + (b + c) * a * d) / (n * n)
        term_r += (b + c) * b * c / (n * n)
    var = term_p / (2 * num * num) + term_q / (2 * num * den) + term_r / (2 * den * den)
    log_or = math.log(num / den)
    se = math.sqrt(var) if var > 0 else float("nan")
    return {
        "OR": float(num / den),
        "OR_low": float(math.exp(log_or - 1.96 * se)) if se == se else float("nan"),
        "OR_high": float(math.exp(log_or + 1.96 * se)) if se == se else float("nan"),
        "n_strata": len(used),
        "tables": tables,
    }


def response_tables(paired: pd.DataFrame) -> dict:
    treated = paired[paired["Responder"].isin(["Responders", "Non-responders"])].copy()
    arm_rows = []
    for stem, g in treated.groupby("stem", sort=True):
        n_r = int((g["Responder"] == "Responders").sum())
        n_nr = int((g["Responder"] == "Non-responders").sum())
        if n_r and n_nr:
            label = "mixed"
        elif n_r:
            label = "R"
        else:
            label = "NR"
        arm_rows.append(
            {
                "stem": stem,
                "label": label,
                "response": 1.0 if label == "R" else (0.0 if label == "NR" else float("nan")),
                "Tacstd2": float(g["Tacstd2"].median()),
                "Cldn4": float(g["Cldn4"].median()),
                "n_mice": int(len(g)),
                "n_R": n_r,
                "n_NR": n_nr,
                "line": g["line"].iloc[0],
                "cancer": g["cancer"].iloc[0],
                "GSE_ID": g["GSE_ID"].iloc[0],
                "Mouse_treatment": g["Mouse_treatment"].iloc[0],
            }
        )
    arms = pd.DataFrame(arm_rows)
    pure = arms[arms["label"].isin(["R", "NR"])].copy()
    study_rows = []
    for gse, g in treated.groupby("GSE_ID", sort=True):
        study_rows.append(
            {
                "unit_id": gse,
                "Tacstd2": float(g["Tacstd2"].median()),
                "Cldn4": float(g["Cldn4"].median()),
                "response_rate": float((g["Responder"] == "Responders").mean()),
                "n_mice": int(len(g)),
                "n_lines": int(g["line"].nunique()),
                "n_arms": int(g["stem"].nunique()),
                "cancer": g["cancer"].iloc[0] if g["cancer"].nunique() == 1 else "mixed",
            }
        )
    studies = pd.DataFrame(study_rows)
    mixed = arms[arms["label"] == "mixed"].copy()
    mixed_detail = []
    for _, arm in mixed.iterrows():
        g = treated[treated["stem"] == arm["stem"]]
        for gene in ("Tacstd2", "Cldn4"):
            r = g.loc[g["Responder"] == "Responders", gene]
            nr = g.loc[g["Responder"] == "Non-responders", gene]
            mixed_detail.append(
                {
                    "gene": gene,
                    "stem": arm["stem"],
                    "GSE_ID": arm["GSE_ID"],
                    "n_R": int(len(r)),
                    "n_NR": int(len(nr)),
                    "median_R": float(r.median()),
                    "median_NR": float(nr.median()),
                    "cliff": cliffs_delta(r.to_numpy(), nr.to_numpy()),
                }
            )
    return {"arms": arms, "pure": pure, "studies": studies, "treated": treated, "mixed_detail": mixed_detail}


def score_response(pure: pd.DataFrame, studies: pd.DataFrame) -> list[dict]:
    rows = []
    # Study unit is the pre-specified ICI unit. Resample series.
    for exposure, control in (("Tacstd2", "Cldn4"), ("Cldn4", "Tacstd2")):
        est = estimate_arrays(studies[exposure], studies["response_rate"], studies[control])
        boot = bootstrap_unit(studies, exposure, control, ("response_rate",), None, f"ici-study-{exposure}")
        b = boot["response_rate"]
        loo_vals = []
        worst_id, worst_partial, worst_gap = None, float("nan"), -1.0
        for cid in studies["unit_id"]:
            sub = studies[studies["unit_id"] != cid]
            one = estimate_arrays(sub[exposure], sub["response_rate"], sub[control])
            if one["partial"] == one["partial"]:
                loo_vals.append(one["partial"])
                gap = abs(one["partial"] - est["partial"]) if est["partial"] == est["partial"] else float("nan")
                if gap == gap and gap > worst_gap:
                    worst_id, worst_partial, worst_gap = cid, one["partial"], gap
        rows.append(
            {
                "unit": "study",
                "exposure": exposure,
                "control": control,
                "n": est["n"],
                "rho": est["rho"],
                "partial": est["partial"],
                "beta": est["beta"],
                "delta": est["rho"] - est["partial"] if est["rho"] == est["rho"] and est["partial"] == est["partial"] else float("nan"),
                "rho_lo": b["rho_lo"],
                "rho_hi": b["rho_hi"],
                "partial_lo": b["partial_lo"],
                "partial_hi": b["partial_hi"],
                "delta_lo": b["delta_lo"],
                "delta_hi": b["delta_hi"],
                "loo_partial_min": float(min(loo_vals)) if loo_vals else float("nan"),
                "loo_partial_max": float(max(loo_vals)) if loo_vals else float("nan"),
                "loo_influence_id": worst_id,
                "loo_partial_at_influence": worst_partial,
                "OR": float("nan"),
                "OR_low": float("nan"),
                "OR_high": float("nan"),
                "OR_table": "",
                "MH_OR": float("nan"),
                "MH_low": float("nan"),
                "MH_high": float("nan"),
                "MH_tables": "",
                "n_from_GSE124821": int((studies["unit_id"] == "GSE124821").sum()),
            }
        )
    # Arm unit. Cluster bootstrap by series, because GSE124821 contributes many arms.
    pure = pure.copy()
    pure["response"] = (pure["label"] == "R").astype(float)
    for exposure, control in (("Tacstd2", "Cldn4"), ("Cldn4", "Tacstd2")):
        est = estimate_arrays(pure[exposure], pure["response"], pure[control])
        boot = bootstrap_unit(pure, exposure, control, ("response",), "GSE_ID", f"ici-arm-{exposure}")
        b = boot["response"]
        odd = median_split_or(pure[exposure], pure["response"])
        mx = float(pure[exposure].median())
        mz = float(pure[control].median())
        high_x = pure[exposure] >= mx
        high_z = pure[control] >= mz
        tables = []
        for flag in (True, False):
            sub = pure[high_z == flag]
            a = int(((sub[exposure] >= mx) & (sub["label"] == "R")).sum())
            b_ = int(((sub[exposure] >= mx) & (sub["label"] == "NR")).sum())
            c = int(((sub[exposure] < mx) & (sub["label"] == "R")).sum())
            d = int(((sub[exposure] < mx) & (sub["label"] == "NR")).sum())
            tables.append((a, b_, c, d))
        mh = mantel_haenszel(tables)
        rows.append(
            {
                "unit": "arm",
                "exposure": exposure,
                "control": control,
                "n": est["n"],
                "rho": est["rho"],
                "partial": est["partial"],
                "beta": est["beta"],
                "delta": est["rho"] - est["partial"] if est["rho"] == est["rho"] and est["partial"] == est["partial"] else float("nan"),
                "rho_lo": b["rho_lo"],
                "rho_hi": b["rho_hi"],
                "partial_lo": b["partial_lo"],
                "partial_hi": b["partial_hi"],
                "delta_lo": b["delta_lo"],
                "delta_hi": b["delta_hi"],
                "loo_partial_min": float("nan"),
                "loo_partial_max": float("nan"),
                "loo_influence_id": None,
                "loo_partial_at_influence": float("nan"),
                "OR": odd["OR"],
                "OR_low": odd["OR_low"],
                "OR_high": odd["OR_high"],
                "OR_table": f"{odd['a']}/{odd['b']}/{odd['c']}/{odd['d']}",
                "fisher_p": odd["fisher_p"],
                "MH_OR": mh["OR"],
                "MH_low": mh["OR_low"],
                "MH_high": mh["OR_high"],
                "MH_tables": ";".join(f"{a}/{b_}/{c}/{d}" for a, b_, c, d in tables),
                "n_from_GSE124821": int((pure["GSE_ID"] == "GSE124821").sum()),
                "median_exposure": mx,
                "median_control": mz,
                "high_z_n": int(high_z.sum()),
                "low_z_n": int((~high_z).sum()),
            }
        )
    return rows


def ci_excludes_zero(lo, hi) -> bool:
    if lo != lo or hi != hi:
        return False
    return hi < 0 or lo > 0


def shrinkage_sentence(rho, partial, delta_lo, delta_hi, partial_lo=float("nan"), partial_hi=float("nan")) -> str:
    if rho != rho or partial != partial:
        return "The partial correlation is undefined."
    closer = abs(partial) < abs(rho) - 0.02
    same_side = rho * partial > 0
    clear = ci_excludes_zero(delta_lo, delta_hi) and delta_lo > 0
    remainder = ci_excludes_zero(partial_lo, partial_hi)
    if abs(rho) < 0.15:
        return "The marginal correlation is already near zero, so there is no link to attenuate."
    if closer and same_side and clear and not remainder:
        return "The shrinkage interval stays above zero. The interval for the partial correlation includes zero."
    if closer and same_side and clear and remainder:
        return "The shrinkage interval stays above zero, and the partial correlation stays on the same side of zero."
    if closer and same_side and not clear:
        return "The point estimate is closer to zero. The interval for the shrinkage includes zero."
    if (not same_side) and abs(partial) >= 0.05:
        return "The partial correlation changes sign relative to the marginal correlation."
    if not closer:
        return "The partial correlation does not shrink toward zero."
    return "The partial correlation moves toward zero."


def lookup(df: pd.DataFrame, unit: str, exposure: str, feature: str) -> pd.Series:
    hit = df[(df["unit"] == unit) & (df["exposure"] == exposure) & (df["feature"] == feature)]
    if len(hit) != 1:
        raise RuntimeError(f"expected one row for {unit} {exposure} {feature}, found {len(hit)}")
    return hit.iloc[0]


def table_block(df: pd.DataFrame, unit: str, family: str) -> list[str]:
    sub = df[(df["unit"] == unit) & (df["family"] == family) & (df["exposure"] == "Tacstd2")].copy()
    rev = df[(df["unit"] == unit) & (df["family"] == family) & (df["exposure"] == "Cldn4")]
    lines = [
        "| Score | n | Tacstd2 ρ | partial given Cldn4 (CI) | shrinkage Δ (CI) | Cldn4 ρ | Cldn4 partial given Tacstd2 |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, row in sub.iterrows():
        other = rev[rev["feature"] == row["feature"]]
        other_rho = float(other.iloc[0]["rho"]) if len(other) else float("nan")
        other_part = float(other.iloc[0]["partial"]) if len(other) else float("nan")
        lines.append(
            f"| {SHORT[row['feature']]} | {int(row['n'])} | {fmt(row['rho'])} | "
            f"{fmt(row['partial'])} ({fmt_ci(row['partial_lo'], row['partial_hi'])}) | "
            f"{fmt(row['delta'])} ({fmt_ci(row['delta_lo'], row['delta_hi'])}) | "
            f"{fmt(other_rho)} | {fmt(other_part)} |"
        )
    return lines


def plot_added_variable(study: pd.DataFrame, path: Path, rho_xz: float, rho: float, partial: float, partial_lo: float, partial_hi: float) -> None:
    from matplotlib.lines import Line2D

    x = study["Tacstd2"].to_numpy(dtype=float)
    y = study[PRIMARY_GRAN].to_numpy(dtype=float)
    z = study["Cldn4"].to_numpy(dtype=float)
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    design = np.column_stack([np.ones(len(rz)), rz])
    bx, *_ = np.linalg.lstsq(design, rx, rcond=None)
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    xr, yr = rx - design @ bx, ry - design @ by
    colors = study["cancer"].map(CANCER_COLOR).fillna("#737373")
    present = [name for name in CANCER_COLOR if name in set(study["cancer"])]
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CANCER_COLOR[name], markeredgecolor="none", markersize=6, label=name)
        for name in present
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8))
    ax = axes[0]
    ax.scatter(z, x, s=28 + 3 * study["n_mice"], c=colors, zorder=3)
    ax.set_xlabel("Cldn4, study median")
    ax.set_ylabel("Tacstd2, study median")
    ax.set_title(f"Tacstd2 and Cldn4\nSpearman ρ = {rho_xz:.3f}, n = {len(study)} studies")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="upper left")

    ax = axes[1]
    ax.axhline(0, color="#d4d4d4", lw=0.8)
    ax.axvline(0, color="#d4d4d4", lw=0.8)
    ax.scatter(xr, yr, s=28 + 3 * study["n_mice"], c=colors, zorder=3)
    slope = np.polyfit(xr, yr, 1)
    grid = np.linspace(xr.min(), xr.max(), 40)
    ax.plot(grid, slope[0] * grid + slope[1], color="#111111", lw=1)
    ax.set_xlabel("Tacstd2 rank residual given Cldn4")
    ax.set_ylabel("Granulocyte mMCP rank residual given Cldn4")
    ax.set_title(
        f"Added-variable, granulocytes\nmarginal ρ = {rho:.3f}, partial ρ = {partial:.3f} "
        f"(CI {partial_lo:.2f} to {partial_hi:.2f})"
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_slope(df: pd.DataFrame, unit: str, path: Path, xlabel: str) -> None:
    sub = df[(df["unit"] == unit) & (df["exposure"] == "Tacstd2")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.2), sharex=True)
    for ax, family, title in (
        (axes[0], "granulocyte", "Granulocyte and neutrophil scores"),
        (axes[1], "CD8", "CD8 scores"),
    ):
        rows = sub[sub["family"] == family].iloc[::-1]
        ypos = np.arange(len(rows))
        for i, (_, row) in enumerate(rows.iterrows()):
            ax.plot([row["rho"], row["partial"]], [i, i], color="#94a3b8", lw=1.4, zorder=2)
            ax.plot(row["rho"], i, "o", color="#1d4ed8", ms=6, zorder=3)
            ax.plot(row["partial"], i, "o", color="#b91c1c", ms=6, zorder=3)
            if row["partial_lo"] == row["partial_lo"]:
                ax.plot([row["partial_lo"], row["partial_hi"]], [i, i], color="#b91c1c", lw=2.4, alpha=0.35, zorder=1)
        ax.axvline(0, color="#111111", lw=0.7)
        ax.set_yticks(ypos)
        ax.set_yticklabels([SHORT[f] for f in rows["feature"]], fontsize=8)
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(-0.55, 1.05)
    axes[0].set_xlabel(xlabel)
    axes[1].set_xlabel(xlabel)
    fig.suptitle("Blue: marginal Spearman. Red: Tacstd2 partial Spearman given Cldn4, with the resample interval.", fontsize=9, y=0.02)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_response(studies: pd.DataFrame, path: Path, rho: float, partial: float) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    sc = ax.scatter(
        studies["Tacstd2"],
        studies["response_rate"],
        s=28 + 3 * studies["n_mice"],
        c=studies["Cldn4"],
        cmap="viridis",
        zorder=3,
    )
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Cldn4, study median of treated mice")
    ax.set_xlabel("Tacstd2, study median of treated mice")
    ax.set_ylabel("Responder fraction")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title(f"ICI response by study\nmarginal ρ = {rho:.3f}, partial given Cldn4 = {partial:.3f}, n = {len(studies)}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_report(summary: dict, path: Path) -> None:
    s = summary
    g = s["primary_granulocyte_study"]
    n = s["primary_neutrophil_mouse"]
    c = s["primary_cd8_study"]
    ici = s["ici_study_tacstd2"]
    gene = s["gene_correlation"]
    neut_s = s["neutrophil_study"]
    gran_m = s["granulocyte_mouse"]
    within_neut = s["within_neutrophil"]
    within_gran = s["within_granulocyte"]
    ge8 = s["granulocyte_study_ge8"]
    cd8_mouse = s["primary_cd8_mouse"]
    rev_g = s["reverse_granulocyte_study"]
    rev_n = s["reverse_neutrophil_mouse"]

    def row_bits(row, label):
        return (
            f"{label}: marginal ρ = {row['rho']:.3f} "
            f"(resample CI {fmt_ci(row['rho_lo'], row['rho_hi'])}), "
            f"partial ρ = {row['partial']:.3f} "
            f"(CI {fmt_ci(row['partial_lo'], row['partial_hi'])}), "
            f"shrinkage Δ = {row['delta']:.3f} "
            f"(CI {fmt_ci(row['delta_lo'], row['delta_hi'])}). "
            f"{shrinkage_sentence(row['rho'], row['partial'], row['delta_lo'], row['delta_hi'], row['partial_lo'], row['partial_hi'])}"
        )

    lines = []
    lines.append("# TISMO: Tacstd2 versus infiltrate and ICI, partialling Cldn4")
    lines.append("")
    lines.append("## Paper sentence")
    lines.append("")
    lines.append(
        f"Tacstd2 and Cldn4 are one bulk axis in the {gene['study_n']} paired TISMO series "
        f"(study Spearman ρ = {gene['study_rho']:.3f}; mouse ρ = {gene['mouse_rho']:.3f}; {gene['mouse_n']} baseline mice). "
        f"On that axis, the study-level Tacstd2–granulocyte link shrinks and the largest Tacstd2–CD8 link does not. "
        f"ICI response has nothing to attenuate."
    )
    lines.append("")
    lines.append(
        f"Study unit, mMCPcounter granulocytes: Tacstd2 ρ = {g['rho']:.3f} falls to a partial ρ of {g['partial']:.3f} "
        f"given Cldn4 (resample CI {fmt_ci(g['partial_lo'], g['partial_hi'])}). "
        f"Shrinkage Δ = {g['delta']:.3f} (CI {fmt_ci(g['delta_lo'], g['delta_hi'])}). "
        f"{shrinkage_sentence(g['rho'], g['partial'], g['delta_lo'], g['delta_hi'], g['partial_lo'], g['partial_hi'])} "
        f"Cldn4 partialled on Tacstd2 moves from {rev_g['rho']:.3f} to {rev_g['partial']:.3f} "
        f"(CI {fmt_ci(rev_g['partial_lo'], rev_g['partial_hi'])}). "
        f"The two leftover partials are {g['partial']:.3f} and {rev_g['partial']:.3f}. "
        f"{'Tacstd2’s interval includes zero. ' if g['partial_lo'] < 0 < g['partial_hi'] else 'Tacstd2’s interval stays off zero. '}"
        f"{'Cldn4’s interval includes zero.' if rev_g['partial_lo'] < 0 < rev_g['partial_hi'] else 'Cldn4’s interval stays off zero.'} "
        f"That is not a clean assignment of the granulocyte link to either gene."
    )
    lines.append("")
    lines.append(
        f"Study unit, CIBERSORT CD8, the largest pre-specified CD8 score: Tacstd2 ρ = {c['rho']:.3f}, "
        f"partial given Cldn4 ρ = {c['partial']:.3f} (CI {fmt_ci(c['partial_lo'], c['partial_hi'])}), "
        f"shrinkage Δ = {c['delta']:.3f} (CI {fmt_ci(c['delta_lo'], c['delta_hi'])}). "
        f"{shrinkage_sentence(c['rho'], c['partial'], c['delta_lo'], c['delta_hi'], c['partial_lo'], c['partial_hi'])} "
        f"The partial stays positive. This is not CD8 exclusion."
    )
    lines.append("")
    lines.append(
        f"Mouse unit, the two myeloid scores diverge. Tacstd2 versus Neutrophils_mMCPcounter falls from {n['rho']:.3f} to {n['partial']:.3f} "
        f"(study-cluster CI {fmt_ci(n['partial_lo'], n['partial_hi'])}; shrinkage CI {fmt_ci(n['delta_lo'], n['delta_hi'])}). "
        f"Cldn4 versus that neutrophil score, partialled on Tacstd2, stays at {rev_n['partial']:.3f} "
        f"(CI {fmt_ci(rev_n['partial_lo'], rev_n['partial_hi'])}). "
        f"The same mice versus Granulocytes_mMCPcounter do the opposite: Tacstd2 moves from {gran_m['rho']:.3f} to a partial of {gran_m['partial']:.3f} "
        f"(CI {fmt_ci(gran_m['partial_lo'], gran_m['partial_hi'])}). "
        f"Within line, the neutrophil residual correlation is already {within_neut['rho']:.3f}. "
        f"ICI, study unit: Tacstd2 versus responder fraction ρ = {ici['rho']:.3f}, partial given Cldn4 ρ = {ici['partial']:.3f} "
        f"(CI {fmt_ci(ici['partial_lo'], ici['partial_hi'])}, n = {int(ici['n'])} series)."
    )
    lines.append("")
    lines.append(
        f"小鼠和 series 分开报，不把两个分数说成同一个答案。{gene['study_n']} 个 series 上 Tacstd2 与 Cldn4 的 ρ = {gene['study_rho']:.3f}。"
        f"series 单位 mMCP 粒细胞：Tacstd2 从 {g['rho']:.3f} 到偏相关 {g['partial']:.3f}"
        f"（区间 {fmt_ci(g['partial_lo'], g['partial_hi'])}），缩小量 {g['delta']:.3f}"
        f"（区间 {fmt_ci(g['delta_lo'], g['delta_hi'])}）。Cldn4 偏掉 Tacstd2 之后是 {rev_g['partial']:.3f}，"
        f"两个剩余相关点估计差不多，不能把粒细胞链单独记到 Cldn4 上。"
        f"CD8 最大的 CIBERSORT 从 {c['rho']:.3f} 到 {c['partial']:.3f}，缩小量的区间含 0，偏相关仍为正，不是排斥。"
        f"小鼠单位 mMCP 中性粒从 {n['rho']:.3f} 到 {n['partial']:.3f}，但 mMCP 粒细胞的小鼠偏相关仍是 {gran_m['partial']:.3f}；"
        f"扣掉细胞系中位数后中性粒边际相关是 {within_neut['rho']:.3f}。"
        f"ICI 应答在 series 单位上是 {ici['rho']:.3f} → {ici['partial']:.3f}。没有重算 49/64。LLC 不是 KL。"
    )
    lines.append("")
    lines.append("## Units")
    lines.append("")
    lines.append(
        f"A mouse is one SRX run. A study is one source series. {s['n_studies_single_line']} of {s['n_studies']} baseline series contain a single cell line. "
        f"GSE124821 pools four mammary lines ({s['gse124821_baseline_mice']} baseline mice) and GSE159344 pools two melanoma lines "
        f"({s['gse159344_baseline_mice']} mice). The study correlation gives every series one vote. The mouse correlation gives GSE124821 "
        f"{s['gse124821_baseline_mice']} of {gene['mouse_n']} votes. Within-line residuals subtract the line median before the partial correlation, "
        f"so that row is the within-model association."
    )
    lines.append("")
    lines.append(
        f"The paired table keeps mice measured for both genes. Tacstd2 baseline mice are a subset of the Cldn4 export: "
        f"{s['dropped_baseline_mice']} Cldn4 baseline mice have no Tacstd2 value "
        f"({s['dropped_baseline_detail']}). RU31562203 (MOC22) is the series that drops out. "
        f"Cldn4 versus Granulocytes_mMCPcounter on the unpaired Cldn4 export is ρ = {s['unpaired_cldn4_granulocyte_rho']:.3f} "
        f"(n = {s['unpaired_cldn4_granulocyte_n']} series). On the paired series it is ρ = {rev_g['rho']:.3f}."
    )
    lines.append("")
    lines.append(
        f"Responder labels are constant inside {s['n_pure_arms']} of {s['n_pure_arms'] + s['n_mixed_arms']} treated arms. "
        "A mouse-level response test would repeat the arm label. "
        "It is not used. The mixed arms are listed below and are too small for a partial correlation."
    )
    lines.append("")
    lines.append("TISMO has no KL or KP lung line. LLC is the only lung carcinoma in this table, and it is not called KL. The locked Tacstd2 49/64 ICB-versus-control count is not recomputed.")
    lines.append("")
    lines.append("## Collinearity")
    lines.append("")
    lines.append(
        f"Study-level Spearman of the two genes is {gene['study_rho']:.3f} (n = {gene['study_n']}; "
        f"resample CI {fmt_ci(gene['study_rho_lo'], gene['study_rho_hi'])}; "
        f"leave-one-study-out {fmt(gene['study_loo_min'])} to {fmt(gene['study_loo_max'])}). "
        f"Mouse-level Spearman is {gene['mouse_rho']:.3f} (n = {gene['mouse_n']}; "
        f"study-cluster CI {fmt_ci(gene['mouse_rho_lo'], gene['mouse_rho_hi'])}). "
        f"With that shared axis, a partial correlation is a noisy split of one program. "
        f"The variance inflation from ρ = {gene['study_rho']:.3f} is {s['vif_study']:.2f}."
    )
    lines.append("")
    lines.append("## Granulocytes")
    lines.append("")
    lines.append(f"Study unit, Tacstd2 partialled on Cldn4. {row_bits(g, 'Granulocytes mMCP')}")
    lines.append("")
    lines.append(
        f"Cldn4 partialled on Tacstd2 is ρ = {rev_g['partial']:.3f} "
        f"(CI {fmt_ci(rev_g['partial_lo'], rev_g['partial_hi'])}), from a marginal ρ of {rev_g['rho']:.3f}. "
        f"Leaving out {g.get('influence_label') or g['loo_influence_id']} moves the Tacstd2 partial from {g['partial']:.3f} to {g['loo_partial_at_influence']:.3f}. "
        f"Leave-one-study-out partials span {fmt(g['loo_partial_min'])} to {fmt(g['loo_partial_max'])}, so no single series creates the positive point estimate. "
        f"The bootstrap interval is wider than that leave-one-out range and includes zero. "
        f"On z-scored ranks the Tacstd2 coefficient in granulocytes ~ Tacstd2 + Cldn4 is {g['beta']:.3f} "
        f"(CI {fmt_ci(g['beta_lo'], g['beta_hi'])}), down from the simple coefficient {g['rho']:.3f}. "
        f"The semipartial correlation, which keeps the granulocyte rank intact, is {g['semipartial']:.3f}."
    )
    lines.append("")
    lines.append(
        f"Studies with at least 8 baseline mice (n = {int(ge8['n'])}): marginal ρ = {ge8['rho']:.3f}, "
        f"partial ρ = {ge8['partial']:.3f} (CI {fmt_ci(ge8['partial_lo'], ge8['partial_hi'])}), "
        f"Δ = {ge8['delta']:.3f} (CI {fmt_ci(ge8['delta_lo'], ge8['delta_hi'])})."
    )
    lines.append("")
    lines.append(
        f"Neutrophil scores are the rest of the pre-specified granulocyte family. They are not a second discovery set. "
        f"At the study unit, Neutrophils_mMCPcounter moves from {neut_s['rho']:.3f} to {neut_s['partial']:.3f} "
        f"(CI {fmt_ci(neut_s['partial_lo'], neut_s['partial_hi'])}; "
        f"shrinkage CI {fmt_ci(neut_s['delta_lo'], neut_s['delta_hi'])}). "
        f"{shrinkage_sentence(neut_s['rho'], neut_s['partial'], neut_s['delta_lo'], neut_s['delta_hi'], neut_s['partial_lo'], neut_s['partial_hi'])}"
    )
    lines.append("")
    lines.extend(table_block(s["specs_df_note"], "study", "granulocyte"))
    lines.append("")
    lines.append(
        f"Mouse unit, neutrophils. {row_bits(n, 'Neutrophils mMCP')} "
        f"Cldn4 versus the same score, partialled on Tacstd2, stays at {rev_n['partial']:.3f} "
        f"(marginal {rev_n['rho']:.3f}; cluster CI {fmt_ci(rev_n['partial_lo'], rev_n['partial_hi'])}). "
        f"Dropping GSE124821 leaves the Tacstd2 neutrophil partial at {s['neutrophil_mouse_no_big']['partial']:.3f} "
        f"(cluster CI {fmt_ci(s['neutrophil_mouse_no_big']['partial_lo'], s['neutrophil_mouse_no_big']['partial_hi'])}; "
        f"marginal {s['neutrophil_mouse_no_big']['rho']:.3f}, n = {int(s['neutrophil_mouse_no_big']['n'])}). "
        f"Within line, the residual Tacstd2–neutrophil correlation is {within_neut['rho']:.3f} and the partial is {within_neut['partial']:.3f} "
        f"(n = {int(within_neut['n'])} mice, {s['n_lines_within']} lines; line-cluster CI for the partial "
        f"{fmt_ci(within_neut['partial_lo'], within_neut['partial_hi'])}). "
        f"The mouse-level attenuation is a between-line fact. There is no within-line Tacstd2–neutrophil link to explain."
    )
    lines.append("")
    lines.append(
        f"Mouse-level granulocytes, the same score as the study result: marginal ρ = {gran_m['rho']:.3f}, "
        f"partial ρ = {gran_m['partial']:.3f} (cluster CI {fmt_ci(gran_m['partial_lo'], gran_m['partial_hi'])}), "
        f"shrinkage Δ = {gran_m['delta']:.3f} (CI {fmt_ci(gran_m['delta_lo'], gran_m['delta_hi'])}). "
        f"{shrinkage_sentence(gran_m['rho'], gran_m['partial'], gran_m['delta_lo'], gran_m['delta_hi'], gran_m['partial_lo'], gran_m['partial_hi'])} "
        f"Cldn4 versus mouse granulocytes, partialled on Tacstd2, is ρ = {s['reverse_granulocyte_mouse']['partial']:.3f} "
        f"from a marginal ρ of {s['reverse_granulocyte_mouse']['rho']:.3f} "
        f"(cluster CI {fmt_ci(s['reverse_granulocyte_mouse']['partial_lo'], s['reverse_granulocyte_mouse']['partial_hi'])}). "
        f"Within line, the residual granulocyte correlation is {within_gran['rho']:.3f} and the partial is {within_gran['partial']:.3f} "
        f"(line-cluster CI {fmt_ci(within_gran['partial_lo'], within_gran['partial_hi'])}; "
        f"shrinkage CI {fmt_ci(within_gran['delta_lo'], within_gran['delta_hi'])}). "
        f"mMCPcounter granulocytes and mMCPcounter neutrophils are not interchangeable in this table."
    )
    lines.append("")
    lines.append("Mouse unit, Tacstd2 partialled on Cldn4. Intervals are study-cluster bootstrap intervals.")
    lines.append("")
    lines.extend(table_block(s["specs_df_note"], "mouse", "granulocyte"))
    lines.append("")
    lines.append("## CD8")
    lines.append("")
    lines.append(
        f"CD8 was pre-specified as the six scores in the earlier TISMO table. "
        f"{row_bits(c, 'Study-level CIBERSORT CD8')} "
        f"Cldn4 versus that score moves from {s['reverse_cd8_study']['rho']:.3f} to {s['reverse_cd8_study']['partial']:.3f} "
        f"after Tacstd2 (CI {fmt_ci(s['reverse_cd8_study']['partial_lo'], s['reverse_cd8_study']['partial_hi'])}; "
        f"shrinkage CI {fmt_ci(s['reverse_cd8_study']['delta_lo'], s['reverse_cd8_study']['delta_hi'])}). "
        f"On CIBERSORT, Cldn4’s CD8 link shrinks and Tacstd2’s shrinkage interval includes zero. "
        f"quanTIseq is the other way around: the Tacstd2 shrinkage interval stays above zero and the partial interval includes zero. "
        f"TIMER, EPIC, xCell, and mMCPcounter CD8 have Tacstd2 partial intervals that include zero. "
        f"No pre-specified study-level CD8 score changes from a positive Tacstd2 correlation to a negative partial correlation. "
        f"This is not CD8 exclusion."
    )
    lines.append("")
    lines.extend(table_block(s["specs_df_note"], "study", "CD8"))
    lines.append("")
    lines.append(
        f"Mouse-level CIBERSORT CD8 is ρ = {cd8_mouse['rho']:.3f}, partial ρ = {cd8_mouse['partial']:.3f} "
        f"(cluster CI {fmt_ci(cd8_mouse['partial_lo'], cd8_mouse['partial_hi'])}). "
        f"The within-line residual marginal correlation for that score is {s['within_cd8']['rho']:.3f}. "
        f"Mouse xCell CD8 is a null marginal correlation (ρ = {s['mouse_cd8_xcell']['rho']:.3f}) whose partial is {s['mouse_cd8_xcell']['partial']:.3f} "
        f"(CI {fmt_ci(s['mouse_cd8_xcell']['partial_lo'], s['mouse_cd8_xcell']['partial_hi'])}). "
        f"That is not a marginal CD8 link that failed to attenuate. "
        f"Mouse CD8 correlations are small once the line is held fixed."
    )
    lines.append("")
    lines.extend(table_block(s["specs_df_note"], "mouse", "CD8"))
    lines.append("")
    lines.append("## Within-line residuals")
    lines.append("")
    lines.append(
        f"After subtracting each line's median, the largest |marginal ρ| in the pre-specified list is "
        f"{s['within_max']['abs_rho']:.3f} ({s['within_max']['exposure']} vs {s['within_max']['feature']}, "
        f"ρ = {s['within_max']['rho']:.3f}, partial {s['within_max']['partial']:.3f}, "
        f"n = {int(s['within_max']['n'])} mice across {s['n_lines_within']} lines). "
        f"Line medians themselves (n = {s['n_lines']} lines) are a smaller between-model table. "
        f"Tacstd2 versus granulocytes at the line unit is ρ = {s['line_granulocyte']['rho']:.3f}, "
        f"partial ρ = {s['line_granulocyte']['partial']:.3f} "
        f"(line-resample CI {fmt_ci(s['line_granulocyte']['partial_lo'], s['line_granulocyte']['partial_hi'])}). "
        f"That interval is wide because seventeen lines are collinear on Tacstd2 and Cldn4 "
        f"(line-level ρ = {s['gene_correlation']['line_rho']:.3f})."
    )
    lines.append("")
    lines.append("## ICI response")
    lines.append("")
    lines.append(
        f"Treated mice with both genes: {s['n_treated_mice']} mice, {s['n_pure_arms']} pure arms, "
        f"{s['n_mixed_arms']} mixed arms, {int(ici['n'])} series. "
        f"GSE124821 contributes {s['n_arms_gse124821']} of the {s['n_pure_arms']} pure arms, so an arm-level correlation "
        f"gives one mammary series nineteen votes. The study-level partial resamples series and is the unit used here."
    )
    lines.append("")
    lines.append(
        f"Study-level Tacstd2 versus responder fraction: ρ = {ici['rho']:.3f} "
        f"(CI {fmt_ci(ici['rho_lo'], ici['rho_hi'])}), partial given Cldn4 ρ = {ici['partial']:.3f} "
        f"(CI {fmt_ci(ici['partial_lo'], ici['partial_hi'])}), Δ = {ici['delta']:.3f} "
        f"(CI {fmt_ci(ici['delta_lo'], ici['delta_hi'])}). "
        f"{shrinkage_sentence(ici['rho'], ici['partial'], ici['delta_lo'], ici['delta_hi'])} "
        f"Cldn4 versus the same responder fraction is ρ = {s['ici_study_cldn4']['rho']:.3f}, "
        f"partial given Tacstd2 ρ = {s['ici_study_cldn4']['partial']:.3f}."
    )
    lines.append("")
    arm = s["ici_arm_tacstd2"]
    lines.append(
        f"Arm-level Tacstd2, descriptive because of the GSE124821 pile-up: ρ = {arm['rho']:.3f}, "
        f"partial ρ = {arm['partial']:.3f} (series-cluster CI {fmt_ci(arm['partial_lo'], arm['partial_hi'])}). "
        f"Median-split OR for response given high Tacstd2 is {fmt(arm['OR'])} "
        f"(Woolf {fmt_ci(arm['OR_low'], arm['OR_high'])}; table high-R/high-NR/low-R/low-NR = {arm['OR_table']}). "
        f"Stratifying that split on the Cldn4 median gives a Mantel–Haenszel OR of {fmt(arm['MH_OR'])} "
        f"(Robins–Breslow–Greenland {fmt_ci(arm['MH_low'], arm['MH_high'])}; stratum tables {arm['MH_tables']}). "
        f"The crude OR and the stratified OR both sit near 1."
    )
    lines.append("")
    lines.append("Mouse-level labels inside the only mixed arms:")
    lines.append("")
    lines.append("| Gene | Arm | Mice R vs NR | Median R | Median NR | Cliff δ |")
    lines.append("|---|---|---|---|---|---|")
    for row in s["mixed_detail"]:
        lines.append(
            f"| {row['gene']} | {row['stem']} | {row['n_R']} vs {row['n_NR']} | {row['median_R']:.3f} | {row['median_NR']:.3f} | {row['cliff']:.3f} |"
        )
    lines.append("")
    lines.append("## What this does not say")
    lines.append("")
    lines.append("- Partialling Cldn4 is an observational split of two correlated bulk measurements. It is not evidence that Cldn4 mediates a Tacstd2 effect, and it is not a knockdown.")
    lines.append("- The other immune scores in the TISMO export were not scanned for the largest attenuation.")
    lines.append("- Mouse-level analytic p-values are not reported as evidence. Mice are nested in lines and series. The mouse intervals above are study-cluster intervals.")
    lines.append("- The locked Tacstd2 49/64 count was not recomputed. LLC is not a KL line. No KP lung line is in this table.")
    lines.append("- Cutpoints for the response odds ratio are the medians. They were not tuned.")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "Expression is the TISMO `value` column from `POST /rtismo/gene/downVivoExprn` for Tacstd2 and Cldn4. "
        "Infiltrate scores are `POST /rtismo/gene/downICBTreated` with `type=3`, joined on the sample id. "
        "The files are the same export used for the study-level granulocyte table. "
        "A partial Spearman correlation is the Pearson correlation of ranks after each rank is residualized on the control rank, with an intercept. "
        "On z-scored ranks the simple coefficient equals the Spearman correlation, and the coefficient in the two-gene model is "
        "`(ρ_xy − ρ_xz ρ_yz) / (1 − ρ_xz²)`. "
        "The semipartial correlation correlates the control-residualized exposure rank with the raw outcome rank. "
        "Bootstrap intervals are percentile intervals from 2000 resamples. "
        "Study rows are resampled as rows. Mouse rows are resampled by series (every mouse from a drawn series), which stops one 91-mouse series from being treated as 91 independent tumors. "
        "Within-line intervals resample lines. "
        "Leave-one-series-out is the influence check. "
        "The response odds ratio is the median split. The stratified odds ratio is the Mantel–Haenszel estimate across the Cldn4 median split, with the Robins–Breslow–Greenland interval. "
        "An analytic partial-correlation p-value is stored in the table. It is not the interval interpreted above. The bootstrap interval is. "
        f"Seed sequence salt is `{SEED}` plus the analysis name."
    )
    lines.append("")
    lines.append("Code: `scripts/tismo_partial_cldn4.py`. Tables: `results/tismo_partial_cldn4/tables/`.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    tac = load_gene(DATA / "vivo" / "Tacstd2.csv", "Tacstd2")
    cld = load_gene(DATA / "vivo" / "Cldn4.csv", "Cldn4")
    meta_cols = ["Samples", "stem", "line", "cancer", "Responder", "Baseline", "GSE_ID", "Mouse_treatment", "cell_line"]
    left = tac[meta_cols + ["Tacstd2"]]
    right = cld[["Samples", "Cldn4", "Baseline", "GSE_ID", "Responder"]].rename(
        columns={"Baseline": "Baseline_cld", "GSE_ID": "GSE_cld", "Responder": "Responder_cld"}
    )
    paired = left.merge(right, on="Samples", how="inner")
    if int((paired["Baseline"] != paired["Baseline_cld"]).sum()) or int((paired["GSE_ID"] != paired["GSE_cld"]).sum()):
        raise RuntimeError("Tacstd2 and Cldn4 metadata disagree on shared samples")
    immune = load_immune(FEATURES)
    paired = paired.merge(immune, left_on="Samples", right_index=True, how="left")

    cld_base = cld[cld["Baseline"] == 1]
    tac_base_ids = set(tac.loc[tac["Baseline"] == 1, "Samples"])
    dropped = cld_base[~cld_base["Samples"].isin(tac_base_ids)]
    dropped_detail = ", ".join(
        f"{gse} n={len(g)}" for gse, g in dropped.groupby("GSE_ID")
    ) or "none"

    frames = build_frames(paired, FEATURES)
    # Gene-gene correlation and its intervals, including the line unit.
    gene_rows = {}
    for unit, cluster in (
        ("study", None),
        ("mouse", "GSE_ID"),
        ("line", None),
    ):
        frame = frames[unit]
        est = estimate_arrays(frame["Tacstd2"], frame["Cldn4"], frame["Tacstd2"])
        # estimate_arrays partials y on z; for the pairwise gene correlation, Spearman is enough.
        rho = safe_corr(rankdata(frame["Tacstd2"]), rankdata(frame["Cldn4"]))
        boot = bootstrap_unit(frame, "Tacstd2", "Cldn4", ("Cldn4",), cluster, f"gene-{unit}")
        # The partial of Cldn4 on itself is undefined. Use the marginal draws, which are the gene-gene correlation.
        loo_min = loo_max = float("nan")
        if unit == "study":
            vals = []
            for cid in frame["unit_id"]:
                sub = frame[frame["unit_id"] != cid]
                vals.append(safe_corr(rankdata(sub["Tacstd2"]), rankdata(sub["Cldn4"])))
            loo_min, loo_max = float(np.nanmin(vals)), float(np.nanmax(vals))
        gene_rows[unit] = {
            "rho": float(rho),
            "n": int(len(frame)),
            "rho_lo": boot["Cldn4"]["rho_lo"],
            "rho_hi": boot["Cldn4"]["rho_hi"],
            "loo_min": loo_min,
            "loo_max": loo_max,
            "unused_partial_guard": est["n"],
        }

    unit_jobs = [
        ("study", "study", None, "unit_id", "study"),
        ("study_n_ge_8", "study_n_ge_8", None, "unit_id", "study-ge8"),
        ("mouse", "mouse", "GSE_ID", "GSE_ID", "mouse-study"),
        ("mouse_without_GSE124821", "mouse_without_GSE124821", "GSE_ID", "GSE_ID", "mouse-nobig"),
        ("line", "line", None, "unit_id", "line"),
        ("within_line", "within_line", "line", "line", "within"),
    ]
    spec_rows = []
    for unit_name, frame_key, boot_cluster, loo_cluster, seed_key in unit_jobs:
        frame = frames[frame_key]
        for exposure, control in (("Tacstd2", "Cldn4"), ("Cldn4", "Tacstd2")):
            scored = score_direction(frame, exposure, control, FEATURES, boot_cluster, loo_cluster, f"{seed_key}-{exposure}")
            for rec in scored:
                rec["unit"] = unit_name
                spec_rows.append(rec)
    specs = pd.DataFrame(spec_rows)
    specs.to_csv(TABLES / "partial_specs.tsv", sep="\t", index=False)

    frames["study"].to_csv(TABLES / "baseline_study_medians.tsv", sep="\t", index=False)
    mouse_keep = ["Samples", "GSE_ID", "line", "cancer", "Tacstd2", "Cldn4", *FEATURES]
    frames["mouse"][mouse_keep].to_csv(TABLES / "baseline_mouse.tsv", sep="\t", index=False)

    resp = response_tables(paired)
    resp["pure"].to_csv(TABLES / "treated_arms.tsv", sep="\t", index=False)
    resp["studies"].to_csv(TABLES / "treated_study.tsv", sep="\t", index=False)
    ici_rows = score_response(resp["pure"], resp["studies"])
    ici = pd.DataFrame(ici_rows)
    ici.to_csv(TABLES / "response_partial.tsv", sep="\t", index=False)
    pd.DataFrame(resp["mixed_detail"]).to_csv(TABLES / "mixed_arm_mice.tsv", sep="\t", index=False)

    def pack(unit, exposure, feature) -> dict:
        return json_ready(lookup(specs, unit, exposure, feature).to_dict())

    # Unpaired Cldn4 granulocyte correlation, to show what the extra series was worth.
    cld_only = cld[cld["Baseline"] == 1].merge(immune, left_on="Samples", right_index=True, how="left")
    unpaired_rows = []
    for gse, g in cld_only.groupby("GSE_ID"):
        unpaired_rows.append({"Cldn4": float(g["Cldn4"].median()), "y": float(g[PRIMARY_GRAN].median())})
    unpaired = pd.DataFrame(unpaired_rows).dropna()
    unpaired_rho = safe_corr(rankdata(unpaired["Cldn4"]), rankdata(unpaired["y"]))

    g = pack("study", "Tacstd2", PRIMARY_GRAN)
    rev_g = pack("study", "Cldn4", PRIMARY_GRAN)
    influence_hit = frames["study"].loc[frames["study"]["unit_id"] == g["loo_influence_id"]]
    if len(influence_hit) == 1:
        hit = influence_hit.iloc[0]
        g["influence_label"] = f"{g['loo_influence_id']} ({hit['line']} {str(hit['cancer']).lower()}, {int(hit['n_mice'])} mice)"
    plot_added_variable(
        frames["study"],
        FIGS / "study_granulocyte_partial.png",
        gene_rows["study"]["rho"],
        g["rho"],
        g["partial"],
        g["partial_lo"],
        g["partial_hi"],
    )
    plot_slope(specs, "study", FIGS / "study_marginal_vs_partial.png", "Study-level Spearman ρ")
    plot_slope(specs, "mouse", FIGS / "mouse_marginal_vs_partial.png", "Mouse-level Spearman ρ")
    ici_tac = next(row for row in ici_rows if row["unit"] == "study" and row["exposure"] == "Tacstd2")
    plot_response(resp["studies"], FIGS / "study_response.png", ici_tac["rho"], ici_tac["partial"])

    within = specs[(specs.unit == "within_line") & (specs.exposure == "Tacstd2")].copy()
    within["abs_rho"] = within["rho"].abs()
    within_max = within.loc[within["abs_rho"].idxmax()]

    base = frames["mouse"]
    n_single = int((frames["study"]["n_lines"] == 1).sum())
    summary = {
        "gene_correlation": {
            "study_rho": gene_rows["study"]["rho"],
            "study_n": gene_rows["study"]["n"],
            "study_rho_lo": gene_rows["study"]["rho_lo"],
            "study_rho_hi": gene_rows["study"]["rho_hi"],
            "study_loo_min": gene_rows["study"]["loo_min"],
            "study_loo_max": gene_rows["study"]["loo_max"],
            "mouse_rho": gene_rows["mouse"]["rho"],
            "mouse_n": gene_rows["mouse"]["n"],
            "mouse_rho_lo": gene_rows["mouse"]["rho_lo"],
            "mouse_rho_hi": gene_rows["mouse"]["rho_hi"],
            "line_rho": gene_rows["line"]["rho"],
        },
        "vif_study": 1.0 / (1.0 - gene_rows["study"]["rho"] ** 2),
        "n_studies": int(frames["study"]["unit_id"].nunique()),
        "n_studies_single_line": n_single,
        "n_lines": int(frames["line"]["unit_id"].nunique()),
        "n_lines_within": int(frames["within_line"]["line"].nunique()),
        "gse124821_baseline_mice": int((base["GSE_ID"] == "GSE124821").sum()),
        "gse159344_baseline_mice": int((base["GSE_ID"] == "GSE159344").sum()),
        "dropped_baseline_mice": int(len(dropped)),
        "dropped_baseline_detail": dropped_detail,
        "unpaired_cldn4_granulocyte_rho": float(unpaired_rho),
        "unpaired_cldn4_granulocyte_n": int(len(unpaired)),
        "n_treated_mice": int(len(resp["treated"])),
        "n_pure_arms": int(len(resp["pure"])),
        "n_mixed_arms": int((resp["arms"]["label"] == "mixed").sum()),
        "n_arms_gse124821": int((resp["pure"]["GSE_ID"] == "GSE124821").sum()),
        "lung_lines": sorted(base.loc[base["cancer"] == "Lung", "line"].unique()),
        "primary_granulocyte_study": g,
        "reverse_granulocyte_study": rev_g,
        "granulocyte_study_ge8": pack("study_n_ge_8", "Tacstd2", PRIMARY_GRAN),
        "neutrophil_study": pack("study", "Tacstd2", PRIMARY_NEUT),
        "primary_neutrophil_mouse": pack("mouse", "Tacstd2", PRIMARY_NEUT),
        "reverse_neutrophil_mouse": pack("mouse", "Cldn4", PRIMARY_NEUT),
        "neutrophil_mouse_no_big": pack("mouse_without_GSE124821", "Tacstd2", PRIMARY_NEUT),
        "granulocyte_mouse": pack("mouse", "Tacstd2", PRIMARY_GRAN),
        "reverse_granulocyte_mouse": pack("mouse", "Cldn4", PRIMARY_GRAN),
        "within_neutrophil": pack("within_line", "Tacstd2", PRIMARY_NEUT),
        "within_granulocyte": pack("within_line", "Tacstd2", PRIMARY_GRAN),
        "primary_cd8_study": pack("study", "Tacstd2", PRIMARY_CD8),
        "reverse_cd8_study": pack("study", "Cldn4", PRIMARY_CD8),
        "primary_cd8_mouse": pack("mouse", "Tacstd2", PRIMARY_CD8),
        "within_cd8": pack("within_line", "Tacstd2", PRIMARY_CD8),
        "mouse_cd8_xcell": pack("mouse", "Tacstd2", "T CD8_xCell"),
        "line_granulocyte": pack("line", "Tacstd2", PRIMARY_GRAN),
        "within_max": json_ready(within_max.to_dict()),
        "ici_study_tacstd2": json_ready(ici_tac),
        "ici_study_cldn4": json_ready(next(row for row in ici_rows if row["unit"] == "study" and row["exposure"] == "Cldn4")),
        "ici_arm_tacstd2": json_ready(next(row for row in ici_rows if row["unit"] == "arm" and row["exposure"] == "Tacstd2")),
        "mixed_detail": resp["mixed_detail"],
        "specs_df_note": specs,
    }
    # The report writer needs the specs frame. JSON gets a cleaned copy without that frame.
    write_report(summary, OUT / "REPORT.md")
    summary.pop("specs_df_note")
    (OUT / "summary.json").write_text(json.dumps(json_ready(summary), indent=2) + "\n", encoding="utf-8")

    # Locked marginals from the earlier TISMO infiltrate table, on this same export.
    cd8_locked = lookup(specs, "study", "Tacstd2", PRIMARY_CD8)
    if abs(cd8_locked["rho"] - 0.738) > 0.01:
        raise RuntimeError(f"study Tacstd2–CIBERSORT CD8 ρ drifted from 0.738: {cd8_locked['rho']}")
    if "KL" in set(base["line"]) or "KP" in set(base["line"]):
        raise RuntimeError("KL or KP labeled in TISMO lines")
    if int((paired["line"] == "LLC").sum()) < 1:
        raise RuntimeError("LLC missing from the paired table")
    if abs(g["partial"] - (g["rho"] - g["rho_xz"] * g["rho_yz"]) / math.sqrt((1 - g["rho_xz"] ** 2) * (1 - g["rho_yz"] ** 2))) > 1e-8:
        raise RuntimeError("partial correlation formula does not match the stored partial")
    if not g.get("loo_influence_id"):
        raise RuntimeError("study granulocyte leave-one-out influence was not recorded")
    for feat in ("T CD8_TIMER", "T CD8_EPIC", "T CD8_xCell", "CD8 T_mMCPcounter"):
        row = lookup(specs, "study", "Tacstd2", feat)
        if not (row["partial_lo"] < 0 < row["partial_hi"]):
            raise RuntimeError(f"CD8 report text assumes {feat} partial CI includes 0")
    quant = lookup(specs, "study", "Tacstd2", "T CD8_quanTIseq")
    if not (quant["delta_lo"] > 0 and quant["partial_lo"] < 0 < quant["partial_hi"]):
        raise RuntimeError("quanTIseq CD8 sentence does not match its intervals")
    if not (cd8_locked["partial_lo"] > 0 and cd8_locked["delta_lo"] < 0 < cd8_locked["delta_hi"]):
        raise RuntimeError("CIBERSORT CD8 sentence does not match its intervals")
    print(json.dumps(json_ready({
        "gene_study": summary["gene_correlation"]["study_rho"],
        "gran_study": [g["rho"], g["partial"], g["delta"], g["partial_lo"], g["partial_hi"], g["delta_lo"], g["delta_hi"]],
        "gran_reverse": [rev_g["rho"], rev_g["partial"], rev_g["partial_lo"], rev_g["partial_hi"]],
        "neut_mouse": [summary["primary_neutrophil_mouse"]["rho"], summary["primary_neutrophil_mouse"]["partial"], summary["primary_neutrophil_mouse"]["partial_lo"], summary["primary_neutrophil_mouse"]["partial_hi"]],
        "cd8_study": [summary["primary_cd8_study"]["rho"], summary["primary_cd8_study"]["partial"], summary["primary_cd8_study"]["partial_lo"], summary["primary_cd8_study"]["partial_hi"]],
        "ici": [ici_tac["rho"], ici_tac["partial"], ici_tac["partial_lo"], ici_tac["partial_hi"]],
    }), indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()

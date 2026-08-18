#!/usr/bin/env python3
"""Concordant-4 CLDN4-only: T/NK infiltrate + malignant IFN/MHC DE.

Datasets: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071 / GSE127465 / GSE207422 / GSE154826. Not dual-high.

1) Patient-level Spearman of malignant CLDN4 %pos vs same-unit T/NK.
   Also Q4 vs Q1. Pooled n and leave-one-dataset-out.
2) Patient-pseudobulk OLS on log2(TMM-CPM+1), malignant Q4 vs Q1.
   IFN (Hallmark α∪γ), MHC-I/APM, chemokine, TJ.

Inferential unit = patient / donor / sample. Do not quote cell counts as n.
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lib_stats import random_effects_dl, spearman, stouffer  # noqa: E402

DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
COMBO = "+".join(COHORTS)
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MIN_N_SINGLE = 4
Q4_MIN_N = 16

# PR #459 / #320 pairwise rows — comparison only, not re-audited.
GIVEN_PAIRS = [
    {
        "source": "PR320_given",
        "combo": "GSE131907+GSE205335",
        "score": "pct",
        "analysis": "q4q1",
        "n": 23,
        "n_q1": 12,
        "n_q4": 11,
        "rho": float("nan"),
        "p": float("nan"),
        "r_rb": -0.705,
        "p_q4q1": 0.000301,
        "note": "PR #320 / #459 given Q4 vs Q1 author %pos; not re-audited",
    },
    {
        "source": "PR459_given",
        "combo": "GSE123902+GSE189357",
        "score": "pct",
        "analysis": "spearman",
        "n": 22,
        "rho": -0.638,
        "p": 0.003,
        "r_rb": float("nan"),
        "p_q4q1": float("nan"),
        "note": "PR #459 pair; not re-audited",
    },
    {
        "source": "PR459_given",
        "combo": "GSE123902+GSE131907",
        "score": "pct",
        "analysis": "spearman",
        "n": 34,
        "rho": -0.575,
        "p": 0.001,
        "r_rb": -0.700,
        "p_q4q1": 0.012,
        "note": "PR #459 pair; not re-audited",
    },
    {
        "source": "PR459_given",
        "combo": "GSE131907+GSE189357",
        "score": "pct",
        "analysis": "spearman",
        "n": 30,
        "rho": -0.542,
        "p": 0.003,
        "r_rb": -0.619,
        "p_q4q1": 0.042,
        "note": "PR #459 pair; not re-audited",
    },
    {
        "source": "PR459_given",
        "combo": "GSE123902+GSE205335",
        "score": "pct",
        "analysis": "spearman",
        "n": 35,
        "rho": -0.522,
        "p": 0.002,
        "r_rb": -0.802,
        "p_q4q1": 0.005,
        "note": "PR #459 pair; not re-audited",
    },
    {
        "source": "PR459_given",
        "combo": "GSE131907+GSE205335",
        "score": "pct",
        "analysis": "spearman",
        "n": 43,
        "rho": -0.479,
        "p": 0.002,
        "r_rb": -0.702,
        "p_q4q1": 0.006,
        "note": "PR #459 pair; not re-audited",
    },
    {
        "source": "PR459_given",
        "combo": "GSE189357+GSE205335",
        "score": "pct",
        "analysis": "spearman",
        "n": 31,
        "rho": -0.478,
        "p": 0.009,
        "r_rb": -0.750,
        "p_q4q1": 0.010,
        "note": "PR #459 pair; not re-audited",
    },
]

CHEMOKINE = [
    "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCL16", "CXCL8", "CX3CL1",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
    "XCL1", "XCL2", "IL15", "IL2", "IL7", "IL21",
    "CXCR3", "CXCR4", "CXCR5", "CCR5", "CCR7",
]
FOCAL_BOX = [
    "CLDN4", "STAT1", "IRF1", "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2",
    "CXCL9", "CXCL10", "CCL5", "OCLN", "TJP1", "KRT8", "KRT18", "KRT19",
]
FAM_COLORS = {
    "IFN": "#d62728",
    "MHC-I/APM": "#1f77b4",
    "TJ": "#2ca02c",
    "chemokine": "#ff7f0e",
}
COHORT_COLORS = {
    "GSE123902": "#4c78a8",
    "GSE131907": "#f58518",
    "GSE205335": "#54a24b",
    "GSE189357": "#b279a2",
    COMBO: "#111111",
}
NEG = "#7a2d0b"
POS = "#4a4a4a"


def load_families() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    tj = (
        set(sets["KEGG_TIGHT_JUNCTION"])
        | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
        | set(g for g in a8["focal_genes"] if g not in sets["KRT_EPITHELIAL"])
    )
    tj.discard("CLDN4")
    return {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC-I/APM": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "TJ": tj,
        "chemokine": set(CHEMOKINE),
    }


def family_of(gene: str, fam: dict[str, set[str]]) -> str:
    hits = [k for k, vs in fam.items() if gene in vs]
    return "|".join(hits) if hits else "other"


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def q4_vs_q1(cldn4, immune, min_n: int = Q4_MIN_N) -> dict | None:
    s = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    n = int(len(s))
    if n < min_n:
        return None
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = s.loc[qs == "Q1", "i"]
    q4 = s.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 3 or n4 < 3:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "r_rb": float(r_rb),
        "p": float(p),
    }


def load_units() -> dict[str, dict]:
    units: dict[str, dict] = {}

    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    units["GSE123902"] = {
        "cohort": "GSE123902",
        "malig_def": "marker_malig",
        "unit": "donor",
        "n": int(len(el)),
        "mean": el["mal_CLDN4_mean"].to_numpy(float),
        "pct": el["mal_CLDN4_pct"].to_numpy(float),
        "tnk": el["frac_tnk"].to_numpy(float),
        "ids": el["patient"].astype(str).tolist(),
        "note": "Laughney 2020 tumor/met donors; marker-malignant; normals dropped",
    }

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(TUMOR_ORIGINS)]
    mal = tumor[tumor["n_malignant"] >= 20].copy()
    units["GSE131907"] = {
        "cohort": "GSE131907",
        "malig_def": "author_malig",
        "unit": "sample",
        "n": int(len(mal)),
        "mean": mal["mal_CLDN4_mean"].to_numpy(float),
        "pct": mal["mal_CLDN4_pct"].to_numpy(float),
        "tnk": mal["frac_tnk"].to_numpy(float),
        "ids": mal["sample"].astype(str).tolist(),
        "note": "Kim 2020; author Malignant; tumor-bearing sites; n_mal>=20; sample-level",
    }

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    units["GSE205335"] = {
        "cohort": "GSE205335",
        "malig_def": "author_malig",
        "unit": "patient",
        "n": int(len(d)),
        "mean": d["mal_CLDN4_mean"].to_numpy(float),
        "pct": d["mal_CLDN4_pct_pos"].to_numpy(float),
        "tnk": d["frac_tnk"].to_numpy(float),
        "ids": d["patient"].astype(str).tolist(),
        "note": "advanced ICI; author malignant; all subtypes; RECIST not required",
    }

    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    units["GSE189357"] = {
        "cohort": "GSE189357",
        "malig_def": "marker_malig",
        "unit": "patient",
        "n": int(len(el)),
        "mean": el["mal_CLDN4_mean"].to_numpy(float),
        "pct": el["mal_CLDN4_pct"].to_numpy(float),
        "tnk": el["frac_tnk"].to_numpy(float),
        "ids": el["patient"].astype(str).tolist(),
        "note": "Zhu/Wang AIS–IAC; 9 patients; marker-malignant",
    }
    return units


def four_meta(units: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for c in COHORTS:
        u = units[c]
        pct = np.asarray(u["pct"], float)
        # GSE123902 / GSE189357 store fraction; others percent.
        if np.nanmax(pct) <= 1.5:
            pct = pct * 100.0
        q = assign_quartiles(pd.Series(pct, index=u["ids"]))
        for i, pid in enumerate(u["ids"]):
            rows.append(
                {
                    "patient": pid,
                    "cohort": c,
                    "unit": u["unit"],
                    "malig_def": u["malig_def"],
                    "cldn4_pct": float(pct[i]),
                    "cldn4_mean": float(u["mean"][i]),
                    "frac_tnk": float(u["tnk"][i]),
                    "quartile": str(q.loc[pid]),
                }
            )
    meta = pd.DataFrame(rows)
    in_mat = set()
    for name in COHORTS:
        path = DATA / f"{name}_malignant_counts.tsv.gz"
        cols = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str)
        in_mat.update(cols.tolist())
    meta["in_count_matrix"] = meta["patient"].isin(in_mat)
    return meta


def nlung_sensitivity() -> dict:
    """nLung cannot enter malignant CLDN4 (n_mal=0). Epithelial %pos vs T/NK only."""
    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    nlung = d[d["origin"] == "nLung"].copy()
    tumor = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    epi_tumor = tumor[tumor["n_epithelial"] >= 20].copy()
    epi_nlung = nlung[nlung["n_epithelial"] >= 20].copy()
    both = pd.concat([epi_tumor, epi_nlung], ignore_index=True)
    rho_t, p_t, n_t = spearman(epi_tumor["epi_CLDN4_pct"], epi_tumor["frac_tnk"])
    rho_b, p_b, n_b = spearman(both["epi_CLDN4_pct"], both["frac_tnk"])
    rho_n, p_n, n_n = spearman(epi_nlung["epi_CLDN4_pct"], epi_nlung["frac_tnk"])
    return {
        "n_nlung": int(len(nlung)),
        "n_nlung_malignant": int((nlung["n_malignant"] > 0).sum()),
        "nlung_note": "all nLung have n_malignant=0; cannot join malignant CLDN4 vs T/NK",
        "epi_tumor_n": n_t,
        "epi_tumor_rho": rho_t,
        "epi_tumor_p": p_t,
        "epi_plus_nlung_n": n_b,
        "epi_plus_nlung_rho": rho_b,
        "epi_plus_nlung_p": p_b,
        "nlung_only_n": n_n,
        "nlung_only_rho": rho_n,
        "nlung_only_p": p_n,
    }


def _within_cohort_ranks(members: list[dict], score: str):
    xs, ys, labels = [], [], []
    for u in members:
        x = np.asarray(u[score], float)
        y = np.asarray(u["tnk"], float)
        m = np.isfinite(x) & np.isfinite(y)
        x, y = x[m], y[m]
        if x.size == 0:
            continue
        r = pd.Series(x).rank(method="average").to_numpy()
        r = (r - 1.0) / max(len(r) - 1, 1)
        xs.append(r)
        ys.append(y)
        labels.extend([u["cohort"]] * len(y))
    if not xs:
        return None, None, None
    return np.concatenate(xs), np.concatenate(ys), labels


def _q4_empty() -> dict:
    return {
        "n_q1": np.nan,
        "n_q4": np.nan,
        "n_compared": np.nan,
        "r_rb": np.nan,
        "p_q4q1": np.nan,
        "delta_median": np.nan,
        "thin_q4": True,
        "median_q1": np.nan,
        "median_q4": np.nan,
    }


def _attach_q4(rec: dict, x, y, min_n: int = Q4_MIN_N) -> dict:
    q = q4_vs_q1(x, y, min_n=min_n)
    if not q:
        rec.update(_q4_empty())
        return rec
    thin = q["n_q1"] < 6 or q["n_q4"] < 6 or abs(q["r_rb"]) >= 0.999
    rec.update(
        {
            "n_q1": q["n_q1"],
            "n_q4": q["n_q4"],
            "n_compared": q["n_compared"],
            "r_rb": q["r_rb"],
            "p_q4q1": q["p"],
            "delta_median": q["delta_median"],
            "thin_q4": thin,
            "median_q1": q["median_q1"],
            "median_q4": q["median_q4"],
        }
    )
    return rec


def stacked_q4q1(members: list[dict], score: str) -> dict | None:
    """Within-cohort quartiles, then MWU on concatenated Q1 vs Q4 T/NK."""
    q1, q4 = [], []
    for u in members:
        x = np.asarray(u[score], float)
        y = np.asarray(u["tnk"], float)
        m = np.isfinite(x) & np.isfinite(y)
        x, y = x[m], y[m]
        if x.size < 4:
            continue
        qs = assign_quartiles(pd.Series(x))
        q1.extend(y[qs.values == "Q1"].tolist())
        q4.extend(y[qs.values == "Q4"].tolist())
    n1, n4 = len(q1), len(q4)
    if n1 < 3 or n4 < 3:
        return None
    u_stat, p = stats.mannwhitneyu(np.asarray(q4), np.asarray(q1), alternative="two-sided")
    r_rb = (2.0 * float(u_stat)) / (n4 * n1) - 1.0
    return {
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "r_rb": float(r_rb),
        "p_q4q1": float(p),
        "median_q1": float(np.median(q1)),
        "median_q4": float(np.median(q4)),
        "delta_median": float(np.median(q4) - np.median(q1)),
        "thin_q4": n1 < 6 or n4 < 6 or abs(r_rb) >= 0.999,
    }


def pool_combo(members: list[dict], score: str, kind: str) -> dict | None:
    names = [u["cohort"] for u in members]
    rhos, ps, ns = [], [], []
    for u in members:
        rho, p, n = spearman(u[score], u["tnk"])
        if n < MIN_N_SINGLE or not math.isfinite(rho):
            return None
        rhos.append(rho)
        ps.append(p)
        ns.append(n)
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    n_total = int(sum(ns))
    rec = {
        "kind": kind,
        "k": len(members),
        "combo": "+".join(names),
        "score": score,
        "n": n_total,
        "rho": re.get("pooled_rho", float("nan")),
        "p": re.get("p", float("nan")),
        "I2": re.get("I2", float("nan")),
        "ci95_lo": re.get("ci95_rho", [float("nan"), float("nan")])[0],
        "ci95_hi": re.get("ci95_rho", [float("nan"), float("nan")])[1],
        "stouffer_z": st.get("z", float("nan")),
        "stouffer_p": st.get("p", float("nan")),
        "member_rhos": ",".join(f"{c}:{r:.3f}" for c, r in zip(names, rhos)),
        "member_ns": ",".join(f"{c}:{n}" for c, n in zip(names, ns)),
        "note": "DerSimonian–Laird Fisher-z of cohort Spearmans; n=sum of units",
        "source": "this_pr",
    }
    xr, yr, _ = _within_cohort_ranks(members, score)
    rec = _attach_q4(rec, xr, yr)
    rec["q4_method"] = "within_cohort_rank_then_qcut"
    stacked = stacked_q4q1(members, score)
    if stacked:
        rec["stacked_n_q1"] = stacked["n_q1"]
        rec["stacked_n_q4"] = stacked["n_q4"]
        rec["stacked_n_compared"] = stacked["n_compared"]
        rec["stacked_r_rb"] = stacked["r_rb"]
        rec["stacked_p_q4q1"] = stacked["p_q4q1"]
        rec["stacked_delta_median"] = stacked["delta_median"]
    return rec


def atomic(u: dict, score: str) -> dict:
    rho, p, n = spearman(u[score], u["tnk"])
    rec = {
        "kind": "single",
        "k": 1,
        "combo": u["cohort"],
        "score": score,
        "n": n,
        "rho": rho,
        "p": p,
        "I2": 0.0,
        "malig_def": u["malig_def"],
        "unit": u["unit"],
        "note": u["note"],
        "source": "this_pr_context",
        "member_rhos": f"{u['cohort']}:{rho:.3f}" if math.isfinite(rho) else "",
        "member_ns": f"{u['cohort']}:{n}",
    }
    rec = _attach_q4(rec, u[score], u["tnk"], min_n=8)
    return rec


def tmm_norm_factors(counts: pd.DataFrame) -> pd.Series:
    lib = counts.sum(axis=0).astype(float).replace(0, np.nan)
    rel = counts.div(lib, axis=1)
    f75 = rel.quantile(0.75, axis=0)
    ref = (f75 - f75.mean()).abs().idxmin()
    ref_c = counts[ref].astype(float)
    ref_lib = float(lib[ref])
    factors = {}
    for col in counts.columns:
        obs = counts[col].astype(float)
        obs_lib = float(lib[col])
        keep = (obs > 0) & (ref_c > 0)
        if keep.sum() < 50:
            factors[col] = 1.0
            continue
        m = np.log2((obs[keep] / obs_lib) / (ref_c[keep] / ref_lib))
        a = 0.5 * np.log2((obs[keep] / obs_lib) * (ref_c[keep] / ref_lib))
        w = (obs_lib - obs[keep]) / (obs_lib * obs[keep]) + (ref_lib - ref_c[keep]) / (
            ref_lib * ref_c[keep]
        )
        ok = np.isfinite(m) & np.isfinite(a) & np.isfinite(w) & (w > 0)
        m, a, w = m[ok], a[ok], w[ok]
        if len(m) < 50:
            factors[col] = 1.0
            continue
        lo_m, hi_m = np.quantile(m, [0.30, 0.70])
        lo_a, hi_a = np.quantile(a, [0.05, 0.95])
        trim = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if trim.sum() < 20:
            factors[col] = 1.0
            continue
        tmm = float(np.average(m[trim], weights=1.0 / w[trim]))
        factors[col] = 2 ** tmm
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def ols_de(logcpm: pd.DataFrame, design: pd.DataFrame, coef: str) -> pd.DataFrame:
    X = design.reindex(logcpm.columns).astype(float)
    if X.isna().any().any():
        keep = ~X.isna().any(axis=1)
        X = X.loc[keep]
        Y = logcpm.loc[:, keep].astype(float).to_numpy()
    else:
        Y = logcpm.astype(float).to_numpy()
    Xv = X.to_numpy()
    n, p = Xv.shape
    df_res = n - p
    if df_res < 1:
        return pd.DataFrame()
    xtx = Xv.T @ Xv
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ Xv.T @ Y.T
    resid = Y.T - Xv @ beta
    sse = np.sum(resid**2, axis=0)
    sigma2 = sse / df_res
    j = list(X.columns).index(coef)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[j, j], 0.0))
    est = beta[j]
    t = np.divide(est, se, out=np.zeros_like(est), where=se > 0)
    pv = 2.0 * stats.t.sf(np.abs(t), df_res)
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "AveExpr": Y.mean(axis=1),
            "t": t,
            "p": pv,
            "se": se,
            "df": df_res,
        }
    )
    out["fdr"] = _bh(out["p"].values)
    return out.sort_values("p")


def _add_cohort_dummies(design: pd.DataFrame, m: pd.DataFrame) -> pd.DataFrame:
    idx = m.set_index("patient")
    for c in COHORTS:
        if c == REF:
            continue
        design[f"cohort_{c}"] = (idx["cohort"] == c).astype(float)
    return design


def _design_q4q1(m: pd.DataFrame, cohort: str | None) -> pd.DataFrame:
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_Q4"] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    if cohort is None and m["cohort"].nunique() > 1:
        design = _add_cohort_dummies(design, m)
    return design


def run_q4q1(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None):
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
    m = m.loc[m["patient"].isin(counts.columns)].copy()
    n1 = int((m["quartile"] == "Q1").sum())
    n4 = int((m["quartile"] == "Q4").sum())
    info = {
        "n_q1": n1,
        "n_q4": n4,
        "n": n1 + n4,
        "patients_q1": ",".join(sorted(m.loc[m["quartile"] == "Q1", "patient"])),
        "patients_q4": ",".join(sorted(m.loc[m["quartile"] == "Q4", "patient"])),
    }
    if n1 < 3 or n4 < 3:
        return pd.DataFrame(), info, pd.DataFrame()
    cts = filter_genes(counts.loc[:, m["patient"]])
    lc = log_cpm(cts, tmm_norm_factors(cts))
    de = ols_de(lc, _design_q4q1(m, cohort), "CLDN4_Q4")
    info["n_genes"] = int(len(de))
    return de, info, lc


def run_continuous(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None):
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    m = m.loc[m["patient"].isin(counts.columns)].copy()
    info = {"n": int(len(m)), "n_q1": np.nan, "n_q4": np.nan}
    if len(m) < 8:
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    lc = log_cpm(cts, tmm_norm_factors(cts))
    z = m.set_index("patient")["cldn4_pct"].astype(float)
    z = (z - z.mean()) / z.std(ddof=1)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_pct_z"] = z
    if cohort is None and m["cohort"].nunique() > 1:
        design = _add_cohort_dummies(design, m)
    de = ols_de(lc, design, "CLDN4_pct_z")
    info["n_genes"] = int(len(de))
    return de, info


def family_score_de(logcpm, meta, fam, cohort, split):
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    if split == "q4q1":
        m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
    m = m.loc[m["patient"].isin(logcpm.columns)].copy()
    if split == "q4q1" and ((m["quartile"] == "Q1").sum() < 3 or (m["quartile"] == "Q4").sum() < 3):
        return pd.DataFrame()
    if split == "continuous" and len(m) < 8:
        return pd.DataFrame()
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    if split == "q4q1":
        coef = "CLDN4_Q4"
        design[coef] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    else:
        coef = "CLDN4_pct_z"
        z = m.set_index("patient")["cldn4_pct"].astype(float)
        design[coef] = (z - z.mean()) / z.std(ddof=1)
    if cohort is None and m["cohort"].nunique() > 1:
        design = _add_cohort_dummies(design, m)
    rows = []
    for name, genes in fam.items():
        present = [g for g in sorted(genes) if g in logcpm.index]
        if len(present) < 3:
            continue
        score = logcpm.loc[present, m["patient"]].astype(float).mean(axis=0).to_frame().T
        score.index = [name]
        de = ols_de(score, design, coef)
        if de.empty:
            continue
        r = de.iloc[0]
        rows.append(
            {
                "family": name,
                "split": split,
                "cohort": cohort or COMBO,
                "n_q1": int((m["quartile"] == "Q1").sum()) if split == "q4q1" else np.nan,
                "n_q4": int((m["quartile"] == "Q4").sum()) if split == "q4q1" else np.nan,
                "n": int(len(m)),
                "n_genes": int(len(present)),
                "logFC": float(r.logFC),
                "p": float(r.p),
                "t": float(r.t),
                "se": float(r.se),
                "df": float(r.df),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = _bh(out["p"].values)
    return out


def _fam_summary(de: pd.DataFrame, fam: dict[str, set[str]], n_q1: int, n_q4: int) -> list[dict]:
    rows = []
    for name, genes in fam.items():
        sub = de[de["gene"].isin(genes)].copy()
        if sub.empty:
            rows.append({"family": name, "n_tested": 0, "n_p05": 0, "n_fdr05": 0, "n_up": 0, "n_down": 0})
            continue
        sig = sub[sub["p"] < 0.05]
        rows.append(
            {
                "family": name,
                "n_tested": int(len(sub)),
                "n_p05": int((sub["p"] < 0.05).sum()),
                "n_fdr05": int((sub["fdr"] < 0.05).sum()),
                "n_up": int((sig["logFC"] > 0).sum()),
                "n_down": int((sig["logFC"] < 0).sum()),
                "median_logFC": float(sub["logFC"].median()),
                "mean_logFC": float(sub["logFC"].mean()),
                "top_gene": str(sub.iloc[0]["gene"]),
                "top_logFC": float(sub.iloc[0]["logFC"]),
                "top_p": float(sub.iloc[0]["p"]),
                "top_fdr": float(sub.iloc[0]["fdr"]),
                "n_q1": n_q1,
                "n_q4": n_q4,
            }
        )
    return rows


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(float(x))):
        return "—"
    return f"{float(x):.{nd}f}"


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and (not math.isfinite(float(p)))):
        return "—"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def forest_effect(rows, title, out, effect_key="rho", p_key="p", n_key="n"):
    rows = [r for r in rows if math.isfinite(float(r.get(effect_key, float("nan"))))]
    if not rows:
        return
    rows = sorted(rows, key=lambda r: float(r[effect_key]))
    fig, ax = plt.subplots(figsize=(9.4, max(2.8, 0.42 * len(rows) + 1.2)))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        e = float(r[effect_key])
        col = NEG if e < 0 else POS
        ax.plot(e, i, "o", color=col, ms=7, zorder=3)
        ax.hlines(i, 0, e, color=col, lw=1.5, zorder=2)
        lab = f"n={r.get(n_key, r.get('n'))}"
        p = r.get(p_key, float("nan"))
        if p == p and math.isfinite(float(p)):
            lab += f"  p={float(p):.3g}"
        ax.text(0.02 if e >= 0 else -0.02, i, lab, va="center", ha="left" if e >= 0 else "right", fontsize=7)
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([r["combo"] for r in rows], fontsize=8)
    ax.set_xlabel(title)
    ax.set_xlim(-1.05, 1.05)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def scatter_tnk(meta: pd.DataFrame, path: Path, pooled: dict):
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for c in COHORTS:
        g = meta[meta["cohort"] == c]
        ax.scatter(g["cldn4_pct"], g["frac_tnk"], s=36, c=COHORT_COLORS[c], label=f"{c} n={len(g)}", alpha=0.9)
    ax.set_xlabel("malignant CLDN4 %pos (within-cohort scale)")
    ax.set_ylabel("same-unit T/NK fraction")
    ax.set_title(
        f"concordant-4  n={pooled['n']}  Fisher-z ρ={pooled['rho']:.3f}  p={_fmt_p(pooled['p'])}"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def box_tnk(meta: pd.DataFrame, path: Path):
    m = meta[meta["quartile"].isin(["Q1", "Q4"])].copy()
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    q1 = m.loc[m["quartile"] == "Q1", "frac_tnk"]
    q4 = m.loc[m["quartile"] == "Q4", "frac_tnk"]
    ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)
    rng = np.random.default_rng(0)
    for vals, x, col, cohort in (
        (q1, 1, "#4c78a8", "Q1"),
        (q4, 2, "#e45756", "Q4"),
    ):
        ax.scatter(x + rng.uniform(-0.08, 0.08, len(vals)), vals, s=18, c=col, zorder=3)
    u, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (len(q4) * len(q1)) - 1.0
    ax.set_ylabel("T/NK fraction")
    ax.set_title(f"within-cohort Q4 vs Q1 T/NK  r={r_rb:.3f}  p={p:.3g}")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar_tnk(units: dict, path: Path):
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    names = COHORTS
    ns = [units[c]["n"] for c in names]
    ax.bar(names, ns, color=[COHORT_COLORS[c] for c in names])
    ax.set_ylabel("honest n (donor / sample / patient)")
    ax.set_title(f"concordant-4 T/NK units  N={sum(ns)}")
    for i, n in enumerate(ns):
        ax.text(i, n + 0.3, str(n), ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def volcano(de, fam, title, path):
    if de.empty:
        return
    fig, ax = plt.subplots(figsize=(7.4, 5.5))
    ax.scatter(de["logFC"], -np.log10(np.clip(de["p"], 1e-300, 1)), s=8, c="#c8c8c8", linewidths=0, alpha=0.7)
    for fam_name, color in FAM_COLORS.items():
        sub = de[de["gene"].isin(fam[fam_name])]
        if sub.empty:
            continue
        ax.scatter(
            sub["logFC"],
            -np.log10(np.clip(sub["p"], 1e-300, 1)),
            s=22,
            c=color,
            linewidths=0.2,
            edgecolors="white",
            label=fam_name,
            zorder=3,
        )
    if "CLDN4" in set(de["gene"]):
        hit = de[de["gene"] == "CLDN4"]
        ax.scatter(hit["logFC"], -np.log10(np.clip(hit["p"], 1e-300, 1)), s=60, c="#111", marker="D", label="CLDN4", zorder=4)
    ax.axhline(-np.log10(0.05), ls="--", c="#666", lw=0.7)
    ax.axvline(0, ls="-", c="#999", lw=0.6)
    ax.set_xlabel("log2FC (CLDN4-high − low)")
    ax.set_ylabel("−log10 p (OLS)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def heatmap_family(logcpm, meta, genes, title, path):
    genes = [g for g in genes if g in logcpm.index]
    if len(genes) < 4 or meta.empty:
        return
    m = meta.loc[meta["patient"].isin(logcpm.columns)].copy().sort_values(["cohort", "cldn4_pct"])
    mat = logcpm.loc[genes, m["patient"]].astype(float)
    z = mat.sub(mat.mean(axis=1), axis=0)
    z = z.div(mat.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(max(8.0, 0.26 * len(m) + 2.4), max(4.8, 0.22 * len(genes) + 1.6)))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.2, vmax=2.2)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xticks(range(len(m)))
    ax.set_xticklabels([f"{r.patient}\n{r.quartile}" for r in m.itertuples()], fontsize=5.5, rotation=90)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="row z")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def box_key_genes(logcpm, meta, title, path):
    genes = [g for g in FOCAL_BOX if g in logcpm.index]
    if not genes:
        return
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(logcpm.columns)].copy()
    ncols = 4
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.6 * nrows), squeeze=False)
    rng = np.random.default_rng(0)
    for i, gene in enumerate(genes):
        ax = axes[i // ncols][i % ncols]
        q1 = logcpm.loc[gene, m.loc[m["quartile"] == "Q1", "patient"]].astype(float)
        q4 = logcpm.loc[gene, m.loc[m["quartile"] == "Q4", "patient"]].astype(float)
        ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(q1)), q1, s=14, c="#4c78a8", zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(q4)), q4, s=14, c="#e45756", zorder=3)
        if len(q1) >= 3 and len(q4) >= 3:
            _, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
            ax.set_title(f"{gene}  p={p:.3g}", fontsize=9)
        else:
            ax.set_title(gene, fontsize=9)
        ax.set_ylabel("log2(CPM+1)")
    for j in range(len(genes), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(title, y=1.01)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def box_family_scores(logcpm, meta, fam, title, path):
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(logcpm.columns)].copy()
    names = [k for k in FAM_COLORS if any(g in logcpm.index for g in fam[k])]
    fig, axes = plt.subplots(1, len(names), figsize=(2.6 * len(names), 3.6), squeeze=False)
    rng = np.random.default_rng(1)
    for i, name in enumerate(names):
        ax = axes[0][i]
        present = [g for g in fam[name] if g in logcpm.index]
        score = logcpm.loc[present, m["patient"]].astype(float).mean(axis=0)
        q1 = score.loc[m.loc[m["quartile"] == "Q1", "patient"]]
        q4 = score.loc[m.loc[m["quartile"] == "Q4", "patient"]]
        ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], widths=0.55)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(q1)), q1, s=16, c="#4c78a8", zorder=3)
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(q4)), q4, s=16, c="#e45756", zorder=3)
        if len(q1) >= 3 and len(q4) >= 3:
            _, p = stats.mannwhitneyu(q4, q1, alternative="two-sided")
            ax.set_title(f"{name}\np={p:.3g}", fontsize=9)
        else:
            ax.set_title(name, fontsize=9)
        ax.set_ylabel("mean log2(CPM+1)")
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def forest_family_scores(fam_de, title, path):
    if fam_de.empty:
        return
    plot = fam_de.sort_values("family").copy()
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    y = np.arange(len(plot))
    ax.axvline(0, c="#666", lw=0.7)
    ax.errorbar(plot["logFC"], y, xerr=1.96 * plot["se"], fmt="o", color="#333", ecolor="#888", elinewidth=1.2, capsize=3)
    ax.scatter(plot["logFC"], y, c=[FAM_COLORS[f] for f in plot["family"]], s=40, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.family}  n={int(r.n)}  p={r.p:.3g}" for r in plot.itertuples()], fontsize=8)
    ax.set_xlabel("family-score log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def forest_family_scores_stacked(fam_de, title, path):
    sub = fam_de[fam_de["split"] == "q4q1"].copy()
    if sub.empty:
        return
    order_c = COHORTS + [COMBO]
    order_f = list(FAM_COLORS)
    sub["cohort"] = pd.Categorical(sub["cohort"], categories=order_c, ordered=True)
    sub["family"] = pd.Categorical(sub["family"], categories=order_f, ordered=True)
    sub = sub.dropna(subset=["cohort", "family"]).sort_values(["family", "cohort"])
    fig, ax = plt.subplots(figsize=(8.4, max(5.4, 0.30 * len(sub) + 1.4)))
    y = np.arange(len(sub))
    ax.axvline(0, c="#666", lw=0.7)
    ax.errorbar(sub["logFC"], y, xerr=1.96 * sub["se"], fmt="none", ecolor="#888", elinewidth=1.1, capsize=2.5)
    ax.scatter(
        sub["logFC"],
        y,
        c=[COHORT_COLORS.get(str(c), "#333") for c in sub["cohort"]],
        s=[42 if c == COMBO else 28 for c in sub["cohort"]],
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.family} · {r.cohort}  n={int(r.n_q1)}/{int(r.n_q4)}  p={r.p:.3g}" for r in sub.itertuples()],
        fontsize=7,
    )
    ax.set_xlabel("family-score log2FC (CLDN4-high − low)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def strip_cldn4(meta, path):
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.6))
    cols = {"Q1": "#4c78a8", "Q2": "#9ecae1", "Q3": "#f4a582", "Q4": "#e45756"}
    for ax, cohort in zip(axes, COHORTS):
        g = meta.loc[meta["cohort"] == cohort].sort_values("cldn4_pct")
        edge = ["#111" if (not bool(v)) else "white" for v in g["in_count_matrix"]]
        ax.scatter(range(len(g)), g["cldn4_pct"], c=[cols[q] for q in g["quartile"]], s=36, edgecolors=edge, linewidths=0.6)
        ax.set_title(f"{cohort}  n={len(g)}  {g['unit'].iloc[0]}")
        ax.set_ylabel("malignant CLDN4 %pos")
        ax.set_xlabel("units (sorted)")
        for q, c in cols.items():
            ax.scatter([], [], c=c, label=q, s=24)
        ax.legend(frameon=False, fontsize=6, ncol=2, loc="best")
    fig.suptitle("concordant-4 CLDN4-only quartiles (within-cohort; black edge = not in malignant matrix)")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def n_bar_de(inv, path):
    sub = inv[inv["split"] == "q4q1"].copy()
    fig, ax = plt.subplots(figsize=(8.8, 3.8))
    x = np.arange(len(sub))
    ax.bar(x - 0.18, sub["n_low"].fillna(0), 0.36, label="CLDN4-low / Q1", color="#4c78a8")
    ax.bar(x + 0.18, sub["n_high"].fillna(0), 0.36, label="CLDN4-high / Q4", color="#e45756")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["contrast"].tolist(), rotation=28, ha="right", fontsize=7)
    ax.set_ylabel("honest n (donor / sample / patient)")
    ax.set_title("Malignant DE sample sizes — concordant-4")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(units, singles, pooled, loos, nlung, meta, inv, fam_rows, fam_score, headline, cldn4_row):
    n_by = {c: units[c]["n"] for c in COHORTS}
    q_by = {c: meta[meta["cohort"] == c]["quartile"].value_counts() for c in COHORTS}
    drop = meta.loc[~meta["in_count_matrix"], "patient"].tolist()
    drop_txt = ", ".join(drop) if drop else "none"
    p_pct = next(r for r in pooled if r["score"] == "pct")
    p_mean = next(r for r in pooled if r["score"] == "mean")

    def loo_md(score):
        lines = [
            "| dropped | remaining | N | ρ (p, I²) | stacked Q4 vs Q1 r (n_Q1/n_Q4, p) |",
            "|---|---|---:|---|---|",
        ]
        for r in loos:
            if r["score"] != score:
                continue
            lines.append(
                f"| {r['loo'].replace('drop_','')} | {r['combo']} | {r['n']} | "
                f"{_fmt(r['rho'])} ({_fmt_p(r['p'])}, I²={_fmt(r['I2'],1)}%) | "
                f"{_fmt(r.get('stacked_r_rb'))} ({r.get('stacked_n_q1','—')}/{r.get('stacked_n_q4','—')}, {_fmt_p(r.get('stacked_p_q4q1'))}) |"
            )
        return "\n".join(lines)

    def single_md(score):
        lines = [
            "| cohort | unit | n | ρ | p | Q4 r (n_Q1/n_Q4, p) |",
            "|---|---|---:|---:|---:|---|",
        ]
        for r in singles:
            if r["score"] != score:
                continue
            lines.append(
                f"| {r['combo']} | {r['unit']} | {r['n']} | {_fmt(r['rho'])} | {_fmt_p(r['p'])} | "
                f"{_fmt(r.get('r_rb'))} ({r.get('n_q1','—')}/{r.get('n_q4','—')}, {_fmt_p(r.get('p_q4q1'))}) |"
            )
        return "\n".join(lines)

    def given_md():
        lines = [
            "| source | combo | analysis | N | effect | p |",
            "|---|---|---|---:|---|---|",
        ]
        for r in GIVEN_PAIRS:
            if r["analysis"] == "q4q1":
                eff = f"r={_fmt(r['r_rb'])}"
                p = _fmt_p(r["p_q4q1"])
            else:
                eff = f"ρ={_fmt(r['rho'])}"
                p = _fmt_p(r["p"])
            lines.append(f"| {r['source']} | {r['combo']} | {r['analysis']} {r['score']} | {r['n']} | {eff} | {p} |")
        return "\n".join(lines)

    def fam_md(contrast):
        sub = fam_rows[fam_rows["contrast"] == contrast]
        if sub.empty:
            return "_No family rows._"
        lines = [
            "| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | mean logFC | top gene (logFC, p, FDR) |",
            "|---|---:|---|---:|---:|---:|---|",
        ]
        for r in sub.itertuples():
            lines.append(
                f"| {r.family} | {r.n_tested} | {r.n_p05} ({r.n_up}/{r.n_down}) | {r.n_fdr05} | "
                f"{r.median_logFC:+.3f} | {r.mean_logFC:+.3f} | {r.top_gene} ({r.top_logFC:+.3f}, {_fmt_p(r.top_p)}, {_fmt_p(r.top_fdr)}) |"
            )
        return "\n".join(lines)

    def score_md(split, cohort=COMBO):
        sub = fam_score[(fam_score["split"] == split) & (fam_score["cohort"] == cohort)]
        if sub.empty:
            return "_No family-score DE._"
        lines = [
            "| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |",
            "|---|---:|---|---:|---:|---:|---:|",
        ]
        for r in sub.itertuples():
            nq = "—" if split != "q4q1" else f"{int(r.n_q1)} / {int(r.n_q4)}"
            lines.append(
                f"| {r.family} | {int(r.n)} | {nq} | {int(r.n_genes)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        return "\n".join(lines)

    def score_md_cohorts(split):
        sub = fam_score[(fam_score["split"] == split) & (fam_score["cohort"] != COMBO)]
        if sub.empty:
            return "_No per-cohort family-score DE._"
        lines = [
            "| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |",
            "|---|---|---:|---|---:|---:|---:|---:|",
        ]
        for r in sub.itertuples():
            nq = "—" if split != "q4q1" else f"{int(r.n_q1)} / {int(r.n_q4)}"
            lines.append(
                f"| {r.cohort} | {r.family} | {int(r.n)} | {nq} | {int(r.n_genes)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        return "\n".join(lines)

    mal = inv[inv["contrast"] == "malignant_q4q1_combined"]
    mal_c = inv[inv["contrast"] == "malignant_continuous_combined"]
    n_de = int(mal.iloc[0]["n"]) if len(mal) else 0
    n1 = int(mal.iloc[0]["n_low"]) if len(mal) else 0
    n4 = int(mal.iloc[0]["n_high"]) if len(mal) else 0
    n_cont = int(mal_c.iloc[0]["n"]) if len(mal_c) else 0

    text = f"""# Concordant-4 CLDN4-only — T/NK + malignant IFN/MHC

ADDITIVE. **CLDN4-only.** No TACSTD2∩CLDN4 dual-high. Not a mega-merge.
The four sets that already point the same way: **GSE123902 + GSE131907 +
GSE205335 + GSE189357**. Not GSE148071, GSE127465, GSE207422, GSE154826,
or CD45+/T-only extracts.

Thesis (already correct; not re-derived): CLDN4-high malignant cells have
lower own IFN/MHC-I, and patients have lower T/NK. CLDN4-low/KD opens IFN/MHC.
Matching extras here are IFN/MHC **DOWN** in CLDN4-high. Do not sell a
DoRothEA IFN-up-in-high as the KD direction.

The **new number is the four-set pooled n**. Pairwise rhos from PR #459 / #320
are comparison rows only and were not re-audited.

Patient / donor / sample is the unit. Do not quote cell counts as n.
GSE123902 = donor. GSE131907 = sample (tumor-bearing). GSE205335 = patient
(RECIST not required). GSE189357 = patient. p-values are descriptive.

## 1. T/NK — four-set pooled n

Primary score = malignant CLDN4 **%pos**. Mean is the matching extra.
Pooling = DerSimonian–Laird on Fisher-z of the four cohort Spearmans.
Q4 vs Q1 primary = **within-cohort quartiles stacked**, then Mann–Whitney
on T/NK (rank-biserial r). That uses the same Q labels as the DE.

| score | k | N | ρ (p, I², 95% CI) | stacked Q4 vs Q1 r (n_Q1/n_Q4, p) |
|---|---:|---:|---|---|
| %pos | 4 | {p_pct['n']} | {_fmt(p_pct['rho'])} ({_fmt_p(p_pct['p'])}, I²={_fmt(p_pct['I2'],1)}%, {_fmt(p_pct.get('ci95_lo'))} to {_fmt(p_pct.get('ci95_hi'))}) | {_fmt(p_pct.get('stacked_r_rb'))} ({p_pct.get('stacked_n_q1')}/{p_pct.get('stacked_n_q4')}, {_fmt_p(p_pct.get('stacked_p_q4q1'))}) |
| mean | 4 | {p_mean['n']} | {_fmt(p_mean['rho'])} ({_fmt_p(p_mean['p'])}, I²={_fmt(p_mean['I2'],1)}%, {_fmt(p_mean.get('ci95_lo'))} to {_fmt(p_mean.get('ci95_hi'))}) | {_fmt(p_mean.get('stacked_r_rb'))} ({p_mean.get('stacked_n_q1')}/{p_mean.get('stacked_n_q4')}, {_fmt_p(p_mean.get('stacked_p_q4q1'))}) |

Honest N = {n_by['GSE123902']} donors + {n_by['GSE131907']} samples + {n_by['GSE205335']} patients + {n_by['GSE189357']} patients = **{p_pct['n']} units**.

Within-cohort rank-then-qcut Q4 (combo-enum method) on %pos: r={_fmt(p_pct.get('r_rb'))}, n_Q1/n_Q4={p_pct.get('n_q1')}/{p_pct.get('n_q4')}, p={_fmt_p(p_pct.get('p_q4q1'))}.

### Leave-one-dataset-out (%pos)

{loo_md("pct")}

LOO that drops GSE189357 recovers the PR #459 triple membership (n=56). That
is a robustness row of this four-set, not a re-audit of the triple rho.

### Singles (context for the pool; not the new number)

{single_md("pct")}

### Comparison rows (PR #459 / #320; not re-audited)

{given_md()}

### GSE131907 nLung sensitivity (not primary)

Primary T/NK uses tumor-bearing sites only ({", ".join(sorted(TUMOR_ORIGINS))};
n_mal≥20). All {nlung['n_nlung']} nLung samples have **n_malignant=0**, so they
cannot join a malignant-CLDN4 vs T/NK test.

Epithelial CLDN4 %pos vs T/NK (sensitivity only):

| cut | n | ρ | p |
|---|---:|---:|---:|
| tumor-bearing epi (n_epi≥20) | {nlung['epi_tumor_n']} | {_fmt(nlung['epi_tumor_rho'])} | {_fmt_p(nlung['epi_tumor_p'])} |
| tumor-bearing epi + nLung | {nlung['epi_plus_nlung_n']} | {_fmt(nlung['epi_plus_nlung_rho'])} | {_fmt_p(nlung['epi_plus_nlung_p'])} |
| nLung epi only | {nlung['nlung_only_n']} | {_fmt(nlung['nlung_only_rho'])} | {_fmt_p(nlung['nlung_only_p'])} |

## 2. Tumor-cell-intrinsic — malignant Q4 vs Q1

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same collapse (patient × cell-type UMI-sum → bulk DE).
Within-cohort malignant CLDN4 %pos Q4 vs Q1, then stacked with cohort
covariates. Positive logFC = higher in CLDN4-high.

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | the four concordant sets | not the 7-pool; not +GSE148071 |
| Units | %pos N={p_pct['n']} | mixed donor / sample / patient |
| Split | within-cohort %pos Q4 vs Q1 | mid quartiles unused in binary DE |
| Malignant | marker-malig UMI-sum (123902, 189357); author-malig (131907, 205335) | marker gate ≠ author CNV |
| Model | ~ cohort + CLDN4_Q4 | small n; no muscat mixed model |
| Families | IFN Hallmark α∪γ; MHC-I/APM custom; chemokine panel; TJ KEGG/GO (CLDN4 held out) | chemokine panel is compact |

Honest DE n (units in the count matrices):

- GSE123902: n={n_by['GSE123902']} donors. Q1={int(q_by['GSE123902'].get('Q1', 0))} Q4={int(q_by['GSE123902'].get('Q4', 0))}.
- GSE131907: n={n_by['GSE131907']} samples. Q1={int(q_by['GSE131907'].get('Q1', 0))} Q4={int(q_by['GSE131907'].get('Q4', 0))}.
- GSE205335: n={n_by['GSE205335']} patients. Q1={int(q_by['GSE205335'].get('Q1', 0))} Q4={int(q_by['GSE205335'].get('Q4', 0))}. Out of malignant DE: {drop_txt}.
- GSE189357: n={n_by['GSE189357']} patients. Q1={int(q_by['GSE189357'].get('Q1', 0))} Q4={int(q_by['GSE189357'].get('Q4', 0))}.

| contrast | n_low / n_high | n_genes | note |
|---|---|---:|---|
| Q4 vs Q1 stacked | {n1} / {n4} | {int(mal.iloc[0]['n_genes']) if len(mal) else 0} | cohort covariates; malignant only |
| continuous stacked | n={n_cont} | {int(mal_c.iloc[0]['n_genes']) if len(mal_c) else 0} | CLDN4 %pos z |

Do not quote N={p_pct['n']} as the DE n: stacked Q4 vs Q1 uses only the tails
that are in the count matrices (n={n_de}).

### Family scores (stacked Q4 vs Q1)

Family score = mean log2(TMM-CPM+1) of family genes present.
Expected under the thesis: IFN down, MHC-I/APM down, chemokine down, TJ up.

{score_md("q4q1")}

Continuous family-score DE (units in the count matrices):

{score_md("continuous")}

Per-cohort family scores (within-cohort Q4 vs Q1; GSE123902 tails are thin).
GSE189357 Q4 n=2: single-cohort binary family DE skipped (need ≥3 each);
TD6 and TD9 stay in the stacked n={n_de}. P4001 is Q1 on the T/NK vector
(stacked T/NK n_Q1=19) but is not in the malignant UMI-sum (DE n_Q1={n1}).

{score_md_cohorts("q4q1")}

CLDN4 itself (held out of TJ) stacked Q4 vs Q1: logFC={cldn4_row.get('logFC', float('nan')):+.3f}, p={_fmt_p(cldn4_row.get('p'))}, n_Q1={cldn4_row.get('n_q1', '—')}, n_Q4={cldn4_row.get('n_q4', '—')} — direction check on the split gene.

### Gene-level family members (stacked Q4 vs Q1)

{fam_md("q4q1_combined")}

Headline genes (family members only, lowest p):

{headline.get("md", "_none_")}

## What this is not

- Not a mega-merge of every lung scRNA set.
- Not GSE148071 / GSE127465 / GSE207422 / GSE154826 / CD45+ or T-only.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a re-audit of PR #459 / #320 pairwise rhos.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not DoRothEA / VIPER IFN activity sold as the KD direction.
- Not evidence that CLDN4 *causes* the T/NK or IFN/MHC change.
- Genome-wide FDR on these n is expected to be thin; family scores and
  gene-count direction (up/down) are the claim.

## Files

- `tables/tnk_pooled.tsv` — **headline T/NK table** (honest n)
- `tables/tnk_loo.tsv` — leave-one-dataset-out
- `tables/tnk_units.tsv` / `tables/tnk_singles.tsv`
- `tables/tnk_nlung_sensitivity.tsv`
- `tables/tnk_comparison_pr459_320.tsv` — given pairwise rows
- `tables/family_de.tsv` / `tables/family_summary.tsv` — **headline IFN/MHC table**
- `tables/de_q4q1_combined_families.tsv` / `tables/de_families.tsv` / `tables/de_all.tsv`
- `tables/n_honest.tsv` / `tables/sample_inventory.tsv`
- `figures/` — T/NK scatter / forests / Q4 box; malignant volcano / heatmap / family forests

Reproduce:

```bash
python3 methods/concordant4_123902_131907_205335_189357_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)


def run_tnk(units):
    singles, pooled, loos = [], [], []
    members = [units[c] for c in COHORTS]
    for score in ("pct", "mean"):
        for c in COHORTS:
            singles.append(atomic(units[c], score))
        rec = pool_combo(members, score, kind="concordant4")
        if rec:
            rec["note"] = "FOUR-SET POOLED n (new number). " + rec["note"]
            pooled.append(rec)
        for drop in members:
            keep = [u for u in members if u["cohort"] != drop["cohort"]]
            loo = pool_combo(keep, score, kind="loo")
            if loo is None:
                continue
            loo["loo"] = f"drop_{drop['cohort']}"
            loo["combo"] = "+".join(u["cohort"] for u in keep)
            loo["note"] = f"leave-one-dataset-out of the concordant-4; dropped {drop['cohort']}"
            loos.append(loo)
    return singles, pooled, loos


def run_de(meta, fam):
    parts = {}
    all_de, inv_rows, fam_rows, fam_score_rows = [], [], [], []
    for cohort in COHORTS:
        path = DATA / f"{cohort}_malignant_counts.tsv.gz"
        if not path.exists():
            raise SystemExit(f"missing {path}")
        cts = read_counts(path)
        parts[cohort] = cts
        de, info, _ = run_q4q1(cts, meta, cohort)
        inv_rows.append(
            {
                "contrast": f"malignant_q4q1_{cohort}",
                "compartment": "malignant",
                "cohort": cohort,
                "split": "q4q1",
                "n_low": info["n_q1"],
                "n_high": info["n_q4"],
                "n": info["n"],
                "n_genes": info.get("n_genes", 0),
                "patients_q1": info.get("patients_q1", ""),
                "patients_q4": info.get("patients_q4", ""),
            }
        )
        if not de.empty:
            de = de.copy()
            de["compartment"] = "malignant"
            de["cohort"] = cohort
            de["contrast"] = f"q4q1_{cohort}"
            de["n_q1"] = info["n_q1"]
            de["n_q4"] = info["n_q4"]
            de["family"] = de["gene"].map(lambda g: family_of(g, fam))
            all_de.append(de)
            for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
                fam_rows.append({**r, "compartment": "malignant", "contrast": f"q4q1_{cohort}"})
        de_c, info_c = run_continuous(cts, meta, cohort)
        inv_rows.append(
            {
                "contrast": f"malignant_continuous_{cohort}",
                "compartment": "malignant",
                "cohort": cohort,
                "split": "continuous",
                "n_low": np.nan,
                "n_high": np.nan,
                "n": info_c["n"],
                "n_genes": info_c.get("n_genes", 0),
                "patients_q1": "",
                "patients_q4": "",
            }
        )

    genes = sorted(set.intersection(*[set(p.index) for p in parts.values()]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts.values()], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()]

    de, info, _ = run_q4q1(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_q4q1_combined",
            "compartment": "malignant",
            "cohort": COMBO,
            "split": "q4q1",
            "n_low": info["n_q1"],
            "n_high": info["n_q4"],
            "n": info["n"],
            "n_genes": info.get("n_genes", 0),
            "patients_q1": info.get("patients_q1", ""),
            "patients_q4": info.get("patients_q4", ""),
        }
    )
    de_c, info_c = run_continuous(combined, meta, None)
    inv_rows.append(
        {
            "contrast": "malignant_continuous_combined",
            "compartment": "malignant",
            "cohort": COMBO,
            "split": "continuous",
            "n_low": np.nan,
            "n_high": np.nan,
            "n": info_c["n"],
            "n_genes": info_c.get("n_genes", 0),
            "patients_q1": "",
            "patients_q4": "",
        }
    )

    cldn4_row = {}
    if not de.empty:
        de = de.copy()
        de["compartment"] = "malignant"
        de["cohort"] = COMBO
        de["contrast"] = "q4q1_combined"
        de["n_q1"] = info["n_q1"]
        de["n_q4"] = info["n_q4"]
        de["family"] = de["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de)
        for r in _fam_summary(de, fam, info["n_q1"], info["n_q4"]):
            fam_rows.append({**r, "compartment": "malignant", "contrast": "q4q1_combined"})
        volcano(
            de,
            fam,
            f"malignant Q4 vs Q1  n={info['n_q4']} vs {info['n_q1']}  (stacked, cohort covariates)",
            FIGS / "volcano_malignant_q4q1_combined",
        )
        hit = de[de["gene"] == "CLDN4"]
        if not hit.empty:
            cldn4_row = {
                "logFC": float(hit.iloc[0]["logFC"]),
                "p": float(hit.iloc[0]["p"]),
                "n_q1": info["n_q1"],
                "n_q4": info["n_q4"],
            }

    if not de_c.empty:
        de_c = de_c.copy()
        de_c["compartment"] = "malignant"
        de_c["cohort"] = COMBO
        de_c["contrast"] = "continuous_combined"
        de_c["n_q1"] = np.nan
        de_c["n_q4"] = np.nan
        de_c["family"] = de_c["gene"].map(lambda g: family_of(g, fam))
        all_de.append(de_c)

    m_all = meta.loc[meta["patient"].isin(combined.columns)].copy()
    cts_all = filter_genes(combined.loc[:, m_all["patient"]])
    lc_all = log_cpm(cts_all, tmm_norm_factors(cts_all))
    for cohort in (None, *COHORTS):
        for split in ("q4q1", "continuous"):
            fs = family_score_de(lc_all, meta, fam, cohort, split)
            if not fs.empty:
                fam_score_rows.append(fs)

    m_q = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(combined.columns)]
    if len(m_q) >= 6:
        cts_q = filter_genes(combined.loc[:, m_q["patient"]])
        lc = log_cpm(cts_q, tmm_norm_factors(cts_q))
        show = []
        prefer = {
            "IFN": ["STAT1", "IRF1", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IFI6", "BST2", "GBP1"],
            "MHC-I/APM": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9", "TAPBP"],
            "TJ": ["CLDN1", "CLDN7", "OCLN", "TJP1", "F11R", "PARD3", "CDH1", "MARVELD2"],
            "chemokine": CHEMOKINE[:12],
        }
        for name in FAM_COLORS:
            show.extend([g for g in prefer[name] if g in lc.index])
        heatmap_family(lc, m_q, show, "malignant family genes, Q1/Q4 only (row z)", FIGS / "heatmap_malignant_families_q4q1")
        box_key_genes(lc, m_q, "malignant key genes, CLDN4 Q4 vs Q1", FIGS / "box_malignant_key_genes")
        box_family_scores(lc, m_q, fam, "malignant family scores, CLDN4 Q4 vs Q1", FIGS / "box_malignant_family_scores")

    inv = pd.DataFrame(inv_rows)
    fam_score = pd.concat(fam_score_rows, ignore_index=True) if fam_score_rows else pd.DataFrame()
    fam_df = pd.DataFrame(fam_rows)
    headline = {"md": "_DE table empty._"}
    if all_de:
        de_all = pd.concat(all_de, ignore_index=True)
        de_all.to_csv(TABLES / "de_all.tsv", sep="\t", index=False)
        fam_de = de_all[de_all["family"] != "other"].copy()
        fam_de.to_csv(TABLES / "de_families.tsv", sep="\t", index=False)
        head = fam_de[fam_de["contrast"] == "q4q1_combined"].sort_values("p")
        head.to_csv(TABLES / "de_q4q1_combined_families.tsv", sep="\t", index=False)
        lines = [
            "| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
        for r in head.head(12).itertuples():
            lines.append(
                f"| {r.family} | {r.gene} | {int(r.n_q1)} | {int(r.n_q4)} | "
                f"{r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        headline = {"md": "\n".join(lines)}
    return inv, fam_df, fam_score, headline, cldn4_row


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    units = load_units()
    print("UNITS", {c: units[c]["n"] for c in COHORTS})
    meta = four_meta(units)
    meta.to_csv(TABLES / "tnk_units.tsv", sep="\t", index=False)
    meta.to_csv(TABLES / "sample_inventory.tsv", sep="\t", index=False)

    singles, pooled, loos = run_tnk(units)
    nlung = nlung_sensitivity()
    pd.DataFrame(singles).to_csv(TABLES / "tnk_singles.tsv", sep="\t", index=False)
    pd.DataFrame(pooled).to_csv(TABLES / "tnk_pooled.tsv", sep="\t", index=False)
    pd.DataFrame(loos).to_csv(TABLES / "tnk_loo.tsv", sep="\t", index=False)
    pd.DataFrame([nlung]).to_csv(TABLES / "tnk_nlung_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(GIVEN_PAIRS).to_csv(TABLES / "tnk_comparison_pr459_320.tsv", sep="\t", index=False)

    p_pct = next(r for r in pooled if r["score"] == "pct")
    scatter_tnk(meta, FIGS / "scatter_tnk_pct", p_pct)
    box_tnk(meta, FIGS / "box_tnk_q4q1")
    n_bar_tnk(units, FIGS / "n_honest_tnk")
    forest_effect(
        [{**r, "combo": f"drop {r['loo'].replace('drop_','')}"} for r in loos if r["score"] == "pct"]
        + [{**p_pct, "combo": f"FOUR-SET n={p_pct['n']}"}],
        "concordant-4 Fisher-z ρ  (malignant CLDN4 %pos vs T/NK)",
        FIGS / "forest_tnk_loo_rho",
    )
    q_rows = []
    for r in loos + pooled:
        if r["score"] != "pct" or not math.isfinite(float(r.get("stacked_r_rb", float("nan")))):
            continue
        q_rows.append(
            {
                "combo": r["combo"] if r["kind"] != "concordant4" else f"FOUR-SET n={r['n']}",
                "r_rb": r["stacked_r_rb"],
                "p_q4q1": r["stacked_p_q4q1"],
                "n": r.get("stacked_n_compared", r["n"]),
            }
        )
    forest_effect(q_rows, "stacked within-cohort Q4 vs Q1 rank-biserial r  (T/NK)", FIGS / "forest_tnk_q4q1", effect_key="r_rb", p_key="p_q4q1")

    fam = load_families()
    strip_cldn4(meta, FIGS / "cldn4_quartile_strip")
    inv, fam_df, fam_score, headline, cldn4_row = run_de(meta, fam)
    inv.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)
    n_bar_de(inv, FIGS / "n_honest_q4q1")
    if not fam_df.empty:
        fam_df.to_csv(TABLES / "family_summary.tsv", sep="\t", index=False)
    if not fam_score.empty:
        fam_score.to_csv(TABLES / "family_de.tsv", sep="\t", index=False)
        head_fs = fam_score[(fam_score["cohort"] == COMBO) & (fam_score["split"] == "q4q1")]
        forest_family_scores(head_fs, "malignant family scores, Q4 vs Q1 stacked", FIGS / "forest_family_scores_q4q1")
        forest_family_scores_stacked(fam_score, "within-cohort Q4 vs Q1 family scores, then stacked", FIGS / "forest_family_scores_stacked")

    write_finding(units, singles, pooled, loos, nlung, meta, inv, fam_df, fam_score, headline, cldn4_row)
    print("WROTE", HERE / "FINDING.md")
    print(pd.DataFrame(pooled)[["score", "n", "rho", "p", "I2", "stacked_r_rb", "stacked_n_compared"]].to_string(index=False))
    if not fam_score.empty:
        print(fam_score[(fam_score["cohort"] == COMBO) & (fam_score["split"] == "q4q1")].to_string(index=False))


if __name__ == "__main__":
    main()

           
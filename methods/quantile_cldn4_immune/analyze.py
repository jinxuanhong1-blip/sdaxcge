#!/usr/bin/env python3
"""ADDITIVE Q4 vs Q1 recut of existing public sample-level tables.

Does not download matrices or invent signatures. Quartiles are cut on the
anchor already stored on each harvested table; endpoints are existing columns.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
HARVEST = HERE / "harvested"
FIG = HERE / "figures"
TAB = HERE / "tables"
RNG = np.random.default_rng(20260817)
N_BOOT = 2000

# Canonical endpoint names used in FINDING / forest
ENDPOINTS = ["CD8A", "ImmuneScore", "GEP18", "IFN", "MHC-I", "OS"]


def _read(name: str, **kwargs) -> pd.DataFrame:
    p = HARVEST / name
    if name.endswith(".csv"):
        return pd.read_csv(p, **kwargs)
    return pd.read_csv(p, sep="\t", **kwargs)


def load_cohorts() -> dict[str, dict]:
    """Return {cohort: {df, n_matrix, anchors, endpoints, notes}}."""
    out: dict[str, dict] = {}

    g218 = _read("GSE218989_per_patient.tsv")
    out["GSE218989"] = {
        "df": g218,
        "n_matrix": 355,
        "layer": "RNA",
        "anchors": {"CLDN4": "CLDN4_log2tpm1", "TACSTD2": "TACSTD2_log2tpm1"},
        "endpoints": {
            "CD8A": "CD8A_log2tpm1",
            "ImmuneScore": "ESTIMATE_Immune_ssgsea",
            "IFN": "IFNG_Ayers6_meanz",
            "MHC-I": "MHCI_meanz",
        },
        "os": {"time": "os_days", "event": "death"},
        "notes": "GEP18 not on table (Ayers6 is IFN, not GEP18). ImmuneScore = ESTIMATE Immune ssGSEA.",
    }

    g285 = _read("GSE285029_sample_scores.csv")
    if "Unnamed: 0" in g285.columns:
        g285 = g285.rename(columns={"Unnamed: 0": "sample"})
    out["GSE285029"] = {
        "df": g285,
        "n_matrix": 234,
        "layer": "RNA",
        "anchors": {"CLDN4": "CLDN4", "TACSTD2": "TACSTD2"},
        "endpoints": {
            "CD8A": "CD8A",
            "ImmuneScore": "immune8",
            "GEP18": "GEP18",
            "IFN": "IFN_compact",
            "MHC-I": "MHC1",
        },
        "os": None,
        "notes": "ImmuneScore column is A1 8-gene (immune8); ESTIMATE ImmuneScore not on table. No OS.",
    }

    luad = _read("TCGA_LUAD_xena_star_tpm.tsv")
    out["TCGA-LUAD"] = {
        "df": luad,
        "n_matrix": 502,
        "layer": "RNA",
        "anchors": {"CLDN4": "CLDN4", "TACSTD2": "TACSTD2"},
        "endpoints": {
            "CD8A": "CD8A",
            "ImmuneScore": "ESTIMATE_Immune_score",
            "GEP18": "GEP18",
            "IFN": "Wolf_IFNgamma",
        },
        "os": None,
        "notes": "Wolf_IFNgamma n=444 (58 NA). MHC-I cassette and OS not on this Xena table.",
    }

    lusc = _read("TCGA_LUSC_xena_sample_table.tsv")
    out["TCGA-LUSC"] = {
        "df": lusc,
        "n_matrix": 502,
        "layer": "RNA",
        "anchors": {"CLDN4": "CLDN4", "TACSTD2": "TACSTD2"},
        "endpoints": {
            "CD8A": "CD8",
            "ImmuneScore": "ESTIMATE_ImmuneScore",
            "GEP18": "GEP18",
        },
        "os": None,
        "notes": "IFN, MHC-I, OS not on this Xena table. CD8 = published CD8 score column.",
    }

    cpl = _read("CPTAC_LUAD_sample_table.tsv")
    out["CPTAC-LUAD-protein"] = {
        "df": cpl,
        "n_matrix": 110,
        "layer": "protein",
        "anchors": {"CLDN4": "CLDN4_protein", "TACSTD2": "TACSTD2_protein"},
        "endpoints": {
            "CD8A": "CD8_Tcell_RNA",
            "ImmuneScore": "ESTIMATE_ImmuneScore",
            "GEP18": "GEP_Tcell_inflamed_RNA",
            "IFN": "IFNG_6gene_RNA",
        },
        "os": {"time": "OS_days", "event": "OS_event"},
        "notes": "Protein anchors. CLDN4 protein n=79. CD8A = CD8_Tcell_RNA (no CD8A gene column). MHC-I absent.",
    }

    lscc_s = _read("CPTAC_LSCC_sample_scores.tsv")
    lscc_f = _read("CPTAC_LSCC_sample_level_features.tsv")
    lscc = lscc_s.merge(
        lscc_f.rename(columns={"Unnamed: 0": "feat_id"})[
            ["feat_id", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "OS_days", "OS_event"]
        ],
        left_on="sample",
        right_on="feat_id",
        how="inner",
    )
    out["CPTAC-LSCC-protein"] = {
        "df": lscc,
        "n_matrix": 108,
        "layer": "protein",
        "anchors": {"CLDN4": "CLDN4_protein", "TACSTD2": "TACSTD2_protein"},
        "endpoints": {
            "CD8A": "CD8A_RNA",
            "ImmuneScore": "ImmuneScore",
            "GEP18": "GEP18_RNA",
            "IFN": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        },
        "os": {"time": "OS_days", "event": "OS_event"},
        "notes": "Protein anchors. CLDN4 protein n=78. MHC-I absent.",
    }

    out["GSE68465"] = {
        "df": _read("GSE68465_samples.tsv"),
        "n_matrix": 443,
        "layer": "RNA",
        "anchors": {"CLDN4": "CLDN4", "TACSTD2": "TACSTD2"},
        "endpoints": {
            "CD8A": "CD8A",
            "ImmuneScore": "ESTIMATE_ImmuneScore",
            "GEP18": "GEP18",
        },
        "os": None,
        "notes": "Director's Challenge U133A. IFN, MHC-I, OS not on this table. GEP18 is 15/18.",
    }

    out["GSE31210"] = {
        "df": _read("GSE31210_samples.tsv"),
        "n_matrix": 226,
        "layer": "RNA",
        "anchors": {"CLDN4": "CLDN4", "TACSTD2": "TACSTD2"},
        "endpoints": {
            "CD8A": "CD8A",
            "ImmuneScore": "ESTIMATE_ImmuneScore",
            "GEP18": "GEP18",
        },
        "os": None,
        "notes": "Okayama/Kohno LUAD. IFN, MHC-I, OS not on this table.",
    }

    onco = _read("OncoSG_sample_table.tsv")
    out["OncoSG"] = {
        "df": onco,
        "n_matrix": 169,
        "layer": "RNA",
        "anchors": {"TACSTD2": "TACSTD2_z"},  # CLDN4 not on public z-score table
        "endpoints": {
            "CD8A": "CD8",
            "ImmuneScore": "immune_tcell_effector",
            "GEP18": "GEP18",
            "IFN": "IMSIG_INTERFERON",
        },
        "os": None,
        "notes": "Public cBioPortal z-scores only; CLDN4 absent. ImmuneScore = A1 8-gene (ESTIMATE not computable). GEP18 17/18.",
    }
    return out


def extreme_quartile(anchor: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Equal-count Q1 / Q4 via rank qcut (ties broken by first appearance)."""
    a = pd.to_numeric(anchor, errors="coerce")
    valid = a.dropna()
    if valid.nunique() < 4:
        raise ValueError("fewer than 4 unique anchor values")
    ranks = valid.rank(method="first")
    q = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    q1 = pd.Series(False, index=a.index)
    q4 = pd.Series(False, index=a.index)
    q1.loc[q.index[q == "Q1"]] = True
    q4.loc[q.index[q == "Q4"]] = True
    return q1, q4


def rank_biserial(x_q4: np.ndarray, x_q1: np.ndarray) -> float:
    """Positive if Q4 > Q1. r_rb = 2U/(n4 n1) - 1 with U = MWU(Q4, Q1)."""
    n4, n1 = len(x_q4), len(x_q1)
    if n4 == 0 or n1 == 0:
        return np.nan
    u = stats.mannwhitneyu(x_q4, x_q1, alternative="two-sided").statistic
    return float(2.0 * u / (n4 * n1) - 1.0)


def bootstrap_rb(x_q4: np.ndarray, x_q1: np.ndarray, n_boot: int = N_BOOT) -> tuple[float, float]:
    if len(x_q4) < 3 or len(x_q1) < 3:
        return np.nan, np.nan
    vals = []
    for _ in range(n_boot):
        a = RNG.choice(x_q4, size=len(x_q4), replace=True)
        b = RNG.choice(x_q1, size=len(x_q1), replace=True)
        vals.append(rank_biserial(a, b))
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def logrank_hr(t1, e1, t4, e4) -> dict:
    """Two-sample log-rank (Q4 vs Q1) and Mantel–Haenszel HR (Q4 / Q1)."""
    t1 = np.asarray(t1, float)
    e1 = np.asarray(e1, float)
    t4 = np.asarray(t4, float)
    e4 = np.asarray(e4, float)
    times = np.unique(np.concatenate([t1[e1 == 1], t4[e4 == 1]]))
    o_minus_e = 0.0
    v = 0.0
    o4 = 0.0
    e4_sum = 0.0
    o1 = 0.0
    e1_sum = 0.0
    for t in times:
        n1 = np.sum(t1 >= t)
        n4 = np.sum(t4 >= t)
        n = n1 + n4
        if n < 2:
            continue
        d1 = np.sum((t1 == t) & (e1 == 1))
        d4 = np.sum((t4 == t) & (e4 == 1))
        d = d1 + d4
        if d == 0:
            continue
        e1_t = d * n1 / n
        e4_t = d * n4 / n
        o_minus_e += d4 - e4_t
        v += (d * (n - d) / (n - 1.0)) * (n1 / n) * (n4 / n)
        o4 += d4
        e4_sum += e4_t
        o1 += d1
        e1_sum += e1_t
    if v <= 0:
        return {"logrank_p": np.nan, "hr": np.nan, "hr_lo": np.nan, "hr_hi": np.nan, "z": np.nan}
    z = o_minus_e / np.sqrt(v)
    p = 2.0 * stats.norm.sf(abs(z))
    hr = np.exp(o_minus_e / v)
    se = 1.0 / np.sqrt(v)
    hr_lo = np.exp((o_minus_e / v) - 1.96 * se)
    hr_hi = np.exp((o_minus_e / v) + 1.96 * se)
    return {
        "logrank_p": float(p),
        "hr": float(hr),
        "hr_lo": float(hr_lo),
        "hr_hi": float(hr_hi),
        "z": float(z),
        "O4": float(o4),
        "E4": float(e4_sum),
        "O1": float(o1),
        "E1": float(e1_sum),
    }


def median_surv(time: np.ndarray, event: np.ndarray) -> float:
    """Kaplan–Meier median; nan if not reached."""
    order = np.argsort(time)
    t = time[order]
    e = event[order]
    n = len(t)
    if n == 0:
        return np.nan
    at_risk = n
    s = 1.0
    last_t = np.nan
    for i in range(n):
        if i == 0 or t[i] != t[i - 1]:
            d = 0
            j = i
            while j < n and t[j] == t[i]:
                d += int(e[j] == 1)
                j += 1
            if at_risk > 0 and d > 0:
                s *= 1.0 - d / at_risk
            last_t = t[i]
            if s <= 0.5:
                return float(t[i])
        at_risk -= 1
    return np.nan


def recut_continuous(df: pd.DataFrame, anchor_col: str, end_col: str) -> dict | None:
    a = pd.to_numeric(df[anchor_col], errors="coerce")
    y = pd.to_numeric(df[end_col], errors="coerce")
    ok = a.notna() & y.notna()
    sub = df.loc[ok].copy()
    if sub.shape[0] < 8:
        return None
    q1, q4 = extreme_quartile(pd.to_numeric(sub[anchor_col], errors="coerce"))
    y1 = pd.to_numeric(sub.loc[q1, end_col], errors="coerce").to_numpy(float)
    y4 = pd.to_numeric(sub.loc[q4, end_col], errors="coerce").to_numpy(float)
    y1 = y1[np.isfinite(y1)]
    y4 = y4[np.isfinite(y4)]
    if len(y1) < 3 or len(y4) < 3:
        return None
    mwu = stats.mannwhitneyu(y4, y1, alternative="two-sided")
    rb = rank_biserial(y4, y1)
    rb_lo, rb_hi = bootstrap_rb(y4, y1)
    return {
        "n_complete": int(ok.sum()),
        "n_q1": int(len(y1)),
        "n_q4": int(len(y4)),
        "median_q1": float(np.median(y1)),
        "median_q4": float(np.median(y4)),
        "delta_median": float(np.median(y4) - np.median(y1)),
        "rank_biserial": rb,
        "rb_lo": rb_lo,
        "rb_hi": rb_hi,
        "mwu_u": float(mwu.statistic),
        "p": float(mwu.pvalue),
        "y1": y1,
        "y4": y4,
    }


def recut_os(df: pd.DataFrame, anchor_col: str, time_col: str, event_col: str) -> dict | None:
    a = pd.to_numeric(df[anchor_col], errors="coerce")
    t = pd.to_numeric(df[time_col], errors="coerce")
    e = pd.to_numeric(df[event_col], errors="coerce")
    ok = a.notna() & t.notna() & e.notna() & (t > 0)
    sub = df.loc[ok].copy()
    if sub.shape[0] < 16:
        return None
    q1, q4 = extreme_quartile(pd.to_numeric(sub[anchor_col], errors="coerce"))
    t1 = pd.to_numeric(sub.loc[q1, time_col], errors="coerce").to_numpy(float)
    e1 = pd.to_numeric(sub.loc[q1, event_col], errors="coerce").to_numpy(float)
    t4 = pd.to_numeric(sub.loc[q4, time_col], errors="coerce").to_numpy(float)
    e4 = pd.to_numeric(sub.loc[q4, event_col], errors="coerce").to_numpy(float)
    if len(t1) < 5 or len(t4) < 5:
        return None
    lr = logrank_hr(t1, e1, t4, e4)
    return {
        "n_complete": int(ok.sum()),
        "n_q1": int(len(t1)),
        "n_q4": int(len(t4)),
        "events_q1": int(e1.sum()),
        "events_q4": int(e4.sum()),
        "median_os_q1": median_surv(t1, e1),
        "median_os_q4": median_surv(t4, e4),
        "delta_median": (
            float(median_surv(t4, e4) - median_surv(t1, e1))
            if np.isfinite(median_surv(t4, e4)) and np.isfinite(median_surv(t1, e1))
            else np.nan
        ),
        "hr": lr["hr"],
        "hr_lo": lr["hr_lo"],
        "hr_hi": lr["hr_hi"],
        "p": lr["logrank_p"],
        "rank_biserial": np.nan,
        "rb_lo": np.nan,
        "rb_hi": np.nan,
    }


def run() -> tuple[pd.DataFrame, dict, dict]:
    cohorts = load_cohorts()
    rows = []
    violin_data: dict[tuple[str, str, str], dict] = {}
    inventory = []

    for cohort, spec in cohorts.items():
        df = spec["df"]
        inventory.append(
            {
                "cohort": cohort,
                "n_matrix": spec["n_matrix"],
                "n_table": int(len(df)),
                "layer": spec["layer"],
                "anchors": ",".join(spec["anchors"]),
                "endpoints_present": ",".join(list(spec["endpoints"]) + (["OS"] if spec["os"] else [])),
                "endpoints_absent": ",".join(
                    e
                    for e in ENDPOINTS
                    if e not in spec["endpoints"] and not (e == "OS" and spec["os"])
                ),
                "notes": spec["notes"],
            }
        )
        for gene, acol in spec["anchors"].items():
            for ep, ecol in spec["endpoints"].items():
                r = recut_continuous(df, acol, ecol)
                if r is None:
                    continue
                violin_data[(cohort, gene, ep)] = r
                rows.append(
                    {
                        "cohort": cohort,
                        "n_matrix": spec["n_matrix"],
                        "layer": spec["layer"],
                        "gene": gene,
                        "endpoint": ep,
                        "endpoint_column": ecol,
                        "n_complete": r["n_complete"],
                        "n_q1": r["n_q1"],
                        "n_q4": r["n_q4"],
                        "median_q1": r["median_q1"],
                        "median_q4": r["median_q4"],
                        "delta_median": r["delta_median"],
                        "rank_biserial": r["rank_biserial"],
                        "rb_lo": r["rb_lo"],
                        "rb_hi": r["rb_hi"],
                        "p": r["p"],
                        "hr": np.nan,
                        "hr_lo": np.nan,
                        "hr_hi": np.nan,
                        "events_q1": np.nan,
                        "events_q4": np.nan,
                    }
                )
            if spec["os"]:
                r = recut_os(df, acol, spec["os"]["time"], spec["os"]["event"])
                if r is None:
                    continue
                rows.append(
                    {
                        "cohort": cohort,
                        "n_matrix": spec["n_matrix"],
                        "layer": spec["layer"],
                        "gene": gene,
                        "endpoint": "OS",
                        "endpoint_column": f"{spec['os']['time']}+{spec['os']['event']}",
                        "n_complete": r["n_complete"],
                        "n_q1": r["n_q1"],
                        "n_q4": r["n_q4"],
                        "median_q1": r["median_os_q1"],
                        "median_q4": r["median_os_q4"],
                        "delta_median": r["delta_median"],
                        "rank_biserial": np.nan,
                        "rb_lo": np.nan,
                        "rb_hi": np.nan,
                        "p": r["p"],
                        "hr": r["hr"],
                        "hr_lo": r["hr_lo"],
                        "hr_hi": r["hr_hi"],
                        "events_q1": r["events_q1"],
                        "events_q4": r["events_q4"],
                    }
                )

    res = pd.DataFrame(rows)
    inv = pd.DataFrame(inventory)
    return res, violin_data, inv


def _fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3f}"


def plot_forest(res: pd.DataFrame, gene: str, path: Path) -> None:
    sub = res[(res["gene"] == gene) & (res["endpoint"] != "OS")].copy()
    if sub.empty:
        return
    ep_order = [e for e in ["CD8A", "ImmuneScore", "GEP18", "IFN", "MHC-I"] if e in set(sub["endpoint"])]
    cohort_order = [
        "GSE218989",
        "GSE285029",
        "TCGA-LUAD",
        "TCGA-LUSC",
        "CPTAC-LUAD-protein",
        "CPTAC-LSCC-protein",
        "GSE68465",
        "GSE31210",
        "OncoSG",
    ]
    rows = []
    for ep in ep_order:
        for c in cohort_order:
            hit = sub[(sub["endpoint"] == ep) & (sub["cohort"] == c)]
            if len(hit) == 1:
                rows.append(hit.iloc[0])
    if not rows:
        return
    plot_df = pd.DataFrame(rows).iloc[::-1]  # top = first
    n = len(plot_df)
    fig_h = max(6.5, 0.38 * n + 1.8)
    fig, ax = plt.subplots(figsize=(9.2, fig_h))
    y = np.arange(n)
    x = plot_df["rank_biserial"].to_numpy(float)
    lo = plot_df["rb_lo"].to_numpy(float)
    hi = plot_df["rb_hi"].to_numpy(float)
    colors = ["#b2182b" if v > 0 else "#2166ac" for v in x]
    ax.axvline(0, color="#444444", lw=0.8)
    for i in range(n):
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=colors[i], lw=1.4, solid_capstyle="round")
        ax.plot(x[i], y[i], "o", color=colors[i], ms=5.5, zorder=3)
    labels = [
        f"{r.cohort}  {r.endpoint}   n={int(r.n_q1)}/{int(r.n_q4)}  p={_fmt_p(r.p)}"
        for r in plot_df.itertuples()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"rank-biserial δ  ({gene} Q4 − Q1; + = higher in Q4)")
    ax.set_title(f"{gene} extreme quartile (Q4 vs Q1) — existing public matrices")
    ax.set_xlim(-1.05, 1.05)
    ax.grid(axis="x", ls=":", alpha=0.4)
    for ep in ep_order:
        idx = [i for i, r in enumerate(plot_df.itertuples()) if r.endpoint == ep]
        if idx:
            ax.axhline(max(idx) + 0.5, color="#dddddd", lw=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_os_forest(res: pd.DataFrame, path: Path) -> None:
    sub = res[res["endpoint"] == "OS"].copy()
    if sub.empty:
        return
    sub = sub.sort_values(["gene", "cohort"])
    n = len(sub)
    fig, ax = plt.subplots(figsize=(8.6, max(3.2, 0.45 * n + 1.4)))
    y = np.arange(n)
    hr = sub["hr"].to_numpy(float)
    lo = sub["hr_lo"].to_numpy(float)
    hi = sub["hr_hi"].to_numpy(float)
    colors = ["#b2182b" if v > 1 else "#2166ac" for v in hr]
    ax.axvline(1, color="#444444", lw=0.8)
    for i in range(n):
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=colors[i], lw=1.5)
        ax.plot(hr[i], y[i], "o", color=colors[i], ms=6)
    labels = [
        f"{r.gene}  {r.cohort}  n={int(r.n_q1)}/{int(r.n_q4)}  ev={int(r.events_q1)}/{int(r.events_q4)}  p={_fmt_p(r.p)}"
        for r in sub.itertuples()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("OS HR (Q4 vs Q1)")
    ax.set_title("OS Q4 vs Q1 — only cohorts with survival on the harvested table")
    ax.set_xscale("log")
    ax.grid(axis="x", ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_violin_grid(
    violin_data: dict,
    gene: str,
    endpoint: str,
    path: Path,
    title: str,
) -> None:
    keys = [(c, gene, endpoint) for (c, g, e) in violin_data if g == gene and e == endpoint]
    cohort_order = [
        "GSE218989",
        "GSE285029",
        "TCGA-LUAD",
        "TCGA-LUSC",
        "CPTAC-LUAD-protein",
        "CPTAC-LSCC-protein",
        "GSE68465",
        "GSE31210",
        "OncoSG",
    ]
    keys = [k for c in cohort_order for k in keys if k[0] == c]
    if not keys:
        return
    n = len(keys)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10.6, 3.15 * nrows), squeeze=False)
    for i, key in enumerate(keys):
        ax = axes[i // ncols][i % ncols]
        r = violin_data[key]
        data = [r["y1"], r["y4"]]
        parts = ax.violinplot(data, positions=[1, 2], showmeans=False, showmedians=True, showextrema=False)
        for j, b in enumerate(parts["bodies"]):
            b.set_facecolor("#6baed6" if j == 0 else "#fc8d59")
            b.set_alpha(0.85)
            b.set_edgecolor("#333333")
        parts["cmedians"].set_color("#222222")
        # jittered points
        for xpos, arr, col in ((1, r["y1"], "#2171b5"), (2, r["y4"], "#d94801")):
            jitter = RNG.normal(0, 0.045, size=len(arr))
            ax.scatter(np.full(len(arr), xpos) + jitter, arr, s=7, c=col, alpha=0.35, linewidths=0, zorder=2)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"Q1\nn={r['n_q1']}", f"Q4\nn={r['n_q4']}"], fontsize=8)
        ax.set_title(f"{key[0]}\nδ={r['delta_median']:+.3g}  p={_fmt_p(r['p'])}", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(title, fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_key_violins(violin_data: dict, path: Path) -> None:
    """One figure: CLDN4 Q1 vs Q4 CD8A across every cohort that has it."""
    plot_violin_grid(
        violin_data,
        "CLDN4",
        "CD8A",
        path,
        "CLDN4 Q4 vs Q1 — CD8A (existing matrices; no continuous cloud)",
    )


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    res, violin_data, inv = run()
    # drop raw arrays from the written table
    res.to_csv(TAB / "q4_vs_q1.tsv", sep="\t", index=False, float_format="%.6g")
    inv.to_csv(TAB / "cohort_inventory.tsv", sep="\t", index=False)

    # compact display table
    show = res.copy()
    show["p"] = show["p"].map(lambda x: f"{x:.3e}" if pd.notna(x) else "")
    show.to_csv(TAB / "q4_vs_q1_display.tsv", sep="\t", index=False)

    plot_forest(res, "CLDN4", FIG / "fig1_forest_cldn4.png")
    plot_forest(res, "TACSTD2", FIG / "fig2_forest_tacstd2.png")
    plot_os_forest(res, FIG / "fig3_forest_os.png")
    plot_key_violins(violin_data, FIG / "fig4_violins_cldn4_cd8a.png")
    plot_violin_grid(
        violin_data,
        "CLDN4",
        "ImmuneScore",
        FIG / "fig5_violins_cldn4_immunescore.png",
        "CLDN4 Q4 vs Q1 — ImmuneScore (existing matrices)",
    )
    plot_violin_grid(
        violin_data,
        "CLDN4",
        "GEP18",
        FIG / "fig6_violins_cldn4_gep18.png",
        "CLDN4 Q4 vs Q1 — GEP18 (existing matrices)",
    )
    plot_violin_grid(
        violin_data,
        "TACSTD2",
        "CD8A",
        FIG / "fig7_violins_tacstd2_cd8a.png",
        "TACSTD2 Q4 vs Q1 — CD8A (existing matrices; no continuous cloud)",
    )
    plot_violin_grid(
        violin_data,
        "TACSTD2",
        "ImmuneScore",
        FIG / "fig8_violins_tacstd2_immunescore.png",
        "TACSTD2 Q4 vs Q1 — ImmuneScore (existing matrices)",
    )
    plot_violin_grid(
        violin_data,
        "CLDN4",
        "IFN",
        FIG / "fig9_violins_cldn4_ifn.png",
        "CLDN4 Q4 vs Q1 — IFN (only cohorts with an IFN column)",
    )
    plot_violin_grid(
        violin_data,
        "CLDN4",
        "MHC-I",
        FIG / "fig10_violins_cldn4_mhci.png",
        "CLDN4 Q4 vs Q1 — MHC-I (only GSE218989 and GSE285029 have this column)",
    )

    summary = {
        "n_tests": int(len(res)),
        "n_cldn4": int((res["gene"] == "CLDN4").sum()),
        "n_tacstd2": int((res["gene"] == "TACSTD2").sum()),
        "cohorts": inv.to_dict(orient="records"),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2))
    print(res.to_string(index=False))
    print("\nWrote", TAB / "q4_vs_q1.tsv")
    print("Figures:", sorted(p.name for p in FIG.glob("*.png")))


if __name__ == "__main__":
    main()

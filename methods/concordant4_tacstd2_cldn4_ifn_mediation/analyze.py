#!/usr/bin/env python3
"""Concordant-4 malignant pseudobulk: TACSTD2, CLDN4, IFN/MHC.

Question
--------
Does TACSTD2-high IFN/MHC downregulation attenuate after residualizing CLDN4,
or inside the CLDN4-low half? Is the patient-level path TACSTD2 → CLDN4 → IFN
supported on the same units?

Pre-specified before the NES and mediation numbers were read back:

* Cohorts are the locked four only: GSE123902, GSE131907, GSE205335, GSE189357.
* Unit is the patient / donor / sample. Expression n is the UMI-sum (P4001 out).
* Expression is log2(TMM-CPM+1). One TMM on genes with count >= 10 in >= 3 units.
* Primary scale is the pseudobulk. TACSTD2 and CLDN4 are the malignant
  log2(TMM-CPM+1) of those two genes. Percent-positive is a second scale.
* TACSTD2-high is the within-cohort Q4 vs Q1 contrast. The continuous
  cohort-adjusted slope is co-primary because residualizing a continuous
  CLDN4 covariate is the mediation model.
* Residualized rank: joint OLS t for TACSTD2 after CLDN4 and cohort dummies.
* CLDN4-low stratum: within-cohort CLDN4 rank at or below the median rank.
  The stratum model does not also residualize CLDN4.
* Positive NES = enriched in TACSTD2-high. GSEA is preranked, weight p=1.
* Mediation is cohort-adjusted product of coefficients on the patient-level
  mean of each gene set. NES attenuation is the change in that NES, not a
  relabeling of a×b.
* Observational. Not a knockdown and not a revision of the locked T/NK rho.

Calibration: CLDN4 %pos Q4 vs Q1 IFN-family score on the Q1/Q4 subset TMM,
the PR #503 contrast (IFN logFC about -0.584, 18 vs 16).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gsea import gsea_on_scores

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "results" / "tables"
FIGS = HERE / "results" / "figures"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
TERMS = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "IFN_HALLMARK_UNION",
    "MHC_I_APM",
]
PRIMARY_TERMS = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "MHC_I_APM",
]
NPERM = 1000
NPERM_BOOT = 200
B_NES = 300
B_MED = 4000
SEED = 42

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 140,
    }
)


def bh(p: np.ndarray) -> np.ndarray:
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
        if int(keep.sum()) < 50:
            factors[col] = 1.0
            continue
        m = np.log2((obs[keep] / obs_lib) / (ref_c[keep] / ref_lib))
        a = 0.5 * np.log2((obs[keep] / obs_lib) * (ref_c[keep] / ref_lib))
        w = (obs_lib - obs[keep]) / (obs_lib * obs[keep]) + (ref_lib - ref_c[keep]) / (
            ref_lib * ref_c[keep]
        )
        ok = np.isfinite(m) & np.isfinite(a) & np.isfinite(w) & (w > 0)
        m, a, w = np.asarray(m)[ok], np.asarray(a)[ok], np.asarray(w)[ok]
        if len(m) < 50:
            factors[col] = 1.0
            continue
        lo_m, hi_m = np.quantile(m, [0.30, 0.70])
        lo_a, hi_a = np.quantile(a, [0.05, 0.95])
        trim = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if int(trim.sum()) < 20:
            factors[col] = 1.0
            continue
        factors[col] = 2 ** float(np.average(m[trim], weights=1.0 / w[trim]))
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


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        # A bootstrap draw can tie enough values that four labeled bins do not exist.
        return pd.Series(["NA"] * len(s), index=s.index)
    return pd.Series(qs.astype(str), index=s.index)


def cohort_dummies(cohort: np.ndarray) -> np.ndarray:
    cols = []
    for c in COHORTS:
        if c == REF:
            continue
        cols.append((cohort == c).astype(float))
    return np.column_stack(cols)


def ols_t(Y: np.ndarray, X: np.ndarray, j: int) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Y is n x g. Returns t, beta, xtx_inv_jj, df."""
    n, p = X.shape
    xtx = X.T @ X
    xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ (X.T @ Y)
    resid = Y - X @ beta
    rank = int(np.linalg.matrix_rank(X, tol=1e-8))
    df = int(n - rank)
    if df < 1:
        nan = np.full(Y.shape[1], np.nan)
        return nan, nan, np.nan, df
    sigma2 = np.sum(resid ** 2, axis=0) / df
    se = np.sqrt(np.maximum(sigma2 * float(xtx_inv[j, j]), 0.0))
    est = beta[j]
    t = np.zeros(Y.shape[1], dtype=float)
    np.divide(est, se, out=t, where=se > 0)
    return t, est, float(xtx_inv[j, j]), df


def ols_coef(y: np.ndarray, X: np.ndarray, j: int) -> dict:
    n, p = X.shape
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ (X.T @ y)
    resid = y - X @ beta
    rank = int(np.linalg.matrix_rank(X, tol=1e-8))
    df = int(n - rank)
    sigma2 = float(np.sum(resid ** 2) / df) if df >= 1 else np.nan
    se = float(np.sqrt(max(sigma2 * float(xtx_inv[j, j]), 0.0))) if df >= 1 else np.nan
    est = float(beta[j])
    if se > 0 and df >= 1:
        t = est / se
        p = float(2 * stats.t.sf(abs(t), df))
    else:
        t, p = np.nan, np.nan
    return {"beta": est, "se": se, "t": t, "p": p, "df": df, "n": n}


def design_matrix(x: np.ndarray, m: np.ndarray | None, cohort: np.ndarray) -> tuple[np.ndarray, int]:
    parts = [np.ones(len(x)), np.asarray(x, float)]
    j = 1
    if m is not None:
        parts.append(np.asarray(m, float))
    parts.append(cohort_dummies(cohort))
    return np.column_stack(parts), j


def spearman(a: np.ndarray, b: np.ndarray) -> tuple[float, float, int]:
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < 5:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(a[m], b[m])
    return float(rho), float(p), n


def residualize(y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    beta = np.linalg.pinv(Z.T @ Z) @ (Z.T @ y)
    return y - Z @ beta


def fisher_dl(rhos: np.ndarray, ns: np.ndarray, k_partial: int = 0) -> dict:
    rhos = np.asarray(rhos, float)
    ns = np.asarray(ns, float)
    ok = np.isfinite(rhos) & (ns - 3 - k_partial > 1)
    rhos, ns = rhos[ok], ns[ok]
    if len(rhos) < 2:
        return {"k": int(len(rhos)), "rho": np.nan, "lo": np.nan, "hi": np.nan, "p": np.nan, "I2": np.nan}
    z = np.arctanh(np.clip(rhos, -0.999, 0.999))
    v = 1.0 / (ns - 3 - k_partial)
    w = 1.0 / v
    zbar = np.sum(w * z) / np.sum(w)
    Q = float(np.sum(w * (z - zbar) ** 2))
    df = len(z) - 1
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wstar = 1.0 / (v + tau2)
    zhat = float(np.sum(wstar * z) / np.sum(wstar))
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zhat / se))) if se > 0 else np.nan
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return {
        "k": int(len(rhos)),
        "rho": float(np.tanh(zhat)),
        "lo": float(np.tanh(zhat - 1.96 * se)),
        "hi": float(np.tanh(zhat + 1.96 * se)),
        "p": p,
        "I2": float(I2),
        "tau2": tau2,
    }


def load_gene_sets() -> dict[str, set[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    ifna = {g.upper() for g in sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]}
    ifng = {g.upper() for g in sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]}
    mhc = {g.upper() for g in sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]}
    for block in (ifna, ifng, mhc):
        block.discard("CLDN4")
        block.discard("TACSTD2")
    return {
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": ifng,
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": ifna,
        "IFN_HALLMARK_UNION": ifna | ifng,
        "MHC_I_APM": mhc,
    }


def load_locked_cldn4_meta() -> pd.DataFrame:
    """PR #503 phenotype, including P4001, with within-cohort CLDN4 % quartiles."""
    rows = []
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"]
    pct = el["mal_CLDN4_pct"].astype(float).to_numpy()
    if np.nanmax(pct) <= 1.5:
        pct = pct * 100.0
    q = assign_quartiles(pd.Series(pct, index=el["patient"].astype(str)))
    for i, pid in enumerate(el["patient"].astype(str)):
        rows.append({"patient": pid, "cohort": "GSE123902", "cldn4_pct_locked": float(pct[i]), "cldn4_q_locked": str(q.loc[pid])})

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    pct = mal["mal_CLDN4_pct"].astype(float)
    if float(np.nanmax(pct)) <= 1.5:
        pct = pct * 100.0
    q = assign_quartiles(pd.Series(pct.to_numpy(), index=mal["sample"].astype(str)))
    for pid, val in zip(mal["sample"].astype(str), pct.astype(float)):
        rows.append({"patient": pid, "cohort": "GSE131907", "cldn4_pct_locked": float(val), "cldn4_q_locked": str(q.loc[pid])})

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    pct = d["mal_CLDN4_pct_pos"].astype(float)
    if float(np.nanmax(pct)) <= 1.5:
        pct = pct * 100.0
    q = assign_quartiles(pd.Series(pct.to_numpy(), index=d["patient"].astype(str)))
    for pid, val in zip(d["patient"].astype(str), pct.astype(float)):
        rows.append({"patient": pid, "cohort": "GSE205335", "cldn4_pct_locked": float(val), "cldn4_q_locked": str(q.loc[pid])})

    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"]
    pct = el["mal_CLDN4_pct"].astype(float).to_numpy()
    if np.nanmax(pct) <= 1.5:
        pct = pct * 100.0
    q = assign_quartiles(pd.Series(pct, index=el["patient"].astype(str)))
    for i, pid in enumerate(el["patient"].astype(str)):
        rows.append({"patient": pid, "cohort": "GSE189357", "cldn4_pct_locked": float(pct[i]), "cldn4_q_locked": str(q.loc[pid])})
    return pd.DataFrame(rows)


def load_pr638() -> pd.DataFrame:
    sc = pd.read_csv(DATA / "patient_scores_pr638.tsv", sep="\t")
    sc = sc.rename(columns={"dataset": "cohort", "unit_id": "patient"})
    keep = ["cohort", "patient", "n_malignant", "pct_TACSTD2", "pct_CLDN4", "mean_TACSTD2", "mean_CLDN4"]
    return sc[keep].copy()


def load_intersect_counts() -> tuple[pd.DataFrame, dict[str, str]]:
    """Genes measured in every cohort, samples side by side.

    Filling genes that exist in only some cohorts changes the TMM factors.
    The intersection is what reproduces the PR #503 IFN family score.
    """
    parts = {}
    cohort_of: dict[str, str] = {}
    for cohort in COHORTS:
        cts = read_counts(DATA / f"{cohort}_malignant_counts.tsv.gz")
        parts[cohort] = cts
        for col in cts.columns:
            cohort_of[str(col)] = cohort
    genes = sorted(set.intersection(*[set(p.index) for p in parts.values()]))
    combined = pd.concat([p.reindex(genes).fillna(0.0) for p in parts.values()], axis=1)
    combined = combined.loc[:, ~combined.columns.duplicated()]
    return combined, cohort_of


def assemble() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    counts, cohort_of = load_intersect_counts()
    counts = filter_genes(counts)
    factors = tmm_norm_factors(counts)
    logc = log_cpm(counts, factors)

    locked = load_locked_cldn4_meta()
    pr = load_pr638()
    rows = []
    for col in logc.columns:
        rows.append({"col": str(col), "cohort": cohort_of[str(col)], "patient": str(col)})
    ph = pd.DataFrame(rows)
    ph = ph.merge(locked, on=["cohort", "patient"], how="left")
    ph = ph.merge(pr, on=["cohort", "patient"], how="left")
    ph["tacstd2_pb"] = [float(logc.loc["TACSTD2", c]) if "TACSTD2" in logc.index else np.nan for c in ph["col"]]
    ph["cldn4_pb"] = [float(logc.loc["CLDN4", c]) if "CLDN4" in logc.index else np.nan for c in ph["col"]]
    ph["tacstd2_q"] = ""
    ph["cldn4_q_pb"] = ""
    for cohort, idx in ph.groupby("cohort").groups.items():
        ph.loc[idx, "tacstd2_q"] = assign_quartiles(ph.loc[idx, "tacstd2_pb"])
        ph.loc[idx, "cldn4_q_pb"] = assign_quartiles(ph.loc[idx, "cldn4_pb"])
        ph.loc[idx, "tacstd2_pct_q"] = assign_quartiles(ph.loc[idx, "pct_TACSTD2"])
    # lower half by within-cohort rank (median rank included in the low half)
    ph["cldn4_low"] = False
    ph["cldn4_pct_low"] = False
    for cohort, idx in ph.groupby("cohort").groups.items():
        r = ph.loc[idx, "cldn4_pb"].rank(method="average")
        ph.loc[idx, "cldn4_low"] = r <= r.median()
        r2 = ph.loc[idx, "pct_CLDN4"].rank(method="average")
        ph.loc[idx, "cldn4_pct_low"] = r2 <= r2.median()
    sets = load_gene_sets()
    return ph.reset_index(drop=True), logc, sets


def score_matrix(logc: pd.DataFrame, ph: pd.DataFrame, sets: dict[str, set[str]]) -> pd.DataFrame:
    out = ph[["col", "cohort", "patient"]].copy()
    for name, genes in sets.items():
        present = [g for g in genes if g in logc.index]
        out[name] = logc.loc[present, ph["col"]].mean(axis=0).to_numpy()
        out[f"{name}__n_genes"] = len(present)
    return out


def Y_of(logc: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    sub = logc.loc[:, cols]
    genes = sub.index.to_numpy()
    Y = sub.to_numpy(dtype=float).T  # n x g
    return Y, genes


def run_gsea_contrast(
    Y: np.ndarray,
    genes: np.ndarray,
    X: np.ndarray,
    j: int,
    sets: dict[str, set[str]],
    nperm: int,
    seed: int,
    with_lead: bool,
) -> list[dict]:
    t, beta, _, df = ols_t(Y, X, j)
    rng = np.random.default_rng(seed)
    rows = gsea_on_scores(t, genes, sets, nperm=nperm, rng=rng, with_lead=with_lead)
    for r in rows:
        r["df"] = df
        r["n"] = int(X.shape[0])
    return rows


def subset_idx(ph: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    return np.flatnonzero(mask)


def contrast_specs(ph: pd.DataFrame) -> list[dict]:
    """Each spec builds X on a row mask. coef j is always column 1 (TACSTD2)."""
    n = len(ph)
    cohort = ph["cohort"].to_numpy()
    specs = []

    def add(scale, contrast, x, m, mask, residualize: bool):
        specs.append(
            {
                "scale": scale,
                "contrast": contrast,
                "x": np.asarray(x, float),
                "m": None if m is None else np.asarray(m, float),
                "mask": np.asarray(mask, bool),
                "residualize": residualize,
                "cohort": cohort,
            }
        )

    q4 = ph["tacstd2_q"].isin(["Q1", "Q4"]).to_numpy()
    x_q = (ph["tacstd2_q"] == "Q4").astype(float).to_numpy()
    add("pseudobulk", "Q4_vs_Q1", x_q, ph["cldn4_pb"].to_numpy(), q4, False)
    add("pseudobulk", "Q4_vs_Q1_resid_CLDN4", x_q, ph["cldn4_pb"].to_numpy(), q4, True)
    add("pseudobulk", "continuous", ph["tacstd2_pb"].to_numpy(), ph["cldn4_pb"].to_numpy(), np.ones(n, bool), False)
    add("pseudobulk", "continuous_resid_CLDN4", ph["tacstd2_pb"].to_numpy(), ph["cldn4_pb"].to_numpy(), np.ones(n, bool), True)
    low = ph["cldn4_low"].to_numpy()
    add("pseudobulk", "CLDN4_low_continuous", ph["tacstd2_pb"].to_numpy(), ph["cldn4_pb"].to_numpy(), low, False)

    q4p = ph["tacstd2_pct_q"].isin(["Q1", "Q4"]).to_numpy()
    x_qp = (ph["tacstd2_pct_q"] == "Q4").astype(float).to_numpy()
    add("percent_positive", "Q4_vs_Q1", x_qp, ph["pct_CLDN4"].to_numpy(), q4p, False)
    add("percent_positive", "Q4_vs_Q1_resid_CLDN4", x_qp, ph["pct_CLDN4"].to_numpy(), q4p, True)
    add("percent_positive", "continuous", ph["pct_TACSTD2"].to_numpy(), ph["pct_CLDN4"].to_numpy(), np.ones(n, bool), False)
    add("percent_positive", "continuous_resid_CLDN4", ph["pct_TACSTD2"].to_numpy(), ph["pct_CLDN4"].to_numpy(), np.ones(n, bool), True)
    add("percent_positive", "CLDN4_low_continuous", ph["pct_TACSTD2"].to_numpy(), ph["pct_CLDN4"].to_numpy(), ph["cldn4_pct_low"].to_numpy(), False)
    return specs


def eval_spec(spec: dict, Y_full: np.ndarray, genes: np.ndarray, sets: dict, nperm: int, seed: int, with_lead: bool) -> list[dict]:
    mask = spec["mask"] & np.isfinite(spec["x"])
    if spec["residualize"]:
        mask = mask & np.isfinite(spec["m"])
    if spec["contrast"].startswith("Q4"):
        # Q1 and Q4 only; the indicator is already 0/1 on that mask
        pass
    idx = np.flatnonzero(mask)
    if idx.size < 12:
        return []
    x = spec["x"][idx]
    cohort = spec["cohort"][idx]
    m = spec["m"][idx] if spec["residualize"] else None
    if spec["contrast"].startswith("Q4") and (np.unique(x).size < 2):
        return []
    if np.unique(cohort).size < 2:
        return []
    X, j = design_matrix(x, m, cohort)
    rows = run_gsea_contrast(Y_full[idx], genes, X, j, sets, nperm, seed, with_lead)
    for r in rows:
        r["scale"] = spec["scale"]
        r["contrast"] = spec["contrast"]
        r["n_high"] = int(np.sum(x == 1)) if spec["contrast"].startswith("Q4") else np.nan
        r["n_low"] = int(np.sum(x == 0)) if spec["contrast"].startswith("Q4") else np.nan
    return rows


def module_models(ph: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cohort = ph["cohort"].to_numpy()
    for scale, xcol, mcol in [
        ("pseudobulk", "tacstd2_pb", "cldn4_pb"),
        ("percent_positive", "pct_TACSTD2", "pct_CLDN4"),
    ]:
        x = ph[xcol].to_numpy(float)
        m = ph[mcol].to_numpy(float)
        for term in TERMS:
            y = scores[term].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(m) & np.isfinite(y)
            Xc, _ = design_matrix(x[ok], None, cohort[ok])
            Xd, _ = design_matrix(x[ok], m[ok], cohort[ok])
            Xm, _ = design_matrix(m[ok], None, cohort[ok])
            # design_matrix puts the exposure in column 1. For M ~ X, exposure is X.
            c = ols_coef(y[ok], Xc, 1)
            direct = ols_coef(y[ok], Xd, 1)
            b = ols_coef(y[ok], Xd, 2)
            a = ols_coef(m[ok], Xc, 1)
            # VIF of X given M + cohort
            Xv, _ = design_matrix(m[ok], None, cohort[ok])
            rx = residualize(x[ok], Xv)
            # R^2 of X ~ M + cohort
            ss_tot = float(np.sum((x[ok] - x[ok].mean()) ** 2))
            ss_res = float(np.sum(rx ** 2))
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            vif = 1 / (1 - r2) if np.isfinite(r2) and r2 < 0.999 else np.nan
            rows.append(
                {
                    "scale": scale,
                    "outcome": term,
                    "n": int(ok.sum()),
                    "a": a["beta"],
                    "a_se": a["se"],
                    "a_p": a["p"],
                    "b": b["beta"],
                    "b_se": b["se"],
                    "b_p": b["p"],
                    "c": c["beta"],
                    "c_se": c["se"],
                    "c_p": c["p"],
                    "c_prime": direct["beta"],
                    "c_prime_se": direct["se"],
                    "c_prime_p": direct["p"],
                    "ab": a["beta"] * b["beta"],
                    "vif_tacstd2": vif,
                    "r2_tacstd2_on_cldn4": r2,
                }
            )
    return pd.DataFrame(rows)


def bootstrap_mediation(ph: pd.DataFrame, scores: pd.DataFrame, B: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    cohort = ph["cohort"].to_numpy()
    groups = {c: np.flatnonzero(cohort == c) for c in COHORTS}
    rows = []
    stores = {
        (scale, term): []
        for scale in ("pseudobulk", "percent_positive")
        for term in TERMS
    }
    xcols = {"pseudobulk": "tacstd2_pb", "percent_positive": "pct_TACSTD2"}
    mcols = {"pseudobulk": "cldn4_pb", "percent_positive": "pct_CLDN4"}
    for b in range(B):
        idx = np.concatenate([rng.choice(groups[c], size=len(groups[c]), replace=True) for c in COHORTS])
        sub_c = cohort[idx]
        for scale in ("pseudobulk", "percent_positive"):
            x = ph[xcols[scale]].to_numpy(float)[idx]
            m = ph[mcols[scale]].to_numpy(float)[idx]
            for term in TERMS:
                y = scores[term].to_numpy(float)[idx]
                ok = np.isfinite(x) & np.isfinite(m) & np.isfinite(y)
                if ok.sum() < 20 or np.unique(sub_c[ok]).size < 4:
                    continue
                Xc, _ = design_matrix(x[ok], None, sub_c[ok])
                Xd, _ = design_matrix(x[ok], m[ok], sub_c[ok])
                a = ols_coef(m[ok], Xc, 1)
                bcoef = ols_coef(y[ok], Xd, 2)
                c = ols_coef(y[ok], Xc, 1)
                cp = ols_coef(y[ok], Xd, 1)
                if not all(np.isfinite([a["beta"], bcoef["beta"], c["beta"], cp["beta"]])):
                    continue
                ab = a["beta"] * bcoef["beta"]
                prop = ab / c["beta"] if c["beta"] != 0 else np.nan
                stores[(scale, term)].append((ab, prop, c["beta"], cp["beta"], a["beta"], bcoef["beta"]))
    for (scale, term), vals in stores.items():
        arr = np.asarray(vals, float) if vals else np.empty((0, 6))
        def pct(col):
            if arr.size == 0:
                return np.nan, np.nan, np.nan
            v = arr[:, col]
            v = v[np.isfinite(v)]
            if v.size == 0:
                return np.nan, np.nan, np.nan
            return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975)), float(2 * min(np.mean(v >= 0), np.mean(v <= 0)))
        ab_lo, ab_hi, ab_p = pct(0)
        pr_lo, pr_hi, pr_p = pct(1)
        rows.append(
            {
                "scale": scale,
                "outcome": term,
                "B_ok": int(arr.shape[0]),
                "ab_lo": ab_lo,
                "ab_hi": ab_hi,
                "ab_boot_p": ab_p,
                "prop_lo": pr_lo,
                "prop_hi": pr_hi,
                "prop_boot_p": pr_p,
            }
        )
    return pd.DataFrame(rows)


def bootstrap_nes(ph: pd.DataFrame, Y: np.ndarray, genes: np.ndarray, sets: dict, B: int, seed: int) -> pd.DataFrame:
    """Patient bootstrap, stratified by cohort, of the primary pseudobulk NES trio plus Q4."""
    rng = np.random.default_rng(seed)
    cohort = ph["cohort"].to_numpy()
    groups = {c: np.flatnonzero(cohort == c) for c in COHORTS}
    x = ph["tacstd2_pb"].to_numpy(float)
    m = ph["cldn4_pb"].to_numpy(float)
    # store nes per contrast/term
    keys = [
        "continuous",
        "continuous_resid_CLDN4",
        "CLDN4_low_continuous",
        "Q4_vs_Q1",
        "Q4_vs_Q1_resid_CLDN4",
        "delta_continuous_resid_minus_total",
        "delta_Q4_resid_minus_total",
    ]
    bag = {(k, t): [] for k in keys for t in PRIMARY_TERMS}
    primary_sets = {t: sets[t] for t in PRIMARY_TERMS}
    for b in range(B):
        idx = np.concatenate([rng.choice(groups[c], size=len(groups[c]), replace=True) for c in COHORTS])
        xb, mb, cb = x[idx], m[idx], cohort[idx]
        Yb = Y[idx]
        # quartiles recomputed within the bootstrap multiset, within cohort
        q = np.array([""] * len(idx), dtype=object)
        low = np.zeros(len(idx), dtype=bool)
        for c, gix in pd.Series(np.arange(len(idx))).groupby(cb).groups.items():
            gix = np.asarray(list(gix))
            q[gix] = assign_quartiles(pd.Series(xb[gix])).to_numpy()
            r = pd.Series(mb[gix]).rank(method="average")
            low[gix] = (r <= r.median()).to_numpy()
        built = {
            "continuous": (np.ones(len(idx), bool), xb, None),
            "continuous_resid_CLDN4": (np.ones(len(idx), bool), xb, mb),
            "CLDN4_low_continuous": (low, xb, None),
            "Q4_vs_Q1": (np.isin(q, ["Q1", "Q4"]), (q == "Q4").astype(float), None),
            "Q4_vs_Q1_resid_CLDN4": (np.isin(q, ["Q1", "Q4"]), (q == "Q4").astype(float), mb),
        }
        this = {}
        for key, (mask, xx, mm) in built.items():
            if mm is not None:
                mask = mask & np.isfinite(mm) & np.isfinite(xx)
            else:
                mask = mask & np.isfinite(xx)
            ix = np.flatnonzero(mask)
            if ix.size < 12 or np.unique(cb[ix]).size < 2:
                continue
            X, j = design_matrix(xx[ix], None if mm is None else mm[ix], cb[ix])
            tstat, _, _, df = ols_t(Yb[ix], X, j)
            if df < 5 or not np.isfinite(tstat).any():
                continue
            subrng = np.random.default_rng(int(rng.integers(0, 2**31 - 1)))
            rows = gsea_on_scores(tstat, genes, primary_sets, nperm=NPERM_BOOT, rng=subrng, with_lead=False)
            this[key] = {r["term"]: r["nes"] for r in rows}
            for r in rows:
                bag[(key, r["term"])].append(r["nes"])
        for src, dst, name in [
            ("continuous", "continuous_resid_CLDN4", "delta_continuous_resid_minus_total"),
            ("Q4_vs_Q1", "Q4_vs_Q1_resid_CLDN4", "delta_Q4_resid_minus_total"),
        ]:
            if src in this and dst in this:
                for term in PRIMARY_TERMS:
                    if term in this[src] and term in this[dst]:
                        bag[(name, term)].append(this[dst][term] - this[src][term])
        if (b + 1) % 50 == 0:
            print(f"  nes bootstrap {b+1}/{B}", flush=True)
    out = []
    for (key, term), vals in bag.items():
        v = np.asarray(vals, float)
        v = v[np.isfinite(v)]
        out.append(
            {
                "scale": "pseudobulk",
                "contrast": key,
                "term": term,
                "B_ok": int(v.size),
                "nes_boot_mean": float(np.mean(v)) if v.size else np.nan,
                "nes_lo": float(np.quantile(v, 0.025)) if v.size else np.nan,
                "nes_hi": float(np.quantile(v, 0.975)) if v.size else np.nan,
            }
        )
    return pd.DataFrame(out)


def leave_one_out(ph, Y, genes, sets, scores) -> pd.DataFrame:
    rows = []
    for drop in COHORTS:
        keep = ph["cohort"].to_numpy() != drop
        sub = ph.loc[keep].reset_index(drop=True)
        Ys = Y[keep]
        specs = contrast_specs(sub)
        want = {
            ("pseudobulk", "continuous"),
            ("pseudobulk", "continuous_resid_CLDN4"),
            ("pseudobulk", "CLDN4_low_continuous"),
            ("pseudobulk", "Q4_vs_Q1"),
            ("pseudobulk", "Q4_vs_Q1_resid_CLDN4"),
        }
        for spec in specs:
            if (spec["scale"], spec["contrast"]) not in want:
                continue
            got = eval_spec(spec, Ys, genes, {t: sets[t] for t in PRIMARY_TERMS}, NPERM, SEED + 17, False)
            for r in got:
                r["dropped"] = drop
                rows.append(r)
        # mediation ab on IFN-gamma
        sc = scores.loc[keep].reset_index(drop=True)
        med = module_models(sub, sc)
        med = med[(med["scale"] == "pseudobulk") & (med["outcome"].isin(PRIMARY_TERMS))]
        for rec in med.itertuples(index=False):
            rows.append(
                {
                    "scale": "pseudobulk",
                    "contrast": "mediation_ab",
                    "term": rec.outcome,
                    "nes": np.nan,
                    "es": np.nan,
                    "nom_p": np.nan,
                    "n_set": np.nan,
                    "mean_stat": np.nan,
                    "lead_genes": "",
                    "n_genes_ranked": np.nan,
                    "df": np.nan,
                    "n": rec.n,
                    "n_high": np.nan,
                    "n_low": np.nan,
                    "dropped": drop,
                    "ab": rec.ab,
                    "c": rec.c,
                    "c_prime": rec.c_prime,
                    "a": rec.a,
                    "b": rec.b,
                }
            )
    return pd.DataFrame(rows)


def partial_and_meta(ph: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scale, xcol, mcol in [
        ("pseudobulk", "tacstd2_pb", "cldn4_pb"),
        ("percent_positive", "pct_TACSTD2", "pct_CLDN4"),
    ]:
        for term in PRIMARY_TERMS:
            rhos, ns, prhos, pns = [], [], [], []
            for cohort, sub in ph.groupby("cohort"):
                y = scores.loc[sub.index, term].to_numpy(float)
                x = sub[xcol].to_numpy(float)
                m = sub[mcol].to_numpy(float)
                rho, p, n = spearman(x, y)
                rows.append({"scale": scale, "outcome": term, "cohort": cohort, "kind": "spearman", "rho": rho, "p": p, "n": n})
                ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(m)
                if ok.sum() >= 6:
                    Z = np.column_stack([np.ones(ok.sum()), m[ok]])
                    rx = residualize(x[ok], Z)
                    ry = residualize(y[ok], Z)
                    pr, pp, pn = spearman(rx, ry)
                else:
                    pr, pp, pn = np.nan, np.nan, int(ok.sum())
                rows.append({"scale": scale, "outcome": term, "cohort": cohort, "kind": "partial_spearman_CLDN4", "rho": pr, "p": pp, "n": pn})
                rhos.append(rho)
                ns.append(n)
                prhos.append(pr)
                pns.append(pn)
            dl = fisher_dl(np.asarray(rhos), np.asarray(ns), 0)
            rows.append({"scale": scale, "outcome": term, "cohort": "DL", "kind": "spearman", **{k: dl[k] for k in ("rho", "p", "I2", "lo", "hi", "k")}, "n": int(np.nansum(ns))})
            dl2 = fisher_dl(np.asarray(prhos, float), np.asarray(pns, float), 1)
            rows.append({"scale": scale, "outcome": term, "cohort": "DL", "kind": "partial_spearman_CLDN4", **{k: dl2[k] for k in ("rho", "p", "I2", "lo", "hi", "k")}, "n": int(np.nansum(pns))})
            # pooled residual spearman
            x = ph[xcol].to_numpy(float)
            m = ph[mcol].to_numpy(float)
            y = scores[term].to_numpy(float)
            cohort = ph["cohort"].to_numpy()
            ok = np.isfinite(x) & np.isfinite(m) & np.isfinite(y)
            Zc = np.column_stack([np.ones(ok.sum()), cohort_dummies(cohort[ok])])
            rho, p, n = spearman(residualize(x[ok], Zc), residualize(y[ok], Zc))
            rows.append({"scale": scale, "outcome": term, "cohort": "pooled_resid_cohort", "kind": "spearman", "rho": rho, "p": p, "n": n})
            Zm = np.column_stack([np.ones(ok.sum()), m[ok], cohort_dummies(cohort[ok])])
            rho, p, n = spearman(residualize(x[ok], Zm), residualize(y[ok], Zm))
            rows.append({"scale": scale, "outcome": term, "cohort": "pooled_resid_cohort", "kind": "partial_spearman_CLDN4", "rho": rho, "p": p, "n": n})
            low = ph["cldn4_low"].to_numpy() if scale == "pseudobulk" else ph["cldn4_pct_low"].to_numpy()
            ok2 = ok & low
            Zc2 = np.column_stack([np.ones(ok2.sum()), cohort_dummies(cohort[ok2])])
            rho, p, n = spearman(residualize(x[ok2], Zc2), residualize(y[ok2], Zc2))
            rows.append({"scale": scale, "outcome": term, "cohort": "CLDN4_low", "kind": "spearman_resid_cohort", "rho": rho, "p": p, "n": n})
    return pd.DataFrame(rows)


def calibration(sets: dict[str, set[str]]) -> pd.DataFrame:
    """PR #503 family-score contrast.

    TMM is fit on every expression unit (n=64), then the OLS uses the
    within-cohort CLDN4 %pos Q1 and Q4 labels. That is the published
    family-score path, not the Q-subset TMM used for the gene-level table.
    """
    locked = load_locked_cldn4_meta()
    counts, _ = load_intersect_counts()
    meta_all = locked[locked["patient"].isin(counts.columns)].drop_duplicates("patient")
    cts = filter_genes(counts.loc[:, meta_all["patient"]])
    lc = log_cpm(cts, tmm_norm_factors(cts))
    meta = meta_all[meta_all["cldn4_q_locked"].isin(["Q1", "Q4"])].copy()
    fam = {
        "IFN": sets["IFN_HALLMARK_UNION"],
        "IFN_gamma": sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        "IFN_alpha": sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        "MHC_I_APM": sets["MHC_I_APM"],
    }
    rows = []
    design = pd.DataFrame({"Intercept": 1.0, "CLDN4_Q4": (meta.set_index("patient")["cldn4_q_locked"] == "Q4").astype(float)})
    design = design.loc[meta["patient"]]
    for c in COHORTS:
        if c == REF:
            continue
        design[f"cohort_{c}"] = (meta.set_index("patient").loc[design.index, "cohort"] == c).astype(float).to_numpy()
    X = design.to_numpy(float)
    n1 = int((meta["cldn4_q_locked"] == "Q1").sum())
    n4 = int((meta["cldn4_q_locked"] == "Q4").sum())
    for name, genes in fam.items():
        present = [g for g in genes if g in lc.index]
        y = lc.loc[present, design.index].mean(axis=0).to_numpy(float)
        coef = ols_coef(y, X, 1)
        rows.append({"family": name, "n_genes": len(present), "n_q1": n1, "n_q4": n4, "logFC": coef["beta"], "se": coef["se"], "p": coef["p"], "df": coef["df"]})
    # CLDN4 itself
    if "CLDN4" in lc.index:
        coef = ols_coef(lc.loc["CLDN4", design.index].to_numpy(float), X, 1)
        rows.append({"family": "CLDN4_gene", "n_genes": 1, "n_q1": n1, "n_q4": n4, "logFC": coef["beta"], "se": coef["se"], "p": coef["p"], "df": coef["df"]})
    return pd.DataFrame(rows)


def stratum_range(ph: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scale, xcol, mcol, lowcol in [
        ("pseudobulk", "tacstd2_pb", "cldn4_pb", "cldn4_low"),
        ("percent_positive", "pct_TACSTD2", "pct_CLDN4", "cldn4_pct_low"),
    ]:
        for label, mask in [("all", np.ones(len(ph), bool)), ("CLDN4_low", ph[lowcol].to_numpy(bool))]:
            sub = ph.loc[mask]
            rho, p, n = spearman(sub[xcol].to_numpy(float), sub[mcol].to_numpy(float))
            rows.append(
                {
                    "scale": scale,
                    "stratum": label,
                    "n": int(mask.sum()),
                    "tacstd2_sd": float(sub[xcol].std(ddof=1)),
                    "cldn4_sd": float(sub[mcol].std(ddof=1)),
                    "spearman_tacstd2_cldn4": rho,
                    "spearman_p": p,
                    "n_by_cohort": ",".join(f"{c}:{int((sub.cohort==c).sum())}" for c in COHORTS),
                }
            )
    return pd.DataFrame(rows)


def plot_nes(nes: pd.DataFrame, boot: pd.DataFrame, path: Path) -> None:
    order_c = [
        ("Q4_vs_Q1", "Q4 vs Q1"),
        ("Q4_vs_Q1_resid_CLDN4", "Q4 vs Q1 | CLDN4"),
        ("continuous", "continuous"),
        ("continuous_resid_CLDN4", "continuous | CLDN4"),
        ("CLDN4_low_continuous", "CLDN4-low half"),
    ]
    terms = PRIMARY_TERMS
    term_lab = {
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α",
        "MHC_I_APM": "MHC-I/APM",
    }
    colors = {
        "Q4_vs_Q1": "#1f4e79",
        "Q4_vs_Q1_resid_CLDN4": "#c45911",
        "continuous": "#1f4e79",
        "continuous_resid_CLDN4": "#c45911",
        "CLDN4_low_continuous": "#548235",
    }
    sub = nes[(nes["scale"] == "pseudobulk") & (nes["contrast"].isin([c for c, _ in order_c]))].copy()
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ypos = []
    y = 0
    ytick = []
    yticklab = []
    for term in terms:
        for contrast, lab in order_c:
            rec = sub[(sub.term == term) & (sub.contrast == contrast)]
            b = boot[(boot.term == term) & (boot.contrast == contrast)]
            if rec.empty:
                continue
            nes_v = float(rec.iloc[0]["nes"])
            if not b.empty and np.isfinite(b.iloc[0]["nes_lo"]):
                lo, hi = float(b.iloc[0]["nes_lo"]), float(b.iloc[0]["nes_hi"])
            else:
                lo = hi = nes_v
            ax.plot([lo, hi], [y, y], color=colors[contrast], lw=1.6, solid_capstyle="round")
            ax.scatter([nes_v], [y], color=colors[contrast], s=28, zorder=3)
            ypos.append(y)
            ytick.append(y)
            yticklab.append(f"{term_lab[term]}  {lab}")
            y += 1
        y += 0.4
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(ytick)
    ax.set_yticklabels(yticklab)
    ax.set_xlabel("NES  (positive = enriched in TACSTD2-high)")
    ax.set_title("Concordant-4 malignant pseudobulk")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"))
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_mediation(med: pd.DataFrame, boot: pd.DataFrame, path: Path) -> None:
    term_lab = {
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α",
        "MHC_I_APM": "MHC-I/APM",
        "IFN_HALLMARK_UNION": "IFN union",
    }
    sub = med[med["scale"] == "pseudobulk"].copy()
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    labels = []
    y = 0
    for term in TERMS:
        rec = sub[sub.outcome == term]
        b = boot[(boot.scale == "pseudobulk") & (boot.outcome == term)]
        if rec.empty:
            continue
        ab = float(rec.iloc[0]["ab"])
        lo = float(b.iloc[0]["ab_lo"]) if not b.empty else ab
        hi = float(b.iloc[0]["ab_hi"]) if not b.empty else ab
        ax.plot([lo, hi], [y, y], color="#1f4e79", lw=1.6)
        ax.scatter([ab], [y], color="#1f4e79", s=28, zorder=3)
        labels.append(term_lab[term])
        y += 1
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Indirect effect a×b  (log2 score units)")
    ax.set_title("TACSTD2 → CLDN4 → module score")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"))
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatter(ph: pd.DataFrame, scores: pd.DataFrame, path: Path) -> None:
    y = scores["HALLMARK_INTERFERON_GAMMA_RESPONSE"].to_numpy(float)
    x = ph["tacstd2_pb"].to_numpy(float)
    m = ph["cldn4_pb"].to_numpy(float)
    cohort = ph["cohort"].to_numpy()
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(m)
    Z = np.column_stack([np.ones(ok.sum()), cohort_dummies(cohort[ok])])
    xr = residualize(x[ok], Z)
    yr = residualize(y[ok], Z)
    mr = residualize(m[ok], Z)
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    sc = ax.scatter(xr, yr, c=mr, cmap="viridis", s=28, edgecolor="none")
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("CLDN4 residual")
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.set_xlabel("TACSTD2 residual (cohort)")
    ax.set_ylabel("IFN-γ score residual (cohort)")
    ax.set_title("Malignant pseudobulk, n=64")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"))
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt(x, nd=3):
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def write_finding(cal, nes, boot, med, medb, partial, stratum, loo, ph) -> None:
    def nes_row(contrast, term):
        r = nes[(nes.scale == "pseudobulk") & (nes.contrast == contrast) & (nes.term == term)]
        b = boot[(boot.contrast == contrast) & (boot.term == term) & (boot.scale == "pseudobulk")]
        if r.empty:
            return None
        rec = r.iloc[0]
        lo = float(b.iloc[0]["nes_lo"]) if not b.empty else np.nan
        hi = float(b.iloc[0]["nes_hi"]) if not b.empty else np.nan
        return rec, lo, hi

    lines = []
    lines.append("# TACSTD2, CLDN4, and IFN/MHC on concordant-4 malignant pseudobulk")
    lines.append("")
    lines.append("Observational. Not a knockdown. The locked T/NK result is not re-fit.")
    lines.append("Cohorts: GSE123902 + GSE131907 + GSE205335 + GSE189357.")
    lines.append(f"Expression n = {len(ph)} (P4001 is in the percent-positive vector and out of the UMI sum).")
    lines.append("Positive NES means the set is enriched in TACSTD2-high.")
    lines.append("The NES interval is a patient bootstrap (stratified by cohort). The GSEA nominal p conditions on the ranking.")
    lines.append("")
    lines.append("## Calibration")
    lines.append("")
    lines.append("CLDN4 %pos within-cohort Q4 vs Q1, family-score OLS on log2(TMM-CPM+1), TMM fit on the Q1/Q4 subset. This is the PR #503 contrast. The published IFN union logFC is −0.584 (p = 0.0024) at 18 vs 16.")
    lines.append("")
    lines.append("| family | n genes | n Q1 / Q4 | logFC | p |")
    lines.append("|---|---:|---|---:|---:|")
    for rec in cal.itertuples(index=False):
        lines.append(f"| {rec.family} | {int(rec.n_genes)} | {int(rec.n_q1)} / {int(rec.n_q4)} | {fmt(rec.logFC)} | {rec.p:.3g} |")
    lines.append("")
    lines.append("## Pseudobulk NES")
    lines.append("")
    lines.append("Primary scale is malignant log2(TMM-CPM+1) of TACSTD2 and CLDN4. Residualized rows are the TACSTD2 coefficient after CLDN4 and cohort. The CLDN4-low row is the lower half by within-cohort CLDN4 rank, without a further CLDN4 covariate.")
    lines.append("")
    lines.append("| contrast | set | n | NES | bootstrap 95% | nominal p | mean rank stat |")
    lines.append("|---|---|---:|---:|---|---:|---:|")
    order = [
        "Q4_vs_Q1",
        "Q4_vs_Q1_resid_CLDN4",
        "continuous",
        "continuous_resid_CLDN4",
        "CLDN4_low_continuous",
    ]
    for contrast in order:
        for term in PRIMARY_TERMS:
            got = nes_row(contrast, term)
            if got is None:
                continue
            rec, lo, hi = got
            nlab = f"{int(rec.n)}"
            if contrast.startswith("Q4") and np.isfinite(rec.n_low):
                nlab = f"{int(rec.n_low)} vs {int(rec.n_high)}"
            lines.append(
                f"| {contrast} | {term} | {nlab} | {fmt(rec.nes)} | {fmt(lo)} to {fmt(hi)} | {rec.nom_p:.3g} | {fmt(rec.mean_stat)} |"
            )
    lines.append("")
    lines.append("Paired bootstrap of residualized NES minus total NES. A negative downregulation that attenuates would move this difference upward through zero from a negative total NES. These rows are the difference itself.")
    lines.append("")
    lines.append("| contrast | set | B | mean | 95% |")
    lines.append("|---|---|---:|---:|---|")
    for contrast in ["delta_Q4_resid_minus_total", "delta_continuous_resid_minus_total"]:
        for term in PRIMARY_TERMS:
            b = boot[(boot.contrast == contrast) & (boot.term == term) & (boot.scale == "pseudobulk")]
            if b.empty:
                continue
            rec = b.iloc[0]
            lines.append(
                f"| {contrast} | {term} | {int(rec.B_ok)} | {fmt(rec.nes_boot_mean)} | {fmt(rec.nes_lo)} to {fmt(rec.nes_hi)} |"
            )
    lines.append("")
    lines.append("## Mediation on the module score")
    lines.append("")
    lines.append("Cohort-adjusted OLS. a is CLDN4 ~ TACSTD2. b is the module ~ CLDN4 given TACSTD2. c is the total TACSTD2 slope. c′ is the slope after CLDN4. Indirect effect is a×b. Units are log2 score per log2 TACSTD2 (pseudobulk) or per percentage point (percent-positive). The bootstrap resamples patients inside cohort.")
    lines.append("")
    lines.append("| scale | outcome | n | a | b | c | c′ | a×b | a×b 95% | VIF |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---|---:|")
    for rec in med.itertuples(index=False):
        b = medb[(medb.scale == rec.scale) & (medb.outcome == rec.outcome)]
        lo = float(b.iloc[0]["ab_lo"]) if not b.empty else np.nan
        hi = float(b.iloc[0]["ab_hi"]) if not b.empty else np.nan
        lines.append(
            f"| {rec.scale} | {rec.outcome} | {int(rec.n)} | {fmt(rec.a)} | {fmt(rec.b)} | {fmt(rec.c)} | {fmt(rec.c_prime)} | {fmt(rec.ab)} | {fmt(lo)} to {fmt(hi)} | {fmt(rec.vif_tacstd2, 2)} |"
        )
    lines.append("")
    lines.append("## Cohort-adjusted Spearman")
    lines.append("")
    lines.append("| scale | outcome | kind | where | rho | p | n |")
    lines.append("|---|---|---|---|---:|---:|---:|")
    show = partial[
        partial["cohort"].isin(["DL", "pooled_resid_cohort", "CLDN4_low"])
        & (partial["scale"] == "pseudobulk")
    ]
    for rec in show.itertuples(index=False):
        lines.append(
            f"| {rec.scale} | {rec.outcome} | {rec.kind} | {rec.cohort} | {fmt(rec.rho)} | {rec.p:.3g} | {int(rec.n) if np.isfinite(rec.n) else 'NA'} |"
        )
    lines.append("")
    lines.append("## Range in the CLDN4-low half")
    lines.append("")
    lines.append("| scale | stratum | n | TACSTD2 SD | CLDN4 SD | Spearman TACSTD2–CLDN4 |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for rec in stratum.itertuples(index=False):
        lines.append(
            f"| {rec.scale} | {rec.stratum} | {int(rec.n)} | {fmt(rec.tacstd2_sd)} | {fmt(rec.cldn4_sd)} | {fmt(rec.spearman_tacstd2_cldn4)} |"
        )
    lines.append("")
    lines.append("## Leave-one-cohort-out (pseudobulk NES and a×b)")
    lines.append("")
    if "dropped" in loo.columns:
        lines.append("| dropped | contrast | set | NES or a×b |")
        lines.append("|---|---|---|---:|")
        for rec in loo.itertuples(index=False):
            val = rec.nes if rec.contrast != "mediation_ab" else rec.ab
            lines.append(f"| {rec.dropped} | {rec.contrast} | {rec.term} | {fmt(val)} |")
    lines.append("")
    lines.append("Numbers are written by `analyze.py` from the tables in `results/tables/`.")
    text = "\n".join(lines) + "\n"
    (HERE / "results" / "AUTO_TABLES.md").write_text(text)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    print("assemble", flush=True)
    ph, logc, sets = assemble()
    ph.to_csv(TABLES / "phenotype.tsv", sep="\t", index=False)
    print("n", len(ph), "genes", logc.shape[0], "P4001 in ph", "P4001" in set(ph.patient))
    print("quartile counts\n", ph.groupby(["cohort", "tacstd2_q"]).size())
    scores = score_matrix(logc, ph, sets)
    scores.to_csv(TABLES / "module_scores.tsv", sep="\t", index=False)
    # align scores to ph row order
    scores = scores.set_index("col").loc[ph["col"]].reset_index()

    print("calibration", flush=True)
    cal = calibration(sets)
    cal.to_csv(TABLES / "calibration_cldn4_q4q1.tsv", sep="\t", index=False)
    print(cal.to_string(index=False))

    Y, genes = Y_of(logc, ph["col"].tolist())
    print("point NES", flush=True)
    nes_rows = []
    for i, spec in enumerate(contrast_specs(ph)):
        got = eval_spec(spec, Y, genes, sets, NPERM, SEED + i * 100, True)
        nes_rows.extend(got)
        print(" ", spec["scale"], spec["contrast"], "rows", len(got), flush=True)
    nes = pd.DataFrame(nes_rows)
    # BH within scale × contrast across the four terms
    nes["fdr"] = np.nan
    for (scale, contrast), idx in nes.groupby(["scale", "contrast"]).groups.items():
        nes.loc[idx, "fdr"] = bh(nes.loc[idx, "nom_p"].to_numpy())
    nes.to_csv(TABLES / "nes.tsv", sep="\t", index=False)

    print("module mediation", flush=True)
    med = module_models(ph, scores)
    print("mediation bootstrap", flush=True)
    medb = bootstrap_mediation(ph, scores, B_MED, SEED + 7)
    med = med.merge(medb, on=["scale", "outcome"], how="left")
    med["prop"] = med["ab"] / med["c"]
    med.to_csv(TABLES / "mediation.tsv", sep="\t", index=False)
    print(med[med.scale == "pseudobulk"][["outcome", "a", "b", "c", "c_prime", "ab", "ab_lo", "ab_hi"]].to_string(index=False))

    print("partial spearman", flush=True)
    partial = partial_and_meta(ph, scores)
    partial.to_csv(TABLES / "spearman_partial.tsv", sep="\t", index=False)

    stratum = stratum_range(ph)
    stratum.to_csv(TABLES / "stratum_range.tsv", sep="\t", index=False)
    print(stratum.to_string(index=False))

    print("leave one out", flush=True)
    loo = leave_one_out(ph, Y, genes, sets, scores)
    loo.to_csv(TABLES / "leave_one_cohort.tsv", sep="\t", index=False)

    print("NES bootstrap", flush=True)
    boot = bootstrap_nes(ph, Y, genes, sets, B_NES, SEED + 99)
    ph_pct = ph.copy()
    ph_pct["tacstd2_pb"] = ph["pct_TACSTD2"].to_numpy(float)
    ph_pct["cldn4_pb"] = ph["pct_CLDN4"].to_numpy(float)
    ph_pct["cldn4_low"] = ph["cldn4_pct_low"].to_numpy(bool)
    for cohort, idx in ph_pct.groupby("cohort").groups.items():
        ph_pct.loc[idx, "tacstd2_q"] = assign_quartiles(ph_pct.loc[idx, "tacstd2_pb"])
    print("NES bootstrap percent-positive", flush=True)
    boot_pct = bootstrap_nes(ph_pct, Y, genes, sets, B_NES, SEED + 199)
    boot_pct["scale"] = "percent_positive"
    boot_pct.to_csv(TABLES / "nes_bootstrap_percent.tsv", sep="\t", index=False)
    boot = pd.concat([boot, boot_pct], ignore_index=True)
    boot.to_csv(TABLES / "nes_bootstrap.tsv", sep="\t", index=False)
    print(boot.to_string(index=False))

    plot_nes(nes, boot, FIGS / "forest_nes")
    plot_mediation(med, medb, FIGS / "forest_mediation_ab")
    plot_scatter(ph, scores, FIGS / "scatter_tacstd2_ifng")
    write_finding(cal, nes, boot, med, medb, partial, stratum, loo, ph)
    print("done", flush=True)


if __name__ == "__main__":
    main()

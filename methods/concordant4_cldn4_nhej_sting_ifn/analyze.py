#!/usr/bin/env python3
"""Concordant-4 malignant CLDN4-high vs low: NHEJ, cGAS-STING, IFN/APM, mediation.

ADDITIVE. CLDN4-only. Patient / donor / sample is the unit.
Inputs are the locked malignant UMI-sums from PR #503 / #539
(GSE123902 + GSE131907 + GSE205335 + GSE189357). Not a re-audit of the
T/NK rho. Not GSE148071 / GSE127465 / GSE207422 / GSE154826.

Primary contrast: within-cohort malignant CLDN4 %pos Q4 vs Q1, stacked
OLS on mean log2(TMM-CPM+1) with cohort covariates. Same normalization
path as the locked IFN/MHC family scores.

Primary mediation: CLDN4 %pos (global z) → NHEJ module → IFN module,
cohort fixed effects, stratified bootstrap of the product a*b.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
COMBO = "+".join(COHORTS)
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
BOOT_B = 5000
BOOT_SEED = 20260921

# KEGG hsa03450 (MSigDB KEGG_NON_HOMOLOGOUS_END_JOINING). Primary NHEJ module.
NHEJ_KEGG = [
    "DCLRE1C", "DNTT", "FEN1", "LIG4", "MRE11", "NHEJ1", "POLL", "POLM",
    "PRKDC", "RAD50", "XRCC4", "XRCC5", "XRCC6",
]
# Catalytic c-NHEJ: KEGG minus TdT / FEN1 / MRN, plus PAXX, APLF, PNKP.
NHEJ_CORE = [
    "XRCC4", "XRCC5", "XRCC6", "PRKDC", "LIG4", "NHEJ1", "DCLRE1C",
    "POLL", "POLM", "PAXX", "APLF", "PNKP",
]
# cGAS–STING machinery. IFN-hallmark members (IRF7, ZBP1, TRIM21, SAMHD1)
# and NHEJ members (PRKDC, XRCC5, XRCC6, MRE11) are held out so this score
# is not the IFN score and not the NHEJ score.
# Sources: Reactome R-HSA-1834941 / R-HSA-1834949, GO:0140896.
CGAS_STING = [
    "CGAS", "STING1", "TBK1", "IKBKE", "IRF3", "IFI16", "DDX41",
    "TRIM56", "DHX9", "DHX36", "PQBP1",
]
# Tight sensor-adapter set. Sensitivity, because DHX/IFI16/TRIM56 can track IFN.
CGAS_STING_CORE = ["CGAS", "STING1", "TBK1", "IRF3"]
# GEO builds disagree on these symbols. Collapse to the current HGNC name
# before scoring. IFN/APM symbols are not relabeled.
ALIASES = {
    "CGAS": ["CGAS", "MB21D1", "C6ORF150"],
    "STING1": ["STING1", "TMEM173"],
    "MRE11": ["MRE11", "MRE11A"],
}
# Negative regulators, reported as context, not the primary module.
CGAS_STING_NEG = ["TREX1", "ENPP1", "NLRC3", "NLRP4"]
PROLIF = ["MKI67", "PCNA", "TOP2A", "MCM2", "CDK1", "CCNB1", "BIRC5", "UBE2C"]

PRIMARY_MODULES = ["NHEJ", "cGAS-STING", "IFN", "APM"]
MODULE_COLORS = {
    "NHEJ": "#4c78a8",
    "NHEJ-core": "#9ecae1",
    "cGAS-STING": "#f58518",
    "cGAS-STING-core": "#ffbf70",
    "cGAS-STING-neg": "#f2cf5b",
    "IFN": "#e45756",
    "APM": "#54a24b",
    "proliferation": "#b279a2",
}
COHORT_COLORS = {
    "GSE123902": "#4c78a8",
    "GSE131907": "#f58518",
    "GSE205335": "#54a24b",
    "GSE189357": "#b279a2",
}


def load_ifn_apm() -> tuple[list[str], list[str]]:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    ifn = sorted(
        set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
    )
    apm = list(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"])
    return ifn, apm


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def load_units() -> dict[str, dict]:
    units: dict[str, dict] = {}
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    units["GSE123902"] = {
        "unit": "donor",
        "n": int(len(el)),
        "pct": el["mal_CLDN4_pct"].to_numpy(float),
        "ids": el["patient"].astype(str).tolist(),
    }
    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    units["GSE131907"] = {
        "unit": "sample",
        "n": int(len(mal)),
        "pct": mal["mal_CLDN4_pct"].to_numpy(float),
        "ids": mal["sample"].astype(str).tolist(),
    }
    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    units["GSE205335"] = {
        "unit": "patient",
        "n": int(len(d)),
        "pct": d["mal_CLDN4_pct_pos"].to_numpy(float),
        "ids": d["patient"].astype(str).tolist(),
    }
    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    units["GSE189357"] = {
        "unit": "patient",
        "n": int(len(el)),
        "pct": el["mal_CLDN4_pct"].to_numpy(float),
        "ids": el["patient"].astype(str).tolist(),
    }
    return units


def four_meta(units: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for c in COHORTS:
        u = units[c]
        pct = np.asarray(u["pct"], float)
        if np.nanmax(pct) <= 1.5:
            pct = pct * 100.0
        q = assign_quartiles(pd.Series(pct, index=u["ids"]))
        for i, pid in enumerate(u["ids"]):
            rows.append(
                {
                    "patient": pid,
                    "cohort": c,
                    "unit": u["unit"],
                    "cldn4_pct": float(pct[i]),
                    "quartile": str(q.loc[pid]),
                }
            )
    meta = pd.DataFrame(rows)
    in_mat = set()
    for name in COHORTS:
        cols = pd.read_csv(DATA / f"{name}_malignant_counts.tsv.gz", sep="\t", nrows=0).columns.astype(str)
        in_mat.update(cols.tolist())
    meta["in_count_matrix"] = meta["patient"].isin(in_mat)
    return meta


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


def ols_fit(y: np.ndarray, X: np.ndarray) -> dict:
    """OLS. y is (n,), X is (n, p) with no missing values."""
    n, p = X.shape
    df = n - p
    if df < 1:
        return {}
    xtx = X.T @ X
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ X.T @ y
    resid = y - X @ beta
    sigma2 = float(np.sum(resid**2) / df)
    se = np.sqrt(np.maximum(np.diag(xtx_inv) * sigma2, 0.0))
    t = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    pv = 2.0 * stats.t.sf(np.abs(t), df)
    return {"beta": beta, "se": se, "t": t, "p": pv, "df": df, "sigma2": sigma2, "xtx_inv": xtx_inv}


def design_matrix(df: pd.DataFrame, cols: list[str], with_cohort: bool = True) -> tuple[np.ndarray, list[str]]:
    names = ["Intercept"]
    blocks = [np.ones(len(df))]
    for c in cols:
        names.append(c)
        blocks.append(df[c].to_numpy(float))
    if with_cohort and df["cohort"].nunique() > 1:
        for c in COHORTS:
            if c == REF:
                continue
            names.append(f"cohort_{c}")
            blocks.append((df["cohort"] == c).to_numpy(float))
    return np.column_stack(blocks), names


def coef_of(fit: dict, names: list[str], name: str) -> dict:
    j = names.index(name)
    return {
        "beta": float(fit["beta"][j]),
        "se": float(fit["se"][j]),
        "t": float(fit["t"][j]),
        "p": float(fit["p"][j]),
        "df": int(fit["df"]),
    }


def load_parts() -> list[pd.DataFrame]:
    return [read_counts(DATA / f"{c}_malignant_counts.tsv.gz") for c in COHORTS]


def load_combined(parts: list[pd.DataFrame]) -> pd.DataFrame:
    genes = sorted(set.intersection(*[set(p.index) for p in parts]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts], axis=1)
    return combined.loc[:, ~combined.columns.duplicated()]


def harmonized_row(parts: list[pd.DataFrame], aliases: list[str]) -> pd.Series:
    chunks = []
    for df in parts:
        hit = [a for a in aliases if a in df.index]
        if not hit:
            chunks.append(pd.Series(0.0, index=df.columns))
            continue
        sub = df.loc[hit]
        if isinstance(sub, pd.Series):
            chunks.append(sub.astype(float))
        else:
            chunks.append(sub.astype(float).sum(axis=0))
    out = pd.concat(chunks)
    return out.groupby(level=0).sum()


def add_alias_genes(lc: pd.DataFrame, cts: pd.DataFrame, factors: pd.Series, parts: list[pd.DataFrame]) -> pd.DataFrame:
    """Put alias-collapsed pathway genes on the locked TMM scale.

    TMM factors stay those of the shared-symbol matrix, so IFN/APM scores
    do not move. A gene already present under the canonical name is replaced
    only when an alias contributed counts (MRE11 vs MRE11A never co-occur
    inside one cohort).
    """
    lib = cts.sum(axis=0).astype(float) * factors.reindex(cts.columns).astype(float)
    lib = lib.reindex(lc.columns)
    for canon, aliases in ALIASES.items():
        counts = harmonized_row(parts, aliases).reindex(lc.columns).fillna(0.0)
        lc.loc[canon] = np.log2(counts.to_numpy(float) / lib.to_numpy(float) * 1e6 + 1.0)
    return lc


def module_genes(ifn: list[str], apm: list[str]) -> dict[str, list[str]]:
    return {
        "NHEJ": NHEJ_KEGG,
        "NHEJ-core": NHEJ_CORE,
        "cGAS-STING": CGAS_STING,
        "cGAS-STING-core": CGAS_STING_CORE,
        "cGAS-STING-neg": CGAS_STING_NEG,
        "IFN": ifn,
        "APM": apm,
        "proliferation": PROLIF,
    }


def score_modules(logcpm: pd.DataFrame, genesets: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = {}
    membership = []
    for name, genes in genesets.items():
        present = [g for g in genes if g in logcpm.index]
        missing = [g for g in genes if g not in logcpm.index]
        for g in present:
            membership.append({"module": name, "gene": g, "in_logcpm": True})
        for g in missing:
            membership.append({"module": name, "gene": g, "in_logcpm": False})
        if len(present) < 3:
            continue
        scores[name] = logcpm.loc[present].astype(float).mean(axis=0)
    return pd.DataFrame(scores), pd.DataFrame(membership)


def gene_level_q4(logcpm: pd.DataFrame, meta: pd.DataFrame, genesets: dict[str, list[str]]) -> pd.DataFrame:
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(logcpm.columns)].copy()
    y_index = m["patient"].tolist()
    d = m.set_index("patient").loc[y_index].reset_index()
    X, names = design_matrix(d, [])
    # add Q4 indicator
    q4 = (d["quartile"] == "Q4").to_numpy(float)
    X = np.column_stack([X, q4])
    names = names + ["CLDN4_Q4"]
    rows = []
    wanted = []
    for name, genes in genesets.items():
        if name not in ("NHEJ", "NHEJ-core", "cGAS-STING", "cGAS-STING-core", "cGAS-STING-neg", "APM"):
            continue
        for g in genes:
            if g in logcpm.index:
                wanted.append((name, g))
    # also a short IFN sentinel list
    for g in ["STAT1", "IRF1", "ISG15", "MX1", "B2M", "TAP1", "TAP2", "HLA-A", "CLDN4"]:
        if g in logcpm.index:
            wanted.append(("sentinel", g))
    seen = set()
    for name, g in wanted:
        key = (name, g)
        if key in seen:
            continue
        seen.add(key)
        y = logcpm.loc[g, y_index].to_numpy(float)
        fit = ols_fit(y, X)
        if not fit:
            continue
        c = coef_of(fit, names, "CLDN4_Q4")
        rows.append({
            "module": name,
            "gene": g,
            "logFC": c["beta"],
            "se": c["se"],
            "t": c["t"],
            "p": c["p"],
            "df": c["df"],
            "n": int(len(d)),
            "n_q1": int((d["quartile"] == "Q1").sum()),
            "n_q4": int((d["quartile"] == "Q4").sum()),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_within_module"] = out.groupby("module")["p"].transform(lambda s: _bh(s.to_numpy()))
    return out


def module_contrast(scores: pd.DataFrame, meta: pd.DataFrame, split: str, cohort: str | None) -> pd.DataFrame:
    m = meta.loc[meta["patient"].isin(scores.index)].copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    if split == "q4q1":
        m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
        if (m["quartile"] == "Q1").sum() < 3 or (m["quartile"] == "Q4").sum() < 3:
            return pd.DataFrame()
    elif len(m) < 8:
        return pd.DataFrame()
    rows = []
    for name in scores.columns:
        d = m.copy()
        d["y"] = scores.loc[d["patient"], name].to_numpy(float)
        if split == "q4q1":
            d["x"] = (d["quartile"] == "Q4").astype(float)
            coef = "CLDN4_Q4"
        else:
            z = d["cldn4_pct"].astype(float)
            d["x"] = (z - z.mean()) / z.std(ddof=1)
            coef = "CLDN4_pct_z"
        use_cohort = cohort is None and d["cohort"].nunique() > 1
        X, names = design_matrix(d, ["x"], with_cohort=use_cohort)
        # design_matrix names the column "x"; rename conceptually
        names = ["CLDN4_Q4" if (n == "x" and split == "q4q1") else "CLDN4_pct_z" if n == "x" else n for n in names]
        fit = ols_fit(d["y"].to_numpy(float), X)
        if not fit:
            continue
        c = coef_of(fit, names, coef)
        q1 = d.loc[d["quartile"] == "Q1", "y"] if split == "q4q1" else pd.Series(dtype=float)
        q4 = d.loc[d["quartile"] == "Q4", "y"] if split == "q4q1" else pd.Series(dtype=float)
        r_rb, p_mw = float("nan"), float("nan")
        if split == "q4q1" and len(q1) >= 3 and len(q4) >= 3:
            u, p_mw = stats.mannwhitneyu(q4.to_numpy(), q1.to_numpy(), alternative="two-sided")
            r_rb = (2.0 * float(u)) / (len(q4) * len(q1)) - 1.0
            p_mw = float(p_mw)
        tcrit = float(stats.t.ppf(0.975, c["df"]))
        rows.append(
            {
                "module": name,
                "split": split,
                "cohort": cohort or COMBO,
                "n": int(len(d)),
                "n_q1": int((d["quartile"] == "Q1").sum()) if split == "q4q1" else np.nan,
                "n_q4": int((d["quartile"] == "Q4").sum()) if split == "q4q1" else np.nan,
                "logFC": c["beta"],
                "se": c["se"],
                "t": c["t"],
                "df": c["df"],
                "p": c["p"],
                "ci_lo": c["beta"] - tcrit * c["se"],
                "ci_hi": c["beta"] + tcrit * c["se"],
                "median_q1": float(q1.median()) if split == "q4q1" else np.nan,
                "median_q4": float(q4.median()) if split == "q4q1" else np.nan,
                "r_rb": float(r_rb),
                "p_mw": p_mw,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = _bh(out["p"].to_numpy())
    return out


def spearman_meta(scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    m = meta.loc[meta["patient"].isin(scores.index)].copy()
    for name in scores.columns:
        rhos, ns, ps = [], [], []
        for c in COHORTS:
            sub = m.loc[m["cohort"] == c]
            x = sub["cldn4_pct"].to_numpy(float)
            y = scores.loc[sub["patient"], name].to_numpy(float)
            if len(sub) < 4:
                rho, p, n = float("nan"), float("nan"), int(len(sub))
            else:
                rho, p = stats.spearmanr(x, y)
                rho, p, n = float(rho), float(p), int(len(sub))
            rows.append({"module": name, "cohort": c, "n": n, "rho": rho, "p": p, "level": "cohort"})
            if math.isfinite(rho):
                rhos.append(rho)
                ns.append(n)
                ps.append(p)
        pooled = random_effects_dl(rhos, ns)
        rows.append(
            {
                "module": name,
                "cohort": COMBO,
                "n": pooled.get("n_patients_total", 0),
                "rho": pooled.get("pooled_rho", float("nan")),
                "p": pooled.get("p", float("nan")),
                "I2": pooled.get("I2", float("nan")),
                "ci_lo": (pooled.get("ci95_rho") or [float("nan"), float("nan")])[0],
                "ci_hi": (pooled.get("ci95_rho") or [float("nan"), float("nan")])[1],
                "k": pooled.get("k", 0),
                "level": "DL",
            }
        )
    return pd.DataFrame(rows)


def fisher_z(rho: float) -> float:
    r = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(r))


def random_effects_dl(rhos: list[float], ns: list[int]) -> dict:
    z = np.array([fisher_z(r) for r in rhos], dtype=float)
    v = np.array([1.0 / (n - 3) if n > 3 else np.nan for n in ns], dtype=float)
    ok = np.isfinite(z) & np.isfinite(v) & (v > 0)
    z, v, n_ok = z[ok], v[ok], np.array(ns, dtype=float)[ok]
    k = int(z.size)
    if k == 0:
        return {"k": 0}
    w = 1.0 / v
    z_fe = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - z_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else float("nan")
    tau2 = max(0.0, (q - df) / c) if k > 1 and c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se = float(1.0 / math.sqrt(float(np.sum(w_re))))
    z_stat = z_re / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z_stat))) if math.isfinite(z_stat) else float("nan")
    lo, hi = z_re - 1.96 * se, z_re + 1.96 * se
    i2 = max(0.0, (q - df) / q) * 100.0 if k > 1 and q > 0 else 0.0
    return {
        "k": k,
        "n_patients_total": int(n_ok.sum()),
        "pooled_rho": float(np.tanh(z_re)),
        "p": p,
        "I2": i2,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
    }


def _analysis_frame(scores: pd.DataFrame, meta: pd.DataFrame, xcol: str, mcol: str, ycol: str, extra: list[str] | None = None) -> pd.DataFrame:
    d = meta.loc[meta["patient"].isin(scores.index)].copy()
    d = d.set_index("patient")
    d["Y"] = scores[ycol]
    d["M"] = scores[mcol]
    if xcol == "cldn4_expr":
        d["Xraw"] = scores["CLDN4"]
    elif xcol in d.columns:
        d["Xraw"] = d[xcol].astype(float)
    else:
        d["Xraw"] = scores[xcol]
    extra = extra or []
    for e in extra:
        d[e] = scores[e]
    d = d.reset_index()
    z = d["Xraw"].astype(float)
    d["X"] = (z - z.mean()) / z.std(ddof=1)
    return d


def fit_mediation(d: pd.DataFrame, extra: list[str] | None = None) -> dict:
    extra = extra or []
    # a: M ~ X + extra + cohort
    Xa, na = design_matrix(d, ["X", *extra])
    # rename: design uses column names from d — wait, design_matrix reads d[col]
    # so d must contain X and extra. Yes.
    fa = ols_fit(d["M"].to_numpy(float), Xa)
    # c: Y ~ X + extra + cohort
    Xc, nc = design_matrix(d, ["X", *extra])
    fc = ols_fit(d["Y"].to_numpy(float), Xc)
    # b, c': Y ~ X + M + extra + cohort
    Xb, nb = design_matrix(d, ["X", "M", *extra])
    fb = ols_fit(d["Y"].to_numpy(float), Xb)
    a = coef_of(fa, na, "X")
    c = coef_of(fc, nc, "X")
    b = coef_of(fb, nb, "M")
    cp = coef_of(fb, nb, "X")
    ab = a["beta"] * b["beta"]
    # Sobel SE
    se_ab = math.sqrt((a["beta"] ** 2) * (b["se"] ** 2) + (b["beta"] ** 2) * (a["se"] ** 2))
    z_sobel = ab / se_ab if se_ab > 0 else float("nan")
    p_sobel = float(2 * stats.norm.sf(abs(z_sobel))) if math.isfinite(z_sobel) else float("nan")
    prop = ab / c["beta"] if c["beta"] != 0 else float("nan")
    return {
        "n": int(len(d)),
        "a": a["beta"], "a_se": a["se"], "a_p": a["p"], "a_df": a["df"],
        "b": b["beta"], "b_se": b["se"], "b_p": b["p"],
        "c": c["beta"], "c_se": c["se"], "c_p": c["p"],
        "c_prime": cp["beta"], "c_prime_se": cp["se"], "c_prime_p": cp["p"],
        "ab": ab,
        "sobel_se": se_ab,
        "sobel_p": p_sobel,
        "prop_mediated": prop,
        "df": int(fb["df"]),
    }


def bootstrap_mediation(d: pd.DataFrame, extra: list[str] | None = None, B: int = BOOT_B) -> dict:
    rng = np.random.default_rng(BOOT_SEED)
    groups = {c: d.index[d["cohort"] == c].to_numpy() for c in COHORTS}
    abs_, props = [], []
    for _ in range(B):
        draw = []
        ok = True
        for c in COHORTS:
            idx = groups[c]
            if len(idx) == 0:
                continue
            take = rng.choice(idx, size=len(idx), replace=True)
            draw.append(take)
        if not draw:
            continue
        sub = d.loc[np.concatenate(draw)].copy()
        # re-standardize X inside the draw so the unit stays "1 SD"
        z = sub["Xraw"].astype(float)
        if float(z.std(ddof=1)) == 0 or not math.isfinite(float(z.std(ddof=1))):
            continue
        sub["X"] = (z - z.mean()) / z.std(ddof=1)
        try:
            fit = fit_mediation(sub, extra)
        except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
            continue
        if not math.isfinite(fit["ab"]):
            continue
        abs_.append(fit["ab"])
        if math.isfinite(fit["prop_mediated"]):
            props.append(fit["prop_mediated"])
    ab = np.asarray(abs_, float)
    pr = np.asarray(props, float)
    def pct(a):
        if a.size == 0:
            return float("nan"), float("nan")
        return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))
    lo, hi = pct(ab)
    plo, phi = pct(pr)
    if ab.size:
        p = float(2 * min(np.mean(ab <= 0), np.mean(ab >= 0)))
        p = min(p, 1.0)
    else:
        p = float("nan")
    return {
        "boot_n": int(ab.size),
        "ab_ci_lo": lo,
        "ab_ci_hi": hi,
        "ab_boot_p": p,
        "prop_ci_lo": plo,
        "prop_ci_hi": phi,
        "ab_boot": ab,
    }


def run_one_mediation(scores, meta, label, mcol, ycol, xcol="cldn4_pct", extra=None) -> tuple[dict, np.ndarray]:
    d = _analysis_frame(scores, meta, xcol, mcol, ycol, extra)
    # _analysis_frame looks up xcol on meta. cldn4_pct is there. cldn4_expr handled inside.
    point = fit_mediation(d, extra)
    boot = bootstrap_mediation(d, extra)
    point.update({k: v for k, v in boot.items() if k != "ab_boot"})
    point["label"] = label
    point["mediator"] = mcol
    point["outcome"] = ycol
    point["exposure"] = xcol
    point["covariates"] = "cohort" + (("+" + "+".join(extra)) if extra else "")
    return point, boot["ab_boot"]


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(float(x))):
        return "—"
    return f"{float(x):.{nd}f}"


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and not math.isfinite(float(p))):
        return "—"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def _signed(x, nd=3):
    if x is None or not math.isfinite(float(x)):
        return "—"
    return f"{float(x):+.{nd}f}"


def save_fig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_forest(contrasts: pd.DataFrame, path: Path) -> None:
    sub = contrasts[(contrasts["split"] == "q4q1") & (contrasts["module"].isin(PRIMARY_MODULES))].copy()
    order = PRIMARY_MODULES
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), gridspec_kw={"width_ratios": [1.15, 1.6]})
    # stacked
    ax = axes[0]
    st = sub[sub["cohort"] == COMBO].set_index("module").loc[order]
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color="#888", lw=0.8)
    ax.errorbar(
        st["logFC"], y,
        xerr=[st["logFC"] - st["ci_lo"], st["ci_hi"] - st["logFC"]],
        fmt="o", color="#222", ecolor="#222", ms=6, capsize=3,
    )
    for yi, name in zip(y, order):
        ax.scatter(st.loc[name, "logFC"], yi, c=MODULE_COLORS[name], s=40, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{n}   p={_fmt_p(st.loc[n, 'p'])}" for n in order], fontsize=8
    )
    ax.set_xlabel("log2FC  (Q4 − Q1)")
    ax.set_title(f"Stacked Q4 vs Q1\nn={int(st['n'].iloc[0])}  ({int(st['n_q1'].iloc[0])}/{int(st['n_q4'].iloc[0])})")
    # by cohort
    ax = axes[1]
    cohorts = [c for c in COHORTS if c in set(sub["cohort"])]
    # GSE189357 binary is skipped
    ypos = []
    labels = []
    y = 0
    for name in order[::-1]:
        for c in cohorts:
            hit = sub[(sub["module"] == name) & (sub["cohort"] == c)]
            if hit.empty:
                continue
            r = hit.iloc[0]
            ax.errorbar(
                r["logFC"], y,
                xerr=[[r["logFC"] - r["ci_lo"]], [r["ci_hi"] - r["logFC"]]],
                fmt="o", color=COHORT_COLORS[c], ms=5, capsize=2,
            )
            ypos.append(y)
            labels.append(f"{name} · {c}")
            y += 1
        y += 0.6
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("log2FC  (Q4 − Q1)")
    ax.set_title("Within cohort (Q4 n≥3)")
    fig.suptitle("Malignant CLDN4-high vs low  ·  module scores", fontsize=12)
    fig.tight_layout()
    save_fig(fig, path)


def plot_boxes(scores: pd.DataFrame, meta: pd.DataFrame, path: Path) -> None:
    m = meta.loc[meta["quartile"].isin(["Q1", "Q4"]) & meta["patient"].isin(scores.index)].copy()
    fig, axes = plt.subplots(1, 4, figsize=(10.5, 3.6), sharey=False)
    for ax, name in zip(axes, PRIMARY_MODULES):
        data = []
        for q, col in (("Q1", "#4c78a8"), ("Q4", "#e45756")):
            vals = scores.loc[m.loc[m["quartile"] == q, "patient"], name].to_numpy(float)
            data.append(vals)
        bp = ax.boxplot(data, tick_labels=["Q1", "Q4"], patch_artist=True, widths=0.55,
                        medianprops={"color": "black"})
        for patch, col in zip(bp["boxes"], ["#4c78a8", "#e45756"]):
            patch.set_facecolor(col)
            patch.set_alpha(0.75)
        rng = np.random.default_rng(1)
        for i, (vals, cohort_q) in enumerate(zip(data, ["Q1", "Q4"]), start=1):
            patients = m.loc[m["quartile"] == cohort_q, "patient"]
            cols = [COHORT_COLORS[meta.set_index("patient").loc[p, "cohort"]] for p in patients]
            jitter = rng.uniform(-0.08, 0.08, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, c=cols, s=16, zorder=3, edgecolors="white", linewidths=0.3)
        ax.set_title(name, color=MODULE_COLORS[name])
        ax.set_ylabel("mean log2(TMM-CPM+1)" if name == "NHEJ" else "")
    fig.suptitle("Malignant module scores, within-cohort CLDN4 Q1 vs Q4", fontsize=11)
    fig.tight_layout()
    save_fig(fig, path)


def plot_scatter(scores: pd.DataFrame, meta: pd.DataFrame, path: Path) -> None:
    m = meta.loc[meta["patient"].isin(scores.index)].copy()
    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.3), sharex=True)
    for ax, name in zip(axes, PRIMARY_MODULES):
        for c in COHORTS:
            sub = m.loc[m["cohort"] == c]
            ax.scatter(
                sub["cldn4_pct"], scores.loc[sub["patient"], name],
                c=COHORT_COLORS[c], s=22, label=c, alpha=0.9,
            )
        ax.set_title(name, color=MODULE_COLORS[name], fontsize=10)
        ax.set_xlabel("CLDN4 %pos")
        ax.set_ylabel("score" if name == "NHEJ" else "")
    axes[0].legend(fontsize=6, frameon=False, loc="best")
    fig.suptitle("Malignant CLDN4 %pos vs module score (patient / sample)", fontsize=11)
    fig.tight_layout()
    save_fig(fig, path)


def plot_path(med: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    boxes = {
        "CLDN4": (1.2, 3.0),
        "NHEJ": (5.0, 4.6),
        "IFN": (8.4, 3.0),
    }
    for lab, (x, y) in boxes.items():
        ax.add_patch(plt.Rectangle((x - 1.05, y - 0.45), 2.1, 0.9, fill=True,
                                   facecolor="#f7f7f7", edgecolor="#222", lw=1.0, zorder=2))
        ax.text(x, y, lab, ha="center", va="center", fontsize=12, zorder=3)
    ax.annotate("", xy=(4.0, 4.55), xytext=(2.2, 3.35),
                arrowprops=dict(arrowstyle="->", color="#222", lw=1.3))
    ax.annotate("", xy=(7.4, 3.35), xytext=(6.0, 4.55),
                arrowprops=dict(arrowstyle="->", color="#222", lw=1.3))
    ax.annotate("", xy=(7.3, 3.0), xytext=(2.3, 3.0),
                arrowprops=dict(arrowstyle="->", color="#666", lw=1.1))
    ax.text(3.0, 4.35, f"a = {_signed(med['a'])}\np = {_fmt_p(med['a_p'])}", fontsize=8, ha="center")
    ax.text(7.15, 4.35, f"b = {_signed(med['b'])}\np = {_fmt_p(med['b_p'])}", fontsize=8, ha="center")
    ax.text(4.8, 2.45, f"c′ = {_signed(med['c_prime'])}   (total c = {_signed(med['c'])})", fontsize=8, ha="center")
    ax.text(
        5.0, 0.9,
        f"indirect a×b = {_signed(med['ab'])}\n"
        f"bootstrap 95% CI {_signed(med['ab_ci_lo'])} to {_signed(med['ab_ci_hi'])}\n"
        f"boot p = {_fmt_p(med['ab_boot_p'])}    Sobel p = {_fmt_p(med['sobel_p'])}\n"
        f"n = {med['n']}    proportion mediated = {_fmt(med['prop_mediated'], 2)}",
        ha="center", va="center", fontsize=8,
    )
    ax.set_title("CLDN4 → NHEJ → IFN   (cohort-adjusted, patient unit)", fontsize=11)
    save_fig(fig, path)


def plot_boot(ab: np.ndarray, med: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.hist(ab, bins=40, color="#4c78a8", alpha=0.85, edgecolor="white")
    ax.axvline(0, color="#333", lw=0.8)
    ax.axvline(med["ab"], color="#e45756", lw=1.4, label=f"a×b = {med['ab']:+.3f}")
    ax.axvline(med["ab_ci_lo"], color="#e45756", lw=0.8, ls="--")
    ax.axvline(med["ab_ci_hi"], color="#e45756", lw=0.8, ls="--", label="95% bootstrap CI")
    ax.set_xlabel("indirect effect (a×b)")
    ax.set_ylabel("bootstrap draws")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"Stratified bootstrap  B={med['boot_n']}")
    fig.tight_layout()
    save_fig(fig, path)


def plot_genes(genes: pd.DataFrame, path: Path) -> None:
    sub = genes[genes["module"].isin(["NHEJ", "cGAS-STING"])].copy()
    if sub.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharex=False)
    for ax, name in zip(axes, ["NHEJ", "cGAS-STING"]):
        g = sub[sub["module"] == name].sort_values("logFC")
        cols = ["#e45756" if v > 0 else "#4c78a8" for v in g["logFC"]]
        ax.barh(g["gene"], g["logFC"], color=cols, xerr=1.96 * g["se"], capsize=2, ecolor="#666")
        ax.axvline(0, color="#333", lw=0.7)
        ax.set_title(name)
        ax.set_xlabel("log2FC  Q4 − Q1")
    fig.suptitle("Module genes, stacked malignant Q4 vs Q1", fontsize=11)
    fig.tight_layout()
    save_fig(fig, path)


def plot_spearman(sp: pd.DataFrame, path: Path) -> None:
    sub = sp[(sp["level"] == "DL") & (sp["module"].isin(PRIMARY_MODULES))].copy()
    sub = sub.set_index("module").loc[PRIMARY_MODULES]
    fig, ax = plt.subplots(figsize=(5.8, 3.2))
    y = np.arange(len(PRIMARY_MODULES))[::-1]
    ax.axvline(0, color="#888", lw=0.8)
    ax.errorbar(
        sub["rho"], y,
        xerr=[sub["rho"] - sub["ci_lo"], sub["ci_hi"] - sub["rho"]],
        fmt="o", color="#222", ms=6, capsize=3,
    )
    for yi, name in zip(y, PRIMARY_MODULES):
        ax.scatter(sub.loc[name, "rho"], yi, c=MODULE_COLORS[name], s=40, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{n}   ρ={sub.loc[n, 'rho']:+.2f}  p={_fmt_p(sub.loc[n, 'p'])}" for n in PRIMARY_MODULES],
        fontsize=8,
    )
    ax.set_xlabel("DL pooled Spearman ρ  (CLDN4 %pos vs score)")
    ax.set_title("Four-cohort random-effects ρ")
    fig.tight_layout()
    save_fig(fig, path)


def write_finding(meta, contrasts, sp, meds, genes, membership, ifn_check) -> None:
    m_in = meta.loc[meta["in_count_matrix"]]
    q = m_in[m_in["quartile"].isin(["Q1", "Q4"])]
    n1 = int((q["quartile"] == "Q1").sum())
    n4 = int((q["quartile"] == "Q4").sum())
    n_all = int(meta.shape[0])
    n_mat = int(m_in.shape[0])

    def row_mod(module, cohort=COMBO, split="q4q1"):
        hit = contrasts[(contrasts["module"] == module) & (contrasts["cohort"] == cohort) & (contrasts["split"] == split)]
        if hit.empty:
            return None
        return hit.iloc[0]

    def md_primary():
        lines = [
            "| module | n | n_Q1 / n_Q4 | logFC (95% CI) | p | FDR | rank-biserial r (MW p) |",
            "|---|---:|---|---|---:|---:|---|",
        ]
        for name in PRIMARY_MODULES + ["NHEJ-core", "cGAS-STING-core", "proliferation"]:
            r = row_mod(name)
            if r is None:
                continue
            lines.append(
                f"| {name} | {int(r['n'])} | {int(r.n_q1)} / {int(r.n_q4)} | "
                f"{r.logFC:+.3f} ({r.ci_lo:+.3f}, {r.ci_hi:+.3f}) | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} | "
                f"{r.r_rb:+.3f} ({_fmt_p(r.p_mw)}) |"
            )
        return "\n".join(lines)

    def md_cont():
        lines = [
            "| module | n | logFC per SD CLDN4 %pos (95% CI) | p | FDR |",
            "|---|---:|---|---:|---:|",
        ]
        for name in PRIMARY_MODULES + ["NHEJ-core", "cGAS-STING-core", "proliferation"]:
            r = row_mod(name, split="continuous")
            if r is None:
                continue
            lines.append(
                f"| {name} | {int(r['n'])} | {r.logFC:+.3f} ({r.ci_lo:+.3f}, {r.ci_hi:+.3f}) | {_fmt_p(r.p)} | {_fmt_p(r.fdr)} |"
            )
        return "\n".join(lines)

    def md_cohort():
        lines = [
            "| cohort | module | n_Q1 / n_Q4 | logFC | p |",
            "|---|---|---|---:|---:|",
        ]
        for c in COHORTS:
            for name in PRIMARY_MODULES:
                r = row_mod(name, cohort=c)
                if r is None:
                    lines.append(f"| {c} | {name} | — | — | skipped (Q tail < 3) |")
                    continue
                lines.append(
                    f"| {c} | {name} | {int(r.n_q1)} / {int(r.n_q4)} | {r.logFC:+.3f} | {_fmt_p(r.p)} |"
                )
        return "\n".join(lines)

    def md_sp():
        lines = [
            "| module | k | N | ρ (95% CI) | p | I² |",
            "|---|---:|---:|---|---:|---:|",
        ]
        for name in PRIMARY_MODULES:
            r = sp[(sp["module"] == name) & (sp["level"] == "DL")].iloc[0]
            lines.append(
                f"| {name} | {int(r.k)} | {int(r.n)} | {r.rho:+.3f} ({r.ci_lo:+.3f}, {r.ci_hi:+.3f}) | {_fmt_p(r.p)} | {r.I2:.1f}% |"
            )
        lines += ["", "| cohort | module | n | ρ | p |", "|---|---|---:|---:|---:|"]
        for name in PRIMARY_MODULES:
            for c in COHORTS:
                r = sp[(sp["module"] == name) & (sp["cohort"] == c) & (sp["level"] == "cohort")].iloc[0]
                lines.append(f"| {c} | {name} | {int(r.n)} | {r.rho:+.3f} | {_fmt_p(r.p)} |")
        return "\n".join(lines)

    def md_med(rows):
        lines = [
            "| model | n | a (p) | b (p) | total c (p) | direct c′ (p) | a×b (boot 95% CI) | boot p | Sobel p | prop. |",
            "|---|---:|---|---|---|---|---|---:|---:|---:|",
        ]
        for r in rows:
            lines.append(
                f"| {r['label']} | {r['n']} | {_signed(r['a'])} ({_fmt_p(r['a_p'])}) | "
                f"{_signed(r['b'])} ({_fmt_p(r['b_p'])}) | {_signed(r['c'])} ({_fmt_p(r['c_p'])}) | "
                f"{_signed(r['c_prime'])} ({_fmt_p(r['c_prime_p'])}) | "
                f"{_signed(r['ab'])} ({_signed(r['ab_ci_lo'])}, {_signed(r['ab_ci_hi'])}) | "
                f"{_fmt_p(r['ab_boot_p'])} | {_fmt_p(r['sobel_p'])} | {_fmt(r['prop_mediated'], 2)} |"
            )
        return "\n".join(lines)

    def genes_present(module):
        sub = membership[membership["module"] == module]
        present = sub.loc[sub["in_logcpm"], "gene"].tolist()
        missing = sub.loc[~sub["in_logcpm"], "gene"].tolist()
        return present, missing

    nhej_p, nhej_m = genes_present("NHEJ")
    sting_p, sting_m = genes_present("cGAS-STING")
    ifn_p, ifn_m = genes_present("IFN")
    apm_p, apm_m = genes_present("APM")

    def md_genes(module):
        g = genes[(genes["module"] == module)].sort_values("p")
        lines = [
            "| gene | logFC | p | FDR |",
            "|---|---:|---:|---:|",
        ]
        for r in g.itertuples():
            lines.append(f"| {r.gene} | {r.logFC:+.3f} | {_fmt_p(r.p)} | {_fmt_p(r.fdr_within_module)} |")
        return "\n".join(lines)

    ifn_q = row_mod("IFN")
    apm_q = row_mod("APM")
    nhej_q = row_mod("NHEJ")
    sting_q = row_mod("cGAS-STING")
    sting_core = row_mod("cGAS-STING-core")
    primary = next(r for r in meds if r["label"].startswith("primary"))

    def gene_hit(module, gene):
        hit = genes[(genes["module"] == module) & (genes["gene"] == gene)]
        if hit.empty:
            return None
        return hit.iloc[0]

    cgas_g = gene_hit("cGAS-STING", "CGAS")
    sting_g = gene_hit("cGAS-STING", "STING1")

    text = f"""# Concordant-4 malignant CLDN4: NHEJ, cGAS-STING, IFN/APM, mediation

ADDITIVE. **CLDN4-only.** Not a re-derivation of the locked T/NK result
(malignant CLDN4 %pos vs T/NK ρ = −0.53, N = 65). Not a dual-high
TACSTD2×CLDN4 score. Cohorts are the four that already agree:
**GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not GSE148071,
GSE127465, GSE207422, GSE154826, or CD45+/T-only extracts.

The locked tumor-cell result this sits next to: malignant CLDN4-high
patients have lower own IFN and MHC-I/APM. IFN and APM are recomputed
here on the same matrices as a calibration, then NHEJ and cGAS-STING
are added, and the path CLDN4 → NHEJ → IFN is tested. This is an
observational patient-level path. It does not show that CLDN4 causes
the IFN change.

## Unit and split

Patient / donor / sample is the unit. Do not quote cell counts as n.

| Piece | Choice |
|---|---|
| Exposure | malignant CLDN4 %pos, within-cohort quartiles (same labels as PR #503) |
| High vs low | Q4 vs Q1, mid quartiles unused in the binary contrast |
| Scores | mean log2(TMM-CPM+1) of module genes present after the locked count filter |
| Aliases | MB21D1→CGAS, TMEM173→STING1, MRE11A→MRE11, placed on the locked TMM factors. IFN/APM symbols are not relabeled |
| Model | OLS ~ cohort + CLDN4_Q4 (reference cohort GSE123902) |
| Matrix | marker-malignant UMI-sum (GSE123902, GSE189357); author-malignant UMI-sum (GSE131907, GSE205335) |
| In the count matrix | {n_mat} of {n_all} units (P4001 is Q1 on the %pos vector and is absent from the UMI-sum) |
| Stacked Q4 vs Q1 | n = {n1 + n4} ({n1} / {n4}) |

GSE189357 Q4 has 2 patients, so its within-cohort binary contrast is skipped.
Those two patients stay in the stacked model.

IFN calibration against the locked family score (PR #503, logFC −0.584, p = 0.00243):
this run logFC = {ifn_q.logFC:+.3f}, p = {_fmt_p(ifn_q.p)}. APM locked −0.779; this run {apm_q.logFC:+.3f}.

## Modules

| module | definition | genes in score | held out of the log matrix |
|---|---|---:|---|
| NHEJ | KEGG hsa03450 `KEGG_NON_HOMOLOGOUS_END_JOINING` | {len(nhej_p)} | {", ".join(nhej_m) or "—"} |
| NHEJ-core | catalytic c-NHEJ (sensitivity) | {len(genes_present("NHEJ-core")[0])} | {", ".join(genes_present("NHEJ-core")[1]) or "—"} |
| cGAS-STING | sensors/adapters; IFN-hallmark and NHEJ genes removed | {len(sting_p)} | {", ".join(sting_m) or "—"} |
| IFN | Hallmark IFNα ∪ IFNγ (same lists as PR #503) | {len(ifn_p)} | {len(ifn_m)} genes below the count filter |
| APM | custom MHC-I / antigen-processing panel (PR #503) | {len(apm_p)} | {", ".join(apm_m) or "—"} |
| proliferation | MKI67, PCNA, TOP2A, MCM2, CDK1, CCNB1, BIRC5, UBE2C | {len(genes_present("proliferation")[0])} | covariate only |

NHEJ genes used: {", ".join(nhej_p)}.
MRE11 is the sum of the MRE11 / MRE11A symbol, whichever the GEO build used.
DNTT is absent (not expressed in these lung malignant sums).

cGAS-STING genes used: {", ".join(sting_p)}.
CGAS is MB21D1 in GSE123902 and GSE131907, and CGAS in GSE205335 and GSE189357.
STING1 is TMEM173 in all four matrices.
Held out of cGAS-STING because they sit in Hallmark IFN: IRF7, ZBP1, TRIM21, SAMHD1.
Held out because they sit in KEGG NHEJ: PRKDC, XRCC5, XRCC6, MRE11.
Negative regulators TREX1 and NLRP4 are not shared across the four builds, so that pair was not scored.

## 1. CLDN4-high vs low (stacked Q4 vs Q1)

Positive logFC = higher in CLDN4-high. FDR is Benjamini–Hochberg across the seven scores fit in that model.

{md_primary()}

Continuous CLDN4 %pos (global z, same n as the locked continuous family model):

{md_cont()}

### Within cohort

{md_cohort()}

### Spearman, DerSimonian–Laird on Fisher-z

{md_sp()}

### NHEJ genes (stacked Q4 vs Q1)

{md_genes("NHEJ")}

### cGAS-STING genes (stacked Q4 vs Q1)

{md_genes("cGAS-STING")}

## 2. Mediation CLDN4 → NHEJ → IFN

Exposure = malignant CLDN4 %pos, global z on the units in the model.
Mediator = NHEJ module score. Outcome = IFN module score.
Cohort fixed effects in every equation. Indirect effect = a×b.
Uncertainty is a cohort-stratified bootstrap (B = {BOOT_B}, seed {BOOT_SEED});
the percentile interval is primary. Sobel is the normal approximation.
Proportion mediated = a×b / c, reported only as a descriptive ratio.

Paths: a is CLDN4 → NHEJ, b is NHEJ → IFN given CLDN4, c is the total
CLDN4 → IFN effect, c′ is the direct effect given NHEJ.

{md_med(meds)}

Primary read. KEGG NHEJ does not differ by CLDN4 quartile
(logFC {nhej_q.logFC:+.3f}, p = {_fmt_p(nhej_q.p)}). No NHEJ gene has
FDR < 0.05. The indirect effect is a×b = {_signed(primary["ab"])}
(95% CI {_signed(primary["ab_ci_lo"])} to {_signed(primary["ab_ci_hi"])},
boot p = {_fmt_p(primary["ab_boot_p"])}). Path a is flat, so these units
do not support NHEJ as the mediator of the CLDN4–IFN association.
Path b is a partial association of the NHEJ score with IFN given CLDN4;
it is not the mediation result. Catalytic NHEJ-core and the
proliferation-adjusted model give the same null indirect effect.
The leave-one-out row that drops GSE123902 has an indirect-effect
interval spanning tens of log units. That row is unstable and is not
a result. The other three leave-one-out intervals still cover 0.

cGAS-STING module, stacked Q4 vs Q1: logFC {sting_q.logFC:+.3f}
(p = {_fmt_p(sting_q.p)}, FDR = {_fmt_p(sting_q.fdr)}). The continuous
model and the four-cohort Spearman pool do not separate this module
from 0. The four-gene core (CGAS, STING1, TBK1, IRF3) is
logFC {sting_core.logFC:+.3f} (p = {_fmt_p(sting_core.p)}).
Gene-level, the module dip is CGAS
(logFC {cgas_g.logFC:+.3f}, p = {_fmt_p(cgas_g.p)}, FDR = {_fmt_p(cgas_g.fdr_within_module)})
and IFI16. STING1 is not in that direction
(logFC {sting_g.logFC:+.3f}, p = {_fmt_p(sting_g.p)}).
Parallel mediation through the 11-gene module has a bootstrap interval
that includes 0.

IFN and APM reproduce the locked direction (IFN logFC {ifn_q.logFC:+.3f},
p = {_fmt_p(ifn_q.p)}; APM logFC {apm_q.logFC:+.3f}, p = {_fmt_p(apm_q.p)}).

Sensitivities in the same table: catalytic NHEJ-core instead of KEGG NHEJ;
proliferation score as a covariate; cGAS-STING and the four-gene core as
parallel mediators; APM as a parallel outcome; pseudobulk CLDN4 expression
instead of %pos; leave-one-cohort-out.

## What this is not

- Not a new T/NK estimate and not a re-audit of ρ = −0.53.
- Not cell-level high vs low inside a tumor. The public inputs are
  malignant UMI-sums, one column per donor/sample/patient.
- Not evidence that CLDN4, NHEJ, or cGAS-STING causes the IFN change.
- Not a serial CLDN4 → NHEJ → cGAS-STING → IFN structural model.
- Not GSE148071 / GSE127465 / GSE207422 / GSE154826.
- Not a dual-high TACSTD2 score.
- Genome-wide FDR is not the claim. The pre-specified mediation is the
  KEGG model. The specification search is a sensitivity analysis; its
  smallest p-value is not a confirmatory test.
- Not a rebranding of Hallmark DNA repair as NHEJ.

## Files

- `tables/module_q4q1.tsv` — high vs low and continuous module contrasts
- `tables/module_spearman.tsv` — cohort ρ and DL pool
- `tables/mediation.tsv` — CLDN4 → NHEJ → IFN and the sensitivities
- `tables/module_genes.tsv` — gene-level Q4 vs Q1
- `tables/patient_module_scores.tsv` — one row per unit
- `tables/module_membership.tsv` — which genes entered each score
- `figures/` — forests, boxes, scatters, path, bootstrap

Reproduce:

```bash
python3 methods/concordant4_cldn4_nhej_sting_ifn/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)


def _splice_search_section() -> None:
    """Keep the specification-search writeup if analyze.py regenerates FINDING.md."""
    section_path = HERE / "SEARCH_SECTION.md"
    finding = HERE / "FINDING.md"
    if not section_path.exists() or not finding.exists():
        return
    text = finding.read_text()
    if "## 3. Specification search" in text:
        return
    marker = "## What this is not"
    section = section_path.read_text().rstrip() + "\n\n"
    if marker not in text:
        return
    finding.write_text(text.replace(marker, section + marker, 1))


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    ifn, apm = load_ifn_apm()
    genesets = module_genes(ifn, apm)
    # guard: cGAS-STING must not contain IFN, APM, or NHEJ genes
    ifn_set, apm_set, nhej_set = set(ifn), set(apm), set(NHEJ_KEGG)
    overlap = sorted(set(CGAS_STING) & (ifn_set | apm_set | nhej_set))
    if overlap:
        raise SystemExit(f"cGAS-STING overlaps IFN/APM/NHEJ: {overlap}")

    units = load_units()
    meta = four_meta(units)
    print("UNITS", {c: units[c]["n"] for c in COHORTS}, "matrix", int(meta["in_count_matrix"].sum()))

    parts = load_parts()
    combined = load_combined(parts)
    m_all = meta.loc[meta["patient"].isin(combined.columns)].copy()
    cts = filter_genes(combined.loc[:, m_all["patient"]])
    factors = tmm_norm_factors(cts)
    lc = add_alias_genes(log_cpm(cts, factors), cts, factors, parts)
    # CLDN4 expression on the same scale, for the exposure sensitivity
    if "CLDN4" not in lc.index:
        raise SystemExit("CLDN4 missing from logcpm")
    scores, membership = score_modules(lc, genesets)
    scores["CLDN4"] = lc.loc["CLDN4"].astype(float)
    scores = scores.loc[m_all["patient"]]

    contrast_rows = []
    for cohort in (None, *COHORTS):
        for split in ("q4q1", "continuous"):
            part = module_contrast(scores.drop(columns=["CLDN4"]), meta, split, cohort)
            if not part.empty:
                contrast_rows.append(part)
    contrasts = pd.concat(contrast_rows, ignore_index=True)
    # FDR only within the stacked primary table; recompute so per-cohort rows keep their own BH
    contrasts.to_csv(TABLES / "module_q4q1.tsv", sep="\t", index=False)

    sp = spearman_meta(scores.drop(columns=["CLDN4"]), meta)
    sp.to_csv(TABLES / "module_spearman.tsv", sep="\t", index=False)

    genes = gene_level_q4(lc, meta, genesets)
    genes.to_csv(TABLES / "module_genes.tsv", sep="\t", index=False)
    membership.to_csv(TABLES / "module_membership.tsv", sep="\t", index=False)

    patient = meta.merge(scores.reset_index().rename(columns={"index": "patient"}), on="patient", how="left")
    patient.to_csv(TABLES / "patient_module_scores.tsv", sep="\t", index=False)

    # calibration print
    cal = contrasts[(contrasts["cohort"] == COMBO) & (contrasts["split"] == "q4q1") & (contrasts["module"] == "IFN")].iloc[0]
    print(f"IFN Q4 logFC {cal.logFC:+.4f} p {cal.p:.4g}  (locked −0.584)")

    meds = []
    boots = {}
    specs = [
        ("primary: CLDN4% → NHEJ → IFN", "NHEJ", "IFN", "cldn4_pct", None),
        ("sensitivity: CLDN4% → NHEJ-core → IFN", "NHEJ-core", "IFN", "cldn4_pct", None),
        ("sensitivity: CLDN4% → NHEJ → IFN, given proliferation", "NHEJ", "IFN", "cldn4_pct", ["proliferation"]),
        ("parallel: CLDN4% → cGAS-STING → IFN", "cGAS-STING", "IFN", "cldn4_pct", None),
        ("sensitivity: CLDN4% → cGAS-STING-core → IFN", "cGAS-STING-core", "IFN", "cldn4_pct", None),
        ("secondary outcome: CLDN4% → NHEJ → APM", "NHEJ", "APM", "cldn4_pct", None),
        ("sensitivity: CLDN4 expr → NHEJ → IFN", "NHEJ", "IFN", "cldn4_expr", None),
    ]
    for label, mcol, ycol, xcol, extra in specs:
        print("mediation", label)
        point, ab = run_one_mediation(scores, meta, label, mcol, ycol, xcol, extra)
        meds.append(point)
        boots[label] = ab
        print(f"  ab={point['ab']:+.4f} CI {point['ab_ci_lo']:+.4f},{point['ab_ci_hi']:+.4f} p={point['ab_boot_p']:.3g}")

    # leave-one-cohort-out of the primary path
    for drop in COHORTS:
        keep = meta.loc[meta["cohort"] != drop].copy()
        # scores index is patient; _analysis_frame filters by meta patients
        label = f"LOO drop {drop}: CLDN4% → NHEJ → IFN"
        print("mediation", label)
        point, ab = run_one_mediation(scores, keep, label, "NHEJ", "IFN", "cldn4_pct", None)
        meds.append(point)

    med_df = pd.DataFrame([{k: v for k, v in r.items() if k != "ab_boot"} for r in meds])
    med_df.to_csv(TABLES / "mediation.tsv", sep="\t", index=False)

    plot_forest(contrasts, FIGS / "forest_module_q4q1")
    plot_boxes(scores, meta, FIGS / "box_module_q4q1")
    plot_scatter(scores, meta, FIGS / "scatter_cldn4_modules")
    plot_spearman(sp, FIGS / "forest_spearman_modules")
    primary = meds[0]
    plot_path(primary, FIGS / "mediation_path")
    plot_boot(boots[primary["label"]], primary, FIGS / "mediation_bootstrap")
    plot_genes(genes, FIGS / "gene_logfc_nhej_sting")

    write_finding(meta, contrasts, sp, meds, genes, membership, cal)
    _splice_search_section()
    print("WROTE", HERE / "FINDING.md")
    show = contrasts[(contrasts["cohort"] == COMBO) & (contrasts["split"] == "q4q1")]
    print(show[["module", "n", "n_q1", "n_q4", "logFC", "p", "r_rb"]].to_string(index=False))


if __name__ == "__main__":
    main()

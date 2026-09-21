#!/usr/bin/env python3
"""Concordant-4 malignant CLDN4-high vs low: decoupler DoRothEA TF activity.

Patient / donor / sample is the unit. Quartile labels are the PR #503
within-cohort malignant CLDN4 %pos labels (not recomputed).

Primary network: DoRothEA A+B+C (decoupler default confidence weights).
Primary statistic: univariate linear model (ULM) t-value per sample,
then OLS of that activity on cohort + CLDN4 Q4.

Also scored on the same matrix: MLM, VIPER (aREA), and weighted mean.
NLRC5 is not in DoRothEA or CollecTRI; its transcript is reported as a gene,
not as a TF activity.
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

import decoupler as dc

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
RES = HERE / "resources"
OUT = HERE / "results"
FIGS = HERE / "figures"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
TMIN = 5

# Expected sign of (CLDN4-high − CLDN4-low) if the observation matched
# CLDN4-loss → IFN/MHC up. Junction TFs have no KD direction.
FOCUS = [
    ("IFN", "STAT1", "down"),
    ("IFN", "STAT2", "down"),
    ("IFN", "IRF1", "down"),
    ("IFN", "IRF2", "down"),
    ("IFN", "IRF3", "down"),
    ("IFN", "IRF7", "down"),
    ("IFN", "IRF8", "down"),
    ("IFN", "IRF9", "down"),
    ("MHC", "RFX5", "down"),
    ("junction", "GRHL2", "none"),
    ("junction", "ELF3", "none"),
    ("junction", "KLF4", "none"),
    ("junction", "TFAP2A", "none"),
]
# Context rows, not the KD test.
CONTEXT = [("context", "STAT3", "none"), ("context", "TP63", "none")]
COLLECTRI_FOCUS = [
    "STAT1", "STAT2", "IRF1", "IRF9", "NLRC5", "CIITA",
    "RFX5", "RFXANK", "RFXAP", "GRHL2", "ELF3", "KLF4", "TFAP2A",
]
GENE_READOUTS = ["NLRC5", "CIITA", "STAT1", "IRF1", "RFX5", "CLDN4", "HLA-A", "B2M", "TAP1", "TAP2"]
METHODS = ["ulm", "mlm", "viper", "wmean"]
PROG_COLOR = {
    "IFN": "#b2182b",
    "MHC": "#2166ac",
    "junction": "#1b7837",
    "context": "#666666",
    "gene": "#7a4e00",
    "geneset": "#4d4d4d",
}


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


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
        m, a, w = m[ok], a[ok], w[ok]
        if len(m) < 50:
            factors[col] = 1.0
            continue
        lo_m, hi_m = np.quantile(m, [0.30, 0.70])
        lo_a, hi_a = np.quantile(a, [0.05, 0.95])
        trim = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if int(trim.sum()) < 20:
            factors[col] = 1.0
            continue
        factors[col] = float(2 ** np.average(m[trim], weights=1.0 / w[trim]))
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


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


def ols_coef(y: pd.Series, design: pd.DataFrame, coef: str) -> dict:
    X = design.reindex(y.index).astype(float)
    keep = y.notna() & ~X.isna().any(axis=1)
    yv = y.loc[keep].astype(float).to_numpy()
    Xv = X.loc[keep].to_numpy()
    n, p = Xv.shape
    df = n - p
    empty = {"n": int(n), "coef": np.nan, "se": np.nan, "t": np.nan, "p": np.nan, "df": int(df)}
    if df < 1 or n < p + 1:
        return empty
    xtx = Xv.T @ Xv
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ Xv.T @ yv
    resid = yv - Xv @ beta
    sigma2 = float(np.sum(resid**2) / df)
    j = list(X.columns).index(coef)
    se = float(np.sqrt(max(sigma2 * xtx_inv[j, j], 0.0)))
    est = float(beta[j])
    t = est / se if se > 0 else np.nan
    pv = float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else np.nan
    return {"n": int(n), "coef": est, "se": se, "t": float(t) if np.isfinite(t) else np.nan, "p": pv, "df": int(df)}


def design_q4(meta: pd.DataFrame, with_cohort: bool) -> pd.DataFrame:
    m = meta.set_index("patient")
    d = pd.DataFrame(index=m.index)
    d["Intercept"] = 1.0
    d["CLDN4_Q4"] = (m["quartile"] == "Q4").astype(float)
    if with_cohort:
        for c in COHORTS:
            if c == REF:
                continue
            d[f"cohort_{c}"] = (m["cohort"] == c).astype(float)
    return d


def design_cont(meta: pd.DataFrame, with_cohort: bool) -> pd.DataFrame:
    m = meta.set_index("patient")
    z = m["cldn4_pct"].astype(float)
    z = (z - z.mean()) / z.std(ddof=1)
    d = pd.DataFrame(index=m.index)
    d["Intercept"] = 1.0
    d["CLDN4_pct_z"] = z
    if with_cohort:
        for c in COHORTS:
            if c == REF:
                continue
            d[f"cohort_{c}"] = (m["cohort"] == c).astype(float)
    return d


def load_net() -> pd.DataFrame:
    net = pd.read_csv(RES / "dorothea_hs_ABC.tsv", sep="\t")
    net["source"] = net["source"].astype(str).str.upper()
    net["target"] = net["target"].astype(str).str.upper()
    net["weight"] = net["weight"].astype(float)
    return net[["source", "target", "weight", "confidence"]]


def load_collectri() -> pd.DataFrame:
    net = pd.read_csv(RES / "collectri.tsv", sep="\t")
    net["source"] = net["source"].astype(str).str.upper()
    net["target"] = net["target"].astype(str).str.upper()
    net["weight"] = net["weight"].astype(float)
    return net[["source", "target", "weight"]]


def load_meta() -> pd.DataFrame:
    meta = pd.read_csv(DATA / "sample_inventory.tsv", sep="\t")
    meta = meta.loc[meta["in_count_matrix"].astype(str).isin(["True", "true", "1"])].copy()
    meta["patient"] = meta["patient"].astype(str)
    meta["cohort"] = meta["cohort"].astype(str)
    meta["quartile"] = meta["quartile"].astype(str)
    meta["cldn4_pct"] = meta["cldn4_pct"].astype(float)
    return meta


def load_counts(meta: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for cohort in COHORTS:
        cts = read_counts(DATA / f"{cohort}_malignant_counts.tsv.gz")
        keep = [c for c in cts.columns if c in set(meta.loc[meta["cohort"] == cohort, "patient"])]
        parts.append(cts.loc[:, keep])
    genes = sorted(set.intersection(*[set(p.index) for p in parts]))
    mat = pd.concat([p.reindex(genes).fillna(0.0) for p in parts], axis=1)
    mat = mat.loc[:, ~mat.columns.duplicated()]
    return mat


def normalize(counts: pd.DataFrame) -> pd.DataFrame:
    cts = filter_genes(counts)
    return log_cpm(cts, tmm_norm_factors(cts))


def score_dorothea(logcpm: pd.DataFrame, net: pd.DataFrame) -> dict[str, pd.DataFrame]:
    data = logcpm.T.astype(float)
    ulm, _ = dc.mt.ulm(data, net, tmin=TMIN)
    mlm, _ = dc.mt.mlm(data, net, tmin=TMIN)
    viper, _ = dc.mt.viper(data, net, tmin=TMIN, pleiotropy=True)
    wmean, _ = dc.mt.waggr(data, net, tmin=TMIN, fun="wmean", times=1, seed=42)
    return {"ulm": ulm, "mlm": mlm, "viper": viper, "wmean": wmean}


def n_targets(net: pd.DataFrame, genes: pd.Index) -> pd.Series:
    sub = net.loc[net["target"].isin(set(genes))]
    return sub.groupby("source")["target"].nunique()


def test_scores(scores: pd.DataFrame, meta: pd.DataFrame, method: str, nt: pd.Series) -> pd.DataFrame:
    rows = []
    q = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    design_q = design_q4(q, with_cohort=True)
    design_c = design_cont(meta, with_cohort=True)
    for tf in scores.columns:
        y_all = scores[tf]
        yq = y_all.reindex(q["patient"])
        fit = ols_coef(yq, design_q, "CLDN4_Q4")
        fit_c = ols_coef(y_all.reindex(meta["patient"]), design_c, "CLDN4_pct_z")
        # unadjusted means
        hi = yq.reindex(q.loc[q["quartile"] == "Q4", "patient"]).astype(float)
        lo = yq.reindex(q.loc[q["quartile"] == "Q1", "patient"]).astype(float)
        rows.append(
            {
                "method": method,
                "tf": tf,
                "n_targets": int(nt.get(tf, 0)),
                "n": fit["n"],
                "n_q1": int((q["quartile"] == "Q1").sum()),
                "n_q4": int((q["quartile"] == "Q4").sum()),
                "mean_q1": float(lo.mean()),
                "mean_q4": float(hi.mean()),
                "delta": float(hi.mean() - lo.mean()),
                "coef": fit["coef"],
                "se": fit["se"],
                "t": fit["t"],
                "p": fit["p"],
                "df": fit["df"],
                "n_continuous": fit_c["n"],
                "coef_continuous": fit_c["coef"],
                "se_continuous": fit_c["se"],
                "p_continuous": fit_c["p"],
            }
        )
    out = pd.DataFrame(rows)
    out["fdr"] = bh(out["p"].to_numpy())
    out["fdr_continuous"] = bh(out["p_continuous"].to_numpy())
    return out


def per_cohort_tests(scores: pd.DataFrame, meta: pd.DataFrame, tfs: list[str]) -> pd.DataFrame:
    rows = []
    for cohort in COHORTS:
        sub = meta.loc[(meta["cohort"] == cohort) & meta["quartile"].isin(["Q1", "Q4"])].copy()
        n1 = int((sub["quartile"] == "Q1").sum())
        n4 = int((sub["quartile"] == "Q4").sum())
        for tf in tfs:
            if tf not in scores.columns:
                rows.append(
                    {
                        "cohort": cohort, "tf": tf, "n_q1": n1, "n_q4": n4,
                        "coef": np.nan, "p": np.nan, "note": "not scored",
                    }
                )
                continue
            if n1 < 3 or n4 < 3:
                rows.append(
                    {
                        "cohort": cohort, "tf": tf, "n_q1": n1, "n_q4": n4,
                        "coef": np.nan, "p": np.nan, "note": "tail n<3; skipped",
                    }
                )
                continue
            y = scores[tf].reindex(sub["patient"])
            fit = ols_coef(y, design_q4(sub, with_cohort=False), "CLDN4_Q4")
            rows.append(
                {
                    "cohort": cohort, "tf": tf, "n_q1": n1, "n_q4": n4,
                    "coef": fit["coef"], "se": fit["se"], "p": fit["p"], "note": "",
                }
            )
    return pd.DataFrame(rows)


def loo_tests(scores: pd.DataFrame, meta: pd.DataFrame, tfs: list[str]) -> pd.DataFrame:
    rows = []
    q = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    for drop in COHORTS:
        sub = q.loc[q["cohort"] != drop].copy()
        design = design_q4(sub, with_cohort=sub["cohort"].nunique() > 1)
        for tf in tfs:
            if tf not in scores.columns:
                continue
            fit = ols_coef(scores[tf].reindex(sub["patient"]), design, "CLDN4_Q4")
            rows.append(
                {
                    "dropped": drop,
                    "tf": tf,
                    "n": fit["n"],
                    "n_q1": int((sub["quartile"] == "Q1").sum()),
                    "n_q4": int((sub["quartile"] == "Q4").sum()),
                    "coef": fit["coef"],
                    "se": fit["se"],
                    "p": fit["p"],
                }
            )
    return pd.DataFrame(rows)


def gene_program_tests(logcpm: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    ifn = set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]) | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
    mhc = set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"])
    tj = set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    tj.discard("CLDN4")
    programs = {
        "IFN_hallmark": ifn,
        "MHC_I_APM": mhc,
        "TJ_no_CLDN4": tj,
    }
    # Hallmark sets from the same a8 file when present.
    for key, name in [
        ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "OXPHOS"),
        ("HALLMARK_APICAL_JUNCTION", "APICAL_JUNCTION"),
    ]:
        if key in sets:
            programs[name] = set(sets[key])
    q = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    design = design_q4(q, with_cohort=True)
    rows = []
    for name, genes in programs.items():
        present = [g for g in sorted(genes) if g in logcpm.index]
        if len(present) < 3:
            continue
        score = logcpm.loc[present, q["patient"]].astype(float).mean(axis=0)
        fit = ols_coef(score, design, "CLDN4_Q4")
        rows.append(
            {
                "program": name,
                "n_genes": len(present),
                "n": fit["n"],
                "coef": fit["coef"],
                "se": fit["se"],
                "p": fit["p"],
                "mean_q1": float(score.reindex(q.loc[q["quartile"] == "Q1", "patient"]).mean()),
                "mean_q4": float(score.reindex(q.loc[q["quartile"] == "Q4", "patient"]).mean()),
            }
        )
    # single genes
    for gene in GENE_READOUTS:
        if gene not in logcpm.index:
            rows.append({"program": f"gene:{gene}", "n_genes": 0, "n": 0, "coef": np.nan, "se": np.nan, "p": np.nan, "note": "absent"})
            continue
        y = logcpm.loc[gene, q["patient"]].astype(float)
        fit = ols_coef(y, design, "CLDN4_Q4")
        rows.append(
            {
                "program": f"gene:{gene}",
                "n_genes": 1,
                "n": fit["n"],
                "coef": fit["coef"],
                "se": fit["se"],
                "p": fit["p"],
                "mean_q1": float(y.reindex(q.loc[q["quartile"] == "Q1", "patient"]).mean()),
                "mean_q4": float(y.reindex(q.loc[q["quartile"] == "Q4", "patient"]).mean()),
            }
        )
    return pd.DataFrame(rows)


def gene_logfc(logcpm: pd.DataFrame, meta: pd.DataFrame) -> pd.Series:
    q = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    design = design_q4(q, with_cohort=True)
    X = design.reindex(q["patient"]).astype(float)
    Y = logcpm.loc[:, q["patient"]].astype(float).to_numpy()
    Xv = X.to_numpy()
    n, p = Xv.shape
    df = n - p
    xtx_inv = np.linalg.pinv(Xv.T @ Xv)
    beta = xtx_inv @ Xv.T @ Y.T
    j = list(X.columns).index("CLDN4_Q4")
    return pd.Series(beta[j], index=logcpm.index, name="logFC")


def stat1_split(net: pd.DataFrame, logfc: pd.Series, ifn_genes: set[str]) -> pd.DataFrame:
    sub = net.loc[net["source"] == "STAT1"].copy()
    sub = sub.loc[sub["target"].isin(logfc.index)]
    sub["logFC"] = sub["target"].map(logfc)
    sub["ifn"] = sub["target"].isin(ifn_genes)
    rows = []
    for label, mask in [
        ("all_STAT1_targets", np.ones(len(sub), dtype=bool)),
        ("positive_weight", sub["weight"] > 0),
        ("negative_weight", sub["weight"] < 0),
        ("IFN_hallmark_overlap", sub["ifn"].to_numpy()),
        ("not_IFN_hallmark", ~sub["ifn"].to_numpy()),
        ("positive_and_IFN", ((sub["weight"] > 0) & sub["ifn"]).to_numpy()),
        ("positive_not_IFN", ((sub["weight"] > 0) & ~sub["ifn"]).to_numpy()),
    ]:
        sl = sub.loc[mask]
        if sl.empty:
            continue
        rows.append(
            {
                "slice": label,
                "n_targets": int(len(sl)),
                "median_logFC": float(sl["logFC"].median()),
                "mean_logFC": float(sl["logFC"].mean()),
                "frac_logFC_neg": float((sl["logFC"] < 0).mean()),
            }
        )
    return pd.DataFrame(rows)


def within_cohort_ulm(counts: pd.DataFrame, meta: pd.DataFrame, net: pd.DataFrame, tfs: list[str]) -> pd.DataFrame:
    rows = []
    for cohort in COHORTS:
        ids = meta.loc[meta["cohort"] == cohort, "patient"]
        sub_counts = counts.loc[:, [c for c in counts.columns if c in set(ids)]]
        lc = normalize(sub_counts)
        scores, _ = dc.mt.ulm(lc.T.astype(float), net, tmin=TMIN)
        sub = meta.loc[(meta["cohort"] == cohort) & meta["quartile"].isin(["Q1", "Q4"])].copy()
        n1 = int((sub["quartile"] == "Q1").sum())
        n4 = int((sub["quartile"] == "Q4").sum())
        for tf in tfs:
            if n1 < 3 or n4 < 3 or tf not in scores.columns:
                rows.append({"cohort": cohort, "tf": tf, "n_q1": n1, "n_q4": n4, "coef": np.nan, "p": np.nan, "note": "skipped"})
                continue
            fit = ols_coef(scores[tf].reindex(sub["patient"]), design_q4(sub, False), "CLDN4_Q4")
            rows.append({"cohort": cohort, "tf": tf, "n_q1": n1, "n_q4": n4, "coef": fit["coef"], "se": fit["se"], "p": fit["p"], "note": "within-cohort TMM+ULM"})
    return pd.DataFrame(rows)


def matches_kd(coef: float, expected: str) -> str:
    if expected == "none" or not np.isfinite(coef):
        return "not a KD prediction" if expected == "none" else "not scored"
    if coef < 0 and expected == "down":
        return "matches KD (lower in CLDN4-high)"
    if coef > 0 and expected == "down":
        return "opposite of KD (higher in CLDN4-high)"
    return "flat"


def forest(focus_tests: pd.DataFrame, path: Path) -> None:
    sub = focus_tests.loc[focus_tests["method"] == "ulm"].copy()
    order = [tf for _, tf, _ in FOCUS + CONTEXT]
    sub = sub.loc[np.isfinite(sub["coef"].astype(float))].copy()
    sub["tf"] = pd.Categorical(sub["tf"], categories=order, ordered=True)
    sub = sub.dropna(subset=["tf"]).sort_values("tf", ascending=False)
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    y = np.arange(len(sub))
    ax.axvline(0, color="#444", lw=0.8)
    ax.errorbar(sub["coef"], y, xerr=1.96 * sub["se"], fmt="none", ecolor="#888", elinewidth=1.1, capsize=2.4)
    ax.scatter(sub["coef"], y, c=[PROG_COLOR.get(p, "#333") for p in sub["program"]], s=42, zorder=3)
    labels = []
    for r in sub.itertuples():
        mark = {"down": "KD expects −", "none": "no KD sign"}[r.expected]
        labels.append(f"{r.program}  {r.tf}   {mark}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("ULM activity, OLS coefficient (CLDN4 Q4 − Q1)")
    ax.set_title("DoRothEA A+B+C on malignant pseudobulk\nconcordant-4, cohort-adjusted, n=18 low / 16 high")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def method_heatmap(focus_tests: pd.DataFrame, path: Path) -> None:
    order = [tf for _, tf, _ in FOCUS + CONTEXT]
    sub = focus_tests.loc[focus_tests["tf"].isin(order)].copy()
    mat = sub.pivot(index="tf", columns="method", values="coef").reindex(index=order, columns=METHODS)
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    vals = mat.to_numpy(dtype=float)
    lim = np.nanmax(np.abs(vals))
    lim = 1.0 if not np.isfinite(lim) or lim == 0 else lim
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels(METHODS)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7, color="black")
            else:
                ax.text(j, i, "NA", ha="center", va="center", fontsize=7, color="#666")
    ax.set_title("Same patients, four decoupler statistics\npositive = higher activity in CLDN4-high")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="OLS coef")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def box_activities(scores: pd.DataFrame, meta: pd.DataFrame, path: Path) -> None:
    tfs = ["STAT1", "IRF1", "RFX5", "GRHL2", "ELF3", "KLF4"]
    q = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    fig, axes = plt.subplots(1, len(tfs), figsize=(12.6, 3.6), sharey=False)
    for ax, tf in zip(axes, tfs):
        if tf not in scores.columns:
            ax.set_visible(False)
            continue
        data = [
            scores.loc[q.loc[q["quartile"] == lab, "patient"], tf].astype(float).dropna()
            for lab in ("Q1", "Q4")
        ]
        bp = ax.boxplot(data, tick_labels=["low", "high"], patch_artist=True, widths=0.6, showfliers=False)
        for patch, color in zip(bp["boxes"], ["#4c78a8", "#e45756"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        rng = np.random.default_rng(0)
        for i, s in enumerate(data, start=1):
            ax.scatter(i + rng.uniform(-0.08, 0.08, len(s)), s, s=10, c="#222", alpha=0.7, zorder=3)
        ax.set_title(tf, fontsize=10)
        ax.set_xlabel("CLDN4")
    axes[0].set_ylabel("ULM activity")
    fig.suptitle("Malignant pseudobulk ULM (DoRothEA), Q1 vs Q4", fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def stat1_bars(split: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    y = np.arange(len(split))[::-1]
    ax.barh(y, split["mean_logFC"], color=["#b2182b" if v < 0 else "#2166ac" for v in split["mean_logFC"]], height=0.7)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.slice}  (n={r.n_targets})" for r in split.itertuples()], fontsize=8)
    ax.set_xlabel("mean gene log2FC (CLDN4 Q4 − Q1), same matrix as ULM")
    ax.set_title("STAT1 DoRothEA targets are not the IFN program")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def cohort_forest(cohort_df: pd.DataFrame, path: Path) -> None:
    tfs = ["STAT1", "IRF1", "RFX5", "GRHL2"]
    fig, axes = plt.subplots(1, len(tfs), figsize=(11.2, 3.8), sharey=True)
    for ax, tf in zip(axes, tfs):
        sub = cohort_df.loc[cohort_df["tf"] == tf].copy()
        sub["cohort"] = pd.Categorical(sub["cohort"], categories=COHORTS, ordered=True)
        sub = sub.sort_values("cohort", ascending=False)
        y = np.arange(len(sub))
        ax.axvline(0, color="#444", lw=0.7)
        ok = np.isfinite(sub["coef"].astype(float))
        ax.errorbar(
            sub.loc[ok, "coef"],
            y[ok.to_numpy()],
            xerr=1.96 * sub.loc[ok, "se"],
            fmt="o",
            color="#333",
            ecolor="#888",
            elinewidth=1.0,
            capsize=2,
            ms=5,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{r.cohort}\n{int(r.n_q1)}/{int(r.n_q4)}" for r in sub.itertuples()],
            fontsize=7,
        )
        ax.set_title(tf, fontsize=10)
        ax.set_xlabel("ULM coef")
    fig.suptitle("Within-cohort Q4 vs Q1 (GSE189357 tails too small to test)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def direction_table(focus_tests: pd.DataFrame, programs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for program, tf, expected in FOCUS + CONTEXT:
        for method in METHODS:
            hit = focus_tests[(focus_tests["method"] == method) & (focus_tests["tf"] == tf)]
            if hit.empty or not np.isfinite(hit.iloc[0]["coef"]):
                rows.append(
                    {
                        "class": "TF activity", "program": program, "name": tf, "method": method,
                        "coef": np.nan, "p": np.nan, "vs_kd": "not scored", "expected": expected,
                    }
                )
                continue
            r = hit.iloc[0]
            rows.append(
                {
                    "class": "TF activity", "program": program, "name": tf, "method": method,
                    "coef": r["coef"], "se": r["se"], "p": r["p"], "fdr": r["fdr"],
                    "delta": r["delta"], "n_targets": r["n_targets"],
                    "vs_kd": matches_kd(float(r["coef"]), expected), "expected": expected,
                }
            )
    for r in programs.itertuples():
        name = r.program
        if name.startswith("gene:"):
            expected = "down" if name.split(":", 1)[1] in {"NLRC5", "CIITA", "STAT1", "IRF1", "HLA-A", "B2M", "TAP1", "TAP2"} else "none"
            if name.endswith("CLDN4"):
                expected = "up_splitter"
            program = "gene"
        elif name in {"IFN_hallmark", "MHC_I_APM"}:
            expected = "down"
            program = "geneset"
        else:
            expected = "none"
            program = "geneset"
        coef = float(r.coef) if np.isfinite(r.coef) else np.nan
        if expected == "up_splitter":
            vs = "splitter check (high should be up)" if coef > 0 else "splitter check failed"
        else:
            vs = matches_kd(coef, expected)
        rows.append(
            {
                "class": program, "program": program, "name": name, "method": "mean_logCPM" if program == "geneset" else "gene_logCPM",
                "coef": coef, "se": getattr(r, "se", np.nan), "p": getattr(r, "p", np.nan),
                "vs_kd": vs, "expected": expected, "n_genes": getattr(r, "n_genes", np.nan),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    meta = load_meta()
    counts = load_counts(meta)
    # align meta to columns
    meta = meta.loc[meta["patient"].isin(counts.columns)].copy()
    logcpm = normalize(counts.loc[:, meta["patient"]])
    meta = meta.loc[meta["patient"].isin(logcpm.columns)].copy()
    print(f"matrix {logcpm.shape[0]} genes x {logcpm.shape[1]} units")
    q = meta.loc[meta["quartile"].isin(["Q1", "Q4"])]
    print("Q1/Q4", int((q.quartile == "Q1").sum()), int((q.quartile == "Q4").sum()))

    net = load_net()
    print("scoring DoRothEA", len(net), "edges", net["source"].nunique(), "TFs")
    scores = score_dorothea(logcpm, net)
    nt = n_targets(net, logcpm.index)
    tests = []
    for method, mat in scores.items():
        print(" testing", method, mat.shape)
        tests.append(test_scores(mat, meta, method, nt))
    tests_df = pd.concat(tests, ignore_index=True)

    focus_names = [tf for _, tf, _ in FOCUS + CONTEXT]
    exp_map = {tf: (prog, exp) for prog, tf, exp in FOCUS + CONTEXT}
    focus = tests_df.loc[tests_df["tf"].isin(focus_names)].copy()
    focus["program"] = focus["tf"].map(lambda t: exp_map.get(t, ("", ""))[0])
    focus["expected"] = focus["tf"].map(lambda t: exp_map.get(t, ("", ""))[1])
    # keep requested TFs that were not scored
    missing = []
    for prog, tf, exp in FOCUS + CONTEXT:
        if tf not in set(focus["tf"]):
            missing.append(
                {
                    "method": "ulm", "tf": tf, "program": prog, "expected": exp,
                    "n_targets": int(nt.get(tf, 0)), "coef": np.nan, "p": np.nan,
                    "note": "below tmin or absent from DoRothEA A+B+C",
                }
            )
    if missing:
        focus = pd.concat([focus, pd.DataFrame(missing)], ignore_index=True)

    programs = gene_program_tests(logcpm, meta)
    a8 = json.loads((DATA / "a8_sets.json").read_text())["sets"]
    ifn_genes = {g.upper() for g in set(a8["HALLMARK_INTERFERON_GAMMA_RESPONSE"]) | set(a8["HALLMARK_INTERFERON_ALPHA_RESPONSE"])}
    logfc = gene_logfc(logcpm, meta)
    split = stat1_split(net, logfc, ifn_genes)

    cohort_ulm = per_cohort_tests(scores["ulm"], meta, focus_names)
    loo = loo_tests(scores["ulm"], meta, ["STAT1", "IRF1", "IRF9", "RFX5", "GRHL2", "ELF3", "KLF4"])
    print("within-cohort ULM")
    within = within_cohort_ulm(counts, meta, net, ["STAT1", "IRF1", "RFX5", "GRHL2", "ELF3"])

    # unscaled ±1 weight sensitivity for ULM
    net_pm = net.copy()
    net_pm["weight"] = np.sign(net_pm["weight"]).astype(float)
    ulm_pm, _ = dc.mt.ulm(logcpm.T.astype(float), net_pm, tmin=TMIN)
    pm_tests = test_scores(ulm_pm, meta, "ulm_sign_only", n_targets(net_pm, logcpm.index))
    pm_focus = pm_tests.loc[pm_tests["tf"].isin(focus_names)].copy()

    print("CollecTRI sensitivity")
    col = load_collectri()
    col_scores, _ = dc.mt.ulm(logcpm.T.astype(float), col, tmin=TMIN)
    col_tests = test_scores(col_scores, meta, "collectri_ulm", n_targets(col, logcpm.index))
    col_focus = col_tests.loc[col_tests["tf"].isin(COLLECTRI_FOCUS)].copy()
    absent_col = [tf for tf in COLLECTRI_FOCUS if tf not in set(col_tests["tf"])]

    direction = direction_table(focus.dropna(subset=["coef"]), programs)

    # save
    tests_df.to_csv(OUT / "tf_tests_all.tsv", sep="\t", index=False)
    focus.to_csv(OUT / "tf_focus.tsv", sep="\t", index=False)
    programs.to_csv(OUT / "gene_programs.tsv", sep="\t", index=False)
    split.to_csv(OUT / "stat1_target_split.tsv", sep="\t", index=False)
    cohort_ulm.to_csv(OUT / "ulm_per_cohort.tsv", sep="\t", index=False)
    loo.to_csv(OUT / "ulm_loo.tsv", sep="\t", index=False)
    within.to_csv(OUT / "ulm_within_cohort_tmm.tsv", sep="\t", index=False)
    pm_focus.to_csv(OUT / "ulm_sign_only_weights.tsv", sep="\t", index=False)
    col_focus.to_csv(OUT / "collectri_ulm_focus.tsv", sep="\t", index=False)
    direction.to_csv(OUT / "direction_vs_kd.tsv", sep="\t", index=False)

    # patient scores for the focus TFs
    long_rows = []
    for method, mat in scores.items():
        keep = [tf for tf in focus_names if tf in mat.columns]
        slim = mat[keep].copy()
        slim["patient"] = slim.index
        long = slim.melt(id_vars="patient", var_name="tf", value_name="activity")
        long["method"] = method
        long_rows.append(long)
    long_df = pd.concat(long_rows, ignore_index=True).merge(
        meta[["patient", "cohort", "quartile", "cldn4_pct"]], on="patient", how="left"
    )
    long_df.to_csv(OUT / "patient_activity_focus.tsv", sep="\t", index=False)

    forest(focus, FIGS / "fig_ulm_forest")
    method_heatmap(focus.dropna(subset=["coef"]), FIGS / "fig_method_coefs")
    cohort_forest(cohort_ulm, FIGS / "fig_ulm_by_cohort")
    box_activities(scores["ulm"], meta, FIGS / "fig_ulm_boxes")
    stat1_bars(split, FIGS / "fig_stat1_targets")

    summary = {
        "decoupler": dc.__version__,
        "n_units": int(meta["patient"].nunique()),
        "n_q1": int((meta["quartile"] == "Q1").sum()),
        "n_q4": int((meta["quartile"] == "Q4").sum()),
        "n_genes": int(logcpm.shape[0]),
        "n_tfs_ulm": int(scores["ulm"].shape[1]),
        "nlrc5_in_dorothea_abc": bool((net["source"] == "NLRC5").any()),
        "ciita_in_dorothea_abc": bool((net["source"] == "CIITA").any()),
        "collectri_absent": absent_col,
        "irf7_n_targets_in_matrix": int(nt.get("IRF7", 0)),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    show = focus.loc[focus["method"].isin(["ulm", "wmean"]), ["method", "program", "tf", "n_targets", "coef", "p", "delta"]]
    print(show.sort_values(["method", "program", "tf"]).to_string(index=False))
    print("--- programs ---")
    print(programs.to_string(index=False))
    print("--- STAT1 split ---")
    print(split.to_string(index=False))
    print("--- LOO STAT1 ---")
    print(loo.loc[loo["tf"] == "STAT1"].to_string(index=False))
    print("--- collectri ---")
    print(col_focus[["tf", "n_targets", "coef", "p"]].to_string(index=False))
    print("--- sign-only STAT1/GRHL2 ---")
    print(pm_focus.loc[pm_focus["tf"].isin(["STAT1", "IRF1", "RFX5", "GRHL2", "ELF3"]), ["tf", "coef", "p"]].to_string(index=False))


if __name__ == "__main__":
    main()

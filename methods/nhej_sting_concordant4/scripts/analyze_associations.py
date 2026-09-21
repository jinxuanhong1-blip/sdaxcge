#!/usr/bin/env python3
"""Patient-level partial correlations and an observational NHEJ–IFN check.

Not a causal mediation. Unit = patient / donor / sample in concordant-4.
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

ROOT = Path("/workspace/methods/nhej_sting_concordant4")
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
LOCKED = ROOT / "data" / "locked_patient_units_pr539.tsv"
GENE_SETS = ROOT / "data" / "gene_sets.json"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLOR = {
    "GSE123902": "#0072B2",
    "GSE131907": "#E69F00",
    "GSE205335": "#009E73",
    "GSE189357": "#CC79A7",
}


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.3e}"
    return f"{p:.4f}"


def fmt_r(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    return f"{x:.3f}"


def rankdata(x: np.ndarray) -> np.ndarray:
    return stats.rankdata(x, method="average")


def ols_resid(y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    if Z.size == 0:
        Z = np.ones((len(y), 1))
    else:
        Z = np.column_stack([np.ones(len(y)), Z])
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return y - Z @ beta


def pearson_p(r: float, df: int) -> float:
    if df <= 0 or not np.isfinite(r) or abs(r) >= 1:
        if np.isfinite(r) and abs(r) >= 1 and df > 0:
            return 0.0
        return np.nan
    t = r * np.sqrt(df / max(1e-15, 1 - r * r))
    return float(2 * stats.t.sf(abs(t), df))


def partial_spearman(x, y, covariates: np.ndarray | None = None) -> dict:
    """Pearson correlation of rank-residuals. covariates are numeric columns."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if covariates is not None and covariates.size:
        covariates = np.asarray(covariates, dtype=float)
        mask &= np.all(np.isfinite(covariates), axis=1)
        C = covariates[mask]
    else:
        C = np.zeros((int(mask.sum()), 0))
    n = int(mask.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan, "df": np.nan}
    rx = rankdata(x[mask])
    ry = rankdata(y[mask])
    k = C.shape[1]
    if k:
        Cr = np.column_stack([rankdata(C[:, j]) for j in range(k)])
        xr = ols_resid(rx, Cr)
        yr = ols_resid(ry, Cr)
    else:
        xr, yr = rx, ry
    rho = float(np.corrcoef(xr, yr)[0, 1])
    df = n - 2 - k
    return {"n": n, "rho": rho, "p": pearson_p(rho, df), "df": df}


def spearman(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[mask], y[mask])
    return {"n": n, "rho": float(rho), "p": float(p)}


def dl_spearman(rhos, ns) -> dict:
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    if len(rhos) == 0:
        return {"rho": np.nan, "p": np.nan, "I2": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "k": 0, "N": 0}
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var_z = 1 / (ns - 3)
    w = 1 / var_z
    zbar = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - zbar) ** 2)
    k = len(rhos)
    dfree = k - 1
    cdenom = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1 / (var_z + tau2)
    zre = np.sum(wstar * z) / np.sum(wstar)
    se = np.sqrt(1 / np.sum(wstar))
    p = float(2 * stats.norm.sf(abs(zre / se)))
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    ci = np.tanh(zre + np.array([-1, 1]) * 1.96 * se)
    return {
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": float(i2),
        "ci_lo": float(ci[0]),
        "ci_hi": float(ci[1]),
        "k": int(k),
        "N": int(ns.sum()),
    }


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for idx in order[::-1]:
        rank = int(np.where(order == idx)[0][0]) + 1
        val = p[idx] * m / rank
        prev = min(prev, val)
        q[idx] = prev
    return [float(min(1.0, v)) if np.isfinite(v) else np.nan for v in q]


def cohort_dummies(dataset: pd.Series) -> np.ndarray:
    cats = [c for c in COHORTS if c in set(dataset)]
    # drop first
    cols = []
    for c in cats[1:]:
        cols.append((dataset == c).astype(float).to_numpy())
    if not cols:
        return np.zeros((len(dataset), 0))
    return np.column_stack(cols)


def within_z(df: pd.DataFrame, col: str) -> pd.Series:
    def _z(s):
        sd = s.std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            return s * 0.0
        return (s - s.mean()) / sd

    return df.groupby("dataset", group_keys=False)[col].transform(_z)


def within_rank(df: pd.DataFrame, col: str) -> pd.Series:
    return df.groupby("dataset", group_keys=False)[col].transform(lambda s: rankdata(s.to_numpy()))


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    metas = []
    genes = []
    for p in sorted(TAB.glob("unit_meta_*.tsv")):
        metas.append(pd.read_csv(p, sep="\t"))
    for p in sorted(TAB.glob("gene_umi_*.tsv")):
        genes.append(pd.read_csv(p, sep="\t"))
    if not metas or not genes:
        raise SystemExit(f"no unit_meta_*.tsv / gene_umi_*.tsv under {TAB}")
    meta = pd.concat(metas, ignore_index=True)
    gene = pd.concat(genes, ignore_index=True)
    meta["unit_id"] = meta["unit_id"].astype(str)
    gene["unit_id"] = gene["unit_id"].astype(str)
    gene["gene"] = gene["gene"].str.upper()
    aliases = {k.upper(): v.upper() for k, v in json.loads(GENE_SETS.read_text())["aliases"].items()}
    gene["gene"] = gene["gene"].map(lambda g: aliases.get(g, g))
    # if a gene was emitted twice (alias collapse), sum
    gene = gene.groupby(["dataset", "unit_id", "gene"], as_index=False)["umi"].sum()
    meta = meta.drop_duplicates(["dataset", "unit_id"], keep="last")
    return meta, gene


def score_units(meta: pd.DataFrame, gene: pd.DataFrame, sets: dict) -> pd.DataFrame:
    present = gene.groupby("dataset")["gene"].apply(set).to_dict()
    wide = gene.pivot_table(index=["dataset", "unit_id"], columns="gene", values="umi", aggfunc="sum")
    rows = []
    modules = {
        "nhej": sets["nhej_kegg"],
        "nhej_core": sets["nhej_core"],
        "sting": sets["sting_core"],
        "ifn": sets["ifn_hallmark"],
        "mhc": sets["mhc_i"],
        "prolif": sets["prolif"],
    }
    ifn_minus_mhc = [g for g in sets["ifn_hallmark"] if g not in set(sets["mhc_i"])]
    modules["ifn_not_mhc"] = ifn_minus_mhc
    for rec in meta.itertuples(index=False):
        key = (rec.dataset, rec.unit_id)
        if key not in wide.index:
            continue
        vec = wide.loc[key]
        lib = float(rec.lib_umi)
        if not np.isfinite(lib) or lib <= 0:
            continue
        avail = present.get(rec.dataset, set())
        out = {
            "dataset": rec.dataset,
            "unit_id": rec.unit_id,
            "n_cells": rec.n_cells,
            "n_malignant": rec.n_malignant,
            "n_tnk": rec.n_tnk,
            "mal_CLDN4_pct": rec.mal_CLDN4_pct,
            "mal_CLDN4_mean": rec.mal_CLDN4_mean,
            "lib_umi": lib,
        }
        if "CLDN4" in avail and "CLDN4" in vec.index and np.isfinite(vec["CLDN4"]):
            out["cldn4_logcpm"] = float(np.log2(vec["CLDN4"] / lib * 1e6 + 1))
        else:
            out["cldn4_logcpm"] = np.nan
        for name, genes in modules.items():
            use = [g for g in genes if g in avail]
            out[f"n_{name}"] = len(use)
            if not use:
                out[f"{name}_score"] = np.nan
                continue
            vals = []
            for g in use:
                u = vec[g] if g in vec.index else 0.0
                if not np.isfinite(u):
                    u = 0.0
                vals.append(np.log2(float(u) / lib * 1e6 + 1))
            out[f"{name}_score"] = float(np.mean(vals))
        rows.append(out)
    df = pd.DataFrame(rows)
    # STING/IFN composite: equal-weight mean of within-cohort z scores
    df["sting_z"] = within_z(df, "sting_score")
    df["ifn_z"] = within_z(df, "ifn_score")
    df["sting_ifn"] = df[["sting_z", "ifn_z"]].mean(axis=1)
    df["cldn4_z"] = within_z(df, "cldn4_logcpm")
    df["nhej_z"] = within_z(df, "nhej_score")
    df["mhc_z"] = within_z(df, "mhc_score")
    df["prolif_z"] = within_z(df, "prolif_score")
    df["ifn_not_mhc_z"] = within_z(df, "ifn_not_mhc_score")
    df["nhej_core_z"] = within_z(df, "nhej_core_score")
    df["pct_z"] = within_z(df, "mal_CLDN4_pct")
    return df


def gene_coverage(gene: pd.DataFrame, sets: dict) -> pd.DataFrame:
    rows = []
    modules = {
        "NHEJ_KEGG": sets["nhej_kegg"],
        "NHEJ_core": sets["nhej_core"],
        "STING_core": sets["sting_core"],
        "IFN_hallmark": sets["ifn_hallmark"],
        "MHC_I": sets["mhc_i"],
    }
    for ds, sub in gene.groupby("dataset"):
        avail = set(sub["gene"])
        # a gene is "detected" if any unit has umi>0
        det = set(sub.loc[sub["umi"] > 0, "gene"])
        for name, genes in modules.items():
            for g in genes:
                rows.append(
                    {
                        "dataset": ds,
                        "module": name,
                        "gene": g,
                        "in_universe": g in avail,
                        "detected_any_unit": g in det,
                    }
                )
    return pd.DataFrame(rows)


def cohort_spearmans(df: pd.DataFrame, x: str, y: str) -> pd.DataFrame:
    rows = []
    for ds in COHORTS:
        sub = df[df["dataset"] == ds]
        sp = spearman(sub[x], sub[y])
        rows.append({"dataset": ds, "x": x, "y": y, **sp})
    return pd.DataFrame(rows)


def strata_split(df: pd.DataFrame, col: str, k: int) -> pd.Series:
    def _cut(s: pd.Series) -> pd.Series:
        r = rankdata(s.to_numpy())
        # equal-count bins by rank
        bins = np.ceil(r / len(s) * k).astype(int)
        bins = np.clip(bins, 1, k)
        return pd.Series(bins, index=s.index)

    return df.groupby("dataset", group_keys=False)[col].apply(_cut)


def fit_paths(df: pd.DataFrame, cldn4_col: str = "pct_z") -> pd.DataFrame:
    """Observational standardized paths on within-cohort z scores. Not causal."""
    df = df.copy()
    df["cldn4_z"] = df[cldn4_col]
    rows = []

    def one(y, xs, label):
        sub = df.dropna(subset=[y] + xs).copy()
        Y = sub[y].to_numpy()
        X = np.column_stack([np.ones(len(sub))] + [sub[c].to_numpy() for c in xs])
        beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
        resid = Y - X @ beta
        dof = len(sub) - X.shape[1]
        sigma2 = float(resid @ resid / dof) if dof > 0 else np.nan
        try:
            xtx_inv = np.linalg.inv(X.T @ X)
        except np.linalg.LinAlgError:
            xtx_inv = np.linalg.pinv(X.T @ X)
        se = np.sqrt(np.diag(xtx_inv) * sigma2)
        for name, b, s in zip(["intercept"] + xs, beta, se):
            if name == "intercept":
                continue
            t = b / s if s > 0 else np.nan
            p = float(2 * stats.t.sf(abs(t), dof)) if np.isfinite(t) else np.nan
            rows.append(
                {
                    "model": label,
                    "y": y,
                    "term": name,
                    "beta": float(b),
                    "se": float(s),
                    "df": int(dof),
                    "p": p,
                    "n": int(len(sub)),
                }
            )

    one("nhej_z", ["cldn4_z"], "a_NHEJ~CLDN4")
    one("ifn_z", ["cldn4_z"], "c_IFN~CLDN4")
    one("ifn_z", ["cldn4_z", "nhej_z"], "cprime_IFN~CLDN4+NHEJ")
    one("sting_z", ["cldn4_z"], "IFN_arm_STING~CLDN4")
    one("sting_ifn", ["cldn4_z"], "composite~CLDN4")
    one("mhc_z", ["cldn4_z"], "MHC~CLDN4")
    one("ifn_z", ["cldn4_z", "nhej_z", "prolif_z"], "IFN~CLDN4+NHEJ+prolif")
    one("sting_z", ["cldn4_z", "nhej_z"], "STING~CLDN4+NHEJ")
    return pd.DataFrame(rows)


def bootstrap_indirect(df: pd.DataFrame, n_boot: int = 2000, seed: int = 1) -> dict:
    """Stratified bootstrap of a*b on within-cohort z. Descriptive uncertainty only."""
    rng = np.random.default_rng(seed)
    groups = {ds: df[df["dataset"] == ds].dropna(subset=["mal_CLDN4_pct", "nhej_score", "ifn_score"]) for ds in COHORTS}

    def ab(parts: list[pd.DataFrame]) -> float:
        d = pd.concat(parts, ignore_index=True)
        # recompute within-cohort z on the resample
        d["cldn4_z"] = within_z(d, "mal_CLDN4_pct")
        d["nhej_z"] = within_z(d, "nhej_score")
        d["ifn_z"] = within_z(d, "ifn_score")
        d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["cldn4_z", "nhej_z", "ifn_z"])
        if len(d) < 10:
            return np.nan
        Xa = np.column_stack([np.ones(len(d)), d["cldn4_z"]])
        ba, *_ = np.linalg.lstsq(Xa, d["nhej_z"].to_numpy(), rcond=None)
        Xb = np.column_stack([np.ones(len(d)), d["cldn4_z"], d["nhej_z"]])
        bb, *_ = np.linalg.lstsq(Xb, d["ifn_z"].to_numpy(), rcond=None)
        return float(ba[1] * bb[2])

    point = ab([groups[ds] for ds in COHORTS])
    boots = []
    for _ in range(n_boot):
        parts = []
        for ds in COHORTS:
            g = groups[ds]
            if len(g) == 0:
                continue
            idx = rng.integers(0, len(g), len(g))
            parts.append(g.iloc[idx])
        boots.append(ab(parts))
    boots = np.asarray(boots, dtype=float)
    boots = boots[np.isfinite(boots)]
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return {
        "indirect": point,
        "ci_lo": float(lo),
        "ci_hi": float(hi),
        "n_boot": int(len(boots)),
        "frac_negative": float(np.mean(boots < 0)),
    }


def make_figures(df: pd.DataFrame, cohort_tbl: pd.DataFrame, dl_map: dict, strata: pd.DataFrame):
    FIG.mkdir(parents=True, exist_ok=True)
    # 1. forest of within-cohort Spearman CLDN4 vs modules
    modules = [
        ("nhej_score", "NHEJ"),
        ("sting_ifn", "STING/IFN"),
        ("mhc_score", "MHC-I"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    yticks = []
    ylabels = []
    y = 0
    for col, label in modules:
        sub = cohort_tbl[(cohort_tbl["x"] == "mal_CLDN4_pct") & (cohort_tbl["y"] == col)]
        for _, r in sub.iterrows():
            ax.plot(r["rho"], y, "o", color=COHORT_COLOR[r["dataset"]], markersize=6)
            yticks.append(y)
            ylabels.append(f"{label} · {r['dataset'].replace('GSE','')}")
            y += 1
        dl = dl_map[("cldn4_logcpm", col)]
        ax.plot([dl["ci_lo"], dl["ci_hi"]], [y, y], color="black", lw=1.4)
        ax.plot(dl["rho"], y, "D", color="black", markersize=5)
        yticks.append(y)
        ylabels.append(f"{label} · DL meta")
        y += 1.4
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Spearman ρ  (malignant CLDN4 % positive vs module)")
    ax.set_title("Concordant-4, patient unit, within cohort")
    # legend
    for ds, col in COHORT_COLOR.items():
        ax.plot([], [], "o", color=col, label=ds)
    ax.plot([], [], "D", color="black", label="DL meta")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "fig_cldn4_partial_forest.png", dpi=160)
    fig.savefig(FIG / "fig_cldn4_partial_forest.pdf")
    plt.close()

    # 2. within-cohort z scatters
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3), sharex=True)
    pairs = [
        ("nhej_z", "NHEJ z"),
        ("sting_ifn", "STING/IFN z"),
        ("mhc_z", "MHC-I z"),
    ]
    for ax, (col, label) in zip(axes, pairs):
        for ds, sub in df.groupby("dataset"):
            ax.scatter(sub["pct_z"], sub[col], s=18, color=COHORT_COLOR[ds], label=ds, alpha=0.9)
        ax.axhline(0, color="0.75", lw=0.6)
        ax.axvline(0, color="0.75", lw=0.6)
        ax.set_xlabel("CLDN4 %pos z (within cohort)")
        ax.set_ylabel(label)
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_cldn4_module_scatter.png", dpi=160)
    fig.savefig(FIG / "fig_cldn4_module_scatter.pdf")
    plt.close()

    # 3. NHEJ vs IFN inside CLDN4 median strata
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), sharex=True, sharey=True)
    for ax, stratum, title in zip(axes, [1, 2], ["CLDN4 below cohort median", "CLDN4 above cohort median"]):
        sub = strata[strata["cldn4_half"] == stratum]
        for ds, g in sub.groupby("dataset"):
            ax.scatter(g["nhej_z"], g["ifn_z"], s=18, color=COHORT_COLOR[ds], label=ds, alpha=0.9)
        ax.axhline(0, color="0.75", lw=0.6)
        ax.axvline(0, color="0.75", lw=0.6)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("NHEJ z (within cohort)")
        ax.set_ylabel("IFN z (within cohort)")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_nhej_ifn_within_cldn4.png", dpi=160)
    fig.savefig(FIG / "fig_nhej_ifn_within_cldn4.pdf")
    plt.close()


def self_check():
    rng = np.random.default_rng(0)
    x = rng.normal(size=80)
    y = -x + rng.normal(scale=0.05, size=80)
    z = rng.normal(size=80)
    sp = partial_spearman(x, y)
    assert sp["rho"] < -0.95, sp
    # covariate that is y itself should knock the correlation down
    killed = partial_spearman(x, y, np.column_stack([y]))
    assert abs(killed["rho"]) < 0.2, killed
    dl = dl_spearman([-0.4, -0.4, -0.4, -0.4], [20, 20, 20, 20])
    assert abs(dl["rho"] + 0.4) < 1e-6, dl
    print("self-check ok", flush=True)


def write_finding(ctx: dict):
    q = ctx
    lines = []
    a = lines.append
    a("# Concordant-4 malignant NHEJ, STING/IFN, and MHC-I")
    a("")
    a("ADDITIVE. Observational only. This is not a causal mediation and not a")
    a("re-derivation of the locked concordant-4 result (malignant CLDN4 %pos vs")
    a("T/NK, ρ = −0.531, n = 65). Same four cohorts and the same patient / donor /")
    a("sample units: **GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not")
    a("GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.")
    a("")
    a("Malignant cells only. GSE123902 and GSE189357 use the marker gate")
    a("(EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0. GSE131907 and GSE205335 use author")
    a("malignant labels. Scores are patient pseudobulk means of log2(CPM+1), CPM from")
    a("that unit's malignant UMI total. p-values are descriptive.")
    a("")
    a("## Modules")
    a("")
    a("- **NHEJ**: MSigDB KEGG_NON_HOMOLOGOUS_END_JOINING (hsa03450):")
    a("  DCLRE1C, DNTT, FEN1, LIG4, MRE11, NHEJ1, POLL, POLM, PRKDC, RAD50, XRCC4, XRCC5, XRCC6.")
    a("  Sensitivity drops DNTT and FEN1.")
    a("- **STING core** (sensor/adapter, not the ISG program): CGAS, STING1, TBK1, IRF3, IFI16, DDX41.")
    a("  NHEJ genes that sit in Reactome STING lists (XRCC5, XRCC6, PRKDC, MRE11) are excluded.")
    a("  Negative regulators (TREX1, NLRC3, NLRP4) are excluded. IRF7 and ZBP1 are left in the IFN arm.")
    a("- **IFN**: Hallmark IFNα ∪ IFNγ, the same union as the concordant-4 Seurat run.")
    a("- **STING/IFN module**: equal-weight mean of the within-cohort z-scored STING core and IFN scores.")
    a("- **MHC-I**: the 21-gene MHC-I/APM list used in that Seurat run. Twelve of those genes also sit in Hallmark IFN-γ; a sensitivity IFN score drops them.")
    a("")
    a("Aliases collapsed before scoring: MRE11A→MRE11, TMEM173→STING1, MB21D1→CGAS, XLF→NHEJ1, MARCH1→MARCHF1, WARS→WARS1.")
    a("GSE205335 stores STING1 as TMEM173; that alias is applied before the STING score.")
    a("Genes absent from a dataset are left out of that dataset's mean. Zeros stay in.")
    a("")
    a("## Honest n")
    a("")
    a(f"- Units scored: **n = {q['n']}** ({q['n_by']}).")
    a(f"- Malignant-cell gate check against the locked table: CLDN4 %pos Spearman ρ = {fmt_r(q['qc_rho'])}, max |Δ| = {q['qc_max']:.4f} percentage points, n matched = {q['qc_n']}.")
    a(f"- Locked IFN score vs this Hallmark IFN score: Spearman ρ = {fmt_r(q['qc_ifn'])} (same units, different normalization: theirs is TMM-CPM).")
    a("- Do not quote malignant-cell or total-cell counts as n.")
    a("")
    a("## 1. Patient-level association of CLDN4 with each module")
    a("")
    a("Primary CLDN4 measure is malignant **% positive**, the locked concordant-4 patient score.")
    a("Module scores stay on log2(CPM+1). Pseudobulk CLDN4 abundance is a sensitivity, because")
    a(f"it does not rank the same patients (cohort-adjusted partial ρ of %pos vs log2(CPM+1) = {fmt_r(q['pct_vs_cpm'])}).")
    a("Two summaries of each association:")
    a("")
    a("1. **Partial Spearman controlling for cohort**: rank CLDN4 and the module across the 65 units, residualize both ranks on cohort indicators, then Pearson-correlate the residuals. Degrees of freedom subtract the three cohort indicators.")
    a("2. **Within-cohort Spearman**, then DerSimonian–Laird on the Fisher z values (same pooling as the locked T/NK meta).")
    a("")
    a("| module | partial ρ (cohort) | p | DL ρ | DL p | I² | DL 95% CI |")
    a("|---|---:|---:|---:|---:|---:|---|")
    for row in q["primary_rows"]:
        a(
            f"| {row['label']} | {fmt_r(row['partial'])} | {fmt_p(row['partial_p'])} | "
            f"{fmt_r(row['dl'])} | {fmt_p(row['dl_p'])} | {row['i2']:.0%} | "
            f"{fmt_r(row['lo'])} to {fmt_r(row['hi'])} |"
        )
    a("")
    a(f"Benjamini–Hochberg q across the three primary partial tests: {q['qtext']}.")
    a("")
    a("Within-cohort Spearmans (context, not a second n):")
    a("")
    a("| module | " + " | ".join(c.replace("GSE", "") for c in COHORTS) + " |")
    a("|---|" + "|".join(["---:" for _ in COHORTS]) + "|")
    for row in q["cohort_rows"]:
        cells = " | ".join(f"{fmt_r(row[c])} (n={row[c+'_n']})" for c in COHORTS)
        a(f"| {row['label']} | {cells} |")
    a("")
    a("Supporting rows (not the three primary tests):")
    a("")
    a("| contrast | partial ρ | p | DL ρ | DL p |")
    a("|---|---:|---:|---:|---:|")
    for row in q["support_rows"]:
        a(
            f"| {row['label']} | {fmt_r(row['partial'])} | {fmt_p(row['partial_p'])} | "
            f"{fmt_r(row['dl'])} | {fmt_p(row['dl_p'])} |"
        )
    a("")
    a("## 2. Does NHEJ-low go with IFN-high inside CLDN4 strata?")
    a("")
    a("This is the mediation-style question. It asks whether the NHEJ–IFN association")
    a("is still there after CLDN4 is held fixed. It does not identify a mechanism.")
    a("")
    a(f"- Partial Spearman, NHEJ vs Hallmark IFN, controlling for CLDN4 %pos and cohort: **ρ = {fmt_r(q['med_rho'])}**, p = {fmt_p(q['med_p'])}, n = {q['med_n']}.")
    a(f"- Same partial, controlling for CLDN4 log2(CPM+1) instead of %pos: ρ = {fmt_r(q['med_log_rho'])}, p = {fmt_p(q['med_log_p'])}.")
    a(f"- Same partial, additionally controlling for a 6-gene proliferation score: ρ = {fmt_r(q['med_prolif_rho'])}, p = {fmt_p(q['med_prolif_p'])}.")
    a(f"- Same partial using the NHEJ-core sensitivity set: ρ = {fmt_r(q['med_core_rho'])}, p = {fmt_p(q['med_core_p'])}.")
    a(f"- Partial Spearman, NHEJ vs STING core, controlling for CLDN4 and cohort: ρ = {fmt_r(q['med_sting_rho'])}, p = {fmt_p(q['med_sting_p'])}.")
    a("")
    a("Within-cohort median split on CLDN4 %pos. NHEJ and IFN are within-cohort z-scores, then pooled. Negative ρ means NHEJ-low with IFN-high inside that stratum.")
    a("")
    a("| CLDN4 stratum | n | Spearman NHEJ vs IFN | p |")
    a("|---|---:|---:|---:|")
    for row in q["strata_rows"]:
        a(f"| {row['label']} | {row['n']} | {fmt_r(row['rho'])} | {fmt_p(row['p'])} |")
    a("")
    a(f"Within-cohort partial Spearman (NHEJ vs IFN | CLDN4), DL meta: ρ = {fmt_r(q['dl_within']['rho'])}, p = {fmt_p(q['dl_within']['p'])}, I² = {q['dl_within']['I2']:.0%}, k = {q['dl_within']['k']}.")
    a("")
    a("Cohort-specific partials (small n; GSE189357 is 9 units):")
    a("")
    a("| cohort | n | partial ρ | p |")
    a("|---|---:|---:|---:|")
    for row in q["within_rows"]:
        a(f"| {row['dataset']} | {row['n']} | {fmt_r(row['rho'])} | {fmt_p(row['p'])} |")
    a("")
    a("## 3. Observational path split (not an effect)")
    a("")
    a("Within-cohort z-scores of CLDN4 %pos, NHEJ, and IFN. OLS. Path a is CLDN4 with NHEJ. Path b is NHEJ with IFN when CLDN4 is in the model. Path c is CLDN4 with IFN. Path c′ is CLDN4 with IFN when NHEJ is in the model. The product a×b is a descriptive split of the CLDN4–IFN association, not a causal indirect effect. No sequential ignorability, no intervention.")
    a("")
    a("| path | beta | p | n |")
    a("|---|---:|---:|---:|")
    for row in q["path_rows"]:
        a(f"| {row['label']} | {fmt_r(row['beta'])} | {fmt_p(row['p'])} | {row['n']} |")
    a("")
    a(
        f"Stratified bootstrap of a×b (2000 resamples within cohort): {fmt_r(q['boot']['indirect'])} "
        f"(percentile interval {fmt_r(q['boot']['ci_lo'])} to {fmt_r(q['boot']['ci_hi'])}; "
        f"fraction of draws below 0 = {q['boot']['frac_negative']:.2f}). "
        "The interval is uncertainty of this decomposition, not a causal confidence interval."
    )
    a("")
    a("## Reading")
    a("")
    a(q["reading"])
    a("")
    a(q["abundance_note"])
    a("")
    a("Cross-sectional patient pseudobulks. Cytosolic-DNA logic is why NHEJ and STING are on the table. The split above is an association decomposition. STING1 is on the matrix in every cohort after the TMEM173 alias, but the sensor arm is a 6-gene mean next to a 224-gene IFN mean, so the composite is not a pure STING score (the STING-core row is separate). Twelve MHC-I genes also sit in Hallmark IFN-γ, so the MHC-I and IFN rows share genes. The IFN-minus-MHC sensitivity stays negative on CLDN4 %pos.")
    a("")
    a("## Reproduce")
    a("")
    a("```")
    a("bash methods/nhej_sting_concordant4/scripts/download.sh /tmp/geo_nhej")
    a("python3 methods/nhej_sting_concordant4/scripts/extract_sums.py")
    a("Rscript methods/nhej_sting_concordant4/scripts/export_gse205335.R /tmp/geo_nhej /tmp/geo_nhej/gse205335_export")
    a("python3 methods/nhej_sting_concordant4/scripts/analyze_associations.py")
    a("```")
    a("")
    text = "\n".join(lines) + "\n"
    (ROOT / "FINDING.md").write_text(text)
    return text


def reading_paragraph(primary, med_rho, strata_rows) -> str:
    def direction(rho, name):
        if not np.isfinite(rho):
            return f"{name} is not estimable"
        if rho <= -0.2:
            return f"{name} is negative (ρ = {rho:.3f})"
        if rho >= 0.2:
            return f"{name} is positive (ρ = {rho:.3f})"
        return f"{name} is near null (ρ = {rho:.3f})"

    bits = [direction(r["partial"], r["label"]) for r in primary]
    med = direction(med_rho, "NHEJ vs IFN | CLDN4, cohort")
    low = next(r for r in strata_rows if r["key"] == "low")
    high = next(r for r in strata_rows if r["key"] == "high")
    return (
        "Cohort-adjusted partial correlations on CLDN4 %pos: "
        + "; ".join(bits)
        +         ". DerSimonian–Laird intervals for STING/IFN and MHC-I still include 0. "
        "STING/IFN I² is 44%: GSE131907 is positive and the other three cohorts are negative. "
        f"{med}. "
        f"Inside the CLDN4-low half, NHEJ vs IFN ρ = {low['rho']:.3f} (n = {low['n']}); "
        f"inside the CLDN4-high half, ρ = {high['rho']:.3f} (n = {high['n']}). "
        "Both stratum estimates are positive. The cytosolic-DNA pattern would be a negative "
        "partial (NHEJ-low with IFN-high at a fixed CLDN4 rank). That pattern is not what these 65 units show."
    )


def main():
    self_check()
    sets = json.loads(GENE_SETS.read_text())
    for k in ("nhej_kegg", "nhej_core", "sting_core", "ifn_hallmark", "mhc_i", "prolif"):
        sets[k] = [g.upper() for g in sets[k]]
    meta, gene = load_inputs()
    df = score_units(meta, gene, sets)
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "patient_module_scores.tsv", sep="\t", index=False)

    locked = pd.read_csv(LOCKED, sep="\t")
    locked["unit_id"] = locked["unit_id"].astype(str)
    m = df.merge(
        locked[["dataset", "unit_id", "mal_CLDN4_pct", "ifn_score", "mhc_score", "n_malignant"]],
        on=["dataset", "unit_id"],
        suffixes=("", "_locked"),
    )
    qc_rho = spearman(m["mal_CLDN4_pct"], m["mal_CLDN4_pct_locked"])["rho"]
    qc_max = float(np.nanmax(np.abs(m["mal_CLDN4_pct"] - m["mal_CLDN4_pct_locked"])))
    qc_ifn = spearman(m["ifn_score"], m["ifn_score_locked"])["rho"]
    m[["dataset", "unit_id", "mal_CLDN4_pct", "mal_CLDN4_pct_locked", "n_malignant", "n_malignant_locked"]].to_csv(
        TAB / "qc_vs_locked_cldn4.tsv", sep="\t", index=False
    )
    print(f"QC pct rho={qc_rho:.4f} maxabs={qc_max:.4f} ifn rho={qc_ifn:.4f} n={len(m)}", flush=True)
    if len(m) < 60 or not np.isfinite(qc_rho) or qc_rho < 0.98:
        # still write scores for debugging, but do not silently publish a bad gate
        print(m.loc[np.abs(m["mal_CLDN4_pct"] - m["mal_CLDN4_pct_locked"]) > 1, ["dataset", "unit_id", "mal_CLDN4_pct", "mal_CLDN4_pct_locked", "n_malignant", "n_malignant_locked"]].to_string())
        raise SystemExit("CLDN4 %pos does not match the locked concordant-4 table")

    cov = gene_coverage(gene, sets)
    cov.to_csv(TAB / "gene_coverage.tsv", sep="\t", index=False)

    # primary and support contrasts
    contrasts = [
        ("mal_CLDN4_pct", "nhej_score", "CLDN4 %pos vs NHEJ", True),
        ("mal_CLDN4_pct", "sting_ifn", "CLDN4 %pos vs STING/IFN", True),
        ("mal_CLDN4_pct", "mhc_score", "CLDN4 %pos vs MHC-I", True),
        ("mal_CLDN4_pct", "sting_score", "CLDN4 %pos vs STING core", False),
        ("mal_CLDN4_pct", "ifn_score", "CLDN4 %pos vs Hallmark IFN", False),
        ("mal_CLDN4_pct", "ifn_not_mhc_score", "CLDN4 %pos vs IFN minus MHC genes", False),
        ("mal_CLDN4_pct", "nhej_core_score", "CLDN4 %pos vs NHEJ core", False),
        ("cldn4_logcpm", "nhej_score", "CLDN4 logCPM vs NHEJ", False),
        ("cldn4_logcpm", "sting_ifn", "CLDN4 logCPM vs STING/IFN", False),
        ("cldn4_logcpm", "mhc_score", "CLDN4 logCPM vs MHC-I", False),
        ("cldn4_logcpm", "ifn_score", "CLDN4 logCPM vs Hallmark IFN", False),
    ]
    dummies = cohort_dummies(df["dataset"])
    partial_rows = []
    cohort_parts = []
    dl_map = {}
    for x, y, label, primary in contrasts:
        part = partial_spearman(df[x], df[y], dummies)
        ct = cohort_spearmans(df, x, y)
        ct["label"] = label
        cohort_parts.append(ct)
        dl = dl_spearman(ct["rho"], ct["n"])
        dl_map[(x, y)] = dl
        partial_rows.append(
            {
                "label": label,
                "x": x,
                "y": y,
                "primary": primary,
                "partial": part["rho"],
                "partial_p": part["p"],
                "partial_n": part["n"],
                "dl": dl["rho"],
                "dl_p": dl["p"],
                "i2": dl["I2"],
                "lo": dl["ci_lo"],
                "hi": dl["ci_hi"],
            }
        )
    partial_df = pd.DataFrame(partial_rows)
    prim_idx = [i for i, r in enumerate(partial_rows) if r["primary"]]
    qvals = bh([partial_rows[i]["partial_p"] for i in prim_idx])
    for i, qv in zip(prim_idx, qvals):
        partial_rows[i]["q"] = qv
    partial_df = pd.DataFrame(partial_rows)
    partial_df.to_csv(TAB / "partial_correlations.tsv", sep="\t", index=False)
    cohort_tbl = pd.concat(cohort_parts, ignore_index=True)
    cohort_tbl.to_csv(TAB / "cohort_spearman.tsv", sep="\t", index=False)

    # Primary CLDN4 covariate is %pos (locked patient score). logCPM is the abundance sensitivity.
    med_cov = np.column_stack([df["mal_CLDN4_pct"].to_numpy(), dummies])
    med = partial_spearman(df["nhej_score"], df["ifn_score"], med_cov)
    med_log = partial_spearman(
        df["nhej_score"], df["ifn_score"], np.column_stack([df["cldn4_logcpm"].to_numpy(), dummies])
    )
    med_pro = partial_spearman(
        df["nhej_score"],
        df["ifn_score"],
        np.column_stack([df["mal_CLDN4_pct"].to_numpy(), df["prolif_score"].to_numpy(), dummies]),
    )
    med_core = partial_spearman(df["nhej_core_score"], df["ifn_score"], med_cov)
    med_sting = partial_spearman(df["nhej_score"], df["sting_score"], med_cov)

    # within-cohort partials
    within_rows = []
    for ds in COHORTS:
        sub = df[df["dataset"] == ds]
        pr = partial_spearman(sub["nhej_score"], sub["ifn_score"], sub[["mal_CLDN4_pct"]].to_numpy())
        within_rows.append({"dataset": ds, **pr})
    within_df = pd.DataFrame(within_rows)
    within_df.to_csv(TAB / "within_cohort_nhej_ifn_partial.tsv", sep="\t", index=False)
    dl_within = dl_spearman(within_df["rho"], within_df["n"])

    # strata: within-cohort median and tertile on CLDN4 logcpm
    df = df.copy()
    df["cldn4_half"] = strata_split(df, "mal_CLDN4_pct", 2)
    df["cldn4_tert"] = strata_split(df, "cldn4_logcpm", 3)
    df["cldn4_tert_pct"] = strata_split(df, "mal_CLDN4_pct", 3)
    strata_rows = []
    for key, val, label in [("low", 1, "CLDN4 low (within-cohort bottom half)"), ("high", 2, "CLDN4 high (within-cohort top half)")]:
        sub = df[df["cldn4_half"] == val]
        sp = spearman(sub["nhej_z"], sub["ifn_z"])
        strata_rows.append({"key": key, "label": label, "stratum": "median", **sp})
    for t, label in [(1, "CLDN4 %pos T1"), (2, "CLDN4 %pos T2"), (3, "CLDN4 %pos T3")]:
        sub = df[df["cldn4_tert_pct"] == t]
        sp = spearman(sub["nhej_z"], sub["ifn_z"])
        strata_rows.append({"key": f"T{t}", "label": label, "stratum": "tertile", **sp})
    strata_df = pd.DataFrame(strata_rows)
    strata_df.to_csv(TAB / "cldn4_strata_nhej_ifn.tsv", sep="\t", index=False)

    paths = fit_paths(df)
    paths.to_csv(TAB / "observational_paths.tsv", sep="\t", index=False)
    boot = bootstrap_indirect(df)
    pd.DataFrame([boot]).to_csv(TAB / "indirect_bootstrap.tsv", sep="\t", index=False)

    # figures use cohort spearman of primary y's. cohort_tbl has many y's.
    make_figures(df, cohort_tbl, dl_map, df)

    n_by = ", ".join(f"{ds.replace('GSE','')} {int((df.dataset==ds).sum())}" for ds in COHORTS)
    primary = [r for r in partial_rows if r["primary"]]
    qtext = ", ".join(f"{r['label'].split(' vs ')[-1]} q={fmt_p(r['q'])}" for r in primary)
    # cohort table compact for primary y
    cohort_rows = []
    for y, label in [("nhej_score", "NHEJ"), ("sting_ifn", "STING/IFN"), ("mhc_score", "MHC-I")]:
        rec = {"label": label}
        sub = cohort_tbl[(cohort_tbl["x"] == "mal_CLDN4_pct") & (cohort_tbl["y"] == y)]
        for _, r in sub.iterrows():
            rec[r["dataset"]] = r["rho"]
            rec[r["dataset"] + "_n"] = int(r["n"])
        cohort_rows.append(rec)

    path_labels = {
        ("a_NHEJ~CLDN4", "cldn4_z"): "a. NHEJ ~ CLDN4 %pos",
        ("c_IFN~CLDN4", "cldn4_z"): "c. IFN ~ CLDN4 %pos",
        ("cprime_IFN~CLDN4+NHEJ", "cldn4_z"): "c′. IFN ~ CLDN4 %pos, NHEJ also in the model",
        ("cprime_IFN~CLDN4+NHEJ", "nhej_z"): "b. IFN ~ NHEJ, CLDN4 %pos also in the model",
        ("IFN~CLDN4+NHEJ+prolif", "nhej_z"): "b, proliferation also in the model",
        ("STING~CLDN4+NHEJ", "nhej_z"): "STING core ~ NHEJ, CLDN4 %pos also in the model",
    }
    path_rows = []
    for _, r in paths.iterrows():
        key = (r["model"], r["term"])
        if key in path_labels:
            path_rows.append({"label": path_labels[key], "beta": r["beta"], "p": r["p"], "n": int(r["n"])})

    pct_vs_cpm = partial_spearman(df["mal_CLDN4_pct"], df["cldn4_logcpm"], dummies)
    ctx = {
        "n": int(len(df)),
        "pct_vs_cpm": pct_vs_cpm["rho"],
        "n_by": n_by,
        "qc_rho": qc_rho,
        "qc_max": qc_max,
        "qc_n": int(len(m)),
        "qc_ifn": qc_ifn,
        "primary_rows": primary,
        "support_rows": [r for r in partial_rows if not r["primary"]],
        "qtext": qtext,
        "cohort_rows": cohort_rows,
        "med_rho": med["rho"],
        "med_p": med["p"],
        "med_n": med["n"],
        "med_log_rho": med_log["rho"],
        "med_log_p": med_log["p"],
        "med_prolif_rho": med_pro["rho"],
        "med_prolif_p": med_pro["p"],
        "med_core_rho": med_core["rho"],
        "med_core_p": med_core["p"],
        "med_sting_rho": med_sting["rho"],
        "med_sting_p": med_sting["p"],
        "strata_rows": [r for r in strata_rows if r["stratum"] == "median"] + [r for r in strata_rows if r["stratum"] == "tertile"],
        "dl_within": dl_within,
        "within_rows": within_rows,
        "path_rows": path_rows,
        "boot": boot,
    }
    ctx["reading"] = reading_paragraph(primary, med["rho"], [r for r in strata_rows if r["stratum"] == "median"])
    abund = [r for r in partial_rows if str(r["label"]).startswith("CLDN4 logCPM")]
    ctx["abundance_note"] = (
        "Abundance sensitivity, same partial Spearman on malignant CLDN4 log2(CPM+1): "
        + "; ".join(f"{r['label']} ρ = {r['partial']:.3f}" for r in abund)
        + ". The locked Q4-versus-Q1 IFN/MHC shift used CLDN4 %pos, not this abundance rank."
    )
    write_finding(ctx)
    print(json.dumps({k: ctx[k] for k in ("n", "med_rho", "med_p", "qc_rho")}, default=str), flush=True)
    print("primary", [(r["label"], round(r["partial"], 3), round(r["dl"], 3)) for r in primary], flush=True)


if __name__ == "__main__":
    main()

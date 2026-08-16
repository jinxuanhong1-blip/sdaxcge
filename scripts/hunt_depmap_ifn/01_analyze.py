#!/usr/bin/env python3
"""DepMap 24Q4 lung lines: TACSTD2 vs IFN/MHC-I/APM and vs CLDN4.

Tests the user claim "TROP2 high suppresses IFN until CLDN4 lost" as
association + interaction in *cell-intrinsic basal* RNA (no immune cells).
Does not claim causation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

# --- curated gene lists (symbols) ---
MHC1_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
APM_GENES = [
    "TAP1", "TAP2", "TAPBP", "TAPBPL",
    "PSMB8", "PSMB9", "PSMB10",
    "NLRC5", "ERAP1", "ERAP2",
    "CALR", "CANX", "PDIA3",
    "PSME1", "PSME2", "PSME3",
]
# Thompson et al. JITC 2020 8-gene APMS
THOMPSON_APMS = ["B2M", "CALR", "NLRC5", "PSMB9", "PSME1", "PSME3", "RFX5", "HSP90AB1"]
IFN_SIGNAL = [
    "STAT1", "STAT2", "IRF1", "IRF9",
    "JAK1", "JAK2", "IFNGR1", "IFNGR2", "IFNAR1", "IFNAR2",
]
EPITHELIAL = ["CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN7"]
MESENCHYMAL = ["VIM", "ZEB1", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2"]
EXTRA_GENES = ["TACSTD2", "CLDN4", "CD274", "NLRC5", "STAT1", "IRF1"]
HALLMARK_KEEP = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
]
MUT_GENES = ["STK11", "KEAP1", "EGFR", "KRAS", "TP53", "SMARCA4", "B2M"]

NSCLC_TYPES = {
    "LUAD", "LUSC", "LCLC", "NSCLC", "LUAS", "GCLC", "SMARCA4-UT",
    "NUTCL", "LUMEC", "NSCLCPD",
}
SCLC_TYPES = {"SCLC"}
EXCLUDE_FROM_MALIGNANT = {"ZIMMEPCL", "ZIMMLUNG", "ZIMMMPLC"}  # immortalized / non-cancer


def parse_gmt(path: Path) -> dict[str, list[str]]:
    out = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            out[parts[0]] = [g for g in parts[2:] if g]
    return out


def symbol_map(columns: list[str]) -> dict[str, str]:
    """'TACSTD2 (4070)' -> TACSTD2."""
    m = {}
    for c in columns:
        if c == "" or c is None:
            continue
        sym = c.split(" (")[0].strip()
        if sym and sym not in m:
            m[sym] = c
    return m


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    if not cols:
        return pd.Series(np.nan, index=df.index)
    z = df[cols].apply(zscore, axis=0)
    return z.mean(axis=1)


def spearman(x: pd.Series, y: pd.Series) -> dict:
    a = pd.concat([x, y], axis=1).dropna()
    n = len(a)
    if n < 8:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def fisher_z_diff(r1: float, n1: int, r2: float, n2: int) -> dict:
    if min(n1, n2) < 6 or not np.isfinite(r1) or not np.isfinite(r2):
        return {"z": np.nan, "p": np.nan}
    r1 = float(np.clip(r1, -0.999999, 0.999999))
    r2 = float(np.clip(r2, -0.999999, 0.999999))
    z = (np.arctanh(r1) - np.arctanh(r2)) / math.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    p = 2 * stats.norm.sf(abs(z))
    return {"z": float(z), "p": float(p)}


def partial_spearman(x: pd.Series, y: pd.Series, cov: pd.DataFrame) -> dict:
    d = pd.concat([x.rename("x"), y.rename("y"), cov], axis=1).dropna()
    n = len(d)
    k = cov.shape[1]
    if n < 10 + k:
        return {"n": n, "rho": np.nan, "p": np.nan}
    C = np.column_stack([np.ones(n), d.iloc[:, 2:].to_numpy()])
    rx = sm.OLS(d["x"].to_numpy(), C).fit().resid
    ry = sm.OLS(d["y"].to_numpy(), C).fit().resid
    rho, p = stats.spearmanr(rx, ry)
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def ols_interaction(y: pd.Series, a: pd.Series, b: pd.Series, extra: pd.DataFrame | None = None) -> dict:
    d = pd.concat([y.rename("y"), a.rename("a"), b.rename("b")], axis=1)
    if extra is not None:
        d = pd.concat([d, extra], axis=1)
    d = d.dropna()
    d["a"] = zscore(d["a"])
    d["b"] = zscore(d["b"])
    d["ab"] = d["a"] * d["b"]
    X = d[["a", "b", "ab"]].copy()
    if extra is not None:
        extras = [c for c in extra.columns if c in d.columns]
        if extras:
            X = pd.concat([X, pd.get_dummies(d[extras], drop_first=True, dtype=float)], axis=1)
    X = sm.add_constant(X.astype(float))
    fit = sm.OLS(d["y"].astype(float), X).fit(cov_type="HC3")
    out = {
        "n": int(len(d)),
        "r2": float(fit.rsquared),
        "params": {k: float(v) for k, v in fit.params.items()},
        "pvalues": {k: float(v) for k, v in fit.pvalues.items()},
        "bse": {k: float(v) for k, v in fit.bse.items()},
    }
    return out


def mw(a: pd.Series, b: pd.Series) -> dict:
    a = a.dropna()
    b = b.dropna()
    if len(a) < 5 or len(b) < 5:
        return {"n_a": int(len(a)), "n_b": int(len(b)), "median_a": np.nan, "median_b": np.nan, "p": np.nan, "rbc": np.nan}
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = 1 - (2 * u) / (len(a) * len(b))  # rank-biserial; + means a < b
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "p": float(p),
        "rbc": float(rbc),
    }


def subtype_group(depmap_type: str) -> str:
    if depmap_type in {"LUAD"}:
        return "LUAD"
    if depmap_type in {"LUSC"}:
        return "LUSC"
    if depmap_type in SCLC_TYPES:
        return "SCLC"
    if depmap_type in {"LUCA"}:
        return "carcinoid"
    if depmap_type in EXCLUDE_FROM_MALIGNANT:
        return "noncancer"
    if depmap_type in NSCLC_TYPES:
        return "other_NSCLC"
    return "other"


def load_matrix(path: Path, wanted_symbols: list[str]) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    cols = list(header.columns)
    smap = symbol_map(cols)
    use = []
    found = {}
    id_col = cols[0]
    use.append(id_col)
    for s in wanted_symbols:
        if s in smap and smap[s] not in use:
            use.append(smap[s])
            found[s] = smap[s]
    df = pd.read_csv(path, usecols=use)
    df = df.rename(columns={id_col: "ModelID"})
    rename = {v: k for k, v in found.items()}
    df = df.rename(columns=rename)
    df = df.set_index("ModelID")
    return df, found


def load_mut(path: Path, genes: list[str]) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    cols = list(header.columns)
    smap = symbol_map(cols)
    id_col = cols[0]
    use = [id_col]
    found = {}
    for g in genes:
        if g in smap:
            use.append(smap[g])
            found[g] = smap[g]
    df = pd.read_csv(path, usecols=use)
    df = df.rename(columns={id_col: "ModelID"})
    df = df.rename(columns={v: k for k, v in found.items()})
    df = df.set_index("ModelID")
    return (df.fillna(0) > 0).astype(int)


def bh(pvals: list[float]) -> list[float]:
    mask = [np.isfinite(p) for p in pvals]
    q = np.full(len(pvals), np.nan)
    if sum(mask) == 0:
        return list(q)
    _, adj, _, _ = multipletests([pvals[i] for i, m in enumerate(mask) if m], method="fdr_bh")
    j = 0
    for i, m in enumerate(mask):
        if m:
            q[i] = adj[j]
            j += 1
    return [float(x) if np.isfinite(x) else np.nan for x in q]


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, float_format="%.6g")


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


COLORS = {
    "LUAD": "#1f77b4",
    "LUSC": "#ff7f0e",
    "other_NSCLC": "#2ca02c",
    "SCLC": "#d62728",
    "carcinoid": "#9467bd",
    "noncancer": "#7f7f7f",
    "other": "#8c564b",
}


def scatter_colored(ax, x, y, groups, xlabel, ylabel, title):
    for g, col in COLORS.items():
        m = groups == g
        if m.any():
            ax.scatter(x[m], y[m], s=18, alpha=0.75, c=col, label=g, edgecolors="none")
    ok = pd.concat([x, y], axis=1).dropna()
    if len(ok) >= 8:
        rho, p = stats.spearmanr(ok.iloc[:, 0], ok.iloc[:, 1])
        ax.set_title(f"{title}\nρ={rho:.2f} p={p:.2e} n={len(ok)}", fontsize=9)
    else:
        ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    style_ax(ax)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/hunt_depmap_ifn")
    ap.add_argument("--outdir", default="results/hunt_depmap_ifn")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    figdir = out / "figures"
    tabdir = out / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    gmt = parse_gmt(data / "h.all.v2024.1.Hs.symbols.gmt")
    hallmark_genes = sorted({g for k in HALLMARK_KEEP for g in gmt.get(k, [])})
    wanted = sorted(
        set(MHC1_GENES + APM_GENES + THOMPSON_APMS + IFN_SIGNAL + EPITHELIAL + MESENCHYMAL + EXTRA_GENES + hallmark_genes)
    )

    expr, found = load_matrix(data / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", wanted)
    model = pd.read_csv(data / "Model.csv")
    lung_ann = model[model["OncotreeLineage"].fillna("").eq("Lung")].copy()
    lung_ann["subtype_group"] = lung_ann["DepmapModelType"].map(subtype_group)
    lung_ann = lung_ann.set_index("ModelID")

    ids = lung_ann.index.intersection(expr.index)
    expr = expr.loc[ids]
    ann = lung_ann.loc[ids].copy()

    hs = load_mut(data / "OmicsSomaticMutationsMatrixHotspot.csv", MUT_GENES)
    dg = load_mut(data / "OmicsSomaticMutationsMatrixDamaging.csv", MUT_GENES)
    mut = pd.DataFrame(index=ids)
    for g in MUT_GENES:
        h = hs[g].reindex(ids).fillna(0) if g in hs.columns else 0
        d = dg[g].reindex(ids).fillna(0) if g in dg.columns else 0
        mut[f"mut_{g}"] = ((h + d) > 0).astype(int)
        mut[f"mut_{g}_hotspot"] = (hs[g].reindex(ids).fillna(0) > 0).astype(int) if g in hs.columns else 0
        mut[f"mut_{g}_damaging"] = (dg[g].reindex(ids).fillna(0) > 0).astype(int) if g in dg.columns else 0

    # signatures computed within each analysis cohort (z across that cohort)
    def add_scores(frame: pd.DataFrame) -> pd.DataFrame:
        s = frame.copy()
        def cols(genes):
            return [g for g in genes if g in s.columns]

        s["sig_MHC1"] = signature(s, cols(MHC1_GENES))
        s["sig_APM"] = signature(s, cols(APM_GENES))
        s["sig_APMS8"] = signature(s, cols(THOMPSON_APMS))
        s["sig_IFN_signal"] = signature(s, cols(IFN_SIGNAL))
        s["sig_IFNG"] = signature(s, cols(gmt.get("HALLMARK_INTERFERON_GAMMA_RESPONSE", [])))
        s["sig_IFNA"] = signature(s, cols(gmt.get("HALLMARK_INTERFERON_ALPHA_RESPONSE", [])))
        s["sig_infl"] = signature(s, cols(gmt.get("HALLMARK_INFLAMMATORY_RESPONSE", [])))
        s["sig_tnfa"] = signature(s, cols(gmt.get("HALLMARK_TNFA_SIGNALING_VIA_NFKB", [])))
        s["sig_EMT_hm"] = signature(s, cols(gmt.get("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", [])))
        s["sig_epi"] = signature(s, cols(EPITHELIAL))
        s["sig_mes"] = signature(s, cols(MESENCHYMAL))
        s["sig_EMT"] = s["sig_mes"] - s["sig_epi"]
        return s

    coverage_rows = []
    for g in wanted:
        coverage_rows.append({
            "symbol": g,
            "in_matrix": g in found,
            "depmap_col": found.get(g, ""),
            "n_nonNA_lung": int(expr[g].notna().sum()) if g in expr.columns else 0,
        })
    cov_df = pd.DataFrame(coverage_rows)
    save_tsv(cov_df, tabdir / "gene_coverage.tsv")

    sig_cov = []
    for name, genes in [
        ("MHC1", MHC1_GENES),
        ("APM", APM_GENES),
        ("APMS8_Thompson", THOMPSON_APMS),
        ("IFN_signal", IFN_SIGNAL),
        ("HALLMARK_IFNG", gmt.get("HALLMARK_INTERFERON_GAMMA_RESPONSE", [])),
        ("HALLMARK_IFNA", gmt.get("HALLMARK_INTERFERON_ALPHA_RESPONSE", [])),
        ("epithelial", EPITHELIAL),
        ("mesenchymal", MESENCHYMAL),
    ]:
        have = [g for g in genes if g in expr.columns]
        sig_cov.append({"signature": name, "n_requested": len(genes), "n_found": len(have), "genes_found": ",".join(have)})
    save_tsv(pd.DataFrame(sig_cov), tabdir / "signature_coverage.tsv")

    cohorts = {
        "NSCLC": ann["subtype_group"].isin(["LUAD", "LUSC", "other_NSCLC"]),
        "LUAD": ann["subtype_group"].eq("LUAD"),
        "LUSC": ann["subtype_group"].eq("LUSC"),
        "SCLC": ann["subtype_group"].eq("SCLC"),
        "all_malignant_lung": ~ann["subtype_group"].isin(["noncancer"]),
        "all_lung_with_expr": pd.Series(True, index=ann.index),
    }

    SIGS = [
        "sig_IFNG", "sig_IFNA", "sig_IFN_signal", "sig_MHC1", "sig_APM", "sig_APMS8",
        "sig_infl", "sig_tnfa", "sig_EMT", "sig_EMT_hm",
    ]
    PRIMARY_SIGS = ["sig_IFNG", "sig_IFNA", "sig_IFN_signal", "sig_MHC1", "sig_APM"]

    spearman_rows = []
    gene_corr_rows = []
    interact_rows = []
    strat_rows = []
    quad_rows = []
    partial_rows = []
    mut_rows = []
    mw_rows = []
    key = {"release": "DepMap Public 24Q4", "doi": "10.25452/figshare.plus.27993248.v1"}

    sample_tables = []

    for cname, mask in cohorts.items():
        sub_ann = ann.loc[mask]
        sub = add_scores(expr.loc[sub_ann.index])
        sub = sub.join(sub_ann[["CellLineName", "StrippedCellLineName", "DepmapModelType",
                                 "OncotreePrimaryDisease", "OncotreeSubtype", "subtype_group", "Sex"]])
        sub = sub.join(mut.reindex(sub.index))
        sub["cohort"] = cname
        mut_cols = [f"mut_{g}" for g in MUT_GENES]
        keep_cols = [
            "cohort", "CellLineName", "StrippedCellLineName", "DepmapModelType",
            "OncotreeSubtype", "subtype_group", "Sex",
            "TACSTD2", "CLDN4",
        ] + SIGS + mut_cols
        sample_tables.append(sub[[c for c in keep_cols if c in sub.columns]].reset_index())

        # TACSTD2 / CLDN4 vs signatures
        for pred in ["TACSTD2", "CLDN4"]:
            for sig in SIGS + ["CLDN4"] if pred == "TACSTD2" else SIGS:
                if pred == sig:
                    continue
                r = spearman(sub[pred], sub[sig] if sig != "CLDN4" else sub["CLDN4"])
                spearman_rows.append({
                    "cohort": cname, "predictor": pred, "target": sig if sig != "CLDN4" else "CLDN4",
                    **r,
                })

        # individual IFN/MHC/APM genes vs TACSTD2 / CLDN4
        for pred in ["TACSTD2", "CLDN4"]:
            for g in MHC1_GENES + APM_GENES + IFN_SIGNAL:
                if g in sub.columns:
                    r = spearman(sub[pred], sub[g])
                    gene_corr_rows.append({"cohort": cname, "predictor": pred, "gene": g, **r})

        # interaction + stratified (need TACSTD2 and CLDN4)
        if "TACSTD2" in sub.columns and "CLDN4" in sub.columns:
            extra = sub[["subtype_group"]] if cname in {"NSCLC", "all_malignant_lung", "all_lung_with_expr"} else None
            for sig in PRIMARY_SIGS:
                inter = ols_interaction(sub[sig], sub["TACSTD2"], sub["CLDN4"], extra)
                interact_rows.append({
                    "cohort": cname, "target": sig, "n": inter["n"], "r2": inter["r2"],
                    "beta_TACSTD2": inter["params"].get("a", np.nan),
                    "p_TACSTD2": inter["pvalues"].get("a", np.nan),
                    "beta_CLDN4": inter["params"].get("b", np.nan),
                    "p_CLDN4": inter["pvalues"].get("b", np.nan),
                    "beta_interaction": inter["params"].get("ab", np.nan),
                    "p_interaction": inter["pvalues"].get("ab", np.nan),
                    "se_interaction": inter["bse"].get("ab", np.nan),
                })

            med_c = sub["CLDN4"].median()
            hi = sub["CLDN4"] >= med_c
            lo = ~hi
            for sig in PRIMARY_SIGS:
                r_hi = spearman(sub.loc[hi, "TACSTD2"], sub.loc[hi, sig])
                r_lo = spearman(sub.loc[lo, "TACSTD2"], sub.loc[lo, sig])
                diff = fisher_z_diff(r_hi["rho"], r_hi["n"], r_lo["rho"], r_lo["n"])
                strat_rows.append({
                    "cohort": cname, "target": sig, "split": "CLDN4_median",
                    "n_high": r_hi["n"], "rho_high": r_hi["rho"], "p_high": r_hi["p"],
                    "n_low": r_lo["n"], "rho_low": r_lo["rho"], "p_low": r_lo["p"],
                    "cldn4_median": float(med_c),
                    "fisher_z": diff["z"], "fisher_p": diff["p"],
                })

            # also TACSTD2-median split, CLDN4 vs IFN (symmetry check)
            med_t = sub["TACSTD2"].median()
            thi = sub["TACSTD2"] >= med_t
            tlo = ~thi
            for sig in PRIMARY_SIGS:
                r_hi = spearman(sub.loc[thi, "CLDN4"], sub.loc[thi, sig])
                r_lo = spearman(sub.loc[tlo, "CLDN4"], sub.loc[tlo, sig])
                diff = fisher_z_diff(r_hi["rho"], r_hi["n"], r_lo["rho"], r_lo["n"])
                strat_rows.append({
                    "cohort": cname, "target": sig, "split": "TACSTD2_median",
                    "n_high": r_hi["n"], "rho_high": r_hi["rho"], "p_high": r_hi["p"],
                    "n_low": r_lo["n"], "rho_low": r_lo["rho"], "p_low": r_lo["p"],
                    "cldn4_median": float(med_t),
                    "fisher_z": diff["z"], "fisher_p": diff["p"],
                })

            q_th = (sub["TACSTD2"] >= med_t)
            q_ch = (sub["CLDN4"] >= med_c)
            labels = pd.Series(index=sub.index, dtype=object)
            labels[q_th & q_ch] = "TROP2hi_CLDN4hi"
            labels[q_th & ~q_ch] = "TROP2hi_CLDN4lo"
            labels[~q_th & q_ch] = "TROP2lo_CLDN4hi"
            labels[~q_th & ~q_ch] = "TROP2lo_CLDN4lo"
            sub["_quad"] = labels
            for sig in PRIMARY_SIGS:
                meds = {q: float(sub.loc[labels == q, sig].median()) for q in labels.dropna().unique()}
                # key contrast: TROP2hi/CLDN4hi vs TROP2hi/CLDN4lo  (claim: IFN rises when CLDN4 lost)
                contrast = mw(sub.loc[labels == "TROP2hi_CLDN4hi", sig], sub.loc[labels == "TROP2hi_CLDN4lo", sig])
                # also TROP2hi/CLDN4hi vs TROP2lo/CLDN4hi
                contrast2 = mw(sub.loc[labels == "TROP2hi_CLDN4hi", sig], sub.loc[labels == "TROP2lo_CLDN4hi", sig])
                # overall Kruskal
                groups = [sub.loc[labels == q, sig].dropna() for q in
                          ["TROP2hi_CLDN4hi", "TROP2hi_CLDN4lo", "TROP2lo_CLDN4hi", "TROP2lo_CLDN4lo"]]
                groups = [g for g in groups if len(g)]
                if len(groups) >= 3:
                    kstat, kp = stats.kruskal(*groups)
                else:
                    kstat, kp = np.nan, np.nan
                quad_rows.append({
                    "cohort": cname, "target": sig,
                    "n_TROP2hi_CLDN4hi": int((labels == "TROP2hi_CLDN4hi").sum()),
                    "n_TROP2hi_CLDN4lo": int((labels == "TROP2hi_CLDN4lo").sum()),
                    "n_TROP2lo_CLDN4hi": int((labels == "TROP2lo_CLDN4hi").sum()),
                    "n_TROP2lo_CLDN4lo": int((labels == "TROP2lo_CLDN4lo").sum()),
                    "median_TROP2hi_CLDN4hi": meds.get("TROP2hi_CLDN4hi", np.nan),
                    "median_TROP2hi_CLDN4lo": meds.get("TROP2hi_CLDN4lo", np.nan),
                    "median_TROP2lo_CLDN4hi": meds.get("TROP2lo_CLDN4hi", np.nan),
                    "median_TROP2lo_CLDN4lo": meds.get("TROP2lo_CLDN4lo", np.nan),
                    "kruskal_H": float(kstat) if np.isfinite(kstat) else np.nan,
                    "kruskal_p": float(kp) if np.isfinite(kp) else np.nan,
                    "MW_hihi_vs_hilo_p": contrast["p"],
                    "MW_hihi_vs_hilo_rbc": contrast["rbc"],
                    "MW_hihi_vs_lohi_p": contrast2["p"],
                    "MW_hihi_vs_lohi_rbc": contrast2["rbc"],
                })

            # partial correlations
            for sig in PRIMARY_SIGS:
                p1 = partial_spearman(sub["TACSTD2"], sub[sig], sub[["CLDN4"]])
                p2 = partial_spearman(sub["CLDN4"], sub[sig], sub[["TACSTD2"]])
                p3 = partial_spearman(sub["TACSTD2"], sub[sig], sub[["CLDN4", "sig_EMT"]])
                p4 = partial_spearman(sub["TACSTD2"], sub[sig], sub[["sig_EMT"]])
                partial_rows.append({"cohort": cname, "predictor": "TACSTD2", "target": sig,
                                     "adjusted_for": "CLDN4", **p1})
                partial_rows.append({"cohort": cname, "predictor": "CLDN4", "target": sig,
                                     "adjusted_for": "TACSTD2", **p2})
                partial_rows.append({"cohort": cname, "predictor": "TACSTD2", "target": sig,
                                     "adjusted_for": "CLDN4+EMT", **p3})
                partial_rows.append({"cohort": cname, "predictor": "TACSTD2", "target": sig,
                                     "adjusted_for": "EMT", **p4})
                if cname in {"NSCLC", "LUAD"} and "mut_STK11" in sub.columns:
                    p5 = partial_spearman(sub["TACSTD2"], sub[sig], sub[["CLDN4", "mut_STK11"]])
                    partial_rows.append({"cohort": cname, "predictor": "TACSTD2", "target": sig,
                                         "adjusted_for": "CLDN4+STK11", **p5})

        # mutation vs signatures / TACSTD2
        for g in MUT_GENES:
            mg = f"mut_{g}"
            if mg not in sub.columns:
                continue
            for target in ["TACSTD2", "CLDN4"] + PRIMARY_SIGS:
                a = sub.loc[sub[mg] == 1, target]
                b = sub.loc[sub[mg] == 0, target]
                r = mw(a, b)
                mut_rows.append({
                    "cohort": cname, "mutation": g, "n_mut": r["n_a"], "n_wt": r["n_b"],
                    "target": target, "median_mut": r["median_a"], "median_wt": r["median_b"],
                    "p": r["p"], "rbc_mut_minus_wt": r["rbc"],
                })

        # median-split TACSTD2 / CLDN4 vs signatures
        for pred in ["TACSTD2", "CLDN4"]:
            med = sub[pred].median()
            for sig in PRIMARY_SIGS:
                r = mw(sub.loc[sub[pred] >= med, sig], sub.loc[sub[pred] < med, sig])
                mw_rows.append({
                    "cohort": cname, "predictor": pred, "target": sig,
                    "n_high": r["n_a"], "n_low": r["n_b"],
                    "median_high": r["median_a"], "median_low": r["median_b"],
                    "p": r["p"], "rbc_high_minus_low": r["rbc"],
                })

        if cname == "NSCLC":
            key["n_NSCLC"] = int(len(sub))
            key["n_NSCLC_by_group"] = sub["subtype_group"].value_counts().to_dict()
            key["TACSTD2_CLDN4"] = spearman(sub["TACSTD2"], sub["CLDN4"])
            key["NSCLC_univariate"] = {}
            for pred in ["TACSTD2", "CLDN4"]:
                key["NSCLC_univariate"][pred] = {sig: spearman(sub[pred], sub[sig]) for sig in PRIMARY_SIGS}
            key["NSCLC_CLDN4_median"] = float(sub["CLDN4"].median())
            key["NSCLC_TACSTD2_median"] = float(sub["TACSTD2"].median())
            key["NSCLC_n_mut"] = {g: int(sub[f"mut_{g}"].sum()) for g in MUT_GENES if f"mut_{g}" in sub.columns}

    samples = pd.concat(sample_tables, ignore_index=True)
    save_tsv(samples, tabdir / "sample_table.tsv")

    sp = pd.DataFrame(spearman_rows)
    # FDR within cohort × predictor for primary signatures
    sp["q_primary"] = np.nan
    for (coh, pred), idx in sp.groupby(["cohort", "predictor"]).groups.items():
        ii = list(idx)
        sub = sp.loc[ii]
        prim = sub["target"].isin(PRIMARY_SIGS + (["CLDN4"] if pred == "TACSTD2" else []))
        qs = bh(sub.loc[prim, "p"].tolist())
        sp.loc[sub.index[prim], "q_primary"] = qs
    save_tsv(sp, tabdir / "spearman.tsv")

    gc = pd.DataFrame(gene_corr_rows)
    gc["q"] = np.nan
    for (coh, pred), idx in gc.groupby(["cohort", "predictor"]).groups.items():
        ii = list(idx)
        gc.loc[ii, "q"] = bh(gc.loc[ii, "p"].tolist())
    save_tsv(gc, tabdir / "gene_level_spearman.tsv")

    inter_df = pd.DataFrame(interact_rows)
    if len(inter_df):
        inter_df["q_interaction"] = np.nan
        for coh, idx in inter_df.groupby("cohort").groups.items():
            inter_df.loc[list(idx), "q_interaction"] = bh(inter_df.loc[list(idx), "p_interaction"].tolist())
    save_tsv(inter_df, tabdir / "interaction_ols.tsv")

    st = pd.DataFrame(strat_rows)
    save_tsv(st, tabdir / "stratified_spearman.tsv")
    qd = pd.DataFrame(quad_rows)
    save_tsv(qd, tabdir / "quadrants.tsv")
    pr = pd.DataFrame(partial_rows)
    save_tsv(pr, tabdir / "partial_spearman.tsv")
    mu = pd.DataFrame(mut_rows)
    if len(mu):
        mu["q"] = np.nan
        for (coh, mutg), idx in mu.groupby(["cohort", "mutation"]).groups.items():
            mu.loc[list(idx), "q"] = bh(mu.loc[list(idx), "p"].tolist())
    save_tsv(mu, tabdir / "mutation_associations.tsv")
    save_tsv(pd.DataFrame(mw_rows), tabdir / "median_split.tsv")

    # cohort counts
    counts = []
    for cname, mask in cohorts.items():
        sub = ann.loc[mask]
        counts.append({
            "cohort": cname,
            "n_annot_lung": int(mask.sum()),
            "n_with_expr": int(len(sub)),
            **{f"n_{g}": int((sub["subtype_group"] == g).sum()) for g in
               ["LUAD", "LUSC", "other_NSCLC", "SCLC", "carcinoid", "noncancer", "other"]},
        })
    save_tsv(pd.DataFrame(counts), tabdir / "cohort_counts.tsv")

    # ---------- figures (primary NSCLC) ----------
    ns_mask = cohorts["NSCLC"]
    ns = add_scores(expr.loc[ann.index[ns_mask]])
    ns = ns.join(ann.loc[ns.index, ["subtype_group", "CellLineName", "DepmapModelType"]])
    ns = ns.join(mut.reindex(ns.index))

    # fig1 TACSTD2 vs CLDN4
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    scatter_colored(ax, ns["TACSTD2"], ns["CLDN4"], ns["subtype_group"],
                    "TACSTD2 log2(TPM+1)", "CLDN4 log2(TPM+1)", "NSCLC lines: TACSTD2 vs CLDN4")
    ax.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(figdir / "fig1_tacstd2_vs_cldn4.png", dpi=160)
    plt.close(fig)

    # fig2 TACSTD2 vs 4 signatures
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2))
    for ax, sig, lab in zip(
        axes.ravel(),
        ["sig_IFNG", "sig_IFNA", "sig_MHC1", "sig_APM"],
        ["Hallmark IFN-γ", "Hallmark IFN-α", "MHC-I (HLA-A/B/C+B2M)", "APM (16 genes)"],
    ):
        scatter_colored(ax, ns["TACSTD2"], ns[sig], ns["subtype_group"],
                        "TACSTD2 log2(TPM+1)", lab, lab)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=7, frameon=False)
    fig.suptitle("NSCLC: TACSTD2 vs IFN / MHC-I / APM (basal RNA)", fontsize=11)
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    fig.savefig(figdir / "fig2_tacstd2_vs_signatures.png", dpi=160)
    plt.close(fig)

    # fig3 heatmap of selected variables
    heat_vars = ["TACSTD2", "CLDN4", "sig_IFNG", "sig_IFNA", "sig_IFN_signal",
                 "sig_MHC1", "sig_APM", "sig_APMS8", "sig_infl", "sig_tnfa",
                 "sig_EMT", "CDH1", "EPCAM", "VIM", "STAT1", "NLRC5", "HLA-A", "B2M"]
    heat_vars = [v for v in heat_vars if v in ns.columns]
    cm = ns[heat_vars].corr(method="spearman")
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    im = ax.imshow(cm.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(heat_vars)))
    ax.set_yticks(range(len(heat_vars)))
    ax.set_xticklabels(heat_vars, rotation=75, ha="right", fontsize=7)
    ax.set_yticklabels(heat_vars, fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman ρ")
    ax.set_title("NSCLC Spearman (expression / signatures)")
    fig.tight_layout()
    fig.savefig(figdir / "fig3_correlation_heatmap.png", dpi=160)
    plt.close(fig)
    cm.reset_index().rename(columns={"index": "var"}).to_csv(tabdir / "nsclc_spearman_matrix.tsv", sep="\t", index=False)

    # fig4 CLDN4-stratified
    med_c = ns["CLDN4"].median()
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), sharey=True)
    for ax, sel, title in [
        (axes[0], ns["CLDN4"] >= med_c, f"CLDN4 high (≥ median {med_c:.2f})"),
        (axes[1], ns["CLDN4"] < med_c, f"CLDN4 low (< median {med_c:.2f})"),
    ]:
        scatter_colored(ax, ns.loc[sel, "TACSTD2"], ns.loc[sel, "sig_IFNG"],
                        ns.loc[sel, "subtype_group"],
                        "TACSTD2 log2(TPM+1)", "Hallmark IFN-γ z-mean", title)
    fig.suptitle("NSCLC: TACSTD2 vs IFN-γ stratified by CLDN4 (median)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig4_cldn4_stratified_ifng.png", dpi=160)
    plt.close(fig)

    # fig5 quadrants
    med_t = ns["TACSTD2"].median()
    q = pd.Series("TROP2lo_CLDN4lo", index=ns.index)
    q[(ns["TACSTD2"] >= med_t) & (ns["CLDN4"] >= med_c)] = "TROP2hi_CLDN4hi"
    q[(ns["TACSTD2"] >= med_t) & (ns["CLDN4"] < med_c)] = "TROP2hi_CLDN4lo"
    q[(ns["TACSTD2"] < med_t) & (ns["CLDN4"] >= med_c)] = "TROP2lo_CLDN4hi"
    order = ["TROP2hi_CLDN4hi", "TROP2hi_CLDN4lo", "TROP2lo_CLDN4hi", "TROP2lo_CLDN4lo"]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2))
    for ax, sig, lab in [(axes[0], "sig_IFNG", "Hallmark IFN-γ"), (axes[1], "sig_MHC1", "MHC-I")]:
        data = [ns.loc[q == name, sig].dropna() for name in order]
        bp = ax.boxplot(data, tick_labels=["T2hi C4hi", "T2hi C4lo", "T2lo C4hi", "T2lo C4lo"],
                        patch_artist=True, widths=0.6)
        cols = ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax.set_ylabel(lab + " z-mean", fontsize=8)
        ax.set_title(lab, fontsize=9)
        style_ax(ax)
        ax.tick_params(axis="x", labelsize=7)
    fig.suptitle("NSCLC median-split quadrants (TROP2 × CLDN4)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig5_quadrants.png", dpi=160)
    plt.close(fig)

    # fig6 CLDN4 vs signatures (for honesty: is CLDN4 the stronger correlate?)
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2))
    for ax, sig, lab in zip(
        axes.ravel(),
        ["sig_IFNG", "sig_IFNA", "sig_MHC1", "sig_APM"],
        ["Hallmark IFN-γ", "Hallmark IFN-α", "MHC-I", "APM"],
    ):
        scatter_colored(ax, ns["CLDN4"], ns[sig], ns["subtype_group"],
                        "CLDN4 log2(TPM+1)", lab, lab)
    fig.suptitle("NSCLC: CLDN4 vs IFN / MHC-I / APM (same lines)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig6_cldn4_vs_signatures.png", dpi=160)
    plt.close(fig)

    # fig7 forest of univariate rhos across cohorts
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    plot_coh = ["NSCLC", "LUAD", "LUSC", "SCLC"]
    plot_sig = ["sig_IFNG", "sig_MHC1", "sig_APM"]
    y = 0
    yticks = []
    ylabels = []
    for sig in plot_sig:
        for coh in plot_coh:
            row = sp[(sp.cohort == coh) & (sp.predictor == "TACSTD2") & (sp.target == sig)]
            if row.empty:
                continue
            rho = row.iloc[0]["rho"]
            n = row.iloc[0]["n"]
            # approx 95% CI via Fisher
            if n > 4 and np.isfinite(rho):
                z = np.arctanh(np.clip(rho, -0.999, 0.999))
                se = 1 / math.sqrt(n - 3)
                lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
            else:
                lo, hi = np.nan, np.nan
            color = {"NSCLC": "#333333", "LUAD": "#1f77b4", "LUSC": "#ff7f0e", "SCLC": "#d62728"}[coh]
            ax.plot([lo, hi], [y, y], color=color, lw=1.4)
            ax.plot(rho, y, "o", color=color, ms=5)
            yticks.append(y)
            ylabels.append(f"{sig.replace('sig_','')}  {coh} n={n}")
            y -= 1
        y -= 0.4
    ax.axvline(0, color="#999", lw=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Spearman ρ  TACSTD2 vs signature (Fisher 95% CI)")
    ax.set_title("TACSTD2–signature correlations by cohort")
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(figdir / "fig7_tacstd2_forest.png", dpi=160)
    plt.close(fig)

    # fig8 distributions
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.6))
    for ax, col, lab in [
        (axes[0], "TACSTD2", "TACSTD2"),
        (axes[1], "CLDN4", "CLDN4"),
        (axes[2], "sig_IFNG", "Hallmark IFN-γ"),
    ]:
        for g, colr in COLORS.items():
            v = ns.loc[ns["subtype_group"] == g, col]
            if len(v):
                ax.hist(v, bins=18, alpha=0.45, color=colr, label=g)
        ax.set_title(lab, fontsize=9)
        ax.set_xlabel(lab, fontsize=8)
        style_ax(ax)
    axes[0].legend(fontsize=6, frameon=False)
    fig.suptitle("NSCLC distributions", fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / "fig8_distributions.png", dpi=160)
    plt.close(fig)

    # pack key stats for writeup
    def rowget(df, **kw):
        q = df
        for k, v in kw.items():
            q = q[q[k] == v]
        return q.iloc[0].to_dict() if len(q) else {}

    key["primary_claim_tests"] = {
        "TACSTD2_vs_CLDN4_NSCLC": rowget(sp, cohort="NSCLC", predictor="TACSTD2", target="CLDN4"),
        "TACSTD2_vs_IFNG_NSCLC": rowget(sp, cohort="NSCLC", predictor="TACSTD2", target="sig_IFNG"),
        "TACSTD2_vs_MHC1_NSCLC": rowget(sp, cohort="NSCLC", predictor="TACSTD2", target="sig_MHC1"),
        "TACSTD2_vs_APM_NSCLC": rowget(sp, cohort="NSCLC", predictor="TACSTD2", target="sig_APM"),
        "CLDN4_vs_IFNG_NSCLC": rowget(sp, cohort="NSCLC", predictor="CLDN4", target="sig_IFNG"),
        "interaction_IFNG_NSCLC": rowget(inter_df, cohort="NSCLC", target="sig_IFNG"),
        "strat_IFNG_CLDN4split_NSCLC": rowget(st[(st.split == "CLDN4_median")], cohort="NSCLC", target="sig_IFNG")
        if len(st) else {},
        "quad_IFNG_NSCLC": rowget(qd, cohort="NSCLC", target="sig_IFNG"),
        "partial_TACSTD2_IFNG_adj_CLDN4": rowget(pr, cohort="NSCLC", predictor="TACSTD2",
                                                 target="sig_IFNG", adjusted_for="CLDN4"),
        "partial_CLDN4_IFNG_adj_TACSTD2": rowget(pr, cohort="NSCLC", predictor="CLDN4",
                                                 target="sig_IFNG", adjusted_for="TACSTD2"),
        "partial_TACSTD2_IFNG_adj_CLDN4_EMT": rowget(pr, cohort="NSCLC", predictor="TACSTD2",
                                                     target="sig_IFNG", adjusted_for="CLDN4+EMT"),
        "LUAD_TACSTD2_IFNG": rowget(sp, cohort="LUAD", predictor="TACSTD2", target="sig_IFNG"),
        "LUSC_TACSTD2_IFNG": rowget(sp, cohort="LUSC", predictor="TACSTD2", target="sig_IFNG"),
        "SCLC_TACSTD2_IFNG": rowget(sp, cohort="SCLC", predictor="TACSTD2", target="sig_IFNG"),
    }
    # fix strat rowget — st filtered may break rowget
    st_ifng = st[(st.cohort == "NSCLC") & (st.target == "sig_IFNG") & (st.split == "CLDN4_median")]
    if len(st_ifng):
        key["primary_claim_tests"]["strat_IFNG_CLDN4split_NSCLC"] = st_ifng.iloc[0].to_dict()

    key["n_lung_annot"] = int((model["OncotreeLineage"].fillna("").eq("Lung")).sum())
    key["n_lung_expr"] = int(len(ann))
    key["missing_genes"] = [g for g in wanted if g not in found]
    key["notes"] = [
        "Expression is log2(TPM+1) protein-coding matrix, DepMap 24Q4.",
        "Signatures are mean of within-cohort z-scored genes.",
        "No immune infiltrate: IFN/MHC/APM are cancer-cell basal transcription.",
        "CLDN4 'loss' is low RNA, not a validated genetic deletion call.",
        "Hotspot OR damaging mutation = mutant for STK11/KEAP1/EGFR/KRAS/TP53/SMARCA4/B2M.",
    ]

    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2, default=str))

    # machine-readable verdict helpers
    print("NSCLC n", key["n_NSCLC"])
    print("TACSTD2 vs CLDN4", key["TACSTD2_CLDN4"])
    print("TACSTD2 vs IFNG", key["primary_claim_tests"]["TACSTD2_vs_IFNG_NSCLC"])
    print("CLDN4 vs IFNG", key["primary_claim_tests"]["CLDN4_vs_IFNG_NSCLC"])
    print("interaction IFNG", key["primary_claim_tests"]["interaction_IFNG_NSCLC"])
    print("strat", key["primary_claim_tests"]["strat_IFNG_CLDN4split_NSCLC"])
    print("quad", key["primary_claim_tests"]["quad_IFNG_NSCLC"])
    print("partial T|C", key["primary_claim_tests"]["partial_TACSTD2_IFNG_adj_CLDN4"])
    print("partial C|T", key["primary_claim_tests"]["partial_CLDN4_IFNG_adj_TACSTD2"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

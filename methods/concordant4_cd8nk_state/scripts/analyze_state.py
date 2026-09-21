#!/usr/bin/env python3
"""Patient-level CD8/NK effector and exhaustion vs locked malignant CLDN4.

Primary question: in concordant-4, does effector intensity (and exhaustion)
inside CD8 cells and inside NK cells differ between patients with high vs low
malignant CLDN4? This is the composition analog of the CosMx neighbor-state
result, not a spatial test and not a re-derivation of the T/NK fraction result.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]

EFFECTOR = ["GZMB", "PRF1", "NKG7", "IFNG"]
EFFECTOR_NO_NKG7 = ["GZMB", "PRF1", "IFNG"]
EXHAUSTION = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"]
EXHAUSTION_EXT = EXHAUSTION + ["CTLA4", "ENTPD1"]
NAIVE = ["IL7R", "TCF7", "CCR7", "SELL", "LEF1"]

PRIMARY = [
    ("cd8", "ALL", "effector_cosmx"),
    ("nk", "ALL", "effector_cosmx"),
    ("cd8", "ALL", "exhaustion"),
    ("nk", "ALL", "exhaustion"),
]

MODULES = {
    "effector_cosmx": EFFECTOR,
    "effector_no_nkg7": EFFECTOR_NO_NKG7,
    "exhaustion": EXHAUSTION,
    "exhaustion_ext": EXHAUSTION_EXT,
    "naive_memory": NAIVE,
}

COHORT_ORDER = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLOR = {
    "GSE123902": "#0072B2",
    "GSE131907": "#E69F00",
    "GSE205335": "#009E73",
    "GSE189357": "#CC79A7",
}

MIN_N = 10
META_MIN_N = 5


def dl_spearman(rhos, ns) -> dict:
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    if len(rhos) == 0:
        return {"rho": np.nan, "p": np.nan, "I2": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "k": 0, "N": 0}
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var_z = 1.0 / (ns - 3.0)
    w = 1.0 / var_z
    zbar = np.sum(w * z) / np.sum(w)
    q = float(np.sum(w * (z - zbar) ** 2))
    k = int(len(rhos))
    dfree = k - 1
    cdenom = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var_z + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = float(np.sqrt(1.0 / np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zre / se))) if se > 0 else np.nan
    i2 = float(max(0.0, (q - dfree) / q)) if q > 0 else 0.0
    ci = np.tanh(zre + np.array([-1.0, 1.0]) * 1.96 * se)
    return {
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": i2,
        "ci_lo": float(ci[0]),
        "ci_hi": float(ci[1]),
        "k": k,
        "N": int(ns.sum()),
        "se_z": se,
    }


def rank_biserial(x4, x1) -> dict:
    x4 = np.asarray(x4, dtype=float)
    x1 = np.asarray(x1, dtype=float)
    x4 = x4[np.isfinite(x4)]
    x1 = x1[np.isfinite(x1)]
    if len(x4) == 0 or len(x1) == 0:
        return {"r": np.nan, "p": np.nan, "n_q1": len(x1), "n_q4": len(x4)}
    res = stats.mannwhitneyu(x4, x1, alternative="two-sided", method="asymptotic", use_continuity=True)
    n4, n1 = len(x4), len(x1)
    u = float(res.statistic)
    return {"r": 2.0 * u / (n4 * n1) - 1.0, "p": float(res.pvalue), "n_q1": n1, "n_q4": n4, "U": u}


def ols_q4(df: pd.DataFrame) -> dict:
    sub = df[df["quartile"].isin(["Q1", "Q4"])].dropna(subset=["score"]).copy()
    both = sub.groupby("dataset")["quartile"].nunique()
    sub = sub[sub["dataset"].isin(both[both == 2].index)]
    if len(sub) < 6:
        return {"delta": np.nan, "se": np.nan, "p": np.nan, "n": int(len(sub)), "n_q1": 0, "n_q4": 0, "df": 0}
    y = sub["score"].to_numpy(dtype=float)
    q4 = (sub["quartile"].to_numpy() == "Q4").astype(float)
    cohorts = [c for c in COHORT_ORDER if c in set(sub["dataset"])]
    dummies = [(sub["dataset"].to_numpy() == c).astype(float) for c in cohorts[1:]]
    cols = [np.ones(len(sub)), q4, *dummies]
    X = np.column_stack(cols)
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df_res = int(len(y) - rank)
    sigma2 = float(np.sum(resid ** 2) / df_res) if df_res > 0 else np.nan
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)) * sigma2)
    tstat = float(beta[1] / se[1]) if se[1] > 0 else np.nan
    p = float(2 * stats.t.sf(abs(tstat), df_res)) if df_res > 0 and np.isfinite(tstat) else np.nan
    return {
        "delta": float(beta[1]),
        "se": float(se[1]),
        "p": p,
        "n": int(len(sub)),
        "n_q1": int((sub["quartile"] == "Q1").sum()),
        "n_q4": int((sub["quartile"] == "Q4").sum()),
        "df": df_res,
    }


def spearman_safe(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 4 or np.unique(x[ok]).size < 2 or np.unique(y[ok]).size < 2:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[ok], y[ok])
    return float(rho), float(p), n


def calibrate(locked: pd.DataFrame) -> dict:
    rows = []
    for ds, sub in locked.groupby("dataset"):
        rho, p, n = spearman_safe(sub["mal_CLDN4_pct"], sub["frac_tnk"])
        rows.append({"dataset": ds, "rho": rho, "p": p, "n": n})
    meta = dl_spearman([r["rho"] for r in rows], [r["n"] for r in rows])
    rb = rank_biserial(
        locked.loc[locked["cldn4_quartile"] == "Q4", "frac_tnk"],
        locked.loc[locked["cldn4_quartile"] == "Q1", "frac_tnk"],
    )
    out = {"cohorts": rows, "meta": meta, "stacked_q4q1": rb}
    print(
        f"CALIBRATE T/NK fraction meta rho={meta['rho']:.4f} p={meta['p']:.3e} "
        f"I2={meta['I2']:.3f} N={meta['N']} stacked r={rb['r']:.3f} p={rb['p']:.4g} "
        f"n={rb['n_q1']}/{rb['n_q4']}",
        flush=True,
    )
    return out


def load_sums(table_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sum_files = sorted(table_dir.glob("compartment_gene_sums*.tsv"))
    qc_files = sorted(table_dir.glob("malignant_qc*.tsv"))
    if not sum_files:
        raise SystemExit(f"no compartment_gene_sums*.tsv in {table_dir}")
    sums = pd.concat([pd.read_csv(p, sep="\t") for p in sum_files], ignore_index=True)
    qc = pd.concat([pd.read_csv(p, sep="\t") for p in qc_files], ignore_index=True) if qc_files else pd.DataFrame()
    if sums.duplicated(["dataset", "unit_id", "compartment", "slice", "gene"]).any():
        raise SystemExit("duplicate gene-sum keys; remove a combined file that overlaps per-dataset files")
    return sums, qc


def module_table(sums: pd.DataFrame, locked: pd.DataFrame) -> pd.DataFrame:
    """One row per dataset, unit, compartment, slice, module."""
    n_key = ["dataset", "unit_id", "compartment", "slice"]
    meta = (
        sums.groupby(n_key, as_index=False)
        .agg(n_cells=("n_cells", "first"), median_lib=("median_lib", "first"), mean_lib=("mean_lib", "first"))
    )
    wide = sums.pivot_table(index=n_key, columns="gene", values="sum_log1p", aggfunc="sum")
    records = []
    modules = dict(MODULES)
    for gene in sorted(set(EFFECTOR + EXHAUSTION + NAIVE + ["CTLA4", "ENTPD1", "GNLY", "GZMA", "GZMH"])):
        modules[f"gene:{gene}"] = [gene]
    for key, sub_idx in meta.groupby(n_key).groups.items():
        pass
    # iterate meta rows
    wide = wide.reset_index()
    merged = meta.merge(wide, on=n_key, how="left")
    out_rows = []
    for rec in merged.itertuples(index=False):
        n = int(rec.n_cells)
        base = {
            "dataset": rec.dataset,
            "unit_id": rec.unit_id,
            "compartment": rec.compartment,
            "slice": rec.slice,
            "n_cells": n,
            "median_lib": rec.median_lib,
            "mean_lib": rec.mean_lib,
        }
        rowmap = rec._asdict() if hasattr(rec, "_asdict") else None
        # itertuples renames columns; use getattr
        for name, genes in modules.items():
            vals = []
            used = []
            for g in genes:
                if g not in merged.columns:
                    continue
                s = getattr(rec, g)
                if pd.notna(s):
                    vals.append(float(s) / n)
                    used.append(g)
            out_rows.append(
                {
                    **base,
                    "score_name": name,
                    "score": float(np.mean(vals)) if vals else np.nan,
                    "n_genes": len(used),
                    "genes_used": ",".join(used),
                }
            )
    scores = pd.DataFrame(out_rows)
    exposure = locked[["dataset", "unit_id", "mal_CLDN4_pct", "cldn4_quartile", "frac_tnk", "n_cells", "n_tnk"]].rename(
        columns={"n_cells": "n_cells_unit"}
    )
    scores = scores.merge(exposure, on=["dataset", "unit_id"], how="left")
    if scores["mal_CLDN4_pct"].isna().any():
        missing = scores.loc[scores["mal_CLDN4_pct"].isna(), ["dataset", "unit_id"]].drop_duplicates()
        raise SystemExit(f"units missing from locked exposure:\n{missing}")
    return scores


def per_cohort_spearman(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ds in COHORT_ORDER:
        sub = df[df["dataset"] == ds]
        rho, p, n = spearman_safe(sub["mal_CLDN4_pct"], sub["score"])
        q4 = sub.loc[sub["quartile"] == "Q4", "score"]
        q1 = sub.loc[sub["quartile"] == "Q1", "score"]
        dlt = float(np.nanmean(q4) - np.nanmean(q1)) if len(q4) and len(q1) else np.nan
        ratio = float(np.nanmean(q4) / np.nanmean(q1)) if len(q4) and len(q1) and np.nanmean(q1) not in (0, np.nan) else np.nan
        rows.append(
            {
                "dataset": ds,
                "n": n,
                "rho": rho,
                "p": p,
                "n_q1": int(np.isfinite(q1).sum()) if len(q1) else 0,
                "n_q4": int(np.isfinite(q4).sum()) if len(q4) else 0,
                "delta_q4_minus_q1": dlt,
                "ratio_q4_over_q1": ratio,
                "mean_q1": float(np.nanmean(q1)) if len(q1) else np.nan,
                "mean_q4": float(np.nanmean(q4)) if len(q4) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_outcome(df: pd.DataFrame, label: str) -> dict:
    cohorts = per_cohort_spearman(df)
    use = cohorts[(cohorts["n"] >= META_MIN_N) & np.isfinite(cohorts["rho"])]
    meta = dl_spearman(use["rho"], use["n"])
    rb = rank_biserial(df.loc[df["quartile"] == "Q4", "score"], df.loc[df["quartile"] == "Q1", "score"])
    ols = ols_q4(df.rename(columns={"quartile": "quartile"}))
    n_down = int((use["delta_q4_minus_q1"] < 0).sum())
    n_rho_neg = int((use["rho"] < 0).sum())
    return {
        "label": label,
        "cohorts": cohorts.to_dict(orient="records"),
        "meta": meta,
        "stacked_q4q1": rb,
        "ols_q4": ols,
        "k_delta_neg": n_down,
        "k_rho_neg": n_rho_neg,
        "k_tested": int(len(use)),
    }


def eligible(scores: pd.DataFrame, compartment: str, slice_name: str, score_name: str, min_n: int) -> pd.DataFrame:
    sub = scores[
        (scores["compartment"] == compartment)
        & (scores["slice"] == slice_name)
        & (scores["score_name"] == score_name)
        & (scores["n_cells"] >= min_n)
    ].copy()
    sub = sub.rename(columns={"cldn4_quartile": "quartile"})
    return sub


def qc_exposure(qc: pd.DataFrame, locked: pd.DataFrame) -> pd.DataFrame:
    m = qc.merge(
        locked[["dataset", "unit_id", "mal_CLDN4_pct", "n_malignant", "n_cells"]],
        on=["dataset", "unit_id"],
        how="outer",
        suffixes=("_re", "_locked"),
    )
    rows = []
    for ds, sub in m.groupby("dataset"):
        rho, p, n = spearman_safe(sub["mal_CLDN4_pct"], sub["mal_CLDN4_pct_recomputed"])
        both = sub.dropna(subset=["mal_CLDN4_pct", "mal_CLDN4_pct_recomputed"])
        pear = float(np.corrcoef(both["mal_CLDN4_pct"], both["mal_CLDN4_pct_recomputed"])[0, 1]) if len(both) > 2 else np.nan
        med_abs = float((both["mal_CLDN4_pct"] - both["mal_CLDN4_pct_recomputed"]).abs().median()) if len(both) else np.nan
        rows.append({"dataset": ds, "n": n, "spearman": rho, "pearson": pear, "median_abs_pp": med_abs, "p": p})
        print(f"QC CLDN4 {ds}: pearson={pear:.4f} median|Δ|={med_abs:.3f} pp n={len(both)}", flush=True)
    return pd.DataFrame(rows)


def composition_table(scores: pd.DataFrame) -> pd.DataFrame:
    all_rows = scores[(scores["slice"] == "ALL") & (scores["score_name"] == "effector_cosmx")][
        ["dataset", "unit_id", "compartment", "n_cells", "mal_CLDN4_pct", "cldn4_quartile"]
    ]
    wide = all_rows.pivot_table(index=["dataset", "unit_id", "mal_CLDN4_pct", "cldn4_quartile"], columns="compartment", values="n_cells")
    wide = wide.reset_index().rename(columns={"cd8": "n_cd8", "nk": "n_nk"})
    for col in ("n_cd8", "n_nk"):
        if col not in wide.columns:
            wide[col] = 0
    wide["n_cd8"] = wide["n_cd8"].fillna(0)
    wide["n_nk"] = wide["n_nk"].fillna(0)
    wide["n_cd8nk"] = wide["n_cd8"] + wide["n_nk"]
    wide["frac_nk"] = np.where(wide["n_cd8nk"] > 0, wide["n_nk"] / wide["n_cd8nk"], np.nan)
    return wide


def pooled_effector(sums: pd.DataFrame, locked: pd.DataFrame) -> pd.DataFrame:
    sub = sums[(sums["slice"] == "ALL") & (sums["compartment"].isin(["cd8", "nk"])) & (sums["gene"].isin(EFFECTOR))]
    rows = []
    for (ds, unit), g in sub.groupby(["dataset", "unit_id"]):
        n_by = g.groupby("compartment")["n_cells"].first()
        n = float(n_by.sum())
        if n <= 0:
            continue
        gene_means = []
        used = []
        for gene, gg in g.groupby("gene"):
            gene_means.append(float(gg["sum_log1p"].sum()) / n)
            used.append(gene)
        rows.append(
            {
                "dataset": ds,
                "unit_id": unit,
                "compartment": "cd8nk",
                "slice": "ALL",
                "score_name": "effector_cosmx",
                "score": float(np.mean(gene_means)) if gene_means else np.nan,
                "n_cells": int(n),
                "n_genes": len(used),
            }
        )
    out = pd.DataFrame(rows)
    out = out.merge(
        locked[["dataset", "unit_id", "mal_CLDN4_pct", "cldn4_quartile"]],
        on=["dataset", "unit_id"],
        how="left",
    )
    return out


def residualize(df: pd.DataFrame, y: str, x: str) -> pd.DataFrame:
    parts = []
    for ds, sub in df.groupby("dataset"):
        sub = sub.dropna(subset=[y, x]).copy()
        if len(sub) < 5 or sub[x].nunique() < 2:
            sub["score"] = np.nan
        else:
            xv = sub[x].to_numpy(dtype=float)
            yv = sub[y].to_numpy(dtype=float)
            X = np.column_stack([np.ones(len(sub)), xv])
            beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
            sub["score"] = yv - X @ beta
        parts.append(sub)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def run_primary_and_sensitivity(scores: pd.DataFrame, sums: pd.DataFrame, locked: pd.DataFrame, comp: pd.DataFrame) -> tuple[list[dict], pd.DataFrame]:
    results = []
    patient_rows = []

    def add(df: pd.DataFrame, label: str, family: str) -> None:
        if df.empty:
            results.append({"label": label, "family": family, "meta": {"k": 0, "N": 0, "rho": np.nan, "p": np.nan}, "note": "empty"})
            return
        rec = summarize_outcome(df, label)
        rec["family"] = family
        results.append(rec)
        tmp = df.copy()
        tmp["label"] = label
        patient_rows.append(tmp)

    for comp_name, slice_name, score_name in PRIMARY:
        df = eligible(scores, comp_name, slice_name, score_name, MIN_N)
        add(df, f"{comp_name}|{slice_name}|{score_name}|min{MIN_N}", "primary")

    # sensitivities on the same machinery
    sens = [
        ("cd8", "ALL", "effector_no_nkg7", MIN_N),
        ("nk", "ALL", "effector_no_nkg7", MIN_N),
        ("cd8", "ALL", "exhaustion_ext", MIN_N),
        ("nk", "ALL", "exhaustion_ext", MIN_N),
        ("cd8", "ALL", "naive_memory", MIN_N),
        ("nk", "ALL", "naive_memory", MIN_N),
        ("cd8", "UMI500", "effector_cosmx", MIN_N),
        ("nk", "UMI500", "effector_cosmx", MIN_N),
        ("cd8", "UMI500", "exhaustion", MIN_N),
        ("nk", "UMI500", "exhaustion", MIN_N),
        ("cd8", "GNLY_POS", "effector_cosmx", MIN_N),
        ("cd8", "ALL", "effector_cosmx", 30),
        ("nk", "ALL", "effector_cosmx", 30),
        ("cd8", "ALL", "exhaustion", 30),
        ("nk", "ALL", "exhaustion", 30),
    ]
    for gene in EFFECTOR + EXHAUSTION:
        sens.append(("cd8", "ALL", f"gene:{gene}", MIN_N))
        sens.append(("nk", "ALL", f"gene:{gene}", MIN_N))
    for comp_name, slice_name, score_name, min_n in sens:
        df = eligible(scores, comp_name, slice_name, score_name, min_n)
        add(df, f"{comp_name}|{slice_name}|{score_name}|min{min_n}", "sensitivity")

    # pooled
    pool = pooled_effector(sums, locked)
    pool = pool[pool["n_cells"] >= MIN_N].rename(columns={"cldn4_quartile": "quartile"})
    add(pool, f"cd8nk|ALL|effector_cosmx|min{MIN_N}", "sensitivity")

    # composition: frac NK
    comp_df = comp[comp["n_cd8nk"] >= MIN_N].rename(columns={"cldn4_quartile": "quartile", "frac_nk": "score"})
    add(comp_df, "composition|frac_nk_among_cd8nk|min10", "composition")

    # naive-adjusted CD8 effector
    cd8_eff = eligible(scores, "cd8", "ALL", "effector_cosmx", MIN_N)[
        ["dataset", "unit_id", "mal_CLDN4_pct", "quartile", "score", "n_cells"]
    ].rename(columns={"score": "effector"})
    cd8_nv = eligible(scores, "cd8", "ALL", "naive_memory", MIN_N)[["dataset", "unit_id", "score"]].rename(
        columns={"score": "naive"}
    )
    both = cd8_eff.merge(cd8_nv, on=["dataset", "unit_id"], how="inner")
    resid = residualize(both, "effector", "naive")
    add(resid, "cd8|ALL|effector_residual_on_naive|min10", "sensitivity")

    # pooled effector residualized on frac_nk
    pool2 = pool.merge(comp[["dataset", "unit_id", "frac_nk"]], on=["dataset", "unit_id"], how="left")
    pool2 = pool2.rename(columns={"score": "effector"})
    # residualize expects y,x and writes score
    resid_p = residualize(pool2, "effector", "frac_nk")
    add(resid_p, "cd8nk|ALL|effector_residual_on_frac_nk|min10", "sensitivity")

    # author slices: within-slice intensity and fraction, per cohort only (no cross-ontology meta)
    author = scores[~scores["slice"].isin(["ALL", "UMI500", "GNLY_POS"])].copy()
    for (ds, comp_name, sl), sub in author.groupby(["dataset", "compartment", "slice"]):
        if sub["score_name"].nunique() == 0:
            continue
        base_n = scores[
            (scores["dataset"] == ds)
            & (scores["compartment"] == comp_name)
            & (scores["slice"] == "ALL")
            & (scores["score_name"] == "effector_cosmx")
        ][["unit_id", "n_cells"]].rename(columns={"n_cells": "n_all"})
        inten = sub[(sub["score_name"] == "effector_cosmx") & (sub["n_cells"] >= MIN_N)].rename(
            columns={"cldn4_quartile": "quartile"}
        )
        if len(inten) >= META_MIN_N:
            add(inten, f"author|{ds}|{comp_name}|{sl}|effector_intensity", "author_subset")
        frac = sub[sub["score_name"] == "effector_cosmx"][["dataset", "unit_id", "n_cells", "mal_CLDN4_pct", "cldn4_quartile"]].drop_duplicates()
        frac = frac.merge(base_n, on="unit_id", how="left")
        frac = frac[frac["n_all"] >= MIN_N].copy()
        frac["score"] = frac["n_cells"] / frac["n_all"]
        frac = frac.rename(columns={"cldn4_quartile": "quartile"})
        if len(frac) >= META_MIN_N:
            # single cohort: still run summarize; meta k=1 is just that cohort
            add(frac, f"author|{ds}|{comp_name}|{sl}|fraction_of_compartment", "author_subset")

    patient = pd.concat(patient_rows, ignore_index=True) if patient_rows else pd.DataFrame()
    return results, patient


def flatten_results(results: list[dict]) -> pd.DataFrame:
    rows = []
    for rec in results:
        meta = rec.get("meta") or {}
        ols = rec.get("ols_q4") or {}
        rb = rec.get("stacked_q4q1") or {}
        rows.append(
            {
                "family": rec.get("family"),
                "label": rec.get("label"),
                "k": meta.get("k"),
                "N": meta.get("N"),
                "rho": meta.get("rho"),
                "p_meta": meta.get("p"),
                "ci_lo": meta.get("ci_lo"),
                "ci_hi": meta.get("ci_hi"),
                "I2": meta.get("I2"),
                "k_rho_neg": rec.get("k_rho_neg"),
                "k_delta_neg": rec.get("k_delta_neg"),
                "k_tested": rec.get("k_tested"),
                "ols_delta": ols.get("delta"),
                "ols_se": ols.get("se"),
                "ols_p": ols.get("p"),
                "ols_n": ols.get("n"),
                "ols_n_q1": ols.get("n_q1"),
                "ols_n_q4": ols.get("n_q4"),
                "stacked_r": rb.get("r"),
                "stacked_p": rb.get("p"),
                "stacked_n_q1": rb.get("n_q1"),
                "stacked_n_q4": rb.get("n_q4"),
            }
        )
    return pd.DataFrame(rows)


def cohort_rows(results: list[dict]) -> pd.DataFrame:
    rows = []
    for rec in results:
        for c in rec.get("cohorts") or []:
            rows.append({"family": rec.get("family"), "label": rec["label"], **c})
    return pd.DataFrame(rows)


def leave_one_out(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for comp_name, slice_name, score_name in PRIMARY:
        df = eligible(scores, comp_name, slice_name, score_name, MIN_N)
        label = f"{comp_name}|{score_name}"
        for drop in COHORT_ORDER:
            sub = df[df["dataset"] != drop]
            cohorts = per_cohort_spearman(sub)
            use = cohorts[(cohorts["n"] >= META_MIN_N) & np.isfinite(cohorts["rho"])]
            meta = dl_spearman(use["rho"], use["n"])
            rows.append({"label": label, "dropped": drop, **{k: meta[k] for k in ("rho", "p", "ci_lo", "ci_hi", "I2", "k", "N")}})
    return pd.DataFrame(rows)


def _forest_panel(ax, rec: dict, title: str) -> None:
    cohorts = pd.DataFrame(rec.get("cohorts") or [])
    meta = rec.get("meta") or {}
    rows = []
    for _, c in cohorts.iterrows():
        if not np.isfinite(c["rho"]) or c["n"] <= 3:
            se = np.nan
            lo = hi = np.nan
        else:
            z = np.arctanh(np.clip(c["rho"], -0.999999, 0.999999))
            se = 1.0 / np.sqrt(c["n"] - 3.0)
            lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
        rows.append((c["dataset"], c["rho"], lo, hi, int(c["n"]), False))
    rows.append(("DL meta", meta.get("rho", np.nan), meta.get("ci_lo", np.nan), meta.get("ci_hi", np.nan), meta.get("N", 0), True))
    rows = rows[::-1]
    ys = np.arange(len(rows))
    for y, (name, rho, lo, hi, n, is_meta) in zip(ys, rows):
        color = "#222222" if is_meta else COHORT_COLOR.get(name, "#555555")
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [y, y], color=color, lw=1.6, solid_capstyle="round")
        if np.isfinite(rho):
            ax.plot(rho, y, "D" if is_meta else "o", color=color, ms=7 if is_meta else 6, zorder=3)
        ax.text(1.02, y, f"n={n}", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=8, color="#333333")
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlim(-1.15, 1.35)
    ax.set_xlabel("Spearman ρ vs malignant CLDN4 %pos", fontsize=8)
    p = meta.get("p", np.nan)
    rho = meta.get("rho", np.nan)
    i2 = meta.get("I2", np.nan)
    subtitle = ""
    if np.isfinite(rho) and np.isfinite(p):
        subtitle = f"\nρ={rho:.2f} ({meta.get('ci_lo', np.nan):.2f}, {meta.get('ci_hi', np.nan):.2f}), p={p:.3g}, I²={100 * i2:.0f}%"
    ax.set_title(title + subtitle, fontsize=9, loc="left")
    ax.tick_params(axis="x", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_figures(results: list[dict], scores: pd.DataFrame, comp: pd.DataFrame, fig_dir: Path) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    by_label = {r["label"]: r for r in results}
    prim = [
        ("cd8|ALL|effector_cosmx|min10", "CD8 effector"),
        ("nk|ALL|effector_cosmx|min10", "NK effector"),
        ("cd8|ALL|exhaustion|min10", "CD8 exhaustion"),
        ("nk|ALL|exhaustion|min10", "NK exhaustion"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.6), constrained_layout=True)
    for ax, (lab, title) in zip(axes.ravel(), prim):
        _forest_panel(ax, by_label.get(lab, {}), title)
    fig.suptitle("Concordant-4 patient CD8/NK state vs malignant CLDN4", fontsize=12)
    fig.savefig(fig_dir / "fig_forest_primary.png", dpi=160, bbox_inches="tight")
    fig.savefig(fig_dir / "fig_forest_primary.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), constrained_layout=True)
    for ax, comp_name, title in zip(axes, ("cd8", "nk"), ("CD8 cells", "NK cells")):
        df = eligible(scores, comp_name, "ALL", "effector_cosmx", MIN_N)
        for ds in COHORT_ORDER:
            sub = df[df["dataset"] == ds]
            if sub.empty:
                continue
            ax.scatter(
                sub["mal_CLDN4_pct"],
                sub["score"],
                s=28,
                c=COHORT_COLOR[ds],
                label=ds,
                alpha=0.9,
                linewidths=0,
            )
        ax.set_xlabel("Malignant CLDN4 % positive")
        ax.set_ylabel("Effector score  mean log1p(CP10k)")
        ax.set_title(title, loc="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02))
    fig.savefig(fig_dir / "fig_scatter_effector.png", dpi=160, bbox_inches="tight")
    fig.savefig(fig_dir / "fig_scatter_effector.pdf", bbox_inches="tight")
    plt.close(fig)

    genes = EFFECTOR + EXHAUSTION
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), constrained_layout=True)
    for ax, comp_name, title in zip(axes, ("cd8", "nk"), ("CD8  Q4−Q1", "NK  Q4−Q1")):
        mat = np.full((len(genes), len(COHORT_ORDER)), np.nan)
        for i, gene in enumerate(genes):
            df = eligible(scores, comp_name, "ALL", f"gene:{gene}", MIN_N)
            for j, ds in enumerate(COHORT_ORDER):
                sub = df[df["dataset"] == ds]
                q4 = sub.loc[sub["quartile"] == "Q4", "score"]
                q1 = sub.loc[sub["quartile"] == "Q1", "score"]
                if len(q4) and len(q1):
                    mat[i, j] = float(np.mean(q4) - np.mean(q1))
        vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
        vmax = max(float(vmax), 0.05)
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(COHORT_ORDER)), [c.replace("GSE", "") for c in COHORT_ORDER], fontsize=8)
        ax.set_yticks(range(len(genes)), genes, fontsize=8)
        ax.set_title(title, loc="left")
        for i in range(len(genes)):
            for j in range(len(COHORT_ORDER)):
                if np.isfinite(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=7, color="#111111")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(fig_dir / "fig_gene_delta.png", dpi=160)
    fig.savefig(fig_dir / "fig_gene_delta.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    # frac NK
    ax = axes[0]
    cdf = comp[comp["n_cd8nk"] >= MIN_N]
    for ds in COHORT_ORDER:
        sub = cdf[cdf["dataset"] == ds]
        ax.scatter(sub["mal_CLDN4_pct"], sub["frac_nk"], s=28, c=COHORT_COLOR[ds], label=ds, linewidths=0)
    ax.set_xlabel("Malignant CLDN4 % positive")
    ax.set_ylabel("NK fraction of CD8+NK")
    ax.set_title("Lineage mix", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    nv = eligible(scores, "cd8", "ALL", "naive_memory", MIN_N)
    for ds in COHORT_ORDER:
        sub = nv[nv["dataset"] == ds]
        ax.scatter(sub["mal_CLDN4_pct"], sub["score"], s=28, c=COHORT_COLOR[ds], linewidths=0)
    ax.set_xlabel("Malignant CLDN4 % positive")
    ax.set_ylabel("CD8 naive/memory score")
    ax.set_title("CD8 naive/memory", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(fig_dir / "fig_composition.png", dpi=160)
    fig.savefig(fig_dir / "fig_composition.pdf")
    plt.close(fig)


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_finding(path: Path, flat: pd.DataFrame, qc: pd.DataFrame, results: list[dict]) -> None:
    by = {r["label"]: r for r in results}

    def line(label: str) -> str:
        r = by[label]
        m = r["meta"]
        o = r.get("ols_q4") or {}
        rb = r.get("stacked_q4q1") or {}
        return (
            f"| {label.split('|')[0]} {label.split('|')[2]} | {m.get('k')} | {m.get('N')} | "
            f"{m.get('rho'):.3f} | {m.get('ci_lo'):.3f} to {m.get('ci_hi'):.3f} | {fmt_p(m.get('p'))} | "
            f"{100 * m.get('I2', float('nan')):.1f}% | {r.get('k_rho_neg')}/{r.get('k_tested')} | "
            f"{o.get('delta', float('nan')):.3f} | {fmt_p(o.get('p'))} | "
            f"{rb.get('n_q1')}/{rb.get('n_q4')} | {rb.get('r', float('nan')):.3f} | {fmt_p(rb.get('p'))} |"
        )

    prim_labels = [f"{c}|ALL|{s}|min10" for c, _, s in ((a, b, c) for a, b, c in PRIMARY)]
    # PRIMARY tuples are (comp, slice, score)
    prim_labels = [f"{c}|ALL|{s}|min10" for c, sl, s in PRIMARY]

    cohort_bits = []
    for lab in prim_labels:
        r = by[lab]
        bits = []
        for c in r["cohorts"]:
            bits.append(f"{c['dataset']} n={c['n']} ρ={c['rho']:.3f} (p={fmt_p(c['p'])}) ΔQ={c['delta_q4_minus_q1']:.3f}")
        cohort_bits.append(f"- {lab}: " + "; ".join(bits))

    qc_lines = []
    for rec in qc.to_dict(orient="records"):
        qc_lines.append(
            f"- {rec['dataset']}: Pearson {rec['pearson']:.4f}, median absolute difference {rec['median_abs_pp']:.3f} percentage points"
        )

    # pull a few sensitivities for the prose table
    sens_labels = [
        "cd8|ALL|naive_memory|min10",
        "cd8|ALL|effector_residual_on_naive|min10",
        "cd8nk|ALL|effector_cosmx|min10",
        "cd8nk|ALL|effector_residual_on_frac_nk|min10",
        "composition|frac_nk_among_cd8nk|min10",
        "cd8|GNLY_POS|effector_cosmx|min10",
        "cd8|UMI500|effector_cosmx|min10",
        "nk|UMI500|effector_cosmx|min10",
        "cd8|ALL|effector_no_nkg7|min10",
        "nk|ALL|effector_no_nkg7|min10",
        "cd8|ALL|effector_cosmx|min30",
        "nk|ALL|effector_cosmx|min30",
    ]

    def sens_line(label: str) -> str:
        r = by.get(label)
        if not r or "meta" not in r or r["meta"].get("k", 0) == 0:
            return f"| {label} |  |  |  |  |"
        m = r["meta"]
        o = r.get("ols_q4") or {}
        return (
            f"| {label} | {m.get('k')} | {m.get('N')} | {m.get('rho'):.3f} | {fmt_p(m.get('p'))} | "
            f"{m.get('ci_lo'):.3f} to {m.get('ci_hi'):.3f} | "
            f"{o.get('delta', float('nan')):.3f} ({fmt_p(o.get('p'))}) |"
        )

    def _m(label: str) -> dict:
        return by[label]["meta"]

    def _o(label: str) -> dict:
        return by[label].get("ols_q4") or {}

    cd8e, nke = _m("cd8|ALL|effector_cosmx|min10"), _m("nk|ALL|effector_cosmx|min10")
    cd8x, nkx = _m("cd8|ALL|exhaustion|min10"), _m("nk|ALL|exhaustion|min10")
    cd8e_o, nke_o = _o("cd8|ALL|effector_cosmx|min10"), _o("nk|ALL|effector_cosmx|min10")
    tox = _m("cd8|ALL|gene:TOX|min10")
    tox_o = _o("cd8|ALL|gene:TOX|min10")
    mt = by.get("author|GSE205335|cd8|CD8+ MT high|fraction_of_compartment", {})
    mt_m = mt.get("meta") or {}

    text = f"""# Concordant-4 CD8/NK state vs malignant CLDN4

Patient-level test of effector intensity and exhaustion inside CD8 cells and inside NK cells, in the locked concordant-4 units (GSE123902, GSE131907, GSE205335, GSE189357).

This is the scRNA composition analog of the CosMx neighbor-state result. CosMx (He 2022) already says CLDN4-high tumor cells have fewer cytotoxic neighbors, and the neighbors that remain are not lower for GZMB, PRF1, NKG7, or IFNG (hi/lo 1.11–1.22, 0/8 down). scRNA has no coordinates. The question here is whether the CD8 and NK cells that are present in CLDN4-high patients have a different effector or exhaustion intensity than those in CLDN4-low patients.

The locked T/NK fraction result is not re-derived. Malignant CLDN4 % positive and the within-cohort quartile are taken from that patient table. A lower T/NK fraction is exclusion. It is not, by itself, a change in the state of the cells that remain.

## Exposure QC

Recomputed malignant CLDN4 % positive versus the locked value (same gates: author malignant in GSE131907 and GSE205335; marker malignant in GSE123902 and GSE189357). The tests use the locked value.

{chr(10).join(qc_lines)}

## Primary result

Unit = patient / donor / sample. Score = mean log1p(CP10k) of the named genes, averaged over cells in the compartment, then averaged over genes. Effector genes are the CosMx set GZMB, PRF1, NKG7, IFNG. Exhaustion genes are PDCD1, HAVCR2, LAG3, TIGIT, TOX. A unit enters a test when that compartment has at least 10 cells. Meta-analysis is DerSimonian–Laird on Fisher z of within-cohort Spearman ρ, cohorts with n≥5. Q4 versus Q1 uses the locked within-cohort quartiles. OLS is score ~ cohort + Q4 among Q1 and Q4 units only (cohort-adjusted difference, Q4 minus Q1). Stacked rank-biserial pools those units across cohorts.

Muzzling, if it were operating at the patient level, would be a negative effector association and a positive exhaustion association. Four tests, two-sided, descriptive.

| compartment | score | k | N | ρ | 95% CI | p | I² | cohorts ρ<0 | OLS Δ | OLS p | Q1/Q4 n | stacked r | stacked p |
|---|---|---:|---:|---:|---|---:|---:|---|---:|---:|---|---:|---:|
{chr(10).join(line(lb) for lb in prim_labels)}

Cohort rows:

{chr(10).join(cohort_bits)}

## Verdict

Patient-level scRNA does not show a drop in effector intensity in CLDN4-high patients, in CD8 cells or in NK cells. It also does not show higher exhaustion. The same meta code, run on the locked T/NK fraction, recovers ρ = −0.531 (p = 1.65×10⁻⁵) and stacked rank-biserial r = −0.724 (n = 19/16). A fraction-sized association would have been visible. These state tests are nulls.

CD8 effector ρ = {cd8e['rho']:.3f} (95% CI {cd8e['ci_lo']:.3f} to {cd8e['ci_hi']:.3f}), p = {fmt_p(cd8e['p'])}, N = {cd8e['N']}, I² = {100 * cd8e['I2']:.0f}%. Cohort-adjusted Q4 minus Q1 = {cd8e_o.get('delta', float('nan')):.3f} (p = {fmt_p(cd8e_o.get('p'))}). The I² is GSE189357 (n = 9, ρ = 0.733, p = 0.025) against GSE205335 (n = 22, ρ = −0.347, p = 0.11). Leave-one-out stays non-significant, including after dropping GSE189357 (ρ = −0.180, p = 0.23). Mean Q4/Q1 ratios for this score are 0.83, 0.95, 0.72, and 1.69. That is not a consistent decrease.

NK effector ρ = {nke['rho']:.3f} (95% CI {nke['ci_lo']:.3f} to {nke['ci_hi']:.3f}), p = {fmt_p(nke['p'])}, N = {nke['N']}, I² = {100 * nke['I2']:.0f}%. The point estimate is slightly higher in CLDN4-high patients, not lower. OLS Δ = {nke_o.get('delta', float('nan')):.3f} (p = {fmt_p(nke_o.get('p'))}).

CD8 exhaustion ρ = {cd8x['rho']:.3f} (95% CI {cd8x['ci_lo']:.3f} to {cd8x['ci_hi']:.3f}), p = {fmt_p(cd8x['p'])}, I² = {100 * cd8x['I2']:.0f}%, and 4/4 cohorts have ρ < 0. The direction is less exhaustion, which is the opposite of a muzzling increase. Dropping GSE205335, the cohort with the largest negative ρ (−0.346, p = 0.11), leaves ρ = −0.096 (p = 0.60). NK exhaustion ρ = {nkx['rho']:.3f}, p = {fmt_p(nkx['p'])}.

Units with fewer than 10 cells in the compartment are out of that test. CD8 keeps 61 of 65. NK keeps 58 of 65. Dropped CD8 units: LX653, LX681, NS_07, NS_12. Dropped NK units: LX653, LX684, EBUS_13, EBUS_28, NS_12, P1090, P1115.

This does not confirm the CosMx neighbor ratios, and it does not contradict them. CosMx scored cells next to a CLDN4-high tumor cell. This scores the CD8 or NK compartment of a CLDN4-high patient. A null here means the patient-level intensity test does not add a muzzling claim on top of the locked exclusion result (fewer T/NK).

Do not quote cell counts as n. Do not add GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, or GSE207422. Do not read a Visium same-spot correlation into this table.

## Sensitivities

Same meta machinery. `effector_residual_on_naive` is the within-cohort residual of the CD8 effector score on the CD8 naive/memory score (IL7R, TCF7, CCR7, SELL, LEF1), so a naive-versus-effector mix is not mistaken for a change inside a fixed program. `GNLY_POS` restricts CD8 to GNLY-positive cells (GNLY is not in the four-gene score) and is the closer analog of intensity inside cells that already look cytotoxic. `frac_nk_among_cd8nk` is lineage composition, not intensity. The Spearman metas all include 0. The closest quartile contrast is the naive-residualized CD8 effector score (OLS p = 0.053); its Spearman meta is ρ = −0.14, p = 0.55.

Gene-level meta-analyses of GZMB, PRF1, NKG7, and IFNG, in CD8 and in NK, all have confidence intervals that include 0. CD8 TOX is the strongest single-gene lean (ρ = {tox['rho']:.3f}, p = {fmt_p(tox['p'])}; Q4−Q1 OLS Δ = {tox_o.get('delta', float('nan')):.3f}, p = {fmt_p(tox_o.get('p'))}). The direction is lower TOX, not higher, and it is one gene in a long sensitivity list. It is not a primary result.

Author cell-state names are not pooled across cohorts. In GSE205335 the fraction of CD8 cells labeled CD8+ MT high tracks malignant CLDN4 (ρ = {mt_m.get('rho', float('nan')):.3f}, p = {fmt_p(mt_m.get('p'))}, n = {mt_m.get('N')}). Effector intensity inside that label does not (ρ = −0.025). That fraction is one label among many author slices, and it is composition inside CD8, not the intensity test. A transitional-NK intensity Spearman in the same cohort (ρ = 0.83, n = 12) is a small subset whose Q4 arm is a single patient; the NK-wide effector score there is ρ = 0.14 (p = 0.56).

| test | k | N | ρ | p | 95% CI | OLS Δ (p) |
|---|---:|---:|---:|---:|---|---|
{chr(10).join(sens_line(lb) for lb in sens_labels)}

## Reproduce

```bash
bash methods/concordant4_cd8nk_state/scripts/download.sh /tmp/concordant4_raw
# GSE205335 RDS is double-gzipped
gzip -dc /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz | gzip -dc > /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds
Rscript methods/concordant4_cd8nk_state/scripts/export_gse205335.R /tmp/concordant4_raw /tmp/concordant4_state/gse205335_cells.tsv.gz
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE123902
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE189357
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE131907
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE205335
python3 methods/concordant4_cd8nk_state/scripts/analyze_state.py
```
"""
    path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", type=Path, default=ROOT / "results" / "tables")
    ap.add_argument("--figures", type=Path, default=ROOT / "results" / "figures")
    ap.add_argument("--calibrate-only", action="store_true")
    args = ap.parse_args()
    locked = pd.read_csv(ROOT / "data" / "locked_patient_units.tsv", sep="\t")
    cal = calibrate(locked)
    if args.calibrate_only:
        return
    sums, qc_raw = load_sums(args.tables)
    qc = qc_exposure(qc_raw, locked)
    scores = module_table(sums, locked)
    comp = composition_table(scores)
    results, patient = run_primary_and_sensitivity(scores, sums, locked, comp)
    flat = flatten_results(results)
    cohorts = cohort_rows(results)
    loo = leave_one_out(scores)
    args.tables.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.tables / "patient_module_scores.tsv", sep="\t", index=False)
    comp.to_csv(args.tables / "cd8_nk_composition.tsv", sep="\t", index=False)
    flat.to_csv(args.tables / "tests.tsv", sep="\t", index=False)
    cohorts.to_csv(args.tables / "tests_by_cohort.tsv", sep="\t", index=False)
    loo.to_csv(args.tables / "leave_one_cohort_out.tsv", sep="\t", index=False)
    qc.to_csv(args.tables / "cldn4_exposure_qc.tsv", sep="\t", index=False)
    # patient-level primary scores, wide enough to audit
    prim = scores[
        (scores["slice"] == "ALL")
        & (scores["score_name"].isin(["effector_cosmx", "exhaustion", "naive_memory", "effector_no_nkg7"]))
    ]
    prim.to_csv(args.tables / "patient_primary_scores.tsv", sep="\t", index=False)
    payload = {
        "calibration_tnk_fraction": cal,
        "qc": qc.to_dict(orient="records"),
        "tests": flat.to_dict(orient="records"),
    }
    (args.tables / "summary.json").write_text(json.dumps(payload, indent=2, default=float))
    draw_figures(results, scores, comp, args.figures)
    write_finding(ROOT / "FINDING.md", flat, qc, results)
    print("wrote", args.tables / "tests.tsv", flush=True)


if __name__ == "__main__":
    main()

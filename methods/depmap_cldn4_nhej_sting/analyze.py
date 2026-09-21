#!/usr/bin/env python3
"""DepMap 24Q4: CLDN4 dependency vs NHEJ / cGAS–STING, and lung RNA coexpression.

Two pre-specified public cuts.

1. Essentiality. Integrated Chronos (CRISPRGeneEffect): CLDN4 gene effect vs
   LIG4, STING1, and CGAS (cGAS). PRKDC is not in that matrix. Its only 24Q4
   gene-effect values are the Humagne-CD Cas12 screens in ScreenGeneEffect
   (n is small and is reported as its own cut). Other classical NHEJ genes
   that are in the integrated matrix (XRCC4/5/6, NHEJ1, DCLRE1C, PAXX) and
   the rest of the cGAS–STING machinery (TBK1, IRF3) are secondary.

2. Expression. Lung cell lines: CLDN4 log2(TPM+1) vs the same NHEJ and
   cGAS–STING genes, including PRKDC, which is present in the RNA matrix.

More negative Chronos = stronger dependency. Spearman is the primary
correlation. BH q is within family × cohort. This is not a genome-wide
discovery screen; the all-gene Spearman files only place the pre-specified
genes on the empirical distribution.
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

NHEJ_GENES = [
    "PRKDC",
    "LIG4",
    "XRCC4",
    "XRCC5",
    "XRCC6",
    "NHEJ1",
    "DCLRE1C",
    "PAXX",
]
STING_GENES = ["CGAS", "STING1", "TBK1", "IRF3"]
NAMED = ["PRKDC", "LIG4", "STING1", "CGAS"]
# Integrated Chronos does not contain PRKDC. Named tests on that matrix:
INTEGRATED_NAMED = ["LIG4", "STING1", "CGAS"]
NHEJ_INTEGRATED = [g for g in NHEJ_GENES if g != "PRKDC"]
STING_OTHER = ["TBK1", "IRF3"]
N_BOOT = 5000
SEED = 0
MIN_N = 8


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * np.nan
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index), []
    block = df[use]
    # Complete cases only, so the mean is not a changing gene set.
    ok = block.notna().all(axis=1)
    z = block.loc[ok].apply(zscore, axis=0)
    out = pd.Series(np.nan, index=df.index, dtype=float)
    out.loc[ok] = z.mean(axis=1).to_numpy()
    return out, use


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt(np.dot(rx, rx) * np.dot(ry, ry))
    if denom == 0:
        return float("nan")
    return float(np.dot(rx, ry) / denom)


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = spearman_rho(x[idx], y[idx])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def pairwise(x: pd.Series, y: pd.Series) -> pd.DataFrame:
    return pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()


def corr_pair(x: pd.Series, y: pd.Series, do_boot: bool = True) -> dict:
    a = pairwise(x, y)
    n = int(len(a))
    if n < MIN_N:
        return {
            "n": n,
            "rho": np.nan,
            "p": np.nan,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
        }
    xv = a["x"].to_numpy(dtype=float)
    yv = a["y"].to_numpy(dtype=float)
    rho, p = stats.spearmanr(xv, yv)
    pr, pp = stats.pearsonr(xv, yv)
    if do_boot:
        lo, hi = bootstrap_spearman_ci(xv, yv)
    else:
        lo, hi = np.nan, np.nan
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci95_low": lo,
        "ci95_high": hi,
        "pearson_r": float(pr),
        "pearson_p": float(pp),
    }


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Pearson correlation of rank residuals after a linear adjustment for z.

    The p-value treats the covariate as fixed. The bootstrap recomputes the
    residualization inside each resample.
    """
    a = pd.concat(
        [x.rename("x"), y.rename("y"), z.rename("z")], axis=1
    ).dropna()
    n = int(len(a))
    if n < MIN_N:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci95_low": np.nan, "ci95_high": np.nan}

    def _partial(xv: np.ndarray, yv: np.ndarray, zv: np.ndarray) -> float:
        rx = stats.rankdata(xv).astype(float)
        ry = stats.rankdata(yv).astype(float)
        rz = stats.rankdata(zv).astype(float)
        design = np.column_stack([np.ones(len(rz)), rz])

        def resid(v: np.ndarray) -> np.ndarray:
            coef, _, _, _ = np.linalg.lstsq(design, v, rcond=None)
            return v - design @ coef

        xr = resid(rx)
        yr = resid(ry)
        if float(np.std(xr)) == 0 or float(np.std(yr)) == 0:
            return float("nan")
        return float(np.corrcoef(xr, yr)[0, 1])

    xv = a["x"].to_numpy(dtype=float)
    yv = a["y"].to_numpy(dtype=float)
    zv = a["z"].to_numpy(dtype=float)
    rho = _partial(xv, yv, zv)
    # Approximate two-sided p from the partial correlation, df = n - 3.
    if not np.isfinite(rho) or abs(rho) >= 1:
        p = np.nan
    else:
        df = n - 3
        t = rho * np.sqrt(df / (1.0 - rho * rho))
        p = float(2 * stats.t.sf(abs(t), df))
    rng = np.random.default_rng(SEED)
    boots = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, n)
        boots[i] = _partial(xv[idx], yv[idx], zv[idx])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return {
        "n": n,
        "rho": float(rho) if np.isfinite(rho) else np.nan,
        "p": p,
        "ci95_low": float(lo),
        "ci95_high": float(hi),
    }


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    q = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return [float("nan")] * n
    pp = p[ok]
    m = len(pp)
    order = np.argsort(pp)
    qq = np.empty(m)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = m - rank + 1
        val = min(prev, pp[i] * m / k)
        qq[i] = val
        prev = val
    q[np.flatnonzero(ok)] = np.clip(qq, 0, 1)
    return [float(v) if np.isfinite(v) else float("nan") for v in q]


def apply_bh(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    out["q_bh"] = np.nan
    if out.empty:
        return out
    for _, idx in out.groupby(group_cols, dropna=False).groups.items():
        idx = list(idx)
        qs = bh_fdr(out.loc[idx, "p"].tolist())
        out.loc[idx, "q_bh"] = qs
    return out


def cohort_tag(row: pd.Series) -> str:
    disease = str(row.get("OncotreePrimaryDisease") or "")
    subtype = str(row.get("OncotreeSubtype") or "")
    if disease == "Non-Small Cell Lung Cancer":
        if subtype == "Lung Adenocarcinoma":
            return "LUAD"
        if subtype == "Lung Squamous Cell Carcinoma":
            return "LUSC"
        return "other_NSCLC"
    if disease == "Lung Neuroendocrine Tumor" or subtype == "Small Cell Lung Cancer":
        return "SCLC_NET"
    return "other_lung"


def chronos_summary(s: pd.Series) -> dict:
    a = s.dropna().astype(float)
    if a.empty:
        return {"n": 0}
    return {
        "n": int(len(a)),
        "median": float(a.median()),
        "q25": float(a.quantile(0.25)),
        "q75": float(a.quantile(0.75)),
        "min": float(a.min()),
        "max": float(a.max()),
        "n_lt_neg0.5": int((a < -0.5).sum()),
        "n_lt_neg1.0": int((a < -1.0).sum()),
        "frac_lt_neg0.5": float((a < -0.5).mean()),
    }


def mannwhitney(high: pd.Series, low: pd.Series) -> dict:
    a = high.dropna()
    b = low.dropna()
    n_high, n_low = int(len(a)), int(len(b))
    empty = {
        "n_high": n_high,
        "n_low": n_low,
        "median_high": np.nan,
        "median_low": np.nan,
        "delta_median": np.nan,
        "U": np.nan,
        "p": np.nan,
        "cliffs_delta": np.nan,
    }
    if n_high < 5 or n_low < 5:
        return empty
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    cliffs = (2.0 * float(u) / (n_high * n_low)) - 1.0
    return {
        "n_high": n_high,
        "n_low": n_low,
        "median_high": float(a.median()),
        "median_low": float(b.median()),
        "delta_median": float(a.median() - b.median()),
        "U": float(u),
        "p": float(p),
        "cliffs_delta": float(cliffs),
    }


def lineage_residual(df: pd.DataFrame, genes: list[str], min_n: int = 8) -> pd.DataFrame:
    """Subtract the OncotreeLineage median. Lineages smaller than min_n become NA."""
    out = pd.DataFrame(index=df.index)
    lineage = df["OncotreeLineage"].astype(str)
    for g in genes:
        if g not in df.columns:
            continue
        vals = df[g]
        counts = vals.groupby(lineage).transform(lambda s: int(s.notna().sum()))
        med = vals.groupby(lineage).transform("median")
        resid = vals - med
        resid = resid.where(counts >= min_n)
        out[g] = resid
    return out


def empirical_place(gw: pd.DataFrame, gene: str) -> dict:
    if gene not in set(gw["gene"]):
        return {"gene": gene, "in_genomewide": False}
    row = gw.loc[gw["gene"] == gene].iloc[0]
    rho = float(row["rho"])
    all_rho = gw["rho"].to_numpy(dtype=float)
    n_g = int(len(all_rho))
    # Includes the gene itself, so the floor is 1/n_genes.
    abs_rho = np.abs(all_rho)
    p_emp = float(np.mean(abs_rho >= abs(rho)))
    abs_rank = int(np.sum(abs_rho > abs(rho)) + 1)
    pct = float(np.mean(all_rho <= rho) * 100.0)
    return {
        "gene": gene,
        "in_genomewide": True,
        "n_pair": int(row["n"]),
        "rho": rho,
        "n_genes_in_null": n_g,
        "abs_rank": abs_rank,
        "empirical_pct_le": pct,
        "empirical_two_sided_p": p_emp,
    }


def rho_quantiles(gw: pd.DataFrame, label: str) -> dict:
    r = gw["rho"].to_numpy(dtype=float)
    qs = np.quantile(r, [0.01, 0.05, 0.5, 0.95, 0.99])
    return {
        "matrix": label,
        "n_genes": int(len(r)),
        "rho_p01": float(qs[0]),
        "rho_p05": float(qs[1]),
        "rho_p50": float(qs[2]),
        "rho_p95": float(qs[3]),
        "rho_p99": float(qs[4]),
        "max_abs_rho": float(np.max(np.abs(r))),
    }


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(rho: float) -> str:
    if rho is None or not np.isfinite(rho):
        return "NA"
    return f"{rho:+.3f}"


def fmt_ci(lo: float, hi: float) -> str:
    if not np.isfinite(lo) or not np.isfinite(hi):
        return "NA"
    return f"[{lo:+.2f}, {hi:+.2f}]"


def add_row(rows: list[dict], **kwargs) -> None:
    rows.append(kwargs)


def quartile_contrast(df: pd.DataFrame, score: str, target: str) -> dict:
    a = df[[score, target]].dropna()
    if len(a) < 20:
        rec = mannwhitney(a[target], a[target].iloc[0:0])
        rec.update({"score": score, "target": target, "n": int(len(a))})
        return rec
    q1 = a[score].quantile(0.25)
    q3 = a[score].quantile(0.75)
    high = a.loc[a[score] >= q3, target]
    low = a.loc[a[score] <= q1, target]
    rec = mannwhitney(high, low)
    rec.update(
        {
            "score": score,
            "target": target,
            "n": int(len(a)),
            "q25_threshold": float(q1),
            "q75_threshold": float(q3),
        }
    )
    return rec


def scatter_panel(ax, df, xcol, ycol, color_col, color_map, xlabel, ylabel, title):
    for key, sub in df.groupby(color_col, dropna=False):
        ax.scatter(
            sub[xcol],
            sub[ycol],
            s=16,
            alpha=0.75,
            c=color_map.get(str(key), "#333333"),
            label=f"{key} n={len(sub)}",
            edgecolors="none",
        )
    ax.axhline(0, color="#dddddd", lw=0.8, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.8, zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_cldn4_nhej_sting")
    ap.add_argument("--outdir", default="methods/depmap_cldn4_nhej_sting")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((data / "download_manifest.json").read_text())
    model = pd.read_csv(data / "Model.csv")
    model["ModelID"] = model["ModelID"].astype(str)
    crispr = pd.read_csv(data / "crispr_gene_effect_slim.csv")
    crispr["ModelID"] = crispr["ModelID"].astype(str)
    expr = pd.read_csv(data / "expression_slim.csv")
    expr["ModelID"] = expr["ModelID"].astype(str)
    screen = pd.read_csv(data / "screen_gene_effect_slim.csv")
    screen["ScreenID"] = screen["ScreenID"].astype(str)
    seq = pd.read_csv(data / "ScreenSequenceMap.csv")
    gw_crispr = pd.read_csv(data / "crispr_cldn4_vs_all_rho.tsv", sep="\t")
    gw_expr = pd.read_csv(data / "expr_lung_cldn4_vs_all_rho.tsv", sep="\t")

    # One library label per screen.
    lib = (
        seq.groupby("ScreenID", as_index=False)["Library"]
        .agg(lambda s: "|".join(sorted(set(map(str, s)))))
    )
    screen = screen.merge(lib, on="ScreenID", how="left")
    screen_model = seq.groupby("ScreenID", as_index=False)["ModelID"].first()
    screen_model["ModelID"] = screen_model["ModelID"].astype(str)
    screen = screen.merge(screen_model, on="ScreenID", how="left")
    screen = screen.merge(
        model[
            [
                "ModelID",
                "CellLineName",
                "OncotreeLineage",
                "OncotreePrimaryDisease",
                "OncotreeSubtype",
                "ModelType",
            ]
        ],
        on="ModelID",
        how="left",
    )

    meta_cols = [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "ModelType",
    ]
    crispr = crispr.merge(model[meta_cols], on="ModelID", how="left")
    expr = expr.merge(model[meta_cols], on="ModelID", how="left")

    integrated_genes = [g for g in ["CLDN4"] + NHEJ_INTEGRATED + STING_GENES if g in crispr.columns]
    missing_integrated = [g for g in ["CLDN4"] + NHEJ_GENES + STING_GENES if g not in crispr.columns]

    # ---------- Chronos summaries ----------
    summary_rows = []
    cohorts_mask = {
        "all_CRISPR": pd.Series(True, index=crispr.index),
        "lung": (crispr["OncotreeLineage"] == "Lung") & (crispr["ModelType"] == "Cell Line"),
        "NSCLC": (crispr["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer")
        & (crispr["ModelType"] == "Cell Line"),
    }
    for cname, mask in cohorts_mask.items():
        sub = crispr.loc[mask]
        for g in integrated_genes:
            rec = chronos_summary(sub[g])
            rec.update({"matrix": "CRISPRGeneEffect", "cohort": cname, "gene": g})
            summary_rows.append(rec)

    cd_mask = screen["ScreenID"].str.contains(r"\.CD\d+$", regex=True) & screen["PRKDC"].notna()
    cd = screen.loc[cd_mask].copy()
    cd_lung = cd[(cd["OncotreeLineage"] == "Lung") & (cd["ModelType"] == "Cell Line")]
    cd_nsclc = cd[
        (cd["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer")
        & (cd["ModelType"] == "Cell Line")
    ]
    for cname, sub in {
        "HumagneCD": cd,
        "HumagneCD_lung": cd_lung,
        "HumagneCD_NSCLC": cd_nsclc,
    }.items():
        for g in [c for c in ["CLDN4", "PRKDC", "LIG4", "STING1", "CGAS"] if c in sub.columns]:
            rec = chronos_summary(sub[g])
            rec.update({"matrix": "ScreenGeneEffect_HumagneCD", "cohort": cname, "gene": g})
            summary_rows.append(rec)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tabdir / "chronos_summary.tsv", sep="\t", index=False)

    # ---------- Dependency correlations ----------
    dep_rows: list[dict] = []

    def add_gene_corrs(frame: pd.DataFrame, cohort: str, matrix: str, genes: list[str], family: str, xcol: str = "CLDN4"):
        for g in genes:
            if g not in frame.columns or xcol not in frame.columns:
                add_row(
                    dep_rows,
                    matrix=matrix,
                    cohort=cohort,
                    family=family,
                    gene=g,
                    n=0,
                    rho=np.nan,
                    p=np.nan,
                    ci95_low=np.nan,
                    ci95_high=np.nan,
                    pearson_r=np.nan,
                    pearson_p=np.nan,
                )
                continue
            rec = corr_pair(frame[xcol], frame[g])
            rec.update(matrix=matrix, cohort=cohort, family=family, gene=g)
            dep_rows.append(rec)

    for cname, mask in cohorts_mask.items():
        sub = crispr.loc[mask]
        add_gene_corrs(sub, cname, "CRISPRGeneEffect", INTEGRATED_NAMED, "named")
        add_gene_corrs(sub, cname, "CRISPRGeneEffect", NHEJ_INTEGRATED, "NHEJ")
        add_gene_corrs(sub, cname, "CRISPRGeneEffect", STING_GENES, "STING")
        # Signatures of integrated genes, within this cohort.
        sig_nhej, used_nhej = signature(sub, NHEJ_INTEGRATED)
        sig_sting, used_sting = signature(sub, STING_GENES)
        for label, series, used in (
            ("sig_NHEJ", sig_nhej, used_nhej),
            ("sig_STING", sig_sting, used_sting),
        ):
            rec = corr_pair(sub["CLDN4"], series)
            rec.update(
                matrix="CRISPRGeneEffect",
                cohort=cname,
                family="signature",
                gene=label,
                n_genes_in_signature=len(used),
                signature_genes=",".join(used),
            )
            dep_rows.append(rec)

    # Lineage-residualized integrated Chronos (pan-cancer models with a lineage).
    resid_genes = [g for g in ["CLDN4"] + NHEJ_INTEGRATED + STING_GENES if g in crispr.columns]
    resid = lineage_residual(crispr, resid_genes, min_n=8)
    resid["CLDN4"] = resid["CLDN4"]
    add_gene_corrs(resid, "all_CRISPR_lineage_resid", "CRISPRGeneEffect", INTEGRATED_NAMED, "named")
    add_gene_corrs(resid, "all_CRISPR_lineage_resid", "CRISPRGeneEffect", NHEJ_INTEGRATED, "NHEJ")
    add_gene_corrs(resid, "all_CRISPR_lineage_resid", "CRISPRGeneEffect", STING_GENES, "STING")
    # Signature of residuals, z-scored inside the residualized rows.
    resid_frame = resid.copy()
    sig_nhej, used_nhej = signature(resid_frame, NHEJ_INTEGRATED)
    sig_sting, used_sting = signature(resid_frame, STING_GENES)
    for label, series, used in (
        ("sig_NHEJ", sig_nhej, used_nhej),
        ("sig_STING", sig_sting, used_sting),
    ):
        rec = corr_pair(resid_frame["CLDN4"], series)
        rec.update(
            matrix="CRISPRGeneEffect",
            cohort="all_CRISPR_lineage_resid",
            family="signature",
            gene=label,
            n_genes_in_signature=len(used),
            signature_genes=",".join(used),
        )
        dep_rows.append(rec)

    # Humagne-CD PRKDC cut. Same screen-level matrix for both axes.
    for cname, sub in {
        "HumagneCD": cd,
        "HumagneCD_lung": cd_lung,
        "HumagneCD_NSCLC": cd_nsclc,
    }.items():
        add_gene_corrs(
            sub,
            cname,
            "ScreenGeneEffect_HumagneCD",
            ["PRKDC", "LIG4", "STING1", "CGAS"],
            "named_HumagneCD",
        )

    dep = pd.DataFrame(dep_rows)
    dep = apply_bh(dep, ["matrix", "cohort", "family"])
    dep.to_csv(tabdir / "dependency_correlations.tsv", sep="\t", index=False)

    # ---------- Expression ----------
    lung_expr = expr[(expr["OncotreeLineage"] == "Lung") & (expr["ModelType"] == "Cell Line")].copy()
    lung_expr = lung_expr.dropna(subset=["CLDN4"]).copy()
    lung_expr["group"] = lung_expr.apply(cohort_tag, axis=1)
    nsclc_expr = lung_expr[lung_expr["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"].copy()
    luad_expr = lung_expr[lung_expr["group"] == "LUAD"].copy()

    expr_rows: list[dict] = []

    def add_expr(frame: pd.DataFrame, cohort: str) -> None:
        families = {
            "NHEJ": NHEJ_GENES,
            "STING": STING_GENES,
        }
        for family, genes in families.items():
            for g in genes:
                if g not in frame.columns:
                    continue
                rec = corr_pair(frame["CLDN4"], frame[g])
                rec.update(cohort=cohort, family=family, gene=g)
                expr_rows.append(rec)
        sig_n, used_n = signature(frame, [g for g in NHEJ_GENES if g in frame.columns])
        sig_s, used_s = signature(frame, [g for g in STING_GENES if g in frame.columns])
        frame = frame.copy()
        frame["sig_NHEJ"] = sig_n
        frame["sig_STING"] = sig_s
        for label, used in (("sig_NHEJ", used_n), ("sig_STING", used_s)):
            rec = corr_pair(frame["CLDN4"], frame[label])
            rec.update(
                cohort=cohort,
                family="signature",
                gene=label,
                n_genes_in_signature=len(used),
                signature_genes=",".join(used),
            )
            expr_rows.append(rec)
        # Store signatures back via attribute for quartile / partial / plots.
        frame.attrs["used_nhej"] = used_n
        frame.attrs["used_sting"] = used_s
        return frame

    lung_scored = add_expr(lung_expr, "lung")
    nsclc_scored = add_expr(nsclc_expr, "NSCLC")
    luad_scored = add_expr(luad_expr, "LUAD")
    expr_corr = pd.DataFrame(expr_rows)
    expr_corr = apply_bh(expr_corr, ["cohort", "family"])
    expr_corr.to_csv(tabdir / "expression_correlations.tsv", sep="\t", index=False)

    # Partials: proliferation (MKI67) and an epithelial marker (EPCAM).
    partial_rows = []
    for cohort, frame in (("lung", lung_scored), ("NSCLC", nsclc_scored)):
        for target in ["sig_NHEJ", "sig_STING", "PRKDC", "LIG4", "STING1", "CGAS"]:
            if target not in frame.columns:
                continue
            for cov in ["MKI67", "EPCAM"]:
                if cov not in frame.columns:
                    continue
                rec = partial_spearman(frame["CLDN4"], frame[target], frame[cov])
                rec.update(cohort=cohort, target=target, covariate=cov)
                partial_rows.append(rec)
    partial = pd.DataFrame(partial_rows)
    partial = apply_bh(partial, ["cohort", "covariate"])
    partial.to_csv(tabdir / "expression_partial.tsv", sep="\t", index=False)

    # CLDN4-high vs low (top vs bottom quartile) on lung and NSCLC RNA.
    q_rows = []
    for cohort, frame in (("lung", lung_scored), ("NSCLC", nsclc_scored)):
        targets = ["sig_NHEJ", "sig_STING"] + [g for g in NAMED if g in frame.columns]
        for target in targets:
            family = (
                "signature"
                if target.startswith("sig_")
                else ("NHEJ" if target in NHEJ_GENES else "STING")
            )
            rec = quartile_contrast(frame, "CLDN4", target)
            rec.update(cohort=cohort, family=family)
            q_rows.append(rec)
    qtab = pd.DataFrame(q_rows)
    qtab = apply_bh(qtab, ["cohort", "family"])
    qtab.to_csv(tabdir / "expression_q4q1.tsv", sep="\t", index=False)

    # ---------- Empirical placement ----------
    place_rows = []
    for gene in INTEGRATED_NAMED + NHEJ_INTEGRATED + STING_GENES:
        rec = empirical_place(gw_crispr, gene)
        rec.update(matrix="CRISPRGeneEffect_all_models", query="CLDN4_Chronos")
        place_rows.append(rec)
    for gene in NHEJ_GENES + STING_GENES:
        rec = empirical_place(gw_expr, gene)
        rec.update(matrix="expression_lung_cell_lines", query="CLDN4_RNA")
        place_rows.append(rec)
    # Also place the signatures approximately by not using genome-wide for them.
    place = pd.DataFrame(place_rows)
    place.to_csv(tabdir / "empirical_placement.tsv", sep="\t", index=False)
    quant = pd.DataFrame(
        [
            rho_quantiles(gw_crispr, "CRISPRGeneEffect_CLDN4_vs_all"),
            rho_quantiles(gw_expr, "expression_lung_CLDN4_vs_all"),
        ]
    )
    quant.to_csv(tabdir / "background_rho_quantiles.tsv", sep="\t", index=False)

    # Context only: genes already used as epithelial / IFN anchors elsewhere.
    # Not a new discovery list. Ranks are |ρ| among the lung RNA genome-wide file.
    context_genes = ["KRT19", "CLDN7", "EPCAM", "TACSTD2", "CD274", "KRT8"]
    ctx_rows = []
    for gene in context_genes:
        rec = empirical_place(gw_expr, gene)
        rec["matrix"] = "expression_lung_cell_lines"
        rec["query"] = "CLDN4_RNA"
        rec["role"] = "context_not_a_new_test"
        ctx_rows.append(rec)
    pd.DataFrame(ctx_rows).to_csv(tabdir / "expression_context_ranks.tsv", sep="\t", index=False)

    # ---------- Line-level tables ----------
    keep_c = [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "ModelType",
    ] + integrated_genes
    crispr[keep_c].to_csv(tabdir / "crispr_model_scores.tsv", sep="\t", index=False)
    expr_keep = [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "group",
        "CLDN4",
        "sig_NHEJ",
        "sig_STING",
    ] + [g for g in NHEJ_GENES + STING_GENES + ["MKI67", "EPCAM"] if g in lung_scored.columns]
    lung_scored[expr_keep].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "lung_expression_lines.tsv", sep="\t", index=False
    )
    cd_cols = [
        c
        for c in [
            "ScreenID",
            "ModelID",
            "Library",
            "CellLineName",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "ModelType",
            "CLDN4",
            "PRKDC",
            "LIG4",
            "STING1",
            "CGAS",
        ]
        if c in cd.columns
    ]
    cd[cd_cols].sort_values("CellLineName").to_csv(
        tabdir / "humagne_cd_screen_scores.tsv", sep="\t", index=False
    )

    # Cohort counts.
    counts = []
    n_lung_models = int(
        ((model["OncotreeLineage"] == "Lung") & (model["ModelType"] == "Cell Line")).sum()
    )
    counts.append({"cohort": "lung_models_Model_csv", "n": n_lung_models})
    counts.append({"cohort": "integrated_CRISPR_models", "n": int(len(crispr))})
    counts.append({"cohort": "integrated_CRISPR_lung", "n": int(cohorts_mask["lung"].sum())})
    counts.append({"cohort": "integrated_CRISPR_NSCLC", "n": int(cohorts_mask["NSCLC"].sum())})
    counts.append({"cohort": "HumagneCD_screens_with_PRKDC", "n": int(len(cd))})
    counts.append({"cohort": "HumagneCD_lung", "n": int(len(cd_lung))})
    counts.append({"cohort": "HumagneCD_NSCLC", "n": int(len(cd_nsclc))})
    counts.append({"cohort": "lung_RNA_with_CLDN4", "n": int(len(lung_scored))})
    counts.append({"cohort": "NSCLC_RNA_with_CLDN4", "n": int(len(nsclc_scored))})
    counts.append({"cohort": "LUAD_RNA_with_CLDN4", "n": int(len(luad_scored))})
    for g, n in lung_scored["group"].value_counts().items():
        counts.append({"cohort": f"lung_RNA_{g}", "n": int(n)})
    pd.DataFrame(counts).to_csv(tabdir / "cohort_counts.tsv", sep="\t", index=False)

    # ---------- Figures ----------
    # 1. Integrated Chronos scatters for the three named genes that exist.
    plot_df = crispr.copy()
    plot_df["panel"] = np.where(cohorts_mask["lung"], "lung", "other")
    color_map = {"lung": "#c2410c", "other": "#94a3b8"}
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8), constrained_layout=True)
    for ax, gene in zip(axes, INTEGRATED_NAMED):
        rec = dep[
            (dep["matrix"] == "CRISPRGeneEffect")
            & (dep["cohort"] == "all_CRISPR")
            & (dep["family"] == "named")
            & (dep["gene"] == gene)
        ].iloc[0]
        scatter_panel(
            ax,
            plot_df.dropna(subset=["CLDN4", gene]),
            "CLDN4",
            gene,
            "panel",
            color_map,
            "CLDN4 Chronos",
            f"{gene} Chronos",
            f"all CRISPR n={int(rec['n'])}\n"
            f"ρ={fmt_rho(rec['rho'])}  p={fmt_p(rec['p'])}  q={fmt_p(rec['q_bh'])}",
        )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle(
        "DepMap 24Q4 integrated Chronos — CLDN4 vs LIG4 / STING1 / CGAS",
        fontsize=11,
    )
    fig.savefig(figdir / "fig_cldn4_chronos_vs_lig4_sting_cgas.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_cldn4_chronos_vs_lig4_sting_cgas.pdf", bbox_inches="tight")
    plt.close(fig)

    # 2. Humagne-CD PRKDC scatter.
    fig, ax = plt.subplots(figsize=(4.6, 4.2), constrained_layout=True)
    if len(cd):
        cd_plot = cd.copy()
        cd_plot["panel"] = np.where(
            (cd_plot["OncotreeLineage"] == "Lung") & (cd_plot["ModelType"] == "Cell Line"),
            "lung",
            "other",
        )
        recs = dep[
            (dep["matrix"] == "ScreenGeneEffect_HumagneCD")
            & (dep["cohort"] == "HumagneCD")
            & (dep["gene"] == "PRKDC")
        ]
        rec = recs.iloc[0] if len(recs) else None
        scatter_panel(
            ax,
            cd_plot.dropna(subset=["CLDN4", "PRKDC"]),
            "CLDN4",
            "PRKDC",
            "panel",
            color_map,
            "CLDN4 screen Chronos",
            "PRKDC screen Chronos",
            (
                f"Humagne-CD Cas12 n={int(rec['n'])}\nρ={fmt_rho(rec['rho'])}  p={fmt_p(rec['p'])}"
                if rec is not None
                else "Humagne-CD"
            ),
        )
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("PRKDC is only in the Humagne-CD screens", fontsize=11)
    fig.savefig(figdir / "fig_humagne_cd_cldn4_vs_prkdc.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_humagne_cd_cldn4_vs_prkdc.pdf", bbox_inches="tight")
    plt.close(fig)

    # 3. Lung expression scatters.
    hist_colors = {
        "LUAD": "#1f77b4",
        "LUSC": "#ff7f0e",
        "other_NSCLC": "#2ca02c",
        "SCLC_NET": "#d62728",
        "other_lung": "#7f7f7f",
    }
    expr_panels = [
        ("sig_NHEJ", "NHEJ mean-z"),
        ("sig_STING", "cGAS–STING mean-z"),
        ("PRKDC", "PRKDC log2(TPM+1)"),
        ("LIG4", "LIG4 log2(TPM+1)"),
        ("STING1", "STING1 log2(TPM+1)"),
        ("CGAS", "CGAS log2(TPM+1)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.6, 7.2), constrained_layout=True)
    for ax, (col, ylab) in zip(axes.ravel(), expr_panels):
        fam = "signature" if col.startswith("sig_") else ("NHEJ" if col in NHEJ_GENES else "STING")
        hit = expr_corr[
            (expr_corr["cohort"] == "lung")
            & (expr_corr["family"] == fam)
            & (expr_corr["gene"] == col)
        ]
        rec = hit.iloc[0]
        for g, sub in lung_scored.groupby("group"):
            ax.scatter(
                sub["CLDN4"],
                sub[col],
                s=16,
                alpha=0.8,
                c=hist_colors.get(g, "#333333"),
                label=f"{g} n={len(sub)}",
                edgecolors="none",
            )
        ax.set_xlabel("CLDN4 log2(TPM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(
            f"lung n={int(rec['n'])}  ρ={fmt_rho(rec['rho'])}\np={fmt_p(rec['p'])}  q={fmt_p(rec['q_bh'])}",
            fontsize=9,
        )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("DepMap 24Q4 lung lines — CLDN4 RNA vs NHEJ / cGAS–STING RNA", fontsize=11)
    fig.savefig(figdir / "fig_lung_cldn4_expr_vs_nhej_sting.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_lung_cldn4_expr_vs_nhej_sting.pdf", bbox_inches="tight")
    plt.close(fig)

    # 4. Forest of primary Spearman rhos.
    forest_specs = [
        ("CRISPRGeneEffect", "all_CRISPR", "named", "LIG4", "Chronos all, LIG4"),
        ("CRISPRGeneEffect", "all_CRISPR", "named", "STING1", "Chronos all, STING1"),
        ("CRISPRGeneEffect", "all_CRISPR", "named", "CGAS", "Chronos all, CGAS"),
        ("CRISPRGeneEffect", "lung", "named", "LIG4", "Chronos lung, LIG4"),
        ("CRISPRGeneEffect", "lung", "named", "STING1", "Chronos lung, STING1"),
        ("CRISPRGeneEffect", "lung", "named", "CGAS", "Chronos lung, CGAS"),
        ("CRISPRGeneEffect", "all_CRISPR_lineage_resid", "named", "LIG4", "Chronos lineage-resid, LIG4"),
        ("CRISPRGeneEffect", "all_CRISPR_lineage_resid", "named", "STING1", "Chronos lineage-resid, STING1"),
        ("CRISPRGeneEffect", "all_CRISPR_lineage_resid", "named", "CGAS", "Chronos lineage-resid, CGAS"),
        ("ScreenGeneEffect_HumagneCD", "HumagneCD", "named_HumagneCD", "PRKDC", "Humagne-CD, PRKDC"),
    ]
    # Expression primary genes on lung.
    forest_expr = [
        ("lung", "NHEJ", "PRKDC", "Lung RNA, PRKDC"),
        ("lung", "NHEJ", "LIG4", "Lung RNA, LIG4"),
        ("lung", "STING", "STING1", "Lung RNA, STING1"),
        ("lung", "STING", "CGAS", "Lung RNA, CGAS"),
        ("lung", "signature", "sig_NHEJ", "Lung RNA, NHEJ score"),
        ("lung", "signature", "sig_STING", "Lung RNA, STING score"),
    ]
    labels_f = []
    rhos_f = []
    lo_f = []
    hi_f = []
    for matrix, cohort, family, gene, label in forest_specs:
        rec = dep[
            (dep["matrix"] == matrix)
            & (dep["cohort"] == cohort)
            & (dep["family"] == family)
            & (dep["gene"] == gene)
        ].iloc[0]
        labels_f.append(f"{label} (n={int(rec['n'])})")
        rhos_f.append(rec["rho"])
        lo_f.append(rec["ci95_low"])
        hi_f.append(rec["ci95_high"])
    for cohort, family, gene, label in forest_expr:
        rec = expr_corr[
            (expr_corr["cohort"] == cohort)
            & (expr_corr["family"] == family)
            & (expr_corr["gene"] == gene)
        ].iloc[0]
        labels_f.append(f"{label} (n={int(rec['n'])})")
        rhos_f.append(rec["rho"])
        lo_f.append(rec["ci95_low"])
        hi_f.append(rec["ci95_high"])
    y = np.arange(len(labels_f))[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 6.4), constrained_layout=True)
    ax.axvline(0, color="#cbd5e1", lw=1)
    for yi, rho, lo, hi in zip(y, rhos_f, lo_f, hi_f):
        ax.plot([lo, hi], [yi, yi], color="#334155", lw=1.4)
        ax.plot(rho, yi, "o", color="#c2410c", ms=5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels_f, fontsize=8)
    ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
    ax.set_title("Pre-specified CLDN4 correlations")
    fig.savefig(figdir / "fig_rho_forest.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_rho_forest.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---------- Key JSON ----------
    def dep_rec(matrix, cohort, family, gene) -> dict:
        rec = dep[
            (dep["matrix"] == matrix)
            & (dep["cohort"] == cohort)
            & (dep["family"] == family)
            & (dep["gene"] == gene)
        ].iloc[0]
        return {
            "n": int(rec["n"]),
            "rho": None if not np.isfinite(rec["rho"]) else float(rec["rho"]),
            "p": None if not np.isfinite(rec["p"]) else float(rec["p"]),
            "q_bh": None if not np.isfinite(rec["q_bh"]) else float(rec["q_bh"]),
            "ci95": [None if not np.isfinite(rec["ci95_low"]) else float(rec["ci95_low"]),
                     None if not np.isfinite(rec["ci95_high"]) else float(rec["ci95_high"])],
            "pearson_r": None if not np.isfinite(rec["pearson_r"]) else float(rec["pearson_r"]),
        }

    def expr_rec(cohort, family, gene) -> dict:
        rec = expr_corr[
            (expr_corr["cohort"] == cohort)
            & (expr_corr["family"] == family)
            & (expr_corr["gene"] == gene)
        ].iloc[0]
        return {
            "n": int(rec["n"]),
            "rho": float(rec["rho"]),
            "p": float(rec["p"]),
            "q_bh": float(rec["q_bh"]),
            "ci95": [float(rec["ci95_low"]), float(rec["ci95_high"])],
            "pearson_r": float(rec["pearson_r"]),
        }

    key = {
        "release": manifest.get("release"),
        "doi": manifest.get("doi"),
        "missing_from_integrated_chronos": missing_integrated,
        "cd_libraries": manifest.get("cd_libraries"),
        "cohort_counts": counts,
        "chronos_cldn4_all": chronos_summary(crispr["CLDN4"]),
        "chronos_cldn4_lung": chronos_summary(crispr.loc[cohorts_mask["lung"], "CLDN4"]),
        "named_dependency": {
            "LIG4_all": dep_rec("CRISPRGeneEffect", "all_CRISPR", "named", "LIG4"),
            "STING1_all": dep_rec("CRISPRGeneEffect", "all_CRISPR", "named", "STING1"),
            "CGAS_all": dep_rec("CRISPRGeneEffect", "all_CRISPR", "named", "CGAS"),
            "LIG4_lung": dep_rec("CRISPRGeneEffect", "lung", "named", "LIG4"),
            "STING1_lung": dep_rec("CRISPRGeneEffect", "lung", "named", "STING1"),
            "CGAS_lung": dep_rec("CRISPRGeneEffect", "lung", "named", "CGAS"),
            "LIG4_resid": dep_rec("CRISPRGeneEffect", "all_CRISPR_lineage_resid", "named", "LIG4"),
            "STING1_resid": dep_rec("CRISPRGeneEffect", "all_CRISPR_lineage_resid", "named", "STING1"),
            "CGAS_resid": dep_rec("CRISPRGeneEffect", "all_CRISPR_lineage_resid", "named", "CGAS"),
            "PRKDC_HumagneCD": dep_rec("ScreenGeneEffect_HumagneCD", "HumagneCD", "named_HumagneCD", "PRKDC"),
            "PRKDC_HumagneCD_lung": dep_rec(
                "ScreenGeneEffect_HumagneCD", "HumagneCD_lung", "named_HumagneCD", "PRKDC"
            ),
        },
        "signatures_dependency": {
            "NHEJ_all": dep_rec("CRISPRGeneEffect", "all_CRISPR", "signature", "sig_NHEJ"),
            "STING_all": dep_rec("CRISPRGeneEffect", "all_CRISPR", "signature", "sig_STING"),
            "NHEJ_lung": dep_rec("CRISPRGeneEffect", "lung", "signature", "sig_NHEJ"),
            "STING_lung": dep_rec("CRISPRGeneEffect", "lung", "signature", "sig_STING"),
        },
        "expression_lung": {
            g: expr_rec("lung", "NHEJ" if g in NHEJ_GENES else "STING", g) for g in NAMED
        },
        "expression_lung_signatures": {
            "NHEJ": expr_rec("lung", "signature", "sig_NHEJ"),
            "STING": expr_rec("lung", "signature", "sig_STING"),
        },
        "expression_nsclc_signatures": {
            "NHEJ": expr_rec("NSCLC", "signature", "sig_NHEJ"),
            "STING": expr_rec("NSCLC", "signature", "sig_STING"),
        },
        "empirical": place.to_dict(orient="records"),
        "notes": [
            "Spearman primary; Pearson stored alongside.",
            "BH q within matrix × cohort × family.",
            "Bootstrap 5000 resamples, seed 0.",
            "Lineage residual = Chronos minus OncotreeLineage median (lineages with n>=8).",
            "Humagne-CD PRKDC is screen-level Chronos, not the integrated model-level matrix.",
            "Expression signatures are the mean of within-cohort gene-wise z-scores, complete cases.",
            "No immune infiltrate and no IFN treatment in these cultured lines.",
        ],
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")

    # Markdown tables for the write-up.
    def md_dep(mask: pd.Series) -> str:
        lines = [
            "| matrix | cohort | family | gene | n | ρ | p | q | 95% CI | Pearson r |",
            "|---|---|---|---|---:|---:|---:|---:|---|---:|",
        ]
        sub = dep.loc[mask]
        for _, r in sub.iterrows():
            lines.append(
                f"| {r['matrix']} | {r['cohort']} | {r['family']} | {r['gene']} | {int(r['n'])} | "
                f"{fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {fmt_p(r['q_bh'])} | "
                f"{fmt_ci(r['ci95_low'], r['ci95_high'])} | {fmt_rho(r['pearson_r'])} |"
            )
        return "\n".join(lines) + "\n"

    focus = dep["gene"].isin(NAMED + ["sig_NHEJ", "sig_STING"]) & dep["cohort"].isin(
        ["all_CRISPR", "lung", "NSCLC", "all_CRISPR_lineage_resid", "HumagneCD", "HumagneCD_lung", "HumagneCD_NSCLC"]
    )
    (tabdir / "dependency_focus.md").write_text(md_dep(focus))

    def md_expr() -> str:
        lines = [
            "| cohort | family | gene | n | ρ | p | q | 95% CI | Pearson r |",
            "|---|---|---|---:|---:|---:|---:|---|---:|",
        ]
        sub = expr_corr[expr_corr["cohort"].isin(["lung", "NSCLC"])]
        for _, r in sub.iterrows():
            lines.append(
                f"| {r['cohort']} | {r['family']} | {r['gene']} | {int(r['n'])} | "
                f"{fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {fmt_p(r['q_bh'])} | "
                f"{fmt_ci(r['ci95_low'], r['ci95_high'])} | {fmt_rho(r['pearson_r'])} |"
            )
        return "\n".join(lines) + "\n"

    (tabdir / "expression_focus.md").write_text(md_expr())

    print(json.dumps({
        "missing_integrated": missing_integrated,
        "n_cd": int(len(cd)),
        "n_cd_lung": int(len(cd_lung)),
        "n_crispr": int(len(crispr)),
        "n_lung_rna": int(len(lung_scored)),
        "named": key["named_dependency"],
        "expr_sig": key["expression_lung_signatures"],
        "expr_named": key["expression_lung"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""CPTAC LSCC protein TACSTD2/CLDN4 vs xCell CD8/immune, purity residual.

Recomputes the xCell slice that PR23 reported marginally, then residualizes
on DNA tumor purity (WES primary; WGS and ESTIMATE-derived as sensitivity).

Treatment-naive surgical LSCC (Satpathy et al. Cell 2021). No ICI labels.
"""

from __future__ import annotations

import json
import math
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "w200" / "CPTAC_LSCC_xcell"
OUT = ROOT / "results" / "w200" / "CPTAC_LSCC_xcell"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

PROTEIN = DATA / "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
PHENO = DATA / "LSCC_phenotype.txt"
META = DATA / "LSCC_meta.txt"

TARGETS = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
}

# Prespecified primary xCell endpoints (user request: CD8 / immune).
PRIMARY_XCELL = ["xCell_T_cell_CD8+", "xCell_immune_score"]
# CD8 subsets: exploratory only.
EXPLORATORY_XCELL = [
    "xCell_T_cell_CD8+_naive",
    "xCell_T_cell_CD8+_central_memory",
    "xCell_T_cell_CD8+_effector_memory",
]
# DNA purity is the honest residual (not RNA-circular with xCell).
PRIMARY_PURITY = "WES_purity"
SENSITIVITY_PURITY = ["WGS_purity", "ESTIMATE_TumorPurity"]

# Yoshihara et al. 2013 ESTIMATE cosine transform.
ESTIMATE_PURITY_A = 0.6049872018
ESTIMATE_PURITY_B = 0.0001467884


def match_ensembl_row(index, ensembl_id: str) -> str | None:
    prefix = ensembl_id.split(".")[0]
    for raw in index:
        if str(raw).split(".")[0] == prefix:
            return str(raw)
    return None


def bh_fdr(pvalues) -> list[float]:
    p = np.asarray(list(pvalues), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty(m, dtype=float)
    restored[order] = q
    out[ok] = restored
    return out.tolist()


def spearman_pair(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(d.shape[0])
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    a = d.iloc[:, 0].to_numpy(float)
    b = d.iloc[:, 1].to_numpy(float)
    if np.nanstd(a) == 0 or np.nanstd(b) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(a, b)
    return {"n": n, "rho": float(rho), "p": float(p)}


def residualize(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(covariate)), covariate])
    coef, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coef


def complete_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for a in arrays:
        mask &= np.isfinite(a)
    return mask


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """One-covariate partial Spearman: Pearson of rank-residuals on ranked z."""
    mask = complete_mask(x, y, z)
    x, y, z = x[mask], y[mask], z[mask]
    n = int(x.size)
    if n < 5:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    if np.nanstd(x) == 0 or np.nanstd(y) == 0 or np.nanstd(z) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    xr = residualize(stats.rankdata(x), stats.rankdata(z))
    yr = residualize(stats.rankdata(y), stats.rankdata(z))
    if np.nanstd(xr) == 0 or np.nanstd(yr) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    rho = float(stats.pearsonr(xr, yr).statistic)
    df = n - 3
    if abs(rho) >= 1:
        p = 0.0
    else:
        t = rho * math.sqrt(df / (1.0 - rho * rho))
        p = float(2 * stats.t.sf(abs(t), df))
    zf = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    # Fisher-z for one controlled covariate uses n-4.
    half = stats.norm.ppf(0.975) / math.sqrt(max(n - 4, 1))
    lo, hi = np.tanh([zf - half, zf + half])
    return {"n": n, "rho": rho, "p": p, "ci_low": float(lo), "ci_high": float(hi)}


def ols_residual_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """Spearman of linear residuals of raw x and y on raw z (literal purity residual)."""
    mask = complete_mask(x, y, z)
    x, y, z = x[mask], y[mask], z[mask]
    n = int(x.size)
    if n < 5 or np.nanstd(z) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    xr = residualize(x, z)
    yr = residualize(y, z)
    if np.nanstd(xr) == 0 or np.nanstd(yr) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(xr, yr)
    return {"n": n, "rho": float(rho), "p": float(p)}


def estimate_tumor_purity(estimate_score: pd.Series) -> pd.Series:
    s = pd.to_numeric(estimate_score, errors="coerce")
    return np.cos(ESTIMATE_PURITY_A + ESTIMATE_PURITY_B * s)


def load_protein_targets() -> tuple[pd.DataFrame, dict]:
    mat = pd.read_csv(PROTEIN, sep="\t").set_index("idx")
    mat.index = mat.index.astype(str)
    info = {}
    cols = {}
    for symbol, ens in TARGETS.items():
        raw = match_ensembl_row(mat.index, ens)
        if raw is None:
            raise RuntimeError(f"{symbol} {ens} not in protein matrix")
        s = pd.to_numeric(mat.loc[raw], errors="coerce")
        if isinstance(s, pd.DataFrame):
            s = s.iloc[0]
        s.name = symbol
        cols[symbol] = s
        info[symbol] = {
            "ensembl": ens,
            "row": raw,
            "n": int(s.notna().sum()),
            "n_na": int(s.isna().sum()),
            "min": float(s.min()) if s.notna().any() else np.nan,
            "median": float(s.median()) if s.notna().any() else np.nan,
            "max": float(s.max()) if s.notna().any() else np.nan,
        }
    return pd.DataFrame(cols), info


def load_pheno() -> pd.DataFrame:
    ph = pd.read_csv(PHENO, sep="\t").set_index("idx")
    ph.index.name = "case_id"
    keep = (
        PRIMARY_XCELL
        + EXPLORATORY_XCELL
        + ["WES_purity", "WGS_purity", "ESTIMATE_ESTIMATEScore", "ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore"]
    )
    missing = [c for c in keep if c not in ph.columns]
    if missing:
        raise RuntimeError(f"missing phenotype columns: {missing}")
    out = ph[keep].apply(pd.to_numeric, errors="coerce")
    out["ESTIMATE_TumorPurity"] = estimate_tumor_purity(out["ESTIMATE_ESTIMATEScore"])
    return out


def load_meta() -> pd.DataFrame:
    meta = pd.read_csv(META, sep="\t")
    meta = meta[meta["case_id"] != "data_type"].copy()
    return meta.set_index("case_id")


def assoc_row(protein: str, score: str, family: str, x, y, purity_name: str | None, z=None) -> dict:
    marg = spearman_pair(pd.Series(x), pd.Series(y))
    rec = {
        "protein": protein,
        "score": score,
        "family": family,
        "purity": purity_name or "none",
        "n_marginal": marg["n"],
        "marginal_rho": marg["rho"],
        "marginal_p": marg["p"],
        "n_partial": np.nan,
        "partial_spearman_rho": np.nan,
        "partial_p": np.nan,
        "partial_ci_low": np.nan,
        "partial_ci_high": np.nan,
        "ols_residual_spearman_rho": np.nan,
        "ols_residual_spearman_p": np.nan,
        "n_ols_residual": np.nan,
    }
    if z is not None:
        part = partial_spearman(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
        ols = ols_residual_spearman(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float))
        rec.update(
            {
                "n_partial": part["n"],
                "partial_spearman_rho": part["rho"],
                "partial_p": part["p"],
                "partial_ci_low": part["ci_low"],
                "partial_ci_high": part["ci_high"],
                "ols_residual_spearman_rho": ols["rho"],
                "ols_residual_spearman_p": ols["p"],
                "n_ols_residual": ols["n"],
            }
        )
    return rec


def fmt(x, digits=3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def plot_scatter_grid(df: pd.DataFrame, xcols: list[str], ycols: list[str], title: str, path: Path, xlabel_suffix=""):
    fig, axes = plt.subplots(len(xcols), len(ycols), figsize=(5.2 * len(ycols), 4.4 * len(xcols)), squeeze=False)
    for i, xc in enumerate(xcols):
        for j, yc in enumerate(ycols):
            ax = axes[i, j]
            sub = df[[xc, yc]].dropna()
            ax.scatter(sub[xc], sub[yc], s=22, alpha=0.75, c="#2c5f8a", edgecolors="none")
            if len(sub) >= 4:
                res = spearman_pair(sub[xc], sub[yc])
                ax.set_title(f"{xc} vs {yc}\nρ={fmt(res['rho'])}  p={fmt_p(res['p'])}  n={res['n']}", fontsize=10)
            else:
                ax.set_title(f"{xc} vs {yc} (n<{len(sub)})", fontsize=10)
            ax.set_xlabel(xc + xlabel_suffix)
            ax.set_ylabel(yc)
    fig.suptitle(title, y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_rho_before_after(assoc: pd.DataFrame, path: Path):
    sub = assoc[(assoc["family"] == "primary") & (assoc["purity"] == PRIMARY_PURITY)].copy()
    labels = [f"{r.protein}\n{r.score.replace('xCell_', '')}" for r in sub.itertuples()]
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    w = 0.35
    ax.bar(x - w / 2, sub["marginal_rho"], w, label="marginal Spearman", color="#4c78a8")
    ax.bar(x + w / 2, sub["partial_spearman_rho"], w, label=f"partial Spearman | {PRIMARY_PURITY}", color="#f58518")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Spearman ρ")
    ax.set_title("CPTAC LSCC protein vs xCell: before vs WES-purity residual")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_purity_context(df: pd.DataFrame, path: Path):
    feats = ["TACSTD2", "CLDN4", "xCell_T_cell_CD8+", "xCell_immune_score"]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
    for ax, feat in zip(axes.ravel(), feats):
        sub = df[[PRIMARY_PURITY, feat]].dropna()
        ax.scatter(sub[PRIMARY_PURITY], sub[feat], s=22, alpha=0.75, c="#3d6b4f", edgecolors="none")
        res = spearman_pair(sub[PRIMARY_PURITY], sub[feat])
        ax.set_title(f"{feat} vs {PRIMARY_PURITY}\nρ={fmt(res['rho'])}  p={fmt_p(res['p'])}  n={res['n']}", fontsize=10)
        ax.set_xlabel(PRIMARY_PURITY)
        ax.set_ylabel(feat)
    fig.suptitle("Purity context (why residualizing matters)", y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    if not PROTEIN.exists() or not PHENO.exists():
        raise SystemExit("Missing data. Run scripts/w200/CPTAC_LSCC_xcell/download.py first.")

    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    protein, prot_info = load_protein_targets()
    pheno = load_pheno()
    meta = load_meta()
    df = protein.join(pheno, how="inner")
    if df.shape[0] != 108:
        raise RuntimeError(f"expected 108 overlapping tumors, got {df.shape[0]}")

    # Sample-level table with WES rank-residuals for the primary pairs.
    sample = df.copy()
    sample = sample.join(meta[["Age", "Sex", "Stage"]], how="left")
    z = sample[PRIMARY_PURITY].to_numpy(float)
    for col in ["TACSTD2", "CLDN4"] + PRIMARY_XCELL:
        x = sample[col].to_numpy(float)
        resid = np.full(len(sample), np.nan)
        mask = complete_mask(x, z)
        if mask.sum() >= 5 and np.nanstd(z[mask]) > 0:
            resid[mask] = residualize(stats.rankdata(x[mask]), stats.rankdata(z[mask]))
        sample[f"{col}__resid_rank_{PRIMARY_PURITY}"] = resid
    sample.index.name = "case_id"
    sample.to_csv(TABLES / "sample_level.tsv", sep="\t")

    rows = []
    # Protein vs purity (confounding check, not a primary endpoint).
    for prot in TARGETS:
        for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
            rows.append(assoc_row(prot, pur, "purity_context", df[prot], df[pur], None))
    for score in PRIMARY_XCELL + EXPLORATORY_XCELL:
        for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
            rows.append(assoc_row("xCell", score if score != pur else score, "purity_context", df[score], df[pur], None))
            rows[-1]["protein"] = score
            rows[-1]["score"] = pur

    # Primary and exploratory protein × xCell, each purity residual.
    for prot in TARGETS:
        for score in PRIMARY_XCELL:
            rows.append(assoc_row(prot, score, "primary", df[prot], df[score], None))
            for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
                rows.append(assoc_row(prot, score, "primary", df[prot], df[score], pur, df[pur]))
        for score in EXPLORATORY_XCELL:
            rows.append(assoc_row(prot, score, "exploratory_cd8_subset", df[prot], df[score], None))
            rows.append(
                assoc_row(prot, score, "exploratory_cd8_subset", df[prot], df[score], PRIMARY_PURITY, df[PRIMARY_PURITY])
            )

    assoc = pd.DataFrame(rows)

    # BH within the 4 primary partial-Spearman tests on WES purity.
    primary_wes = (assoc["family"] == "primary") & (assoc["purity"] == PRIMARY_PURITY)
    assoc.loc[primary_wes, "partial_fdr_primary_wes"] = bh_fdr(assoc.loc[primary_wes, "partial_p"])
    # BH within the 4 primary marginal tests (for comparison with PR23).
    primary_marg = (assoc["family"] == "primary") & (assoc["purity"] == "none")
    assoc.loc[primary_marg, "marginal_fdr_primary"] = bh_fdr(assoc.loc[primary_marg, "marginal_p"])

    assoc.to_csv(TABLES / "associations.tsv", sep="\t", index=False)

    # Compact primary table.
    compact = assoc[(assoc["family"] == "primary") & (assoc["purity"].isin(["none", PRIMARY_PURITY]))].copy()
    compact.to_csv(TABLES / "primary_wes_partial.tsv", sep="\t", index=False)

    # Residual scatter frame.
    resid_df = pd.DataFrame(
        {
            "TACSTD2": sample[f"TACSTD2__resid_rank_{PRIMARY_PURITY}"],
            "CLDN4": sample[f"CLDN4__resid_rank_{PRIMARY_PURITY}"],
            "xCell_T_cell_CD8+": sample[f"xCell_T_cell_CD8+__resid_rank_{PRIMARY_PURITY}"],
            "xCell_immune_score": sample[f"xCell_immune_score__resid_rank_{PRIMARY_PURITY}"],
        }
    )
    plot_scatter_grid(
        df,
        ["TACSTD2", "CLDN4"],
        PRIMARY_XCELL,
        "Marginal: protein vs xCell (no purity residual)",
        FIGS / "fig1_marginal_scatter.png",
    )
    plot_scatter_grid(
        resid_df,
        ["TACSTD2", "CLDN4"],
        PRIMARY_XCELL,
        f"Purity residual: rank-residuals after {PRIMARY_PURITY}",
        FIGS / "fig2_wes_residual_scatter.png",
        xlabel_suffix=" residual",
    )
    plot_rho_before_after(assoc, FIGS / "fig3_rho_before_after.png")
    plot_purity_context(df, FIGS / "fig4_purity_context.png")

    def pick(protein, score, purity):
        hit = assoc[(assoc["protein"] == protein) & (assoc["score"] == score) & (assoc["purity"] == purity)]
        if hit.empty:
            return {}
        return hit.iloc[0].to_dict()

    summary = {
        "cohort": "CPTAC LSCC / LUSC",
        "citation": "Satpathy et al. Cell 2021 PMID 34358469; freeze data_freeze_v1.2_reorganized",
        "treatment": "treatment-naive resected LSCC; no prior chemo/RT; no ICI labels",
        "n_tumors": int(df.shape[0]),
        "targets": prot_info,
        "purity_primary": PRIMARY_PURITY,
        "purity_n": {
            "WES_purity": int(df["WES_purity"].notna().sum()),
            "WGS_purity": int(df["WGS_purity"].notna().sum()),
            "ESTIMATE_TumorPurity": int(df["ESTIMATE_TumorPurity"].notna().sum()),
        },
        "primary_tests": {
            f"{p}__{s}": {
                "marginal": pick(p, s, "none"),
                "partial_WES": pick(p, s, PRIMARY_PURITY),
                "partial_WGS": pick(p, s, "WGS_purity"),
                "partial_ESTIMATE": pick(p, s, "ESTIMATE_TumorPurity"),
            }
            for p in TARGETS
            for s in PRIMARY_XCELL
        },
        "purity_context": {
            f"{p}__{pur}": pick(p, pur, "none")
            for p in list(TARGETS) + PRIMARY_XCELL
            for pur in [PRIMARY_PURITY]
        },
        "notes": [
            "xCell scores are RNA-derived; proteins are TMT. Cross-layer correlation.",
            "WES/WGS purity are DNA-based and are the honest residual covariates.",
            "ESTIMATE TumorPurity is a cosine transform of ESTIMATEScore and is circular with immune content.",
            "CLDN4 protein missingness is TMT dropout (30/108), not a join error.",
            "This folder is a focused recompute of the PR23 xCell CD8/immune slice plus purity residual.",
        ],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    # JSON cannot serialize numpy types from Series; coerce.
    def _clean(obj):
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_clean(v) for v in obj]
        if isinstance(obj, (np.floating,)):
            return float(obj) if np.isfinite(obj) else None
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if obj is pd.NA or (isinstance(obj, float) and not np.isfinite(obj)):
            return None
        return obj

    (TABLES / "run_summary.json").write_text(json.dumps(_clean(summary), indent=2) + "\n")

    # Seed data dir from /tmp cache if download.py was not used in this workspace.
    DATA.mkdir(parents=True, exist_ok=True)
    for src in (PROTEIN, PHENO, META):
        dest = DATA / src.name
        if src.exists() and not dest.exists():
            shutil.copy2(src, dest)

    print("Wrote", TABLES)
    print("Wrote", FIGS)
    print(assoc[primary_wes][["protein", "score", "n_partial", "marginal_rho", "partial_spearman_rho", "partial_p"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

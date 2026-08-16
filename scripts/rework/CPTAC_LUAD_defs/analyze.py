#!/usr/bin/env python3
"""CPTAC LUAD TACSTD2 protein vs MCP-counter / ESTIMATE / CIBERSORT, WES residual.

Rework of the PR6 LUAD protein–immune slice. PR6's headline was TACSTD2
protein vs xCell immune ρ = −0.309. This script keeps that as a
reproduction check and asks whether the same protein associates with
public ESTIMATE, public CIBERSORT, and Becht-2016 MCP-counter after a
WES purity residual.

Treatment-naive surgical LUAD (Gillette et al. Cell 2020). No ICI labels.
"""

from __future__ import annotations

import json
import math
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
HERE = Path(__file__).resolve().parent
DATA = ROOT / "data" / "rework" / "CPTAC_LUAD_defs"
OUT = ROOT / "results" / "rework" / "CPTAC_LUAD_defs"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

PROTEIN = DATA / "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
RNA = DATA / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
PHENO = DATA / "LUAD_phenotype.txt"
META = DATA / "LUAD_meta.txt"
MCP_GENES = HERE / "mcp_counter_genes.tsv"

TARGETS = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
}

# Rework primary immune definitions (not xCell).
PRIMARY_SCORES = [
    "ESTIMATE_ImmuneScore",
    "MCPcounter_T_cells",
    "MCPcounter_Cytotoxic_lymphocytes",
    "CIBERSORT_T_cell_CD8+",
]
# Reproduction of PR6 (not in the rework FDR family).
REPRO_XCELL = ["xCell_immune_score", "xCell_T_cell_CD8+"]
# Extra public / computed scores (exploratory).
EXPLORATORY_SCORES = [
    "ESTIMATE_StromalScore",
    "MCPcounter_CD8_T_cells",
    "MCPcounter_NK_cells",
    "MCPcounter_B_lineage",
    "MCPcounter_Monocytic_lineage",
    "MCPcounter_Neutrophils",
    "MCPcounter_Fibroblasts",
    "CIBERSORT_leukocyte_sum",
    "CIBERSORT_T_cell_regulatory_(Tregs)",
    "CIBERSORT_Macrophage_M2",
    "CIBERSORT_Neutrophil",
]

PRIMARY_PURITY = "WES_purity"
SENSITIVITY_PURITY = ["WGS_purity", "ESTIMATE_TumorPurity"]

ESTIMATE_PURITY_A = 0.6049872018
ESTIMATE_PURITY_B = 0.0001467884

CIBERSORT_COLS = [
    "CIBERSORT_B_cell_naive",
    "CIBERSORT_B_cell_memory",
    "CIBERSORT_B_cell_plasma",
    "CIBERSORT_T_cell_CD8+",
    "CIBERSORT_T_cell_CD4+_naive",
    "CIBERSORT_T_cell_CD4+_memory_resting",
    "CIBERSORT_T_cell_CD4+_memory_activated",
    "CIBERSORT_T_cell_follicular_helper",
    "CIBERSORT_T_cell_regulatory_(Tregs)",
    "CIBERSORT_T_cell_gamma_delta",
    "CIBERSORT_NK_cell_resting",
    "CIBERSORT_NK_cell_activated",
    "CIBERSORT_Monocyte",
    "CIBERSORT_Macrophage_M0",
    "CIBERSORT_Macrophage_M1",
    "CIBERSORT_Macrophage_M2",
    "CIBERSORT_Myeloid_dendritic_cell_resting",
    "CIBERSORT_Myeloid_dendritic_cell_activated",
    "CIBERSORT_Mast_cell_activated",
    "CIBERSORT_Mast_cell_resting",
    "CIBERSORT_Eosinophil",
    "CIBERSORT_Neutrophil",
]


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
    half = stats.norm.ppf(0.975) / math.sqrt(max(n - 4, 1))
    lo, hi = np.tanh([zf - half, zf + half])
    return {"n": n, "rho": rho, "p": p, "ci_low": float(lo), "ci_high": float(hi)}


def ols_residual_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
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


def slug_population(name: str) -> str:
    return "MCPcounter_" + name.replace(" ", "_")


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


def compute_mcp_counter() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Becht 2016 MCP-counter on log2 RNA: mean of log2 = log2 geometric mean."""
    genes = pd.read_csv(MCP_GENES, sep="\t")
    rna = pd.read_csv(RNA, sep="\t").set_index("idx")
    rna.index = rna.index.astype(str)
    bare_to_row = {str(i).split(".")[0]: str(i) for i in rna.index}
    coverage_rows = []
    scores = {}
    for pop, sub in genes.groupby("population", sort=False):
        official = [str(e) for e in sub["ensembl"] if pd.notna(e)]
        present = [e for e in official if e in bare_to_row]
        missing = [e for e in official if e not in bare_to_row]
        if not present:
            raise RuntimeError(f"no MCP-counter genes present for {pop}")
        rows = [bare_to_row[e] for e in present]
        block = rna.loc[rows].apply(pd.to_numeric, errors="coerce")
        # Official score is geometric mean of linear expression.
        # Freeze RNA is already log2; arithmetic mean of log2 is that transform.
        score = block.mean(axis=0, skipna=True)
        col = slug_population(pop)
        scores[col] = score
        coverage_rows.append(
            {
                "population": pop,
                "score_column": col,
                "n_official_ensembl": len(official),
                "n_present": len(present),
                "n_missing": len(missing),
                "missing_ensembl": ",".join(missing),
                "present_ensembl": ",".join(present),
                "method": "mean_of_log2_RSEM_UQ1500 (== log2 geometric mean of linear)",
            }
        )
    return pd.DataFrame(scores), pd.DataFrame(coverage_rows)


def load_pheno() -> tuple[pd.DataFrame, dict]:
    ph = pd.read_csv(PHENO, sep="\t").set_index("idx")
    ph.index.name = "case_id"
    need = (
        [
            "ESTIMATE_ImmuneScore",
            "ESTIMATE_StromalScore",
            "ESTIMATE_ESTIMATEScore",
            "CIBERSORT_T_cell_CD8+",
            "CIBERSORT_T_cell_regulatory_(Tregs)",
            "CIBERSORT_Macrophage_M2",
            "CIBERSORT_Neutrophil",
            "xCell_immune_score",
            "xCell_T_cell_CD8+",
            "WES_purity",
            "WGS_purity",
        ]
        + CIBERSORT_COLS
    )
    missing = [c for c in need if c not in ph.columns]
    if missing:
        raise RuntimeError(f"missing phenotype columns: {missing}")
    out = ph[need].apply(pd.to_numeric, errors="coerce")
    out = out.loc[:, ~out.columns.duplicated()]
    out["ESTIMATE_TumorPurity"] = estimate_tumor_purity(out["ESTIMATE_ESTIMATEScore"])
    cib = out[CIBERSORT_COLS]
    out["CIBERSORT_leukocyte_sum"] = cib.sum(axis=1, skipna=False)
    cib_info = {
        "n_populations": len(CIBERSORT_COLS),
        "row_sum_min": float(out["CIBERSORT_leukocyte_sum"].min()),
        "row_sum_median": float(out["CIBERSORT_leukocyte_sum"].median()),
        "row_sum_max": float(out["CIBERSORT_leukocyte_sum"].max()),
        "note": (
            "Row sums are not 1 (median ~1.8). Freeze CIBERSORT is not "
            "relative-mode fractions. Treat as immunedeconv-style scores."
        ),
    }
    return out, cib_info


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
    fig, axes = plt.subplots(len(xcols), len(ycols), figsize=(4.6 * len(ycols), 4.0 * len(xcols)), squeeze=False)
    for i, xc in enumerate(xcols):
        for j, yc in enumerate(ycols):
            ax = axes[i, j]
            sub = df[[xc, yc]].dropna()
            ax.scatter(sub[xc], sub[yc], s=20, alpha=0.75, c="#2c5f8a", edgecolors="none")
            if len(sub) >= 4:
                res = spearman_pair(sub[xc], sub[yc])
                ax.set_title(f"{yc}\nρ={fmt(res['rho'])}  p={fmt_p(res['p'])}  n={res['n']}", fontsize=9)
            ax.set_xlabel(xc + xlabel_suffix, fontsize=8)
            ax.set_ylabel(yc, fontsize=8)
    fig.suptitle(title, y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_rho_before_after(assoc: pd.DataFrame, path: Path):
    sub = assoc[
        (assoc["family"] == "primary")
        & (assoc["purity"] == PRIMARY_PURITY)
        & (assoc["protein"] == "TACSTD2")
    ].copy()
    labels = [r.score.replace("MCPcounter_", "MCP ").replace("CIBERSORT_", "CIB ").replace("ESTIMATE_", "") for r in sub.itertuples()]
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    w = 0.35
    ax.bar(x - w / 2, sub["marginal_rho"], w, label="marginal Spearman", color="#4c78a8")
    ax.bar(x + w / 2, sub["partial_spearman_rho"], w, label=f"partial Spearman | {PRIMARY_PURITY}", color="#f58518")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Spearman ρ")
    ax.set_title("LUAD TACSTD2 protein vs public immune defs: before vs WES residual")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_purity_context(df: pd.DataFrame, path: Path):
    feats = ["TACSTD2"] + PRIMARY_SCORES
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2))
    for ax, feat in zip(axes.ravel(), feats):
        sub = df[[PRIMARY_PURITY, feat]].dropna()
        ax.scatter(sub[PRIMARY_PURITY], sub[feat], s=18, alpha=0.75, c="#3d6b4f", edgecolors="none")
        res = spearman_pair(sub[PRIMARY_PURITY], sub[feat])
        ax.set_title(f"{feat}\nρ={fmt(res['rho'])}  p={fmt_p(res['p'])}  n={res['n']}", fontsize=9)
        ax.set_xlabel(PRIMARY_PURITY, fontsize=8)
        ax.set_ylabel(feat, fontsize=8)
    axes.ravel()[-1].axis("off")
    fig.suptitle("Purity context (WES vs TACSTD2 protein and primary immune defs)", y=1.01, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_method_compare(assoc: pd.DataFrame, path: Path):
    """TACSTD2 protein: xCell reproduction vs rework defs, marginal + WES partial."""
    scores = REPRO_XCELL + PRIMARY_SCORES
    rows = []
    for score in scores:
        for purity, label in [("none", "marginal"), (PRIMARY_PURITY, "WES partial")]:
            hit = assoc[
                (assoc["protein"] == "TACSTD2")
                & (assoc["score"] == score)
                & (assoc["purity"] == purity)
            ]
            if hit.empty:
                continue
            r = hit.iloc[0]
            rows.append(
                {
                    "score": score,
                    "adjust": label,
                    "rho": r["marginal_rho"] if purity == "none" else r["partial_spearman_rho"],
                    "p": r["marginal_p"] if purity == "none" else r["partial_p"],
                }
            )
    plot_df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    labels = [s.replace("MCPcounter_", "MCP ").replace("CIBERSORT_", "CIB ").replace("ESTIMATE_", "").replace("xCell_", "xCell ") for s in scores]
    x = np.arange(len(scores))
    w = 0.35
    marg = plot_df[plot_df.adjust == "marginal"].set_index("score").reindex(scores)
    part = plot_df[plot_df.adjust == "WES partial"].set_index("score").reindex(scores)
    ax.bar(x - w / 2, marg["rho"], w, label="marginal", color="#4c78a8")
    ax.bar(x + w / 2, part["rho"], w, label="WES partial", color="#f58518")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("Spearman ρ (TACSTD2 protein)")
    ax.set_title("PR6 xCell reproduction vs MCP-counter / ESTIMATE / CIBERSORT")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


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


def main() -> int:
    for path in (PROTEIN, RNA, PHENO):
        if not path.exists():
            raise SystemExit(f"Missing {path}. Run scripts/rework/CPTAC_LUAD_defs/download.py first.")

    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    protein, prot_info = load_protein_targets()
    mcp, mcp_cov = compute_mcp_counter()
    pheno, cib_info = load_pheno()
    meta = load_meta()

    df = protein.join(pheno, how="inner").join(mcp, how="inner")
    if df.shape[0] != 110:
        raise RuntimeError(f"expected 110 overlapping tumors, got {df.shape[0]}")

    mcp_cov.to_csv(TABLES / "mcp_counter_coverage.tsv", sep="\t", index=False)

    sample = df.copy()
    keep_meta = [c for c in ["Age", "Sex", "Stage"] if c in meta.columns]
    if keep_meta:
        sample = sample.join(meta[keep_meta], how="left")
    z = sample[PRIMARY_PURITY].to_numpy(float)
    resid_cols = ["TACSTD2", "CLDN4"] + PRIMARY_SCORES + REPRO_XCELL
    for col in resid_cols:
        x = sample[col].to_numpy(float)
        resid = np.full(len(sample), np.nan)
        mask = complete_mask(x, z)
        if mask.sum() >= 5 and np.nanstd(z[mask]) > 0:
            resid[mask] = residualize(stats.rankdata(x[mask]), stats.rankdata(z[mask]))
        sample[f"{col}__resid_rank_{PRIMARY_PURITY}"] = resid
    sample.index.name = "case_id"
    sample.to_csv(TABLES / "sample_level.tsv", sep="\t")

    rows = []
    score_universe = PRIMARY_SCORES + REPRO_XCELL + EXPLORATORY_SCORES
    for prot in TARGETS:
        for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
            rows.append(assoc_row(prot, pur, "purity_context", df[prot], df[pur], None))
    for score in score_universe:
        for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
            rec = assoc_row(score, pur, "purity_context", df[score], df[pur], None)
            rec["protein"] = score
            rec["score"] = pur
            rows.append(rec)

    for prot in TARGETS:
        for score in PRIMARY_SCORES:
            rows.append(assoc_row(prot, score, "primary", df[prot], df[score], None))
            for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
                rows.append(assoc_row(prot, score, "primary", df[prot], df[score], pur, df[pur]))
        for score in REPRO_XCELL:
            rows.append(assoc_row(prot, score, "xcell_reproduction", df[prot], df[score], None))
            for pur in [PRIMARY_PURITY] + SENSITIVITY_PURITY:
                rows.append(assoc_row(prot, score, "xcell_reproduction", df[prot], df[score], pur, df[pur]))
        for score in EXPLORATORY_SCORES:
            rows.append(assoc_row(prot, score, "exploratory", df[prot], df[score], None))
            rows.append(assoc_row(prot, score, "exploratory", df[prot], df[score], PRIMARY_PURITY, df[PRIMARY_PURITY]))

    assoc = pd.DataFrame(rows)

    # Primary FDR: TACSTD2 protein × 4 rework defs, WES partial only.
    primary_wes = (
        (assoc["family"] == "primary")
        & (assoc["purity"] == PRIMARY_PURITY)
        & (assoc["protein"] == "TACSTD2")
    )
    assoc.loc[primary_wes, "partial_fdr_primary_wes"] = bh_fdr(assoc.loc[primary_wes, "partial_p"])
    primary_marg = (
        (assoc["family"] == "primary")
        & (assoc["purity"] == "none")
        & (assoc["protein"] == "TACSTD2")
    )
    assoc.loc[primary_marg, "marginal_fdr_primary"] = bh_fdr(assoc.loc[primary_marg, "marginal_p"])

    # CLDN4 comparator FDR (same 4 defs, WES partial) — reported separately.
    cldn4_wes = (
        (assoc["family"] == "primary")
        & (assoc["purity"] == PRIMARY_PURITY)
        & (assoc["protein"] == "CLDN4")
    )
    assoc.loc[cldn4_wes, "partial_fdr_cldn4_wes"] = bh_fdr(assoc.loc[cldn4_wes, "partial_p"])

    assoc.to_csv(TABLES / "associations.tsv", sep="\t", index=False)
    compact = assoc[
        (assoc["family"].isin(["primary", "xcell_reproduction"]))
        & (assoc["purity"].isin(["none", PRIMARY_PURITY]))
    ].copy()
    compact.to_csv(TABLES / "primary_wes_partial.tsv", sep="\t", index=False)

    resid_df = pd.DataFrame({c: sample[f"{c}__resid_rank_{PRIMARY_PURITY}"] for c in resid_cols})
    plot_scatter_grid(
        df,
        ["TACSTD2"],
        PRIMARY_SCORES,
        "Marginal: TACSTD2 protein vs MCP-counter / ESTIMATE / CIBERSORT",
        FIGS / "fig1_marginal_scatter.png",
    )
    plot_scatter_grid(
        resid_df,
        ["TACSTD2"],
        PRIMARY_SCORES,
        f"WES rank-residual: TACSTD2 protein vs primary immune defs",
        FIGS / "fig2_wes_residual_scatter.png",
        xlabel_suffix=" residual",
    )
    plot_rho_before_after(assoc, FIGS / "fig3_rho_before_after.png")
    plot_purity_context(df, FIGS / "fig4_purity_context.png")
    plot_method_compare(assoc, FIGS / "fig5_xcell_vs_rework_defs.png")

    def pick(protein, score, purity):
        hit = assoc[(assoc["protein"] == protein) & (assoc["score"] == score) & (assoc["purity"] == purity)]
        if hit.empty:
            return {}
        return hit.iloc[0].to_dict()

    key_rows = []
    for prot in TARGETS:
        for score in PRIMARY_SCORES + REPRO_XCELL:
            m = pick(prot, score, "none")
            w = pick(prot, score, PRIMARY_PURITY)
            key_rows.append(
                {
                    "protein": prot,
                    "score": score,
                    "family": "primary" if score in PRIMARY_SCORES else "xcell_reproduction",
                    "n_marginal": m.get("n_marginal"),
                    "marginal_rho": m.get("marginal_rho"),
                    "marginal_p": m.get("marginal_p"),
                    "marginal_fdr_primary": m.get("marginal_fdr_primary"),
                    "n_WES": w.get("n_partial"),
                    "partial_rho_WES": w.get("partial_spearman_rho"),
                    "partial_p_WES": w.get("partial_p"),
                    "partial_fdr_WES": w.get("partial_fdr_primary_wes") if prot == "TACSTD2" else w.get("partial_fdr_cldn4_wes"),
                    "partial_ci_low": w.get("partial_ci_low"),
                    "partial_ci_high": w.get("partial_ci_high"),
                    "ols_residual_rho_WES": w.get("ols_residual_spearman_rho"),
                    "ols_residual_p_WES": w.get("ols_residual_spearman_p"),
                }
            )
    key = pd.DataFrame(key_rows)
    key.to_csv(TABLES / "key_findings.tsv", sep="\t", index=False)

    summary = {
        "cohort": "CPTAC LUAD",
        "citation": "Gillette et al. Cell 2020 PMID 32649874; freeze data_freeze_v1.2_reorganized",
        "treatment": "treatment-naive resected LUAD; no ICI labels",
        "n_tumors": int(df.shape[0]),
        "targets": prot_info,
        "question": (
            "Does TACSTD2 protein vs immune anti-correlation (PR6 xCell ρ=-0.309) "
            "hold for public ESTIMATE, public CIBERSORT, and Becht MCP-counter "
            "after WES purity residual?"
        ),
        "purity_primary": PRIMARY_PURITY,
        "purity_n": {
            "WES_purity": int(df["WES_purity"].notna().sum()),
            "WGS_purity": int(df["WGS_purity"].notna().sum()),
            "ESTIMATE_TumorPurity": int(df["ESTIMATE_TumorPurity"].notna().sum()),
        },
        "cibersort": cib_info,
        "mcp_counter": mcp_cov.to_dict(orient="records"),
        "primary_scores": PRIMARY_SCORES,
        "primary_tests_TACSTD2": {
            s: {
                "marginal": pick("TACSTD2", s, "none"),
                "partial_WES": pick("TACSTD2", s, PRIMARY_PURITY),
                "partial_WGS": pick("TACSTD2", s, "WGS_purity"),
                "partial_ESTIMATE": pick("TACSTD2", s, "ESTIMATE_TumorPurity"),
            }
            for s in PRIMARY_SCORES
        },
        "xcell_reproduction_TACSTD2": {
            s: {
                "marginal": pick("TACSTD2", s, "none"),
                "partial_WES": pick("TACSTD2", s, PRIMARY_PURITY),
            }
            for s in REPRO_XCELL
        },
        "purity_context": {
            f"{p}__{PRIMARY_PURITY}": pick(p, PRIMARY_PURITY, "none")
            for p in list(TARGETS) + PRIMARY_SCORES + REPRO_XCELL
        },
        "notes": [
            "xCell is the PR6 definition; it is reproduced but is not in the rework FDR family.",
            "ESTIMATE and CIBERSORT are freeze phenotype columns (public).",
            "MCP-counter is not a freeze column (standalone files HTTP 403). Computed here from public RNA + Becht 2016 genes.",
            "WES/WGS purity are DNA-based and are the honest residual covariates.",
            "ESTIMATE TumorPurity is a cosine of ESTIMATEScore and is circular with ImmuneScore.",
            "CIBERSORT row sums are not 1; do not treat as relative fractions.",
            "CLDN4 protein missingness is TMT dropout (31/110), not a join error.",
            "No ICI labels. Do not read correlations as immunotherapy outcomes.",
            "LSCC TACSTD2 protein vs xCell was null (PR23/PR84); this folder is LUAD only.",
        ],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (TABLES / "run_summary.json").write_text(json.dumps(_clean(summary), indent=2) + "\n")

    print("Wrote", TABLES)
    print("Wrote", FIGS)
    show = assoc[primary_wes][
        ["protein", "score", "n_partial", "marginal_rho", "partial_spearman_rho", "partial_p", "partial_fdr_primary_wes"]
    ]
    print(show.to_string(index=False))
    print("--- xCell reproduction TACSTD2 ---")
    repro = assoc[
        (assoc["protein"] == "TACSTD2")
        & (assoc["family"] == "xcell_reproduction")
        & (assoc["purity"].isin(["none", PRIMARY_PURITY]))
    ]
    print(repro[["score", "purity", "n_marginal", "marginal_rho", "marginal_p", "partial_spearman_rho", "partial_p"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

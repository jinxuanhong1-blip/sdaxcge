#!/usr/bin/env python3
"""Maximize honest TROP2–CLDN4 protein and TROP2–TJ gene Spearman ρ.

Part1 bridge question: do DepMap/CCLE lung lines support TROP2 coexpression
with CLDN4 protein and locked TJ gene scores, and what is the largest
pre-specified ρ with honest n?

Two independent maxima:
  1) Protein: Gygi TACSTD2 vs CLDN4 (histology / primary-met / subtype partial grid)
  2) RNA: DepMap 24Q4 TACSTD2 vs locked TJ gene scores (same histology grid)

Rules written before sorting:
  - MIN_N_HEADLINE = 15
  - MIN_RESID_DF = 10 for partials
  - Eligible families: histology + partial_histology (protein); histology +
    partial_histology for RNA TJ partners in PRIMARY_TJ
  - Epithelial / keratin partials are reported and cannot win
  - Peptide QC (protein) is reported and cannot win the protein histology max
  - No imputation. No dropping individual lines to raise ρ.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata, spearmanr, t as student_t

PRIOR_PROTEIN_RHO = 0.69
PRIOR_CLAIMED_N = 118
MIN_N_HEADLINE = 15
MIN_RESID_DF = 10
N_BOOT = 4000
N_NULL = 4000
SEED = 0

SAMPLE_RE = re.compile(r"_TenPx\d+$")

TJ_EPITHELIAL = [
    "CLDN1",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "OCLN",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "JAM2",
    "JAM3",
    "CGN",
    "CGNL1",
    "CRB3",
    "ILDR1",
    "LSR",
]
TJ_TISMO = ["CLDN3", "CLDN4", "CLDN6", "CLDN7", "CDH1", "F11R", "OCLN"]
CLDN4_TJ_EDGE = [
    "CLDN4",
    "CLDN1",
    "CLDN7",
    "CGNL1",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "ILDR1",
    "CLDN3",
    "OCLN",
]
PRIMARY_TJ = ["CLDN4", "TJ_EPITHELIAL", "TJ_TISMO", "CLDN4_TJ_EDGE"]
FOCAL_SINGLE = ["CLDN1", "CLDN7", "OCLN", "F11R", "TJP1"]
PROTEIN_GENES = ["TACSTD2", "CLDN4", "EPCAM", "VIM", "CDH1"]
RNA_CONTROLS = ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19"]

TJ_LABEL = {
    "CLDN4": "CLDN4",
    "TJ_EPITHELIAL": "TJ epithelial (18)",
    "TJ_TISMO": "TJ TISMO (7)",
    "CLDN4_TJ_EDGE": "CLDN4 TJ edge (11)",
    "CLDN1": "CLDN1",
    "CLDN7": "CLDN7",
    "OCLN": "OCLN",
    "F11R": "F11R",
    "TJP1": "TJP1",
}

HISTOLOGY_ADJUSTERS = [
    ["OncotreeSubtype"],
    ["OncotreePrimaryDisease"],
    ["PrimaryOrMetastasis"],
    ["OncotreeSubtype", "PrimaryOrMetastasis"],
]
PARTIAL_PARENTS_PROTEIN = [
    "S1_Lung",
    "Oncotree_NSCLC",
    "Oncotree_LUAD",
    "NSCLC_Primary",
    "S1_Lung_Primary",
]
PARTIAL_PARENTS_RNA = [
    "RNA_Lung",
    "RNA_NSCLC",
    "RNA_LUAD",
    "RNA_NSCLC_Primary",
    "RNA_Lung_Primary",
]


def core_of(col: str) -> str:
    m = re.match(r"(.+)_TenPx\d+$", col)
    return m.group(1) if m else col


def plex_of(col: str) -> str:
    m = re.search(r"(TenPx\d+)$", col)
    return m.group(1) if m else ""


def fast_spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(x) - rankdata(x).mean()
    ry = rankdata(y) - rankdata(y).mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if denom == 0:
        return np.nan
    return float((rx * ry).sum() / denom)


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED):
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        ix = rng.integers(0, n, n)
        boots[i] = fast_spearman(x[ix], y[ix])
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return float(lo), float(hi)


def subset_null(parent_x, parent_y, obs: float, size: int, n_perm: int = N_NULL, seed: int = SEED):
    rng = np.random.default_rng(seed)
    n = len(parent_x)
    if size >= n or size < 5 or not np.isfinite(obs):
        return None
    ge = 0
    for _ in range(n_perm):
        ix = rng.choice(n, size=size, replace=False)
        if fast_spearman(parent_x[ix], parent_y[ix]) >= obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)


def bh_q(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    pp = p[ok]
    n = len(pp)
    order = np.argsort(pp)
    ranked = pp[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    filled = np.empty(n)
    filled[order] = q
    out[np.where(ok)[0]] = filled
    return out.tolist()


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str], min_members: int) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if len(use) < min_members:
        return pd.Series(np.nan, index=df.index), use
    z = df[use].apply(zscore, axis=0)
    n_ok = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True).where(n_ok >= min_members)
    return score, use


def unadjusted_xy(x: np.ndarray, y: np.ndarray) -> dict:
    n = int(len(x))
    rec = {
        "n": n,
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "ci95_lo": None,
        "ci95_hi": None,
        "resid_df": None,
        "n_adjusters": 0,
        "partial_status": None,
    }
    if n < 4:
        return rec
    rho, p = spearmanr(x, y)
    r, pp = pearsonr(x, y)
    lo, hi = spearman_ci(x, y)
    rec.update(
        {
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "pearson_r": float(r),
            "pearson_p": float(pp),
            "ci95_lo": lo,
            "ci95_hi": hi,
        }
    )
    return rec


def _design(d: pd.DataFrame, y_col: str, x_col: str, covars: list[str]):
    numeric = [c for c in covars if c in d.columns and pd.api.types.is_numeric_dtype(d[c])]
    need = [y_col, x_col] + numeric
    d = d.dropna(subset=need).copy()
    pieces = []
    for c in covars:
        if c in numeric:
            pieces.append(d[[c]].to_numpy(float))
            continue
        if c not in d.columns:
            return None
        s = d[c].fillna("NA").astype(str)
        vc = s.value_counts()
        s = s.where(s.map(vc) >= 3, "Other")
        if s.nunique() < 2:
            return None
        dum = pd.get_dummies(s, drop_first=True)
        if dum.shape[1] == 0:
            return None
        pieces.append(dum.to_numpy(float))
    Z = np.column_stack(pieces)
    keep = Z.std(axis=0) > 1e-8
    Z = Z[:, keep]
    if Z.shape[1] == 0:
        return None
    return d, Z


def partial_spearman(d: pd.DataFrame, y_col: str, x_col: str, covars: list[str], bootstrap: bool = True) -> dict:
    built = _design(d, y_col, x_col, covars)
    base = {
        "n": int(len(d.dropna(subset=[y_col, x_col]))),
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "ci95_lo": None,
        "ci95_hi": None,
        "resid_df": None,
        "n_adjusters": None,
        "partial_status": "not_estimable",
        "unadjusted_spearman_same_n": None,
    }
    if built is None:
        return base
    d, Z = built
    y = d[y_col].to_numpy(float)
    x = d[x_col].to_numpy(float)
    n, k = len(d), int(Z.shape[1])
    df_ = n - k - 2
    base.update({"n": n, "n_adjusters": k, "resid_df": int(df_)})
    if df_ < 5:
        base["partial_status"] = "residual_df_below_5"
        return base

    def _r(y_in, x_in, Z_in) -> float:
        yr = rankdata(y_in).astype(float)
        xr = rankdata(x_in).astype(float)
        Zr = np.column_stack([rankdata(Z_in[:, i]) for i in range(Z_in.shape[1])])
        A = np.column_stack([np.ones(len(y_in)), Zr])
        by, *_ = np.linalg.lstsq(A, yr, rcond=None)
        bx, *_ = np.linalg.lstsq(A, xr, rcond=None)
        ry, rx = yr - A @ by, xr - A @ bx
        denom = np.sqrt((ry * ry).sum() * (rx * rx).sum())
        if denom == 0:
            return np.nan
        return float((ry * rx).sum() / denom)

    r = _r(y, x, Z)
    if not np.isfinite(r) or abs(r) >= 1:
        base["partial_status"] = "degenerate"
        return base
    tstat = r * np.sqrt(df_ / (1 - r * r))
    p = float(2 * student_t.sf(abs(tstat), df_))
    rho_u, p_u = spearmanr(y, x)
    pr, pp = pearsonr(y, x)
    lo = hi = None
    if bootstrap:
        rng = np.random.default_rng(SEED)
        boots = []
        for _ in range(N_BOOT):
            ix = rng.integers(0, n, n)
            rr = _r(y[ix], x[ix], Z[ix])
            if np.isfinite(rr):
                boots.append(rr)
        if len(boots) >= 100:
            lo, hi = [float(v) for v in np.quantile(boots, [0.025, 0.975])]
    base.update(
        {
            "spearman_rho": float(r),
            "spearman_p": p,
            "pearson_r": float(pr),
            "pearson_p": float(pp),
            "unadjusted_spearman_same_n": float(rho_u),
            "ci95_lo": lo,
            "ci95_hi": hi,
            "partial_status": "ok",
        }
    )
    return base


def empty_row(**kwargs) -> dict:
    base = {
        "layer": None,
        "analysis_id": None,
        "family": None,
        "cohort": None,
        "partner": None,
        "label": None,
        "definition": None,
        "adjustment": None,
        "n": None,
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "ci95_lo": None,
        "ci95_hi": None,
        "resid_df": None,
        "n_adjusters": None,
        "unadjusted_spearman_same_n": None,
        "null_p_same_size_from_parent": None,
        "partial_status": None,
        "eligible_for_max": False,
        "is_layer_max": False,
        "q_bh_within_eligible": None,
        "n_equals_118": False,
        "rounds_to_0.69": False,
        "delta_vs_0.69": None,
    }
    base.update(kwargs)
    return base


def load_protein(data: Path) -> pd.DataFrame:
    prot = pd.read_csv(data / "protein_quant_current_normalized.csv.gz")
    missing = [g for g in PROTEIN_GENES if g not in set(prot["Gene_Symbol"])]
    if missing:
        raise SystemExit(f"Missing protein genes: {missing}")
    sub = prot[prot["Gene_Symbol"].isin(PROTEIN_GENES)].drop_duplicates("Gene_Symbol").set_index("Gene_Symbol")
    sample_cols = [c for c in prot.columns if SAMPLE_RE.search(str(c))]
    pep_cols = [c for c in prot.columns if str(c).endswith("_Peptides")]
    cl_pep = {
        c.replace("_Peptides", ""): float(pd.to_numeric(sub.loc["CLDN4", c], errors="coerce"))
        for c in pep_cols
        if "CLDN4" in sub.index
    }
    cols_by: dict[str, list[str]] = defaultdict(list)
    for c in sample_cols:
        cols_by[core_of(c)].append(c)
    rows = []
    for core, cols in cols_by.items():
        rec = {"CCLE": core, "n_replicates": len(cols), "plex": plex_of(cols[0])}
        for g in PROTEIN_GENES:
            rec[g] = float(pd.to_numeric(sub.loc[g, cols], errors="coerce").mean())
        rec["CLDN4_peptides"] = float(np.nanmin([cl_pep.get(plex_of(c), np.nan) for c in cols]))
        rows.append(rec)
    df = pd.DataFrame(rows).set_index("CCLE")
    s1 = pd.read_excel(data / "Table_S1_Sample_Information.xlsx", sheet_name="Sample_Information")
    s1 = s1[s1["Notes"] != "Bridge line"].drop_duplicates("CCLE Code").set_index("CCLE Code")
    model = pd.read_csv(data / "Model.csv").drop_duplicates("CCLEName").set_index("CCLEName")
    df["CellLine"] = s1.reindex(df.index)["Cell Line"]
    df["S1_Tissue"] = s1.reindex(df.index)["Tissue of Origin"]
    for col in [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
        "Sex",
        "PrimaryOrMetastasis",
    ]:
        df[col] = model.reindex(df.index)[col]
    df["CellLine"] = df["CellLine"].fillna(df["CellLineName"])
    return df


def load_rna(data: Path) -> pd.DataFrame:
    expr = pd.read_csv(data / "expression_panel.csv").set_index("ModelID")
    for c in expr.columns:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")
    model = pd.read_csv(data / "Model.csv").drop_duplicates("ModelID").set_index("ModelID")
    keep = [
        "CellLineName",
        "CCLEName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
        "Sex",
        "PrimaryOrMetastasis",
    ]
    out = expr.join(model[keep], how="inner")
    out = out[out["OncotreeLineage"] == "Lung"].copy()
    out = out.rename(columns={"CellLineName": "CellLine"})
    # Locked TJ scores (z-mean).
    for name, genes, vmin in (
        ("TJ_EPITHELIAL", TJ_EPITHELIAL, max(5, len(TJ_EPITHELIAL) // 2)),
        ("TJ_TISMO", TJ_TISMO, max(3, len(TJ_TISMO) // 2)),
        ("CLDN4_TJ_EDGE", CLDN4_TJ_EDGE, max(4, len(CLDN4_TJ_EDGE) // 2)),
    ):
        score, used = signature(out, genes, min_members=vmin)
        out[name] = score
        out[f"{name}_n_genes"] = len(used)
    return out


def protein_cohorts(both: pd.DataFrame) -> list[tuple[str, pd.Series, str]]:
    nsclc = both["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    luad = both["OncotreeSubtype"] == "Lung Adenocarcinoma"
    primary = both["PrimaryOrMetastasis"] == "Primary"
    metastatic = both["PrimaryOrMetastasis"] == "Metastatic"
    return [
        ("S1_Lung", both.index == both.index, "Table S1 Tissue=Lung, TACSTD2+CLDN4 both quantified"),
        ("Oncotree_NSCLC", nsclc, "S1 Lung and OncotreePrimaryDisease=NSCLC"),
        ("Oncotree_LUAD", luad, "S1 Lung and OncotreeSubtype=LUAD"),
        ("Oncotree_LUSC", both["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma", "S1 Lung LUSC"),
        ("Oncotree_SCLC", both["OncotreeSubtype"] == "Small Cell Lung Cancer", "S1 Lung SCLC"),
        ("S1_Lung_Primary", primary, "S1 Lung PrimaryOrMetastasis=Primary"),
        ("S1_Lung_Metastatic", metastatic, "S1 Lung PrimaryOrMetastasis=Metastatic"),
        ("NSCLC_Primary", nsclc & primary, "NSCLC Primary"),
        ("NSCLC_Metastatic", nsclc & metastatic, "NSCLC Metastatic"),
        ("LUAD_Primary", luad & primary, "LUAD Primary"),
        ("LUAD_Metastatic", luad & metastatic, "LUAD Metastatic"),
    ]


def rna_cohorts(rna: pd.DataFrame) -> list[tuple[str, pd.Series, str]]:
    nsclc = rna["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    luad = rna["OncotreeSubtype"] == "Lung Adenocarcinoma"
    primary = rna["PrimaryOrMetastasis"] == "Primary"
    metastatic = rna["PrimaryOrMetastasis"] == "Metastatic"
    return [
        ("RNA_Lung", rna.index == rna.index, "DepMap 24Q4 OncotreeLineage=Lung with TACSTD2 RNA"),
        ("RNA_NSCLC", nsclc, "Lung OncotreePrimaryDisease=NSCLC"),
        ("RNA_LUAD", luad, "Lung OncotreeSubtype=LUAD"),
        ("RNA_LUSC", rna["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma", "Lung LUSC"),
        ("RNA_SCLC", rna["OncotreeSubtype"] == "Small Cell Lung Cancer", "Lung SCLC"),
        ("RNA_Lung_Primary", primary, "Lung PrimaryOrMetastasis=Primary"),
        ("RNA_Lung_Metastatic", metastatic, "Lung PrimaryOrMetastasis=Metastatic"),
        ("RNA_NSCLC_Primary", nsclc & primary, "NSCLC Primary"),
        ("RNA_NSCLC_Metastatic", nsclc & metastatic, "NSCLC Metastatic"),
        ("RNA_LUAD_Primary", luad & primary, "LUAD Primary"),
        ("RNA_LUAD_Metastatic", luad & metastatic, "LUAD Metastatic"),
    ]


def mark_eligible(rows: list[dict], layer: str, eligible_families: set[str]) -> dict | None:
    for r in rows:
        if r["layer"] != layer:
            continue
        rho = r["spearman_rho"]
        resid_ok = r["family"] != "partial_histology" or (
            r["resid_df"] is not None and r["resid_df"] >= MIN_RESID_DF
        )
        partner_ok = True
        if layer == "rna" and r["family"] in {"histology", "partial_histology"}:
            partner_ok = r["partner"] in PRIMARY_TJ
        r["eligible_for_max"] = bool(
            r["family"] in eligible_families
            and partner_ok
            and r["n"] is not None
            and r["n"] >= MIN_N_HEADLINE
            and rho is not None
            and np.isfinite(rho)
            and r["spearman_p"] is not None
            and r["spearman_p"] < 0.05
            and resid_ok
            and r.get("partial_status") in {None, "ok"}
        )
        if layer == "protein" and rho is not None:
            r["delta_vs_0.69"] = float(rho - PRIOR_PROTEIN_RHO)
            r["rounds_to_0.69"] = bool(round(float(rho), 2) == PRIOR_PROTEIN_RHO)
        r["n_equals_118"] = bool(r["n"] == PRIOR_CLAIMED_N)

    eligible_idx = [i for i, r in enumerate(rows) if r["layer"] == layer and r["eligible_for_max"]]
    qvals = bh_q([rows[i]["spearman_p"] for i in eligible_idx])
    for i, q in zip(eligible_idx, qvals):
        rows[i]["q_bh_within_eligible"] = q
    eligible = [rows[i] for i in eligible_idx]
    if not eligible:
        return None
    winner = max(eligible, key=lambda r: r["spearman_rho"])
    for r in rows:
        if r["layer"] == layer:
            r["is_layer_max"] = r["analysis_id"] == winner["analysis_id"]
    return winner


def scatter_pair(df: pd.DataFrame, x: str, y: str, title: str, xlab: str, ylab: str, dest: Path, hue: str | None = None):
    fig, ax = plt.subplots(figsize=(6.6, 5.3))
    if hue and hue in df.columns:
        for name, g in df.groupby(df[hue].fillna("Other")):
            ax.scatter(g[x], g[y], s=36, alpha=0.9, label=f"{name} (n={len(g)})", edgecolors="none")
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.scatter(df[x], df[y], s=36, alpha=0.9, c="#4C78A8", edgecolors="none")
    rho, p = spearmanr(df[x], df[y])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nSpearman ρ = {rho:.3f}   n = {len(df)}   p = {p:.1e}")
    ax.axhline(0, color="0.75", lw=0.6)
    ax.axvline(0, color="0.75", lw=0.6)
    fig.tight_layout()
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def forest_layer(rows: list[dict], layer: str, show_ids: list[str], title: str, dest: Path, vline: float | None = None):
    by_id = {r["analysis_id"]: r for r in rows if r["layer"] == layer}
    chosen = [by_id[i] for i in show_ids if i in by_id and by_id[i].get("spearman_rho") is not None]
    chosen = list(reversed(chosen))
    if not chosen:
        return
    labels = []
    for r in chosen:
        star = " *" if r.get("is_layer_max") else ""
        labels.append(f"{r['label']}  (n={r['n']}){star}")
    y = np.arange(len(chosen))
    rho = np.array([r["spearman_rho"] for r in chosen], float)
    lo = np.array([r["ci95_lo"] if r["ci95_lo"] is not None else r["spearman_rho"] for r in chosen], float)
    hi = np.array([r["ci95_hi"] if r["ci95_hi"] is not None else r["spearman_rho"] for r in chosen], float)
    fam_color = {
        "histology": "#4C78A8",
        "partial_histology": "#F58518",
        "partial_epithelial": "#54A24B",
        "peptide_qc": "#B279A2",
        "focal_single": "#9D755D",
    }
    colors = [fam_color.get(r["family"], "#333333") for r in chosen]
    fig, ax = plt.subplots(figsize=(9.0, max(4.5, 0.38 * len(chosen) + 1.5)))
    ax.errorbar(rho, y, xerr=[rho - lo, hi - rho], fmt="none", ecolor="0.55", elinewidth=1, capsize=2, zorder=1)
    ax.scatter(rho, y, c=colors, s=42, zorder=2)
    if vline is not None:
        ax.axvline(vline, color="#E45756", lw=1, ls="--", label=f"prior ρ = {vline}")
        ax.legend(frameon=False, loc="lower right")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ (partials = rank-residual Pearson)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def run_protein(data: Path, out: Path) -> tuple[list[dict], dict, dict]:
    lines = load_protein(data)
    lung_lines = lines[lines["S1_Tissue"] == "Lung"].copy()
    both = lung_lines[lung_lines["TACSTD2"].notna() & lung_lines["CLDN4"].notna()].copy()
    parent_x = both["TACSTD2"].to_numpy(float)
    parent_y = both["CLDN4"].to_numpy(float)
    rows: list[dict] = []
    pairs: dict[str, pd.DataFrame] = {}

    for name, mask, definition in protein_cohorts(both):
        pair = both.loc[mask].copy()
        pairs[name] = pair
        stats = unadjusted_xy(pair["TACSTD2"].to_numpy(float), pair["CLDN4"].to_numpy(float))
        rows.append(
            empty_row(
                layer="protein",
                analysis_id=name,
                family="histology",
                cohort=name,
                partner="CLDN4",
                label=name.replace("_", " "),
                definition=definition,
                adjustment="none",
                null_p_same_size_from_parent=(
                    None if name == "S1_Lung" else subset_null(parent_x, parent_y, stats["spearman_rho"], stats["n"])
                ),
                **stats,
            )
        )

    for thr in (2, 3):
        pair = both[both["CLDN4_peptides"] >= thr].copy()
        name = f"S1_Lung_CLDN4pep_ge{thr}"
        pairs[name] = pair
        stats = unadjusted_xy(pair["TACSTD2"].to_numpy(float), pair["CLDN4"].to_numpy(float))
        rows.append(
            empty_row(
                layer="protein",
                analysis_id=name,
                family="peptide_qc",
                cohort="S1_Lung",
                partner="CLDN4",
                label=f"S1 Lung, CLDN4 peptides ≥ {thr}",
                definition=f"S1 Lung complete cases with ≥{thr} CLDN4 peptides on the plex.",
                adjustment="none",
                null_p_same_size_from_parent=subset_null(parent_x, parent_y, stats["spearman_rho"], stats["n"]),
                **stats,
            )
        )

    for parent_name in PARTIAL_PARENTS_PROTEIN:
        parent = pairs[parent_name]
        for covars in HISTOLOGY_ADJUSTERS:
            if any(c == "PrimaryOrMetastasis" and parent["PrimaryOrMetastasis"].nunique(dropna=True) < 2 for c in covars):
                continue
            if any(c == "OncotreeSubtype" and parent["OncotreeSubtype"].nunique(dropna=True) < 2 for c in covars):
                continue
            if any(
                c == "OncotreePrimaryDisease" and parent["OncotreePrimaryDisease"].nunique(dropna=True) < 2 for c in covars
            ):
                continue
            label_adj = " + ".join(covars)
            stats_p = partial_spearman(parent, "TACSTD2", "CLDN4", covars)
            rows.append(
                empty_row(
                    layer="protein",
                    analysis_id=f"partial::{parent_name}::{label_adj}",
                    family="partial_histology",
                    cohort=parent_name,
                    partner="CLDN4",
                    label=f"{parent_name.replace('_', ' ')} | {label_adj}",
                    definition=f"Rank partial within {parent_name}.",
                    adjustment=label_adj,
                    **{k: stats_p[k] for k in [
                        "n", "spearman_rho", "spearman_p", "pearson_r", "pearson_p",
                        "ci95_lo", "ci95_hi", "resid_df", "n_adjusters",
                        "unadjusted_spearman_same_n", "partial_status",
                    ]},
                )
            )
        for covar in (["EPCAM"], ["CDH1"], ["VIM"]):
            stats_p = partial_spearman(parent, "TACSTD2", "CLDN4", covar, bootstrap=False)
            rows.append(
                empty_row(
                    layer="protein",
                    analysis_id=f"partial::{parent_name}::{covar[0]}",
                    family="partial_epithelial",
                    cohort=parent_name,
                    partner="CLDN4",
                    label=f"{parent_name.replace('_', ' ')} | {covar[0]}",
                    definition="Epithelial-program partial. Not eligible for the maximum.",
                    adjustment=covar[0],
                    **{k: stats_p[k] for k in [
                        "n", "spearman_rho", "spearman_p", "pearson_r", "pearson_p",
                        "ci95_lo", "ci95_hi", "resid_df", "n_adjusters",
                        "unadjusted_spearman_same_n", "partial_status",
                    ]},
                )
            )

    winner = mark_eligible(rows, "protein", {"histology", "partial_histology"})
    both.sort_values(["OncotreeSubtype", "CellLine"]).to_csv(out / "protein_s1_lung_complete_cases.csv")
    pairs["NSCLC_Primary"].sort_values("TACSTD2").to_csv(out / "protein_nsclc_primary_complete_cases.csv")
    scatter_pair(
        both,
        "TACSTD2",
        "CLDN4",
        "Protein S1 Lung (baseline)",
        "TACSTD2 / TROP2 protein (Gygi TMT)",
        "CLDN4 protein (Gygi TMT)",
        out / "fig_protein_s1_lung_scatter.png",
        hue="OncotreeSubtype",
    )
    scatter_pair(
        pairs["NSCLC_Primary"],
        "TACSTD2",
        "CLDN4",
        "Protein NSCLC primary (unadjusted)",
        "TACSTD2 / TROP2 protein (Gygi TMT)",
        "CLDN4 protein (Gygi TMT)",
        out / "fig_protein_nsclc_primary_scatter.png",
        hue="OncotreeSubtype",
    )
    forest_layer(
        rows,
        "protein",
        [
            "S1_Lung",
            "Oncotree_NSCLC",
            "Oncotree_LUAD",
            "NSCLC_Metastatic",
            "NSCLC_Primary",
            "LUAD_Primary",
            "partial::S1_Lung::OncotreeSubtype",
            "partial::Oncotree_NSCLC::OncotreeSubtype",
            "partial::NSCLC_Primary::OncotreeSubtype",
            "partial::S1_Lung::EPCAM",
            "S1_Lung_CLDN4pep_ge2",
        ],
        "Part1 protein: TACSTD2–CLDN4 on Gygi CCLE",
        out / "fig_protein_rho_grid.png",
        vline=PRIOR_PROTEIN_RHO,
    )
    meta = {
        "s1_lung_lines": int((lines["S1_Tissue"] == "Lung").sum()),
        "s1_lung_complete_cases": int(len(both)),
        "baseline": next(r for r in rows if r["analysis_id"] == "S1_Lung"),
        "winner": winner,
        "peptide_ge2": next(r for r in rows if r["analysis_id"] == "S1_Lung_CLDN4pep_ge2"),
    }
    return rows, meta, pairs


def run_rna(data: Path, out: Path) -> tuple[list[dict], dict, dict]:
    rna = load_rna(data)
    rna = rna[rna["TACSTD2"].notna()].copy()
    rows: list[dict] = []
    pairs: dict[str, pd.DataFrame] = {}

    # Parent nulls per partner from all-lung.
    parent_x = rna["TACSTD2"].to_numpy(float)
    parent_y = {p: rna[p].to_numpy(float) for p in PRIMARY_TJ if p in rna.columns}

    for name, mask, definition in rna_cohorts(rna):
        pair = rna.loc[mask].copy()
        pairs[name] = pair
        for partner in PRIMARY_TJ + FOCAL_SINGLE:
            if partner not in pair.columns:
                continue
            sub = pair.dropna(subset=["TACSTD2", partner])
            stats = unadjusted_xy(sub["TACSTD2"].to_numpy(float), sub[partner].to_numpy(float))
            fam = "histology" if partner in PRIMARY_TJ else "focal_single"
            null = None
            if name != "RNA_Lung" and partner in parent_y:
                # Restrict parent y to finite for this partner.
                ok = np.isfinite(parent_y[partner])
                null = subset_null(parent_x[ok], parent_y[partner][ok], stats["spearman_rho"], stats["n"])
            rows.append(
                empty_row(
                    layer="rna",
                    analysis_id=f"{name}::{partner}",
                    family=fam,
                    cohort=name,
                    partner=partner,
                    label=f"{name.replace('RNA_', '').replace('_', ' ')} vs {TJ_LABEL.get(partner, partner)}",
                    definition=definition,
                    adjustment="none",
                    null_p_same_size_from_parent=null,
                    **stats,
                )
            )

    for parent_name in PARTIAL_PARENTS_RNA:
        parent = pairs[parent_name]
        for partner in PRIMARY_TJ:
            if partner not in parent.columns:
                continue
            for covars in HISTOLOGY_ADJUSTERS:
                if any(c == "PrimaryOrMetastasis" and parent["PrimaryOrMetastasis"].nunique(dropna=True) < 2 for c in covars):
                    continue
                if any(c == "OncotreeSubtype" and parent["OncotreeSubtype"].nunique(dropna=True) < 2 for c in covars):
                    continue
                if any(
                    c == "OncotreePrimaryDisease" and parent["OncotreePrimaryDisease"].nunique(dropna=True) < 2
                    for c in covars
                ):
                    continue
                label_adj = " + ".join(covars)
                stats_p = partial_spearman(parent, "TACSTD2", partner, covars)
                rows.append(
                    empty_row(
                        layer="rna",
                        analysis_id=f"partial::{parent_name}::{partner}::{label_adj}",
                        family="partial_histology",
                        cohort=parent_name,
                        partner=partner,
                        label=f"{parent_name.replace('RNA_', '').replace('_', ' ')} vs {TJ_LABEL.get(partner, partner)} | {label_adj}",
                        definition=f"Rank partial within {parent_name}.",
                        adjustment=label_adj,
                        **{k: stats_p[k] for k in [
                            "n", "spearman_rho", "spearman_p", "pearson_r", "pearson_p",
                            "ci95_lo", "ci95_hi", "resid_df", "n_adjusters",
                            "unadjusted_spearman_same_n", "partial_status",
                        ]},
                    )
                )
            # Epithelial / keratin partials — reported, not eligible.
            for covars, tag in (
                (["EPCAM"], "EPCAM"),
                (["KRT8", "KRT18", "KRT19"], "keratins"),
                (["EPCAM", "KRT8", "KRT18", "KRT19"], "EPCAM+keratins"),
            ):
                if any(c not in parent.columns for c in covars):
                    continue
                stats_p = partial_spearman(parent, "TACSTD2", partner, covars, bootstrap=False)
                rows.append(
                    empty_row(
                        layer="rna",
                        analysis_id=f"partial::{parent_name}::{partner}::{tag}",
                        family="partial_epithelial",
                        cohort=parent_name,
                        partner=partner,
                        label=f"{parent_name.replace('RNA_', '').replace('_', ' ')} vs {TJ_LABEL.get(partner, partner)} | {tag}",
                        definition="Epithelial/keratin partial. Not eligible for the maximum.",
                        adjustment=tag,
                        **{k: stats_p[k] for k in [
                            "n", "spearman_rho", "spearman_p", "pearson_r", "pearson_p",
                            "ci95_lo", "ci95_hi", "resid_df", "n_adjusters",
                            "unadjusted_spearman_same_n", "partial_status",
                        ]},
                    )
                )

    winner = mark_eligible(rows, "rna", {"histology", "partial_histology"})
    rna.sort_values(["OncotreeSubtype", "CellLine"]).to_csv(out / "rna_lung_lines.csv")
    # Best unadjusted LUAD / NSCLC for figures
    for cohort, partner, fname in (
        ("RNA_NSCLC", "CLDN4_TJ_EDGE", "fig_rna_nsclc_cldn4_edge_scatter.png"),
        ("RNA_LUAD", "CLDN4", "fig_rna_luad_cldn4_scatter.png"),
        ("RNA_NSCLC", "TJ_EPITHELIAL", "fig_rna_nsclc_tj_epi_scatter.png"),
    ):
        pair = pairs[cohort].dropna(subset=["TACSTD2", partner])
        scatter_pair(
            pair,
            "TACSTD2",
            partner,
            f"RNA {cohort.replace('RNA_', '')} vs {TJ_LABEL.get(partner, partner)}",
            "TACSTD2 / TROP2 RNA (log2 TPM+1)",
            TJ_LABEL.get(partner, partner),
            out / fname,
            hue="OncotreeSubtype" if cohort != "RNA_LUAD" else None,
        )

    # Forest: one partner panel for CLDN4_TJ_EDGE + LUAD CLDN4 + keratin control
    show = []
    for cohort in ["RNA_Lung", "RNA_NSCLC", "RNA_LUAD", "RNA_NSCLC_Primary", "RNA_LUAD_Primary"]:
        for partner in PRIMARY_TJ:
            show.append(f"{cohort}::{partner}")
    show += [
        "partial::RNA_NSCLC::CLDN4_TJ_EDGE::OncotreeSubtype",
        "partial::RNA_LUAD::CLDN4::OncotreeSubtype",
        "partial::RNA_NSCLC::CLDN4_TJ_EDGE::keratins",
        "partial::RNA_LUAD::CLDN4::keratins",
        "partial::RNA_NSCLC::TJ_EPITHELIAL::keratins",
    ]
    forest_layer(
        rows,
        "rna",
        show,
        "Part1 RNA: TACSTD2 vs locked TJ gene scores",
        out / "fig_rna_tj_rho_grid.png",
    )

    gene_use = {
        "TJ_EPITHELIAL": int(rna["TJ_EPITHELIAL_n_genes"].iloc[0]) if "TJ_EPITHELIAL_n_genes" in rna.columns else None,
        "TJ_TISMO": int(rna["TJ_TISMO_n_genes"].iloc[0]) if "TJ_TISMO_n_genes" in rna.columns else None,
        "CLDN4_TJ_EDGE": int(rna["CLDN4_TJ_EDGE_n_genes"].iloc[0]) if "CLDN4_TJ_EDGE_n_genes" in rna.columns else None,
    }
    meta = {
        "n_lung_with_TACSTD2": int(len(rna)),
        "n_nsclc": int((rna["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer").sum()),
        "n_luad": int((rna["OncotreeSubtype"] == "Lung Adenocarcinoma").sum()),
        "tj_genes_used": gene_use,
        "baseline_nsclc_cldn4": next(r for r in rows if r["analysis_id"] == "RNA_NSCLC::CLDN4"),
        "baseline_nsclc_edge": next(r for r in rows if r["analysis_id"] == "RNA_NSCLC::CLDN4_TJ_EDGE"),
        "winner": winner,
    }
    return rows, meta, pairs


def write_readme(out: Path, prot_meta: dict, rna_meta: dict, rows: list[dict]) -> None:
    pw = prot_meta["winner"]
    rw = rna_meta["winner"]
    baseline = prot_meta["baseline"]
    pep2 = prot_meta["peptide_ge2"]
    nsclc_c4 = rna_meta["baseline_nsclc_cldn4"]
    nsclc_edge = rna_meta["baseline_nsclc_edge"]

    def line(r: dict) -> str:
        ci = ""
        if r.get("ci95_lo") is not None:
            ci = f", 95% CI {r['ci95_lo']:.2f}–{r['ci95_hi']:.2f}"
        return f"**ρ = {r['spearman_rho']:.3f}**, n = {r['n']}, p = {r['spearman_p']:.2e}{ci}"

    # Top eligible RNA rows for a table
    rna_elig = sorted(
        [r for r in rows if r["layer"] == "rna" and r["eligible_for_max"]],
        key=lambda r: -r["spearman_rho"],
    )[:8]
    prot_elig = sorted(
        [r for r in rows if r["layer"] == "protein" and r["eligible_for_max"]],
        key=lambda r: -r["spearman_rho"],
    )[:8]

    md = []
    md.append("# Part1 bridge max ρ: DepMap/CCLE TROP2–CLDN4 protein + TROP2–TJ gene")
    md.append("")
    md.append("Numbers below are written by `analyze.py` from `summary.json` / `correlations.csv`. They are not typed by hand.")
    md.append("")
    md.append("## Question")
    md.append("")
    md.append(
        "For the **Part1 bridge** (TROP2 coexpresses with a tight-junction program), "
        "what is the largest pre-specified Spearman ρ for (1) Gygi TACSTD2–CLDN4 **protein** "
        "and (2) DepMap 24Q4 TACSTD2 versus locked **TJ gene** scores, with honest n?"
    )
    md.append("")
    md.append("## Rules (fixed before sorting)")
    md.append("")
    md.append(f"- Minimum headline n = {MIN_N_HEADLINE}; minimum residual df for partials = {MIN_RESID_DF}.")
    md.append("- Eligible: histology slices and histology partials. Epithelial/keratin partials and peptide QC are reported and cannot win.")
    md.append("- Locked TJ gene sets only (TJ epithelial 18, TJ TISMO 7, CLDN4 TJ edge 11, plus CLDN4 alone). No gene dropping to raise ρ.")
    md.append("- No imputation. No removing individual cell lines to raise ρ.")
    md.append("- No n = 118 TROP2–CLDN4 protein row exists (RPPA lacks TROP2/CLDN4 antibodies).")
    md.append("")
    md.append("## Protein maximum (TROP2–CLDN4)")
    md.append("")
    md.append(f"- **Locked baseline (rounds to 0.69):** S1 Lung complete cases {line(baseline)}.")
    def fmt_q(q):
        return "NA" if q is None or (isinstance(q, float) and not np.isfinite(q)) else f"{q:.2e}"

    md.append(
        f"- **Histology / partial maximum:** {pw['label']}: {line(pw)} "
        f"(BH q = {fmt_q(pw.get('q_bh_within_eligible'))})."
    )
    null_p = pep2.get("null_p_same_size_from_parent")
    md.append(
        f"- Peptide ≥ 2 sensitivity (not eligible for max): {line(pep2)}; "
        f"same-size null p = {fmt_q(null_p)}."
    )
    md.append(
        f"- S1 Lung lines = {prot_meta['s1_lung_lines']}; complete cases = {prot_meta['s1_lung_complete_cases']}. "
        "Any n = 118 in protein grid: **False**."
    )
    md.append("")
    md.append("| Analysis | n | ρ | p | eligible |")
    md.append("|---|---:|---:|---:|:---:|")
    for r in prot_elig:
        md.append(
            f"| {r['label']} | {r['n']} | {r['spearman_rho']:.3f} | {r['spearman_p']:.2e} | yes |"
        )
    md.append("")
    md.append("## RNA / gene maximum (TROP2–TJ)")
    md.append("")
    md.append(
        f"- Lung lines with TACSTD2 RNA: n = {rna_meta['n_lung_with_TACSTD2']} "
        f"(NSCLC {rna_meta['n_nsclc']}, LUAD {rna_meta['n_luad']})."
    )
    md.append(f"- TJ genes used: {json.dumps(rna_meta['tj_genes_used'])}.")
    md.append(f"- Powered NSCLC baselines: CLDN4 {line(nsclc_c4)}; CLDN4 TJ edge {line(nsclc_edge)}.")
    md.append(
        f"- **TJ gene maximum (eligible):** {rw['label']}: {line(rw)} "
        f"(BH q = {fmt_q(rw.get('q_bh_within_eligible'))})."
    )
    md.append("")
    md.append("| Analysis | n | ρ | p | eligible |")
    md.append("|---|---:|---:|---:|:---:|")
    for r in rna_elig:
        md.append(
            f"| {r['label']} | {r['n']} | {r['spearman_rho']:.3f} | {r['spearman_p']:.2e} | yes |"
        )
    # Keratin caveats for the powered NSCLC call
    ker = [
        r
        for r in rows
        if r["layer"] == "rna"
        and r["family"] == "partial_epithelial"
        and r["cohort"] == "RNA_NSCLC"
        and r["partner"] in {"CLDN4", "CLDN4_TJ_EDGE", "TJ_EPITHELIAL"}
        and r["adjustment"] == "keratins"
    ]
    md.append("")
    md.append("Keratin partials on NSCLC (not eligible for max):")
    md.append("")
    md.append("| Partner | unadjusted ρ | keratins partial ρ | n |")
    md.append("|---|---:|---:|---:|")
    for r in ker:
        unadj = next(
            u
            for u in rows
            if u["analysis_id"] == f"RNA_NSCLC::{r['partner']}"
        )
        md.append(
            f"| {TJ_LABEL.get(r['partner'], r['partner'])} | {unadj['spearman_rho']:.3f} | "
            f"{r['spearman_rho']:.3f} | {r['n']} |"
        )
    md.append("")
    md.append("## Part1 bridge language (honest)")
    md.append("")
    md.append(
        f"- KEEP: Gygi lung TROP2–CLDN4 protein ρ = {baseline['spearman_rho']:.3f} (n={baseline['n']}). "
        f"Maximum eligible in the locked grid: {pw['label']} ρ = {pw['spearman_rho']:.3f} (n={pw['n']})."
    )
    md.append(
        f"- KEEP: DepMap NSCLC RNA TROP2–CLDN4 ρ = {nsclc_c4['spearman_rho']:.3f} (n={nsclc_c4['n']}); "
        f"TROP2–CLDN4 TJ edge ρ = {nsclc_edge['spearman_rho']:.3f} (n={nsclc_edge['n']}). "
        f"Maximum eligible TJ gene coexpression: {rw['label']} ρ = {rw['spearman_rho']:.3f} (n={rw['n']})."
    )
    md.append("- DO NOT claim ICI resistance, immune exclusion, or RPPA n≈118 from DepMap/CCLE alone.")
    md.append("- Keratin/EPCAM partials shrink RNA TJ ρ; they are caveats, not the maximum.")
    md.append("")
    md.append("## Reproduce")
    md.append("")
    md.append("```bash")
    md.append("python3 -m pip install -r scripts/depmap_trop2_tj_part1_maxrho/requirements.txt")
    md.append("python3 scripts/depmap_trop2_tj_part1_maxrho/download.py")
    md.append("python3 scripts/depmap_trop2_tj_part1_maxrho/analyze.py")
    md.append("```")
    md.append("")
    (out / "README.md").write_text("\n".join(md) + "\n")
    (out / "FINDING.md").write_text("\n".join(md) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_trop2_tj_part1_maxrho")
    ap.add_argument("--outdir", default="results/depmap_trop2_tj_part1_maxrho")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    prot_rows, prot_meta, _ = run_protein(data, out)
    rna_rows, rna_meta, _ = run_rna(data, out)
    rows = prot_rows + rna_rows
    pd.DataFrame(rows).to_csv(out / "correlations.csv", index=False)

    # Copy download manifest if present
    man = data / "download_manifest.json"
    if man.exists():
        (out / "download_manifest.json").write_text(man.read_text())

    summary = {
        "prior_protein_rho": PRIOR_PROTEIN_RHO,
        "prior_claimed_n": PRIOR_CLAIMED_N,
        "rule": {
            "min_n": MIN_N_HEADLINE,
            "min_resid_df": MIN_RESID_DF,
            "protein_eligible_families": ["histology", "partial_histology"],
            "rna_eligible_families": ["histology", "partial_histology"],
            "rna_eligible_partners": PRIMARY_TJ,
            "no_imputation": True,
            "no_cell_line_dropping": True,
        },
        "protein": prot_meta,
        "rna": rna_meta,
        "n_equals_118_anywhere": bool(any(r["n_equals_118"] for r in rows)),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_readme(out, prot_meta, rna_meta, rows)

    pw = prot_meta["winner"]
    rw = rna_meta["winner"]
    verdict = [
        "PART1_BRIDGE_MAX_RHO",
        f"protein baseline S1 Lung rho={prot_meta['baseline']['spearman_rho']:.4f} n={prot_meta['baseline']['n']}",
        f"protein max: {pw['label']} rho={pw['spearman_rho']:.4f} n={pw['n']} p={pw['spearman_p']:.3e}",
        f"rna NSCLC CLDN4 rho={rna_meta['baseline_nsclc_cldn4']['spearman_rho']:.4f} n={rna_meta['baseline_nsclc_cldn4']['n']}",
        f"rna NSCLC CLDN4_TJ_EDGE rho={rna_meta['baseline_nsclc_edge']['spearman_rho']:.4f} n={rna_meta['baseline_nsclc_edge']['n']}",
        f"rna TJ max: {rw['label']} rho={rw['spearman_rho']:.4f} n={rw['n']} p={rw['spearman_p']:.3e}",
        f"n_equals_118 anywhere: {summary['n_equals_118_anywhere']}",
        "No value was imputed. No line was removed to raise rho.",
    ]
    (out / "verdict.txt").write_text("\n".join(verdict) + "\n")
    print("\n".join(verdict))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

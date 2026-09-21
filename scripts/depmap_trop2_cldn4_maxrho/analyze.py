#!/usr/bin/env python3
"""Maximum TACSTD2–CLDN4 protein correlation on the Gygi CCLE table.

The cohort list and the partial-correlation adjusters are fixed below.
The maximum is the largest estimate in that list. It is not a search over
individual cell lines, and it is not forced toward n=118.

Prior result, recomputed here as the baseline: Table S1 Lung, pairwise
complete TACSTD2 and CLDN4, Spearman ρ ≈ 0.693 on n=45.
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

USER_RHO = 0.69
USER_N = 118
MIN_N_HEADLINE = 15
MIN_RESID_DF = 10
N_BOOT = 5000
N_NULL = 5000
SEED = 0

SAMPLE_RE = re.compile(r"_TenPx\d+$")
GENES = ["TACSTD2", "CLDN4", "EPCAM", "VIM", "CDH1"]
# Histology / disease / specimen-site adjusters. Epithelial proteins are a
# separate family and cannot win the maximum.
HISTOLOGY_ADJUSTERS = [
    ["OncotreeSubtype"],
    ["OncotreePrimaryDisease"],
    ["PrimaryOrMetastasis"],
    ["OncotreeSubtype", "PrimaryOrMetastasis"],
]
EPITHELIAL_ADJUSTERS = [["EPCAM"], ["CDH1"], ["VIM"]]
PARTIAL_PARENTS = [
    "S1_Lung",
    "Oncotree_NSCLC",
    "Oncotree_LUAD",
    "NSCLC_Primary",
    "S1_Lung_Primary",
]


def core_of(col: str) -> str:
    m = re.match(r"(.+)_TenPx\d+$", col)
    return m.group(1) if m else col


def plex_of(col: str) -> str:
    m = re.search(r"(TenPx\d+)$", col)
    return m.group(1) if m else ""


def fast_spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(x)
    ry = rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
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
    """Fraction of random same-size subsets of the parent whose ρ is at least obs."""
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


def load_lines(data: Path) -> pd.DataFrame:
    prot = pd.read_csv(data / "protein_quant_current_normalized.csv.gz")
    missing = [g for g in GENES if g not in set(prot["Gene_Symbol"])]
    if missing:
        raise SystemExit(f"Missing genes: {missing}")
    sub = prot[prot["Gene_Symbol"].isin(GENES)].drop_duplicates("Gene_Symbol").set_index("Gene_Symbol")
    sample_cols = [c for c in prot.columns if SAMPLE_RE.search(str(c))]
    pep_cols = [c for c in prot.columns if str(c).endswith("_Peptides")]
    cl_pep = {
        c.replace("_Peptides", ""): float(pd.to_numeric(sub.loc["CLDN4", c], errors="coerce"))
        for c in pep_cols
    }
    cols_by: dict[str, list[str]] = defaultdict(list)
    for c in sample_cols:
        cols_by[core_of(c)].append(c)

    rows = []
    for core, cols in cols_by.items():
        rec = {"CCLE": core, "n_replicates": len(cols), "plex": plex_of(cols[0])}
        for g in GENES:
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


def unadjusted(pair: pd.DataFrame) -> dict:
    x = pair["TACSTD2"].to_numpy(float)
    y = pair["CLDN4"].to_numpy(float)
    n = int(len(pair))
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


def _design(d: pd.DataFrame, covars: list[str]) -> tuple[pd.DataFrame, np.ndarray] | None:
    numeric = [c for c in covars if c in GENES or c == "CLDN4_peptides"]
    d = d.dropna(subset=["TACSTD2", "CLDN4"] + numeric).copy()
    pieces = []
    for c in covars:
        if c in numeric:
            pieces.append(d[[c]].to_numpy(float))
            continue
        s = d[c].fillna("NA").astype(str)
        vc = s.value_counts()
        # Levels with fewer than 3 lines are pooled so a single rare code
        # cannot become its own dummy.
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


def partial_spearman(d: pd.DataFrame, covars: list[str], bootstrap: bool = True) -> dict:
    built = _design(d, covars)
    base = {
        "n": int(len(d.dropna(subset=["TACSTD2", "CLDN4"]))),
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "ci95_lo": None,
        "ci95_hi": None,
        "resid_df": None,
        "n_adjusters": None,
        "partial_status": "not_estimable",
    }
    if built is None:
        return base
    d, Z = built
    y = d["TACSTD2"].to_numpy(float)
    x = d["CLDN4"].to_numpy(float)
    n, k = len(d), int(Z.shape[1])
    df_ = n - k - 2
    base.update({"n": n, "n_adjusters": k, "resid_df": int(df_)})
    if df_ < 5:
        base["partial_status"] = "residual_df_below_5"
        return base

    def _r(y_in, x_in, Z_in) -> float:
        def rank(a):
            return rankdata(a).astype(float)

        yr, xr = rank(y_in), rank(x_in)
        Zr = np.column_stack([rank(Z_in[:, i]) for i in range(Z_in.shape[1])])
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
            # Keep the original design codes, re-estimated inside _r via ranks.
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
            "unadjusted_spearman_p_same_n": float(p_u),
            "ci95_lo": lo,
            "ci95_hi": hi,
            "partial_status": "ok",
        }
    )
    return base


def coarse_class(row: pd.Series) -> str:
    if row["OncotreeSubtype"] == "Lung Adenocarcinoma":
        return "LUAD"
    if row["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma":
        return "LUSC"
    if row["OncotreeSubtype"] == "Small Cell Lung Cancer":
        return "SCLC"
    if row["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer":
        return "Other NSCLC"
    return "Other lung-labeled"


def scatter_baseline(lung: pd.DataFrame, dest: Path) -> None:
    colors = {
        "LUAD": "#F58518",
        "LUSC": "#E45756",
        "SCLC": "#54A24B",
        "Other NSCLC": "#4C78A8",
        "Other lung-labeled": "#9D755D",
    }
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    plot = lung.copy()
    plot["cls"] = plot.apply(coarse_class, axis=1)
    for name, g in plot.groupby("cls"):
        ax.scatter(
            g["TACSTD2"],
            g["CLDN4"],
            s=36,
            alpha=0.9,
            c=colors.get(name, "#333333"),
            label=f"{name} (n={len(g)})",
            edgecolors="none",
        )
    rho, p = spearmanr(plot["TACSTD2"], plot["CLDN4"])
    ax.set_xlabel("TACSTD2 / TROP2 protein (Gygi TMT, normalized)")
    ax.set_ylabel("CLDN4 protein (Gygi TMT, normalized)")
    ax.set_title(f"S1 Lung, both proteins\nSpearman ρ = {rho:.3f}   n = {len(plot)}   p = {p:.1e}")
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0, color="0.75", lw=0.6)
    ax.axvline(0, color="0.75", lw=0.6)
    fig.tight_layout()
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def scatter_nsclc_primary(pair: pd.DataFrame, rho: float, p: float, dest: Path) -> None:
    colors = {
        "Lung Adenocarcinoma": "#F58518",
        "Lung Squamous Cell Carcinoma": "#E45756",
        "Non-Small Cell Lung Cancer": "#4C78A8",
    }
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    for name, g in pair.groupby(pair["OncotreeSubtype"].fillna("Other")):
        ax.scatter(
            g["TACSTD2"],
            g["CLDN4"],
            s=42,
            alpha=0.9,
            c=colors.get(name, "#9D755D"),
            label=f"{name} (n={len(g)})",
            edgecolors="none",
        )
        for _, row in g.iterrows():
            ax.annotate(
                str(row["CellLine"]),
                (row["TACSTD2"], row["CLDN4"]),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=7,
                color="0.25",
            )
    ax.set_xlabel("TACSTD2 / TROP2 protein (Gygi TMT, normalized)")
    ax.set_ylabel("CLDN4 protein (Gygi TMT, normalized)")
    ax.set_title(
        "NSCLC lines from primary tumors, unadjusted\n"
        f"Spearman ρ = {rho:.3f}   n = {len(pair)}   p = {p:.1e}"
    )
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0, color="0.75", lw=0.6)
    ax.axvline(0, color="0.75", lw=0.6)
    fig.tight_layout()
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def forest(rows: list[dict], dest: Path) -> None:
    show_ids = [
        "S1_Lung",
        "Oncotree_NSCLC",
        "Oncotree_LUAD",
        "Oncotree_LUSC",
        "Oncotree_SCLC",
        "NSCLC_Metastatic",
        "NSCLC_Primary",
        "LUAD_Primary",
        "partial::S1_Lung::OncotreeSubtype",
        "partial::Oncotree_NSCLC::OncotreeSubtype",
        "partial::NSCLC_Primary::OncotreeSubtype",
        "partial::S1_Lung::EPCAM",
        "partial::NSCLC_Primary::EPCAM",
        "S1_Lung_CLDN4pep_ge2",
        "S1_Lung_CLDN4pep_ge3",
    ]
    by_id = {r["analysis_id"]: r for r in rows}
    chosen = [by_id[i] for i in show_ids if i in by_id and r_ok(by_id[i])]
    chosen = list(reversed(chosen))
    labels = []
    for r in chosen:
        star = "  *" if r.get("is_histology_max") else ""
        labels.append(f"{r['label']}  (n={r['n']}){star}")
    y = np.arange(len(chosen))
    rho = np.array([r["spearman_rho"] for r in chosen], float)
    lo = np.array([r["ci95_lo"] if r["ci95_lo"] is not None else r["spearman_rho"] for r in chosen], float)
    hi = np.array([r["ci95_hi"] if r["ci95_hi"] is not None else r["spearman_rho"] for r in chosen], float)
    family_color = {
        "histology": "#4C78A8",
        "partial_histology": "#F58518",
        "partial_epithelial": "#54A24B",
        "peptide_qc": "#B279A2",
    }
    colors = [family_color.get(r["family"], "#333333") for r in chosen]
    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    ax.errorbar(
        rho,
        y,
        xerr=[rho - lo, hi - rho],
        fmt="none",
        ecolor="0.55",
        elinewidth=1,
        capsize=2,
        zorder=1,
    )
    ax.scatter(rho, y, c=colors, s=42, zorder=2)
    ax.axvline(USER_RHO, color="#E45756", lw=1, ls="--", label="Prior lung ρ = 0.69")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ  (partial rows are rank-residual Pearson)")
    ax.set_xlim(0.2, 1.02)
    ax.set_title("TACSTD2–CLDN4 protein on the Gygi CCLE table")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def r_ok(r: dict) -> bool:
    return r.get("spearman_rho") is not None and np.isfinite(r["spearman_rho"])


def empty_row(**kwargs) -> dict:
    base = {
        "analysis_id": None,
        "family": None,
        "cohort": None,
        "label": None,
        "definition": None,
        "adjustment": None,
        "n": None,
        "n_lines_in_filter": None,
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "ci95_lo": None,
        "ci95_hi": None,
        "resid_df": None,
        "n_adjusters": None,
        "unadjusted_spearman_same_n": None,
        "delta_vs_0.69": None,
        "rounds_to_0.69": False,
        "n_equals_118": False,
        "eligible_for_histology_max": False,
        "null_p_same_size_from_S1_lung": None,
        "partial_status": None,
    }
    base.update(kwargs)
    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_trop2_cldn4_maxrho")
    ap.add_argument("--outdir", default="results/depmap_trop2_cldn4_maxrho")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    lines = load_lines(data)
    lung_lines = lines[lines["S1_Tissue"] == "Lung"].copy()
    both = lung_lines[lung_lines["TACSTD2"].notna() & lung_lines["CLDN4"].notna()].copy()
    parent_x = both["TACSTD2"].to_numpy(float)
    parent_y = both["CLDN4"].to_numpy(float)

    nsclc = both["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    luad = both["OncotreeSubtype"] == "Lung Adenocarcinoma"
    primary = both["PrimaryOrMetastasis"] == "Primary"
    metastatic = both["PrimaryOrMetastasis"] == "Metastatic"

    cohorts = [
        ("S1_Lung", both.index == both.index, "Table S1 Tissue of Origin = Lung, TACSTD2 and CLDN4 both quantified"),
        ("Oncotree_Lung", both["OncotreeLineage"] == "Lung", "S1 Lung and DepMap 24Q4 OncotreeLineage = Lung"),
        ("Oncotree_NSCLC", nsclc, "S1 Lung and OncotreePrimaryDisease = Non-Small Cell Lung Cancer"),
        ("Oncotree_LUAD", luad, "S1 Lung and OncotreeSubtype = Lung Adenocarcinoma"),
        ("Oncotree_LUSC", both["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma", "S1 Lung and OncotreeSubtype = Lung Squamous Cell Carcinoma"),
        ("Oncotree_SCLC", both["OncotreeSubtype"] == "Small Cell Lung Cancer", "S1 Lung and OncotreeSubtype = Small Cell Lung Cancer"),
        ("NSCLC_not_LUSC", nsclc & (both["OncotreeSubtype"] != "Lung Squamous Cell Carcinoma"), "NSCLC excluding lung squamous cell carcinoma"),
        ("NSCLC_not_LUAD", nsclc & ~luad, "NSCLC excluding lung adenocarcinoma"),
        ("S1_Lung_Primary", primary, "S1 Lung and PrimaryOrMetastasis = Primary"),
        ("S1_Lung_Metastatic", metastatic, "S1 Lung and PrimaryOrMetastasis = Metastatic"),
        ("NSCLC_Primary", nsclc & primary, "NSCLC and PrimaryOrMetastasis = Primary"),
        ("NSCLC_Metastatic", nsclc & metastatic, "NSCLC and PrimaryOrMetastasis = Metastatic"),
        ("LUAD_Primary", luad & primary, "LUAD and PrimaryOrMetastasis = Primary"),
        ("LUAD_Metastatic", luad & metastatic, "LUAD and PrimaryOrMetastasis = Metastatic"),
    ]

    rows: list[dict] = []
    pairs: dict[str, pd.DataFrame] = {}
    for name, mask, definition in cohorts:
        pair = both.loc[mask].copy()
        pairs[name] = pair
        stats = unadjusted(pair)
        rows.append(
            empty_row(
                analysis_id=name,
                family="histology",
                cohort=name,
                label=name.replace("_", " "),
                definition=definition,
                adjustment="none",
                n_lines_in_filter=int(mask.sum()) if not isinstance(mask, bool) else int(len(pair)),
                null_p_same_size_from_S1_lung=(
                    None
                    if name == "S1_Lung"
                    else subset_null(parent_x, parent_y, stats["spearman_rho"], stats["n"])
                ),
                **stats,
            )
        )

    # Peptide-depth sensitivities. Thresholds are the observed integer peptide
    # counts above the minimum of 1 (every quantified CLDN4 value in this
    # lung set already has at least one peptide).
    for thr in (2, 3):
        pair = both[both["CLDN4_peptides"] >= thr].copy()
        name = f"S1_Lung_CLDN4pep_ge{thr}"
        pairs[name] = pair
        stats = unadjusted(pair)
        rows.append(
            empty_row(
                analysis_id=name,
                family="peptide_qc",
                cohort="S1_Lung",
                label=f"S1 Lung, CLDN4 peptides ≥ {thr}",
                definition=(
                    "S1 Lung complete cases whose TMT plex identified CLDN4 with "
                    f"at least {thr} peptides. Peptide count is a plex-level depth, "
                    "not a per-line abundance filter."
                ),
                adjustment="none",
                n_lines_in_filter=int(len(pair)),
                null_p_same_size_from_S1_lung=subset_null(
                    parent_x, parent_y, stats["spearman_rho"], stats["n"]
                ),
                **stats,
            )
        )
        # Histology partial inside the ≥2-peptide set only (pre-set, one row).
        if thr == 2:
            stats_p = partial_spearman(pair, ["OncotreeSubtype"])
            rows.append(
                empty_row(
                    analysis_id="partial::S1_Lung_CLDN4pep_ge2::OncotreeSubtype",
                    family="peptide_qc",
                    cohort="S1_Lung_CLDN4pep_ge2",
                    label="S1 Lung, CLDN4 peptides ≥ 2 | subtype",
                    definition="Peptide ≥ 2 set, Spearman partial adjusted for Oncotree subtype.",
                    adjustment="OncotreeSubtype",
                    n_lines_in_filter=int(len(pair)),
                    partial_status=stats_p.get("partial_status"),
                    unadjusted_spearman_same_n=stats_p.get("unadjusted_spearman_same_n"),
                    **{k: stats_p[k] for k in ["n", "spearman_rho", "spearman_p", "pearson_r", "pearson_p", "ci95_lo", "ci95_hi", "resid_df", "n_adjusters"]},
                )
            )

    for parent_name in PARTIAL_PARENTS:
        parent = pairs[parent_name]
        for covars, family in (
            (c, "partial_histology") for c in HISTOLOGY_ADJUSTERS
        ):
            # Skip adjusters that cannot vary inside this cohort.
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
            stats_p = partial_spearman(parent, covars)
            rows.append(
                empty_row(
                    analysis_id=f"partial::{parent_name}::{label_adj}",
                    family=family,
                    cohort=parent_name,
                    label=f"{parent_name.replace('_', ' ')} | {label_adj}",
                    definition=f"Spearman partial (rank, then residual Pearson) within {parent_name}.",
                    adjustment=label_adj,
                    n_lines_in_filter=int(len(parent)),
                    partial_status=stats_p.get("partial_status"),
                    unadjusted_spearman_same_n=stats_p.get("unadjusted_spearman_same_n"),
                    **{
                        k: stats_p[k]
                        for k in [
                            "n",
                            "spearman_rho",
                            "spearman_p",
                            "pearson_r",
                            "pearson_p",
                            "ci95_lo",
                            "ci95_hi",
                            "resid_df",
                            "n_adjusters",
                        ]
                    },
                )
            )
        for covars in EPITHELIAL_ADJUSTERS:
            label_adj = covars[0]
            stats_p = partial_spearman(parent, covars, bootstrap=False)
            rows.append(
                empty_row(
                    analysis_id=f"partial::{parent_name}::{label_adj}",
                    family="partial_epithelial",
                    cohort=parent_name,
                    label=f"{parent_name.replace('_', ' ')} | {label_adj}",
                    definition=(
                        "Epithelial-program partial. Reported because it shrinks ρ. "
                        "Not eligible for the maximum."
                    ),
                    adjustment=label_adj,
                    n_lines_in_filter=int(len(parent)),
                    partial_status=stats_p.get("partial_status"),
                    unadjusted_spearman_same_n=stats_p.get("unadjusted_spearman_same_n"),
                    **{
                        k: stats_p[k]
                        for k in [
                            "n",
                            "spearman_rho",
                            "spearman_p",
                            "pearson_r",
                            "pearson_p",
                            "ci95_lo",
                            "ci95_hi",
                            "resid_df",
                            "n_adjusters",
                        ]
                    },
                )
            )

    for r in rows:
        rho = r["spearman_rho"]
        r["delta_vs_0.69"] = None if rho is None else float(rho - USER_RHO)
        r["rounds_to_0.69"] = bool(rho is not None and round(float(rho), 2) == USER_RHO)
        r["n_equals_118"] = bool(r["n"] == USER_N)
        resid_ok = r["family"] != "partial_histology" or (
            r["resid_df"] is not None and r["resid_df"] >= MIN_RESID_DF
        )
        r["eligible_for_histology_max"] = bool(
            r["family"] in {"histology", "partial_histology"}
            and r["n"] is not None
            and r["n"] >= MIN_N_HEADLINE
            and rho is not None
            and np.isfinite(rho)
            and r["spearman_p"] is not None
            and r["spearman_p"] < 0.05
            and resid_ok
            and r.get("partial_status") in {None, "ok"}
        )

    eligible_idx = [i for i, r in enumerate(rows) if r["eligible_for_histology_max"]]
    qvals = bh_q([rows[i]["spearman_p"] for i in eligible_idx])
    for i, q in zip(eligible_idx, qvals):
        rows[i]["q_bh_within_eligible"] = q
    for r in rows:
        r.setdefault("q_bh_within_eligible", None)

    eligible = [r for r in rows if r["eligible_for_histology_max"]]
    winner = max(eligible, key=lambda r: r["spearman_rho"])
    for r in rows:
        r["is_histology_max"] = r["analysis_id"] == winner["analysis_id"]

    # Leave-one-out for the winning unadjusted cohort sitting under a partial,
    # and for NSCLC primary itself.
    loo_rows = []
    for cohort_name in ["NSCLC_Primary", "Oncotree_LUAD", "S1_Lung_CLDN4pep_ge2"]:
        pair = pairs[cohort_name]
        for cc in pair.index:
            held = pair.drop(index=cc)
            rho, p = spearmanr(held["TACSTD2"], held["CLDN4"])
            loo_rows.append(
                {
                    "cohort": cohort_name,
                    "held_out_ccle": cc,
                    "held_out_cell_line": pair.loc[cc, "CellLine"],
                    "n": int(len(held)),
                    "spearman_rho": float(rho),
                    "spearman_p": float(p),
                }
            )

    corr = pd.DataFrame(rows)
    corr.to_csv(out / "correlations.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(out / "leave_one_out.csv", index=False)

    export_cols = [
        "CellLine",
        "ModelID",
        "S1_Tissue",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
        "PrimaryOrMetastasis",
        "Sex",
        "plex",
        "CLDN4_peptides",
        "TACSTD2",
        "CLDN4",
        "EPCAM",
        "CDH1",
        "VIM",
    ]
    both.sort_values(["OncotreeSubtype", "CellLine"]).to_csv(out / "s1_lung_complete_cases.csv", columns=export_cols)
    pairs["NSCLC_Primary"].sort_values("TACSTD2").to_csv(
        out / "nsclc_primary_complete_cases.csv", columns=export_cols
    )

    scatter_baseline(both, out / "fig_scatter_s1_lung.png")
    npair = pairs["NSCLC_Primary"]
    nstats = next(r for r in rows if r["analysis_id"] == "NSCLC_Primary")
    scatter_nsclc_primary(npair, nstats["spearman_rho"], nstats["spearman_p"], out / "fig_scatter_nsclc_primary.png")
    forest(rows, out / "fig_rho_grid.png")

    any_118 = bool(corr["n_equals_118"].any())
    baseline = next(r for r in rows if r["analysis_id"] == "S1_Lung")
    pep2 = next(r for r in rows if r["analysis_id"] == "S1_Lung_CLDN4pep_ge2")
    pep_vs_tac, pep_vs_tac_p = spearmanr(both["CLDN4_peptides"], both["TACSTD2"])
    pep_vs_cld, pep_vs_cld_p = spearmanr(both["CLDN4_peptides"], both["CLDN4"])
    peptide_vs_abundance = {
        "n": int(len(both)),
        "spearman_CLDN4_peptides_vs_TACSTD2": float(pep_vs_tac),
        "p_vs_TACSTD2": float(pep_vs_tac_p),
        "spearman_CLDN4_peptides_vs_CLDN4": float(pep_vs_cld),
        "p_vs_CLDN4": float(pep_vs_cld_p),
        "note": "Plex-level CLDN4 peptide count versus abundance in the 45 S1 lung complete cases.",
    }
    pd.DataFrame([peptide_vs_abundance]).to_csv(out / "cldn4_peptide_vs_abundance.csv", index=False)
    summary = {
        "prior_rho": USER_RHO,
        "prior_claimed_n": USER_N,
        "n_118_found": any_118,
        "s1_lung_lines": int((lines["S1_Tissue"] == "Lung").sum()),
        "s1_lung_complete_cases": int(len(both)),
        "baseline": baseline,
        "histology_max": winner,
        "peptide_ge2": pep2,
        "peptide_vs_abundance": peptide_vs_abundance,
        "rule": {
            "min_n": MIN_N_HEADLINE,
            "min_resid_df": MIN_RESID_DF,
            "families_eligible": ["histology", "partial_histology"],
            "families_reported_not_eligible": ["partial_epithelial", "peptide_qc"],
            "no_cell_line_dropping": True,
            "no_imputation": True,
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    loo_np = [r for r in loo_rows if r["cohort"] == "NSCLC_Primary"]
    loo_min = min(loo_np, key=lambda r: r["spearman_rho"])
    loo_max = max(loo_np, key=lambda r: r["spearman_rho"])
    verdict = [
        "MAX_HISTOLOGY_PARTIAL",
        f"baseline S1 Lung Spearman {baseline['spearman_rho']:.4f} n={baseline['n']} (rounds to 0.69)",
        f"histology max: {winner['label']} rho={winner['spearman_rho']:.4f} n={winner['n']} p={winner['spearman_p']:.3e}",
        f"n_equals_118 anywhere in the grid: {any_118}",
        f"NSCLC primary leave-one-out rho range: {loo_min['spearman_rho']:.4f} to {loo_max['spearman_rho']:.4f}",
        f"CLDN4 peptides>=2: rho={pep2['spearman_rho']:.4f} n={pep2['n']} same-size null p={pep2['null_p_same_size_from_S1_lung']:.4f}",
        "No value was imputed. No line was removed to raise rho.",
    ]
    (out / "verdict.txt").write_text("\n".join(verdict) + "\n")
    print("\n".join(verdict))
    show = corr.loc[
        corr["spearman_rho"].notna(),
        ["analysis_id", "family", "n", "spearman_rho", "spearman_p", "ci95_lo", "ci95_hi", "null_p_same_size_from_S1_lung", "eligible_for_histology_max"],
    ]
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

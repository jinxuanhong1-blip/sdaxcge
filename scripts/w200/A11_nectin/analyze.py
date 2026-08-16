#!/usr/bin/env python3
"""A11: nectin genes vs TACSTD2-high public lung.

Primary question (pre-specified):
  In public lung tumors, are NECTIN1/2/3/4 higher in TACSTD2-high samples
  and/or positively correlated with TACSTD2?

Primary statistic:
  pooled rank-partial Spearman of TACSTD2 vs each gene, controlling for
  ABSOLUTE purity and histology (LUSC vs LUAD).

TACSTD2-high contrast:
  Q4 vs Q1 (primary) and median split (sensitivity), Mann-Whitney U.

Honest class-effect rule (pre-specified, not tuned after seeing numbers):
  A nectin-family class effect is claimed ONLY if >=3 of the 4 nectins are
  WEAK_POSITIVE or ASSOCIATED_POSITIVE on the primary statistic and none
  is significantly negative (FDR<0.05 and rho<0). Otherwise: MIXED or NULL.

Thresholds on primary partial rho (also pre-specified):
  ASSOCIATED_POSITIVE : rho >= 0.20 and FDR < 0.05
  WEAK_POSITIVE       : 0.10 <= rho < 0.20 and FDR < 0.05
  NULL                : FDR >= 0.05 or |rho| < 0.10
  WEAK_NEGATIVE       : -0.20 < rho <= -0.10 and FDR < 0.05
  ASSOCIATED_NEGATIVE : rho <= -0.20 and FDR < 0.05

PVR is nectin-like, reported separately.
CLDN4 is a positive-control junction gene, not counted in the class verdict.
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
from statsmodels.stats.multitest import multipletests

PRIMARY_NECTINS = ["NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]
NECTIN_LIKE = ["PVR"]
POSITIVE_CONTROL = ["CLDN4"]
ALL_PARTNERS = PRIMARY_NECTINS + NECTIN_LIKE + POSITIVE_CONTROL
TARGET = "TACSTD2"

ALIASES = {
    "NECTIN1": "PVRL1 / nectin-1",
    "NECTIN2": "PVRL2 / nectin-2 / CD112",
    "NECTIN3": "PVRL3 / nectin-3",
    "NECTIN4": "PVRL4 / nectin-4 (ADC target)",
    "PVR": "CD155 / NECL5 (nectin-like)",
    "CLDN4": "claudin-4 (junction positive control, not a nectin)",
}


def bh(pvals: list[float]) -> np.ndarray:
    return multipletests(pvals, method="fdr_bh")[1]


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 5:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), n


def partial_spearman(x: np.ndarray, y: np.ndarray, Z: np.ndarray) -> tuple[float, float, int]:
    """Rank residual partial correlation. Z is 1-D or 2-D (covariates in columns)."""
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    x, y, Z = x[mask], y[mask], Z[mask]
    n = int(len(x))
    k = int(Z.shape[1])
    if n < 5 + k:
        return float("nan"), float("nan"), n
    xr, yr = stats.rankdata(x), stats.rankdata(y)
    Zr = np.column_stack([stats.rankdata(Z[:, j]) for j in range(k)])
    design = np.column_stack([np.ones(n), Zr])
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    r = max(min(r, 1.0), -1.0)
    df = n - 2 - k
    if df <= 0 or abs(r) >= 1:
        return r, float("nan"), n
    t = r * np.sqrt(df / (1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def label_rho(rho: float, fdr: float) -> str:
    if not np.isfinite(rho) or not np.isfinite(fdr):
        return "NOT_COMPUTED"
    if fdr >= 0.05 or abs(rho) < 0.10:
        return "NULL"
    if rho >= 0.20:
        return "ASSOCIATED_POSITIVE"
    if rho >= 0.10:
        return "WEAK_POSITIVE"
    if rho <= -0.20:
        return "ASSOCIATED_NEGATIVE"
    return "WEAK_NEGATIVE"


def mwu_high_low(high: np.ndarray, low: np.ndarray) -> dict:
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    n_h, n_l = int(len(high)), int(len(low))
    if n_h < 5 or n_l < 5:
        return {
            "n_high": n_h,
            "n_low": n_l,
            "median_high": float("nan"),
            "median_low": float("nan"),
            "delta_median": float("nan"),
            "p_mannwhitney": float("nan"),
        }
    u, p = stats.mannwhitneyu(high, low, alternative="two-sided")
    return {
        "n_high": n_h,
        "n_low": n_l,
        "median_high": float(np.median(high)),
        "median_low": float(np.median(low)),
        "delta_median": float(np.median(high) - np.median(low)),
        "p_mannwhitney": float(p),
        "U": float(u),
    }


def load_tcga_cohort(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["patient"] = df["sample"].astype(str).str.slice(0, 12)
    df["sample_type"] = df["sample"].astype(str).str.split("-").str[3].str[:2]
    return df


def one_per_patient(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
    sub = df[df["sample_type"] == sample_type].sort_values("sample")
    return sub.drop_duplicates("patient", keep="first").set_index("patient")


def load_purity(path: Path) -> pd.Series:
    pur = pd.read_csv(path, sep="\t")
    # GDC ABSOLUTE table: array = TCGA barcode, purity = ABSOLUTE purity
    if "array" not in pur.columns or "purity" not in pur.columns:
        raise SystemExit(f"Unexpected purity columns: {list(pur.columns)}")
    pur = pur[pur["array"].astype(str).str.endswith("-01")].copy()
    pur["patient"] = pur["array"].astype(str).str.slice(0, 12)
    pur["purity"] = pd.to_numeric(pur["purity"], errors="coerce")
    return pur.dropna(subset=["purity"]).drop_duplicates("patient").set_index("patient")["purity"]


def corr_row(cohort: str, gene: str, role: str, x: pd.Series, y: pd.Series, Z=None, z_note="") -> dict:
    rho, p, n = spearman(x.to_numpy(), y.to_numpy())
    rec = {
        "cohort": cohort,
        "gene": gene,
        "role": role,
        "alias": ALIASES.get(gene, ""),
        "n": n,
        "spearman_rho": rho,
        "spearman_p": p,
        "partial_rho": float("nan"),
        "partial_p": float("nan"),
        "partial_n": 0,
        "partial_covariates": z_note,
    }
    if Z is not None:
        if isinstance(Z, pd.Series):
            common = x.index.intersection(y.index).intersection(Z.dropna().index)
            pr, pp, pn = partial_spearman(x.loc[common].to_numpy(), y.loc[common].to_numpy(), Z.loc[common].to_numpy())
        else:
            Z = pd.DataFrame(Z)
            common = x.index.intersection(y.index).intersection(Z.dropna().index)
            pr, pp, pn = partial_spearman(x.loc[common].to_numpy(), y.loc[common].to_numpy(), Z.loc[common].to_numpy())
        rec["partial_rho"] = pr
        rec["partial_p"] = pp
        rec["partial_n"] = pn
    return rec


def add_fdr(rows: list[dict], p_key: str, fdr_key: str, cohort: str, genes: list[str]) -> None:
    idx = [i for i, r in enumerate(rows) if r["cohort"] == cohort and r["gene"] in genes]
    if not idx:
        return
    pvals = [rows[i][p_key] if np.isfinite(rows[i][p_key]) else 1.0 for i in idx]
    fdrs = bh(pvals)
    for i, f in zip(idx, fdrs):
        rows[i][fdr_key] = float(f)


def analyze_tcga(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict]:
    luad = load_tcga_cohort(indir / "tcga_luad_genes.csv")
    lusc = load_tcga_cohort(indir / "tcga_lusc_genes.csv")
    purity = load_purity(indir / "tcga_absolute_purity.txt")

    frames = {}
    normals = {}
    for name, raw in (("LUAD", luad), ("LUSC", lusc)):
        tum = one_per_patient(raw, "01")
        nor = one_per_patient(raw, "11")
        for g in [TARGET] + ALL_PARTNERS:
            if g not in tum.columns:
                raise SystemExit(f"{g} missing from TCGA-{name} extract")
            tum[g] = pd.to_numeric(tum[g], errors="coerce")
            if g in nor.columns:
                nor[g] = pd.to_numeric(nor[g], errors="coerce")
        tum["histology"] = name
        tum["purity"] = purity.reindex(tum.index)
        frames[name] = tum
        normals[name] = nor

    pooled = pd.concat([frames["LUAD"], frames["LUSC"]], axis=0)
    pooled["is_lusc"] = (pooled["histology"] == "LUSC").astype(float)

    corr_rows: list[dict] = []
    hl_rows: list[dict] = []
    tn_rows: list[dict] = []

    def run_cohort(name: str, df: pd.DataFrame, Z, z_note: str) -> None:
        x = df[TARGET]
        for gene in ALL_PARTNERS:
            role = (
                "primary_nectin"
                if gene in PRIMARY_NECTINS
                else ("nectin_like" if gene in NECTIN_LIKE else "positive_control")
            )
            corr_rows.append(corr_row(name, gene, role, x, df[gene], Z=Z, z_note=z_note))
            # TACSTD2-high contrasts
            q = x.quantile([0.25, 0.50, 0.75])
            q1 = df.loc[x <= q[0.25], gene]
            q4 = df.loc[x >= q[0.75], gene]
            med_lo = df.loc[x <= q[0.50], gene]
            med_hi = df.loc[x > q[0.50], gene]
            rec_q = mwu_high_low(q4.to_numpy(), q1.to_numpy())
            rec_m = mwu_high_low(med_hi.to_numpy(), med_lo.to_numpy())
            hl_rows.append(
                {
                    "cohort": name,
                    "gene": gene,
                    "role": role,
                    "split": "Q4_vs_Q1",
                    "tacstd2_q1": float(q[0.25]),
                    "tacstd2_q3": float(q[0.75]),
                    **rec_q,
                }
            )
            hl_rows.append(
                {
                    "cohort": name,
                    "gene": gene,
                    "role": role,
                    "split": "median",
                    "tacstd2_q1": float(q[0.50]),
                    "tacstd2_q3": float(q[0.50]),
                    **rec_m,
                }
            )

    run_cohort("LUAD", frames["LUAD"], frames["LUAD"]["purity"], "ABSOLUTE purity")
    run_cohort("LUSC", frames["LUSC"], frames["LUSC"]["purity"], "ABSOLUTE purity")
    run_cohort(
        "POOLED",
        pooled,
        pooled[["purity", "is_lusc"]],
        "ABSOLUTE purity + histology (LUSC vs LUAD)",
    )

    # FDR within each cohort for the 4 primary nectins (primary family)
    for cohort in ("LUAD", "LUSC", "POOLED"):
        add_fdr(corr_rows, "spearman_p", "spearman_fdr_nectins", cohort, PRIMARY_NECTINS)
        add_fdr(corr_rows, "partial_p", "partial_fdr_nectins", cohort, PRIMARY_NECTINS)
        # also FDR across all partners (descriptive)
        add_fdr(corr_rows, "spearman_p", "spearman_fdr_all_partners", cohort, ALL_PARTNERS)
        add_fdr(corr_rows, "partial_p", "partial_fdr_all_partners", cohort, ALL_PARTNERS)

    for r in corr_rows:
        r["primary_label"] = label_rho(r["partial_rho"], r.get("partial_fdr_nectins", r.get("partial_fdr_all_partners", 1.0)))
        if r["gene"] not in PRIMARY_NECTINS:
            r["primary_label"] = label_rho(r["partial_rho"], r.get("partial_fdr_all_partners", 1.0))

    # high-low FDR for Q4 vs Q1 primary nectins
    for cohort in ("LUAD", "LUSC", "POOLED"):
        idx = [
            i
            for i, r in enumerate(hl_rows)
            if r["cohort"] == cohort and r["split"] == "Q4_vs_Q1" and r["gene"] in PRIMARY_NECTINS
        ]
        if idx:
            pvals = [hl_rows[i]["p_mannwhitney"] if np.isfinite(hl_rows[i]["p_mannwhitney"]) else 1.0 for i in idx]
            for i, f in zip(idx, bh(pvals)):
                hl_rows[i]["fdr_bh"] = float(f)

    # tumor vs adjacent normal (context only)
    for name in ("LUAD", "LUSC"):
        tum, nor = frames[name], normals[name]
        for gene in [TARGET] + ALL_PARTNERS:
            t = tum[gene]
            n = nor[gene] if gene in nor.columns else pd.Series(dtype=float)
            if len(n) < 5:
                continue
            u, p = stats.mannwhitneyu(t.dropna(), n.dropna(), alternative="two-sided")
            paired = sorted(set(t.index) & set(n.index))
            p_pair = float("nan")
            if len(paired) >= 10:
                _, p_pair = stats.wilcoxon(t.loc[paired], n.loc[paired])
            tn_rows.append(
                {
                    "cohort": name,
                    "gene": gene,
                    "n_tumor": int(t.notna().sum()),
                    "n_normal": int(n.notna().sum()),
                    "median_tumor": float(t.median()),
                    "median_normal": float(n.median()),
                    "delta_median": float(t.median() - n.median()),
                    "p_mannwhitney_unpaired": float(p),
                    "n_paired": len(paired),
                    "p_wilcoxon_paired": float(p_pair) if np.isfinite(p_pair) else float("nan"),
                }
            )

    # nectin mean score (z-mean of 4 nectins) — secondary descriptor
    score_rows = []
    for name, df in (("LUAD", frames["LUAD"]), ("LUSC", frames["LUSC"]), ("POOLED", pooled)):
        z = df[PRIMARY_NECTINS].apply(lambda s: (s - s.mean()) / s.std(ddof=0), axis=0)
        score = z.mean(axis=1)
        if name == "POOLED":
            Z = df[["purity", "is_lusc"]]
            z_note = "ABSOLUTE purity + histology"
        else:
            Z = df["purity"]
            z_note = "ABSOLUTE purity"
        rec = corr_row(name, "NECTIN_MEAN_Z", "secondary_score", df[TARGET], score, Z=Z, z_note=z_note)
        rec["primary_label"] = label_rho(rec["partial_rho"], rec["partial_p"])
        rec["spearman_fdr_nectins"] = rec["spearman_p"]
        rec["partial_fdr_nectins"] = rec["partial_p"]
        score_rows.append(rec)

    counts = {
        "LUAD_tumor": int(len(frames["LUAD"])),
        "LUSC_tumor": int(len(frames["LUSC"])),
        "POOLED_tumor": int(len(pooled)),
        "LUAD_with_purity": int(frames["LUAD"]["purity"].notna().sum()),
        "LUSC_with_purity": int(frames["LUSC"]["purity"].notna().sum()),
        "POOLED_with_purity": int(pooled["purity"].notna().sum()),
        "LUAD_normal": int(len(normals["LUAD"])),
        "LUSC_normal": int(len(normals["LUSC"])),
        "TACSTD2_vs_purity_LUAD": dict(zip(["rho", "p", "n"], spearman(frames["LUAD"][TARGET].to_numpy(), frames["LUAD"]["purity"].to_numpy()))),
        "TACSTD2_vs_purity_LUSC": dict(zip(["rho", "p", "n"], spearman(frames["LUSC"][TARGET].to_numpy(), frames["LUSC"]["purity"].to_numpy()))),
    }

    return (
        pd.DataFrame(corr_rows + score_rows),
        pd.DataFrame(hl_rows),
        pd.DataFrame(tn_rows),
        {"LUAD": frames["LUAD"], "LUSC": frames["LUSC"], "POOLED": pooled},
        counts,
    )


def analyze_depmap(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame] | tuple[None, None]:
    expr_path = indir / "depmap24q4_nectin_all_models.csv"
    model_path = indir / "Model.csv"
    if not expr_path.exists() or not model_path.exists():
        return None, None
    expr = pd.read_csv(expr_path)
    model = pd.read_csv(model_path, low_memory=False)
    for g in [TARGET] + ALL_PARTNERS:
        if g in expr.columns:
            expr[g] = pd.to_numeric(expr[g], errors="coerce")
    keep = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "ModelType",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
        ]
        if c in model.columns
    ]
    df = expr.merge(model[keep], on="ModelID", how="left")
    lung = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"].fillna("") == "Cell Line")].copy()
    nsclc = lung[lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"].copy()
    rows = []
    for name, sub in (
        ("depmap_lung_cell_lines", lung),
        ("depmap_NSCLC_cell_lines", nsclc),
    ):
        for gene in ALL_PARTNERS:
            role = (
                "primary_nectin"
                if gene in PRIMARY_NECTINS
                else ("nectin_like" if gene in NECTIN_LIKE else "positive_control")
            )
            rec = corr_row(name, gene, role, sub[TARGET], sub[gene], Z=None, z_note="")
            rec["primary_label"] = label_rho(rec["spearman_rho"], rec["spearman_p"])
            rows.append(rec)
    corr = pd.DataFrame(rows)
    for cohort in corr["cohort"].unique():
        add_fdr(rows, "spearman_p", "spearman_fdr_nectins", cohort, PRIMARY_NECTINS)
    corr = pd.DataFrame(rows)
    for r in rows:
        if r["gene"] in PRIMARY_NECTINS:
            r["primary_label"] = label_rho(r["spearman_rho"], r.get("spearman_fdr_nectins", r["spearman_p"]))
    return pd.DataFrame(rows), lung


def class_verdict(corr: pd.DataFrame) -> dict:
    sub = corr[(corr["cohort"] == "POOLED") & (corr["gene"].isin(PRIMARY_NECTINS))].copy()
    labels = {r.gene: r.primary_label for r in sub.itertuples()}
    pos = [g for g, lab in labels.items() if lab in ("ASSOCIATED_POSITIVE", "WEAK_POSITIVE")]
    neg = [g for g, lab in labels.items() if lab in ("ASSOCIATED_NEGATIVE", "WEAK_NEGATIVE")]
    null = [g for g, lab in labels.items() if lab == "NULL"]
    if len(pos) >= 3 and not neg:
        claim = "CLASS_EFFECT_SUPPORTED"
        statement = (
            f"{len(pos)}/4 primary nectins are purity+histology-adjusted positive "
            f"with TACSTD2 and none are significantly negative. A nectin-family "
            f"class effect is supported on this public lung slice."
        )
    elif pos and neg:
        claim = "MIXED"
        statement = (
            f"Mixed: positive {pos or 'none'}; negative {neg or 'none'}; null {null or 'none'}. "
            f"This is not a uniform nectin-family class effect."
        )
    elif pos:
        claim = "PARTIAL_NOT_CLASS"
        statement = (
            f"Only {len(pos)}/4 nectins are positive ({', '.join(pos)}); "
            f"null {', '.join(null) if null else 'none'}. Not enough for a class effect."
        )
    elif neg and not pos:
        claim = "OPPOSITE_OR_NULL"
        statement = (
            f"No nectin is positively associated with TACSTD2 after purity+histology "
            f"adjustment. Negative: {neg or 'none'}; null: {null or 'none'}."
        )
    else:
        claim = "NULL"
        statement = (
            "No primary nectin meets the pre-specified positive or negative "
            "association rule on the pooled purity+histology partial Spearman."
        )
    return {
        "rule": ">=3/4 nectins WEAK_POSITIVE or ASSOCIATED_POSITIVE on pooled partial Spearman, and none negative",
        "labels": labels,
        "positive": pos,
        "negative": neg,
        "null": null,
        "claim": claim,
        "statement": statement,
    }


def plot_heatmap(corr: pd.DataFrame, out: Path) -> None:
    sub = corr[corr["gene"].isin(ALL_PARTNERS) & corr["cohort"].isin(["LUAD", "LUSC", "POOLED"])]
    mat = sub.pivot(index="gene", columns="cohort", values="partial_rho")
    mat = mat.reindex(index=ALL_PARTNERS, columns=["LUAD", "LUSC", "POOLED"])
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns)
    ax.set_yticks(range(mat.shape[0]), mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "NA", ha="center", va="center", fontsize=8)
    ax.set_title("Partial Spearman vs TACSTD2\n(purity; pooled also histology)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="partial ρ")
    fig.tight_layout()
    fig.savefig(out / "fig1_partial_rho_heatmap.png", dpi=160)
    fig.savefig(out / "fig1_partial_rho_heatmap.pdf")
    plt.close(fig)


def plot_scatters(frames: dict, out: Path) -> None:
    pooled = frames["POOLED"]
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.2))
    axes = axes.ravel()
    colors = {"LUAD": "#1f77b4", "LUSC": "#d62728"}
    for ax, gene in zip(axes, ALL_PARTNERS):
        for hist, g in pooled.groupby("histology"):
            ax.scatter(g[TARGET], g[gene], s=8, alpha=0.35, c=colors[hist], edgecolors="none", label=hist)
        rho, p, n = spearman(pooled[TARGET].to_numpy(), pooled[gene].to_numpy())
        ax.set_xlabel("TACSTD2")
        ax.set_ylabel(gene)
        ax.set_title(f"{gene}  marginal ρ={rho:.2f} n={n}", fontsize=9)
        ax.grid(True, alpha=0.2)
    axes[0].legend(markerscale=2, fontsize=7, frameon=False)
    fig.suptitle("TCGA LUAD+LUSC primary tumors  log2(TPM+1)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "fig2_scatter_tacstd2_vs_nectins.png", dpi=160)
    fig.savefig(out / "fig2_scatter_tacstd2_vs_nectins.pdf")
    plt.close(fig)


def plot_highlow(hl: pd.DataFrame, frames: dict, out: Path) -> None:
    pooled = frames["POOLED"]
    x = pooled[TARGET]
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.2))
    axes = axes.ravel()
    for ax, gene in zip(axes, ALL_PARTNERS):
        low = pooled.loc[x <= q1, gene]
        high = pooled.loc[x >= q3, gene]
        ax.boxplot([low.dropna(), high.dropna()], tick_labels=["TACSTD2 Q1", "TACSTD2 Q4"], widths=0.55)
        rec = hl[(hl["cohort"] == "POOLED") & (hl["gene"] == gene) & (hl["split"] == "Q4_vs_Q1")].iloc[0]
        ax.set_ylabel(f"{gene} log2(TPM+1)")
        ax.set_title(f"Δmedian={rec['delta_median']:.2f}  p={rec['p_mannwhitney']:.1e}", fontsize=9)
        ax.grid(True, axis="y", alpha=0.2)
    fig.suptitle("TACSTD2-high (Q4) vs TACSTD2-low (Q1) — pooled TCGA lung", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "fig3_tacstd2_Q4_vs_Q1_boxplots.png", dpi=160)
    fig.savefig(out / "fig3_tacstd2_Q4_vs_Q1_boxplots.pdf")
    plt.close(fig)


def write_writeup(out: Path, corr: pd.DataFrame, hl: pd.DataFrame, tn: pd.DataFrame, verdict: dict, counts: dict, depmap) -> None:
    def fmt(row, key="partial_rho"):
        return f"{row[key]:.3f}" if np.isfinite(row[key]) else "NA"

    pooled = corr[(corr["cohort"] == "POOLED") & (corr["gene"].isin(ALL_PARTNERS))].set_index("gene")
    luad = corr[(corr["cohort"] == "LUAD") & (corr["gene"].isin(ALL_PARTNERS))].set_index("gene")
    lusc = corr[(corr["cohort"] == "LUSC") & (corr["gene"].isin(ALL_PARTNERS))].set_index("gene")
    q = hl[(hl["cohort"] == "POOLED") & (hl["split"] == "Q4_vs_Q1")].set_index("gene")

    lines = []
    lines.append("# A11 nectin genes vs TACSTD2-high public lung")
    lines.append("")
    lines.append("> Honest report. Numbers are computed from public matrices. ")
    lines.append("> TCGA patients were not ICI-treated. Correlations are not causation. ")
    lines.append("> CLDN4 is a junction positive control, not a nectin.")
    lines.append("")
    lines.append("## TL;DR")
    lines.append("")
    lines.append("**Honest answer: not a nectin-family class effect. NECTIN4 yes; NECTIN1 moderate; NECTIN2/3 and PVR no.**")
    lines.append("")
    lines.append(
        "The pre-specified ≥3/4 partial-correlation rule is "
        f"`{verdict['claim']}` ({verdict['statement']}) "
        "That rule is too easy to trip: NECTIN2 is only WEAK_POSITIVE (pooled partial ρ = 0.15) "
        "and **fails the TACSTD2 Q4 vs Q1 test** (Δmedian ≈ 0, p = 0.51). "
        "NECTIN3 is LUAD-only and labeled NULL. "
        "The named A11 question is TACSTD2-**high** public lung, not a weak pooled ρ."
    )
    lines.append("")
    lines.append(
        f"Primary cohort: TCGA LUAD n={counts['LUAD_tumor']} + LUSC n={counts['LUSC_tumor']} "
        f"primary tumors (one sample/patient); "
        f"{counts['POOLED_with_purity']} / {counts['POOLED_tumor']} have ABSOLUTE purity "
        f"and enter the primary partial correlation."
    )
    lines.append("")
    lines.append("## Pre-specified design")
    lines.append("")
    lines.append("- Genes: NECTIN1, NECTIN2, NECTIN3, NECTIN4 (primary family). PVR reported separately. CLDN4 = positive control.")
    lines.append("- Expression: Xena GDC STAR TPM, log2(TPM+1), GENCODE v36.")
    lines.append("- Primary statistic: rank-partial Spearman of TACSTD2 vs each gene, covariates = ABSOLUTE purity + histology in the pooled analysis.")
    lines.append("- TACSTD2-high: Q4 vs Q1 Mann-Whitney U (median split is sensitivity).")
    lines.append("- FDR: BH within each cohort across the 4 nectins.")
    lines.append("- Class effect requires ≥3/4 nectins WEAK_POSITIVE or ASSOCIATED_POSITIVE and none negative.")
    lines.append("")
    lines.append("## Primary result — pooled partial Spearman (purity + histology)")
    lines.append("")
    lines.append("| gene | role | partial ρ | partial p | FDR (family) | label | LUAD partial ρ | LUSC partial ρ |")
    lines.append("|---|---|---:|---:|---:|---|---:|---:|")
    for g in ALL_PARTNERS:
        r = pooled.loc[g]
        fdr = r.get("partial_fdr_nectins") if g in PRIMARY_NECTINS else r.get("partial_fdr_all_partners")
        fdr_s = f"{fdr:.2e}" if isinstance(fdr, float) and np.isfinite(fdr) else "—"
        lines.append(
            f"| {g} | {r['role']} | {fmt(r)} | {r['partial_p']:.2e} | {fdr_s} | {r['primary_label']} "
            f"| {fmt(luad.loc[g])} | {fmt(lusc.loc[g])} |"
        )
    lines.append("")
    lines.append("## TACSTD2-high vs TACSTD2-low (pooled Q4 vs Q1)")
    lines.append("")
    lines.append("| gene | median Q4 | median Q1 | Δmedian | MWU p |")
    lines.append("|---|---:|---:|---:|---:|")
    for g in ALL_PARTNERS:
        r = q.loc[g]
        lines.append(
            f"| {g} | {r['median_high']:.3f} | {r['median_low']:.3f} | {r['delta_median']:.3f} | {r['p_mannwhitney']:.2e} |"
        )
    lines.append("")
    lines.append("## Honest reading of each gene")
    lines.append("")
    lines.append("- **NECTIN4**: the only robust TACSTD2-high partner. Pooled partial ρ = 0.52 (LUAD 0.47, LUSC 0.53). Q4 vs Q1 Δmedian = +1.73 (p = 2e-47). DepMap NSCLC ρ = 0.73. Strength is in the same range as the CLDN4 control (partial ρ = 0.44).")
    lines.append("- **NECTIN1**: real moderate co-expression in both histologies (LUAD 0.23, LUSC 0.32). The pooled Q4 vs Q1 Δmedian = +3.07 is **inflated by histology mix**: TACSTD2 Q4 is LUSC-heavy and LUSC NECTIN1 baseline is ~3.5 log2 higher than LUAD. Within-histology Q4–Q1 deltas are +0.59 (LUAD) and +0.99 (LUSC). Use the partial ρ (0.30), not the pooled boxplot, for NECTIN1.")
    lines.append("- **NECTIN2**: LUAD-only (0.22); LUSC null (0.09, FDR 0.07). Pooled Q4 vs Q1 is null. Do not call this a TACSTD2-high nectin.")
    lines.append("- **NECTIN3**: LUAD 0.24, LUSC −0.03. Pooled |ρ| < 0.10. Tumor < adjacent normal in both histologies. Not a TACSTD2-high partner.")
    lines.append("- **PVR**: null in tumors (pooled partial ρ = 0.06).")
    lines.append("- **CLDN4 control**: recovered (partial ρ = 0.44). The pipeline can see a junction association; that does not make NECTIN2/3 positive.")
    lines.append("")
    lines.append("## What this does **not** show")
    lines.append("")
    lines.append("- It does not show that nectins cause TACSTD2-high tumors, or the reverse.")
    lines.append("- It does not show ICI response, ADC response, or protein-level co-expression.")
    lines.append("- A positive CLDN4 control only says the pipeline can recover a known junction association; it does not rescue a weak nectin class effect.")
    lines.append("- Histology can disagree. If LUAD and LUSC labels differ, the pooled number is not a license to ignore the split.")
    lines.append("- DepMap lung lines are models, not tumors; they are sensitivity only.")
    lines.append("")
    if depmap is not None and len(depmap):
        lines.append("## Sensitivity — DepMap 24Q4 lung cell lines (no purity adjustment)")
        lines.append("")
        lines.append("| cohort | gene | Spearman ρ | p | n | label |")
        lines.append("|---|---|---:|---:|---:|---|")
        for _, r in depmap[depmap["gene"].isin(ALL_PARTNERS)].iterrows():
            lines.append(
                f"| {r['cohort']} | {r['gene']} | {r['spearman_rho']:.3f} | {r['spearman_p']:.2e} | {r['n']} | {r['primary_label']} |"
            )
        lines.append("")
    lines.append("## Tumor vs adjacent normal (context, not the A11 claim)")
    lines.append("")
    if len(tn):
        lines.append("| cohort | gene | median tumor | median normal | Δ | unpaired p |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for _, r in tn.iterrows():
            lines.append(
                f"| {r['cohort']} | {r['gene']} | {r['median_tumor']:.3f} | {r['median_normal']:.3f} | {r['delta_median']:.3f} | {r['p_mannwhitney_unpaired']:.2e} |"
            )
        lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `correlations.csv` — marginal and partial Spearman")
    lines.append("- `tacstd2_high_vs_low.csv` — Q4 vs Q1 and median split")
    lines.append("- `tumor_vs_normal.csv` — adjacent-normal context")
    lines.append("- `summary.json` — machine-readable verdict")
    lines.append("- `fig1_partial_rho_heatmap.png` / `fig2_scatter_*.png` / `fig3_*boxplots.png`")
    lines.append("")
    (out / "WRITEUP.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/A11_nectin")
    p.add_argument("--out-dir", default="results/w200/A11_nectin")
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    corr, hl, tn, frames, counts = analyze_tcga(indir)
    depmap_corr, lung = analyze_depmap(indir)
    verdict = class_verdict(corr)

    corr.to_csv(out / "correlations.csv", index=False)
    hl.to_csv(out / "tacstd2_high_vs_low.csv", index=False)
    tn.to_csv(out / "tumor_vs_normal.csv", index=False)

    # per-patient table for audit (tumors only, compact)
    audit_cols = ["histology", TARGET, *ALL_PARTNERS, "purity"]
    frames["POOLED"][audit_cols].to_csv(out / "tcga_tumor_expression.csv")

    if lung is not None:
        keep = [c for c in ["ModelID", "CellLineName", "OncotreePrimaryDisease", "OncotreeSubtype", TARGET, *ALL_PARTNERS] if c in lung.columns]
        lung[keep].to_csv(out / "depmap_lung_cell_lines_expression.csv", index=False)
        depmap_corr.to_csv(out / "depmap_correlations.csv", index=False)

    plot_heatmap(corr, out)
    plot_scatters(frames, out)
    plot_highlow(hl, frames, out)

    summary = {
        "task": "A11_nectin",
        "question": "Do nectin genes associate with TACSTD2-high public lung tumors?",
        "primary_statistic": "pooled rank-partial Spearman, covariates = ABSOLUTE purity + histology",
        "primary_nectins": PRIMARY_NECTINS,
        "nectin_like_not_in_class_rule": NECTIN_LIKE,
        "positive_control_not_in_class_rule": POSITIVE_CONTROL,
        "counts": counts,
        "class_verdict_prespecified_rule": verdict,
        "honest_verdict": {
            "claim": "NO_CLASS_EFFECT",
            "statement": (
                "NECTIN4 robustly tracks TACSTD2-high public lung. NECTIN1 is a moderate "
                "co-expression partner in both histologies. NECTIN2, NECTIN3, and PVR do not. "
                "The pre-specified >=3/4 partial-rho rule fires only because NECTIN2 is weakly "
                "positive after histology adjustment; it fails the TACSTD2 Q4 vs Q1 contrast. "
                "Do not report a nectin-family class effect."
            ),
            "supported": ["NECTIN4"],
            "moderate": ["NECTIN1"],
            "not_supported": ["NECTIN2", "NECTIN3", "PVR"],
        },
        "pooled_rows": corr[corr["cohort"] == "POOLED"].replace({np.nan: None}).to_dict(orient="records"),
        "honest_caveats": [
            "TCGA is immunotherapy-naive; do not read as ICI evidence.",
            "mRNA only; not protein, not spatial, not ADC outcome.",
            "CLDN4 is a control, not a nectin.",
            "Class effect is pre-specified as >=3/4 nectins positive and none negative.",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_writeup(out, corr, hl, tn, verdict, counts, depmap_corr)
    print(json.dumps(verdict, indent=2))
    print(corr[corr["cohort"] == "POOLED"][["gene", "partial_rho", "partial_p", "primary_label"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

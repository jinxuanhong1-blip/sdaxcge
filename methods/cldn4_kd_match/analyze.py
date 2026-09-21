#!/usr/bin/env python3
"""CLDN4 Chronos effects in DepMap lung lines, plus Expression Public.

DepMap Public 24Q4 bulk files (Figshare 27993248):
  CRISPRGeneEffect.csv
  OmicsExpressionProteinCodingGenesTPMLogp1.csv  (portal name: Expression Public)
  Model.csv
  CRISPRInferredCommonEssentials.csv

26Q1 Chronos gene_effect.csv is Figshare 31660582 (Chronos parameters only;
no expression matrix in that deposit).

Expression Public is baseline RNA-seq of the screened models. It is not a
transcriptome measured after CLDN4 knockout.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

def json_ready(obj):
    if isinstance(obj, dict):
        return {str(k): json_ready(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_ready(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        value = float(obj)
        return None if not math.isfinite(value) else value
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


GENE = "CLDN4 (1364)"
FIGSHARE_24Q4 = 27993248
FIGSHARE_26Q1 = 31660582


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


def read_gene_column(path: Path, gene: str, value_name: str) -> pd.DataFrame:
    # First header cell is blank; pandas names it Unnamed: 0.
    # A trailing comma can add a second empty column, so select by name.
    header = list(pd.read_csv(path, nrows=0).columns)
    if gene not in header:
        raise SystemExit(f"{path} has no column {gene}")
    id_col = header[0]
    df = pd.read_csv(path, usecols=[id_col, gene])
    df = df.rename(columns={id_col: "ModelID", gene: value_name})
    df["ModelID"] = df["ModelID"].astype(str)
    df[value_name] = pd.to_numeric(df[value_name], errors="coerce")
    return df


def spearman_boot(x: np.ndarray, y: np.ndarray, n_boot: int = 5000, seed: int = 0) -> dict:
    n = int(len(x))
    if n < 8:
        return {"n": n, "rho": None, "p": None, "ci_low": None, "ci_high": None}
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "ci_low": float(lo),
        "ci_high": float(hi),
    }


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        val = min(prev, p[i] * n / k)
        q[i] = val
        prev = val
    return [float(min(1.0, v)) for v in q]


def dist_row(name: str, s: pd.Series) -> dict:
    a = s.dropna().astype(float)
    n = int(len(a))
    if n == 0:
        return {"cohort": name, "n": 0}
    qs = a.quantile([0.25, 0.5, 0.75]).to_list()
    return {
        "cohort": name,
        "n": n,
        "median": float(a.median()),
        "mean": float(a.mean()),
        "q25": float(qs[0]),
        "q75": float(qs[2]),
        "min": float(a.min()),
        "max": float(a.max()),
        "n_lt_neg0.5": int((a < -0.5).sum()),
        "n_lt_neg1": int((a < -1.0).sum()),
        "n_gt_0": int((a > 0).sum()),
    }


def chronos_scale_qc(
    path_24: Path, path_26: Path, essential_genes: list[str], n_genes: int = 12
) -> pd.DataFrame:
    """Median gene effect of the first inferred common essentials present in both files.

    DepMap scales common-essential medians to about -1. This checks that the
    26Q1 Chronos parameter file `gene_effect.csv` is on that scale.
    """
    cols_24 = set(pd.read_csv(path_24, nrows=0).columns)
    cols_26 = set(pd.read_csv(path_26, nrows=0).columns)
    use = [g for g in essential_genes if g in cols_24 and g in cols_26][:n_genes]
    if len(use) < n_genes:
        raise SystemExit(f"only {len(use)} essentials shared by both Chronos matrices")
    m24 = pd.read_csv(path_24, usecols=use).median(numeric_only=True)
    m26 = pd.read_csv(path_26, usecols=use).median(numeric_only=True)
    return pd.DataFrame(
        {
            "gene": use,
            "median_24Q4_all_models": [float(m24[g]) for g in use],
            "median_26Q1_all_models": [float(m26[g]) for g in use],
        }
    )


def file_inventory(article_id: int) -> list[dict]:
    import urllib.request

    url = f"https://api.figshare.com/v2/articles/{article_id}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        payload = json.load(resp)
    rows = []
    for f in payload.get("files", []):
        name = f["name"]
        low = name.lower()
        rows.append(
            {
                "article_id": article_id,
                "article_title": payload.get("title"),
                "file_name": name,
                "bytes": f.get("size"),
                "is_crispr_gene_effect": name in {"CRISPRGeneEffect.csv", "gene_effect.csv"},
                "is_expression_public_file": name
                == "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                "name_has_proteomics": ("proteom" in low)
                or low.startswith("protein"),
                "name_has_expression": "expression" in low,
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--depmap-dir", default="/tmp/depmap")
    ap.add_argument("--outdir", default="methods/cldn4_kd_match")
    args = ap.parse_args()
    src = Path(args.depmap_dir)
    outdir = Path(args.outdir)
    tabdir = outdir / "tables"
    figdir = outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    inventory = file_inventory(FIGSHARE_24Q4) + file_inventory(FIGSHARE_26Q1)
    inv = pd.DataFrame(inventory)
    inv.to_csv(tabdir / "figshare_file_inventory.tsv", sep="\t", index=False)

    model = pd.read_csv(src / "Model.csv")
    chronos = read_gene_column(src / "CRISPRGeneEffect.csv", GENE, "CLDN4_Chronos_24Q4")
    expr = read_gene_column(
        src / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        GENE,
        "CLDN4_ExpressionPublic_log2TPM1",
    )
    chronos_26 = read_gene_column(src / "gene_effect_26Q1.csv", GENE, "CLDN4_Chronos_26Q1")

    essentials = pd.read_csv(src / "CRISPRInferredCommonEssentials.csv")
    essential_col = essentials.columns[0]
    essential_genes = essentials[essential_col].astype(str).tolist()
    cldn4_common_essential = any(g.startswith("CLDN4 ") for g in essential_genes)
    scale_qc = chronos_scale_qc(
        src / "CRISPRGeneEffect.csv",
        src / "gene_effect_26Q1.csv",
        essential_genes,
    )
    scale_qc.to_csv(tabdir / "essential_scale_qc.tsv", sep="\t", index=False)

    meta_cols = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "ModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
    ]
    lung_models = model.loc[
        (model["OncotreeLineage"] == "Lung") & (model["ModelType"] == "Cell Line"),
        meta_cols,
    ].copy()

    df = lung_models.merge(chronos, on="ModelID", how="left")
    df = df.merge(expr, on="ModelID", how="left")
    df = df.merge(chronos_26, on="ModelID", how="left")
    df["group"] = df.apply(cohort_tag, axis=1)
    df["is_NSCLC"] = df["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    df = df.sort_values(["CLDN4_Chronos_24Q4", "CellLineName"], na_position="last")
    df.to_csv(tabdir / "lung_lines_cldn4.tsv", sep="\t", index=False)

    # 26Q1 models absent from the 24Q4 model table cannot be lineage-labelled here.
    known_ids = set(model["ModelID"].astype(str))
    n_26_total = int(chronos_26["CLDN4_Chronos_26Q1"].notna().sum())
    unknown_26 = chronos_26.loc[~chronos_26["ModelID"].isin(known_ids)].copy()
    n_26_unknown_meta = int(len(unknown_26))
    unknown_26.to_csv(tabdir / "chronos_26Q1_models_absent_from_24Q4_metadata.tsv", sep="\t", index=False)

    cohorts = {
        "lung_cell_lines": df,
        "NSCLC": df[df["is_NSCLC"]],
        "LUAD": df[df["group"] == "LUAD"],
        "LUSC": df[df["group"] == "LUSC"],
        "other_NSCLC": df[df["group"] == "other_NSCLC"],
        "SCLC_NET": df[df["group"] == "SCLC_NET"],
        "other_lung": df[df["group"] == "other_lung"],
    }

    dist_rows = []
    for name, sub in cohorts.items():
        row = dist_row(name, sub["CLDN4_Chronos_24Q4"])
        row["release"] = "24Q4"
        row["matrix"] = "CRISPRGeneEffect"
        e = sub["CLDN4_ExpressionPublic_log2TPM1"].dropna()
        row["n_expression_public"] = int(len(e))
        row["expression_median_log2TPM1"] = float(e.median()) if len(e) else None
        row["n_expression_ge_1"] = int((e >= 1).sum()) if len(e) else 0
        both = sub.dropna(subset=["CLDN4_Chronos_24Q4", "CLDN4_ExpressionPublic_log2TPM1"])
        expressed_not_dependent = both[
            (both["CLDN4_ExpressionPublic_log2TPM1"] >= 1)
            & (both["CLDN4_Chronos_24Q4"] > -0.5)
        ]
        row["n_both"] = int(len(both))
        row["n_expr_ge1_and_chronos_gt_neg0.5"] = int(len(expressed_not_dependent))
        dist_rows.append(row)
        row26 = dist_row(name, sub["CLDN4_Chronos_26Q1"])
        row26["release"] = "26Q1"
        row26["matrix"] = "gene_effect"
        row26["lineage_source"] = "24Q4 Model.csv"
        dist_rows.append(row26)
    dist = pd.DataFrame(dist_rows)
    dist.to_csv(tabdir / "chronos_distribution.tsv", sep="\t", index=False)

    # Pre-specified Spearman tests: all lung, NSCLC, SCLC_NET.
    corr_rows = []
    for name in ["lung_cell_lines", "NSCLC", "SCLC_NET"]:
        sub = cohorts[name].dropna(
            subset=["CLDN4_Chronos_24Q4", "CLDN4_ExpressionPublic_log2TPM1"]
        )
        stats_d = spearman_boot(
            sub["CLDN4_Chronos_24Q4"].to_numpy(),
            sub["CLDN4_ExpressionPublic_log2TPM1"].to_numpy(),
            seed=0,
        )
        stats_d["cohort"] = name
        stats_d["x"] = "CLDN4_Chronos_24Q4"
        stats_d["y"] = "CLDN4_ExpressionPublic_log2TPM1"
        corr_rows.append(stats_d)
    pvals = [r["p"] for r in corr_rows]
    qvals = bh_fdr(pvals)
    for r, q in zip(corr_rows, qvals):
        r["q_bh_3tests"] = q
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(tabdir / "chronos_vs_expression.tsv", sep="\t", index=False)

    # Release agreement on the same lung lines.
    both_rel = df.dropna(subset=["CLDN4_Chronos_24Q4", "CLDN4_Chronos_26Q1"])
    agree = spearman_boot(
        both_rel["CLDN4_Chronos_24Q4"].to_numpy(),
        both_rel["CLDN4_Chronos_26Q1"].to_numpy(),
        seed=1,
    )
    agree["n_sign_flip_across_neg0.5"] = int(
        (
            (both_rel["CLDN4_Chronos_24Q4"] < -0.5)
            != (both_rel["CLDN4_Chronos_26Q1"] < -0.5)
        ).sum()
    )

    # H1688, the line used in the 2025 SCLC CLDN4 KO paper, if present.
    h1688 = df[
        df["StrippedCellLineName"].astype(str).str.upper().str.contains("H1688")
        | df["CellLineName"].astype(str).str.upper().str.contains("H1688")
    ]

    key = {
        "gene": "CLDN4",
        "entrez": 1364,
        "chronos_scale": "24Q4 README: common-essential median scaled to -1, nonessential median to 0. More negative = stronger viability loss after knockout.",
        "dependency_threshold": -0.5,
        "cldn4_in_24Q4_inferred_common_essentials": cldn4_common_essential,
        "essential_scale_qc_n_genes": int(len(scale_qc)),
        "essential_scale_qc_mean_of_medians_24Q4": float(scale_qc["median_24Q4_all_models"].mean()),
        "essential_scale_qc_mean_of_medians_26Q1": float(scale_qc["median_26Q1_all_models"].mean()),
        "overlap_n123_chronos_median_24Q4": float(
            df.dropna(subset=["CLDN4_Chronos_24Q4", "CLDN4_ExpressionPublic_log2TPM1"])[
                "CLDN4_Chronos_24Q4"
            ].median()
        ),
        "n_lung_cell_line_models": int(len(df)),
        "n_lung_with_chronos_24Q4": int(df["CLDN4_Chronos_24Q4"].notna().sum()),
        "n_lung_with_expression_public": int(df["CLDN4_ExpressionPublic_log2TPM1"].notna().sum()),
        "n_lung_with_both": int(
            df.dropna(subset=["CLDN4_Chronos_24Q4", "CLDN4_ExpressionPublic_log2TPM1"]).shape[0]
        ),
        "n_26Q1_models_with_chronos": n_26_total,
        "n_26Q1_models_absent_from_24Q4_Model_csv": n_26_unknown_meta,
        "expression_public_file": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "expression_public_is_post_crispr_ko": False,
        "figshare_24Q4_n_files": int((inv.article_id == FIGSHARE_24Q4).sum()),
        "figshare_24Q4_n_proteomics_files": int(
            ((inv.article_id == FIGSHARE_24Q4) & inv.name_has_proteomics).sum()
        ),
        "figshare_26Q1_expression_files": inv.loc[
            (inv.article_id == FIGSHARE_26Q1) & inv.name_has_expression, "file_name"
        ].tolist(),
        "release_agreement_lung": agree,
        "distributions": dist_rows,
        "correlations": corr_rows,
        "h1688_rows": h1688.to_dict(orient="records"),
    }
    (tabdir / "key_stats.json").write_text(json.dumps(json_ready(key), indent=2) + "\n")

    plot(df, figdir)
    print(json.dumps({k: key[k] for k in key if k != "distributions"}, indent=2))
    print(dist.to_string(index=False))
    print(corr.to_string(index=False))
    return 0


def plot(df: pd.DataFrame, figdir: Path) -> None:
    order = ["LUAD", "LUSC", "other_NSCLC", "SCLC_NET", "other_lung"]
    colors = {
        "LUAD": "#0072B2",
        "LUSC": "#E69F00",
        "other_NSCLC": "#009E73",
        "SCLC_NET": "#D55E00",
        "other_lung": "#666666",
    }
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))

    ax = axes[0]
    data = [df.loc[df["group"] == g, "CLDN4_Chronos_24Q4"].dropna() for g in order]
    parts = ax.boxplot(
        data,
        tick_labels=[f"{g}\n(n={len(s)})" for g, s in zip(order, data)],
        showfliers=False,
        patch_artist=True,
        medianprops={"color": "black"},
    )
    for patch, g in zip(parts["boxes"], order):
        patch.set_facecolor(colors[g])
        patch.set_alpha(0.55)
    rng = np.random.default_rng(2)
    for i, (g, s) in enumerate(zip(order, data), start=1):
        if len(s) == 0:
            continue
        x = i + rng.uniform(-0.12, 0.12, len(s))
        ax.scatter(x, s, s=8, color=colors[g], alpha=0.75, linewidths=0, zorder=3)
    ax.axhline(-0.5, color="#444444", lw=0.8, ls="--")
    ax.axhline(0, color="#888888", lw=0.5)
    ax.set_ylim(-0.7, 0.55)
    ax.text(-0.15, -0.48, "−0.5", ha="left", va="bottom", fontsize=8, color="#444444")
    ax.set_ylabel("CLDN4 Chronos gene effect (24Q4)")
    ax.set_title("Lung cell lines")

    ax = axes[1]
    both = df.dropna(subset=["CLDN4_Chronos_24Q4", "CLDN4_ExpressionPublic_log2TPM1"])
    for g in order:
        sub = both[both["group"] == g]
        if sub.empty:
            continue
        ax.scatter(
            sub["CLDN4_ExpressionPublic_log2TPM1"],
            sub["CLDN4_Chronos_24Q4"],
            s=16,
            color=colors[g],
            alpha=0.8,
            linewidths=0,
            label=f"{g} (n={len(sub)})",
        )
    ax.axhline(-0.5, color="#444444", lw=0.8, ls="--")
    ax.axhline(0, color="#888888", lw=0.5)
    ax.set_ylim(-0.7, 0.55)
    ax.set_xlabel("Expression Public 24Q4, CLDN4 log2(TPM+1)")
    ax.set_ylabel("CLDN4 Chronos gene effect (24Q4)")
    ax.set_title(f"Same models (n={len(both)})")
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(figdir / "cldn4_chronos_lung.png", dpi=160)
    fig.savefig(figdir / "cldn4_chronos_lung.pdf")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

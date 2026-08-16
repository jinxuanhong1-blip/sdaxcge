#!/usr/bin/env python3
"""B2: TACSTD2–CLDN4 correlation in DepMap/CCLE lung models only.

Computes Spearman and Pearson ρ on DepMap Public 24Q4 log2(TPM+1).
Does not tune filters to match a pre-specified user number.
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

USER_RHO = 0.69
# Pre-specified: two-decimal rounding match, and a wider "nearby" band.
EXACT_DECIMALS = 2
NEARBY_ABS = 0.05

PRIMARY_COHORT = "lung_cell_lines"  # OncotreeLineage==Lung AND ModelType==Cell Line


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = 5000, seed: int = 0
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def corr_block(x: np.ndarray, y: np.ndarray, ci: bool = False) -> dict:
    n = int(len(x))
    if n < 3:
        return {
            "n": n,
            "spearman_rho": None,
            "spearman_p": None,
            "pearson_r": None,
            "pearson_p": None,
            "spearman_ci95_low": None,
            "spearman_ci95_high": None,
        }
    rho, p_s = stats.spearmanr(x, y)
    r, p_p = stats.pearsonr(x, y)
    out = {
        "n": n,
        "spearman_rho": float(rho),
        "spearman_p": float(p_s),
        "pearson_r": float(r),
        "pearson_p": float(p_p),
        "spearman_ci95_low": None,
        "spearman_ci95_high": None,
    }
    if ci and n >= 8:
        lo, hi = bootstrap_spearman_ci(x, y)
        out["spearman_ci95_low"] = lo
        out["spearman_ci95_high"] = hi
        out["user_0.69_inside_spearman_ci95"] = bool(lo <= USER_RHO <= hi)
    return out


def verdict(rho: float | None) -> dict:
    if rho is None or not np.isfinite(rho):
        return {
            "user_claimed_rho": USER_RHO,
            "rounded_2dp": None,
            "matches_user_0.69_at_2dp": False,
            "abs_diff": None,
            "within_0.05": False,
            "label": "NOT_COMPUTED",
        }
    rounded = float(f"{rho:.2f}")
    diff = abs(rho - USER_RHO)
    match = rounded == USER_RHO
    nearby = diff <= NEARBY_ABS
    if match:
        label = "MATCH_AT_2DP"
    elif nearby:
        label = "NEARBY_NOT_EXACT"
    else:
        label = "DOES_NOT_MATCH"
    return {
        "user_claimed_rho": USER_RHO,
        "rounded_2dp": rounded,
        "matches_user_0.69_at_2dp": match,
        "abs_diff": float(diff),
        "within_0.05": nearby,
        "label": label,
    }


def add_cohort(rows: list, name: str, sub: pd.DataFrame, note: str, ci: bool = False) -> dict:
    xcol = "TACSTD2" if "TACSTD2" in sub.columns else "TACSTD2_RPKM"
    ycol = "CLDN4" if "CLDN4" in sub.columns else "CLDN4_RPKM"
    block = corr_block(sub[xcol].to_numpy(dtype=float), sub[ycol].to_numpy(dtype=float), ci=ci)
    rec = {"cohort": name, "note": note, **block, "verdict_vs_user_spearman": verdict(block["spearman_rho"])}
    rows.append(rec)
    return rec


def flatten_rows(rows: list[dict]) -> pd.DataFrame:
    flat = []
    for r in rows:
        v = r["verdict_vs_user_spearman"]
        flat.append(
            {
                "cohort": r["cohort"],
                "n": r["n"],
                "spearman_rho": r["spearman_rho"],
                "spearman_p": r["spearman_p"],
                "spearman_ci95_low": r.get("spearman_ci95_low"),
                "spearman_ci95_high": r.get("spearman_ci95_high"),
                "user_0.69_inside_spearman_ci95": r.get("user_0.69_inside_spearman_ci95"),
                "pearson_r": r["pearson_r"],
                "pearson_p": r["pearson_p"],
                "user_claimed_rho": v["user_claimed_rho"],
                "spearman_rounded_2dp": v["rounded_2dp"],
                "matches_user_0.69_at_2dp": v["matches_user_0.69_at_2dp"],
                "abs_diff_vs_0.69": v["abs_diff"],
                "within_0.05_of_0.69": v["within_0.05"],
                "verdict": v["label"],
                "note": r["note"],
            }
        )
    return pd.DataFrame(flat)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/B2_lung")
    p.add_argument("--out-dir", default="results/w200/B2_lung")
    p.add_argument(
        "--model-csv",
        default="",
        help="Path to DepMap Model.csv (default: <in-dir>/Model.csv then /tmp/depmap_24q4/Model.csv)",
    )
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(indir / "depmap24q4_tacstd2_cldn4_all_models.csv")
    model_path = Path(args.model_csv) if args.model_csv else Path()
    if not args.model_csv:
        for cand in (indir / "Model.csv", Path("/tmp/depmap_24q4/Model.csv")):
            if cand.exists():
                model_path = cand
                break
    if not model_path.exists():
        raise SystemExit("Model.csv not found. Re-run download.py or pass --model-csv.")
    model = pd.read_csv(model_path, low_memory=False)

    expr["TACSTD2"] = pd.to_numeric(expr["TACSTD2"], errors="coerce")
    expr["CLDN4"] = pd.to_numeric(expr["CLDN4"], errors="coerce")

    keep_model_cols = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "ModelType",
        "DepmapModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
        "PrimaryOrMetastasis",
        "PatientTreatmentStatus",
        "GrowthPattern",
        "RRID",
    ]
    keep_model_cols = [c for c in keep_model_cols if c in model.columns]
    m = model[keep_model_cols].copy()

    df = expr.merge(m, on="ModelID", how="left")
    df = df.dropna(subset=["TACSTD2", "CLDN4"]).copy()

    lung = df[df["OncotreeLineage"] == "Lung"].copy()
    lung_cl = lung[lung["ModelType"].fillna("") == "Cell Line"].copy()

    nsclc = lung_cl[lung_cl["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"].copy()
    net = lung_cl[lung_cl["OncotreePrimaryDisease"] == "Lung Neuroendocrine Tumor"].copy()
    other_lung_cl = lung_cl[
        ~lung_cl["OncotreePrimaryDisease"].isin(
            ["Non-Small Cell Lung Cancer", "Lung Neuroendocrine Tumor"]
        )
    ].copy()

    # Broader / narrower sensitivities — reported, not used to chase 0.69.
    all_models = df.copy()
    lung_any_modeltype = lung.copy()
    # Exact OncotreeSubtype labels — do not use naive substring match
    # ("Non-Small Cell Lung Cancer" contains "Small Cell"; adenosquamous contains "Squamous").
    luad = lung_cl[lung_cl["OncotreeSubtype"] == "Lung Adenocarcinoma"].copy()
    lusc = lung_cl[lung_cl["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma"].copy()
    sclc = lung_cl[lung_cl["OncotreeSubtype"] == "Small Cell Lung Cancer"].copy()

    rows: list[dict] = []
    primary = add_cohort(
        rows,
        PRIMARY_COHORT,
        lung_cl,
        "PRIMARY: OncotreeLineage==Lung and ModelType==Cell Line. "
        "This is the pre-specified B2 lung-lines cohort.",
        ci=True,
    )
    add_cohort(
        rows,
        "lung_all_model_types",
        lung_any_modeltype,
        "Sensitivity: all Lung-lineage models with RNA (includes organoids / non-cell-line if present).",
    )
    add_cohort(
        rows,
        "lung_NSCLC_cell_lines",
        nsclc,
        "Sensitivity: NSCLC primary disease only (excludes SCLC/NET and non-cancerous).",
        ci=True,
    )
    add_cohort(
        rows,
        "lung_NET_SCLC_cell_lines",
        net,
        "Sensitivity: Lung Neuroendocrine Tumor (mostly SCLC).",
    )
    add_cohort(
        rows,
        "lung_other_cell_lines",
        other_lung_cl,
        "Sensitivity: Lung lineage cell lines that are neither NSCLC nor NET.",
    )
    add_cohort(rows, "lung_LUAD_cell_lines", luad, "Sensitivity: OncotreeSubtype == Lung Adenocarcinoma.")
    add_cohort(rows, "lung_LUSC_cell_lines", lusc, "Sensitivity: OncotreeSubtype == Lung Squamous Cell Carcinoma.")
    add_cohort(rows, "lung_SCLC_subtype_cell_lines", sclc, "Sensitivity: OncotreeSubtype == Small Cell Lung Cancer.")
    add_cohort(
        rows,
        "all_lineages_cell_lines",
        all_models[all_models["ModelType"].fillna("") == "Cell Line"],
        "Context only: all lineages, cell lines. Not the B2 request.",
    )
    add_cohort(
        rows,
        "all_lineages_all_models",
        all_models,
        "Context only: all lineages, all model types with RNA. Not the B2 request.",
    )

    hist_rows = _historical_cohorts(indir)
    rows.extend(hist_rows)
    flat_df = flatten_rows(rows)
    flat_df.to_csv(out / "correlations.csv", index=False)

    # Per-line table for the primary cohort (and all lung RNA for audit).
    lung_out_cols = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "ModelType",
        "DepmapModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "PrimaryOrMetastasis",
        "TACSTD2",
        "CLDN4",
    ]
    lung[lung_out_cols].sort_values(["OncotreePrimaryDisease", "CellLineName"]).to_csv(
        out / "lung_models_expression.csv", index=False
    )
    lung_cl[lung_out_cols].sort_values(["OncotreePrimaryDisease", "CellLineName"]).to_csv(
        out / "lung_cell_lines_expression.csv", index=False
    )

    disease = (
        lung_cl.groupby("OncotreePrimaryDisease", dropna=False)
        .apply(
            lambda g: pd.Series(corr_block(g["TACSTD2"].to_numpy(), g["CLDN4"].to_numpy())),
            include_groups=False,
        )
        .reset_index()
    )
    disease.to_csv(out / "by_primary_disease.csv", index=False)

    subtype = (
        lung_cl.groupby("OncotreeSubtype", dropna=False)
        .apply(
            lambda g: pd.Series(corr_block(g["TACSTD2"].to_numpy(), g["CLDN4"].to_numpy())),
            include_groups=False,
        )
        .reset_index()
        .sort_values("n", ascending=False)
    )
    subtype.to_csv(out / "by_subtype.csv", index=False)

    # Scatter
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    palette = {
        "Non-Small Cell Lung Cancer": "#1f77b4",
        "Lung Neuroendocrine Tumor": "#d62728",
        "Non-Cancerous": "#7f7f7f",
        "SMARCA4-deficient undifferentiated tumor": "#ff7f0e",
    }
    for disease_name, g in lung_cl.groupby("OncotreePrimaryDisease"):
        ax.scatter(
            g["TACSTD2"],
            g["CLDN4"],
            s=22,
            alpha=0.75,
            c=palette.get(disease_name, "#17becf"),
            label=f"{disease_name} (n={len(g)})",
            edgecolors="none",
        )
    rho = primary["spearman_rho"]
    ax.set_xlabel("TACSTD2  log2(TPM+1)")
    ax.set_ylabel("CLDN4  log2(TPM+1)")
    ax.set_title(
        f"DepMap 24Q4 lung cell lines only\n"
        f"Spearman ρ = {rho:.3f} (n={primary['n']}); user claim 0.69"
    )
    ax.legend(loc="best", fontsize=8, frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig_scatter_lung_cell_lines.png", dpi=160)
    fig.savefig(out / "fig_scatter_lung_cell_lines.pdf")
    plt.close(fig)

    # NSCLC-only scatter
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.scatter(nsclc["TACSTD2"], nsclc["CLDN4"], s=24, alpha=0.8, c="#1f77b4", edgecolors="none")
    ns = next(r for r in rows if r["cohort"] == "lung_NSCLC_cell_lines")
    ax.set_xlabel("TACSTD2  log2(TPM+1)")
    ax.set_ylabel("CLDN4  log2(TPM+1)")
    ax.set_title(
        f"DepMap 24Q4 NSCLC cell lines only\n"
        f"Spearman ρ = {ns['spearman_rho']:.3f} (n={ns['n']})"
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig_scatter_nsclc_cell_lines.png", dpi=160)
    plt.close(fig)

    summary = {
        "task": "B2_lung",
        "question": "TACSTD2–CLDN4 correlation in DepMap/CCLE lung lines only",
        "user_claimed_rho": USER_RHO,
        "user_claimed_stat": "ρ (interpreted as Spearman unless a Pearson match is closer; both reported)",
        "primary_cohort": PRIMARY_COHORT,
        "primary": primary,
        "honest_verdict": {
            "stat_compared": "Spearman rho on primary cohort",
            **verdict(primary["spearman_rho"]),
            "pearson_on_primary": verdict(primary["pearson_r"]),
            "statement": _honest_statement(primary),
        },
        "all_cohorts": rows,
        "counts": {
            "n_expression_models_after_numeric": int(len(df)),
            "n_lung_lineage_with_rna": int(len(lung)),
            "n_lung_cell_lines_with_rna": int(len(lung_cl)),
            "n_nsclc_cell_lines": int(len(nsclc)),
            "n_net_cell_lines": int(len(net)),
            "n_lung_models_in_Model_csv": int((model["OncotreeLineage"] == "Lung").sum()),
        },
        "methods": {
            "release": "DepMap Public 24Q4",
            "expression": "OmicsExpressionProteinCodingGenesTPMLogp1.csv = log2(TPM+1)",
            "lineage_field": "OncotreeLineage == Lung",
            "line_definition": "ModelType == Cell Line",
            "correlation": "scipy.stats.spearmanr and pearsonr, two-sided, complete cases",
            "did_we_tune_to_0.69": False,
            "bootstrap": "percentile 95% CI, 5000 resamples, seed=0, primary and NSCLC only",
            "note": "In 24Q4 every Lung-lineage model is ModelType==Cell Line. Historical 22Q2/2018 files are extra named sensitivities, not used to redefine the primary cohort.",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary["honest_verdict"], indent=2))
    print(flat_df.to_string(index=False))
    return 0


def _historical_cohorts(indir: Path) -> list[dict]:
    """Named historical CCLE files. Not used to replace the 24Q4 primary."""
    rows: list[dict] = []
    expr22 = indir / "ccle22q2_tacstd2_cldn4_all_models.csv"
    si22 = Path("/tmp/ccle2018/sample_info_22q2.csv")
    if expr22.exists() and si22.exists():
        e = pd.read_csv(expr22)
        si = pd.read_csv(si22)
        d = e.merge(si, left_on="ModelID", right_on="DepMap_ID", how="left")
        d["TACSTD2"] = pd.to_numeric(d["TACSTD2"], errors="coerce")
        d["CLDN4"] = pd.to_numeric(d["CLDN4"], errors="coerce")
        d = d.dropna(subset=["TACSTD2", "CLDN4"])
        lung = d[d["lineage"] == "lung"]
        add_cohort(
            rows,
            "hist_22q2_lineage_lung",
            lung,
            "Historical: DepMap 22Q2 CCLE_expression.csv, sample_info lineage==lung "
            "(this older annotation includes mesothelioma).",
            ci=True,
        )
        add_cohort(
            rows,
            "hist_22q2_lineage_lung_no_meso",
            lung[lung["lineage_subtype"] != "mesothelioma"],
            "Historical: 22Q2 lineage==lung excluding mesothelioma.",
        )
        add_cohort(
            rows,
            "hist_22q2_NSCLC",
            lung[lung["lineage_subtype"] == "NSCLC"],
            "Historical: 22Q2 lineage_subtype==NSCLC.",
        )
        add_cohort(
            rows,
            "hist_22q2_all_lines",
            d,
            "Context only: 22Q2 all lineages. Not the B2 lung-only request.",
        )
        lung[
            [
                "ModelID",
                "cell_line_name",
                "CCLE_Name",
                "lineage",
                "lineage_subtype",
                "primary_disease",
                "TACSTD2",
                "CLDN4",
            ]
        ].to_csv(indir / "ccle22q2_lung_expression.csv", index=False)

    gct = indir / "ccle2018_rpkm_tacstd2_cldn4.csv"
    ann = Path("/tmp/ccle2018/Cell_lines_annotations_20181226.txt")
    if gct.exists():
        g = pd.read_csv(gct)
        g["TACSTD2_RPKM"] = pd.to_numeric(g["TACSTD2_RPKM"], errors="coerce")
        g["CLDN4_RPKM"] = pd.to_numeric(g["CLDN4_RPKM"], errors="coerce")
        g["tissue_suffix"] = g["CCLE_ID"].astype(str).str.split("_").str[-1]
        add_cohort(
            rows,
            "hist_2018_suffix_LUNG_RPKM",
            g[g["tissue_suffix"] == "LUNG"].dropna(subset=["TACSTD2_RPKM", "CLDN4_RPKM"]),
            "Historical: CCLE RNAseq RPKM 20180929, CCLE_ID suffix == LUNG. Spearman is rank-based so log transform does not change ρ.",
            ci=True,
        )
        if ann.exists():
            a = pd.read_csv(ann, sep="\t")
            m = g.merge(a, on="CCLE_ID", how="left")
            site = m[m["Site_Primary"].fillna("").str.lower() == "lung"].dropna(
                subset=["TACSTD2_RPKM", "CLDN4_RPKM"]
            )
            add_cohort(
                rows,
                "hist_2018_SitePrimary_lung_RPKM",
                site,
                "Historical: CCLE 2018 RPKM, annotation Site_Primary==lung.",
                ci=True,
            )
            site[["CCLE_ID", "Name", "Site_Primary", "Histology", "TACSTD2_RPKM", "CLDN4_RPKM"]].to_csv(
                indir / "ccle2018_lung_expression.csv", index=False
            )
        add_cohort(
            rows,
            "hist_2018_all_lines_RPKM",
            g.dropna(subset=["TACSTD2_RPKM", "CLDN4_RPKM"]),
            "Context only: CCLE 2018 RPKM all tissues. Not the B2 lung-only request.",
        )
    return rows


def _honest_statement(primary: dict) -> str:
    rho = primary["spearman_rho"]
    n = primary["n"]
    v = verdict(rho)
    r = primary["pearson_r"]
    ci = ""
    if primary.get("spearman_ci95_low") is not None:
        ci = (
            f" Bootstrap 95% CI [{primary['spearman_ci95_low']:.3f}, "
            f"{primary['spearman_ci95_high']:.3f}]; 0.69 "
            f"{'is' if primary.get('user_0.69_inside_spearman_ci95') else 'is not'} inside this CI."
        )
    if v["label"] == "MATCH_AT_2DP":
        return (
            f"Primary lung cell-line Spearman ρ = {rho:.4f} (n={n}), which rounds to 0.69 "
            f"and matches the user claim at 2 decimal places. Pearson r = {r:.4f}.{ci}"
        )
    if v["label"] == "NEARBY_NOT_EXACT":
        return (
            f"Primary lung cell-line Spearman ρ = {rho:.4f} (n={n}), within 0.05 of the "
            f"user claim 0.69 but does not match at 2 decimal places "
            f"(rounds to {v['rounded_2dp']:.2f}). Pearson r = {r:.4f}. "
            f"We do not treat this as confirmation of 0.69.{ci}"
        )
    return (
        f"Primary lung cell-line Spearman ρ = {rho:.4f} (n={n}), which does not match "
        f"the user claim 0.69 (absolute difference {v['abs_diff']:.4f}; rounds to "
        f"{v['rounded_2dp']:.2f}). Pearson r = {r:.4f}. The user number is not reproduced "
        f"on this public 24Q4 lung-cell-line slice.{ci}"
    )


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""New DepMap cut: CLDN4 Chronos (and CRISPR-overlap RNA) vs IFN / MHC-I / CD274.

Additive only. The all-lung CLDN4-RNA Spearman (n=214, ρ=+0.28 / +0.33 vs
Hallmark IFN-γ / CD274) is treated as already reported in
methods/depmap_cldn4_ifn and is not restated as the finding.
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

MHC1_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
AXES = ["sig_IFNG", "sig_MHC1", "CD274"]
AXIS_LABELS = {
    "sig_IFNG": "Hallmark IFN-γ (mean-z)",
    "sig_MHC1": "MHC-I mean-z (HLA-A/B/C + B2M)",
    "CD274": "CD274 log2(TPM+1)",
}


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index), []
    z = df[use].apply(zscore, axis=0)
    return z.mean(axis=1), use


def spearman(x: pd.Series, y: pd.Series) -> dict:
    a = pd.concat([x, y], axis=1).dropna()
    n = int(len(a))
    if n < 8:
        return {"n": n, "rho": None, "p": None}
    rho, p = stats.spearmanr(a.iloc[:, 0], a.iloc[:, 1])
    return {"n": n, "rho": float(rho), "p": float(p)}


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


def mannwhitney(high: pd.Series, low: pd.Series) -> dict:
    a = high.dropna()
    b = low.dropna()
    n_high, n_low = int(len(a)), int(len(b))
    if n_high < 5 or n_low < 5:
        return {
            "n_high": n_high,
            "n_low": n_low,
            "median_high": None,
            "median_low": None,
            "delta_median": None,
            "U": None,
            "p": None,
            "cliffs_delta": None,
        }
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Cliff's δ = (2U)/(n1 n2) − 1. Positive: Q4 tends larger than Q1.
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


def add_scores(frame: pd.DataFrame, ifng_genes: list[str]) -> tuple[pd.DataFrame, dict]:
    out = frame.copy()
    sig_ifng, used_ifng = signature(out, ifng_genes)
    sig_mhc, used_mhc = signature(out, MHC1_GENES)
    out["sig_IFNG"] = sig_ifng
    out["sig_MHC1"] = sig_mhc
    cov = {
        "ifng_genes_used": used_ifng,
        "ifng_n_used": len(used_ifng),
        "ifng_n_gmt": len(ifng_genes),
        "mhc1_genes_used": used_mhc,
        "mhc1_n_used": len(used_mhc),
    }
    return out, cov


def fmt_p(p: float | None) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3g}"


def fmt_rho(rho: float | None) -> str:
    if rho is None or (isinstance(rho, float) and not np.isfinite(rho)):
        return "NA"
    return f"{rho:+.3f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_cldn4_extra")
    ap.add_argument("--outdir", default="methods/depmap_cldn4_extra")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    model = pd.read_csv(data / "Model.csv")
    expr = pd.read_csv(data / "expression_ifn_panel.csv")
    crispr = pd.read_csv(data / "crispr_cldn4.csv")
    ifng_genes = [
        g.strip()
        for g in (data / "hallmark_ifng_genes.txt").read_text().splitlines()
        if g.strip()
    ]
    manifest = json.loads((data / "download_manifest.json").read_text())

    crispr = crispr.rename(columns={"CLDN4": "CLDN4_Chronos"})
    # Chronos files can store blanks as empty strings.
    crispr["CLDN4_Chronos"] = pd.to_numeric(crispr["CLDN4_Chronos"], errors="coerce")

    df = expr.merge(model, on="ModelID", how="left")
    lung_models = model[
        (model["OncotreeLineage"] == "Lung") & (model["ModelType"] == "Cell Line")
    ].copy()
    n_lung_models = int(len(lung_models))

    rna_need = ["CLDN4", "CD274"] + [g for g in MHC1_GENES if g in df.columns]
    lung_rna = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"] == "Cell Line")].copy()
    n_lung_rna_any = int(len(lung_rna))
    lung_rna = lung_rna.dropna(subset=[c for c in rna_need if c in lung_rna.columns]).copy()
    n_lung_rna = int(len(lung_rna))

    # New cut: lung cell lines with complete RNA *and* a finite CLDN4 Chronos.
    overlap = lung_rna.merge(crispr[["ModelID", "CLDN4_Chronos"]], on="ModelID", how="left")
    n_rna_no_crispr = int(overlap["CLDN4_Chronos"].isna().sum())
    overlap = overlap.dropna(subset=["CLDN4_Chronos"]).copy()
    overlap["group"] = overlap.apply(cohort_tag, axis=1)
    overlap["is_NSCLC"] = overlap["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    overlap["is_SCLC_NET"] = overlap["group"] == "SCLC_NET"

    # Signatures z-scored *within the CRISPR∩RNA cohort* (not the n=214 RNA-only set).
    overlap, cov_all = add_scores(overlap, ifng_genes)
    nsclc, cov_nsclc = add_scores(overlap[overlap["is_NSCLC"]].copy(), ifng_genes)
    sclc, _ = add_scores(overlap[overlap["is_SCLC_NET"]].copy(), ifng_genes)
    luad, _ = add_scores(overlap[overlap["group"] == "LUAD"].copy(), ifng_genes)
    lusc, _ = add_scores(overlap[overlap["group"] == "LUSC"].copy(), ifng_genes)

    cohorts = {
        "CRISPR_RNA_lung": overlap,
        "CRISPR_RNA_NSCLC": nsclc,
        "CRISPR_RNA_SCLC_NET": sclc,
        "CRISPR_RNA_LUAD": luad,
        "CRISPR_RNA_LUSC": lusc,
    }

    corr_rows = []
    for cname, frame in cohorts.items():
        for predictor, pred_col in [
            ("CLDN4_Chronos", "CLDN4_Chronos"),
            ("CLDN4_RNA", "CLDN4"),
        ]:
            for target in AXES:
                rec = spearman(frame[pred_col], frame[target])
                rec.update(
                    {
                        "cut": cname,
                        "predictor": predictor,
                        "target": target,
                    }
                )
                do_ci = cname in {"CRISPR_RNA_lung", "CRISPR_RNA_NSCLC"} and rec["n"] and rec["n"] >= 8
                if do_ci:
                    a = pd.concat([frame[pred_col], frame[target]], axis=1).dropna()
                    lo, hi = bootstrap_spearman_ci(
                        a.iloc[:, 0].to_numpy(), a.iloc[:, 1].to_numpy()
                    )
                    rec["ci95_low"] = lo
                    rec["ci95_high"] = hi
                else:
                    rec["ci95_low"] = None
                    rec["ci95_high"] = None
                corr_rows.append(rec)

    corr = pd.DataFrame(corr_rows)
    # BH within the three new axes, separately per cut × predictor family.
    for cname in corr["cut"].unique():
        for predictor in ["CLDN4_Chronos", "CLDN4_RNA"]:
            mask = (corr["cut"] == cname) & (corr["predictor"] == predictor)
            p = corr.loc[mask, "p"]
            if p.notna().sum() == 0:
                continue
            q = pd.Series(index=p.index, dtype=float)
            finite = p.notna()
            q.loc[finite] = bh_fdr(p.loc[finite].astype(float).tolist())
            corr.loc[mask, "q_bh"] = q

    corr.to_csv(tabdir / "extra_correlations.tsv", sep="\t", index=False)

    # Q4 vs Q1 of CLDN4 RNA on the CRISPR∩RNA lung set (new statistic, new n).
    q_rows = []
    q1 = overlap["CLDN4"].quantile(0.25)
    q3 = overlap["CLDN4"].quantile(0.75)
    high = overlap[overlap["CLDN4"] >= q3]
    low = overlap[overlap["CLDN4"] <= q1]
    for target in AXES:
        rec = mannwhitney(high[target], low[target])
        rec.update(
            {
                "cut": "CRISPR_RNA_lung_CLDN4_Q4_vs_Q1",
                "predictor": "CLDN4_RNA_Q4_vs_Q1",
                "target": target,
                "q1_threshold": float(q1),
                "q3_threshold": float(q3),
            }
        )
        q_rows.append(rec)
    qtab = pd.DataFrame(q_rows)
    p = qtab["p"].astype(float).tolist()
    qtab["q_bh"] = bh_fdr(p)
    qtab.to_csv(tabdir / "extra_q4q1.tsv", sep="\t", index=False)

    # Chronos QC on the new-cut cohort.
    ch = overlap["CLDN4_Chronos"]
    chronos_qc = {
        "n": int(ch.notna().sum()),
        "median": float(ch.median()),
        "mean": float(ch.mean()),
        "q25": float(ch.quantile(0.25)),
        "q75": float(ch.quantile(0.75)),
        "min": float(ch.min()),
        "max": float(ch.max()),
        "n_chronos_lt_-0.5": int((ch < -0.5).sum()),
        "n_chronos_lt_-1.0": int((ch < -1.0).sum()),
        "frac_chronos_lt_-0.5": float((ch < -0.5).mean()),
        "note": "CLDN4 is not a common essential; Chronos near 0 is expected.",
    }

    sample_cols = [
        "ModelID",
        "CellLineName",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "group",
        "CLDN4",
        "CLDN4_Chronos",
        "CD274",
        "sig_IFNG",
        "sig_MHC1",
    ]
    overlap[sample_cols].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "crispr_rna_lung_lines.tsv", sep="\t", index=False
    )

    counts = (
        overlap.groupby("group", dropna=False)
        .size()
        .rename("n")
        .reset_index()
        .sort_values("n", ascending=False)
    )
    counts.to_csv(tabdir / "cohort_counts.tsv", sep="\t", index=False)

    n_funnel = {
        "n_lung_models_Model_csv": n_lung_models,
        "n_lung_cell_lines_in_RNA_matrix": n_lung_rna_any,
        "n_lung_RNA_complete_CLDN4_CD274_MHC1": n_lung_rna,
        "n_RNA_complete_without_CLDN4_Chronos": n_rna_no_crispr,
        "n_CRISPR_RNA_overlap": int(len(overlap)),
        "n_CRISPR_RNA_NSCLC": int(len(nsclc)),
        "n_CRISPR_RNA_LUAD": int(len(luad)),
        "n_CRISPR_RNA_LUSC": int(len(lusc)),
        "n_CRISPR_RNA_SCLC_NET": int(len(sclc)),
        "n_CRISPR_RNA_other_NSCLC": int((overlap["group"] == "other_NSCLC").sum()),
        "n_CRISPR_RNA_other_lung": int((overlap["group"] == "other_lung").sum()),
    }

    def pick(cname: str, predictor: str, target: str) -> dict:
        r = corr[
            (corr["cut"] == cname)
            & (corr["predictor"] == predictor)
            & (corr["target"] == target)
        ].iloc[0]
        return {
            "n": int(r["n"]),
            "rho": None if pd.isna(r["rho"]) else float(r["rho"]),
            "p": None if pd.isna(r["p"]) else float(r["p"]),
            "q_bh": None if "q_bh" not in r or pd.isna(r["q_bh"]) else float(r["q_bh"]),
            "ci95_low": None if pd.isna(r["ci95_low"]) else float(r["ci95_low"]),
            "ci95_high": None if pd.isna(r["ci95_high"]) else float(r["ci95_high"]),
        }

    # Extra table (the PR deliverable): Chronos vs the three axes, plus
    # CRISPR-overlap RNA (new n, not the reported n=214).
    extra_rows = []
    for cname in ["CRISPR_RNA_lung", "CRISPR_RNA_NSCLC"]:
        for predictor in ["CLDN4_Chronos", "CLDN4_RNA"]:
            for target in AXES:
                rec = pick(cname, predictor, target)
                extra_rows.append(
                    {
                        "cut": cname,
                        "predictor": predictor,
                        "target": target,
                        "n": rec["n"],
                        "rho": rec["rho"],
                        "p": rec["p"],
                        "q_bh": rec["q_bh"],
                        "ci95_low": rec["ci95_low"],
                        "ci95_high": rec["ci95_high"],
                    }
                )
    extra = pd.DataFrame(extra_rows)
    extra.to_csv(tabdir / "extra_table.tsv", sep="\t", index=False)

    # Markdown extra table for FINDING.md assembly.
    md_lines = [
        "| cut | predictor | target | n | ρ | p | q | 95% CI |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    for _, r in extra.iterrows():
        ci = "NA"
        if r["ci95_low"] is not None and pd.notna(r["ci95_low"]):
            ci = f"[{r['ci95_low']:+.2f}, {r['ci95_high']:+.2f}]"
        md_lines.append(
            f"| {r['cut']} | {r['predictor']} | {r['target']} | {int(r['n'])} | "
            f"{fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {fmt_p(r['q_bh'])} | {ci} |"
        )
    (tabdir / "extra_table.md").write_text("\n".join(md_lines) + "\n")

    key = {
        "release": manifest.get("release"),
        "doi": manifest.get("doi"),
        "already_reported_not_this_cut": {
            "source": "methods/depmap_cldn4_ifn",
            "n": 214,
            "CLDN4_RNA_vs_Hallmark_IFNG_rho": 0.284,
            "CLDN4_RNA_vs_CD274_rho": 0.330,
            "note": "Not recomputed as the finding. This page is the CRISPR∩RNA cut.",
        },
        "n_funnel": n_funnel,
        "chronos_qc": chronos_qc,
        "signature_coverage_CRISPR_RNA_lung": cov_all,
        "signature_coverage_CRISPR_RNA_NSCLC": cov_nsclc,
        "correlations": {
            cname: {
                pred: {t: pick(cname, pred, t) for t in AXES}
                for pred in ["CLDN4_Chronos", "CLDN4_RNA"]
            }
            for cname in cohorts
        },
        "q4q1": qtab.to_dict(orient="records"),
        "notes": [
            "New cut = lung cell lines with complete RNA and a finite CLDN4 Chronos.",
            "Signatures are mean of within-cut gene-wise z-scores.",
            "More negative Chronos = stronger CLDN4 dependency.",
            "No immune infiltrate and no IFN treatment.",
        ],
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")

    color = {
        "LUAD": "#1f77b4",
        "LUSC": "#ff7f0e",
        "other_NSCLC": "#2ca02c",
        "SCLC_NET": "#d62728",
        "other_lung": "#7f7f7f",
    }

    # Primary extra scatter: Chronos vs the three axes.
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0), constrained_layout=True)
    panels = [
        ("sig_IFNG", "Hallmark IFN-γ score (mean-z)"),
        ("sig_MHC1", "MHC-I score (HLA-A/B/C + B2M, mean-z)"),
        ("CD274", "CD274  log2(TPM+1)"),
    ]
    for ax, (col, ylab) in zip(axes, panels):
        for g, sub in overlap.groupby("group"):
            ax.scatter(
                sub["CLDN4_Chronos"],
                sub[col],
                s=18,
                alpha=0.75,
                c=color.get(g, "#333333"),
                label=f"{g} n={len(sub)}",
                edgecolors="none",
            )
        rec = pick("CRISPR_RNA_lung", "CLDN4_Chronos", col)
        ax.set_xlabel("CLDN4 Chronos (more negative = more dependent)")
        ax.set_ylabel(ylab)
        ax.set_title(
            f"CRISPR∩RNA lung n={rec['n']}\n"
            f"ρ={rec['rho']:.2f}  p={fmt_p(rec['p'])}",
            fontsize=9,
        )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.08),
    )
    fig.suptitle(
        "DepMap 24Q4 extra cut — CLDN4 Chronos vs IFN / MHC-I / CD274",
        fontsize=11,
        y=1.08,
    )
    fig.savefig(figdir / "fig_cldn4_chronos_vs_ifn_mhc_cd274.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_cldn4_chronos_vs_ifn_mhc_cd274.pdf", bbox_inches="tight")
    plt.close(fig)

    print(
        json.dumps(
            {
                "n_funnel": n_funnel,
                "chronos_qc": chronos_qc,
                "CLDN4_Chronos_vs_IFNG": pick("CRISPR_RNA_lung", "CLDN4_Chronos", "sig_IFNG"),
                "CLDN4_Chronos_vs_MHC1": pick("CRISPR_RNA_lung", "CLDN4_Chronos", "sig_MHC1"),
                "CLDN4_Chronos_vs_CD274": pick("CRISPR_RNA_lung", "CLDN4_Chronos", "CD274"),
                "CLDN4_RNA_overlap_vs_IFNG": pick("CRISPR_RNA_lung", "CLDN4_RNA", "sig_IFNG"),
                "CLDN4_RNA_overlap_vs_CD274": pick("CRISPR_RNA_lung", "CLDN4_RNA", "CD274"),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

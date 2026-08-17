#!/usr/bin/env python3
"""Score CLDN4 vs Hallmark IFN-γ, MHC-I mean, CD274, TACSTD2 in DepMap lung lines.

Additive: TACSTD2–CLDN4 co-expression is treated as already known and is not
the finding. The new axis is IFN / MHC-I / CD274 (cancer-cell basal RNA).
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
TARGETS = ["sig_IFNG", "sig_MHC1", "CD274", "TACSTD2"]
TARGET_LABELS = {
    "sig_IFNG": "Hallmark IFN-γ (mean-z)",
    "sig_MHC1": "MHC-I mean-z (HLA-A/B/C + B2M)",
    "CD274": "CD274 log2(TPM+1)",
    "TACSTD2": "TACSTD2 log2(TPM+1)",
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


def bootstrap_spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 5000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


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
        # rank from the end: i is original index of the largest remaining p
        k = n - rank + 1
        val = min(prev, p[i] * n / k)
        q[i] = val
        prev = val
    return [float(min(1.0, v)) for v in q]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_cldn4_ifn")
    ap.add_argument("--outdir", default="methods/depmap_cldn4_ifn")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    model = pd.read_csv(data / "Model.csv")
    expr = pd.read_csv(data / "expression_ifn_panel.csv")
    ifng_genes = [g.strip() for g in (data / "hallmark_ifng_genes.txt").read_text().splitlines() if g.strip()]
    manifest = json.loads((data / "download_manifest.json").read_text())

    df = expr.merge(model, on="ModelID", how="left")
    lung = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"] == "Cell Line")].copy()
    # Complete cases for CLDN4 (required) and the four scored partners.
    need = ["CLDN4", "TACSTD2", "CD274"] + [g for g in MHC1_GENES if g in lung.columns]
    lung = lung.dropna(subset=[c for c in need if c in lung.columns]).copy()
    lung["group"] = lung.apply(cohort_tag, axis=1)
    lung["is_NSCLC"] = lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    lung["is_SCLC_NET"] = lung["group"] == "SCLC_NET"

    # Signatures are z-scored *within the analysis cohort* so they are not
    # contaminated by non-lung lines. Primary scores use all-lung; NSCLC
    # scores are recomputed within NSCLC for the NSCLC-only table.
    def add_scores(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
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

    lung, cov_all = add_scores(lung)
    nsclc, cov_nsclc = add_scores(lung[lung["is_NSCLC"]].copy())
    sclc, _ = add_scores(lung[lung["is_SCLC_NET"]].copy())
    luad, _ = add_scores(lung[lung["group"] == "LUAD"].copy())
    lusc, _ = add_scores(lung[lung["group"] == "LUSC"].copy())

    cohorts = {
        "lung_cell_lines": lung,
        "NSCLC": nsclc,
        "SCLC_NET": sclc,
        "LUAD": luad,
        "LUSC": lusc,
    }

    rows = []
    for cname, frame in cohorts.items():
        for target in TARGETS:
            rec = spearman(frame["CLDN4"], frame[target])
            rec.update({"cohort": cname, "predictor": "CLDN4", "target": target})
            if cname in {"lung_cell_lines", "NSCLC"} and rec["n"] and rec["n"] >= 8:
                a = pd.concat([frame["CLDN4"], frame[target]], axis=1).dropna()
                lo, hi = bootstrap_spearman_ci(a.iloc[:, 0].to_numpy(), a.iloc[:, 1].to_numpy())
                rec["ci95_low"] = lo
                rec["ci95_high"] = hi
            else:
                rec["ci95_low"] = None
                rec["ci95_high"] = None
            rows.append(rec)

    corr = pd.DataFrame(rows)
    # BH within the three *new* axes on the primary (all-lung) and NSCLC tables.
    # TACSTD2 is given and is excluded from FDR.
    for cname in ["lung_cell_lines", "NSCLC"]:
        mask = (corr["cohort"] == cname) & (corr["target"].isin(["sig_IFNG", "sig_MHC1", "CD274"]))
        p = corr.loc[mask, "p"].astype(float).tolist()
        corr.loc[mask, "q_bh"] = bh_fdr(p)

    corr.to_csv(tabdir / "correlations.tsv", sep="\t", index=False)

    sample_cols = [
        "ModelID",
        "CellLineName",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "group",
        "CLDN4",
        "TACSTD2",
        "CD274",
        "sig_IFNG",
        "sig_MHC1",
    ]
    lung[sample_cols].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "lung_cell_lines.tsv", sep="\t", index=False
    )

    counts = (
        lung.groupby("group", dropna=False)
        .size()
        .rename("n")
        .reset_index()
        .sort_values("n", ascending=False)
    )
    counts.to_csv(tabdir / "cohort_counts.tsv", sep="\t", index=False)

    def pick(cname: str, target: str) -> dict:
        r = corr[(corr["cohort"] == cname) & (corr["target"] == target)].iloc[0]
        return {
            "n": int(r["n"]),
            "rho": None if pd.isna(r["rho"]) else float(r["rho"]),
            "p": None if pd.isna(r["p"]) else float(r["p"]),
            "q_bh": None if "q_bh" not in r or pd.isna(r["q_bh"]) else float(r["q_bh"]),
            "ci95_low": None if pd.isna(r["ci95_low"]) else float(r["ci95_low"]),
            "ci95_high": None if pd.isna(r["ci95_high"]) else float(r["ci95_high"]),
        }

    key = {
        "release": manifest.get("release"),
        "doi": manifest.get("doi"),
        "expression": "OmicsExpressionProteinCodingGenesTPMLogp1.csv = log2(TPM+1)",
        "lineage": "OncotreeLineage == Lung AND ModelType == Cell Line",
        "n_lung_cell_lines": int(len(lung)),
        "n_NSCLC": int(len(nsclc)),
        "n_SCLC_NET": int(len(sclc)),
        "n_LUAD": int(len(luad)),
        "n_LUSC": int(len(lusc)),
        "n_other_lung": int((lung["group"] == "other_lung").sum()),
        "n_other_NSCLC": int((lung["group"] == "other_NSCLC").sum()),
        "signature_coverage_all_lung": cov_all,
        "signature_coverage_NSCLC": cov_nsclc,
        "CLDN4_vs": {
            cname: {t: pick(cname, t) for t in TARGETS} for cname in cohorts
        },
        "notes": [
            "TACSTD2–CLDN4 is a given companion correlation; not the finding.",
            "Hallmark IFN-γ and MHC-I are mean of within-cohort gene-wise z-scores.",
            "No immune infiltrate: scores are cancer-cell basal transcription.",
            "Primary n is complete lung cell lines with CLDN4 + scored genes.",
        ],
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2) + "\n")

    # --- scatter: CLDN4 vs the three additive axes (IFN / MHC / CD274) ---
    color = {
        "LUAD": "#1f77b4",
        "LUSC": "#ff7f0e",
        "other_NSCLC": "#2ca02c",
        "SCLC_NET": "#d62728",
        "other_lung": "#7f7f7f",
    }
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0), constrained_layout=True)
    panels = [
        ("sig_IFNG", "Hallmark IFN-γ score (mean-z)"),
        ("sig_MHC1", "MHC-I score (HLA-A/B/C + B2M, mean-z)"),
        ("CD274", "CD274  log2(TPM+1)"),
    ]
    for ax, (col, ylab) in zip(axes, panels):
        for g, sub in lung.groupby("group"):
            ax.scatter(
                sub["CLDN4"],
                sub[col],
                s=18,
                alpha=0.75,
                c=color.get(g, "#333333"),
                label=f"{g} n={len(sub)}",
                edgecolors="none",
            )
        rec = pick("lung_cell_lines", col)
        ax.set_xlabel("CLDN4  log2(TPM+1)")
        ax.set_ylabel(ylab)
        ax.set_title(
            f"all lung n={rec['n']}\n"
            f"ρ={rec['rho']:.2f}  p={rec['p']:.1e}\n"
            f"NSCLC n={pick('NSCLC', col)['n']}  ρ={pick('NSCLC', col)['rho']:.2f}",
            fontsize=9,
        )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle(
        "DepMap 24Q4 lung cell lines — CLDN4 vs IFN / MHC-I / CD274 (basal RNA)",
        fontsize=11,
        y=1.08,
    )
    fig.savefig(figdir / "fig_cldn4_vs_ifn_mhc_cd274.png", dpi=160, bbox_inches="tight")
    fig.savefig(figdir / "fig_cldn4_vs_ifn_mhc_cd274.pdf", bbox_inches="tight")
    plt.close(fig)

    # Single primary scatter (CLDN4 vs Hallmark IFN-γ) for the "a scatter" deliverable.
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for g, sub in lung.groupby("group"):
        ax.scatter(
            sub["CLDN4"],
            sub["sig_IFNG"],
            s=22,
            alpha=0.8,
            c=color.get(g, "#333333"),
            label=f"{g} n={len(sub)}",
            edgecolors="none",
        )
    rec = pick("lung_cell_lines", "sig_IFNG")
    ns = pick("NSCLC", "sig_IFNG")
    ax.set_xlabel("CLDN4  log2(TPM+1)")
    ax.set_ylabel("Hallmark IFN-γ response (within-lung mean-z)")
    ax.set_title(
        "DepMap 24Q4 lung cell lines\n"
        f"CLDN4 vs Hallmark IFN-γ   all-lung n={rec['n']}  ρ={rec['rho']:.2f}  p={rec['p']:.1e}\n"
        f"NSCLC n={ns['n']}  ρ={ns['rho']:.2f}  p={ns['p']:.1e}"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fig_cldn4_vs_ifng.png", dpi=160)
    fig.savefig(figdir / "fig_cldn4_vs_ifng.pdf")
    plt.close(fig)

    print(json.dumps({
        "n_lung": key["n_lung_cell_lines"],
        "n_NSCLC": key["n_NSCLC"],
        "n_SCLC_NET": key["n_SCLC_NET"],
        "CLDN4_vs_IFNG_lung": key["CLDN4_vs"]["lung_cell_lines"]["sig_IFNG"],
        "CLDN4_vs_MHC1_lung": key["CLDN4_vs"]["lung_cell_lines"]["sig_MHC1"],
        "CLDN4_vs_CD274_lung": key["CLDN4_vs"]["lung_cell_lines"]["CD274"],
        "CLDN4_vs_IFNG_NSCLC": key["CLDN4_vs"]["NSCLC"]["sig_IFNG"],
        "CLDN4_vs_MHC1_NSCLC": key["CLDN4_vs"]["NSCLC"]["sig_MHC1"],
        "CLDN4_vs_CD274_NSCLC": key["CLDN4_vs"]["NSCLC"]["CD274"],
        "CLDN4_vs_TACSTD2_lung": key["CLDN4_vs"]["lung_cell_lines"]["TACSTD2"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

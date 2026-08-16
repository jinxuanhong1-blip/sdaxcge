#!/usr/bin/env python3
"""
Claim A2 verification: TACSTD2 (TROP2) vs immune score in durvalumab-treated NSCLC.

User claim: Spearman rho = -0.65 (raw), -0.46 (tumor-purity-adjusted), p = 2e-4.

Data (only public human durvalumab NSCLC bulk RNA-seq we could locate):
  - GSE253564: pre-treatment tumors, randomized phase II neoadjuvant durvalumab
    +/- SBRT trial (Altorki et al.), n=32, FPKM.
  - GSE248378: post-treatment resected tumors from the same trial, n=29, FPKM.

Immune score / purity: exact Python port of the ESTIMATE R package v1.0.13
(rank-normalized ssGSEA, weight exponent 0.25; purity from the published
Affymetrix-calibrated formula cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)).
ESTIMATE depends only on within-sample ranks, so FPKM vs log-FPKM is irrelevant.

Partial Spearman = Pearson partial correlation on ranks; p from t with n-3 df.
"""

import json
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "results", "claim_A2")
os.makedirs(OUT, exist_ok=True)

CLAIM = {"rho_raw": -0.65, "rho_purity_adj": -0.46, "p": 2e-4}


def load_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def estimate_scores(expr, gene_sets):
    """Exact port of estimateScore() from the ESTIMATE R package (v1.0.13).

    expr: DataFrame genes x samples (any monotone-equivalent scale).
    Returns DataFrame samples x [StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity].
    """
    genes = expr.index.to_numpy()
    n_genes = expr.shape[0]
    # Per-sample rank normalization (average ties), scaled to 10000*rank/Ng
    m = expr.rank(axis=0, method="average").to_numpy() * (10000.0 / n_genes)

    rows = {}
    overlaps = {}
    for name in ("StromalSignature", "ImmuneSignature"):
        gs = set(gene_sets[name])
        in_set = np.isin(genes, list(gs))
        overlaps[name] = int(in_set.sum())
        es = np.empty(expr.shape[1])
        for j in range(expr.shape[1]):
            order = np.argsort(-m[:, j], kind="stable")
            tag = in_set[order].astype(float)
            w = np.abs(m[order, j]) ** 0.25
            nh = tag.sum()
            nm = n_genes - nh
            pn = np.cumsum(tag * w) / (tag * w).sum()
            p0 = np.cumsum((1.0 - tag) / nm)
            es[j] = float(np.sum(pn - p0))  # ES = sum(RES), as in the R code
        rows["Stromal" if name.startswith("Str") else "Immune"] = es

    df = pd.DataFrame(
        {"StromalScore": rows["Stromal"], "ImmuneScore": rows["Immune"]},
        index=expr.columns,
    )
    df["ESTIMATEScore"] = df["StromalScore"] + df["ImmuneScore"]
    purity = np.cos(0.6049872018 + 0.0001467884 * df["ESTIMATEScore"].to_numpy())
    purity[purity < 0] = np.nan
    df["TumorPurity"] = purity
    return df, overlaps


def spearman(x, y):
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def partial_spearman(x, y, z):
    """Partial Spearman rho of x,y controlling z; p via t with n-3 df."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    r = (rxy - rxz * ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2))
    n = len(x)
    dof = n - 3
    t = r * np.sqrt(dof / (1 - r**2))
    p = 2 * stats.t.sf(abs(t), dof)
    return float(r), float(p)


def analyze(tag, expr, gene_sets):
    est, overlaps = estimate_scores(expr, gene_sets)
    if "TACSTD2" not in expr.index:
        raise KeyError(f"TACSTD2 missing from {tag}")
    tac = np.log2(expr.loc["TACSTD2"].to_numpy(dtype=float) + 1.0)
    # ESTIMATE TumorPurity = cos(a + b*ESTIMATEScore) is strictly monotone
    # decreasing over the valid range, and the Affymetrix-calibrated formula
    # goes out of bounds (NaN) for many FPKM samples here. Partial Spearman is
    # rank-based and invariant to monotone transforms of the covariate, so
    # controlling for ESTIMATEScore is exactly the purity adjustment while
    # keeping all samples.
    purity_covariate = est["ESTIMATEScore"].to_numpy()
    n = len(tac)

    out = {
        "dataset": tag,
        "n_samples": n,
        "geneset_overlap": overlaps,
        "n_purity_formula_out_of_bounds": int(np.isnan(est["TumorPurity"]).sum()),
        "purity_adjustment_note": (
            "partial Spearman controls for ESTIMATEScore ranks, which are "
            "rank-identical to ESTIMATE TumorPurity (monotone transform); "
            "avoids dropping samples where the Affymetrix purity formula is "
            "out of bounds on FPKM data"
        ),
        "metrics": {},
    }

    immune_defs = {"ESTIMATE_ImmuneScore": est["ImmuneScore"].to_numpy()}
    for g in ("PTPRC", "CD8A"):
        if g in expr.index:
            immune_defs[f"log2_{g}"] = np.log2(expr.loc[g].to_numpy(dtype=float) + 1.0)
    if {"GZMA", "PRF1"}.issubset(expr.index):
        immune_defs["CYT_gzma_prf1"] = 0.5 * (
            np.log2(expr.loc["GZMA"].to_numpy(dtype=float) + 1.0)
            + np.log2(expr.loc["PRF1"].to_numpy(dtype=float) + 1.0)
        )

    for name, vec in immune_defs.items():
        rho, p = spearman(tac, vec)
        pr, pp = partial_spearman(tac, vec, purity_covariate)
        out["metrics"][name] = {
            "spearman_rho_raw": round(rho, 4),
            "p_raw": float(f"{p:.3e}"),
            "spearman_rho_purity_adjusted": round(pr, 4),
            "p_purity_adjusted": float(f"{pp:.3e}"),
        }

    # per-sample table
    tbl = est.copy()
    tbl.insert(0, "TACSTD2_log2FPKMp1", tac)
    tbl.to_csv(os.path.join(OUT, f"scores_{tag}.csv"))
    return out, est, tac


def main():
    gene_sets = load_gmt(os.path.join(DATA, "estimate", "inst", "extdata", "SI_geneset.gmt"))

    pre = pd.read_csv(os.path.join(DATA, "GSE253564_Pre_FPKMs.txt.gz"), sep="\t")
    pre = pre.drop(columns=["Entrez.ID"]).groupby("gene").max()

    post = pd.read_csv(os.path.join(DATA, "GSE248378_Durva_Post_FPKMs.txt.gz"), sep="\t")
    post = post.groupby("gene").max()

    results = {"claim": CLAIM, "datasets": []}
    plots = []
    for tag, expr in (("GSE253564_pretreatment", pre), ("GSE248378_posttreatment", post)):
        res, est, tac = analyze(tag, expr, gene_sets)
        results["datasets"].append(res)
        plots.append((tag, est, tac))
        m = res["metrics"]["ESTIMATE_ImmuneScore"]
        print(
            f"{tag} (n={res['n_samples']}): raw rho={m['spearman_rho_raw']} "
            f"(p={m['p_raw']:.2e}); purity-adj rho={m['spearman_rho_purity_adjusted']} "
            f"(p={m['p_purity_adjusted']:.2e})"
        )

    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for ax, (tag, est, tac) in zip(axes, plots):
            imm = est["ImmuneScore"].to_numpy()
            rho, p = spearman(tac, imm)
            sc = ax.scatter(
                tac, imm, c=est["ESTIMATEScore"], cmap="viridis_r",
                s=45, edgecolor="k", lw=0.4,
            )
            ax.set_xlabel("TACSTD2 log2(FPKM+1)")
            ax.set_ylabel("ESTIMATE ImmuneScore")
            ax.set_title(f"{tag} (n={len(tac)})\nSpearman rho={rho:.2f}, p={p:.1e}")
            fig.colorbar(sc, ax=ax, label="ESTIMATEScore (purity rank proxy)")
        fig.suptitle("Claim A2: TACSTD2 vs immune score, durvalumab NSCLC trial (GSE253564/GSE248378)")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "tacstd2_vs_immune.png"), dpi=150)
        print("wrote figure")
    except ImportError:
        print("matplotlib unavailable; skipping figure")


if __name__ == "__main__":
    main()

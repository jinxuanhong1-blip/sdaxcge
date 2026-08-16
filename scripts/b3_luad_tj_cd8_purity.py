#!/usr/bin/env python3
"""
B3 - TCGA-LUAD: Tight-Junction (TJ) score vs CD8 T-cell score, before and after
adjusting for ESTIMATE-derived tumor purity.

Scientific question
-------------------
Epithelial tight-junction gene expression (TJ) and CD8 T-cell abundance are both
strongly confounded by tumor purity: a sample with more tumor epithelium tends to
have a higher TJ signal and a lower relative immune signal simply because of its
cellular composition. This script asks whether any TJ-vs-CD8 association survives
after controlling for ESTIMATE tumor purity (a partial Spearman correlation).

The point of this analysis is to report the association *honestly*: both the raw
correlation and the purity-adjusted one, with sample sizes, p-values and bootstrap
confidence intervals, and with the exact gene sets and formulas written out so the
result can be reproduced and criticized.

Inputs
------
data/LUAD_HiSeqV2.gz          UCSC Xena TCGA-LUAD RNA-seq (IlluminaHiSeq, gene-level
                              log2(norm_count+1)). Rows = genes, cols = samples.
data/LUAD_ESTIMATE_RNAseqV2.txt
                              MD Anderson precomputed ESTIMATE scores for TCGA-LUAD
                              RNAseqV2 (Stromal_score, Immune_score, ESTIMATE_score).

Outputs (results/w200/B3_purity/)
---------------------------------
sample_scores.csv             per-sample TJ score, CD8 score, purity, ESTIMATE scores
tj_genes_used.csv             TJ signature genes and whether each was found
cd8_genes_used.csv            CD8 signature genes and whether each was found
correlation_results.json      all correlation statistics (machine readable)
correlation_results.csv       the same statistics as a flat table
scatter_tj_cd8.png            TJ vs CD8 scatter, colored by purity
scatter_partial_residuals.png purity-residualized TJ vs CD8 scatter
purity_associations.png       TJ~purity and CD8~purity scatters
README.md                     methods + results narrative
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
OUT = REPO / "results" / "w200" / "B3_purity"
OUT.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20260816
N_BOOT = 2000

# ---------------------------------------------------------------------------
# Gene sets (curated, explicit, and written to disk so they can be audited).
# ---------------------------------------------------------------------------
# TJ = core structural tight-junction components (claudins that form epithelial
# TJ strands, occludin, the ZO scaffold proteins, JAM adhesion molecules, the
# tricellular-junction MARVEL proteins, and cingulin). This is a structural set
# deliberately kept away from broad signaling genes so the score reflects
# epithelial tight-junction content.
TJ_GENES = [
    "CLDN1", "CLDN3", "CLDN4", "CLDN7",   # epithelial claudins
    "OCLN",                                 # occludin
    "TJP1", "TJP2", "TJP3",                 # ZO-1/2/3 scaffolds
    "F11R", "JAM2", "JAM3",                 # JAM-A/B/C adhesion molecules
    "MARVELD2", "MARVELD3",                 # tricellulin / MarvelD3
    "CGN", "CGNL1",                         # cingulin / cingulin-like
]

# CD8 = canonical CD8 T-cell coreceptor chains. Kept intentionally minimal and
# lineage-specific (CD8A/CD8B) so the readout is CD8 T-cell abundance rather than
# a broad cytotoxicity program.
CD8_GENES = ["CD8A", "CD8B"]


def load_expression() -> pd.DataFrame:
    """Return genes x samples expression matrix (log2 norm_count+1)."""
    df = pd.read_csv(DATA / "LUAD_HiSeqV2.gz", sep="\t", index_col=0)
    df.index.name = "gene"
    return df


def load_estimate() -> pd.DataFrame:
    est = pd.read_csv(DATA / "LUAD_ESTIMATE_RNAseqV2.txt", sep="\t")
    est = est.rename(columns={"ID": "sample"}).set_index("sample")
    return est


def estimate_purity(estimate_score: pd.Series) -> pd.Series:
    """ESTIMATE tumor purity (Yoshihara et al. 2013, Nat Commun).

    purity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)

    The formula was calibrated on Affymetrix data; applying it to RNAseqV2
    ESTIMATE scores is common practice but approximate (see README caveats).
    Values are only defined where the cosine argument lies in [0, pi]; outside
    that range purity is set to NaN.
    """
    arg = 0.6049872018 + 0.0001467884 * estimate_score
    purity = np.cos(arg)
    purity[(arg < 0) | (arg > np.pi)] = np.nan
    return purity


def zscore_rows(sub: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene (row) across samples."""
    mu = sub.mean(axis=1)
    sd = sub.std(axis=1, ddof=0).replace(0, np.nan)
    return sub.sub(mu, axis=0).div(sd, axis=0)


def signature_score(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, pd.DataFrame]:
    """Mean z-score signature. Returns (score per sample, gene-presence table)."""
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    z = zscore_rows(expr.loc[present])
    score = z.mean(axis=0)
    table = pd.DataFrame(
        {"gene": genes, "found": [g in expr.index for g in genes]}
    )
    return score, table, present, missing


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    """Partial Spearman correlation of x and y controlling for z.

    Rank-transform all three variables, regress the ranks of x and y on the
    ranks of z (with intercept) via least squares, and Pearson-correlate the
    residuals. This is the standard rank-based partial correlation. The p-value
    uses the Student-t approximation with n-3 degrees of freedom.
    """
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    Z = np.column_stack([np.ones_like(rz), rz])

    def resid(v):
        beta, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ beta

    ex, ey = resid(rx), resid(ry)
    r, _ = stats.pearsonr(ex, ey)
    n = len(x)
    dof = n - 3
    if dof > 0 and abs(r) < 1.0:
        t = r * np.sqrt(dof / (1 - r ** 2))
        p = 2 * stats.t.sf(abs(t), dof)
    else:
        p = np.nan
    return r, p, n, ex, ey


def boot_ci(func, *arrays, n_boot=N_BOOT, seed=RNG_SEED):
    """Percentile bootstrap 95% CI for a statistic returned by func(*arrays)."""
    rng = np.random.default_rng(seed)
    n = len(arrays[0])
    stats_out = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            stats_out.append(func(*[a[idx] for a in arrays]))
        except Exception:
            continue
    lo, hi = np.nanpercentile(stats_out, [2.5, 97.5])
    return float(lo), float(hi)


def main():
    expr = load_expression()
    est = load_estimate()

    # Match samples between expression and ESTIMATE on the 15-char TCGA barcode
    # (patient + sample-type), and keep only primary solid tumors (code -01).
    expr_cols = {c[:15]: c for c in expr.columns}
    common = [s for s in est.index if s in expr_cols and s.endswith("-01")]
    common = sorted(set(common))

    expr_m = expr[[expr_cols[s] for s in common]]
    expr_m.columns = common
    est_m = est.loc[common]

    tj_score, tj_tbl, tj_present, tj_missing = signature_score(expr_m, TJ_GENES)
    cd8_score, cd8_tbl, cd8_present, cd8_missing = signature_score(expr_m, CD8_GENES)

    purity = estimate_purity(est_m["ESTIMATE_score"])

    scores = pd.DataFrame({
        "sample": common,
        "TJ_score": tj_score.values,
        "CD8_score": cd8_score.values,
        "purity": purity.values,
        "Immune_score": est_m["Immune_score"].values,
        "Stromal_score": est_m["Stromal_score"].values,
        "ESTIMATE_score": est_m["ESTIMATE_score"].values,
    }).set_index("sample")

    # Analysis set: complete cases across TJ, CD8 and purity.
    a = scores.dropna(subset=["TJ_score", "CD8_score", "purity"]).copy()
    tj = a["TJ_score"].to_numpy()
    cd8 = a["CD8_score"].to_numpy()
    pur = a["purity"].to_numpy()

    # 1) Unadjusted Spearman TJ vs CD8
    rho_raw, p_raw = stats.spearmanr(tj, cd8)
    ci_raw = boot_ci(lambda x, y: stats.spearmanr(x, y)[0], tj, cd8)

    # 2) Partial Spearman TJ vs CD8 controlling for purity
    rho_par, p_par, n_par, ex, ey = partial_spearman(tj, cd8, pur)
    ci_par = boot_ci(lambda x, y, z: partial_spearman(x, y, z)[0], tj, cd8, pur)

    # Context: how strongly each score tracks purity
    rho_tj_pur, p_tj_pur = stats.spearmanr(tj, pur)
    rho_cd8_pur, p_cd8_pur = stats.spearmanr(cd8, pur)

    # Sensitivity: partial correlations controlling for Immune / Stromal score
    imm = a["Immune_score"].to_numpy()
    strm = a["Stromal_score"].to_numpy()
    rho_par_imm, p_par_imm, _, _, _ = partial_spearman(tj, cd8, imm)
    rho_par_strm, p_par_strm, _, _, _ = partial_spearman(tj, cd8, strm)

    results = {
        "dataset": "TCGA-LUAD (UCSC Xena IlluminaHiSeq log2 norm_count+1)",
        "n_samples_expression": int(expr.shape[1]),
        "n_samples_estimate": int(est.shape[0]),
        "n_matched_primary_tumor": int(len(common)),
        "n_analysis_complete_cases": int(len(a)),
        "tj_genes_requested": TJ_GENES,
        "tj_genes_found": tj_present,
        "tj_genes_missing": tj_missing,
        "cd8_genes_requested": CD8_GENES,
        "cd8_genes_found": cd8_present,
        "cd8_genes_missing": cd8_missing,
        "unadjusted_spearman_TJ_vs_CD8": {
            "rho": float(rho_raw), "p_value": float(p_raw),
            "ci95_low": ci_raw[0], "ci95_high": ci_raw[1],
        },
        "partial_spearman_TJ_vs_CD8_given_purity": {
            "rho": float(rho_par), "p_value": float(p_par), "n": int(n_par),
            "ci95_low": ci_par[0], "ci95_high": ci_par[1],
        },
        "context_TJ_vs_purity_spearman": {"rho": float(rho_tj_pur), "p_value": float(p_tj_pur)},
        "context_CD8_vs_purity_spearman": {"rho": float(rho_cd8_pur), "p_value": float(p_cd8_pur)},
        "sensitivity_partial_given_Immune_score": {"rho": float(rho_par_imm), "p_value": float(p_par_imm)},
        "sensitivity_partial_given_Stromal_score": {"rho": float(rho_par_strm), "p_value": float(p_par_strm)},
        "purity_formula": "cos(0.6049872018 + 0.0001467884 * ESTIMATE_score)",
        "score_method": "mean of per-gene z-scores across samples",
        "correlation_method": "Spearman; partial via rank-residual regression; 95% CI = 2000x percentile bootstrap",
        "rng_seed": RNG_SEED,
    }

    # ---- write tables ----
    scores.to_csv(OUT / "sample_scores.csv")
    tj_tbl.to_csv(OUT / "tj_genes_used.csv", index=False)
    cd8_tbl.to_csv(OUT / "cd8_genes_used.csv", index=False)
    with open(OUT / "correlation_results.json", "w") as f:
        json.dump(results, f, indent=2)

    flat = {
        "n_analysis": len(a),
        "unadj_rho": rho_raw, "unadj_p": p_raw,
        "unadj_ci_low": ci_raw[0], "unadj_ci_high": ci_raw[1],
        "partial_purity_rho": rho_par, "partial_purity_p": p_par,
        "partial_purity_ci_low": ci_par[0], "partial_purity_ci_high": ci_par[1],
        "TJ_vs_purity_rho": rho_tj_pur, "CD8_vs_purity_rho": rho_cd8_pur,
        "partial_Immune_rho": rho_par_imm, "partial_Stromal_rho": rho_par_strm,
    }
    pd.DataFrame([flat]).to_csv(OUT / "correlation_results.csv", index=False)

    # ---- plots ----
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(tj, cd8, c=pur, cmap="viridis", s=18, alpha=0.85, edgecolor="none")
    ax.set_xlabel("Tight-junction score (mean z)")
    ax.set_ylabel("CD8 T-cell score (mean z)")
    ax.set_title(f"TCGA-LUAD TJ vs CD8 (n={len(a)})\n"
                 f"Spearman rho={rho_raw:.3f} (p={p_raw:.1e})")
    fig.colorbar(sc, ax=ax, label="ESTIMATE purity")
    fig.tight_layout()
    fig.savefig(OUT / "scatter_tj_cd8.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(ex, ey, s=18, alpha=0.85, color="#c0392b", edgecolor="none")
    ax.set_xlabel("TJ rank residual (purity removed)")
    ax.set_ylabel("CD8 rank residual (purity removed)")
    ax.set_title(f"Purity-adjusted TJ vs CD8 (n={len(a)})\n"
                 f"partial rho={rho_par:.3f} (p={p_par:.1e})")
    fig.tight_layout()
    fig.savefig(OUT / "scatter_partial_residuals.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].scatter(pur, tj, s=16, alpha=0.8, color="#2c7fb8", edgecolor="none")
    axes[0].set_xlabel("ESTIMATE purity"); axes[0].set_ylabel("TJ score")
    axes[0].set_title(f"TJ vs purity (rho={rho_tj_pur:.3f})")
    axes[1].scatter(pur, cd8, s=16, alpha=0.8, color="#d95f0e", edgecolor="none")
    axes[1].set_xlabel("ESTIMATE purity"); axes[1].set_ylabel("CD8 score")
    axes[1].set_title(f"CD8 vs purity (rho={rho_cd8_pur:.3f})")
    fig.tight_layout()
    fig.savefig(OUT / "purity_associations.png", dpi=150)
    plt.close(fig)

    # ---- console summary ----
    print("=== B3 TCGA-LUAD TJ vs CD8 (ESTIMATE purity) ===")
    print(f"matched primary tumors: {len(common)}; complete cases: {len(a)}")
    print(f"TJ genes found:  {len(tj_present)}/{len(TJ_GENES)} missing={tj_missing}")
    print(f"CD8 genes found: {len(cd8_present)}/{len(CD8_GENES)} missing={cd8_missing}")
    print(f"Unadjusted   Spearman rho={rho_raw:+.3f}  p={p_raw:.2e}  CI95=({ci_raw[0]:+.3f},{ci_raw[1]:+.3f})")
    print(f"Partial|pur  Spearman rho={rho_par:+.3f}  p={p_par:.2e}  CI95=({ci_par[0]:+.3f},{ci_par[1]:+.3f})")
    print(f"TJ~purity rho={rho_tj_pur:+.3f}; CD8~purity rho={rho_cd8_pur:+.3f}")
    print(f"Sensitivity partial|Immune rho={rho_par_imm:+.3f}; partial|Stromal rho={rho_par_strm:+.3f}")
    return results


if __name__ == "__main__":
    main()

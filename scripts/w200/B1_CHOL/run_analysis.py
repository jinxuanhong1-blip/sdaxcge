#!/usr/bin/env python3
"""B1 analog in TCGA-CHOL: TACSTD2 (TROP2) × CLDN4 surface rank + co-expression.

Mirrors two sibling B1 slices on the same public-data stack:

  * B1_BRCA  — honest rank of CLDN4 among the Bausch-Fluck 2018 in-silico
    surfaceome as a TACSTD2 co-expression partner (window w200).
  * B1_PAAD  — TACSTD2–CLDN4 Spearman/Pearson with bootstrap CI, purity
    partial correlation, and transcriptome-wide partner rank.

Cohort is small (n=36 primary tumours). Every number written under
results/w200/B1_CHOL/ is computed here from the files fetched by
download_data.py. Nothing is hard-coded.

Analyses
  1. Primary tumours only (barcode sample-type 01), one sample per patient.
  2. Surface-abundance rank of TACSTD2 and CLDN4 among surfaceome genes
     (median log2(norm_count+1) in primary tumours).
  3. Pearson + Spearman of TACSTD2 vs CLDN4, with 2,000-fold case-resampling
     bootstrap 95% CI on Spearman rho (seed 20260816).
  4. Honest surfaceome co-expression ranking of every surface gene vs
     TACSTD2; report where CLDN4 actually lands. Write the full ranking
     plus the top-200 window.
  5. Transcriptome-wide Spearman rank of CLDN4 among expressed genes.
  6. Sensitivity: rank-based partial Spearman given ABSOLUTE purity
     (all 36 tumours have a called ABSOLUTE solution).
  7. Sensitivity: intrahepatic-only (n=32); the other four cases are
     extrahepatic / liver / gallbladder.
  8. Tumour vs adjacent-normal (9 matched pairs): Wilcoxon signed-rank
     plus unpaired Mann–Whitney on 36 vs 9.

Inputs:  $B1_CHOL_DATA  (default /tmp/b1_chol_data)
Outputs: results/w200/B1_CHOL/
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.environ.get("B1_CHOL_DATA", "/tmp/b1_chol_data")
OUT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "results", "w200", "B1_CHOL")
)
os.makedirs(OUT_DIR, exist_ok=True)

RNG = np.random.default_rng(20260816)
ANCHOR = "TACSTD2"
FOCUS = "CLDN4"
WINDOW = 200
BOOT = 2000


def load_expression(path: str) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    return expr


def load_surfaceome(path: str) -> list[str]:
    surf = pd.read_excel(
        path, sheet_name="in silico surfaceome only", engine="openpyxl", header=1
    )
    genes = (
        surf["UniProt gene"]
        .dropna()
        .astype(str)
        .replace({"nan": np.nan, "None": np.nan})
        .dropna()
        .tolist()
    )
    return sorted(set(genes))


def sample_type(barcode: str) -> str:
    parts = str(barcode).split("-")
    return parts[3][:2] if len(parts) >= 4 else ""


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = BOOT) -> tuple[float, float]:
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def partial_spearman(x, y, z) -> tuple[float, float]:
    """Spearman partial correlation of x,y given z (residualised ranks)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))

    def resid(a, b):
        beta = np.polyfit(b, a, 1)
        return a - np.polyval(beta, b)

    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def corr_block(xv, yv, label: str) -> dict:
    pr, pp = stats.pearsonr(xv, yv)
    sr, sp = stats.spearmanr(xv, yv)
    lo, hi = spearman_ci(np.asarray(xv, float), np.asarray(yv, float))
    return {
        "subset": label,
        "n": int(len(xv)),
        "pearson_r": round(float(pr), 4),
        "pearson_p": float(pp),
        "spearman_rho": round(float(sr), 4),
        "spearman_p": float(sp),
        "spearman_rho_ci95_lo": round(lo, 4),
        "spearman_rho_ci95_hi": round(hi, 4),
    }


def _corr_matrix_vs_vector(mat: np.ndarray, vec: np.ndarray) -> np.ndarray:
    vc = vec - vec.mean()
    mc = mat - mat.mean(axis=1, keepdims=True)
    num = mc @ vc
    den = np.sqrt((mc**2).sum(axis=1) * (vc**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def _pvalues(r: np.ndarray, n: int) -> np.ndarray:
    r = np.clip(r, -0.999999999, 0.999999999)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / (1.0 - r**2))
    return 2 * stats.t.sf(np.abs(t), df=n - 2)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    q = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def rank_partners(expr: pd.DataFrame, universe: list[str], samples: list[str]) -> pd.DataFrame:
    sub = expr.loc[universe, samples]
    anchor_vec = expr.loc[ANCHOR, samples].to_numpy(dtype=float)
    mat = sub.to_numpy(dtype=float)
    n = len(samples)
    pear_r = _corr_matrix_vs_vector(mat, anchor_vec)
    anchor_rank = stats.rankdata(anchor_vec)
    mat_rank = np.apply_along_axis(stats.rankdata, 1, mat)
    spear_r = _corr_matrix_vs_vector(mat_rank, anchor_rank)
    df = pd.DataFrame(
        {
            "gene": universe,
            "spearman_r": spear_r,
            "spearman_p": _pvalues(spear_r, n),
            "pearson_r": pear_r,
            "pearson_p": _pvalues(pear_r, n),
            "mean_log2_expr": mat.mean(axis=1),
            "median_log2_expr": np.median(mat, axis=1),
            "pct_expressed": (mat > 0).mean(axis=1) * 100.0,
        }
    )
    df["spearman_q"] = _bh_fdr(df["spearman_p"].to_numpy())
    df["pearson_q"] = _bh_fdr(df["pearson_p"].to_numpy())
    # Zero-variance surface genes yield NaN correlations and cannot be ranked.
    n_dropped = int(df["spearman_r"].isna().sum())
    df = df.dropna(subset=["spearman_r"]).copy()
    df.attrs["n_dropped_zero_variance"] = n_dropped
    df = df.sort_values("spearman_r", ascending=False, kind="mergesort").reset_index(drop=True)
    df.insert(1, "spearman_rank", np.arange(1, len(df) + 1))
    df["pearson_rank"] = (
        df["pearson_r"].rank(ascending=False, method="min").astype("Int64")
    )
    return df


def focus_report(df: pd.DataFrame, n_universe: int) -> dict:
    row = df[df["gene"] == FOCUS]
    if row.empty:
        return {"gene": FOCUS, "in_universe": False}
    r = row.iloc[0]
    s_rank = int(r["spearman_rank"])
    p_rank = int(r["pearson_rank"])
    return {
        "gene": FOCUS,
        "in_universe": True,
        "is_top_spearman": s_rank == 1,
        "is_top_pearson": p_rank == 1,
        "spearman_rank": s_rank,
        "spearman_r": float(r["spearman_r"]),
        "spearman_q": float(r["spearman_q"]),
        "spearman_percentile": round(100 * (1 - (s_rank - 1) / n_universe), 3),
        "pearson_rank": p_rank,
        "pearson_r": float(r["pearson_r"]),
        "pearson_q": float(r["pearson_q"]),
        "pearson_percentile": round(100 * (1 - (p_rank - 1) / n_universe), 3),
    }


def main() -> None:
    expr = load_expression(os.path.join(DATA_DIR, "CHOL.HiSeqV2.gz"))
    surf_genes = load_surfaceome(os.path.join(DATA_DIR, "table_S3_surfaceome.xlsx"))

    tumor_cols = sorted(c for c in expr.columns if sample_type(c) == "01")
    normal_cols = sorted(c for c in expr.columns if sample_type(c) == "11")
    # one sample per patient (lexicographically first vial)
    patients: dict[str, str] = {}
    for c in tumor_cols:
        patients.setdefault(c[:12], c)
    tumor_cols = sorted(patients.values())
    n_tumor = len(tumor_cols)

    if ANCHOR not in expr.index or FOCUS not in expr.index:
        raise SystemExit(f"{ANCHOR}/{FOCUS} missing from expression matrix")

    x = expr.loc[ANCHOR, tumor_cols].to_numpy(dtype=float)
    y = expr.loc[FOCUS, tumor_cols].to_numpy(dtype=float)

    # ---- clinical site-of-origin ------------------------------------------
    cl = pd.read_csv(os.path.join(DATA_DIR, "TCGA-CHOL.clinical.tsv.gz"), sep="\t", low_memory=False)
    site = (
        cl.drop_duplicates("submitter_id")
        .set_index("submitter_id")["tissue_or_organ_of_origin.diagnoses"]
    )
    site_series = pd.Series({c: site.get(c[:12], np.nan) for c in tumor_cols})
    ih_mask = site_series.eq("Intrahepatic bile duct").to_numpy()

    # ---- ABSOLUTE purity ---------------------------------------------------
    pur = pd.read_csv(os.path.join(DATA_DIR, "tcga_absolute_purity.txt"), sep="\t")
    pur["sample15"] = pur["array"].astype(str).str[:15]
    purity = pur.drop_duplicates("sample15").set_index("sample15")["purity"]
    pvec = pd.Series(purity).reindex([c[:15] for c in tumor_cols]).to_numpy(dtype=float)
    p_mask = ~np.isnan(pvec)

    # ---- 3. TACSTD2 vs CLDN4 co-expression --------------------------------
    results = [
        corr_block(x, y, "all_primary_tumors"),
        corr_block(x[ih_mask], y[ih_mask], "intrahepatic_only"),
    ]
    pd.DataFrame(results).to_csv(os.path.join(OUT_DIR, "coexpression_main.csv"), index=False)

    r_part, p_part = partial_spearman(x[p_mask], y[p_mask], pvec[p_mask])
    sr_sub, sp_sub = stats.spearmanr(x[p_mask], y[p_mask])
    purity_rows = [
        {
            "quantity": "spearman_unadjusted_purity_subset",
            "n": int(p_mask.sum()),
            "rho": round(float(sr_sub), 4),
            "p": float(sp_sub),
        },
        {
            "quantity": "spearman_partial_given_ABSOLUTE_purity",
            "n": int(p_mask.sum()),
            "rho": round(r_part, 4),
            "p": p_part,
        },
        {
            "quantity": "spearman_TACSTD2_vs_purity",
            "n": int(p_mask.sum()),
            "rho": round(float(stats.spearmanr(x[p_mask], pvec[p_mask]).statistic), 4),
            "p": float(stats.spearmanr(x[p_mask], pvec[p_mask]).pvalue),
        },
        {
            "quantity": "spearman_CLDN4_vs_purity",
            "n": int(p_mask.sum()),
            "rho": round(float(stats.spearmanr(y[p_mask], pvec[p_mask]).statistic), 4),
            "p": float(stats.spearmanr(y[p_mask], pvec[p_mask]).pvalue),
        },
    ]
    pd.DataFrame(purity_rows).to_csv(os.path.join(OUT_DIR, "purity_adjusted.csv"), index=False)

    # ---- 2. surface abundance rank ----------------------------------------
    surf_in_matrix = [g for g in surf_genes if g in expr.index]
    med = expr.loc[surf_in_matrix, tumor_cols].median(axis=1).sort_values(ascending=False)
    abund = med.reset_index()
    abund.columns = ["gene", "median_log2_normcount1"]
    abund.insert(0, "abundance_rank", range(1, len(abund) + 1))
    abund["mean_log2_normcount1"] = expr.loc[abund["gene"], tumor_cols].mean(axis=1).to_numpy()
    abund.to_csv(os.path.join(OUT_DIR, "surface_abundance_rank.csv"), index=False)

    def abund_row(gene: str) -> dict:
        rnk = int(abund.index[abund["gene"] == gene][0]) + 1
        return {
            "gene": gene,
            "abundance_rank": rnk,
            "n_surface_genes_ranked": int(len(abund)),
            "percentile": round(100 * (1 - (rnk - 1) / len(abund)), 3),
            "median_log2_normcount1": float(med[gene]),
            "mean_log2_normcount1": float(expr.loc[gene, tumor_cols].mean()),
        }

    abund_focus = {g: abund_row(g) for g in (ANCHOR, FOCUS)}

    # ---- 4. surfaceome co-expression ranking ------------------------------
    universe = [g for g in surf_in_matrix if g != ANCHOR]
    rank_df = rank_partners(expr, universe, tumor_cols)
    n_dropped_zv = int(rank_df.attrs.get("n_dropped_zero_variance", 0))
    rank_df.to_csv(os.path.join(OUT_DIR, f"coexpression_{ANCHOR}_surfaceome.csv"), index=False)
    rank_df.head(WINDOW).to_csv(os.path.join(OUT_DIR, f"top{WINDOW}.csv"), index=False)
    report = focus_report(rank_df, int(len(rank_df)))

    # ---- 5. transcriptome-wide rank ---------------------------------------
    mat = expr.loc[:, tumor_cols]
    mat = mat[~mat.index.isna()]
    expressed = mat[(mat > 1.0).mean(axis=1) >= 0.20]
    expressed = expressed[expressed.index != ANCHOR]
    expressed = (
        expressed.assign(_mean=expressed.mean(axis=1))
        .sort_values("_mean", ascending=False)
        .drop(columns="_mean")
    )
    expressed = expressed[~expressed.index.duplicated()]
    xr = stats.rankdata(x)
    xr_c = xr - xr.mean()
    m_ranks = np.apply_along_axis(stats.rankdata, 1, expressed.to_numpy(dtype=float))
    m_c = m_ranks - m_ranks.mean(axis=1, keepdims=True)
    num = m_c @ xr_c
    den = np.sqrt((m_c**2).sum(axis=1) * (xr_c**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        rhos = pd.Series(num / den, index=expressed.index).sort_values(ascending=False)
    cldn4_tw_rank = int(rhos.index.get_loc(FOCUS)) + 1
    top25 = rhos.head(25).round(4).reset_index()
    top25.columns = ["gene", "spearman_rho_vs_TACSTD2"]
    top25.insert(0, "rank", range(1, len(top25) + 1))
    top25.to_csv(os.path.join(OUT_DIR, "top25_tacstd2_partners.csv"), index=False)

    # ---- 8. tumour vs normal ----------------------------------------------
    paired_pats = sorted(set(c[:12] for c in tumor_cols) & set(c[:12] for c in normal_cols))
    tn_rows = []
    for gene in (ANCHOR, FOCUS):
        t_all = expr.loc[gene, tumor_cols].to_numpy(dtype=float)
        n_all = expr.loc[gene, normal_cols].to_numpy(dtype=float)
        u, up = stats.mannwhitneyu(t_all, n_all, alternative="two-sided")
        t_pair = expr.loc[gene, [patients[p] for p in paired_pats]].to_numpy(dtype=float)
        n_pair = expr.loc[gene, [c for c in normal_cols if c[:12] in set(paired_pats)]].to_numpy(
            dtype=float
        )
        # align pairs by patient
        n_map = {c[:12]: c for c in normal_cols}
        t_pair = np.array([expr.loc[gene, patients[p]] for p in paired_pats], dtype=float)
        n_pair = np.array([expr.loc[gene, n_map[p]] for p in paired_pats], dtype=float)
        w, wp = stats.wilcoxon(t_pair, n_pair, alternative="two-sided", zero_method="wilcox")
        tn_rows.append(
            {
                "gene": gene,
                "n_tumor": int(len(t_all)),
                "n_normal": int(len(n_all)),
                "median_tumor": round(float(np.median(t_all)), 4),
                "median_normal": round(float(np.median(n_all)), 4),
                "median_delta_unpaired": round(float(np.median(t_all) - np.median(n_all)), 4),
                "mannwhitney_U": float(u),
                "mannwhitney_p": float(up),
                "n_paired": int(len(paired_pats)),
                "median_paired_delta": round(float(np.median(t_pair - n_pair)), 4),
                "wilcoxon_W": float(w),
                "wilcoxon_p": float(wp),
            }
        )
    pd.DataFrame(tn_rows).to_csv(os.path.join(OUT_DIR, "tumor_vs_normal.csv"), index=False)

    # ---- per-sample table --------------------------------------------------
    sample_df = pd.DataFrame(
        {
            "barcode": tumor_cols,
            "patient": [c[:12] for c in tumor_cols],
            "TACSTD2": x,
            "CLDN4": y,
            "ABSOLUTE_purity": pvec,
            "site_of_origin": site_series.to_numpy(),
            "intrahepatic": ih_mask,
        }
    )
    sample_df.to_csv(os.path.join(OUT_DIR, "sample_level.csv"), index=False)

    # ---- figures -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.4, 5.1), dpi=150)
    ax.scatter(
        x[ih_mask],
        y[ih_mask],
        s=28,
        alpha=0.75,
        color="#2b6a99",
        label=f"intrahepatic (n={int(ih_mask.sum())})",
    )
    ax.scatter(
        x[~ih_mask],
        y[~ih_mask],
        s=36,
        alpha=0.9,
        color="#d1495b",
        marker="^",
        label=f"other site (n={int((~ih_mask).sum())})",
    )
    r0 = results[0]
    ax.set_xlabel("TACSTD2  log2(norm_count+1)")
    ax.set_ylabel("CLDN4  log2(norm_count+1)")
    ax.set_title(
        f"TCGA-CHOL primary tumors (n={n_tumor})\n"
        f"Spearman rho = {r0['spearman_rho']:.2f} "
        f"[{r0['spearman_rho_ci95_lo']:.2f}, {r0['spearman_rho_ci95_hi']:.2f}], "
        f"p = {r0['spearman_p']:.1e}",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig1_scatter_tacstd2_cldn4.png"))
    fig.savefig(os.path.join(OUT_DIR, "scatter_TACSTD2_vs_CLDN4.png"))
    plt.close(fig)

    topn = rank_df.head(20).iloc[::-1]
    colors = ["#d62728" if g == FOCUS else "#4c78a8" for g in topn["gene"]]
    fig, ax = plt.subplots(figsize=(6.2, 7.0), dpi=150)
    ax.barh(topn["gene"], topn["spearman_r"], color=colors)
    ax.set_xlabel("Spearman rho with TACSTD2")
    ax.set_title(
        f"Top surface-gene co-expression with TACSTD2\n"
        f"TCGA-CHOL (n={n_tumor}; focus CLDN4 highlighted)"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig2_top_surface_partners.png"))
    fig.savefig(os.path.join(OUT_DIR, "top_partners_TACSTD2.png"))
    plt.close(fig)

    # abundance context: TACSTD2 / CLDN4 vs surfaceome median distribution
    fig, ax = plt.subplots(figsize=(6.4, 3.6), dpi=150)
    ax.hist(med.values, bins=40, color="#9bb7d4", edgecolor="white")
    ax.axvline(med[ANCHOR], color="#2b6a99", lw=2, label=f"{ANCHOR} rank #{abund_focus[ANCHOR]['abundance_rank']}")
    ax.axvline(med[FOCUS], color="#d1495b", lw=2, label=f"{FOCUS} rank #{abund_focus[FOCUS]['abundance_rank']}")
    ax.set_xlabel("median log2(norm_count+1) in TCGA-CHOL primary tumors")
    ax.set_ylabel("surfaceome genes")
    ax.set_title("Surface-abundance rank (Bausch-Fluck 2018 ∩ HiSeqV2)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig3_surface_abundance.png"))
    plt.close(fig)

    # ---- summary json ------------------------------------------------------
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cohort": "TCGA-CHOL",
        "expression_dataset": "UCSC Xena TCGA.CHOL.sampleMap/HiSeqV2 (log2 norm_count+1)",
        "surfaceome_source": "Bausch-Fluck et al. 2018 in-silico surfaceome (table S3)",
        "purity_source": "PanCanAtlas ABSOLUTE (GDC 4f277128-f793-4354-a13d-30cc7fe9f6b5)",
        "anchor_gene": ANCHOR,
        "focus_gene": FOCUS,
        "sample_type_codes": ["01"],
        "n_samples": n_tumor,
        "n_normals": len(normal_cols),
        "n_paired_normals": len(paired_pats),
        "n_intrahepatic": int(ih_mask.sum()),
        "n_with_ABSOLUTE_purity": int(p_mask.sum()),
        "n_surfaceome_table": len(surf_genes),
        "n_surfaceome_in_matrix": len(surf_in_matrix),
        "surface_gene_universe": int(len(rank_df)),
        "n_surface_zero_variance_dropped": n_dropped_zv,
        "correlation_primary": "spearman",
        "window": WINDOW,
        "n_expressed_genes_ranked": int(len(rhos)),
        "CLDN4_rank_among_TACSTD2_partners_transcriptome": cldn4_tw_rank,
        "CLDN4_percentile_transcriptome": round(100 * (1 - cldn4_tw_rank / len(rhos)), 2),
        "spearman_rho_main": results[0]["spearman_rho"],
        "spearman_p_main": results[0]["spearman_p"],
        "spearman_rho_ci95": [
            results[0]["spearman_rho_ci95_lo"],
            results[0]["spearman_rho_ci95_hi"],
        ],
        "pearson_r_main": results[0]["pearson_r"],
        "median_log2_TACSTD2": round(float(np.median(x)), 3),
        "median_log2_CLDN4": round(float(np.median(y)), 3),
        "surface_abundance": abund_focus,
        "focus_result": report,
        "top10_spearman": rank_df.head(10)[
            ["gene", "spearman_rank", "spearman_r", "spearman_q"]
        ].to_dict("records"),
        "purity": purity_rows,
        "tumor_vs_normal": tn_rows,
        "coexpression_subsets": results,
    }
    with open(os.path.join(OUT_DIR, "cohort_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(results).to_string(index=False))
    print(pd.DataFrame(purity_rows).to_string(index=False))
    print(pd.DataFrame(tn_rows).to_string(index=False))
    print("\nTop 10 surface partners:")
    print(rank_df.head(10)[["spearman_rank", "gene", "spearman_r", "spearman_q"]].to_string(index=False))
    print(f"\n[done] outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()

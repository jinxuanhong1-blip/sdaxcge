#!/usr/bin/env python3
"""B1 analog in TCGA-PAAD: TACSTD2 (TROP2) x CLDN4 co-expression.

Mirrors the B1 co-expression readout from the LUAD/LUSC slice
(results/fable_tcga: rho = 0.53 LUAD, 0.39 LUSC) in TCGA-PAAD.

Analyses
  1. Primary tumors only (sample code -01), one sample per patient.
  2. Pearson + Spearman correlation of log2(TPM+1) TACSTD2 vs CLDN4,
     with bootstrap 95% CI on Spearman rho.
  3. Transcriptome-wide Spearman of every gene vs TACSTD2; report the
     rank/percentile of CLDN4 among expressed genes.
  4. Sensitivity: exclude neuroendocrine carcinomas (morphology 8246/3),
     which are known non-PDAC contaminants of TCGA-PAAD.
  5. Sensitivity: rank-based partial correlation adjusting for ABSOLUTE
     tumor purity (PAAD is stroma-rich; shared epithelial content can
     inflate co-expression of two epithelial genes).

Inputs come from scripts/w200/B1_PAAD/download_data.py (default
/tmp/b1_paad_data). Outputs land in results/w200/B1_PAAD/.
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.environ.get("B1_PAAD_DATA", "/tmp/b1_paad_data")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                       "results", "w200", "B1_PAAD")
OUT_DIR = os.path.abspath(OUT_DIR)
os.makedirs(OUT_DIR, exist_ok=True)

RNG = np.random.default_rng(20260816)
GENES = ["TACSTD2", "CLDN4"]


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000) -> tuple:
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y given z (residualised ranks)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    def resid(a, b):
        beta = np.polyfit(b, a, 1)
        return a - np.polyval(beta, b)
    r, p = stats.pearsonr(resid(rx, rz), resid(ry, rz))
    return float(r), float(p)


def main() -> None:
    # ---- expression matrix -------------------------------------------------
    expr = pd.read_csv(os.path.join(DATA_DIR, "TCGA-PAAD.star_tpm.tsv.gz"),
                       sep="\t", index_col=0)
    probemap = pd.read_csv(os.path.join(DATA_DIR, "gencode.v36.probemap"),
                           sep="\t")
    id2sym = probemap.set_index("id")["gene"]
    expr.index = id2sym.reindex(expr.index).values

    # primary tumors only (-01), one per patient (lexicographically first vial)
    tumor_cols = sorted(c for c in expr.columns if c[13:15] == "01")
    patients = {}
    for c in tumor_cols:
        patients.setdefault(c[:12], c)
    tumor_cols = sorted(patients.values())
    tum = expr[tumor_cols]

    # ---- clinical: neuroendocrine flag ------------------------------------
    cl = pd.read_csv(os.path.join(DATA_DIR, "TCGA-PAAD.clinical.tsv.gz"),
                     sep="\t", low_memory=False)
    morph = cl.set_index(cl["sample"].str[:12])["morphology.diagnoses"]
    morph = morph[~morph.index.duplicated()]
    ne_patients = set(morph[morph == "8246/3"].index)

    # ---- ABSOLUTE purity ---------------------------------------------------
    pur = pd.read_csv(os.path.join(DATA_DIR, "tcga_absolute_purity.txt"),
                      sep="\t")
    pur["sample15"] = pur["array"].str[:15]
    purity = pur.set_index("sample15")["purity"]
    purity = purity[~purity.index.duplicated()]

    # ---- gene vectors ------------------------------------------------------
    g = {}
    for gene in GENES:
        rows = tum.loc[tum.index == gene]
        assert len(rows) == 1, f"{gene}: {len(rows)} rows found"
        g[gene] = rows.iloc[0]

    x = g["TACSTD2"].values
    y = g["CLDN4"].values
    n_all = len(x)

    def corr_block(xv, yv, label, n_boot=2000):
        pr, pp = stats.pearsonr(xv, yv)
        sr, sp = stats.spearmanr(xv, yv)
        lo, hi = spearman_ci(np.asarray(xv), np.asarray(yv), n_boot)
        return {
            "subset": label, "n": int(len(xv)),
            "pearson_r": round(float(pr), 4), "pearson_p": float(pp),
            "spearman_rho": round(float(sr), 4), "spearman_p": float(sp),
            "spearman_rho_ci95_lo": round(lo, 4),
            "spearman_rho_ci95_hi": round(hi, 4),
        }

    results = [corr_block(x, y, "all_primary_tumors")]

    # sensitivity: exclude neuroendocrine carcinomas
    keep = [i for i, c in enumerate(tumor_cols) if c[:12] not in ne_patients]
    results.append(corr_block(x[keep], y[keep],
                              "excl_neuroendocrine_8246_3"))

    # sensitivity: purity-adjusted partial Spearman
    pvec = pd.Series(purity).reindex([c[:15] for c in tumor_cols]).values
    mask = ~np.isnan(pvec)
    r_part, p_part = partial_spearman(x[mask], y[mask], pvec[mask])
    sr_sub, sp_sub = stats.spearmanr(x[mask], y[mask])
    purity_rows = [
        {"quantity": "spearman_unadjusted_purity_subset",
         "n": int(mask.sum()), "rho": round(float(sr_sub), 4),
         "p": float(sp_sub)},
        {"quantity": "spearman_partial_given_ABSOLUTE_purity",
         "n": int(mask.sum()), "rho": round(r_part, 4), "p": p_part},
        {"quantity": "spearman_TACSTD2_vs_purity", "n": int(mask.sum()),
         "rho": round(float(stats.spearmanr(x[mask], pvec[mask]).statistic), 4),
         "p": float(stats.spearmanr(x[mask], pvec[mask]).pvalue)},
        {"quantity": "spearman_CLDN4_vs_purity", "n": int(mask.sum()),
         "rho": round(float(stats.spearmanr(y[mask], pvec[mask]).statistic), 4),
         "p": float(stats.spearmanr(y[mask], pvec[mask]).pvalue)},
    ]

    pd.DataFrame(results).to_csv(
        os.path.join(OUT_DIR, "coexpression_main.csv"), index=False)
    pd.DataFrame(purity_rows).to_csv(
        os.path.join(OUT_DIR, "purity_adjusted.csv"), index=False)

    # ---- transcriptome-wide rank of CLDN4 among TACSTD2 partners ----------
    mat = tum[tumor_cols]
    mat = mat[~mat.index.isna()]
    expressed = mat[(mat > 1.0).mean(axis=1) >= 0.20]          # log2(TPM+1) > 1
    expressed = expressed[expressed.index != "TACSTD2"]
    # collapse duplicate symbols by keeping the highest-mean row
    expressed = (expressed.assign(_mean=expressed.mean(axis=1))
                 .sort_values("_mean", ascending=False)
                 .drop(columns="_mean"))
    expressed = expressed[~expressed.index.duplicated()]

    xr = stats.rankdata(x)
    xr_c = xr - xr.mean()
    m_ranks = np.apply_along_axis(stats.rankdata, 1, expressed.values)
    m_c = m_ranks - m_ranks.mean(axis=1, keepdims=True)
    num = m_c @ xr_c
    den = np.sqrt((m_c ** 2).sum(axis=1) * (xr_c ** 2).sum())
    rhos = pd.Series(num / den, index=expressed.index).sort_values(
        ascending=False)

    cldn4_rank = int(rhos.index.get_loc("CLDN4")) + 1
    top = rhos.head(25).round(4)
    top_df = top.reset_index()
    top_df.columns = ["gene", "spearman_rho_vs_TACSTD2"]
    top_df.insert(0, "rank", range(1, len(top_df) + 1))
    top_df.to_csv(os.path.join(OUT_DIR, "top25_tacstd2_partners.csv"),
                  index=False)

    summary = {
        "cohort": "TCGA-PAAD",
        "expression": "Xena GDC hub STAR TPM, log2(TPM+1), GENCODE v36",
        "n_primary_tumors": n_all,
        "n_excl_neuroendocrine": len(keep),
        "n_with_ABSOLUTE_purity": int(mask.sum()),
        "n_expressed_genes_ranked": int(len(rhos)),
        "CLDN4_rank_among_TACSTD2_partners": cldn4_rank,
        "CLDN4_percentile": round(100 * (1 - cldn4_rank / len(rhos)), 2),
        "spearman_rho_main": results[0]["spearman_rho"],
        "median_log2tpm1_TACSTD2": round(float(np.median(x)), 3),
        "median_log2tpm1_CLDN4": round(float(np.median(y)), 3),
    }
    with open(os.path.join(OUT_DIR, "cohort_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- scatter figure ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.2, 5.0), dpi=150)
    ne_mask = np.array([c[:12] in ne_patients for c in tumor_cols])
    ax.scatter(x[~ne_mask], y[~ne_mask], s=18, alpha=0.65,
               color="#2b6a99", label=f"other histology (n={int((~ne_mask).sum())})")
    ax.scatter(x[ne_mask], y[ne_mask], s=26, alpha=0.9, color="#d1495b",
               marker="^", label=f"neuroendocrine (n={int(ne_mask.sum())})")
    r0 = results[0]
    ax.set_xlabel("TACSTD2  log2(TPM+1)")
    ax.set_ylabel("CLDN4  log2(TPM+1)")
    ax.set_title(f"TCGA-PAAD primary tumors (n={n_all})\n"
                 f"Spearman rho = {r0['spearman_rho']:.2f} "
                 f"[{r0['spearman_rho_ci95_lo']:.2f}, "
                 f"{r0['spearman_rho_ci95_hi']:.2f}], "
                 f"p = {r0['spearman_p']:.1e}", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig1_scatter_tacstd2_cldn4.png"))

    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(results).to_string(index=False))
    print(pd.DataFrame(purity_rows).to_string(index=False))
    print("\nTop 10 TACSTD2 partners:")
    print(top_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()

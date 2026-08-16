"""Bulk RNA-seq: TACSTD2/CLDN4 vs approximate SCLC subtype and immune context.

Cohorts:
  - George et al. 2015 (cBioPortal sclc_ucologne_2015, RNA-seq FPKM panel)
  - GSE60052 (Jiang et al. 2016, 79 SCLC tumors + 7 normal lung, log2 matrix)

Subtypes are assigned by NMF (k=4, top 1250 variable genes), replicating the
approach of Gay et al. 2021 (Cancer Cell); factors are labelled by their
ASCL1/NEUROD1/POU2F3 loadings and the TF-low factor is SCLC-I. A simpler
marker-threshold rule is written alongside as a sensitivity annotation.

Statistics (all two-sided, alpha = 0.05, BH FDR within each family):
  - Kruskal-Wallis across the four NMF subtypes (SCLC-A/N/P/I)
  - Mann-Whitney U: SCLC-I vs all other subtypes pooled
  - Spearman correlation of TACSTD2/CLDN4 with immune signatures
  - GSE60052 only: tumor vs normal Mann-Whitney U
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from common import (DATA_DIR, TABLES_DIR, FIGURES_DIR, TARGET_GENES,
                    TEFF_GENES, CYT_GENES, APM_GENES, ensure_dirs,
                    assign_sclc_subtype, nmf_subtypes, signature_score, bh_fdr)

SUBTYPE_ORDER = ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I"]
IMMUNE_SCORES = {
    "Teff_score": TEFF_GENES,
    "Cytolytic_score": CYT_GENES,
    "AntigenPresentation_score": APM_GENES,
}
IMMUNE_SINGLE_GENES = ["CD8A", "CD274", "CXCL9", "STAT1", "CD68", "CD163"]


def load_george():
    fpkm = pd.read_csv(DATA_DIR / "george2015_fpkm_full.csv.gz", index_col=0)
    log = np.log2(fpkm.clip(lower=0) + 1)
    return log.dropna(axis=1, how="any")


def load_gse60052():
    df = pd.read_csv(
        DATA_DIR / "GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv.gz",
        sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df = df.drop_duplicates(subset="gene").set_index("gene")
    tumors = [c for c in df.columns if not c.endswith(".normal")]
    normals = [c for c in df.columns if c.endswith(".normal")]
    assert len(tumors) == 79 and len(normals) == 7, (len(tumors), len(normals))
    return df[tumors].T, df[normals].T  # samples x genes, already log2


def analyze_cohort(name, log_expr, rows_subtype, rows_corr):
    subtype, weights = nmf_subtypes(log_expr)
    marker_subtype = assign_sclc_subtype(log_expr)
    agree = (subtype == marker_subtype).mean()
    print(f"[{name}] NMF vs marker-rule agreement: {agree:.0%}")
    scores = {k: signature_score(log_expr, v) for k, v in IMMUNE_SCORES.items()}

    # sanity check: SCLC-I should be the most immune-infiltrated class
    teff = scores["Teff_score"]
    u, p = stats.mannwhitneyu(teff[subtype == "SCLC-I"],
                              teff[subtype != "SCLC-I"], alternative="two-sided")
    print(f"[{name}] subtype counts: {subtype.value_counts().to_dict()}")
    print(f"[{name}] Teff score SCLC-I vs rest: U={u:.1f}, p={p:.2e} "
          f"(median {teff[subtype == 'SCLC-I'].median():.2f} vs "
          f"{teff[subtype != 'SCLC-I'].median():.2f})")

    for gene in TARGET_GENES:
        if gene not in log_expr.columns:
            continue
        v = log_expr[gene]
        groups = [v[subtype == s] for s in SUBTYPE_ORDER if (subtype == s).any()]
        kw_h, kw_p = stats.kruskal(*groups)
        in_i, out_i = v[subtype == "SCLC-I"], v[subtype != "SCLC-I"]
        mw_u, mw_p = stats.mannwhitneyu(in_i, out_i, alternative="two-sided")
        rows_subtype.append({
            "cohort": name, "gene": gene, "n": len(v),
            **{f"median_{s}": v[subtype == s].median() for s in SUBTYPE_ORDER},
            "kruskal_H": kw_h, "kruskal_p": kw_p,
            "median_SCLC-I": in_i.median(), "median_others": out_i.median(),
            "mwu_U_IvsRest": mw_u, "mwu_p_IvsRest": mw_p,
        })
        for score_name, score in scores.items():
            rho, sp = stats.spearmanr(v, score)
            rows_corr.append({"cohort": name, "gene": gene,
                              "covariate": score_name, "n": len(v),
                              "spearman_rho": rho, "spearman_p": sp})
        for ig in IMMUNE_SINGLE_GENES:
            if ig in log_expr.columns:
                rho, sp = stats.spearmanr(v, log_expr[ig])
                rows_corr.append({"cohort": name, "gene": gene,
                                  "covariate": ig, "n": len(v),
                                  "spearman_rho": rho, "spearman_p": sp})
    return subtype, marker_subtype, scores


def plot_cohort(name, log_expr, subtype, teff):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    df = log_expr.copy()
    df["subtype"] = subtype.values
    for ax, gene in zip(axes[:2], TARGET_GENES):
        sns.boxplot(data=df, x="subtype", y=gene, order=SUBTYPE_ORDER,
                    hue="subtype", legend=False,
                    palette="Set2", showfliers=False, ax=ax)
        sns.stripplot(data=df, x="subtype", y=gene, order=SUBTYPE_ORDER,
                      color="k", size=2.5, alpha=0.6, ax=ax)
        ax.set_ylabel(f"{gene} (log2 expr)")
        ax.set_xlabel("NMF subtype")
        ax.set_title(f"{name}: {gene}")
    ax = axes[2]
    ax.scatter(teff, df["TACSTD2"], s=12, alpha=0.7)
    rho, p = stats.spearmanr(teff, df["TACSTD2"])
    ax.set_xlabel("T-effector score (z)")
    ax.set_ylabel("TACSTD2 (log2 expr)")
    ax.set_title(f"{name}: rho={rho:.2f}, p={p:.1e}")
    fig.tight_layout()
    out = FIGURES_DIR / f"bulk_{name}_tacstd2_cldn4.png"
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print(f"[fig ] {out}")


def main():
    ensure_dirs()
    rows_subtype, rows_corr = [], []

    george = load_george()
    subtype_g, marker_g, scores_g = analyze_cohort("George2015", george,
                                                   rows_subtype, rows_corr)
    plot_cohort("George2015", george, subtype_g, scores_g["Teff_score"])

    tum, norm = load_gse60052()
    subtype_j, marker_j, scores_j = analyze_cohort("GSE60052", tum,
                                                   rows_subtype, rows_corr)
    plot_cohort("GSE60052", tum, subtype_j, scores_j["Teff_score"])

    # tumor vs normal in GSE60052
    tvn = []
    for gene in TARGET_GENES:
        u, p = stats.mannwhitneyu(tum[gene], norm[gene], alternative="two-sided")
        tvn.append({"gene": gene, "median_tumor": tum[gene].median(),
                    "median_normal": norm[gene].median(),
                    "n_tumor": len(tum), "n_normal": len(norm),
                    "mwu_U": u, "mwu_p": p})
    tvn = pd.DataFrame(tvn)
    tvn["mwu_p_BH"] = bh_fdr(tvn["mwu_p"])
    tvn.to_csv(TABLES_DIR / "bulk_gse60052_tumor_vs_normal.csv", index=False)

    sub = pd.DataFrame(rows_subtype)
    sub["kruskal_p_BH"] = bh_fdr(sub["kruskal_p"])
    sub["mwu_p_IvsRest_BH"] = bh_fdr(sub["mwu_p_IvsRest"])
    sub.to_csv(TABLES_DIR / "bulk_subtype_stats.csv", index=False)

    corr = pd.DataFrame(rows_corr)
    corr["spearman_p_BH"] = bh_fdr(corr["spearman_p"])
    corr.to_csv(TABLES_DIR / "bulk_immune_correlations.csv", index=False)

    # per-sample annotation table for reproducibility
    ann = []
    for name, expr, st, mk, sc in [
            ("George2015", george, subtype_g, marker_g, scores_g),
            ("GSE60052", tum, subtype_j, marker_j, scores_j)]:
        d = pd.DataFrame({
            "cohort": name, "sample": expr.index,
            "subtype_nmf": st.values,
            "subtype_marker_rule": mk.values,
            "TACSTD2_log2": expr["TACSTD2"].values,
            "CLDN4_log2": expr["CLDN4"].values,
            "Teff_score": sc["Teff_score"].values,
        })
        ann.append(d)
    pd.concat(ann).to_csv(TABLES_DIR / "bulk_sample_annotations.csv", index=False)
    print("[done] bulk tables written to", TABLES_DIR)


if __name__ == "__main__":
    main()

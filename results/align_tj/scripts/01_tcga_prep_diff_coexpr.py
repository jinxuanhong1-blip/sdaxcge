"""TCGA NSCLC (LUAD+LUSC): build TROP2-high/low groups, differential expression,
and TACSTD2 coexpression with user genes / TFs.

Input : data/TCGA.{LUAD,LUSC}.HiSeqV2.gz  (log2(norm_count+1), gene x sample)
Output: tables/tcga_diff_high_vs_low.tsv     (per-gene high-vs-low stats + Spearman r vs TACSTD2)
        tables/tcga_rank_tstat.rnk           (preranked list for GSEA, by high-vs-low t-stat)
        tables/tcga_coexpr_targets.tsv       (TACSTD2 vs user TJ genes + TFs, Spearman & Pearson)
        tables/tcga_group_summary.tsv        (n samples per cohort / group)
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import common as C

LOWQ, HIGHQ = 1 / 3, 2 / 3  # tertile split on TACSTD2


def load_nsclc():
    luad = C.read_xena_matrix(f"{C.DATA}/TCGA.LUAD.HiSeqV2.gz")
    lusc = C.read_xena_matrix(f"{C.DATA}/TCGA.LUSC.HiSeqV2.gz")
    luad = luad[C.tumor_samples(luad.columns)]
    lusc = lusc[C.tumor_samples(lusc.columns)]
    genes = luad.index.intersection(lusc.index)
    luad, lusc = luad.loc[genes], lusc.loc[genes]
    cohort = pd.Series(["LUAD"] * luad.shape[1] + ["LUSC"] * lusc.shape[1],
                       index=list(luad.columns) + list(lusc.columns))
    expr = pd.concat([luad, lusc], axis=1)
    return expr, cohort


def main():
    expr, cohort = load_nsclc()
    print(f"NSCLC tumors: {expr.shape[1]} (LUAD {int((cohort=='LUAD').sum())}, "
          f"LUSC {int((cohort=='LUSC').sum())}); genes {expr.shape[0]}")

    trop2 = expr.loc[C.TROP2]
    lo_thr, hi_thr = trop2.quantile(LOWQ), trop2.quantile(HIGHQ)
    high = trop2[trop2 >= hi_thr].index
    low = trop2[trop2 <= lo_thr].index
    print(f"TACSTD2 high n={len(high)}  low n={len(low)}  (thr low={lo_thr:.3f} high={hi_thr:.3f})")

    # group summary
    gs = []
    for grp, idx in [("high", high), ("low", low)]:
        c = cohort.loc[idx].value_counts()
        gs.append({"group": grp, "n": len(idx),
                   "LUAD": int(c.get("LUAD", 0)), "LUSC": int(c.get("LUSC", 0)),
                   "mean_TACSTD2": float(trop2.loc[idx].mean())})
    pd.DataFrame(gs).to_csv(f"{C.TABLES}/tcga_group_summary.tsv", sep="\t", index=False)

    Xh = expr[high].values
    Xl = expr[low].values
    # Welch t-test high vs low per gene (data already log2)
    t, p = stats.ttest_ind(Xh, Xl, axis=1, equal_var=False)
    log2fc = Xh.mean(axis=1) - Xl.mean(axis=1)
    # Spearman corr of each gene with continuous TACSTD2 across all tumors
    ranks_trop2 = stats.rankdata(trop2.values)
    exprvals = expr.values
    rho = np.array([stats.spearmanr(exprvals[i], trop2.values).correlation
                    for i in range(exprvals.shape[0])])

    res = pd.DataFrame({
        "gene": expr.index,
        "log2FC_high_minus_low": log2fc,
        "t_stat": t,
        "p_ttest": p,
        "spearman_r_vs_TACSTD2": rho,
    }).set_index("gene")
    res["padj_ttest"] = multipletests(res["p_ttest"].fillna(1), method="fdr_bh")[1]
    res = res.sort_values("t_stat", ascending=False)
    res.to_csv(f"{C.TABLES}/tcga_diff_high_vs_low.tsv", sep="\t")

    # preranked list by t-stat (drop TACSTD2 itself to avoid trivial self-drive not needed for GSEA)
    rnk = res[["t_stat"]].dropna().sort_values("t_stat", ascending=False)
    rnk.to_csv(f"{C.TABLES}/tcga_rank_tstat.rnk", sep="\t", header=False)

    # coexpression of TACSTD2 with target genes
    targets = C.USER_TJ_GENES + C.USER_TFS
    rows = []
    for g in targets:
        if g not in expr.index:
            rows.append({"gene": g, "present": False}); continue
        sr, sp = stats.spearmanr(expr.loc[g].values, trop2.values)
        pr, pp = stats.pearsonr(expr.loc[g].values, trop2.values)
        rows.append({"gene": g, "present": True, "class":
                     "TJ" if g in C.USER_TJ_GENES else "TF",
                     "spearman_r": sr, "spearman_p": sp,
                     "pearson_r": pr, "pearson_p": pp})
    co = pd.DataFrame(rows)
    co.to_csv(f"{C.TABLES}/tcga_coexpr_targets.tsv", sep="\t", index=False)
    print(co.to_string(index=False))

    # Also: coexpression among the TFs and CLDN4 (for the ELF3/GRHL1/KLF4/TFAP2A vs CLDN4 part)
    panel = [C.TROP2, "CLDN4"] + C.USER_TFS
    panel = [g for g in panel if g in expr.index]
    cm = expr.loc[panel].T.corr(method="spearman")
    cm.to_csv(f"{C.TABLES}/tcga_tf_cldn4_trop2_spearman_matrix.tsv", sep="\t")
    print("\nTF/CLDN4/TROP2 Spearman matrix:\n", cm.round(3).to_string())


if __name__ == "__main__":
    main()

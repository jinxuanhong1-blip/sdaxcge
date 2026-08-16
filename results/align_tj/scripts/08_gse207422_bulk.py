"""GSE207422 bulk RNA-seq (log2 TPM): TACSTD2 vs user TJ/TF genes.

Small n (~29 samples). Reported as a same-study bulk check, not an
independent cohort.

Outputs: tables/gse207422_bulk_coexpr_targets.tsv
"""
import pandas as pd
from scipy import stats
import common as C

MAT = f"{C.DATA}/GSE207422_bulk_log2TPM.txt.gz"


def main():
    expr = pd.read_csv(MAT, sep="\t", index_col=0)
    print("bulk", expr.shape)
    trop = expr.loc["TACSTD2"] if "TACSTD2" in expr.index else None
    if trop is None:
        raise SystemExit("TACSTD2 missing")
    rows = []
    for g in C.USER_TJ_GENES + C.USER_TFS:
        if g not in expr.index:
            rows.append({"gene": g, "present": False})
            continue
        r, p = stats.spearmanr(expr.loc[g], trop)
        r4 = p4 = None
        if "CLDN4" in expr.index and g != "CLDN4":
            r4, p4 = stats.spearmanr(expr.loc[g], expr.loc["CLDN4"])
        rows.append({"gene": g, "present": True,
                     "class": "TJ" if g in C.USER_TJ_GENES else "TF",
                     "n": expr.shape[1],
                     "spearman_r_vs_TACSTD2": r, "spearman_p_vs_TACSTD2": p,
                     "spearman_r_vs_CLDN4": r4, "spearman_p_vs_CLDN4": p4})
    out = pd.DataFrame(rows)
    out.to_csv(f"{C.TABLES}/gse207422_bulk_coexpr_targets.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

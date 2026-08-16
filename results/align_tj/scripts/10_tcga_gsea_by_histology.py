"""LUAD-only and LUSC-only TROP2-high vs low preranked GSEA (same libraries
as the pooled run). Tests whether keratinization/barrier/TJ are only a
squamous-mix effect.

Outputs: tables/tcga_{luad,lusc}_gsea_themes.tsv
         tables/tcga_{luad,lusc}_rank_tstat.rnk
"""
import re
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp
import common as C

LIBS = ["GO_Biological_Process_2021", "GO_Cellular_Component_2021",
        "MSigDB_Hallmark_2020", "KEGG_2021_Human"]


def load_one(path):
    expr = C.read_xena_matrix(path)
    expr = expr[C.tumor_samples(expr.columns)]
    return expr


def rank_and_gsea(expr, tag, gene_sets):
    trop2 = expr.loc[C.TROP2]
    lo, hi = trop2.quantile(1 / 3), trop2.quantile(2 / 3)
    high, low = trop2[trop2 >= hi].index, trop2[trop2 <= lo].index
    print(f"{tag}: n={expr.shape[1]} high={len(high)} low={len(low)}")
    t, _p = stats.ttest_ind(expr[high].values, expr[low].values,
                            axis=1, equal_var=False)
    rnk = pd.Series(t, index=expr.index, name="t").dropna().sort_values(ascending=False)
    rnk.to_csv(f"{C.TABLES}/tcga_{tag}_rank_tstat.rnk", sep="\t", header=False)
    pre = gp.prerank(rnk=rnk, gene_sets=gene_sets, min_size=5, max_size=1500,
                     permutation_num=1000, seed=7, threads=2,
                     outdir=None, no_plot=True)
    res = pre.res2d.copy()
    res.columns = [c.strip() for c in res.columns]
    text = res["Term"].str.lower()
    rows = []
    for theme, kws in C.THEME_KEYWORDS.items():
        mask = pd.Series(False, index=res.index)
        for kw in kws:
            mask |= text.str.contains(re.escape(kw))
        sub = res[mask].copy()
        sub.insert(0, "theme", theme)
        sub.insert(0, "cohort", tag.upper())
        rows.append(sub)
    themes = pd.concat(rows, ignore_index=True)
    themes.to_csv(f"{C.TABLES}/tcga_{tag}_gsea_themes.tsv", sep="\t", index=False)
    keep = [c for c in ["cohort", "theme", "Term", "NES", "NOM p-val", "FDR q-val"]
            if c in themes.columns]
    top = themes.sort_values("NES", ascending=False).groupby("theme").head(2)
    print(top[keep].to_string(index=False))
    return themes


def main():
    gene_sets = {}
    for lib in LIBS:
        d = gp.get_library(name=lib, organism="Human")
        for term, genes in d.items():
            gene_sets[f"{lib} :: {term}"] = genes
        print("loaded", lib, len(d))
    luad = load_one(f"{C.DATA}/TCGA.LUAD.HiSeqV2.gz")
    lusc = load_one(f"{C.DATA}/TCGA.LUSC.HiSeqV2.gz")
    rank_and_gsea(luad, "luad", gene_sets)
    rank_and_gsea(lusc, "lusc", gene_sets)


if __name__ == "__main__":
    main()

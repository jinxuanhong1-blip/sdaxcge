"""Preranked GSEA on TCGA NSCLC TROP2-high vs low (ranking = Welch t-stat).

Gene-set libraries (Enrichr / MSigDB via gseapy):
  GO_Biological_Process_2021, GO_Cellular_Component_2021,
  MSigDB_Hallmark_2020, KEGG_2021_Human
We report all terms matching the user's four themes:
  keratinization, skin barrier, tight junction, EMT.

Output: tables/tcga_gsea_full.tsv        (all terms)
        tables/tcga_gsea_themes.tsv       (theme-matched terms)
"""
import re
import pandas as pd
import gseapy as gp
import common as C

LIBS = ["GO_Biological_Process_2021", "GO_Cellular_Component_2021",
        "MSigDB_Hallmark_2020", "KEGG_2021_Human"]


def load_libraries():
    sets = {}
    for lib in LIBS:
        try:
            d = gp.get_library(name=lib, organism="Human")
            for term, genes in d.items():
                sets[f"{lib} :: {term}"] = genes
            print(f"loaded {lib}: {len(d)} terms")
        except Exception as e:
            print(f"WARN could not load {lib}: {e}")
    return sets


def main():
    rnk = pd.read_csv(f"{C.TABLES}/tcga_rank_tstat.rnk", sep="\t", header=None,
                      index_col=0)[1]
    rnk = rnk.sort_values(ascending=False)
    gene_sets = load_libraries()
    print(f"total gene sets: {len(gene_sets)}")

    pre = gp.prerank(rnk=rnk, gene_sets=gene_sets, min_size=5, max_size=1500,
                     permutation_num=1000, seed=7, threads=4,
                     outdir=None, no_plot=True)
    res = pre.res2d.copy()
    # normalize column names
    res.columns = [c.strip() for c in res.columns]
    res.to_csv(f"{C.TABLES}/tcga_gsea_full.tsv", sep="\t", index=False)

    term_col = "Term"
    text = res[term_col].str.lower()
    theme_rows = []
    for theme, kws in C.THEME_KEYWORDS.items():
        mask = pd.Series(False, index=res.index)
        for kw in kws:
            mask |= text.str.contains(re.escape(kw))
        sub = res[mask].copy()
        sub.insert(0, "theme", theme)
        theme_rows.append(sub)
    themes = pd.concat(theme_rows, ignore_index=True)
    themes.to_csv(f"{C.TABLES}/tcga_gsea_themes.tsv", sep="\t", index=False)

    cols = [c for c in ["theme", "Term", "NES", "NOM p-val", "FDR q-val",
                        "Lead_genes"] if c in themes.columns]
    print(themes[cols].to_string(index=False)[:6000])


if __name__ == "__main__":
    main()

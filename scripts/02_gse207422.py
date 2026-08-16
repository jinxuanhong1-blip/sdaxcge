#!/usr/bin/env python3
"""GSE207422 NSCLC pre-treatment bulk RNA-seq (author-provided log2TPM).

n=24 pre-treatment biopsies (adeno + squamous + NOS). Small-n is reported
as-is; no histology subsetting (that would drop power below a usable tertile).
"""
import os
import pandas as pd
import lib_analysis as L

RAW = "data/raw"
OUT = "results/claim_A8A9/gse207422"
os.makedirs(OUT, exist_ok=True)


def main():
    expr = pd.read_csv(f"{RAW}/GSE207422_bulk_log2TPM.txt.gz", sep="\t", index_col=0)
    expr = expr.groupby(level=0).max()
    meta = pd.read_excel(f"{RAW}/GSE207422_bulk_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"]).copy()
    meta["Sample"] = meta["Sample"].astype(str)
    keep = [c for c in expr.columns if c in set(meta["Sample"])]
    expr = expr[keep]
    print("GSE207422 bulk matrix:", expr.shape)
    print("  pathology counts:")
    print(meta.set_index("Sample").loc[keep, "Pathology"].value_counts().to_string())

    if L.TROP2 not in expr.index:
        raise SystemExit(f"{L.TROP2} missing")
    trop2 = expr.loc[L.TROP2]
    trop2.to_csv(f"{OUT}/tacstd2_per_sample.csv")

    stats = L.per_gene_stats(expr, trop2, min_frac_expressed=0.10)
    stats.to_csv(f"{OUT}/per_gene_stats.csv")
    print("  genes tested:", stats.shape[0], "n_high/n_low/n_total =",
          int(stats["n_high"].iloc[0]), int(stats["n_low"].iloc[0]), int(stats["n_total"].iloc[0]))

    panel = stats.reindex(L.PANEL)
    panel.to_csv(f"{OUT}/panel_stats.csv")
    print(panel[["log2FC_high_vs_low", "spearman_rho", "welch_fdr", "up_in_trop2_high"]])

    gsea = L.run_prerank_gsea(stats, OUT)
    print(gsea[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))
    print("GSE207422 done.")


if __name__ == "__main__":
    main()

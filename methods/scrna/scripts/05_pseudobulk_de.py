"""PRIMARY differential expression = PSEUDOBULK (decoupler + pydeseq2).

Aggregate raw counts per sample within a chosen population, then fit a
sample-level model. This is the correct unit of replication for a
patient-level label (MPR/RECIST). NOT cells-as-replicates.
主要差异分析 = 伪bulk：在选定细胞群内按样本聚合原始counts，再做样本层建模。

Usage:
  python 05_pseudobulk_de.py annotated.h5ad OUTDIR \
      --population Epithelial --group_key path_response \
      --group_a MPR --group_b NMPR --covariates histology
"""
import argparse, os
import numpy as np
import pandas as pd
import scanpy as sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad"); ap.add_argument("outdir")
    ap.add_argument("--sample_key", default="sample")
    ap.add_argument("--population_key", default="lineage")
    ap.add_argument("--population", default="Epithelial")
    ap.add_argument("--group_key", default="path_response")
    ap.add_argument("--group_a", default="MPR")
    ap.add_argument("--group_b", default="NMPR")
    ap.add_argument("--covariates", nargs="*", default=[])
    ap.add_argument("--min_cells", type=int, default=10)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    import decoupler as dc
    adata = sc.read_h5ad(args.h5ad)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]      # decoupler needs raw counts

    sub = adata[adata.obs[args.population_key] == args.population].copy()
    print(f"{sub.n_obs} {args.population} cells")

    # sum counts per sample -> pseudobulk AnnData (samples x genes)
    pdata = dc.get_pseudobulk(
        sub, sample_col=args.sample_key, groups_col=None,
        mode="sum", min_cells=args.min_cells, min_counts=1000,
    )
    pdata.obs = pdata.obs.join(
        sub.obs.drop_duplicates(args.sample_key).set_index(args.sample_key)[
            [args.group_key] + args.covariates],
        on=args.sample_key)

    # keep only the two groups of interest
    keep = pdata.obs[args.group_key].isin([args.group_a, args.group_b])
    pdata = pdata[keep].copy()
    print(pdata.obs[[args.group_key] + args.covariates])

    # edgeR filterByExpr-equivalent
    dc.pp.filter_by_expr(pdata, group=pdata.obs[args.group_key]) \
        if hasattr(dc, "pp") else None

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    design = "~ " + " + ".join(args.covariates + [args.group_key]) \
        if args.covariates else f"~ {args.group_key}"
    dds = DeseqDataSet(adata=pdata, design=design, refit_cooks=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=[args.group_key, args.group_a, args.group_b])
    stat.summary()
    res = stat.results_df.sort_values("padj")
    out = os.path.join(args.outdir,
                       f"de_{args.population}_{args.group_a}_vs_{args.group_b}.csv")
    res.to_csv(out)
    print(f"DE table -> {out}\n", res.head(20))

    # export the pseudobulk matrix + design for edgeR/DESeq2 in R too
    pd.DataFrame(pdata.X, index=pdata.obs_names, columns=pdata.var_names)\
      .to_csv(os.path.join(args.outdir, "pseudobulk_counts.csv"))
    pdata.obs.to_csv(os.path.join(args.outdir, "pseudobulk_design.csv"))

if __name__ == "__main__":
    main()

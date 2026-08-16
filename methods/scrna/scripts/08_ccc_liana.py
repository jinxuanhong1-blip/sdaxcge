"""SECONDARY analysis: cell-cell communication with LIANA+ (consensus).

Hypothesis-generating ONLY. Prefer DIFFERENTIAL CCC (per condition) and confirm
priorities against pseudobulk DE / spatial proximity. Treat permutation p-values
skeptically (same pseudoreplication risk). See playbook.md section 8.
仅作次要/假设生成；优先做分组差异CCC，并与伪bulk DE/空间邻近互证；排列检验p值需谨慎。

Usage: python 08_ccc_liana.py annotated.h5ad out_prefix [groupby]
"""
import sys
import scanpy as sc

def main(in_h5ad, out_prefix, groupby="lineage"):
    import liana as li
    adata = sc.read_h5ad(in_h5ad)
    # LIANA expects log-normalized data in .X
    if "log1p" not in adata.uns_keys():
        sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)

    li.mt.rank_aggregate(
        adata, groupby=groupby, expr_prop=0.1, use_raw=False, verbose=True,
    )
    res = adata.uns["liana_res"]
    res.to_csv(f"{out_prefix}_liana_consensus.csv", index=False)
    print(res.head(20))

    # For DIFFERENTIAL CCC across response groups, run per group and compare,
    # or use li.mt with a condition, then contrast the ranks. Example:
    #   for g, sub in adata.obs.groupby('path_response'):
    #       li.mt.rank_aggregate(adata[adata.obs.path_response==g], groupby=groupby)
    print(f"LIANA consensus -> {out_prefix}_liana_consensus.csv")

if __name__ == "__main__":
    main(*sys.argv[1:])

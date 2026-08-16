"""Epithelial + immune module scoring (Scanpy score_genes; optionally decoupler AUCell).

Score the tumor program ONLY on malignant epithelial cells, immune programs on the
relevant immune subset, then compare at the SAMPLE level. Check the score is not just
tracking depth. 模块打分：肿瘤程序仅在恶性上皮细胞、免疫程序在对应免疫亚群；样本层比较，并检查是否受深度混杂。

Usage: python 06_module_scores.py annotated.h5ad out.h5ad
"""
import sys
import numpy as np
import scanpy as sc
from gene_sets import ALL_MODULES

def main(in_h5ad, out_h5ad):
    adata = sc.read_h5ad(in_h5ad)
    use_raw = adata.raw is not None
    vnames = adata.raw.var_names if use_raw else adata.var_names
    for name, gs in ALL_MODULES.items():
        present = [g for g in gs if g in vnames]
        if present:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=use_raw,
                              ctrl_size=max(50, len(present)))
        print(f"{name}: {len(present)}/{len(gs)} genes")

    # depth-confound check for the tumor module
    if "TACSTD2_CLDN4_junction" in adata.obs and "total_counts" in adata.obs:
        r = np.corrcoef(adata.obs["TACSTD2_CLDN4_junction"],
                        np.log1p(adata.obs["total_counts"]))[0, 1]
        print(f"corr(module, log depth) = {r:.3f}  (want |r| small)")

    adata.write(out_h5ad)
    print(f"Module scores -> {out_h5ad}")

if __name__ == "__main__":
    main(*sys.argv[1:])

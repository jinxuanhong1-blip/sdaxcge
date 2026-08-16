"""Per-sample QC + doublet detection (Scanpy + scrublet/scDblFinder-style).

Run per sample BEFORE integration. Prefer R's scDblFinder for the doublet call;
this Python template uses Scanpy's built-in scrublet as a portable fallback.
整合前、每样本单独运行；双细胞检测优先用 R 的 scDblFinder，此处用 scrublet 作可移植回退。

Usage: python 02_qc_doublets.py sample.h5ad out.h5ad
"""
import sys
import numpy as np
import scanpy as sc

def mad_outlier(x, nmads=3.0, log=False, upper=True, lower=True):
    v = np.log1p(x) if log else np.asarray(x, float)
    med, mad = np.median(v), np.median(np.abs(v - np.median(v))) * 1.4826
    hi = v > med + nmads * mad if upper else np.zeros(len(v), bool)
    lo = v < med - nmads * mad if lower else np.zeros(len(v), bool)
    return hi | lo

def main(in_h5ad, out_h5ad):
    adata = sc.read_h5ad(in_h5ad)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True,
                               percent_top=None, log1p=False)

    # per-sample MAD thresholds (this file already is one sample)
    bad = (mad_outlier(adata.obs["total_counts"], log=True, upper=False)
           | mad_outlier(adata.obs["n_genes_by_counts"], log=True, upper=False)
           | mad_outlier(adata.obs["pct_counts_mt"], lower=False))
    print(f"QC discards: {int(bad.sum())} / {adata.n_obs}")
    adata = adata[~bad].copy()

    # doublets (scrublet); for production use R scDblFinder
    try:
        sc.pp.scrublet(adata, random_state=0)
        n_db = int(adata.obs["predicted_doublet"].sum())
        print(f"Doublets: {n_db} / {adata.n_obs}")
        adata = adata[~adata.obs["predicted_doublet"]].copy()
    except Exception as exc:  # noqa: BLE001
        print(f"scrublet skipped: {exc}")

    adata.write(out_h5ad)
    print(f"Clean singlets -> {out_h5ad}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

"""Join MPR / RECIST clinical labels at the SAMPLE level and test associations.

Keep MPR and RECIST as SEPARATE columns (they are not interchangeable). Tests are at
the SAMPLE level (never cells-as-replicates). Correct across signatures. Be honest with n.
在样本层对接 MPR 与 RECIST（二者不可互换、分开保存），样本层检验，多重校正，诚实报告小样本。

Usage:
  python 09_join_clinical.py annotated.h5ad sample_metadata.csv out_prefix \
      --readouts TACSTD2_CLDN4_junction CD8_cytotoxic TLS_core CXCL13 \
      --population Epithelial --group_key path_response
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad"); ap.add_argument("clinical_csv"); ap.add_argument("out_prefix")
    ap.add_argument("--sample_key", default="sample")
    ap.add_argument("--readouts", nargs="+", required=True)
    ap.add_argument("--population_key", default="lineage")
    ap.add_argument("--population", default=None,
                    help="restrict readout averaging to this lineage (e.g. Epithelial)")
    ap.add_argument("--group_key", default="path_response")
    ap.add_argument("--group_a", default="MPR"); ap.add_argument("--group_b", default="NMPR")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    obs = adata.obs
    if args.population and args.population_key in obs:
        obs = obs[obs[args.population_key] == args.population]

    # sample-level means of each readout
    per_sample = obs.groupby(args.sample_key)[args.readouts].mean()

    clin = pd.read_csv(args.clinical_csv).set_index(args.sample_key)
    tab = per_sample.join(clin, how="inner")
    tab.to_csv(f"{args.out_prefix}_sample_readouts.csv")

    # sample-level Wilcoxon between the two response groups, BH across readouts
    a = tab[tab[args.group_key] == args.group_a]
    b = tab[tab[args.group_key] == args.group_b]
    rows = []
    for r in args.readouts:
        x, y = a[r].dropna(), b[r].dropna()
        if len(x) >= 2 and len(y) >= 2:
            u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            rows.append({"readout": r, "n_a": len(x), "n_b": len(y),
                         "median_a": x.median(), "median_b": y.median(),
                         "U": u, "p": p})
    res = pd.DataFrame(rows)
    if len(res):
        # Benjamini-Hochberg
        res = res.sort_values("p").reset_index(drop=True)
        m = len(res)
        res["padj"] = (res["p"] * m / (res.index + 1)).clip(upper=1.0)[::-1].cummin()[::-1]
    res.to_csv(f"{args.out_prefix}_association_tests.csv", index=False)
    print(res.to_string(index=False))
    print(f"\nn={len(a)} {args.group_a} vs n={len(b)} {args.group_b} samples "
          f"-> underpowered; treat as EXPLORATORY and replicate externally.")

if __name__ == "__main__":
    main()

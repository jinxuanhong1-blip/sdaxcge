"""Chan et al. 2021 SCLC single-cell atlas: TACSTD2/CLDN4 in tumor cells.

Data: cellxgene "Combined samples" h5ad (147,137 cells; HTAN MSK,
dbGaP phs002371). Raw counts are taken from .raw.X and normalized to
log1p(CP10K) for the extracted gene panel only (backed mode keeps memory low).

Analyses (donor-level to avoid pseudoreplication):
  - Expression summary of TACSTD2/CLDN4 by coarse cell type (SCLC samples)
  - Donor pseudobulk (mean log1p CP10K in epithelial/tumor cells of SCLC
    samples); approximate subtype per donor via ASCL1/NEUROD1/POU2F3 z-max
    (triple-low = SCLC-I-like, as in the bulk analysis)
  - Kruskal-Wallis across subtypes; Mann-Whitney U ICI-exposed vs ICI-naive
"""

import numpy as np
import pandas as pd
from scipy import stats
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from common import (DATA_DIR, TABLES_DIR, FIGURES_DIR, SUBTYPE_TFS,
                    ensure_dirs, bh_fdr, zscore)

GENES = ["TACSTD2", "CLDN4", "ASCL1", "NEUROD1", "POU2F3", "YAP1",
         "EPCAM", "CD8A"]
MIN_TUMOR_CELLS = 50


def extract(h5ad_path):
    a = ad.read_h5ad(h5ad_path, backed="r")
    fn = a.raw.var["feature_name"].astype(str)
    idx = [int(np.where(fn == g)[0][0]) for g in GENES]
    raw = a.raw.X[:, idx]
    raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)
    # per-cell total raw counts, streamed in chunks (obs['libsize'] is not
    # the raw count total in this export)
    totals = np.concatenate([
        np.asarray(a.raw.X[i:i + 20000].sum(axis=1)).ravel()
        for i in range(0, a.n_obs, 20000)])
    cp10k = np.log1p(raw / totals[:, None] * 1e4)
    df = pd.DataFrame(cp10k, columns=GENES, index=a.obs_names)
    for c in ["donor_id", "histo", "treatment", "cell_type_coarse",
              "cell_type_fine", "tissue", "procedure"]:
        df[c] = a.obs[c].values
    return df


def main():
    ensure_dirs()
    df = extract(DATA_DIR / "chan2021_combined.h5ad")
    sclc = df[df["histo"] == "SCLC"].copy()
    print(f"SCLC cells: {len(sclc)}, donors: {sclc['donor_id'].nunique()}")

    # cell-type-level summary within SCLC samples
    summ = (sclc.groupby("cell_type_coarse", observed=True)
            .agg(n_cells=("TACSTD2", "size"),
                 TACSTD2_mean=("TACSTD2", "mean"),
                 TACSTD2_pct_pos=("TACSTD2", lambda x: 100 * (x > 0).mean()),
                 CLDN4_mean=("CLDN4", "mean"),
                 CLDN4_pct_pos=("CLDN4", lambda x: 100 * (x > 0).mean())))
    summ.to_csv(TABLES_DIR / "scrna_celltype_summary.csv")
    print(summ)

    # donor pseudobulk in epithelial (tumor) cells
    epi = sclc[sclc["cell_type_coarse"] == "Epithelial"]
    counts = epi.groupby("donor_id", observed=True).size()
    keep = counts[counts >= MIN_TUMOR_CELLS].index
    pb = (epi[epi["donor_id"].isin(keep)]
          .groupby("donor_id", observed=True)[GENES].mean())
    treat = (sclc.drop_duplicates("donor_id").set_index("donor_id")["treatment"]
             .reindex(pb.index).astype(str))
    pb["ici_exposed"] = treat.str.contains("Immunotherapy")
    pb["n_tumor_cells"] = counts.reindex(pb.index)

    z = zscore(pb[SUBTYPE_TFS])
    mapping = {"ASCL1": "SCLC-A", "NEUROD1": "SCLC-N", "POU2F3": "SCLC-P"}
    pb["subtype_approx"] = [
        mapping[r.idxmax()] if r.max() > 0 else "SCLC-I-like"
        for _, r in z.iterrows()]
    pb.to_csv(TABLES_DIR / "scrna_donor_pseudobulk.csv")
    print(pb["subtype_approx"].value_counts())

    rows = []
    for gene in ["TACSTD2", "CLDN4"]:
        groups = [g[gene].values for _, g in pb.groupby("subtype_approx")
                  if len(g) >= 2]
        if len(groups) >= 3:
            h, p = stats.kruskal(*groups)
            rows.append({"gene": gene, "test": "Kruskal-Wallis subtype",
                         "n": len(pb), "stat": h, "p": p})
        ici = pb[pb["ici_exposed"]][gene]
        naive = pb[~pb["ici_exposed"]][gene]
        u, p = stats.mannwhitneyu(ici, naive, alternative="two-sided")
        rows.append({"gene": gene, "test": "MWU ICI-exposed vs naive",
                     "n": len(pb), "stat": u, "p": p,
                     "median_ici": ici.median(), "median_naive": naive.median(),
                     "n_ici": len(ici), "n_naive": len(naive)})
    res = pd.DataFrame(rows)
    res["p_BH"] = bh_fdr(res["p"])
    res.to_csv(TABLES_DIR / "scrna_donor_stats.csv", index=False)
    print(res)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    order = ["SCLC-A", "SCLC-N", "SCLC-P", "SCLC-I-like"]
    order = [o for o in order if o in set(pb["subtype_approx"])]
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
        sns.boxplot(data=pb, x="subtype_approx", y=gene, order=order,
                    hue="subtype_approx", legend=False, palette="Set2",
                    showfliers=False, ax=ax)
        sns.stripplot(data=pb, x="subtype_approx", y=gene, order=order,
                      color="k", size=4, ax=ax)
        ax.set_ylabel(f"{gene} (donor mean log1p CP10K, tumor cells)")
        ax.set_title(f"Chan 2021 atlas: {gene}")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "scrna_donor_pseudobulk.png", dpi=170)
    plt.close(fig)
    print("[done] scRNA tables/figures written")


if __name__ == "__main__":
    main()

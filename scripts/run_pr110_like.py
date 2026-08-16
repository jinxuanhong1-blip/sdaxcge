"""PR #110-like sensitivity on the same Visium objects.

This is NOT a re-implementation of PR #110 (that lives on a separate branch).
It is a one-table comparison using THIS repo's marker lists:

  index spots = epithelial-broad ≥ section median  AND  TB score < section Q3
  stat        = partial Spearman of index CLDN4 vs ring-k mean TB,
                controlling for the index epithelial-broad score
  sections    = E-MTAB-13530 tumor (P*_T*) and GSE189487 (all 6)

If this also collapses toward 0, the raw negative neighborhood ρ is
composition / spatial autocorrelation, not CLDN4-specific exclusion.
"""
import os
import glob
import warnings
import numpy as np
import pandas as pd
import scanpy as sc

import common
from neighborhood import hex_adj, rings_from_adj, neighbor_means, partial_spearman, spearman
from run_visium_neighborhood import load_emtab, load_gse189, emtab_meta

warnings.simplefilter("ignore")
sc.settings.verbosity = 0
OUT = "results/claim_B6/tables"


def one_section(a, sample, dataset):
    rows = a.obs["array_row"].values
    cols = a.obs["array_col"].values
    keys = list(zip(rows.astype(int), cols.astype(int)))
    adj = hex_adj(rows, cols)
    rmap = rings_from_adj(adj, max_k=3)
    epi = a.obs["epi_broad_score"].values
    tb = a.obs["tb_score"].values
    cldn4 = a.obs["cldn4"].values
    epi_thr = np.nanmedian(epi)
    tb_q3 = np.nanquantile(tb, 0.75)
    index = (epi >= epi_thr) & (tb < tb_q3)
    recs = []
    for k in range(1, 4):
        tb_nb = neighbor_means(keys, tb, rmap, k)
        # restrict to index spots that have a neighborhood value
        m = index & np.isfinite(tb_nb) & np.isfinite(cldn4) & np.isfinite(epi)
        rho, p, n = spearman(cldn4[m], tb_nb[m])
        pr, pp, pn = partial_spearman(cldn4[m], tb_nb[m], epi[m])
        recs.append(dict(
            sample=sample, dataset=dataset, ring=k,
            n_qc=int(a.n_obs), n_index=int(index.sum()), n_used=n,
            rho_raw=rho, p_raw=p,
            rho_partialEpi=pr, p_partialEpi=pp,
        ))
    return recs


def main():
    meta = emtab_meta()
    rows = []
    for h5 in sorted(glob.glob("data/emtab13530/*-filtered_feature_bc_matrix.h5")):
        sid = os.path.basename(h5).replace("-filtered_feature_bc_matrix.h5", "")
        if meta.get(sid) != "tumor":
            continue
        a = load_emtab(sid)
        recs = one_section(a, sid, "E-MTAB-13530")
        rows.extend(recs)
        print(f"{sid} ring1 n_index={recs[0]['n_index']} "
              f"raw={recs[0]['rho_raw']:+.3f} partial={recs[0]['rho_partialEpi']:+.3f}")

    for mtx in sorted(glob.glob("data/gse189487/raw/*_matrix.mtx.gz")):
        pref = os.path.basename(mtx).rsplit("_", 1)[0]
        sid = pref.split("_", 1)[1]
        a = load_gse189(pref)
        recs = one_section(a, sid, "GSE189487")
        rows.extend(recs)
        print(f"{sid} ring1 n_index={recs[0]['n_index']} "
              f"raw={recs[0]['rho_raw']:+.3f} partial={recs[0]['rho_partialEpi']:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/visium_pr110like_rings.csv", index=False)
    print("wrote", f"{OUT}/visium_pr110like_rings.csv")


if __name__ == "__main__":
    main()

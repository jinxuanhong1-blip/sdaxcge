"""CLAIM B6 on E-MTAB-13530 (10x Visium, human NSCLC)."""
import os
import glob
import warnings
import numpy as np
import pandas as pd
import scanpy as sc

import common

warnings.simplefilter("ignore")
sc.settings.verbosity = 0

DATA = "data/emtab13530"
OUT = "results/claim_B6"
os.makedirs(f"{OUT}/tables", exist_ok=True)


def load_sdrf():
    df = pd.read_csv(f"{DATA}/E-MTAB-13530.sdrf.txt", sep="\t")
    m = {}
    for _, r in df.iterrows():
        s = r["Source Name"]
        site = str(r.get("Characteristics[sampling site]", ""))
        dis = str(r.get("Characteristics[disease]", ""))
        indiv = str(r.get("Characteristics[individual]", ""))
        if "tumor" == site.strip():
            grp = "tumor"
        elif "adjacent" in site:
            grp = "adjacent_normal"
        else:
            grp = "healthy"
        m[s] = dict(sampling_site=site, disease=dis, patient=indiv, group=grp)
    return m


def main():
    meta = load_sdrf()
    h5s = sorted(glob.glob(f"{DATA}/*-filtered_feature_bc_matrix.h5"))
    rows = []
    for h5 in h5s:
        sid = os.path.basename(h5).replace("-filtered_feature_bc_matrix.h5", "")
        a = sc.read_10x_h5(h5)
        a.var_names_make_unique()
        sc.pp.filter_cells(a, min_genes=200)
        sc.pp.filter_genes(a, min_cells=3)
        if a.n_obs < 50:
            print(f"skip {sid}: only {a.n_obs} spots pass QC")
            continue
        a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        info = meta.get(sid, {})
        extra = dict(dataset="E-MTAB-13530", platform="Visium",
                     group=info.get("group", "NA"),
                     patient=info.get("patient", "NA"),
                     disease=info.get("disease", "NA"))
        row = common.analyze_sample(a, sid, extra)
        rows.append(row)
        print(f"{sid:10s} [{extra['group']:15s}] n={row['n_spots']:5d} "
              f"rho(epi,TB)={row['rho_epi_vs_TB']:+.3f} "
              f"rho(epi,T)={row['rho_epi_vs_T']:+.3f} "
              f"rho(epi,B)={row['rho_epi_vs_B']:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/tables/emtab13530_per_section.csv", index=False)
    print(f"\nwrote {OUT}/tables/emtab13530_per_section.csv ({len(df)} sections)")
    return df


if __name__ == "__main__":
    main()

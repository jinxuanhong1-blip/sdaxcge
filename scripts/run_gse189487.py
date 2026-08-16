"""CLAIM B6 same-spot test on GSE189487 (10x Visium, early LUAD AIS/MIA/IAC)."""
import os
import glob
import warnings
import pandas as pd
import scanpy as sc

import common

warnings.simplefilter("ignore")
sc.settings.verbosity = 0

DATA = "data/gse189487/raw"
OUT = "results/claim_B6"
os.makedirs(f"{OUT}/tables", exist_ok=True)

STAGE = {
    "TD1": ("IAC", "GSM5702473"),
    "TD2": ("IAC", "GSM5702474"),
    "TD3": ("MIA", "GSM5702475"),
    "TD5": ("AIS", "GSM5702476"),
    "TD6": ("MIA", "GSM5702477"),
    "TD8": ("AIS", "GSM5702478"),
}


def load_section(prefix):
    """prefix like GSM5702473_TD1"""
    mtx = f"{DATA}/{prefix}_matrix.mtx.gz"
    feat = f"{DATA}/{prefix}_features.tsv.gz"
    bc = f"{DATA}/{prefix}_barcodes.tsv.gz"
    pos = f"{DATA}/{prefix}_tissue_positions_list.csv.gz"
    a = sc.read_mtx(mtx).T
    genes = pd.read_csv(feat, sep="\t", header=None)
    barcodes = pd.read_csv(bc, header=None)[0].astype(str).values
    a.var_names = genes[1].astype(str).values
    a.obs_names = barcodes
    a.var_names_make_unique()
    p = pd.read_csv(pos, header=None)
    p.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    p["barcode"] = p["barcode"].astype(str)
    p = p.set_index("barcode")
    shared = a.obs_names.intersection(p.index)
    a = a[shared].copy()
    a.obs["in_tissue"] = p.loc[shared, "in_tissue"].values
    a.obs["array_row"] = p.loc[shared, "array_row"].values
    a.obs["array_col"] = p.loc[shared, "array_col"].values
    a.obsm["spatial"] = p.loc[shared, ["pxl_col", "pxl_row"]].to_numpy()
    # keep in-tissue when the flag is present; some files mark all 0/1
    if (a.obs["in_tissue"] == 1).any():
        a = a[a.obs["in_tissue"] == 1].copy()
    return a


def main():
    rows = []
    prefixes = sorted({os.path.basename(p).rsplit("_", 1)[0]
                       for p in glob.glob(f"{DATA}/*_matrix.mtx.gz")})
    for pref in prefixes:
        sid = pref.split("_", 1)[1]  # TD1
        stage, gsm = STAGE[sid]
        a = load_section(pref)
        sc.pp.filter_cells(a, min_genes=200)
        sc.pp.filter_genes(a, min_cells=3)
        if a.n_obs < 50:
            print(f"skip {sid}: n={a.n_obs}")
            continue
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        extra = dict(dataset="GSE189487", platform="Visium", group="tumor",
                     patient=sid, disease=f"LUAD_{stage}", stage=stage, gsm=gsm)
        row = common.analyze_sample(a, sid, extra)
        rows.append(row)
        print(f"{sid} [{stage}] n={row['n_spots']} "
              f"rho(epi,TB)={row['rho_epi_vs_TB']:+.3f} "
              f"rho(epi,T)={row['rho_epi_vs_T']:+.3f} "
              f"rho(epi,B)={row['rho_epi_vs_B']:+.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/tables/gse189487_per_section.csv", index=False)
    print(f"wrote {OUT}/tables/gse189487_per_section.csv ({len(df)} sections)")


if __name__ == "__main__":
    main()

"""Hex-ring neighborhood analysis for E-MTAB-13530 and GSE189487.

Pre-specified (not tuned after seeing results):
  * Rings 1/2/3 on the Visium honeycomb (self excluded, ≥3 neighbors).
  * Index = per-spot CLDN4 (log-norm) and, separately, CLDN4+TACSTD2 epi score.
  * Outcome = mean T / B / TB signature in the ring.
  * Primary cross-section summary: median Spearman and two-sided sign test
    on tumor sections only.
  * Sensitivity: partial Spearman of CLDN4 vs ring-TB controlling for the
    index epithelial score (EPCAM/KRT panel). This is closer to PR #110
    but does NOT drop immune-rich index spots (that filter is PR #110 only).
"""
import os
import glob
import tarfile
import tempfile
import warnings
import numpy as np
import pandas as pd
import scanpy as sc

import common
from neighborhood import neighborhood_rows

warnings.simplefilter("ignore")
sc.settings.verbosity = 0

OUT = "results/claim_B6"
os.makedirs(f"{OUT}/tables", exist_ok=True)


def gene_vec(adata, gene):
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    x = adata[:, gene].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    return x.astype(float)


def attach_scores(a):
    common.score_signature(a, common.EPI_MARKERS, "epi_score")
    common.score_signature(a, common.EPI_PANEL_BROAD, "epi_broad_score")
    common.score_signature(a, common.T_MARKERS, "t_score")
    common.score_signature(a, common.B_MARKERS, "b_score")
    tb = common.present(common.T_MARKERS + common.B_MARKERS, a.var_names)
    sc.tl.score_genes(a, tb, score_name="tb_score", use_raw=False)
    a.obs["cldn4"] = gene_vec(a, "CLDN4")
    a.obs["tacstd2"] = gene_vec(a, "TACSTD2")
    return a


def run_section(a, sample, group, dataset):
    if "array_row" not in a.obs.columns:
        return []
    return neighborhood_rows(
        sample, group, dataset,
        a.obs["array_row"].values, a.obs["array_col"].values,
        a.obs["cldn4"].values, a.obs["epi_broad_score"].values,
        a.obs["t_score"].values, a.obs["b_score"].values, a.obs["tb_score"].values,
    )


def emtab_meta():
    df = pd.read_csv("data/emtab13530/E-MTAB-13530.sdrf.txt", sep="\t")
    m = {}
    for _, r in df.iterrows():
        s = r["Source Name"]
        site = str(r.get("Characteristics[sampling site]", ""))
        if site.strip() == "tumor":
            grp = "tumor"
        elif "adjacent" in site:
            grp = "adjacent_normal"
        else:
            grp = "healthy"
        m[s] = grp
    return m


def load_emtab(sid):
    h5 = f"data/emtab13530/{sid}-filtered_feature_bc_matrix.h5"
    star = f"data/emtab13530/{sid}-spatial.tar"
    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    with tarfile.open(star) as tf:
        names = tf.getnames()
        pos_name = [n for n in names if n.endswith("tissue_positions_list.csv")
                    or n.endswith("tissue_positions.csv")][0]
        with tempfile.TemporaryDirectory() as td:
            tf.extract(pos_name, path=td)
            p = pd.read_csv(os.path.join(td, pos_name), header=None)
    # 10x tissue_positions_list: barcode, in_tissue, row, col, pxl_row, pxl_col
    # newer: header present
    if str(p.iloc[0, 0]).lower() in ("barcode", "barcodes"):
        p.columns = p.iloc[0]
        p = p.iloc[1:].reset_index(drop=True)
        cols = {c.lower(): c for c in p.columns}
        p = p.rename(columns={
            cols.get("barcode", p.columns[0]): "barcode",
            [c for c in p.columns if "tissue" in str(c).lower()][0]: "in_tissue",
            [c for c in p.columns if "row" in str(c).lower() and "pxl" not in str(c).lower()][0]: "array_row",
            [c for c in p.columns if "col" in str(c).lower() and "pxl" not in str(c).lower()][0]: "array_col",
        })
    else:
        p.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    p["barcode"] = p["barcode"].astype(str)
    p = p.set_index("barcode")
    shared = a.obs_names.intersection(p.index)
    a = a[shared].copy()
    a.obs["in_tissue"] = pd.to_numeric(p.loc[shared, "in_tissue"], errors="coerce").values
    a.obs["array_row"] = pd.to_numeric(p.loc[shared, "array_row"], errors="coerce").values
    a.obs["array_col"] = pd.to_numeric(p.loc[shared, "array_col"], errors="coerce").values
    if (a.obs["in_tissue"] == 1).any():
        a = a[a.obs["in_tissue"] == 1].copy()
    sc.pp.filter_cells(a, min_genes=200)
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return attach_scores(a)


def load_gse189(prefix):
    from run_gse189487 import load_section
    a = load_section(prefix)
    sc.pp.filter_cells(a, min_genes=200)
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return attach_scores(a)


def main():
    rows = []
    meta = emtab_meta()
    for h5 in sorted(glob.glob("data/emtab13530/*-filtered_feature_bc_matrix.h5")):
        sid = os.path.basename(h5).replace("-filtered_feature_bc_matrix.h5", "")
        grp = meta.get(sid, "NA")
        try:
            a = load_emtab(sid)
        except Exception as e:
            print(f"EMTAB {sid} FAIL {e}")
            continue
        if a.n_obs < 80:
            print(f"skip {sid} n={a.n_obs}")
            continue
        recs = run_section(a, sid, grp, "E-MTAB-13530")
        rows.extend(recs)
        r1 = recs[0]
        print(f"{sid:10s} [{grp:15s}] ring1 CLDN4~TBnb ρ={r1['rho_CLDN4_vs_TBnb']:+.3f} "
              f"partial={r1['rho_CLDN4_vs_TBnb_partialEpi']:+.3f}")

    for mtx in sorted(glob.glob("data/gse189487/raw/*_matrix.mtx.gz")):
        pref = os.path.basename(mtx).rsplit("_", 1)[0]
        sid = pref.split("_", 1)[1]
        a = load_gse189(pref)
        if a.n_obs < 80:
            print(f"skip {sid} n={a.n_obs}")
            continue
        recs = run_section(a, sid, "tumor", "GSE189487")
        rows.extend(recs)
        r1 = recs[0]
        print(f"{sid:10s} [GSE189487 tumor] ring1 CLDN4~TBnb ρ={r1['rho_CLDN4_vs_TBnb']:+.3f} "
              f"partial={r1['rho_CLDN4_vs_TBnb_partialEpi']:+.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/tables/visium_neighborhood_rings.csv", index=False)
    print(f"wrote {OUT}/tables/visium_neighborhood_rings.csv ({len(df)} rows)")


if __name__ == "__main__":
    main()

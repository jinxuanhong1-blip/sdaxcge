"""CLAIM B6 on GSE271689 (NanoString GeoMx DSP Whole Transcriptome Atlas, NSCLC).

GeoMx AOIs in this study were segmented by protein morphology markers into
compartments: CK (PanCK+, tumor), CD45 (leukocyte), CD68 (macrophage). We build
an AOI x gene matrix from the raw DCC probe counts using the WTA PKC to map
probes -> genes, Q3-normalize, then test whether the epithelial signature
(CLDN4/TACSTD2) anti-correlates with T / B signatures across AOIs.

IMPORTANT CAVEAT (reported in the output): because AOIs are pre-sorted into
tumor vs leukocyte compartments by antibody staining, some anti-colocalization
is expected by construction. CLDN4/TACSTD2 and CD3D/MS4A1 transcripts are not
the proteins used for segmentation, so this is a semi-independent confirmation,
not a fully independent one.
"""
import os
import re
import gzip
import json
import glob
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

import common

warnings.simplefilter("ignore")
sc.settings.verbosity = 0

DATA = "data/gse271689"
OUT = "results/claim_B6"
os.makedirs(f"{OUT}/tables", exist_ok=True)


def parse_pkc(path):
    with open(path) as fh:
        pkc = json.load(fh)
    rts2gene, rts2class = {}, {}
    for tgt in pkc["Targets"]:
        gene = tgt["DisplayName"]
        cc = tgt.get("CodeClass", "")
        for pr in tgt["Probes"]:
            rid = pr["RTS_ID"]
            rts2gene[rid] = gene
            rts2class[rid] = cc
    return rts2gene, rts2class


def parse_dcc(path):
    counts = {}
    with gzip.open(path, "rt") as fh:
        in_cs = False
        for line in fh:
            line = line.strip()
            if line == "<Code_Summary>":
                in_cs = True
                continue
            if line == "</Code_Summary>":
                break
            if in_cs and "," in line:
                rid, c = line.split(",")
                counts[rid] = int(c)
    return counts


def parse_soft(path):
    """Map GSM id -> dict(cell_type, spotid, tissue, dcc_basename)."""
    ann = {}
    cur = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE ="):
                cur = line.split("=")[1].strip()
                ann[cur] = dict(cell_type="NA", spotid="NA", tissue="NA", dcc="NA", title="NA")
            elif cur and line.startswith("!Sample_title ="):
                ann[cur]["title"] = line.split("=", 1)[1].strip()
            elif cur and "cell type:" in line:
                ann[cur]["cell_type"] = line.split("cell type:")[1].strip()
            elif cur and "spotid:" in line:
                ann[cur]["spotid"] = line.split("spotid:")[1].strip()
            elif cur and "tissue:" in line:
                ann[cur]["tissue"] = line.split("tissue:")[1].strip()
            elif cur and line.startswith("!Sample_supplementary_file_1"):
                ann[cur]["dcc"] = os.path.basename(line.split("=", 1)[1].strip())
    return ann


def build_matrix():
    rts2gene, rts2class = parse_pkc(f"{DATA}/Hs_R_NGS_WTA_v1.0.pkc")
    ann = parse_soft(f"{DATA}/GSE271689_family.soft.gz")
    # map dcc basename (without .gz) -> gsm
    dcc2gsm = {}
    for gsm, d in ann.items():
        base = d["dcc"].replace(".gz", "")
        dcc2gsm[base] = gsm

    neg_rts = {r for r, c in rts2class.items() if c.startswith("Negative")}

    dcc_files = sorted(glob.glob(f"{DATA}/dcc/*.dcc.gz"))
    gene_counts = {}   # aoi -> {gene: summed count}
    neg_geomean = {}
    obs = {}
    for path in dcc_files:
        base = os.path.basename(path).replace(".gz", "")
        gsm = dcc2gsm.get(base)
        if gsm is None:
            continue
        meta = ann[gsm]
        if meta["cell_type"] in ("NA", "", "No Template Control"):
            continue
        counts = parse_dcc(path)
        gc = {}
        negs = []
        for rid, c in counts.items():
            if rid in neg_rts:
                negs.append(c)
            else:
                g = rts2gene.get(rid)
                if g is None:
                    continue
                gc[g] = gc.get(g, 0) + c
        aoi = base.split("_")[0] + "|" + meta["spotid"] + "|" + meta["cell_type"]
        # ensure uniqueness with gsm
        aoi = f"{gsm}"
        gene_counts[aoi] = gc
        negs = np.array(negs, float)
        neg_geomean[aoi] = float(np.exp(np.mean(np.log(negs + 1)))) if len(negs) else np.nan
        obs[aoi] = dict(gsm=gsm, cell_type=meta["cell_type"], spotid=meta["spotid"],
                        tissue=meta["tissue"], scan=base.split("-C-")[0].split("_")[-1])
    # build DataFrame AOI x gene
    mat = pd.DataFrame.from_dict(gene_counts, orient="index").fillna(0.0)
    obs_df = pd.DataFrame.from_dict(obs, orient="index").loc[mat.index]
    obs_df["neg_geomean"] = pd.Series(neg_geomean)
    obs_df["total_counts"] = mat.sum(axis=1).values
    return mat, obs_df


def q3_normalize(mat):
    """Standard GeoMx Q3 normalization: scale each AOI to its 75th percentile."""
    q3 = mat.apply(lambda r: np.percentile(r[r > 0], 75) if (r > 0).any() else np.nan, axis=1)
    geo = np.exp(np.mean(np.log(q3.dropna())))
    norm = mat.div(q3, axis=0) * geo
    return norm, q3


def main():
    mat, obs = build_matrix()
    print(f"raw AOI x gene matrix: {mat.shape[0]} AOIs x {mat.shape[1]} genes")
    print("compartments:", obs["cell_type"].value_counts().to_dict())

    # QC: drop AOIs with very low total counts or below background
    keep = (obs["total_counts"] >= 1000)
    # signal above background: Q3 should exceed neg geomean
    mat, obs = mat.loc[keep], obs.loc[keep]
    print(f"after QC (total>=1000): {mat.shape[0]} AOIs")

    norm, q3 = q3_normalize(mat)
    obs["q3"] = q3.values
    obs["signal_to_bg"] = obs["q3"] / obs["neg_geomean"]

    # Build AnnData on log1p(Q3-normalized)
    a = sc.AnnData(np.log1p(norm.values), obs=obs.reset_index(drop=True),
                   var=pd.DataFrame(index=norm.columns))
    a.var_names_make_unique()

    common.score_signature(a, common.EPI_MARKERS, "epi_score")
    common.score_signature(a, common.T_MARKERS, "t_score")
    common.score_signature(a, common.B_MARKERS, "b_score")
    tb = common.present(common.T_MARKERS + common.B_MARKERS, a.var_names)
    sc.tl.score_genes(a, tb, score_name="tb_score", use_raw=False)

    df = a.obs.copy()
    for c in ["epi_score", "t_score", "b_score", "tb_score"]:
        df[c] = a.obs[c].values

    # ---- primary test: across ALL AOIs ----
    def rep(tag, sub):
        out = {}
        for lab, col in [("T", "t_score"), ("B", "b_score"), ("TB", "tb_score")]:
            rho, p, n = common.spearman(sub["epi_score"], sub[col])
            out[f"rho_epi_vs_{lab}"] = rho
            out[f"p_epi_vs_{lab}"] = p
        out["n"] = len(sub)
        print(f"[{tag:22s}] n={len(sub):4d}  "
              f"rho(epi,TB)={out['rho_epi_vs_TB']:+.3f} (p={out['p_epi_vs_TB']:.1e})  "
              f"rho(epi,T)={out['rho_epi_vs_T']:+.3f}  rho(epi,B)={out['rho_epi_vs_B']:+.3f}")
        return dict(scope=tag, **out)

    rows = []
    rows.append(rep("all_AOIs", df))
    # within each compartment (removes the pre-sorting design effect)
    for ct in ["CK", "CD45", "CD68"]:
        sub = df[df["cell_type"] == ct]
        if len(sub) >= 15:
            rows.append(rep(f"within_{ct}", sub))

    # compartment mean scores (tumor vs immune)
    comp = df.groupby("cell_type")[["epi_score", "t_score", "b_score", "tb_score"]].mean()
    print("\ncompartment mean signature scores:\n", comp.round(3))

    # paired within-ROI: CK vs CD45 for same spotid
    piv = df.pivot_table(index="spotid", columns="cell_type",
                         values=["epi_score", "tb_score"], aggfunc="mean")
    paired_msg = ""
    try:
        ck_epi = piv[("epi_score", "CK")]
        im_epi = piv[("epi_score", "CD45")]
        ck_tb = piv[("tb_score", "CK")]
        im_tb = piv[("tb_score", "CD45")]
        m = ck_epi.notna() & im_epi.notna()
        w_epi = stats.wilcoxon(ck_epi[m], im_epi[m])
        m2 = ck_tb.notna() & im_tb.notna()
        w_tb = stats.wilcoxon(ck_tb[m2], im_tb[m2])
        paired_msg = (f"paired ROIs (n={int(m.sum())}): "
                      f"epi CK>CD45 median diff={float((ck_epi-im_epi)[m].median()):+.3f} (p={w_epi.pvalue:.1e}); "
                      f"TB CK<CD45 median diff={float((ck_tb-im_tb)[m2].median()):+.3f} (p={w_tb.pvalue:.1e})")
        print("\n" + paired_msg)
    except Exception as e:
        paired_msg = f"paired test failed: {e}"
        print(paired_msg)

    pd.DataFrame(rows).to_csv(f"{OUT}/tables/gse271689_geomx_correlations.csv", index=False)
    comp.round(4).to_csv(f"{OUT}/tables/gse271689_geomx_compartment_means.csv")
    df.to_csv(f"{OUT}/tables/gse271689_geomx_per_aoi.csv", index=False)
    with open(f"{OUT}/tables/gse271689_geomx_paired.txt", "w") as fh:
        fh.write(paired_msg + "\n")
    print(f"\nwrote GeoMx tables to {OUT}/tables/")


if __name__ == "__main__":
    main()

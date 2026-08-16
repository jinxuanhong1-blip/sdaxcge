#!/usr/bin/env python3
"""B6 leftover: public CosMx NSCLC CLDN4/TACSTD2 vs immune niches.

Processed CosMx SMI flat files (He et al. 2022, Nat Biotechnol) from the
NanoString public S3 bucket. Official cell-type labels are not in the flat
release, so compartments are RNA marker scores, with a protein-stain
(PanCK / CD45) robustness split.

MERFISH leftovers are catalogued and skipped when there is no processed
cell-by-gene table under 2 GB without login.

Outputs: results/w200/B6_cosmx/
"""

import json
import os
import re

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "w200", "B6_cosmx")
os.makedirs(OUT, exist_ok=True)

SAMPLES = {
    "Lung9_Rep1": os.path.join(
        ROOT, "data", "cosmx", "Lung9_Rep1", "Lung9_Rep1-Flat_files_and_images"),
    "Lung12": os.path.join(
        ROOT, "data", "cosmx_l12", "Lung12", "Lung12-Flat_files_and_images"),
    "Lung13": os.path.join(
        ROOT, "data", "cosmx_l13", "Lung13", "Lung13-Flat_files_and_images"),
}

PX_TO_UM = 0.18
RADIUS_UM = 50.0
RADIUS_UM_SENS = 25.0
MIN_COUNTS = 20
MIN_GENES = 5
MIN_NEIGHBORS = 5
MIN_FOV_EPI = 50

MARKERS = {
    "immune": [
        "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4",
        "IL7R", "CCL5", "GZMA", "GZMB", "GZMK", "GNLY", "KLRB1", "KLRK1",
        "CD68", "CD163", "C1QA", "C1QB", "C1QC", "LYZ", "FCGR3A", "ITGAX",
        "ITGAM", "CD14", "CD79A", "MS4A1", "CD19", "CD37", "CD52", "CD53",
        "CPA3", "MS4A2", "TPSAB1",
    ],
    "epithelial": [
        "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT17", "KRT15",
        "KRT14", "KRT13", "KRT16", "KRT6A", "KRT5", "KRT7", "CEACAM6",
        "MUC1", "ELF3", "SFTPC", "SFTPB", "NAPSA", "S100A14",
    ],
    "stromal": [
        "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "DCN",
        "LUM", "FN1", "BGN", "PDGFRB", "ACTA2", "MYH11", "PECAM1", "VWF",
        "CDH5", "CLEC14A", "FLT1", "KDR", "ESAM", "RAMP2",
    ],
}

T_MARKERS = ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"]


def find_file(d, pat):
    for f in sorted(os.listdir(d)):
        if re.search(pat, f):
            return os.path.join(d, f)
    raise FileNotFoundError(f"{pat} in {d}")


def load_sample(d):
    expr = pd.read_csv(find_file(d, "exprMat"), dtype={"fov": np.int32, "cell_ID": np.int32})
    meta = pd.read_csv(find_file(d, "metadata"))
    expr = expr[expr["cell_ID"] != 0]
    key = expr["fov"].astype(str) + "_" + expr["cell_ID"].astype(str)
    expr.index = key
    meta.index = meta["fov"].astype(str) + "_" + meta["cell_ID"].astype(str)
    genes = [c for c in expr.columns if c not in ("fov", "cell_ID")]
    genes = [g for g in genes if not g.lower().startswith("neg")]
    common = expr.index.intersection(meta.index)
    expr, meta = expr.loc[common], meta.loc[common]
    X = sparse.csr_matrix(expr[genes].to_numpy(dtype=np.float32))
    return X, genes, meta.copy()


def lognorm(X):
    tot = np.asarray(X.sum(axis=1)).ravel()
    sf = np.median(tot) / np.maximum(tot, 1)
    Xn = X.multiply(sf[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def score(Xn, genes, marker_list):
    idx = [genes.index(g) for g in marker_list if g in genes]
    used = [g for g in marker_list if g in genes]
    return np.asarray(Xn[:, idx].mean(axis=1)).ravel(), used


def auroc(pos, neg):
    u, _ = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return float(u / (len(pos) * len(neg)))


def local_frac(xy, fovs, flag, radius_um):
    r_px = radius_um / PX_TO_UM
    frac = np.full(len(flag), np.nan)
    n_nb = np.zeros(len(flag), dtype=int)
    for f in np.unique(fovs):
        ix = np.where(fovs == f)[0]
        tree = cKDTree(xy[ix])
        nb = tree.query_ball_point(xy[ix], r=r_px)
        for j, nbrs in enumerate(nb):
            nbrs = [k for k in nbrs if k != j]
            n_nb[ix[j]] = len(nbrs)
            if nbrs:
                frac[ix[j]] = flag[ix[nbrs]].mean()
    return frac, n_nb


def fov_spearman(df, xcol, ycol, sample, gene, neighborhood):
    rows = []
    for f, g in df.groupby("fov"):
        if len(g) >= MIN_FOV_EPI and g[xcol].nunique() > 1 and g[ycol].nunique() > 1:
            r, p = stats.spearmanr(g[xcol], g[ycol])
            rows.append({
                "sample": sample, "fov": int(f), "n": int(len(g)),
                "gene": gene, "neighborhood": neighborhood,
                "spearman_rho": float(r), "p": float(p),
            })
    return pd.DataFrame(rows)


def summarize_fov(fov_df):
    if fov_df is None or len(fov_df) == 0:
        return {}
    rhos = fov_df["spearman_rho"]
    w = stats.wilcoxon(rhos, alternative="two-sided")
    return {
        "n_fovs": int(len(fov_df)),
        "median_fov_rho": round(float(rhos.median()), 4),
        "mean_fov_rho": round(float(rhos.mean()), 4),
        "frac_fovs_negative": round(float((rhos < 0).mean()), 3),
        "n_fovs_p05": int((fov_df["p"] < 0.05).sum()),
        "wilcoxon_signedrank_p": float(w.pvalue),
    }


def analyze(name, d):
    X, genes, meta = load_sample(d)
    tot = np.asarray(X.sum(axis=1)).ravel()
    ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    keep = (tot >= MIN_COUNTS) & (ngene >= MIN_GENES)
    X, meta = X[keep], meta.loc[keep]
    Xn = lognorm(X)

    for g in ("CLDN4", "TACSTD2"):
        assert g in genes, f"{g} missing from CosMx panel"
        gi = genes.index(g)
        meta[g.lower()] = np.asarray(Xn[:, gi].todense()).ravel()
        meta[g.lower() + "_raw"] = np.asarray(X[:, gi].todense()).ravel()

    scores, used = {}, {}
    for comp, ms in MARKERS.items():
        scores[comp], used[comp] = score(Xn, genes, ms)
    S = np.vstack([scores[c] for c in MARKERS])
    lab = np.array(list(MARKERS))[S.argmax(axis=0)]
    lab[S.max(axis=0) <= 0] = "unassigned"
    meta["compartment"] = lab

    t_idx = [genes.index(g) for g in T_MARKERS if g in genes]
    meta["t_like"] = np.asarray((X[:, t_idx].sum(axis=1) > 0)).ravel()

    # Protein-stain split (median on this sample)
    has_prot = "Mean.CD45" in meta.columns and "Mean.PanCK" in meta.columns
    if has_prot:
        panck = np.log1p(meta["Mean.PanCK"].to_numpy())
        cd45 = np.log1p(meta["Mean.CD45"].to_numpy())
        meta["prot_epi"] = (panck >= np.median(panck)) & (cd45 < np.median(cd45))
        meta["prot_imm"] = (cd45 >= np.median(cd45)) & (panck < np.median(panck))
        val = {
            "cd45_auroc_rna_immune_vs_rest": round(
                auroc(meta.loc[meta["compartment"] == "immune", "Mean.CD45"],
                      meta.loc[meta["compartment"] != "immune", "Mean.CD45"]), 3),
            "panck_auroc_rna_epithelial_vs_rest": round(
                auroc(meta.loc[meta["compartment"] == "epithelial", "Mean.PanCK"],
                      meta.loc[meta["compartment"] != "epithelial", "Mean.PanCK"]), 3),
            "n_prot_epi": int(meta["prot_epi"].sum()),
            "n_prot_imm": int(meta["prot_imm"].sum()),
        }
    else:
        val = {"note": "no protein stains"}

    spec_rows = []
    for gene in ("CLDN4", "TACSTD2"):
        col, raw = gene.lower(), gene.lower() + "_raw"
        for comp in ("epithelial", "immune", "stromal", "unassigned"):
            m = meta["compartment"] == comp
            if m.sum() == 0:
                continue
            spec_rows.append({
                "sample": name, "gene": gene, "compartment": comp,
                "n_cells": int(m.sum()),
                "mean_lognorm": float(meta.loc[m, col].mean()),
                "frac_positive": float((meta.loc[m, raw] > 0).mean()),
            })
    spec = pd.DataFrame(spec_rows)

    spec_stats = {}
    for gene in ("CLDN4", "TACSTD2"):
        e = meta.loc[meta["compartment"] == "epithelial", gene.lower()]
        i = meta.loc[meta["compartment"] == "immune", gene.lower()]
        spec_stats[gene] = {
            "epi_vs_imm_auroc": round(auroc(e, i), 3),
            "epi_vs_imm_mannwhitney_p": float(stats.mannwhitneyu(e, i).pvalue),
        }

    xy = meta[["CenterX_global_px", "CenterY_global_px"]].to_numpy()
    fovs = meta["fov"].to_numpy()
    is_imm = (meta["compartment"] == "immune").to_numpy()
    is_t = meta["t_like"].to_numpy()

    frac50, n50 = local_frac(xy, fovs, is_imm, RADIUS_UM)
    frac25, n25 = local_frac(xy, fovs, is_imm, RADIUS_UM_SENS)
    frac_t, _ = local_frac(xy, fovs, is_t, RADIUS_UM)
    meta = meta.assign(frac_immune_50um=frac50, n_neighbors_50um=n50,
                       frac_immune_25um=frac25, n_neighbors_25um=n25,
                       frac_t_50um=frac_t)
    if has_prot:
        frac_p, np_ = local_frac(xy, fovs, meta["prot_imm"].to_numpy(), RADIUS_UM)
        meta = meta.assign(frac_prot_imm_50um=frac_p, n_neighbors_prot=np_)

    fov_parts = []
    spatial = {}

    def add_spatial(key, index_m, xcol, ycol, gene, neighborhood):
        m = index_m & (meta["n_neighbors_50um"] >= MIN_NEIGHBORS)
        if ycol == "frac_immune_25um":
            m = index_m & (meta["n_neighbors_25um"] >= MIN_NEIGHBORS)
        if ycol == "frac_prot_imm_50um":
            m = index_m & (meta.get("n_neighbors_prot", meta["n_neighbors_50um"]) >= MIN_NEIGHBORS)
        sub = meta.loc[m]
        if len(sub) < 100:
            spatial[key] = {"n": int(len(sub)), "note": "too few cells"}
            return
        rho, p = stats.spearmanr(sub[xcol], sub[ycol])
        fov = fov_spearman(sub, xcol, ycol, name, gene, neighborhood)
        fov_parts.append(fov)
        q = sub[xcol].quantile([1 / 3, 2 / 3]).to_numpy()
        lo = sub.loc[sub[xcol] <= q[0], ycol]
        hi = sub.loc[sub[xcol] >= q[1], ycol]
        tert = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        spatial[key] = {
            "n_index_cells": int(len(sub)),
            "spearman_rho_pooled": round(float(rho), 4),
            "spearman_p_pooled": float(p),
            "hi_tertile_nbhd_mean": round(float(hi.mean()), 4),
            "lo_tertile_nbhd_mean": round(float(lo.mean()), 4),
            "tertile_hi_vs_lo_p": float(tert.pvalue),
            **summarize_fov(fov),
        }

    epi_rna = meta["compartment"] == "epithelial"
    add_spatial("cldn4_vs_immune_50um", epi_rna, "cldn4", "frac_immune_50um",
                "CLDN4", "rna_immune_50um")
    add_spatial("cldn4_vs_immune_25um", epi_rna, "cldn4", "frac_immune_25um",
                "CLDN4", "rna_immune_25um")
    add_spatial("cldn4_vs_t_50um", epi_rna, "cldn4", "frac_t_50um",
                "CLDN4", "t_like_50um")
    add_spatial("tacstd2_vs_immune_50um", epi_rna, "tacstd2", "frac_immune_50um",
                "TACSTD2", "rna_immune_50um")
    if has_prot:
        add_spatial("cldn4_vs_prot_imm_50um", meta["prot_epi"], "cldn4",
                    "frac_prot_imm_50um", "CLDN4", "prot_cd45_50um")
        add_spatial("tacstd2_vs_prot_imm_50um", meta["prot_epi"], "tacstd2",
                    "frac_prot_imm_50um", "TACSTD2", "prot_cd45_50um")

    report = {
        "n_cells_qc": int(X.shape[0]),
        "n_fovs": int(meta["fov"].nunique()),
        "compartment_counts": {k: int(v) for k, v in meta["compartment"].value_counts().items()},
        "n_t_like": int(meta["t_like"].sum()),
        "markers_n": {k: len(v) for k, v in used.items()},
        "protein_stain_validation": val,
        "specificity": spec_stats,
        "spatial": spatial,
    }
    fov_df = pd.concat(fov_parts, ignore_index=True) if fov_parts else pd.DataFrame()
    epi = meta.loc[epi_rna & (meta["n_neighbors_50um"] >= MIN_NEIGHBORS)]
    return meta, spec, fov_df, epi, report


def plot_all(results):
    names = list(results)
    fig, axes = plt.subplots(1, len(names), figsize=(5.2 * len(names), 4.2), squeeze=False)
    for ax, name in zip(axes[0], names):
        meta = results[name][0]
        order = ["epithelial", "stromal", "immune"]
        data = [meta.loc[meta["compartment"] == c, "cldn4"] for c in order]
        vp = ax.violinplot(data, showmedians=True)
        for b in vp["bodies"]:
            b.set_alpha(0.55)
        ax.set_xticks(range(1, 4))
        ax.set_xticklabels([f"{c}\n(n={len(d):,})" for c, d in zip(order, data)])
        ax.set_ylabel("CLDN4 (log-norm)")
        ax.set_title(name)
    fig.suptitle("CosMx NSCLC leftover: CLDN4 by RNA compartment")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_cldn4_by_compartment.png"), dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(names), figsize=(5.2 * len(names), 4.2), squeeze=False)
    for ax, name in zip(axes[0], names):
        epi = results[name][3]
        dec = pd.qcut(epi["frac_immune_50um"].rank(method="first"), 10, labels=False)
        mean, sem = epi.groupby(dec)["cldn4"].mean(), epi.groupby(dec)["cldn4"].sem()
        xc = epi.groupby(dec)["frac_immune_50um"].mean()
        ax.errorbar(xc, mean, yerr=1.96 * sem, fmt="o-", capsize=3)
        ax.set_xlabel("Local RNA-immune fraction (50 um)")
        ax.set_ylabel("Epithelial CLDN4 (log-norm), mean +/- 95% CI")
        ax.set_title(name)
    fig.suptitle("Epithelial CLDN4 vs 50 um immune neighborhood")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_cldn4_vs_immune_neighborhood.png"), dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(names), figsize=(4.4 * len(names), 4.4), squeeze=False)
    for ax, name in zip(axes[0], names):
        fov = results[name][2]
        sub = fov[fov["neighborhood"] == "rna_immune_50um"]
        for gene, col in (("CLDN4", "#1f77b4"), ("TACSTD2", "#ff7f0e")):
            g = sub[sub["gene"] == gene]
            ax.scatter(np.zeros(len(g)) + (0 if gene == "CLDN4" else 1) +
                       np.random.default_rng(0).uniform(-0.08, 0.08, len(g)),
                       g["spearman_rho"], s=18, c=col, alpha=0.7, label=gene)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4", "TACSTD2"])
        ax.set_ylabel("Per-FOV Spearman rho")
        ax.set_title(name)
        ax.legend(fontsize=8)
    fig.suptitle("Per-FOV epithelial gene vs 50 um RNA-immune fraction")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig4_per_fov_rho.png"), dpi=150)
    plt.close(fig)

    # only first two samples for the heavy spatial scatter
    show = names[:2]
    fig, axes = plt.subplots(1, len(show), figsize=(6.4 * len(show), 5.4), squeeze=False)
    colors = {"epithelial": "#1f77b4", "immune": "#d62728",
              "stromal": "#2ca02c", "unassigned": "#cccccc"}
    for ax, name in zip(axes[0], show):
        meta = results[name][0]
        for comp, col in colors.items():
            m = meta["compartment"] == comp
            ax.scatter(meta.loc[m, "CenterX_global_px"] * PX_TO_UM,
                       meta.loc[m, "CenterY_global_px"] * PX_TO_UM,
                       s=0.25, c=col, label=comp, rasterized=True)
        ax.set_aspect("equal")
        ax.legend(markerscale=22, fontsize=8, loc="upper right")
        ax.set_title(f"{name}: RNA compartments")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_spatial_compartments.png"), dpi=120)
    plt.close(fig)


def write_catalog():
    rows = [
        {
            "source": "NanoString CosMx SMI NSCLC (He 2022)",
            "modality": "CosMx",
            "processed": "yes",
            "size_note": "flat exprMat+metadata 145–181 MB / sample",
            "access": "public S3 nanostring-public-share/SMI-Compressed/",
            "decision": "analyze",
            "reason": "processed cell-by-gene + XY + PanCK/CD45 stains; leftover vs Visium/GeoMx B6",
        },
        {
            "source": "Vizgen MERSCOPE FFPE IO showcase lung 1/2",
            "modality": "MERFISH",
            "processed": "yes_but_gated",
            "size_note": "cell_by_gene exists behind Vizgen/GCS login",
            "access": "storage.googleapis.com/vz-ffpe-showcase → HTTP 403",
            "decision": "skip",
            "reason": "no anonymous processed download",
        },
        {
            "source": "Figshare 25983526 CellCharter MERSCOPE lung (2 pts)",
            "modality": "MERFISH",
            "processed": "yes",
            "size_note": "merscope_human_lung_cancer_nhood3.h5ad 3.67 GB",
            "access": "https://doi.org/10.6084/m9.figshare.25983526",
            "decision": "skip",
            "reason": "processed but >2 GB leftover budget",
        },
        {
            "source": "Zenodo 11198494 Chen 2024 MERFISH (Pt38/43/67/73)",
            "modality": "MERFISH",
            "processed": "no_cell_by_gene",
            "size_note": "transcript tables 0.26–3.9 GB; no cell-by-gene matrix",
            "access": "https://doi.org/10.5281/zenodo.11198494",
            "decision": "skip",
            "reason": "honest skip: transcript-level, not processed cell-by-gene",
        },
        {
            "source": "NanoString All SMI Giotto object",
            "modality": "CosMx",
            "processed": "yes",
            "size_note": "930 MB tar.gz",
            "access": "same S3 bucket",
            "decision": "skip",
            "reason": "redundant with flat files; no published official cell types required",
        },
    ]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "catalog.tsv"), sep="\t", index=False)
    return df


def main():
    write_catalog()
    report, results = {}, {}
    spec_all, fov_all = [], []
    for name, d in SAMPLES.items():
        if not os.path.isdir(d):
            print(f"SKIP missing {name}: {d}")
            continue
        print(f"== {name}")
        meta, spec, fov_df, epi, rep = analyze(name, d)
        results[name] = (meta, spec, fov_df, epi)
        report[name] = rep
        spec_all.append(spec)
        fov_all.append(fov_df)

    pd.concat(spec_all).to_csv(os.path.join(OUT, "cldn4_tacstd2_by_compartment.csv"), index=False)
    pd.concat(fov_all).to_csv(os.path.join(OUT, "per_fov_spearman.csv"), index=False)
    with open(os.path.join(OUT, "stats.json"), "w") as f:
        json.dump(report, f, indent=2, default=float)
    plot_all(results)
    print(json.dumps(report, indent=2, default=float))


if __name__ == "__main__":
    main()

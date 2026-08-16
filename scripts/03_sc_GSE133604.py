#!/usr/bin/env python3
"""GSE133604 (mouse KP lung tumors, 10x scRNA-seq, 4 samples): Tacstd2-high vs
-low epithelial/tumor cells -> preranked GSEA on keratin / tight-junction sets.

The four samples share one barcodes.tsv / genes.tsv and per-sample matrices
(GSM..._matrix.mtx.gz). Cells are pooled, QC-filtered, normalized (CP10k+log1p),
epithelial cells identified by an epithelial marker score (Epcam/Krt8/Krt18/
Sftpc/...), and epithelial cells split into Tacstd2-high vs -low (Tacstd2>0 vs
==0 among epithelial cells). Ranking metric = Wilcoxon score (scanpy
rank_genes_groups) on all genes; preranked GSEA via gseapy.
"""
from pathlib import Path

import gseapy as gp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp
from scipy import stats

DATA = Path("/workspace/data/GSE133604")
GMT = Path("/workspace/results/w200/A8_mouse/genesets/mouse_keratin_tj.gmt")
OUT = Path("/workspace/results/w200/A8_mouse/GSE133604_sc")
OUT.mkdir(parents=True, exist_ok=True)
SEED = 42
sc.settings.verbosity = 1

SAMPLES = {
    "GSM3912860_1_Ctrl_matrix.mtx.gz": "Ctrl",
    "GSM3912861_2_CtrlplusPD1_matrix.mtx.gz": "Ctrl_PD1",
    "GSM3912862_3_ko_matrix.mtx.gz": "Asf1aKO",
    "GSM3912863_4_KOplusPD1_matrix.mtx.gz": "Asf1aKO_PD1",
}
EPI_MARKERS = ["Epcam", "Krt8", "Krt18", "Krt19", "Cdh1", "Sftpc", "Sftpb",
               "Scgb1a1", "Sftpa1", "Nkx2-1"]


def load():
    genes = pd.read_csv(DATA / "GSE133604_genes.tsv.gz", sep="\t", header=None)
    barcodes = pd.read_csv(DATA / "GSE133604_barcodes.tsv.gz", sep="\t", header=None)[0].tolist()
    symbols = genes[1].astype(str).to_numpy()
    adatas = []
    for fname, label in SAMPLES.items():
        m = sio.mmread(DATA / fname).tocsc()  # genes x cells
        a = sc.AnnData(sp.csr_matrix(m.T.astype(np.float32)))
        a.var_names = symbols
        a.var_names_make_unique()
        a.obs_names = [f"{label}_{b}" for b in barcodes[: a.n_obs]]
        a.obs["sample"] = label
        adatas.append(a)
    adata = sc.concat(adatas, join="outer")
    adata.obs_names_make_unique()
    return adata


def main():
    adata = load()
    print(f"loaded {adata.n_obs} cells x {adata.n_vars} genes")

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[(adata.obs.n_genes_by_counts < 6000) & (adata.obs.pct_counts_mt < 20)].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"after QC: {adata.n_obs} cells x {adata.n_vars} genes")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    markers = [g for g in EPI_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, markers, score_name="epi_score", random_state=SEED)
    # epithelial = positive epithelial score AND detectable Epcam
    epcam = np.asarray(adata[:, "Epcam"].X.todense()).ravel() if "Epcam" in adata.var_names else np.zeros(adata.n_obs)
    adata.obs["is_epi"] = (adata.obs["epi_score"] > 0.1) & (epcam > 0)
    n_epi = int(adata.obs["is_epi"].sum())
    print(f"epithelial cells: {n_epi}")

    epi = adata[adata.obs["is_epi"]].copy()
    tac = np.asarray(epi[:, "Tacstd2"].X.todense()).ravel() if "Tacstd2" in epi.var_names else np.zeros(epi.n_obs)
    cld = np.asarray(epi[:, "Cldn4"].X.todense()).ravel() if "Cldn4" in epi.var_names else np.zeros(epi.n_obs)
    epi.obs["Tacstd2_expr"] = tac
    epi.obs["Cldn4_expr"] = cld
    epi.obs["tac_group"] = np.where(tac > 0, "Tacstd2_high", "Tacstd2_low")
    epi.obs["cld_group"] = np.where(cld > 0, "Cldn4_high", "Cldn4_low")
    print("epi Tacstd2 groups:\n", epi.obs["tac_group"].value_counts())
    print("epi Cldn4 groups:\n", epi.obs["cld_group"].value_counts())
    both = int(((tac > 0) & (cld > 0)).sum())
    print(f"epi Tacstd2+ and Cldn4+: {both} / {n_epi}")

    summary = pd.DataFrame({
        "metric": [
            "total_cells_after_QC", "epithelial_cells",
            "epi_Tacstd2_high", "epi_Tacstd2_low", "pct_epi_Tacstd2_pos",
            "epi_Cldn4_high", "epi_Cldn4_low", "pct_epi_Cldn4_pos",
            "epi_double_pos_Tacstd2_Cldn4",
        ],
        "value": [
            adata.n_obs, n_epi,
            int((tac > 0).sum()), int((tac == 0).sum()), round(100 * (tac > 0).mean(), 2),
            int((cld > 0).sum()), int((cld == 0).sum()), round(100 * (cld > 0).mean(), 2),
            both,
        ],
    })
    summary.to_csv(OUT / "cell_summary.tsv", sep="\t", index=False)
    epi.obs.groupby(["sample", "tac_group"]).size().unstack(fill_value=0).to_csv(
        OUT / "epi_Tacstd2_by_sample.tsv", sep="\t")
    epi.obs.groupby(["sample", "cld_group"]).size().unstack(fill_value=0).to_csv(
        OUT / "epi_Cldn4_by_sample.tsv", sep="\t")

    # cell-level co-expression is descriptive only (pseudoreplication)
    r_tc, p_tc = stats.spearmanr(tac, cld)
    pd.DataFrame([{"unit": "epithelial_cell", "n": n_epi,
                   "spearman_r": r_tc, "spearman_p": p_tc,
                   "note": "cell-level p is pseudoreplication; inferential unit is the 4 libraries"}]
                 ).to_csv(OUT / "tacstd2_vs_cldn4_epi.tsv", sep="\t", index=False)

    all_gsea = []
    for label, col, pos, neg, drop in (
        ("Tacstd2_high_vs_low_epi", "tac_group", "Tacstd2_high", "Tacstd2_low", "Tacstd2"),
        ("Cldn4_high_vs_low_epi", "cld_group", "Cldn4_high", "Cldn4_low", "Cldn4"),
    ):
        counts = epi.obs[col].value_counts()
        if counts.min() < 20:
            print(f"skip {label}: too few cells ({counts.to_dict()})")
            continue
        sc.tl.rank_genes_groups(epi, col, groups=[pos], reference=neg, method="wilcoxon")
        de = sc.get.rank_genes_groups_df(epi, group=pos)
        de.to_csv(OUT / f"de_{label}.tsv", sep="\t", index=False)
        rnk = de.set_index("names")["scores"].dropna().drop(index=drop, errors="ignore")
        rnk.sort_values(ascending=False).to_csv(OUT / f"rank_{label}.rnk", sep="\t", header=False)
        res = gp.prerank(
            rnk=rnk.sort_values(ascending=False).reset_index(),
            gene_sets=str(GMT), min_size=5, max_size=1000,
            permutation_num=10000, seed=SEED, threads=4, outdir=None, no_plot=True,
        )
        gsea = res.res2d.copy()
        gsea.insert(0, "contrast", label)
        all_gsea.append(gsea)
        g = gsea.sort_values("NES")
        colors = ["#c0392b" if float(f) < 0.05 else "#f1948a" if float(f) < 0.25 else "#95a5a6"
                  for f in g["FDR q-val"]]
        fig, ax = plt.subplots(figsize=(7, 0.42 * len(g) + 1.4))
        ax.barh(g["Term"], g["NES"].astype(float), color=colors)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel(f"NES (positive = enriched in {pos})")
        ax.set_title(f"GSE133604 scRNA-seq: {label}", fontsize=10)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#c0392b", "#f1948a", "#95a5a6")]
        ax.legend(handles, ["FDR<0.05", "FDR<0.25", "n.s."], fontsize=8, loc="lower right")
        fig.tight_layout()
        fig.savefig(OUT / f"nes_{label}.png", dpi=200)
        plt.close(fig)
        if label.startswith("Tacstd2"):
            for term in ("CURATED_KERATIN_FAMILY", "CURATED_TIGHT_JUNCTION_CORE",
                         "GO_KERATIN_FILAMENT", "GO_BICELLULAR_TIGHT_JUNCTION",
                         "BARRIER_CORE_TACSTD2_CLDN4"):
                try:
                    res.plot(terms=term)
                    plt.savefig(OUT / f"enrichment_{term}.png", dpi=200, bbox_inches="tight")
                    plt.close("all")
                except Exception as e:
                    print(f"plot failed for {term}: {e}")

    if all_gsea:
        out = pd.concat(all_gsea, ignore_index=True)
        out.to_csv(OUT / "gsea_results.tsv", sep="\t", index=False)
        print(out[["contrast", "Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))


if __name__ == "__main__":
    main()

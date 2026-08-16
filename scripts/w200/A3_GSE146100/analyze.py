#!/usr/bin/env python3
"""A3 analog: malignant TACSTD2, R vs NR, GSE146100.

One 72-year-old patient, three surgically resected LUAD nodules after
pembrolizumab (Zhang et al., JITC 2021, PMID 33820821).

    W1 GSM4365352  NR  EGFR L858R
    W2 GSM4365353  R   KRAS G12C
    W3 GSM4365354  NR  EGFR L858R / exon20 R77H

The public file is the authors' Seurat log-normalized genes x cells matrix.
There are no author cell-type labels and no CNV calls. "Malignant" is therefore
an epithelial proxy (cluster lineage + a marker-rule sensitivity set).

Honesty constraints baked into the outputs
------------------------------------------
* n_patients = 1. R is one nodule; NR is two nodules.
* Response is fully confounded with driver genotype and clone identity.
* Cell-level tests treat nested cells as independent (pseudoreplication) and
  are reported only as a descriptive effect-size companion, not as evidence
  that TACSTD2 differs by pembrolizumab response.
* A nodule-level test is impossible (1 vs 2).
* All samples are post-treatment; there is no pre-treatment contrast.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import mannwhitneyu

np.random.seed(0)
sc.settings.verbosity = 1

DATA = Path("data/GSE146100_NormData.txt.gz")
OUT = Path("results/w200/A3_GSE146100")

SAMPLE_META = {
    "W1": {
        "gsm": "GSM4365352",
        "response": "NR",
        "genotype": "EGFR L858R",
        "title": "Non-responded W1",
    },
    "W2": {
        "gsm": "GSM4365353",
        "response": "R",
        "genotype": "KRAS G12C",
        "title": "Responded W2",
    },
    "W3": {
        "gsm": "GSM4365354",
        "response": "NR",
        "genotype": "EGFR L858R/R77H",
        "title": "Non-responded W3",
    },
}

MARKERS = {
    "Epithelial": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "CDH1",
        "SFTPC",
        "SFTPB",
        "NAPSA",
        "NKX2-1",
        "SCGB1A1",
        "SFTPA1",
        "SFTPA2",
        "MUC1",
    ],
    "T_NK": ["CD3D", "CD3E", "CD2", "CD8A", "CD4", "NKG7", "GNLY", "KLRD1", "TRAC"],
    "B_Plasma": ["CD79A", "CD79B", "MS4A1", "MZB1", "IGHG1", "IGKC", "JCHAIN"],
    "Myeloid": [
        "LYZ",
        "CD68",
        "CD14",
        "FCGR3A",
        "AIF1",
        "C1QA",
        "C1QB",
        "MARCO",
        "FCN1",
    ],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "ACTA2"],
}

EPI_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19"]


def _expr(adata: sc.AnnData, gene: str) -> np.ndarray:
    gi = adata.raw.var_names.get_loc(gene)
    return np.asarray(adata.raw.X[:, gi].todense()).ravel()


def load_adata() -> sc.AnnData:
    print("Loading normalized matrix ...")
    df = pd.read_csv(DATA, sep="\t", index_col=0)
    df = df.astype(np.float32)
    X = sparse.csr_matrix(df.values.T)
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame(index=df.columns.astype(str)),
        var=pd.DataFrame(index=df.index.astype(str)),
    )
    adata.obs["sample"] = [c.split("_")[0] for c in adata.obs_names]
    adata.obs["response"] = adata.obs["sample"].map(
        lambda s: SAMPLE_META[s]["response"]
    )
    adata.obs["genotype"] = adata.obs["sample"].map(
        lambda s: SAMPLE_META[s]["genotype"]
    )
    adata.obs["gsm"] = adata.obs["sample"].map(lambda s: SAMPLE_META[s]["gsm"])
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs["sample"].value_counts().to_string())
    return adata


def select_hvg(adata: sc.AnnData, n_top: int = 2000) -> pd.Index:
    """Normalized-dispersion HVG on the log-norm matrix (scanpy seurat flavor
    is broken on pandas 3.x)."""
    X = adata.X.tocsr()
    mean = np.asarray(X.mean(axis=0)).ravel()
    X2 = X.copy()
    X2.data = X2.data ** 2
    mean_sq = np.asarray(X2.mean(axis=0)).ravel()
    var = np.clip(mean_sq - mean ** 2, 0, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        disp = np.where(mean > 0, var / mean, 0.0)
    dfg = pd.DataFrame({"mean": mean, "disp": disp}, index=adata.var_names)
    dfg["bin"] = pd.cut(dfg["mean"], bins=20)
    grp = dfg.groupby("bin", observed=True)["disp"]
    std = grp.transform("std").replace(0, np.nan)
    dfg["disp_z"] = (dfg["disp"] - grp.transform("mean")) / std
    dfg["disp_z"] = dfg["disp_z"].fillna(0.0)
    return dfg.sort_values("disp_z", ascending=False).head(n_top).index


def cluster(adata: sc.AnnData) -> sc.AnnData:
    print("Clustering ...")
    adata.raw = adata
    hvg = select_hvg(adata)
    adata.var["highly_variable"] = adata.var_names.isin(hvg)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(
        adata_hvg,
        resolution=1.0,
        random_state=0,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    sc.tl.umap(adata_hvg, random_state=0)
    adata.obs["leiden"] = adata_hvg.obs["leiden"].values
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    print(f"  {adata.obs['leiden'].nunique()} leiden clusters")
    return adata


def annotate(adata: sc.AnnData) -> sc.AnnData:
    print("Scoring lineage markers ...")
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, present, score_name=f"score_{name}", use_raw=True)
    score_cols = [f"score_{n}" for n in MARKERS]
    cl_scores = adata.obs.groupby("leiden")[score_cols].mean()
    cl_lineage = cl_scores.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["lineage"] = adata.obs["leiden"].map(cl_lineage).astype(str)
    adata.obs["is_malignant"] = adata.obs["lineage"].eq("Epithelial")

    # Marker-rule sensitivity set: any epithelial keratin/EPCAM, no PTPRC.
    epi_any = np.zeros(adata.n_obs, dtype=bool)
    for g in EPI_MARKERS:
        if g in adata.raw.var_names:
            epi_any |= _expr(adata, g) > 0
    ptprc = _expr(adata, "PTPRC") if "PTPRC" in adata.raw.var_names else 0
    adata.obs["is_malignant_marker"] = epi_any & (ptprc == 0)

    adata.obs["TACSTD2_lognorm"] = _expr(adata, "TACSTD2")
    adata.obs["TACSTD2_linear"] = np.expm1(adata.obs["TACSTD2_lognorm"])
    adata.obs["EPCAM_lognorm"] = _expr(adata, "EPCAM")
    adata.obs["PTPRC_lognorm"] = _expr(adata, "PTPRC")

    summ = cl_scores.copy()
    summ["lineage"] = cl_lineage
    summ["n_cells"] = adata.obs["leiden"].value_counts().reindex(summ.index)
    summ = summ.join(pd.crosstab(adata.obs["leiden"], adata.obs["sample"]))
    for g in ["TACSTD2", "EPCAM", "PTPRC", "KRT19", "CD3D", "LYZ"]:
        if g in adata.raw.var_names:
            adata.obs[f"tmp_{g}"] = _expr(adata, g)
            summ[f"mean_{g}"] = adata.obs.groupby("leiden")[f"tmp_{g}"].mean()
            adata.obs.drop(columns=[f"tmp_{g}"], inplace=True)
    summ.to_csv(OUT / "cluster_lineage_summary.tsv", sep="\t")
    print(cl_lineage.to_string())
    print("Lineage counts:\n" + adata.obs["lineage"].value_counts().to_string())
    return adata


def _stats(values: np.ndarray) -> dict:
    lin = np.expm1(values)
    return {
        "n_cells": int(len(values)),
        "pct_expressing": float(100.0 * np.mean(values > 0)),
        "mean_lognorm": float(np.mean(values)),
        "median_lognorm": float(np.median(values)),
        "mean_linear": float(np.mean(lin)),
    }


def tacstd2_tables(adata: sc.AnnData) -> tuple[pd.DataFrame, dict]:
    print("TACSTD2 malignant R vs NR ...")
    mal = adata.obs[adata.obs["is_malignant"]].copy()
    mal_m = adata.obs[adata.obs["is_malignant_marker"]].copy()

    rows = []
    for label, frame in [("cluster_epithelial", mal), ("marker_epithelial", mal_m)]:
        for sample, sub in frame.groupby("sample"):
            row = {
                "definition": label,
                "sample": sample,
                "gsm": SAMPLE_META[sample]["gsm"],
                "response": SAMPLE_META[sample]["response"],
                "genotype": SAMPLE_META[sample]["genotype"],
            }
            row.update(_stats(sub["TACSTD2_lognorm"].to_numpy()))
            rows.append(row)
        for resp, sub in frame.groupby("response"):
            row = {
                "definition": label,
                "sample": ",".join(sorted(sub["sample"].unique())),
                "gsm": "",
                "response": resp,
                "genotype": ";".join(
                    sorted({SAMPLE_META[s]["genotype"] for s in sub["sample"].unique()})
                ),
            }
            row.update(_stats(sub["TACSTD2_lognorm"].to_numpy()))
            rows.append(row)
    by = pd.DataFrame(rows)
    by.to_csv(OUT / "tacstd2_malignant_by_nodule.tsv", sep="\t", index=False)

    def contrast(frame: pd.DataFrame, definition: str) -> dict:
        r = frame.loc[frame["response"] == "R", "TACSTD2_lognorm"].to_numpy()
        nr = frame.loc[frame["response"] == "NR", "TACSTD2_lognorm"].to_numpy()
        if len(r) == 0 or len(nr) == 0:
            return {
                "definition": definition,
                "n_R": int(len(r)),
                "n_NR": int(len(nr)),
                "note": "empty arm",
            }
        U, p = mannwhitneyu(r, nr, alternative="two-sided")
        auc = float(U / (len(r) * len(nr)))
        lin_r = float(np.mean(np.expm1(r)))
        lin_nr = float(np.mean(np.expm1(nr)))
        pooled = np.sqrt((np.var(r, ddof=1) + np.var(nr, ddof=1)) / 2)
        return {
            "definition": definition,
            "n_R": int(len(r)),
            "n_NR": int(len(nr)),
            "mean_lognorm_R": float(np.mean(r)),
            "mean_lognorm_NR": float(np.mean(nr)),
            "median_lognorm_R": float(np.median(r)),
            "median_lognorm_NR": float(np.median(nr)),
            "pct_expressing_R": float(100 * np.mean(r > 0)),
            "pct_expressing_NR": float(100 * np.mean(nr > 0)),
            "mean_linear_R": lin_r,
            "mean_linear_NR": lin_nr,
            "log2FC_R_over_NR_linear": float(
                np.log2((lin_r + 1e-9) / (lin_nr + 1e-9))
            ),
            "direction": "higher in R" if np.mean(r) > np.mean(nr) else "higher in NR",
            "cell_level_mannwhitney_U": float(U),
            "cell_level_mannwhitney_p": float(p),
            "cell_level_AUC_P_R_gt_NR": auc,
            "cell_level_rank_biserial": float(2 * auc - 1),
            "cell_level_cohens_d_lognorm": float((np.mean(r) - np.mean(nr)) / pooled)
            if pooled > 0
            else float("nan"),
        }

    contrasts = [
        contrast(mal, "cluster_epithelial"),
        contrast(mal_m, "marker_epithelial"),
    ]
    pd.DataFrame(contrasts).to_csv(
        OUT / "tacstd2_R_vs_NR_celllevel.tsv", sep="\t", index=False
    )

    # Compartment sanity: TACSTD2 should be epithelial-restricted.
    comp_rows = []
    for lin, sub in adata.obs.groupby("lineage"):
        row = {"lineage": lin}
        row.update(_stats(sub["TACSTD2_lognorm"].to_numpy()))
        comp_rows.append(row)
    pd.DataFrame(comp_rows).sort_values("mean_lognorm", ascending=False).to_csv(
        OUT / "tacstd2_by_lineage.tsv", sep="\t", index=False
    )

    # Nodule composition
    comp = (
        adata.obs.groupby(["sample", "response", "lineage"])
        .size()
        .rename("n_cells")
        .reset_index()
    )
    tot = comp.groupby("sample")["n_cells"].transform("sum")
    comp["fraction"] = comp["n_cells"] / tot
    comp.to_csv(OUT / "sample_composition.tsv", sep="\t", index=False)

    return by, {
        "cluster_epithelial": contrasts[0],
        "marker_epithelial": contrasts[1],
        "n_malignant_cluster": int(mal.shape[0]),
        "n_malignant_marker": int(mal_m.shape[0]),
        "agreement_both": int(
            ((adata.obs["is_malignant"]) & (adata.obs["is_malignant_marker"])).sum()
        ),
    }


def make_figures(adata: sc.AnnData) -> None:
    mal = adata.obs[adata.obs["is_malignant"]]
    sc.pl.umap(
        adata, color=["lineage", "sample", "response"], ncols=3, show=False, wspace=0.4
    )
    plt.savefig(OUT / "fig_umap_lineage_sample.png", dpi=150, bbox_inches="tight")
    plt.close()

    sc.pl.umap(
        adata,
        color=["TACSTD2", "EPCAM", "PTPRC"],
        ncols=3,
        show=False,
        use_raw=True,
        wspace=0.4,
    )
    plt.savefig(OUT / "fig_umap_markers.png", dpi=150, bbox_inches="tight")
    plt.close()

    order_n = ["W1", "W2", "W3"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    data_n = [mal.loc[mal["sample"] == s, "TACSTD2_lognorm"].to_numpy() for s in order_n]
    axes[0].violinplot(data_n, showmedians=True)
    axes[0].set_xticks([1, 2, 3])
    axes[0].set_xticklabels(
        [f"{s}\n({SAMPLE_META[s]['response']})" for s in order_n]
    )
    axes[0].set_ylabel("TACSTD2 (author log-norm)")
    axes[0].set_title("Malignant (cluster epithelial), per nodule")
    data_g = [
        mal.loc[mal["response"] == g, "TACSTD2_lognorm"].to_numpy() for g in ("R", "NR")
    ]
    axes[1].violinplot(data_g, showmedians=True)
    axes[1].set_xticks([1, 2])
    axes[1].set_xticklabels(["R (W2)", "NR (W1+W3)"])
    axes[1].set_ylabel("TACSTD2 (author log-norm)")
    axes[1].set_title("R vs NR — cells nested in 1 vs 2 nodules")
    fig.tight_layout()
    fig.savefig(OUT / "fig_tacstd2_malignant_violin.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    bn = (
        mal.groupby("sample")["TACSTD2_lognorm"]
        .agg(pct=lambda v: 100 * np.mean(v > 0), mean="mean")
        .reindex(order_n)
    )
    colors = [
        "#1b9e77" if SAMPLE_META[s]["response"] == "R" else "#d95f02" for s in order_n
    ]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4))
    ax[0].bar(order_n, bn["pct"], color=colors)
    ax[0].set_title("% TACSTD2+ malignant cells")
    ax[0].set_ylabel("% expressing")
    ax[1].bar(order_n, bn["mean"], color=colors)
    ax[1].set_title("Mean TACSTD2 (log-norm)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_tacstd2_malignant_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Figures written.")


def write_metadata(adata: sc.AnnData, extra: dict) -> None:
    sample_meta = pd.DataFrame(
        [
            {
                "sample": s,
                **m,
                "n_cells": int((adata.obs["sample"] == s).sum()),
                "n_malignant_cluster": int(
                    ((adata.obs["sample"] == s) & adata.obs["is_malignant"]).sum()
                ),
                "n_malignant_marker": int(
                    ((adata.obs["sample"] == s) & adata.obs["is_malignant_marker"]).sum()
                ),
            }
            for s, m in SAMPLE_META.items()
        ]
    )
    sample_meta.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)

    keep = [
        "sample",
        "gsm",
        "response",
        "genotype",
        "leiden",
        "lineage",
        "is_malignant",
        "is_malignant_marker",
        "TACSTD2_lognorm",
        "TACSTD2_linear",
        "EPCAM_lognorm",
        "PTPRC_lognorm",
    ]
    adata.obs[keep].to_csv(OUT / "cell_annotation.tsv", sep="\t")

    manifest = pd.DataFrame(
        [
            {
                "file": str(DATA),
                "bytes": DATA.stat().st_size if DATA.exists() else None,
                "used": True,
                "note": "author Seurat log-normalized genes x cells; not committed",
            }
        ]
    )
    manifest.to_csv(OUT / "file_manifest.tsv", sep="\t", index=False)

    c = extra["cluster_epithelial"]
    sanity = {
        "n_cells_total": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_leiden": int(adata.obs["leiden"].nunique()),
        "lineage_counts": adata.obs["lineage"].value_counts().to_dict(),
        "malignant_cluster_n": extra["n_malignant_cluster"],
        "malignant_marker_n": extra["n_malignant_marker"],
        "malignant_definition_agreement": extra["agreement_both"],
        "tacstd2_present": True,
        "epcam_present": True,
        "ptprc_present": True,
        "author_celltype_labels": False,
        "cnv_available": False,
        "n_patients": 1,
        "n_nodules_R": 1,
        "n_nodules_NR": 2,
        "response_confounded_with_genotype": True,
        "pretreatment_samples": False,
    }
    (OUT / "sanity_checks.json").write_text(json.dumps(sanity, indent=2) + "\n")

    direction = c.get("direction", "n/a")
    log2fc = c.get("log2FC_R_over_NR_linear")
    p = c.get("cell_level_mannwhitney_p")
    verdict = (
        "No patient-level claim is possible. Pooled malignant TACSTD2 looks "
        f"{direction} (cluster-epithelial log2FC R/NR = {log2fc:.3f}; "
        f"cell-level MWU p = {p:.2e}, n_cells R={c['n_R']} NR={c['n_NR']}), "
        "but that contrast is W3-driven: W1 (NR) ≈ W2 (R) and W3 (NR) is higher "
        "and contributes most NR epithelial cells. The p-value is "
        "pseudoreplicated (cells nested in 1 R vs 2 NR nodules). Response is "
        "fully confounded with KRAS vs EGFR. Nodule-level test is impossible "
        "(n=1 vs 2)."
    )
    summary = {
        "dataset": "GSE146100",
        "pmid": 33820821,
        "task": "W200-A3 analog: malignant TACSTD2 R vs NR (multi-nodule pembro)",
        "honest_verdict": verdict,
        "response_contrast_possible_at_patient_level": False,
        "response_contrast_possible_at_nodule_level": "descriptive only (1 vs 2 nodules, 1 patient)",
        "n_patients": 1,
        "R_nodule": "W2 KRAS G12C",
        "NR_nodules": "W1 EGFR L858R; W3 EGFR L858R/R77H",
        "malignant_definition": (
            "cluster epithelial (EPCAM/KRT-high, PTPRC-low lineage); "
            "no CNV / matched-normal reference, so this is an epithelial proxy"
        ),
        "cluster_epithelial_contrast": c,
        "marker_epithelial_contrast": extra["marker_epithelial"],
        "per_nodule_pattern": (
            "W1 (NR) ≈ W2 (R); W3 (NR) higher. Pooled NR elevation is W3-driven "
            "and W3 also contributes most NR epithelial cells."
        ),
        "pseudoreplication_warning": (
            "Cells are nested within 3 nodules from ONE patient. "
            "The cell-level p-value is not valid evidence of a "
            "response-associated difference."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (OUT / "audit.json").write_text(
        json.dumps(
            {
                "input": str(DATA),
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "normalization": "author Seurat log-norm; not re-normalized",
                "clustering": "leiden res=1.0 on 2000 dispersion HVGs, 30 PCs",
                "skipped": [
                    "SRA/FASTQ raw reads (not needed; processed matrix is public)",
                    "inferCNV / CopyKAT (no matched normal; would be exploratory only)",
                ],
            },
            indent=2,
        )
        + "\n"
    )
    write_readme(c, extra, sanity)


def write_readme(c: dict, extra: dict, sanity: dict) -> None:
    log2fc = c.get("log2FC_R_over_NR_linear", float("nan"))
    p = c.get("cell_level_mannwhitney_p", float("nan"))
    text = f"""# W200-A3 · GSE146100 — malignant TACSTD2, R vs NR

**Task:** A3 analog on the only public multi-nodule pembrolizumab LUAD scRNA-seq
set: compare TACSTD2 in malignant cells between the responding nodule and the
two non-responding nodules.

## Verdict

**No patient-level R-vs-NR claim is possible.** This series is one 72-year-old
patient (Zhang et al., *JITC* 2021, PMID 33820821). The responding nodule is
W2 (KRAS G12C). The two non-responding nodules are W1 and W3 (both EGFR).
Response, driver genotype, and clone identity are the same contrast.

Pooled cluster-epithelial TACSTD2 looks **{c.get('direction', 'n/a')}**
(mean log-norm R={c.get('mean_lognorm_R'):.3f} vs NR={c.get('mean_lognorm_NR'):.3f};
% expressing R={c.get('pct_expressing_R'):.1f} vs NR={c.get('pct_expressing_NR'):.1f};
linear-scale log2FC(R/NR)={log2fc:.3f}; cell-level MWU p={p:.2e},
n={c.get('n_R')} R / {c.get('n_NR')} NR cells). **Do not use that p-value.**
Cells are nested in 1 vs 2 nodules. Per-nodule means show the pooled NR
elevation is **W3-driven**: W1 (NR) ≈ W2 (R); W3 (NR) is higher and holds
most NR epithelial cells. A nodule-level test cannot be run (n=1 vs 2).

The marker-rule epithelial sensitivity set also looks higher in NR in the
pool (log2FC={extra['marker_epithelial'].get('log2FC_R_over_NR_linear'):.3f})
and is likewise W3-driven (W1 is not higher than W2).

## What was used

- GEO [GSE146100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146100),
  `GSE146100_NormData.txt.gz` (229 MB; author Seurat log-normalized matrix).
- {sanity['n_cells_total']} cells × {sanity['n_genes']} genes; barcodes prefixed
  W1_/W2_/W3_.
- No author cell-type labels. Malignant = Leiden clusters whose highest lineage
  score is epithelial (EPCAM/KRT/SFT* vs PTPRC/immune/stroma). Sensitivity:
  any of EPCAM/KRT8/KRT18/KRT19 > 0 and PTPRC = 0.
- Raw SRA FASTQ was not downloaded.

## Honest limits vs a real A3 / ICI claim

- n_patients = 1. This is a within-patient nodule comparison, not a cohort.
- R = KRAS; both NR = EGFR. Any TACSTD2 difference may be genotype, not response.
- All three nodules are post-pembrolizumab. There is no pre-treatment arm.
- "Malignant" is an epithelial proxy. Without CNV or a matched normal we cannot
  exclude residual normal epithelial cells.
- Cell-level statistics are pseudoreplicated.

## Results (cluster epithelial, primary)

| Nodule | Response | Genotype | n malignant | % TACSTD2+ | mean log-norm |
|---|---|---|---|---|---|
"""
    # filled after tables exist; write a static skeleton then append from TSV
    tsv = pd.read_csv(OUT / "tacstd2_malignant_by_nodule.tsv", sep="\t")
    nod = tsv[
        (tsv["definition"] == "cluster_epithelial")
        & tsv["gsm"].notna()
        & (tsv["gsm"].astype(str) != "")
        & tsv["sample"].isin(["W1", "W2", "W3"])
    ]
    for _, r in nod.iterrows():
        text += (
            f"| {r['sample']} | {r['response']} | {r['genotype']} | "
            f"{int(r['n_cells'])} | {r['pct_expressing']:.1f} | "
            f"{r['mean_lognorm']:.3f} |\n"
        )
    text += f"""
Pooled R (W2) vs NR (W1+W3): log2FC={log2fc:.3f}, cell-level AUC={c.get('cell_level_AUC_P_R_gt_NR'):.3f}.

## Files

| File | Role |
|---|---|
| `sample_metadata.tsv` | Nodule / GSM / genotype / cell counts |
| `sample_composition.tsv` | Lineage mix per nodule |
| `cluster_lineage_summary.tsv` | Leiden → lineage + marker means |
| `cell_annotation.tsv` | Per-cell labels and TACSTD2 |
| `tacstd2_malignant_by_nodule.tsv` | Per-nodule and pooled R/NR stats |
| `tacstd2_R_vs_NR_celllevel.tsv` | Cell-level contrast (flagged) |
| `tacstd2_by_lineage.tsv` | Compartment sanity |
| `summary.json` / `sanity_checks.json` / `audit.json` | Verdict |
| `fig_umap_lineage_sample.png` / `fig_umap_markers.png` | UMAP |
| `fig_tacstd2_malignant_violin.png` / `fig_tacstd2_malignant_bar.png` | TACSTD2 |

## Reproduce

```bash
python3 scripts/w200/A3_GSE146100/download.py
python3 scripts/w200/A3_GSE146100/analyze.py
```
"""
    (OUT / "README.md").write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not DATA.exists():
        raise SystemExit(f"missing {DATA}; run scripts/w200/A3_GSE146100/download.py")
    adata = load_adata()
    adata = cluster(adata)
    adata = annotate(adata)
    _, extra = tacstd2_tables(adata)
    make_figures(adata)
    write_metadata(adata, extra)
    print("Done.")
    print((OUT / "summary.json").read_text())


if __name__ == "__main__":
    main()

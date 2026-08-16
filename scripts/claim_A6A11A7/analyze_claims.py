#!/usr/bin/env python3
"""Honest public-data tests for claims A6, A11, and the A7 paired-ICI hunt."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/tmp/data")
OUT = Path("/workspace/results/claim_A6A11A7")
OUT.mkdir(parents=True, exist_ok=True)

A11_GENES = [
    "LGALS1",
    "LGALS3",
    "LGALS9",
    "NECTIN1",
    "NECTIN2",
    "NECTIN3",
    "NECTIN4",
    "PVR",
    "TGFB1",
    "TGFB2",
    "TGFB3",
    "CD47",
    "SIRPA",
    "CD274",
    "EPCAM",
    "CLDN4",
]
CD8_SUBTYPES = {
    "Exhausted CD8+ T",
    "Naive CD8+ T",
    "Cytotoxic CD8+ T",
    "CD8 low T",
}
NK_SUBTYPES = {"NK"}
MALIGNANT = {"Malignant cells"}
EPITHELIAL_TYPES = {"Epithelial cells"}


def spearman(x, y) -> dict:
    x = pd.Series(x).astype(float)
    y = pd.Series(y).astype(float)
    mask = x.notna() & y.notna()
    if mask.sum() < 5:
        return {"n": int(mask.sum()), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[mask], y[mask])
    return {"n": int(mask.sum()), "rho": float(rho), "p": float(p)}


def mw_and_median(high, low) -> dict:
    high = pd.Series(high).dropna().astype(float)
    low = pd.Series(low).dropna().astype(float)
    if len(high) < 3 or len(low) < 3:
        return {
            "n_high": int(len(high)),
            "n_low": int(len(low)),
            "median_high": float(high.median()) if len(high) else np.nan,
            "median_low": float(low.median()) if len(low) else np.nan,
            "p": np.nan,
        }
    p = float(stats.mannwhitneyu(high, low, alternative="two-sided").pvalue)
    return {
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "median_high": float(high.median()),
        "median_low": float(low.median()),
        "p": p,
    }


def partial_spearman(x, y, z) -> dict:
    """Spearman of residuals after ranking and linear residualization on z."""
    df = pd.DataFrame({"x": x, "y": y, "z": z}).dropna().astype(float)
    if len(df) < 8:
        return {"n": int(len(df)), "rho": np.nan, "p": np.nan}
    rx = df["x"].rank()
    ry = df["y"].rank()
    rz = df["z"].rank()
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    resx = rx - (bx[0] * rz + bx[1])
    resy = ry - (by[0] * rz + by[1])
    rho, p = stats.spearmanr(resx, resy)
    return {"n": int(len(df)), "rho": float(rho), "p": float(p)}


def lognorm(umi: pd.Series, n_umi: pd.Series) -> pd.Series:
    scale = 1e4 / n_umi.replace(0, np.nan)
    return np.log1p(umi * scale)


def load_extracted(prefix: str) -> tuple[pd.DataFrame, pd.Series]:
    expr = pd.read_csv(DATA / "extracted" / f"{prefix}.genes.csv.gz", index_col=0)
    lib = pd.read_csv(DATA / "extracted" / f"{prefix}.lib.csv.gz", index_col=0)["n_umi"]
    return expr, lib


def analyze_gse131907(expr: pd.DataFrame, lib: pd.Series) -> dict:
    ann = pd.read_csv(DATA / "GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep="\t")
    ann = ann.set_index("Index")
    common = expr.index.intersection(ann.index)
    ann = ann.loc[common].copy()
    expr = expr.loc[common]
    lib = lib.loc[common]
    tac = lognorm(expr["TACSTD2"], lib)
    ann["tacstd2"] = tac
    ann["is_malignant"] = ann["Cell_subtype"].isin(MALIGNANT)
    ann["is_epithelial"] = ann["Cell_type"].isin(EPITHELIAL_TYPES)
    ann["is_cd8"] = ann["Cell_subtype"].isin(CD8_SUBTYPES)
    ann["is_nk"] = ann["Cell_subtype"].isin(NK_SUBTYPES)
    ann["is_t"] = ann["Cell_type"].eq("T lymphocytes")
    ann["is_immune"] = ann["Cell_type"].isin(
        ["T lymphocytes", "NK cells", "B lymphocytes", "Myeloid cells", "MAST cells"]
    )

    tumor_origins = ["tLung", "tL/B", "mLN", "mBrain", "PE"]
    rows = []
    for sample, g in ann.groupby("Sample"):
        origin = g["Sample_Origin"].iloc[0]
        n = len(g)
        malig = g.loc[g["is_malignant"], "tacstd2"]
        epi = g.loc[g["is_epithelial"], "tacstd2"]
        rows.append(
            {
                "sample": sample,
                "origin": origin,
                "n_cells": n,
                "tumor_frac": float(g["is_malignant"].mean()),
                "epithelial_frac": float(g["is_epithelial"].mean()),
                "cd8_frac": float(g["is_cd8"].mean()),
                "nk_frac": float(g["is_nk"].mean()),
                "t_frac": float(g["is_t"].mean()),
                "immune_frac": float(g["is_immune"].mean()),
                "cd8_among_immune": float(g.loc[g["is_immune"], "is_cd8"].mean())
                if g["is_immune"].any()
                else np.nan,
                "nk_among_immune": float(g.loc[g["is_immune"], "is_nk"].mean())
                if g["is_immune"].any()
                else np.nan,
                "tacstd2_all": float(g["tacstd2"].mean()),
                "tacstd2_malignant": float(malig.mean()) if len(malig) >= 10 else np.nan,
                "tacstd2_epithelial": float(epi.mean()) if len(epi) >= 10 else np.nan,
                "n_malignant": int(g["is_malignant"].sum()),
                "n_epithelial": int(g["is_epithelial"].sum()),
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(OUT / "A6_GSE131907_sample_composition.csv", index=False)

    results = {"dataset": "GSE131907", "n_cells": int(len(ann)), "n_samples": int(sample_df["sample"].nunique())}
    contrasts = {
        "primary_tumor_tLung": sample_df[sample_df["origin"].eq("tLung")],
        "all_tumor_origins": sample_df[sample_df["origin"].isin(tumor_origins)],
    }
    claim_rows = []
    for label, sdf in contrasts.items():
        sdf = sdf.copy()
        if len(sdf) < 6:
            continue
        for score_name in ["tacstd2_all", "tacstd2_malignant", "tacstd2_epithelial"]:
            scored = sdf.dropna(subset=[score_name])
            if len(scored) < 6:
                continue
            med = scored[score_name].median()
            high = scored[scored[score_name] >= med]
            low = scored[scored[score_name] < med]
            for metric in [
                "tumor_frac",
                "epithelial_frac",
                "cd8_frac",
                "nk_frac",
                "immune_frac",
                "cd8_among_immune",
                "nk_among_immune",
            ]:
                mw = mw_and_median(high[metric], low[metric])
                sp = spearman(scored[score_name], scored[metric])
                claim_rows.append(
                    {
                        "dataset": "GSE131907",
                        "subset": label,
                        "trop2_definition": score_name,
                        "metric": metric,
                        "n_samples": int(len(scored)),
                        "median_split": float(med),
                        **{f"mw_{k}": v for k, v in mw.items()},
                        **{f"spearman_{k}": v for k, v in sp.items()},
                    }
                )
        # composition-adjusted tests
        if sdf["tacstd2_all"].notna().sum() >= 8:
            for metric in ["cd8_frac", "nk_frac"]:
                zcol = "epithelial_frac" if label.startswith("primary") else "tumor_frac"
                # tLung has no author malignant calls; use epithelium as the composition covariate
                if sdf[zcol].nunique() < 2:
                    zcol = "epithelial_frac"
                pc = partial_spearman(sdf["tacstd2_all"], sdf[metric], sdf[zcol])
                claim_rows.append(
                    {
                        "dataset": "GSE131907",
                        "subset": label,
                        "trop2_definition": f"tacstd2_all_partial_{zcol}",
                        "metric": metric,
                        "n_samples": pc["n"],
                        "median_split": np.nan,
                        "mw_n_high": np.nan,
                        "mw_n_low": np.nan,
                        "mw_median_high": np.nan,
                        "mw_median_low": np.nan,
                        "mw_p": np.nan,
                        "spearman_n": pc["n"],
                        "spearman_rho": pc["rho"],
                        "spearman_p": pc["p"],
                    }
                )
    claim_tab = pd.DataFrame(claim_rows)
    claim_tab.to_csv(OUT / "A6_GSE131907_tests.csv", index=False)
    results["tests"] = claim_rows

    # A11 in malignant cells and sample-level
    a11_rows = []
    mal = ann[ann["is_malignant"]].copy()
    epi_cells = ann[ann["is_epithelial"]].copy()
    for gene in A11_GENES:
        if gene not in expr.columns:
            continue
        mal[gene] = lognorm(expr.loc[mal.index, gene], lib.loc[mal.index])
        epi_cells[gene] = lognorm(expr.loc[epi_cells.index, gene], lib.loc[epi_cells.index])
        a11_rows.append(
            {
                "level": "malignant_cells",
                "dataset": "GSE131907",
                "gene": gene,
                **spearman(mal["tacstd2"], mal[gene]),
            }
        )
        a11_rows.append(
            {
                "level": "epithelial_cells",
                "dataset": "GSE131907",
                "gene": gene,
                **spearman(epi_cells["tacstd2"], epi_cells[gene]),
            }
        )
    # sample means within malignant cells, primary tumors
    mal_sample = (
        mal.groupby(["Sample", "Sample_Origin"], observed=True)
        .agg({c: "mean" for c in ["tacstd2"] + [g for g in A11_GENES if g in mal.columns]})
        .reset_index()
    )
    mal_sample.to_csv(OUT / "A11_GSE131907_malignant_sample_means.csv", index=False)
    epi = ann[ann["is_epithelial"]].copy()
    epi_sample = (
        epi.groupby(["Sample", "Sample_Origin"], observed=True)
        .agg({c: "mean" for c in ["tacstd2"] + [g for g in A11_GENES if g in epi.columns]})
        .reset_index()
    )
    epi_sample.to_csv(OUT / "A11_GSE131907_epithelial_sample_means.csv", index=False)
    for origin_label, sub in [
        ("tLung_epithelial_sample", epi_sample[epi_sample["Sample_Origin"].eq("tLung")]),
        ("all_tumor_epithelial_sample", epi_sample[epi_sample["Sample_Origin"].isin(tumor_origins)]),
        ("tLung_malignant_sample", mal_sample[mal_sample["Sample_Origin"].eq("tLung")]),
        ("all_tumor_malignant_sample", mal_sample[mal_sample["Sample_Origin"].isin(tumor_origins)]),
    ]:
        for gene in A11_GENES:
            if gene not in sub.columns:
                continue
            a11_rows.append({"level": origin_label, "dataset": "GSE131907", "gene": gene, **spearman(sub["tacstd2"], sub[gene])})
    # all-cell sample means (composition-confounded)
    for gene in A11_GENES:
        if gene not in expr.columns:
            continue
        ann[gene] = lognorm(expr[gene], lib)
    all_sample = ann.groupby("Sample").agg({g: "mean" for g in ["tacstd2"] + [g for g in A11_GENES if g in ann.columns]})
    all_sample = all_sample.join(sample_df.set_index("sample")[["origin"]], how="left")
    for origin_label, sub in [
        ("tLung_allcells_sample", all_sample[all_sample["origin"].eq("tLung")]),
        ("all_tumor_allcells_sample", all_sample[all_sample["origin"].isin(tumor_origins)]),
    ]:
        for gene in A11_GENES:
            if gene not in sub.columns:
                continue
            a11_rows.append({"level": origin_label, "dataset": "GSE131907", "gene": gene, **spearman(sub["tacstd2"], sub[gene])})
    a11_tab = pd.DataFrame(a11_rows)
    a11_tab.to_csv(OUT / "A11_GSE131907_correlations.csv", index=False)
    results["a11_n_malignant_cells"] = int(len(mal))
    return results, sample_df, claim_tab, a11_tab


def annotate_gse207422(expr: pd.DataFrame, lib: pd.Series) -> pd.DataFrame:
    scores = pd.DataFrame(index=expr.index)
    scores["epithelial"] = lognorm(expr.get("EPCAM", 0), lib) + lognorm(expr.get("KRT19", 0), lib) + lognorm(
        expr.get("KRT8", 0), lib
    )
    scores["t"] = lognorm(expr.get("CD3D", 0), lib) + lognorm(expr.get("CD3E", 0), lib)
    scores["nk"] = lognorm(expr.get("NKG7", 0), lib) + lognorm(expr.get("KLRD1", 0), lib) + lognorm(
        expr.get("GNLY", 0), lib
    )
    scores["b"] = lognorm(expr.get("MS4A1", 0), lib) + lognorm(expr.get("CD79A", 0), lib)
    scores["myeloid"] = lognorm(expr.get("LYZ", 0), lib)
    scores["fibroblast"] = lognorm(expr.get("COL1A1", 0), lib)
    scores["endothelial"] = lognorm(expr.get("PECAM1", 0), lib)
    lineage = scores.idxmax(axis=1)
    # NK vs T: if T score is competitive, prefer T
    t_vs_nk = scores["t"] - scores["nk"]
    lineage = lineage.mask((lineage.eq("nk")) & (t_vs_nk > 0), "t")
    out = pd.DataFrame({"lineage": lineage, "sample": expr.index.to_series().str.split("_").str[:2].str.join("_")})
    out["tacstd2"] = lognorm(expr["TACSTD2"], lib)
    out["cd8a"] = lognorm(expr.get("CD8A", 0), lib)
    out["is_epithelial"] = out["lineage"].eq("epithelial")
    out["is_cd8"] = out["lineage"].eq("t") & (out["cd8a"] > 0.5)
    out["is_nk"] = out["lineage"].eq("nk")
    return out


def analyze_gse207422(expr: pd.DataFrame, lib: pd.Series) -> dict:
    meta = pd.read_excel(DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    meta = meta.dropna(subset=["Sample"]).copy()
    cells = annotate_gse207422(expr, lib)
    cells = cells.merge(meta, left_on="sample", right_on="Sample", how="left")
    cells.index = expr.index
    lineage_counts = (
        cells.groupby(["sample", "lineage"]).size().unstack(fill_value=0).reset_index()
    )
    lineage_counts.to_csv(OUT / "A6_GSE207422_lineage_counts.csv", index=False)
    rows = []
    for sample, g in cells.groupby("sample"):
        m = g.iloc[0]
        epi = g.loc[g["is_epithelial"], "tacstd2"]
        rows.append(
            {
                "sample": sample,
                "patient": m.get("Patient"),
                "resource": m.get("Resource"),
                "pathology": m.get("Pathology"),
                "response": m.get("Pathologic Response"),
                "n_cells": int(len(g)),
                "tumor_frac": float(g["is_epithelial"].mean()),
                "cd8_frac": float(g["is_cd8"].mean()),
                "nk_frac": float(g["is_nk"].mean()),
                "tacstd2_all": float(g["tacstd2"].mean()),
                "tacstd2_epithelial": float(epi.mean()) if len(epi) >= 10 else np.nan,
                "n_epithelial": int(g["is_epithelial"].sum()),
            }
        )
    sdf = pd.DataFrame(rows)
    sdf.to_csv(OUT / "A6_GSE207422_sample_composition.csv", index=False)

    tests = []
    scored = sdf.dropna(subset=["tacstd2_all"])
    if len(scored) >= 6:
        med = scored["tacstd2_all"].median()
        high = scored[scored["tacstd2_all"] >= med]
        low = scored[scored["tacstd2_all"] < med]
        for metric in ["tumor_frac", "cd8_frac", "nk_frac"]:
            mw = mw_and_median(high[metric], low[metric])
            sp = spearman(scored["tacstd2_all"], scored[metric])
            tests.append(
                {
                    "dataset": "GSE207422",
                    "trop2_definition": "tacstd2_all",
                    "metric": metric,
                    **mw,
                    **{f"spearman_{k}": v for k, v in sp.items()},
                }
            )
        pc_cd8 = partial_spearman(scored["tacstd2_all"], scored["cd8_frac"], scored["tumor_frac"])
        pc_nk = partial_spearman(scored["tacstd2_all"], scored["nk_frac"], scored["tumor_frac"])
        tests.append({"dataset": "GSE207422", "trop2_definition": "partial_tumor_frac", "metric": "cd8_frac", **pc_cd8})
        tests.append({"dataset": "GSE207422", "trop2_definition": "partial_tumor_frac", "metric": "nk_frac", **pc_nk})
    epi_scored = sdf.dropna(subset=["tacstd2_epithelial"])
    if len(epi_scored) >= 6:
        med = epi_scored["tacstd2_epithelial"].median()
        high = epi_scored[epi_scored["tacstd2_epithelial"] >= med]
        low = epi_scored[epi_scored["tacstd2_epithelial"] < med]
        for metric in ["tumor_frac", "cd8_frac", "nk_frac"]:
            mw = mw_and_median(high[metric], low[metric])
            sp = spearman(epi_scored["tacstd2_epithelial"], epi_scored[metric])
            tests.append(
                {
                    "dataset": "GSE207422",
                    "trop2_definition": "tacstd2_epithelial",
                    "metric": metric,
                    **mw,
                    **{f"spearman_{k}": v for k, v in sp.items()},
                }
            )
    pd.DataFrame(tests).to_csv(OUT / "A6_GSE207422_tests.csv", index=False)

    # A7 analog: unpaired pre vs post TACSTD2
    a7 = []
    for score in ["tacstd2_all", "tacstd2_epithelial", "tumor_frac", "cd8_frac", "nk_frac"]:
        pre = sdf.loc[sdf["resource"].astype(str).str.contains("Pre", na=False), score]
        post = sdf.loc[sdf["resource"].astype(str).str.contains("Post", na=False), score]
        a7.append(
            {
                "score": score,
                "n_pre": int(pre.notna().sum()),
                "n_post": int(post.notna().sum()),
                "median_pre": float(pre.median()) if pre.notna().any() else np.nan,
                "median_post": float(post.median()) if post.notna().any() else np.nan,
                "p_mw": float(stats.mannwhitneyu(pre.dropna(), post.dropna(), alternative="two-sided").pvalue)
                if pre.dropna().size >= 2 and post.dropna().size >= 2
                else np.nan,
                "paired_same_patient": False,
            }
        )
    # response among post-treatment
    post = sdf[sdf["resource"].astype(str).str.contains("Post", na=False)]
    mpr = post[post["response"].isin(["MPR", "pCR", "MPR (pCR)"])]
    nmpr = post[post["response"].eq("NMPR")]
    for score in ["tacstd2_all", "tacstd2_epithelial"]:
        a7.append(
            {
                "score": f"{score}_post_MPR_vs_NMPR",
                "n_pre": int(mpr[score].notna().sum()),
                "n_post": int(nmpr[score].notna().sum()),
                "median_pre": float(mpr[score].median()) if mpr[score].notna().any() else np.nan,
                "median_post": float(nmpr[score].median()) if nmpr[score].notna().any() else np.nan,
                "p_mw": float(
                    stats.mannwhitneyu(mpr[score].dropna(), nmpr[score].dropna(), alternative="two-sided").pvalue
                )
                if mpr[score].dropna().size >= 2 and nmpr[score].dropna().size >= 2
                else np.nan,
                "paired_same_patient": False,
            }
        )
    pd.DataFrame(a7).to_csv(OUT / "A7_GSE207422_prepost_unpaired.csv", index=False)

    # A11 epithelial cells
    a11 = []
    epi_cells = cells[cells["is_epithelial"]].copy()
    for gene in A11_GENES:
        if gene not in expr.columns:
            continue
        epi_cells[gene] = lognorm(expr.loc[epi_cells.index, gene], lib.loc[epi_cells.index])
        a11.append({"level": "epithelial_cells", "dataset": "GSE207422", "gene": gene, **spearman(epi_cells["tacstd2"], epi_cells[gene])})
    pd.DataFrame(a11).to_csv(OUT / "A11_GSE207422_epithelial_correlations.csv", index=False)

    # bulk pre-treatment TACSTD2 vs response
    bulk_meta = pd.read_excel(DATA / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx").dropna(subset=["Sample", "Patient"])
    bulk = pd.read_csv(DATA / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t")
    bulk = bulk.set_index(bulk.columns[0])
    keep = [c for c in bulk.columns if c in set(bulk_meta["Sample"])]
    genes_bulk = ["TACSTD2", "EPCAM", "CD8A", "NKG7"] + [g for g in A11_GENES if g in bulk.index]
    gene_rows = [g for g in dict.fromkeys(genes_bulk) if g in bulk.index]
    bmat = bulk.loc[gene_rows, keep]
    bmat = bmat[~bmat.index.duplicated(keep="first")]
    bsub = bmat.T
    bsub.index.name = "Sample"
    bsub = bsub.merge(bulk_meta.set_index("Sample"), left_index=True, right_index=True, how="left")
    bsub.to_csv(OUT / "A7_GSE207422_bulk_pretreatment.csv")
    bulk_tests = []
    if "TACSTD2" in bsub.columns:
        mpr = bsub[bsub["Pathologic Response"].isin(["MPR", "pCR", "MPR (pCR)"])]["TACSTD2"]
        nmpr = bsub[bsub["Pathologic Response"].eq("NMPR")]["TACSTD2"]
        bulk_tests.append(
            {
                "comparison": "pre_bulk_TACSTD2_MPR_vs_NMPR",
                "n_mpr": int(mpr.notna().sum()),
                "n_nmpr": int(nmpr.notna().sum()),
                "median_mpr": float(mpr.median()) if len(mpr) else np.nan,
                "median_nmpr": float(nmpr.median()) if len(nmpr) else np.nan,
                "p": float(stats.mannwhitneyu(mpr.dropna(), nmpr.dropna(), alternative="two-sided").pvalue)
                if mpr.dropna().size >= 2 and nmpr.dropna().size >= 2
                else np.nan,
            }
        )
        for gene in A11_GENES:
            if gene in bsub.columns:
                bulk_tests.append({"comparison": f"pre_bulk_TACSTD2_vs_{gene}", **spearman(bsub["TACSTD2"], bsub[gene])})
    pd.DataFrame(bulk_tests).to_csv(OUT / "A7_GSE207422_bulk_tests.csv", index=False)
    return {"n_sc_samples": int(len(sdf)), "n_bulk_samples": int(len(bsub)), "tests": tests, "a7": a7}


def analyze_tcga() -> pd.DataFrame:
    rows = []
    for cohort, fname in [("LUAD", "TCGA.LUAD.HiSeqV2.gz"), ("LUSC", "TCGA.LUSC.HiSeqV2.gz")]:
        mat = pd.read_csv(DATA / fname, sep="\t", index_col=0)
        tumors = [c for c in mat.columns if str(c).endswith("-01")]
        genes = ["TACSTD2", "CD8A", "NKG7", "NCAM1"] + A11_GENES
        present = [g for g in genes if g in mat.index]
        sub = mat.loc[present, tumors].T
        sub.to_csv(OUT / f"A11_TCGA_{cohort}_tumor_gene_panel.csv")
        for gene in present:
            if gene == "TACSTD2":
                continue
            rows.append({"dataset": f"TCGA-{cohort}", "level": "bulk_tumor", "gene": gene, **spearman(sub["TACSTD2"], sub[gene])})
        if "EPCAM" in sub.columns:
            for gene in ["CD8A", "NKG7", "LGALS9", "NECTIN2", "PVR", "TGFB1", "CD47"]:
                if gene in sub.columns:
                    pc = partial_spearman(sub["TACSTD2"], sub[gene], sub["EPCAM"])
                    rows.append(
                        {
                            "dataset": f"TCGA-{cohort}",
                            "level": "bulk_tumor_partial_EPCAM",
                            "gene": gene,
                            **pc,
                        }
                    )
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "A11_TCGA_correlations.csv", index=False)
    return tab


def make_plots(gse131907_samples: pd.DataFrame) -> None:
    tlung = gse131907_samples[gse131907_samples["origin"].eq("tLung")].copy()
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, y, title in [
        (axes[0], "epithelial_frac", "Epithelial fraction"),
        (axes[1], "cd8_frac", "CD8 fraction"),
        (axes[2], "nk_frac", "NK fraction"),
    ]:
        ax.scatter(tlung["tacstd2_all"], tlung[y], c="#1f4e79", s=40, label="all-cell TACSTD2")
        ax.scatter(tlung["tacstd2_epithelial"], tlung[y], c="#c44e52", s=40, marker="D", label="epithelial TACSTD2")
        ax.set_xlabel("TACSTD2 (mean log1p CPM)")
        ax.set_ylabel(title)
        ax.set_title(f"GSE131907 tLung: {title}")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "A6_GSE131907_tLung_scatter.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6))
    for ax, score, title in [
        (axes[0], "tacstd2_all", "Split on all-cell TACSTD2"),
        (axes[1], "tacstd2_epithelial", "Split on epithelial TACSTD2"),
    ]:
        sdf = tlung.dropna(subset=[score])
        med = sdf[score].median()
        sdf = sdf.assign(grp=np.where(sdf[score] >= med, "TROP2-high", "TROP2-low"))
        data = [sdf.loc[sdf["grp"].eq(g), m] for g in ["TROP2-low", "TROP2-high"] for m in ["epithelial_frac", "cd8_frac", "nk_frac"]]
        positions = [1, 2, 3, 5, 6, 7]
        bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True)
        colors = ["#9ecae1", "#fcbba1", "#c7e9c0"] * 2
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
        ax.set_xticks([2, 6])
        ax.set_xticklabels(["TROP2-low", "TROP2-high"])
        ax.set_ylabel("Sample fraction")
        ax.set_title(f"GSE131907 tLung\n{title}")
    fig.tight_layout()
    fig.savefig(OUT / "A6_GSE131907_tLung_boxplots.png", dpi=160)
    plt.close(fig)


def verdicts(a6_131907: pd.DataFrame, a11_131907: pd.DataFrame, a11_tcga: pd.DataFrame) -> dict:
    def grab(df, **kwargs):
        q = df
        for k, v in kwargs.items():
            q = q[q[k].eq(v)]
        return q.iloc[0].to_dict() if len(q) else {}

    tlung_all_tumor = grab(
        a6_131907, subset="primary_tumor_tLung", trop2_definition="tacstd2_all", metric="epithelial_frac"
    )
    tlung_all_cd8 = grab(a6_131907, subset="primary_tumor_tLung", trop2_definition="tacstd2_all", metric="cd8_frac")
    tlung_epi_cd8 = grab(
        a6_131907, subset="primary_tumor_tLung", trop2_definition="tacstd2_epithelial", metric="cd8_frac"
    )
    tlung_partial_cd8 = grab(
        a6_131907, subset="primary_tumor_tLung", trop2_definition="tacstd2_all_partial_epithelial_frac", metric="cd8_frac"
    )
    a11_mal = a11_131907[a11_131907["level"].eq("malignant_cells")]
    a11_tlung_epi = a11_131907[a11_131907["level"].eq("tLung_epithelial_sample")]
    return {
        "A6_tLung_allcell_TACSTD2_vs_tumor_frac_rho": tlung_all_tumor.get("spearman_rho"),
        "A6_tLung_allcell_TACSTD2_vs_tumor_frac_p": tlung_all_tumor.get("spearman_p"),
        "A6_tLung_allcell_TACSTD2_vs_cd8_rho": tlung_all_cd8.get("spearman_rho"),
        "A6_tLung_allcell_TACSTD2_vs_cd8_p": tlung_all_cd8.get("spearman_p"),
        "A6_tLung_epithelial_TACSTD2_vs_cd8_rho": tlung_epi_cd8.get("spearman_rho"),
        "A6_tLung_epithelial_TACSTD2_vs_cd8_p": tlung_epi_cd8.get("spearman_p"),
        "A6_tLung_partial_TACSTD2_cd8_given_tumorfrac_rho": tlung_partial_cd8.get("spearman_rho"),
        "A6_tLung_partial_TACSTD2_cd8_given_tumorfrac_p": tlung_partial_cd8.get("spearman_p"),
        "A11_malignant_cell_top": a11_mal.sort_values("rho", ascending=False).head(8).to_dict(orient="records"),
        "A11_tLung_epithelial_sample": a11_tlung_epi.to_dict(orient="records"),
        "A11_TCGA_LUAD_bulk": a11_tcga[a11_tcga["dataset"].eq("TCGA-LUAD")].to_dict(orient="records"),
    }


def analyze_pretreatment_ici_bulk() -> None:
    """Public ICI RNA that is NOT paired pre/post — recorded so A7 is not overstated."""
    rows = []
    # GSE126044: Cho et al., pretreatment anti-PD-1, 16 NSCLC
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    titles = resp = None
    with gzip.open(DATA / "GSE126044_series_matrix.txt.gz", "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "patient response" in line:
                resp = [x.strip('"').split(":", 1)[-1].strip() for x in line.rstrip().split("\t")[1:]]
    titles = [t.replace("RNA-seq_", "") for t in titles]
    m = pd.DataFrame({"sample": titles, "response": resp})
    if "TACSTD2" in counts.index:
        s = counts.loc["TACSTD2"]
        df = m.merge(s.rename("TACSTD2"), left_on="sample", right_index=True, how="left")
        df.to_csv(OUT / "A7_GSE126044_pretreatment_TACSTD2.csv", index=False)
        r = df[df["response"].eq("responder")]["TACSTD2"]
        nr = df[df["response"].eq("non-responder")]["TACSTD2"]
        rows.append(
            {
                "dataset": "GSE126044",
                "design": "pretreatment_only_not_paired",
                "n_responder": int(r.notna().sum()),
                "n_nonresponder": int(nr.notna().sum()),
                "median_responder": float(r.median()) if len(r) else np.nan,
                "median_nonresponder": float(nr.median()) if len(nr) else np.nan,
                "p": float(stats.mannwhitneyu(r.dropna(), nr.dropna(), alternative="two-sided").pvalue)
                if r.dropna().size >= 2 and nr.dropna().size >= 2
                else np.nan,
            }
        )
    # GSE135222 Jung et al. pretreatment, PFS event
    exp = pd.read_csv(DATA / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t")
    tac = exp[exp.iloc[:, 0].astype(str).str.startswith("ENSG00000184292")]
    if len(tac):
        ser = tac.iloc[0, 1:].astype(float)
        ser.name = "TACSTD2"
        with gzip.open(DATA / "GSE135222_series_matrix.txt.gz", "rt") as handle:
            titles = pfs = None
            for line in handle:
                if line.startswith("!Sample_title"):
                    titles = [x.strip('"').replace(" ", "") for x in line.rstrip().split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1") and "progression-free survival (pfs)" in line and "pfs.time" not in line:
                    pfs = [x.strip('"').split(":", 1)[-1].strip() for x in line.rstrip().split("\t")[1:]]
        m2 = pd.DataFrame({"sample": titles, "pfs_event": pfs})
        m2 = m2.merge(ser, left_on="sample", right_index=True, how="left")
        m2.to_csv(OUT / "A7_GSE135222_pretreatment_TACSTD2.csv", index=False)
        # pfs_event 1 = progressor typically in this series; keep raw
        rows.append(
            {
                "dataset": "GSE135222",
                "design": "pretreatment_only_not_paired",
                "n_samples": int(m2["TACSTD2"].notna().sum()),
                "tacstd2_median": float(m2["TACSTD2"].median()),
                "note": "Ensembl ENSG00000184292; no paired post-ICI RNA",
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "A7_public_ICI_RNA_not_paired.csv", index=False)


def main() -> None:
    expr1319, lib1319 = load_extracted("GSE131907")
    expr2074, lib2074 = load_extracted("GSE207422")
    r1319, samples, a6tab, a11tab = analyze_gse131907(expr1319, lib1319)
    r2074 = analyze_gse207422(expr2074, lib2074)
    a11_tcga = analyze_tcga()
    analyze_pretreatment_ici_bulk()
    make_plots(samples)
    summary = {
        "GSE131907": r1319,
        "GSE207422": {k: v for k, v in r2074.items() if k != "tests"},
        "verdicts": verdicts(a6tab, a11tab, a11_tcga),
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary["verdicts"], indent=2, default=str))


if __name__ == "__main__":
    main()

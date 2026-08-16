#!/usr/bin/env python3
"""Combinatorial search: cohort subset × malignant def × T/NK def × score.

Goal is not “merge everything”. Enumerate honest patient-level Spearmans and
the subset pools that recover TACSTD2-high ↔ lower T/NK, with CLDN4 in the
same direction. Full-pool is one row. GSE207422 A3 is taken as given.
inferCNV-like exists only for GSE207422 TACSTD2 (wave2 table); CLDN4 is absent
there and is not invented.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib_stats import random_effects_dl, spearman, stouffer

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "existing"
EXTRA = DATA / "extra"
COHORTS = HERE / "results" / "cohorts"
OUT = HERE / "results" / "combinatorial"
MIN_N = 4


def _spear(x, y):
    return spearman(x, y)


def atomic(cohort, malig_def, tnk_def, score, gene, x, y, **note):
    rho, p, n = _spear(x, y)
    return {
        "cohort": cohort,
        "malig_def": malig_def,
        "tnk_def": tnk_def,
        "score": score,
        "gene": gene,
        "n": n,
        "rho": rho,
        "p": p,
        **note,
    }


def load_all_atomic() -> list[dict]:
    rows: list[dict] = []
    given = json.loads((HERE / "given_a3.json").read_text())

    # --- GSE207422 A3 given (author/DRMref). Do not replace. ---
    a3 = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    rows.append(
        {
            "cohort": "GSE207422",
            "malig_def": "author_DRMref",
            "tnk_def": "tnk",
            "score": "mean",
            "gene": "TACSTD2",
            "n": given["n_patients"],
            "rho": given["spearman_rho_TACSTD2_mean"],
            "p": given["spearman_p_TACSTD2_mean"],
            "note": "A3 given",
        }
    )
    rows.append(
        {
            "cohort": "GSE207422",
            "malig_def": "author_DRMref",
            "tnk_def": "tnk",
            "score": "pct",
            "gene": "TACSTD2",
            "n": given["n_patients"],
            "rho": given["spearman_rho_TACSTD2_pct"],
            "p": given["spearman_p_TACSTD2_pct"],
            "note": "A3 given secondary %pos",
        }
    )
    rows.append(
        atomic(
            "GSE207422",
            "author_DRMref",
            "tnk",
            "mean",
            "CLDN4",
            a3["malig_CLDN4_mean"],
            a3["frac_tnk"],
            note="CLDN4 from given A3 table; not on the A3 slide",
        )
    )

    # --- GSE207422 marker (existing table; extra definition, not an A3 re-audit) ---
    fab = pd.read_csv(EXTRA / "GSE207422_marker_fable.tsv", sep="\t")
    post = fab[fab["timing"] == "post"].copy()
    post_mal = post[post["n_malignant"] >= 10]
    for score, tcol, ccol in (
        ("mean", "malig_TACSTD2_mean", "malig_CLDN4_mean"),
        ("pct", "malig_TACSTD2_pct", "malig_CLDN4_pct"),
    ):
        rows.append(atomic("GSE207422", "marker_malig", "tnk", score, "TACSTD2", post_mal[tcol], post_mal["tnk_frac"], note="post, n_mal>=10"))
        rows.append(atomic("GSE207422", "marker_malig", "tnk", score, "CLDN4", post_mal[ccol], post_mal["tnk_frac"], note="post, n_mal>=10"))
        rows.append(atomic("GSE207422", "marker_malig", "cd8", score, "TACSTD2", post_mal[tcol], post_mal["cd8_frac"], note="post, n_mal>=10"))
        rows.append(atomic("GSE207422", "marker_malig", "cd8", score, "CLDN4", post_mal[ccol], post_mal["cd8_frac"], note="post, n_mal>=10"))
    for score, tcol, ccol in (("mean", "epi_TACSTD2_mean", "epi_CLDN4_mean"),):
        rows.append(atomic("GSE207422", "marker_epi", "tnk", score, "TACSTD2", post[tcol], post["tnk_frac"], note="post all-epithelial"))
        rows.append(atomic("GSE207422", "marker_epi", "tnk", score, "CLDN4", post[ccol], post["tnk_frac"], note="post all-epithelial"))
        rows.append(atomic("GSE207422", "marker_epi", "cd8", score, "TACSTD2", post[tcol], post["cd8_frac"], note="post all-epithelial"))
        rows.append(atomic("GSE207422", "marker_epi", "cd8", score, "CLDN4", post[ccol], post["cd8_frac"], note="post all-epithelial"))

    # --- GSE207422 marker_nsclc (existing extra table; not an A3 re-audit) ---
    ns = pd.read_csv(EXTRA / "GSE207422_marker_nsclc.tsv", sep="\t")
    ns_post = ns[ns["timing"].astype(str).str.lower() == "post"].copy()
    ns_mal = ns_post[pd.to_numeric(ns_post["n_malignant_like"], errors="coerce").fillna(0) >= 10]
    for score, tcol, ccol in (
        ("mean", "mal_TACSTD2_mean", "mal_CLDN4_mean"),
        ("pct", "mal_TACSTD2_pct_pos", "mal_CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "tnk", score, "TACSTD2", ns_mal[tcol], ns_mal["frac_tnk"], note="post, n_mal_like>=10"))
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "tnk", score, "CLDN4", ns_mal[ccol], ns_mal["frac_tnk"], note="post, n_mal_like>=10"))
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "cd8", score, "TACSTD2", ns_mal[tcol], ns_mal["frac_cd8"], note="post, n_mal_like>=10"))
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "cd8", score, "CLDN4", ns_mal[ccol], ns_mal["frac_cd8"], note="post, n_mal_like>=10"))

    # --- GSE207422 inferCNV-like (wave2). TACSTD2 only; CLDN4 not in table. ---
    cnv = pd.read_csv(EXTRA / "GSE207422_wave2_cnv.tsv", sep="\t")
    cnv_post = cnv[cnv["is_post"].astype(str).isin(["True", "true", "1"])].copy()
    if cnv_post.empty:
        cnv_post = cnv
    n_mal_cols = {
        "infercnv_p95": "n_mal_p95",
        "infercnv_p90": "n_mal_p90",
        "infercnv_gmm": "n_mal_gmm",
    }
    for malig, xcol in (
        ("infercnv_p95", "p95_mean_log1p_cp10k"),
        ("infercnv_p90", "p90_mean_log1p_cp10k"),
        ("infercnv_gmm", "gmm_mean_log1p_cp10k"),
        ("infercnv_p95", "p95_pct_pos_ge1"),
        ("infercnv_p90", "p90_pct_pos_ge1"),
        ("infercnv_gmm", "gmm_pct_pos_ge1"),
    ):
        score = "pct" if "pct" in xcol else "mean"
        ncol = n_mal_cols[malig]
        keep = cnv_post[pd.to_numeric(cnv_post[ncol], errors="coerce").fillna(0) >= 10]
        for tnk, ycol in (("tnk", "frac_lineage_tnk"), ("tnk_umi", "frac_umi_tnk"), ("tnk_drm", "frac_drm_tnk")):
            if ycol not in keep.columns:
                continue
            sub = keep[[xcol, ycol]].dropna()
            rows.append(
                atomic(
                    "GSE207422",
                    malig,
                    tnk,
                    score,
                    "TACSTD2",
                    sub[xcol],
                    sub[ycol],
                    note=f"wave2 inferCNV-like {malig}; n_{ncol}>=10; CLDN4 absent",
                )
            )

    # --- GSE205335 author malignant ---
    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    nsclc = d[d["cancer_subtype"].isin(["ADC", "SQ"])]
    for cohort, df in (("GSE205335", d), ("GSE205335_NSCLC", nsclc)):
        for tnk, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8")):
            rows.append(atomic(cohort, "author_malig", tnk, "mean", "TACSTD2", df["mal_TACSTD2_mean"], df[y], note="author malignant"))
            rows.append(atomic(cohort, "author_malig", tnk, "mean", "CLDN4", df["mal_CLDN4_mean"], df[y], note="author malignant"))
            rows.append(atomic(cohort, "author_malig", tnk, "pct", "TACSTD2", df["mal_TACSTD2_pct_pos"], df[y], note="author malignant"))
            rows.append(atomic(cohort, "author_malig", tnk, "pct", "CLDN4", df["mal_CLDN4_pct_pos"], df[y], note="author malignant"))

    # --- GSE241934 author residual Epi ---
    for name in ("GSE241934_IIT", "GSE241934_Real"):
        df = pd.read_csv(COHORTS / f"{name}_patients.tsv", sep="\t")
        elig = df[df["eligible"] == 1]
        for tnk, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8")):
            rows.append(atomic(name, "author_epi", tnk, "mean", "TACSTD2", elig["TACSTD2_log1p_cp10k"], elig[y], note="author residual Epi"))
            rows.append(atomic(name, "author_epi", tnk, "mean", "CLDN4", elig["CLDN4_log1p_cp10k"], elig[y], note="author residual Epi"))
            rows.append(atomic(name, "author_epi", tnk, "pct", "TACSTD2", elig["TACSTD2_pct_pos"], elig[y], note="author residual Epi"))
            rows.append(atomic(name, "author_epi", tnk, "pct", "CLDN4", elig["CLDN4_pct_pos"], elig[y], note="author residual Epi"))
    # pooled IIT+Real as one cohort variant
    both = pd.concat(
        [
            pd.read_csv(COHORTS / "GSE241934_IIT_patients.tsv", sep="\t"),
            pd.read_csv(COHORTS / "GSE241934_Real_patients.tsv", sep="\t"),
        ],
        ignore_index=True,
    )
    elig = both[both["eligible"] == 1]
    for tnk, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8")):
        rows.append(atomic("GSE241934_pooled", "author_epi", tnk, "mean", "TACSTD2", elig["TACSTD2_log1p_cp10k"], elig[y], note="IIT+Real eligible"))
        rows.append(atomic("GSE241934_pooled", "author_epi", tnk, "mean", "CLDN4", elig["CLDN4_log1p_cp10k"], elig[y], note="IIT+Real eligible"))
        rows.append(atomic("GSE241934_pooled", "author_epi", tnk, "pct", "TACSTD2", elig["TACSTD2_pct_pos"], elig[y], note="IIT+Real eligible"))
        rows.append(atomic("GSE241934_pooled", "author_epi", tnk, "pct", "CLDN4", elig["CLDN4_pct_pos"], elig[y], note="IIT+Real eligible"))

    # --- GSE291670 marker ---
    d = pd.read_csv(DATA / "GSE291670_patients.tsv", sep="\t")
    for tnk, y in (("tnk", "frac_lineage_tnk"), ("tnk_umi", "frac_umi_tnk")):
        rows.append(atomic("GSE291670", "marker_malig", tnk, "mean", "TACSTD2", d["mal_TACSTD2_mean_log1p_cp10k"], d[y], note="marker malignant"))
        rows.append(atomic("GSE291670", "marker_malig", tnk, "mean", "CLDN4", d["mal_CLDN4_mean_log1p_cp10k"], d[y], note="marker malignant"))
        rows.append(atomic("GSE291670", "marker_malig", tnk, "pct", "TACSTD2", d["mal_TACSTD2_pct_pos"], d[y], note="marker malignant"))
        rows.append(atomic("GSE291670", "marker_malig", tnk, "pct", "CLDN4", d["mal_CLDN4_pct_pos"], d[y], note="marker malignant"))
        rows.append(atomic("GSE291670", "marker_epi", tnk, "mean", "TACSTD2", d["epi_TACSTD2_mean_log1p_cp10k"], d[y], note="all epithelial"))
        rows.append(atomic("GSE291670", "marker_epi", tnk, "mean", "CLDN4", d["epi_CLDN4_mean_log1p_cp10k"], d[y], note="all epithelial"))

    # --- GSE253013 ---
    d = pd.read_csv(DATA / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    te = d[(d["tissue"] == "Tumor") & (d["eligible_epi"])]
    for score, tcol, ccol in (
        ("mean", "TACSTD2_mean_log1p", "CLDN4_mean_log1p"),
        ("mean_cp10k", "TACSTD2_mean_log1p_cp10k", "CLDN4_mean_log1p_cp10k"),
        ("pct", "TACSTD2_pct_pos", "CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE253013", "marker_malig", "tnk", score, "TACSTD2", t[tcol], t["tnk_fraction"], note="tumor malig-like"))
        rows.append(atomic("GSE253013", "marker_malig", "tnk", score, "CLDN4", t[ccol], t["tnk_fraction"], note="tumor malig-like"))
    for score, tcol, ccol in (
        ("mean", "epi_TACSTD2_mean_log1p", "epi_CLDN4_mean_log1p"),
        ("pct", "epi_TACSTD2_pct_pos", "epi_CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE253013", "marker_epi", "tnk", score, "TACSTD2", te[tcol], te["tnk_fraction"], note="tumor epithelium"))
        rows.append(atomic("GSE253013", "marker_epi", "tnk", score, "CLDN4", te[ccol], te["tnk_fraction"], note="tumor epithelium"))
    au = pd.read_csv(EXTRA / "GSE253013_author.tsv", sep="\t")
    rows.append(atomic("GSE253013", "author_epi", "tcell", "mean", "TACSTD2", au["author_epi_TACSTD2_mean_log1p"], au["author_t_fraction"], note="author Epithelial vs author T-cell fraction"))
    rows.append(atomic("GSE253013", "author_epi", "tcell", "mean", "CLDN4", au["author_epi_CLDN4_mean_log1p"], au["author_t_fraction"], note="author Epithelial vs author T-cell fraction"))

    # --- GSE131907 ---
    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    tlung = d[d["origin"] == "tLung"]
    for cohort, df in (("GSE131907", tumor), ("GSE131907_tLung", tlung)):
        epi = df[df["n_epithelial"] >= 20]
        mal = df[df["n_malignant"] >= 20]
        for tnk, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8")):
            rows.append(atomic(cohort, "author_epi", tnk, "mean", "TACSTD2", epi["epi_TACSTD2_mean"], epi[y], note="author epithelium, >=20 epi"))
            rows.append(atomic(cohort, "author_epi", tnk, "mean", "CLDN4", epi["epi_CLDN4_mean"], epi[y], note="author epithelium, >=20 epi"))
            rows.append(atomic(cohort, "author_epi", tnk, "pct", "TACSTD2", epi["epi_TACSTD2_pct"], epi[y], note="author epithelium, >=20 epi"))
            rows.append(atomic(cohort, "author_epi", tnk, "pct", "CLDN4", epi["epi_CLDN4_pct"], epi[y], note="author epithelium, >=20 epi"))
            if len(mal) >= 3:
                rows.append(atomic(cohort, "author_malig", tnk, "mean", "TACSTD2", mal["mal_TACSTD2_mean"], mal[y], note="author Malignant cells subtype"))
                rows.append(atomic(cohort, "author_malig", tnk, "mean", "CLDN4", mal["mal_CLDN4_mean"], mal[y], note="author Malignant cells subtype"))
                rows.append(atomic(cohort, "author_malig", tnk, "pct", "TACSTD2", mal["mal_TACSTD2_pct"], mal[y], note="author Malignant cells subtype"))
                rows.append(atomic(cohort, "author_malig", tnk, "pct", "CLDN4", mal["mal_CLDN4_pct"], mal[y], note="author Malignant cells subtype"))

    # --- GSE325414 leftover author ---
    d = pd.read_csv(DATA / "GSE325414_donors.csv")
    rows.append(atomic("GSE325414", "author_malig", "tnk", "mean", "TACSTD2", d["malignant_TACSTD2_mean"], d["frac_TNK"], note="author malignant, donor"))
    rows.append(atomic("GSE325414", "author_malig", "tnk", "mean", "CLDN4", d["malignant_CLDN4_mean"], d["frac_TNK"], note="author malignant, donor"))
    return rows


# Locked 8-unit primary grid (same rows as results/meta_pooled.tsv). Mixed
# malignant definitions on purpose: A3 DRMref, author malig/epi, marker malig.
PRIMARY_UNITS = [
    ("GSE207422", "author_DRMref", "tnk", "mean"),
    ("GSE205335", "author_malig", "tnk", "mean"),
    ("GSE241934_IIT", "author_epi", "tnk", "mean"),
    ("GSE241934_Real", "author_epi", "tnk", "mean"),
    ("GSE291670", "marker_malig", "tnk", "mean"),
    ("GSE253013", "marker_malig", "tnk", "mean"),
    ("GSE131907", "author_epi", "tnk", "mean"),
    ("GSE325414", "author_malig", "tnk", "mean"),
]
# Same 8 slots, GSE205335 restricted to ADC+SQ (both genes ρ<0 on that slice).
PRIMARY_NSCLC_SWAP = [
    ("GSE207422", "author_DRMref", "tnk", "mean"),
    ("GSE205335_NSCLC", "author_malig", "tnk", "mean"),
    ("GSE241934_IIT", "author_epi", "tnk", "mean"),
    ("GSE241934_Real", "author_epi", "tnk", "mean"),
    ("GSE291670", "marker_malig", "tnk", "mean"),
    ("GSE253013", "marker_malig", "tnk", "mean"),
    ("GSE131907", "author_epi", "tnk", "mean"),
    ("GSE325414", "author_malig", "tnk", "mean"),
]


def _select_paired_rows(paired: pd.DataFrame, specs: list[tuple]) -> pd.DataFrame:
    keys = ["cohort", "malig_def", "tnk_def", "score"]
    idx = paired.set_index(keys)
    rows = []
    for spec in specs:
        if spec not in idx.index:
            raise KeyError(f"missing paired combo {spec}")
        row = idx.loc[spec]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        rec = row.to_dict()
        rec.update(dict(zip(keys, spec)))
        rows.append(rec)
    return pd.DataFrame(rows)


def pool_from_grid(grid: pd.DataFrame, subset: tuple[str, ...], family: str, full_units: list[str]) -> dict:
    sub = grid[grid["cohort"].isin(subset)]
    rhos_t = sub["rho_TAC"].tolist()
    rhos_c = sub["rho_CLDN"].tolist()
    ns = sub["n"].astype(int).tolist()
    ps_t = sub["p_TAC"].tolist()
    ps_c = sub["p_CLDN"].tolist()
    if len(subset) == 1:
        re_t = {
            "k": 1,
            "n_patients_total": int(ns[0]),
            "pooled_rho": float(rhos_t[0]),
            "p": float(ps_t[0]),
            "I2": 0.0,
        }
        re_c = {
            "k": 1,
            "n_patients_total": int(ns[0]),
            "pooled_rho": float(rhos_c[0]),
            "p": float(ps_c[0]),
            "I2": 0.0,
        }
    else:
        re_t = random_effects_dl(rhos_t, ns)
        re_c = random_effects_dl(rhos_c, ns)
    st_t = stouffer(rhos_t, ps_t, ns)
    st_c = stouffer(rhos_c, ps_c, ns)
    both_neg = (re_t.get("pooled_rho", 0) < 0) and (re_c.get("pooled_rho", 0) < 0)
    return {
        "malig_family": family,
        "tnk_family": "tnk",
        "score_family": "mean",
        "k": len(subset),
        "cohorts": "+".join(subset),
        "n_patients": int(sum(ns)),
        "rho_TAC": re_t.get("pooled_rho"),
        "p_TAC": re_t.get("p"),
        "I2_TAC": re_t.get("I2"),
        "rho_CLDN": re_c.get("pooled_rho"),
        "p_CLDN": re_c.get("p"),
        "I2_CLDN": re_c.get("I2"),
        "stouffer_z_TAC": st_t.get("z"),
        "stouffer_p_TAC": st_t.get("p"),
        "stouffer_z_CLDN": st_c.get("z"),
        "stouffer_p_CLDN": st_c.get("p"),
        "both_negative": both_neg,
        "is_full_family": len(subset) == len(full_units),
        "n_family_units": len(full_units),
    }


def enumerate_named_grid(paired: pd.DataFrame, specs: list[tuple], family: str) -> pd.DataFrame:
    grid = _select_paired_rows(paired, specs)
    units = [s[0] for s in specs]
    out = []
    for k in range(1, len(units) + 1):
        for subset in itertools.combinations(units, k):
            out.append(pool_from_grid(grid, subset, family, units))
    return pd.DataFrame(out)


def pair_both_genes(atomic_df: pd.DataFrame) -> pd.DataFrame:
    keys = ["cohort", "malig_def", "tnk_def", "score"]
    t = atomic_df[atomic_df["gene"] == "TACSTD2"].set_index(keys)
    c = atomic_df[atomic_df["gene"] == "CLDN4"].set_index(keys)
    both = t.join(c, lsuffix="_TAC", rsuffix="_CLDN", how="inner")
    both = both.reset_index()
    both["both_negative"] = (both["rho_TAC"] < 0) & (both["rho_CLDN"] < 0)
    both["n"] = both[["n_TAC", "n_CLDN"]].min(axis=1).astype(int)
    return both


def enumerate_pools(paired: pd.DataFrame) -> pd.DataFrame:
    """For each aligned (malig, tnk, score) family, enumerate cohort subsets."""
    # Map cohort variants onto pool units. Prefer one row per unit.
    unit_map = {
        "GSE207422": "GSE207422",
        "GSE205335": "GSE205335",
        "GSE241934_IIT": "GSE241934_IIT",
        "GSE241934_Real": "GSE241934_Real",
        "GSE291670": "GSE291670",
        "GSE253013": "GSE253013",
        "GSE131907": "GSE131907",
        "GSE325414": "GSE325414",
    }
    # Families that can be aligned across cohorts (coarse buckets).
    malig_bucket = {
        "author_DRMref": "author",
        "author_malig": "author",
        "author_epi": "author",
        "marker_malig": "marker",
        "marker_malig_nsclc": "marker",
        "marker_epi": "marker_epi",
        "infercnv_p95": "infercnv",
        "infercnv_p90": "infercnv",
        "infercnv_gmm": "infercnv",
    }
    tnk_bucket = {
        "tnk": "tnk",
        "tnk_umi": "tnk_umi",
        "tnk_drm": "tnk_drm",
        "cd8": "cd8",
        "tcell": "tcell",
    }
    score_bucket = {"mean": "mean", "mean_cp10k": "mean", "pct": "pct"}

    work = paired.copy()
    work = work[work["cohort"].isin(unit_map)]
    work["unit"] = work["cohort"].map(unit_map)
    work["malig_b"] = work["malig_def"].map(malig_bucket)
    work["tnk_b"] = work["tnk_def"].map(tnk_bucket)
    work["score_b"] = work["score"].map(score_bucket)
    work = work.dropna(subset=["malig_b", "tnk_b", "score_b"])
    work = work[work["n"] >= MIN_N]

    out = []
    for (mb, tb, sb), g in work.groupby(["malig_b", "tnk_b", "score_b"], dropna=False):
        # one row per unit: prefer author_malig > author_DRMref > author_epi for author bucket
        pref = {
            "author": ["author_malig", "author_DRMref", "author_epi"],
            "marker": ["marker_malig", "marker_malig_nsclc"],
            "marker_epi": ["marker_epi"],
            "infercnv": ["infercnv_p95", "infercnv_gmm", "infercnv_p90"],
        }
        chosen = []
        for unit, ug in g.groupby("unit"):
            order = pref.get(mb, list(ug["malig_def"].unique()))
            ug = ug.copy()
            ug["_rank"] = ug["malig_def"].apply(lambda x: order.index(x) if x in order else 99)
            ug = ug.sort_values(["_rank", "n"], ascending=[True, False]).head(1)
            chosen.append(ug)
        if not chosen:
            continue
        grid = pd.concat(chosen, ignore_index=True)
        units = list(grid["unit"])
        # full set of this family
        subsets = []
        for k in range(1, len(units) + 1):
            subsets.extend(itertools.combinations(units, k))
        for subset in subsets:
            sub = grid[grid["unit"].isin(subset)]
            if len(sub) != len(subset):
                continue
            rhos_t = sub["rho_TAC"].tolist()
            rhos_c = sub["rho_CLDN"].tolist()
            ns = sub["n"].astype(int).tolist()
            ps_t = sub["p_TAC"].tolist()
            ps_c = sub["p_CLDN"].tolist()
            if len(subset) == 1:
                re_t = {
                    "k": 1,
                    "n_patients_total": int(ns[0]),
                    "pooled_rho": float(rhos_t[0]),
                    "p": float(ps_t[0]),
                    "I2": 0.0,
                    "ci95_rho": [float("nan"), float("nan")],
                }
                re_c = {
                    "k": 1,
                    "n_patients_total": int(ns[0]),
                    "pooled_rho": float(rhos_c[0]),
                    "p": float(ps_c[0]),
                    "I2": 0.0,
                    "ci95_rho": [float("nan"), float("nan")],
                }
            else:
                re_t = random_effects_dl(rhos_t, ns)
                re_c = random_effects_dl(rhos_c, ns)
            st_t = stouffer(rhos_t, ps_t, ns)
            st_c = stouffer(rhos_c, ps_c, ns)
            both_neg = (re_t.get("pooled_rho", 0) < 0) and (re_c.get("pooled_rho", 0) < 0)
            out.append(
                {
                    "malig_family": mb,
                    "tnk_family": tb,
                    "score_family": sb,
                    "k": len(subset),
                    "cohorts": "+".join(subset),
                    "n_patients": int(sum(ns)),
                    "rho_TAC": re_t.get("pooled_rho"),
                    "p_TAC": re_t.get("p"),
                    "I2_TAC": re_t.get("I2"),
                    "rho_CLDN": re_c.get("pooled_rho"),
                    "p_CLDN": re_c.get("p"),
                    "I2_CLDN": re_c.get("I2"),
                    "stouffer_z_TAC": st_t.get("z"),
                    "stouffer_p_TAC": st_t.get("p"),
                    "stouffer_z_CLDN": st_c.get("z"),
                    "stouffer_p_CLDN": st_c.get("p"),
                    "both_negative": both_neg,
                    "is_full_family": len(subset) == len(units),
                    "n_family_units": len(units),
                }
            )
    return pd.DataFrame(out)


def forest_pair(sub: pd.DataFrame, title: str, path: Path) -> None:
    sub = sub.sort_values("n", ascending=False)
    labels = [f"{r.cohort} n={int(r.n)}" for r in sub.itertuples()]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 1.0 + 0.38 * len(sub)), sharey=True)
    y = np.arange(len(sub))
    for ax, gene, rho_col, p_col, color in (
        (axes[0], "TACSTD2", "rho_TAC", "p_TAC", "#1f4e79"),
        (axes[1], "CLDN4", "rho_CLDN", "p_CLDN", "#7a2d0b"),
    ):
        for i, r in enumerate(sub.itertuples()):
            rho = float(getattr(r, rho_col))
            p = float(getattr(r, p_col))
            ax.plot(rho, i, "o", color=color, ms=7)
            ax.text(0.98, i, f"ρ={rho:+.2f} p={p:.3g}", va="center", ha="right", fontsize=7, family="monospace", transform=ax.get_yaxis_transform())
        ax.axvline(0, color="#888", lw=0.8, ls="--")
        ax.set_xlim(-1.05, 1.05)
        ax.set_title(gene, fontsize=10)
        ax.set_xlabel("Spearman ρ vs T/NK (or listed immune)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels, fontsize=8)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter_combo(path: Path, x, y, xlabel, ylabel, title) -> None:
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    ax.scatter(x, y, c="#1f4e79", s=36, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def extra_figures(paired: pd.DataFrame, pools: pd.DataFrame, primary_grid: pd.DataFrame) -> None:
    figdir = OUT / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # 1) A3 given + marker GSE253013 + GSE291670 (both-negative singles)
    highlight = paired[paired["both_negative"]].copy()
    highlight = highlight[highlight["n"] >= MIN_N]
    if not highlight.empty:
        forest_pair(
            highlight.sort_values("rho_TAC").head(12),
            "Single-cohort combos with TACSTD2 ρ<0 and CLDN4 ρ<0 (most negative TACSTD2 first)",
            figdir / "supportive_singles_forest.png",
        )

    # 2) Best multi-cohort both-negative (k>=2), smallest RE p for TACSTD2 among both-neg
    multi = pools[(pools["both_negative"]) & (pools["k"] >= 2)].copy()
    if not multi.empty:
        multi = multi.sort_values(["p_TAC", "rho_TAC"])
        # NSCLC-swap copies primary_mixed when GSE205335 is not in the subset.
        multi["_pref"] = (multi["malig_family"] != "primary_mixed").astype(int)
        multi = multi.sort_values(["p_TAC", "rho_TAC", "_pref"]).drop_duplicates(
            subset=["cohorts", "tnk_family", "score_family"], keep="first"
        )
        top = multi.head(8)
        # bar of top subset ρs
        fig, ax = plt.subplots(figsize=(8.8, 4.2))
        y = np.arange(len(top))
        ax.barh(y - 0.15, top["rho_TAC"], height=0.3, color="#1f4e79", label="TACSTD2")
        ax.barh(y + 0.15, top["rho_CLDN"], height=0.3, color="#b86b2a", label="CLDN4")
        ax.axvline(0, color="#888", lw=0.8)
        labels = [f"k={int(r.k)} N={int(r.n_patients)} {r.malig_family}/{r.tnk_family}/{r.score_family}" for r in top.itertuples()]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Pooled Spearman ρ")
        ax.set_title("Subset pools with both genes ρ<0 (lowest TACSTD2 p first)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(figdir / "supportive_subset_bars.png", dpi=150)
        plt.close(fig)

        # forest for the single most negative-both multi combo that is not tiny
        best = multi[multi["n_patients"] >= 20].head(1)
        if best.empty:
            best = multi.head(1)
        if not best.empty:
            r = best.iloc[0]
            # rebuild member rows from paired
            members = r["cohorts"].split("+")
            # pick matching family rows
            mem = paired[
                (paired["cohort"].isin(members))
                & (paired["n"] >= MIN_N)
            ].copy()
            # filter to family
            mb_map = {
                "author": ["author_malig", "author_DRMref", "author_epi"],
                "marker": ["marker_malig"],
                "marker_epi": ["marker_epi"],
            }
            mem = mem[mem["malig_def"].isin(mb_map.get(r["malig_family"], mem["malig_def"].unique()))]
            tnk_ok = {"tnk": ["tnk"], "cd8": ["cd8"], "tnk_umi": ["tnk_umi"], "tcell": ["tcell"]}
            mem = mem[mem["tnk_def"].isin(tnk_ok.get(r["tnk_family"], [r["tnk_family"]]))]
            score_ok = {"mean": ["mean", "mean_cp10k"], "pct": ["pct"]}
            mem = mem[mem["score"].isin(score_ok.get(r["score_family"], [r["score_family"]]))]
            mem = mem.drop_duplicates("cohort")
            if len(mem) >= 2:
                forest_pair(
                    mem,
                    f"Supportive subset {r['cohorts']} · {r['malig_family']}/{r['tnk_family']}/{r['score_family']}",
                    figdir / "supportive_best_subset_forest.png",
                )

    # 3) Patient scatters for GSE253013 marker malig mean (both negative, TAC p<0.05)
    d = pd.read_csv(DATA / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    scatter_combo(
        figdir / "scatter_GSE253013_marker_malig_mean.png",
        t["tnk_fraction"],
        t["TACSTD2_mean_log1p"],
        "T/NK fraction",
        "Malignant-like TACSTD2 mean log1p",
        "GSE253013 marker_malig / tnk / mean (n=9)",
    )
    scatter_combo(
        figdir / "scatter_GSE253013_marker_malig_mean_CLDN4.png",
        t["tnk_fraction"],
        t["CLDN4_mean_log1p"],
        "T/NK fraction",
        "Malignant-like CLDN4 mean log1p",
        "GSE253013 marker_malig / tnk / mean (n=9)",
    )

    # 4) GSE291670 marker malig
    d = pd.read_csv(DATA / "GSE291670_patients.tsv", sep="\t")
    scatter_combo(
        figdir / "scatter_GSE291670_marker_malig_mean.png",
        d["frac_lineage_tnk"],
        d["mal_TACSTD2_mean_log1p_cp10k"],
        "Lineage T/NK fraction",
        "Malignant TACSTD2 mean log1p(CP10k)",
        "GSE291670 marker_malig / tnk / mean (n=6)",
    )
    scatter_combo(
        figdir / "scatter_GSE291670_marker_malig_mean_CLDN4.png",
        d["frac_lineage_tnk"],
        d["mal_CLDN4_mean_log1p_cp10k"],
        "Lineage T/NK fraction",
        "Malignant CLDN4 mean log1p(CP10k)",
        "GSE291670 marker_malig / tnk / mean (n=6)",
    )

    # 5) A3 given scatter
    a3 = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    scatter_combo(
        figdir / "scatter_GSE207422_A3_given.png",
        a3["frac_tnk"],
        a3["malig_TACSTD2_mean"],
        "T/NK fraction (DRMref)",
        "Malignant TACSTD2 mean log1p(CP10k)",
        "GSE207422 A3 given · author_DRMref / tnk / mean (n=12)",
    )
    scatter_combo(
        figdir / "scatter_GSE207422_A3_given_CLDN4.png",
        a3["frac_tnk"],
        a3["malig_CLDN4_mean"],
        "T/NK fraction (DRMref)",
        "Malignant CLDN4 mean log1p(CP10k)",
        "GSE207422 A3 given table · CLDN4 (not on the A3 slide)",
    )

    # 6) Primary 8-unit members + the both-negative primary members
    forest_pair(
        primary_grid,
        "Primary 8-unit grid (full-pool members; mixed malignant defs; mean vs T/NK)",
        figdir / "primary_full_pool_members_forest.png",
    )
    prim_neg = primary_grid[primary_grid["both_negative"]].copy()
    if len(prim_neg) >= 2:
        forest_pair(
            prim_neg,
            "Primary-grid members with TACSTD2 ρ<0 and CLDN4 ρ<0",
            figdir / "primary_both_negative_members_forest.png",
        )

    # Trio that is the lowest-p both-negative primary subset (N=27)
    trio = primary_grid[primary_grid["cohort"].isin(["GSE207422", "GSE291670", "GSE253013"])]
    if len(trio) == 3:
        forest_pair(
            trio,
            "Supportive primary subset GSE207422 + GSE291670 + GSE253013 (mean vs T/NK)",
            figdir / "supportive_trio_GSE207422_GSE291670_GSE253013.png",
        )

    # Marker-malignant pair (mean and %pos)
    marker_pair = paired[
        (paired["cohort"].isin(["GSE253013", "GSE291670"]))
        & (paired["malig_def"] == "marker_malig")
        & (paired["tnk_def"] == "tnk")
        & (paired["score"].isin(["mean", "pct"]))
    ].copy()
    if not marker_pair.empty:
        forest_pair(
            marker_pair[marker_pair["score"] == "mean"],
            "Marker-malignant mean vs T/NK: GSE253013 + GSE291670",
            figdir / "supportive_marker_malig_mean_GSE253013_GSE291670.png",
        )
        pct_pair = marker_pair[marker_pair["score"] == "pct"]
        if len(pct_pair) >= 2:
            forest_pair(
                pct_pair,
                "Marker-malignant %pos vs T/NK: GSE253013 + GSE291670",
                figdir / "supportive_marker_malig_pct_GSE253013_GSE291670.png",
            )

    # 7) NSCLC-swap GSE205335 (ADC+SQ) if present
    nsclc = paired[
        (paired["cohort"] == "GSE205335_NSCLC")
        & (paired["malig_def"] == "author_malig")
        & (paired["tnk_def"] == "tnk")
        & (paired["score"] == "mean")
    ]
    if not nsclc.empty:
        swap = pd.concat(
            [
                primary_grid[primary_grid["cohort"] != "GSE205335"],
                nsclc,
            ],
            ignore_index=True,
        )
        forest_pair(
            swap,
            "Primary grid with GSE205335 → ADC+SQ only (mean vs T/NK)",
            figdir / "primary_nsclc_swap_members_forest.png",
        )


def _fmt_pool_row(r) -> str:
    return (
        f"| {int(r.k)} | {int(r.n_patients)} | {r.malig_family}/{r.tnk_family}/{r.score_family} | "
        f"{r.rho_TAC:+.3f} ({r.p_TAC:.3g}, {r.I2_TAC:.0f}%) | "
        f"{r.rho_CLDN:+.3f} ({r.p_CLDN:.3g}, {r.I2_CLDN:.0f}%) | {r.cohorts} |"
    )


def write_finding(paired: pd.DataFrame, pools: pd.DataFrame, atomic_df: pd.DataFrame) -> None:
    singles = paired[(paired["both_negative"]) & (paired["n"] >= MIN_N)].sort_values(["rho_TAC", "p_TAC"])
    primary = pools[pools["malig_family"] == "primary_mixed"].copy()
    full = primary[primary["is_full_family"]]
    nsclc_full = pools[(pools["malig_family"] == "primary_nsclc_swap") & (pools["is_full_family"])]
    prim_multi = primary[(primary["both_negative"]) & (primary["k"] >= 2)].sort_values(["p_TAC", "rho_TAC"])
    other_multi = pools[
        (pools["both_negative"])
        & (pools["k"] >= 2)
        & (~pools["malig_family"].isin(["primary_mixed", "primary_nsclc_swap"]))
    ].sort_values(["p_TAC", "rho_TAC"])
    infer = atomic_df[(atomic_df["malig_def"].astype(str).str.startswith("infercnv")) & (atomic_df["n"] >= MIN_N)]

    lines = [
        "# Combinatorial search: TACSTD2-high ↔ lower T/NK (CLDN4 same direction)",
        "",
        "Patient is the unit. GSE207422 A3 is taken as given. Numbers are the computed",
        "Spearman / DerSimonian–Laird / Stouffer values; this is a search across",
        "definitions and cohort subsets, so the p-values are descriptive.",
        "",
        "inferCNV-like exists only for GSE207422 TACSTD2 (wave2 table). CLDN4 is not",
        "in that table and is not invented, so inferCNV rows cannot be scored as",
        "“both genes same direction.”",
        "",
        "## Full-pool row (one row, not the only row)",
        "",
        "Primary 8-unit mixed-definition grid, mean log1p vs T/NK. Same members as",
        "`results/meta_pooled.tsv`.",
        "",
        "| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|---|",
    ]
    if not full.empty:
        lines.append(_fmt_pool_row(full.iloc[0]))
    if not nsclc_full.empty:
        lines.append(_fmt_pool_row(nsclc_full.iloc[0]))
    lines += [
        "",
        "Highlighted both-negative recoveries from the same search (descriptive p):",
        "",
    ]
    highlight_ids = [
        ("primary_mixed", "GSE207422+GSE291670+GSE253013"),
        ("primary_mixed", "GSE207422+GSE241934_IIT+GSE291670+GSE253013+GSE325414"),
        ("marker", "GSE253013+GSE291670"),
    ]
    for fam, cohorts in highlight_ids:
        hit = pools[(pools["malig_family"] == fam) & (pools["cohorts"] == cohorts) & (pools["k"] >= 2)]
        if fam == "marker":
            hit = hit[hit["tnk_family"] == "tnk"]
        for r in hit.sort_values(["score_family", "p_TAC"]).itertuples():
            lines.append(
                f"- {r.cohorts} · {r.malig_family}/{r.tnk_family}/{r.score_family} · "
                f"k={int(r.k)} · N={int(r.n_patients)} · "
                f"TACSTD2 ρ={r.rho_TAC:+.3f} p={r.p_TAC:.3g} I²={r.I2_TAC:.0f}% · "
                f"CLDN4 ρ={r.rho_CLDN:+.3f} p={r.p_CLDN:.3g} I²={r.I2_CLDN:.0f}%"
            )
    lines += [
        "",
        "## Primary-grid members (the 8 units behind the full-pool row)",
        "",
        "| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) | both ρ<0 |",
        "|---|---|---|---|---:|---|---|---|",
    ]
    grid = _select_paired_rows(paired, PRIMARY_UNITS)
    for r in grid.itertuples():
        flag = "yes" if r.both_negative else "no"
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.tnk_def} | {r.score} | {int(r.n)} | "
            f"{r.rho_TAC:+.3f} ({r.p_TAC:.3g}) | {r.rho_CLDN:+.3f} ({r.p_CLDN:.3g}) | {flag} |"
        )
    lines += [
        "",
        "## Primary-grid subset pools with both genes pooled ρ < 0 (k≥2)",
        "",
        "All 2⁸−1 = 255 nonempty subsets of the 8 primary units were enumerated.",
        "Rows below are both-negative subsets with N≥20, lowest TACSTD2 RE p first",
        "(top 20).",
        "",
        "| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|---|",
    ]
    show = prim_multi[prim_multi["n_patients"] >= 20].head(20)
    for r in show.itertuples():
        lines.append(_fmt_pool_row(r))
    n_prim_both = int(((primary["both_negative"]) & (primary["k"] >= 2)).sum())
    lines += [
        "",
        f"Primary-grid both-negative subset pools (k≥2): {n_prim_both} of "
        f"{int((primary['k']>=2).sum())} enumerated k≥2 subsets.",
        "",
        "## Single-cohort combinations with both genes ρ < 0 (n≥4)",
        "",
        "Every available (cohort × malignant def × T/NK def × score) that has both",
        "genes. inferCNV is absent here because CLDN4 is missing.",
        "",
        "| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) |",
        "|---|---|---|---|---:|---|---|",
    ]
    for r in singles.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.tnk_def} | {r.score} | {int(r.n)} | "
            f"{r.rho_TAC:+.3f} ({r.p_TAC:.3g}) | {r.rho_CLDN:+.3f} ({r.p_CLDN:.3g}) |"
        )
    lines += [
        "",
        f"Count: {len(singles)} both-negative single-cohort combos "
        f"(of {len(paired[paired['n']>=MIN_N])} paired combos with n≥4).",
        "",
        "## Aligned-family subset pools with both genes pooled ρ < 0 (k≥2)",
        "",
        "Families hold malignant / T/NK / score buckets fixed, then enumerate cohort",
        "subsets. Showing the 15 lowest TACSTD2 RE p among both-negative subsets with N≥15.",
        "",
        "| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|---|",
    ]
    show = other_multi[other_multi["n_patients"] >= 15].head(15)
    for r in show.itertuples():
        lines.append(_fmt_pool_row(r))
    lines += [
        "",
        f"Aligned-family both-negative subset pools (k≥2): "
        f"{int(len(other_multi))} of "
        f"{int(((~pools['malig_family'].isin(['primary_mixed','primary_nsclc_swap'])) & (pools['k']>=2)).sum())} "
        f"enumerated k≥2 aligned pools.",
        "",
        "## inferCNV-like (GSE207422 TACSTD2 only)",
        "",
        "CLDN4 is not in the wave2 table. These rows are TACSTD2 vs T/NK only.",
        "",
        "| malig | T/NK | score | n | TACSTD2 ρ (p) |",
        "|---|---|---|---:|---|",
    ]
    infer = infer.sort_values(["rho", "p"])
    for r in infer.itertuples():
        lines.append(
            f"| {r.malig_def} | {r.tnk_def} | {r.score} | {int(r.n)} | {r.rho:+.3f} ({r.p:.3g}) |"
        )
    lines += [
        "",
        "Figures: `figures/supportive_singles_forest.png`, `figures/supportive_subset_bars.png`,",
        "`figures/supportive_trio_GSE207422_GSE291670_GSE253013.png`,",
        "`figures/supportive_marker_malig_*_GSE253013_GSE291670.png`,",
        "`figures/primary_full_pool_members_forest.png`, `figures/primary_both_negative_members_forest.png`,",
        "`figures/primary_nsclc_swap_members_forest.png`, `figures/scatter_GSE253013_*.png`,",
        "`figures/scatter_GSE291670_*.png`, `figures/scatter_GSE207422_A3_given*.png`.",
        "",
    ]
    (OUT / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    atomic_rows = load_all_atomic()
    atomic_df = pd.DataFrame(atomic_rows)
    atomic_df.to_csv(OUT / "atomic_effects.tsv", sep="\t", index=False)
    paired = pair_both_genes(atomic_df)
    paired.to_csv(OUT / "paired_combos.tsv", sep="\t", index=False)
    both = paired[(paired["both_negative"]) & (paired["n"] >= MIN_N)].sort_values("rho_TAC")
    both.to_csv(OUT / "both_negative_singles.tsv", sep="\t", index=False)
    primary_grid = _select_paired_rows(paired, PRIMARY_UNITS)
    primary_grid.to_csv(OUT / "primary_grid.tsv", sep="\t", index=False)
    aligned = enumerate_pools(paired)
    prim_pools = enumerate_named_grid(paired, PRIMARY_UNITS, "primary_mixed")
    nsclc_pools = enumerate_named_grid(paired, PRIMARY_NSCLC_SWAP, "primary_nsclc_swap")
    pools = pd.concat([prim_pools, nsclc_pools, aligned], ignore_index=True)
    pools.to_csv(OUT / "subset_pools.tsv", sep="\t", index=False)
    pools[pools["both_negative"]].sort_values(["k", "p_TAC"]).to_csv(
        OUT / "both_negative_pools.tsv", sep="\t", index=False
    )
    extra_figures(paired, pools, primary_grid)
    write_finding(paired, pools, atomic_df)
    print("atomic", len(atomic_df), "paired", len(paired), "both-neg singles", len(both))
    print("pools", len(pools), "both-neg pools", int(pools["both_negative"].sum()) if len(pools) else 0)
    full = pools[(pools["malig_family"] == "primary_mixed") & (pools["is_full_family"])]
    if not full.empty:
        r = full.iloc[0]
        print(
            "FULL-POOL",
            f"k={int(r.k)} N={int(r.n_patients)}",
            f"TAC {r.rho_TAC:+.3f} p={r.p_TAC:.3g}",
            f"CLDN4 {r.rho_CLDN:+.3f} p={r.p_CLDN:.3g}",
        )
    print(both[["cohort", "malig_def", "tnk_def", "score", "n", "rho_TAC", "p_TAC", "rho_CLDN", "p_CLDN"]].to_string(index=False))


if __name__ == "__main__":
    main()

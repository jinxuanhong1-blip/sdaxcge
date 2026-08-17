#!/usr/bin/env python3
"""CLDN4-first combinatorial search on existing public lung scRNA patient scores.

Additive to PR #279 / #271 / #274. Those TACSTD2-primary combinations are taken
as given and are not re-audited. GSE207422 A3 TACSTD2 is taken as given.
CLDN4 on the same given 12-patient table is one honest row.

Patient is the unit. Full-pool is one row, not the answer.
GSE253013 uses the existing extract; the 9 GB RDS is not downloaded.
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
TNK = HERE / "data" / "tnk"
TLS = HERE / "data" / "tls"
EXH = HERE / "data" / "exh"
LEFT = HERE / "data" / "leftover"
OUT = HERE
TABLES = HERE / "tables"
FIGS = HERE / "figures"
MIN_N = 4


def _spear(x, y):
    return spearman(x, y)


def atomic(cohort, malig_def, immune_def, score, gene, x, y, **note):
    rho, p, n = _spear(x, y)
    rec = {
        "cohort": cohort,
        "malig_def": malig_def,
        "immune_def": immune_def,
        "score": score,
        "gene": gene,
        "n": n,
        "rho": rho,
        "p": p,
    }
    rec.update(note)
    return rec


def load_tnk_atomic() -> list[dict]:
    """CLDN4 (and TACSTD2 companion) vs T/NK from the PR #279 extracts."""
    rows: list[dict] = []
    given = json.loads((HERE / "data" / "given_a3.json").read_text())

    a3 = pd.read_csv(TNK / "GSE207422_drmref_patients.tsv", sep="\t")
    rows.append(
        {
            "cohort": "GSE207422",
            "malig_def": "author_DRMref",
            "immune_def": "tnk",
            "score": "mean",
            "gene": "TACSTD2",
            "n": given["n_patients"],
            "rho": given["spearman_rho_TACSTD2_mean"],
            "p": given["spearman_p_TACSTD2_mean"],
            "note": "A3 given TACSTD2; not re-audited",
            "source": "a3_given",
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
            note="CLDN4 from the given A3 12-patient table; not an A3 re-cut",
            source="a3_given_table",
        )
    )
    rows.append(
        atomic(
            "GSE207422",
            "author_DRMref",
            "cyto",
            "mean",
            "CLDN4",
            a3["malig_CLDN4_mean"],
            a3["tnk_cyto_mean"],
            note="CLDN4 vs T/NK cytotoxicity on the given A3 table",
            source="a3_given_table",
        )
    )

    fab = pd.read_csv(TNK / "GSE207422_marker_fable.tsv", sep="\t")
    post = fab[fab["timing"] == "post"].copy()
    post_mal = post[post["n_malignant"] >= 10]
    for score, tcol, ccol in (
        ("mean", "malig_TACSTD2_mean", "malig_CLDN4_mean"),
        ("pct", "malig_TACSTD2_pct", "malig_CLDN4_pct"),
    ):
        rows.append(atomic("GSE207422", "marker_malig", "tnk", score, "TACSTD2", post_mal[tcol], post_mal["tnk_frac"], note="post, n_mal>=10", source="tnk_extract"))
        rows.append(atomic("GSE207422", "marker_malig", "tnk", score, "CLDN4", post_mal[ccol], post_mal["tnk_frac"], note="post, n_mal>=10", source="tnk_extract"))
        rows.append(atomic("GSE207422", "marker_malig", "cd8", score, "CLDN4", post_mal[ccol], post_mal["cd8_frac"], note="post, n_mal>=10", source="tnk_extract"))
    if "cxcl13_log2cpm" in post_mal.columns:
        rows.append(
            atomic(
                "GSE207422",
                "marker_malig",
                "cxcl13_mean",
                "mean",
                "CLDN4",
                post_mal["malig_CLDN4_mean"],
                post_mal["cxcl13_log2cpm"],
                note="post, n_mal>=10; CXCL13 log2cpm from marker table",
                source="tnk_extract",
            )
        )

    ns = pd.read_csv(TNK / "GSE207422_marker_nsclc.tsv", sep="\t")
    ns_post = ns[ns["timing"].astype(str).str.lower() == "post"].copy()
    ns_mal = ns_post[pd.to_numeric(ns_post["n_malignant_like"], errors="coerce").fillna(0) >= 10]
    for score, tcol, ccol in (
        ("mean", "mal_TACSTD2_mean", "mal_CLDN4_mean"),
        ("pct", "mal_TACSTD2_pct_pos", "mal_CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "tnk", score, "TACSTD2", ns_mal[tcol], ns_mal["frac_tnk"], note="post, n_mal_like>=10", source="tnk_extract"))
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "tnk", score, "CLDN4", ns_mal[ccol], ns_mal["frac_tnk"], note="post, n_mal_like>=10", source="tnk_extract"))
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "cd8", score, "CLDN4", ns_mal[ccol], ns_mal["frac_cd8"], note="post, n_mal_like>=10", source="tnk_extract"))

    d = pd.read_csv(TNK / "GSE205335_patients.tsv", sep="\t")
    nsclc = d[d["cancer_subtype"].isin(["ADC", "SQ"])]
    for cohort, df in (("GSE205335", d), ("GSE205335_NSCLC", nsclc)):
        for immune, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8"), ("b_plasma", "frac_b_plasma")):
            rows.append(atomic(cohort, "author_malig", immune, "mean", "TACSTD2", df["mal_TACSTD2_mean"], df[y], note="author malignant", source="tnk_extract"))
            rows.append(atomic(cohort, "author_malig", immune, "mean", "CLDN4", df["mal_CLDN4_mean"], df[y], note="author malignant", source="tnk_extract"))
            rows.append(atomic(cohort, "author_malig", immune, "pct", "CLDN4", df["mal_CLDN4_pct_pos"], df[y], note="author malignant", source="tnk_extract"))

    for name in ("GSE241934_IIT", "GSE241934_Real"):
        df = pd.read_csv(TNK / f"{name}_patients.tsv", sep="\t")
        elig = df[df["eligible"] == 1]
        for immune, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8")):
            rows.append(atomic(name, "author_epi", immune, "mean", "TACSTD2", elig["TACSTD2_log1p_cp10k"], elig[y], note="author residual Epi", source="tnk_extract"))
            rows.append(atomic(name, "author_epi", immune, "mean", "CLDN4", elig["CLDN4_log1p_cp10k"], elig[y], note="author residual Epi", source="tnk_extract"))
            rows.append(atomic(name, "author_epi", immune, "pct", "CLDN4", elig["CLDN4_pct_pos"], elig[y], note="author residual Epi", source="tnk_extract"))

    d = pd.read_csv(TNK / "GSE291670_patients.tsv", sep="\t")
    for immune, y in (("tnk", "frac_lineage_tnk"), ("tnk_umi", "frac_umi_tnk")):
        rows.append(atomic("GSE291670", "marker_malig", immune, "mean", "TACSTD2", d["mal_TACSTD2_mean_log1p_cp10k"], d[y], note="marker malignant", source="tnk_extract"))
        rows.append(atomic("GSE291670", "marker_malig", immune, "mean", "CLDN4", d["mal_CLDN4_mean_log1p_cp10k"], d[y], note="marker malignant", source="tnk_extract"))
        rows.append(atomic("GSE291670", "marker_malig", immune, "pct", "CLDN4", d["mal_CLDN4_pct_pos"], d[y], note="marker malignant", source="tnk_extract"))
        rows.append(atomic("GSE291670", "marker_epi", immune, "mean", "CLDN4", d["epi_CLDN4_mean_log1p_cp10k"], d[y], note="all epithelial", source="tnk_extract"))

    d = pd.read_csv(TNK / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    te = d[(d["tissue"] == "Tumor") & (d["eligible_epi"])]
    for score, tcol, ccol in (
        ("mean", "TACSTD2_mean_log1p", "CLDN4_mean_log1p"),
        ("mean_cp10k", "TACSTD2_mean_log1p_cp10k", "CLDN4_mean_log1p_cp10k"),
        ("pct", "TACSTD2_pct_pos", "CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE253013", "marker_malig", "tnk", score, "TACSTD2", t[tcol], t["tnk_fraction"], note="tumor malig-like; existing extract, no RDS", source="tnk_extract"))
        rows.append(atomic("GSE253013", "marker_malig", "tnk", score, "CLDN4", t[ccol], t["tnk_fraction"], note="tumor malig-like; existing extract, no RDS", source="tnk_extract"))
    for score, tcol, ccol in (
        ("mean", "epi_TACSTD2_mean_log1p", "epi_CLDN4_mean_log1p"),
        ("pct", "epi_TACSTD2_pct_pos", "epi_CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE253013", "marker_epi", "tnk", score, "CLDN4", te[ccol], te["tnk_fraction"], note="tumor epithelium; existing extract", source="tnk_extract"))

    d = pd.read_csv(TNK / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    tlung = d[d["origin"] == "tLung"]
    for cohort, df in (("GSE131907", tumor), ("GSE131907_tLung", tlung)):
        epi = df[df["n_epithelial"] >= 20]
        mal = df[df["n_malignant"] >= 20]
        for immune, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8"), ("b", "frac_b")):
            rows.append(atomic(cohort, "author_epi", immune, "mean", "TACSTD2", epi["epi_TACSTD2_mean"], epi[y], note="author epithelium, >=20 epi", source="tnk_extract"))
            rows.append(atomic(cohort, "author_epi", immune, "mean", "CLDN4", epi["epi_CLDN4_mean"], epi[y], note="author epithelium, >=20 epi", source="tnk_extract"))
            rows.append(atomic(cohort, "author_epi", immune, "pct", "CLDN4", epi["epi_CLDN4_pct"], epi[y], note="author epithelium, >=20 epi", source="tnk_extract"))
            if len(mal) >= 3:
                rows.append(atomic(cohort, "author_malig", immune, "mean", "CLDN4", mal["mal_CLDN4_mean"], mal[y], note="author Malignant cells subtype", source="tnk_extract"))
                rows.append(atomic(cohort, "author_malig", immune, "pct", "CLDN4", mal["mal_CLDN4_pct"], mal[y], note="author Malignant cells subtype", source="tnk_extract"))

    d = pd.read_csv(TNK / "GSE325414_donors.csv")
    rows.append(atomic("GSE325414", "author_malig", "tnk", "mean", "TACSTD2", d["malignant_TACSTD2_mean"], d["frac_TNK"], note="author malignant, donor", source="tnk_extract"))
    rows.append(atomic("GSE325414", "author_malig", "tnk", "mean", "CLDN4", d["malignant_CLDN4_mean"], d["frac_TNK"], note="author malignant, donor", source="tnk_extract"))
    return rows


def load_tls_atomic() -> list[dict]:
    """CLDN4 vs B and CXCL13+ from the PR #274 patient table. TACSTD2-vs-B meta is given."""
    rows: list[dict] = []
    d = pd.read_csv(TLS / "patient_scores.tsv", sep="\t")
    elig = d[d["eligible"] == True].copy()
    # GSE241934_RWC in TLS == GSE241934_Real in the T/NK extract
    elig["cohort"] = elig["cohort"].replace({"GSE241934_RWC": "GSE241934_Real"})
    for cohort, g in elig.groupby("cohort"):
        malig = "author_epi" if "241934" in cohort else "malig_or_epi"
        if cohort == "GSE207422":
            malig = "author_epi_tls"
        rows.append(
            atomic(
                cohort,
                malig,
                "b",
                "mean",
                "CLDN4",
                g["cldn4_mal_mean_log1p_cp10k"],
                g["frac_B"],
                note="PR #274 extract; TACSTD2-vs-B meta taken as given",
                source="tls_extract",
            )
        )
        rows.append(
            atomic(
                cohort,
                malig,
                "cxcl13pos",
                "mean",
                "CLDN4",
                g["cldn4_mal_mean_log1p_cp10k"],
                g["frac_CXCL13pos_T"],
                note="PR #274 extract; frac CXCL13+ among T",
                source="tls_extract",
            )
        )
        rows.append(
            atomic(
                cohort,
                malig,
                "cxcl13_mean",
                "mean",
                "CLDN4",
                g["cldn4_mal_mean_log1p_cp10k"],
                g["cxcl13_mean_log1p_cp10k"],
                note="PR #274 extract; CXCL13 mean log1p(CP10k)",
                source="tls_extract",
            )
        )
    return rows


def load_exh_atomic() -> list[dict]:
    """CLDN4 vs T/NK cytotoxicity (cyto-high T score) from PR #260."""
    rows: list[dict] = []
    d = pd.read_csv(EXH / "patient_scores.tsv", sep="\t")
    d["cohort"] = d["cohort"].replace({"GSE241934_RWC": "GSE241934_Real"})
    for cohort, g in d.groupby("cohort"):
        keep = g[g["mal_CLDN4"].notna() & g["tnk_cyto"].notna() & (pd.to_numeric(g["n_malignant"], errors="coerce").fillna(0) >= 10)]
        if len(keep) < 3:
            continue
        malig = "author_epi" if "241934" in cohort else "marker_or_author_malig"
        if cohort == "GSE205335":
            malig = "author_malig"
        if cohort == "GSE291670":
            malig = "marker_malig"
        if cohort == "GSE207422":
            malig = "marker_malig"
        rows.append(
            atomic(
                cohort,
                malig,
                "cyto",
                "mean",
                "CLDN4",
                keep["mal_CLDN4"],
                keep["tnk_cyto"],
                note="PR #260 extract; T/NK cytotoxicity mean log1p(CP10k); n_mal>=10",
                source="exh_extract",
            )
        )
    return rows


def load_leftover_atomic() -> list[dict]:
    """Extra n from PR #280 leftover extracts. Not in the primary 8-unit grid."""
    rows: list[dict] = []
    g = pd.read_csv(LEFT / "GSE267108_per_sample.tsv", sep="\t")
    keep = g[(g["n_epi"] >= 10) & (g["n_tnk"] >= 20)]
    rows.append(
        atomic(
            "GSE267108",
            "marker_epi",
            "tnk",
            "mean",
            "CLDN4",
            keep["epi_CLDN4"],
            keep["frac_tnk"],
            note="leftover PR #280; treatment-naive LUAD; extra n",
            source="leftover",
        )
    )
    g = pd.read_csv(LEFT / "GSE274595_per_sample.tsv", sep="\t")
    keep = g[(g["n_epi"] >= 10) & (g["n_tnk"] >= 20)]
    rows.append(
        atomic(
            "GSE274595",
            "marker_epi",
            "tnk",
            "mean",
            "CLDN4",
            keep["epi_CLDN4"],
            keep["frac_tnk"],
            note="leftover PR #280; surgical, not ICI; extra n",
            source="leftover",
        )
    )
    g = pd.read_csv(LEFT / "EMTAB13526_per_patient.tsv", sep="\t")
    keep = g[(g["subset"] == "CD235a_minus") & (g["n_epi"] >= 20) & (g["n_tnk"] >= 20)]
    rows.append(
        atomic(
            "EMTAB13526",
            "author_epi",
            "tnk",
            "mean",
            "CLDN4",
            keep["epi_CLDN4"],
            keep["frac_tnk"],
            note="leftover PR #280 / #259; Cvejic atlas, not ICI; extra n",
            source="leftover",
        )
    )
    return rows


# Locked 8-unit primary T/NK grid (same members as PR #279). CLDN4-first ranking.
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


def _select_rows(df: pd.DataFrame, specs: list[tuple], gene: str = "CLDN4") -> pd.DataFrame:
    keys = ["cohort", "malig_def", "immune_def", "score"]
    sub = df[df["gene"] == gene].set_index(keys)
    rows = []
    for spec in specs:
        if spec not in sub.index:
            raise KeyError(f"missing {gene} combo {spec}")
        row = sub.loc[spec]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        rec = row.to_dict()
        rec.update(dict(zip(keys, spec)))
        rec["gene"] = gene
        rows.append(rec)
    return pd.DataFrame(rows)


def _pool_one(rhos, ns, ps) -> dict:
    rhos = [float(r) for r in rhos]
    ns = [int(n) for n in ns]
    ps = [float(p) for p in ps]
    if len(rhos) == 1:
        re = {
            "k": 1,
            "n_patients_total": ns[0],
            "pooled_rho": rhos[0],
            "p": ps[0],
            "I2": 0.0,
            "ci95_rho": [float("nan"), float("nan")],
        }
    else:
        re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    return re, st


def enumerate_named_grid(cldn: pd.DataFrame, specs: list[tuple], family: str) -> pd.DataFrame:
    grid = _select_rows(cldn, specs, "CLDN4")
    units = [s[0] for s in specs]
    out = []
    for k in range(1, len(units) + 1):
        for subset in itertools.combinations(units, k):
            sub = grid[grid["cohort"].isin(subset)]
            re, st = _pool_one(sub["rho"], sub["n"], sub["p"])
            out.append(
                {
                    "family": family,
                    "immune_family": "tnk",
                    "score_family": "mean",
                    "k": len(subset),
                    "cohorts": "+".join(subset),
                    "n_patients": int(sub["n"].sum()),
                    "rho_CLDN": re.get("pooled_rho"),
                    "p_CLDN": re.get("p"),
                    "I2_CLDN": re.get("I2"),
                    "stouffer_z_CLDN": st.get("z"),
                    "stouffer_p_CLDN": st.get("p"),
                    "cldn4_negative": float(re.get("pooled_rho", 0) or 0) < 0,
                    "is_full_family": len(subset) == len(units),
                    "n_family_units": len(units),
                }
            )
    return pd.DataFrame(out)


def enumerate_aligned(cldn: pd.DataFrame) -> pd.DataFrame:
    """Hold immune + score bucket fixed; enumerate cohort subsets (CLDN4-first)."""
    unit_ok = {
        "GSE207422",
        "GSE205335",
        "GSE241934_IIT",
        "GSE241934_Real",
        "GSE291670",
        "GSE253013",
        "GSE131907",
        "GSE325414",
        "GSE148071",
        "GSE154826",
        "GSE267108",
        "GSE233203",
    }
    malig_bucket = {
        "author_DRMref": "author",
        "author_malig": "author",
        "author_epi": "author",
        "author_epi_tls": "tls",
        "malig_or_epi": "tls",
        "marker_malig": "marker",
        "marker_malig_nsclc": "marker",
        "marker_epi": "marker_epi",
        "marker_or_author_malig": "exh",
    }
    immune_bucket = {
        "tnk": "tnk",
        "tnk_umi": "tnk_umi",
        "cd8": "cd8",
        "b": "b",
        "b_plasma": "b",
        "cxcl13pos": "cxcl13pos",
        "cxcl13_mean": "cxcl13_mean",
        "cyto": "cyto",
    }
    score_bucket = {"mean": "mean", "mean_cp10k": "mean", "pct": "pct"}
    work = cldn[cldn["cohort"].isin(unit_ok)].copy()
    work["malig_b"] = work["malig_def"].map(malig_bucket)
    work["imm_b"] = work["immune_def"].map(immune_bucket)
    work["score_b"] = work["score"].map(score_bucket)
    work = work.dropna(subset=["malig_b", "imm_b", "score_b"])
    work = work[work["n"] >= MIN_N]
    pref = {
        "author": ["author_malig", "author_DRMref", "author_epi"],
        "marker": ["marker_malig", "marker_malig_nsclc"],
        "marker_epi": ["marker_epi"],
        "tls": ["malig_or_epi", "author_epi_tls", "author_epi"],
        "exh": ["marker_or_author_malig", "author_malig", "marker_malig", "author_epi"],
    }
    out = []
    for (mb, ib, sb), g in work.groupby(["malig_b", "imm_b", "score_b"], dropna=False):
        chosen = []
        for unit, ug in g.groupby("cohort"):
            order = pref.get(mb, list(ug["malig_def"].unique()))
            ug = ug.copy()
            ug["_rank"] = ug["malig_def"].apply(lambda x: order.index(x) if x in order else 99)
            chosen.append(ug.sort_values(["_rank", "n"], ascending=[True, False]).head(1))
        if not chosen:
            continue
        grid = pd.concat(chosen, ignore_index=True)
        units = list(grid["cohort"])
        for k in range(1, len(units) + 1):
            for subset in itertools.combinations(units, k):
                sub = grid[grid["cohort"].isin(subset)]
                if len(sub) != len(subset):
                    continue
                re, st = _pool_one(sub["rho"], sub["n"], sub["p"])
                out.append(
                    {
                        "family": f"{mb}/{ib}/{sb}",
                        "immune_family": ib,
                        "score_family": sb,
                        "k": len(subset),
                        "cohorts": "+".join(subset),
                        "n_patients": int(sub["n"].sum()),
                        "rho_CLDN": re.get("pooled_rho"),
                        "p_CLDN": re.get("p"),
                        "I2_CLDN": re.get("I2"),
                        "stouffer_z_CLDN": st.get("z"),
                        "stouffer_p_CLDN": st.get("p"),
                        "cldn4_negative": float(re.get("pooled_rho", 0) or 0) < 0,
                        "is_full_family": len(subset) == len(units),
                        "n_family_units": len(units),
                    }
                )
    return pd.DataFrame(out)


def forest_cldn(sub: pd.DataFrame, title: str, path: Path, immune_label: str = "listed immune") -> None:
    sub = sub.sort_values("n", ascending=False).copy()
    labels = [f"{r.cohort} n={int(r.n)}" for r in sub.itertuples()]
    fig, ax = plt.subplots(figsize=(8.4, 1.05 + 0.40 * len(sub)))
    y = np.arange(len(sub))
    for i, r in enumerate(sub.itertuples()):
        rho = float(r.rho)
        p = float(r.p)
        color = "#7a2d0b" if rho < 0 else "#4a4a4a"
        ax.plot(rho, i, "o", color=color, ms=7)
        ax.text(
            0.98,
            i,
            f"ρ={rho:+.2f} p={p:.3g}",
            va="center",
            ha="right",
            fontsize=7,
            family="monospace",
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_xlim(-1.05, 1.05)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"Spearman ρ  (malignant/epithelial CLDN4 vs {immune_label})")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter_combo(path: Path, x, y, xlabel, ylabel, title, color="#7a2d0b") -> None:
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    ax.scatter(list(x), list(y), c=color, s=36, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def extra_figures(cldn: pd.DataFrame, pools: pd.DataFrame, primary_grid: pd.DataFrame) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)

    forest_cldn(
        primary_grid,
        "Primary 8-unit grid · CLDN4 mean vs T/NK (full-pool members)",
        FIGS / "primary_full_pool_members_forest.png",
        "T/NK",
    )
    prim_neg = primary_grid[primary_grid["rho"] < 0]
    if len(prim_neg) >= 2:
        forest_cldn(
            prim_neg,
            "Primary-grid members with CLDN4 ρ<0 vs T/NK",
            FIGS / "primary_cldn4_negative_members_forest.png",
            "T/NK",
        )

    singles = cldn[(cldn["immune_def"] == "tnk") & (cldn["n"] >= MIN_N) & (cldn["rho"] < 0)].sort_values(["p", "rho"])
    if not singles.empty:
        forest_cldn(
            singles.head(14),
            "Single-cohort CLDN4 ρ<0 vs T/NK (lowest p first)",
            FIGS / "cldn4_tnk_neg_singles_forest.png",
            "T/NK",
        )

    # Best T/NK multi-cohort CLDN4-negative (not the full pool)
    multi = pools[
        (pools["cldn4_negative"])
        & (pools["k"] >= 2)
        & (pools["immune_family"] == "tnk")
        & (~pools["is_full_family"])
    ].copy()
    if not multi.empty:
        multi = multi.sort_values(["p_CLDN", "rho_CLDN"])
        top = multi[multi["n_patients"] >= 15].head(10)
        if top.empty:
            top = multi.head(10)
        fig, ax = plt.subplots(figsize=(9.0, 4.6))
        y = np.arange(len(top))
        ax.barh(y, top["rho_CLDN"], color="#7a2d0b")
        ax.axvline(0, color="#888", lw=0.8)
        labels = [
            f"k={int(r.k)} N={int(r.n_patients)} p={r.p_CLDN:.3g}  {r.family}"
            for r in top.itertuples()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Pooled CLDN4 Spearman ρ")
        ax.set_title("CLDN4-negative subset pools vs T/NK (lowest CLDN4 p; full-pool excluded)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIGS / "cldn4_tnk_subset_bars.png", dpi=150)
        fig.savefig(FIGS / "cldn4_tnk_subset_bars.pdf")
        plt.close(fig)

    # Marker pair — strongest CLDN4 multi on T/NK
    marker_pair = cldn[
        (cldn["cohort"].isin(["GSE253013", "GSE291670"]))
        & (cldn["malig_def"] == "marker_malig")
        & (cldn["immune_def"] == "tnk")
        & (cldn["score"].isin(["mean", "pct"]))
    ]
    if len(marker_pair) >= 2:
        forest_cldn(
            marker_pair[marker_pair["score"] == "mean"],
            "Best CLDN4 T/NK cut · marker-malignant mean: GSE253013 + GSE291670",
            FIGS / "best_marker_malig_mean_GSE253013_GSE291670.png",
            "T/NK",
        )
        pct = marker_pair[marker_pair["score"] == "pct"]
        if len(pct) >= 2:
            forest_cldn(
                pct,
                "Best CLDN4 T/NK cut · marker-malignant %pos: GSE253013 + GSE291670",
                FIGS / "best_marker_malig_pct_GSE253013_GSE291670.png",
                "T/NK",
            )

    # GSE131907 author malignant %pos — largest-n CLDN4-negative single with p<0.05
    g131 = cldn[
        (cldn["cohort"] == "GSE131907")
        & (cldn["malig_def"] == "author_malig")
        & (cldn["immune_def"].isin(["tnk", "cd8"]))
        & (cldn["score"] == "pct")
    ]
    if not g131.empty:
        forest_cldn(
            g131,
            "GSE131907 author-malignant CLDN4 %pos vs T/NK and CD8",
            FIGS / "best_GSE131907_author_malig_pct.png",
            "T/NK or CD8",
        )

    # Best aligned T/NK %pos pair (CLDN4-first; GSE205335 TACSTD2 is the other sign)
    pct_pair = cldn[
        (cldn["cohort"].isin(["GSE131907", "GSE205335"]))
        & (cldn["malig_def"].isin(["author_malig", "author_epi"]))
        & (cldn["immune_def"] == "tnk")
        & (cldn["score"] == "pct")
    ].copy()
    if "author_malig" in set(pct_pair["malig_def"]):
        pct_pair = pct_pair[pct_pair["malig_def"] == "author_malig"]
    if len(pct_pair) >= 2:
        forest_cldn(
            pct_pair,
            "Best CLDN4 T/NK cut · author %pos: GSE131907 + GSE205335",
            FIGS / "best_author_pct_GSE131907_GSE205335.png",
            "T/NK",
        )

    cx_trio = cldn[
        (cldn["cohort"].isin(["GSE148071", "GSE207422", "GSE253013"]))
        & (cldn["immune_def"] == "cxcl13pos")
        & (cldn["score"] == "mean")
    ]
    if len(cx_trio) >= 2:
        forest_cldn(
            cx_trio,
            "Best CLDN4 CXCL13+ cut · GSE148071 + GSE207422 + GSE253013",
            FIGS / "best_cxcl13pos_GSE148071_GSE207422_GSE253013.png",
            "CXCL13+ T fraction",
        )

    # CXCL13+ CLDN4-negative singles
    cx = cldn[(cldn["immune_def"] == "cxcl13pos") & (cldn["n"] >= MIN_N) & (cldn["rho"] < 0)]
    if not cx.empty:
        forest_cldn(
            cx.sort_values("p"),
            "CLDN4 ρ<0 vs CXCL13+ T fraction (PR #274 extract)",
            FIGS / "cldn4_cxcl13pos_neg_forest.png",
            "CXCL13+ T fraction",
        )
    bneg = cldn[(cldn["immune_def"].isin(["b", "b_plasma"])) & (cldn["n"] >= MIN_N) & (cldn["rho"] < 0)]
    if not bneg.empty:
        forest_cldn(
            bneg.sort_values("p"),
            "CLDN4 ρ<0 vs B fraction",
            FIGS / "cldn4_B_neg_forest.png",
            "B fraction",
        )
    cy = cldn[(cldn["immune_def"] == "cyto") & (cldn["n"] >= MIN_N)]
    if not cy.empty:
        forest_cldn(
            cy.sort_values("rho"),
            "CLDN4 vs T/NK cytotoxicity (all scored cohorts)",
            FIGS / "cldn4_cyto_forest.png",
            "T/NK cytotoxicity",
        )

    # Patient scatters for the best CLDN4 T/NK cuts
    d = pd.read_csv(TNK / "GSE291670_patients.tsv", sep="\t")
    scatter_combo(
        FIGS / "scatter_GSE291670_CLDN4_tnk.png",
        d["frac_lineage_tnk"],
        d["mal_CLDN4_mean_log1p_cp10k"],
        "Lineage T/NK fraction",
        "Malignant CLDN4 mean log1p(CP10k)",
        "GSE291670 marker_malig / tnk / mean (n=6)",
    )
    d = pd.read_csv(TNK / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    scatter_combo(
        FIGS / "scatter_GSE253013_CLDN4_tnk.png",
        t["tnk_fraction"],
        t["CLDN4_mean_log1p"],
        "T/NK fraction",
        "Malignant-like CLDN4 mean log1p",
        "GSE253013 marker_malig / tnk / mean (n=9; existing extract)",
    )
    d = pd.read_csv(TNK / "GSE131907_samples.tsv", sep="\t")
    mal = d[(d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])) & (d["n_malignant"] >= 20)]
    scatter_combo(
        FIGS / "scatter_GSE131907_CLDN4_tnk_pct.png",
        mal["frac_tnk"],
        mal["mal_CLDN4_pct"],
        "T/NK fraction",
        "Malignant CLDN4 %pos",
        "GSE131907 author_malig / tnk / %pos",
    )
    a3 = pd.read_csv(TNK / "GSE207422_drmref_patients.tsv", sep="\t")
    scatter_combo(
        FIGS / "scatter_GSE207422_CLDN4_tnk.png",
        a3["frac_tnk"],
        a3["malig_CLDN4_mean"],
        "T/NK fraction (DRMref)",
        "Malignant CLDN4 mean log1p(CP10k)",
        "GSE207422 given table · CLDN4 vs T/NK (n=12)",
    )


def _fmt(r) -> str:
    return (
        f"| {int(r.k)} | {int(r.n_patients)} | {r.family} | "
        f"{r.rho_CLDN:+.3f} ({r.p_CLDN:.3g}, {r.I2_CLDN:.0f}%) | {r.cohorts} |"
    )


def write_finding(cldn: pd.DataFrame, pools: pd.DataFrame, primary_grid: pd.DataFrame) -> None:
    full = pools[(pools["family"] == "primary_mixed") & (pools["is_full_family"])]
    prim = pools[pools["family"] == "primary_mixed"].copy()
    prim_neg = prim[(prim["cldn4_negative"]) & (prim["k"] >= 2)].sort_values(["p_CLDN", "rho_CLDN"])
    tnk_singles = cldn[(cldn["immune_def"] == "tnk") & (cldn["n"] >= MIN_N)].sort_values(["p", "rho"])
    tnk_neg = tnk_singles[tnk_singles["rho"] < 0]
    cx_neg = cldn[(cldn["immune_def"].isin(["cxcl13pos", "cxcl13_mean"])) & (cldn["n"] >= MIN_N) & (cldn["rho"] < 0)].sort_values(["p", "rho"])
    b_neg = cldn[(cldn["immune_def"].isin(["b", "b_plasma"])) & (cldn["n"] >= MIN_N) & (cldn["rho"] < 0)].sort_values(["p", "rho"])
    cy = cldn[(cldn["immune_def"] == "cyto") & (cldn["n"] >= MIN_N)].sort_values(["rho", "p"])
    leftover = cldn[cldn["source"] == "leftover"]
    aligned_neg = pools[
        (pools["cldn4_negative"])
        & (pools["k"] >= 2)
        & (pools["family"] != "primary_mixed")
        & (pools["n_patients"] >= 15)
    ].sort_values(["p_CLDN", "rho_CLDN"])

    lines = [
        "# Combinatorial search: malignant CLDN4 vs T/NK (and CXCL13+ / B / cyto-high T)",
        "",
        "CLDN4-first. Patient is the unit. Numbers are computed Spearman /",
        "DerSimonian–Laird / Stouffer values across definitions and cohort subsets;",
        "p-values are descriptive.",
        "",
        "TACSTD2-primary combinations in PR #279 / #271 / #274 are **taken as given**",
        "and are not re-audited. GSE207422 A3 TACSTD2 is taken as given. CLDN4 on the",
        "same given 12-patient table is one honest row (not an A3 re-cut).",
        "",
        "GSE253013 uses the existing 9-patient tumor extract. The 9 GB GEO RDS was",
        "not downloaded.",
        "",
        "## Full-pool row (one row, not the answer)",
        "",
        "Primary 8-unit mixed-definition grid, mean log1p vs T/NK. Same members as",
        "PR #279. This is the merge-everything row.",
        "",
        "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    if not full.empty:
        lines.append(_fmt(full.iloc[0]))
    lines += [
        "",
        "The search below is the answer: combinations where **CLDN4 ρ < 0** vs T/NK",
        "or CXCL13+, with honest n and p.",
        "",
        "## Combinations that recover CLDN4 ρ < 0 vs T/NK",
        "",
        "Highlighted cuts (descriptive p; ranked by CLDN4, not TACSTD2):",
        "",
    ]
    highlight = [
        ("author/tnk/pct", "GSE131907+GSE205335"),
        ("marker/tnk/pct", "GSE253013+GSE291670"),
        ("marker/tnk/mean", "GSE253013+GSE291670"),
        ("primary_mixed", "GSE291670"),
        ("primary_mixed", "GSE291670+GSE253013"),
        ("primary_mixed", "GSE205335+GSE291670+GSE253013+GSE131907"),
        ("primary_mixed", "GSE207422+GSE291670+GSE253013"),
        ("tls/cxcl13pos/mean", "GSE148071+GSE207422+GSE253013"),
        ("tls/cxcl13pos/mean", "GSE148071+GSE207422"),
        ("tls/cxcl13_mean/mean", "GSE207422+GSE253013"),
    ]
    seen = set()
    for fam, cohorts in highlight:
        hit = pools[(pools["family"] == fam) & (pools["cohorts"] == cohorts)]
        for r in hit.itertuples():
            key = (r.family, r.cohorts)
            if key in seen:
                continue
            seen.add(key)
            lines.append(
                f"- {r.cohorts} · {r.family} · k={int(r.k)} · N={int(r.n_patients)} · "
                f"CLDN4 ρ={r.rho_CLDN:+.3f} p={r.p_CLDN:.3g} I²={r.I2_CLDN:.0f}%"
            )
    # seven CLDN4-negative primary members (drop Real, which is ρ>0)
    seven = "GSE207422+GSE205335+GSE241934_IIT+GSE291670+GSE253013+GSE131907+GSE325414"
    hit = pools[(pools["family"] == "primary_mixed") & (pools["cohorts"] == seven)]
    if not hit.empty:
        r = hit.iloc[0]
        lines.append(
            f"- seven primary members with single-cohort CLDN4 ρ<0 (drop GSE241934_Real) · "
            f"k={int(r.k)} · N={int(r.n_patients)} · "
            f"CLDN4 ρ={r.rho_CLDN:+.3f} p={r.p_CLDN:.3g} I²={r.I2_CLDN:.0f}%"
        )
    lines += [
        "",
        "## Primary-grid members (the 8 units behind the full-pool row)",
        "",
        "| cohort | malig | immune | score | n | CLDN4 ρ (p) | CLDN4 ρ<0 |",
        "|---|---|---|---|---:|---|---|",
    ]
    for r in primary_grid.itertuples():
        flag = "yes" if r.rho < 0 else "no"
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) | {flag} |"
        )
    lines += [
        "",
        "## Primary-grid subset pools with CLDN4 pooled ρ < 0 (k≥2)",
        "",
        "All 2⁸−1 = 255 nonempty subsets of the 8 primary units were enumerated.",
        "Rows below are CLDN4-negative subsets with N≥20, lowest CLDN4 RE p first",
        "(top 20). Full-pool is not repeated here.",
        "",
        "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    show = prim_neg[(prim_neg["n_patients"] >= 20) & (~prim_neg["is_full_family"])].head(20)
    for r in show.itertuples():
        lines.append(_fmt(r))
    n_prim_neg = int(((prim["cldn4_negative"]) & (prim["k"] >= 2)).sum())
    lines += [
        "",
        f"Primary-grid CLDN4-negative subset pools (k≥2): {n_prim_neg} of "
        f"{int((prim['k']>=2).sum())} enumerated k≥2 subsets.",
        "",
        "## Single-cohort CLDN4 ρ < 0 vs T/NK (n≥4)",
        "",
        "Every available (cohort × malignant def × T/NK × score) with CLDN4 ρ<0.",
        "",
        "| cohort | malig | immune | score | n | CLDN4 ρ (p) | source |",
        "|---|---|---|---|---:|---|---|",
    ]
    for r in tnk_neg.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) | {r.source} |"
        )
    lines += [
        "",
        f"Count: {len(tnk_neg)} CLDN4-negative T/NK singles "
        f"(of {len(tnk_singles)} T/NK singles with n≥4).",
        "",
        "## CLDN4 ρ < 0 vs CXCL13+ (n≥4)",
        "",
        "CXCL13+ T fraction and CXCL13 mean from the PR #274 extract. The TACSTD2",
        "UCell-vs-CXCL13+ grid in PR #271 is TACSTD2-primary and is not re-audited.",
        "",
        "| cohort | malig | immune | score | n | CLDN4 ρ (p) |",
        "|---|---|---|---|---:|---|",
    ]
    for r in cx_neg.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) |"
        )
    if cx_neg.empty:
        lines.append("| — | — | — | — | — | no CLDN4 ρ<0 row with n≥4 |")
    cx_pools = aligned_neg[aligned_neg["immune_family"].isin(["cxcl13pos", "cxcl13_mean"])].head(12)
    if not cx_pools.empty:
        lines += [
            "",
            "CXCL13+ / CXCL13-mean subset pools with CLDN4 ρ<0, N≥15, lowest p first:",
            "",
            "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in cx_pools.itertuples():
            lines.append(_fmt(r))
    lines += [
        "",
        "## CLDN4 ρ < 0 vs B (n≥4)",
        "",
        "B fraction from the PR #274 extract plus B/plasma columns already on the",
        "T/NK tables. PR #274 TACSTD2-vs-B RE meta is taken as given.",
        "",
        "| cohort | malig | immune | score | n | CLDN4 ρ (p) |",
        "|---|---|---|---|---:|---|",
    ]
    for r in b_neg.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) |"
        )
    b_pools = aligned_neg[aligned_neg["immune_family"] == "b"].head(10)
    if not b_pools.empty:
        lines += [
            "",
            "B-fraction subset pools with CLDN4 ρ<0, N≥15, lowest p first:",
            "",
            "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in b_pools.itertuples():
            lines.append(_fmt(r))
    lines += [
        "",
        "## CLDN4 vs cyto-high T / T/NK cytotoxicity",
        "",
        "Continuous T/NK cytotoxicity score (GZMB/PRF1/GNLY/NKG7 mean log1p) from",
        "the PR #260 extract and the given GSE207422 table. This is the public",
        "cyto-high T score that already exists; PR #271 cyto-high T fractions are",
        "TACSTD2-primary (UCell module) and are not re-audited.",
        "",
        "| cohort | malig | immune | score | n | CLDN4 ρ (p) |",
        "|---|---|---|---|---:|---|",
    ]
    for r in cy.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) |"
        )
    cy_neg_pools = aligned_neg[aligned_neg["immune_family"] == "cyto"].head(8)
    if not cy_neg_pools.empty:
        lines += [
            "",
            "Cytotoxicity subset pools with CLDN4 ρ<0, N≥15:",
            "",
            "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in cy_neg_pools.itertuples():
            lines.append(_fmt(r))
    lines += [
        "",
        "## Extra leftover n (not in the primary 8-unit grid)",
        "",
        "| cohort | malig | immune | n | CLDN4 ρ (p) | note |",
        "|---|---|---|---:|---|---|",
    ]
    for r in leftover.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) | {r.note} |"
        )
    lines += [
        "",
        "## Aligned-family CLDN4-negative pools (k≥2, N≥15, lowest CLDN4 p)",
        "",
        "Families hold malignant / immune / score buckets fixed. Top 15.",
        "",
        "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    for r in aligned_neg.head(15).itertuples():
        lines.append(_fmt(r))
    lines += [
        "",
        "## Given TACSTD2-primary combos (not re-audited)",
        "",
        "From PR #279 (TACSTD2-first): GSE207422+GSE291670+GSE253013 mean/T/NK",
        "N=27 TACSTD2 ρ=−0.585 p=0.00452; marker-malignant %pos GSE253013+GSE291670",
        "N=15 TACSTD2 ρ=−0.736 p=0.00766. From PR #271: GSE207422 UCell",
        "{TACSTD2, CLDN4, EPCAM} vs CXCL13+ ρ=−0.657 p=0.020 (n=12). From PR #274:",
        "TACSTD2 vs B RE ρ=−0.110 p=0.47 (k=7, n=162). Those rows stay there.",
        "",
        "Figures: `figures/primary_full_pool_members_forest.png`,",
        "`figures/primary_cldn4_negative_members_forest.png`,",
        "`figures/cldn4_tnk_neg_singles_forest.png`, `figures/cldn4_tnk_subset_bars.png`,",
        "`figures/best_marker_malig_*_GSE253013_GSE291670.png`,",
        "`figures/best_author_pct_GSE131907_GSE205335.png`,",
        "`figures/best_cxcl13pos_GSE148071_GSE207422_GSE253013.png`,",
        "`figures/best_GSE131907_author_malig_pct.png`,",
        "`figures/cldn4_cxcl13pos_neg_forest.png`, `figures/cldn4_B_neg_forest.png`,",
        "`figures/cldn4_cyto_forest.png`, `figures/scatter_*.png`.",
        "",
        "Reproduce: `python3 methods/scrna_cldn4_combo/combinatorial_search.py`",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rows = load_tnk_atomic() + load_tls_atomic() + load_exh_atomic() + load_leftover_atomic()
    atomic_df = pd.DataFrame(rows)
    atomic_df.to_csv(TABLES / "atomic_effects.tsv", sep="\t", index=False)
    cldn = atomic_df[atomic_df["gene"] == "CLDN4"].copy()
    cldn.to_csv(TABLES / "cldn4_atomic.tsv", sep="\t", index=False)
    cldn_neg = cldn[(cldn["rho"] < 0) & (cldn["n"] >= MIN_N)].sort_values(["immune_def", "p", "rho"])
    cldn_neg.to_csv(TABLES / "cldn4_negative_singles.tsv", sep="\t", index=False)
    primary_grid = _select_rows(cldn, PRIMARY_UNITS, "CLDN4")
    primary_grid.to_csv(TABLES / "primary_grid_cldn4.tsv", sep="\t", index=False)
    prim_pools = enumerate_named_grid(cldn, PRIMARY_UNITS, "primary_mixed")
    aligned = enumerate_aligned(cldn)
    pools = pd.concat([prim_pools, aligned], ignore_index=True)
    pools.to_csv(TABLES / "subset_pools.tsv", sep="\t", index=False)
    pools[pools["cldn4_negative"]].sort_values(["immune_family", "p_CLDN"]).to_csv(
        TABLES / "cldn4_negative_pools.tsv", sep="\t", index=False
    )
    extra_figures(cldn, pools, primary_grid)
    write_finding(cldn, pools, primary_grid)
    highlight_keys = [
        ("author/tnk/pct", "GSE131907+GSE205335"),
        ("marker/tnk/pct", "GSE253013+GSE291670"),
        ("marker/tnk/mean", "GSE253013+GSE291670"),
        ("primary_mixed", "GSE291670"),
        ("primary_mixed", "GSE291670+GSE253013"),
        ("primary_mixed", "GSE205335+GSE291670+GSE253013+GSE131907"),
        ("primary_mixed", "GSE207422+GSE291670+GSE253013"),
        ("primary_mixed", "GSE207422+GSE205335+GSE241934_IIT+GSE291670+GSE253013+GSE131907+GSE325414"),
        ("tls/cxcl13pos/mean", "GSE148071+GSE207422+GSE253013"),
        ("tls/cxcl13pos/mean", "GSE148071+GSE207422"),
        ("tls/cxcl13_mean/mean", "GSE207422+GSE253013"),
        ("author/b/mean", "GSE131907+GSE241934_IIT"),
    ]
    hrows = []
    for fam, cohorts in highlight_keys:
        hit = pools[(pools["family"] == fam) & (pools["cohorts"] == cohorts)]
        if not hit.empty:
            hrows.append(hit.iloc[0])
    if hrows:
        pd.DataFrame(hrows).to_csv(TABLES / "highlighted_combos.tsv", sep="\t", index=False)
    print("atomic", len(atomic_df), "CLDN4", len(cldn), "CLDN4-neg singles", len(cldn_neg))
    print("pools", len(pools), "CLDN4-neg pools", int(pools["cldn4_negative"].sum()))
    full = pools[(pools["family"] == "primary_mixed") & (pools["is_full_family"])]
    if not full.empty:
        r = full.iloc[0]
        print(f"FULL-POOL k={int(r.k)} N={int(r.n_patients)} CLDN4 {r.rho_CLDN:+.3f} p={r.p_CLDN:.3g}")
    print(
        cldn_neg[cldn_neg["immune_def"] == "tnk"][
            ["cohort", "malig_def", "score", "n", "rho", "p"]
        ]
        .head(15)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()

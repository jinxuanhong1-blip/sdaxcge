#!/usr/bin/env python3
"""Strict-malignant combinatorial search: CLDN4 and TACSTD2 vs T/NK.

ADDITIVE. Reuses existing scRNA patient tables. Keeps only author-malignant,
marker-malignant, or DRMref rows. Drops all-epithelial. Enumerates cohort
subsets until a cut has ρ ≤ −0.35, with honest n and p.
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
DATA = HERE / "data"
LEFT = DATA / "leftover"
OUT = HERE / "tables"
FIG = HERE / "figures"
MIN_N = 4
RHO_CUT = -0.35

KEEP_MALIG = {
    "author_malig",
    "marker_malig",
    "author_DRMref",
    "marker_malig_nsclc",
}
DROP_MALIG = {"author_epi", "marker_epi", "infercnv_p95", "infercnv_p90", "infercnv_gmm"}
TNK_DEFS = {"tnk", "tnk_umi"}  # T/NK only; CD8 / tcell / CXCL13 are not this search


def atomic(cohort, malig_def, tnk_def, score, gene, x, y, **note):
    rho, p, n = spearman(x, y)
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


def load_kept_atomic() -> list[dict]:
    """Score only author-malignant / marker-malignant / DRMref vs T/NK."""
    rows: list[dict] = []
    given = json.loads((HERE / "given_a3.json").read_text())

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
            "note": "A3 given DRMref; TACSTD2 not re-audited",
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
            "note": "A3 given secondary %pos; same 12-patient table",
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
            note="CLDN4 from given A3 DRMref table; not on the A3 slide",
        )
    )

    fab = pd.read_csv(DATA / "GSE207422_marker_fable.tsv", sep="\t")
    post = fab[fab["timing"] == "post"].copy()
    post_mal = post[post["n_malignant"] >= 10]
    for score, tcol, ccol in (
        ("mean", "malig_TACSTD2_mean", "malig_CLDN4_mean"),
        ("pct", "malig_TACSTD2_pct", "malig_CLDN4_pct"),
    ):
        rows.append(atomic("GSE207422", "marker_malig", "tnk", score, "TACSTD2", post_mal[tcol], post_mal["tnk_frac"], note="post, n_mal>=10"))
        rows.append(atomic("GSE207422", "marker_malig", "tnk", score, "CLDN4", post_mal[ccol], post_mal["tnk_frac"], note="post, n_mal>=10"))

    ns = pd.read_csv(DATA / "GSE207422_marker_nsclc.tsv", sep="\t")
    ns_post = ns[ns["timing"].astype(str).str.lower() == "post"].copy()
    ns_mal = ns_post[pd.to_numeric(ns_post["n_malignant_like"], errors="coerce").fillna(0) >= 10]
    for score, tcol, ccol in (
        ("mean", "mal_TACSTD2_mean", "mal_CLDN4_mean"),
        ("pct", "mal_TACSTD2_pct_pos", "mal_CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "tnk", score, "TACSTD2", ns_mal[tcol], ns_mal["frac_tnk"], note="post, n_mal_like>=10"))
        rows.append(atomic("GSE207422", "marker_malig_nsclc", "tnk", score, "CLDN4", ns_mal[ccol], ns_mal["frac_tnk"], note="post, n_mal_like>=10"))

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    nsclc = d[d["cancer_subtype"].isin(["ADC", "SQ"])]
    for cohort, df in (("GSE205335", d), ("GSE205335_NSCLC", nsclc)):
        for score, tcol, ccol in (
            ("mean", "mal_TACSTD2_mean", "mal_CLDN4_mean"),
            ("pct", "mal_TACSTD2_pct_pos", "mal_CLDN4_pct_pos"),
        ):
            rows.append(atomic(cohort, "author_malig", "tnk", score, "TACSTD2", df[tcol], df["frac_tnk"], note="author malignant"))
            rows.append(atomic(cohort, "author_malig", "tnk", score, "CLDN4", df[ccol], df["frac_tnk"], note="author malignant"))

    d = pd.read_csv(DATA / "GSE291670_patients.tsv", sep="\t")
    for tnk, y in (("tnk", "frac_lineage_tnk"), ("tnk_umi", "frac_umi_tnk")):
        for score, tcol, ccol in (
            ("mean", "mal_TACSTD2_mean_log1p_cp10k", "mal_CLDN4_mean_log1p_cp10k"),
            ("pct", "mal_TACSTD2_pct_pos", "mal_CLDN4_pct_pos"),
        ):
            rows.append(atomic("GSE291670", "marker_malig", tnk, score, "TACSTD2", d[tcol], d[y], note="marker malignant"))
            rows.append(atomic("GSE291670", "marker_malig", tnk, score, "CLDN4", d[ccol], d[y], note="marker malignant"))

    d = pd.read_csv(DATA / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    for score, tcol, ccol in (
        ("mean", "TACSTD2_mean_log1p", "CLDN4_mean_log1p"),
        ("mean_cp10k", "TACSTD2_mean_log1p_cp10k", "CLDN4_mean_log1p_cp10k"),
        ("pct", "TACSTD2_pct_pos", "CLDN4_pct_pos"),
    ):
        rows.append(atomic("GSE253013", "marker_malig", "tnk", score, "TACSTD2", t[tcol], t["tnk_fraction"], note="tumor malig-like"))
        rows.append(atomic("GSE253013", "marker_malig", "tnk", score, "CLDN4", t[ccol], t["tnk_fraction"], note="tumor malig-like"))

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    tlung = d[d["origin"] == "tLung"]
    for cohort, df in (("GSE131907", tumor), ("GSE131907_tLung", tlung)):
        mal = df[df["n_malignant"] >= 20]
        if len(mal) < 3:
            continue
        for score, tcol, ccol in (
            ("mean", "mal_TACSTD2_mean", "mal_CLDN4_mean"),
            ("pct", "mal_TACSTD2_pct", "mal_CLDN4_pct"),
        ):
            rows.append(atomic(cohort, "author_malig", "tnk", score, "TACSTD2", mal[tcol], mal["frac_tnk"], note="author Malignant cells subtype, n_mal>=20"))
            rows.append(atomic(cohort, "author_malig", "tnk", score, "CLDN4", mal[ccol], mal["frac_tnk"], note="author Malignant cells subtype, n_mal>=20"))

    d = pd.read_csv(DATA / "GSE325414_donors.csv")
    rows.append(atomic("GSE325414", "author_malig", "tnk", "mean", "TACSTD2", d["malignant_TACSTD2_mean"], d["frac_TNK"], note="author malignant, donor"))
    rows.append(atomic("GSE325414", "author_malig", "tnk", "mean", "CLDN4", d["malignant_CLDN4_mean"], d["frac_TNK"], note="author malignant, donor"))
    return rows


def drop_audit() -> pd.DataFrame:
    """Document all-epithelial / empty-malignant rows that are not scored."""
    recs = [
        {
            "source": "GSE241934_IIT",
            "reason": "all-epithelial",
            "detail": "author residual Epi only; no author-malignant / marker-malignant / DRMref column",
            "n_rows": 11,
        },
        {
            "source": "GSE241934_Real",
            "reason": "all-epithelial",
            "detail": "author residual Epi only",
            "n_rows": 24,
        },
        {
            "source": "GSE131907 author_epi",
            "reason": "all-epithelial",
            "detail": "author epithelium columns dropped; author_malig n_mal>=20 is kept",
            "n_rows": 36,
        },
        {
            "source": "GSE253013 author_epi / marker_epi",
            "reason": "all-epithelial",
            "detail": "author Epithelial and tumor epithelium columns dropped; marker_malig kept",
            "n_rows": 9,
        },
        {
            "source": "GSE207422 marker_epi",
            "reason": "all-epithelial",
            "detail": "post all-epithelial extra dropped; marker_malig and DRMref kept",
            "n_rows": 12,
        },
        {
            "source": "GSE291670 marker_epi",
            "reason": "all-epithelial",
            "detail": "epi_* columns dropped; marker_malig kept",
            "n_rows": 6,
        },
        {
            "source": "E-MTAB-13526",
            "reason": "all-epithelial",
            "detail": "author Cell types epithelium only; no malignant split",
            "n_rows": 13,
        },
        {
            "source": "GSE267108 malignant",
            "reason": "malignant n<4",
            "detail": "leftover stats: mal_TACSTD2 n=1; epithelial-only would be dropped anyway",
            "n_rows": 1,
        },
        {
            "source": "GSE274595 malignant",
            "reason": "malignant n<4",
            "detail": "leftover stats: mal_TACSTD2 n=0",
            "n_rows": 0,
        },
        {
            "source": "GSE207422 inferCNV-like",
            "reason": "not author/marker/DRMref",
            "detail": "wave2 inferCNV defs dropped; CLDN4 absent in that table",
            "n_rows": 12,
        },
    ]
    return pd.DataFrame(recs)


def pair_both_genes(atomic_df: pd.DataFrame) -> pd.DataFrame:
    keys = ["cohort", "malig_def", "tnk_def", "score"]
    t = atomic_df[atomic_df["gene"] == "TACSTD2"].set_index(keys)
    c = atomic_df[atomic_df["gene"] == "CLDN4"].set_index(keys)
    both = t.join(c, lsuffix="_TAC", rsuffix="_CLDN", how="inner")
    both = both.reset_index()
    both["n"] = both[["n_TAC", "n_CLDN"]].min(axis=1).astype(int)
    both["both_negative"] = (both["rho_TAC"] < 0) & (both["rho_CLDN"] < 0)
    both["tac_le_cut"] = both["rho_TAC"] <= RHO_CUT
    both["cldn_le_cut"] = both["rho_CLDN"] <= RHO_CUT
    both["both_le_cut"] = both["tac_le_cut"] & both["cldn_le_cut"]
    both["weaker_rho"] = both[["rho_TAC", "rho_CLDN"]].max(axis=1)
    return both


def _re_or_single(rhos, ns, ps):
    if len(rhos) == 1:
        return {
            "k": 1,
            "n_patients_total": int(ns[0]),
            "pooled_rho": float(rhos[0]),
            "p": float(ps[0]),
            "I2": 0.0,
            "ci95_rho": [float("nan"), float("nan")],
        }
    return random_effects_dl(rhos, ns)


def pool_rows(sub: pd.DataFrame, family: str, full_n: int) -> dict:
    rhos_t = sub["rho_TAC"].tolist()
    rhos_c = sub["rho_CLDN"].tolist()
    ns = sub["n"].astype(int).tolist()
    ps_t = sub["p_TAC"].tolist()
    ps_c = sub["p_CLDN"].tolist()
    re_t = _re_or_single(rhos_t, ns, ps_t)
    re_c = _re_or_single(rhos_c, ns, ps_c)
    st_t = stouffer(rhos_t, ps_t, ns)
    st_c = stouffer(rhos_c, ps_c, ns)
    rho_t = re_t.get("pooled_rho")
    rho_c = re_c.get("pooled_rho")
    return {
        "family": family,
        "k": int(len(sub)),
        "cohorts": "+".join(sub["cohort"].tolist()),
        "malig_defs": "+".join(sub["malig_def"].tolist()),
        "tnk_def": "+".join(sorted(set(sub["tnk_def"]))),
        "score": "+".join(sorted(set(sub["score"]))),
        "n_patients": int(sum(ns)),
        "rho_TAC": rho_t,
        "p_TAC": re_t.get("p"),
        "I2_TAC": re_t.get("I2"),
        "rho_CLDN": rho_c,
        "p_CLDN": re_c.get("p"),
        "I2_CLDN": re_c.get("I2"),
        "stouffer_p_TAC": st_t.get("p"),
        "stouffer_p_CLDN": st_c.get("p"),
        "both_negative": bool(rho_t is not None and rho_c is not None and rho_t < 0 and rho_c < 0),
        "tac_le_cut": bool(rho_t is not None and rho_t <= RHO_CUT),
        "cldn_le_cut": bool(rho_c is not None and rho_c <= RHO_CUT),
        "both_le_cut": bool(
            rho_t is not None and rho_c is not None and rho_t <= RHO_CUT and rho_c <= RHO_CUT
        ),
        "weaker_rho": max(rho_t, rho_c) if rho_t is not None and rho_c is not None else float("nan"),
        "is_full_family": int(len(sub) == full_n),
        "n_family_units": full_n,
    }


def _select(paired: pd.DataFrame, specs: list[tuple]) -> pd.DataFrame:
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


def enumerate_grid(paired: pd.DataFrame, specs: list[tuple], family: str) -> pd.DataFrame:
    grid = _select(paired, specs)
    units = list(range(len(specs)))
    out = []
    for k in range(1, len(units) + 1):
        for subset in itertools.combinations(units, k):
            sub = grid.iloc[list(subset)]
            out.append(pool_rows(sub, family, len(specs)))
    return pd.DataFrame(out)


def enumerate_aligned(paired: pd.DataFrame) -> pd.DataFrame:
    """Aligned families: one malignant bucket × T/NK × score, then cohort subsets."""
    malig_bucket = {
        "author_DRMref": "author_or_DRMref",
        "author_malig": "author_or_DRMref",
        "marker_malig": "marker_malig",
        "marker_malig_nsclc": "marker_malig",
    }
    score_bucket = {"mean": "mean", "mean_cp10k": "mean", "pct": "pct"}
    work = paired.copy()
    work = work[work["malig_def"].isin(KEEP_MALIG)]
    work = work[work["tnk_def"].isin(TNK_DEFS)]
    work = work[work["n"] >= MIN_N]
    # one row per base cohort (drop NSCLC / tLung variants from aligned families)
    work = work[~work["cohort"].str.contains("_NSCLC|_tLung|_pooled", regex=True)]
    work["malig_b"] = work["malig_def"].map(malig_bucket)
    work["score_b"] = work["score"].map(score_bucket)
    work = work.dropna(subset=["malig_b", "score_b"])

    pref = {
        "author_or_DRMref": ["author_malig", "author_DRMref"],
        "marker_malig": ["marker_malig", "marker_malig_nsclc"],
    }
    out = []
    for (mb, tb, sb), g in work.groupby(["malig_b", "tnk_def", "score_b"], dropna=False):
        chosen = []
        for _unit, ug in g.groupby("cohort"):
            order = pref.get(mb, list(ug["malig_def"].unique()))
            ug = ug.copy()
            ug["_rank"] = ug["malig_def"].apply(lambda x: order.index(x) if x in order else 99)
            chosen.append(ug.sort_values(["_rank", "n"], ascending=[True, False]).head(1))
        grid = pd.concat(chosen, ignore_index=True)
        family = f"{mb}/{tb}/{sb}"
        units = list(range(len(grid)))
        for k in range(1, len(units) + 1):
            for subset in itertools.combinations(units, k):
                sub = grid.iloc[list(subset)]
                out.append(pool_rows(sub, family, len(grid)))
    return pd.DataFrame(out)


# Strict 6-unit primary grid: no all-epithelial. Mean vs T/NK.
PRIMARY_UNITS = [
    ("GSE207422", "author_DRMref", "tnk", "mean"),
    ("GSE205335", "author_malig", "tnk", "mean"),
    ("GSE291670", "marker_malig", "tnk", "mean"),
    ("GSE253013", "marker_malig", "tnk", "mean"),
    ("GSE131907", "author_malig", "tnk", "mean"),
    ("GSE325414", "author_malig", "tnk", "mean"),
]
PRIMARY_NSCLC_SWAP = [
    ("GSE207422", "author_DRMref", "tnk", "mean"),
    ("GSE205335_NSCLC", "author_malig", "tnk", "mean"),
    ("GSE291670", "marker_malig", "tnk", "mean"),
    ("GSE253013", "marker_malig", "tnk", "mean"),
    ("GSE131907", "author_malig", "tnk", "mean"),
    ("GSE325414", "author_malig", "tnk", "mean"),
]
PRIMARY_PCT = [
    ("GSE207422", "author_DRMref", "tnk", "pct"),
    ("GSE205335", "author_malig", "tnk", "pct"),
    ("GSE291670", "marker_malig", "tnk", "pct"),
    ("GSE253013", "marker_malig", "tnk", "pct"),
    ("GSE131907", "author_malig", "tnk", "pct"),
]


def forest_pair(sub: pd.DataFrame, title: str, path: Path) -> None:
    sub = sub.sort_values("n", ascending=False).reset_index(drop=True)
    labels = [f"{r.cohort} n={int(r.n)}" for r in sub.itertuples()]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 1.15 + 0.40 * max(len(sub), 1)), sharey=True)
    y = np.arange(len(sub))
    for ax, gene, rho_col, p_col, color in (
        (axes[0], "TACSTD2", "rho_TAC", "p_TAC", "#1f4e79"),
        (axes[1], "CLDN4", "rho_CLDN", "p_CLDN", "#7a2d0b"),
    ):
        for i, r in enumerate(sub.itertuples()):
            rho = float(getattr(r, rho_col))
            p = float(getattr(r, p_col))
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
        ax.axvline(RHO_CUT, color="#c44", lw=0.7, ls=":")
        ax.set_xlim(-1.05, 1.05)
        ax.set_title(gene, fontsize=10)
        ax.set_xlabel("Spearman ρ vs T/NK")
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
    ax.scatter(list(x), list(y), c="#1f4e79", s=36, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt_rho(rho, p, i2=None) -> str:
    if rho is None or not np.isfinite(rho):
        return "NA"
    s = f"{rho:+.3f} (p={p:.3g}"
    if i2 is not None and np.isfinite(i2):
        s += f", I²={i2:.0f}%"
    return s + ")"


def write_finding(
    primary_grid: pd.DataFrame,
    primary_pools: pd.DataFrame,
    nsclc_pools: pd.DataFrame,
    pct_pools: pd.DataFrame,
    aligned: pd.DataFrame,
    paired: pd.DataFrame,
    first_cut: dict,
    dropped: pd.DataFrame,
) -> None:
    full = primary_pools[primary_pools["is_full_family"] == 1].iloc[0]
    nsclc_full = nsclc_pools[nsclc_pools["is_full_family"] == 1].iloc[0]

    both_cut_pools = pd.concat(
        [
            primary_pools[primary_pools["both_le_cut"]],
            aligned[aligned["both_le_cut"]],
            pct_pools[pct_pools["both_le_cut"]] if not pct_pools.empty else pct_pools,
        ],
        ignore_index=True,
    )
    # One row per (cohort set, score). Prefer primary-grid family labels.
    both_cut_pools["_fam_rank"] = both_cut_pools["family"].map(
        lambda s: 0 if str(s).startswith("strict_primary") else 1
    )
    both_cut_pools["_set"] = both_cut_pools["cohorts"].map(lambda s: "+".join(sorted(str(s).split("+"))))
    both_cut_pools = both_cut_pools.sort_values(["_fam_rank", "k", "n_patients", "weaker_rho"])
    both_cut_pools = both_cut_pools.drop_duplicates(subset=["_set", "score", "n_patients"], keep="first")
    both_cut_pools = both_cut_pools.sort_values(["k", "n_patients", "weaker_rho"], ascending=[True, False, True])

    singles_cut = paired[(paired["both_le_cut"]) & (paired["n"] >= MIN_N)].sort_values("weaker_rho")
    singles_neg = paired[(paired["both_negative"]) & (paired["n"] >= MIN_N)].sort_values("rho_TAC")

    prim_cut = primary_pools[primary_pools["both_le_cut"]].sort_values(["k", "p_TAC"])
    prim_neg = primary_pools[(primary_pools["both_negative"]) & (primary_pools["k"] >= 2)]
    prim_neg = prim_neg.sort_values("p_TAC")

    aligned_cut = aligned[aligned["both_le_cut"]].sort_values(["k", "weaker_rho"])

    lines = []
    lines.append("# Strict malignant T/NK recut: combinatorial CLDN4 and TACSTD2")
    lines.append("")
    lines.append("ADDITIVE. Patient is the unit. Existing scRNA patient tables are reused;")
    lines.append("GEO is not re-downloaded. **Keep only author-malignant, marker-malignant,")
    lines.append("or DRMref.** All-epithelial rows are dropped (GSE241934 residual Epi,")
    lines.append("GSE131907/GSE253013/GSE207422/GSE291670 epithelium columns, E-MTAB-13526).")
    lines.append("Leftover GSE267108/GSE274595 malignant splits have n<4 and cannot enter.")
    lines.append("inferCNV-like is not author/marker/DRMref and is dropped (CLDN4 also absent).")
    lines.append("")
    lines.append("GSE207422 A3 TACSTD2 is taken as given (n=12, mean ρ=−0.490, p=0.106).")
    lines.append("CLDN4 is from the same given DRMref table. p-values are descriptive")
    lines.append("(many subsets). The ρ ≤ −0.35 rule is **both genes**, not a p-value gate.")
    lines.append("")
    lines.append("Locked search order: (1) the 6 primary-grid mean-vs-T/NK members,")
    lines.append("(2) other mean singles, (3) %pos singles, (4) k=2 aligned mean,")
    lines.append("(5) larger subsets. First hit is reported with its honest n and p;")
    lines.append("every other both-≤−0.35 cut is in `tables/combo_table.tsv`.")
    lines.append("")
    lines.append("## First cut with both genes ρ ≤ −0.35")
    lines.append("")
    lines.append(
        f"- **{first_cut['label']}** · k={first_cut['k']} · N={first_cut['n']} · "
        f"TACSTD2 ρ={first_cut['rho_TAC']:+.3f} p={first_cut['p_TAC']:.3g} · "
        f"CLDN4 ρ={first_cut['rho_CLDN']:+.3f} p={first_cut['p_CLDN']:.3g}"
    )
    lines.append("")
    lines.append("GSE253013 marker_malig / tnk / mean is a near-miss: TACSTD2 ρ=−0.717")
    lines.append("p=0.0298, CLDN4 ρ=−0.333 p=0.381 (CLDN4 just above −0.35). The same")
    lines.append("two marker-malignant cohorts pooled (mean) are N=15, TACSTD2 ρ=−0.666")
    lines.append("p=0.016 I²=0%, CLDN4 ρ=−0.582 p=0.102 I²=29%.")
    lines.append("")
    lines.append("## Strict full-pool row (one row, not the answer)")
    lines.append("")
    lines.append("6-unit grid after dropping all-epithelial. Mean log1p vs T/NK.")
    lines.append("")
    lines.append("| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |")
    lines.append("|---:|---:|---|---|---|---|")
    lines.append(
        f"| {int(full.k)} | {int(full.n_patients)} | {full.family} | "
        f"{fmt_rho(full.rho_TAC, full.p_TAC, full.I2_TAC)} | "
        f"{fmt_rho(full.rho_CLDN, full.p_CLDN, full.I2_CLDN)} | {full.cohorts} |"
    )
    lines.append(
        f"| {int(nsclc_full.k)} | {int(nsclc_full.n_patients)} | {nsclc_full.family} | "
        f"{fmt_rho(nsclc_full.rho_TAC, nsclc_full.p_TAC, nsclc_full.I2_TAC)} | "
        f"{fmt_rho(nsclc_full.rho_CLDN, nsclc_full.p_CLDN, nsclc_full.I2_CLDN)} | {nsclc_full.cohorts} |"
    )
    lines.append("")
    lines.append("## Strict primary-grid members")
    lines.append("")
    lines.append("| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) | both ρ<0 | both ρ≤−0.35 |")
    lines.append("|---|---|---|---|---:|---|---|---|---|")
    for r in primary_grid.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.tnk_def} | {r.score} | {int(r.n)} | "
            f"{r.rho_TAC:+.3f} ({r.p_TAC:.3g}) | {r.rho_CLDN:+.3f} ({r.p_CLDN:.3g}) | "
            f"{'yes' if r.both_negative else 'no'} | {'yes' if r.both_le_cut else 'no'} |"
        )
    lines.append("")
    lines.append("## Cuts with both genes ρ ≤ −0.35 (honest n and p)")
    lines.append("")
    lines.append("Singles first, then pooled subsets. Full combo table:")
    lines.append("`tables/combo_table.tsv`.")
    lines.append("")
    lines.append("| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |")
    lines.append("|---:|---:|---|---|---|---|")
    shown = both_cut_pools.head(20)
    for r in shown.itertuples():
        lines.append(
            f"| {int(r.k)} | {int(r.n_patients)} | {r.family} | "
            f"{fmt_rho(r.rho_TAC, r.p_TAC, r.I2_TAC)} | "
            f"{fmt_rho(r.rho_CLDN, r.p_CLDN, r.I2_CLDN)} | {r.cohorts} |"
        )
    lines.append("")
    lines.append(f"Count: {int(both_cut_pools.shape[0])} enumerated cuts with both genes ρ ≤ −0.35.")
    lines.append("")
    lines.append("## Single-cohort strict T/NK (n≥4, both genes present)")
    lines.append("")
    lines.append("| cohort | malig | T/NK | score | n | TACSTD2 ρ (p) | CLDN4 ρ (p) | both ρ≤−0.35 |")
    lines.append("|---|---|---|---|---:|---|---|---|")
    for r in singles_neg.itertuples():
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.tnk_def} | {r.score} | {int(r.n)} | "
            f"{r.rho_TAC:+.3f} ({r.p_TAC:.3g}) | {r.rho_CLDN:+.3f} ({r.p_CLDN:.3g}) | "
            f"{'yes' if r.both_le_cut else 'no'} |"
        )
    lines.append("")
    lines.append(
        f"Both-negative singles: {int(singles_neg.shape[0])} of "
        f"{int(((paired['n'] >= MIN_N)).sum())} paired strict T/NK combos with n≥4. "
        f"Both ρ≤−0.35: {int(singles_cut.shape[0])}."
    )
    lines.append("")
    lines.append("## Primary-grid subset pools with both genes ρ < 0 (k≥2)")
    lines.append("")
    lines.append("All 2⁶−1 = 63 nonempty subsets of the 6 strict units were enumerated.")
    lines.append("Rows below: both-negative, N≥15, lowest TACSTD2 RE p first (top 15).")
    lines.append("")
    lines.append("| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | both ρ≤−0.35 | cohorts |")
    lines.append("|---:|---:|---|---|---|---|---|")
    showp = prim_neg[prim_neg["n_patients"] >= 15].head(15)
    for r in showp.itertuples():
        lines.append(
            f"| {int(r.k)} | {int(r.n_patients)} | {r.family} | "
            f"{fmt_rho(r.rho_TAC, r.p_TAC, r.I2_TAC)} | "
            f"{fmt_rho(r.rho_CLDN, r.p_CLDN, r.I2_CLDN)} | "
            f"{'yes' if r.both_le_cut else 'no'} | {r.cohorts} |"
        )
    lines.append("")
    lines.append(
        f"Primary-grid both-negative k≥2 subsets: "
        f"{int(prim_neg.shape[0])} of {int((primary_pools['k'] >= 2).sum())} enumerated. "
        f"Both ρ≤−0.35: {int(prim_cut.shape[0])}."
    )
    lines.append("")
    lines.append("## Aligned-family cuts with both genes ρ ≤ −0.35")
    lines.append("")
    lines.append("| k | N | family | TACSTD2 ρ (p, I²) | CLDN4 ρ (p, I²) | cohorts |")
    lines.append("|---:|---:|---|---|---|---|")
    for r in aligned_cut.head(15).itertuples():
        lines.append(
            f"| {int(r.k)} | {int(r.n_patients)} | {r.family} | "
            f"{fmt_rho(r.rho_TAC, r.p_TAC, r.I2_TAC)} | "
            f"{fmt_rho(r.rho_CLDN, r.p_CLDN, r.I2_CLDN)} | {r.cohorts} |"
        )
    lines.append("")
    lines.append("## Dropped rows")
    lines.append("")
    lines.append("| source | reason | n | detail |")
    lines.append("|---|---|---:|---|")
    for r in dropped.itertuples():
        lines.append(f"| {r.source} | {r.reason} | {int(r.n_rows)} | {r.detail} |")
    lines.append("")
    lines.append("## Extra figures")
    lines.append("")
    lines.append("See `figures/`: strict primary forest, NSCLC swap, both-≤−0.35 singles,")
    lines.append("marker pair, trio, subset bars, and patient scatters.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/strict_malig_tnk/combinatorial_search.py")
    lines.append("```")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def extra_figures(paired: pd.DataFrame, primary_grid: pd.DataFrame, nsclc_grid: pd.DataFrame, pools: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)

    forest_pair(
        primary_grid,
        "Strict 6-unit primary grid (author-malig / marker-malig / DRMref; mean vs T/NK)",
        FIG / "strict_primary_members_forest.png",
    )
    forest_pair(
        nsclc_grid,
        "Strict 6-unit grid, GSE205335 → ADC+SQ only",
        FIG / "strict_nsclc_swap_members_forest.png",
    )

    prim_neg = primary_grid[primary_grid["both_negative"]]
    if len(prim_neg) >= 2:
        forest_pair(
            prim_neg,
            "Strict primary members with TACSTD2 ρ<0 and CLDN4 ρ<0",
            FIG / "strict_primary_both_negative_forest.png",
        )

    singles_cut = paired[(paired["both_le_cut"]) & (paired["n"] >= MIN_N) & (paired["tnk_def"] == "tnk")]
    if not singles_cut.empty:
        forest_pair(
            singles_cut.sort_values("weaker_rho"),
            "Single-cohort strict T/NK cuts with both genes ρ ≤ −0.35",
            FIG / "both_le_cut_singles_forest.png",
        )

    trio = primary_grid[primary_grid["cohort"].isin(["GSE207422", "GSE291670", "GSE253013"])]
    if len(trio) == 3:
        forest_pair(
            trio,
            "Strict trio GSE207422 DRMref + GSE291670 + GSE253013 (mean vs T/NK)",
            FIG / "strict_trio_GSE207422_GSE291670_GSE253013.png",
        )

    pair = paired[
        (paired["cohort"].isin(["GSE253013", "GSE291670"]))
        & (paired["malig_def"] == "marker_malig")
        & (paired["tnk_def"] == "tnk")
        & (paired["score"] == "pct")
    ]
    if len(pair) == 2:
        forest_pair(
            pair,
            "Marker-malignant %pos pair GSE253013 + GSE291670 vs T/NK",
            FIG / "marker_malig_pct_GSE253013_GSE291670.png",
        )
    pair_m = paired[
        (paired["cohort"].isin(["GSE253013", "GSE291670"]))
        & (paired["malig_def"] == "marker_malig")
        & (paired["tnk_def"] == "tnk")
        & (paired["score"] == "mean")
    ]
    if len(pair_m) == 2:
        forest_pair(
            pair_m,
            "Marker-malignant mean pair GSE253013 + GSE291670 vs T/NK",
            FIG / "marker_malig_mean_GSE253013_GSE291670.png",
        )

    cuts = pools[pools["both_le_cut"] & (pools["k"] >= 2)].copy()
    if not cuts.empty:
        cuts = cuts.sort_values(["k", "weaker_rho"]).drop_duplicates(subset=["cohorts", "family"]).head(8)
        fig, ax = plt.subplots(figsize=(9.2, 4.4))
        y = np.arange(len(cuts))
        ax.barh(y - 0.15, cuts["rho_TAC"], height=0.3, color="#1f4e79", label="TACSTD2")
        ax.barh(y + 0.15, cuts["rho_CLDN"], height=0.3, color="#b86b2a", label="CLDN4")
        ax.axvline(0, color="#888", lw=0.8)
        ax.axvline(RHO_CUT, color="#c44", lw=0.8, ls=":")
        labels = [
            f"k={int(r.k)} N={int(r.n_patients)} {r.family}"
            for r in cuts.itertuples()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Pooled Spearman ρ vs T/NK")
        ax.set_title("Subset pools with both genes ρ ≤ −0.35 (dotted line)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIG / "both_le_cut_subset_bars.png", dpi=150)
        fig.savefig(FIG / "both_le_cut_subset_bars.pdf")
        plt.close(fig)

    # Patient scatters
    d = pd.read_csv(DATA / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    scatter_combo(
        FIG / "scatter_GSE253013_marker_malig_mean.png",
        t["tnk_fraction"],
        t["TACSTD2_mean_log1p"],
        "T/NK fraction",
        "Malignant-like TACSTD2 mean log1p",
        "GSE253013 marker_malig / tnk / mean (n=9)",
    )
    scatter_combo(
        FIG / "scatter_GSE253013_marker_malig_mean_CLDN4.png",
        t["tnk_fraction"],
        t["CLDN4_mean_log1p"],
        "T/NK fraction",
        "Malignant-like CLDN4 mean log1p",
        "GSE253013 marker_malig / tnk / mean (n=9)",
    )
    scatter_combo(
        FIG / "scatter_GSE253013_marker_malig_pct.png",
        t["tnk_fraction"],
        t["TACSTD2_pct_pos"],
        "T/NK fraction",
        "Malignant-like TACSTD2 %pos",
        "GSE253013 marker_malig / tnk / pct (n=9)",
    )
    scatter_combo(
        FIG / "scatter_GSE253013_marker_malig_pct_CLDN4.png",
        t["tnk_fraction"],
        t["CLDN4_pct_pos"],
        "T/NK fraction",
        "Malignant-like CLDN4 %pos",
        "GSE253013 marker_malig / tnk / pct (n=9)",
    )

    d = pd.read_csv(DATA / "GSE291670_patients.tsv", sep="\t")
    scatter_combo(
        FIG / "scatter_GSE291670_marker_malig_mean.png",
        d["frac_lineage_tnk"],
        d["mal_TACSTD2_mean_log1p_cp10k"],
        "Lineage T/NK fraction",
        "Malignant TACSTD2 mean log1p(CP10k)",
        "GSE291670 marker_malig / tnk / mean (n=6)",
    )
    scatter_combo(
        FIG / "scatter_GSE291670_marker_malig_mean_CLDN4.png",
        d["frac_lineage_tnk"],
        d["mal_CLDN4_mean_log1p_cp10k"],
        "Lineage T/NK fraction",
        "Malignant CLDN4 mean log1p(CP10k)",
        "GSE291670 marker_malig / tnk / mean (n=6)",
    )

    a3 = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    scatter_combo(
        FIG / "scatter_GSE207422_DRMref_TACSTD2.png",
        a3["frac_tnk"],
        a3["malig_TACSTD2_mean"],
        "T/NK fraction (DRMref)",
        "Malignant TACSTD2 mean log1p(CP10k)",
        "GSE207422 A3 given · author_DRMref / tnk / mean (n=12)",
    )
    scatter_combo(
        FIG / "scatter_GSE207422_DRMref_CLDN4.png",
        a3["frac_tnk"],
        a3["malig_CLDN4_mean"],
        "T/NK fraction (DRMref)",
        "Malignant CLDN4 mean log1p(CP10k)",
        "GSE207422 A3 given table · CLDN4 (not on the A3 slide)",
    )

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    mal = tumor[tumor["n_malignant"] >= 20]
    scatter_combo(
        FIG / "scatter_GSE131907_author_malig_pct_CLDN4.png",
        mal["frac_tnk"],
        mal["mal_CLDN4_pct"],
        "T/NK fraction",
        "Author-malignant CLDN4 %pos",
        f"GSE131907 author_malig / tnk / pct (n={len(mal)})",
    )
    scatter_combo(
        FIG / "scatter_GSE131907_author_malig_mean_TACSTD2.png",
        mal["frac_tnk"],
        mal["mal_TACSTD2_mean"],
        "T/NK fraction",
        "Author-malignant TACSTD2 mean",
        f"GSE131907 author_malig / tnk / mean (n={len(mal)})",
    )

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    scatter_combo(
        FIG / "scatter_GSE205335_author_malig_pct_CLDN4.png",
        d["frac_tnk"],
        d["mal_CLDN4_pct_pos"],
        "T/NK fraction",
        "Author-malignant CLDN4 %pos",
        f"GSE205335 author_malig / tnk / pct (n={len(d)})",
    )


def _cut_record(label, k, n, rho_t, p_t, rho_c, p_c, stage) -> dict:
    return {
        "label": label,
        "k": int(k),
        "n": int(n),
        "rho_TAC": float(rho_t),
        "p_TAC": float(p_t),
        "rho_CLDN": float(rho_c),
        "p_CLDN": float(p_c),
        "stage": stage,
    }


def first_cut_search(paired: pd.DataFrame, aligned: pd.DataFrame, primary_pools: pd.DataFrame) -> dict:
    """Locked order: primary-grid mean → other mean → %pos → k=2 mean → larger."""
    # 1) The 6 locked primary-grid members (mean vs T/NK).
    prim = paired.merge(
        pd.DataFrame(PRIMARY_UNITS, columns=["cohort", "malig_def", "tnk_def", "score"]),
        on=["cohort", "malig_def", "tnk_def", "score"],
        how="inner",
    )
    hit = prim[prim["both_le_cut"]]
    if not hit.empty:
        r = hit.iloc[0]
        return _cut_record(
            f"{r.cohort} · {r.malig_def}/{r.tnk_def}/{r.score} (primary-grid member)",
            1,
            r.n,
            r.rho_TAC,
            r.p_TAC,
            r.rho_CLDN,
            r.p_CLDN,
            "primary_grid_mean",
        )
    # 2) Any other mean vs T/NK single.
    mean_s = paired[
        (paired["n"] >= MIN_N)
        & (paired["tnk_def"] == "tnk")
        & (paired["score"].isin(["mean", "mean_cp10k"]))
        & (paired["both_le_cut"])
    ].sort_values(["n", "weaker_rho"], ascending=[False, True])
    if not mean_s.empty:
        r = mean_s.iloc[0]
        return _cut_record(
            f"{r.cohort} · {r.malig_def}/{r.tnk_def}/{r.score}",
            1,
            r.n,
            r.rho_TAC,
            r.p_TAC,
            r.rho_CLDN,
            r.p_CLDN,
            "single_mean",
        )
    # 3) %pos vs T/NK singles.
    pct_s = paired[
        (paired["n"] >= MIN_N) & (paired["tnk_def"] == "tnk") & (paired["score"] == "pct") & (paired["both_le_cut"])
    ].sort_values(["n", "weaker_rho"], ascending=[False, True])
    if not pct_s.empty:
        r = pct_s.iloc[0]
        return _cut_record(
            f"{r.cohort} · {r.malig_def}/{r.tnk_def}/{r.score}",
            1,
            r.n,
            r.rho_TAC,
            r.p_TAC,
            r.rho_CLDN,
            r.p_CLDN,
            "single_pct",
        )
    # 4) k=2 aligned mean, then any k=2.
    k2_mean = aligned[(aligned["k"] == 2) & (aligned["both_le_cut"]) & (aligned["family"].str.endswith("/mean"))]
    k2 = k2_mean if not k2_mean.empty else aligned[(aligned["k"] == 2) & (aligned["both_le_cut"])]
    k2 = k2.sort_values("n_patients", ascending=False)
    if not k2.empty:
        r = k2.iloc[0]
        return _cut_record(f"{r.cohorts} · {r.family}", 2, r.n_patients, r.rho_TAC, r.p_TAC, r.rho_CLDN, r.p_CLDN, "k2")
    later = primary_pools[primary_pools["both_le_cut"]].sort_values(["k", "n_patients"])
    if not later.empty:
        r = later.iloc[0]
        return _cut_record(
            f"{r.cohorts} · {r.family}",
            r.k,
            r.n_patients,
            r.rho_TAC,
            r.p_TAC,
            r.rho_CLDN,
            r.p_CLDN,
            "primary_subset",
        )
    return _cut_record("NO CUT", 0, 0, float("nan"), float("nan"), float("nan"), float("nan"), "none")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    raw = pd.DataFrame(load_kept_atomic())
    assert set(raw["malig_def"]).issubset(KEEP_MALIG), set(raw["malig_def"])
    assert set(raw["tnk_def"]).issubset(TNK_DEFS), set(raw["tnk_def"])
    assert not set(raw["malig_def"]) & DROP_MALIG

    atomic_df = raw[raw["n"] >= MIN_N].copy()
    atomic_df.to_csv(OUT / "atomic_effects.tsv", sep="\t", index=False)

    paired = pair_both_genes(atomic_df)
    paired.to_csv(OUT / "paired_combos.tsv", sep="\t", index=False)

    dropped = drop_audit()
    dropped.to_csv(OUT / "dropped_rows.tsv", sep="\t", index=False)

    primary_grid = _select(paired, PRIMARY_UNITS)
    nsclc_grid = _select(paired, PRIMARY_NSCLC_SWAP)
    primary_pools = enumerate_grid(paired, PRIMARY_UNITS, "strict_primary/tnk/mean")
    nsclc_pools = enumerate_grid(paired, PRIMARY_NSCLC_SWAP, "strict_nsclc_swap/tnk/mean")
    pct_specs = [s for s in PRIMARY_PCT if s in set(map(tuple, paired[["cohort", "malig_def", "tnk_def", "score"]].to_numpy()))]
    pct_pools = enumerate_grid(paired, pct_specs, "strict_primary/tnk/pct") if pct_specs else pd.DataFrame()
    aligned = enumerate_aligned(paired)

    primary_grid.to_csv(OUT / "primary_grid.tsv", sep="\t", index=False)
    primary_pools.to_csv(OUT / "primary_subset_pools.tsv", sep="\t", index=False)
    nsclc_pools.to_csv(OUT / "nsclc_swap_pools.tsv", sep="\t", index=False)
    if not pct_pools.empty:
        pct_pools.to_csv(OUT / "pct_subset_pools.tsv", sep="\t", index=False)
    aligned.to_csv(OUT / "aligned_family_pools.tsv", sep="\t", index=False)

    combo = pd.concat(
        [primary_pools, nsclc_pools, pct_pools, aligned],
        ignore_index=True,
        sort=False,
    )
    # also add paired singles as k=1 combo rows
    single_rows = []
    for r in paired[paired["n"] >= MIN_N].itertuples():
        single_rows.append(
            {
                "family": f"single/{r.malig_def}/{r.tnk_def}/{r.score}",
                "k": 1,
                "cohorts": r.cohort,
                "malig_defs": r.malig_def,
                "tnk_def": r.tnk_def,
                "score": r.score,
                "n_patients": int(r.n),
                "rho_TAC": r.rho_TAC,
                "p_TAC": r.p_TAC,
                "I2_TAC": 0.0,
                "rho_CLDN": r.rho_CLDN,
                "p_CLDN": r.p_CLDN,
                "I2_CLDN": 0.0,
                "stouffer_p_TAC": r.p_TAC,
                "stouffer_p_CLDN": r.p_CLDN,
                "both_negative": bool(r.both_negative),
                "tac_le_cut": bool(r.tac_le_cut),
                "cldn_le_cut": bool(r.cldn_le_cut),
                "both_le_cut": bool(r.both_le_cut),
                "weaker_rho": r.weaker_rho,
                "is_full_family": 0,
                "n_family_units": 1,
            }
        )
    combo = pd.concat([pd.DataFrame(single_rows), combo], ignore_index=True)
    combo = combo.drop_duplicates(subset=["family", "cohorts", "n_patients", "score"], keep="first")
    combo = combo.sort_values(
        ["both_le_cut", "both_negative", "k", "n_patients"],
        ascending=[False, False, True, False],
    )
    combo.to_csv(OUT / "combo_table.tsv", sep="\t", index=False)
    combo[combo["both_le_cut"]].to_csv(OUT / "cuts_both_le_neg035.tsv", sep="\t", index=False)

    first = first_cut_search(paired, aligned, primary_pools)
    pd.DataFrame([first]).to_csv(OUT / "first_cut.tsv", sep="\t", index=False)

    extra_figures(paired, primary_grid, nsclc_grid, combo)
    write_finding(
        primary_grid,
        primary_pools,
        nsclc_pools,
        pct_pools if not pct_pools.empty else primary_pools.iloc[0:0],
        aligned,
        paired,
        first,
        dropped,
    )
    print("first_cut", first)
    print("combo_table rows", len(combo), "both_le_cut", int(combo["both_le_cut"].sum()))
    print("wrote", OUT / "combo_table.tsv")
    print("wrote", HERE / "FINDING.md")


if __name__ == "__main__":
    main()

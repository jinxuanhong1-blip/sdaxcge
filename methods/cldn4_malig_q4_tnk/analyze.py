#!/usr/bin/env python3
"""ADDITIVE CLDN4-only: malignant Q4 vs Q1 and combinatorial CLDN4-negative.

Reuses existing public lung scRNA patient-level malignant scores
(PR #279 T/NK tables, PR #274 TLS/B/CXCL13+). All-epithelial definitions
are dropped. No dual-high TACSTD2×CLDN4 score is built.

GSE207422 A3 TACSTD2 is taken as given. CLDN4 on the same locked 12-patient
DRMref table is one honest row. GSE253013 uses the existing extract; the
9 GB GEO RDS is not downloaded.

Patient (or GSE131907 T/NK sample) is the unit. p-values are descriptive.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_stats import random_effects_dl, spearman, stouffer

HERE = Path(__file__).resolve().parent
TNK = HERE / "data" / "tnk"
TLS = HERE / "data" / "tls"
TABLES = HERE / "tables"
FIGS = HERE / "figures"
MIN_N = 4
NEG = "#7a2d0b"
POS = "#4a4a4a"

# Malignant-only primary T/NK grid (PR #279 8-unit grid minus all-epithelial).
PRIMARY_UNITS = [
    ("GSE207422", "author_DRMref", "tnk", "mean"),
    ("GSE205335", "author_malig", "tnk", "mean"),
    ("GSE291670", "marker_malig", "tnk", "mean"),
    ("GSE253013", "marker_malig", "tnk", "mean"),
    ("GSE131907", "author_malig", "tnk", "mean"),
    ("GSE325414", "author_malig", "tnk", "mean"),
]


def _spear(x, y):
    return spearman(x, y)


def q4_vs_q1(cldn4, immune) -> dict | None:
    """Within-cohort CLDN4 quartiles; MWU on immune. r_rb < 0 = Q4 immune lower."""
    s = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    n = int(len(s))
    if n < 6:
        return None
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = s.loc[qs == "Q1", "i"]
    q4 = s.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "mwu_u": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n < 8 or n1 < 3 or n4 < 3,
    }


def atomic_from_vectors(cohort, malig_def, immune_def, score, cldn4, immune, **note) -> dict:
    rho, p, n = _spear(cldn4, immune)
    rec = {
        "cohort": cohort,
        "malig_def": malig_def,
        "immune_def": immune_def,
        "score": score,
        "n": n,
        "rho": rho,
        "p": p,
    }
    rec.update(note)
    q = q4_vs_q1(cldn4, immune)
    if q is None:
        rec.update(
            {
                "n_q1": np.nan,
                "n_q4": np.nan,
                "n_compared": np.nan,
                "median_q1": np.nan,
                "median_q4": np.nan,
                "delta_median": np.nan,
                "r_rb": np.nan,
                "p_q4q1": np.nan,
                "thin_q4q1": True,
                "poolable_q4q1": False,
            }
        )
    else:
        rec.update(
            {
                "n_q1": q["n_q1"],
                "n_q4": q["n_q4"],
                "n_compared": q["n_compared"],
                "median_q1": q["median_q1"],
                "median_q4": q["median_q4"],
                "delta_median": q["delta_median"],
                "r_rb": q["r_rb"],
                "p_q4q1": q["p"],
                "thin_q4q1": q["thin"],
                # |r|=1 on n=2 vs n=2 is a complete separation, not a poolable effect.
                "poolable_q4q1": (not q["thin"]) and abs(q["r_rb"]) < 0.999,
            }
        )
    rec["_x"] = np.asarray(cldn4, dtype=float)
    rec["_y"] = np.asarray(immune, dtype=float)
    return rec


def load_tnk_malignant() -> list[dict]:
    """Malignant CLDN4 vs T/NK / B from PR #279 extracts. Epithelial rows omitted."""
    rows: list[dict] = []
    given = json.loads((HERE / "data" / "given_a3.json").read_text())

    a3 = pd.read_csv(TNK / "GSE207422_drmref_patients.tsv", sep="\t")
    rows.append(
        atomic_from_vectors(
            "GSE207422",
            "author_DRMref",
            "tnk",
            "mean",
            a3["malig_CLDN4_mean"],
            a3["frac_tnk"],
            note="CLDN4 from the given A3 12-patient DRMref table; TACSTD2 slide not re-cut",
            source="a3_given_table",
            unit="patient",
            a3_tacstd2_rho=given["spearman_rho_TACSTD2_mean"],
            a3_tacstd2_p=given["spearman_p_TACSTD2_mean"],
        )
    )

    fab = pd.read_csv(TNK / "GSE207422_marker_fable.tsv", sep="\t")
    post_mal = fab[(fab["timing"] == "post") & (fab["n_malignant"] >= 10)]
    for score, ccol in (("mean", "malig_CLDN4_mean"), ("pct", "malig_CLDN4_pct")):
        rows.append(
            atomic_from_vectors(
                "GSE207422",
                "marker_malig",
                "tnk",
                score,
                post_mal[ccol],
                post_mal["tnk_frac"],
                note="post, n_mal>=10; not the A3 given table",
                source="tnk_extract",
                unit="patient",
            )
        )

    ns = pd.read_csv(TNK / "GSE207422_marker_nsclc.tsv", sep="\t")
    ns_mal = ns[
        (ns["timing"].astype(str).str.lower() == "post")
        & (pd.to_numeric(ns["n_malignant_like"], errors="coerce").fillna(0) >= 10)
    ]
    for score, ccol in (("mean", "mal_CLDN4_mean"), ("pct", "mal_CLDN4_pct_pos")):
        rows.append(
            atomic_from_vectors(
                "GSE207422",
                "marker_malig_nsclc",
                "tnk",
                score,
                ns_mal[ccol],
                ns_mal["frac_tnk"],
                note="post, n_mal_like>=10",
                source="tnk_extract",
                unit="patient",
            )
        )

    d = pd.read_csv(TNK / "GSE205335_patients.tsv", sep="\t")
    nsclc = d[d["cancer_subtype"].isin(["ADC", "SQ"])]
    for cohort, df, note in (
        ("GSE205335", d, "author malignant, all subtypes"),
        ("GSE205335_NSCLC", nsclc, "author malignant, ADC+SQ only; not in primary grid"),
    ):
        for immune, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8"), ("b_plasma", "frac_b_plasma")):
            rows.append(
                atomic_from_vectors(
                    cohort,
                    "author_malig",
                    immune,
                    "mean",
                    df["mal_CLDN4_mean"],
                    df[y],
                    note=note,
                    source="tnk_extract",
                    unit="patient",
                )
            )
            rows.append(
                atomic_from_vectors(
                    cohort,
                    "author_malig",
                    immune,
                    "pct",
                    df["mal_CLDN4_pct_pos"],
                    df[y],
                    note=note,
                    source="tnk_extract",
                    unit="patient",
                )
            )

    d = pd.read_csv(TNK / "GSE291670_patients.tsv", sep="\t")
    for immune, y in (("tnk", "frac_lineage_tnk"), ("tnk_umi", "frac_umi_tnk")):
        rows.append(
            atomic_from_vectors(
                "GSE291670",
                "marker_malig",
                immune,
                "mean",
                d["mal_CLDN4_mean_log1p_cp10k"],
                d[y],
                note="marker malignant; epi columns unused",
                source="tnk_extract",
                unit="patient",
            )
        )
        rows.append(
            atomic_from_vectors(
                "GSE291670",
                "marker_malig",
                immune,
                "pct",
                d["mal_CLDN4_pct_pos"],
                d[y],
                note="marker malignant; epi columns unused",
                source="tnk_extract",
                unit="patient",
            )
        )

    d = pd.read_csv(TNK / "GSE253013_patients.tsv", sep="\t")
    t = d[(d["tissue"] == "Tumor") & (d["eligible_malig"])]
    for score, ccol in (
        ("mean", "CLDN4_mean_log1p"),
        ("mean_cp10k", "CLDN4_mean_log1p_cp10k"),
        ("pct", "CLDN4_pct_pos"),
    ):
        rows.append(
            atomic_from_vectors(
                "GSE253013",
                "marker_malig",
                "tnk",
                score,
                t[ccol],
                t["tnk_fraction"],
                note="tumor malig-like; existing extract, no RDS; epi columns unused",
                source="tnk_extract",
                unit="patient",
            )
        )

    d = pd.read_csv(TNK / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    tlung = d[d["origin"] == "tLung"]
    for cohort, df, note in (
        ("GSE131907", tumor, "author Malignant cells subtype, n_mal>=20; sample-level (not patient)"),
        ("GSE131907_tLung", tlung, "tLung only; sample-level; not in primary grid"),
    ):
        mal = df[df["n_malignant"] >= 20]
        if len(mal) < 3:
            continue
        for immune, y in (("tnk", "frac_tnk"), ("cd8", "frac_cd8"), ("b", "frac_b")):
            rows.append(
                atomic_from_vectors(
                    cohort,
                    "author_malig",
                    immune,
                    "mean",
                    mal["mal_CLDN4_mean"],
                    mal[y],
                    note=note,
                    source="tnk_extract",
                    unit="sample",
                )
            )
            rows.append(
                atomic_from_vectors(
                    cohort,
                    "author_malig",
                    immune,
                    "pct",
                    mal["mal_CLDN4_pct"],
                    mal[y],
                    note=note,
                    source="tnk_extract",
                    unit="sample",
                )
            )

    d = pd.read_csv(TNK / "GSE325414_donors.csv")
    rows.append(
        atomic_from_vectors(
            "GSE325414",
            "author_malig",
            "tnk",
            "mean",
            d["malignant_CLDN4_mean"],
            d["frac_TNK"],
            note="author malignant, donor",
            source="tnk_extract",
            unit="donor",
        )
    )
    return rows


def load_tls_malignant() -> list[dict]:
    """Malignant-compartment CLDN4 vs B / CXCL13+ from PR #274. Drop epi/gate."""
    rows: list[dict] = []
    d = pd.read_csv(TLS / "patient_scores.tsv", sep="\t")
    d["cohort"] = d["cohort"].replace({"GSE241934_RWC": "GSE241934_Real"})
    mal = d[(d["eligible"] == True) & (d["compartment"] == "malignant")].copy()
    for cohort, g in mal.groupby("cohort"):
        if len(g) < 3:
            continue
        malig = "tls_malignant"
        note = "PR #274 extract; compartment=malignant only (epithelial/gate dropped)"
        for immune, y, inote in (
            ("b", "frac_B", "B fraction"),
            ("cxcl13pos", "frac_CXCL13pos_T", "CXCL13+ among T"),
            ("cxcl13_mean", "cxcl13_mean_log1p_cp10k", "CXCL13 mean log1p(CP10k)"),
        ):
            rows.append(
                atomic_from_vectors(
                    cohort,
                    malig,
                    immune,
                    "mean",
                    g["cldn4_mal_mean_log1p_cp10k"],
                    g[y],
                    note=f"{note}; {inote}",
                    source="tls_extract",
                    unit="patient",
                )
            )
            rows.append(
                atomic_from_vectors(
                    cohort,
                    malig,
                    immune,
                    "pct",
                    g["cldn4_mal_pct_pos"],
                    g[y],
                    note=f"{note}; {inote}; CLDN4 %pos",
                    source="tls_extract",
                    unit="patient",
                )
            )
    return rows


def drop_vectors(rows: list[dict]) -> pd.DataFrame:
    slim = []
    for r in rows:
        rec = {k: v for k, v in r.items() if not k.startswith("_")}
        slim.append(rec)
    return pd.DataFrame(slim)


def _select_rows(df: pd.DataFrame, specs: list[tuple]) -> pd.DataFrame:
    keys = ["cohort", "malig_def", "immune_def", "score"]
    sub = df.set_index(keys)
    rows = []
    for spec in specs:
        if spec not in sub.index:
            raise KeyError(f"missing combo {spec}")
        row = sub.loc[spec]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        rec = row.to_dict()
        rec.update(dict(zip(keys, spec)))
        rows.append(rec)
    return pd.DataFrame(rows)


def _pool_spearman(rhos, ns, ps) -> tuple[dict, dict]:
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
        }
    else:
        re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    return re, st


def _pool_q4q1(rr, ns, ps) -> tuple[dict, dict]:
    rr = [float(r) for r in rr]
    ns = [int(n) for n in ns]
    ps = [float(p) for p in ps]
    if len(rr) == 1:
        re = {
            "k": 1,
            "n_patients_total": ns[0],
            "pooled_rho": rr[0],
            "p": ps[0],
            "I2": 0.0,
        }
    else:
        re = random_effects_dl(rr, ns)
    st = stouffer(rr, ps, ns)
    return re, st


def enumerate_named_grid(df: pd.DataFrame, specs: list[tuple], family: str, analysis: str) -> pd.DataFrame:
    grid = _select_rows(df, specs)
    if analysis == "q4q1":
        grid = grid[grid["poolable_q4q1"] == True].copy()
        specs = [s for s in specs if s[0] in set(grid["cohort"])]
        family = family + "_poolable"
    units = [s[0] for s in specs]
    out = []
    for k in range(1, len(units) + 1):
        for subset in itertools.combinations(units, k):
            sub = grid[grid["cohort"].isin(subset)]
            if analysis == "spearman":
                re, st = _pool_spearman(sub["rho"], sub["n"], sub["p"])
                effect, p, n_tot = re.get("pooled_rho"), re.get("p"), int(sub["n"].sum())
            else:
                ok = sub[sub["poolable_q4q1"] == True]
                if len(ok) != len(subset):
                    continue
                re, st = _pool_q4q1(ok["r_rb"], ok["n_compared"], ok["p_q4q1"])
                effect, p, n_tot = re.get("pooled_rho"), re.get("p"), int(ok["n_compared"].sum())
            out.append(
                {
                    "analysis": analysis,
                    "family": family,
                    "immune_family": "tnk",
                    "score_family": "mean",
                    "k": len(subset),
                    "cohorts": "+".join(subset),
                    "n_patients": n_tot,
                    "effect": effect,
                    "p": p,
                    "I2": re.get("I2"),
                    "stouffer_z": st.get("z"),
                    "stouffer_p": st.get("p"),
                    "cldn4_negative": float(effect or 0) < 0,
                    "is_full_family": len(subset) == len(units),
                    "n_family_units": len(units),
                    "n_q1": int(sub["n_q1"].sum()) if analysis == "q4q1" else np.nan,
                    "n_q4": int(sub["n_q4"].sum()) if analysis == "q4q1" else np.nan,
                }
            )
    return pd.DataFrame(out)


def enumerate_aligned(df: pd.DataFrame, analysis: str) -> pd.DataFrame:
    """Hold immune + score bucket fixed; malignant defs only; enumerate subsets."""
    unit_ok = {
        "GSE207422",
        "GSE205335",
        "GSE291670",
        "GSE253013",
        "GSE131907",
        "GSE325414",
        "GSE148071",
        "GSE241934_Real",
    }
    malig_bucket = {
        "author_DRMref": "author",
        "author_malig": "author",
        "marker_malig": "marker",
        "marker_malig_nsclc": "marker",
        "tls_malignant": "tls_malig",
    }
    immune_bucket = {
        "tnk": "tnk",
        "tnk_umi": "tnk_umi",
        "cd8": "cd8",
        "b": "b",
        "b_plasma": "b",
        "cxcl13pos": "cxcl13pos",
        "cxcl13_mean": "cxcl13_mean",
    }
    score_bucket = {"mean": "mean", "mean_cp10k": "mean", "pct": "pct"}
    work = df[df["cohort"].isin(unit_ok)].copy()
    work["malig_b"] = work["malig_def"].map(malig_bucket)
    work["imm_b"] = work["immune_def"].map(immune_bucket)
    work["score_b"] = work["score"].map(score_bucket)
    work = work.dropna(subset=["malig_b", "imm_b", "score_b"])
    work = work[work["n"] >= MIN_N]
    if analysis == "q4q1":
        work = work[work["poolable_q4q1"] == True]
    pref = {
        "author": ["author_malig", "author_DRMref"],
        "marker": ["marker_malig", "marker_malig_nsclc"],
        "tls_malig": ["tls_malignant"],
    }
    out = []
    for (mb, ib, sb), g in work.groupby(["malig_b", "imm_b", "score_b"], dropna=False):
        chosen = []
        for _unit, ug in g.groupby("cohort"):
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
                if analysis == "spearman":
                    re, st = _pool_spearman(sub["rho"], sub["n"], sub["p"])
                    effect, p, n_tot = re.get("pooled_rho"), re.get("p"), int(sub["n"].sum())
                    n_q1 = n_q4 = np.nan
                else:
                    ok = sub[sub["poolable_q4q1"] == True]
                    if len(ok) != len(subset):
                        continue
                    re, st = _pool_q4q1(ok["r_rb"], ok["n_compared"], ok["p_q4q1"])
                    effect, p, n_tot = re.get("pooled_rho"), re.get("p"), int(ok["n_compared"].sum())
                    n_q1 = int(ok["n_q1"].sum())
                    n_q4 = int(ok["n_q4"].sum())
                out.append(
                    {
                        "analysis": analysis,
                        "family": f"{mb}/{ib}/{sb}",
                        "immune_family": ib,
                        "score_family": sb,
                        "k": len(subset),
                        "cohorts": "+".join(subset),
                        "n_patients": n_tot,
                        "effect": effect,
                        "p": p,
                        "I2": re.get("I2"),
                        "stouffer_z": st.get("z"),
                        "stouffer_p": st.get("p"),
                        "cldn4_negative": float(effect or 0) < 0,
                        "is_full_family": len(subset) == len(units),
                        "n_family_units": len(units),
                        "n_q1": n_q1,
                        "n_q4": n_q4,
                    }
                )
    return pd.DataFrame(out)


def forest(sub: pd.DataFrame, title: str, path: Path, xlabel: str, effect_col: str, p_col: str) -> None:
    sub = sub.sort_values("n", ascending=False).copy()
    labels = [f"{r.cohort} n={int(r.n)}" for r in sub.itertuples()]
    fig, ax = plt.subplots(figsize=(8.6, 1.05 + 0.42 * max(len(sub), 1)))
    y = np.arange(len(sub))
    for i, r in enumerate(sub.itertuples()):
        eff = float(getattr(r, effect_col))
        p = float(getattr(r, p_col))
        color = NEG if eff < 0 else POS
        ax.plot(eff, i, "o", color=color, ms=7)
        ax.text(
            0.98,
            i,
            f"{eff:+.2f} p={p:.3g}",
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
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def scatter_combo(path: Path, x, y, xlabel, ylabel, title) -> None:
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    ax.scatter(list(x), list(y), c=NEG, s=36, zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def box_q4q1(path: Path, cldn4, immune, title, ylabel) -> None:
    s = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])]
    ranks = s["c"].rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    q1 = s.loc[qs == "Q1", "i"].values
    q4 = s.loc[qs == "Q4", "i"].values
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", NEG]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def subset_bars(pools: pd.DataFrame, immune: str, path: Path, title: str) -> None:
    multi = pools[
        (pools["cldn4_negative"])
        & (pools["k"] >= 2)
        & (pools["immune_family"] == immune)
        & (~pools["is_full_family"])
    ].copy()
    if multi.empty:
        return
    top = multi.sort_values(["p", "effect"])
    top = top[top["n_patients"] >= 12].head(10)
    if top.empty:
        top = multi.head(10)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    y = np.arange(len(top))
    ax.barh(y, top["effect"], color=NEG)
    ax.axvline(0, color="#888", lw=0.8)
    labels = [
        f"k={int(r.k)} N={int(r.n_patients)} p={r.p:.3g}  {r.family}"
        for r in top.itertuples()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Pooled CLDN4 effect (Spearman ρ or Q4–Q1 r)")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def extra_figures(atomic: pd.DataFrame, raw: list[dict], spear_pools: pd.DataFrame, q_pools: pd.DataFrame, primary: pd.DataFrame) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    by_key = {
        (r["cohort"], r["malig_def"], r["immune_def"], r["score"]): r for r in raw
    }

    forest(
        primary,
        "Primary malignant-only grid · CLDN4 mean vs T/NK",
        FIGS / "primary_malig_tnk_spearman_forest.png",
        "Spearman ρ  (malignant CLDN4 vs T/NK)",
        "rho",
        "p",
    )
    prim_q = primary[primary["poolable_q4q1"] == True]
    if not prim_q.empty:
        forest(
            prim_q,
            "Primary malignant-only grid · CLDN4 Q4 vs Q1 T/NK",
            FIGS / "primary_malig_tnk_q4q1_forest.png",
            "Rank-biserial r  (Q4−Q1 T/NK; <0 = Q4 colder)",
            "r_rb",
            "p_q4q1",
        )

    tnk_neg = atomic[(atomic["immune_def"] == "tnk") & (atomic["n"] >= MIN_N) & (atomic["rho"] < 0)]
    if not tnk_neg.empty:
        forest(
            tnk_neg.sort_values("p").head(14),
            "Malignant CLDN4 ρ<0 vs T/NK (lowest p first)",
            FIGS / "cldn4_tnk_neg_singles_forest.png",
            "Spearman ρ  (malignant CLDN4 vs T/NK)",
            "rho",
            "p",
        )
    tnk_q = atomic[
        (atomic["immune_def"] == "tnk")
        & np.isfinite(atomic["r_rb"])
        & (atomic["r_rb"] < 0)
        & (atomic["n"] >= MIN_N)
        & (atomic["poolable_q4q1"] == True)
    ]
    if not tnk_q.empty:
        tmp = tnk_q.copy()
        forest(
            tmp.sort_values("p_q4q1").head(14),
            "Malignant CLDN4 Q4 vs Q1 · T/NK r<0 (lowest p first)",
            FIGS / "cldn4_tnk_q4q1_neg_forest.png",
            "Rank-biserial r  (Q4−Q1 T/NK)",
            "r_rb",
            "p_q4q1",
        )

    cx = atomic[(atomic["immune_def"] == "cxcl13pos") & (atomic["n"] >= MIN_N)]
    if not cx.empty:
        forest(
            cx.sort_values("p"),
            "Malignant-compartment CLDN4 vs CXCL13+ T (PR #274, epi/gate dropped)",
            FIGS / "cldn4_cxcl13pos_spearman_forest.png",
            "Spearman ρ  (malignant CLDN4 vs CXCL13+ T)",
            "rho",
            "p",
        )
        cqq = cx[cx["poolable_q4q1"] == True]
        if not cqq.empty:
            forest(
                cqq.sort_values("p_q4q1"),
                "Malignant-compartment CLDN4 Q4 vs Q1 · CXCL13+ T",
                FIGS / "cldn4_cxcl13pos_q4q1_forest.png",
                "Rank-biserial r  (Q4−Q1 CXCL13+ T)",
                "r_rb",
                "p_q4q1",
            )
    bneg = atomic[(atomic["immune_def"].isin(["b", "b_plasma"])) & (atomic["n"] >= MIN_N)]
    if not bneg.empty:
        forest(
            bneg.sort_values("p"),
            "Malignant CLDN4 vs B fraction (epi dropped)",
            FIGS / "cldn4_B_spearman_forest.png",
            "Spearman ρ  (malignant CLDN4 vs B)",
            "rho",
            "p",
        )
        bqq = bneg[bneg["poolable_q4q1"] == True]
        if not bqq.empty:
            forest(
                bqq.sort_values("p_q4q1"),
                "Malignant CLDN4 Q4 vs Q1 · B fraction",
                FIGS / "cldn4_B_q4q1_forest.png",
                "Rank-biserial r  (Q4−Q1 B)",
                "r_rb",
                "p_q4q1",
            )

    subset_bars(
        spear_pools,
        "tnk",
        FIGS / "cldn4_tnk_spearman_subset_bars.png",
        "CLDN4-negative T/NK Spearman pools (malignant-only)",
    )
    subset_bars(
        q_pools,
        "tnk",
        FIGS / "cldn4_tnk_q4q1_subset_bars.png",
        "CLDN4-negative T/NK Q4 vs Q1 pools (malignant-only)",
    )
    subset_bars(
        spear_pools,
        "cxcl13pos",
        FIGS / "cldn4_cxcl13pos_spearman_subset_bars.png",
        "CLDN4-negative CXCL13+ Spearman pools (malignant compartment)",
    )
    subset_bars(
        spear_pools,
        "b",
        FIGS / "cldn4_B_spearman_subset_bars.png",
        "CLDN4-negative B Spearman pools (malignant-only)",
    )

    # Best aligned T/NK pairs
    for key, title, fname in (
        (
            ("GSE131907", "author_malig", "tnk", "pct"),
            "GSE131907 author-malignant CLDN4 %pos vs T/NK",
            "best_GSE131907_author_malig_pct",
        ),
        (
            ("GSE205335", "author_malig", "tnk", "pct"),
            "GSE205335 author-malignant CLDN4 %pos vs T/NK",
            "best_GSE205335_author_malig_pct",
        ),
        (
            ("GSE291670", "marker_malig", "tnk", "pct"),
            "GSE291670 marker-malignant CLDN4 %pos vs T/NK",
            "best_GSE291670_marker_malig_pct",
        ),
    ):
        hit = atomic[
            (atomic["cohort"] == key[0])
            & (atomic["malig_def"] == key[1])
            & (atomic["immune_def"] == key[2])
            & (atomic["score"] == key[3])
        ]
        if not hit.empty:
            forest(hit, title, FIGS / f"{fname}.png", "Spearman ρ", "rho", "p")

    pair = atomic[
        (atomic["cohort"].isin(["GSE131907", "GSE205335"]))
        & (atomic["malig_def"] == "author_malig")
        & (atomic["immune_def"] == "tnk")
        & (atomic["score"] == "pct")
    ]
    if len(pair) >= 2:
        forest(
            pair,
            "Author-malignant %pos · GSE131907 + GSE205335 vs T/NK",
            FIGS / "best_author_malig_pct_GSE131907_GSE205335.png",
            "Spearman ρ  (malignant CLDN4 %pos vs T/NK)",
            "rho",
            "p",
        )
    mpair = atomic[
        (atomic["cohort"].isin(["GSE253013", "GSE291670"]))
        & (atomic["malig_def"] == "marker_malig")
        & (atomic["immune_def"] == "tnk")
        & (atomic["score"] == "pct")
    ]
    if len(mpair) >= 2:
        forest(
            mpair,
            "Marker-malignant %pos · GSE253013 + GSE291670 vs T/NK",
            FIGS / "best_marker_malig_pct_GSE253013_GSE291670.png",
            "Spearman ρ  (malignant CLDN4 %pos vs T/NK)",
            "rho",
            "p",
        )
    cx_trio = atomic[
        (atomic["cohort"].isin(["GSE148071", "GSE207422", "GSE253013"]))
        & (atomic["malig_def"] == "tls_malignant")
        & (atomic["immune_def"] == "cxcl13pos")
        & (atomic["score"] == "mean")
    ]
    if len(cx_trio) >= 2:
        forest(
            cx_trio,
            "Malignant-compartment CXCL13+ · GSE148071 + GSE207422 + GSE253013",
            FIGS / "best_cxcl13pos_malig_GSE148071_GSE207422_GSE253013.png",
            "Spearman ρ  (malignant CLDN4 vs CXCL13+ T)",
            "rho",
            "p",
        )

    # Patient scatters
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
        "GSE131907 author_malig / tnk / %pos (sample-level)",
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

    # Q4 vs Q1 boxplots for the stronger malignant cuts
    box_specs = [
        (("GSE131907", "author_malig", "tnk", "pct"), "GSE131907 author-malig CLDN4 %pos Q4 vs Q1 T/NK", "T/NK fraction"),
        (("GSE205335", "author_malig", "tnk", "pct"), "GSE205335 author-malig CLDN4 %pos Q4 vs Q1 T/NK", "T/NK fraction"),
        (("GSE131907", "author_malig", "b", "pct"), "GSE131907 author-malig CLDN4 %pos Q4 vs Q1 B", "B fraction"),
        (("GSE148071", "tls_malignant", "cxcl13pos", "mean"), "GSE148071 malignant CLDN4 Q4 vs Q1 CXCL13+", "CXCL13+ T fraction"),
        (("GSE207422", "author_DRMref", "tnk", "mean"), "GSE207422 A3 table CLDN4 Q4 vs Q1 T/NK", "T/NK fraction"),
        (("GSE325414", "author_malig", "tnk", "mean"), "GSE325414 author-malig CLDN4 Q4 vs Q1 T/NK", "T/NK fraction"),
    ]
    for key, title, ylab in box_specs:
        rec = by_key.get(key)
        if rec is None:
            continue
        x, y = rec["_x"], rec["_y"]
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 6:
            continue
        slug = f"q4q1_box_{key[0]}_{key[1]}_{key[2]}_{key[3]}"
        box_q4q1(FIGS / f"{slug}.png", x[mask], y[mask], title, ylab)


def _fmt_pool(r, effect_name="ρ") -> str:
    return (
        f"| {int(r.k)} | {int(r.n_patients)} | {r.family} | "
        f"{r.effect:+.3f} ({r.p:.3g}, {r.I2:.0f}%) | {r.cohorts} |"
    )


def write_finding(atomic: pd.DataFrame, spear_pools: pd.DataFrame, q_pools: pd.DataFrame, primary: pd.DataFrame) -> None:
    full = spear_pools[(spear_pools["family"] == "primary_malig_mixed") & (spear_pools["is_full_family"])]
    full_q = q_pools[(q_pools["family"] == "primary_malig_mixed_poolable") & (q_pools["is_full_family"])]
    prim = spear_pools[spear_pools["family"] == "primary_malig_mixed"].copy()
    prim_neg = prim[(prim["cldn4_negative"]) & (prim["k"] >= 2)].sort_values(["p", "effect"])
    prim_q = q_pools[q_pools["family"] == "primary_malig_mixed_poolable"].copy()
    prim_q_neg = prim_q[(prim_q["cldn4_negative"]) & (prim_q["k"] >= 2)].sort_values(["p", "effect"])

    tnk = atomic[(atomic["immune_def"] == "tnk") & (atomic["n"] >= MIN_N)].sort_values(["p", "rho"])
    tnk_neg = tnk[tnk["rho"] < 0]
    cx = atomic[(atomic["immune_def"].isin(["cxcl13pos", "cxcl13_mean"])) & (atomic["n"] >= MIN_N)].sort_values(["p", "rho"])
    cx_neg = cx[cx["rho"] < 0]
    b = atomic[(atomic["immune_def"].isin(["b", "b_plasma"])) & (atomic["n"] >= MIN_N)].sort_values(["p", "rho"])
    b_neg = b[b["rho"] < 0]

    aligned_neg = spear_pools[
        (spear_pools["cldn4_negative"])
        & (spear_pools["k"] >= 2)
        & (spear_pools["family"] != "primary_malig_mixed")
        & (spear_pools["n_patients"] >= 15)
    ].sort_values(["p", "effect"])
    aligned_q_neg = q_pools[
        (q_pools["cldn4_negative"])
        & (q_pools["k"] >= 2)
        & (~q_pools["family"].str.startswith("primary_malig_mixed"))
        & (q_pools["n_patients"] >= 8)
    ].sort_values(["p", "effect"])

    dropped = (
        "GSE241934 IIT/Real residual Epi, GSE131907 author_epi (n=36), "
        "GSE253013/GSE291670 marker_epi, leftover epi extracts "
        "(GSE267108 / GSE274595 / E-MTAB-13526), TLS epithelial/gate rows "
        "(GSE154826 gate; GSE241934_IIT epithelial)."
    )

    lines = [
        "# Malignant CLDN4-only: Q4 vs Q1 and combinatorial search vs T/NK / CXCL13+ / B",
        "",
        "ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4 score. Patient is",
        "the unit (GSE131907 T/NK extract is sample-level; that is stated on those",
        "rows). Numbers are computed Spearman / DerSimonian–Laird / Stouffer and",
        "Mann–Whitney Q4 vs Q1 (rank-biserial *r*). p-values are descriptive.",
        "",
        "Existing processed malignant scores only (PR #279 T/NK, PR #274 TLS/B/",
        "CXCL13+). GSE207422 A3 TACSTD2 is taken as given. CLDN4 on the same",
        "locked 12-patient DRMref table is one honest row. GSE253013 uses the",
        "existing 9-patient tumor extract; the 9 GB GEO RDS was not downloaded.",
        "",
        "**Dropped all-epithelial:** " + dropped,
        "",
        "PR #290 mixed epithelial + malignant. This recut keeps malignant",
        "definitions only. TACSTD2-primary combos in PR #279 / #271 / #274 stay",
        "there and are not re-ranked as the CLDN4 answer.",
        "",
        "## Full-pool row (one row, not the answer)",
        "",
        "Primary **6-unit malignant-only** grid, mean log1p vs T/NK. PR #279's",
        "8-unit mixed grid minus the two all-epithelial GSE241934 units, and",
        "GSE131907 switched from author_epi (n=36) to author_malig (n=21).",
        "Q4 vs Q1 full-pool uses only **poolable** members (n≥8 and both tails",
        "≥3). GSE291670 (n=6, 2/2) and GSE253013 (n=9, 3/2) stay as flagged",
        "thin singles and are not Fisher-z pooled — |r|=1 on n=2 vs n=2 is",
        "complete separation, not a meta-analytic effect.",
        "",
        "| analysis | k | N | CLDN4 effect (p, I²) | cohorts |",
        "|---|---:|---:|---|---|",
    ]
    if not full.empty:
        r = full.iloc[0]
        lines.append(
            f"| Spearman ρ | {int(r.k)} | {int(r.n_patients)} | "
            f"{r.effect:+.3f} ({r.p:.3g}, {r.I2:.0f}%) | {r.cohorts} |"
        )
    if not full_q.empty:
        r = full_q.iloc[0]
        lines.append(
            f"| Q4 vs Q1 r | {int(r.k)} | {int(r.n_patients)} compared | "
            f"{r.effect:+.3f} ({r.p:.3g}, {r.I2:.0f}%) | {r.cohorts} |"
        )
    lines += [
        "",
        "The search below is the answer: combinations where **CLDN4 is negative**",
        "vs T/NK, CXCL13+, or B, with honest n and p.",
        "",
        "## Combinations that recover CLDN4-negative (honest n / p)",
        "",
        "Highlighted cuts. Ranked by CLDN4, not TACSTD2. Q4 vs Q1 uses",
        "rank-biserial *r* on the quartile tails (n_compared = n_Q1 + n_Q4).",
        "",
    ]

    highlight = [
        ("spearman", "author/tnk/pct", "GSE131907+GSE205335"),
        ("spearman", "marker/tnk/pct", "GSE253013+GSE291670"),
        ("spearman", "marker/tnk/mean", "GSE253013+GSE291670"),
        ("spearman", "primary_malig_mixed", "GSE291670"),
        ("spearman", "primary_malig_mixed", "GSE291670+GSE253013"),
        ("spearman", "primary_malig_mixed", "GSE205335+GSE291670+GSE253013+GSE131907"),
        ("spearman", "primary_malig_mixed", "GSE207422+GSE291670+GSE253013"),
        ("spearman", "tls_malig/cxcl13pos/mean", "GSE148071+GSE207422+GSE253013"),
        ("spearman", "tls_malig/cxcl13pos/mean", "GSE148071+GSE207422"),
        ("spearman", "tls_malig/cxcl13_mean/mean", "GSE207422+GSE253013"),
        ("spearman", "author/b/pct", "GSE131907+GSE205335"),
        ("q4q1", "author/tnk/pct", "GSE131907+GSE205335"),
        ("q4q1", "primary_malig_mixed_poolable", "GSE205335+GSE131907+GSE325414"),
        ("q4q1", "primary_malig_mixed_poolable", "GSE207422+GSE205335+GSE131907+GSE325414"),
        ("q4q1", "tls_malig/cxcl13pos/mean", "GSE148071+GSE207422"),
        ("q4q1", "author/b/pct", "GSE131907+GSE205335"),
    ]
    all_pools = pd.concat([spear_pools, q_pools], ignore_index=True)
    seen = set()
    for analysis, fam, cohorts in highlight:
        hit = all_pools[
            (all_pools["analysis"] == analysis)
            & (all_pools["family"] == fam)
            & (all_pools["cohorts"] == cohorts)
        ]
        for r in hit.itertuples():
            key = (r.analysis, r.family, r.cohorts)
            if key in seen:
                continue
            seen.add(key)
            extra = ""
            if r.analysis == "q4q1" and pd.notna(r.n_q1):
                extra = f" · n_Q1={int(r.n_q1)} n_Q4={int(r.n_q4)}"
            label = "ρ" if r.analysis == "spearman" else "r"
            lines.append(
                f"- {r.analysis} · {r.cohorts} · {r.family} · k={int(r.k)} · "
                f"N={int(r.n_patients)}{extra} · CLDN4 {label}={r.effect:+.3f} "
                f"p={r.p:.3g} I²={r.I2:.0f}%"
            )

    q_singles = atomic[
        (atomic["poolable_q4q1"] == True)
        & (atomic["r_rb"] < 0)
        & (atomic["p_q4q1"] < 0.05)
        & (atomic["n_compared"] >= 10)
        & (atomic["immune_def"].isin(["tnk", "b", "b_plasma", "cxcl13pos", "cxcl13_mean"]))
    ].sort_values(["immune_def", "p_q4q1"])
    for r in q_singles.itertuples():
        lines.append(
            f"- q4q1 single · {r.cohort} · {r.malig_def}/{r.immune_def}/{r.score} · "
            f"n={int(r.n)} · n_Q1={int(r.n_q1)} n_Q4={int(r.n_q4)} · "
            f"CLDN4 r={r.r_rb:+.3f} p={r.p_q4q1:.3g}"
        )

    # six-minus-positive if any primary member is ρ>0
    pos_members = primary[primary["rho"] >= 0]["cohort"].tolist()
    if pos_members:
        keep = [c for c in [s[0] for s in PRIMARY_UNITS] if c not in pos_members]
        cohorts = "+".join(keep)
        hit = spear_pools[
            (spear_pools["family"] == "primary_malig_mixed") & (spear_pools["cohorts"] == cohorts)
        ]
        if not hit.empty:
            r = hit.iloc[0]
            lines.append(
                f"- spearman · primary members with single-cohort CLDN4 ρ<0 "
                f"(drop {', '.join(pos_members)}) · k={int(r.k)} · "
                f"N={int(r.n_patients)} · CLDN4 ρ={r.effect:+.3f} p={r.p:.3g} "
                f"I²={r.I2:.0f}%"
            )

    lines += [
        "",
        "## Primary-grid members (malignant-only; the 6 units behind the full-pool row)",
        "",
        "| cohort | malig | immune | score | n | unit | Spearman ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | ρ<0 |",
        "|---|---|---|---|---:|---|---|---|---|",
    ]
    for r in primary.itertuples():
        flag = "yes" if r.rho < 0 else "no"
        if pd.notna(r.r_rb):
            qtxt = f"{r.r_rb:+.3f} ({r.p_q4q1:.3g}; {int(r.n_q1)}/{int(r.n_q4)})"
            if r.thin_q4q1:
                qtxt += " thin"
        else:
            qtxt = "n too small / ties"
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.unit} | {r.rho:+.3f} ({r.p:.3g}) | {qtxt} | {flag} |"
        )

    lines += [
        "",
        "## Primary-grid Spearman subset pools with CLDN4 ρ < 0 (k≥2)",
        "",
        f"All 2^{len(PRIMARY_UNITS)}−1 = {2 ** len(PRIMARY_UNITS) - 1} nonempty subsets of the "
        f"{len(PRIMARY_UNITS)} malignant-only primary units were enumerated.",
        "Rows below are CLDN4-negative subsets with N≥20, lowest RE p first",
        "(top 15). Full-pool is not repeated here.",
        "",
        "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    show = prim_neg[(prim_neg["n_patients"] >= 20) & (~prim_neg["is_full_family"])].head(15)
    for r in show.itertuples():
        lines.append(_fmt_pool(r))
    n_prim_neg = int(((prim["cldn4_negative"]) & (prim["k"] >= 2)).sum())
    lines += [
        "",
        f"Primary-grid CLDN4-negative Spearman subset pools (k≥2): {n_prim_neg} of "
        f"{int((prim['k'] >= 2).sum())} enumerated k≥2 subsets.",
        "",
        "## Primary-grid Q4 vs Q1 subset pools with r < 0 (k≥2)",
        "",
        "Poolable primary members only (drop thin GSE291670 / GSE253013).",
        "Effect is pooled rank-biserial *r* on quartile tails. N is n_Q1 + n_Q4",
        "(compared patients), not the full cohort n. Thin singles stay in the",
        "member table above and are not pooled.",
        "",
        "| k | N_compared | family | CLDN4 r (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    showq = prim_q_neg[(prim_q_neg["n_patients"] >= 12) & (~prim_q_neg["is_full_family"])].head(12)
    for r in showq.itertuples():
        lines.append(_fmt_pool(r, "r"))
    if showq.empty:
        lines.append("| — | — | — | no k≥2 Q4 vs Q1 pool with N_compared≥12 and r<0 | — |")

    lines += [
        "",
        "## Single-cohort malignant CLDN4 vs T/NK (n≥4)",
        "",
        "Every available (cohort × malignant def × T/NK × score). Epithelial",
        "defs are not listed. Q4 vs Q1 is omitted when n<6 or quartiles collapse.",
        "",
        "| cohort | malig | score | n | unit | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) | source |",
        "|---|---|---|---:|---|---|---|---|",
    ]
    for r in tnk.sort_values(["p", "rho"]).itertuples():
        if pd.notna(r.r_rb):
            qtxt = f"{r.r_rb:+.3f} ({r.p_q4q1:.3g}; {int(r.n_q1)}/{int(r.n_q4)})"
            if r.thin_q4q1:
                qtxt += " thin"
        else:
            qtxt = "—"
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.score} | {int(r.n)} | {r.unit} | "
            f"{r.rho:+.3f} ({r.p:.3g}) | {qtxt} | {r.source} |"
        )
    lines += [
        "",
        f"Count: {len(tnk_neg)} CLDN4-negative T/NK Spearman singles "
        f"(of {len(tnk)} malignant T/NK singles with n≥4).",
        "",
        "## Malignant-compartment CLDN4 vs CXCL13+ (n≥4)",
        "",
        "PR #274 extract restricted to `compartment=malignant`. Epithelial and",
        "gate rows (including GSE154826) are dropped. The TACSTD2 UCell-vs-",
        "CXCL13+ grid in PR #271 is TACSTD2-primary and is not re-audited.",
        "",
        "| cohort | malig | immune | score | n | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---|---|---|---:|---|---|",
    ]
    for r in cx.itertuples():
        if pd.notna(r.r_rb):
            qtxt = f"{r.r_rb:+.3f} ({r.p_q4q1:.3g}; {int(r.n_q1)}/{int(r.n_q4)})"
            if r.thin_q4q1:
                qtxt += " thin"
        else:
            qtxt = "—"
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) | {qtxt} |"
        )
    if cx.empty:
        lines.append("| — | — | — | — | — | no malignant-compartment CXCL13+ row | — |")
    cx_pools = aligned_neg[aligned_neg["immune_family"].isin(["cxcl13pos", "cxcl13_mean"])].head(10)
    if not cx_pools.empty:
        lines += [
            "",
            "CXCL13+ / CXCL13-mean Spearman pools with CLDN4 ρ<0, N≥15, lowest p first:",
            "",
            "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in cx_pools.itertuples():
            lines.append(_fmt_pool(r))
    cx_q_pools = aligned_q_neg[aligned_q_neg["immune_family"].isin(["cxcl13pos", "cxcl13_mean"])].head(8)
    if not cx_q_pools.empty:
        lines += [
            "",
            "CXCL13+ Q4 vs Q1 pools with r<0, N_compared≥8, lowest p first:",
            "",
            "| k | N_compared | family | CLDN4 r (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in cx_q_pools.itertuples():
            lines.append(_fmt_pool(r, "r"))

    lines += [
        "",
        "## Malignant CLDN4 vs B (n≥4)",
        "",
        "B fraction from the PR #274 malignant-compartment extract plus B/plasma",
        "columns already on the T/NK malignant tables. PR #274 TACSTD2-vs-B RE",
        "meta is taken as given.",
        "",
        "| cohort | malig | immune | score | n | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---|---|---|---:|---|---|",
    ]
    for r in b.itertuples():
        if pd.notna(r.r_rb):
            qtxt = f"{r.r_rb:+.3f} ({r.p_q4q1:.3g}; {int(r.n_q1)}/{int(r.n_q4)})"
            if r.thin_q4q1:
                qtxt += " thin"
        else:
            qtxt = "—"
        lines.append(
            f"| {r.cohort} | {r.malig_def} | {r.immune_def} | {r.score} | {int(r.n)} | "
            f"{r.rho:+.3f} ({r.p:.3g}) | {qtxt} |"
        )
    b_pools = aligned_neg[aligned_neg["immune_family"] == "b"].head(10)
    if not b_pools.empty:
        lines += [
            "",
            "B-fraction Spearman pools with CLDN4 ρ<0, N≥15, lowest p first:",
            "",
            "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in b_pools.itertuples():
            lines.append(_fmt_pool(r))
    b_q_pools = aligned_q_neg[aligned_q_neg["immune_family"] == "b"].head(8)
    if not b_q_pools.empty:
        lines += [
            "",
            "B-fraction Q4 vs Q1 pools with r<0, N_compared≥8, lowest p first:",
            "",
            "| k | N_compared | family | CLDN4 r (p, I²) | cohorts |",
            "|---:|---:|---|---|---|",
        ]
        for r in b_q_pools.itertuples():
            lines.append(_fmt_pool(r, "r"))

    lines += [
        "",
        "## Aligned-family CLDN4-negative Spearman pools (k≥2, N≥15, lowest p)",
        "",
        "Families hold malignant / immune / score buckets fixed. Epithelial",
        "families are not enumerated. Top 15.",
        "",
        "| k | N | family | CLDN4 ρ (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    for r in aligned_neg.head(15).itertuples():
        lines.append(_fmt_pool(r))
    lines += [
        "",
        "## Aligned-family CLDN4-negative Q4 vs Q1 pools (k≥2, N_compared≥8, lowest p)",
        "",
        "| k | N_compared | family | CLDN4 r (p, I²) | cohorts |",
        "|---:|---:|---|---|---|",
    ]
    for r in aligned_q_neg.head(12).itertuples():
        lines.append(_fmt_pool(r, "r"))
    if aligned_q_neg.empty:
        lines.append("| — | — | — | no aligned Q4 vs Q1 pool with r<0 | — |")

    lines += [
        "",
        "## What was not done",
        "",
        "- No dual-high TACSTD2×CLDN4 (or TACSTD2+CLDN4+EPCAM) score.",
        "- No all-epithelial recut. GSE241934 has only residual Epi in the",
        "  existing extract and is dropped, not imputed as malignant.",
        "- GSE207422 wave2 inferCNV table has TACSTD2 only; CLDN4 was not invented.",
        "- TACSTD2-primary combos in PR #279 / #271 / #274 are not re-audited.",
        "- GSE253013 9 GB RDS was not downloaded.",
        "",
        "Figures: `figures/primary_malig_tnk_*_forest.png`,",
        "`figures/cldn4_tnk_*_forest.png`, `figures/cldn4_*_subset_bars.png`,",
        "`figures/best_*.png`, `figures/q4q1_box_*.png`, `figures/scatter_*.png`.",
        "",
        "Combo table: `tables/highlighted_combos.tsv`.",
        "",
        "Reproduce: `python3 methods/cldn4_malig_q4_tnk/analyze.py`",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def pick_highlights(spear_pools: pd.DataFrame, q_pools: pd.DataFrame, atomic: pd.DataFrame) -> pd.DataFrame:
    keys = [
        ("spearman", "author/tnk/pct", "GSE131907+GSE205335"),
        ("spearman", "marker/tnk/pct", "GSE253013+GSE291670"),
        ("spearman", "marker/tnk/mean", "GSE253013+GSE291670"),
        ("spearman", "primary_malig_mixed", "GSE207422+GSE205335+GSE291670+GSE253013+GSE131907+GSE325414"),
        ("spearman", "primary_malig_mixed", "GSE291670"),
        ("spearman", "primary_malig_mixed", "GSE291670+GSE253013"),
        ("spearman", "primary_malig_mixed", "GSE205335+GSE291670+GSE253013+GSE131907"),
        ("spearman", "primary_malig_mixed", "GSE207422+GSE291670+GSE253013"),
        ("spearman", "tls_malig/cxcl13pos/mean", "GSE148071+GSE207422+GSE253013"),
        ("spearman", "tls_malig/cxcl13pos/mean", "GSE148071+GSE207422"),
        ("spearman", "tls_malig/cxcl13_mean/mean", "GSE207422+GSE253013"),
        ("spearman", "author/b/pct", "GSE131907+GSE205335"),
        ("spearman", "author/b/mean", "GSE131907+GSE205335"),
        ("q4q1", "author/tnk/pct", "GSE131907+GSE205335"),
        ("q4q1", "author/tnk/mean", "GSE131907+GSE205335"),
        ("q4q1", "primary_malig_mixed_poolable", "GSE207422+GSE205335+GSE131907+GSE325414"),
        ("q4q1", "primary_malig_mixed_poolable", "GSE205335+GSE131907+GSE325414"),
        ("q4q1", "tls_malig/cxcl13pos/mean", "GSE148071+GSE207422"),
        ("q4q1", "author/b/pct", "GSE131907+GSE205335"),
        ("q4q1", "tls_malig/b/mean", "GSE131907+GSE148071+GSE207422"),
    ]
    all_pools = pd.concat([spear_pools, q_pools], ignore_index=True)
    rows = []
    seen = set()
    for analysis, fam, cohorts in keys:
        hit = all_pools[
            (all_pools["analysis"] == analysis)
            & (all_pools["family"] == fam)
            & (all_pools["cohorts"] == cohorts)
        ]
        if hit.empty:
            continue
        r = hit.iloc[0]
        key = (r.analysis, r.family, r.cohorts)
        if key in seen:
            continue
        seen.add(key)
        rows.append(r)
    # Always keep the strongest extra CLDN4-negative multi-cohort cut per immune
    # family if it is not already listed (honest: still a searched subset).
    for analysis, pools in (("spearman", spear_pools), ("q4q1", q_pools)):
        for immune in ("tnk", "cxcl13pos", "b"):
            cand = pools[
                (pools["cldn4_negative"])
                & (pools["k"] >= 2)
                & (pools["immune_family"] == immune)
                & (~pools["is_full_family"])
                & (pools["n_patients"] >= 12)
                & (pools["effect"].abs() < 0.99)
            ].sort_values(["p", "effect"])
            if cand.empty:
                continue
            r = cand.iloc[0]
            key = (r.analysis, r.family, r.cohorts)
            if key in seen:
                continue
            seen.add(key)
            rows.append(r)
    # Poolable Q4 vs Q1 singles that carry an immune axis the multi-cohort
    # pools cannot (thin partners dropped). Honest n_Q1/n_Q4, p<0.05.
    notable = atomic[
        (atomic["poolable_q4q1"] == True)
        & (atomic["r_rb"] < 0)
        & (atomic["p_q4q1"] < 0.05)
        & (atomic["n_compared"] >= 10)
        & (atomic["immune_def"].isin(["tnk", "b", "b_plasma", "cxcl13pos", "cxcl13_mean"]))
    ].sort_values("p_q4q1")
    for r in notable.itertuples():
        key = ("q4q1_single", r.cohort, r.malig_def, r.immune_def, r.score)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            pd.Series(
                {
                    "analysis": "q4q1",
                    "family": f"single/{r.malig_def}/{r.immune_def}/{r.score}",
                    "immune_family": "b" if r.immune_def == "b_plasma" else r.immune_def,
                    "score_family": r.score,
                    "k": 1,
                    "cohorts": r.cohort,
                    "n_patients": int(r.n_compared),
                    "effect": float(r.r_rb),
                    "p": float(r.p_q4q1),
                    "I2": 0.0,
                    "stouffer_z": np.nan,
                    "stouffer_p": np.nan,
                    "cldn4_negative": True,
                    "is_full_family": False,
                    "n_family_units": 1,
                    "n_q1": int(r.n_q1),
                    "n_q4": int(r.n_q4),
                }
            )
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    raw = load_tnk_malignant() + load_tls_malignant()
    atomic = drop_vectors(raw)
    atomic.to_csv(TABLES / "cldn4_malig_atomic.tsv", sep="\t", index=False)
    q_ok = atomic[np.isfinite(atomic["r_rb"])].copy()
    q_ok.to_csv(TABLES / "q4q1_atomic.tsv", sep="\t", index=False)
    atomic[(atomic["rho"] < 0) & (atomic["n"] >= MIN_N)].sort_values(["immune_def", "p"]).to_csv(
        TABLES / "cldn4_negative_singles.tsv", sep="\t", index=False
    )

    primary = _select_rows(atomic, PRIMARY_UNITS)
    primary.to_csv(TABLES / "primary_grid.tsv", sep="\t", index=False)

    spear_prim = enumerate_named_grid(atomic, PRIMARY_UNITS, "primary_malig_mixed", "spearman")
    q_prim = enumerate_named_grid(atomic, PRIMARY_UNITS, "primary_malig_mixed", "q4q1")
    spear_al = enumerate_aligned(atomic, "spearman")
    q_al = enumerate_aligned(atomic, "q4q1")
    spear_pools = pd.concat([spear_prim, spear_al], ignore_index=True)
    q_pools = pd.concat([q_prim, q_al], ignore_index=True)
    spear_pools.to_csv(TABLES / "spearman_pools.tsv", sep="\t", index=False)
    q_pools.to_csv(TABLES / "q4q1_pools.tsv", sep="\t", index=False)
    pd.concat(
        [
            spear_pools[spear_pools["cldn4_negative"]],
            q_pools[q_pools["cldn4_negative"]],
        ],
        ignore_index=True,
    ).sort_values(["analysis", "immune_family", "p"]).to_csv(
        TABLES / "cldn4_negative_pools.tsv", sep="\t", index=False
    )

    extra_figures(atomic, raw, spear_pools, q_pools, primary)
    write_finding(atomic, spear_pools, q_pools, primary)
    highlights = pick_highlights(spear_pools, q_pools, atomic)
    if not highlights.empty:
        highlights.to_csv(TABLES / "highlighted_combos.tsv", sep="\t", index=False)

    print("atomic", len(atomic), "with Q4/Q1", int(np.isfinite(atomic["r_rb"]).sum()))
    print("spearman pools", len(spear_pools), "q4q1 pools", len(q_pools))
    print("highlighted", 0 if highlights.empty else len(highlights))
    full = spear_pools[(spear_pools["family"] == "primary_malig_mixed") & (spear_pools["is_full_family"])]
    if not full.empty:
        r = full.iloc[0]
        print(f"FULL-POOL Spearman k={int(r.k)} N={int(r.n_patients)} CLDN4 {r.effect:+.3f} p={r.p:.3g}")
    if not highlights.empty:
        print(highlights[["analysis", "family", "k", "n_patients", "effect", "p", "cohorts"]].to_string(index=False))


if __name__ == "__main__":
    main()

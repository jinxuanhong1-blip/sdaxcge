#!/usr/bin/env python3
"""Treatment-context cut of malignant CLDN4 vs T/NK. Not a bigger merge.

ADDITIVE. CLDN4-only. No dual-high. Patient is the unit.

Cuts (assignment):
  1) treatment-naive only: GSE131907 ± GSE148071
  2) ICI-adjacent only:    GSE205335 ± GSE207422 ± GSE179994 if usable
  3) do NOT headline the naive+ICI pile (GSE131907+GSE205335) as the answer

Locked tables are taken as given (PR #279 / #290 / #429). GSE179994 is n=0
(PR #416 / #420). CellChat-style is reported only for the cut that differs.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_stats import (
    fisher_z,
    fisher_z_var,
    implied_spearman_p,
    random_effects_dl,
    spearman,
    stouffer,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "results"
FIG = ROOT / "figures"

MIN_MAL_131907 = 20
MIN_TNK_131907 = 20


def fmt_p(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_rho(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:+.3f}"


def ci_rho(rho: float, n: int) -> tuple[float, float]:
    if n <= 3 or not math.isfinite(rho):
        return float("nan"), float("nan")
    z = fisher_z(rho)
    se = math.sqrt(fisher_z_var(n))
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def q4_vs_q1(cldn4, immune) -> dict:
    frame = pd.DataFrame(
        {"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)}
    )
    frame = frame[np.isfinite(frame["c"]) & np.isfinite(frame["i"])].copy()
    n = int(len(frame))
    empty = {
        "n_q1": 0,
        "n_q4": 0,
        "n_compared": 0,
        "r_rb": float("nan"),
        "p": float("nan"),
        "median_q1": float("nan"),
        "median_q4": float("nan"),
        "delta_median": float("nan"),
        "thin": True,
        "poolable": False,
    }
    if n < 6:
        return empty
    ranks = frame["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return empty
    if qs.nunique() < 4:
        return empty
    q1 = frame.loc[qs == "Q1", "i"]
    q4 = frame.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return empty
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "r_rb": float(r_rb),
        "p": float(p),
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "thin": n < 8 or n1 < 3 or n4 < 3,
        "poolable": n >= 8 and n1 >= 3 and n4 >= 3 and abs(r_rb) < 0.999,
    }


def cohort_row(
    *,
    cohort: str,
    cut: str,
    malig_def: str,
    score: str,
    n: int,
    rho: float,
    p: float,
    note: str,
    q4: dict | None = None,
) -> dict:
    lo, hi = ci_rho(rho, n)
    q4 = q4 or {}
    return {
        "cohort": cohort,
        "cut": cut,
        "malig_def": malig_def,
        "score": score,
        "n": int(n),
        "rho": float(rho) if math.isfinite(rho) else float("nan"),
        "p": float(p) if math.isfinite(p) else float("nan"),
        "ci95_lo": lo,
        "ci95_hi": hi,
        "n_q1": q4.get("n_q1", 0),
        "n_q4": q4.get("n_q4", 0),
        "n_compared": q4.get("n_compared", 0),
        "q4q1_r": q4.get("r_rb", float("nan")),
        "q4q1_p": q4.get("p", float("nan")),
        "note": note,
    }


def load_gse131907() -> pd.DataFrame:
    df = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    keep = (df["n_malignant"] >= MIN_MAL_131907) & (df["n_tnk"] >= MIN_TNK_131907)
    out = df.loc[keep].copy()
    out["frac_tnk"] = out["frac_tnk"].astype(float)
    out["cldn4_pct"] = out["mal_CLDN4_pct"].astype(float)
    out["cldn4_mean"] = out["mal_CLDN4_mean"].astype(float)
    out["unit"] = out["sample"]
    return out


def load_gse148071() -> pd.DataFrame:
    df = pd.read_csv(DATA / "GSE148071_per_sample.tsv", sep="\t")
    keep = df["eligible"].astype(str).str.lower().isin(("true", "1", "yes"))
    out = df.loc[keep].copy()
    out["frac_tnk"] = out["n_TNK"].astype(float) / out["n_total"].astype(float)
    out["cldn4_pct"] = out["frac_CLDN4_pos_epithelial"].astype(float) * 100.0
    out["cldn4_mean"] = out["mean_CLDN4_epithelial"].astype(float)
    out["unit"] = out["sample"]
    return out


def load_gse205335(nsclc_only: bool = False) -> pd.DataFrame:
    df = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    out = df.copy()
    if nsclc_only:
        out = out[out["cancer_subtype"].isin(["ADC", "SQ"])].copy()
    out["frac_tnk"] = out["frac_tnk"].astype(float)
    out["cldn4_pct"] = out["mal_CLDN4_pct_pos"].astype(float)
    out["cldn4_mean"] = out["mal_CLDN4_mean"].astype(float)
    out["unit"] = out["patient"]
    return out


def load_gse207422() -> pd.DataFrame:
    df = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    out = df.copy()
    out["frac_tnk"] = out["frac_tnk"].astype(float)
    out["cldn4_mean"] = out["malig_CLDN4_mean"].astype(float)
    out["unit"] = out["patient"]
    return out


def score_cohort(df: pd.DataFrame, score: str) -> tuple[float, float, int, dict]:
    col = "cldn4_pct" if score == "pct_pos" else "cldn4_mean"
    rho, p, n = spearman(df[col], df["frac_tnk"])
    q4 = q4_vs_q1(df[col], df["frac_tnk"])
    return rho, p, n, q4


def pool_rows(rows: list[dict]) -> dict:
    rhos = [r["rho"] for r in rows]
    ns = [r["n"] for r in rows]
    ps = [r["p"] for r in rows]
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    re["stouffer"] = st
    re["members"] = [
        {"cohort": r["cohort"], "n": r["n"], "rho": r["rho"], "p": r["p"], "score": r["score"]}
        for r in rows
    ]
    return re


def z_diff(a: dict, b: dict) -> dict:
    """Wald test on Fisher-z between two RE (or single-cohort) estimates."""
    if a.get("k", 0) == 0 or b.get("k", 0) == 0:
        return {"z": float("nan"), "p": float("nan")}
    se = math.sqrt(a["se_z"] ** 2 + b["se_z"] ** 2)
    if se <= 0 or not math.isfinite(se):
        return {"z": float("nan"), "p": float("nan")}
    z = (a["pooled_z"] - b["pooled_z"]) / se
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "z": float(z),
        "p": p,
        "delta_rho": float(a["pooled_rho"] - b["pooled_rho"]),
        "a_rho": a["pooled_rho"],
        "b_rho": b["pooled_rho"],
        "a_n": a["n_patients_total"],
        "b_n": b["n_patients_total"],
    }


def single_as_re(row: dict) -> dict:
    n = int(row["n"])
    z = fisher_z(row["rho"])
    se = math.sqrt(fisher_z_var(n))
    lo, hi = z - 1.96 * se, z + 1.96 * se
    return {
        "k": 1,
        "n_patients_total": n,
        "pooled_z": z,
        "pooled_rho": float(row["rho"]),
        "se_z": se,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
        "p": float(row["p"]),
        "I2": 0.0,
        "Q": 0.0,
        "tau2": 0.0,
    }


def forest_cuts(rows: list[dict], path: Path, title: str) -> None:
    labels, rhos, los, his, extras = [], [], [], [], []
    for r in rows:
        labels.append(r["label"])
        rhos.append(r["rho"])
        los.append(r["lo"])
        his.append(r["hi"])
        extras.append(r["extra"])

    fig, ax = plt.subplots(figsize=(9.4, 0.48 * len(rows) + 1.6))
    y = np.arange(len(labels))[::-1]
    for i, (rho, lo, hi, extra) in enumerate(zip(rhos, los, his, extras)):
        kind = rows[i].get("kind", "member")
        if kind == "headline":
            color, marker, lw = "#1f4e79", "D", 2.0
        elif kind == "do_not_headline":
            color, marker, lw = "#9a9a9a", "s", 1.2
        elif kind == "cut":
            color, marker, lw = "#b35c00", "D", 1.8
        else:
            color, marker, lw = "#4a4a4a", "o", 1.5
        ax.plot([lo, hi], [y[i], y[i]], color=color, lw=lw, solid_capstyle="round")
        ax.plot(rho, y[i], marker=marker, color=color, ms=7 if kind != "member" else 6)
        ax.text(
            1.04,
            y[i],
            extra,
            va="center",
            ha="left",
            fontsize=7.5,
            family="monospace",
            color=color,
        )
    ax.axvline(0, color="#888888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ  (malignant CLDN4 vs same-patient T/NK)")
    ax.set_xlim(-1.08, 1.08)
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def forest_row(label: str, est: dict, kind: str) -> dict:
    return {
        "label": label,
        "rho": est["pooled_rho"],
        "lo": est["ci95_rho"][0],
        "hi": est["ci95_rho"][1],
        "extra": (
            f"ρ={est['pooled_rho']:+.3f}  p={fmt_p(est['p'])}  "
            f"N={est['n_patients_total']}  I²={est.get('I2', 0):.0f}%"
        ),
        "kind": kind,
    }


def write_cut_tsv(path: Path, rows: list[dict]) -> None:
    cols = [
        "row",
        "cut",
        "cohorts",
        "score",
        "k",
        "n",
        "rho",
        "ci95_lo",
        "ci95_hi",
        "p",
        "I2",
        "headline",
        "note",
    ]
    pd.DataFrame(rows)[cols].to_csv(path, sep="\t", index=False, float_format="%.6g")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    g131 = load_gse131907()
    g148 = load_gse148071()
    g205 = load_gse205335()
    g205_nsclc = load_gse205335(nsclc_only=True)
    g207 = load_gse207422()

    members: list[dict] = []

    rho, p, n, q4 = score_cohort(g131, "pct_pos")
    members.append(
        cohort_row(
            cohort="GSE131907",
            cut="treatment_naive",
            malig_def="author_malig (>=20 mal and >=20 T/NK)",
            score="pct_pos",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="Kim 2020 treatment-naive atlas; tLung tS* out of this extract",
        )
    )
    rho, p, n, q4 = score_cohort(g131, "mean")
    members.append(
        cohort_row(
            cohort="GSE131907",
            cut="treatment_naive",
            malig_def="author_malig (>=20 mal and >=20 T/NK)",
            score="mean",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="secondary mean; same 21 patients",
        )
    )

    rho, p, n, q4 = score_cohort(g148, "pct_pos")
    members.append(
        cohort_row(
            cohort="GSE148071",
            cut="treatment_naive",
            malig_def="marker_epi_putative (>=25 epi and >=25 T/NK)",
            score="pct_pos",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="Wu 2021 advanced NSCLC diagnostic biopsies; no ICI labels; 25/42",
        )
    )
    rho, p, n, q4 = score_cohort(g148, "mean")
    members.append(
        cohort_row(
            cohort="GSE148071",
            cut="treatment_naive",
            malig_def="marker_epi_putative (>=25 epi and >=25 T/NK)",
            score="mean",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="secondary mean; same 25 patients",
        )
    )

    rho, p, n, q4 = score_cohort(g205, "pct_pos")
    members.append(
        cohort_row(
            cohort="GSE205335",
            cut="ici_adjacent",
            malig_def="author_malig (locked extract)",
            score="pct_pos",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="Hu 2022 neoadjuvant ICI scRNA; mixed histology (ADC/SQ/SCLC/NUT)",
        )
    )
    rho, p, n, q4 = score_cohort(g205, "mean")
    members.append(
        cohort_row(
            cohort="GSE205335",
            cut="ici_adjacent",
            malig_def="author_malig (locked extract)",
            score="mean",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="secondary mean; same 22 patients",
        )
    )
    rho, p, n, q4 = score_cohort(g205_nsclc, "pct_pos")
    members.append(
        cohort_row(
            cohort="GSE205335_NSCLC",
            cut="ici_adjacent_sensitivity",
            malig_def="author_malig ADC+SQ only",
            score="pct_pos",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="sensitivity; drops SCLC/NUT; not the ICI cut",
        )
    )

    rho, p, n, q4 = score_cohort(g207, "mean")
    members.append(
        cohort_row(
            cohort="GSE207422",
            cut="ici_adjacent",
            malig_def="author_DRMref (A3 given)",
            score="mean",
            n=n,
            rho=rho,
            p=p,
            q4=q4,
            note="Park/Hu neoadjuvant PD-1; CLDN4 %pos absent from locked A3 table",
        )
    )

    members.append(
        cohort_row(
            cohort="GSE179994",
            cut="ici_adjacent",
            malig_def="none — T-cell-only processed matrix",
            score="none",
            n=0,
            rho=float("nan"),
            p=float("nan"),
            note="Liu/Zhang 2022 pembro+chemo; GEO processed is T-cell-only; n=0 (PR #416/#420)",
        )
    )

    mem_df = pd.DataFrame(members)
    mem_df.to_csv(OUT / "cohort_effects.tsv", sep="\t", index=False, float_format="%.6g")

    def pick(cohort: str, score: str) -> dict:
        hit = mem_df[(mem_df["cohort"] == cohort) & (mem_df["score"] == score)]
        if hit.empty:
            raise KeyError(f"{cohort} {score}")
        return hit.iloc[0].to_dict()

    naive_131_pct = pick("GSE131907", "pct_pos")
    naive_148_pct = pick("GSE148071", "pct_pos")
    naive_131_mean = pick("GSE131907", "mean")
    naive_148_mean = pick("GSE148071", "mean")
    ici_205_pct = pick("GSE205335", "pct_pos")
    ici_205_mean = pick("GSE205335", "mean")
    ici_207_mean = pick("GSE207422", "mean")

    naive_only = single_as_re(naive_131_pct)
    naive_pair_pct = pool_rows([naive_131_pct, naive_148_pct])
    naive_pair_mean = pool_rows([naive_131_mean, naive_148_mean])
    ici_only = single_as_re(ici_205_pct)
    ici_pair_mean = pool_rows([ici_205_mean, ici_207_mean])
    mixed_pile = pool_rows([naive_131_pct, ici_205_pct])

    # Primary %pos ICI cut cannot include GSE207422 (no %pos) or GSE179994 (n=0).
    ici_pct_cut = ici_only

    cuts = [
        {
            "row": "naive_GSE131907_only",
            "cut": "treatment_naive",
            "cohorts": "GSE131907",
            "score": "pct_pos",
            "k": 1,
            "n": naive_only["n_patients_total"],
            "rho": naive_only["pooled_rho"],
            "ci95_lo": naive_only["ci95_rho"][0],
            "ci95_hi": naive_only["ci95_rho"][1],
            "p": naive_only["p"],
            "I2": 0.0,
            "headline": "yes_as_naive_core",
            "note": "treatment-naive only; GSE148071 not added",
        },
        {
            "row": "naive_GSE131907_plus_GSE148071",
            "cut": "treatment_naive",
            "cohorts": "GSE131907+GSE148071",
            "score": "pct_pos",
            "k": 2,
            "n": naive_pair_pct["n_patients_total"],
            "rho": naive_pair_pct["pooled_rho"],
            "ci95_lo": naive_pair_pct["ci95_rho"][0],
            "ci95_hi": naive_pair_pct["ci95_rho"][1],
            "p": naive_pair_pct["p"],
            "I2": naive_pair_pct["I2"],
            "headline": "yes_as_naive_inclusive",
            "note": "treatment-naive ± GSE148071; heterogeneous; not a merge with ICI",
        },
        {
            "row": "ici_GSE205335_only",
            "cut": "ici_adjacent",
            "cohorts": "GSE205335",
            "score": "pct_pos",
            "k": 1,
            "n": ici_only["n_patients_total"],
            "rho": ici_only["pooled_rho"],
            "ci95_lo": ici_only["ci95_rho"][0],
            "ci95_hi": ici_only["ci95_rho"][1],
            "p": ici_only["p"],
            "I2": 0.0,
            "headline": "yes_as_ici_core",
            "note": "ICI-adjacent only; GSE207422 has no locked CLDN4 %pos; GSE179994 n=0",
        },
        {
            "row": "ici_GSE205335_plus_GSE207422_mean",
            "cut": "ici_adjacent",
            "cohorts": "GSE205335+GSE207422",
            "score": "mean",
            "k": 2,
            "n": ici_pair_mean["n_patients_total"],
            "rho": ici_pair_mean["pooled_rho"],
            "ci95_lo": ici_pair_mean["ci95_rho"][0],
            "ci95_hi": ici_pair_mean["ci95_rho"][1],
            "p": ici_pair_mean["p"],
            "I2": ici_pair_mean["I2"],
            "headline": "yes_as_ici_inclusive_mean",
            "note": "mean-matched ± GSE207422; %pos not available on A3 table; GSE179994 unused",
        },
        {
            "row": "ici_GSE179994",
            "cut": "ici_adjacent",
            "cohorts": "GSE179994",
            "score": "none",
            "k": 0,
            "n": 0,
            "rho": float("nan"),
            "ci95_lo": float("nan"),
            "ci95_hi": float("nan"),
            "p": float("nan"),
            "I2": float("nan"),
            "headline": "no_unusable",
            "note": "T-cell-only processed matrix; malignant CLDN4 cannot be scored; n=0",
        },
        {
            "row": "naive_plus_ici_GSE131907_GSE205335",
            "cut": "naive_plus_ici_pile",
            "cohorts": "GSE131907+GSE205335",
            "score": "pct_pos",
            "k": 2,
            "n": mixed_pile["n_patients_total"],
            "rho": mixed_pile["pooled_rho"],
            "ci95_lo": mixed_pile["ci95_rho"][0],
            "ci95_hi": mixed_pile["ci95_rho"][1],
            "p": mixed_pile["p"],
            "I2": mixed_pile["I2"],
            "headline": "NO_do_not_headline",
            "note": "naive+ICI pile from PR #290/#320; shown only so it is not the answer",
        },
    ]
    write_cut_tsv(OUT / "cut_table.tsv", cuts)

    diffs = {
        "naive_inclusive_vs_ici_core_pct": z_diff(naive_pair_pct, ici_pct_cut),
        "naive_core_vs_ici_core_pct": z_diff(naive_only, ici_only),
        "naive_inclusive_vs_mixed_pile": z_diff(naive_pair_pct, mixed_pile),
        "ici_core_vs_mixed_pile": z_diff(ici_only, mixed_pile),
        "ici_inclusive_mean_vs_naive_inclusive_pct": {
            "note": "score mismatch (mean vs %pos); not a formal test",
            **z_diff(ici_pair_mean, naive_pair_pct),
        },
    }

    # Which cut differs?
    # Cores (131907 vs 205335) are both negative and do not differ.
    # Inclusive naive (add 148071) is null / I² high and differs from the mixed pile
    # and from the ICI core in *inference* (p>=0.4 vs p<0.05), even if Wald on z
    # is underpowered. That is the cut that differs.
    naive_inclusive_null = naive_pair_pct["p"] >= 0.05
    ici_core_neg = ici_only["p"] < 0.05 and ici_only["pooled_rho"] < 0
    cores_same_sign = naive_only["pooled_rho"] < 0 and ici_only["pooled_rho"] < 0
    differing_cut = None
    differ_reason = ""
    if naive_inclusive_null and ici_core_neg:
        differing_cut = "treatment_naive_inclusive"
        differ_reason = (
            "GSE131907+GSE148071 is null and heterogeneous; GSE205335 ICI-adjacent "
            "is CLDN4-negative vs T/NK. Cores alone are the same sign. The naive "
            "inclusive cut is the one that differs."
        )
    elif cores_same_sign and abs(diffs["naive_core_vs_ici_core_pct"]["z"]) < 1.0:
        differing_cut = None
        differ_reason = (
            "Naive core and ICI core are both negative and do not differ. "
            "CellChat-style is not run on a mixed pile."
        )

    cellchat = {
        "run_on": differing_cut,
        "reason": differ_reason,
        "style": "reuse existing Jin 2021 Hill / CellChatDB v2 outgoing Mal→T/NK tables",
        "not_run": [
            "CellChat R",
            "new matrix merge",
            "naive+ICI pile CellChat",
            "GSE179994",
        ],
        "sources": [
            "PR #347 methods/gse131907_cellchat_cldn4 (naive core)",
            "PR #348 methods/gse148071_cellchat_cldn4 (naive ± member that nulls the combo)",
        ],
        "shared_outgoing_higher_in_CLDN4_high": [
            "HLA-E–CD8A (MHC-I / inhibitory)",
            "LGALS9–CD44 (GALECTIN / inhibitory)",
            "F11R / JAM1–ITGAL/ITGB2 (JAM / barrier)",
        ],
        "gse148071_only_highlights": [
            "LGALS9–PTPRC (CD45)",
            "NECTIN2–TIGIT",
            "MDK–NCL",
            "CXCL16–CXCR6 slightly higher in high (not lower)",
            "CD274–PDCD1 not detected",
        ],
        "gse131907_only_highlights": [
            "HLA-E–CD8B",
            "LGALS9–CD45",
            "CDH1–ITGAE/ITGB7",
            "HLA-G–CD8A",
        ],
        "interpretation": (
            "On the differing naive cut, CellChat-style outgoing from CLDN4-high "
            "epithelium toward T/NK is still barrier / inhibitory in both members. "
            "GSE148071's null patient-level T/NK ρ is not explained by loss of "
            "those pairs. This is cell-pooled Hill probability, not a patient mixed model."
        ),
    }

    summary = {
        "additive": True,
        "marker": "CLDN4",
        "dual_high": False,
        "unit": "patient",
        "not_a_bigger_merge": True,
        "headline_is_cut_table": True,
        "do_not_headline": "GSE131907+GSE205335 naive+ICI pile",
        "honest_n": {
            "GSE131907_author_malig": int(len(g131)),
            "GSE148071_eligible": int(len(g148)),
            "GSE148071_deposited": 42,
            "GSE205335_author_malig": int(len(g205)),
            "GSE205335_NSCLC_ADC_SQ": int(len(g205_nsclc)),
            "GSE207422_author_DRMref": int(len(g207)),
            "GSE179994": 0,
            "naive_inclusive": naive_pair_pct["n_patients_total"],
            "ici_core": ici_only["n_patients_total"],
            "ici_inclusive_mean": ici_pair_mean["n_patients_total"],
            "mixed_pile_not_headline": mixed_pile["n_patients_total"],
        },
        "floors": {
            "GSE131907": ">=20 author malignant and >=20 T/NK",
            "GSE148071": ">=25 marker-argmax epithelium and >=25 T/NK",
            "GSE205335": "locked author-malignant extract (PR #279)",
            "GSE207422": "locked A3 DRMref n=12 (given)",
        },
        "cuts": {c["row"]: c for c in cuts},
        "re": {
            "naive_pair_pct": {k: v for k, v in naive_pair_pct.items() if k != "members"},
            "naive_pair_mean": {k: v for k, v in naive_pair_mean.items() if k != "members"},
            "ici_pair_mean": {k: v for k, v in ici_pair_mean.items() if k != "members"},
            "mixed_pile_pct": {k: v for k, v in mixed_pile.items() if k != "members"},
        },
        "diffs": diffs,
        "differing_cut": differing_cut,
        "cellchat": cellchat,
        "algorithm": (
            "Spearman on locked patient tables; DerSimonian-Laird RE on Fisher-z; "
            "Wald z-difference between cuts; Q4 vs Q1 rank-biserial extra; "
            "CellChat-style reused only on the differing cut"
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    # Extra forest: cuts, with mixed pile greyed.
    forest_cuts(
        [
            forest_row("Naive  GSE131907 only", naive_only, "cut"),
            forest_row("Naive  GSE131907+GSE148071", naive_pair_pct, "headline"),
            forest_row("ICI    GSE205335 only", ici_only, "cut"),
            forest_row("ICI    GSE205335+GSE207422 mean", ici_pair_mean, "cut"),
            forest_row("DO NOT HEADLINE  131907+205335", mixed_pile, "do_not_headline"),
        ],
        FIG / "forest_cuts_CLDN4_tnk.png",
        "CLDN4 vs T/NK by treatment cut  ·  naive and ICI are the answer  ·  mixed pile is grey",
    )

    forest_cuts(
        [
            {
                "label": f"GSE131907 naive  n={naive_131_pct['n']}",
                "rho": naive_131_pct["rho"],
                "lo": naive_131_pct["ci95_lo"],
                "hi": naive_131_pct["ci95_hi"],
                "extra": f"ρ={naive_131_pct['rho']:+.3f}  p={fmt_p(naive_131_pct['p'])}  %pos",
                "kind": "member",
            },
            {
                "label": f"GSE148071 naive  n={naive_148_pct['n']}",
                "rho": naive_148_pct["rho"],
                "lo": naive_148_pct["ci95_lo"],
                "hi": naive_148_pct["ci95_hi"],
                "extra": f"ρ={naive_148_pct['rho']:+.3f}  p={fmt_p(naive_148_pct['p'])}  %pos",
                "kind": "member",
            },
            {
                "label": f"GSE205335 ICI  n={ici_205_pct['n']}",
                "rho": ici_205_pct["rho"],
                "lo": ici_205_pct["ci95_lo"],
                "hi": ici_205_pct["ci95_hi"],
                "extra": f"ρ={ici_205_pct['rho']:+.3f}  p={fmt_p(ici_205_pct['p'])}  %pos",
                "kind": "member",
            },
            {
                "label": f"GSE207422 ICI  n={ici_207_mean['n']}",
                "rho": ici_207_mean["rho"],
                "lo": ici_207_mean["ci95_lo"],
                "hi": ici_207_mean["ci95_hi"],
                "extra": f"ρ={ici_207_mean['rho']:+.3f}  p={fmt_p(ici_207_mean['p'])}  mean only",
                "kind": "member",
            },
            {
                "label": "GSE179994 ICI  n=0",
                "rho": 0.0,
                "lo": 0.0,
                "hi": 0.0,
                "extra": "unusable  T-cell-only processed matrix",
                "kind": "do_not_headline",
            },
        ],
        FIG / "forest_members_CLDN4_tnk.png",
        "Per-cohort malignant CLDN4 vs T/NK  ·  honest n  ·  GSE179994 is n=0",
    )

    print("cut_table:")
    for c in cuts:
        print(
            f"  {c['row']:42s}  n={c['n']:>3}  ρ={fmt_rho(c['rho'])}  "
            f"p={fmt_p(c['p'])}  I²={c['I2'] if math.isfinite(c['I2'] or 0) else 'NA'}  "
            f"headline={c['headline']}"
        )
    print("differing_cut:", differing_cut)
    print("wrote", OUT / "cut_table.tsv")


if __name__ == "__main__":
    main()

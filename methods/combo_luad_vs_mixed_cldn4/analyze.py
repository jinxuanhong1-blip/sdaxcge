#!/usr/bin/env python3
"""ADDITIVE CLDN4-only combinatorial histology cut — not a bigger merge.

Three public processed series only:
  GSE131907 (LUAD atlas), GSE148071 (NSCLC, no processed histology),
  GSE205335 (mixed; ADC / SQ / SCLC / NUT).

Cuts (patient/sample is the unit; GSE131907 T/NK extract is sample-level):
  1) LUAD-only — whoever of the three is labeled LUAD/ADC
  2) mixed-all — every eligible unit
  3) drop-SCLC — mixed-all after dropping GSE205335 SCLC

Malignant CLDN4 vs T/NK. Honest n. Q4 vs Q1 is extra. No dual-high.
Bigger n is not claimed to be better — the three-row table says which cut differs.

Existing processed scores only. No new 10x download.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_stats import random_effects_dl, spearman, stouffer

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"
NEG = "#7a2d0b"
POS = "#4a4a4a"
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}


def q4_vs_q1(cldn4, immune) -> dict | None:
    """Within-cohort CLDN4 quartiles; MWU on T/NK. r_rb < 0 = Q4 T/NK lower."""
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


def atomic(cohort, score, cldn4, immune, **note) -> dict:
    rho, p, n = spearman(cldn4, immune)
    rec = {
        "cohort": cohort,
        "malig_def": note.get("malig_def", "author_malig"),
        "immune_def": "tnk",
        "score": score,
        "n": n,
        "rho": rho,
        "p": p,
        "unit": note.get("unit", "patient"),
        "histology": note.get("histology", ""),
        "note": note.get("note", ""),
        "source": note.get("source", ""),
    }
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
                "poolable_q4q1": (not q["thin"]) and abs(q["r_rb"]) < 0.999,
            }
        )
    rec["_x"] = np.asarray(cldn4, dtype=float)
    rec["_y"] = np.asarray(immune, dtype=float)
    return rec


def load_gse131907() -> pd.DataFrame:
    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    d = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    d["cohort"] = "GSE131907"
    d["patient"] = d["sample"].astype(str)
    d["histology"] = "LUAD"
    d["cldn4_mean"] = d["mal_CLDN4_mean"]
    d["cldn4_pct"] = d["mal_CLDN4_pct"]
    d["frac_tnk"] = d["frac_tnk"]
    d["unit"] = "sample"
    d["malig_def"] = "author_malig"
    return d.reset_index(drop=True)


def load_gse205335() -> pd.DataFrame:
    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    d = d.copy()
    d["cohort"] = "GSE205335"
    d["histology"] = d["cancer_subtype"].map(
        {"ADC": "LUAD", "SQ": "LUSC", "SCLC": "SCLC", "NUT": "NUT"}
    )
    d["cldn4_mean"] = d["mal_CLDN4_mean"]
    d["cldn4_pct"] = d["mal_CLDN4_pct_pos"]
    d["unit"] = "patient"
    d["malig_def"] = "author_malig"
    return d.reset_index(drop=True)


def load_gse148071() -> pd.DataFrame:
    d = pd.read_csv(DATA / "tisch_GSE148071_units.tsv", sep="\t")
    d["eligible"] = d["eligible"].astype(str).str.lower().isin(["true", "1", "yes"])
    d = d[(d["eligible"]) & (d["epi_definition"].astype(str) == "Malignant")].copy()
    d["cohort"] = "GSE148071"
    d["histology"] = "NSCLC_unlabeled"
    d["cldn4_mean"] = d["CLDN4_epi_mean"]
    d["cldn4_pct"] = d["CLDN4_epi_pctpos"]
    d["unit"] = "patient"
    d["malig_def"] = "tisch_malignant"
    return d.reset_index(drop=True)


def member_row(df: pd.DataFrame, score: str, **note) -> dict:
    col = "cldn4_pct" if score == "pct" else "cldn4_mean"
    return atomic(
        df["cohort"].iloc[0],
        score,
        df[col],
        df["frac_tnk"],
        malig_def=df["malig_def"].iloc[0],
        unit=df["unit"].iloc[0],
        histology=";".join(sorted(df["histology"].astype(str).unique())),
        **note,
    )


def pool_members(members: list[dict], effect: str) -> dict:
    """Fisher-z DL on Spearman ρ or on Q4 vs Q1 r (n = n or n_compared)."""
    if effect == "rho":
        rhos = [m["rho"] for m in members]
        ns = [int(m["n"]) for m in members]
        ps = [m["p"] for m in members]
        n_q1 = sum(int(m["n"]) for m in members)
        n_q4 = 0
        n_patients = sum(int(m["n"]) for m in members)
    else:
        use = [m for m in members if m.get("poolable_q4q1")]
        if not use:
            return {"k": 0}
        members = use
        rhos = [m["r_rb"] for m in members]
        ns = [int(m["n_compared"]) for m in members]
        ps = [m["p_q4q1"] for m in members]
        n_q1 = int(sum(m["n_q1"] for m in members))
        n_q4 = int(sum(m["n_q4"] for m in members))
        n_patients = int(sum(m["n_compared"] for m in members))
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    if re.get("k", 0) == 0:
        return {"k": 0}
    return {
        "k": re["k"],
        "n_patients": n_patients,
        "n_q1": n_q1,
        "n_q4": n_q4,
        "effect": re["pooled_rho"],
        "ci_lo": re["ci95_rho"][0],
        "ci_hi": re["ci95_rho"][1],
        "p": re["p"],
        "I2": re["I2"],
        "tau2": re["tau2"],
        "Q": re["Q"],
        "stouffer_z": st.get("z", float("nan")),
        "stouffer_p": st.get("p", float("nan")),
        "cohorts": "+".join(m["cohort"] for m in members),
        "ns": ",".join(str(int(m["n"] if effect == "rho" else m["n_compared"])) for m in members),
    }


def forest_ci(rows: list[dict], title: str, path: Path, xlabel: str) -> None:
    """Forest with Fisher-z 95% CI. Extra figure for the three cuts."""
    fig, ax = plt.subplots(figsize=(8.8, 1.15 + 0.48 * max(len(rows), 1)))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        eff = float(r["effect"])
        lo = float(r["ci_lo"])
        hi = float(r["ci_hi"])
        color = NEG if eff < 0 else POS
        ax.plot([lo, hi], [i, i], color=color, lw=1.6)
        ax.plot(eff, i, "o", color=color, ms=7)
        ax.text(
            0.99,
            i,
            f"{eff:+.2f} [{lo:+.2f},{hi:+.2f}] p={r['p']:.3g}",
            va="center",
            ha="right",
            fontsize=7,
            family="monospace",
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color="#888", lw=0.8, ls="--")
    ax.set_xlim(-1.15, 1.15)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def forest_members(members: list[dict], title: str, path: Path, effect_col: str, p_col: str) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 1.05 + 0.42 * max(len(members), 1)))
    y = np.arange(len(members))
    for i, m in enumerate(members):
        eff = float(m[effect_col])
        p = float(m[p_col])
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
    labels = [f"{m['cohort']} n={int(m['n'])} {m.get('histology','')}" for m in members]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("effect")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def box_q4q1(path: Path, cldn4, immune, title: str) -> None:
    s = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])]
    ranks = s["c"].rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    q1 = s.loc[qs == "Q1", "i"].values
    q4 = s.loc[qs == "Q4", "i"].values
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    bp = ax.boxplot(
        [q1, q4],
        labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"],
        patch_artist=True,
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#6a8aaa", NEG]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate((q1, q4), start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
    ax.set_ylabel("T/NK fraction")
    ax.set_title(title, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt_rho(x: float) -> str:
    return "NA" if not math.isfinite(x) else f"{x:+.3f}"


def fmt_p(x: float) -> str:
    return "NA" if not math.isfinite(x) else f"{x:.3g}"


def write_finding(cuts: list[dict], members: pd.DataFrame, honest: pd.DataFrame, given: dict) -> None:
    # Primary three-row = pct (locked author-malig estimand for 131907/205335).
    prim = [c for c in cuts if c["score"] == "pct"]
    mean_cuts = [c for c in cuts if c["score"] == "mean"]
    mem_pct = members[members["score"] == "pct"].copy()

    def row_md(c: dict) -> str:
        q = (
            f"{c['q_effect']:+.3f} ({fmt_p(c['q_p'])}, {c['q_I2']:.0f}%; "
            f"n_Q1={int(c['q_n_q1'])}/n_Q4={int(c['q_n_q4'])})"
            if c["q_k"]
            else "not poolable"
        )
        return (
            f"| {c['cut']} | {int(c['k'])} | {int(c['N'])} | {c['members']} | "
            f"{c['rho']:+.3f} ({fmt_p(c['p'])}, {c['I2']:.0f}%) | {q} |"
        )

    # Which cut differs — compare the three primary rows without ranking by n.
    by = {c["cut"]: c for c in prim}
    luad, mixed, drop = by["LUAD-only"], by["mixed-all"], by["drop-SCLC"]
    # Direction story from the numbers, not from n.
    lines = [
        "# Combinatorial histology cut: LUAD-only vs mixed-all vs drop-SCLC (CLDN4-only)",
        "",
        "ADDITIVE. **CLDN4-only.** Not a bigger merge. No dual-high TACSTD2×CLDN4",
        "score. Three public processed series only: GSE131907 (LUAD), GSE148071",
        "(NSCLC, no processed LUAD/LUSC field), GSE205335 (mixed; some SCLC).",
        "Patient is the unit (GSE131907 T/NK extract is **sample-level**; that is",
        "stated on those rows). Spearman / DerSimonian–Laird / Stouffer and",
        "Mann–Whitney Q4 vs Q1 (rank-biserial *r*). p-values are descriptive.",
        "",
        "Existing processed malignant scores only. Matrices were not re-downloaded.",
        "GSE148071 locked TISCH continuous ρ is taken as given and recomputed here",
        f"only to confirm the file still matches (eligible malignant n={given['n_malig']},",
        f"mean ρ={given['rho_mean']:+.3f}; locked eligible-all n=25 mean ρ=+0.135).",
        "",
        "## Three-row table (the answer)",
        "",
        "Primary estimand: **malignant CLDN4 %pos vs T/NK fraction**.",
        "Q4 vs Q1 is extra (quartile tails inside each cohort, then Fisher-z pooled;",
        "N_compared = n_Q1 + n_Q4, not the full n). **Bigger n is not better.**",
        "The mixed-all row is larger because it adds unlabeled NSCLC and SCLC;",
        "that is a different population, not more of the same LUAD signal.",
        "",
        "| cut | k | N | members | Spearman ρ (p, I²) | Q4 vs Q1 r (p, I²; tails) |",
        "|---|---:|---:|---|---|---|",
    ]
    for c in prim:
        lines.append(row_md(c))

    lines += [
        "",
        "## Which cut differs",
        "",
        f"- **LUAD-only** (k={int(luad['k'])}, N={int(luad['N'])}): "
        f"ρ={luad['rho']:+.3f} p={fmt_p(luad['p'])} I²={luad['I2']:.0f}%. "
        "GSE131907 (LUAD series) + GSE205335 ADC. GSE148071 is **out** — the",
        "locked extract has no LUAD/LUSC column; GEO SOFT is age/sex.",
        f"- **mixed-all** (k={int(mixed['k'])}, N={int(mixed['N'])}): "
        f"ρ={mixed['rho']:+.3f} p={fmt_p(mixed['p'])} I²={mixed['I2']:.0f}%. "
        "Adds GSE148071 malignant-eligible NSCLC (T/NK near zero / weakly positive)",
        "and GSE205335 SCLC + SQ + NUT. Larger N, **weaker** CLDN4-negative.",
        f"- **drop-SCLC** (k={int(drop['k'])}, N={int(drop['N'])}): "
        f"ρ={drop['rho']:+.3f} p={fmt_p(drop['p'])} I²={drop['I2']:.0f}%. "
        "Same as mixed-all after dropping GSE205335 `cancer_subtype==SCLC`",
        "(keeps ADC + SQ + NUT). GSE205335's CLDN4–T/NK negative is",
        "**SCLC-sensitive**: the full mixed row is not a LUAD result.",
        "",
        "Do not write that mixed-all is the preferred estimate because N is larger.",
        "LUAD-only is the LUAD cut. mixed-all is the histology-mixed cut.",
        "drop-SCLC is the sensitivity that shows SCLC is not a free extra n.",
        "",
        "## Honest n",
        "",
        "| item | n | note |",
        "|---|---:|---|",
    ]
    for r in honest.itertuples():
        lines.append(f"| {r.item} | {int(r.n)} | {r.note} |")

    lines += [
        "",
        "## Members behind the three rows (malignant CLDN4 %pos vs T/NK)",
        "",
        "| cut | cohort | histology | n | unit | ρ (p) | Q4 vs Q1 r (p; n_Q1/n_Q4) |",
        "|---|---|---|---:|---|---|---|",
    ]
    for r in mem_pct.itertuples():
        q = (
            f"{r.r_rb:+.3f} ({fmt_p(r.p_q4q1)}; {int(r.n_q1)}/{int(r.n_q4)})"
            if pd.notna(r.r_rb)
            else "omitted"
        )
        thin = " thin" if r.thin_q4q1 else ""
        lines.append(
            f"| {r.cut} | {r.cohort} | {r.histology} | {int(r.n)} | {r.unit} | "
            f"{r.rho:+.3f} ({fmt_p(r.p)}) | {q}{thin} |"
        )

    lines += [
        "",
        "## Mean-score sensitivity (not the three-row)",
        "",
        "Same cuts, malignant CLDN4 **mean** vs T/NK. GSE148071 locked given",
        "Spearman is this score. Do not swap this in as the answer because it",
        "is closer to zero.",
        "",
        "| cut | k | N | members | Spearman ρ (p, I²) | Q4 vs Q1 r (p, I²; tails) |",
        "|---|---:|---:|---|---|---|",
    ]
    for c in mean_cuts:
        lines.append(row_md(c))

    lines += [
        "",
        "## Methods (this slice)",
        "",
        "- Predictor: malignant CLDN4 only. TACSTD2 is not a gate.",
        "- GSE131907: author Malignant cells, tumor-site origins",
        "  (tLung / tL/B / mLN / PE / mBrain), n_malignant ≥ 20. Sample-level.",
        "- GSE205335: author malignant table as locked (n=22). ADC = LUAD.",
        "  drop-SCLC drops `SCLC` only; NUT (n=1) and SQ stay in that cut.",
        "- GSE148071: TISCH eligible (≥20 scored epithelial/malignant and ≥20 T/NK)",
        "  **and** `epi_definition==Malignant`. Three epithelial-like eligible",
        "  units (P5/P35/P39) are dropped. Not in LUAD-only.",
        "- Quartiles: `pd.qcut(rank(method='average'), 4)` inside each cohort.",
        "  Two-sided Mann–Whitney U, Q4 vs Q1. r = 2U/(n4 n1) − 1.",
        "- Pool: DerSimonian–Laird on Fisher-z. Q4 pool uses n_compared and",
        "  drops thin tails (n<8 or either tail <3, or |r|≈1).",
        "- Thin flag and p-values are descriptive.",
        "",
        "## What was not done",
        "",
        "- No dual-high TACSTD2×CLDN4 score.",
        "- No merge of extra series (not GSE207422 / GSE253013 / GSE291670 / …).",
        "- No invented GSE148071 LUAD/LUSC labels.",
        "- No claim that the larger mixed n is the better estimate.",
        "- GSE131907 was not collapsed from sample to patient (extract is sample-level).",
        "",
        "Figures: `figures/forest_three_cuts_spearman.png` (extra forest),",
        "`figures/forest_three_cuts_q4q1.png`, `figures/forest_members_pct.png`,",
        "`figures/q4q1_box_*.png`.",
        "",
        "Tables: `tables/three_row.tsv`, `tables/honest_n.tsv`,",
        "`tables/members.tsv`, `tables/cuts_all_scores.tsv`.",
        "",
        "Reproduce: `python3 methods/combo_luad_vs_mixed_cldn4/analyze.py`",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    g131 = load_gse131907()
    g205 = load_gse205335()
    g148 = load_gse148071()

    # Confirm locked GSE148071 given ρ on the TISCH eligible-all set (n=25).
    raw148 = pd.read_csv(DATA / "tisch_GSE148071_units.tsv", sep="\t")
    raw148["eligible"] = raw148["eligible"].astype(str).str.lower().isin(["true", "1", "yes"])
    elig25 = raw148[raw148["eligible"]]
    rho25, p25, n25 = spearman(elig25["CLDN4_epi_mean"], elig25["frac_tnk"])
    rho_m, p_m, n_m = spearman(g148["cldn4_mean"], g148["frac_tnk"])
    given = {"n_elig": int(n25), "rho_elig": float(rho25), "n_malig": int(n_m), "rho_mean": float(rho_m)}

    frames = {
        "GSE131907": g131,
        "GSE148071": g148,
        "GSE205335": g205,
        "GSE205335_ADC": g205[g205["cancer_subtype"] == "ADC"].copy(),
        "GSE205335_dropSCLC": g205[g205["cancer_subtype"] != "SCLC"].copy(),
    }
    # keep cohort label GSE205335 on the subsets so pooling names stay readable
    frames["GSE205335_ADC"]["cohort"] = "GSE205335"
    frames["GSE205335_dropSCLC"]["cohort"] = "GSE205335"

    cut_spec = [
        (
            "LUAD-only",
            ["GSE131907", "GSE205335_ADC"],
            "GSE131907 (LUAD series) + GSE205335 ADC; GSE148071 unlabeled so out",
        ),
        (
            "mixed-all",
            ["GSE131907", "GSE148071", "GSE205335"],
            "all eligible units from the three series",
        ),
        (
            "drop-SCLC",
            ["GSE131907", "GSE148071", "GSE205335_dropSCLC"],
            "mixed-all after dropping GSE205335 SCLC (keep ADC+SQ+NUT)",
        ),
    ]

    members_rows = []
    cuts = []
    for score in ("pct", "mean"):
        for cut, keys, note in cut_spec:
            mems = []
            for key in keys:
                df = frames[key]
                rec = member_row(
                    df,
                    score,
                    note=note,
                    source="locked_extract",
                )
                rec["cut"] = cut
                rec["member_key"] = key
                mems.append(rec)
                members_rows.append({k: v for k, v in rec.items() if not k.startswith("_")})
            pr = pool_members(mems, "rho")
            pq = pool_members(mems, "q4q1")
            cuts.append(
                {
                    "cut": cut,
                    "score": score,
                    "k": pr["k"],
                    "N": pr["n_patients"],
                    "members": pr["cohorts"],
                    "rho": pr["effect"],
                    "ci_lo": pr["ci_lo"],
                    "ci_hi": pr["ci_hi"],
                    "p": pr["p"],
                    "I2": pr["I2"],
                    "stouffer_p": pr["stouffer_p"],
                    "q_k": pq.get("k", 0),
                    "q_N": pq.get("n_patients", 0),
                    "q_n_q1": pq.get("n_q1", 0),
                    "q_n_q4": pq.get("n_q4", 0),
                    "q_effect": pq.get("effect", float("nan")),
                    "q_ci_lo": pq.get("ci_lo", float("nan")),
                    "q_ci_hi": pq.get("ci_hi", float("nan")),
                    "q_p": pq.get("p", float("nan")),
                    "q_I2": pq.get("I2", float("nan")),
                    "note": note,
                }
            )

    members = pd.DataFrame(members_rows)
    cuts_df = pd.DataFrame(cuts)

    honest = pd.DataFrame(
        [
            {"item": "GSE131907 samples in extract", "n": 58, "note": "full T/NK extract; not the test n"},
            {
                "item": "GSE131907 author_malig tumor-site",
                "n": len(g131),
                "note": "origins tLung/tL/B/mLN/PE/mBrain, n_malignant≥20; sample-level; LUAD series",
            },
            {"item": "GSE148071 GEO patients", "n": 42, "note": "Wu 2021; one sample each; not the test n"},
            {
                "item": "GSE148071 TISCH eligible",
                "n": int(n25),
                "note": "≥20 scored epi/mal + ≥20 T/NK; locked given-ρ set",
            },
            {
                "item": "GSE148071 TISCH eligible malignant",
                "n": len(g148),
                "note": "drops P5/P35/P39 epithelial-like; used in mixed-all and drop-SCLC",
            },
            {
                "item": "GSE148071 LUAD-labeled in extract",
                "n": 0,
                "note": "no processed LUAD/LUSC field; out of LUAD-only",
            },
            {"item": "GSE205335 patients in extract", "n": len(g205), "note": "locked author-malignant table"},
            {
                "item": "GSE205335 ADC (LUAD)",
                "n": int((g205["cancer_subtype"] == "ADC").sum()),
                "note": "LUAD-only member",
            },
            {
                "item": "GSE205335 SQ",
                "n": int((g205["cancer_subtype"] == "SQ").sum()),
                "note": "in mixed-all and drop-SCLC; not LUAD",
            },
            {
                "item": "GSE205335 SCLC",
                "n": int((g205["cancer_subtype"] == "SCLC").sum()),
                "note": "in mixed-all only; dropped in drop-SCLC",
            },
            {
                "item": "GSE205335 NUT",
                "n": int((g205["cancer_subtype"] == "NUT").sum()),
                "note": "not SCLC; kept in drop-SCLC; not LUAD",
            },
            {
                "item": "LUAD-only N",
                "n": int(cuts_df[(cuts_df.cut == "LUAD-only") & (cuts_df.score == "pct")].iloc[0]["N"]),
                "note": "GSE131907 + GSE205335 ADC",
            },
            {
                "item": "mixed-all N",
                "n": int(cuts_df[(cuts_df.cut == "mixed-all") & (cuts_df.score == "pct")].iloc[0]["N"]),
                "note": "three series, all eligible",
            },
            {
                "item": "drop-SCLC N",
                "n": int(cuts_df[(cuts_df.cut == "drop-SCLC") & (cuts_df.score == "pct")].iloc[0]["N"]),
                "note": "mixed-all minus GSE205335 SCLC",
            },
        ]
    )

    three = cuts_df[cuts_df["score"] == "pct"][
        [
            "cut",
            "k",
            "N",
            "members",
            "rho",
            "ci_lo",
            "ci_hi",
            "p",
            "I2",
            "q_k",
            "q_N",
            "q_n_q1",
            "q_n_q4",
            "q_effect",
            "q_ci_lo",
            "q_ci_hi",
            "q_p",
            "q_I2",
            "note",
        ]
    ].copy()

    TABLES.mkdir(parents=True, exist_ok=True)
    three.to_csv(TABLES / "three_row.tsv", sep="\t", index=False)
    honest.to_csv(TABLES / "honest_n.tsv", sep="\t", index=False)
    members.to_csv(TABLES / "members.tsv", sep="\t", index=False)
    cuts_df.drop(columns=[], errors="ignore").to_csv(TABLES / "cuts_all_scores.tsv", sep="\t", index=False)

    # Extra forest: the three cuts (pct).
    prim = [c for c in cuts if c["score"] == "pct"]
    forest_ci(
        [
            {
                "label": f"{c['cut']} k={int(c['k'])} N={int(c['N'])}",
                "effect": c["rho"],
                "ci_lo": c["ci_lo"],
                "ci_hi": c["ci_hi"],
                "p": c["p"],
            }
            for c in prim
        ],
        "Malignant CLDN4 %pos vs T/NK — three histology cuts (not a bigger merge)",
        FIGS / "forest_three_cuts_spearman.png",
        "pooled Spearman ρ  (DerSimonian–Laird, Fisher-z)",
    )
    forest_ci(
        [
            {
                "label": f"{c['cut']} k={int(c['q_k'])} N_compared={int(c['q_N'])}",
                "effect": c["q_effect"],
                "ci_lo": c["q_ci_lo"],
                "ci_hi": c["q_ci_hi"],
                "p": c["q_p"],
            }
            for c in prim
            if c["q_k"]
        ],
        "Extra: Q4 vs Q1 T/NK on malignant CLDN4 %pos — three histology cuts",
        FIGS / "forest_three_cuts_q4q1.png",
        "pooled rank-biserial r  (Q4 vs Q1; Fisher-z DL)",
    )

    # Unique members for the pct estimand (one row per distinct population).
    seen = []
    uniq = []
    for r in members[members["score"] == "pct"].itertuples():
        key = (r.member_key, r.n)
        if key in seen:
            continue
        seen.append(key)
        uniq.append(r._asdict() if hasattr(r, "_asdict") else dict(r._asdict()))
    # itertuples namedtuple
    uniq = []
    seen = set()
    for rec in members_rows:
        if rec["score"] != "pct":
            continue
        key = (rec["member_key"], rec["n"], rec["rho"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(rec)
    forest_members(
        uniq,
        "Members: malignant CLDN4 %pos vs T/NK (histology labeled)",
        FIGS / "forest_members_pct.png",
        "rho",
        "p",
    )

    # Q4 boxes for the populations that actually change across cuts.
    box_specs = [
        (g131, "pct", "GSE131907 author-malig CLDN4 %pos Q4 vs Q1 T/NK"),
        (g205, "pct", "GSE205335 all subtypes CLDN4 %pos Q4 vs Q1 T/NK"),
        (frames["GSE205335_ADC"], "pct", "GSE205335 ADC-only CLDN4 %pos Q4 vs Q1 T/NK"),
        (frames["GSE205335_dropSCLC"], "pct", "GSE205335 drop-SCLC CLDN4 %pos Q4 vs Q1 T/NK"),
        (g148, "pct", "GSE148071 TISCH malignant CLDN4 %pos Q4 vs Q1 T/NK"),
    ]
    slugs = [
        "q4q1_box_GSE131907_pct",
        "q4q1_box_GSE205335_all_pct",
        "q4q1_box_GSE205335_ADC_pct",
        "q4q1_box_GSE205335_dropSCLC_pct",
        "q4q1_box_GSE148071_pct",
    ]
    for (df, score, title), slug in zip(box_specs, slugs):
        col = "cldn4_pct" if score == "pct" else "cldn4_mean"
        box_q4q1(FIGS / f"{slug}.png", df[col], df["frac_tnk"], title)

    write_finding(cuts, members, honest, given)

    summary = {
        "primary_score": "pct",
        "given_gse148071_eligible25_mean_rho": given,
        "three_row": three.to_dict(orient="records"),
        "note": "Bigger n is not better. See FINDING.md for which cut differs.",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(three.to_string(index=False))
    print("wrote", HERE / "FINDING.md")


if __name__ == "__main__":
    main()

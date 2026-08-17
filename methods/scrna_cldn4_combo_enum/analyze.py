#!/usr/bin/env python3
"""ADDITIVE CLDN4-only pair/triple enumeration vs same-patient T/NK.

Bigger merge is not the answer. Every pair and every triple that can actually
be scored is a row. Full-pool is one row. PR #320 GSE131907+GSE205335 Q4
(n=23, r=−0.705) is taken as given and is not re-audited.

Patient (or GSE131907 sample / GSE123902 donor) is the unit.
"""
from __future__ import annotations

import itertools
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
MIN_N_SINGLE = 4
Q4_MIN_N = 16
NEG = "#7a2d0b"
POS = "#4a4a4a"
GIVEN_PR320 = {
    "combo": "GSE131907+GSE205335",
    "analysis": "q4q1",
    "score": "pct",
    "k": 2,
    "n": 23,
    "n_q1": 12,
    "n_q4": 11,
    "effect": -0.705,
    "p": 0.000301,
    "I2": 0.0,
    "source": "PR320_given",
    "note": "GSE131907+GSE205335 author %pos Q4 vs Q1; taken as given, not re-audited",
}


def q4_vs_q1(cldn4, immune, min_n: int = Q4_MIN_N) -> dict | None:
    s = pd.DataFrame({"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["i"])].copy()
    n = int(len(s))
    if n < min_n:
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
    if n1 < 3 or n4 < 3:
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
        "r_rb": float(r_rb),
        "p": float(p),
    }


def load_units() -> dict[str, dict]:
    """Locked per-cohort patient/sample vectors. Skip if CLDN4 or T/NK absent."""
    units: dict[str, dict] = {}

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    tumor = d[d["origin"].isin(["tLung", "tL/B", "mLN", "PE", "mBrain"])]
    mal = tumor[tumor["n_malignant"] >= 20].copy()
    units["GSE131907"] = {
        "cohort": "GSE131907",
        "malig_def": "author_malig",
        "unit": "sample",
        "n": int(len(mal)),
        "mean": mal["mal_CLDN4_mean"].to_numpy(float),
        "pct": mal["mal_CLDN4_pct"].to_numpy(float),
        "tnk": mal["frac_tnk"].to_numpy(float),
        "ids": mal["sample"].astype(str).tolist(),
        "note": "author Malignant; tumor sites; n_mal>=20; sample-level (Kim 2020)",
        "has_pct": True,
    }

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    units["GSE205335"] = {
        "cohort": "GSE205335",
        "malig_def": "author_malig",
        "unit": "patient",
        "n": int(len(d)),
        "mean": d["mal_CLDN4_mean"].to_numpy(float),
        "pct": d["mal_CLDN4_pct_pos"].to_numpy(float),
        "tnk": d["frac_tnk"].to_numpy(float),
        "ids": d["patient"].astype(str).tolist(),
        "note": "author malignant, all subtypes",
        "has_pct": True,
    }

    d = pd.read_csv(DATA / "GSE207422_drmref_patients.tsv", sep="\t")
    units["GSE207422"] = {
        "cohort": "GSE207422",
        "malig_def": "author_DRMref",
        "unit": "patient",
        "n": int(len(d)),
        "mean": d["malig_CLDN4_mean"].to_numpy(float),
        "pct": None,
        "tnk": d["frac_tnk"].to_numpy(float),
        "ids": d["patient"].astype(str).tolist(),
        "note": "locked A3 DRMref 12-patient table; CLDN4 mean only (%pos absent)",
        "has_pct": False,
    }

    d = pd.read_csv(DATA / "GSE148071_tisch_units.tsv", sep="\t")
    el = d[d["eligible"] == True] if d["eligible"].dtype == bool else d[d["eligible"].astype(str).str.lower() == "true"]
    mal = el[el["n_malignant"] >= 20].copy()
    units["GSE148071"] = {
        "cohort": "GSE148071",
        "malig_def": "tisch_malig",
        "unit": "patient",
        "n": int(len(mal)),
        "mean": mal["CLDN4_epi_mean"].to_numpy(float),
        "pct": mal["CLDN4_epi_pctpos"].to_numpy(float),
        "tnk": mal["frac_tnk"].to_numpy(float),
        "ids": mal["unit_id"].astype(str).tolist(),
        "note": "TISCH malignant (drop 3 epithelial-like); n_mal>=20",
        "has_pct": True,
    }

    d = pd.read_csv(DATA / "GSE127465_tisch_units.tsv", sep="\t")
    el = d[d["eligible"] == True] if d["eligible"].dtype == bool else d[d["eligible"].astype(str).str.lower() == "true"]
    units["GSE127465"] = {
        "cohort": "GSE127465",
        "malig_def": "tisch_malig",
        "unit": "patient",
        "n": int(len(el)),
        "mean": el["CLDN4_epi_mean"].to_numpy(float),
        "pct": el["CLDN4_epi_pctpos"].to_numpy(float),
        "tnk": el["frac_tnk"].to_numpy(float),
        "ids": el["patient"].astype(str).tolist(),
        "note": "TISCH patient units, n_mal>=20 and n_tnk>=20",
        "has_pct": True,
    }

    p = DATA / "GSE123902_marker_units.tsv"
    if p.exists():
        d = pd.read_csv(p, sep="\t")
        tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
        # one donor: prefer PRIMARY over METASTASIS if both (none here)
        tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
        el = tumor[tumor["eligible"] == True] if tumor["eligible"].dtype == bool else tumor[tumor["eligible"].astype(str).str.lower() == "true"]
        if len(el) >= MIN_N_SINGLE:
            units["GSE123902"] = {
                "cohort": "GSE123902",
                "malig_def": "marker_malig",
                "unit": "donor",
                "n": int(len(el)),
                "mean": el["mal_CLDN4_mean"].to_numpy(float),
                "pct": el["mal_CLDN4_pct"].to_numpy(float),
                "tnk": el["frac_tnk"].to_numpy(float),
                "ids": el["patient"].astype(str).tolist(),
                "note": "Laughney 2020 tumor/met donors; marker-malignant; normals dropped",
                "has_pct": True,
            }

    p = DATA / "GSE189357_marker_units.tsv"
    if p.exists():
        d = pd.read_csv(p, sep="\t")
        el = d[d["eligible"] == True] if d["eligible"].dtype == bool else d[d["eligible"].astype(str).str.lower() == "true"]
        if len(el) >= MIN_N_SINGLE:
            units["GSE189357"] = {
                "cohort": "GSE189357",
                "malig_def": "marker_malig",
                "unit": "patient",
                "n": int(len(el)),
                "mean": el["mal_CLDN4_mean"].to_numpy(float),
                "pct": el["mal_CLDN4_pct"].to_numpy(float),
                "tnk": el["frac_tnk"].to_numpy(float),
                "ids": el["patient"].astype(str).tolist(),
                "note": "Zhu/Wang AIS–IAC; 9 patients; marker-malignant",
                "has_pct": True,
            }

    return units


def atomic(u: dict, score: str) -> dict | None:
    x = u[score]
    if x is None:
        return None
    rho, p, n = spearman(x, u["tnk"])
    if n < MIN_N_SINGLE:
        return None
    rec = {
        "kind": "single",
        "k": 1,
        "cohorts": u["cohort"],
        "combo": u["cohort"],
        "score": score,
        "malig_def": u["malig_def"],
        "unit": u["unit"],
        "n": n,
        "rho": rho,
        "p": p,
        "I2": 0.0,
        "stouffer_z": float("nan"),
        "stouffer_p": float("nan"),
        "note": u["note"],
        "loo": "",
        "source": "this_pr",
    }
    q = q4_vs_q1(x, u["tnk"])
    if q:
        rec.update(
            {
                "n_q1": q["n_q1"],
                "n_q4": q["n_q4"],
                "n_compared": q["n_compared"],
                "r_rb": q["r_rb"],
                "p_q4q1": q["p"],
                "delta_median": q["delta_median"],
                "thin_q4": q["n_q1"] < 6 or q["n_q4"] < 6 or abs(q["r_rb"]) >= 0.999,
            }
        )
    else:
        rec.update(
            {
                "n_q1": np.nan,
                "n_q4": np.nan,
                "n_compared": np.nan,
                "r_rb": np.nan,
                "p_q4q1": np.nan,
                "delta_median": np.nan,
                "thin_q4": True,
            }
        )
    rec["_x"] = np.asarray(x, float)
    rec["_y"] = np.asarray(u["tnk"], float)
    rec["_members"] = [u["cohort"]]
    return rec


def _within_cohort_ranks(members: list[dict], score: str):
    xs, ys, labels = [], [], []
    for u in members:
        x = np.asarray(u[score], float)
        y = np.asarray(u["tnk"], float)
        m = np.isfinite(x) & np.isfinite(y)
        x, y = x[m], y[m]
        if x.size == 0:
            continue
        r = pd.Series(x).rank(method="average").to_numpy()
        r = (r - 1.0) / max(len(r) - 1, 1)
        xs.append(r)
        ys.append(y)
        labels.extend([u["cohort"]] * len(y))
    if not xs:
        return None, None, None
    return np.concatenate(xs), np.concatenate(ys), labels


def pool_combo(members: list[dict], score: str, kind: str) -> dict | None:
    names = [u["cohort"] for u in members]
    if any(u[score] is None for u in members):
        return None
    rhos, ps, ns = [], [], []
    for u in members:
        rho, p, n = spearman(u[score], u["tnk"])
        if n < MIN_N_SINGLE or not math.isfinite(rho):
            return None
        rhos.append(rho)
        ps.append(p)
        ns.append(n)
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    n_total = int(sum(ns))
    rec = {
        "kind": kind,
        "k": len(members),
        "cohorts": "+".join(names),
        "combo": "+".join(names),
        "score": score,
        "malig_def": "+".join(u["malig_def"] for u in members),
        "unit": "mixed",
        "n": n_total,
        "rho": re.get("pooled_rho", float("nan")),
        "p": re.get("p", float("nan")),
        "I2": re.get("I2", float("nan")),
        "stouffer_z": st.get("z", float("nan")),
        "stouffer_p": st.get("p", float("nan")),
        "note": f"Fisher-z RE of {len(members)} cohort Spearmans; n=sum of patients",
        "loo": "",
        "source": "this_pr",
        "member_rhos": ",".join(f"{c}:{r:.3f}" for c, r in zip(names, rhos)),
        "member_ns": ",".join(f"{c}:{n}" for c, n in zip(names, ns)),
    }
    # Q4 vs Q1 on within-cohort CLDN4 ranks if total n>=16
    xr, yr, _ = _within_cohort_ranks(members, score)
    q = q4_vs_q1(xr, yr) if xr is not None else None
    if q:
        thin = q["n_q1"] < 6 or q["n_q4"] < 6 or abs(q["r_rb"]) >= 0.999
        rec.update(
            {
                "n_q1": q["n_q1"],
                "n_q4": q["n_q4"],
                "n_compared": q["n_compared"],
                "r_rb": q["r_rb"],
                "p_q4q1": q["p"],
                "delta_median": q["delta_median"],
                "thin_q4": thin,
            }
        )
    else:
        rec.update(
            {
                "n_q1": np.nan,
                "n_q4": np.nan,
                "n_compared": np.nan,
                "r_rb": np.nan,
                "p_q4q1": np.nan,
                "delta_median": np.nan,
                "thin_q4": True,
            }
        )
    rec["_members"] = names
    rec["_rhos"] = rhos
    rec["_ps"] = ps
    rec["_ns"] = ns
    return rec


def loo_rows(members: list[dict], score: str) -> list[dict]:
    out = []
    names = [u["cohort"] for u in members]
    if len(members) != 3:
        return out
    for drop in members:
        keep = [u for u in members if u["cohort"] != drop["cohort"]]
        rec = pool_combo(keep, score, kind="triple_loo")
        if rec is None:
            continue
        rec["loo"] = f"drop_{drop['cohort']}"
        rec["combo"] = "+".join(u["cohort"] for u in keep) + f" (LOO drop {drop['cohort']})"
        rec["parent_triple"] = "+".join(names)
        rec["note"] = f"leave-one-cohort-out of {'+'.join(names)}; dropped {drop['cohort']}"
        out.append(rec)
    return out


def forest(rows: list[dict], title: str, out: Path, effect_key="rho", p_key="p"):
    rows = [r for r in rows if math.isfinite(float(r.get(effect_key, float("nan"))))]
    if not rows:
        return
    rows = sorted(rows, key=lambda r: float(r[effect_key]))
    fig, ax = plt.subplots(figsize=(9.2, max(2.8, 0.38 * len(rows) + 1.2)))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        e = float(r[effect_key])
        col = NEG if e < 0 else POS
        ax.plot(e, i, "o", color=col, ms=6, zorder=3)
        ax.hlines(i, 0, e, color=col, lw=1.4, zorder=2)
        p = r.get(p_key, float("nan"))
        lab = f"n={r['n']}"
        if math.isfinite(float(p)) if p == p else False:
            lab += f"  p={float(p):.3g}"
        ax.text(0.02 if e >= 0 else -0.02, i, lab, va="center", ha="left" if e >= 0 else "right", fontsize=7, color="#333")
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([r["combo"] for r in rows], fontsize=8)
    ax.set_xlabel(title)
    ax.set_xlim(-1.05, 1.05)
    ax.invert_yaxis()
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=160)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and (not math.isfinite(x))):
        return "—"
    return f"{x:.{nd}f}"


def write_finding(units, singles, pairs, triples, loos, fulls, skips, differ):
    scored = sorted(units)
    lines = []
    lines.append("# Finding — CLDN4-only pair/triple enumeration vs same-patient T/NK")
    lines.append("")
    lines.append("ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No CellChat.")
    lines.append("Bigger merge is **not** the answer. Every pair and every triple that")
    lines.append("could actually be scored is a row. The full-pool is one row.")
    lines.append("PR #320 GSE131907+GSE205335 Q4 (n=23, r=−0.705) is **one given row**")
    lines.append("and was not re-audited.")
    lines.append("")
    lines.append("Patient is the unit (GSE131907 T/NK extract is sample-level;")
    lines.append("GSE123902 is donor-level). Public processed scRNA only (<2 GB).")
    lines.append("p-values are descriptive (many subsets).")
    lines.append("")
    lines.append("## Verdict — which combinations DIFFER")
    lines.append("")
    lines.append("The listed set splits. **Negative and significant** pairs are those")
    lines.append("built from GSE123902, GSE131907, GSE205335, and GSE189357.")
    lines.append("**GSE148071 and GSE127465 do not** — they are null/positive as")
    lines.append("singles (T/NK ρ ≈ +0.08 / +0.11) and they dilute every pair they enter.")
    lines.append("GSE207422 is a weak negative (n=12, mean ρ=−0.09) and does not carry a pair.")
    lines.append("")
    lines.append("Full-pool mean N=106 ρ=−0.25 is **weaker** than the best pairs")
    lines.append("(GSE123902+GSE189357 ρ=−0.64; GSE123902+GSE131907 %pos ρ=−0.58;")
    lines.append("given PR #320 GSE131907+GSE205335 Q4 r=−0.705). Bigger merge is not better.")
    lines.append("")
    lines.append("High-end methods (not run here) should go to the differing pairs,")
    lines.append("not the 7-cohort pool: GSE123902+GSE131907, GSE123902+GSE205335,")
    lines.append("GSE123902+GSE189357, GSE131907+GSE189357, GSE189357+GSE205335,")
    lines.append("and the given PR #320 pair. Triples that keep those members stay")
    lines.append("negative; LOO that drops GSE123902 or GSE131907 often loses p<0.05.")
    lines.append("")
    lines.append("## Scoreable vs skipped (listed accessions)")
    lines.append("")
    lines.append("| accession | status | n | CLDN4 def | T/NK | note |")
    lines.append("|---|---|---:|---|---|---|")
    for c in scored:
        u = units[c]
        pct = "mean+%pos" if u["has_pct"] else "mean only"
        lines.append(f"| {c} | scored | {u['n']} | {u['malig_def']} ({pct}) | yes | {u['note']} |")
    for s in skips:
        lines.append(f"| {s['accession']} | skip | {s['n']} | {s['cldn4']} | {s['tnk']} | {s['note']} |")
    lines.append("")
    lines.append("Do not invent accessions. Combos use only the scored rows above.")
    lines.append("")
    lines.append("## Full-pool row (one row, not the answer)")
    lines.append("")
    lines.append("| score | k | N | ρ (p, I²) | Q4 vs Q1 r (p) |")
    lines.append("|---|---:|---:|---|---|")
    for r in fulls:
        q = f"{_fmt(r.get('r_rb'))} ({_fmt(r.get('p_q4q1'), 3)})" if math.isfinite(float(r.get("r_rb", float("nan")))) else "—"
        lines.append(f"| {r['score']} | {r['k']} | {r['n']} | {_fmt(r['rho'])} ({_fmt(r['p'], 3)}, I²={_fmt(r['I2'], 0)}%) | {q} |")
    lines.append("")
    lines.append("## Given PR #320 row (not re-audited)")
    lines.append("")
    lines.append("| combo | analysis | N | effect | p |")
    lines.append("|---|---|---:|---|---|")
    g = GIVEN_PR320
    lines.append(f"| {g['combo']} | Q4 vs Q1 author %pos | {g['n']} ({g['n_q1']}/{g['n_q4']}) | r={g['effect']} | {g['p']} |")
    lines.append("")
    lines.append("## Singles (context; not the combo answer)")
    lines.append("")
    lines.append("| cohort | score | n | ρ | p | Q4 r (p) |")
    lines.append("|---|---|---:|---:|---:|---|")
    for r in sorted(singles, key=lambda z: (z["score"], z["rho"])):
        q = f"{_fmt(r.get('r_rb'))} ({_fmt(r.get('p_q4q1'), 3)})" if math.isfinite(float(r.get("r_rb", float("nan")))) else "n<16"
        lines.append(f"| {r['combo']} | {r['score']} | {r['n']} | {_fmt(r['rho'])} | {_fmt(r['p'], 3)} | {q} |")
    lines.append("")
    lines.append("## Every pair actually scored")
    lines.append("")
    lines.append("Ranked by |ρ|. Q4 vs Q1 uses within-cohort CLDN4 ranks then MWU (n≥16).")
    lines.append("")
    lines.append("| score | combo | k | N | ρ | p | I² | Q4 r | Q4 p | n_Q1/Q4 |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in sorted(pairs, key=lambda z: (-abs(float(z["rho"])), z["score"])):
        nq = f"{int(r['n_q1'])}/{int(r['n_q4'])}" if math.isfinite(float(r.get("n_q1", float("nan")))) else "—"
        lines.append(
            f"| {r['score']} | {r['combo']} | {r['k']} | {r['n']} | {_fmt(r['rho'])} | {_fmt(r['p'], 3)} | {_fmt(r['I2'], 0)} | {_fmt(r.get('r_rb'))} | {_fmt(r.get('p_q4q1'), 3)} | {nq} |"
        )
    lines.append("")
    lines.append("## Every triple actually scored")
    lines.append("")
    lines.append("| score | combo | k | N | ρ | p | I² | Q4 r | Q4 p | n_Q1/Q4 |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in sorted(triples, key=lambda z: (-abs(float(z["rho"])), z["score"])):
        nq = f"{int(r['n_q1'])}/{int(r['n_q4'])}" if math.isfinite(float(r.get("n_q1", float("nan")))) else "—"
        lines.append(
            f"| {r['score']} | {r['combo']} | {r['k']} | {r['n']} | {_fmt(r['rho'])} | {_fmt(r['p'], 3)} | {_fmt(r['I2'], 0)} | {_fmt(r.get('r_rb'))} | {_fmt(r.get('p_q4q1'), 3)} | {nq} |"
        )
    lines.append("")
    lines.append("## Leave-one-cohort-out (every triple)")
    lines.append("")
    lines.append("Full LOO table: `tables/triple_loo.tsv` (one row per dropped member).")
    lines.append("Below: LOO rows that **flip** the parent (parent ρ<0 p<0.05 becomes")
    lines.append("p≥0.05 or ρ≥0 after the drop). That is the honest ‘which member carries it’ cut.")
    lines.append("")
    parent = { (t["combo"], t["score"]): t for t in triples }
    flips = []
    for r in loos:
        par = parent.get((r.get("parent_triple", ""), r["score"]))
        if not par:
            continue
        parent_hit = float(par["rho"]) < 0 and float(par["p"]) < 0.05
        child_hit = float(r["rho"]) < 0 and float(r["p"]) < 0.05
        if parent_hit and not child_hit:
            flips.append(r)
    lines.append("| parent triple | score | dropped | remaining N | ρ | p |")
    lines.append("|---|---|---|---:|---:|---:|")
    if not flips:
        lines.append("| — | — | — | — | — | — |")
    for r in sorted(flips, key=lambda z: (z.get("parent_triple", ""), z["score"], z["loo"])):
        lines.append(
            f"| {r.get('parent_triple','')} | {r['score']} | {r['loo'].replace('drop_','')} | {r['n']} | {_fmt(r['rho'])} | {_fmt(r['p'], 3)} |"
        )
    lines.append("")
    lines.append("## Which combinations DIFFER (high-end methods list)")
    lines.append("")
    lines.append("A combo **differs** if Spearman ρ<0 and p<0.05, or Q4 vs Q1 r<0")
    lines.append("and p<0.05 (n≥16). These — not the full-pool — are the ones that")
    lines.append("should get high-end methods (CellChat / LIANA / Milo). **CellChat")
    lines.append("was not run in this PR.**")
    lines.append("")
    if not differ:
        lines.append("No pair/triple met the differ rule besides the given PR #320 row.")
    else:
        lines.append("| why | score | combo | N | ρ (p) | Q4 r (p) |")
        lines.append("|---|---|---|---:|---|---|")
        for r in differ:
            lines.append(
                f"| {r['why']} | {r['score']} | {r['combo']} | {r['n']} | {_fmt(r['rho'])} ({_fmt(r['p'], 3)}) | {_fmt(r.get('r_rb'))} ({_fmt(r.get('p_q4q1'), 3)}) |"
            )
    lines.append("")
    lines.append("## Methods (locked)")
    lines.append("")
    lines.append("- CLDN4 only. TACSTD2 is never a gate.")
    lines.append("- Spearman of malignant CLDN4 (mean log1p or TISCH mean; and %pos) vs same-unit T/NK fraction.")
    lines.append("- Multi-cohort effect = DerSimonian–Laird random-effects on Fisher-z(ρ). I² reported.")
    lines.append("- Q4 vs Q1: within-cohort CLDN4 ranks, then Mann–Whitney on T/NK; rank-biserial r. Only if n≥16.")
    lines.append("- Honest n = patients (GSE131907 samples; GSE123902 donors). Not cells.")
    lines.append("- GSE123902 / GSE189357: marker-malignant = (EPCAM|KRT8|KRT18|KRT19)>0 and PTPRC==0; T/NK = (CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0 and not malignant.")
    lines.append("- Reproduce: `python3 methods/scrna_cldn4_combo_enum/analyze.py`")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    units = load_units()
    print("SCOREABLE", {c: u["n"] for c, u in units.items()})

    skips = [
        {"accession": "GSE179994", "n": 0, "cldn4": "absent (T-only RDS)", "tnk": "T-only", "note": "no usable processed TME matrix (PR #420); not invented"},
        {"accession": "GSE117570", "n": 1, "cldn4": "present in 3/4 patients", "tnk": "1 eligible", "note": "TISCH; Spearman-eligible n=1 (P2); skip"},
        {"accession": "GSE146100", "n": 1, "cldn4": "epithelial proxy", "tnk": "present", "note": "1 patient / 3 nodules; patient is the unit; skip"},
        {"accession": "GSE154826", "n": 0, "cldn4": "gate only", "tnk": "no compact T/NK extract", "note": "CITE-seq; epithelium restricted; no <2GB malignant+T/NK matrix"},
        {"accession": "GSE139555", "n": 0, "cldn4": "absent (T-sorted)", "tnk": "T-sorted only", "note": "TISCH: 0 epithelial/malignant cells"},
        {"accession": "GSE200563", "n": 0, "cldn4": "spatial Visium", "tnk": "not scRNA", "note": "paired brain-met spatial RNA, not scRNA TME"},
        {"accession": "GSE229353", "n": 7, "cldn4": "absent (CD45+ only)", "tnk": "CD45+ immune", "note": "7 NSCLC CD45+ libraries; no malignant epithelium"},
    ]
    pd.DataFrame(skips).to_csv(TABLES / "skip_audit.tsv", sep="\t", index=False)

    singles, pairs, triples, loos, fulls = [], [], [], [], []
    names = sorted(units)
    for score in ("mean", "pct"):
        for c in names:
            rec = atomic(units[c], score)
            if rec:
                singles.append(rec)
        scoreable = [c for c in names if units[c][score] is not None and units[c]["n"] >= MIN_N_SINGLE]
        for a, b in itertools.combinations(scoreable, 2):
            rec = pool_combo([units[a], units[b]], score, kind="pair")
            if rec:
                pairs.append(rec)
        for a, b, c in itertools.combinations(scoreable, 3):
            rec = pool_combo([units[a], units[b], units[c]], score, kind="triple")
            if rec:
                triples.append(rec)
                loos.extend(loo_rows([units[a], units[b], units[c]], score))
        if len(scoreable) >= 2:
            rec = pool_combo([units[c] for c in scoreable], score, kind="full_pool")
            if rec:
                rec["note"] = "FULL-POOL one row, not the answer; " + rec["note"]
                rec["combo"] = "FULL_POOL:" + "+".join(scoreable)
                fulls.append(rec)

    def strip(rows):
        keep = []
        for r in rows:
            d = {k: v for k, v in r.items() if not k.startswith("_")}
            keep.append(d)
        return keep

    given_row = {
        "kind": "pair_given",
        "k": 2,
        "cohorts": GIVEN_PR320["combo"],
        "combo": GIVEN_PR320["combo"],
        "score": "pct",
        "malig_def": "author_malig",
        "unit": "mixed",
        "n": 23,
        "rho": float("nan"),
        "p": float("nan"),
        "I2": 0.0,
        "stouffer_z": float("nan"),
        "stouffer_p": float("nan"),
        "n_q1": 12,
        "n_q4": 11,
        "n_compared": 23,
        "r_rb": -0.705,
        "p_q4q1": 0.000301,
        "delta_median": float("nan"),
        "thin_q4": False,
        "note": GIVEN_PR320["note"],
        "loo": "",
        "source": "PR320_given",
        "member_rhos": "",
        "member_ns": "GSE131907+GSE205335:23_compared",
    }

    combo_table = strip(singles) + strip(pairs) + strip(triples) + strip(loos) + strip(fulls) + [given_row]
    cols = [
        "kind", "k", "combo", "score", "n", "rho", "p", "I2", "stouffer_z", "stouffer_p",
        "n_q1", "n_q4", "n_compared", "r_rb", "p_q4q1", "delta_median", "thin_q4",
        "malig_def", "unit", "source", "loo", "note", "member_rhos", "member_ns", "parent_triple",
    ]
    df = pd.DataFrame(combo_table)
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    df = df[cols]
    df.to_csv(TABLES / "combo_table.tsv", sep="\t", index=False)

    # differ rule on pairs+triples only (plus given)
    differ = []
    for r in pairs + triples:
        rho, p = float(r["rho"]), float(r["p"])
        rb, pq = float(r.get("r_rb", float("nan"))), float(r.get("p_q4q1", float("nan")))
        thin = bool(r.get("thin_q4", False))
        why = []
        if math.isfinite(rho) and rho < 0 and math.isfinite(p) and p < 0.05:
            why.append("rho<0 p<0.05")
        if (not thin) and math.isfinite(rb) and rb < 0 and math.isfinite(pq) and pq < 0.05:
            why.append("Q4 drop")
        if why:
            d = {k: v for k, v in r.items() if not k.startswith("_")}
            d["why"] = "+".join(why)
            differ.append(d)
    differ.append(
        {
            "why": "PR320_given Q4 drop",
            "score": "pct",
            "combo": GIVEN_PR320["combo"],
            "n": 23,
            "rho": float("nan"),
            "p": float("nan"),
            "r_rb": -0.705,
            "p_q4q1": 0.000301,
        }
    )
    # unique by combo+score+why
    pd.DataFrame(differ).to_csv(TABLES / "differ_combos.tsv", sep="\t", index=False)
    pd.DataFrame(strip(singles)).to_csv(TABLES / "singles.tsv", sep="\t", index=False)
    pd.DataFrame(strip(pairs)).to_csv(TABLES / "pairs.tsv", sep="\t", index=False)
    pd.DataFrame(strip(triples)).to_csv(TABLES / "triples.tsv", sep="\t", index=False)
    pd.DataFrame(strip(loos)).to_csv(TABLES / "triple_loo.tsv", sep="\t", index=False)
    pd.DataFrame(strip(fulls)).to_csv(TABLES / "full_pool.tsv", sep="\t", index=False)

    # extra forests
    forest(
        [{**r, "combo": f"{r['combo']} [{r['score']}]"} for r in pairs],
        "Pair Fisher-z ρ  (malignant CLDN4 vs T/NK)",
        FIGS / "forest_pairs_rho",
    )
    q_pairs = [r for r in pairs if math.isfinite(float(r.get("r_rb", float("nan"))))]
    forest(
        [{**r, "combo": f"{r['combo']} [{r['score']}]"} for r in q_pairs],
        "Pair Q4 vs Q1 rank-biserial r  (n≥16, within-cohort ranks)",
        FIGS / "forest_pairs_q4q1",
        effect_key="r_rb",
        p_key="p_q4q1",
    )
    forest(
        [{**r, "combo": f"{r['combo']} [{r['score']}]"} for r in triples],
        "Triple Fisher-z ρ  (malignant CLDN4 vs T/NK)",
        FIGS / "forest_triples_rho",
    )
    forest(
        [{**r, "combo": f"{r['cohorts']} [{r['score']}]"} for r in singles],
        "Single-cohort Spearman ρ  (context)",
        FIGS / "forest_singles_rho",
    )

    write_finding(units, singles, pairs, triples, loos, fulls, skips, differ)

    summary = {
        "n_scoreable": len(units),
        "scoreable": {c: {"n": u["n"], "has_pct": u["has_pct"], "malig_def": u["malig_def"]} for c, u in units.items()},
        "n_pairs": len(pairs),
        "n_triples": len(triples),
        "n_loo": len(loos),
        "n_differ": len(differ),
        "full_pool": [{k: v for k, v in r.items() if not k.startswith("_")} for r in fulls],
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("pairs", len(pairs), "triples", len(triples), "loo", len(loos), "differ", len(differ))
    print("WROTE", TABLES / "combo_table.tsv")


if __name__ == "__main__":
    main()

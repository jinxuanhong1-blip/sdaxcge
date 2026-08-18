#!/usr/bin/env python3
"""ADDITIVE CLDN4-only: can GSE200563 join the concordant-4 pool?

Solo: tumor-core CLDN4 vs same-patient TIME T/NK on the public
GSE200563_processed_data.txt.gz matrix. Also tumor-core Q4 vs Q1
IFN/MHC/TJ if n_units >= 8.

Join if the T/NK Spearman sign is negative (same as the concordant
sets) OR tumor-cell IFN/MHC is down in CLDN4-high. If the sign flips
and IFN/MHC is not down, or n is unusable, stop after the solo table.

Not Wu 2021 scRNA (that is GSE148071; excluded). This series is
Zhang/Abdo 2022 GeoMx DSP. No dual-high. No GSE148071.
"""
from __future__ import annotations

import json
import math
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_stats import random_effects_dl, spearman

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200563/"
    "suppl/GSE200563_processed_data.txt.gz"
)
MATRIX_NAME = "GSE200563_processed_data.txt.gz"
COL_RE = re.compile(r"^(TIME-L|TIME-B|TBME|mLN|LB|BC|L)(\d+)([a-z])?$")

TUMOR_COMPS = ("L", "LB")  # PanCK/tumor cores: primary lung + brain met
TIME_COMPS = ("TIME-L", "TIME-B")
TNK_GENES = ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "NKG7", "GNLY", "KLRD1"]
AUDIT_GENES = ["CLDN4", "TACSTD2", "EPCAM", "PTPRC"]
CONCORDANT = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
MIN_N_SPEARMAN = 4
MIN_N_Q4 = 8
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    return f"{x:.{nd}f}"


def ensure_matrix() -> Path:
    dest = DATA / MATRIX_NAME
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print("DOWNLOAD", GEO_URL)
    urllib.request.urlretrieve(GEO_URL, dest)
    return dest


def parse_roi(name: str) -> dict:
    m = COL_RE.match(str(name))
    if not m:
        raise ValueError(f"unparsed ROI column: {name}")
    return {
        "roi": str(name),
        "compartment": m.group(1),
        "patient": m.group(2),
        "rep": m.group(3) or "",
    }


def log2p1(s: pd.Series) -> pd.Series:
    return np.log2(pd.to_numeric(s, errors="coerce") + 1.0)


def mean_log_genes(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in mat.index]
    if not present:
        return pd.Series(np.nan, index=mat.columns)
    return mat.loc[present].apply(log2p1).mean(axis=0)


def q4_vs_q1(cldn4, y, min_n: int = MIN_N_Q4) -> dict:
    s = pd.DataFrame({"c": np.asarray(cldn4, float), "y": np.asarray(y, float)})
    s = s[np.isfinite(s["c"]) & np.isfinite(s["y"])].copy()
    n = int(len(s))
    out = {
        "n": n,
        "n_q1": np.nan,
        "n_q4": np.nan,
        "median_q1": np.nan,
        "median_q4": np.nan,
        "delta_median": np.nan,
        "r_rb": np.nan,
        "p": np.nan,
        "usable": False,
    }
    if n < min_n:
        return out
    ranks = s["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return out
    if qs.nunique() < 4:
        return out
    q1 = s.loc[qs == "Q1", "y"]
    q4 = s.loc[qs == "Q4", "y"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return out
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    out.update(
        {
            "n_q1": n1,
            "n_q4": n4,
            "median_q1": float(q1.median()),
            "median_q4": float(q4.median()),
            "delta_median": float(q4.median() - q1.median()),
            "r_rb": float(r_rb),
            "p": float(p),
            "usable": True,
        }
    )
    return out


def load_families() -> dict[str, list[str]]:
    a8 = json.loads((DATA / "a8_families.json").read_text())
    sets = a8["sets"]
    tj = set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    for g in a8.get("focal_genes", []):
        if g not in sets.get("KRT_EPITHELIAL", []):
            tj.add(g)
    tj.discard("CLDN4")
    return {
        "IFN": sorted(
            set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
            | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
        ),
        "MHC-I/APM": sorted(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "TJ": sorted(tj),
    }


def load_concordant_units() -> dict[str, dict]:
    units: dict[str, dict] = {}
    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    units["GSE131907"] = {
        "cohort": "GSE131907",
        "n": int(len(mal)),
        "mean": mal["mal_CLDN4_mean"].to_numpy(float),
        "tnk": mal["frac_tnk"].to_numpy(float),
        "note": "author Malignant; tumor sites; n_mal>=20; sample-level (Kim 2020)",
    }
    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    units["GSE205335"] = {
        "cohort": "GSE205335",
        "n": int(len(d)),
        "mean": d["mal_CLDN4_mean"].to_numpy(float),
        "tnk": d["frac_tnk"].to_numpy(float),
        "note": "author malignant, all subtypes",
    }
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"]
    units["GSE123902"] = {
        "cohort": "GSE123902",
        "n": int(len(el)),
        "mean": el["mal_CLDN4_mean"].to_numpy(float),
        "tnk": el["frac_tnk"].to_numpy(float),
        "note": "Laughney 2020 tumor/met donors; marker-malignant",
    }
    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"]
    units["GSE189357"] = {
        "cohort": "GSE189357",
        "n": int(len(el)),
        "mean": el["mal_CLDN4_mean"].to_numpy(float),
        "tnk": el["frac_tnk"].to_numpy(float),
        "note": "Zhu/Wang AIS–IAC; 9 patients; marker-malignant",
    }
    return units


def pool_units(members: list[dict]) -> dict:
    rhos, ps, ns = [], [], []
    for u in members:
        rho, p, n = spearman(u["mean"], u["tnk"])
        rhos.append(rho)
        ps.append(p)
        ns.append(n)
    re = random_effects_dl(rhos, ns)
    return {
        "k": len(members),
        "n": int(sum(ns)),
        "rho": re.get("pooled_rho", float("nan")),
        "p": re.get("p", float("nan")),
        "I2": re.get("I2", float("nan")),
        "ci95": re.get("ci95_rho", [float("nan"), float("nan")]),
        "member_rhos": {u["cohort"]: r for u, r in zip(members, rhos)},
        "member_ns": {u["cohort"]: n for u, n in zip(members, ns)},
        "member_ps": {u["cohort"]: p for u, p in zip(members, ps)},
    }


def load_geomx() -> tuple[pd.DataFrame, pd.DataFrame]:
    mat = pd.read_csv(ensure_matrix(), sep="\t", index_col=0)
    mat.index = mat.index.astype(str)
    mat = mat.groupby(mat.index).mean()
    meta = pd.DataFrame([parse_roi(c) for c in mat.columns])
    meta["is_tumor"] = meta["compartment"].isin(TUMOR_COMPS)
    meta["is_time"] = meta["compartment"].isin(TIME_COMPS)
    meta["is_mln"] = meta["compartment"] == "mLN"
    meta["is_control"] = meta["compartment"] == "BC"
    return mat, meta


def patient_table(mat: pd.DataFrame, meta: pd.DataFrame, fam: dict[str, list[str]]) -> pd.DataFrame:
    cldn4 = log2p1(mat.loc["CLDN4"])
    tacstd2 = log2p1(mat.loc["TACSTD2"]) if "TACSTD2" in mat.index else pd.Series(np.nan, index=mat.columns)
    tnk = mean_log_genes(mat, TNK_GENES)
    fam_scores = {name: mean_log_genes(mat, genes) for name, genes in fam.items()}
    present_tnk = [g for g in TNK_GENES if g in mat.index]
    rows = []
    for pid, g in meta.groupby("patient"):
        if g["is_control"].all():
            continue
        tumor = g[g["is_tumor"]]
        time = g[g["is_time"]]
        mln = g[g["is_mln"]]
        lung_t = g[g["compartment"] == "L"]
        lung_i = g[g["compartment"] == "TIME-L"]
        br_t = g[g["compartment"] == "LB"]
        br_i = g[g["compartment"] == "TIME-B"]
        rec = {
            "patient": pid,
            "n_tumor_roi": int(len(tumor)),
            "n_time_roi": int(len(time)),
            "n_mln_roi": int(len(mln)),
            "has_L": int(len(lung_t) > 0),
            "has_LB": int(len(br_t) > 0),
            "has_TIME_L": int(len(lung_i) > 0),
            "has_TIME_B": int(len(br_i) > 0),
            "cldn4_tumor": float(cldn4[tumor["roi"]].mean()) if len(tumor) else np.nan,
            "tacstd2_tumor": float(tacstd2[tumor["roi"]].mean()) if len(tumor) else np.nan,
            "tnk_tumor": float(tnk[tumor["roi"]].mean()) if len(tumor) else np.nan,
            "tnk_time": float(tnk[time["roi"]].mean()) if len(time) else np.nan,
            "cldn4_L": float(cldn4[lung_t["roi"]].mean()) if len(lung_t) else np.nan,
            "tnk_TIME_L": float(tnk[lung_i["roi"]].mean()) if len(lung_i) else np.nan,
            "cldn4_LB": float(cldn4[br_t["roi"]].mean()) if len(br_t) else np.nan,
            "tnk_TIME_B": float(tnk[br_i["roi"]].mean()) if len(br_i) else np.nan,
            "cldn4_mln": float(cldn4[mln["roi"]].mean()) if len(mln) else np.nan,
        }
        for name, sc in fam_scores.items():
            rec[f"{name}_tumor"] = float(sc[tumor["roi"]].mean()) if len(tumor) else np.nan
        rec["eligible_tnk"] = bool(np.isfinite(rec["cldn4_tumor"]) and np.isfinite(rec["tnk_time"]))
        rec["eligible_ifn"] = bool(np.isfinite(rec["cldn4_tumor"]))
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values("patient")
    out.attrs["tnk_genes_present"] = present_tnk
    return out


def contrast_row(name: str, x, y) -> dict:
    rho, p, n = spearman(x, y)
    q = q4_vs_q1(x, y)
    return {
        "contrast": name,
        "n": n,
        "rho": rho,
        "p": p,
        "n_q1": q["n_q1"],
        "n_q4": q["n_q4"],
        "q4_delta": q["delta_median"],
        "q4_r": q["r_rb"],
        "q4_p": q["p"],
        "q4_usable": q["usable"],
    }


def figures(pat: pd.DataFrame, solo: dict, join: bool) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    el = pat[pat["eligible_tnk"]].copy()
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.scatter(el["cldn4_tumor"], el["tnk_time"], c="#7a2d0b", s=42, zorder=3)
    for r in el.itertuples():
        ax.text(r.cldn4_tumor, r.tnk_time, r.patient, fontsize=7, ha="left", va="bottom", color="#444")
    ax.set_xlabel("Tumor-core CLDN4  log2(Q3+1)")
    ax.set_ylabel("Same-patient TIME T/NK score")
    ax.set_title(f"GSE200563 solo  n={solo['n']}  ρ={solo['rho']:.3f}  p={solo['p']:.3g}")
    ax.axhline(el["tnk_time"].median(), color="#bbb", lw=0.6)
    fig.tight_layout()
    fig.savefig(FIGS / "solo_cldn4_vs_tnk.png", dpi=160)
    fig.savefig(FIGS / "solo_cldn4_vs_tnk.pdf")
    plt.close(fig)

    fam = pat[pat["eligible_ifn"]].copy()
    ranks = fam["cldn4_tumor"].rank(method="average")
    fam["q"] = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.4), sharex=True)
    for ax, name in zip(axes, ["IFN", "MHC-I/APM", "TJ"]):
        col = f"{name}_tumor"
        parts = [fam.loc[fam["q"] == q, col].dropna().to_numpy() for q in ["Q1", "Q4"]]
        ax.boxplot(parts, tick_labels=["Q1", "Q4"], widths=0.55)
        ax.set_title(name)
        ax.set_ylabel("tumor-core family score" if name == "IFN" else "")
    fig.suptitle("GSE200563 tumor-core CLDN4 Q4 vs Q1 family scores", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "q4q1_ifn_mhc_tj.png", dpi=160)
    fig.savefig(FIGS / "q4q1_ifn_mhc_tj.pdf")
    plt.close(fig)


def write_finding(ctx: dict) -> None:
    solo = ctx["solo"]
    ifn = ctx["ifn_q"]
    mhc = ctx["mhc_q"]
    tj = ctx["tj_q"]
    base = ctx["base_pool"]
    new = ctx.get("new_pool")
    join = ctx["join"]
    reason = ctx["reason"]
    pat = ctx["patients"]
    n_time = int(pat["eligible_tnk"].sum())
    n_ifn = int(pat["eligible_ifn"].sum())
    n_catalog = int(len(pat))

    decision = "JOINS" if join else "does not join"
    lines = [
        "# Finding — GSE200563 CLDN4-only add-if-same-sign",
        "",
        "ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE148071.",
        "Public matrix only: `GSE200563_processed_data.txt.gz`.",
        "",
        "This accession is **not** Wu 2021 scRNA (that is GSE148071, excluded).",
        "GSE200563 is Zhang/Abdo 2022 GeoMx DSP of NSCLC primary + paired brain",
        "metastasis (PMID 36216799). PR #459 skipped it as “spatial, not scRNA”.",
        "This folder tests the **compartment analog**: tumor-core (L/LB) CLDN4",
        "vs same-patient TIME (TIME-L / TIME-B) T/NK. Honest n = patients.",
        "",
        f"## Decision: GSE200563 **{decision}** the concordant pool",
        "",
        reason,
        "",
        "## Solo table (this accession)",
        "",
        "| item | n / effect | note |",
        "|---|---|---|",
        f"| catalog patients (exclude BC) | **{n_catalog}** | numeric IDs on ROI names |",
        f"| patients with L or LB tumor core | **{n_ifn}** | epithelial/malignant analog |",
        f"| patients with tumor core **and** TIME ROI | **{n_time}** | Spearman unit |",
        f"| CLDN4 present | yes | Q3 range {ctx['cldn4_q3_min']:.1f}–{ctx['cldn4_q3_max']:.1f} |",
        f"| TACSTD2 present | yes | audit only; never a gate |",
        f"| T/NK genes used | {', '.join(ctx['tnk_present'])} | mean log2(Q3+1) |",
        f"| tumor CLDN4 vs TIME T/NK Spearman | n={solo['n']} ρ={_fmt(solo['rho'])} p={_fmt(solo['p'], 3)} | primary solo |",
        (
            f"| same, Q4 vs Q1 T/NK | n_Q1/Q4={int(solo['n_q1'])}/{int(solo['n_q4'])} "
            f"r={_fmt(solo['q4_r'])} p={_fmt(solo['q4_p'], 3)} | same 16 patients; tails thin |"
            if solo["q4_usable"]
            else f"| same, Q4 vs Q1 T/NK | n={solo['n']} < 8 or thin tails | not used |"
        ),
        "",
        "Patient is the unit. Tumor-core = L (primary) and/or LB (brain met).",
        "TIME = CD45/immune AOIs labeled TIME-L / TIME-B. mLN and TBME are",
        "sensitivity only. BC controls are dropped. Values are depositor Q3.",
        "",
        "## Per-patient TIME-paired units",
        "",
        "| patient | L | LB | TIME-L | TIME-B | CLDN4 tumor | T/NK TIME | IFN tumor | MHC tumor | TJ tumor |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    el = pat[pat["eligible_tnk"]].sort_values("patient")
    for _, r in el.iterrows():
        lines.append(
            f"| {r.patient} | {int(r.has_L)} | {int(r.has_LB)} | {int(r.has_TIME_L)} | {int(r.has_TIME_B)} | "
            f"{r.cldn4_tumor:.3f} | {r.tnk_time:.3f} | {r['IFN_tumor']:.3f} | "
            f"{r['MHC-I/APM_tumor']:.3f} | {r['TJ_tumor']:.3f} |"
        )
    lines += [
        "",
        f"Honest Spearman n = **{n_time}** (not 35 cases, not 120 ROIs, not 109 paper ROIs).",
        "",
        "## Sensitivity T/NK (not the join key unless primary n is unusable)",
        "",
        "| contrast | n | ρ | p |",
        "|---|---:|---:|---:|",
    ]
    for r in ctx["sens"]:
        lines.append(f"| {r['contrast']} | {r['n']} | {_fmt(r['rho'])} | {_fmt(r['p'], 3)} |")
    lines += [
        "",
        "## Tumor-core IFN / MHC / TJ  (Q4 vs Q1 if n_units ≥ 8)",
        "",
        f"Units = patients with L or LB (n={n_ifn}). Family score = mean log2(Q3+1)",
        "of A8 genes present on the matrix. TJ holds CLDN4 out. Positive delta =",
        "higher in CLDN4-high.",
        "",
        "| family | n | n_Q1 / n_Q4 | Δ median (Q4−Q1) | rank-biserial r | p | Spearman vs CLDN4 ρ (p) |",
        "|---|---:|---|---:|---:|---:|---|",
    ]
    for name, q, cont in (
        ("IFN", ifn, ctx["ifn_s"]),
        ("MHC-I/APM", mhc, ctx["mhc_s"]),
        ("TJ (CLDN4 held out)", tj, ctx["tj_s"]),
    ):
        nq = f"{int(q['n_q1'])}/{int(q['n_q4'])}" if q["usable"] else "—"
        lines.append(
            f"| {name} | {q['n']} | {nq} | {_fmt(q['delta_median'])} | {_fmt(q['r_rb'])} | "
            f"{_fmt(q['p'], 3)} | {_fmt(cont['rho'])} ({_fmt(cont['p'], 3)}) |"
        )
    lines += [
        "",
        "IFN/MHC **down** in CLDN4-high means Q4−Q1 delta < 0 or continuous ρ < 0.",
        "",
        "## Concordant-4 locked pool (PR #459; not re-audited as singles)",
        "",
        "Members: GSE123902 + GSE131907 + GSE205335 + GSE189357.",
        "Score = malignant CLDN4 **mean** vs same-unit T/NK fraction.",
        "Pooled ρ = DerSimonian–Laird Fisher-z (same as PR #459).",
        "",
        "| cohort | n | ρ | p |",
        "|---|---:|---:|---:|",
    ]
    for c in CONCORDANT:
        lines.append(
            f"| {c} | {base['member_ns'][c]} | {_fmt(base['member_rhos'][c])} | {_fmt(base['member_ps'][c], 3)} |"
        )
    lines.append(
        f"| **concordant-4** | **{base['n']}** | **{_fmt(base['rho'])}** | {_fmt(base['p'], 3)} "
        f"(I²={_fmt(base['I2'], 0)}%) |"
    )
    if join and new is not None:
        lines += [
            "",
            "## Stacked pool after adding GSE200563",
            "",
            f"| pool | k | N | ρ | p | I² |",
            "|---|---:|---:|---:|---:|---:|",
            f"| concordant-4 | 4 | {base['n']} | {_fmt(base['rho'])} | {_fmt(base['p'], 3)} | {_fmt(base['I2'], 0)}% |",
            f"| concordant-4 **+ GSE200563** | 5 | **{new['n']}** | **{_fmt(new['rho'])}** | {_fmt(new['p'], 3)} | {_fmt(new['I2'], 0)}% |",
            "",
            f"GSE200563 member: n={solo['n']} ρ={_fmt(solo['rho'])} p={_fmt(solo['p'], 3)} "
            "(tumor-core CLDN4 vs TIME T/NK; GeoMx analog, not a cell fraction).",
        ]
    else:
        lines += [
            "",
            "## Stacked pool",
            "",
            "Not computed. Join rule failed (sign flip / IFN-MHC not down / n unusable).",
        ]
    lines += [
        "",
        "## Methods (locked)",
        "",
        "- CLDN4 only. TACSTD2 is an audit gene, never a gate.",
        "- Matrix: GEO `GSE200563_processed_data.txt.gz` only (depositor Q3).",
        "- Tumor-core ROIs = columns `L##` and `LB##`. TIME ROIs = `TIME-L##`, `TIME-B##`.",
        "- Patient ID = the numeric suffix. Replicates (a/b) are averaged.",
        "- CLDN4 score = mean log2(Q3+1) across a patient's tumor-core ROIs.",
        "- T/NK score = mean log2(Q3+1) of CD3D/CD3E/CD3G/CD2/CD8A/CD8B/NKG7/GNLY/KLRD1",
        "  on TIME ROIs (same genes as the PR #459 marker gate, scored as a mean).",
        "- Spearman requires n≥4 finite pairs. Q4 vs Q1 requires n_units≥8 and two tails.",
        "- Family scores use A8 Hallmark IFN-α/γ, custom MHC-I/APM, KEGG+GOBP TJ (CLDN4 held out).",
        "- Join if solo T/NK ρ < 0 **or** tumor IFN or MHC Q4−Q1 delta < 0 (n≥8).",
        "- Stop if n<4 for T/NK **and** IFN/MHC Q4 is unusable, or if T/NK ρ>0 and IFN/MHC not down.",
        "- Stacked ρ = Fisher-z DerSimonian–Laird of cohort Spearmans. Honest n = sum of units.",
        "- Concordant-4 vectors are the PR #459 locked tables. GSE148071 is not used.",
        "- Reproduce: `python3 methods/gse200563_cldn4_add/analyze.py`",
        "",
        "## How to read this",
        "",
        "- This is a **GeoMx patient analog**, not a single-cell T/NK fraction.",
        "- Do not quote n=35 cases or n=120 ROIs as the Spearman n.",
        "- Do not call this Wu scRNA. Do not add GSE148071.",
        "- No dual-high. No CellChat.",
        "",
    ]
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    fam = load_families()
    mat, meta = load_geomx()
    present = {g: g in mat.index for g in AUDIT_GENES + TNK_GENES}
    if not present["CLDN4"] or not present["TACSTD2"]:
        raise SystemExit("CLDN4 or TACSTD2 missing from public matrix")
    tnk_present = [g for g in TNK_GENES if g in mat.index]
    if len(tnk_present) < 3:
        raise SystemExit("too few T/NK genes on the matrix")

    pat = patient_table(mat, meta, fam)
    el = pat[pat["eligible_tnk"]]
    ifn_units = pat[pat["eligible_ifn"]]

    solo_rho, solo_p, solo_n = spearman(el["cldn4_tumor"], el["tnk_time"])
    solo_q = q4_vs_q1(el["cldn4_tumor"], el["tnk_time"])
    solo = {
        "n": solo_n,
        "rho": solo_rho,
        "p": solo_p,
        "n_q1": solo_q["n_q1"],
        "n_q4": solo_q["n_q4"],
        "q4_r": solo_q["r_rb"],
        "q4_p": solo_q["p"],
        "q4_usable": solo_q["usable"],
    }

    c_mln, y_mln = [], []
    for r in pat.itertuples():
        if not r.eligible_tnk:
            continue
        vals = [v for v in [r.cldn4_tumor, r.cldn4_mln] if np.isfinite(v)]
        if vals:
            c_mln.append(float(np.mean(vals)))
            y_mln.append(r.tnk_time)
    sens = [
        contrast_row("L vs TIME-L (lung pair)", pat["cldn4_L"], pat["tnk_TIME_L"]),
        contrast_row("LB vs TIME-B (brain-met pair)", pat["cldn4_LB"], pat["tnk_TIME_B"]),
        contrast_row("tumor-core CLDN4 vs tumor-core T/NK (same ROI)", ifn_units["cldn4_tumor"], ifn_units["tnk_tumor"]),
        contrast_row("tumor+mLN CLDN4 vs TIME T/NK", c_mln, y_mln),
    ]

    ifn_q = q4_vs_q1(ifn_units["cldn4_tumor"], ifn_units["IFN_tumor"])
    mhc_q = q4_vs_q1(ifn_units["cldn4_tumor"], ifn_units["MHC-I/APM_tumor"])
    tj_q = q4_vs_q1(ifn_units["cldn4_tumor"], ifn_units["TJ_tumor"])
    ifn_s = contrast_row("CLDN4 vs IFN", ifn_units["cldn4_tumor"], ifn_units["IFN_tumor"])
    mhc_s = contrast_row("CLDN4 vs MHC", ifn_units["cldn4_tumor"], ifn_units["MHC-I/APM_tumor"])
    tj_s = contrast_row("CLDN4 vs TJ", ifn_units["cldn4_tumor"], ifn_units["TJ_tumor"])

    tnk_usable = solo_n >= MIN_N_SPEARMAN and math.isfinite(solo_rho)
    ifn_down = (ifn_q["usable"] and ifn_q["delta_median"] < 0) or (
        math.isfinite(ifn_s["rho"]) and ifn_s["rho"] < 0 and ifn_s["n"] >= MIN_N_Q4
    )
    mhc_down = (mhc_q["usable"] and mhc_q["delta_median"] < 0) or (
        math.isfinite(mhc_s["rho"]) and mhc_s["rho"] < 0 and mhc_s["n"] >= MIN_N_Q4
    )
    ifn_mhc_down = bool(ifn_down or mhc_down)
    tnk_neg = bool(tnk_usable and solo_rho < 0)

    units = load_concordant_units()
    base = pool_units([units[c] for c in CONCORDANT])

    if not tnk_usable and not (ifn_q["usable"] or mhc_q["usable"]):
        join = False
        reason = (
            f"n is unusable: TIME-paired Spearman n={solo_n} (need ≥{MIN_N_SPEARMAN}) "
            f"and tumor-core Q4 vs Q1 n={ifn_q['n']} (need ≥{MIN_N_Q4}). "
            "Stop after the solo table. GSE200563 does not join."
        )
        new = None
    elif tnk_usable and (not tnk_neg) and (not ifn_mhc_down):
        join = False
        reason = (
            f"T/NK Spearman sign **flips** relative to the concordant sets "
            f"(ρ={solo_rho:+.3f}, n={solo_n}, p={solo_p:.3g}) and tumor-cell "
            f"IFN/MHC is **not** down in CLDN4-high "
            f"(IFN Δ={_fmt(ifn_q['delta_median'])}, MHC Δ={_fmt(mhc_q['delta_median'])}). "
            "Stop after the solo table. GSE200563 does not join."
        )
        new = None
    else:
        join = True
        why = []
        if tnk_neg:
            why.append(f"T/NK Spearman is negative (ρ={solo_rho:.3f}, n={solo_n}, same sign as concordant-4)")
        if ifn_mhc_down:
            why.append(
                f"tumor-cell IFN/MHC is down in CLDN4-high "
                f"(IFN Δ={_fmt(ifn_q['delta_median'])}, MHC Δ={_fmt(mhc_q['delta_median'])})"
            )
        reason = (
            "Join rule hit: " + "; ".join(why) + ". "
            "Stack GSE200563 onto GSE123902+GSE131907+GSE205335+GSE189357."
        )
        geo_unit = {
            "cohort": "GSE200563",
            "n": solo_n,
            "mean": el["cldn4_tumor"].to_numpy(float),
            "tnk": el["tnk_time"].to_numpy(float),
            "note": "GeoMx tumor-core CLDN4 vs TIME T/NK",
        }
        new = pool_units([units[c] for c in CONCORDANT] + [geo_unit])

    cldn4_raw = pd.to_numeric(mat.loc["CLDN4"], errors="coerce")
    ctx = {
        "solo": solo,
        "sens": sens,
        "ifn_q": ifn_q,
        "mhc_q": mhc_q,
        "tj_q": tj_q,
        "ifn_s": ifn_s,
        "mhc_s": mhc_s,
        "tj_s": tj_s,
        "base_pool": base,
        "new_pool": new,
        "join": join,
        "reason": reason,
        "patients": pat,
        "tnk_present": tnk_present,
        "cldn4_q3_min": float(cldn4_raw.min()),
        "cldn4_q3_max": float(cldn4_raw.max()),
        "present": present,
    }
    figures(pat, solo, join)
    write_finding(ctx)

    pat.to_csv(TABLES / "patient_units.tsv", sep="\t", index=False)
    meta.to_csv(TABLES / "roi_annotation.tsv", sep="\t", index=False)
    pd.DataFrame([solo]).to_csv(TABLES / "solo_tnk.tsv", sep="\t", index=False)
    pd.DataFrame(sens).to_csv(TABLES / "sensitivity_tnk.tsv", sep="\t", index=False)
    fam_rows = []
    for name, q, cont in (
        ("IFN", ifn_q, ifn_s),
        ("MHC-I/APM", mhc_q, mhc_s),
        ("TJ", tj_q, tj_s),
    ):
        fam_rows.append({"family": name, **q, "spearman_rho": cont["rho"], "spearman_p": cont["p"]})
    pd.DataFrame(fam_rows).to_csv(TABLES / "family_q4q1.tsv", sep="\t", index=False)
    summary = {
        "dataset": "GSE200563",
        "platform": "GeoMx DSP Q3 (Zhang/Abdo 2022), not Wu scRNA",
        "join": join,
        "reason": reason,
        "solo": solo,
        "concordant4": {k: v for k, v in base.items() if k != "ci95"} | {"ci95": base["ci95"]},
        "stacked": None
        if new is None
        else {k: v for k, v in new.items() if k != "ci95"} | {"ci95": new["ci95"]},
        "n_catalog": int(len(pat)),
        "n_tumor": int(pat["eligible_ifn"].sum()),
        "n_time_paired": int(pat["eligible_tnk"].sum()),
        "tnk_genes": tnk_present,
        "cldn4_present": True,
        "tacstd2_present": True,
        "dual_high": False,
        "gse148071": False,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print("JOIN" if join else "NO-JOIN", reason)
    print("SOLO", solo)
    print("BASE", base["n"], base["rho"], base["p"])
    if new:
        print("NEW", new["n"], new["rho"], new["p"])
    print("WROTE", HERE / "FINDING.md")


if __name__ == "__main__":
    main()

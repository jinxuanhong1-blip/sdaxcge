#!/usr/bin/env python3
"""ADDITIVE extra — prerank GSEA on malignant-cell patient pseudobulk.

Does NOT re-run TCGA-LUAD/LUSC A8. User A8 (keratin/TJ up, Hallmark EMT down)
is taken as given for TCGA. This slice: public lung tumor scRNA → malignant
UMI-sum per patient → TACSTD2-high vs low → same GSEA engine → meta NES.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from gsea_core import (  # noqa: E402
    bh_fdr,
    gsea_prerank,
    median_split,
    quartile_split,
    rank_high_vs_low,
    rank_spearman_vs_target,
)

OUT = ROOT / "results" / "scrna_pseudobulk_gsea_meta"
FIG = OUT / "figures"
TAB = OUT / "tables"
PB = OUT / "pseudobulk"
SET_JSON = ROOT / "data" / "genesets" / "a8_sets.json"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

TARGET = "TACSTD2"
MIN_ARM_GSEA = 3
MIN_ARM_META = 5
MIN_ARM_QUARTILE = 6
PRIMARY_FDR = 0.05

# A8 freeze primary 12 + IFN + MHC-I (this extra only)
EXTRA_PRIMARY = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
]
HEADLINE = [
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "KEGG_TIGHT_JUNCTION",
    "HALLMARK_APICAL_JUNCTION",
    "GOBP_KERATINIZATION",
    "KRT_EPITHELIAL",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
]
COHORT_ORDER = [
    "GSE207422",
    "GSE241934",
    "GSE131907",
    "GSE205335",
    "GSE291670",
    "GSE253013",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1)


def verdict_from_primary(gsea: pd.DataFrame) -> str:
    if gsea is None or gsea.empty:
        return "not_run"
    prim = gsea[gsea["primary"] == True].copy()  # noqa: E712
    if prim.empty:
        return "not_run"

    def sig(term: str, want_up: bool) -> bool:
        row = prim[prim["term"] == term]
        if row.empty:
            return False
        nes = float(row.iloc[0]["nes"])
        q = float(row.iloc[0]["fdr_primary"])
        if q >= PRIMARY_FDR or nes != nes:
            return False
        return (nes > 0) if want_up else (nes < 0)

    def any_bucket(bucket: str, want_up: bool) -> bool:
        ok = False
        for _, r in prim[prim["bucket"] == bucket].iterrows():
            if r["term"] == "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION":
                continue
            if r["fdr_primary"] < PRIMARY_FDR and ((r["nes"] > 0) if want_up else (r["nes"] < 0)):
                ok = True
        return ok

    emt_down = sig("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", False)
    emt_up = sig("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", True)
    tj_up = any_bucket("TJ", True)
    krt_up = any_bucket("KERATIN_BARRIER", True)
    tj_down = any_bucket("TJ", False) and not tj_up
    krt_down = any_bucket("KERATIN_BARRIER", False) and not krt_up
    if emt_down and tj_up and krt_up:
        return "supportive"
    if tj_up and krt_up and emt_up:
        return "keratin_TJ_up_Hallmark_EMT_opposite"
    if tj_up and krt_up and not emt_down and not emt_up:
        return "keratin_TJ_up_Hallmark_EMT_null"
    if emt_down and (tj_up or krt_up):
        return "partial"
    if tj_down and krt_down:
        return "contradicts_TJ_KRT"
    if tj_up or krt_up or emt_down:
        return "mixed"
    return "null"


def annotate(gsea: pd.DataFrame, set_meta: dict, contrast: str, cohort: str, n_high: int, n_low: int, primary_names: set[str]) -> pd.DataFrame:
    if gsea.empty:
        return gsea
    gsea = gsea.copy()
    gsea["cohort"] = cohort
    gsea["contrast"] = contrast
    gsea["n_high"] = n_high
    gsea["n_low"] = n_low
    meta_df = pd.DataFrame.from_dict(set_meta, orient="index").rename_axis("term").reset_index()
    gsea = gsea.merge(meta_df, on="term", how="left")
    gsea["primary"] = gsea["term"].isin(primary_names)
    gsea["fdr_all"] = bh_fdr(gsea["nom_p"])
    gsea["fdr_primary"] = np.nan
    prim = gsea["primary"]
    if prim.any():
        gsea.loc[prim, "fdr_primary"] = bh_fdr(gsea.loc[prim, "nom_p"])
    return gsea


def fmt_p(p) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_nes(x) -> str:
    if x != x:
        return "NA"
    return f"{x:+.3f}"


def stouffer_signs(nes: pd.Series) -> tuple[float, float]:
    s = nes.dropna()
    if s.empty:
        return np.nan, np.nan
    z = float(np.sign(s).sum() / np.sqrt(len(s)))
    p = float(2 * stats.norm.sf(abs(z)))
    return z, p


def load_expr(cohort: str) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    cpath = PB / f"{cohort}_counts.tsv.gz"
    mpath = PB / f"{cohort}_patient_meta.tsv"
    if not cpath.exists():
        return None, pd.read_csv(mpath, sep="\t") if mpath.exists() else None
    counts = pd.read_csv(cpath, sep="\t", index_col=0)
    meta = pd.read_csv(mpath, sep="\t") if mpath.exists() else None
    expr = log2_cpm(counts)
    return expr, meta


def main() -> None:
    payload = json.loads(SET_JSON.read_text())
    all_sets = payload["sets"]
    set_meta = payload["meta"]
    a8_primary = [k for k, v in set_meta.items() if v.get("primary")]
    primary_names = list(dict.fromkeys(a8_primary + EXTRA_PRIMARY))
    primary_sets = {k: all_sets[k] for k in primary_names if k in all_sets}
    for term in primary_names:
        members = set(primary_sets.get(term, []))
        if TARGET in members:
            raise SystemExit(f"{TARGET} is inside {term} — ranking is circular")

    gsea_rows = []
    verdicts = []
    sample_rows = []
    inventory = []

    for cohort in COHORT_ORDER:
        log(f"=== {cohort} ===")
        expr, meta = load_expr(cohort)
        if expr is None or expr.empty or TARGET not in expr.index:
            inventory.append(
                {
                    "cohort": cohort,
                    "usable": False,
                    "reason": "no_pseudobulk_or_no_TACSTD2" if expr is None or expr.empty else "TACSTD2_absent",
                    "n": 0 if expr is None or expr.empty else int(expr.shape[1]),
                }
            )
            verdicts.append({"cohort": cohort, "contrast": "tacstd2_median", "verdict": "not_run", "n_high": 0, "n_low": 0})
            continue
        tac = expr.loc[TARGET].dropna()
        expr = expr.loc[:, tac.index]
        n = int(tac.shape[0])
        for sid in tac.index:
            sample_rows.append({"cohort": cohort, "patient": sid, "TACSTD2": float(tac.loc[sid])})
        high, low = median_split(tac)
        n_h, n_l = len(high), len(low)
        log(f"  n={n} genes={expr.shape[0]} TACSTD2 median={float(tac.median()):.3f} high={n_h} low={n_l}")
        inventory.append(
            {
                "cohort": cohort,
                "usable": n_h >= MIN_ARM_GSEA and n_l >= MIN_ARM_GSEA,
                "reason": "ok" if n_h >= MIN_ARM_GSEA and n_l >= MIN_ARM_GSEA else "n_arm<3",
                "n": n,
                "n_high": n_h,
                "n_low": n_l,
                "n_genes": int(expr.shape[0]),
                "tacstd2_median": float(tac.median()),
                "meta_eligible": n_h >= MIN_ARM_META and n_l >= MIN_ARM_META,
            }
        )
        if n_h < MIN_ARM_GSEA or n_l < MIN_ARM_GSEA:
            verdicts.append(
                {
                    "cohort": cohort,
                    "contrast": "tacstd2_median",
                    "verdict": "n_too_small",
                    "n_high": n_h,
                    "n_low": n_l,
                    "underpowered": True,
                }
            )
            continue
        g = annotate(
            gsea_prerank(rank_high_vs_low(expr, high, low), primary_sets),
            set_meta,
            "tacstd2_median",
            cohort,
            n_h,
            n_l,
            set(primary_names),
        )
        gsea_rows.append(g)
        v = verdict_from_primary(g)
        if n_h < MIN_ARM_META or n_l < MIN_ARM_META:
            v = "n_too_small"
        verdicts.append(
            {
                "cohort": cohort,
                "contrast": "tacstd2_median",
                "verdict": v,
                "n_high": n_h,
                "n_low": n_l,
                "underpowered": n_h < MIN_ARM_META or n_l < MIN_ARM_META,
            }
        )
        log(f"  median verdict={v}")

        qh, ql = quartile_split(tac)
        if len(qh) >= MIN_ARM_QUARTILE and len(ql) >= MIN_ARM_QUARTILE:
            gq = annotate(
                gsea_prerank(rank_high_vs_low(expr, qh, ql), primary_sets),
                set_meta,
                "tacstd2_quartile",
                cohort,
                len(qh),
                len(ql),
                set(primary_names),
            )
            gsea_rows.append(gq)
            verdicts.append(
                {
                    "cohort": cohort,
                    "contrast": "tacstd2_quartile",
                    "verdict": verdict_from_primary(gq),
                    "n_high": len(qh),
                    "n_low": len(ql),
                    "underpowered": False,
                }
            )

        rank_s = rank_spearman_vs_target(expr, tac)
        if len(rank_s) >= 50:
            gs = annotate(gsea_prerank(rank_s, primary_sets), set_meta, "tacstd2_spearman", cohort, n, n, set(primary_names))
            gsea_rows.append(gs)

    gsea = pd.concat(gsea_rows, ignore_index=True) if gsea_rows else pd.DataFrame()
    verd = pd.DataFrame(verdicts)
    inv = pd.DataFrame(inventory)
    samp = pd.DataFrame(sample_rows)
    if not gsea.empty:
        gsea.to_csv(TAB / "gsea_prerank_all.tsv", sep="\t", index=False)
        gsea.loc[gsea.contrast == "tacstd2_median"].to_csv(TAB / "gsea_tacstd2_median.tsv", sep="\t", index=False)
    verd.to_csv(TAB / "verdicts.tsv", sep="\t", index=False)
    inv.to_csv(TAB / "gsea_cohort_n.tsv", sep="\t", index=False)
    samp.to_csv(TAB / "patient_tacstd2.tsv", sep="\t", index=False)

    # Meta on median-split, meta-eligible cohorts
    meta_rows = []
    if not gsea.empty:
        med = gsea.loc[(gsea.contrast == "tacstd2_median") & (gsea.primary == True)].copy()  # noqa: E712
        eligible = set(inv.loc[inv.get("meta_eligible", False) == True, "cohort"])  # noqa: E712
        med_e = med.loc[med.cohort.isin(eligible)]
        for term, sub in med_e.groupby("term"):
            nes = sub["nes"].astype(float)
            z, p = stouffer_signs(nes)
            meta_rows.append(
                {
                    "term": term,
                    "bucket": sub["bucket"].iloc[0] if "bucket" in sub else "",
                    "n_cohorts": int(nes.notna().sum()),
                    "cohorts": ",".join(sub["cohort"].tolist()),
                    "median_NES": float(nes.median()) if nes.notna().any() else np.nan,
                    "mean_NES": float(nes.mean()) if nes.notna().any() else np.nan,
                    "n_pos": int((nes > 0).sum()),
                    "n_neg": int((nes < 0).sum()),
                    "stouffer_z_signs": z,
                    "stouffer_p_signs": p,
                    "n_high_sum": int(sub["n_high"].sum()),
                    "n_low_sum": int(sub["n_low"].sum()),
                }
            )
    meta = pd.DataFrame(meta_rows)
    if not meta.empty:
        meta["stouffer_fdr"] = bh_fdr(meta["stouffer_p_signs"])
        meta = meta.sort_values("term")
        meta.to_csv(TAB / "meta_nes.tsv", sep="\t", index=False)

    _figures(gsea, meta, inv, verd)
    _report(gsea, meta, inv, verd, primary_names)
    summary = {
        "n_cohorts_gsea": int(inv["usable"].sum()) if not inv.empty and "usable" in inv else 0,
        "n_cohorts_meta": int(inv["meta_eligible"].sum()) if not inv.empty and "meta_eligible" in inv else 0,
        "primary_sets": primary_names,
        "min_arm_gsea": MIN_ARM_GSEA,
        "min_arm_meta": MIN_ARM_META,
        "note": "TCGA A8 was not re-run. Additive scRNA patient-pseudobulk extra only.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    log("done")


def _figures(gsea: pd.DataFrame, meta: pd.DataFrame, inv: pd.DataFrame, verd: pd.DataFrame) -> None:
    if gsea.empty:
        return
    med = gsea.loc[(gsea.contrast == "tacstd2_median") & (gsea.primary == True)].copy()  # noqa: E712
    if med.empty:
        return
    terms = [t for t in HEADLINE if t in set(med.term)]
    extra = [t for t in sorted(med.term.unique()) if t not in terms]
    terms = terms + extra
    cohorts = [c for c in COHORT_ORDER if c in set(med.cohort)]
    mat = pd.DataFrame(np.nan, index=terms, columns=cohorts)
    qmat = pd.DataFrame(np.nan, index=terms, columns=cohorts)
    for _, r in med.iterrows():
        if r["term"] in mat.index and r["cohort"] in mat.columns:
            mat.loc[r["term"], r["cohort"]] = r["nes"]
            qmat.loc[r["term"], r["cohort"]] = r["fdr_primary"]

    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    vmax = np.nanmax(np.abs(mat.to_numpy())) if np.isfinite(mat.to_numpy()).any() else 1
    vmax = max(float(vmax), 1.0)
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(cohorts, rotation=35, ha="right")
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms, fontsize=8)
    nlab = []
    for c in cohorts:
        row = inv.loc[inv.cohort == c]
        if row.empty:
            nlab.append(c)
        else:
            nh = int(row.iloc[0].get("n_high", 0) or 0)
            nl = int(row.iloc[0].get("n_low", 0) or 0)
            nlab.append(f"{c}\n{nh} vs {nl}")
    ax.set_xticklabels(nlab, rotation=0, ha="center", fontsize=8)
    for i, t in enumerate(terms):
        for j, c in enumerate(cohorts):
            nes = mat.loc[t, c]
            q = qmat.loc[t, c]
            if nes != nes:
                txt = "NA"
            else:
                star = "*" if q == q and q < 0.05 else ""
                txt = f"{nes:+.2f}{star}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.5, color="black")
    ax.set_title("Malignant-cell patient pseudobulk GSEA (TACSTD2 median)\nNES; * BH-FDR<0.05 within 15 primary sets. TCGA A8 not re-run.")
    fig.colorbar(im, ax=ax, fraction=0.03, label="NES (high vs low)")
    fig.tight_layout()
    fig.savefig(FIG / "nes_heatmap_median.png", dpi=160)
    fig.savefig(FIG / "nes_heatmap_median.pdf")
    plt.close(fig)

    # Headline bars + meta
    if meta.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    show = [t for t in HEADLINE if t in set(meta.term)]
    y = np.arange(len(show))
    medn = [float(meta.loc[meta.term == t, "median_NES"].iloc[0]) for t in show]
    colors = ["#b2182b" if v > 0 else "#2166ac" for v in medn]
    ax.barh(y, medn, color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(show, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Median NES across meta-eligible cohorts")
    ncoh = int(meta["n_cohorts"].iloc[0]) if not meta.empty else 0
    ax.set_title(f"Meta median NES (k={ncoh} cohorts, both arms ≥5)")
    for i, t in enumerate(show):
        r = meta.loc[meta.term == t].iloc[0]
        ax.text(
            medn[i] + (0.04 if medn[i] >= 0 else -0.04),
            i,
            f" {r['n_pos']}/{r['n_cohorts']} up  Stouffer p={fmt_p(r['stouffer_p_signs'])}",
            va="center",
            ha="left" if medn[i] >= 0 else "right",
            fontsize=7,
        )

    ax = axes[1]
    # per-cohort Hallmark EMT / KEGG TJ / KRT / IFNG / MHC-I
    focus = [
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "KEGG_TIGHT_JUNCTION",
        "GOBP_KERATINIZATION",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    ]
    x = np.arange(len(cohorts))
    width = 0.16
    for i, term in enumerate(focus):
        vals = [mat.loc[term, c] if term in mat.index else np.nan for c in cohorts]
        ax.bar(x + (i - 2) * width, vals, width=width, label=term.replace("HALLMARK_", "H_").replace("CUSTOM_", ""), linewidth=0.3, edgecolor="black")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(nlab, fontsize=7)
    ax.set_ylabel("NES")
    ax.set_title("Per-cohort NES (honest n on x-axis)")
    ax.legend(fontsize=6, loc="best")
    fig.tight_layout()
    fig.savefig(FIG / "nes_meta_and_headline.png", dpi=160)
    fig.savefig(FIG / "nes_meta_and_headline.pdf")
    plt.close(fig)


def _report(gsea: pd.DataFrame, meta: pd.DataFrame, inv: pd.DataFrame, verd: pd.DataFrame, primary_names: list[str]) -> None:
    lines = []
    lines.append("# ADDITIVE extra — malignant-cell patient pseudobulk GSEA (TACSTD2-high vs low)")
    lines.append("")
    lines.append("**Public data only. Additive to the TCGA A8 slide. That TCGA audit was not re-run.**")
    lines.append("")
    lines.append(
        "User A8 (TROP2-high enrich keratin/TJ, Hallmark EMT down) is **taken as given** for TCGA-LUAD/LUSC. "
        "This file asks the same prerank question after collapsing **malignant cells to one patient** in public lung tumor scRNA."
    )
    lines.append("")
    lines.append("## 一句话结论 / TL;DR")
    lines.append("")
    tldr = _tldr(gsea, meta, inv, verd)
    lines.append(tldr)
    lines.append("")
    lines.append("## Why this extra exists")
    lines.append("")
    lines.append(
        "A8 is TCGA bulk. Conservation of keratin/TJ-high (and Hallmark EMT-low) inside **malignant epithelium**, "
        "after removing stroma/immune mixture, is a different question. Patient is the unit. Cell-level tests are not reported as GSEA."
    )
    lines.append("")
    lines.append("## Cohorts")
    lines.append("")
    lines.append("| Cohort | n patients (gated) | n high vs low | Meta-eligible | Note |")
    lines.append("|---|---:|---|---|---|")
    notes = {
        "GSE207422": "Hu 2023 neoadj PD-1+chemo; marker malignant-like (no public CopyKAT table)",
        "GSE241934": "NEOTIDE + real-world neoadj IO; author tumor Epi (IIT+RWC)",
        "GSE253013": "Sze/Xiang 2024 treatment-naive LUAD; 9.3 GB RDS not loaded; n would be < meta gate",
        "GSE131907": "Kim 2020 atlas; author Malignant cells; drop nLung/nLN",
        "GSE291670": "Neoadj anlotinib+camrelizumab; 6 tumors; marker malignant-like",
        "GSE205335": "Park/Ahn/Lee ICI atlas; author Malignant cells; patient from SOFT; drop normal tissues",
    }
    for c in COHORT_ORDER:
        row = inv.loc[inv.cohort == c]
        if row.empty:
            lines.append(f"| {c} | 0 | — | no | {notes.get(c, '')} |")
            continue
        r = row.iloc[0]
        n = int(r.get("n", 0) or 0)
        nh = r.get("n_high", "")
        nl = r.get("n_low", "")
        elig = "yes" if r.get("meta_eligible") else "no"
        hl = f"{nh} vs {nl}" if nh == nh and nh != "" else "—"
        lines.append(f"| {c} | {n} | {hl} | {elig} | {notes.get(c, '')} |")
    lines.append("")
    lines.append("## Pre-specified design")
    lines.append("")
    lines.append("| Piece | Choice | Honest limitation |")
    lines.append("|---|---|---|")
    lines.append("| Split | TACSTD2 median, Welch t | Quartile only if both arms ≥6 |")
    lines.append("| Unit | Patient malignant UMI-sum, log2(CPM+1) | Not cell-level; not mixed-cell bulk |")
    lines.append("| Min cells | 30 malignant/Epi per patient | Marker gates ≠ CopyKAT |")
    lines.append("| GSEA | A8 engine, 1000 gene-set perm, seed=42 | Gene-set permutation |")
    lines.append("| FDR | BH within 15 primary sets | Not nested Broad FDR |")
    lines.append("| EMT rule | Hallmark EMT decides | GOBP EMT is recorded, not Hallmark |")
    lines.append("| IFN / MHC-I | Reported, undirected | Do not enter `supportive` |")
    lines.append("| Meta | Median NES + Stouffer of NES signs; both arms ≥5 | Sign Stouffer ignores magnitude; not a TCGA NES |")
    lines.append("")
    lines.append("Positive NES = enriched in TACSTD2-high.")
    lines.append("")
    lines.append("## Primary NES — TACSTD2 median split")
    lines.append("")
    if gsea.empty:
        lines.append("No GSEA rows.")
    else:
        med = gsea.loc[(gsea.contrast == "tacstd2_median") & (gsea.primary == True)].copy()  # noqa: E712
        cohorts = [c for c in COHORT_ORDER if c in set(med.cohort)]
        show = [t for t in HEADLINE if t in set(med.term)]
        header = "| Gene set |" + "".join(f" {c} NES | {c} FDR |" for c in cohorts)
        lines.append(header)
        lines.append("|---|" + "---:|---:|" * len(cohorts))
        for term in show:
            bits = [term]
            for c in cohorts:
                sub = med[(med.term == term) & (med.cohort == c)]
                if sub.empty:
                    bits += ["NA", "NA"]
                else:
                    bits += [fmt_nes(sub.iloc[0]["nes"]), fmt_p(sub.iloc[0]["fdr_primary"])]
            lines.append("| " + " | ".join(bits) + " |")
        lines.append("")
        lines.append("| Cohort | n_high | n_low | Verdict | Underpowered |")
        lines.append("|---|---:|---:|---|---|")
        for _, r in verd.loc[verd.contrast == "tacstd2_median"].iterrows():
            lines.append(
                f"| {r['cohort']} | {r['n_high']} | {r['n_low']} | {r['verdict']} | {r.get('underpowered', '')} |"
            )
    lines.append("")
    lines.append("## Meta NES (cohorts with both arms ≥ 5)")
    lines.append("")
    if meta.empty:
        lines.append("No meta-eligible cohort.")
    else:
        lines.append("| Gene set | k | median NES | n up / k | Stouffer z (signs) | Stouffer p | BH-FDR |")
        lines.append("|---|---:|---:|---|---:|---:|---:|")
        show = [t for t in HEADLINE if t in set(meta.term)]
        for term in show:
            r = meta.loc[meta.term == term].iloc[0]
            lines.append(
                f"| {term} | {int(r['n_cohorts'])} | {fmt_nes(r['median_NES'])} | "
                f"{int(r['n_pos'])}/{int(r['n_cohorts'])} | {fmt_nes(r['stouffer_z_signs'])} | "
                f"{fmt_p(r['stouffer_p_signs'])} | {fmt_p(r['stouffer_fdr'])} |"
            )
        lines.append("")
        lines.append(
            "Stouffer here is \(z=\\sum\\mathrm{sign}(\\mathrm{NES}_i)/\\sqrt{k}\). "
            "It tests sign concordance, not a weighted pooled NES. Do not quote it as a TCGA number."
        )
    lines.append("")
    lines.append("## What this is not")
    lines.append("")
    lines.append("- Not a re-audit of TCGA A8.")
    lines.append("- Not evidence that TACSTD2 *causes* keratin/TJ programs.")
    lines.append("- Not a TROP2 × ICI interaction test.")
    lines.append("- GSE207422 / GSE291670 malignant calls are marker proxies.")
    lines.append("- GSE241934 uses author Epi on resected tumors, not a CNV malignant call.")
    lines.append("- GSE253013 was not GSEA'd (file size / RAM / n).")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/gsea_prerank_all.tsv` — all contrasts")
    lines.append("- `tables/gsea_tacstd2_median.tsv` — locked split")
    lines.append("- `tables/meta_nes.tsv`")
    lines.append("- `tables/verdicts.tsv`")
    lines.append("- `figures/nes_heatmap_median.png`")
    lines.append("- `figures/nes_meta_and_headline.png`")
    lines.append("- `pseudobulk/*_counts.tsv.gz`")
    lines.append("")
    lines.append("Reproduce: `bash scripts/scrna_pseudobulk_gsea_meta/download.sh` then the two Python scripts. Methods: `methods/scrna_pseudobulk_gsea_meta/playbook.md`.")
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")
    (ROOT / "notes" / "scrna_pseudobulk_gsea_meta").mkdir(parents=True, exist_ok=True)
    (ROOT / "notes" / "scrna_pseudobulk_gsea_meta" / "WRITEUP.md").write_text(
        tldr + "\n\nFull report: results/scrna_pseudobulk_gsea_meta/REPORT.md\n"
    )


def _tldr(gsea, meta, inv, verd) -> str:
    if gsea.empty:
        return "No GSEA completed. See inventory."
    parts = []
    medv = verd.loc[verd.contrast == "tacstd2_median"]
    for _, r in medv.iterrows():
        parts.append(f"{r['cohort']}={r['verdict']} (n_high={r['n_high']}, n_low={r['n_low']})")
    meta_bits = []
    if not meta.empty:
        for term in [
            "KEGG_TIGHT_JUNCTION",
            "GOBP_KERATINIZATION",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
        ]:
            sub = meta.loc[meta.term == term]
            if sub.empty:
                continue
            r = sub.iloc[0]
            meta_bits.append(
                f"{term} median NES {fmt_nes(r['median_NES'])} "
                f"({int(r['n_pos'])}/{int(r['n_cohorts'])} up; Stouffer p={fmt_p(r['stouffer_p_signs'])})"
            )
    body = "Per-cohort median verdicts: " + "; ".join(parts) + "."
    if meta_bits:
        body += " Meta (arms ≥5): " + "; ".join(meta_bits) + "."
    body += (
        " This does not retract the TCGA A8 slide. Do not write “A8 conserved in scRNA” unless a gated cohort is `supportive`. "
        "n is hypothesis-generating. Hallmark EMT, not GOBP EMT, decides the EMT arm."
    )
    return body


if __name__ == "__main__":
    main()

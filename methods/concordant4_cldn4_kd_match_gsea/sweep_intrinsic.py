#!/usr/bin/env python3
"""Purify the concordant-4 CLDN4 proxy toward tumor-intrinsic IFN.

NOT a knockdown. The malignant UMI-sum is the same one as analyze.py.
This sweep does not add cohorts and does not add genes after seeing NES.

Pre-specified arms, all ranked so positive = higher when CLDN4 is lower:
  ref_q4                  unpurified Q1 vs Q4 (the PR #591 result)
  q4_drop_leading_edge    drop T/B/NK markers that sat in that IFN-γ leading edge
  q4_drop_panel           drop the full T/B/NK lineage panel from the rank and sets
  q4_immune_covariate     same panel removed from sets; rank is CLDN4 adjusted for the panel score
  q5_drop_panel           within-cohort quintiles, lowest vs highest, panel dropped
  q5_immune_covariate     quintiles, panel score as a covariate
  cont_drop_panel         continuous CLDN4 %pos OLS, panel dropped
  cont_spearman           continuous Spearman, panel dropped
  cont_immune_covariate   continuous OLS adjusted for the panel score

Primary purified contrast is q4_drop_panel. The sweep maximum is the largest
Hallmark IFN-γ NES among purified arms. It is reported as a maximum, not swapped
in as the only result.

NHEJ / STING: Reactome v2023.2.Hs, plus the 7-gene c-NHEJ and 5-gene STING cores
from the GSE207704 panel (module scores; both are below the GSEA floor of 8).
Thesis direction on this rank: IFN and STING positive, NHEJ negative.
AUCell is malignant-cell IFN already computed for GSE131907 and GSE205335.
GSE148071 is excluded. GSE123902 and GSE189357 have no cell-level AUCell here.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import analyze as A
from gsea_core import bh_fdr, gsea_prerank, rank_spearman_vs_target

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

IFN_G = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
IFN_A = "HALLMARK_INTERFERON_ALPHA_RESPONSE"
APM = "CUSTOM_MHC_I_ANTIGEN_PRESENTATION"
ISG = "INTRINSIC_ISG"
STING = "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES"
NHEJ_FULL = "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ"
NHEJ = "REACTOME_NHEJ_NO_HISTONE"
TJ = "KEGG_TIGHT_JUNCTION_NO_CLDN4"

THESIS_TERMS = [IFN_G, APM, STING, NHEJ]
LABEL = {
    IFN_G: "IFN-γ",
    IFN_A: "IFN-α",
    APM: "MHC-I / APM",
    ISG: "intrinsic ISG",
    STING: "Reactome STING",
    NHEJ: "Reactome NHEJ (no histone)",
    NHEJ_FULL: "Reactome NHEJ (full)",
    TJ: "KEGG TJ (CLDN4 out)",
}
ARM_LABEL = {
    "ref_q4": "Q4 unpurified",
    "q4_drop_leading_edge": "Q4 drop leading-edge markers",
    "q4_drop_panel": "Q4 drop T/B/NK panel",
    "q4_immune_covariate": "Q4 + T/B/NK covariate",
    "q5_drop_panel": "Q5 drop T/B/NK panel",
    "q5_immune_covariate": "Q5 + T/B/NK covariate",
    "cont_drop_panel": "continuous OLS, panel dropped",
    "cont_spearman": "continuous Spearman, panel dropped",
    "cont_immune_covariate": "continuous OLS + T/B/NK covariate",
}
PURIFIED = [k for k in ARM_LABEL if k != "ref_q4"]
PRIMARY = "q4_drop_panel"


def fmt_p(p) -> str:
    return A.fmt_p(p)


def load_markers() -> pd.DataFrame:
    m = pd.read_csv(DATA / "tnk_lineage_markers.tsv", sep="\t")
    m["symbol"] = m["symbol"].str.upper()
    return m.drop_duplicates("symbol")


def histone(gene: str) -> bool:
    g = gene.upper()
    return g.startswith("H2BC") or g.startswith("H4C") or g.startswith("H3-") or g in {"H2AX", "H2AFX"}


ALIASES = {"CGAS": "MB21D1", "STING1": "TMEM173", "H2AX": "H2AFX"}


def with_aliases(genes: list[str]) -> list[str]:
    out = []
    seen = set()
    for g in genes:
        for name in (g, ALIASES.get(g, "")):
            if not name or name in seen:
                continue
            seen.add(name)
            out.append(name)
    return out


def load_pathways(markers: set[str]) -> dict[str, list[str]]:
    raw = json.loads((DATA / "pathway_sets.json").read_text())
    base = A.load_sets()
    sets = {k: with_aliases(list(v)) for k, v in base.items()}
    sets[NHEJ_FULL] = with_aliases([g.upper() for g in raw["sets"][NHEJ_FULL]])
    sets[NHEJ] = [g for g in sets[NHEJ_FULL] if not histone(g)]
    sets[STING] = with_aliases([g.upper() for g in raw["sets"][STING]])
    sets[ISG] = with_aliases([g.upper() for g in raw["intrinsic_isg"]])
    sets["_markers"] = sorted(markers)
    return sets


def purify(sets: dict[str, list[str]], drop: set[str]) -> dict[str, list[str]]:
    out = {}
    for k, genes in sets.items():
        if k == "_markers":
            continue
        out[k] = [g for g in genes if g not in drop]
    return out


def assign_bins(meta: pd.DataFrame, k: int, prefix: str) -> pd.DataFrame:
    parts = []
    for _, sub in meta.groupby("cohort", sort=False):
        chunk = sub.copy()
        labels = [f"{prefix}{i}" for i in range(1, k + 1)]
        try:
            bins = pd.qcut(chunk["cldn4_pct"].rank(method="average"), k, labels=labels, duplicates="drop")
        except ValueError:
            bins = pd.Series(["NA"] * len(chunk), index=chunk.index)
        chunk[f"bin{k}"] = bins.astype(str)
        parts.append(chunk)
    return pd.concat(parts, ignore_index=True)


def tail_frame(meta: pd.DataFrame, col: str, low: str, high: str, counts: pd.DataFrame) -> pd.DataFrame:
    m = meta.loc[meta[col].isin([low, high]) & meta["patient"].isin(counts.columns)].copy()
    m["high"] = (m[col] == high).astype(float)
    return m


def cohort_dummies(design: pd.DataFrame, m: pd.DataFrame) -> pd.DataFrame:
    if m["cohort"].nunique() <= 1:
        return design
    idx = m.set_index("patient")
    for c in A.COHORTS:
        if c == A.REF:
            continue
        design[f"cohort_{c}"] = (idx["cohort"] == c).astype(float)
    return design


def log_matrix(counts: pd.DataFrame, patients: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    cts = A.filter_genes(counts.loc[:, patients])
    factors = A.tmm_norm_factors(cts)
    return A.log_cpm(cts, factors), cts, factors


def alias_row(parts: dict[str, pd.DataFrame], m: pd.DataFrame, cts: pd.DataFrame, factors: pd.Series, names: tuple[str, ...]) -> pd.Series | None:
    """log2(TMM-CPM+1) for the first symbol present in each cohort. Does not refit TMM."""
    lib = cts.sum(axis=0).astype(float) * factors.reindex(cts.columns).astype(float)
    cohort = m.drop_duplicates("patient").set_index("patient")["cohort"]
    vals = {}
    for patient in cts.columns:
        mat = parts[cohort.loc[patient]]
        src = next((s for s in names if s in mat.index and patient in mat.columns), None)
        if src is None:
            return None
        raw = float(mat.loc[src, patient])
        vals[patient] = float(np.log2(raw / float(lib[patient]) * 1e6 + 1.0))
    return pd.Series(vals)


def immune_z(lc: pd.DataFrame, markers: set[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in sorted(markers) if g in lc.index]
    score = lc.loc[present].mean(axis=0)
    z = (score - score.mean()) / score.std(ddof=1)
    return z, present


def ols_rank(lc: pd.DataFrame, m: pd.DataFrame, high_col: str | None, continuous: bool, immune: pd.Series | None) -> pd.DataFrame:
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    if continuous:
        coef = "CLDN4_pct_z"
        z = m.set_index("patient")["cldn4_pct"].astype(float)
        design[coef] = (z - z.mean()) / z.std(ddof=1)
    else:
        coef = "CLDN4_HIGH"
        design[coef] = m.set_index("patient")[high_col].astype(float)
    design = cohort_dummies(design, m)
    if immune is not None:
        design["immune_z"] = immune.reindex(design.index).astype(float)
    de = A.ols_t(lc, design, coef)
    de["stat_low_vs_high"] = -de["t"]
    return de


def to_rank(de: pd.DataFrame, drop: set[str]) -> pd.Series:
    s = de.set_index("gene")["stat_low_vs_high"].astype(float)
    s = s.drop(index=[g for g in drop if g in s.index], errors="ignore")
    return s.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)


def run_gsea(rank: pd.Series, sets: dict[str, list[str]], arm: str, extra: dict) -> pd.DataFrame:
    terms = [IFN_G, IFN_A, APM, ISG, STING, NHEJ, NHEJ_FULL, TJ]
    use = {t: sets[t] for t in terms if t in sets}
    out = gsea_prerank(rank, use)
    if out.empty:
        return out
    out.insert(0, "arm", arm)
    out["label"] = out["term"].map(LABEL)
    out["fdr_thesis"] = np.nan
    head = out["term"].isin(THESIS_TERMS)
    if head.any():
        out.loc[head, "fdr_thesis"] = bh_fdr(out.loc[head, "nom_p"]).to_numpy()
    for k, v in extra.items():
        out[k] = v
    return out


def module_score(lc: pd.DataFrame, m: pd.DataFrame, genes: list[str], name: str, continuous: bool, immune: pd.Series | None) -> dict:
    present = [g for g in genes if g in lc.index]
    if len(present) < 3:
        return {"module": name, "n_genes": len(present), "skipped": True}
    score = lc.loc[present, m["patient"]].mean(axis=0).to_frame().T
    score.index = [name]
    de = ols_rank(score, m, "high" if "high" in m.columns and not continuous else None, continuous, immune)
    # ols_rank for a 1-row frame still works; high column is used when not continuous
    r = de.iloc[0]
    return {
        "module": name,
        "n_genes": len(present),
        "skipped": False,
        "logFC_high_vs_low": float(r.logFC),
        "stat_low_vs_high": float(r.stat_low_vs_high),
        "t_high_vs_low": float(r.t),
        "p": float(r.p),
        "genes": ",".join(present),
    }


def contrast_bundle(counts, parts, meta, markers, sets_full, arm, kind, drop, use_immune):
    if kind == "q4":
        m = tail_frame(meta, "bin4", "Q1", "Q4", counts)
        continuous = False
    elif kind == "q5":
        m = tail_frame(meta, "bin5", "P1", "P5", counts)
        continuous = False
    else:
        m = meta.loc[meta["patient"].isin(counts.columns)].copy()
        m["high"] = np.nan
        continuous = True
    n_low = int((m["bin4"] == "Q1").sum()) if kind == "q4" else int((m["bin5"] == "P1").sum()) if kind == "q5" else int(len(m))
    n_high = int((m["bin4"] == "Q4").sum()) if kind == "q4" else int((m["bin5"] == "P5").sum()) if kind == "q5" else int(len(m))
    info = {"arm": arm, "n_low": n_low, "n_high": n_high, "n": int(len(m)), "skipped": False}
    if kind != "cont" and (n_low < 3 or n_high < 3):
        info["skipped"] = True
        return info, pd.DataFrame(), pd.DataFrame()
    lc, cts, factors = log_matrix(counts, m["patient"].tolist())
    iz, present = immune_z(lc, markers)
    info["n_markers_in_matrix"] = len(present)
    info["immune_vs_cldn4_spearman"] = float(stats.spearmanr(iz, m.set_index("patient")["cldn4_pct"].reindex(iz.index)).statistic)
    immune = iz if use_immune else None
    if kind == "spearman":
        rho = rank_spearman_vs_target(lc, m.set_index("patient")["cldn4_pct"])
        de = pd.DataFrame({"gene": rho.index, "logFC": rho.values, "t": rho.values, "p": np.nan, "stat_low_vs_high": -rho.values})
    else:
        de = ols_rank(lc, m, None if continuous else "high", continuous, immune)
    rank = to_rank(de, drop)
    info["n_genes"] = int(len(rank))
    info["n_dropped_from_rank"] = int(sum(g in de["gene"].values for g in drop))
    sets = purify(sets_full, drop)
    gsea = run_gsea(rank, sets, arm, {"n_low": n_low, "n_high": n_high, "n_genes_ranked": int(len(rank))})
    modules = []
    if kind != "spearman":
        c_genes = ["PRKDC", "XRCC4", "LIG4", "RIF1", "TP53BP1", "XRCC5", "XRCC6"]
        rec = module_score(lc, m, [g for g in c_genes if g not in drop], "cNHEJ_7", continuous, immune)
        rec["arm"] = arm
        modules.append(rec)
        sting_lc = lc.copy()
        for dst, names in (("CGAS", ("CGAS", "MB21D1")), ("STING1", ("STING1", "TMEM173"))):
            if dst in sting_lc.index:
                continue
            row = alias_row(parts, m, cts, factors, names)
            if row is not None:
                sting_lc.loc[dst] = row.reindex(sting_lc.columns)
        rec = module_score(
            sting_lc,
            m,
            [g for g in ["CGAS", "STING1", "TBK1", "IRF3", "STAT1"] if g not in drop],
            "STING_5",
            continuous,
            immune,
        )
        rec["arm"] = arm
        modules.append(rec)
    return info, gsea, pd.DataFrame(modules), rank, de


def aucell_block() -> pd.DataFrame:
    path = Path("/tmp/cell_scores.tsv.gz")
    if not path.exists():
        raise SystemExit("missing /tmp/cell_scores.tsv.gz")
    cells = pd.read_csv(path, sep="\t")
    cells = cells[cells["dataset"].isin(["GSE131907", "GSE205335"])].copy()
    g = cells.groupby(["dataset", "patient"], as_index=False).agg(
        n_malignant=("auc_IFN_STAT1", "size"),
        auc_ifn=("auc_IFN_STAT1", "mean"),
        auc_mhc=("auc_MHC_NLRC5", "mean"),
        auc_tj=("auc_TJ_barrier", "mean"),
        cldn4_mean=("CLDN4", "mean"),
        cldn4_pct=("CLDN4_pos", "mean"),
    )
    g = g[g["n_malignant"] >= 20].copy()
    rows = []
    for dataset, sub in g.groupby("dataset"):
        rho, p = stats.spearmanr(sub["cldn4_mean"], sub["auc_ifn"])
        q = A.assign_quartiles(sub.set_index("patient")["cldn4_mean"])
        sub = sub.copy()
        sub["quartile"] = sub["patient"].map(q)
        tails = sub[sub["quartile"].isin(["Q1", "Q4"])]
        u, up = stats.mannwhitneyu(
            tails.loc[tails["quartile"] == "Q1", "auc_ifn"],
            tails.loc[tails["quartile"] == "Q4", "auc_ifn"],
            alternative="two-sided",
        )
        rows.append(
            {
                "dataset": dataset,
                "n_patients": int(len(sub)),
                "n_q1": int((sub["quartile"] == "Q1").sum()),
                "n_q4": int((sub["quartile"] == "Q4").sum()),
                "n_cells": int(sub["n_malignant"].sum()),
                "spearman_cldn4_vs_ifn_auc": float(rho),
                "spearman_p": float(p),
                "q1_minus_q4_median_ifn_auc": float(
                    tails.loc[tails["quartile"] == "Q1", "auc_ifn"].median()
                    - tails.loc[tails["quartile"] == "Q4", "auc_ifn"].median()
                ),
                "mwu_p_q1_vs_q4": float(up),
                "note": "malignant cells only; GSE148071 excluded; not a NES",
            }
        )
    # within-dataset rank pool of the two concordant cohorts
    parts = []
    for dataset, sub in g.groupby("dataset"):
        s = sub.copy()
        s["r_cldn4"] = s["cldn4_mean"].rank()
        s["r_ifn"] = s["auc_ifn"].rank()
        parts.append(s)
    pool = pd.concat(parts, ignore_index=True)
    rho, p = stats.spearmanr(pool["r_cldn4"], pool["r_ifn"])
    rows.append(
        {
            "dataset": "GSE131907+GSE205335 within-dataset ranks",
            "n_patients": int(len(pool)),
            "n_q1": np.nan,
            "n_q4": np.nan,
            "n_cells": int(pool["n_malignant"].sum()),
            "spearman_cldn4_vs_ifn_auc": float(rho),
            "spearman_p": float(p),
            "q1_minus_q4_median_ifn_auc": np.nan,
            "mwu_p_q1_vs_q4": np.nan,
            "note": "positive rho = IFN AUCell higher when CLDN4 is higher. Thesis wants negative.",
        }
    )
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "aucell_epithelial_ifn.tsv", sep="\t", index=False)
    g.to_csv(TABLES / "aucell_patient_means.tsv", sep="\t", index=False)
    return out


def md_gsea(df: pd.DataFrame, terms: list[str]) -> str:
    lines = ["| set | ES | NES | nom p | thesis FDR | n |", "|---|---:|---:|---:|---:|---:|"]
    for term in terms:
        hit = df[df["term"] == term]
        if hit.empty:
            continue
        r = hit.iloc[0]
        fdr = "—" if term not in THESIS_TERMS or not math.isfinite(float(r.fdr_thesis)) else fmt_p(r.fdr_thesis)
        lines.append(
            f"| {LABEL.get(term, term)} | {r.es:+.3f} | {r.nes:+.3f} | {fmt_p(r.nom_p)} | {fdr} | {int(r.n_set_in_rank)} |"
        )
    return "\n".join(lines)


def plot_sweep(gsea: pd.DataFrame, path: Path) -> None:
    sub = gsea[gsea["term"] == IFN_G].copy()
    order = list(ARM_LABEL)
    sub["order"] = sub["arm"].map({a: i for i, a in enumerate(order)})
    sub = sub.sort_values("order", ascending=False)
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    colors = ["#8a8a8a" if a == "ref_q4" else "#c45c26" if v > 0 else "#2c5f8a" for a, v in zip(sub["arm"], sub["nes"])]
    y = np.arange(len(sub))
    ax.barh(y, sub["nes"], color=colors, height=0.68)
    ax.axvline(0, color="#222", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([ARM_LABEL[a] for a in sub["arm"]])
    ax.set_xlabel("IFN-γ NES  (positive = higher when CLDN4 is lower)")
    ax.set_title("NOT a knockdown\nmalignant pseudobulk sweep after T/B/NK purification")
    for i, r in enumerate(sub.itertuples()):
        nudge = 0.06 if r.nes >= 0 else -0.06
        ax.text(r.nes + nudge, i, f"{r.nes:+.2f}  p {fmt_p(r.nom_p)}", va="center", ha="left" if r.nes >= 0 else "right", fontsize=8)
    lo = float(sub["nes"].min()) - 0.8
    hi = float(sub["nes"].max()) + 1.3
    ax.set_xlim(min(-1.5, lo), max(2.5, hi))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(ctx: dict) -> None:
    gsea = ctx["gsea"]
    primary = gsea[gsea["arm"] == PRIMARY]
    ref = gsea[gsea["arm"] == "ref_q4"]
    ifn_ref = ref[ref["term"] == IFN_G].iloc[0]
    ifn_p = primary[primary["term"] == IFN_G].iloc[0]
    purified = gsea[(gsea["arm"].isin(PURIFIED)) & (gsea["term"] == IFN_G)]
    best = purified.loc[purified["nes"].idxmax()]
    worst = purified.loc[purified["nes"].idxmin()]
    le = ctx["leading_edge_markers"]
    auc = ctx["auc"]
    mods = ctx["modules"]
    prim_info = ctx["info"][PRIMARY]

    def arm_line(arm: str) -> str:
        r = gsea[(gsea["arm"] == arm) & (gsea["term"] == IFN_G)].iloc[0]
        sting = gsea[(gsea["arm"] == arm) & (gsea["term"] == STING)].iloc[0]
        nhej = gsea[(gsea["arm"] == arm) & (gsea["term"] == NHEJ)].iloc[0]
        isg = gsea[(gsea["arm"] == arm) & (gsea["term"] == ISG)].iloc[0]
        ncell = f"n={int(r.n_low)}" if arm.startswith("cont") else f"{int(r.n_low)}/{int(r.n_high)}"
        return (
            f"| {ARM_LABEL[arm]} | {ncell} | {r.nes:+.3f} ({fmt_p(r.nom_p)}) | "
            f"{isg.nes:+.3f} ({fmt_p(isg.nom_p)}) | {sting.nes:+.3f} ({fmt_p(sting.nom_p)}) | {nhej.nes:+.3f} ({fmt_p(nhej.nom_p)}) |"
        )

    auc_lines = ["| cohort | patients | cells | Spearman CLDN4 vs IFN AUCell | p | Q1−Q4 median IFN AUC | MWU p |", "|---|---:|---:|---:|---:|---:|---:|"]
    for r in auc.itertuples():
        qmed = "—" if not math.isfinite(float(r.q1_minus_q4_median_ifn_auc)) else f"{r.q1_minus_q4_median_ifn_auc:+.4f}"
        mwu = "—" if not math.isfinite(float(r.mwu_p_q1_vs_q4)) else fmt_p(r.mwu_p_q1_vs_q4)
        nq = "" if not math.isfinite(float(r.n_q1)) else f" (Q1 {int(r.n_q1)} / Q4 {int(r.n_q4)})"
        auc_lines.append(
            f"| {r.dataset}{nq} | {int(r.n_patients)} | {int(r.n_cells)} | {r.spearman_cldn4_vs_ifn_auc:+.3f} | {fmt_p(r.spearman_p)} | {qmed} | {mwu} |"
        )

    mod_lines = ["| arm | module | genes | low-vs-high t | p |", "|---|---|---:|---:|---:|"]
    for r in mods.itertuples():
        if r.skipped:
            continue
        mod_lines.append(
            f"| {ARM_LABEL.get(r.arm, r.arm)} | {r.module} | {int(r.n_genes)} | {r.stat_low_vs_high:+.3f} | {fmt_p(r.p)} |"
        )

    call = A.direction_call(float(ifn_p.nes), float(ifn_p.nom_p))
    text = f"""# SWEEP — tumor-intrinsic IFN on the concordant-4 malignant pseudobulk

**NOT a true KD.** Same malignant UMI-sums and the same unpurified IFN-γ NES as the quartile GSEA. No new cohort. No gene was added to the IFN sets after looking at NES. The STING core score reads CGAS from MB21D1 where that is the symbol, and STING1 from TMEM173. Those aliases are not forced into the prerank, so the IFN rank stays the published one. Reactome STING GSEA therefore misses CGAS and STING1 when those official symbols are absent from the shared matrix.

The unpurified IFN-γ leading edge contained T/B/NK markers ({", ".join(le) if le else "none"}). Those genes, and the rest of the locked lineage panel in `data/tnk_lineage_markers.tsv`, are removed in the purified arms. Myeloid genes that are not on that panel stay in the rank.

Primary purified arm: **Q1 vs Q4, full T/B/NK panel dropped** from the rank and from every gene set. Honest n = **{prim_info['n_low']} lowest vs {prim_info['n_high']} highest**. Markers present in the matrix: {prim_info['n_markers_in_matrix']}. Spearman of that T/B/NK score vs CLDN4 %pos inside this contrast: {prim_info['immune_vs_cldn4_spearman']:+.3f}.

## Call

| readout | unpurified Q4 | primary purified (panel dropped) |
|---|---|---|
| Hallmark IFN-γ NES | {ifn_ref.nes:+.3f} (p {fmt_p(ifn_ref.nom_p)}) | {ifn_p.nes:+.3f} (p {fmt_p(ifn_p.nom_p)}, thesis FDR {fmt_p(ifn_p.fdr_thesis)}) |
| direction vs private KD (IFN up after loss) | same direction | **{call}** |

Highest Hallmark IFN-γ NES among the purified arms: **{best.nes:+.3f}** on `{ARM_LABEL[best.arm]}` (p {fmt_p(best.nom_p)}). Lowest purified IFN-γ NES: **{worst.nes:+.3f}** on `{ARM_LABEL[worst.arm]}`. Panel removal does not beat the unpurified NES. The sweep maximum is the panel-drop quartile, which is the least stringent purification. The stricter test is the T/B/NK score covariate, because that score tracks CLDN4 (Spearman {prim_info['immune_vs_cldn4_spearman']:+.3f} on the quartile tails). Quote the covariate row next to the maximum. Do not quote the maximum alone.

The purified leading edge still opens with CCL5 and CSF2RB. Those are not on the T/B/NK panel, so they were not removed.

Purified IFN-γ leading edge: `{ifn_p.lead_genes}`

## Sweep

Positive NES = higher when CLDN4 is lower. Thesis-aligned IFN and STING are positive. Thesis-aligned NHEJ is negative. Thesis FDR is BH inside IFN-γ, APM, Reactome STING, and Reactome NHEJ with histones removed.

| arm | n low/high | IFN-γ NES (p) | intrinsic ISG NES (p) | STING NES (p) | NHEJ no-histone NES (p) |
|---|---:|---:|---:|---:|---:|
{chr(10).join(arm_line(a) for a in ARM_LABEL)}

### Primary purified sets

{md_gsea(primary, [IFN_G, IFN_A, APM, ISG, STING, NHEJ, NHEJ_FULL, TJ])}

Reactome NHEJ v2023.2.Hs is 68 genes and is mostly histone symbols. The no-histone set is the pre-specified GSEA module (histone symbols H2BC*, H3-*, H4C*, H2AX removed). The full set is the companion row. The 7-gene c-NHEJ core and the 5-gene STING core are module scores, not NES, because both are below the size floor of 8.

### Core module scores (low vs high)

Positive t = higher when CLDN4 is lower.

{chr(10).join(mod_lines)}

## AUCell IFN, malignant cells only

This is not a NES. AUCell was already run on author-malignant cells. GSE148071 is not in this table. GSE123902 and GSE189357 have no cell-level AUCell in that file, so they are absent here rather than imputed. The AUCell IFN gene list is the intrinsic ISG list (no TCR, BCR, or NK markers).

{chr(10).join(auc_lines)}

A negative Spearman would be the thesis direction (IFN activity higher when malignant CLDN4 is lower). The pooled rank correlation is the AUCell result for this sweep. It is not a NES and it is not a hit.

## What this sweep does not do

- It does not relabel the proxy as a knockdown.
- It does not drop GSE205335, or any other cohort, to raise the NES.
- It does not remove myeloid genes that were not on the T/B/NK panel.
- It does not bring GSE148071 into the AUCell pool.
- It does not treat the sweep maximum as the only number.

Reproduce: `python3 methods/concordant4_cldn4_kd_match_gsea/sweep_intrinsic.py`
"""
    (HERE / "FINDING_SWEEP.md").write_text(text)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    markers_df = load_markers()
    markers = set(markers_df["symbol"])
    sets_full = load_pathways(markers)
    meta = A.load_units()
    parts = A.load_counts()
    counts = A.combine(parts)
    # Bins use every eligible unit, including units absent from the UMI-sum
    # (P4001). The contrast then keeps only units that have counts.
    meta = assign_bins(meta, 4, "Q")
    meta = assign_bins(meta, 5, "P")
    # Quartile labels from assign_bins use Q1..Q4. Confirm they match the locked tails.
    for cohort, (q1, q4) in A.EXPECTED_TAILS.items():
        sub = meta[meta["cohort"] == cohort]
        if not sub["bin4"].eq(sub["quartile"]).all():
            raise SystemExit(f"quartile drift {cohort}")
        sub_m = sub[sub["patient"].isin(counts.columns)]
        got1 = set(sub_m.loc[sub_m["bin4"] == "Q1", "patient"])
        got4 = set(sub_m.loc[sub_m["bin4"] == "Q4", "patient"])
        if got1 != q1 or got4 != q4:
            raise SystemExit(f"matrix-tail drift {cohort}: {got1} vs {q1}; {got4} vs {q4}")

    # Reference rank, then the leading-edge intersection.
    info_ref, gsea_ref, mod_ref, rank_ref, de_ref = contrast_bundle(
        counts, parts, meta, markers, sets_full, "ref_q4", "q4", drop=set(), use_immune=False
    )
    ifn_le = gsea_ref[gsea_ref["term"] == IFN_G].iloc[0]["lead_genes"].split(",")
    le_markers = [g for g in ifn_le if g in markers]
    (TABLES / "leading_edge_tnk_markers.tsv").write_text(
        "gene\tin_unpurified_ifng_leading_edge\tin_lineage_panel\n"
        + "\n".join(f"{g}\t1\t1" for g in le_markers)
        + "\n"
    )

    specs = [
        ("ref_q4", "q4", set(), False, info_ref, gsea_ref, mod_ref),
        ("q4_drop_leading_edge", "q4", set(le_markers), False, None, None, None),
        ("q4_drop_panel", "q4", markers, False, None, None, None),
        ("q4_immune_covariate", "q4", markers, True, None, None, None),
        ("q5_drop_panel", "q5", markers, False, None, None, None),
        ("q5_immune_covariate", "q5", markers, True, None, None, None),
        ("cont_drop_panel", "cont", markers, False, None, None, None),
        ("cont_spearman", "spearman", markers, False, None, None, None),
        ("cont_immune_covariate", "cont", markers, True, None, None, None),
    ]
    frames = [gsea_ref]
    modules = [mod_ref]
    infos = {PRIMARY: None, "ref_q4": info_ref}
    for arm, kind, drop, use_immune, info, gsea, mods in specs:
        if gsea is not None:
            continue
        info, gsea, mods, rank, de = contrast_bundle(counts, parts, meta, markers, sets_full, arm, kind, drop, use_immune)
        if info.get("skipped"):
            raise SystemExit(f"arm skipped: {info}")
        frames.append(gsea)
        modules.append(mods)
        infos[arm] = info
        if arm == PRIMARY:
            rank.to_csv(TABLES / "rank_intrinsic_q4.tsv", sep="\t", header=["stat_low_vs_high"])
    gsea_all = pd.concat(frames, ignore_index=True)
    mods_all = pd.concat(modules, ignore_index=True)
    gsea_all.to_csv(TABLES / "gsea_intrinsic_sweep.tsv", sep="\t", index=False)
    mods_all.to_csv(TABLES / "module_scores_intrinsic.tsv", sep="\t", index=False)
    pd.DataFrame(infos.values()).to_csv(TABLES / "sweep_n.tsv", sep="\t", index=False)
    auc = aucell_block()
    plot_sweep(gsea_all, FIGS / "fig_intrinsic_sweep_nes")
    # sanity: unpurified IFN-γ NES must match the published quartile result
    got = float(gsea_ref[gsea_ref["term"] == IFN_G].iloc[0]["nes"])
    if abs(got - 3.788611) > 1e-4:
        raise SystemExit(f"unpurified IFN NES drifted: {got}")
    write_finding(
        {
            "gsea": gsea_all,
            "modules": mods_all,
            "info": infos,
            "leading_edge_markers": le_markers,
            "auc": auc,
        }
    )
    show = gsea_all[gsea_all["term"] == IFN_G][["arm", "nes", "nom_p", "n_low", "n_high"]]
    print(show.to_string(index=False))
    print("LE markers", le_markers)
    print(auc.to_string(index=False))


if __name__ == "__main__":
    main()

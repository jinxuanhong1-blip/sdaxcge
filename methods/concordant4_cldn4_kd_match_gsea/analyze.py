#!/usr/bin/env python3
"""Observational proxy: concordant-4 malignant CLDN4 Q1 vs Q4 prerank GSEA.

NOT a knockdown. Lowest vs highest malignant CLDN4 %pos quartile, then
IFN / APM / tight-junction GSEA. The rank is oriented so a positive score
means higher in the lowest CLDN4 quartile (the observational analog of
"after loss"). The private KD comparator is a direction only: IFN up after
CLDN4 loss. That matrix is not in this repository and is not reanalyzed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gsea_core import bh_fdr, enrichment_walk, gsea_prerank

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}

# Locked tails from the concordant-4 malignant stack (units in the UMI-sum).
EXPECTED_TAILS = {
    "GSE123902": (
        {"LX675", "LX682", "LX699", "LX701"},
        {"LX653", "LX680", "LX684"},
    ),
    "GSE131907": (
        {"EBUS_13", "EBUS_15", "EBUS_49", "NS_02", "NS_06", "NS_16"},
        {"EBUS_19", "EBUS_28", "NS_03", "NS_04", "NS_07"},
    ),
    "GSE205335": (
        {"P1015", "P1062", "P1063", "P1090", "P1119"},
        {"P1016", "P1025", "P1037", "P1084", "P1089", "P1115"},
    ),
    "GSE189357": ({"TD2", "TD4", "TD7"}, {"TD6", "TD9"}),
}

HEADLINE = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    "KEGG_TIGHT_JUNCTION_NO_CLDN4",
]
LABEL = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "MHC-I / APM",
    "KEGG_TIGHT_JUNCTION_NO_CLDN4": "KEGG TJ (CLDN4 out)",
    "KEGG_TIGHT_JUNCTION": "KEGG TJ (CLDN4 in)",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION_NO_CLDN4": "GO TJ organization (CLDN4 out)",
}
PRIVATE_IFN = "up after CLDN4 loss"


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.astype(float).rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=values.index)


def load_units() -> pd.DataFrame:
    rows = []

    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"]
    for r in el.itertuples(index=False):
        rows.append(("GSE123902", "donor", str(r.patient), float(r.mal_CLDN4_pct), float(r.mal_CLDN4_mean)))

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)]
    for r in mal.itertuples(index=False):
        rows.append(("GSE131907", "sample", str(r.sample), float(r.mal_CLDN4_pct), float(r.mal_CLDN4_mean)))

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    for r in d.itertuples(index=False):
        rows.append(("GSE205335", "patient", str(r.patient), float(r.mal_CLDN4_pct_pos), float(r.mal_CLDN4_mean)))

    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"]
    for r in el.itertuples(index=False):
        rows.append(("GSE189357", "patient", str(r.patient), float(r.mal_CLDN4_pct), float(r.mal_CLDN4_mean)))

    meta = pd.DataFrame(rows, columns=["cohort", "unit", "patient", "cldn4_pct", "cldn4_mean"])
    scaled = []
    for cohort, sub in meta.groupby("cohort", sort=False):
        pct = sub["cldn4_pct"].to_numpy(float)
        if np.nanmax(pct) <= 1.5:
            pct = pct * 100.0
        q = assign_quartiles(pd.Series(pct, index=sub.index))
        chunk = sub.copy()
        chunk["cldn4_pct"] = pct
        chunk["quartile"] = q.to_numpy()
        scaled.append(chunk)
    return pd.concat(scaled, ignore_index=True)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


def load_counts() -> dict[str, pd.DataFrame]:
    return {c: read_counts(DATA / f"{c}_malignant_counts.tsv.gz") for c in COHORTS}


def combine(parts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    genes = sorted(set.intersection(*[set(p.index) for p in parts.values()]))
    combined = pd.concat([p.reindex(genes).fillna(0) for p in parts.values()], axis=1)
    return combined.loc[:, ~combined.columns.duplicated()]


def tmm_norm_factors(counts: pd.DataFrame) -> pd.Series:
    lib = counts.sum(axis=0).astype(float).replace(0, np.nan)
    rel = counts.div(lib, axis=1)
    f75 = rel.quantile(0.75, axis=0)
    ref = (f75 - f75.mean()).abs().idxmin()
    ref_c = counts[ref].astype(float)
    ref_lib = float(lib[ref])
    factors = {}
    for col in counts.columns:
        obs = counts[col].astype(float)
        obs_lib = float(lib[col])
        keep = (obs > 0) & (ref_c > 0)
        if int(keep.sum()) < 50:
            factors[col] = 1.0
            continue
        m = np.log2((obs[keep] / obs_lib) / (ref_c[keep] / ref_lib))
        a = 0.5 * np.log2((obs[keep] / obs_lib) * (ref_c[keep] / ref_lib))
        w = (obs_lib - obs[keep]) / (obs_lib * obs[keep]) + (ref_lib - ref_c[keep]) / (ref_lib * ref_c[keep])
        ok = np.isfinite(m) & np.isfinite(a) & np.isfinite(w) & (w > 0)
        m, a, w = m[ok], a[ok], w[ok]
        if len(m) < 50:
            factors[col] = 1.0
            continue
        lo_m, hi_m = np.quantile(m, [0.30, 0.70])
        lo_a, hi_a = np.quantile(a, [0.05, 0.95])
        trim = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if int(trim.sum()) < 20:
            factors[col] = 1.0
            continue
        factors[col] = 2 ** float(np.average(m[trim], weights=1.0 / w[trim]))
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def ols_t(logcpm: pd.DataFrame, design: pd.DataFrame, coef: str) -> pd.DataFrame:
    X = design.reindex(logcpm.columns).astype(float)
    if X.isna().any().any():
        keep = ~X.isna().any(axis=1)
        X = X.loc[keep]
        Y = logcpm.loc[:, keep].astype(float).to_numpy()
    else:
        Y = logcpm.astype(float).to_numpy()
    Xv = X.to_numpy()
    n, p = Xv.shape
    df_res = n - p
    xtx_inv = np.linalg.pinv(Xv.T @ Xv)
    beta = xtx_inv @ Xv.T @ Y.T
    resid = Y.T - Xv @ beta
    sigma2 = np.sum(resid**2, axis=0) / df_res
    j = list(X.columns).index(coef)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[j, j], 0.0))
    est = beta[j]
    t = np.divide(est, se, out=np.zeros_like(est), where=se > 0)
    pv = 2.0 * stats.t.sf(np.abs(t), df_res)
    return pd.DataFrame({"gene": logcpm.index.astype(str), "logFC": est, "t": t, "se": se, "p": pv, "df": df_res})


def add_cohort_dummies(design: pd.DataFrame, m: pd.DataFrame) -> pd.DataFrame:
    idx = m.set_index("patient")
    for c in COHORTS:
        if c == REF:
            continue
        design[f"cohort_{c}"] = (idx["cohort"] == c).astype(float)
    return design


def design_q4(m: pd.DataFrame, cohort: str | None) -> pd.DataFrame:
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_Q4"] = (m.set_index("patient")["quartile"] == "Q4").astype(float)
    if cohort is None and m["cohort"].nunique() > 1:
        design = add_cohort_dummies(design, m)
    return design


def prepare_q4(counts: pd.DataFrame, meta: pd.DataFrame, cohort: str | None):
    m = meta.copy()
    if cohort:
        m = m.loc[m["cohort"] == cohort]
    m = m.loc[m["quartile"].isin(["Q1", "Q4"])]
    m = m.loc[m["patient"].isin(counts.columns)].copy()
    n1 = int((m["quartile"] == "Q1").sum())
    n4 = int((m["quartile"] == "Q4").sum())
    info = {
        "cohort": cohort or "+".join(COHORTS),
        "n_q1": n1,
        "n_q4": n4,
        "n": n1 + n4,
        "patients_q1": ",".join(sorted(m.loc[m["quartile"] == "Q1", "patient"])),
        "patients_q4": ",".join(sorted(m.loc[m["quartile"] == "Q4", "patient"])),
        "skipped": n1 < 3 or n4 < 3,
    }
    if info["skipped"]:
        return pd.DataFrame(), info
    cts = filter_genes(counts.loc[:, m["patient"]])
    de = ols_t(log_cpm(cts, tmm_norm_factors(cts)), design_q4(m, cohort), "CLDN4_Q4")
    info["n_genes"] = int(len(de))
    return de, info


def prepare_continuous(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    m = meta.loc[meta["patient"].isin(counts.columns)].copy()
    cts = filter_genes(counts.loc[:, m["patient"]])
    lc = log_cpm(cts, tmm_norm_factors(cts))
    z = m.set_index("patient")["cldn4_pct"].astype(float)
    z = (z - z.mean()) / z.std(ddof=1)
    design = pd.DataFrame(index=m["patient"])
    design["Intercept"] = 1.0
    design["CLDN4_pct_z"] = z
    design = add_cohort_dummies(design, m)
    de = ols_t(lc, design, "CLDN4_pct_z")
    return de, {"n": int(len(m)), "n_genes": int(len(de))}


def low_vs_high_rank(de: pd.DataFrame, drop_genes: set[str] | None = None) -> pd.Series:
    """Positive = higher in the lowest CLDN4 quartile (Q1 minus Q4)."""
    s = de.set_index("gene")["t"].astype(float) * -1.0
    if drop_genes:
        s = s.drop(index=[g for g in drop_genes if g in s.index], errors="ignore")
    return s.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)


def load_sets() -> dict[str, list[str]]:
    raw = json.loads((DATA / "gene_sets.json").read_text())["sets"]
    sets = {k: [g.upper() for g in v] for k, v in raw.items()}
    sets["KEGG_TIGHT_JUNCTION_NO_CLDN4"] = [g for g in sets["KEGG_TIGHT_JUNCTION"] if g != "CLDN4"]
    sets["GOBP_TIGHT_JUNCTION_ORGANIZATION_NO_CLDN4"] = [
        g for g in sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"] if g != "CLDN4"
    ]
    return sets


def run_sets(rank: pd.Series, sets: dict[str, list[str]], terms: list[str], contrast: str) -> pd.DataFrame:
    use = {t: sets[t] for t in terms}
    out = gsea_prerank(rank, use)
    if out.empty:
        return out
    out.insert(0, "contrast", contrast)
    out["label"] = out["term"].map(LABEL)
    out["fdr_headline"] = np.nan
    head = out["term"].isin(HEADLINE)
    out.loc[head, "fdr_headline"] = bh_fdr(out.loc[head, "nom_p"]).to_numpy()
    return out


def direction_call(nes: float, p: float) -> str:
    if not math.isfinite(nes) or not math.isfinite(p):
        return "not resolved"
    if p >= 0.05:
        return "not resolved"
    if nes > 0:
        return "same direction"
    if nes < 0:
        return "opposite direction"
    return "not resolved"


def fmt_p(p) -> str:
    p = float(p)
    if not math.isfinite(p):
        return "—"
    # 1000 permutations cannot resolve below 1/1001.
    if p <= (1.0 / 1001.0) + 1e-12:
        return "0.001"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt(x, nd=3) -> str:
    x = float(x)
    if not math.isfinite(x):
        return "—"
    return f"{x:.{nd}f}"


def plot_headline(df: pd.DataFrame, path: Path) -> None:
    sub = df[df["term"].isin(HEADLINE)].copy()
    sub["order"] = sub["term"].map({t: i for i, t in enumerate(HEADLINE)})
    sub = sub.sort_values("order", ascending=False)
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    colors = ["#c45c26" if v > 0 else "#2c5f8a" for v in sub["nes"]]
    y = np.arange(len(sub))
    ax.barh(y, sub["nes"], color=colors, height=0.62, zorder=2)
    ax.axvline(0, color="#222", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"])
    ax.set_xlabel("NES  (positive = higher in lowest CLDN4 quartile)")
    ax.set_title("NOT a knockdown\nconcordant-4 malignant Q1 vs Q4")
    for i, r in enumerate(sub.itertuples()):
        x = r.nes + (0.04 if r.nes >= 0 else -0.04)
        ha = "left" if r.nes >= 0 else "right"
        ax.text(x, i, f"{r.nes:+.2f}   FDR {fmt_p(r.fdr_headline)}", va="center", ha=ha, fontsize=8)
    ax.set_xlim(min(-2.4, float(sub["nes"].min()) - 0.8), max(2.4, float(sub["nes"].max()) + 0.9))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_curves(rank: pd.Series, sets: dict[str, list[str]], path: Path) -> None:
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(float)
    abs_s = np.abs(scores)
    pos = {g: i for i, g in enumerate(genes)}
    terms = [
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
        "KEGG_TIGHT_JUNCTION_NO_CLDN4",
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2), sharey=True)
    for ax, term in zip(axes, terms):
        idx = np.unique(np.array([pos[g] for g in sets[term] if g in pos], dtype=int))
        hit = np.zeros(len(genes), dtype=bool)
        hit[idx] = True
        es, walk, _ = enrichment_walk(abs_s, hit)
        ax.plot(np.arange(len(walk)), walk, color="#c45c26" if es >= 0 else "#2c5f8a", lw=1.4)
        ax.axhline(0, color="#222", lw=0.6)
        ax.set_title(f"{LABEL[term]}\nES {es:+.3f}", fontsize=9)
        ax.set_xlabel("rank (low CLDN4 →)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("enrichment score")
    fig.suptitle("Running enrichment, lowest vs highest CLDN4 quartile", fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_cohorts(df: pd.DataFrame, path: Path) -> None:
    terms = [
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
        "KEGG_TIGHT_JUNCTION_NO_CLDN4",
    ]
    cohorts = [c for c in COHORTS if c != "GSE189357"] + ["stacked"]
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    width = 0.24
    x = np.arange(len(cohorts))
    colors = ["#c45c26", "#e0a106", "#2c5f8a"]
    for i, term in enumerate(terms):
        vals = []
        for c in cohorts:
            key = "stacked_q1_vs_q4" if c == "stacked" else f"q1_vs_q4_{c}"
            hit = df[(df["contrast"] == key) & (df["term"] == term)]
            vals.append(float(hit.iloc[0]["nes"]) if len(hit) else np.nan)
        ax.bar(x + (i - 1) * width, vals, width=width, color=colors[i], label=LABEL[term])
    ax.axhline(0, color="#222", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(["GSE123902\n4 vs 3", "GSE131907\n6 vs 5", "GSE205335\n5 vs 6", "stacked\n18 vs 16"])
    ax.set_ylabel("NES (positive = higher in Q1)")
    ax.set_title("Per-cohort NES. GSE189357 Q4 n=2 is in the stack only.")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def md_table(df: pd.DataFrame, terms: list[str]) -> str:
    lines = [
        "| set | ES | NES | nom p | FDR | n in rank |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for term in terms:
        hit = df[df["term"] == term]
        if hit.empty:
            continue
        r = hit.iloc[0]
        fdr = "—" if term not in HEADLINE or not math.isfinite(float(r.fdr_headline)) else fmt_p(r.fdr_headline)
        lines.append(
            f"| {LABEL.get(term, term)} | {r.es:+.3f} | {r.nes:+.3f} | {fmt_p(r.nom_p)} | {fdr} | {int(r.n_set_in_rank)} |"
        )
    return "\n".join(lines)


def write_finding(ctx: dict) -> None:
    h = ctx["headline"]
    ifn = h[h["term"] == "HALLMARK_INTERFERON_GAMMA_RESPONSE"].iloc[0]
    apm = h[h["term"] == "CUSTOM_MHC_I_ANTIGEN_PRESENTATION"].iloc[0]
    tj = h[h["term"] == "KEGG_TIGHT_JUNCTION_NO_CLDN4"].iloc[0]
    ifna = h[h["term"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE"].iloc[0]
    call = direction_call(float(ifn.nes), float(ifn.nom_p))
    cldn4 = ctx["cldn4"]
    info = ctx["info"]
    lead_ifn = ifn.lead_genes
    lead_apm = apm.lead_genes
    lead_tj = tj.lead_genes
    top_genes = ", ".join(ctx["top_genes"])

    def sens_block(name: str) -> str:
        sub = ctx["all"]
        sub = sub[sub["contrast"] == name]
        return md_table(sub, HEADLINE + ["KEGG_TIGHT_JUNCTION", "GOBP_TIGHT_JUNCTION_ORGANIZATION_NO_CLDN4"])

    cohort_rows = [
        "| cohort | n Q1 / n Q4 | IFN-γ NES (p) | APM NES (p) | TJ NES (p) |",
        "|---|---|---:|---:|---:|",
    ]
    for cohort in COHORTS:
        key = f"q1_vs_q4_{cohort}"
        sub = ctx["all"][ctx["all"]["contrast"] == key]
        if sub.empty:
            cohort_rows.append(f"| {cohort} | Q4 n<3 | skipped | skipped | skipped |")
            continue
        def cell(term: str) -> str:
            r = sub[sub["term"] == term].iloc[0]
            return f"{r.nes:+.3f} ({fmt_p(r.nom_p)})"
        n1 = int(sub.iloc[0]["n_q1"]) if "n_q1" in sub.columns else ""
        # n stored on contrast table
        meta_c = ctx["cohort_n"][cohort]
        cohort_rows.append(
            f"| {cohort} | {meta_c['n_q1']} / {meta_c['n_q4']} | "
            f"{cell('HALLMARK_INTERFERON_GAMMA_RESPONSE')} | "
            f"{cell('CUSTOM_MHC_I_ANTIGEN_PRESENTATION')} | "
            f"{cell('KEGG_TIGHT_JUNCTION_NO_CLDN4')} |"
        )
    loo_rows = [
        "| dropped | IFN-γ NES (p) | APM NES (p) | TJ NES (p) |",
        "|---|---:|---:|---:|",
    ]
    for cohort in COHORTS:
        sub = ctx["all"][ctx["all"]["contrast"] == f"loo_drop_{cohort}"]
        def cell(term: str) -> str:
            r = sub[sub["term"] == term].iloc[0]
            return f"{r.nes:+.3f} ({fmt_p(r.nom_p)})"
        loo_rows.append(
            f"| {cohort} | {cell('HALLMARK_INTERFERON_GAMMA_RESPONSE')} | "
            f"{cell('CUSTOM_MHC_I_ANTIGEN_PRESENTATION')} | "
            f"{cell('KEGG_TIGHT_JUNCTION_NO_CLDN4')} |"
        )

    text = f"""# FINDING — concordant-4 CLDN4 quartile GSEA is an observational proxy, not a knockdown

**NOT a true KD.** No CRISPR, no siRNA, no CLDN4-loss culture. The split is the endogenous lowest versus highest malignant CLDN4 %pos quartile in the concordant-4 cohorts (GSE123902, GSE131907, GSE205335, GSE189357). Patient, donor, or sample is the unit. Positive NES means the gene set is higher in the **lowest** CLDN4 quartile.

The private comparator is a direction supplied for this wave, not a statistic recomputed here: **IFN up after CLDN4 loss**. The private knockdown matrix is not in this repository.

Pre-specified match rule, IFN-γ only: on this lowest-versus-highest rank, NES > 0 and nominal p < 0.05 is the **same direction** as that private result. NES < 0 and nominal p < 0.05 is the **opposite direction**. Otherwise the proxy does not resolve the direction. APM and tight junction are reported beside it. They were not given a private direction.

---

## Call

| readout | this proxy | private KD |
|---|---|---|
| IFN-γ | NES {ifn.nes:+.3f}, nominal p {fmt_p(ifn.nom_p)}, headline FDR {fmt_p(ifn.fdr_headline)} | IFN up after loss |
| direction | **{call}** | — |
| IFN-α | NES {ifna.nes:+.3f}, p {fmt_p(ifna.nom_p)}, FDR {fmt_p(ifna.fdr_headline)} | not supplied |
| MHC-I / APM | NES {apm.nes:+.3f}, p {fmt_p(apm.nom_p)}, FDR {fmt_p(apm.fdr_headline)} | not supplied |
| KEGG TJ, CLDN4 removed | NES {tj.nes:+.3f}, p {fmt_p(tj.nom_p)}, FDR {fmt_p(tj.fdr_headline)} — not resolved | not supplied |

The IFN leading edge is not a clean tumor-cell ISG list. The top of the lowest-versus-highest rank is `{top_genes}`. IFN-γ's leading edge includes `GZMA`, `CD69`, `CD86`, and `LCP2` along with `STAT1`, `TAP1`, and `PSMB9`. That mix fits immune transcripts inside the malignant gate, in line with the locked concordant-4 result that CLDN4-high tumors carry fewer T/NK cells. A tumor-cell knockdown is a different experiment. The IFN sign matches the private KD. This proxy does not identify the mechanism.

Honest GSEA n is **{info['n_q1']} lowest vs {info['n_q4']} highest** malignant pseudobulks (not the T/NK N=65). Genes ranked: **{info['n_genes']}**. CLDN4 itself, on the Q4-minus-Q1 coefficient this rank negates, is logFC {cldn4['logFC']:+.3f}, p {fmt_p(cldn4['p'])}. That is the split-gene check, not a pathway result.

Same direction means the low-CLDN4 quartile is the IFN-high side, which is the sign the private knockdown has. It does not mean this quartile split replicates the knockdown, estimates a knockdown effect, or replaces it.

---

## Locked design

| Item | Choice |
|---|---|
| What this is | Observational proxy |
| What this is not | A CLDN4 knockdown, knockout, or siRNA |
| Cohorts | GSE123902 + GSE131907 + GSE205335 + GSE189357 |
| Not included | GSE148071, GSE127465, GSE207422, GSE154826, CD45+ or T-only sets |
| Malignant matrix | Existing UMI-sums (marker-malignant in 123902 and 189357; author-malignant in 131907 and 205335) |
| Quartile | Within-cohort malignant CLDN4 %pos, rank then qcut. Q2 and Q3 are not in the binary contrast |
| Model | log2(TMM-CPM+1), OLS, cohort covariates, coefficient CLDN4_Q4 |
| Rank | −t of that coefficient. Positive = higher in Q1 |
| Engine | Weighted KS p=1, 1000 gene-set permutations, seed 42. Same `gsea_core.py` as the GSE68465 prerank |
| Headline sets | Hallmark IFN-γ, Hallmark IFN-α, custom MHC-I/APM (21 genes), KEGG tight junction with CLDN4 removed |
| FDR | BH inside those four sets, per rank |
| GSE189357 | Q4 n=2, so no single-cohort GSEA. TD6 and TD9 stay in the stack |
| P4001 | Not in the GSE205335 malignant UMI-sum |

---

## Headline NES (lowest vs highest)

{md_table(h, HEADLINE)}

Leading edge (first genes on the low-CLDN4 side of the rank):

- IFN-γ: `{lead_ifn}`
- MHC-I / APM: `{lead_apm}`
- KEGG TJ without CLDN4: `{lead_tj}`

---

## Per cohort

GSE189357 is not tested alone. The other three are thin; they are context for the stack, not separate claims. IFN-γ stays positive in each testable cohort. GSE205335 has the largest IFN-γ NES. Dropping it leaves IFN-γ NES positive (see the leave-one-out table). GSE205335 is also the cohort where KEGG TJ flips positive on this rank (higher TJ in the low-CLDN4 arm). The stacked TJ result does not survive that, and it is not a TJ claim.

{chr(10).join(cohort_rows)}

Leave-one-cohort-out of the stacked rank:

{chr(10).join(loo_rows)}

---

## Sensitivities

CLDN4 removed from the ranked list. The primary TJ set already excluded CLDN4; this also stops the split gene from occupying the bottom of the rank.

{sens_block("stacked_q1_vs_q4_drop_CLDN4")}

Same primary rank, with the circular KEGG TJ set that still contains CLDN4, and with GO tight-junction organization (CLDN4 removed). GO organization is not in the four-set FDR. On this quartile rank its nominal p is 0.054. That is not a TJ hit.

{sens_block("stacked_q1_vs_q4")}

Continuous malignant CLDN4 %pos (n={ctx['cont_n']} units in the count matrices, cohort covariates). Rank = −t of the %pos z-score, so positive still means higher when CLDN4 is lower. IFN-γ and APM stay positive. KEGG TJ stays unresolved. GO organization reaches nominal p < 0.05 on this continuous rank only; it was not a pre-specified headline set.

{sens_block("continuous_low_when_CLDN4_high")}

---

## What this is not

- Not a true CLDN4 knockdown, and not a substitute for the private KD.
- Not the public CLDN4/TACSTD2 KD experiments (those remain a separate result: they do not cleanly open IFN/APM).
- Not a re-audit of the concordant-4 T/NK ρ = −0.53.
- Not evidence that low CLDN4 causes IFN. The quartile is endogenous. Malignant UMI-sums can still carry immune transcripts. Histology and cohort are only partly held by the cohort covariate.
- Not a claim about APM or tight junction in the private KD. Only the IFN direction was supplied.
- Not GSE148071, GSE127465, GSE207422, or GSE154826.

## Files

- `tables/gsea_headline.tsv` — primary four sets
- `tables/gsea_all.tsv` — primary, per cohort, leave-one-out, drop-CLDN4, continuous
- `tables/direction_vs_private_kd.tsv`
- `tables/rank_low_vs_high.tsv` — gene, Q4-vs-Q1 t, lowest-vs-highest stat
- `tables/quartile_membership.tsv`
- `figures/fig_headline_nes.png` — NES bar
- `figures/fig_enrichment_curves.png`
- `figures/fig_per_cohort_nes.png`

Reproduce:

```bash
python3 methods/concordant4_cldn4_kd_match_gsea/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(text)
    (HERE / "README.md").write_text(
        "Observational proxy only: concordant-4 malignant lowest vs highest "
        "CLDN4 quartile, prerank GSEA of IFN / APM / tight junction. "
        "Not a knockdown. See FINDING.md.\n\n"
        "```bash\npython3 methods/concordant4_cldn4_kd_match_gsea/analyze.py\n```\n"
    )


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    meta = load_units()
    parts = load_counts()
    in_mat = set()
    for cts in parts.values():
        in_mat.update(cts.columns.astype(str))
    meta["in_count_matrix"] = meta["patient"].isin(in_mat)
    for cohort, (q1, q4) in EXPECTED_TAILS.items():
        sub = meta[(meta["cohort"] == cohort) & meta["in_count_matrix"]]
        got1 = set(sub.loc[sub["quartile"] == "Q1", "patient"])
        got4 = set(sub.loc[sub["quartile"] == "Q4", "patient"])
        if got1 != q1 or got4 != q4:
            raise SystemExit(f"quartile mismatch {cohort}: Q1 {got1} Q4 {got4}")
    meta.to_csv(TABLES / "quartile_membership.tsv", sep="\t", index=False)

    counts = combine(parts)
    sets = load_sets()
    terms = HEADLINE + ["KEGG_TIGHT_JUNCTION", "GOBP_TIGHT_JUNCTION_ORGANIZATION_NO_CLDN4"]

    de, info = prepare_q4(counts, meta, None)
    if info["skipped"] or de.empty:
        raise SystemExit(f"stacked contrast skipped: {info}")
    if info["n_q1"] != 18 or info["n_q4"] != 16:
        raise SystemExit(f"unexpected stacked n: {info}")
    cldn4 = de[de["gene"] == "CLDN4"].iloc[0]
    if not (cldn4["logFC"] > 1.0 and cldn4["p"] < 0.05):
        raise SystemExit(f"CLDN4 split check failed: logFC={cldn4['logFC']} p={cldn4['p']}")

    rank = low_vs_high_rank(de)
    de_out = de.copy()
    de_out["stat_low_vs_high"] = -de_out["t"]
    de_out.sort_values("stat_low_vs_high", ascending=False).to_csv(
        TABLES / "rank_low_vs_high.tsv", sep="\t", index=False
    )

    frames = [run_sets(rank, sets, terms, "stacked_q1_vs_q4")]
    frames.append(run_sets(low_vs_high_rank(de, {"CLDN4"}), sets, terms, "stacked_q1_vs_q4_drop_CLDN4"))

    cohort_n = {}
    for cohort in COHORTS:
        de_c, info_c = prepare_q4(parts[cohort], meta, cohort)
        cohort_n[cohort] = info_c
        if info_c["skipped"]:
            continue
        g = run_sets(low_vs_high_rank(de_c), sets, HEADLINE, f"q1_vs_q4_{cohort}")
        g["n_q1"] = info_c["n_q1"]
        g["n_q4"] = info_c["n_q4"]
        frames.append(g)

    for drop in COHORTS:
        keep = {c: parts[c] for c in COHORTS if c != drop}
        de_l, info_l = prepare_q4(combine(keep), meta, None)
        g = run_sets(low_vs_high_rank(de_l), sets, HEADLINE, f"loo_drop_{drop}")
        g["n_q1"] = info_l["n_q1"]
        g["n_q4"] = info_l["n_q4"]
        frames.append(g)

    de_z, info_z = prepare_continuous(counts, meta)
    frames.append(run_sets(low_vs_high_rank(de_z), sets, terms, "continuous_low_when_CLDN4_high"))

    all_g = pd.concat(frames, ignore_index=True)
    all_g.to_csv(TABLES / "gsea_all.tsv", sep="\t", index=False)
    headline = all_g[all_g["contrast"] == "stacked_q1_vs_q4"].copy()
    headline[headline["term"].isin(HEADLINE)].to_csv(TABLES / "gsea_headline.tsv", sep="\t", index=False)

    ifn = headline[headline["term"] == "HALLMARK_INTERFERON_GAMMA_RESPONSE"].iloc[0]
    call = direction_call(float(ifn.nes), float(ifn.nom_p))
    pd.DataFrame(
        [
            {
                "proxy": "concordant-4 malignant lowest vs highest CLDN4 quartile",
                "label": "NOT a true KD",
                "private_kd_ifn": PRIVATE_IFN,
                "private_matrix_in_repo": False,
                "ifn_gamma_nes_low_vs_high": float(ifn.nes),
                "ifn_gamma_nom_p": float(ifn.nom_p),
                "ifn_gamma_fdr": float(ifn.fdr_headline),
                "direction_vs_private_kd": call,
                "n_q1": info["n_q1"],
                "n_q4": info["n_q4"],
                "match_rule": "IFN-γ NES>0 and nom p<0.05 on the Q1-vs-Q4 rank = same direction as IFN up after loss",
            }
        ]
    ).to_csv(TABLES / "direction_vs_private_kd.tsv", sep="\t", index=False)

    plot_headline(headline, FIGS / "fig_headline_nes")
    plot_curves(rank, sets, FIGS / "fig_enrichment_curves")
    plot_cohorts(all_g, FIGS / "fig_per_cohort_nes")
    write_finding(
        {
            "headline": headline,
            "all": all_g,
            "info": info,
            "cldn4": {"logFC": float(cldn4["logFC"]), "p": float(cldn4["p"])},
            "cohort_n": cohort_n,
            "cont_n": info_z["n"],
            "top_genes": de_out.sort_values("stat_low_vs_high", ascending=False)["gene"].head(8).tolist(),
        }
    )
    print("CLDN4 Q4 logFC", float(cldn4["logFC"]), "p", float(cldn4["p"]))
    print(headline[headline["term"].isin(HEADLINE)][["label", "nes", "nom_p", "fdr_headline", "n_set_in_rank"]].to_string(index=False))
    print("DIRECTION", call)


if __name__ == "__main__":
    main()

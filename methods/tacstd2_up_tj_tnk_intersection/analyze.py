#!/usr/bin/env python3
"""Story gap after PR #741.

Among genes UP in TACSTD2-high malignant cells on the #741 list, which are
tight-junction family members, and which of those individually predict low T/NK?

Primary UP list = the #741 ORA query (p<0.01 and logFC>0.25; TACSTD2 held out).
TJ family = the #741 family (KEGG tight junction ∪ GO tight-junction
organization ∪ the custom epithelial-adhesion panel), TACSTD2 held out.
"Predicts low T/NK" = DerSimonian–Laird meta Spearman ρ < 0 and p < 0.05
of malignant log2(TMM-CPM+1) vs the same-unit frac_tnk. That is the #741
meta, applied gene by gene. BH-FDR is computed inside the screened set.

Concordant-4 only. P4001 has no count column, so expression n = 64.
Never fabricates a %pos column: only CLDN4 has one in the locked unit tables.
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

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
COMBO = "+".join(COHORTS)
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
CORE = ["CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "TJP1"]

# Anchors copied from PR #741 tables. The script stops if a recompute drifts.
ANCHOR_TACSTD2_LOGFC = 3.9040230854753495
ANCHOR_TACSTD2_P = 7.708786908785916e-07
ANCHOR_CLDN4_LOGFC = 1.6323196517230478
ANCHOR_CLDN4_P = 0.012564175272030588
ANCHOR_N_UP = 274
ANCHOR_UP_TJ = {
    "SYNPO",
    "EPHA2",
    "NECTIN4",
    "ERBB2",
    "CLDN23",
    "DSG2",
    "GRHL2",
    "LSR",
    "AMOTL2",
    "JUP",
    "PLEC",
}
ANCHOR_CLDN4_PCT_RHO = -0.5228600495638758
ANCHOR_CLDN4_PCT_P = 2.859572051598e-05


def load_units() -> pd.DataFrame:
    rows = []
    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    for _, r in el.iterrows():
        pct = float(r["mal_CLDN4_pct"])
        if pct <= 1.5:
            pct *= 100.0
        rows.append(
            {
                "patient": str(r["patient"]),
                "cohort": "GSE123902",
                "unit": "donor",
                "cldn4_pct": pct,
                "frac_tnk": float(r["frac_tnk"]),
            }
        )

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    for _, r in mal.iterrows():
        pct = float(r["mal_CLDN4_pct"])
        if pct <= 1.5:
            pct *= 100.0
        rows.append(
            {
                "patient": str(r["sample"]),
                "cohort": "GSE131907",
                "unit": "sample",
                "cldn4_pct": pct,
                "frac_tnk": float(r["frac_tnk"]),
            }
        )

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    for _, r in d.iterrows():
        pct = float(r["mal_CLDN4_pct_pos"])
        if pct <= 1.5:
            pct *= 100.0
        rows.append(
            {
                "patient": str(r["patient"]),
                "cohort": "GSE205335",
                "unit": "patient",
                "cldn4_pct": pct,
                "frac_tnk": float(r["frac_tnk"]),
            }
        )

    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    for _, r in el.iterrows():
        pct = float(r["mal_CLDN4_pct"])
        if pct <= 1.5:
            pct *= 100.0
        rows.append(
            {
                "patient": str(r["patient"]),
                "cohort": "GSE189357",
                "unit": "patient",
                "cldn4_pct": pct,
                "frac_tnk": float(r["frac_tnk"]),
            }
        )
    return pd.DataFrame(rows)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str)
    return df.groupby(df.index).sum()


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
        w = (obs_lib - obs[keep]) / (obs_lib * obs[keep]) + (ref_lib - ref_c[keep]) / (
            ref_lib * ref_c[keep]
        )
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
        tmm = float(np.average(m[trim], weights=1.0 / w[trim]))
        factors[col] = 2**tmm
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def load_merged_counts(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mats = []
    for c in COHORTS:
        mats.append(read_counts(DATA / f"{c}_malignant_counts.tsv.gz"))
    all_genes = sorted(set().union(*[set(m.index) for m in mats]))
    pieces = []
    for c, m in zip(COHORTS, mats):
        want = [p for p in meta.loc[meta["cohort"] == c, "patient"] if p in m.columns]
        if not want:
            continue
        pieces.append(m.reindex(all_genes).fillna(0.0)[want])
    counts = pd.concat(pieces, axis=1)
    counts = counts.loc[:, ~counts.columns.duplicated()]
    meta2 = meta.loc[meta["patient"].isin(counts.columns)].copy()
    return counts, meta2


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    ranked = np.nan_to_num(p[order], nan=1.0)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def ols_de(logcpm: pd.DataFrame, design: pd.DataFrame, coef: str) -> pd.DataFrame:
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
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "t": t,
            "p": pv,
            "se": se,
            "df": df_res,
        }
    )
    out["fdr"] = bh(out["p"].values)
    return out.sort_values("p")


def spearman(x, y) -> tuple[float, float, int]:
    xa = np.asarray(list(x), dtype=float)
    ya = np.asarray(list(y), dtype=float)
    mask = np.isfinite(xa) & np.isfinite(ya)
    xa, ya = xa[mask], ya[mask]
    n = int(xa.size)
    if n < 4 or np.unique(xa).size < 2 or np.unique(ya).size < 2:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def fisher_z(rho: float) -> float:
    r = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(r))


def random_effects_dl(rhos: list[float], ns: list[int]) -> dict:
    z = np.array([fisher_z(r) for r in rhos], dtype=float)
    v = np.array([1.0 / (n - 3) for n in ns], dtype=float)
    ok = np.isfinite(z) & np.isfinite(v) & (v > 0)
    z, v = z[ok], v[ok]
    n_ok = np.array(ns, dtype=float)[ok]
    k = int(z.size)
    if k == 0:
        return {"k": 0, "pooled_rho": float("nan"), "p": float("nan"), "I2": float("nan"),
                "ci95_rho": [float("nan"), float("nan")], "n_patients_total": 0}
    w = 1.0 / v
    z_fe = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - z_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else float("nan")
    tau2 = max(0.0, (q - df) / c) if k > 1 and c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se = float(1.0 / math.sqrt(float(np.sum(w_re))))
    z_stat = z_re / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z_stat))) if math.isfinite(z_stat) else float("nan")
    lo, hi = z_re - 1.96 * se, z_re + 1.96 * se
    i2 = max(0.0, (q - df) / q) * 100.0 if k > 1 and q > 0 else 0.0
    return {
        "k": k,
        "n_patients_total": int(n_ok.sum()),
        "pooled_rho": float(np.tanh(z_re)),
        "p": p,
        "I2": i2,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
    }


def build_sets(a8: dict) -> dict[str, set[str]]:
    kegg = {str(g).upper() for g in a8["sets"]["KEGG_TIGHT_JUNCTION"]}
    gobp = {str(g).upper() for g in a8["sets"]["GOBP_TIGHT_JUNCTION_ORGANIZATION"]}
    custom = {
        "CDH1", "EPCAM", "F11R", "OCLN", "TJP1", "TJP2", "TJP3", "CGN", "MARVELD2",
        "PVRL1", "NECTIN1", "NECTIN2", "NECTIN4", "JUP", "DSP", "DSG2", "PKP2",
        "CTNNA1", "CTNNB1", "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    }
    cldn = {g for g in kegg if g.startswith("CLDN")} | {
        "CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLDN12", "CLDN18",
    }
    tj = (kegg | gobp | custom) - {"TACSTD2"}
    return {"kegg": kegg, "gobp": gobp, "custom": custom, "cldn": cldn, "tj": tj}


def source_label(gene: str, sets: dict[str, set[str]]) -> str:
    parts = []
    if gene in sets["kegg"]:
        parts.append("KEGG_TJ")
    if gene in sets["gobp"]:
        parts.append("GOBP_TJ_ORG")
    if gene in sets["custom"]:
        parts.append("CUSTOM_ADHESION")
    return "|".join(parts)


def fmt_p(p: float) -> str:
    if not math.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r: float) -> str:
    if not math.isfinite(r):
        return "NA"
    return f"{r:+.3f}"


def gene_vs_tnk(gene: str, lc: pd.DataFrame, meta: pd.DataFrame) -> dict | None:
    if gene not in lc.index:
        return None
    x = lc.loc[gene]
    rhos, ns, members = [], [], []
    cohort_rho = {}
    for c in COHORTS:
        m = meta.loc[meta["cohort"] == c]
        rho, p, n = spearman(x.reindex(m["patient"]), m["frac_tnk"])
        cohort_rho[c] = {"rho": rho, "p": p, "n": n}
        if n >= 4 and math.isfinite(rho):
            rhos.append(rho)
            ns.append(n)
            members.append(c)
    if len(rhos) < 2:
        return None
    dl = random_effects_dl(rhos, ns)
    n_neg = int(sum(1 for r in rhos if r < 0))
    return {
        "rho": dl["pooled_rho"],
        "p": dl["p"],
        "I2": dl["I2"],
        "ci_lo": dl["ci95_rho"][0],
        "ci_hi": dl["ci95_rho"][1],
        "k": dl["k"],
        "n": dl["n_patients_total"],
        "n_neg": n_neg,
        "n_cohorts": len(rhos),
        "member_rhos": ",".join(f"{c}:{r:.3f}" for c, r in zip(members, rhos)),
        "cohorts": cohort_rho,
        "predicts_low_tnk": bool(dl["pooled_rho"] < 0 and dl["p"] < 0.05),
    }


def stacked_q(meta: pd.DataFrame, values: pd.Series) -> dict:
    q1, q4 = [], []
    for c in COHORTS:
        m = meta.loc[meta["cohort"] == c].set_index("patient")
        xv = values.reindex(m.index).astype(float)
        if int(xv.nunique(dropna=True)) < 4 or len(xv.dropna()) < 8:
            continue
        try:
            qs = assign_quartiles(xv.dropna())
        except ValueError:
            continue
        q1.extend(m.loc[qs.index[qs == "Q1"], "frac_tnk"].tolist())
        q4.extend(m.loc[qs.index[qs == "Q4"], "frac_tnk"].tolist())
    if len(q1) < 3 or len(q4) < 3:
        return {"mwu_p": float("nan"), "r_rb": float("nan"), "delta_median": float("nan"),
                "n_q1": len(q1), "n_q4": len(q4)}
    u, p = stats.mannwhitneyu(np.asarray(q4), np.asarray(q1), alternative="two-sided")
    r_rb = (2.0 * float(u)) / (len(q4) * len(q1)) - 1.0
    return {
        "mwu_p": float(p),
        "r_rb": float(r_rb),
        "delta_median": float(np.median(q4) - np.median(q1)),
        "n_q1": len(q1),
        "n_q4": len(q4),
    }


def screen_genes(genes: list[str], lc: pd.DataFrame, meta: pd.DataFrame, de: pd.DataFrame,
                 sets: dict[str, set[str]]) -> pd.DataFrame:
    de_ix = de.set_index("gene")
    rows = []
    for gene in genes:
        hit = gene_vs_tnk(gene, lc, meta)
        if hit is None:
            continue
        st = stacked_q(meta, lc.loc[gene])
        rec = {
            "gene": gene,
            "tj_source": source_label(gene, sets),
            "classical_tj": gene in (sets["kegg"] | sets["gobp"]),
            "logFC": float(de_ix.loc[gene, "logFC"]) if gene in de_ix.index else float("nan"),
            "de_p": float(de_ix.loc[gene, "p"]) if gene in de_ix.index else float("nan"),
            "de_fdr": float(de_ix.loc[gene, "fdr"]) if gene in de_ix.index else float("nan"),
        }
        rec.update({k: hit[k] for k in
                    ["rho", "p", "I2", "ci_lo", "ci_hi", "k", "n", "n_neg", "n_cohorts",
                     "member_rhos", "predicts_low_tnk"]})
        rec.update(st)
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr_screen"] = bh(out["p"].values)
    out["predicts_low_tnk_fdr"] = (out["rho"] < 0) & (out["fdr_screen"] < 0.05)
    return out.sort_values(["predicts_low_tnk", "rho"], ascending=[False, True])


def assert_close(name: str, got: float, exp: float, tol: float = 1e-6) -> None:
    if not math.isfinite(got) or abs(got - exp) > tol * max(1.0, abs(exp)):
        raise SystemExit(f"anchor drift {name}: got {got} expected {exp}")


def write_finding(summary: dict, up_tj: pd.DataFrame, nominal: pd.DataFrame,
                  tj_low: pd.DataFrame, cldn: pd.DataFrame) -> None:
    def row_md(r: pd.Series) -> str:
        call = "yes" if bool(r["predicts_low_tnk"]) else "no"
        return (
            f"| {r['gene']} | {r['logFC']:+.3f} | {fmt_p(r['de_p'])} | {r['tj_source']} | "
            f"{fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {fmt_rho(r['ci_lo'])} to {fmt_rho(r['ci_hi'])} | "
            f"{int(r['n_neg'])}/{int(r['n_cohorts'])} | {r['I2']:.1f}% | {call} |"
        )

    lines = []
    lines.append("# TACSTD2-high UP ∩ TJ ∩ low T/NK")
    lines.append("")
    lines.append("Story gap after PR #741. Question: among genes up in TACSTD2-high")
    lines.append("malignant cells on the #741 list, which are tight-junction members")
    lines.append("**and** individually predict low T/NK?")
    lines.append("")
    lines.append("Locked cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.")
    lines.append("Split and the UP list are TACSTD2, not CLDN4. This page does not")
    lines.append("replace the locked CLDN4 %pos vs T/NK result (ρ = −0.531, n = 65).")
    lines.append("")
    lines.append("## Call")
    lines.append("")
    lines.append("| piece | rule |")
    lines.append("|---|---|")
    lines.append("| UP list | #741 ORA query: p<0.01 and logFC>0.25 on the TACSTD2 Q4 vs Q1 OLS. TACSTD2 held out. FDR arm of that DEG is empty, which is why #741 used this threshold. |")
    lines.append("| TJ member | #741 TJ family: KEGG_TIGHT_JUNCTION ∪ GOBP_TIGHT_JUNCTION_ORGANIZATION ∪ CUSTOM_EPITHELIAL_ADHESION. TACSTD2 held out. |")
    lines.append("| Predicts low T/NK | Within-cohort Spearman of malignant log2(TMM-CPM+1) vs frac_tnk, DerSimonian–Laird meta, **ρ < 0 and p < 0.05**. BH-FDR is reported inside the screened set and is not required for the nominal call. |")
    lines.append("| Unit | Patient / donor / sample. Expression n = 64 (P4001 has T/NK and CLDN4 %pos but no count column). |")
    lines.append("")
    lines.append("A sensitivity UP list (p<0.05 and logFC>0) is reported so CLDN4, which")
    lines.append("misses the #741 threshold at p = 0.013, is not dropped in silence.")
    lines.append("There is no per-cell matrix here, so this screen is pseudobulk expression.")
    lines.append("The only %pos column in the locked tables is CLDN4.")
    lines.append("")
    lines.append("## Set sizes")
    lines.append("")
    lines.append("| set | n |")
    lines.append("|---|---:|")
    for k, v in summary["counts"].items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("## Primary table: #741 UP ∩ TJ, each gene vs T/NK")
    lines.append("")
    lines.append("Positive logFC = higher in TACSTD2 Q4. Positive ρ = higher malignant expression with higher T/NK.")
    lines.append("")
    lines.append("| gene | logFC | DE p | TJ source | meta ρ | meta p | 95% CI | cohorts ρ<0 | I² | predicts low T/NK |")
    lines.append("|---|---:|---:|---|---:|---:|---|---:|---:|---|")
    show = up_tj.sort_values("rho")
    for _, r in show.iterrows():
        lines.append(row_md(r))
    lines.append("")
    n_call = int(up_tj["predicts_low_tnk"].sum()) if len(up_tj) else 0
    n_fdr = int(up_tj["predicts_low_tnk_fdr"].sum()) if len(up_tj) else 0
    lines.append(
        f"**Three-way intersection (#741 UP ∩ TJ ∩ predicts low T/NK) = {n_call} genes.** "
        f"FDR < 0.05 inside this screen: {n_fdr}."
    )
    lines.append("")
    if n_call == 0:
        lines.append(
            "No gene that is on the #741 UP list and in the TJ family individually "
            "predicts low T/NK. The 11 meta ρ values sit between "
            f"{fmt_rho(float(show['rho'].min()))} and {fmt_rho(float(show['rho'].max()))}; "
            f"the smallest meta p is {fmt_p(float(show['p'].min()))}."
        )
        lines.append("")
    classical = up_tj[up_tj["classical_tj"]]
    lines.append(
        f"Restricting TJ to KEGG ∪ GO organization (dropping adhesion-only members) "
        f"leaves {int(classical.shape[0])} genes and "
        f"{int(classical['predicts_low_tnk'].sum())} low-T/NK calls. "
        "The empty intersection does not depend on the custom adhesion panel."
    )
    lines.append("")
    lines.append(
        f"The same call on the entire #741 UP list, TJ or not: "
        f"{summary['up_list_low_tnk']} / {summary['counts']['#741 UP list (p<0.01, logFC>0.25, TACSTD2 out)']} genes "
        f"have meta ρ < 0 and p < 0.05. "
        f"The most negative meta ρ inside the UP list is {summary['up_most_neg_gene']} "
        f"(ρ = {fmt_rho(summary['up_most_neg_rho'])}, p = {fmt_p(summary['up_most_neg_p'])}). "
        f"Significant UP-list associations that do exist are in the other direction "
        f"({summary['up_list_high_tnk']} genes with ρ > 0 and p < 0.05; "
        f"best FDR among those is {fmt_p(summary['up_pos_best_fdr'])})."
    )
    lines.append("")
    lines.append(
        "NECTIN4 (and NECTIN2 on the nominal list) have no rank correlation in "
        "GSE131907 because expression does not vary there, so those metas use 3 cohorts."
    )
    lines.append("")
    lines.append("## Sensitivity: nominal UP (p<0.05, logFC>0) ∩ TJ")
    lines.append("")
    lines.append(
        f"n = {len(nominal)}. Genes with meta ρ < 0 and p < 0.05: "
        f"**{int(nominal['predicts_low_tnk'].sum())}**. "
        f"FDR < 0.05 and ρ < 0: **{int(nominal['predicts_low_tnk_fdr'].sum())}**."
    )
    lines.append("CLDN4 is in this nominal list and is absent from the #741 UP list.")
    lines.append("")
    lines.append("| gene | on #741 UP list | logFC | DE p | meta ρ | meta p | cohorts ρ<0 | predicts low T/NK |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|")
    nom_show = nominal.sort_values("rho")
    up_names = set(up_tj["gene"])
    for _, r in nom_show.iterrows():
        on = "yes" if r["gene"] in up_names else "no"
        call = "yes" if bool(r["predicts_low_tnk"]) else "no"
        lines.append(
            f"| {r['gene']} | {on} | {r['logFC']:+.3f} | {fmt_p(r['de_p'])} | "
            f"{fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {int(r['n_neg'])}/{int(r['n_cohorts'])} | {call} |"
        )
    lines.append("")
    c16 = nominal.loc[nominal["gene"] == "CLDN16"]
    if len(c16):
        r = c16.iloc[0]
        lines.append(
            f"CLDN16 is the one gene in this nominal table with meta p < 0.05 "
            f"(ρ = {fmt_rho(r['rho'])}, p = {fmt_p(r['p'])}). The direction is higher T/NK, "
            "so it stays out of the low-T/NK intersection."
        )
        lines.append("")
    lines.append("## Claudin rows, and the CLDN4 %pos reference")
    lines.append("")
    lines.append(
        "Every claudin in the #741 claudin panel that is present in the expression "
        "matrix is listed. CLDN4 %pos is a different measurement (cell-level percent "
        "positive), recomputed on these 64 units with the same meta. It is not a "
        "member of the UP list."
    )
    lines.append("")
    lines.append("| gene | on #741 UP list | logFC | DE p | meta ρ | meta p | cohorts ρ<0 | readout |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|")
    for _, r in cldn.sort_values(["readout", "rho"]).iterrows():
        on = "yes" if bool(r["on_741_up"]) else "no"
        if r["readout"] == "pct_pos":
            lines.append(
                f"| {r['gene']} | {on} | — | — | {fmt_rho(r['rho'])} | {fmt_p(r['p'])} | "
                f"{int(r['n_neg'])}/{int(r['n_cohorts'])} | %pos |"
            )
        else:
            lines.append(
                f"| {r['gene']} | {on} | {r['logFC']:+.3f} | {fmt_p(r['de_p'])} | "
                f"{fmt_rho(r['rho'])} | {fmt_p(r['p'])} | {int(r['n_neg'])}/{int(r['n_cohorts'])} | pseudobulk |"
            )
    lines.append("")
    lines.append("## TJ genes with a nominal low-T/NK call are not the UP list")
    lines.append("")
    lines.append(
        f"Across {summary['counts']['TJ family genes in the expression matrix']} TJ-family genes "
        f"in the expression matrix, {len(tj_low)} have meta ρ < 0 and p < 0.05. "
        f"Of those, {int(tj_low['on_741_up'].sum()) if len(tj_low) else 0} are on the #741 UP list "
        f"and {int(tj_low['on_nominal_up'].sum()) if len(tj_low) else 0} are on the nominal UP list. "
        f"None survives BH-FDR across the TJ family "
        f"(smallest FDR = {fmt_p(summary['tj_low_best_fdr'])}). "
        "They are recorded so the empty intersection is visible as disjoint sets, "
        "not as an untested corner. They are not an exclusion signature."
    )
    lines.append("")
    if len(tj_low):
        lines.append("| gene | logFC | DE p | on #741 UP | meta ρ | meta p | FDR in TJ family | cohorts ρ<0 |")
        lines.append("|---|---:|---:|---|---:|---:|---:|---:|")
        for _, r in tj_low.sort_values("rho").iterrows():
            lines.append(
                f"| {r['gene']} | {r['logFC']:+.3f} | {fmt_p(r['de_p'])} | "
                f"{'yes' if r['on_741_up'] else 'no'} | {fmt_rho(r['rho'])} | {fmt_p(r['p'])} | "
                f"{fmt_p(r['fdr_tj'])} | {int(r['n_neg'])}/{int(r['n_cohorts'])} |"
            )
        lines.append("")
    lines.append("## Lead-in")
    lines.append("")
    lines.append(
        "The TACSTD2-high malignant program and the low-T/NK direction do not meet "
        "on a TJ gene. The claudin that clears the #741 UP list is CLDN23 "
        f"(logFC {fmt_rho(summary['cldn23_logfc']).replace('+','+')}, DE p = {fmt_p(summary['cldn23_p'])}), "
        f"and its meta ρ vs T/NK is {fmt_rho(summary['cldn23_rho'])} "
        f"(p = {fmt_p(summary['cldn23_meta_p'])}). "
        "CLDN4 misses that UP list (logFC "
        f"{summary['cldn4_logfc']:+.3f}, DE p = {fmt_p(summary['cldn4_de_p'])}). "
        "On the same pseudobulk scale CLDN4 is negative in "
        f"{summary['cldn4_n_neg']}/{summary['cldn4_n_cohorts']} cohorts "
        f"(ρ = {fmt_rho(summary['cldn4_rho'])}, p = {fmt_p(summary['cldn4_meta_p'])}, I² = {summary['cldn4_I2']:.1f}%) "
        "and does not meet p < 0.05. "
        "The measurement that predicts low T/NK is CLDN4 %pos on these units "
        f"(ρ = {fmt_rho(summary['cldn4_pct_rho'])}, p = {fmt_p(summary['cldn4_pct_p'])}, "
        f"negative in {summary['cldn4_pct_n_neg']}/4, I² = {summary['cldn4_pct_I2']:.1f}%). "
        "The screen this intersection points at is cell-level CLDN4 %pos."
    )
    lines.append("")
    lines.append("CLDN1 is nominally up (DE p = "
                 f"{fmt_p(summary['cldn1_de_p'])}) and its meta ρ vs T/NK is "
                 f"{fmt_rho(summary['cldn1_rho'])} "
                 f"({summary['cldn1_n_neg']}/{summary['cldn1_n_cohorts']} cohorts negative). "
                 "A claudin that is up with TACSTD2 is not automatically a low-T/NK gene.")
    lines.append("")
    lines.append("## Not claimed")
    lines.append("")
    lines.append("- A TJ gene from the #741 UP list excludes T/NK. None does, at this rule.")
    lines.append("- Pseudobulk CLDN4 expression is a significant T/NK predictor. It is directionally negative in 4/4 and not significant.")
    lines.append("- The nine nominal TJ hits above are an exclusion list. They miss the UP list and miss FDR.")
    lines.append("- Cell-level %pos for any gene other than CLDN4. Those columns are not in the locked tables.")
    lines.append("- Private 8KL, KD coculture, or a Visium spatial claim.")
    lines.append("- Pooling non-concordant accessions.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/tacstd2_up_tj_tnk_intersection/analyze.py")
    lines.append("```")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def forest_plot(up_tj: pd.DataFrame, cldn4_expr: dict, cldn4_pct: dict) -> None:
    block = up_tj.sort_values("rho", ascending=True).copy()
    rows = []
    for _, r in block.iterrows():
        rows.append((r["gene"], float(r["rho"]), float(r["ci_lo"]), float(r["ci_hi"]), float(r["p"]), "up_tj"))
    rows.append((
        "CLDN4 expression", float(cldn4_expr["rho"]), float(cldn4_expr["ci_lo"]),
        float(cldn4_expr["ci_hi"]), float(cldn4_expr["p"]), "cldn4_expr",
    ))
    rows.append((
        "CLDN4 %pos", float(cldn4_pct["rho"]), float(cldn4_pct["ci_lo"]),
        float(cldn4_pct["ci_hi"]), float(cldn4_pct["p"]), "cldn4_pct",
    ))
    # Most-negative UP gene at the top; CLDN4 reference rows at the bottom.
    y = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(8.6, 7.4))
    colors = {"up_tj": "#4C78A8", "cldn4_expr": "#F58518", "cldn4_pct": "#E45756"}
    handles = {}
    for yy, (lab, rho, lo, hi, p, k) in zip(y, rows):
        ax.plot([lo, hi], [yy, yy], color=colors[k], lw=1.6, solid_capstyle="round")
        handles[k], = ax.plot(rho, yy, "o", color=colors[k], ms=6, zorder=3)
        ax.text(0.70, yy, f"{rho:+.2f}  p={p:.3g}", va="center", ha="left", fontsize=7, color="#333333")
    ax.axvline(0, color="#333333", lw=0.8)
    ax.axhline(1.5, color="#bbbbbb", lw=0.6, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel("DL meta Spearman ρ versus frac T/NK (95% CI)")
    ax.set_xlim(-0.90, 1.15)
    ax.set_title("TACSTD2-high UP ∩ TJ versus T/NK")
    fig.subplots_adjust(left=0.24, right=0.98, top=0.92, bottom=0.24)
    fig.legend(
        [handles["up_tj"], handles["cldn4_expr"], handles["cldn4_pct"]],
        [
            "#741 UP ∩ TJ (0/11 with ρ<0 and p<0.05)",
            "CLDN4 pseudobulk (nominal UP only)",
            "CLDN4 %pos on the same 64 units",
        ],
        frameon=False,
        fontsize=8,
        loc="lower center",
        bbox_to_anchor=(0.58, 0.01),
    )
    fig.savefig(FIGS / "forest_up_tj_vs_tnk.png", dpi=160)
    fig.savefig(FIGS / "forest_up_tj_vs_tnk.pdf")
    plt.close(fig)


def main() -> None:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = build_sets(a8)
    units = load_units()
    if len(units) != 65:
        raise SystemExit(f"expected 65 locked units, got {len(units)}")
    counts_raw, meta = load_merged_counts(units)
    missing = sorted(set(units["patient"]) - set(meta["patient"]))
    if missing != ["P4001"]:
        raise SystemExit(f"unexpected units missing from counts: {missing}")
    counts = filter_genes(counts_raw)
    lc = log_cpm(counts, tmm_norm_factors(counts))

    meta = meta.set_index("patient")
    meta["tacstd2"] = lc.loc["TACSTD2", meta.index].astype(float)
    qs = []
    for c in COHORTS:
        idx = meta.index[meta["cohort"] == c]
        qs.append(assign_quartiles(meta.loc[idx, "tacstd2"]))
    meta["quartile"] = pd.concat(qs)
    meta = meta.reset_index()

    m_de = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    design = pd.DataFrame(index=m_de["patient"])
    design["Intercept"] = 1.0
    design["TACSTD2_Q4"] = (m_de.set_index("patient")["quartile"] == "Q4").astype(float)
    for c in COHORTS:
        if c == REF:
            continue
        design[f"cohort_{c}"] = (m_de.set_index("patient")["cohort"] == c).astype(float)
    de = ols_de(lc.loc[:, m_de["patient"]], design, "TACSTD2_Q4")
    de.to_csv(TABLES / "de_q4q1_stacked.tsv.gz", sep="\t", index=False, compression="gzip")

    tac = de.loc[de["gene"] == "TACSTD2"].iloc[0]
    c4 = de.loc[de["gene"] == "CLDN4"].iloc[0]
    assert_close("TACSTD2 logFC", float(tac["logFC"]), ANCHOR_TACSTD2_LOGFC, tol=1e-5)
    assert_close("TACSTD2 p", float(tac["p"]), ANCHOR_TACSTD2_P, tol=1e-4)
    assert_close("CLDN4 logFC", float(c4["logFC"]), ANCHOR_CLDN4_LOGFC, tol=1e-5)
    assert_close("CLDN4 p", float(c4["p"]), ANCHOR_CLDN4_P, tol=1e-4)

    up_mask = (de["p"] < 0.01) & (de["logFC"] > 0.25) & (de["gene"] != "TACSTD2")
    nom_mask = (de["p"] < 0.05) & (de["logFC"] > 0) & (de["gene"] != "TACSTD2")
    up_genes = de.loc[up_mask, "gene"].tolist()
    nom_genes = de.loc[nom_mask, "gene"].tolist()
    if len(up_genes) != ANCHOR_N_UP:
        raise SystemExit(f"UP list n drifted: {len(up_genes)} vs {ANCHOR_N_UP}")

    tj_present = [g for g in sorted(sets["tj"]) if g in lc.index]
    up_tj_genes = [g for g in up_genes if g in sets["tj"]]
    if set(up_tj_genes) != ANCHOR_UP_TJ:
        raise SystemExit(f"UP ∩ TJ drifted: {sorted(up_tj_genes)}")

    up_tj = screen_genes(up_tj_genes, lc, meta, de, sets)
    # FDR must be within this 11, which screen_genes already does.
    nominal_genes = [g for g in nom_genes if g in sets["tj"]]
    nominal = screen_genes(nominal_genes, lc, meta, de, sets)
    tj_all = screen_genes(tj_present, lc, meta, de, sets)
    tj_all["on_741_up"] = tj_all["gene"].isin(up_tj_genes)
    tj_all["on_nominal_up"] = tj_all["gene"].isin(nominal_genes)
    tj_all["fdr_tj"] = tj_all["fdr_screen"]
    tj_low = tj_all[tj_all["predicts_low_tnk"]].copy()

    # Whole UP list, for the "zero in 274" statement. Expression only.
    up_screen = screen_genes(up_genes, lc, meta, de, sets)
    up_low = up_screen[up_screen["predicts_low_tnk"]]
    up_high = up_screen[(up_screen["rho"] > 0) & (up_screen["p"] < 0.05)]
    most_neg = up_screen.loc[up_screen["rho"].idxmin()]

    # Claudin panel + CLDN4 %pos reference.
    cldn_genes = [g for g in sorted(sets["cldn"]) if g in lc.index and g != "TACSTD2"]
    cldn_expr = screen_genes(cldn_genes, lc, meta, de, sets)
    cldn_rows = []
    for _, r in cldn_expr.iterrows():
        cldn_rows.append(
            {
                "gene": r["gene"],
                "readout": "pseudobulk",
                "on_741_up": r["gene"] in set(up_tj_genes),
                "on_nominal_up": r["gene"] in set(nominal_genes),
                "logFC": r["logFC"],
                "de_p": r["de_p"],
                "rho": r["rho"],
                "p": r["p"],
                "I2": r["I2"],
                "ci_lo": r["ci_lo"],
                "ci_hi": r["ci_hi"],
                "n_neg": r["n_neg"],
                "n_cohorts": r["n_cohorts"],
                "n": r["n"],
                "member_rhos": r["member_rhos"],
                "predicts_low_tnk": r["predicts_low_tnk"],
            }
        )
    pct = gene_vs_tnk_from_column(meta, "cldn4_pct")
    assert_close("CLDN4 %pos rho", pct["rho"], ANCHOR_CLDN4_PCT_RHO, tol=1e-5)
    assert_close("CLDN4 %pos p", pct["p"], ANCHOR_CLDN4_PCT_P, tol=1e-4)
    cldn_rows.append(
        {
            "gene": "CLDN4",
            "readout": "pct_pos",
            "on_741_up": False,
            "on_nominal_up": False,
            "logFC": float("nan"),
            "de_p": float("nan"),
            "rho": pct["rho"],
            "p": pct["p"],
            "I2": pct["I2"],
            "ci_lo": pct["ci_lo"],
            "ci_hi": pct["ci_hi"],
            "n_neg": pct["n_neg"],
            "n_cohorts": pct["n_cohorts"],
            "n": pct["n"],
            "member_rhos": pct["member_rhos"],
            "predicts_low_tnk": pct["predicts_low_tnk"],
        }
    )
    cldn = pd.DataFrame(cldn_rows)

    c4_hit = gene_vs_tnk("CLDN4", lc, meta)
    c23 = up_tj.loc[up_tj["gene"] == "CLDN23"].iloc[0]
    c1 = cldn_expr.loc[cldn_expr["gene"] == "CLDN1"].iloc[0]

    counts = {
        "#741 UP list (p<0.01, logFC>0.25, TACSTD2 out)": len(up_genes),
        "nominal UP (p<0.05, logFC>0, TACSTD2 out)": len(nom_genes),
        "TJ family defined (TACSTD2 out)": len(sets["tj"]),
        "TJ family genes in the expression matrix": len(tj_present),
        "#741 UP ∩ TJ": len(up_tj_genes),
        "#741 UP ∩ classical TJ (KEGG ∪ GO)": int(up_tj["classical_tj"].sum()),
        "#741 UP ∩ TJ ∩ predicts low T/NK": int(up_tj["predicts_low_tnk"].sum()),
        "nominal UP ∩ TJ": len(nominal_genes),
        "nominal UP ∩ TJ ∩ predicts low T/NK": int(nominal["predicts_low_tnk"].sum()),
        "TJ family with nominal low T/NK (ρ<0, p<0.05)": int(tj_all["predicts_low_tnk"].sum()),
        "of those, also on the #741 UP list": int(tj_low["on_741_up"].sum()) if len(tj_low) else 0,
    }

    up_tj.to_csv(TABLES / "up_tj_vs_tnk.tsv", sep="\t", index=False)
    # Explicit three-way table: header plus the rows that pass. Empty is a result.
    three = up_tj[up_tj["predicts_low_tnk"]].copy()
    three.to_csv(TABLES / "intersection_three_way.tsv", sep="\t", index=False)
    nominal.to_csv(TABLES / "nominal_up_tj_vs_tnk.tsv", sep="\t", index=False)
    tj_low.sort_values("rho").to_csv(TABLES / "tj_nominal_low_tnk_not_the_up_list.tsv", sep="\t", index=False)
    cldn.to_csv(TABLES / "claudin_panel_vs_tnk.tsv", sep="\t", index=False)
    pd.DataFrame([{"set": k, "n": v} for k, v in counts.items()]).to_csv(
        TABLES / "intersection_counts.tsv", sep="\t", index=False
    )

    summary = {
        "counts": counts,
        "up_list_low_tnk": int(len(up_low)),
        "up_list_high_tnk": int(len(up_high)),
        "up_most_neg_gene": str(most_neg["gene"]),
        "up_most_neg_rho": float(most_neg["rho"]),
        "up_most_neg_p": float(most_neg["p"]),
        "up_pos_best_fdr": float(up_high["fdr_screen"].min()) if len(up_high) else float("nan"),
        "tj_low_best_fdr": float(tj_low["fdr_tj"].min()) if len(tj_low) else float("nan"),
        "cldn23_logfc": float(c23["logFC"]),
        "cldn23_p": float(c23["de_p"]),
        "cldn23_rho": float(c23["rho"]),
        "cldn23_meta_p": float(c23["p"]),
        "cldn4_logfc": float(c4["logFC"]),
        "cldn4_de_p": float(c4["p"]),
        "cldn4_rho": float(c4_hit["rho"]),
        "cldn4_meta_p": float(c4_hit["p"]),
        "cldn4_I2": float(c4_hit["I2"]),
        "cldn4_n_neg": int(c4_hit["n_neg"]),
        "cldn4_n_cohorts": int(c4_hit["n_cohorts"]),
        "cldn4_pct_rho": float(pct["rho"]),
        "cldn4_pct_p": float(pct["p"]),
        "cldn4_pct_I2": float(pct["I2"]),
        "cldn4_pct_n_neg": int(pct["n_neg"]),
        "cldn1_de_p": float(c1["de_p"]),
        "cldn1_rho": float(c1["rho"]),
        "cldn1_n_neg": int(c1["n_neg"]),
        "cldn1_n_cohorts": int(c1["n_cohorts"]),
        "de_n_q1": int((m_de["quartile"] == "Q1").sum()),
        "de_n_q4": int((m_de["quartile"] == "Q4").sum()),
        "expression_n": int(meta["patient"].nunique()),
        "locked_units_n": int(len(units)),
        "missing_count_unit": missing,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_finding(summary, up_tj, nominal, tj_low, cldn)
    forest_plot(up_tj, c4_hit, pct)
    print(json.dumps(counts, indent=2))
    print("three-way rows", len(three))
    print("wrote", HERE / "FINDING.md")


def gene_vs_tnk_from_column(meta: pd.DataFrame, col: str) -> dict:
    rhos, ns, members = [], [], []
    for c in COHORTS:
        m = meta.loc[meta["cohort"] == c]
        rho, p, n = spearman(m[col], m["frac_tnk"])
        if n >= 4 and math.isfinite(rho):
            rhos.append(rho)
            ns.append(n)
            members.append(c)
    dl = random_effects_dl(rhos, ns)
    return {
        "rho": dl["pooled_rho"],
        "p": dl["p"],
        "I2": dl["I2"],
        "ci_lo": dl["ci95_rho"][0],
        "ci_hi": dl["ci95_rho"][1],
        "k": dl["k"],
        "n": dl["n_patients_total"],
        "n_neg": int(sum(1 for r in rhos if r < 0)),
        "n_cohorts": len(rhos),
        "member_rhos": ",".join(f"{c}:{r:.3f}" for c, r in zip(members, rhos)),
        "predicts_low_tnk": bool(dl["pooled_rho"] < 0 and dl["p"] < 0.05),
    }


if __name__ == "__main__":
    main()

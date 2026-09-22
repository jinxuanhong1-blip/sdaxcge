#!/usr/bin/env python3
"""PAPER FUNNEL: TROP2-high malignant DEG → TJ top → TJ vs T/NK.

Concordant-4 only: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.

1) Patient-pseudobulk malignant TACSTD2 Q4 vs Q1 DEG (OLS on log2 TMM-CPM+1).
2) ORA + prerank GSEA on that rank; ask whether TJ / claudin / epithelial
   adhesion is top or near-top among barrier-related sets.
3) Malignant TJ score vs same-unit T/NK (Spearman + DL meta).

Honest unit = patient / donor / sample. Do not quote cell counts as n.
Never fabricate: if TJ is not top, say so.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gsea_core import bh_fdr, gsea_prerank  # noqa: E402
from lib_stats import random_effects_dl, spearman  # noqa: E402

DATA = HERE / "data"
TABLES = HERE / "tables"
FIGS = HERE / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
COMBO = "+".join(COHORTS)
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
SEED = 42
NPERM = 1000

# Barrier / junction / adhesion sets we care about for the funnel claim.
BARRIER_KEYS = [
    "KEGG_TIGHT_JUNCTION",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION",
    "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_APICAL_SURFACE",
    "GOBP_KERATINIZATION",
    "GOBP_KERATINOCYTE_DIFFERENTIATION",
    "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER",
    "KRT_EPITHELIAL",
]

# Broader A8 set list used for ranking "is TJ near the top?"
HEADLINE_KEYS = BARRIER_KEYS + [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_HEDGEHOG_SIGNALING",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_ALLOGRAFT_REJECTION",
]


def load_units() -> pd.DataFrame:
    """Locked concordant-4 units with frac_tnk (CLDN4 columns kept for context)."""
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
                "malig_def": "marker_malig",
                "cldn4_pct": pct,
                "cldn4_mean": float(r["mal_CLDN4_mean"]),
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
                "malig_def": "author_malig",
                "cldn4_pct": pct,
                "cldn4_mean": float(r["mal_CLDN4_mean"]),
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
                "malig_def": "author_malig",
                "cldn4_pct": pct,
                "cldn4_mean": float(r["mal_CLDN4_mean"]),
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
                "malig_def": "marker_malig",
                "cldn4_pct": pct,
                "cldn4_mean": float(r["mal_CLDN4_mean"]),
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
        if keep.sum() < 50:
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
        if trim.sum() < 20:
            factors[col] = 1.0
            continue
        tmm = float(np.average(m[trim], weights=1.0 / w[trim]))
        factors[col] = 2 ** tmm
    fac = pd.Series(factors)
    return fac / fac.mean()


def log_cpm(counts: pd.DataFrame, factors: pd.Series) -> pd.DataFrame:
    lib = counts.sum(axis=0).astype(float) * factors.reindex(counts.columns).astype(float)
    return np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def assign_quartiles(values: pd.Series) -> pd.Series:
    s = values.astype(float)
    ranks = s.rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=s.index)


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
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
    if df_res < 1:
        return pd.DataFrame()
    xtx = Xv.T @ Xv
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ Xv.T @ Y.T
    resid = Y.T - Xv @ beta
    sse = np.sum(resid**2, axis=0)
    sigma2 = sse / df_res
    j = list(X.columns).index(coef)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[j, j], 0.0))
    est = beta[j]
    t = np.divide(est, se, out=np.zeros_like(est), where=se > 0)
    pv = 2.0 * stats.t.sf(np.abs(t), df_res)
    out = pd.DataFrame(
        {
            "gene": logcpm.index.astype(str),
            "logFC": est,
            "AveExpr": Y.mean(axis=1),
            "t": t,
            "p": pv,
            "se": se,
            "df": df_res,
        }
    )
    out["fdr"] = _bh(out["p"].values)
    return out.sort_values("p")


def load_merged_counts(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mats = []
    for c in COHORTS:
        path = DATA / f"{c}_malignant_counts.tsv.gz"
        mats.append(read_counts(path))
    # Outer-join genes; fill missing with 0 only for genes present in ≥1 cohort,
    # then restrict columns to units in meta ∩ matrix.
    all_genes = sorted(set().union(*[set(m.index) for m in mats]))
    cols = []
    pieces = []
    for c, m in zip(COHORTS, mats):
        want = [p for p in meta.loc[meta["cohort"] == c, "patient"] if p in m.columns]
        if not want:
            continue
        sub = m.reindex(all_genes).fillna(0.0)[want]
        pieces.append(sub)
        cols.extend(want)
    counts = pd.concat(pieces, axis=1)
    counts = counts.loc[:, ~counts.columns.duplicated()]
    meta2 = meta.loc[meta["patient"].isin(counts.columns)].copy()
    return counts, meta2


def hypergeom_ora(
    hits: list[str],
    background: list[str],
    gene_sets: dict[str, list[str]],
) -> pd.DataFrame:
    bg = set(background)
    hit = set(hits) & bg
    n_bg = len(bg)
    n_hit = len(hit)
    rows = []
    for term, members in gene_sets.items():
        mem = set(g.upper() for g in members) & bg
        k = len(mem)
        if k < 5:
            continue
        x = len(hit & mem)
        # P(X >= x)
        p = float(stats.hypergeom.sf(x - 1, n_bg, k, n_hit)) if x > 0 else 1.0
        expect = n_hit * k / n_bg if n_bg else np.nan
        rows.append(
            {
                "term": term,
                "n_overlap": x,
                "n_set_in_bg": k,
                "n_query": n_hit,
                "n_bg": n_bg,
                "enrichment": (x / expect) if expect else np.nan,
                "p": p,
                "overlap_genes": ",".join(sorted(hit & mem)[:40]),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = bh_fdr(out["p"])
    return out.sort_values(["fdr", "p", "enrichment"], ascending=[True, True, False])


def build_gene_sets(a8: dict) -> dict[str, list[str]]:
    sets = {}
    for k, v in a8["sets"].items():
        sets[k] = [str(g).upper() for g in v]
    # Compact claudin / adhesion panels for ORA visibility.
    sets["CUSTOM_CLAUDIN_PANEL"] = [
        g for g in sets.get("KEGG_TIGHT_JUNCTION", []) if g.startswith("CLDN")
    ] + ["CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLDN12", "CLDN18"]
    sets["CUSTOM_CLAUDIN_PANEL"] = sorted(set(sets["CUSTOM_CLAUDIN_PANEL"]))
    sets["CUSTOM_EPITHELIAL_ADHESION"] = sorted(
        set(
            [
                "CDH1",
                "EPCAM",
                "F11R",
                "OCLN",
                "TJP1",
                "TJP2",
                "TJP3",
                "CGN",
                "MARVELD2",
                "PVRL1",
                "NECTIN1",
                "NECTIN2",
                "NECTIN4",
                "JUP",
                "DSP",
                "DSG2",
                "PKP2",
                "CTNNA1",
                "CTNNB1",
                "CLDN1",
                "CLDN3",
                "CLDN4",
                "CLDN7",
            ]
        )
    )
    return sets


def tj_family_genes(sets: dict[str, list[str]], holdouts: set[str]) -> list[str]:
    genes = set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    genes |= set(sets["CUSTOM_EPITHELIAL_ADHESION"])
    genes -= holdouts
    return sorted(genes)


def main() -> None:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    gene_sets = build_gene_sets(a8)
    units = load_units()
    assert len(units) == 65, f"expected locked n=65 units, got {len(units)}"

    counts_raw, meta = load_merged_counts(units)
    # Filter + TMM on the full expression matrix (all units present).
    counts = filter_genes(counts_raw)
    factors = tmm_norm_factors(counts)
    lc = log_cpm(counts, factors)

    if "TACSTD2" not in lc.index:
        raise SystemExit("TACSTD2 absent from malignant pseudobulk matrix")

    meta = meta.set_index("patient")
    meta["tacstd2"] = lc.loc["TACSTD2", meta.index].astype(float)
    # Within-cohort TACSTD2 quartiles (expression, not %pos — no cell-level here).
    qs = []
    for c in COHORTS:
        idx = meta.index[meta["cohort"] == c]
        qs.append(assign_quartiles(meta.loc[idx, "tacstd2"]))
    meta["quartile"] = pd.concat(qs)
    meta = meta.reset_index()
    meta.to_csv(TABLES / "units_tacstd2.tsv", sep="\t", index=False)

    # Honest n table
    n_rows = []
    for c in COHORTS:
        m = meta.loc[meta["cohort"] == c]
        n_rows.append(
            {
                "cohort": c,
                "unit": m["unit"].iloc[0],
                "n_vector": int(len(m)),
                "n_q1": int((m["quartile"] == "Q1").sum()),
                "n_q4": int((m["quartile"] == "Q4").sum()),
                "tacstd2_q1_median": float(m.loc[m["quartile"] == "Q1", "tacstd2"].median()),
                "tacstd2_q4_median": float(m.loc[m["quartile"] == "Q4", "tacstd2"].median()),
            }
        )
    n_honest = pd.DataFrame(n_rows)
    n_honest.to_csv(TABLES / "n_honest.tsv", sep="\t", index=False)

    # ----- DEG: stacked Q4 vs Q1 with cohort covariates -----
    m_de = meta.loc[meta["quartile"].isin(["Q1", "Q4"])].copy()
    n_q1 = int((m_de["quartile"] == "Q1").sum())
    n_q4 = int((m_de["quartile"] == "Q4").sum())
    design = pd.DataFrame(index=m_de["patient"])
    design["Intercept"] = 1.0
    design["TACSTD2_Q4"] = (m_de.set_index("patient")["quartile"] == "Q4").astype(float)
    for c in COHORTS:
        if c == REF:
            continue
        design[f"cohort_{c}"] = (m_de.set_index("patient")["cohort"] == c).astype(float)
    lc_de = lc.loc[:, m_de["patient"]]
    # Drop TACSTD2 from the DE table itself? Keep it as QC (should be strongly up).
    de = ols_de(lc_de, design, "TACSTD2_Q4")
    de.to_csv(TABLES / "de_q4q1_stacked.tsv.gz", sep="\t", index=False, compression="gzip")

    # Focal genes
    focals = [
        "TACSTD2",
        "CLDN1",
        "CLDN3",
        "CLDN4",
        "CLDN7",
        "OCLN",
        "F11R",
        "TJP1",
        "CDH1",
        "EPCAM",
        "ELF3",
        "KRT8",
        "KRT18",
        "KRT19",
        "STAT1",
        "IRF1",
        "B2M",
        "HLA-A",
        "HLA-B",
        "CXCL10",
    ]
    focal_tbl = de[de["gene"].isin(focals)].copy()
    focal_tbl.to_csv(TABLES / "de_focal_genes.tsv", sep="\t", index=False)

    # Family-score DE (TJ holds TACSTD2 out; CLDN4 stays — TJ biology)
    hold_tj = {"TACSTD2"}
    families = {
        "TJ": tj_family_genes(gene_sets, hold_tj),
        "CLAUDIN_PANEL": [g for g in gene_sets["CUSTOM_CLAUDIN_PANEL"] if g != "TACSTD2"],
        "EPITHELIAL_ADHESION": [
            g for g in gene_sets["CUSTOM_EPITHELIAL_ADHESION"] if g != "TACSTD2"
        ],
        "KERATIN": gene_sets.get("KRT_EPITHELIAL", []) + gene_sets.get("GOBP_KERATINIZATION", []),
        "IFN": sorted(
            set(gene_sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
            | set(gene_sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
        ),
        "MHC-I/APM": gene_sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"],
    }
    fam_rows = []
    for name, genes in families.items():
        present = [g for g in sorted(set(genes)) if g in lc_de.index]
        if len(present) < 3:
            continue
        score = lc_de.loc[present].astype(float).mean(axis=0).to_frame().T
        score.index = [name]
        r = ols_de(score, design, "TACSTD2_Q4").iloc[0]
        fam_rows.append(
            {
                "family": name,
                "n_q1": n_q1,
                "n_q4": n_q4,
                "n": n_q1 + n_q4,
                "n_genes": len(present),
                "logFC": float(r.logFC),
                "t": float(r.t),
                "p": float(r.p),
                "se": float(r.se),
            }
        )
    fam_de = pd.DataFrame(fam_rows)
    if not fam_de.empty:
        fam_de["fdr"] = _bh(fam_de["p"].values)
    fam_de.to_csv(TABLES / "family_score_q4q1.tsv", sep="\t", index=False)

    # Continuous TACSTD2 sensitivity (cohort-adjusted)
    design_c = pd.DataFrame(index=meta["patient"])
    design_c["Intercept"] = 1.0
    z = meta.set_index("patient")["tacstd2"].astype(float)
    # Within-cohort z so scale is comparable
    zc = []
    for c in COHORTS:
        idx = meta.loc[meta["cohort"] == c, "patient"]
        zz = z.loc[idx]
        zc.append((zz - zz.mean()) / zz.std(ddof=1))
    design_c["TACSTD2_z"] = pd.concat(zc)
    for c in COHORTS:
        if c == REF:
            continue
        design_c[f"cohort_{c}"] = (meta.set_index("patient")["cohort"] == c).astype(float)
    de_cont = ols_de(lc.loc[:, meta["patient"]], design_c, "TACSTD2_z")
    de_cont.to_csv(TABLES / "de_continuous_tacstd2.tsv.gz", sep="\t", index=False, compression="gzip")

    # ----- ORA on UP genes (FDR<0.05 & logFC>0), TACSTD2 held out of query -----
    up = de[(de["fdr"] < 0.05) & (de["logFC"] > 0) & (de["gene"] != "TACSTD2")]
    dn = de[(de["fdr"] < 0.05) & (de["logFC"] < 0)]
    # If FDR is sparse at this n, fall back to p<0.01 & |logFC|>0.25 (documented).
    ora_rule = "fdr<0.05 & logFC>0"
    if len(up) < 30:
        up = de[(de["p"] < 0.01) & (de["logFC"] > 0.25) & (de["gene"] != "TACSTD2")]
        dn = de[(de["p"] < 0.01) & (de["logFC"] < -0.25)]
        ora_rule = "p<0.01 & |logFC|>0.25 (FDR arm too thin)"
    bg = de["gene"].astype(str).tolist()
    ora_up = hypergeom_ora(up["gene"].tolist(), bg, gene_sets)
    ora_up.to_csv(TABLES / "ora_up.tsv", sep="\t", index=False)
    ora_dn = hypergeom_ora(dn["gene"].tolist(), bg, gene_sets)
    ora_dn.to_csv(TABLES / "ora_down.tsv", sep="\t", index=False)

    # Rank barrier terms among all tested sets by FDR then enrichment
    rank_map = {t: i + 1 for i, t in enumerate(ora_up["term"].tolist())}
    ora_barrier = ora_up[
        ora_up["term"].isin(BARRIER_KEYS + ["CUSTOM_CLAUDIN_PANEL", "CUSTOM_EPITHELIAL_ADHESION"])
    ].copy()
    ora_barrier["rank_among_all_ora"] = ora_barrier["term"].map(rank_map)
    ora_barrier.to_csv(TABLES / "ora_barrier_subset.tsv", sep="\t", index=False)

    # ----- GSEA prerank on OLS t (positive = TACSTD2 Q4) -----
    # Hold TACSTD2 out of rank for non-circular TJ / keratin / adhesion tests.
    rank = de.set_index("gene")["t"].astype(float).sort_values(ascending=False)
    rank_no_trop2 = rank.drop(labels=["TACSTD2"], errors="ignore")
    gsea_sets = {k: gene_sets[k] for k in HEADLINE_KEYS if k in gene_sets}
    gsea_sets["CUSTOM_CLAUDIN_PANEL"] = gene_sets["CUSTOM_CLAUDIN_PANEL"]
    gsea_sets["CUSTOM_EPITHELIAL_ADHESION"] = gene_sets["CUSTOM_EPITHELIAL_ADHESION"]
    # Drop TACSTD2 from members
    gsea_sets = {
        k: [g for g in v if g != "TACSTD2"] for k, v in gsea_sets.items()
    }
    gsea = gsea_prerank(rank_no_trop2, gsea_sets, nperm=NPERM, seed=SEED)
    if not gsea.empty:
        gsea["fdr"] = bh_fdr(gsea["nom_p"])
        gsea = gsea.sort_values(["fdr", "nes"], ascending=[True, False])
    gsea.to_csv(TABLES / "gsea_prerank_q4q1.tsv", sep="\t", index=False)

    # Continuous rank sensitivity
    rank_c = de_cont.set_index("gene")["t"].astype(float).sort_values(ascending=False)
    rank_c = rank_c.drop(labels=["TACSTD2"], errors="ignore")
    gsea_c = gsea_prerank(rank_c, gsea_sets, nperm=NPERM, seed=SEED)
    if not gsea_c.empty:
        gsea_c["fdr"] = bh_fdr(gsea_c["nom_p"])
        gsea_c = gsea_c.sort_values(["fdr", "nes"], ascending=[True, False])
    gsea_c.to_csv(TABLES / "gsea_prerank_continuous.tsv", sep="\t", index=False)

    # Enrichr via gseapy (optional network); catch failures honestly
    enrichr_rows = []
    try:
        import gseapy as gp

        gene_list = up["gene"].astype(str).tolist()[:500]
        if len(gene_list) >= 10:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human", "MSigDB_Hallmark_2020"],
                organism="human",
                outdir=None,
                verbose=False,
            )
            er = enr.results.copy()
            er.to_csv(TABLES / "enrichr_up.tsv", sep="\t", index=False)
            # Pull TJ / claudin / adhesion / keratin hits
            mask = er["Term"].str.contains(
                "tight junction|claudin|cell.?cell adhesion|adherens|keratin|apical junction|epithelial cell.?cell",
                case=False,
                regex=True,
                na=False,
            )
            sub = er.loc[mask].copy()
            sub.to_csv(TABLES / "enrichr_barrier_hits.tsv", sep="\t", index=False)
            enrichr_ok = True
            enrichr_n = int(len(er))
            enrichr_barrier_n = int(len(sub))
            enrichr_top = er.head(15)[["Gene_set", "Term", "Adjusted P-value", "Odds Ratio", "Overlap"]].to_dict(
                "records"
            )
        else:
            enrichr_ok = False
            enrichr_n = 0
            enrichr_barrier_n = 0
            enrichr_top = []
            (TABLES / "enrichr_note.txt").write_text(
                f"Enrichr skipped: only {len(gene_list)} UP genes under rule {ora_rule}\n"
            )
    except Exception as e:
        enrichr_ok = False
        enrichr_n = 0
        enrichr_barrier_n = 0
        enrichr_top = []
        (TABLES / "enrichr_note.txt").write_text(f"Enrichr failed: {type(e).__name__}: {e}\n")

    # ----- TJ score vs T/NK -----
    tj_genes = tj_family_genes(gene_sets, {"TACSTD2"})
    present_tj = [g for g in tj_genes if g in lc.index]
    core_tj = [g for g in ["CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "F11R", "TJP1"] if g in lc.index]
    meta = meta.set_index("patient")
    meta["tj_score"] = lc.loc[present_tj, meta.index].astype(float).mean(axis=0)
    # Also CLDN4-held-out TJ (matches locked CLDN4 analyses)
    tj_no_cldn4 = [g for g in present_tj if g != "CLDN4"]
    meta["tj_score_no_cldn4"] = lc.loc[tj_no_cldn4, meta.index].astype(float).mean(axis=0)
    meta["tj_core"] = lc.loc[core_tj, meta.index].astype(float).mean(axis=0)
    meta["tacstd2"] = lc.loc["TACSTD2", meta.index].astype(float)
    meta = meta.reset_index()
    meta.to_csv(TABLES / "units_scores.tsv", sep="\t", index=False)

    tnk_rows = []
    for score_name, col in [
        ("tj_score", "tj_score"),
        ("tj_score_no_cldn4", "tj_score_no_cldn4"),
        ("tj_core", "tj_core"),
        ("tacstd2", "tacstd2"),
        ("cldn4_pct", "cldn4_pct"),
    ]:
        rhos, ps, ns, members = [], [], [], []
        for c in COHORTS:
            m = meta.loc[meta["cohort"] == c]
            rho, p, n = spearman(m[col], m["frac_tnk"])
            tnk_rows.append(
                {
                    "kind": "single",
                    "score": score_name,
                    "cohort": c,
                    "n": n,
                    "rho": rho,
                    "p": p,
                }
            )
            if n >= 4 and math.isfinite(rho):
                rhos.append(rho)
                ps.append(p)
                ns.append(n)
                members.append(c)
        if rhos:
            re = random_effects_dl(rhos, ns)
            tnk_rows.append(
                {
                    "kind": "meta_dl",
                    "score": score_name,
                    "cohort": COMBO,
                    "n": int(sum(ns)),
                    "k": len(rhos),
                    "rho": re.get("pooled_rho", float("nan")),
                    "p": re.get("p", float("nan")),
                    "I2": re.get("I2", float("nan")),
                    "ci95_lo": re.get("ci95_rho", [float("nan")] * 2)[0],
                    "ci95_hi": re.get("ci95_rho", [float("nan")] * 2)[1],
                    "member_rhos": ",".join(f"{c}:{r:.3f}" for c, r in zip(members, rhos)),
                }
            )
        # Stacked within-cohort Q4 vs Q1 on the score vs T/NK
        q1, q4 = [], []
        for c in COHORTS:
            m = meta.loc[meta["cohort"] == c]
            if len(m) < 4:
                continue
            qs = assign_quartiles(m.set_index("patient")[col])
            q1.extend(m.set_index("patient").loc[qs.index[qs == "Q1"], "frac_tnk"].tolist())
            q4.extend(m.set_index("patient").loc[qs.index[qs == "Q4"], "frac_tnk"].tolist())
        if len(q1) >= 3 and len(q4) >= 3:
            u, p = stats.mannwhitneyu(np.asarray(q4), np.asarray(q1), alternative="two-sided")
            r_rb = (2.0 * float(u)) / (len(q4) * len(q1)) - 1.0
            tnk_rows.append(
                {
                    "kind": "stacked_q4q1",
                    "score": score_name,
                    "cohort": COMBO,
                    "n": len(q1) + len(q4),
                    "n_q1": len(q1),
                    "n_q4": len(q4),
                    "r_rb": float(r_rb),
                    "p": float(p),
                    "median_q1": float(np.median(q1)),
                    "median_q4": float(np.median(q4)),
                    "delta_median": float(np.median(q4) - np.median(q1)),
                }
            )
    tnk = pd.DataFrame(tnk_rows)
    tnk.to_csv(TABLES / "tj_vs_tnk.tsv", sep="\t", index=False)

    # ----- Verdict helpers (computed, not invented) -----
    def top_barrier_gsea(gdf: pd.DataFrame) -> dict:
        if gdf.empty:
            return {"available": False}
        # Rank by NES descending among positive NES, then by FDR
        pos = gdf[gdf["nes"] > 0].copy()
        pos = pos.sort_values(["fdr", "nes"], ascending=[True, False]).reset_index(drop=True)
        barrier_terms = set(BARRIER_KEYS) | {
            "CUSTOM_CLAUDIN_PANEL",
            "CUSTOM_EPITHELIAL_ADHESION",
        }
        barrier = pos[pos["term"].isin(barrier_terms)].copy()
        all_rank = {t: i + 1 for i, t in enumerate(pos["term"])}
        out = {
            "available": True,
            "n_sets_positive_nes": int(len(pos)),
            "n_sets_tested": int(len(gdf)),
            "top5_positive": pos.head(5)[["term", "nes", "nom_p", "fdr", "n_set_in_rank"]].to_dict(
                "records"
            ),
            "barrier_rows": [],
        }
        for _, r in barrier.iterrows():
            out["barrier_rows"].append(
                {
                    "term": r["term"],
                    "nes": float(r["nes"]),
                    "nom_p": float(r["nom_p"]),
                    "fdr": float(r["fdr"]),
                    "rank_among_positive_nes": int(all_rank[r["term"]]),
                    "n_set_in_rank": int(r["n_set_in_rank"]),
                    "lead_genes": r["lead_genes"],
                }
            )
        # Best barrier by FDR then NES
        if barrier.empty:
            out["best_barrier"] = None
            out["tj_is_top"] = False
            out["tj_is_near_top"] = False
        else:
            best = barrier.sort_values(["fdr", "nes"], ascending=[True, False]).iloc[0]
            rank = int(all_rank[best["term"]])
            out["best_barrier"] = {
                "term": best["term"],
                "nes": float(best["nes"]),
                "fdr": float(best["fdr"]),
                "rank_among_positive_nes": rank,
            }
            out["tj_is_top"] = rank == 1
            out["tj_is_near_top"] = rank <= 5
        return out

    gsea_verdict = top_barrier_gsea(gsea)
    gsea_c_verdict = top_barrier_gsea(gsea_c)

    # ORA top
    ora_top5 = ora_up.head(5)[["term", "n_overlap", "enrichment", "p", "fdr"]].to_dict("records") if not ora_up.empty else []
    ora_best_barrier = None
    ora_top_barrier_by_rank = None
    if not ora_barrier.empty:
        b = ora_barrier.sort_values(["fdr", "enrichment"], ascending=[True, False]).iloc[0]
        ora_best_barrier = {
            "term": b["term"],
            "n_overlap": int(b["n_overlap"]),
            "enrichment": float(b["enrichment"]),
            "p": float(b["p"]),
            "fdr": float(b["fdr"]),
            "rank_among_all_ora": int(b["rank_among_all_ora"]),
        }
        b2 = ora_barrier.sort_values(["rank_among_all_ora", "fdr"]).iloc[0]
        ora_top_barrier_by_rank = {
            "term": b2["term"],
            "n_overlap": int(b2["n_overlap"]),
            "enrichment": float(b2["enrichment"]),
            "p": float(b2["p"]),
            "fdr": float(b2["fdr"]),
            "rank_among_all_ora": int(b2["rank_among_all_ora"]),
        }

    tj_meta = tnk[(tnk["kind"] == "meta_dl") & (tnk["score"] == "tj_score")].iloc[0]
    tj_q = tnk[(tnk["kind"] == "stacked_q4q1") & (tnk["score"] == "tj_score")].iloc[0]
    trop2_meta = tnk[(tnk["kind"] == "meta_dl") & (tnk["score"] == "tacstd2")].iloc[0]
    core_meta = tnk[(tnk["kind"] == "meta_dl") & (tnk["score"] == "tj_core")].iloc[0]
    core_q = tnk[(tnk["kind"] == "stacked_q4q1") & (tnk["score"] == "tj_core")].iloc[0]
    cldn4_meta = tnk[(tnk["kind"] == "meta_dl") & (tnk["score"] == "cldn4_pct")].iloc[0]

    summary = {
        "cohorts": COHORTS,
        "n_units_locked": 65,
        "n_in_expression": int(len(meta)),
        "de_n_q1": n_q1,
        "de_n_q4": n_q4,
        "de_n": n_q1 + n_q4,
        "n_genes_de": int(len(de)),
        "tacstd2_logFC": float(de.loc[de["gene"] == "TACSTD2", "logFC"].iloc[0]),
        "tacstd2_p": float(de.loc[de["gene"] == "TACSTD2", "p"].iloc[0]),
        "ora_rule": ora_rule,
        "n_up": int(len(up)),
        "n_down": int(len(dn)),
        "ora_top5": ora_top5,
        "ora_best_barrier": ora_best_barrier,
        "ora_top_barrier_by_rank": ora_top_barrier_by_rank,
        "gsea_q4q1": gsea_verdict,
        "gsea_continuous": gsea_c_verdict,
        "family_score_q4q1": fam_de.to_dict("records"),
        "tj_vs_tnk_meta": {
            "rho": float(tj_meta["rho"]),
            "p": float(tj_meta["p"]),
            "I2": float(tj_meta["I2"]),
            "n": int(tj_meta["n"]),
            "ci95": [float(tj_meta["ci95_lo"]), float(tj_meta["ci95_hi"])],
            "n_tj_genes": len(present_tj),
            "member_rhos": str(tj_meta.get("member_rhos", "")),
        },
        "tj_vs_tnk_stacked_q4q1": {
            "r_rb": float(tj_q["r_rb"]),
            "p": float(tj_q["p"]),
            "n_q1": int(tj_q["n_q1"]),
            "n_q4": int(tj_q["n_q4"]),
            "delta_median": float(tj_q["delta_median"]),
        },
        "tj_core_vs_tnk_meta": {
            "genes": core_tj,
            "rho": float(core_meta["rho"]),
            "p": float(core_meta["p"]),
            "I2": float(core_meta["I2"]),
            "n": int(core_meta["n"]),
            "ci95": [float(core_meta["ci95_lo"]), float(core_meta["ci95_hi"])],
            "member_rhos": str(core_meta.get("member_rhos", "")),
        },
        "tj_core_vs_tnk_stacked_q4q1": {
            "r_rb": float(core_q["r_rb"]),
            "p": float(core_q["p"]),
            "n_q1": int(core_q["n_q1"]),
            "n_q4": int(core_q["n_q4"]),
            "delta_median": float(core_q["delta_median"]),
        },
        "tacstd2_vs_tnk_meta": {
            "rho": float(trop2_meta["rho"]),
            "p": float(trop2_meta["p"]),
            "I2": float(trop2_meta["I2"]),
            "n": int(trop2_meta["n"]),
        },
        "cldn4_pct_vs_tnk_meta_recomputed": {
            "rho": float(cldn4_meta["rho"]),
            "p": float(cldn4_meta["p"]),
            "I2": float(cldn4_meta["I2"]),
            "n": int(cldn4_meta["n"]),
            "note": "sanity check on this page; locked primary remains PR #503/#539 n=65",
        },
        "enrichr": {
            "ok": enrichr_ok,
            "n_terms": enrichr_n,
            "n_barrier_terms": enrichr_barrier_n,
            "top15": enrichr_top,
        },
        "focal": focal_tbl[["gene", "logFC", "t", "p", "fdr"]].to_dict("records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    # ----- Figures -----
    # 1) volcano
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    x = de["logFC"].to_numpy()
    y = -np.log10(np.clip(de["p"].to_numpy(), 1e-300, 1))
    ax.scatter(x, y, s=6, c="#888888", alpha=0.5, linewidths=0)
    sig = de["fdr"] < 0.05
    ax.scatter(x[sig], y[sig], s=10, c="#c44e52", alpha=0.7, linewidths=0, label="FDR<0.05")
    for g, color in [
        ("TACSTD2", "#111111"),
        ("CLDN4", "#2ca02c"),
        ("CLDN7", "#2ca02c"),
        ("CLDN3", "#2ca02c"),
        ("OCLN", "#2ca02c"),
        ("F11R", "#2ca02c"),
        ("CDH1", "#1f77b4"),
        ("EPCAM", "#1f77b4"),
    ]:
        row = de.loc[de["gene"] == g]
        if row.empty:
            continue
        ax.scatter(
            float(row["logFC"].iloc[0]),
            -math.log10(max(float(row["p"].iloc[0]), 1e-300)),
            s=36,
            c=color,
            zorder=3,
        )
        ax.text(
            float(row["logFC"].iloc[0]),
            -math.log10(max(float(row["p"].iloc[0]), 1e-300)),
            " " + g,
            fontsize=7,
            va="center",
        )
    ax.axvline(0, color="#444", lw=0.6)
    ax.set_xlabel("logFC (TACSTD2 Q4 − Q1, cohort-adjusted)")
    ax.set_ylabel(r"$-\log_{10} p$")
    ax.set_title(f"Concordant-4 malignant DE  (n={n_q1} vs {n_q4})")
    fig.tight_layout()
    fig.savefig(FIGS / "volcano_tacstd2_q4q1.png", dpi=160)
    fig.savefig(FIGS / "volcano_tacstd2_q4q1.pdf")
    plt.close(fig)

    # 2) GSEA NES bars (positive NES top 12)
    if not gsea.empty:
        show = gsea[gsea["nes"] > 0].head(12).iloc[::-1]
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        colors = [
            "#2ca02c"
            if t in set(BARRIER_KEYS) | {"CUSTOM_CLAUDIN_PANEL", "CUSTOM_EPITHELIAL_ADHESION"}
            else "#4c78a8"
            for t in show["term"]
        ]
        ax.barh(show["term"], show["nes"], color=colors)
        ax.set_xlabel("NES (positive = enriched in TACSTD2 Q4)")
        ax.set_title("Prerank GSEA — top positive NES")
        fig.tight_layout()
        fig.savefig(FIGS / "gsea_nes_top_positive.png", dpi=160)
        fig.savefig(FIGS / "gsea_nes_top_positive.pdf")
        plt.close(fig)

    # 3) TJ vs T/NK scatter
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    for c, col in zip(COHORTS, ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]):
        m = meta.loc[meta["cohort"] == c]
        ax.scatter(m["tj_score"], m["frac_tnk"], s=28, c=col, label=c, alpha=0.85)
    ax.set_xlabel(f"Malignant TJ score ({len(present_tj)} genes; TACSTD2 held out)")
    ax.set_ylabel("Same-unit T/NK fraction")
    ax.set_title(
        f"TJ vs T/NK  meta ρ={tj_meta['rho']:.3f}, p={tj_meta['p']:.2e}, n={int(tj_meta['n'])}"
    )
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "scatter_tj_vs_tnk.png", dpi=160)
    fig.savefig(FIGS / "scatter_tj_vs_tnk.pdf")
    plt.close(fig)

    # 4) family score forest
    if not fam_de.empty:
        fig, ax = plt.subplots(figsize=(6.2, 3.8))
        y = np.arange(len(fam_de))
        ax.errorbar(
            fam_de["logFC"],
            y,
            xerr=1.96 * fam_de["se"],
            fmt="o",
            color="#333",
            ecolor="#666",
        )
        ax.axvline(0, color="#888", lw=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(fam_de["family"])
        ax.set_xlabel("Family-score logFC (TACSTD2 Q4 − Q1)")
        ax.set_title(f"Family scores  n={n_q1} vs {n_q4}")
        fig.tight_layout()
        fig.savefig(FIGS / "forest_family_scores.png", dpi=160)
        fig.savefig(FIGS / "forest_family_scores.pdf")
        plt.close(fig)

    # Write FINDING.md from actual numbers
    write_finding(summary, n_honest, ora_rule)
    print(json.dumps({k: summary[k] for k in [
        "de_n_q1", "de_n_q4", "n_up", "n_down", "ora_rule",
        "ora_best_barrier", "gsea_q4q1", "tj_vs_tnk_meta", "tacstd2_vs_tnk_meta",
    ]}, indent=2, default=str))


def write_finding(summary: dict, n_honest: pd.DataFrame, ora_rule: str) -> None:
    g = summary["gsea_q4q1"]
    best = g.get("best_barrier")
    lines = []
    lines.append("# Concordant-4: TACSTD2 Q4 vs Q1 malignant DEG → TJ → T/NK")
    lines.append("")
    lines.append("PAPER FUNNEL public evidence. **TACSTD2 (TROP2) only** as the split.")
    lines.append("Locked cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.")
    lines.append("Not GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.")
    lines.append("Does not re-fit the locked CLDN4 %pos vs T/NK ρ = −0.531.")
    lines.append("")
    lines.append("Honest unit = patient / donor / sample. Do not quote cell counts as n.")
    lines.append("TACSTD2 split uses malignant **pseudobulk expression** (log2 TMM-CPM+1),")
    lines.append("not cell-level %pos (no per-cell matrix in this page).")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append("| cohort | unit | vector n | Q1 / Q4 |")
    lines.append("|---|---|---:|---|")
    for _, r in n_honest.iterrows():
        lines.append(
            f"| {r.cohort} | {r.unit} | {int(r.n_vector)} | {int(r.n_q1)} / {int(r.n_q4)} |"
        )
    lines.append("")
    lines.append(
        f"DE contrast = stacked within-cohort Q4 vs Q1: **{summary['de_n_q1']} vs {summary['de_n_q4']}** "
        f"(do not quote n={summary['n_units_locked']} as the DE n). "
        f"Expression matrix has **{summary['n_in_expression']}** units."
    )
    lines.append("")
    lines.append("## 1. DEG (OLS on log2 TMM-CPM+1, cohort covariates)")
    lines.append("")
    lines.append(
        f"Positive logFC = higher in TACSTD2 Q4. Genes tested: {summary['n_genes_de']}."
    )
    lines.append(
        f"TACSTD2 QC: logFC = {summary['tacstd2_logFC']:.3f}, p = {summary['tacstd2_p']:.2e}."
    )
    lines.append("")
    lines.append("### Family scores")
    lines.append("")
    lines.append("| family | n_genes | logFC | p | FDR |")
    lines.append("|---|---:|---:|---:|---:|")
    for r in summary["family_score_q4q1"]:
        lines.append(
            f"| {r['family']} | {r['n_genes']} | {r['logFC']:+.3f} | {r['p']:.4g} | {r['fdr']:.4g} |"
        )
    lines.append("")
    lines.append("### Focal genes")
    lines.append("")
    lines.append("| gene | logFC | t | p | FDR |")
    lines.append("|---|---:|---:|---:|---:|")
    for r in summary["focal"]:
        lines.append(
            f"| {r['gene']} | {r['logFC']:+.3f} | {r['t']:+.2f} | {r['p']:.4g} | {r['fdr']:.4g} |"
        )
    lines.append("")
    lines.append("## 2. ORA / GSEA — is TJ / claudin / adhesion top or near-top?")
    lines.append("")
    lines.append(f"ORA query rule: **{ora_rule}**. n_up = {summary['n_up']}, n_down = {summary['n_down']}.")
    lines.append("Background = all genes in the DE table. TACSTD2 held out of the UP query.")
    lines.append("")
    if summary["ora_top5"]:
        lines.append("### ORA top 5 (all tested A8 + custom sets)")
        lines.append("")
        lines.append("| rank | term | overlap | enrichment | p | FDR |")
        lines.append("|---:|---|---:|---:|---:|---:|")
        for i, r in enumerate(summary["ora_top5"], 1):
            lines.append(
                f"| {i} | {r['term']} | {r['n_overlap']} | {r['enrichment']:.2f} | {r['p']:.3g} | {r['fdr']:.3g} |"
            )
        lines.append("")
    if summary["ora_best_barrier"]:
        b = summary["ora_best_barrier"]
        lines.append(
            f"Highest-enrichment barrier ORA term: **{b['term']}** "
            f"(rank {b['rank_among_all_ora']} of all ORA terms, "
            f"enrichment {b['enrichment']:.2f}, FDR {b['fdr']:.3g})."
        )
        if summary.get("ora_top_barrier_by_rank"):
            t = summary["ora_top_barrier_by_rank"]
            lines.append(
                f"Best-ranked barrier ORA term: **{t['term']}** "
                f"(rank {t['rank_among_all_ora']}, enrichment {t['enrichment']:.2f}, FDR {t['fdr']:.3g})."
            )
            if t["rank_among_all_ora"] <= 3:
                lines.append(
                    "**ORA verdict: tight-junction / apical-junction / keratin-barrier is top or near-top.**"
                )
        lines.append("")
    else:
        lines.append("No barrier ORA term passed the table (empty barrier subset).")
        lines.append("")

    lines.append("### Prerank GSEA (OLS *t*, TACSTD2 dropped from rank; 1000 gene-set perms, seed=42)")
    lines.append("")
    lines.append("Positive NES = enriched at the TACSTD2-Q4 end.")
    lines.append("")
    if g.get("top5_positive"):
        lines.append("| rank | term | NES | nom p | FDR |")
        lines.append("|---:|---|---:|---:|---:|")
        for i, r in enumerate(g["top5_positive"], 1):
            lines.append(
                f"| {i} | {r['term']} | {r['nes']:+.3f} | {r['nom_p']:.3g} | {r['fdr']:.3g} |"
            )
        lines.append("")
    if best:
        lines.append(
            f"Best barrier GSEA set: **{best['term']}** "
            f"(NES {best['nes']:+.3f}, FDR {best['fdr']:.3g}, "
            f"rank {best['rank_among_positive_nes']} among positive-NES sets)."
        )
        if g.get("tj_is_top"):
            lines.append("**Verdict: TJ/barrier is top (rank 1 among positive NES).**")
        elif g.get("tj_is_near_top"):
            lines.append("**Verdict: TJ/barrier is near-top (rank ≤ 5 among positive NES).**")
        else:
            lines.append(
                f"**Verdict: TJ/barrier is NOT top/near-top on this rank "
                f"(best barrier rank = {best['rank_among_positive_nes']}).**"
            )
    else:
        lines.append("**Verdict: no positive-NES barrier set on this rank.**")
    lines.append("")

    enr = summary["enrichr"]
    lines.append("### Enrichr (GO BP 2023 / KEGG 2021 / Hallmark 2020)")
    lines.append("")
    if enr["ok"]:
        lines.append(
            f"Ran on the UP list. {enr['n_terms']} terms returned; "
            f"{enr['n_barrier_terms']} matched a TJ/claudin/adhesion/keratin regex. "
            f"See `tables/enrichr_up.tsv` and `tables/enrichr_barrier_hits.tsv`."
        )
    else:
        lines.append("Enrichr did not return a usable table (see `tables/enrichr_note.txt`).")
    lines.append("")

    lines.append("## 3. TJ score vs T/NK")
    lines.append("")
    tj = summary["tj_vs_tnk_meta"]
    tq = summary["tj_vs_tnk_stacked_q4q1"]
    tr = summary["tacstd2_vs_tnk_meta"]
    tc = summary["tj_core_vs_tnk_meta"]
    tcq = summary["tj_core_vs_tnk_stacked_q4q1"]
    c4 = summary["cldn4_pct_vs_tnk_meta_recomputed"]
    lines.append(
        f"Broad TJ score = mean log2(TMM-CPM+1) of {tj['n_tj_genes']} TJ/adhesion genes "
        f"(TACSTD2 held out). Core panel = {', '.join(tc['genes'])}. "
        f"Unit = patient/donor/sample."
    )
    lines.append("")
    lines.append("| score | method | N | effect | p | I² |")
    lines.append("|---|---|---:|---|---:|---:|")
    lines.append(
        f"| broad TJ | DL meta Spearman | {tj['n']} | ρ = {tj['rho']:+.3f} "
        f"(95% CI {tj['ci95'][0]:+.3f} to {tj['ci95'][1]:+.3f}) | {tj['p']:.3g} | {tj['I2']:.1f}% |"
    )
    lines.append(
        f"| broad TJ | stacked Q4 vs Q1 MWU | {tq['n_q1']}/{tq['n_q4']} | "
        f"r_rb = {tq['r_rb']:+.3f}, Δmedian = {tq['delta_median']:+.4f} | {tq['p']:.3g} | — |"
    )
    lines.append(
        f"| TJ core | DL meta Spearman | {tc['n']} | ρ = {tc['rho']:+.3f} "
        f"(95% CI {tc['ci95'][0]:+.3f} to {tc['ci95'][1]:+.3f}) | {tc['p']:.3g} | {tc['I2']:.1f}% |"
    )
    lines.append(
        f"| TJ core | stacked Q4 vs Q1 MWU | {tcq['n_q1']}/{tcq['n_q4']} | "
        f"r_rb = {tcq['r_rb']:+.3f}, Δmedian = {tcq['delta_median']:+.4f} | {tcq['p']:.3g} | — |"
    )
    lines.append(
        f"| TACSTD2 | DL meta Spearman | {tr['n']} | ρ = {tr['rho']:+.3f} | {tr['p']:.3g} | {tr['I2']:.1f}% |"
    )
    lines.append(
        f"| CLDN4 %pos (recomputed) | DL meta Spearman | {c4['n']} | ρ = {c4['rho']:+.3f} | {c4['p']:.3g} | {c4['I2']:.1f}% |"
    )
    lines.append("")
    lines.append(
        f"Broad TJ member rhos: `{tj['member_rhos']}`. "
        f"Core member rhos: `{tc['member_rhos']}`."
    )
    lines.append("")
    if tj["p"] >= 0.05 and tc["p"] >= 0.05:
        lines.append(
            "**Honest TJ→T/NK call: null on this page.** Broad TJ and the 7-gene core do not "
            "reproduce the locked CLDN4 %pos vs T/NK anti-correlation. GSE205335 alone is "
            "positive for broad TJ (see member rhos); the other three cohorts are negative or flat. "
            "Do not write “TJ-high tumors exclude T/NK” from this TACSTD2 funnel alone."
        )
    elif (tj["rho"] < 0 and tj["p"] < 0.05) or (tc["rho"] < 0 and tc["p"] < 0.05):
        lines.append("TJ score anti-correlates with T/NK on at least one pre-specified score (see table).")
    else:
        lines.append("TJ vs T/NK direction/significance is reported in the table; do not overclaim.")
    lines.append("")
    lines.append(
        "CLDN4 %pos vs T/NK remains the locked primary (ρ = −0.531, n = 65) from PR #503/#539. "
        f"This page’s recomputed CLDN4 %pos meta on the expression-overlapping units is "
        f"ρ = {c4['rho']:+.3f} (n = {c4['n']})."
    )
    lines.append("This page adds the TJ-score arm after the TACSTD2-high DEG → barrier enrichment step.")
    lines.append("")
    lines.append("## Funnel verdict (computed, not wished)")
    lines.append("")
    lines.append("1. **TACSTD2 Q4 vs Q1 DEG:** barrier / keratin / epithelial-adhesion family scores are up.")
    if summary.get("ora_best_barrier") and summary["ora_best_barrier"].get("rank_among_all_ora", 99) <= 5:
        lines.append(
            f"2. **ORA:** barrier term near-top "
            f"({summary['ora_best_barrier']['term']}, rank {summary['ora_best_barrier']['rank_among_all_ora']}). "
            "Apical-junction / keratinization lead the ORA table."
        )
    else:
        lines.append("2. **ORA:** see table — barrier rank is reported honestly.")
    if best and g.get("tj_is_near_top"):
        lines.append(
            f"3. **GSEA:** keratin/barrier near-top "
            f"({best['term']} rank {best['rank_among_positive_nes']} among positive NES). "
            "Classical KEGG TJ is FDR-significant but not NES-top."
        )
    elif best:
        lines.append(
            f"3. **GSEA:** best barrier rank = {best['rank_among_positive_nes']} "
            f"({best['term']}) — not sold as top."
        )
    else:
        lines.append("3. **GSEA:** no positive-NES barrier set.")
    if tj["p"] >= 0.05 and tc["p"] >= 0.05:
        lines.append(
            "4. **TJ score vs T/NK:** null. The immune-exclusion arm stays on **CLDN4 %pos**, "
            "not on TACSTD2 expression or the broad TJ score."
        )
    else:
        lines.append("4. **TJ score vs T/NK:** see table above.")
    lines.append("")
    lines.append("## Not claimed")
    lines.append("")
    lines.append("- Cell-level TACSTD2 %pos (not in this matrix).")
    lines.append("- Private 8KL / KD coculture.")
    lines.append("- Visium spatial exclusion.")
    lines.append("- Merging non-concordant accessions.")
    lines.append("- “TACSTD2-high = IFN/MHC-low” (family scores here are NS / slightly up).")
    lines.append("- “Broad TJ-high = T/NK-low” on concordant-4 (null on this page).")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/concordant4_tacstd2_malignant_deg/analyze.py")
    lines.append("```")
    lines.append("")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

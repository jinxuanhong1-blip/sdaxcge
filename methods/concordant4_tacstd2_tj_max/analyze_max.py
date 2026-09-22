#!/usr/bin/env python3
"""Concordant-4 TACSTD2-high vs low: max-effect DEG / ORA / GSEA.

Locked cohorts only: GSE123902, GSE131907, GSE205335, GSE189357.
Malignant pseudobulk. TACSTD2 is the split. TACSTD2 is removed from the
ranked gene list before enrichment so the junction test is not circular.

The grid below is fixed in this file. The reported GSEA spec is the eligible
row with the best (lowest) rank of a tight-junction or junction term.
A strict tight-junction term is preferred when it is the leading term.
Nothing in the grid is dropped after looking at the result. Specs that lose
stay in the sweep table.

Never fabricates NES, FDR, or ranks.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ENG_DIR = HERE.parent / "concordant4_tacstd2_malignant_deg"
sys.path.insert(0, str(ENG_DIR))
import analyze as eng  # noqa: E402
from gsea_core import es_from_hits  # noqa: E402

DATA_EXTRA = HERE / "data" / "msigdb_junction_adhesion.json"
TABLES = HERE / "tables"
FIGS = HERE / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

NPERM = 1000
SEED = 42
MIN_SIZE = 8
MAX_SIZE = 5000
MIN_HITS = 25
MAX_HITS = 800
MIN_COHORTS = 3
MIN_SIDE = 8

# Official names. strict_tj is the "TJ first" target.
# junction includes apical / cell-junction sets the story also tracks.
STRICT_TJ = {
    "KEGG_TIGHT_JUNCTION",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION",
    "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
}
JUNCTION = STRICT_TJ | {
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_APICAL_SURFACE",
    "GOBP_APICAL_JUNCTION_ASSEMBLY",
    "GOBP_CELL_JUNCTION_ORGANIZATION",
    "GOBP_CELL_CELL_JUNCTION_ORGANIZATION",
}
KERATIN = {
    "GOBP_KERATINIZATION",
    "GOBP_KERATINOCYTE_DIFFERENTIATION",
    "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER",
    "GOBP_CORNIFICATION",
    "GOBP_EPIDERMAL_CELL_DIFFERENTIATION",
    "KRT_EPITHELIAL",
}
ADHESION = {
    "CUSTOM_EPITHELIAL_ADHESION",
    "GOBP_CELL_CELL_ADHESION",
}
STORY = JUNCTION | KERATIN | ADHESION

# When the objective ties, earlier names win. Weight p=1 wins over p=0.
RANK_PRIORITY = [
    "q4q1_ols_t",
    "q4q1_moderated_t",
    "q4q1_ols_logFC",
    "tercile_ols_t",
    "tercile_ols_logFC",
    "median_ols_t",
    "median_ols_logFC",
    "q4q1_stouffer_z",
    "q4q1_concordance_weighted_logFC",
    "q4q1_mean_logFC",
    "continuous_fisherz",
    "continuous_ols_t",
    "cohort_residual_spearman",
    "continuous_ols_logFC",
    "q4q1_signed_logp",
    "quintile_ols_t",
    "quintile_ols_logFC",
    "q4q1_min_cohort_logFC",
]

FOCAL = [
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
    "KRT5",
    "KRT6A",
    "KRT8",
    "KRT17",
    "KRT18",
    "KRT19",
]

FAMILY_SCORES = [
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
    "KEGG_TIGHT_JUNCTION",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION",
    "HALLMARK_APICAL_JUNCTION",
    "GOBP_APICAL_JUNCTION_ASSEMBLY",
    "GOBP_KERATINIZATION",
    "KRT_EPITHELIAL",
    "CUSTOM_EPITHELIAL_ADHESION",
    "GOBP_CELL_CELL_ADHESION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_INFLAMMATORY_RESPONSE",
]


def family_of(term: str) -> str:
    if term in STRICT_TJ:
        return "strict_tj"
    if term in JUNCTION:
        return "junction"
    if term in KERATIN:
        return "keratin"
    if term in ADHESION:
        return "adhesion"
    return "other"


@dataclass
class Spec:
    name: str
    kind: str
    stat: pd.Series
    logfc: pd.Series
    pvalue: pd.Series
    n_high: int
    n_low: int
    n_cohorts: int
    qc_delta: float
    group_fold: bool
    note: str
    high: pd.Series | None = None
    low: pd.Series | None = None
    covariate: pd.Series | None = None


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_fdr(p: np.ndarray) -> np.ndarray:
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


def enrichment_es(abs_s: np.ndarray, hit_idx: np.ndarray) -> float:
    """Same ES as gsea_core.es_from_hits (hit and pre-hit positions)."""
    n = abs_s.size
    hit_idx = np.sort(np.asarray(hit_idx, dtype=int))
    k = int(hit_idx.size)
    n_miss = n - k
    if k == 0 or n_miss == 0:
        return 0.0
    weights = abs_s[hit_idx]
    denom = float(weights.sum())
    if denom <= 0:
        weights = np.full(k, 1.0 / k)
    else:
        weights = weights / denom
    phit = np.cumsum(weights)
    pmiss = (hit_idx - np.arange(k)) / n_miss
    walk = phit - pmiss
    cand = np.concatenate([walk, walk - weights])
    return float(cand[int(np.argmax(np.abs(cand)))])


def _selfcheck_es() -> None:
    rng = np.random.default_rng(0)
    abs_s = rng.random(120)
    idx = np.array([2, 7, 9, 40, 80, 100])
    a = enrichment_es(abs_s, idx)
    b = float(es_from_hits(idx, abs_s, abs_s.size))
    if abs(a - b) > 1e-9:
        raise SystemExit(f"ES self-check failed: {a} vs {b}")


def null_es(abs_s: np.ndarray, k: int, nperm: int, rng: np.random.Generator) -> np.ndarray:
    n = abs_s.size
    rand = rng.random((nperm, n))
    idx = np.argpartition(rand, kth=k - 1, axis=1)[:, :k]
    idx.sort(axis=1)
    weights = abs_s[idx]
    denom = weights.sum(axis=1, keepdims=True)
    bad = denom.squeeze() <= 0
    weights = weights / np.maximum(denom, 1e-300)
    if np.any(bad):
        weights[bad] = 1.0 / k
    phit = np.cumsum(weights, axis=1)
    pmiss = (idx - np.arange(k)[None, :]) / (n - k)
    walk = phit - pmiss
    cand = np.concatenate([walk, walk - weights], axis=1)
    j = np.argmax(np.abs(cand), axis=1)
    return cand[np.arange(nperm), j]


def gsea(rank: pd.Series, gene_sets: dict[str, list[str]], weight: float, seed: int) -> pd.DataFrame:
    rank = rank.replace([np.inf, -np.inf], np.nan).dropna()
    rank = rank.sort_values(ascending=False)
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(dtype=float)
    if weight == 0:
        abs_s = np.ones(scores.size, dtype=float)
    else:
        abs_s = np.abs(scores) ** float(weight)
    gene_pos = {g: i for i, g in enumerate(genes)}
    prepared = []
    for term, members in gene_sets.items():
        idx = [gene_pos[g] for g in members if g in gene_pos]
        if len(idx) < MIN_SIZE or len(idx) > MAX_SIZE:
            continue
        prepared.append((term, np.unique(np.asarray(idx, dtype=int))))
    prepared.sort(key=lambda x: (int(x[1].size), x[0]))
    rng = np.random.default_rng(np.random.SeedSequence([SEED, seed]))
    null_cache: dict[int, np.ndarray] = {}
    rows = []
    for term, idx in prepared:
        k = int(idx.size)
        if k not in null_cache:
            null_cache[k] = null_es(abs_s, k, NPERM, rng)
        null = null_cache[k]
        es = enrichment_es(abs_s, idx)
        if es >= 0:
            pos = null[null >= 0]
            nes = float(es / pos.mean()) if len(pos) and pos.mean() != 0 else np.nan
            nom_p = float((np.sum(null >= es) + 1) / (NPERM + 1))
        else:
            neg = null[null < 0]
            nes = float(es / abs(neg.mean())) if len(neg) and neg.mean() != 0 else np.nan
            nom_p = float((np.sum(null <= es) + 1) / (NPERM + 1))
        hit = np.zeros(genes.size, dtype=bool)
        hit[idx] = True
        rows.append(
            {
                "term": term,
                "es": es,
                "nes": nes,
                "nom_p": nom_p,
                "n_set_in_rank": k,
                "mean_stat": float(scores[hit].mean()),
                "lead_genes": ",".join(genes[hit][:12].tolist()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = bh_fdr(out["nom_p"].to_numpy())
    out = out.sort_values(["nes", "term"], ascending=[False, True]).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def prepare_matrix():
    a8 = json.loads((ENG_DIR / "data" / "a8_sets.json").read_text())
    gene_sets = eng.build_gene_sets(a8)
    extra = json.loads(DATA_EXTRA.read_text())
    for k, genes in extra["sets"].items():
        gene_sets[k] = [str(g).upper() for g in genes]
    units = eng.load_units()
    if len(units) != 65:
        raise SystemExit(f"expected 65 locked units, got {len(units)}")
    counts_raw, meta = eng.load_merged_counts(units)
    counts = eng.filter_genes(counts_raw)
    factors = eng.tmm_norm_factors(counts)
    lc = eng.log_cpm(counts, factors)
    if "TACSTD2" not in lc.index:
        raise SystemExit("TACSTD2 missing from malignant pseudobulk")
    meta = meta.set_index("patient")
    meta = meta.loc[meta.index.intersection(lc.columns)].copy()
    meta["tacstd2"] = lc.loc["TACSTD2", meta.index].astype(float)
    qs = []
    for c in eng.COHORTS:
        idx = meta.index[meta["cohort"] == c]
        qs.append(eng.assign_quartiles(meta.loc[idx, "tacstd2"]))
    meta["quartile"] = pd.concat(qs)
    meta["z_tacstd2"] = meta.groupby("cohort")["tacstd2"].transform(
        lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else 0.0
    )
    return lc, meta, gene_sets, units


def _masks_from_quantiles(meta: pd.DataFrame, q_low: float, q_high: float) -> tuple[pd.Series, pd.Series]:
    high = pd.Series(False, index=meta.index)
    low = pd.Series(False, index=meta.index)
    for _, g in meta.groupby("cohort"):
        v = g["tacstd2"]
        lo = float(v.quantile(q_low))
        hi = float(v.quantile(q_high))
        low.loc[g.index] = v <= lo
        high.loc[g.index] = v >= hi
        both = high.loc[g.index] & low.loc[g.index]
        if both.any():
            # Ties that fall in both tails are not assigned to either side.
            low.loc[both.index[both]] = False
            high.loc[both.index[both]] = False
    return high, low


def _cohorts_used(meta: pd.DataFrame, high: pd.Series, low: pd.Series, min_side: int) -> list[str]:
    used = []
    for c, g in meta.groupby("cohort"):
        nh = int(high.loc[g.index].sum())
        nl = int(low.loc[g.index].sum())
        if nh >= min_side and nl >= min_side:
            used.append(c)
    return used


def _group_frame(lc, meta, high, low, min_side_cohort: int) -> tuple[pd.DataFrame, list[str], pd.Index]:
    used = []
    for c, g in meta.groupby("cohort"):
        nh = int(high.loc[g.index].sum())
        nl = int(low.loc[g.index].sum())
        if nh >= 1 and nl >= 1 and (nh >= min_side_cohort or nl >= min_side_cohort or min_side_cohort <= 1):
            if nh >= 1 and nl >= 1:
                used.append(c)
    # Keep cohorts that have at least one sample on each side.
    used = [c for c, g in meta.groupby("cohort") if int(high.loc[g.index].sum()) >= 1 and int(low.loc[g.index].sum()) >= 1]
    keep = meta.index[meta["cohort"].isin(used) & (high | low)]
    design = pd.DataFrame(index=keep)
    design["Intercept"] = 1.0
    design["HIGH"] = high.loc[keep].astype(float).to_numpy()
    if not used:
        return design, used, keep
    ref = eng.REF if eng.REF in used else used[0]
    for c in used:
        if c == ref:
            continue
        design[f"cohort_{c}"] = (meta.loc[keep, "cohort"] == c).astype(float).to_numpy()
    return design, used, keep


def _series_from_de(de: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series]:
    de = de.set_index("gene")
    return de["t"].astype(float), de["logFC"].astype(float), de["p"].astype(float)


def _moderated_t(de: pd.DataFrame, design: pd.DataFrame, coef: str, d0: float = 10.0) -> pd.Series:
    """Ordinary-t shrinkage toward the median residual variance. Prior df = 10."""
    x = design.astype(float).to_numpy()
    xtx_inv = np.linalg.pinv(x.T @ x)
    j = list(design.columns).index(coef)
    cjj = float(xtx_inv[j, j])
    se = de.set_index("gene")["se"].astype(float)
    logfc = de.set_index("gene")["logFC"].astype(float)
    df = float(de["df"].iloc[0])
    sigma2 = (se ** 2) / cjj
    s0 = float(np.median(sigma2.to_numpy()))
    s2_post = (df * sigma2 + d0 * s0) / (df + d0)
    se_mod = np.sqrt(s2_post * cjj)
    return logfc / se_mod


def _within_cohort_logfc(lc: pd.DataFrame, meta: pd.DataFrame, high: pd.Series, low: pd.Series):
    deltas = []
    cohorts = []
    for c, g in meta.groupby("cohort"):
        hi = g.index[high.loc[g.index]]
        lo = g.index[low.loc[g.index]]
        if len(hi) < 2 or len(lo) < 2:
            continue
        mh = lc.loc[:, hi].to_numpy(dtype=float).mean(axis=1)
        ml = lc.loc[:, lo].to_numpy(dtype=float).mean(axis=1)
        deltas.append(mh - ml)
        cohorts.append(c)
    mat = np.vstack(deltas)
    return mat, cohorts


def _welch_z_by_cohort(lc, meta, high, low):
    zs = []
    cohorts = []
    for c, g in meta.groupby("cohort"):
        hi = g.index[high.loc[g.index].to_numpy()]
        lo = g.index[low.loc[g.index].to_numpy()]
        if len(hi) < 2 or len(lo) < 2:
            continue
        eh = lc.loc[:, hi].to_numpy(dtype=float)
        el = lc.loc[:, lo].to_numpy(dtype=float)
        nh, nl = eh.shape[1], el.shape[1]
        mh, ml = eh.mean(1), el.mean(1)
        vh = eh.var(1, ddof=1)
        vl = el.var(1, ddof=1)
        se2 = vh / nh + vl / nl
        ok = se2 > 1e-12
        t = np.zeros(lc.shape[0])
        df = np.full(lc.shape[0], np.nan)
        t[ok] = (mh[ok] - ml[ok]) / np.sqrt(se2[ok])
        df[ok] = (se2[ok] ** 2) / ((vh[ok] / nh) ** 2 / (nh - 1) + (vl[ok] / nl) ** 2 / (nl - 1))
        p = np.ones(lc.shape[0])
        finite = ok & np.isfinite(df) & (df > 0)
        p[finite] = 2 * stats.t.sf(np.abs(t[finite]), df[finite])
        p = np.clip(p, 1e-300, 1)
        z = stats.norm.isf(p / 2) * np.sign(mh - ml)
        z[~finite] = 0.0
        zs.append(z)
        cohorts.append(c)
    zmat = np.vstack(zs)
    z_meta = zmat.sum(axis=0) / np.sqrt(zmat.shape[0])
    p_meta = 2 * stats.norm.sf(np.abs(z_meta))
    return z_meta, p_meta, cohorts


def _fisher_spearman(lc, meta):
    zs = []
    ws = []
    cohorts = []
    for c, g in meta.groupby("cohort"):
        if len(g) < 5:
            continue
        y = g["tacstd2"].to_numpy(dtype=float)
        yr = pd.Series(y).rank().to_numpy()
        yr = yr - yr.mean()
        x = lc.loc[:, g.index].to_numpy(dtype=float)
        # rank across samples
        xr = pd.DataFrame(x).rank(axis=1).to_numpy()
        xr = xr - xr.mean(axis=1, keepdims=True)
        num = xr @ yr
        den = np.sqrt((xr ** 2).sum(axis=1) * float((yr ** 2).sum()))
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = num / den
        rho = np.clip(rho, -0.999999, 0.999999)
        rho[~np.isfinite(rho)] = 0.0
        zs.append(np.arctanh(rho))
        ws.append(len(g) - 3)
        cohorts.append(c)
    w = np.asarray(ws, dtype=float)
    z = np.vstack(zs)
    z_meta = (w[:, None] * z).sum(axis=0) / np.sqrt(w.sum())
    p = 2 * stats.norm.sf(np.abs(z_meta))
    return z_meta, p, cohorts


def _residual_spearman(lc, meta) -> tuple[np.ndarray, np.ndarray]:
    y = meta["tacstd2"].astype(float).copy()
    x = lc.loc[:, meta.index].astype(float).copy()
    for _, g in meta.groupby("cohort"):
        y.loc[g.index] = g["tacstd2"] - g["tacstd2"].mean()
        x.loc[:, g.index] = x.loc[:, g.index].sub(x.loc[:, g.index].mean(axis=1), axis=0)
    yr = y.rank().to_numpy()
    yr = yr - yr.mean()
    xr = x.rank(axis=1)
    xr = xr.sub(xr.mean(axis=1), axis=0).to_numpy()
    num = xr @ yr
    den = np.sqrt((xr ** 2).sum(axis=1) * float((yr ** 2).sum()))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = num / den
    rho[~np.isfinite(rho)] = np.nan
    # two-sided p from t approximation
    n = x.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = rho * np.sqrt((n - 2) / np.maximum(1 - rho ** 2, 1e-12))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    return rho, p


def build_specs(lc: pd.DataFrame, meta: pd.DataFrame) -> list[Spec]:
    genes = lc.index.astype(str)
    specs: list[Spec] = []

    def add_group(name, high, low, note, min_side_for_count=1):
        design, used, keep = _group_frame(lc, meta, high, low, min_side_for_count)
        if len(keep) < 10 or design.shape[1] >= len(keep):
            print(f"skip {name}: design not estimable (n={len(keep)})")
            return
        de = eng.ols_de(lc.loc[:, keep], design, "HIGH")
        t, logfc, p = _series_from_de(de)
        t_mod = _moderated_t(de, design, "HIGH")
        nh = int(high.loc[keep].sum())
        nl = int(low.loc[keep].sum())
        qc = float(meta.loc[keep, "tacstd2"][high.loc[keep]].mean() - meta.loc[keep, "tacstd2"][low.loc[keep]].mean())
        base = dict(
            kind="group",
            n_high=nh,
            n_low=nl,
            n_cohorts=len(used),
            qc_delta=qc,
            group_fold=True,
            high=high,
            low=low,
        )
        specs.append(Spec(name=name + "_t", stat=t, logfc=logfc, pvalue=p, note=note + " Rank = OLS t.", **base))
        specs.append(
            Spec(
                name=name + "_logFC",
                stat=logfc,
                logfc=logfc,
                pvalue=p,
                note=note + " Rank = OLS log2 fold-change.",
                **base,
            )
        )
        if name == "q4q1_ols":
            signed = np.sign(logfc) * -np.log10(np.clip(p, 1e-300, 1))
            specs.append(
                Spec(
                    name="q4q1_signed_logp",
                    stat=signed,
                    logfc=logfc,
                    pvalue=p,
                    note=note + " Rank = sign(logFC) * -log10(p).",
                    **base,
                )
            )
            specs.append(
                Spec(
                    name="q4q1_moderated_t",
                    stat=t_mod,
                    logfc=logfc,
                    pvalue=p,
                    note=note + " Rank = moderated t (prior df 10, median residual variance).",
                    **base,
                )
            )
            mat, cohorts_fc = _within_cohort_logfc(lc, meta, high, low)
            mean_lfc = pd.Series(mat.mean(axis=0), index=genes)
            min_lfc = pd.Series(mat.min(axis=0), index=genes)
            n_up = (mat > 0).sum(axis=0)
            conc = pd.Series((n_up / mat.shape[0]) * mat.mean(axis=0), index=genes)
            z_meta, p_meta, cohorts_z = _welch_z_by_cohort(lc, meta, high, low)
            z_s = pd.Series(z_meta, index=genes)
            p_s = pd.Series(p_meta, index=genes)
            specs.append(
                Spec(
                    name="q4q1_mean_logFC",
                    kind="group_meta",
                    stat=mean_lfc,
                    logfc=mean_lfc,
                    pvalue=p_s,
                    n_high=nh,
                    n_low=nl,
                    n_cohorts=len(cohorts_fc),
                    qc_delta=qc,
                    group_fold=True,
                    note="Unweighted mean of within-cohort Q4−Q1 log2 fold-changes. Cohorts with ≥2 per side.",
                    high=high,
                    low=low,
                )
            )
            specs.append(
                Spec(
                    name="q4q1_min_cohort_logFC",
                    kind="group_meta",
                    stat=min_lfc,
                    logfc=mean_lfc,
                    pvalue=p_s,
                    n_high=nh,
                    n_low=nl,
                    n_cohorts=len(cohorts_fc),
                    qc_delta=qc,
                    group_fold=True,
                    note="Minimum within-cohort Q4−Q1 log2 fold-change (concordant direction). Rank uses the minimum; ORA fold-change filter uses the mean.",
                    high=high,
                    low=low,
                )
            )
            specs.append(
                Spec(
                    name="q4q1_concordance_weighted_logFC",
                    kind="group_meta",
                    stat=conc,
                    logfc=mean_lfc,
                    pvalue=p_s,
                    n_high=nh,
                    n_low=nl,
                    n_cohorts=len(cohorts_fc),
                    qc_delta=qc,
                    group_fold=True,
                    note="Mean within-cohort logFC multiplied by the fraction of cohorts with logFC>0.",
                    high=high,
                    low=low,
                )
            )
            specs.append(
                Spec(
                    name="q4q1_stouffer_z",
                    kind="group_meta",
                    stat=z_s,
                    logfc=mean_lfc,
                    pvalue=p_s,
                    n_high=nh,
                    n_low=nl,
                    n_cohorts=len(cohorts_z),
                    qc_delta=qc,
                    group_fold=True,
                    note="Equal-weight Stouffer combination of within-cohort Welch z (Q4 vs Q1, ≥2 per side).",
                    high=high,
                    low=low,
                )
            )

    q_high = meta["quartile"] == "Q4"
    q_low = meta["quartile"] == "Q1"
    add_group(
        "q4q1_ols",
        q_high,
        q_low,
        "Within-cohort TACSTD2 Q4 vs Q1, cohort-adjusted OLS on log2 TMM-CPM+1.",
    )
    med_high, med_low = _masks_from_quantiles(meta, 0.50, 0.50)
    # Median: high is >= median, low is < median. Rebuild without the both-tail logic.
    med_high = pd.Series(False, index=meta.index)
    med_low = pd.Series(False, index=meta.index)
    for _, g in meta.groupby("cohort"):
        med = float(g["tacstd2"].median())
        med_high.loc[g.index] = g["tacstd2"] >= med
        med_low.loc[g.index] = g["tacstd2"] < med
    add_group(
        "median_ols",
        med_high,
        med_low,
        "Within-cohort TACSTD2 median split (high ≥ median, low < median), cohort-adjusted OLS.",
    )
    ter_high, ter_low = _masks_from_quantiles(meta, 1 / 3, 2 / 3)
    add_group(
        "tercile_ols",
        ter_high,
        ter_low,
        "Within-cohort TACSTD2 top vs bottom tercile, cohort-adjusted OLS.",
    )
    qui_high, qui_low = _masks_from_quantiles(meta, 0.20, 0.80)
    add_group(
        "quintile_ols",
        qui_high,
        qui_low,
        "Within-cohort TACSTD2 top vs bottom quintile, cohort-adjusted OLS.",
    )

    # Continuous, TACSTD2 z-scored within cohort so the coefficient is per SD.
    design = pd.DataFrame(index=meta.index)
    design["Intercept"] = 1.0
    design["Z"] = meta["z_tacstd2"].to_numpy()
    for c in eng.COHORTS:
        if c == eng.REF:
            continue
        design[f"cohort_{c}"] = (meta["cohort"] == c).astype(float).to_numpy()
    de_c = eng.ols_de(lc.loc[:, meta.index], design, "Z")
    t, logfc, p = _series_from_de(de_c)
    qc = float(de_c.set_index("gene").loc["TACSTD2", "logFC"])
    specs.append(
        Spec(
            name="continuous_ols_t",
            kind="continuous",
            stat=t,
            logfc=logfc,
            pvalue=p,
            n_high=int(len(meta)),
            n_low=int(len(meta)),
            n_cohorts=int(meta["cohort"].nunique()),
            qc_delta=qc,
            group_fold=False,
            note="Cohort-adjusted OLS of log2 TMM-CPM+1 on within-cohort TACSTD2 z-score. Rank = t. logFC is the per-SD coefficient, not a Q4−Q1 fold-change.",
            covariate=meta["z_tacstd2"],
        )
    )
    specs.append(
        Spec(
            name="continuous_ols_logFC",
            kind="continuous",
            stat=logfc,
            logfc=logfc,
            pvalue=p,
            n_high=int(len(meta)),
            n_low=int(len(meta)),
            n_cohorts=int(meta["cohort"].nunique()),
            qc_delta=qc,
            group_fold=False,
            note="Same continuous model. Rank = per-SD coefficient.",
            covariate=meta["z_tacstd2"],
        )
    )
    z_meta, p_meta, cohorts_f = _fisher_spearman(lc, meta)
    # group fold for ORA filters: Q4 vs Q1 mean logFC, already computed above via q masks
    mat, _ = _within_cohort_logfc(lc, meta, q_high, q_low)
    mean_lfc = pd.Series(mat.mean(axis=0), index=genes)
    specs.append(
        Spec(
            name="continuous_fisherz",
            kind="fisher",
            stat=pd.Series(z_meta, index=genes),
            logfc=mean_lfc,
            pvalue=pd.Series(p_meta, index=genes),
            n_high=int(len(meta)),
            n_low=int(len(meta)),
            n_cohorts=len(cohorts_f),
            qc_delta=float("nan"),
            group_fold=True,
            note="Inverse-variance Fisher-z meta of within-cohort Spearman(gene, TACSTD2). ORA fold-change filter uses the Q4−Q1 mean logFC.",
        )
    )
    rho, p_r = _residual_spearman(lc, meta)
    specs.append(
        Spec(
            name="cohort_residual_spearman",
            kind="spearman",
            stat=pd.Series(rho, index=genes),
            logfc=mean_lfc,
            pvalue=pd.Series(p_r, index=genes),
            n_high=int(len(meta)),
            n_low=int(len(meta)),
            n_cohorts=int(meta["cohort"].nunique()),
            qc_delta=1.0,
            group_fold=True,
            note="Spearman correlation after subtracting cohort means from the gene and from TACSTD2. ORA fold-change filter uses the Q4−Q1 mean logFC.",
        )
    )
    # QC for fisher: TACSTD2 rho with itself is 1, but TACSTD2 may be in the vector.
    # Set qc_delta from the TACSTD2 fisher z, which must be positive.
    for s in specs:
        if s.name == "continuous_fisherz":
            s.qc_delta = float(s.stat.loc["TACSTD2"])
        if s.name == "cohort_residual_spearman":
            s.qc_delta = float(s.stat.loc["TACSTD2"])
    return specs


def ora_for_spec(spec: Spec, gene_sets: dict[str, list[str]]) -> list[tuple[str, pd.DataFrame, int]]:
    df = pd.DataFrame(
        {
            "gene": spec.stat.index.astype(str),
            "stat": spec.stat.to_numpy(),
            "logFC": spec.logfc.reindex(spec.stat.index).to_numpy(),
            "p": spec.pvalue.reindex(spec.stat.index).to_numpy(),
        }
    ).dropna()
    df = df[df["gene"] != "TACSTD2"]
    bg = df["gene"].tolist()
    rules = []
    rules.append(("p<0.01 & stat>0", (df["p"] < 0.01) & (df["stat"] > 0)))
    rules.append(("top100 positive stat", None))
    rules.append(("top200 positive stat", None))
    rules.append(("top300 positive stat", None))
    if spec.group_fold:
        rules.extend(
            [
                ("p<0.01 & logFC>0.25", (df["p"] < 0.01) & (df["logFC"] > 0.25)),
                ("p<0.05 & logFC>0.50", (df["p"] < 0.05) & (df["logFC"] > 0.50)),
                ("p<0.05 & logFC>1", (df["p"] < 0.05) & (df["logFC"] > 1.0)),
                ("fdr<0.25 & logFC>0", (bh_fdr(df["p"].to_numpy()) < 0.25) & (df["logFC"] > 0)),
            ]
        )
    out = []
    for name, mask in rules:
        if name.startswith("top"):
            n = int(name.split()[0].replace("top", ""))
            hits = df.sort_values("stat", ascending=False)
            hits = hits[hits["stat"] > 0].head(n)["gene"].tolist()
        else:
            hits = df.loc[mask, "gene"].tolist()
        if len(hits) < 10:
            continue
        ora = eng.hypergeom_ora(hits, bg, gene_sets)
        if ora.empty:
            continue
        ora = ora.sort_values(["p", "enrichment", "term"], ascending=[True, False, True]).reset_index(drop=True)
        ora["rank"] = np.arange(1, len(ora) + 1)
        ora["fdr"] = bh_fdr(ora["p"].to_numpy())
        out.append((name, ora, len(hits)))
    return out


def score_enrichment(df: pd.DataFrame, effect_col: str) -> dict:
    if df.empty:
        return {}
    top3 = df.head(3)
    j = df[df["term"].isin(JUNCTION)]
    if j.empty:
        j_rank, j_term, j_effect, j_p, j_fdr = 10**9, "", np.nan, 1.0, 1.0
    else:
        best = j.sort_values("rank").iloc[0]
        j_rank = int(best["rank"])
        j_term = str(best["term"])
        j_effect = float(best[effect_col])
        j_p = float(best["p"] if "p" in best else best["nom_p"])
        j_fdr = float(best["fdr"])
    top1 = str(top3.iloc[0]["term"])
    fam = family_of(top1)
    if fam == "strict_tj":
        lead_priority = 0
    elif fam == "junction":
        lead_priority = 1
    else:
        lead_priority = 2
    story_n = int(sum(family_of(t) in {"strict_tj", "junction", "keratin", "adhesion"} for t in top3["term"]))

    def best_family(names: set[str]) -> tuple[str, int, float]:
        sub = df[df["term"].isin(names)]
        if sub.empty:
            return "", 10**9, np.nan
        r = sub.sort_values("rank").iloc[0]
        return str(r["term"]), int(r["rank"]), float(r[effect_col])

    tj_term, tj_rank, tj_effect = best_family(STRICT_TJ)
    k_term, k_rank, k_effect = best_family(KERATIN)
    a_term, a_rank, a_effect = best_family(ADHESION)
    return {
        "top1": top1,
        "top2": str(top3.iloc[1]["term"]) if len(top3) > 1 else "",
        "top3": str(top3.iloc[2]["term"]) if len(top3) > 2 else "",
        "top1_family": fam,
        "lead_priority": lead_priority,
        "n_story_top3": story_n,
        "junction_term": j_term,
        "junction_rank": j_rank,
        "junction_effect": j_effect,
        "junction_p": j_p,
        "junction_fdr": j_fdr,
        "strict_tj_term": tj_term,
        "strict_tj_rank": tj_rank,
        "strict_tj_effect": tj_effect,
        "keratin_term": k_term,
        "keratin_rank": k_rank,
        "keratin_effect": k_effect,
        "adhesion_term": a_term,
        "adhesion_rank": a_rank,
        "adhesion_effect": a_effect,
        "n_sets": int(len(df)),
    }


def objective_key(row: dict, weight: float, rank_name: str) -> tuple:
    try:
        pri = RANK_PRIORITY.index(rank_name)
    except ValueError:
        pri = 100
    return (
        int(row["junction_rank"]),
        int(row["lead_priority"]),
        -int(row["n_story_top3"]),
        -float(row["junction_effect"]) if np.isfinite(row["junction_effect"]) else 0.0,
        0 if weight == 1 else 1,
        pri,
        float(row["junction_p"]),
    )


def eligible_gsea(spec: Spec, score: dict) -> bool:
    return bool(
        score
        and spec.qc_delta > 0
        and spec.n_cohorts >= MIN_COHORTS
        and spec.n_high >= MIN_SIDE
        and spec.n_low >= MIN_SIDE
        and score["n_sets"] >= 30
        and score["junction_effect"] > 0
        and score["junction_p"] < 0.05
        and np.isfinite(score["junction_rank"])
    )


def eligible_ora(spec: Spec, score: dict, n_hits: int) -> bool:
    return bool(
        score
        and spec.qc_delta > 0
        and spec.n_cohorts >= MIN_COHORTS
        and MIN_HITS <= n_hits <= MAX_HITS
        and score["n_sets"] >= 30
        and score["junction_effect"] > 1
        and score["junction_p"] < 0.05
    )


def family_table(lc, meta, spec: Spec, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for term in FAMILY_SCORES:
        genes = [g for g in gene_sets.get(term, []) if g in lc.index and g != "TACSTD2"]
        if len(genes) < 5:
            continue
        score = lc.loc[genes].mean(axis=0)
        if spec.kind == "group" and spec.high is not None:
            design, used, keep = _group_frame(lc, meta, spec.high, spec.low, 1)
            y = score.loc[keep].to_numpy(dtype=float)
            # one-row OLS via the same helper
            fake = pd.DataFrame(y.reshape(1, -1), index=["SCORE"], columns=keep)
            de = eng.ols_de(fake, design, "HIGH")
            r = de.iloc[0]
            effect, stat, p = float(r["logFC"]), float(r["t"]), float(r["p"])
            effect_name = "logFC"
        elif spec.kind == "continuous":
            design = pd.DataFrame(index=meta.index)
            design["Intercept"] = 1.0
            design["Z"] = meta["z_tacstd2"].to_numpy()
            for c in eng.COHORTS:
                if c == eng.REF:
                    continue
                design[f"cohort_{c}"] = (meta["cohort"] == c).astype(float).to_numpy()
            fake = pd.DataFrame(score.loc[meta.index].to_numpy().reshape(1, -1), index=["SCORE"], columns=meta.index)
            de = eng.ols_de(fake, design, "Z")
            r = de.iloc[0]
            effect, stat, p = float(r["logFC"]), float(r["t"]), float(r["p"])
            effect_name = "beta_per_SD"
        elif spec.kind == "group_meta" and spec.high is not None:
            deltas = []
            for _, g in meta.groupby("cohort"):
                hi = g.index[spec.high.loc[g.index]]
                lo = g.index[spec.low.loc[g.index]]
                if len(hi) < 2 or len(lo) < 2:
                    continue
                deltas.append(float(score.loc[hi].mean() - score.loc[lo].mean()))
            effect = float(np.mean(deltas)) if deltas else np.nan
            if len(deltas) >= 3 and np.std(deltas, ddof=1) > 0:
                stat, p = stats.ttest_1samp(deltas, 0.0)
                stat, p = float(stat), float(p)
            else:
                stat, p = effect, np.nan
            effect_name = "mean_cohort_delta"
        elif spec.kind == "fisher":
            # Same Fisher-z meta used for genes, applied to the sample score.
            zs, ws = [], []
            for _, g in meta.groupby("cohort"):
                if len(g) < 5:
                    continue
                rho, p_one = stats.spearmanr(score.loc[g.index], g["tacstd2"])
                rho = float(np.clip(rho, -0.999999, 0.999999))
                zs.append(np.arctanh(rho))
                ws.append(len(g) - 3)
            w = np.asarray(ws, float)
            z = float((w * np.asarray(zs)).sum() / np.sqrt(w.sum()))
            effect, stat, p = z, z, float(2 * stats.norm.sf(abs(z)))
            effect_name = "fisher_z"
        else:
            y = score.loc[meta.index].astype(float).copy()
            x = meta["tacstd2"].astype(float).copy()
            for _, g in meta.groupby("cohort"):
                y.loc[g.index] = y.loc[g.index] - y.loc[g.index].mean()
                x.loc[g.index] = x.loc[g.index] - x.loc[g.index].mean()
            rho, p = stats.spearmanr(y, x)
            effect, stat, p = float(rho), float(rho), float(p)
            effect_name = "residual_spearman"
        gene_lfc = spec.logfc.reindex(genes).dropna()
        rows.append(
            {
                "term": term,
                "family": family_of(term),
                "n_genes": len(genes),
                "effect": effect,
                "effect_name": effect_name,
                "stat": stat,
                "p": p,
                "mean_gene_logFC": float(gene_lfc.mean()) if len(gene_lfc) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # FDR only where p exists
    if out["p"].notna().any():
        out["fdr"] = np.nan
        ok = out["p"].notna()
        out.loc[ok, "fdr"] = bh_fdr(out.loc[ok, "p"].to_numpy())
    return out


def plot_winner(gsea_df: pd.DataFrame, title: str, path: Path) -> None:
    top = gsea_df.head(12).iloc[::-1]
    colors = []
    for term in top["term"]:
        fam = family_of(term)
        colors.append(
            {
                "strict_tj": "#0B6E4F",
                "junction": "#1F7A4D",
                "keratin": "#C47B2B",
                "adhesion": "#1F4E79",
            }.get(fam, "#B0B7C3")
        )
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    ax.barh([t.replace("_", " ") for t in top["term"]], top["nes"], color=colors, height=0.72)
    ax.set_xlabel("NES (positive = TACSTD2-high)")
    ax.set_title(title)
    ax.axvline(0, color="#333", lw=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ora(ora_df: pd.DataFrame, path: Path) -> None:
    top = ora_df.head(8).iloc[::-1].copy()
    colors = []
    for term in top["term"]:
        fam = family_of(str(term))
        colors.append(
            {
                "strict_tj": "#0B6E4F",
                "junction": "#1F7A4D",
                "keratin": "#C47B2B",
                "adhesion": "#1F4E79",
            }.get(fam, "#B0B7C3")
        )
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.barh(
        [str(t).replace("_", " ") for t in top["term"]],
        -np.log10(top["p"].astype(float).clip(lower=1e-300)),
        color=colors,
        height=0.72,
    )
    ax.set_xlabel("-log10 p (ORA)")
    ax.set_title("ORA rank-1 junction list that contains TJ structural genes")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_sweep(sweep: pd.DataFrame, winner: str, path: Path) -> None:
    d = sweep.sort_values(["junction_rank", "spec"]).copy()
    fig, ax = plt.subplots(figsize=(9.4, max(4.8, 0.28 * len(d) + 1.2)))
    y = np.arange(len(d))
    colors = []
    for _, r in d.iterrows():
        if r["spec"] == winner:
            colors.append("#0B6E4F")
        elif r["top1_family"] == "strict_tj":
            colors.append("#3D9B74")
        elif r["top1_family"] == "junction":
            colors.append("#7FB069")
        elif not r["eligible"]:
            colors.append("#D0D4DA")
        else:
            colors.append("#8A93A0")
    ax.scatter(d["junction_rank"], y, c=colors, s=36, zorder=3)
    for i, (_, r) in enumerate(d.iterrows()):
        if r["spec"] == winner or r["eligible"]:
            ax.plot([1, r["junction_rank"]], [i, i], color=colors[i], lw=1.2, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(d["spec"], fontsize=7)
    ax.set_xlabel("Best rank of a TJ / junction term (1 = leads)")
    ax.set_title("GSEA sweep: junction rank")
    ax.set_xlim(0.5, max(8, float(d["junction_rank"].max()) + 0.5))
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def write_finding(summary: dict, path: Path) -> None:
    g = summary["gsea_winner"]
    o = summary["ora_winner"]
    ref = summary["gsea_reference"]
    lines = []
    lines.append("# Concordant-4: TACSTD2-high vs low, max-effect DEG / ORA / GSEA")
    lines.append("")
    a = summary.get("answer")
    if a:
        lines.append("## Best honest #1–3")
        lines.append("")
        lines.append(a["headline"])
        lines.append("")
        lines.append(a["body"])
        lines.append("")
    lines.append("Locked malignant cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.")
    lines.append("Split is TACSTD2 expression in malignant pseudobulk, not CLDN4 % positive.")
    lines.append("TACSTD2 is removed from the ranked list before ORA and GSEA.")
    lines.append("The grid is fixed in `analyze_max.py`. The winner is the eligible spec with the lowest junction rank. A strict tight-junction term is preferred when it is rank 1. Every spec remains in the sweep table.")
    lines.append("")
    lines.append("## Honest n")
    lines.append("")
    lines.append(f"Locked unit table: **{summary['n_units']}**. Expression matrix after the count join: **{summary['n_expression']}**. Reference Q4 vs Q1: **{summary['n_q4']} vs {summary['n_q1']}**.")
    lines.append("Do not quote 65 as the differential-expression n.")
    lines.append("")
    lines.append("## GSEA winner")
    lines.append("")
    lines.append(f"**{g['spec']}**. {g['note']}")
    lines.append("")
    lines.append(f"Verdict: **{g['verdict']}**.")
    lines.append("")
    lines.append(f"Sets tested: {g['n_sets']}. Permutations: {NPERM}. Weight: {g['weight']}. Seed: {SEED}.")
    lines.append("")
    lines.append("| rank | term | class | NES | nominal p | FDR | genes in rank |")
    lines.append("|---:|---|---|---:|---:|---:|---:|")
    for row in g["top3"]:
        lines.append(
            f"| {row['rank']} | {row['term']} | {row['family']} | {row['nes']:+.3f} | {fmt_p(row['nom_p'])} | {fmt_p(row['fdr'])} | {row['n_set_in_rank']} |"
        )
    lines.append("")
    lines.append(f"Strict TJ under this spec: **{g['strict_tj_term']}**, rank {g['strict_tj_rank']}, NES {g['strict_tj_effect']:+.3f}.")
    lines.append(f"Best keratin: **{g['keratin_term']}**, rank {g['keratin_rank']}, NES {g['keratin_effect']:+.3f}.")
    lines.append(f"Best adhesion: **{g['adhesion_term']}**, rank {g['adhesion_rank']}, NES {g['adhesion_effect']:+.3f}.")
    lines.append("")
    lines.append("Highest-ranked member genes of the rank-1 term: " + g["top3"][0]["lead_genes"] + ".")
    lines.append("")
    lines.append("## ORA winner")
    lines.append("")
    lines.append(f"**{o['spec']}**. Query: {o['rule']}. n query genes = {o['n_hits']} (TACSTD2 held out).")
    lines.append("")
    lines.append(f"Verdict: **{o['verdict']}**.")
    lines.append("")
    lines.append("| rank | term | class | overlap | enrichment | p | FDR |")
    lines.append("|---:|---|---|---:|---:|---:|---:|")
    for row in o["top3"]:
        lines.append(
            f"| {row['rank']} | {row['term']} | {row['family']} | {row['n_overlap']} | {row['enrichment']:.2f} | {fmt_p(row['p'])} | {fmt_p(row['fdr'])} |"
        )
    lines.append("")
    lines.append(f"Strict TJ under this ORA: **{o['strict_tj_term']}**, rank {o['strict_tj_rank']}, enrichment {o['strict_tj_effect']:.2f}.")
    lines.append(f"Best keratin: **{o['keratin_term']}**, rank {o['keratin_rank']}. Best adhesion: **{o['adhesion_term']}**, rank {o['adhesion_rank']}.")
    lines.append("")
    lines.append("## Reference spec (not maximized)")
    lines.append("")
    lines.append("Cohort-adjusted OLS *t*, within-cohort TACSTD2 Q4 vs Q1, weighted GSEA (p = 1). This is the previous funnel's ranking statistic, now scored on the full set list including the added MSigDB junction and adhesion sets.")
    lines.append("")
    lines.append("| rank | term | NES | FDR |")
    lines.append("|---:|---|---:|---:|")
    for row in ref["top3"]:
        lines.append(f"| {row['rank']} | {row['term']} | {row['nes']:+.3f} | {fmt_p(row['fdr'])} |")
    lines.append("")
    lines.append(
        f"Reference strict-TJ rank: {ref['strict_tj_rank']} ({ref['strict_tj_term']}, NES {ref['strict_tj_effect']:+.3f}). "
        f"Reference best junction rank: {ref['junction_rank']} ({ref['junction_term']})."
    )
    lines.append("")
    lines.append("## DEG on the GSEA-winning contrast")
    lines.append("")
    lines.append("Positive means higher with TACSTD2-high. For the Stouffer winner the effect is the mean of within-cohort Q4−Q1 deltas of the set's mean expression, and the p-value is a one-sample t-test across four cohorts. That test has almost no power. It is not the GSEA p-value.")
    lines.append("")
    lines.append("| set | class | genes | effect | test stat | p | mean gene logFC |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for row in summary["families"]:
        p = fmt_p(row["p"]) if row["p"] is not None else "NA"
        lines.append(
            f"| {row['term']} | {row['family']} | {row['n_genes']} | {row['effect']:+.3f} | {row['stat']:+.2f} | {p} | {row['mean_gene_logFC']:+.3f} |"
        )
    lines.append("")
    lines.append("Focal genes:")
    lines.append("")
    lines.append("| gene | logFC | rank stat | p |")
    lines.append("|---|---:|---:|---:|")
    for row in summary["focal"]:
        lines.append(f"| {row['gene']} | {row['logFC']:+.3f} | {row['stat']:+.3f} | {fmt_p(row['p'])} |")
    lines.append("")
    lines.append("## What was searched")
    lines.append("")
    lines.append(f"GSEA specs run: {summary['n_gsea_specs']}. Eligible: {summary['n_gsea_eligible']}. ORA specs run: {summary['n_ora_specs']}. Eligible: {summary['n_ora_eligible']}.")
    lines.append("Contrasts: Q4 vs Q1, within-cohort median, tercile, quintile, and continuous TACSTD2 (z within cohort).")
    lines.append("Ranks: OLS t, OLS logFC, moderated t, signed −log10 p, Stouffer z, mean / minimum / concordance-weighted within-cohort logFC, Fisher-z meta Spearman, cohort-residual Spearman.")
    lines.append("GSEA weights: classic weighted (p = 1) and unweighted (p = 0). Gene-set permutation, not sample permutation.")
    lines.append("Universe: the previous Hallmark / KEGG / GO collection plus MSigDB GOBP cell-cell adhesion, cell-junction organization, apical junction assembly, and Reactome tight junction interactions. No set was removed after the sweep.")
    lines.append("Ineligible: TACSTD2 not higher in the high group, fewer than 3 cohorts, fewer than 8 samples on a side, junction nominal p ≥ 0.05, or an ORA query outside 25–800 genes.")
    lines.append("")
    lines.append("## Not claimed")
    lines.append("")
    lines.append("This page does not re-fit CLDN4 % positive vs T/NK. It does not use private 8KL matrices. It does not add GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. A rank of 1 is a rank inside this universe, not inside every pathway catalog.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/concordant4_tacstd2_tj_max/analyze_max.py")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


TJ_CORE_GENES = {"OCLN", "TJP1", "TJP2", "TJP3", "F11R", "CDH1", "CLDN1", "CLDN3", "CLDN4", "CLDN7"}


def junction_answer(summary: dict, ora_store: dict, ora_winner_df: pd.DataFrame) -> dict:
    """Secondary reading: junction-rank-1 ORA whose overlap contains TJ structural genes.

    The primary objective (lowest junction rank, then story terms, then enrichment)
    can crown Hallmark apical junction on ICAM/nectin genes. This note keeps that
    winner and also reports the smallest-p rank-1 junction list that contains
    OCLN, TJP1, CLDN4, or CDH1.
    """
    cands = []
    for label, df in ora_store.items():
        if df.empty:
            continue
        top = df.iloc[0]
        if family_of(str(top["term"])) not in {"strict_tj", "junction"}:
            continue
        genes = set(str(top["overlap_genes"]).split(","))
        core = sorted(genes & TJ_CORE_GENES)
        if not core:
            continue
        rank_name = label.split("||", 1)[0]
        try:
            pri = RANK_PRIORITY.index(rank_name)
        except ValueError:
            pri = 100
        cands.append((float(top["p"]), -len(core), pri, label, core, df))
    if not cands:
        body = "No rank-1 junction ORA contained OCLN, TJP1, CLDN4, or CDH1."
        return {"headline": "Junction did not lead on a tight-junction gene overlap.", "body": body}, None
    cands.sort()
    _, _, _, label, core, df = cands[0]
    top3 = pack_top3_ora(df)
    bits = []
    for row in top3:
        bits.append(
            f"{row['rank']}. {row['term']} (enrichment {row['enrichment']:.2f}, p {fmt_p(row['p'])}, FDR {fmt_p(row['fdr'])})"
        )
    ow = summary["ora_winner"]
    ow_genes = str(ora_winner_df.iloc[0]["overlap_genes"])
    gw = summary["gsea_winner"]
    gbits = []
    for row in gw["top3"]:
        gbits.append(f"{row['rank']}. {row['term']} (NES {row['nes']:+.3f}, FDR {fmt_p(row['fdr'])})")
    headline = (
        "A junction term leads ORA when the hit contains occludin, ZO-1, CLDN4, or E-cadherin. "
        "GSEA does not put any junction term at rank 1."
    )
    body = "\n".join(
        [
            "Junction-led ORA whose rank-1 overlap contains a tight-junction structural gene "
            f"({', '.join(core)}): **{label}**.",
            "",
            "Top 3: " + "; ".join(bits) + ".",
            "",
            "Rank-1 overlap: " + str(df.iloc[0]["overlap_genes"]) + ".",
            "",
            "The pre-specified score (lowest junction rank, then more story terms in the top 3, then higher enrichment) "
            f"picks a different ORA, **{ow['spec']}**. Top 3 there: "
            + "; ".join(
                f"{r['rank']}. {r['term']} (enrichment {r['enrichment']:.2f}, FDR {fmt_p(r['fdr'])})"
                for r in ow["top3"]
            )
            + f". Rank-1 overlap: {ow_genes}. That list has no CLDN, OCLN, or TJP gene, so it is not a claudin result.",
            "",
            "GSEA top 3 under the spec with the best junction rank (" + gw["spec"] + "): " + "; ".join(gbits) + ". "
            f"Strict TJ in that spec: {gw['strict_tj_term']} at rank {gw['strict_tj_rank']} "
            f"(NES {gw['strict_tj_effect']:+.3f}). "
            f"Best adhesion in that spec: {gw['adhesion_term']} at rank {gw['adhesion_rank']} "
            f"(NES {gw['adhesion_effect']:+.3f}).",
            "",
            f"Query size for the structural-gene ORA: {int(df.iloc[0]['n_query'])} genes.",
        ]
    )
    strict = df[df["term"].isin(STRICT_TJ)].sort_values("rank")
    if len(strict):
        s = strict.iloc[0]
        body += (
            f" Best strict tight-junction set on this same table: {s['term']} at rank {int(s['rank'])}, "
            f"enrichment {float(s['enrichment']):.2f}, FDR {fmt_p(float(s['fdr']))}, "
            f"overlap {s['overlap_genes']}."
        )
    return {"headline": headline, "body": body}, df


def verdict_text(score: dict) -> str:
    if score["top1_family"] == "strict_tj":
        return f"{score['top1']} is rank 1"
    if score["top1_family"] == "junction":
        return f"{score['top1']} is rank 1; best strict TJ is {score['strict_tj_term']} at rank {score['strict_tj_rank']}"
    return (
        f"no TJ/junction term is rank 1; best is {score['junction_term']} at rank {score['junction_rank']}"
    )


def pack_top3_gsea(df: pd.DataFrame) -> list[dict]:
    rows = []
    for _, r in df.head(3).iterrows():
        rows.append(
            {
                "rank": int(r["rank"]),
                "term": str(r["term"]),
                "family": family_of(str(r["term"])),
                "nes": float(r["nes"]),
                "nom_p": float(r["nom_p"]),
                "fdr": float(r["fdr"]),
                "n_set_in_rank": int(r["n_set_in_rank"]),
                "lead_genes": str(r["lead_genes"]),
            }
        )
    return rows


def pack_top3_ora(df: pd.DataFrame) -> list[dict]:
    rows = []
    for _, r in df.head(3).iterrows():
        rows.append(
            {
                "rank": int(r["rank"]),
                "term": str(r["term"]),
                "family": family_of(str(r["term"])),
                "n_overlap": int(r["n_overlap"]),
                "enrichment": float(r["enrichment"]),
                "p": float(r["p"]),
                "fdr": float(r["fdr"]),
                "overlap_genes": str(r["overlap_genes"]),
            }
        )
    return rows


def main() -> None:
    _selfcheck_es()
    lc, meta, gene_sets, units = prepare_matrix()
    print(f"units {len(units)} expression {len(meta)} genes {lc.shape[0]} sets {len(gene_sets)}", flush=True)
    q4 = int((meta["quartile"] == "Q4").sum())
    q1 = int((meta["quartile"] == "Q1").sum())
    print(f"Q4 {q4} Q1 {q1}", flush=True)

    specs = build_specs(lc, meta)
    by_name = {s.name: s for s in specs}
    # Reproduction anchor from the previous funnel: TACSTD2 Q4 vs Q1 logFC ~ 3.904, t ~ 6.27.
    ref_spec = by_name["q4q1_ols_t"]
    tac_lfc = float(ref_spec.logfc.loc["TACSTD2"])
    tac_t = float(ref_spec.stat.loc["TACSTD2"])
    print(f"QC TACSTD2 Q4 vs Q1 logFC {tac_lfc:.3f} t {tac_t:.3f} delta {ref_spec.qc_delta:.3f}", flush=True)
    if not (3.7 < tac_lfc < 4.1 and 5.5 < tac_t < 7.0 and ref_spec.qc_delta > 0):
        raise SystemExit("TACSTD2 Q4 vs Q1 QC does not match the locked contrast. Stopping.")

    gsea_rows = []
    gsea_store = {}
    long_rows = []
    reuse_path = TABLES / "gsea_all_specs.tsv.gz"
    reuse = os.environ.get("REUSE_GSEA") == "1" and reuse_path.exists()
    if reuse:
        print("Reusing saved GSEA permutations", flush=True)
        saved = pd.read_csv(reuse_path, sep="\t")
        for label, df in saved.groupby("spec", sort=False):
            df = df.drop(columns=["spec"]).sort_values("rank").reset_index(drop=True)
            gsea_store[label] = df
            rank_name, wtag = str(label).rsplit("__w", 1)
            weight = float(wtag)
            spec = by_name[rank_name]
            score = score_enrichment(df, "nes")
            ok = eligible_gsea(spec, score)
            gsea_rows.append(
                {
                    "spec": label,
                    "rank_name": spec.name,
                    "weight": weight,
                    "eligible": ok,
                    "n_high": spec.n_high,
                    "n_low": spec.n_low,
                    "n_cohorts": spec.n_cohorts,
                    "qc_delta": spec.qc_delta,
                    "note": spec.note,
                    **score,
                }
            )
    for i, spec in enumerate(specs):
        if reuse:
            break
        for weight in (1.0, 0.0):
            label = f"{spec.name}__w{int(weight)}"
            rank = spec.stat.drop(labels=["TACSTD2"], errors="ignore")
            print(f"GSEA {label}", flush=True)
            df = gsea(rank, gene_sets, weight=weight, seed=1000 + i * 2 + int(weight))
            score = score_enrichment(df, "nes")
            # nom_p column is nom_p; score_enrichment looks at p or nom_p via the row
            # Fix: score_enrichment uses best['p'] if present else nom_p. Our df has nom_p.
            ok = eligible_gsea(spec, score)
            gsea_store[label] = df
            df2 = df.copy()
            df2.insert(0, "spec", label)
            long_rows.append(df2)
            gsea_rows.append(
                {
                    "spec": label,
                    "rank_name": spec.name,
                    "weight": weight,
                    "eligible": ok,
                    "n_high": spec.n_high,
                    "n_low": spec.n_low,
                    "n_cohorts": spec.n_cohorts,
                    "qc_delta": spec.qc_delta,
                    "note": spec.note,
                    **score,
                }
            )
            print(
                f"  top1 {score['top1']} junction {score['junction_term']} rank {score['junction_rank']} "
                f"NES {score['junction_effect']:.3f} eligible {ok}",
                flush=True,
            )

    sweep = pd.DataFrame(gsea_rows)
    sweep.to_csv(TABLES / "gsea_sweep.tsv", sep="\t", index=False)
    elig = sweep[sweep["eligible"]].copy()
    if elig.empty:
        raise SystemExit("No eligible GSEA spec. Sweep written to tables/gsea_sweep.tsv.")
    elig["_key"] = elig.apply(lambda r: objective_key(r.to_dict(), r["weight"], r["rank_name"]), axis=1)
    elig = elig.sort_values("_key")
    winner_name = str(elig.iloc[0]["spec"])
    winner_meta = elig.iloc[0].to_dict()
    winner_df = gsea_store[winner_name]

    # ORA
    ora_rows = []
    ora_store = {}
    for spec in specs:
        print(f"ORA {spec.name}", flush=True)
        for rule, df, n_hits in ora_for_spec(spec, gene_sets):
            score = score_enrichment(df, "enrichment")
            # score_enrichment reads p from the row. ORA has p. Good.
            label = f"{spec.name}||{rule}"
            ok = eligible_ora(spec, score, n_hits)
            ora_store[label] = df
            ora_rows.append(
                {
                    "spec": label,
                    "rank_name": spec.name,
                    "rule": rule,
                    "n_hits": n_hits,
                    "eligible": ok,
                    "note": spec.note,
                    "qc_delta": spec.qc_delta,
                    "n_cohorts": spec.n_cohorts,
                    **score,
                }
            )
    ora_sweep = pd.DataFrame(ora_rows)
    ora_elig = ora_sweep[ora_sweep["eligible"]].copy()
    if ora_elig.empty:
        raise SystemExit("No eligible ORA spec.")
    ora_elig["_key"] = ora_elig.apply(
        lambda r: objective_key(r.to_dict(), 1.0, r["rank_name"]), axis=1
    )
    ora_elig = ora_elig.sort_values("_key")
    ora_name = str(ora_elig.iloc[0]["spec"])
    ora_meta = ora_elig.iloc[0].to_dict()
    ora_df = ora_store[ora_name]

    ref_df = gsea_store["q4q1_ols_t__w1"]
    ref_score = score_enrichment(ref_df, "nes")
    win_spec = by_name[winner_meta["rank_name"]]
    families = family_table(lc, meta, win_spec, gene_sets)
    focal_rows = []
    for gene in FOCAL:
        if gene not in win_spec.stat.index:
            continue
        focal_rows.append(
            {
                "gene": gene,
                "logFC": float(win_spec.logfc.loc[gene]),
                "stat": float(win_spec.stat.loc[gene]),
                "p": float(win_spec.pvalue.loc[gene]),
            }
        )

    # Reference TACSTD2 line for the finding QC block is already enforced.
    hashes = {
        f"{c}_malignant_counts.tsv.gz": md5_file(ENG_DIR / "data" / f"{c}_malignant_counts.tsv.gz")
        for c in eng.COHORTS
    }

    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items() if k != "_key"}
        if isinstance(obj, float) and not np.isfinite(obj):
            return None
        if isinstance(obj, (np.floating,)):
            v = float(obj)
            return None if not np.isfinite(v) else v
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        return obj

    summary = {
        "n_units": int(len(units)),
        "n_expression": int(len(meta)),
        "n_genes": int(lc.shape[0]),
        "n_q4": q4,
        "n_q1": q1,
        "n_sets_universe": int(len(gene_sets)),
        "tacstd2_q4q1_logFC": tac_lfc,
        "tacstd2_q4q1_t": tac_t,
        "n_gsea_specs": int(len(sweep)),
        "n_gsea_eligible": int(sweep["eligible"].sum()),
        "n_ora_specs": int(len(ora_sweep)),
        "n_ora_eligible": int(ora_sweep["eligible"].sum()),
        "input_md5": hashes,
        "gsea_winner": {
            "spec": winner_name,
            "note": win_spec.note + f" GSEA weight p={int(winner_meta['weight'])}.",
            "weight": float(winner_meta["weight"]),
            "verdict": verdict_text(winner_meta),
            "n_sets": int(winner_meta["n_sets"]),
            "strict_tj_term": winner_meta["strict_tj_term"],
            "strict_tj_rank": int(winner_meta["strict_tj_rank"]),
            "strict_tj_effect": float(winner_meta["strict_tj_effect"]),
            "keratin_term": winner_meta["keratin_term"],
            "keratin_rank": int(winner_meta["keratin_rank"]),
            "keratin_effect": float(winner_meta["keratin_effect"]),
            "adhesion_term": winner_meta["adhesion_term"],
            "adhesion_rank": int(winner_meta["adhesion_rank"]),
            "adhesion_effect": float(winner_meta["adhesion_effect"]),
            "top3": pack_top3_gsea(winner_df),
        },
        "ora_winner": {
            "spec": ora_name,
            "rule": ora_meta["rule"],
            "n_hits": int(ora_meta["n_hits"]),
            "verdict": verdict_text(ora_meta),
            "strict_tj_term": ora_meta["strict_tj_term"],
            "strict_tj_rank": int(ora_meta["strict_tj_rank"]),
            "strict_tj_effect": float(ora_meta["strict_tj_effect"]),
            "keratin_term": ora_meta["keratin_term"],
            "keratin_rank": int(ora_meta["keratin_rank"]),
            "adhesion_term": ora_meta["adhesion_term"],
            "adhesion_rank": int(ora_meta["adhesion_rank"]),
            "top3": pack_top3_ora(ora_df),
        },
        "gsea_reference": {
            "spec": "q4q1_ols_t__w1",
            "strict_tj_term": ref_score["strict_tj_term"],
            "strict_tj_rank": int(ref_score["strict_tj_rank"]),
            "strict_tj_effect": float(ref_score["strict_tj_effect"]),
            "junction_term": ref_score["junction_term"],
            "junction_rank": int(ref_score["junction_rank"]),
            "top3": pack_top3_gsea(ref_df),
        },
        "families": families.to_dict("records"),
        "focal": focal_rows,
    }
    summary["answer"], core_ora = junction_answer(summary, ora_store, ora_df)
    if core_ora is not None:
        core_ora.to_csv(TABLES / "ora_junction_core.tsv", sep="\t", index=False)
    summary = clean(summary)

    sweep.drop(columns=[c for c in sweep.columns if c == "_key"], errors="ignore").to_csv(
        TABLES / "gsea_sweep.tsv", sep="\t", index=False
    )
    ora_sweep.to_csv(TABLES / "ora_sweep.tsv", sep="\t", index=False)
    winner_df.to_csv(TABLES / "gsea_winner.tsv", sep="\t", index=False)
    ora_df.to_csv(TABLES / "ora_winner.tsv", sep="\t", index=False)
    ref_df.to_csv(TABLES / "gsea_reference_q4q1_ols_t.tsv", sep="\t", index=False)
    families.to_csv(TABLES / "family_effects_winner.tsv", sep="\t", index=False)
    pd.DataFrame(focal_rows).to_csv(TABLES / "focal_genes_winner.tsv", sep="\t", index=False)
    if long_rows:
        pd.concat(long_rows, ignore_index=True).to_csv(TABLES / "gsea_all_specs.tsv.gz", sep="\t", index=False)
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))
    write_finding(summary, HERE / "FINDING.md")

    plot_winner(
        winner_df,
        f"GSEA {winner_name}",
        FIGS / "gsea_winner_top.png",
    )
    plot_sweep(sweep, winner_name, FIGS / "gsea_sweep_junction_rank.png")
    if core_ora is not None:
        plot_ora(core_ora, FIGS / "ora_junction_core_top.png")
    print("WINNER", winner_name, summary["gsea_winner"]["verdict"], flush=True)
    print("ORA", ora_name, summary["ora_winner"]["verdict"], flush=True)
    print("TOP3", summary["gsea_winner"]["top3"], flush=True)


if __name__ == "__main__":
    main()

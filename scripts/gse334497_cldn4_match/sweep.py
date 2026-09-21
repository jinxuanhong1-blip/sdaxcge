#!/usr/bin/env python3
"""Exploratory sweep on GSE334497. Trop2 KO, not a CLDN4 knockdown.

The pre-specified 5-vs-5 call lives in analyze.py. This script searches
DE methods, gene sets, sample filters, and Cldn4 contrasts for the
strongest thesis-aligned signals:

  barrier / epithelial program down, IFN / APM / immune / STING up, NHEJ down.

Every row names its contrast, filter, and method. Minimum p across the
grid is not a confirmatory p. BH-FDR is reported inside the grid.
"""
from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "gse334497_cldn4_match"
DATA = OUT / "data"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

KO = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
WT = ["RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R", "control170"]

# Thesis direction for delta = groupA − groupB, with A = KO or Cldn4-low.
THESIS_DOWN = {
    "barrier",
    "epithelial_program",
    "nhej",
    "cldn4_gene",
}
THESIS_UP = {
    "ifn_isg",
    "ifn_alpha",
    "ifn_gamma",
    "apm",
    "immune",
    "inflammatory",
    "sting",
}

BARRIER_STRUCT = [
    "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Cldn15", "Cldn12", "Ocln",
    "Tjp1", "Tjp2", "Tjp3", "F11r", "Cdh1", "Epcam", "Marveld2", "Cgn",
    "Jam2", "Jam3", "Pard3", "Crb3", "Ctnna1", "Ctnnd1",
]
APM = [
    "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q6", "H2-Q7", "H2-T23",
    "Tap1", "Tap2", "Tapbp", "Tapbpl", "Psmb8", "Psmb9", "Psmb10",
    "Psme1", "Psme2", "Nlrc5", "Erap1", "Irf1", "Calr", "Canx", "Pdia3",
]
ISG = [
    "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas2", "Oas3",
    "Oasl1", "Rsad2", "Stat1", "Irf7", "Ifi27", "Usp18", "Bst2", "Ifih1",
]
IMMUNE = [
    "Cd3e", "Cd3d", "Cd8a", "Cd8b1", "Cd4", "Gzmb", "Gzma", "Prf1",
    "Ifng", "Nkg7", "Klrd1", "Ncr1", "Cxcl9", "Cxcl10", "Ptprc",
]
# Sensors and adapters, not the downstream ISG cloud.
STING = [
    "Cgas", "Sting1", "Tbk1", "Irf3", "Ifnb1", "Ifnar1", "Ifnar2",
    "Ddx41", "Zbp1", "Aim2", "Mavs", "Ddx58", "Ifih1", "Trex1", "Enpp1",
    "Ripk3", "Mlkl",
]
NHEJ = [
    "Xrcc4", "Xrcc5", "Xrcc6", "Lig4", "Prkdc", "Nhej1", "Paxx",
    "Dclre1c", "Poll", "Polm", "Trp53bp1", "Rnf8", "Rnf168", "Mad2l2",
    "Xrcc1", "Polb", "Apex1",
]
FOCAL_GENES = [
    "Tacstd2", "Cldn4", "Cldn1", "Cldn7", "Ocln", "Epcam", "Cxcl9",
    "Cd8a", "Prf1", "B2m", "Tap1", "Psmb8", "Nlrc5", "Stat1", "Isg15",
    "Cgas", "Sting1", "Tbk1", "Irf3", "Ifnb1", "Xrcc5", "Xrcc6", "Lig4",
    "Prkdc", "Nhej1",
]


def load_log2() -> pd.DataFrame:
    raw = pd.read_csv(DATA / "GSE334497_normalized_counts.csv.gz", index_col=0)
    sym = {}
    with open(DATA / "ensembl_to_symbol.tsv") as f:
        next(f)
        for line in f:
            ens, s = line.rstrip("\n").split("\t")
            sym[ens] = s
    raw = raw.copy()
    raw["symbol"] = pd.Series(raw.index.astype(str), index=raw.index).map(sym)
    raw = raw.dropna(subset=["symbol"])
    cols = KO + WT
    raw["_mean"] = raw[cols].mean(axis=1)
    raw = raw.sort_values("_mean", ascending=False)
    raw = raw.loc[~raw["symbol"].duplicated()].set_index("symbol")
    mat = raw[cols].astype(float)
    return np.log2(mat + 1.0)


def catalog() -> list[dict]:
    with open(DATA / "sweep_sets_humanfold.json") as f:
        shipped = json.load(f)
    specs = [
        ("BARRIER_STRUCT", BARRIER_STRUCT, "barrier"),
        ("KEGG_TIGHT_JUNCTION", shipped["KEGG_TIGHT_JUNCTION"], "barrier"),
        ("GOBP_TJ_ORGANIZATION", shipped["GOBP_TIGHT_JUNCTION_ORGANIZATION"], "barrier"),
        ("GOBP_TJ_ASSEMBLY", shipped["GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY"], "barrier"),
        ("HALLMARK_APICAL_JUNCTION", shipped["HALLMARK_APICAL_JUNCTION"], "barrier"),
        ("GOBP_KERATINIZATION", shipped["GOBP_KERATINIZATION"], "epithelial_program"),
        ("GOBP_SKIN_BARRIER", shipped["GOBP_ESTABLISHMENT_OF_SKIN_BARRIER"], "epithelial_program"),
        ("KRT_EPITHELIAL", shipped["KRT_EPITHELIAL"], "epithelial_program"),
        ("EPITHELIAL_ISG", ISG, "ifn_isg"),
        ("HALLMARK_IFNA", shipped["HALLMARK_INTERFERON_ALPHA_RESPONSE"], "ifn_alpha"),
        ("HALLMARK_IFNG", shipped["HALLMARK_INTERFERON_GAMMA_RESPONSE"], "ifn_gamma"),
        ("APM_MHCI", APM, "apm"),
        ("BULK_IMMUNE", IMMUNE, "immune"),
        ("HALLMARK_ALLOGRAFT", shipped["HALLMARK_ALLOGRAFT_REJECTION"], "immune"),
        ("HALLMARK_INFLAMMATORY", shipped["HALLMARK_INFLAMMATORY_RESPONSE"], "inflammatory"),
        ("HALLMARK_IL6_STAT3", shipped["HALLMARK_IL6_JAK_STAT3_SIGNALING"], "inflammatory"),
        ("STING_CORE", STING, "sting"),
        ("NHEJ_CORE", NHEJ, "nhej"),
        ("HALLMARK_DNA_REPAIR", shipped["HALLMARK_DNA_REPAIR"], "nhej"),
        ("HALLMARK_OXPHOS", shipped["HALLMARK_OXIDATIVE_PHOSPHORYLATION"], "control"),
    ]
    return [{"name": n, "genes": g, "arm": a} for n, g, a in specs]


def ebayes_t(diff, s2, df_pool, n1, n2) -> tuple[np.ndarray, float]:
    """Limma-style moderated t. Returns t and posterior df."""
    ok = np.isfinite(s2) & (s2 > 0)
    eg = np.log(s2[ok]) - digamma(df_pool / 2) + np.log(df_pool / 2)
    e = float(np.mean(eg))
    v = float(np.var(eg, ddof=1))
    lo, hi = 1e-4, 1e6
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if polygamma(1, mid / 2) > v:
            lo = mid
        else:
            hi = mid
    d0 = 0.5 * (lo + hi)
    s02 = math.exp(e + digamma(d0 / 2) - math.log(d0 / 2))
    d_post = df_pool + d0
    s2_post = (df_pool * s2 + d0 * s02) / d_post
    se = np.sqrt(s2_post * (1 / n1 + 1 / n2))
    t = diff / np.where(se > 0, se, np.nan)
    return t, float(d_post)


def group_effects(log2: pd.DataFrame, a_cols: list[str], b_cols: list[str]) -> dict:
    A = log2[a_cols].to_numpy(dtype=float)
    B = log2[b_cols].to_numpy(dtype=float)
    n1, n2 = A.shape[1], B.shape[1]
    m1 = A.mean(axis=1)
    m2 = B.mean(axis=1)
    v1 = A.var(axis=1, ddof=1)
    v2 = B.var(axis=1, ddof=1)
    diff = m1 - m2
    se_w = np.sqrt(v1 / n1 + v2 / n2)
    with np.errstate(divide="ignore", invalid="ignore"):
        tw = diff / se_w
        df_w = (v1 / n1 + v2 / n2) ** 2 / (
            (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
        )
    df_pool = n1 + n2 - 2
    s2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / df_pool
    tm, d_post = ebayes_t(diff, s2, df_pool, n1, n2)
    return {
        "diff": diff,
        "welch_t": tw,
        "df_w": df_w,
        "mod_t": tm,
        "df_mod": d_post,
        "n1": n1,
        "n2": n2,
        "index": log2.index.to_numpy(),
    }


def welch_p_gene(eff: dict, gene: str, alternative: str) -> float:
    i = int(np.where(eff["index"] == gene)[0][0])
    t = float(eff["welch_t"][i])
    df = float(eff["df_w"][i])
    if not np.isfinite(t):
        return 1.0
    if alternative == "less":
        return float(stats.t.cdf(t, df))
    return float(stats.t.sf(t, df))


def mod_p_gene(eff: dict, gene: str, alternative: str) -> float:
    i = int(np.where(eff["index"] == gene)[0][0])
    t = float(eff["mod_t"][i])
    if not np.isfinite(t):
        return 1.0
    if alternative == "less":
        return float(stats.t.cdf(t, eff["df_mod"]))
    return float(stats.t.sf(t, eff["df_mod"]))


def mwu_p_gene(log2: pd.DataFrame, gene: str, a_cols, b_cols, alternative: str) -> float:
    a = log2.loc[gene, a_cols].to_numpy()
    b = log2.loc[gene, b_cols].to_numpy()
    try:
        _u, p = stats.mannwhitneyu(a, b, alternative=alternative)
    except ValueError:
        return 1.0
    return float(p)


def exact_perm_p(scores: np.ndarray, k: int, observed: float, alternative: str) -> tuple[float, int]:
    n = len(scores)
    n_ext = 0
    n_tot = 0
    for comb in combinations(range(n), k):
        mask = np.zeros(n, dtype=bool)
        mask[list(comb)] = True
        delta = float(scores[mask].mean() - scores[~mask].mean())
        n_tot += 1
        if alternative == "greater" and delta >= observed - 1e-12:
            n_ext += 1
        elif alternative == "less" and delta <= observed + 1e-12:
            n_ext += 1
    return n_ext / n_tot, n_tot


def sample_scores(log2: pd.DataFrame, genes: list[str], cols: list[str], how: str) -> np.ndarray | None:
    present = [g for g in genes if g in log2.index]
    if len(present) < 5:
        return None
    sub = log2.loc[present, cols]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)
    if how == "mean":
        return z.mean(axis=0).to_numpy()
    return z.median(axis=0).to_numpy()


def gene_rank_p(stat: np.ndarray, index: np.ndarray, genes: list[str], alternative: str) -> tuple[float, int]:
    present = [g for g in genes if g in set(index)]
    if len(present) < 5:
        return float("nan"), len(present)
    inset = np.isin(index, present)
    a = stat[inset]
    b = stat[~inset]
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return float("nan"), len(present)
    _u, p = stats.mannwhitneyu(a, b, alternative=alternative)
    return float(p), int(inset.sum())


def thesis_alt(arm: str) -> str | None:
    if arm in THESIS_DOWN:
        return "less"
    if arm in THESIS_UP:
        return "greater"
    return None


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 0.001:
        return f"{p:.3g}"
    return f"{p:.3f}"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    log2 = load_log2()
    sets = catalog()
    cldn4 = log2.loc["Cldn4"]

    filters = [("all10", list(KO), list(WT))]
    for s in KO + WT:
        filters.append((f"drop_{s}", [x for x in KO if x != s], [x for x in WT if x != s]))
    ko_hi = str(cldn4[KO].idxmax())
    wt_lo = str(cldn4[WT].idxmin())
    filters.append(
        (
            f"drop_{ko_hi}_and_{wt_lo}_Cldn4gap",
            [x for x in KO if x != ko_hi],
            [x for x in WT if x != wt_lo],
        )
    )

    rows: list[dict] = []
    gene_rows: list[dict] = []

    def add_set_row(**kwargs):
        rows.append(kwargs)

    # --- Genotype contrasts under each sample filter ---
    for fname, a_cols, b_cols in filters:
        cols = a_cols + b_cols
        eff = group_effects(log2, a_cols, b_cols)
        # Focal genes
        for gene in FOCAL_GENES:
            if gene not in log2.index:
                continue
            i = int(np.where(eff["index"] == gene)[0][0])
            diff = float(eff["diff"][i])
            gene_rows.append(
                {
                    "contrast": "Trop2KO_minus_WT",
                    "filter": fname,
                    "gene": gene,
                    "n_A": len(a_cols),
                    "n_B": len(b_cols),
                    "log2FC": diff,
                    "welch_p_down": welch_p_gene(eff, gene, "less"),
                    "welch_p_up": welch_p_gene(eff, gene, "greater"),
                    "modt_p_down": mod_p_gene(eff, gene, "less"),
                    "modt_p_up": mod_p_gene(eff, gene, "greater"),
                    "mwu_p_down": mwu_p_gene(log2, gene, a_cols, b_cols, "less"),
                    "mwu_p_up": mwu_p_gene(log2, gene, a_cols, b_cols, "greater"),
                }
            )
        for spec in sets:
            genes = spec["genes"]
            arm = spec["arm"]
            alt = thesis_alt(arm)
            # Control is scored in both directions but is not a thesis arm.
            alts = [alt] if alt else ["greater", "less"]
            present = [g for g in genes if g in log2.index]
            for alternative in alts:
                # sample mean / median z, exact perm
                for how, method in (("mean", "sample_mean_z_perm"), ("median", "sample_median_z_perm")):
                    sc = sample_scores(log2, genes, cols, how)
                    if sc is None:
                        continue
                    observed = float(sc[: len(a_cols)].mean() - sc[len(a_cols) :].mean())
                    p, nperm = exact_perm_p(sc, len(a_cols), observed, alternative)
                    add_set_row(
                        contrast="Trop2KO_minus_WT",
                        filter=fname,
                        set=spec["name"],
                        arm=arm,
                        method=method,
                        n_A=len(a_cols),
                        n_B=len(b_cols),
                        n_genes=len(present),
                        effect=observed,
                        alternative=alternative,
                        p=p,
                        n_perm=nperm,
                        thesis=(arm != "control" and alternative == alt),
                    )
                # gene-rank competitive tests (no sample perm)
                for stat, method in (
                    (eff["welch_t"], "gene_rank_welch"),
                    (eff["mod_t"], "gene_rank_modt"),
                ):
                    p, n_g = gene_rank_p(stat, eff["index"], genes, alternative)
                    if p != p:
                        continue
                    inset = np.isin(eff["index"], present)
                    effect = float(np.nanmean(eff["diff"][inset]))
                    add_set_row(
                        contrast="Trop2KO_minus_WT",
                        filter=fname,
                        set=spec["name"],
                        arm=arm,
                        method=method,
                        n_A=len(a_cols),
                        n_B=len(b_cols),
                        n_genes=n_g,
                        effect=effect,
                        alternative=alternative,
                        p=p,
                        n_perm=0,
                        thesis=(arm != "control" and alternative == alt),
                    )

    # --- Cldn4 low vs high, all 10 samples. Cldn4 removed from every set. ---
    order = list(cldn4.sort_values().index)  # low to high
    # median split: 5 lowest vs 5 highest
    cldn_contrasts = [
        ("Cldn4low_minus_Cldn4high_median", order[:5], order[5:]),
        ("Cldn4low_minus_Cldn4high_top3", order[:3], order[-3:]),
    ]
    for cname, a_cols, b_cols in cldn_contrasts:
        n_ko_in_low = sum(s in KO for s in a_cols)
        eff = group_effects(log2, a_cols, b_cols)
        for spec in sets:
            genes = [g for g in spec["genes"] if g != "Cldn4"]
            arm = spec["arm"]
            alt = thesis_alt(arm)
            alts = [alt] if alt else ["greater", "less"]
            present = [g for g in genes if g in log2.index]
            for alternative in alts:
                for how, method in (("mean", "sample_mean_z_perm"), ("median", "sample_median_z_perm")):
                    sc = sample_scores(log2, genes, a_cols + b_cols, how)
                    if sc is None:
                        continue
                    observed = float(sc[: len(a_cols)].mean() - sc[len(a_cols) :].mean())
                    p, nperm = exact_perm_p(sc, len(a_cols), observed, alternative)
                    add_set_row(
                        contrast=cname,
                        filter="all10",
                        set=spec["name"],
                        arm=arm,
                        method=method,
                        n_A=len(a_cols),
                        n_B=len(b_cols),
                        n_genes=len(present),
                        effect=observed,
                        alternative=alternative,
                        p=p,
                        n_perm=nperm,
                        thesis=(arm != "control" and alternative == alt),
                        n_KO_in_low=n_ko_in_low,
                    )
                for stat, method in (
                    (eff["welch_t"], "gene_rank_welch"),
                    (eff["mod_t"], "gene_rank_modt"),
                ):
                    p, n_g = gene_rank_p(stat, eff["index"], genes, alternative)
                    if p != p:
                        continue
                    inset = np.isin(eff["index"], present)
                    effect = float(np.nanmean(eff["diff"][inset]))
                    add_set_row(
                        contrast=cname,
                        filter="all10",
                        set=spec["name"],
                        arm=arm,
                        method=method,
                        n_A=len(a_cols),
                        n_B=len(b_cols),
                        n_genes=n_g,
                        effect=effect,
                        alternative=alternative,
                        p=p,
                        n_perm=0,
                        thesis=(arm != "control" and alternative == alt),
                        n_KO_in_low=n_ko_in_low,
                    )

    # Spearman of Cldn4 vs sample score, Cldn4 removed from the set.
    cols = KO + WT
    cldn_vec = log2.loc["Cldn4", cols].to_numpy()
    for spec in sets:
        genes = [g for g in spec["genes"] if g != "Cldn4"]
        arm = spec["arm"]
        alt = thesis_alt(arm)
        sc = sample_scores(log2, genes, cols, "mean")
        if sc is None or alt is None:
            continue
        rho, p_two = stats.spearmanr(cldn_vec, sc)
        rho = float(rho)
        p_two = float(p_two)
        # Thesis: barrier/NHEJ move with Cldn4 (rho>0); IFN/immune move against Cldn4 (rho<0).
        want_pos = arm in THESIS_DOWN
        if want_pos:
            p = p_two / 2 if rho > 0 else 1 - p_two / 2
            alternative = "rho_positive"
        else:
            p = p_two / 2 if rho < 0 else 1 - p_two / 2
            alternative = "rho_negative"
        add_set_row(
            contrast="Spearman_Cldn4",
            filter="all10",
            set=spec["name"],
            arm=arm,
            method="spearman_score",
            n_A=10,
            n_B=0,
            n_genes=sum(g in log2.index for g in genes),
            effect=rho,
            alternative=alternative,
            p=p,
            n_perm=0,
            thesis=True,
        )

    sweep = pd.DataFrame(rows)
    genes = pd.DataFrame(gene_rows)

    # BH within thesis-directed tests, two scopes.
    sweep["fdr_sweep"] = np.nan
    thesis_mask = sweep["thesis"] == True  # noqa: E712
    if thesis_mask.any():
        sweep.loc[thesis_mask, "fdr_sweep"] = multipletests(
            sweep.loc[thesis_mask, "p"], method="fdr_bh"
        )[1]
    grid = thesis_mask & (sweep["contrast"] == "Trop2KO_minus_WT") & (sweep["filter"] == "all10")
    sweep["fdr_full10"] = np.nan
    if grid.any():
        sweep.loc[grid, "fdr_full10"] = multipletests(
            sweep.loc[grid, "p"], method="fdr_bh"
        )[1]

    sweep = sweep.sort_values(["thesis", "p"], ascending=[False, True])
    sweep.to_csv(TABLES / "sweep_all.tsv", sep="\t", index=False)
    genes.to_csv(TABLES / "sweep_genes.tsv", sep="\t", index=False)

    best = best_table(sweep)
    best.to_csv(TABLES / "sweep_best.tsv", sep="\t", index=False)
    joint = joint_table(sweep)
    joint.to_csv(TABLES / "sweep_joint.tsv", sep="\t", index=False)

    plot_grid(sweep, log2)
    plot_best(sweep, log2, best)
    write_sweep_md(sweep, genes, best, joint, ko_hi, wt_lo)
    print(f"sweep rows {len(sweep)} thesis {int(thesis_mask.sum())}")
    print(best.to_string(index=False))
    print("--- joint head ---")
    print(joint.head(12).to_string(index=False))


def best_table(sweep: pd.DataFrame) -> pd.DataFrame:
    """Best thesis-aligned row per arm, full 10-sample genotype, and across filters."""
    out = []
    sub = sweep[(sweep["thesis"] == True) & (sweep["contrast"] == "Trop2KO_minus_WT")]  # noqa: E712
    for arm, g in sub.groupby("arm"):
        full = g[g["filter"] == "all10"]
        for scope, block in (("full10", full), ("any_filter", g)):
            if block.empty:
                continue
            r = block.sort_values("p").iloc[0]
            out.append(
                {
                    "scope": scope,
                    "arm": arm,
                    "set": r["set"],
                    "method": r["method"],
                    "filter": r["filter"],
                    "contrast": r["contrast"],
                    "n_A": r["n_A"],
                    "n_B": r["n_B"],
                    "n_genes": r["n_genes"],
                    "effect": r["effect"],
                    "p": r["p"],
                    "fdr_full10": r["fdr_full10"],
                    "fdr_sweep": r["fdr_sweep"],
                    "alternative": r["alternative"],
                }
            )
    # Also the best Cldn4-split and Spearman row per arm.
    other = sweep[(sweep["thesis"] == True) & (sweep["contrast"] != "Trop2KO_minus_WT")]  # noqa: E712
    for arm, g in other.groupby("arm"):
        r = g.sort_values("p").iloc[0]
        out.append(
            {
                "scope": "cldn4_contrast",
                "arm": arm,
                "set": r["set"],
                "method": r["method"],
                "filter": r["filter"],
                "contrast": r["contrast"],
                "n_A": r["n_A"],
                "n_B": r["n_B"],
                "n_genes": r["n_genes"],
                "effect": r["effect"],
                "p": r["p"],
                "fdr_full10": r["fdr_full10"],
                "fdr_sweep": r["fdr_sweep"],
                "alternative": r["alternative"],
            }
        )
    return pd.DataFrame(out).sort_values(["scope", "p"])


def joint_table(sweep: pd.DataFrame) -> pd.DataFrame:
    """Same contrast, filter, and method: barrier down AND an immune-side arm up."""
    sub = sweep[(sweep["thesis"] == True) & (sweep["contrast"] == "Trop2KO_minus_WT")].copy()  # noqa: E712
    keys = ["filter", "method", "contrast"]
    barrier_sets = ["BARRIER_STRUCT", "KEGG_TIGHT_JUNCTION", "GOBP_TJ_ASSEMBLY"]
    up_sets = {
        "APM_MHCI": "apm",
        "BULK_IMMUNE": "immune",
        "HALLMARK_IFNG": "ifn_gamma",
        "HALLMARK_IFNA": "ifn_alpha",
        "EPITHELIAL_ISG": "ifn_isg",
        "STING_CORE": "sting",
        "HALLMARK_ALLOGRAFT": "immune",
    }
    rows = []
    for key, g in sub.groupby(keys):
        b = g[g["set"].isin(barrier_sets)].sort_values("p")
        if b.empty:
            continue
        b0 = b.iloc[0]
        ups = []
        for sname in up_sets:
            hit = g[g["set"] == sname]
            if hit.empty:
                continue
            ups.append(hit.sort_values("p").iloc[0])
        if not ups:
            continue
        up = min(ups, key=lambda r: r["p"])
        both = float(b0["p"]) < 0.05 and float(up["p"]) < 0.05
        rows.append(
            {
                "filter": key[0],
                "method": key[1],
                "barrier_set": b0["set"],
                "barrier_effect": b0["effect"],
                "barrier_p": b0["p"],
                "up_set": up["set"],
                "up_arm": up["arm"],
                "up_effect": up["effect"],
                "up_p": up["p"],
                "both_p_lt_0.05": both,
                "max_p": max(float(b0["p"]), float(up["p"])),
                "n_A": b0["n_A"],
                "n_B": b0["n_B"],
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["both_p_lt_0.05", "max_p"], ascending=[False, True])


def plot_grid(sweep: pd.DataFrame, log2: pd.DataFrame) -> None:
    sub = sweep[
        (sweep["contrast"] == "Trop2KO_minus_WT")
        & (sweep["filter"] == "all10")
        & (sweep["arm"] != "control")
        & (sweep["thesis"] == True)  # noqa: E712
    ].copy()
    methods = ["sample_mean_z_perm", "sample_median_z_perm", "gene_rank_welch", "gene_rank_modt"]
    sets = list(dict.fromkeys(sub["set"]))
    mat = np.full((len(sets), len(methods)), np.nan)
    for i, s in enumerate(sets):
        for j, m in enumerate(methods):
            hit = sub[(sub["set"] == s) & (sub["method"] == m)]
            if not hit.empty:
                mat[i, j] = -np.log10(max(float(hit.iloc[0]["p"]), 1e-6))
    fig, ax = plt.subplots(figsize=(8.4, 7.2), dpi=160)
    im = ax.imshow(mat, cmap="YlOrRd", aspect="auto", vmin=0, vmax=4)
    ax.set_xticks(range(len(methods)), ["sample mean z", "sample median z", "gene-rank Welch", "gene-rank mod-t"], rotation=20, ha="right")
    ax.set_yticks(range(len(sets)), sets, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                p = 10 ** (-mat[i, j])
                ax.text(j, i, f"{p:.3g}" if p < 0.001 else f"{p:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if mat[i, j] > 2.2 else "black")
    ax.set_title("Trop2 KO vs WT, all 10 tumors\n−log10 p in the thesis direction (exploratory)")
    fig.colorbar(im, ax=ax, fraction=0.03, label="−log10 p")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_sweep_grid.png", bbox_inches="tight")
    fig.savefig(FIGS / "fig_sweep_grid.pdf", bbox_inches="tight")
    plt.close(fig)


def _row(sweep: pd.DataFrame, set_name: str, method: str) -> pd.Series | None:
    hit = sweep[
        (sweep["contrast"] == "Trop2KO_minus_WT")
        & (sweep["filter"] == "all10")
        & (sweep["set"] == set_name)
        & (sweep["method"] == method)
        & (sweep["thesis"] == True)  # noqa: E712
    ]
    if hit.empty:
        return None
    return hit.iloc[0]


def plot_best(sweep: pd.DataFrame, log2: pd.DataFrame, best: pd.DataFrame) -> None:
    """Sample scores whose p matches the sample permutation, not a gene-rank p."""
    panels = [
        ("KEGG_TIGHT_JUNCTION", "KEGG tight junction\nbarrier down"),
        ("BARRIER_STRUCT", "Claudin / TJ core\nbarrier down"),
        ("GOBP_KERATINIZATION", "Keratinization\nepithelial program down"),
        ("BULK_IMMUNE", "Bulk immune\nup"),
        ("STING_CORE", "STING core\nup"),
        ("NHEJ_CORE", "NHEJ core\ndown"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.6), dpi=160)
    rng = np.random.default_rng(1)
    for ax, (set_name, title) in zip(axes.ravel(), panels):
        r = _row(sweep, set_name, "sample_mean_z_perm")
        sc = sample_scores(log2, _genes_for(set_name), KO + WT, "mean")
        if r is None or sc is None:
            ax.set_axis_off()
            continue
        ax.scatter(rng.normal(0, 0.04, 5), sc[5:], c="#4c78a8", s=28, label="WT", zorder=3)
        ax.scatter(rng.normal(1, 0.04, 5), sc[:5], c="#e45756", s=28, label="Trop2 KO", zorder=3)
        ax.hlines(float(sc[5:].mean()), -0.18, 0.18, color="#4c78a8", lw=1.6)
        ax.hlines(float(sc[:5].mean()), 0.82, 1.18, color="#e45756", lw=1.6)
        ax.axhline(0, color="#ddd", lw=0.6)
        ax.set_xticks([0, 1], ["WT", "KO"])
        ax.set_xlim(-0.55, 1.55)
        ax.set_title(f"{title}\nsample perm p={fmt_p(r['p'])}", fontsize=8.5)
        if set_name == "KEGG_TIGHT_JUNCTION":
            ax.legend(frameon=False, fontsize=7)
    fig.suptitle(
        "GSE334497 Trop2 KO, not CLDN4 KD — sample scores on all 10 tumors",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(FIGS / "fig_sweep_best.png", bbox_inches="tight")
    fig.savefig(FIGS / "fig_sweep_best.pdf", bbox_inches="tight")
    plt.close(fig)


def _genes_for(name: str) -> list[str]:
    with open(DATA / "sweep_sets_humanfold.json") as f:
        shipped = json.load(f)
    custom = {
        "BARRIER_STRUCT": BARRIER_STRUCT,
        "APM_MHCI": APM,
        "EPITHELIAL_ISG": ISG,
        "BULK_IMMUNE": IMMUNE,
        "STING_CORE": STING,
        "NHEJ_CORE": NHEJ,
    }
    alias = {
        "GOBP_TJ_ORGANIZATION": "GOBP_TIGHT_JUNCTION_ORGANIZATION",
        "GOBP_TJ_ASSEMBLY": "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
        "GOBP_SKIN_BARRIER": "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER",
        "HALLMARK_IFNA": "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_IFNG": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_ALLOGRAFT": "HALLMARK_ALLOGRAFT_REJECTION",
        "HALLMARK_INFLAMMATORY": "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_IL6_STAT3": "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_OXPHOS": "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    }
    if name in custom:
        return custom[name]
    key = alias.get(name, name)
    return shipped[key]


def write_sweep_md(sweep, genes, best, joint, ko_hi, wt_lo) -> None:
    n = len(sweep)
    n_th = int((sweep["thesis"] == True).sum())  # noqa: E712
    full = best[best["scope"] == "full10"].sort_values("p")
    anyf = best[best["scope"] == "any_filter"].sort_values("p")
    cld = best[best["scope"] == "cldn4_contrast"].sort_values("p")
    both = joint[joint["both_p_lt_0.05"] == True] if len(joint) else joint  # noqa: E712

    def block(df: pd.DataFrame) -> str:
        lines = [
            "| Arm | Set | Method | Filter / contrast | Effect | p | FDR in full 5v5 grid | FDR in whole sweep |",
            "|---|---|---|---|---:|---:|---:|---:|",
        ]
        for _, r in df.iterrows():
            where = r["filter"] if r["scope"] != "cldn4_contrast" else r["contrast"]
            fdr10 = "—" if pd.isna(r["fdr_full10"]) else fmt_p(r["fdr_full10"])
            fdrs = "—" if pd.isna(r["fdr_sweep"]) else fmt_p(r["fdr_sweep"])
            lines.append(
                f"| {r['arm']} | {r['set']} | {r['method']} | {where} | {r['effect']:+.3f} | "
                f"{fmt_p(r['p'])} | {fdr10} | {fdrs} |"
            )
        return "\n".join(lines)

    # Cldn4 gene across methods on all10
    g10 = genes[(genes["filter"] == "all10") & (genes["gene"] == "Cldn4")].iloc[0]
    # best Cldn4 down p across filters and methods
    cldn_long = []
    for _, r in genes[genes["gene"] == "Cldn4"].iterrows():
        for col, method in (
            ("welch_p_down", "welch"),
            ("modt_p_down", "mod_t"),
            ("mwu_p_down", "mwu"),
        ):
            cldn_long.append((r["filter"], method, r["log2FC"], r[col]))
    cldn_long.sort(key=lambda x: x[3])
    best_cldn = cldn_long[0]
    full_cldn = [x for x in cldn_long if x[0] == "all10"]
    full_cldn.sort(key=lambda x: x[3])

    # STING / NHEJ focal genes on all10
    def gene_blurb(name: str, direction: str) -> str:
        r = genes[(genes["filter"] == "all10") & (genes["gene"] == name)]
        if r.empty:
            return f"{name} absent"
        r = r.iloc[0]
        if direction == "up":
            return f"{name} log2FC {r['log2FC']:+.3f}, mod-t up p {fmt_p(r['modt_p_up'])}, Welch up p {fmt_p(r['welch_p_up'])}"
        return f"{name} log2FC {r['log2FC']:+.3f}, mod-t down p {fmt_p(r['modt_p_down'])}, Welch down p {fmt_p(r['welch_p_down'])}"

    joint_lines = [
        "| Filter | Method | Barrier set | Barrier p | Up set | Up p | Both < 0.05 |",
        "|---|---|---|---:|---|---:|---|",
    ]
    show = joint.head(8)
    for _, r in show.iterrows():
        joint_lines.append(
            f"| {r['filter']} | {r['method']} | {r['barrier_set']} | {fmt_p(r['barrier_p'])} | "
            f"{r['up_set']} | {fmt_p(r['up_p'])} | {'yes' if r['both_p_lt_0.05'] else 'no'} |"
        )

    # How many full10 thesis tests pass raw 0.05 and FDR 0.05
    grid = sweep[(sweep["contrast"] == "Trop2KO_minus_WT") & (sweep["filter"] == "all10") & (sweep["thesis"] == True)]  # noqa: E712
    n_grid = len(grid)
    n_raw = int((grid["p"] < 0.05).sum())
    n_fdr = int((grid["fdr_full10"] < 0.05).sum())

    headline_sets = [
        ("Barrier", "KEGG_TIGHT_JUNCTION", "down"),
        ("Claudin/TJ core", "BARRIER_STRUCT", "down"),
        ("Keratinization", "GOBP_KERATINIZATION", "down"),
        ("Bulk immune", "BULK_IMMUNE", "up"),
        ("Allograft / immune", "HALLMARK_ALLOGRAFT", "up"),
        ("Hallmark IFN-γ", "HALLMARK_IFNG", "up"),
        ("Hallmark IFN-α", "HALLMARK_IFNA", "up"),
        ("MHC-I / APM", "APM_MHCI", "up"),
        ("STING core", "STING_CORE", "up"),
        ("Epithelial ISG", "EPITHELIAL_ISG", "up"),
        ("NHEJ core", "NHEJ_CORE", "down"),
        ("DNA repair hallmark", "HALLMARK_DNA_REPAIR", "down"),
    ]

    def _p(set_name: str, method: str) -> str:
        r = _row(sweep, set_name, method)
        if r is None:
            return "NA"
        return fmt_p(float(r["p"]))

    def _effect(set_name: str) -> str:
        r = _row(sweep, set_name, "gene_rank_modt")
        if r is None:
            r = _row(sweep, set_name, "gene_rank_welch")
        if r is None:
            return "NA"
        return f"{float(r['effect']):+.3f}"

    hlines = []
    for label, set_name, direction in headline_sets:
        sp = _p(set_name, "sample_mean_z_perm")
        gp = _p(set_name, "gene_rank_modt")
        # Call from the sample test first; gene-rank can support but does not override a null sample test.
        try:
            spv = float(_row(sweep, set_name, "sample_mean_z_perm")["p"])
            gpv = float(_row(sweep, set_name, "gene_rank_modt")["p"])
        except (TypeError, KeyError):
            spv, gpv = 1.0, 1.0
        if spv < 0.05 and gpv < 0.05:
            call = "both"
        elif spv < 0.05:
            call = "samples"
        elif gpv < 0.05:
            call = "gene-rank only"
        else:
            call = "no"
        hlines.append(
            f"| {label} | {set_name} | {_effect(set_name)} | {sp} | {gp} | {call} |"
        )
    headline = "\n".join(hlines)

    text = f"""# SWEEP — GSE334497 Trop2 KO, thesis-aligned signals

**This is a Trop2 (Tacstd2) knockout, not a CLDN4 knockdown.** Wu *et al.*, JITC 2026, GEO GSE334497, 4T1 tumors, frozen sections. The pre-specified 5-vs-5 match (Cldn4 Welch and epithelial-ISG sample score) is unchanged in FINDING.md. This file is the exploratory grid run after that call.

Delta for genotype tests is KO minus WT. Thesis direction: barrier and NHEJ **down**, IFN / APM / immune / STING **up**. Cldn4-split delta is Cldn4-low minus Cldn4-high, with *Cldn4* removed from every set. Human GMT sets were case-folded to mouse symbols (`Cldn4` from `CLDN4`); genes absent from the matrix are dropped.

## Grid

| Item | Count |
|---|---:|
| Rows in `tables/sweep_all.tsv` | {n} |
| Thesis-directed tests (the BH family) | {n_th} |
| Full 5 vs 5 genotype × thesis tests | {n_grid} |
| Of those, raw p < 0.05 | {n_raw} |
| Of those, BH-FDR < 0.05 inside the 5 vs 5 grid | {n_fdr} |

Methods: exact permutation of the sample mean z-score, exact permutation of the sample median z-score, competitive Mann–Whitney of gene-level Welch *t* versus the rest of the genome, and the same competitive test on limma-style moderated *t*. Filters: all 10 samples, each leave-one-out, and one named pair that widens the *Cldn4* gap (drop KO `{ko_hi}`, the KO with the highest *Cldn4*, and WT `{wt_lo}`, the WT with the lowest *Cldn4*).

BH-FDR in the full 5 vs 5 grid is a descriptor inside that grid. The 76 tests reuse the same 10 tumors and overlapping genes, so they are not 76 independent experiments. Gene-rank p-values treat genes as exchangeable; they are smaller than sample-permutation p-values because the null is not “these 10 tumors could have swapped labels.” Sample permutation is the test that matches n = 5 vs 5. Leave-one-out minima are a search, reported below the full-matrix table.

## Headline on all 10 tumors (KO − WT)

| Program | Set | Mean log2FC | Sample mean-z perm p | Gene-rank mod-t p | Call |
|---|---|---:|---:|---:|---|
{headline}

“Both” means raw *p* < 0.05 on the sample permutation **and** on the gene-rank test. It is not a claim that the within-grid BH-FDR is < 0.05. STING is the borderline case: sample perm *p* = 0.048, gene-rank Welch *p* = 0.034, and the BH-FDR inside the 76-test grid is **0.059–0.072**. Its mean log2FC is only **+0.036**. *Cgas* moves the other way (log2FC −0.381). *Sting1* (Welch up *p* = 0.035) and *Irf3* (Welch up *p* = 0.009) are the STING genes that rise. *Ifnb1* is absent from the matrix. NHEJ core does not fall (mean log2FC +0.032). Epithelial ISG does not rise. *Cldn4* as one gene does not clear 0.05 on the full 5 vs 5 (moderated-*t* down *p* = 0.094; Welch 0.126). Dropping the KO with the highest *Cldn4* and the WT with the lowest *Cldn4* (`KO165` and `RESUB-170R`) is what brings the *Cldn4* moderated-*t* down *p* to 0.016. That filter was chosen because it widens the *Cldn4* gap.

On the *Cldn4* median split, 3 of the 5 low-*Cldn4* tumors are Trop2 KO, so that contrast is mixed with genotype. *Cldn4* was removed from the sets before testing.

## Best single method per arm on all 10 tumors

{block(full)}

## Best row per arm after sample filters

{block(anyf)}

## Best *Cldn4*-defined contrast per arm

These contrasts are not the knockout. On the median split, the number of Trop2-KO samples inside the *Cldn4*-low half is in `sweep_all.tsv` (`n_KO_in_low`). *Cldn4* itself is removed from the sets.

{block(cld)}

## Same specification, barrier down and an up-arm together

{chr(10).join(joint_lines)}

## *Cldn4* gene, three tests

On all 10 tumors the smallest down-sided p is **{full_cldn[0][1]} p = {fmt_p(full_cldn[0][3])}** at log2FC {g10['log2FC']:+.3f}. Welch down p = {fmt_p(g10['welch_p_down'])}, moderated-t down p = {fmt_p(g10['modt_p_down'])}, MWU down p = {fmt_p(g10['mwu_p_down'])}.

Across filters and those three tests, the smallest *Cldn4* down p is **{best_cldn[0]} / {best_cldn[1]} / p = {fmt_p(best_cldn[3])}** (log2FC {best_cldn[2]:+.3f}). That row is a sensitivity, not a replacement for the 5 vs 5 test.

## STING and NHEJ genes (all 10, KO − WT)

STING, up side: {gene_blurb('Cgas','up')}; {gene_blurb('Sting1','up')}; {gene_blurb('Tbk1','up')}; {gene_blurb('Irf3','up')}; {gene_blurb('Ifnb1','up')}.

NHEJ, down side: {gene_blurb('Xrcc5','down')}; {gene_blurb('Xrcc6','down')}; {gene_blurb('Lig4','down')}; {gene_blurb('Prkdc','down')}; {gene_blurb('Nhej1','down')}.

Figures: `figures/fig_sweep_grid.png` (full 5 vs 5 grid), `figures/fig_sweep_best.png` (sample scores of the winning full-matrix set in each arm).

Reproduce: `python3 scripts/gse334497_cldn4_match/sweep.py`
"""
    (OUT / "SWEEP.md").write_text(text)
    print("wrote", OUT / "SWEEP.md")


if __name__ == "__main__":
    main()

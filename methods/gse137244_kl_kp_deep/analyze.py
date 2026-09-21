#!/usr/bin/env python3
"""GSE137244 KL vs KP: beyond mean Tacstd2 / Cldn4.

Public cell-line RNA-seq only (Deng et al., Cancer Discovery 2021, PMID 34142094).
Five KrasG12D;Lkb1 (KL) libraries versus five KrasG12D;Trp53 (KP) libraries.
Normal lung is excluded. No private 8-KL matrices. No TISMO LLC.

Primary scale is mean log2(FPKM+1), the scale on which Tacstd2 Δ=+3.24 and
Cldn4 Δ=+5.57 were locked. Pathway tests are preranked GSEA (gene permutation)
plus per-line ssGSEA tested with the exact Mann-Whitney floor for n=5 vs n=5.
"""

from __future__ import annotations

import json
import math
import urllib.request
from collections import Counter
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "cache"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
GENESETS = ROOT / "gene_sets"
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/"
    "GSE137244_counts.fpkm.csv.gz"
)
FPKM_PATH = CACHE / "GSE137244_counts.fpkm.csv.gz"

KP_LIBS = [
    "B6AL10-1-RNA",
    "B6AL10-2-RNA",
    "B6AL10-3-RNA",
    "B6AL10-4-RNA",
    "B6AL10-5-RNA",
]
KL_LIBS = [
    "KL155mix-control-2-RNA",
    "KL47-1-untreated-1-RNA",
    "KLC-RNA",
    "KLD-RNA",
    "KLE-RNA",
]
SHORT = {
    "B6AL10-1-RNA": "KP B6AL10-1",
    "B6AL10-2-RNA": "KP B6AL10-2",
    "B6AL10-3-RNA": "KP B6AL10-3",
    "B6AL10-4-RNA": "KP B6AL10-4",
    "B6AL10-5-RNA": "KP B6AL10-5",
    "KL155mix-control-2-RNA": "KL KL155",
    "KL47-1-untreated-1-RNA": "KL KL47",
    "KLC-RNA": "KL KLC",
    "KLD-RNA": "KL KLD",
    "KLE-RNA": "KL KLE",
}
GSM = {
    "B6AL10-1-RNA": "GSM4073816",
    "B6AL10-2-RNA": "GSM4073817",
    "B6AL10-3-RNA": "GSM4073818",
    "B6AL10-4-RNA": "GSM4073819",
    "B6AL10-5-RNA": "GSM4073820",
    "KL155mix-control-2-RNA": "GSM4073821",
    "KL47-1-untreated-1-RNA": "GSM4073822",
    "KLC-RNA": "GSM4073823",
    "KLD-RNA": "GSM4073824",
    "KLE-RNA": "GSM4073825",
    "normal-lung-RNA": "GSM4073826",
}

# Near-replicate collapse. Single-linkage on Pearson r of log2(FPKM+1)
# among genes with mean FPKM >= 1. Threshold is reported, not hidden.
CORR_THR = 0.98
N_PERM = 1000
SEED = 3033
FLOOR_P = 2 / math.comb(10, 5)  # 2/252

KP_COLOR = "#4C78A8"
KL_COLOR = "#E45756"
TJ_COLOR = "#7B4B94"

# ssGSEA / mean-score sets, in figure order.
SS_ORDER = [
    "EPITHELIAL_TJ_CORE",
    "GO_BICELLULAR_TIGHT_JUNCTION",
    "KEGG_TIGHT_JUNCTION",
    "HALLMARK_APICAL_JUNCTION",
    "GOCC_GAP_JUNCTION",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "COMPACT_ISG",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
]
SS_LABEL = {
    "EPITHELIAL_TJ_CORE": "Epithelial TJ core",
    "GO_BICELLULAR_TIGHT_JUNCTION": "GO bicellular TJ",
    "KEGG_TIGHT_JUNCTION": "KEGG tight junction",
    "HALLMARK_APICAL_JUNCTION": "Hallmark apical junction",
    "GOCC_GAP_JUNCTION": "GO gap junction",
    "GOCC_ADHERENS_JUNCTION": "GO adherens junction",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Hallmark IFN-α",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Hallmark IFN-γ",
    "COMPACT_ISG": "Compact ISG",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "Hallmark EMT",
}

CALLOUT_EXTRA = [
    "Tacstd2",
    "Ocln",
    "Marveld2",
    "Marveld3",
    "Tjp1",
    "Tjp2",
    "Tjp3",
    "F11r",
    "Jam2",
    "Jam3",
    "Cgn",
    "Cgnl1",
    "Crb3",
    "Ildr1",
    "Lsr",
    "Esam",
]


def ensure_fpkm() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    if FPKM_PATH.exists() and FPKM_PATH.stat().st_size > 1000:
        return
    urllib.request.urlretrieve(FPKM_URL, FPKM_PATH)


def load_gmt(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            genes = []
            for gene in parts[2:]:
                if gene and gene not in genes:
                    genes.append(gene)
            out[parts[0]] = genes
    return out


def load_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(FPKM_PATH)
    # Two symbols are duplicated Excel date conversions (1-Mar, 2-Mar).
    n_dup_rows = int(raw["gene"].duplicated().sum())
    expr = raw.groupby("gene", as_index=True).mean(numeric_only=True)
    libs = KP_LIBS + KL_LIBS
    missing = [c for c in libs + ["normal-lung-RNA"] if c not in expr.columns]
    if missing:
        raise SystemExit(f"Missing columns: {missing}")
    meta = pd.DataFrame(
        {
            "library": libs,
            "short": [SHORT[c] for c in libs],
            "gsm": [GSM[c] for c in libs],
            "arm": ["KP"] * 5 + ["KL"] * 5,
            "genotype": ["KrasG12D;Trp53"] * 5 + ["KrasG12D;Lkb1"] * 5,
        }
    )
    meta.attrs["n_duplicate_symbol_rows"] = n_dup_rows
    return expr, meta


def welch_t(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    va = a.var(ddof=1)
    vb = b.var(ddof=1)
    se = math.sqrt(va / len(a) + vb / len(b))
    diff = float(a.mean() - b.mean())
    if se == 0:
        if diff == 0:
            return 0.0
        return math.copysign(math.inf, diff)
    return diff / se


def pairwise_cliff(a: np.ndarray, b: np.ndarray) -> tuple[int, int, int, float]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n_gt = n_lt = n_eq = 0
    for x in a:
        for y in b:
            if x > y:
                n_gt += 1
            elif x < y:
                n_lt += 1
            else:
                n_eq += 1
    denom = len(a) * len(b)
    cliff = (n_gt - n_lt) / denom
    return n_gt, n_lt, n_eq, float(cliff)


def exact_mw(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    res = mannwhitneyu(a, b, alternative="two-sided", method="exact")
    n_gt, n_lt, n_eq, cliff = pairwise_cliff(a, b)
    complete = bool(a.min() > b.max() or b.min() > a.max())
    # The 2/252 floor is the 5-vs-5 library test only.
    hits_floor = bool(
        len(a) == 5 and len(b) == 5 and np.isclose(res.pvalue, FLOOR_P, rtol=0, atol=1e-12)
    )
    return {
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "n_kl_gt_kp": n_gt,
        "n_kl_lt_kp": n_lt,
        "n_tie_pairs": n_eq,
        "cliffs_delta": cliff,
        "complete_separation": complete,
        "hits_floor": hits_floor,
        "delta_mean": float(a.mean() - b.mean()),
        "mean_kl": float(a.mean()),
        "mean_kp": float(b.mean()),
    }


def mw_ladder() -> pd.DataFrame:
    """Exact two-sided MW p for every U under n=5 vs n=5, no ties."""
    seen = {}
    for comb in combinations(range(10), 5):
        other = [i for i in range(10) if i not in comb]
        wins = sum(x > y for x in comb for y in other)
        if wins not in seen:
            seen[wins] = (list(comb), other)
    rows = []
    counts = Counter()
    for comb in combinations(range(10), 5):
        other = [i for i in range(10) if i not in comb]
        wins = sum(x > y for x in comb for y in other)
        counts[wins] += 1
    for u in range(0, 26):
        a, b = seen[u]
        p = float(mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue)
        cliff = (2 * u - 25) / 25
        rows.append(
            {
                "U_kl_wins": u,
                "n_configurations": int(counts[u]),
                "two_sided_exact_p": p,
                "cliffs_delta": cliff,
                "rejects_alpha_0.05": p < 0.05,
                "is_floor": u in (0, 25),
            }
        )
    return pd.DataFrame(rows)


def union_find_clusters(names: list[str], corr: pd.DataFrame, thr: float) -> dict[str, str]:
    parent = {n: n for n in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            if float(corr.loc[a, b]) >= thr:
                union(a, b)
    groups: dict[str, list[str]] = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    # stable labels
    labels = {}
    for members in groups.values():
        members = sorted(members, key=lambda x: names.index(x))
        label = "+".join(SHORT[m].split(" ", 1)[1] for m in members)
        for m in members:
            labels[m] = label
    return labels


def mean_score(log_expr: pd.DataFrame, genes: list[str], cols: list[str]) -> pd.Series:
    present = [g for g in genes if g in log_expr.index]
    if not present:
        return pd.Series(np.nan, index=cols)
    return log_expr.loc[present, cols].mean(axis=0)


def gene_level(log_expr: pd.DataFrame, fpkm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    libs = KP_LIBS + KL_LIBS
    for gene in log_expr.index:
        kl = log_expr.loc[gene, KL_LIBS].to_numpy(dtype=float)
        kp = log_expr.loc[gene, KP_LIBS].to_numpy(dtype=float)
        stats = exact_mw(kl, kp)
        t = welch_t(kl, kp)
        fp = fpkm.loc[gene, libs].to_numpy(dtype=float)
        rows.append(
            {
                "gene": gene,
                "mean_fpkm": float(fp.mean()),
                "max_fpkm": float(fp.max()),
                "n_libs_fpkm_ge1": int((fp >= 1).sum()),
                "mean_fpkm_kl": float(fpkm.loc[gene, KL_LIBS].mean()),
                "mean_fpkm_kp": float(fpkm.loc[gene, KP_LIBS].mean()),
                "delta_log2": stats["delta_mean"],
                "welch_t": t if math.isfinite(t) else np.nan,
                "U": stats["U"],
                "mw_p": stats["p"],
                "cliffs_delta": stats["cliffs_delta"],
                "n_tie_pairs": stats["n_tie_pairs"],
                "complete_separation": stats["complete_separation"],
                "hits_floor": stats["hits_floor"],
            }
        )
    return pd.DataFrame(rows)


def run_prerank(rank: pd.Series, gene_sets: dict[str, list[str]], min_size: int, tag: str) -> pd.DataFrame:
    import gseapy as gp

    rnk = rank.replace([np.inf, -np.inf], np.nan).dropna()
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        organism="human",  # ignored when gene_sets is an explicit dict/gmt of mouse symbols
        min_size=min_size,
        max_size=500,
        permutation_num=N_PERM,
        weight=1.0,
        ascending=False,
        threads=4,
        seed=SEED,
        no_plot=True,
        outdir=None,
        verbose=False,
        pheno_pos="KL",
        pheno_neg="KP",
    )
    df = pre.res2d.copy()
    df.insert(0, "rank_metric", tag)
    df.insert(1, "n_genes_ranked", int(rnk.shape[0]))
    return df


def run_ssgsea(log_expr: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    import gseapy as gp

    data = log_expr[KP_LIBS + KL_LIBS].copy()
    ss = gp.ssgsea(
        data=data,
        gene_sets=gene_sets,
        min_size=10,
        max_size=500,
        sample_norm_method="rank",
        correl_norm_type="rank",
        permutation_num=0,
        weight=0.25,
        threads=4,
        seed=SEED,
        no_plot=True,
        outdir=None,
        verbose=False,
    )
    return ss.res2d.copy()


def contrast_scores(scores: pd.DataFrame, value_col: str, set_col: str, sample_col: str) -> pd.DataFrame:
    rows = []
    for term, sub in scores.groupby(set_col):
        wide = sub.set_index(sample_col)[value_col]
        kl = wide.reindex(KL_LIBS).to_numpy(dtype=float)
        kp = wide.reindex(KP_LIBS).to_numpy(dtype=float)
        stats = exact_mw(kl, kp)
        row = {"set": term, **stats}
        for lib in KP_LIBS + KL_LIBS:
            row[SHORT[lib]] = float(wide.loc[lib])
        rows.append(row)
    return pd.DataFrame(rows)


def collapsed_units(log_expr: pd.DataFrame, labels: dict[str, str]) -> pd.DataFrame:
    """Equal weight per near-replicate cluster: mean of member libraries on log2(FPKM+1)."""
    frame = log_expr[KP_LIBS + KL_LIBS].copy()
    clusters = []
    for arm, libs in (("KP", KP_LIBS), ("KL", KL_LIBS)):
        for label in dict.fromkeys(labels[lib] for lib in libs):
            members = [lib for lib in libs if labels[lib] == label]
            col = frame[members].mean(axis=1)
            col.name = f"{arm}:{label}"
            clusters.append(col)
    return pd.concat(clusters, axis=1)


def style_ax(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)


def plot_volcano(stats: pd.DataFrame, callouts: pd.DataFrame) -> None:
    show = stats.loc[stats["max_fpkm"] >= 1].copy()
    show["nlp"] = -np.log10(show["mw_p"].clip(lower=1e-300))
    floor_y = -math.log10(FLOOR_P)
    step_y = -math.log10(8 / 252)  # U>=23, p=0.031746

    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    background = show.loc[~show["gene"].isin(callouts["gene"])]
    ax.scatter(
        background["delta_log2"],
        background["nlp"],
        s=8,
        c="#C8C8C8",
        linewidths=0,
        alpha=0.45,
        rasterized=True,
        zorder=1,
    )
    on_floor = background.loc[background["hits_floor"]]
    ax.scatter(
        on_floor["delta_log2"],
        on_floor["nlp"],
        s=10,
        c="#E6C35C",
        linewidths=0,
        alpha=0.7,
        zorder=2,
        label="Other genes at MW floor",
    )
    tj = show.loc[show["gene"].isin(callouts["gene"]) & (show["gene"] != "Tacstd2")]
    ax.scatter(
        tj["delta_log2"],
        tj["nlp"],
        s=28,
        c=TJ_COLOR,
        linewidths=0,
        zorder=3,
        label="TJ callout",
    )
    tac = show.loc[show["gene"] == "Tacstd2"]
    ax.scatter(
        tac["delta_log2"],
        tac["nlp"],
        s=42,
        c=KL_COLOR,
        linewidths=0,
        zorder=4,
        label="Tacstd2",
    )
    ax.axhline(floor_y, color="#222222", lw=1.0, zorder=0)
    ax.axhline(step_y, color="#222222", lw=0.8, ls=":", zorder=0)
    ax.axvline(0, color="#888888", lw=0.6, zorder=0)
    ax.text(
        0.99,
        floor_y,
        " exact MW floor  2/252 = 0.00794",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8,
        color="#222222",
    )
    ax.text(
        0.99,
        step_y,
        " largest p still < 0.05  (8/252 = 0.032)",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8,
        color="#444444",
    )
    label_genes = callouts.loc[callouts["label"]].copy()
    texts = []
    for _, row in label_genes.iterrows():
        hit = show.loc[show["gene"] == row["gene"]]
        if hit.empty:
            continue
        x = float(hit["delta_log2"].iloc[0])
        y = float(hit["nlp"].iloc[0])
        texts.append(ax.text(x, y, row["gene"], fontsize=8, color="#222222"))
    if texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="#666666", lw=0.5),
            expand=(1.15, 1.25),
            force_text=(0.4, 0.6),
        )
    ax.set_xlabel("KL − KP   mean log2(FPKM + 1)")
    ax.set_ylabel("−log10 exact Mann–Whitney p")
    ax.set_title("GSE137244 libraries  ·  5 KL vs 5 KP")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(FIGS / "volcano_tj_callouts.png", dpi=160)
    fig.savefig(FIGS / "volcano_tj_callouts.pdf")
    plt.close(fig)


def plot_ssgsea(wide_scores: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    sets = [s for s in SS_ORDER if s in set(wide_scores["set"])]
    n = len(sets)
    ncols = 3
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(9.6, 2.55 * nrows), sharey=False)
    axes = np.atleast_1d(axes).ravel()
    y_pos = {lib: i for i, lib in enumerate(KP_LIBS + KL_LIBS)}
    for ax, term in zip(axes, sets):
        sub = wide_scores.loc[wide_scores["set"] == term].set_index("library")
        con = contrasts.loc[contrasts["set"] == term].iloc[0]
        for lib in KP_LIBS + KL_LIBS:
            arm = "KP" if lib in KP_LIBS else "KL"
            ax.scatter(
                float(sub.loc[lib, "ES"]),
                y_pos[lib],
                s=28,
                c=KP_COLOR if arm == "KP" else KL_COLOR,
                zorder=3,
            )
        ax.axvline(con["mean_kp"], color=KP_COLOR, lw=0.7, ls="--", alpha=0.8)
        ax.axvline(con["mean_kl"], color=KL_COLOR, lw=0.7, ls="--", alpha=0.8)
        p = con["p"]
        flag = "floor" if con["hits_floor"] else f"p={p:.3f}"
        ax.set_title(
            f"{SS_LABEL.get(term, term)}\nΔ={con['delta_mean']:+.3f}  {flag}",
            fontsize=9,
        )
        ax.set_yticks(list(y_pos.values()))
        ax.set_yticklabels([SHORT[lib] for lib in y_pos], fontsize=7)
        ax.set_xlabel("ssGSEA ES", fontsize=8)
        style_ax(ax)
        ax.tick_params(axis="x", labelsize=8)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("ssGSEA enrichment score per library  ·  KL red, KP blue", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(FIGS / "ssgsea_per_line.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGS / "ssgsea_per_line.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_gsea(gsea: pd.DataFrame) -> None:
    """Dumbbell of NES under Welch-t rank vs log2FC rank."""
    keep_terms = [
        "EPITHELIAL_TJ_CORE",
        "GO_BICELLULAR_TIGHT_JUNCTION",
        "KEGG_TIGHT_JUNCTION",
        "GOCC_ADHERENS_JUNCTION",
        "GOCC_GAP_JUNCTION",
        "HALLMARK_APICAL_JUNCTION",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ]
    welch = gsea.loc[gsea["rank_metric"] == "welch_t"].set_index("Term")
    delta = gsea.loc[gsea["rank_metric"] == "delta_log2"].set_index("Term")
    terms = [t for t in keep_terms if t in welch.index and t in delta.index]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    y = np.arange(len(terms))
    for i, term in enumerate(terms):
        x1 = float(welch.loc[term, "NES"])
        x2 = float(delta.loc[term, "NES"])
        ax.plot([x1, x2], [i, i], color="#B0B0B0", lw=1.4, zorder=1)
        ax.scatter([x1], [i], s=46, c="#222222", zorder=2, label="Welch t" if i == 0 else None)
        ax.scatter(
            [x2],
            [i],
            s=46,
            facecolors="none",
            edgecolors="#222222",
            linewidths=1.2,
            zorder=2,
            label="mean log2FC" if i == 0 else None,
        )
        fdr = float(welch.loc[term, "FDR q-val"])
        nom = float(welch.loc[term, "NOM p-val"])
        ax.text(
            max(x1, x2) + 0.06,
            i,
            f"FDR {fdr:.3f}   nom {nom:.3f}",
            va="center",
            fontsize=8,
            color="#333333",
        )
    ax.axvline(0, color="#888888", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([SS_LABEL.get(t, t) for t in terms])
    ax.set_xlabel("Preranked NES   (positive = enriched in KL)")
    ax.set_title("GSEA  ·  filled = Welch t rank, open = mean log2FC rank\nFDR: hallmark sets among 50; junction sets among 5")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style_ax(ax)
    # leave room for FDR labels
    xmax = max(float(welch.loc[t, "NES"]) for t in terms)
    xmin = min(min(float(welch.loc[t, "NES"]), float(delta.loc[t, "NES"])) for t in terms)
    ax.set_xlim(xmin - 0.35, xmax + 1.55)
    fig.tight_layout()
    fig.savefig(FIGS / "gsea_nes.png", dpi=160)
    fig.savefig(FIGS / "gsea_nes.pdf")
    plt.close(fig)


def plot_corr(corr: pd.DataFrame, labels: dict[str, str]) -> None:
    libs = KP_LIBS + KL_LIBS
    mat = corr.loc[libs, libs].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(mat, vmin=0.55, vmax=1.0, cmap="Blues")
    ax.set_xticks(range(10))
    ax.set_yticks(range(10))
    names = [SHORT[c].replace(" ", "\n") for c in libs]
    ax.set_xticklabels(names, fontsize=7)
    ax.set_yticklabels([SHORT[c] for c in libs], fontsize=8)
    for i in range(10):
        for j in range(10):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color="#111111")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r  log2(FPKM+1)")
    ax.set_title(f"Library correlation  ·  clusters join at r ≥ {CORR_THR:.2f}")
    fig.tight_layout()
    fig.savefig(FIGS / "library_correlation.png", dpi=160)
    fig.savefig(FIGS / "library_correlation.pdf")
    plt.close(fig)
    _ = labels


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    ensure_fpkm()
    fpkm, meta = load_expression()
    libs = KP_LIBS + KL_LIBS
    log_expr = np.log2(fpkm[libs + ["normal-lung-RNA"]] + 1)

    hallmark = load_gmt(GENESETS / "mh.all.v2024.1.Mm.symbols.gmt")
    junction = load_gmt(GENESETS / "junction_controls.gmt")
    curated = load_gmt(GENESETS / "epithelial_and_isg.gmt")
    # Hallmark names in the 2024.1 GMT.
    needed_h = [
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_APICAL_JUNCTION",
    ]
    for name in needed_h:
        if name not in hallmark:
            raise SystemExit(f"Missing hallmark set {name}")

    focused = {
        "EPITHELIAL_TJ_CORE": curated["EPITHELIAL_TJ_CORE"],
        "GO_BICELLULAR_TIGHT_JUNCTION": junction["GO_BICELLULAR_TIGHT_JUNCTION"],
        "KEGG_TIGHT_JUNCTION": junction["KEGG_TIGHT_JUNCTION"],
        "GOCC_GAP_JUNCTION": junction["GOCC_GAP_JUNCTION"],
        "GOCC_ADHERENS_JUNCTION": junction["GOCC_ADHERENS_JUNCTION"],
    }
    ssgsea_sets = {
        **focused,
        **{k: hallmark[k] for k in needed_h},
        "COMPACT_ISG": curated["COMPACT_ISG"],
    }

    print("Gene-level exact MW ...", flush=True)
    stats = gene_level(log_expr[libs], fpkm)
    expressed = stats.loc[stats["max_fpkm"] >= 1].copy()
    reject, qvals, _, _ = multipletests(expressed["mw_p"].to_numpy(), method="fdr_bh")
    expressed["bh_q"] = qvals
    expressed["bh_reject_0.05"] = reject
    n_tested = int(len(expressed))
    n_floor = int(expressed["hits_floor"].sum())
    n_complete = int(expressed["complete_separation"].sum())
    k_needed = math.ceil(FLOOR_P * n_tested / 0.05)
    n_bh = int(reject.sum())

    # Locked single genes.
    locked_rows = []
    for gene in ["Tacstd2", "Cldn4", "Stk11", "Trp53"]:
        row = stats.loc[stats["gene"] == gene].iloc[0]
        locked_rows.append(row)
    locked = pd.DataFrame(locked_rows)

    # Callouts: every claudin plus structural TJ genes and Tacstd2.
    callout_genes = set(CALLOUT_EXTRA)
    callout_genes.update(g for g in stats["gene"] if str(g).startswith("Cldn"))
    callouts = stats.loc[stats["gene"].isin(callout_genes)].copy()
    callouts["label"] = (callouts["max_fpkm"] >= 1) & (callouts["delta_log2"].abs() >= 1.0)
    # Always label the two locked genes.
    callouts.loc[callouts["gene"].isin(["Tacstd2", "Cldn4"]), "label"] = True
    callouts = callouts.sort_values("delta_log2", ascending=False)

    print("Library correlation and clusters ...", flush=True)
    detected = log_expr.loc[fpkm[libs].mean(axis=1) >= 1, libs]
    corr = detected.corr()
    labels = {}
    labels.update(union_find_clusters(KP_LIBS, corr, CORR_THR))
    labels.update(union_find_clusters(KL_LIBS, corr, CORR_THR))
    meta["cluster"] = meta["library"].map(labels)
    cluster_expr = collapsed_units(log_expr, labels)
    kl_clusters = [c for c in cluster_expr.columns if c.startswith("KL:")]
    kp_clusters = [c for c in cluster_expr.columns if c.startswith("KP:")]

    # Mean signature scores on the locked scale, library and cluster.
    score_sets = ssgsea_sets
    mean_rows = []
    for name, genes in score_sets.items():
        present = [g for g in genes if g in log_expr.index]
        sc = mean_score(log_expr, present, libs)
        stats_lib = exact_mw(sc.reindex(KL_LIBS), sc.reindex(KP_LIBS))
        # cluster-level: mean of member-library scores
        cl_scores = {}
        for col in cluster_expr.columns:
            arm, label = col.split(":", 1)
            members = [lib for lib in (KL_LIBS if arm == "KL" else KP_LIBS) if labels[lib] == label]
            cl_scores[col] = float(sc.reindex(members).mean())
        stats_cl = exact_mw(
            [cl_scores[c] for c in kl_clusters],
            [cl_scores[c] for c in kp_clusters],
        )
        mean_rows.append(
            {
                "set": name,
                "n_genes_in_matrix": len(present),
                "n_genes_listed": len(genes),
                "library_delta": stats_lib["delta_mean"],
                "library_p": stats_lib["p"],
                "library_U": stats_lib["U"],
                "library_cliffs": stats_lib["cliffs_delta"],
                "library_complete": stats_lib["complete_separation"],
                "library_hits_floor": stats_lib["hits_floor"],
                "cluster_delta": stats_cl["delta_mean"],
                "cluster_p": stats_cl["p"],
                "cluster_U": stats_cl["U"],
                "cluster_cliffs": stats_cl["cliffs_delta"],
                "cluster_complete": stats_cl["complete_separation"],
                "n_kl_clusters": len(kl_clusters),
                "n_kp_clusters": len(kp_clusters),
                **{SHORT[lib]: float(sc.loc[lib]) for lib in libs},
            }
        )
    mean_tbl = pd.DataFrame(mean_rows)

    # Cluster-level locked genes.
    cluster_gene_rows = []
    for gene in ["Tacstd2", "Cldn4", "Stk11", "Trp53"]:
        stats_cl = exact_mw(
            cluster_expr.loc[gene, kl_clusters],
            cluster_expr.loc[gene, kp_clusters],
        )
        cluster_gene_rows.append({"gene": gene, **stats_cl})
    cluster_genes = pd.DataFrame(cluster_gene_rows)

    print("Preranked GSEA ...", flush=True)
    gsea_universe = stats.loc[stats["n_libs_fpkm_ge1"] >= 2].copy()
    welch_rank = gsea_universe.set_index("gene")["welch_t"].astype(float)
    delta_rank = gsea_universe.set_index("gene")["delta_log2"].astype(float)
    # Drop non-finite Welch (identical expression within both arms).
    n_nonfinite = int((~np.isfinite(welch_rank.to_numpy())).sum())
    welch_rank = welch_rank[np.isfinite(welch_rank)]

    gsea_hall_w = run_prerank(welch_rank, {k: hallmark[k] for k in needed_h} | {
        k: hallmark[k] for k in hallmark
    }, min_size=15, tag="welch_t")
    # The dict union above passes the full hallmark. Tag the collection.
    gsea_hall_w.insert(0, "collection", "hallmark_50")
    gsea_hall_d = run_prerank(delta_rank, hallmark, min_size=15, tag="delta_log2")
    gsea_hall_d.insert(0, "collection", "hallmark_50")
    gsea_foc_w = run_prerank(welch_rank, focused, min_size=10, tag="welch_t")
    gsea_foc_w.insert(0, "collection", "junction_focused")
    gsea_foc_d = run_prerank(delta_rank, focused, min_size=10, tag="delta_log2")
    gsea_foc_d.insert(0, "collection", "junction_focused")
    gsea = pd.concat([gsea_hall_w, gsea_hall_d, gsea_foc_w, gsea_foc_d], ignore_index=True)

    # Cluster-equal rank: mean log2 difference across clusters, same gene universe.
    cl_delta = cluster_expr.loc[delta_rank.index, kl_clusters].mean(axis=1) - cluster_expr.loc[
        delta_rank.index, kp_clusters
    ].mean(axis=1)
    gsea_cl = run_prerank(cl_delta, {**focused, **{k: hallmark[k] for k in needed_h}}, min_size=10, tag="cluster_delta")
    gsea_cl.insert(0, "collection", "cluster_equal_weight")

    print("ssGSEA ...", flush=True)
    ss = run_ssgsea(log_expr, ssgsea_sets)
    # gseapy names the sample column "Name" and the set "Term".
    ss_use = ss.rename(columns={"Name": "library", "Term": "set"})
    # Name may be the sample. Confirm.
    if "library" not in ss_use.columns or ss_use["library"].iloc[0] not in libs:
        raise SystemExit(f"Unexpected ssGSEA columns: {list(ss.columns)} head={ss.head(2).to_dict()}")
    ss_contrast = contrast_scores(ss_use, "ES", "set", "library")

    # IFN-cold rule: both hallmark IFN sets completely lower in every KL library.
    def ifn_call(term: str) -> dict:
        row = ss_contrast.loc[ss_contrast["set"] == term].iloc[0]
        mean_row = mean_tbl.loc[mean_tbl["set"] == term].iloc[0]
        kl_max = max(float(row[SHORT[lib]]) for lib in KL_LIBS)
        kp_min = min(float(row[SHORT[lib]]) for lib in KP_LIBS)
        return {
            "set": term,
            "ssgsea_delta_kl_minus_kp": float(row["delta_mean"]),
            "ssgsea_p": float(row["p"]),
            "ssgsea_complete_kl_lower": bool(kl_max < kp_min),
            "ssgsea_hits_floor": bool(row["hits_floor"]),
            "meanlog_delta": float(mean_row["library_delta"]),
            "meanlog_p": float(mean_row["library_p"]),
            "meanlog_complete_kl_lower": bool(mean_row["library_complete"] and mean_row["library_delta"] < 0),
            "kl_max_es": kl_max,
            "kp_min_es": kp_min,
        }

    ifn_alpha = ifn_call("HALLMARK_INTERFERON_ALPHA_RESPONSE")
    ifn_gamma = ifn_call("HALLMARK_INTERFERON_GAMMA_RESPONSE")
    ifn_cold = bool(ifn_alpha["ssgsea_complete_kl_lower"] and ifn_gamma["ssgsea_complete_kl_lower"])

    # Ladder and soft-power summary.
    ladder = mw_ladder()
    soft = {
        "n_kl_libraries": 5,
        "n_kp_libraries": 5,
        "n_assignments": math.comb(10, 5),
        "two_sided_floor_p": FLOOR_P,
        "one_sided_floor_p": 1 / math.comb(10, 5),
        "attainable_p_still_below_0.05": 8 / 252,
        "min_cliffs_delta_to_reject_0.05": 21 / 25,
        "u_equals_22_p": 14 / 252,
        "n_genes_max_fpkm_ge1": n_tested,
        "n_genes_at_mw_floor": n_floor,
        "n_genes_complete_separation": n_complete,
        "bh_k_needed_for_floor_to_pass_fdr_0.05": k_needed,
        "n_genes_bh_fdr_0.05": n_bh,
        "floor_passes_bh": bool(n_floor >= k_needed and n_bh > 0),
        "corr_threshold": CORR_THR,
        "kp_clusters": sorted(set(labels[c] for c in KP_LIBS)),
        "kl_clusters": sorted(set(labels[c] for c in KL_LIBS)),
        "n_kp_clusters": len(kp_clusters),
        "n_kl_clusters": len(kl_clusters),
        "cluster_two_sided_floor_p": 2 / math.comb(len(kl_clusters) + len(kp_clusters), len(kl_clusters)),
        "n_welch_nonfinite_dropped": n_nonfinite,
        "ifn_cold_complete_separation": ifn_cold,
    }

    # Provenance of gene-set overlap actually used.
    universe = set(welch_rank.index)
    overlap_rows = []
    for name, genes in {**focused, **{k: hallmark[k] for k in needed_h}, "COMPACT_ISG": curated["COMPACT_ISG"]}.items():
        present = [g for g in genes if g in fpkm.index]
        in_rank = [g for g in present if g in universe]
        overlap_rows.append(
            {
                "set": name,
                "n_listed": len(genes),
                "n_in_fpkm": len(present),
                "n_in_gsea_universe": len(in_rank),
                "missing_from_fpkm": ",".join(g for g in genes if g not in fpkm.index),
            }
        )
    overlap = pd.DataFrame(overlap_rows)

    # Leading-edge extract for the primary (Welch, both collections).
    lead_terms = [
        "EPITHELIAL_TJ_CORE",
        "GO_BICELLULAR_TIGHT_JUNCTION",
        "KEGG_TIGHT_JUNCTION",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_APICAL_JUNCTION",
        "GOCC_GAP_JUNCTION",
    ]
    lead = gsea.loc[
        (gsea["rank_metric"] == "welch_t") & (gsea["Term"].isin(lead_terms)),
        ["collection", "Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes"],
    ].copy()

    # Write tables.
    meta.to_csv(TABLES / "samples.tsv", sep="\t", index=False)
    locked.to_csv(TABLES / "locked_genes.tsv", sep="\t", index=False)
    callouts.to_csv(TABLES / "tj_gene_callouts.tsv", sep="\t", index=False)
    expressed.sort_values(["mw_p", "delta_log2"], ascending=[True, False]).to_csv(
        TABLES / "gene_stats_expressed.tsv.gz", sep="\t", index=False
    )
    ladder.to_csv(TABLES / "mw_exact_ladder.tsv", sep="\t", index=False)
    mean_tbl.to_csv(TABLES / "mean_signature_scores.tsv", sep="\t", index=False)
    cluster_genes.to_csv(TABLES / "cluster_locked_genes.tsv", sep="\t", index=False)
    gsea.to_csv(TABLES / "gsea_prerank.tsv", sep="\t", index=False)
    gsea_cl.to_csv(TABLES / "gsea_prerank_cluster.tsv", sep="\t", index=False)
    ss_use.to_csv(TABLES / "ssgsea_long.tsv", sep="\t", index=False)
    ss_contrast.to_csv(TABLES / "ssgsea_contrasts.tsv", sep="\t", index=False)
    pd.DataFrame([ifn_alpha, ifn_gamma]).to_csv(TABLES / "ifn_cold.tsv", sep="\t", index=False)
    overlap.to_csv(TABLES / "geneset_overlap.tsv", sep="\t", index=False)
    lead.to_csv(TABLES / "gsea_leading_edge.tsv", sep="\t", index=False)
    corr.loc[libs, libs].to_csv(TABLES / "library_correlation.tsv", sep="\t")
    pd.DataFrame(
        [{"library": lib, "short": SHORT[lib], "arm": "KL" if lib in KL_LIBS else "KP", "cluster": labels[lib]} for lib in libs]
    ).to_csv(TABLES / "library_clusters.tsv", sep="\t", index=False)

    # Per-library log2 for the genes a reader will want next to the volcano.
    watch = ["Tacstd2", "Cldn4", "Cldn1", "Cldn3", "Cldn7", "Cldn6", "Cldn2", "Cldn5", "Cldn18", "Ocln", "Tjp1", "Stk11", "Trp53"]
    per_lib = log_expr.loc[[g for g in watch if g in log_expr.index], libs].T
    per_lib.index = [SHORT[i] for i in per_lib.index]
    per_lib.to_csv(TABLES / "watch_genes_log2.tsv", sep="\t")

    summary = {
        "soft_power": soft,
        "locked": locked[["gene", "delta_log2", "mw_p", "cliffs_delta", "complete_separation", "hits_floor"]].to_dict(
            orient="records"
        ),
        "cluster_genes": cluster_genes.to_dict(orient="records"),
        "ifn_alpha": ifn_alpha,
        "ifn_gamma": ifn_gamma,
        "ifn_cold": ifn_cold,
        "mean_signatures": mean_tbl[
            ["set", "n_genes_in_matrix", "library_delta", "library_p", "library_hits_floor", "library_complete", "cluster_delta", "cluster_p", "cluster_complete"]
        ].to_dict(orient="records"),
    }
    # GSEA headline rows
    headline_terms = lead_terms
    gsea_head = gsea.loc[gsea["Term"].isin(headline_terms), ["collection", "rank_metric", "Term", "NES", "NOM p-val", "FDR q-val", "ES", "Tag %"]]
    summary["gsea_headline"] = gsea_head.to_dict(orient="records")
    gsea_cl_head = gsea_cl.loc[gsea_cl["Term"].isin(headline_terms), ["rank_metric", "Term", "NES", "NOM p-val", "FDR q-val"]]
    summary["gsea_cluster"] = gsea_cl_head.to_dict(orient="records")
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    print("Figures ...", flush=True)
    plot_volcano(stats, callouts)
    # wide frame for the strip plot: one row per set-library from ss_use
    plot_frame = ss_use.loc[ss_use["set"].isin(SS_ORDER), ["set", "library", "ES"]].copy()
    plot_ssgsea(plot_frame, ss_contrast)
    plot_gsea(gsea)
    plot_corr(corr, labels)

    # Reproduction guards. Handoff rounds are +3.24 and +5.57.
    tac = float(locked.loc[locked["gene"] == "Tacstd2", "delta_log2"].iloc[0])
    cld = float(locked.loc[locked["gene"] == "Cldn4", "delta_log2"].iloc[0])
    tac_p = float(locked.loc[locked["gene"] == "Tacstd2", "mw_p"].iloc[0])
    if abs(tac - 3.2382) > 0.005 or abs(cld - 5.5696) > 0.005:
        raise SystemExit(f"Locked deltas did not reproduce: Tacstd2 {tac}, Cldn4 {cld}")
    if not np.isclose(tac_p, FLOOR_P, rtol=0, atol=1e-12):
        raise SystemExit(f"Tacstd2 p {tac_p} is not the MW floor {FLOOR_P}")
    print(json.dumps({"Tacstd2": tac, "Cldn4": cld, "floor": FLOOR_P, "ifn_cold": ifn_cold, "n_bh": n_bh}, indent=2))
    print("Done", flush=True)


if __name__ == "__main__":
    main()

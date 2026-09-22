#!/usr/bin/env python3
"""PAPER FUNNEL: Tacstd2-high vs low malignant enrichment GSEA.

Cohorts (public only; no private 8KL):
  1. GSE137244 — KL/KP GEMM lung-tumor nodule–derived epithelial cell lines (n=10).
  2. GSE164758 — public GEMM primary NSCLC tumors, untreated KL+KP (n=17).
  3. GSE137396 — public GEMM lung tumor nodules, KL+KP (n=10).

Design
  - Split samples by Tacstd2 median (drop Tacstd2 from the ranked list).
  - Rank metric = Welch t (Tacstd2-high minus Tacstd2-low) on log2(expr+1)
    when the deposited matrix is FPKM/raw-like; deposited log-scale matrices
    are used as-is.
  - Preranked GSEA (gseapy), gene-set permutation, seed fixed.
  - Ranking universe = Hallmark + KEGG_2019_Mouse + WikiPathways_2019_Mouse
    + GO CC junction/adhesion terms + curated Claudin family.
  - Force-report every set with NES / NOM p / FDR; highlight whether
    tight-junction / cell-adhesion / Claudin-family terms fall in ranks 1–3
    among positive-NES sets (sorted by NES descending, ties by FDR).

Honesty
  - Never invent NES/FDR.
  - On GSE137244 Tacstd2 completely separates KL from KP, so the Tacstd2
    median split is genotype-confounded; that is reported, not hidden.
  - LLC / TISMO are not used as KL substitutes.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "cache"
GENESETS = ROOT / "gene_sets"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
for d in (CACHE, GENESETS, TABLES, FIGS):
    d.mkdir(parents=True, exist_ok=True)

SEED = 49901
N_PERM = 2000
MIN_SIZE = 10
MAX_SIZE = 500

# Highlight keywords for TJ / adhesion / claudin family.
HIGHLIGHT_PATTERNS = [
    re.compile(r"tight.?junction", re.I),
    re.compile(r"bicellular.?tight", re.I),
    re.compile(r"cell.?adhesion", re.I),
    re.compile(r"adherens.?junction", re.I),
    re.compile(r"claudin", re.I),
]

CLAUDIN_FAMILY = [
    f"Cldn{i}" for i in list(range(1, 28)) + [34]
] + [
    "Cldn34a",
    "Cldn34b1",
    "Cldn34b2",
    "Cldn34b3",
    "Cldn34b4",
    "Cldn34c1",
    "Cldn34c2",
    "Cldn34c4",
    "Cldn34c6",
    "Cldn34d",
]


def load_gmt(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            genes = []
            for g in parts[2:]:
                if g and g not in genes:
                    genes.append(g)
            out[parts[0]] = genes
    return out


def load_enrichr_json(path: Path) -> dict[str, list[str]]:
    raw = json.loads(path.read_text())
    # Enrichr mouse libraries ship UPPERCASE symbols; keep as deposited and
    # match case-insensitively against the ranked list later.
    return {k: list(dict.fromkeys(v)) for k, v in raw.items()}


def titlecase_mouse(sym: str) -> str:
    """Best-effort Enrichr UPPER -> mouse Title case (Cldn4, H2-K1, etc.)."""
    if not sym:
        return sym
    parts = []
    for chunk in sym.split("-"):
        if not chunk:
            parts.append(chunk)
            continue
        if chunk.isdigit():
            parts.append(chunk)
        else:
            parts.append(chunk[0].upper() + chunk[1:].lower())
    return "-".join(parts)


def to_mouse_symbols(sets: dict[str, list[str]]) -> dict[str, list[str]]:
    out = {}
    for name, genes in sets.items():
        mapped = []
        for g in genes:
            mg = titlecase_mouse(g)
            if mg not in mapped:
                mapped.append(mg)
        out[name] = mapped
    return out


def build_ranking_universe() -> dict[str, list[str]]:
    hallmark = load_gmt(GENESETS / "mh.all.v2024.1.Mm.symbols.gmt")
    go_cc = load_gmt(GENESETS / "m5.go.cc.v2024.1.Mm.symbols.gmt")
    go_bp = load_gmt(GENESETS / "m5.go.bp.v2024.1.Mm.symbols.gmt")
    reactome = load_gmt(GENESETS / "m2.cp.reactome.v2024.1.Mm.symbols.gmt")
    kegg = to_mouse_symbols(load_enrichr_json(GENESETS / "KEGG_2019_Mouse.json"))
    wiki = to_mouse_symbols(load_enrichr_json(GENESETS / "WikiPathways_2019_Mouse.json"))

    # GO CC / BP / Reactome terms that are TJ, adhesion, or Claudin-named.
    go_focus = {}
    for coll in (go_cc, go_bp, reactome):
        for name, genes in coll.items():
            if any(p.search(name) for p in HIGHLIGHT_PATTERNS):
                go_focus[name] = genes

    universe: dict[str, list[str]] = {}
    universe.update({f"HALLMARK__{k}": v for k, v in hallmark.items()})
    universe.update({f"KEGG__{k}": v for k, v in kegg.items()})
    universe.update({f"WIKI__{k}": v for k, v in wiki.items()})
    universe.update({f"GOFOCUS__{k}": v for k, v in go_focus.items()})
    universe["CURATED__CLAUDIN_FAMILY"] = CLAUDIN_FAMILY
    return universe


def welch_t_matrix(expr: pd.DataFrame, hi: list[str], lo: list[str]) -> pd.Series:
    a = expr[hi].to_numpy(dtype=float)
    b = expr[lo].to_numpy(dtype=float)
    t, _ = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    return pd.Series(t, index=expr.index)


def median_split(series: pd.Series) -> tuple[list[str], list[str], float]:
    order = series.sort_values(ascending=False)
    n = len(order)
    hi = list(order.index[: n // 2])
    lo = list(order.index[n // 2 :])
    # split threshold = midpoint between last high and first low
    thr = float((order.iloc[n // 2 - 1] + order.iloc[n // 2]) / 2.0)
    return hi, lo, thr


def run_prerank(rank: pd.Series, gene_sets: dict[str, list[str]], tag: str) -> pd.DataFrame:
    import gseapy as gp

    rnk = rank.replace([np.inf, -np.inf], np.nan).dropna()
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    # drop non-finite
    rnk = rnk[np.isfinite(rnk.to_numpy())]
    pre = gp.prerank(
        rnk=rnk.sort_values(ascending=False).rename("score").reset_index(),
        gene_sets=gene_sets,
        min_size=MIN_SIZE,
        max_size=MAX_SIZE,
        permutation_num=N_PERM,
        seed=SEED,
        threads=4,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    df = pre.res2d.copy()
    df.insert(0, "cohort", tag)
    df.insert(1, "n_genes_ranked", int(rnk.shape[0]))
    # standardize numeric columns
    for col in ("ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def is_highlight(term: str) -> bool:
    return any(p.search(term) for p in HIGHLIGHT_PATTERNS) or term.endswith("CLAUDIN_FAMILY")


def highlight_family(term: str) -> str:
    t = term.lower()
    if "claudin" in t:
        return "claudin_family"
    if "tight" in t and "junction" in t:
        return "tight_junction"
    if "bicellular" in t and "tight" in t:
        return "tight_junction"
    if "cell adhesion" in t or "cell_adhesion" in t:
        return "cell_adhesion"
    if "adherens" in t:
        return "cell_adhesion"
    return "other_junction_like"


def rank_positive(gsea: pd.DataFrame) -> pd.DataFrame:
    pos = gsea.loc[gsea["NES"] > 0].copy()
    pos = pos.sort_values(["NES", "FDR q-val"], ascending=[False, True], kind="mergesort")
    pos.insert(0, "rank_pos_NES", np.arange(1, len(pos) + 1))
    return pos


def load_gse137244() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = CACHE / "GSE137244_counts.fpkm.csv.gz"
    raw = pd.read_csv(path)
    expr = raw.groupby("gene", as_index=True).mean(numeric_only=True)
    libs = [
        "B6AL10-1-RNA",
        "B6AL10-2-RNA",
        "B6AL10-3-RNA",
        "B6AL10-4-RNA",
        "B6AL10-5-RNA",
        "KL155mix-control-2-RNA",
        "KL47-1-untreated-1-RNA",
        "KLC-RNA",
        "KLD-RNA",
        "KLE-RNA",
    ]
    arm = ["KP"] * 5 + ["KL"] * 5
    meta = pd.DataFrame({"sample": libs, "arm": arm, "setting": "GEMM_derived_epithelial_cell_line"})
    log = np.log2(expr[libs].clip(lower=0) + 1.0)
    return log, meta


def load_gse137396() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(CACHE / "GSE137396_KL_vs_KP_nodules.expr.tsv.gz", sep="\t", index_col=0)
    design = pd.read_csv(CACHE / "GSE137396_KL_vs_KP_nodules.design.tsv", sep="\t")
    cols = list(design["column"])
    meta = design.rename(columns={"column": "sample", "group": "arm"})[
        ["sample", "arm"]
    ].copy()
    meta["setting"] = "GEMM_in_vivo_nodule"
    # Matrix is already log-like (values ~0–15); do not double-log.
    return expr[cols], meta


def load_gse164758() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = CACHE / "GSE164758_primary_tumors_fpkm.txt.gz"
    df = pd.read_csv(path, sep="\t")
    ann_i = list(df.columns).index("Annotation/Divergence")
    fpkm_cols = [c for c in df.columns[ann_i + 1 :] if str(c).endswith("FPKM")]
    symbols = df["Annotation/Divergence"].astype(str).str.split("|").str[0]
    mat = df[fpkm_cols].copy()
    mat.index = symbols
    # collapse duplicate symbols by mean
    mat = mat.groupby(level=0).mean()
    shorts = [re.sub(r"Aligned\.out\.sam FPKM$", "", c) for c in fpkm_cols]
    mat.columns = shorts

    def bucket(s: str) -> str:
        if s.startswith("KL"):
            return "KL"
        if s.startswith("KPL"):
            return "KPL"
        if s.startswith("KP"):
            return "KP"
        if s.startswith("Kras"):
            return "K"
        return "other"

    # Primary public-GEMM contrast: untreated KL + KP primary tumors only
    # (exclude K, KPL; KPa columns are absent from this FPKM deposit).
    keep = [s for s in shorts if bucket(s) in ("KL", "KP")]
    meta = pd.DataFrame(
        {
            "sample": keep,
            "arm": [bucket(s) for s in keep],
            "setting": "GEMM_in_vivo_primary_tumor",
        }
    )
    log = np.log2(mat[keep].clip(lower=0) + 1.0)
    return log, meta


def cohort_runners():
    return [
        ("GSE137244_epithelium_cell_lines", load_gse137244),
        ("GSE164758_GEMM_KL_KP_primary", load_gse164758),
        ("GSE137396_GEMM_KL_KP_nodules", load_gse137396),
    ]


def analyze_cohort(
    tag: str,
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    universe: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    if "Tacstd2" not in expr.index:
        raise SystemExit(f"{tag}: Tacstd2 missing from matrix")
    tac = expr.loc["Tacstd2"].astype(float)
    meta = meta.copy()
    meta["Tacstd2"] = meta["sample"].map(tac.to_dict())
    hi, lo, thr = median_split(tac)
    meta["tac_group"] = meta["sample"].map(
        lambda s: "Tacstd2_high" if s in hi else "Tacstd2_low"
    )
    # genotype confound table
    xtab = pd.crosstab(meta["arm"], meta["tac_group"])
    # rank
    rnk = welch_t_matrix(expr, hi, lo).dropna()
    rnk = rnk.drop(index="Tacstd2", errors="ignore")
    # filter genes expressed somewhere
    keep = expr.max(axis=1) > 0
    rnk = rnk.loc[rnk.index.intersection(expr.index[keep])]
    gsea = run_prerank(rnk, universe, tag)
    pos = rank_positive(gsea)
    hi_all = gsea.loc[gsea["Term"].map(is_highlight)].copy()
    hi_all["family"] = hi_all["Term"].map(highlight_family)
    hi_pos = pos.loc[pos["Term"].map(is_highlight)].copy()
    hi_pos["family"] = hi_pos["Term"].map(highlight_family)
    # top-3 check per family among positive-NES sets; also force-report
    # the strongest NES (can be negative) so nulls are never invented.
    checks = []
    for fam in ("tight_junction", "cell_adhesion", "claudin_family"):
        sub_pos = hi_pos.loc[hi_pos["family"] == fam]
        sub_all = hi_all.loc[hi_all["family"] == fam].sort_values(
            "NES", ascending=False, kind="mergesort"
        )
        if sub_all.empty:
            checks.append(
                {
                    "cohort": tag,
                    "family": fam,
                    "best_term": None,
                    "rank_pos_NES": None,
                    "NES": None,
                    "FDR_q": None,
                    "NOM_p": None,
                    "in_top3": False,
                    "direction": "absent",
                    "n_positive_sets_total": int(pos.shape[0]),
                }
            )
            continue
        best_any = sub_all.iloc[0]
        if sub_pos.empty:
            checks.append(
                {
                    "cohort": tag,
                    "family": fam,
                    "best_term": best_any["Term"],
                    "rank_pos_NES": None,
                    "NES": float(best_any["NES"]),
                    "FDR_q": float(best_any["FDR q-val"])
                    if pd.notna(best_any["FDR q-val"])
                    else None,
                    "NOM_p": float(best_any["NOM p-val"])
                    if pd.notna(best_any["NOM p-val"])
                    else None,
                    "in_top3": False,
                    "direction": "not_positive",
                    "n_positive_sets_total": int(pos.shape[0]),
                }
            )
        else:
            best = sub_pos.iloc[0]
            checks.append(
                {
                    "cohort": tag,
                    "family": fam,
                    "best_term": best["Term"],
                    "rank_pos_NES": int(best["rank_pos_NES"]),
                    "NES": float(best["NES"]),
                    "FDR_q": float(best["FDR q-val"]) if pd.notna(best["FDR q-val"]) else None,
                    "NOM_p": float(best["NOM p-val"]) if pd.notna(best["NOM p-val"]) else None,
                    "in_top3": bool(int(best["rank_pos_NES"]) <= 3),
                    "direction": "positive",
                    "n_positive_sets_total": int(pos.shape[0]),
                }
            )
    summary = {
        "cohort": tag,
        "n_samples": int(meta.shape[0]),
        "n_high": len(hi),
        "n_low": len(lo),
        "tacstd2_split_threshold": thr,
        "tacstd2_high_min": float(tac[hi].min()),
        "tacstd2_low_max": float(tac[lo].max()),
        "tacstd2_vs_arm_complete_separation": bool(
            set(meta.loc[meta["tac_group"] == "Tacstd2_high", "arm"])
            .isdisjoint(set(meta.loc[meta["tac_group"] == "Tacstd2_low", "arm"]))
        ),
        "n_sets_tested": int(gsea.shape[0]),
        "n_positive_NES": int(pos.shape[0]),
        "any_highlight_in_top3": any(c["in_top3"] for c in checks),
        "cross_tab": xtab.to_dict(),
    }
    return meta, gsea, pos, {"checks": checks, "summary": summary, "xtab": xtab}


def plot_top(pos: pd.DataFrame, tag: str, out: Path, top_n: int = 25) -> None:
    df = pos.head(top_n).iloc[::-1].copy()
    colors = []
    for term in df["Term"]:
        if is_highlight(term):
            colors.append("#C44E52")
        elif float(df.loc[df["Term"] == term, "FDR q-val"].iloc[0]) < 0.05:
            colors.append("#4C78A8")
        elif float(df.loc[df["Term"] == term, "FDR q-val"].iloc[0]) < 0.25:
            colors.append("#72B7B2")
        else:
            colors.append("#B0B0B0")
    fig, ax = plt.subplots(figsize=(10, 0.34 * len(df) + 1.8))
    ax.barh(range(len(df)), df["NES"].astype(float), color=colors)
    labels = []
    for _, row in df.iterrows():
        lab = str(row["Term"])
        if len(lab) > 70:
            lab = lab[:67] + "..."
        labels.append(f"#{int(row['rank_pos_NES'])}  {lab}")
    ax.set_yticks(range(len(df)), labels, fontsize=7)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("NES (positive = enriched in Tacstd2-high)")
    ax.set_title(
        f"{tag}\nTop {top_n} positive-NES pathways (red = TJ / adhesion / Claudin)",
        fontsize=10,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main() -> None:
    print("Building ranking universe ...", flush=True)
    universe = build_ranking_universe()
    sizes = pd.DataFrame(
        [{"Term": k, "n_genes": len(v)} for k, v in universe.items()]
    ).sort_values("Term")
    sizes.to_csv(TABLES / "geneset_universe_sizes.tsv", sep="\t", index=False)
    print(f"  {len(universe)} sets", flush=True)

    all_meta = []
    all_gsea = []
    all_pos = []
    all_checks = []
    summaries = []

    for tag, loader in cohort_runners():
        print(f"\n=== {tag} ===", flush=True)
        expr, meta = loader()
        meta_out, gsea, pos, pack = analyze_cohort(tag, expr, meta, universe)
        meta_out.insert(0, "cohort", tag)
        all_meta.append(meta_out)
        all_gsea.append(gsea)
        all_pos.append(pos)
        all_checks.extend(pack["checks"])
        summaries.append(pack["summary"])
        # write per-cohort
        meta_out.to_csv(TABLES / f"{tag}__samples.tsv", sep="\t", index=False)
        gsea.sort_values(["NES"], ascending=False).to_csv(
            TABLES / f"{tag}__gsea_all.tsv", sep="\t", index=False
        )
        pos.to_csv(TABLES / f"{tag}__gsea_positive_ranked.tsv", sep="\t", index=False)
        hi = pos.loc[pos["Term"].map(is_highlight)].copy()
        hi["family"] = hi["Term"].map(highlight_family)
        hi.to_csv(TABLES / f"{tag}__gsea_highlight_terms.tsv", sep="\t", index=False)
        pack["xtab"].to_csv(TABLES / f"{tag}__tacstd2_by_arm.tsv", sep="\t")
        plot_top(pos, tag, FIGS / f"{tag}__top25_positive_NES.png")
        # print forced report for highlights
        print("Highlight ranks (positive NES):", flush=True)
        for c in pack["checks"]:
            print(
                f"  {c['family']}: rank={c['rank_pos_NES']} NES={c['NES']} FDR={c['FDR_q']} "
                f"top3={c['in_top3']} term={c['best_term']}",
                flush=True,
            )
        print("Top 10 positive NES:", flush=True)
        cols = ["rank_pos_NES", "Term", "NES", "NOM p-val", "FDR q-val"]
        print(pos[cols].head(10).to_string(index=False), flush=True)

    pd.concat(all_meta, ignore_index=True).to_csv(
        TABLES / "samples_all.tsv", sep="\t", index=False
    )
    pd.concat(all_gsea, ignore_index=True).to_csv(
        TABLES / "gsea_all_cohorts.tsv", sep="\t", index=False
    )
    pd.concat(all_pos, ignore_index=True).to_csv(
        TABLES / "gsea_positive_ranked_all.tsv", sep="\t", index=False
    )
    checks = pd.DataFrame(all_checks)
    checks.to_csv(TABLES / "highlight_rank_summary.tsv", sep="\t", index=False)
    summary = {
        "seed": SEED,
        "n_perm": N_PERM,
        "min_size": MIN_SIZE,
        "max_size": MAX_SIZE,
        "n_sets_in_universe": len(universe),
        "cohorts": summaries,
        "highlight_checks": all_checks,
        "any_cohort_highlight_in_top3": any(c["in_top3"] for c in all_checks),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\nDONE", json.dumps({
        "any_top3": summary["any_cohort_highlight_in_top3"],
        "checks": all_checks,
    }, indent=2, default=str), flush=True)


if __name__ == "__main__":
    main()

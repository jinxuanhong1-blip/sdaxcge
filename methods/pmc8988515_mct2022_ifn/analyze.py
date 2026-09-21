#!/usr/bin/env python3
"""Score IFN / MHC-I / tight-junction on the expression tables tied to PMC8988515."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from gsea_core import bh_fdr, gsea_prerank, median_split, rank_high_vs_low, rank_spearman_vs_target

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

HEADLINE = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    "KEGG_TIGHT_JUNCTION",
]
KERATIN = ["KRT8", "KRT18", "KRT19"]
NON_OVCA = {"ME180", "ME180C13", "IOSE397", "IOSE7576"}
ANCHORS = ["CLDN4", "XRCC1", "TP53BP1"]


def load_sets() -> dict[str, list[str]]:
    raw = json.loads((ROOT / "gene_sets.json").read_text())
    return {k: list(v) for k, v in raw["sets"].items()}


def mean_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    use = [g for g in genes if g in expr.index]
    if len(use) < 3:
        return pd.Series(dtype=float)
    sub = expr.loc[use]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    df = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    if len(df) < 5:
        return float("nan"), float("nan"), int(len(df))
    rho, p = stats.spearmanr(df["x"], df["y"])
    return float(rho), float(p), int(len(df))


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> tuple[float, float, int]:
    df = pd.concat([x.rename("x"), y.rename("y"), z.rename("z")], axis=1).dropna()
    if len(df) < 8:
        return float("nan"), float("nan"), int(len(df))
    xr = df["x"].rank()
    yr = df["y"].rank()
    zr = df["z"].rank()
    bx = np.polyfit(zr, xr, 1)
    by = np.polyfit(zr, yr, 1)
    rx = xr - (bx[0] * zr + bx[1])
    ry = yr - (by[0] * zr + by[1])
    rho, p = stats.spearmanr(rx, ry)
    return float(rho), float(p), int(len(df))


def run_gsea(rank: pd.Series, sets: dict[str, list[str]], contrast: str) -> pd.DataFrame:
    out = gsea_prerank(rank, sets)
    if out.empty:
        return out
    out.insert(0, "contrast", contrast)
    out["fdr"] = bh_fdr(out["nom_p"])
    return out.sort_values("term")


def load_s1() -> pd.DataFrame:
    df = pd.read_excel(DATA / "table_s1.xlsx", header=1)
    df.columns = [str(c).strip() for c in df.columns]
    gene_col = df.columns[0]
    rename = {gene_col: "gene"}
    for col in df.columns:
        key = col.lower()
        if key.startswith("log2 of ratio"):
            rename[col] = "log2_high_over_low"
        elif key in {"p-value", "p value"}:
            rename[col] = "p"
        elif key in {"q-value", "q value"}:
            rename[col] = "q"
        elif key.startswith("higher expression"):
            rename[col] = "higher_in"
    df = df.rename(columns=rename)
    df["gene"] = df["gene"].astype(str).str.strip()
    for col in ["log2_high_over_low", "p", "q"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["gene", "log2_high_over_low"])
    df = df[df["gene"].str.match(r"^[A-Za-z0-9][A-Za-z0-9\-\.]*$")]
    df = df.drop_duplicates("gene", keep="first")
    return df.set_index("gene")


def load_xena() -> pd.DataFrame:
    expr = pd.read_csv(DATA / "tcga_ov_hiseqv2.gz", sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.split("|").str[0]
    expr = expr[~expr.index.duplicated(keep="first")]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr


def load_coscia() -> pd.DataFrame:
    raw = pd.read_excel(DATA / "coscia_moesm1558.xlsx", header=1)
    raw.columns = [str(c).strip() for c in raw.columns]
    lines = [c for c in raw.columns[4:34]]
    piece = raw[["Gene names"] + lines].copy()
    records = []
    for _, row in piece.iterrows():
        names = str(row["Gene names"] if pd.notna(row["Gene names"]) else "")
        genes = [g.strip() for g in names.replace(",", ";").split(";") if g.strip() and g.strip() != "nan"]
        if not genes:
            continue
        vals = pd.to_numeric(row[lines], errors="coerce")
        for g in genes:
            records.append((g, vals))
    if not records:
        raise RuntimeError("Coscia table had no gene names")
    mat = pd.concat([v.rename(g) for g, v in records], axis=1).T
    mat = mat.groupby(level=0).mean()
    mat = mat.apply(pd.to_numeric, errors="coerce")
    return mat


def set_overlap(s1: pd.DataFrame, sets: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    sig = s1[s1["q"] < 0.05]
    for name, genes in sets.items():
        in_table = [g for g in genes if g in s1.index]
        in_sig = [g for g in in_table if g in sig.index]
        direction = s1.loc[in_sig, "higher_in"].astype(str).str.contains("High", case=False)
        rows.append(
            {
                "set": name,
                "n_in_table": len(in_table),
                "n_q_lt_0.05": len(in_sig),
                "n_higher_in_cldn4_high": int(direction.sum()) if len(in_sig) else 0,
                "n_higher_in_cldn4_low": int((~direction).sum()) if len(in_sig) else 0,
                "median_log2_high_over_low_all_in_set": float(s1.loc[in_table, "log2_high_over_low"].median())
                if in_table
                else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def score_rows(expr: pd.DataFrame, cldn4: pd.Series, sets: dict[str, list[str]], cohort: str) -> pd.DataFrame:
    rows = []
    keratin = mean_score(expr, KERATIN)
    for name, genes in sets.items():
        score = mean_score(expr, genes)
        rho, p, n = spearman(cldn4, score)
        pr, pp, pn = partial_spearman(cldn4, score, keratin)
        high, low = median_split(cldn4.dropna())
        sh = score.reindex(high).dropna()
        sl = score.reindex(low).dropna()
        if len(sh) > 2 and len(sl) > 2:
            t, tp = stats.ttest_ind(sh, sl, equal_var=False)
            delta = float(sh.mean() - sl.mean())
        else:
            t, tp, delta = float("nan"), float("nan"), float("nan")
        rows.append(
            {
                "cohort": cohort,
                "set": name,
                "n_genes_in_matrix": int(len([g for g in genes if g in expr.index])),
                "n": n,
                "spearman_rho_vs_CLDN4": rho,
                "spearman_p": p,
                "partial_spearman_rho_given_KRT8_18_19": pr,
                "partial_spearman_p": pp,
                "partial_n": pn,
                "median_split_high_minus_low": delta,
                "welch_t": float(t) if np.isfinite(t) else float("nan"),
                "welch_p": float(tp) if np.isfinite(tp) else float("nan"),
                "n_high": int(len(sh)),
                "n_low": int(len(sl)),
            }
        )
    return pd.DataFrame(rows)


def coscia_gene_table(mat: pd.DataFrame, genes: list[str], cohort: str) -> pd.DataFrame:
    sub = mat.drop(columns=[c for c in NON_OVCA if c in mat.columns], errors="ignore") if cohort == "ovca26" else mat
    if "CLDN4" not in sub.index:
        return pd.DataFrame()
    cldn4 = sub.loc["CLDN4"]
    rows = []
    for g in genes:
        if g not in sub.index:
            rows.append({"cohort": cohort, "gene": g, "present": 0, "n": 0, "spearman_rho_vs_CLDN4": np.nan, "spearman_p": np.nan})
            continue
        rho, p, n = spearman(cldn4, sub.loc[g])
        rows.append({"cohort": cohort, "gene": g, "present": 1, "n": n, "spearman_rho_vs_CLDN4": rho, "spearman_p": p})
    return pd.DataFrame(rows)


def plot_effects(overlap: pd.DataFrame, gsea: pd.DataFrame) -> None:
    """Plot the mean shift, not NES. NES is large here for a very small gene-level effect."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    order = HEADLINE
    labels = ["IFN-γ", "IFN-α", "MHC-I", "Tight junction"]
    med = overlap.set_index("set")
    rho = gsea[gsea["contrast"] == "Xena_TCGA_OV_Spearman"].set_index("term")
    left = [float(med.loc[t, "median_log2_high_over_low_all_in_set"]) for t in order]
    right = [float(rho.loc[t, "mean_stat"]) for t in order]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5))
    x = np.arange(len(order))
    axes[0].bar(x, left, color="#4C78A8")
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=15)
    axes[0].set_ylabel("Median log2(CLDN4 high / low)")
    axes[0].set_title("Table S1, genes inside the set")
    axes[1].bar(x, right, color="#F58518")
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=15)
    axes[1].set_ylabel("Mean Spearman ρ vs CLDN4")
    axes[1].set_title("Xena TCGA-OV primary, n=304")
    fig.suptitle("IFN side of PMC8988515 is a small positive shift", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_ifn_effect.png", dpi=160)
    fig.savefig(FIGURES / "fig_ifn_effect.pdf")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    sets = load_sets()
    s1 = load_s1()
    rank_s1 = s1["log2_high_over_low"].sort_values(ascending=False)
    gsea_s1 = run_gsea(rank_s1, sets, "S1_published_log2_ratio")
    overlap = set_overlap(s1, sets)

    expr = load_xena()
    if "CLDN4" not in expr.index:
        raise SystemExit("CLDN4 missing from Xena TCGA-OV")
    cldn4 = expr.loc["CLDN4"]
    # Xena barcodes are TCGA-XX-XXXX-01; field 4 is the sample type.
    parts_ok = []
    for c in cldn4.index:
        bits = str(c).split("-")
        if len(bits) >= 4 and bits[3][:2] == "01":
            parts_ok.append(c)
    primary = cldn4.loc[parts_ok] if parts_ok else cldn4
    expr_p = expr.loc[:, primary.index]
    high, low = median_split(primary.dropna())
    rank_x = rank_high_vs_low(expr_p, high, low)
    gsea_x = run_gsea(rank_x, sets, "Xena_TCGA_OV_Welch_t")
    rank_rho = rank_spearman_vs_target(expr_p, primary)
    gsea_rho = run_gsea(rank_rho, sets, "Xena_TCGA_OV_Spearman")

    tj_drop = {"KEGG_TIGHT_JUNCTION_drop_CLDN4": [g for g in sets["KEGG_TIGHT_JUNCTION"] if g != "CLDN4"]}
    gsea_drop = pd.concat(
        [
            run_gsea(rank_s1, tj_drop, "S1_published_log2_ratio_TJdrop"),
            run_gsea(rank_x, tj_drop, "Xena_TCGA_OV_Welch_t_TJdrop"),
        ],
        ignore_index=True,
    )

    scores = pd.concat(
        [
            score_rows(expr_p, primary, sets, "Xena_TCGA_OV_primary"),
        ],
        ignore_index=True,
    )
    kscore = mean_score(expr_p, KERATIN)
    k_rho, k_p, k_n = spearman(primary, kscore)

    coscia = load_coscia()
    focus_genes = sorted(set(ANCHORS + sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"] + sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"] + ["CLDN3", "CLDN7", "STAT1", "IRF1", "MX1", "ISG15", "IFI6", "B2M", "TAP1", "PSMB8", "PSMB9"]))
    gene_tbl = pd.concat(
        [coscia_gene_table(coscia, focus_genes, "all30"), coscia_gene_table(coscia, focus_genes, "ovca26")],
        ignore_index=True,
    )
    coscia_scores = pd.concat(
        [
            score_rows(coscia, coscia.loc["CLDN4"], sets, "Coscia_all30"),
            score_rows(
                coscia.drop(columns=list(NON_OVCA & set(coscia.columns))),
                coscia.drop(columns=list(NON_OVCA & set(coscia.columns))).loc["CLDN4"],
                sets,
                "Coscia_ovca26",
            ),
        ],
        ignore_index=True,
    )

    gsea = pd.concat([gsea_s1, gsea_x, gsea_rho], ignore_index=True)
    gsea.to_csv(TABLES / "gsea_prerank.tsv", sep="\t", index=False)
    gsea_drop.to_csv(TABLES / "gsea_tj_drop_cldn4.tsv", sep="\t", index=False)
    overlap.to_csv(TABLES / "s1_set_overlap.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / "score_spearman.tsv", sep="\t", index=False)
    coscia_scores.to_csv(TABLES / "coscia_score_spearman.tsv", sep="\t", index=False)
    gene_tbl.to_csv(TABLES / "coscia_gene_spearman.tsv", sep="\t", index=False)

    headline = gsea[gsea["contrast"].isin(["S1_published_log2_ratio", "Xena_TCGA_OV_Welch_t"])][
        ["contrast", "term", "nes", "fdr", "nom_p", "n_set_in_rank", "mean_stat"]
    ]
    headline.to_csv(TABLES / "gsea_headline.tsv", sep="\t", index=False)

    n_q = int((s1["q"] < 0.05).sum())
    summary = {
        "s1_n_genes": int(s1.shape[0]),
        "s1_n_q_lt_0.05": n_q,
        "s1_cldn4_log2_ratio": float(s1.loc["CLDN4", "log2_high_over_low"]) if "CLDN4" in s1.index else None,
        "xena_n_primary": int(primary.dropna().shape[0]),
        "xena_n_high": int(len(high)),
        "xena_n_low": int(len(low)),
        "xena_cldn4_vs_krt8_18_19_rho": k_rho,
        "xena_cldn4_vs_krt8_18_19_p": k_p,
        "xena_cldn4_vs_krt8_18_19_n": k_n,
        "coscia_n_lines_all": int(coscia.shape[1]),
        "coscia_n_lines_ovca26": int(coscia.shape[1] - len(NON_OVCA & set(coscia.columns))),
        "coscia_cldn4_present": bool("CLDN4" in coscia.index),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    non01 = [c for c in expr.columns if c not in set(primary.index)]
    pd.DataFrame({"sample": non01, "reason": "not sample-type 01"}).to_csv(
        TABLES / "xena_excluded_samples.tsv", sep="\t", index=False
    )
    cldn4_lines = coscia.loc["CLDN4"].dropna().rename("CLDN4_log2_LFQ").reset_index()
    cldn4_lines.columns = ["cell_line", "CLDN4_log2_LFQ"]
    cldn4_lines["ovca_line"] = ~cldn4_lines["cell_line"].isin(NON_OVCA)
    cldn4_lines.to_csv(TABLES / "coscia_cldn4_detected_lines.tsv", sep="\t", index=False)
    plot_effects(overlap, gsea)
    print(json.dumps(summary, indent=2))
    print(headline.to_string(index=False))
    print(scores.to_string(index=False))
    print(coscia_scores.to_string(index=False))
    anchors = gene_tbl[gene_tbl["gene"].isin(ANCHORS)]
    print(anchors.to_string(index=False))


if __name__ == "__main__":
    main()

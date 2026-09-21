#!/usr/bin/env python3
"""Recheck CLDN4, IFN, and MHC-I/APM on GSE289287 Trop-2 KO T-47D xenografts.

Primary numbers are the author DESeq2 columns in the deposited table.
The recheck adds a Welch t-test and an exact label permutation on
log2(author-normalized counts + 1), because an earlier note only said
the CLDN4 padj might be non-significant.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz"
SETS = ROOT / "gene_sets.json"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

WT_ANIMALS = ["2808", "2810", "2812"]
KO_ANIMALS = ["2807", "2815", "2817", "2818"]
WT_GSM = {"2808": "GSM8788420", "2810": "GSM8788421", "2812": "GSM8788419"}
KO_GSM = {"2807": "GSM8788425", "2815": "GSM8788422", "2817": "GSM8788423", "2818": "GSM8788424"}


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    filled = np.empty(n)
    filled[order] = q
    out[ok] = filled
    return out


def load_table() -> pd.DataFrame:
    df = pd.read_csv(DATA, sep="\t")
    for col in ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "stat"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    norm_cols = [c for c in df.columns if c.endswith("_normCounts")]
    raw_cols = [c for c in df.columns if c.endswith("_rawCounts")]
    df[norm_cols + raw_cols] = df[norm_cols + raw_cols].apply(pd.to_numeric, errors="coerce")
    return df


def animal_of(col: str) -> str:
    # OV.2808.RNA_normCounts -> 2808
    return col.split(".")[1]


def group_cols(df: pd.DataFrame, kind: str) -> tuple[list[str], list[str]]:
    cols = [c for c in df.columns if c.endswith(f"_{kind}")]
    wt = [c for c in cols if animal_of(c) in WT_ANIMALS]
    ko = [c for c in cols if animal_of(c) in KO_ANIMALS]
    if len(wt) != 3 or len(ko) != 4:
        raise SystemExit(f"unexpected {kind} columns: WT={wt} KO={ko}")
    return wt, ko


def exact_perm_p(wt: np.ndarray, ko: np.ndarray) -> tuple[float, float, int]:
    """Two-sided exact p for mean(KO) - mean(WT). All C(7,3) labelings."""
    values = np.concatenate([wt, ko])
    n_wt = len(wt)
    obs = float(ko.mean() - wt.mean())
    n = 0
    extreme = 0
    for idx in itertools.combinations(range(len(values)), n_wt):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(idx)] = True
        diff = float(values[~mask].mean() - values[mask].mean())
        n += 1
        if abs(diff) + 1e-12 >= abs(obs):
            extreme += 1
    return obs, extreme / n, n


def welch_row(wt: np.ndarray, ko: np.ndarray) -> tuple[float, float, float]:
    res = stats.ttest_ind(ko, wt, equal_var=False, alternative="two-sided")
    # Cohen's d, KO minus WT, pooled SD
    n1, n2 = len(ko), len(wt)
    v1, v2 = ko.var(ddof=1), wt.var(ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    d = (ko.mean() - wt.mean()) / sp if sp > 0 else np.nan
    return float(res.statistic), float(res.pvalue), float(d)


def one_symbol(df: pd.DataFrame, symbol: str) -> pd.Series | None:
    hit = df[df["Feature_name"] == symbol]
    if hit.empty:
        return None
    hit = hit.sort_values(["baseMean"], ascending=False)
    pc = hit[hit["biotype"] == "protein_coding"]
    return (pc if len(pc) else hit).iloc[0]


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    df = load_table()
    wt_norm, ko_norm = group_cols(df, "normCounts")
    wt_raw, ko_raw = group_cols(df, "rawCounts")
    sets = json.loads(SETS.read_text())

    inventory = []
    for animal, gsm in {**WT_GSM, **KO_GSM}.items():
        inventory.append(
            {
                "animal": animal,
                "gsm": gsm,
                "genotype": "WT" if animal in WT_ANIMALS else "Trop2_KO",
                "context": "xenograft",
                "norm_col": f"OV.{animal}.RNA_normCounts",
                "raw_libsize": int(df[f"OV.{animal}.RNA_rawCounts"].sum()),
            }
        )
    inv = pd.DataFrame(inventory).sort_values(["genotype", "animal"])
    inv.to_csv(TABLES / "sample_inventory.tsv", sep="\t", index=False)

    # Independent tests on log2(norm + 1) for every protein-coding symbol
    # that has a finite author p-value (the genes DESeq2 actually tested
    # before independent filtering removed padj).
    pc = df[df["biotype"] == "protein_coding"].copy()
    pc = pc.sort_values("baseMean", ascending=False).drop_duplicates("Feature_name")
    log_wt = np.log2(pc[wt_norm].to_numpy(dtype=float) + 1.0)
    log_ko = np.log2(pc[ko_norm].to_numpy(dtype=float) + 1.0)
    welch_p = np.empty(len(pc))
    welch_t = np.empty(len(pc))
    cohens_d = np.empty(len(pc))
    mean_diff = log_ko.mean(axis=1) - log_wt.mean(axis=1)
    for i in range(len(pc)):
        t, p, d = welch_row(log_wt[i], log_ko[i])
        welch_t[i] = t
        welch_p[i] = p
        cohens_d[i] = d
    pc = pc.copy()
    pc["welch_t"] = welch_t
    pc["welch_p"] = welch_p
    pc["welch_padj"] = bh(welch_p)
    pc["cohens_d_log2norm"] = cohens_d
    pc["mean_log2norm_diff"] = mean_diff

    focus_symbols = [
        "TACSTD2",
        "DSG2",
        "DSG3",
        "MMP14",
        "CLDN1",
        "CLDN3",
        "CLDN4",
        "CLDN7",
        "EPCAM",
        "ISG15",
        "IFI44L",
        "IFI44",
        "IFIT1",
        "MX1",
        "OAS2",
        "STAT1",
        "IRF1",
        "CXCL10",
        "CD274",
        "B2M",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "TAP1",
        "PSMB8",
        "PSMB9",
        "NLRC5",
    ]
    # CLDN4 exact permutation is the recheck the note asked for.
    recheck_rows = []
    for symbol in ["TACSTD2", "CLDN4", "CLDN7", "ISG15", "B2M"]:
        row = one_symbol(df, symbol)
        if row is None:
            continue
        wt = np.log2(row[wt_norm].to_numpy(dtype=float) + 1.0)
        ko = np.log2(row[ko_norm].to_numpy(dtype=float) + 1.0)
        obs, perm_p, nperm = exact_perm_p(wt, ko)
        w = pc.loc[pc["Feature_name"] == symbol].iloc[0]
        recheck_rows.append(
            {
                "symbol": symbol,
                "ensembl": row["Ensembl_Id"],
                "baseMean": row["baseMean"],
                "author_log2FC": row["log2FoldChange"],
                "author_lfcSE": row["lfcSE"],
                "author_stat": row["stat"],
                "author_pvalue": row["pvalue"],
                "author_padj": row["padj"],
                "author_significant_DE": row["significant_DE"],
                "welch_t": w["welch_t"],
                "welch_p": w["welch_p"],
                "welch_padj_protein_coding": w["welch_padj"],
                "cohens_d_log2norm": w["cohens_d_log2norm"],
                "mean_log2norm_KO_minus_WT": obs,
                "exact_perm_p": perm_p,
                "n_permutations": nperm,
                "wt_log2norm": ";".join(f"{v:.4f}" for v in wt),
                "ko_log2norm": ";".join(f"{v:.4f}" for v in ko),
            }
        )
    recheck = pd.DataFrame(recheck_rows)
    recheck.to_csv(TABLES / "cldn4_recheck.tsv", sep="\t", index=False)

    key = []
    for symbol in focus_symbols:
        row = one_symbol(df, symbol)
        if row is None:
            key.append({"symbol": symbol, "present": False})
            continue
        w = pc.loc[pc["Feature_name"] == symbol]
        key.append(
            {
                "symbol": symbol,
                "present": True,
                "ensembl": row["Ensembl_Id"],
                "baseMean": row["baseMean"],
                "log2FC": row["log2FoldChange"],
                "lfcSE": row["lfcSE"],
                "stat": row["stat"],
                "pvalue": row["pvalue"],
                "padj": row["padj"],
                "welch_p": float(w["welch_p"].iloc[0]) if len(w) else np.nan,
                "welch_padj": float(w["welch_padj"].iloc[0]) if len(w) else np.nan,
            }
        )
    pd.DataFrame(key).to_csv(TABLES / "key_genes.tsv", sep="\t", index=False)

    def set_table(name: str, symbols: list[str]) -> pd.DataFrame:
        rows = []
        for symbol in symbols:
            row = one_symbol(df, symbol)
            if row is None:
                rows.append({"set": name, "symbol": symbol, "present": False})
                continue
            w = pc.loc[pc["Feature_name"] == symbol]
            rows.append(
                {
                    "set": name,
                    "symbol": symbol,
                    "present": True,
                    "ensembl": row["Ensembl_Id"],
                    "baseMean": row["baseMean"],
                    "log2FC": row["log2FoldChange"],
                    "lfcSE": row["lfcSE"],
                    "stat": row["stat"],
                    "pvalue": row["pvalue"],
                    "padj": row["padj"],
                    "welch_p": float(w["welch_p"].iloc[0]) if len(w) else np.nan,
                    "welch_padj": float(w["welch_padj"].iloc[0]) if len(w) else np.nan,
                }
            )
        return pd.DataFrame(rows)

    ifn = pd.concat(
        [
            set_table("HALLMARK_INTERFERON_GAMMA_RESPONSE", sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]),
            set_table("HALLMARK_INTERFERON_ALPHA_RESPONSE", sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        ],
        ignore_index=True,
    )
    apm = set_table("CUSTOM_MHC_I_ANTIGEN_PRESENTATION", sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"])
    ifn.to_csv(TABLES / "ifn_genes.tsv", sep="\t", index=False)
    apm.to_csv(TABLES / "apm_genes.tsv", sep="\t", index=False)

    # Competitive Mann-Whitney on author Wald stat.
    # Background = other protein-coding genes with a finite stat.
    bg = pc.dropna(subset=["stat"])
    summaries = []
    for name, symbols in sets.items():
        sub = set_table(name, symbols)
        present = sub[sub["present"] == True].dropna(subset=["stat"])  # noqa: E712
        in_set = set(present["symbol"])
        background = bg[~bg["Feature_name"].isin(in_set)]["stat"].to_numpy(dtype=float)
        inset = present["stat"].to_numpy(dtype=float)
        mw = stats.mannwhitneyu(inset, background, alternative="greater")
        tested = present.dropna(subset=["padj"])
        # Sample-level mean of log2(norm+1) across present genes, exact permutation.
        mat_wt = []
        mat_ko = []
        for symbol in present["symbol"]:
            row = one_symbol(df, symbol)
            mat_wt.append(np.log2(row[wt_norm].to_numpy(dtype=float) + 1.0))
            mat_ko.append(np.log2(row[ko_norm].to_numpy(dtype=float) + 1.0))
        score_wt = np.vstack(mat_wt).mean(axis=0)
        score_ko = np.vstack(mat_ko).mean(axis=0)
        obs, perm_p, nperm = exact_perm_p(score_wt, score_ko)
        summaries.append(
            {
                "set": name,
                "n_symbols": len(symbols),
                "n_present": int(present.shape[0]),
                "n_with_padj": int(tested.shape[0]),
                "n_up_log2FC": int((present["log2FC"] > 0).sum()),
                "n_down_log2FC": int((present["log2FC"] < 0).sum()),
                "n_padj_lt_0.05_up": int(((tested["padj"] < 0.05) & (tested["log2FC"] > 0)).sum()),
                "n_padj_lt_0.05_down": int(((tested["padj"] < 0.05) & (tested["log2FC"] < 0)).sum()),
                "median_log2FC": float(present["log2FC"].median()),
                "mean_log2FC": float(present["log2FC"].mean()),
                "competitive_MW_greater_p": float(mw.pvalue),
                "competitive_MW_U": float(mw.statistic),
                "sample_score_KO_minus_WT": obs,
                "sample_score_exact_perm_p": perm_p,
                "n_permutations": nperm,
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(TABLES / "geneset_summary.tsv", sep="\t", index=False)

    universe = {
        "n_rows": int(len(df)),
        "n_protein_coding_symbols_tested": int(len(pc)),
        "n_author_padj_lt_0.05": int((df["padj"] < 0.05).sum()),
        "n_author_padj_na": int(df["padj"].isna().sum()),
        "raw_libsize": {a: int(df[f"OV.{a}.RNA_rawCounts"].sum()) for a in WT_ANIMALS + KO_ANIMALS},
    }
    (TABLES / "universe.json").write_text(json.dumps(universe, indent=2))

    plot(df, wt_norm, ko_norm, pc, apm, ifn, summary)
    print(recheck.to_string(index=False))
    print(summary.to_string(index=False))
    print(json.dumps(universe, indent=2))


def plot(df, wt_norm, ko_norm, pc, apm, ifn, summary) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(11.2, 8.4), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.35], width_ratios=[1.05, 1.15])

    # A. Per-sample log2(norm+1) for the KO control gene and CLDN4.
    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(1)
    for i, symbol in enumerate(["TACSTD2", "CLDN4"]):
        row = one_symbol(df, symbol)
        wt = np.log2(row[wt_norm].to_numpy(dtype=float) + 1.0)
        ko = np.log2(row[ko_norm].to_numpy(dtype=float) + 1.0)
        xw = i - 0.16 + rng.uniform(-0.03, 0.03, size=len(wt))
        xk = i + 0.16 + rng.uniform(-0.03, 0.03, size=len(ko))
        ax.scatter(xw, wt, s=36, c="#4C78A8", zorder=3, label="WT xenograft" if i == 0 else None)
        ax.scatter(xk, ko, s=36, c="#E45756", zorder=3, label="Trop-2 KO xenograft" if i == 0 else None)
        ax.plot([i - 0.16, i + 0.16], [wt.mean(), ko.mean()], color="#222222", lw=1.2, zorder=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["TACSTD2", "CLDN4"])
    ax.set_ylabel("log2(normalized count + 1)")
    cldn = one_symbol(df, "CLDN4")
    ax.set_title(
        "A   Xenograft samples, 3 WT vs 4 KO\n"
        f"CLDN4 author log2FC {cldn['log2FoldChange']:+.3f}, p {cldn['pvalue']:.3f}, padj {cldn['padj']:.3f}"
    )
    ax.legend(frameon=False, loc="upper right")

    # B. Set-level counts.
    ax = fig.add_subplot(gs[0, 1])
    order = [
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
    ]
    labels = ["IFN-γ\nHallmark", "IFN-α\nHallmark", "MHC-I\nAPM"]
    sub = summary.set_index("set").loc[order]
    x = np.arange(len(order))
    ax.bar(x - 0.18, sub["n_padj_lt_0.05_up"], width=0.36, color="#E45756", label="padj < 0.05, up in KO")
    ax.bar(x + 0.18, sub["n_padj_lt_0.05_down"], width=0.36, color="#4C78A8", label="padj < 0.05, down in KO")
    ymax = 0
    for i, name in enumerate(order):
        r = sub.loc[name]
        top = max(r["n_padj_lt_0.05_up"], r["n_padj_lt_0.05_down"])
        ax.text(
            i,
            top + 0.35,
            f"{int(r['n_up_log2FC'])}/{int(r['n_present'])} log2FC>0",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )
        ymax = max(ymax, top)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            f"{lab}\nrank p={sub.loc[name, 'competitive_MW_greater_p']:.1e}"
            for lab, name in zip(labels, order)
        ]
    )
    ax.set_ylabel("Genes with author padj < 0.05")
    ax.set_ylim(0, ymax + 4.5)
    ax.set_title("B   Author calls. Rank p is a gene-level Mann–Whitney")
    ax.legend(frameon=False, loc="upper right", fontsize=8)

    # C. Forest of CLDN4, claudins, APM, and a short IFN panel.
    ax = fig.add_subplot(gs[1, :])
    forest_symbols = [
        "TACSTD2",
        "CLDN1",
        "CLDN3",
        "CLDN4",
        "CLDN7",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "HLA-F",
        "B2M",
        "TAP1",
        "TAP2",
        "PSMB8",
        "PSMB9",
        "NLRC5",
        "IRF1",
        "ISG15",
        "IFI44L",
        "IFIT1",
        "MX1",
        "OAS2",
        "STAT1",
        "CXCL10",
    ]
    rows = []
    for symbol in forest_symbols:
        row = one_symbol(df, symbol)
        if row is None or not np.isfinite(row["log2FoldChange"]):
            continue
        rows.append(row)
    forest = pd.DataFrame(rows).iloc[::-1]
    y = np.arange(len(forest))
    sig = forest["padj"].to_numpy(dtype=float) < 0.05
    colors = np.where(sig, "#E45756", "#7A7A7A")
    ax.axvline(0, color="#888888", lw=0.8)
    ax.errorbar(
        forest["log2FoldChange"],
        y,
        xerr=1.96 * forest["lfcSE"],
        fmt="none",
        ecolor="#B0B0B0",
        elinewidth=0.8,
        capsize=0,
        zorder=1,
    )
    ax.scatter(forest["log2FoldChange"], y, c=colors, s=28, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(forest["Feature_name"])
    ax.set_xlabel("Author DESeq2 log2FC (Trop-2 KO − WT), bars = 1.96 × lfcSE")
    ax.set_title("C   Collateral genes. Red = author padj < 0.05. CLDN4 interval crosses 0.")
    fig.savefig(FIGURES / "fig_cldn4_ifn_apm.png", dpi=160)
    fig.savefig(FIGURES / "fig_cldn4_ifn_apm.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()

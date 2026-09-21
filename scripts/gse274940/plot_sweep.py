#!/usr/bin/env python3
"""Merge the limma/edgeR/DESeq2 sweep with preranked GSEA and draw the figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "gse274940"


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order]
    n = len(ranked)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def main() -> None:
    base = pd.read_csv(RES / "sweep_all_tests.tsv", sep="\t")
    gsea = pd.read_csv(RES / "sweep_prerank_gsea.tsv", sep="\t")
    for col in ("sweep_q_within_contrast", "primary_family", "primary_family_q", "sweep_id"):
        if col not in gsea.columns:
            gsea[col] = np.nan
    gsea["primary_family"] = False
    merged = pd.concat([base, gsea], ignore_index=True, sort=False)
    # Gene-label GSEA p of 0 means 0/2000 permutations. Use the plus-one
    # estimator so a later BH is not a zero.
    zero = (merged["p_kind"] == "prerank_gsea_two_sided_absES") & (merged["p_method"] == 0)
    merged.loc[zero, "p_method"] = 1 / 2001
    merged.loc[zero, "extra"] = merged.loc[zero, "extra"].astype(str) + ";p_reported_as_lt_1/2000"

    merged["sweep_q_within_contrast"] = np.nan
    key = merged["contrast"].astype(str) + "|" + merged["filter"].astype(str)
    usable = (
        np.isfinite(merged["thesis_p"])
        & ~merged["variant"].astype(str).str.startswith("loo_gene")
        & (merged["p_kind"] != "effect_only")
    )
    for k in key[usable].unique():
        ix = np.where((key == k) & usable)[0]
        merged.loc[merged.index[ix], "sweep_q_within_contrast"] = bh(merged["thesis_p"].to_numpy()[ix])

    primary_methods = {
        "fgsea:limma_t", "fgsea:deseq2_wald", "fgsea:signal_to_noise",
        "prerank_gsea:limma_t", "prerank_gsea:deseq2_wald", "prerank_gsea:signal_to_noise",
        "camera_voom:estimated", "camera_voom:mild", "roast_voom_mean", "fry_voom",
        "zmean_welch", "gsva", "ssgsea",
    }
    prim = (
        (merged["contrast"] == "all_null_vs_WT")
        & (merged["filter"] == "cpmfilter_ge10_in2")
        & (merged["variant"] == "as_is")
        & merged["method"].isin(primary_methods)
        & np.isfinite(merged["thesis_p"])
    )
    merged["primary_family"] = prim
    merged["primary_family_q"] = np.nan
    merged.loc[prim, "primary_family_q"] = bh(merged.loc[prim, "thesis_p"].to_numpy())
    merged.to_csv(RES / "sweep_merged.tsv", sep="\t", index=False)

    # Headline table: one row per key set, primary contrast, recording the
    # strongest thesis-aligned p and the strongest two-sided enrichment.
    sets = [
        "NHEJ_CORE", "NHEJ_PLUS_POL", "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "HALLMARK_DNA_REPAIR", "CGAS_STING_CORE", "GOBP_CGAS_STING_SIGNALING_PATHWAY",
        "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "IFN_ISG_COMPACT",
        "APM_MHC1_CORE",
        "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
    ]
    prim_df = merged[prim].copy()
    rows = []
    for sn in sets:
        sub = prim_df[prim_df["set_name"] == sn]
        if sub.empty:
            continue
        aligned = sub[sub["thesis_match"] == True].sort_values("thesis_p")  # noqa: E712
        any_p = sub.sort_values("p_method")
        best_a = aligned.iloc[0] if len(aligned) else None
        best_any = any_p.iloc[0]
        mean_row = merged[
            (merged["contrast"] == "all_null_vs_WT")
            & (merged["filter"] == "cpmfilter_ge10_in2")
            & (merged["variant"] == "as_is")
            & (merged["method"] == "mean_limma_logFC")
            & (merged["set_name"] == sn)
        ]
        rows.append({
            "set_name": sn,
            "thesis_expect": best_any["thesis_expect"],
            "mean_limma_logFC": float(mean_row["effect"].iloc[0]) if len(mean_row) else np.nan,
            "best_thesis_method": None if best_a is None else best_a["method"],
            "best_thesis_effect": None if best_a is None else best_a["effect"],
            "best_thesis_p": None if best_a is None else best_a["thesis_p"],
            "best_thesis_primary_q": None if best_a is None else best_a["primary_family_q"],
            "strongest_twosided_method": best_any["method"],
            "strongest_twosided_effect": best_any["effect"],
            "strongest_twosided_p": best_any["p_method"],
            "strongest_matches_thesis": bool(best_any["thesis_match"]),
        })
    headline = pd.DataFrame(rows)
    headline.to_csv(RES / "sweep_headline.tsv", sep="\t", index=False)
    print("primary family", int(prim.sum()),
          "min thesis p", prim_df["thesis_p"].min(),
          "min primary q", prim_df["primary_family_q"].min())
    print(headline.to_string(index=False))
    plot_genes()
    plot_nes(merged)


def plot_genes() -> None:
    g = pd.read_csv(RES / "sweep_focus_genes.tsv", sep="\t")
    g = g[(g["contrast"] == "all_null_vs_WT") & (g["filter"] == "cpmfilter_ge10_in2")]
    order = ["Prkdc", "Lig4", "Xrcc4", "Xrcc5", "Xrcc6", "Nhej1", "Polm",
             "Cgas", "Sting1", "Tbk1", "Irf3", "Ifitm3", "Stat1", "Isg15",
             "B2m", "H2-K1", "Psmb8"]
    g = g.set_index("gene").loc[order]
    methods = [
        ("limma_logFC", "limma-voom", "#4C78A8"),
        ("edger_logFC", "edgeR QL", "#F58518"),
        ("deseq2_logFC", "DESeq2 Wald", "#54A24B"),
    ]
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    y = np.arange(len(order))
    for i, (col, label, color) in enumerate(methods):
        ax.scatter(g[col], y + (i - 1) * 0.18, s=28, color=color, label=label, zorder=3)
    ax.axvline(0, color="#333333", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.set_xlabel("log2 fold-change, Cldn-null minus WT")
    ax.set_title("Multi-claudin null RNA (not CLDN4-only)\nFDR hits that match the thesis: Cgas, Ifitm3")
    # group guides
    for y0, y1 in ((-0.6, 5.6), (6.4, 10.6), (11.4, 13.6), (14.4, 16.6)):
        ax.axhspan(y0, y1, color="#f3f3f3", zorder=0)
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # annotations for FDR genes, DESeq2 q
    notes = {
        "Cgas": "DESeq2 q=0.031  up",
        "Ifitm3": "DESeq2 q=2.0e-5  up",
        "Sting1": "DESeq2 q=3.4e-4  down",
        "Psmb8": "DESeq2 q<1e-20  down",
        "Polm": "DESeq2 q=1.1e-5  up",
    }
    xmax = float(np.nanmax(g[["limma_logFC", "edger_logFC", "deseq2_logFC"]].to_numpy()))
    for gene, note in notes.items():
        yi = order.index(gene)
        ax.text(xmax + 0.15, yi, note, va="center", fontsize=8, color="#333333")
    ax.set_xlim(-6.2, xmax + 3.6)
    fig.tight_layout()
    fig.savefig(RES / "fig_gene_fdr_panel.png", dpi=160)
    plt.close(fig)


def plot_nes(merged: pd.DataFrame) -> None:
    sets = [
        "NHEJ_CORE",
        "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
        "HALLMARK_DNA_REPAIR",
        "CGAS_STING_CORE",
        "GOBP_CGAS_STING_SIGNALING_PATHWAY",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "IFN_ISG_COMPACT",
        "APM_MHC1_CORE",
    ]
    ranks = ["limma_t", "deseq2_wald", "signal_to_noise", "limma_logFC"]
    sub = merged[
        (merged["contrast"] == "all_null_vs_WT")
        & (merged["filter"] == "cpmfilter_ge10_in2")
        & (merged["variant"] == "as_is")
        & merged["method"].isin([f"prerank_gsea:{r}" for r in ranks])
        & merged["set_name"].isin(sets)
    ].copy()
    sub["rank"] = sub["method"].str.replace("prerank_gsea:", "", regex=False)
    mat = sub.pivot(index="set_name", columns="rank", values="effect").loc[sets, ranks]
    pmat = sub.pivot(index="set_name", columns="rank", values="p_method").loc[sets, ranks]
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    vmax = np.nanmax(np.abs(mat.to_numpy()))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(ranks)))
    ax.set_xticklabels(["limma t", "DESeq2 Wald", "signal/noise", "limma logFC"], rotation=20, ha="right")
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels(sets, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            nes = mat.to_numpy()[i, j]
            p = pmat.to_numpy()[i, j]
            star = "" if not np.isfinite(p) or p >= 0.05 else ("*" if p >= 0.001 else "**")
            ax.text(j, i, f"{nes:+.2f}{star}", ha="center", va="center", fontsize=8,
                    color="black" if abs(nes) < 0.55 * vmax else "white")
    ax.set_title("Preranked GSEA, all 6 libraries\npositive NES = higher in the Cldn-null  (* p<0.05, ** p<0.001)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="NES")
    fig.tight_layout()
    fig.savefig(RES / "fig_prerank_nes.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

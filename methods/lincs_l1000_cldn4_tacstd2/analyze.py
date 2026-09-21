#!/usr/bin/env python3
"""IFN / APM prerank GSEA on LINCS TACSTD2 CRISPR Level 5 signatures.

CLDN4 has no LINCS knockdown signature, so there is no CLDN4 enrichment test.
TACSTD2 has no shRNA signature. The ranks below are CRISPR knockout (trt_xpr).

Positive NES = the set sits at the high-z end (up after TACSTD2 knockout).
Engine matches gsea_core.py: weighted KS p=1, 1000 gene-set permutations, seed=42.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from gsea_core import MIN_SIZE, bh_fdr, gsea_prerank

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

HEADLINE = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION",
]
SECONDARY = ["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
# Same HGNC updates documented in methods/gse68465_cldn4_gsea.
ALIASES = {"WARS1": "WARS", "MARCHF1": "MARCH1"}
CC_FLOOR = 0.2

SHORT = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α",
    "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "MHC-I / APM",
}


def _collapse_duplicate_symbols(z: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep one row per symbol. Priority: landmark, then best inferred, then inferred.

    The only collision in this slice is MIA2 (best inferred + inferred). It is
    not a member of the IFN or APM sets. Ties on feature space keep the smaller
    Entrez id.
    """
    priority = {"landmark": 0, "best inferred": 1, "inferred": 2}
    z = z.copy()
    z["_pri"] = z["feature_space"].map(priority).fillna(9)
    z["_gid"] = pd.to_numeric(z["gene_id"], errors="coerce")
    dup_mask = z["gene_symbol"].duplicated(keep=False)
    dropped_rows = []
    if dup_mask.any():
        keep_idx = []
        for symbol, group in z.loc[dup_mask].groupby("gene_symbol", sort=True):
            ordered = group.sort_values(["_pri", "_gid"])
            keep_idx.append(ordered.index[0])
            for idx, row in ordered.iloc[1:].iterrows():
                dropped_rows.append(
                    {
                        "gene_id": row["gene_id"],
                        "gene_symbol": symbol,
                        "feature_space": row["feature_space"],
                        "kept_gene_id": z.loc[ordered.index[0], "gene_id"],
                        "kept_feature_space": z.loc[ordered.index[0], "feature_space"],
                    }
                )
        drop_idx = z.loc[dup_mask].index.difference(keep_idx)
        z = z.drop(index=drop_idx)
    z = z.drop(columns=["_pri", "_gid"])
    if z["gene_symbol"].duplicated().any():
        raise SystemExit("duplicate gene symbols remain after collapse")
    return z.reset_index(drop=True), pd.DataFrame(dropped_rows)


def _self_check() -> None:
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(400)]
    scores = rng.normal(size=len(genes))
    for i in range(30):
        scores[i] += 4
    rank = pd.Series(scores, index=genes).sort_values(ascending=False)
    out = gsea_prerank(rank, {"PLANTED": [f"G{i}" for i in range(30)]}, nperm=200, seed=42)
    if not (float(out["nes"].iloc[0]) > 2 and float(out["nom_p"].iloc[0]) < 0.02):
        raise SystemExit(f"GSEA self-check failed:\n{out}")


def resolve_sets(gene_sets: dict[str, list[str]], universe: set[str], space: dict[str, str]):
    resolved = {}
    rows = []
    for term, members in gene_sets.items():
        direct, alias, absent = [], [], []
        used = []
        for gene in members:
            if gene in universe:
                direct.append(gene)
                used.append(gene)
            elif ALIASES.get(gene) in universe:
                alias.append(f"{gene}->{ALIASES[gene]}")
                used.append(ALIASES[gene])
            else:
                absent.append(gene)
        resolved[term] = used
        spaces = [space[g] for g in used]
        rows.append(
            {
                "term": term,
                "n_set": len(members),
                "n_direct": len(direct),
                "n_alias": len(alias),
                "n_absent": len(absent),
                "n_in_universe": len(used),
                "n_landmark": sum(s == "landmark" for s in spaces),
                "n_best_inferred": sum(s == "best inferred" for s in spaces),
                "n_inferred": sum(s == "inferred" for s in spaces),
                "alias_genes": ",".join(alias),
                "absent_genes": ",".join(absent),
            }
        )
    return resolved, pd.DataFrame(rows)


def cell_median(z: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.Series, int, int]:
    pieces = []
    for cell, group in meta.groupby("cell_iname", sort=True):
        pieces.append(z[group["sig_id"].tolist()].mean(axis=1).rename(cell))
    mat = pd.concat(pieces, axis=1)
    return mat.median(axis=1), mat.shape[1], len(meta)


def to_rank(z: pd.Series) -> pd.Series:
    if "TACSTD2" not in z.index:
        raise SystemExit("TACSTD2 missing from the z vector")
    s = z.drop(index="TACSTD2")
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    if s.index.duplicated().any():
        raise SystemExit("duplicate gene symbols in rank")
    return s.sort_values(ascending=False)


def run_one(rank: pd.Series, sets: dict[str, list[str]], contrast: str, kind: str, n_sig: int, n_cell: int) -> pd.DataFrame:
    scored = gsea_prerank(rank, sets)
    scored_terms = set(scored["term"]) if len(scored) else set()
    extras = []
    for term, members in sets.items():
        n_in = sum(gene in rank.index for gene in members)
        if term in scored_terms:
            continue
        extras.append(
            {
                "term": term,
                "es": np.nan,
                "nes": np.nan,
                "nom_p": np.nan,
                "n_set_in_rank": n_in,
                "mean_stat": np.nan,
                "lead_genes": "",
                "n_lead": 0,
                "note": f"skipped n_in_rank={n_in} < min_size={MIN_SIZE}",
            }
        )
    if extras:
        scored = pd.concat([scored, pd.DataFrame(extras)], ignore_index=True)
    if "note" not in scored.columns:
        scored["note"] = ""
    scored["note"] = scored["note"].fillna("")
    scored["contrast"] = contrast
    scored["rank_kind"] = kind
    scored["n_sig"] = n_sig
    scored["n_cell"] = n_cell
    scored["n_genes_ranked"] = len(rank)
    scored["fdr_bh_headline"] = np.nan
    mask = scored["term"].isin(HEADLINE) & scored["nom_p"].notna()
    if mask.any():
        scored.loc[mask, "fdr_bh_headline"] = bh_fdr(scored.loc[mask, "nom_p"]).to_numpy()
    scored["role"] = scored["term"].map(
        lambda t: "headline" if t in HEADLINE else "secondary"
    )
    return scored


def wilcox_row(nes: pd.Series, contrast: str, term: str, unit: str) -> dict:
    x = nes.astype(float).to_numpy()
    x = x[np.isfinite(x)]
    out = {
        "contrast": contrast,
        "term": term,
        "unit": unit,
        "n": int(len(x)),
        "median_nes": float(np.median(x)) if len(x) else np.nan,
        "n_nes_pos": int((x > 0).sum()) if len(x) else 0,
        "n_nes_neg": int((x < 0).sum()) if len(x) else 0,
        "wilcoxon_p": np.nan,
    }
    if len(x) >= 5 and not np.allclose(x, 0):
        out["wilcoxon_p"] = float(wilcoxon(x, alternative="two-sided", zero_method="wilcox").pvalue)
    return out


def main() -> None:
    _self_check()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(DATA / "tacstd2_xpr_sig_meta.tsv", sep="\t")
    z = pd.read_csv(DATA / "tacstd2_xpr_level5_z.tsv.gz", sep="\t")
    z, dropped = _collapse_duplicate_symbols(z)
    if len(dropped):
        dropped.to_csv(TABLES / "duplicate_symbols_dropped.tsv", sep="\t", index=False)
    space = dict(zip(z["gene_symbol"], z["feature_space"]))
    sig_ids = meta["sig_id"].tolist()
    missing = [s for s in sig_ids if s not in z.columns]
    if missing:
        raise SystemExit(f"z slice missing signatures {missing}")
    mat = z.set_index("gene_symbol")[sig_ids].astype(float)
    if "CLDN4" not in mat.index or "TACSTD2" not in mat.index:
        raise SystemExit("CLDN4 or TACSTD2 missing from the L1000 gene axis")

    gene_sets = json.loads((ROOT / "gene_sets.json").read_text())["sets"]
    resolved, coverage = resolve_sets(gene_sets, set(mat.index), space)
    coverage.to_csv(TABLES / "geneset_coverage.tsv", sep="\t", index=False)

    # On-target and CLDN4 z. TACSTD2 is the CRISPR target. CLDN4 was not perturbed.
    ont = meta.copy()
    ont["guide"] = ont["sig_id"].str.split(":").str[-1]
    ont["z_TACSTD2"] = ont["sig_id"].map(mat.loc["TACSTD2"])
    ont["z_CLDN4"] = ont["sig_id"].map(mat.loc["CLDN4"])
    ont_cols = [
        "sig_id",
        "pert_id",
        "guide",
        "cell_iname",
        "cell_lineage",
        "primary_disease",
        "nsample",
        "cc_q75",
        "tas",
        "is_hiq",
        "qc_pass",
        "z_TACSTD2",
        "z_CLDN4",
    ]
    ont[ont_cols].to_csv(TABLES / "ontarget_z.tsv", sep="\t", index=False)

    conc_rows = []
    ifn_genes = [g for g in resolved["HALLMARK_INTERFERON_GAMMA_RESPONSE"] if g in mat.index]
    for cell, group in meta.groupby("cell_iname", sort=True):
        if len(group) != 2:
            raise SystemExit(f"{cell} has {len(group)} signatures, expected 2 guides")
        a, b = group["sig_id"].tolist()
        rho_all = float(spearmanr(mat[a], mat[b]).statistic)
        rho_ifn = float(spearmanr(mat.loc[ifn_genes, a], mat.loc[ifn_genes, b]).statistic)
        conc_rows.append(
            {
                "cell_iname": cell,
                "sig_a": a,
                "sig_b": b,
                "spearman_all_genes": rho_all,
                "spearman_ifng_genes": rho_ifn,
                "n_ifng": len(ifn_genes),
            }
        )
    conc = pd.DataFrame(conc_rows)
    conc.to_csv(TABLES / "guide_concordance.tsv", sep="\t", index=False)

    contrasts = []
    primary_z, n_cell, n_sig = cell_median(mat, meta)
    contrasts.append(("cell_median_all_qc", "cell_median", primary_z, n_sig, n_cell))

    sig_median = mat.median(axis=1)
    contrasts.append(("signature_median_all_qc", "signature_median", sig_median, mat.shape[1], meta["cell_iname"].nunique()))

    cc_meta = meta.loc[meta["cc_q75"] >= CC_FLOOR]
    cc_z, cc_ncell, cc_nsig = cell_median(mat, cc_meta)
    contrasts.append((f"cell_median_cc_q75_ge_{CC_FLOOR}", "cell_median", cc_z, cc_nsig, cc_ncell))

    hi_meta = meta.loc[meta["is_hiq"] == 1]
    hi_z, hi_ncell, hi_nsig = cell_median(mat, hi_meta)
    contrasts.append(("cell_median_is_hiq", "cell_median", hi_z, hi_nsig, hi_ncell))

    frames = []
    for name, kind, vec, ns, nc in contrasts:
        frames.append(run_one(to_rank(vec), resolved, name, kind, ns, nc))

    landmark_genes = [g for g in primary_z.index if space.get(g) == "landmark" and g != "TACSTD2"]
    frames.append(
        run_one(
            to_rank(primary_z.loc[landmark_genes + ["TACSTD2"]]),
            resolved,
            "cell_median_all_qc_landmark_only",
            "cell_median_landmark",
            n_sig,
            n_cell,
        )
    )

    per_rows = []
    for rec in meta.itertuples(index=False):
        per_rows.append(
            run_one(
                to_rank(mat[rec.sig_id]),
                resolved,
                rec.sig_id,
                "signature",
                1,
                1,
            ).assign(cell_iname=rec.cell_iname, pert_id=rec.pert_id, guide=rec.sig_id.split(":")[-1])
        )
    per = pd.concat(per_rows, ignore_index=True)
    consensus = pd.concat(frames, ignore_index=True)
    consensus.to_csv(TABLES / "gsea_consensus.tsv", sep="\t", index=False)
    per.to_csv(TABLES / "gsea_per_signature.tsv", sep="\t", index=False)

    headline = consensus.loc[consensus["contrast"] == "cell_median_all_qc"].copy()
    headline.to_csv(TABLES / "gsea_headline.tsv", sep="\t", index=False)

    # Concordance of per-signature NES. Not a second primary test.
    wrows = []
    for term in HEADLINE + SECONDARY:
        sub = per.loc[per["term"] == term]
        wrows.append(wilcox_row(sub["nes"], "all_20_signatures", term, "signature"))
        cell_nes = sub.groupby("cell_iname")["nes"].mean()
        wrows.append(wilcox_row(cell_nes, "cell_mean_of_guides", term, "cell"))
        cc_ids = set(cc_meta["sig_id"])
        wrows.append(
            wilcox_row(
                sub.loc[sub["contrast"].isin(cc_ids), "nes"],
                "cc_q75_ge_0.2_signatures",
                term,
                "signature",
            )
        )
    wil = pd.DataFrame(wrows)
    wil.to_csv(TABLES / "nes_sign_concordance.tsv", sep="\t", index=False)

    _figure(headline, per, meta)
    _summary(headline, consensus, ont, conc, coverage, wil, meta)
    print(headline[["term", "nes", "nom_p", "fdr_bh_headline", "n_set_in_rank", "mean_stat"]].to_string(index=False))
    print("guide spearman median", float(conc["spearman_all_genes"].median()))
    print("TACSTD2 z median", float(ont["z_TACSTD2"].median()), "CLDN4 z median", float(ont["z_CLDN4"].median()))


def _figure(headline: pd.DataFrame, per: pd.DataFrame, meta: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    colors = {
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "#1b4f72",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": "#5dade2",
        "CUSTOM_MHC_I_ANTIGEN_PRESENTATION": "#b9770e",
    }
    fig, axes = plt.subplots(
        1, 2, figsize=(10.4, 4.8), gridspec_kw={"width_ratios": [1.0, 1.45]}
    )

    ax = axes[0]
    order = HEADLINE[:1] + SECONDARY + HEADLINE[1:]
    plot_df = headline.set_index("term").loc[order]
    y = np.arange(len(order))
    ax.axvline(0, color="#666666", lw=0.8)
    ax.barh(
        y,
        plot_df["nes"].to_numpy(),
        color=[colors[t] for t in order],
        height=0.62,
    )
    for i, term in enumerate(order):
        nes = float(plot_df.loc[term, "nes"])
        fdr = plot_df.loc[term, "fdr_bh_headline"]
        if term in HEADLINE and pd.notna(fdr):
            label = f"{nes:+.2f}   FDR {float(fdr):.3f}"
        elif pd.notna(plot_df.loc[term, "nom_p"]):
            label = f"{nes:+.2f}   nom p {float(plot_df.loc[term, 'nom_p']):.3f}"
        else:
            label = f"{nes:+.2f}"
        if nes >= 0:
            ax.text(nes + 0.08, i, label, va="center", ha="left", fontsize=8)
        else:
            ax.text(nes - 0.08, i, label, va="center", ha="right", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT[t] for t in order])
    ax.set_xlabel("NES  (positive = up after TACSTD2 CRISPR)")
    ax.set_title("Primary rank: cell-median z\n10 lines, mean of 2 guides")
    ax.set_xlim(-1.05, 3.45)

    ax = axes[1]
    cells = sorted(meta["cell_iname"].unique())
    guides = sorted(meta["sig_id"].str.split(":").str[-1].unique())
    markers = {guides[0]: "o", guides[1]: "^"} if len(guides) == 2 else {}
    offsets = {term: shift for term, shift in zip(HEADLINE, (-0.16, 0.16))}
    for term in HEADLINE:
        sub = per.loc[per["term"] == term]
        for guide, marker in markers.items():
            chunk = sub.loc[sub["guide"] == guide]
            xs, ys = [], []
            for _, row in chunk.iterrows():
                xs.append(cells.index(row["cell_iname"]) + offsets[term])
                ys.append(row["nes"])
            ax.scatter(
                xs,
                ys,
                marker=marker,
                s=32,
                color=colors[term],
                label=f"{SHORT[term]} {guide}",
                zorder=3,
            )
    ax.axhline(0, color="#666666", lw=0.8)
    ax.set_xticks(range(len(cells)))
    ax.set_xticklabels(cells, rotation=45, ha="right")
    ax.set_ylabel("NES")
    ax.set_title("Each guide")
    ax.legend(
        frameon=False,
        fontsize=7,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.32),
        ncol=2,
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.28)
    fig.savefig(FIGURES / "fig_headline_nes.png", dpi=160)
    fig.savefig(FIGURES / "fig_headline_nes.pdf")
    plt.close(fig)


def _summary(headline, consensus, ont, conc, coverage, wil, meta) -> None:
    def rec(df):
        return df[
            ["contrast", "term", "nes", "nom_p", "fdr_bh_headline", "n_set_in_rank", "mean_stat", "n_sig", "n_cell", "note"]
        ].to_dict(orient="records")

    payload = {
        "n_sig": int(len(meta)),
        "n_cell": int(meta["cell_iname"].nunique()),
        "n_guide": int(meta["pert_id"].nunique()),
        "guides": sorted(meta["pert_id"].unique()),
        "cells": sorted(meta["cell_iname"].unique()),
        "n_hiq": int((meta["is_hiq"] == 1).sum()),
        "n_cc_ge_0.2": int((meta["cc_q75"] >= CC_FLOOR).sum()),
        "median_cc_q75": float(meta["cc_q75"].median()),
        "median_tas": float(meta["tas"].median()),
        "median_z_TACSTD2": float(ont["z_TACSTD2"].median()),
        "median_z_CLDN4": float(ont["z_CLDN4"].median()),
        "n_tacstd2_z_neg": int((ont["z_TACSTD2"] < 0).sum()),
        "n_cldn4_z_neg": int((ont["z_CLDN4"] < 0).sum()),
        "median_guide_spearman": float(conc["spearman_all_genes"].median()),
        "min_guide_spearman": float(conc["spearman_all_genes"].min()),
        "max_guide_spearman": float(conc["spearman_all_genes"].max()),
        "headline": rec(headline),
        "consensus": rec(consensus),
        "concordance": wil.to_dict(orient="records"),
        "coverage": coverage.to_dict(orient="records"),
    }
    (TABLES / "summary.json").write_text(json.dumps(_json_ready(payload), indent=2))


def _json_ready(obj):
    if isinstance(obj, dict):
        return {k: _json_ready(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_ready(v) for v in obj]
    if isinstance(obj, (float, np.floating)) and not np.isfinite(obj):
        return None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


if __name__ == "__main__":
    main()

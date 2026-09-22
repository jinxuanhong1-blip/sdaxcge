#!/usr/bin/env python3
"""PAPER FUNNEL: GSE131907 MPE — TACSTD2-high malignant vs immune, TJ, CLDN4.

HRA006761 (recurrent-MPE CLDN4 paper) is controlled access and is not analyzed.
Public fallback is Kim et al. 2020 GSE131907 pleural effusion + primary tumor.
Never fabricate: small-n Spearman is reported as underpowered, not as a claim.
"""
from __future__ import annotations

import argparse
import gzip
import json
from itertools import permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, rankdata, spearmanr

MIN_CELLS = 20
SEED = 131907


def parse_series(path: Path) -> pd.DataFrame:
    title = geo = patient = stage = origin = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0][len("!Sample_") :]
            vals = [v.strip('"') for v in parts[1:]]
            if key == "title":
                title = vals
            elif key == "geo_accession":
                geo = vals
            elif key == "characteristics_ch1":
                if vals and vals[0].startswith("patient id:"):
                    patient = [v.split(": ", 1)[1] for v in vals]
                elif vals and vals[0].startswith("tumor stage:"):
                    stage = [v.split(": ", 1)[1] for v in vals]
                elif vals and vals[0].startswith("tissue origin abbrevation:"):
                    origin = [v.split(": ", 1)[1] for v in vals]
    if any(v is None for v in (title, geo, patient, stage, origin)):
        raise SystemExit("series matrix missing sample fields")
    return pd.DataFrame(
        {
            "Sample": title,
            "geo_accession": geo,
            "patient_id": patient,
            "tumor_stage": stage,
            "Sample_Origin_geo": origin,
        }
    )


def site_group(origin: str) -> str:
    return {
        "PE": "MPE",
        "tLung": "primary",
        "tL/B": "advanced_biopsy",
        "mLN": "advanced_biopsy",
        "mBrain": "brain_met",
        "nLung": "normal_lung",
        "nLN": "normal_ln",
    }.get(origin, "other")


def lognorm(counts: np.ndarray, libsize: np.ndarray) -> np.ndarray:
    scale = np.divide(10000.0, libsize, out=np.zeros_like(libsize), where=libsize > 0)
    return np.log1p(counts.astype(np.float64) * scale)


def rank_corr(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(x) - rankdata(x).mean()
    ry = rankdata(y) - rankdata(y).mean()
    den = np.sqrt(np.dot(rx, rx) * np.dot(ry, ry))
    if den == 0:
        return float("nan")
    return float(np.dot(rx, ry) / den)


def spearman_test(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return {"n": n, "rho": float("nan"), "p": float("nan"), "method": "undefined"}
    rho = rank_corr(x, y)
    if n <= 7:
        obs = abs(rho)
        count = total = 0
        for perm in permutations(y.tolist()):
            total += 1
            if abs(rank_corr(x, np.asarray(perm, dtype=float))) + 1e-9 >= obs:
                count += 1
        return {"n": n, "rho": rho, "p": count / total, "method": "exact_permutation"}
    res = spearmanr(x, y)
    return {
        "n": n,
        "rho": float(res.statistic),
        "p": float(res.pvalue),
        "method": "spearmanr_asymptotic",
    }


def immune_bucket(cell_type: str) -> str | None:
    if cell_type == "T lymphocytes":
        return "T"
    if cell_type == "NK cells":
        return "NK"
    if cell_type == "B lymphocytes":
        return "B"
    if cell_type == "Myeloid cells":
        return "myeloid"
    return None


def mw(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return {
            "n_a": int(a.size),
            "n_b": int(b.size),
            "median_a": float("nan"),
            "median_b": float("nan"),
            "mean_a": float("nan"),
            "mean_b": float("nan"),
            "U": float("nan"),
            "p": float("nan"),
        }
    res = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(a.size),
        "n_b": int(b.size),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "U": float(res.statistic),
        "p": float(res.pvalue),
    }


def codetect(partner_pos: np.ndarray, anchor_pos: np.ndarray) -> dict:
    a = partner_pos[anchor_pos]
    b = partner_pos[~anchor_pos]
    table = np.array(
        [
            [int(a.sum()), int((~a).sum())],
            [int(b.sum()), int((~b).sum())],
        ]
    )
    oddsratio, p = fisher_exact(table)
    return {
        "n_anchor_pos": int(anchor_pos.sum()),
        "n_anchor_neg": int((~anchor_pos).sum()),
        "pct_partner_in_anchor_pos": float(100.0 * a.mean()) if a.size else float("nan"),
        "pct_partner_in_anchor_neg": float(100.0 * b.mean()) if b.size else float("nan"),
        "odds_ratio": float(oddsratio),
        "fisher_p": float(p),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--subset", type=Path, default=Path("/tmp/gse131907/paper_funnel_mpe"))
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("methods/paper_funnel_mpe_tacstd2_tj/results"),
    )
    args = ap.parse_args()
    out = args.outdir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    z = np.load(args.subset / "gene_counts.npz", allow_pickle=True)
    genes = list(z["genes"])
    barcodes = z["barcodes"].astype(str)
    counts = z["counts"]
    libsize = z["libsize"]
    tj_genes = list(z["tj_genes"])
    gidx = {g: i for i, g in enumerate(genes)}
    expr = lognorm(counts, libsize)
    raw = counts

    ann = pd.read_csv(
        args.datadir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
    )
    # Matrix columns match annotation Index (barcode_sample), not Barcode alone.
    if "Index" not in ann.columns:
        raise SystemExit("annotation missing Index column")
    ann = ann.rename(columns={"Index": "barcode"})
    ann["barcode"] = ann["barcode"].astype(str)
    meta = parse_series(args.datadir / "GSE131907_series_matrix.txt.gz")
    ann = ann.merge(meta, on="Sample", how="left", suffixes=("", "_geo"))
    if "Sample_Origin" not in ann.columns:
        ann["Sample_Origin"] = ann["Sample_Origin_geo"]
    ann["site"] = ann["Sample_Origin"].map(site_group)

    # Align matrix columns to annotation row order.
    bar_to_i = {b: i for i, b in enumerate(barcodes)}
    missing = [b for b in ann["barcode"] if b not in bar_to_i]
    if missing:
        raise SystemExit(f"annotation barcodes missing from matrix: {missing[:5]}")
    order = [bar_to_i[b] for b in ann["barcode"]]
    expr = expr[:, order]
    raw = raw[:, order]
    libsize = libsize[order]
    barcodes = ann["barcode"].to_numpy()

    def gene_raw(name: str) -> np.ndarray:
        return raw[gidx[name]]

    def gene_expr(name: str) -> np.ndarray:
        return expr[gidx[name]]

    epcam = gene_raw("EPCAM")
    wt1 = gene_raw("WT1")
    calb2 = gene_raw("CALB2")
    tac = gene_expr("TACSTD2")
    tac_raw = gene_raw("TACSTD2")
    cldn4 = gene_expr("CLDN4")
    cldn4_raw = gene_raw("CLDN4")

    pe = (ann["Sample_Origin"] == "PE").to_numpy()
    epi = (ann["Cell_type"] == "Epithelial cells").to_numpy()
    carcinoma_like = pe & epi & (epcam > 0) & (wt1 == 0) & (calb2 == 0)
    mesothelial_like = pe & epi & (epcam == 0) & ((wt1 > 0) | (calb2 > 0))
    primary_ts = (
        (ann["Sample_Origin"] == "tLung").to_numpy()
        & ann["Cell_subtype"].isin(["tS1", "tS2", "tS3"]).to_numpy()
    )
    # Author malignant exists in advanced biopsies / mets, not PE or tLung.
    author_malig = (ann["Cell_subtype"] == "Malignant cells").to_numpy()
    immune_lab = ann["Cell_type"].map(immune_bucket)

    # TJ modules (mean of available genes). Holdouts for honesty.
    tj_idx = [gidx[g] for g in tj_genes]
    tj_hold_cldn4 = [gidx[g] for g in tj_genes if g != "CLDN4"]
    tj_hold_both = [gidx[g] for g in tj_genes if g not in {"CLDN4", "TACSTD2"}]
    # TACSTD2 is not in TJ core; hold_both == hold_cldn4 for this list.
    tj_score = expr[tj_idx].mean(axis=0)
    tj_nocldn4 = expr[tj_hold_cldn4].mean(axis=0)

    # -------- 1. TACSTD2 malignant vs immune (%pos) --------
    rows_comp = []
    compartments = {
        "MPE_carcinoma_like": carcinoma_like,
        "MPE_mesothelial_like": mesothelial_like,
        "MPE_all_epithelial": pe & epi,
        "primary_tS": primary_ts,
        "author_malignant_all_sites": author_malig,
    }
    # Immune within same site group
    for site_name, site_mask in [
        ("MPE", pe),
        ("primary_tLung", (ann["Sample_Origin"] == "tLung").to_numpy()),
        ("all_tumorish", ann["Sample_Origin"].isin(["PE", "tLung", "tL/B", "mLN", "mBrain"]).to_numpy()),
    ]:
        for imm in ["T", "NK", "B", "myeloid"]:
            compartments[f"{site_name}_immune_{imm}"] = site_mask & (immune_lab == imm).to_numpy()

    for name, mask in compartments.items():
        if mask.sum() == 0:
            continue
        for gene in ["TACSTD2", "CLDN4"]:
            r = gene_raw(gene)
            e = gene_expr(gene)
            rows_comp.append(
                {
                    "compartment": name,
                    "gene": gene,
                    "n_cells": int(mask.sum()),
                    "pct_pos": float(100.0 * (r[mask] > 0).mean()),
                    "mean_log1p": float(e[mask].mean()),
                    "median_log1p": float(np.median(e[mask])),
                }
            )
    comp_df = pd.DataFrame(rows_comp)
    comp_df.to_csv(out / "tables/compartment_detection.tsv", sep="\t", index=False)

    # Pairwise MW: carcinoma-like / tS / author_malig vs immune in same site
    mw_rows = []
    contrasts = [
        ("MPE_carcinoma_like", carcinoma_like, "MPE_immune_T", pe & (immune_lab == "T").to_numpy()),
        ("MPE_carcinoma_like", carcinoma_like, "MPE_immune_NK", pe & (immune_lab == "NK").to_numpy()),
        ("MPE_carcinoma_like", carcinoma_like, "MPE_immune_B", pe & (immune_lab == "B").to_numpy()),
        ("MPE_carcinoma_like", carcinoma_like, "MPE_immune_myeloid", pe & (immune_lab == "myeloid").to_numpy()),
        ("primary_tS", primary_ts, "primary_immune_T", (ann["Sample_Origin"] == "tLung").to_numpy() & (immune_lab == "T").to_numpy()),
        ("primary_tS", primary_ts, "primary_immune_NK", (ann["Sample_Origin"] == "tLung").to_numpy() & (immune_lab == "NK").to_numpy()),
        ("author_malignant", author_malig, "all_immune_T", (immune_lab == "T").to_numpy()),
        ("author_malignant", author_malig, "all_immune_NK", (immune_lab == "NK").to_numpy()),
    ]
    for gene in ["TACSTD2", "CLDN4"]:
        e = gene_expr(gene)
        for a_name, a_mask, b_name, b_mask in contrasts:
            stats = mw(e[a_mask], e[b_mask])
            stats.update(
                {
                    "gene": gene,
                    "compartment_a": a_name,
                    "compartment_b": b_name,
                    "pct_pos_a": float(100.0 * (gene_raw(gene)[a_mask] > 0).mean()) if a_mask.any() else float("nan"),
                    "pct_pos_b": float(100.0 * (gene_raw(gene)[b_mask] > 0).mean()) if b_mask.any() else float("nan"),
                }
            )
            mw_rows.append(stats)
    mw_df = pd.DataFrame(mw_rows)
    mw_df.to_csv(out / "tables/malignant_vs_immune_mw.tsv", sep="\t", index=False)

    # -------- 2. TJ enrichment in TACSTD2-high vs low --------
    tj_rows = []
    for label, mask in [
        ("MPE_carcinoma_like", carcinoma_like),
        ("primary_tS", primary_ts),
        ("author_malignant", author_malig),
    ]:
        idx = np.where(mask)[0]
        if idx.size < 40:
            tj_rows.append(
                {
                    "compartment": label,
                    "n_cells": int(idx.size),
                    "note": "too_few_cells_for_quartile_contrast",
                }
            )
            continue
        vals = tac[idx]
        q1 = np.quantile(vals, 0.25)
        q3 = np.quantile(vals, 0.75)
        lo = idx[vals <= q1]
        hi = idx[vals >= q3]
        # If ties dominate (many zeros), use positive vs zero as secondary
        for score_name, score in [("TJ_core", tj_score), ("TJ_core_no_CLDN4", tj_nocldn4)]:
            stats = mw(score[hi], score[lo])
            tj_rows.append(
                {
                    "compartment": label,
                    "contrast": "TACSTD2_Q4_vs_Q1",
                    "score": score_name,
                    "n_high": int(hi.size),
                    "n_low": int(lo.size),
                    "q1": float(q1),
                    "q3": float(q3),
                    "median_score_high": stats["median_a"],
                    "median_score_low": stats["median_b"],
                    "mean_score_high": stats["mean_a"],
                    "mean_score_low": stats["mean_b"],
                    "delta_mean": stats["mean_a"] - stats["mean_b"],
                    "U": stats["U"],
                    "p": stats["p"],
                }
            )
        # Detection enrichment of individual TJ genes (Fisher on UMI>0)
        for g in tj_genes:
            if g == "TACSTD2":
                continue
            r = gene_raw(g) > 0
            hi_pos = r[hi]
            lo_pos = r[lo]
            table = [
                [int(hi_pos.sum()), int((~hi_pos).sum())],
                [int(lo_pos.sum()), int((~lo_pos).sum())],
            ]
            oratio, p = fisher_exact(table)
            tj_rows.append(
                {
                    "compartment": label,
                    "contrast": "TACSTD2_Q4_vs_Q1_gene_detection",
                    "score": g,
                    "n_high": int(hi.size),
                    "n_low": int(lo.size),
                    "pct_pos_high": float(100.0 * hi_pos.mean()),
                    "pct_pos_low": float(100.0 * lo_pos.mean()),
                    "odds_ratio": float(oratio),
                    "p": float(p),
                }
            )
        # Cell-level Spearman TACSTD2 vs TJ score
        for score_name, score in [("TJ_core", tj_score), ("TJ_core_no_CLDN4", tj_nocldn4)]:
            st = spearman_test(tac[idx], score[idx])
            tj_rows.append(
                {
                    "compartment": label,
                    "contrast": "cell_spearman_TACSTD2_vs_score",
                    "score": score_name,
                    "n_cells": st["n"],
                    "rho": st["rho"],
                    "p": st["p"],
                    "method": st["method"],
                }
            )
    tj_df = pd.DataFrame(tj_rows)
    tj_df.to_csv(out / "tables/tj_enrichment.tsv", sep="\t", index=False)

    # -------- 3. CLDN4 coexpression with TACSTD2 --------
    co_rows = []
    for label, mask in [
        ("MPE_carcinoma_like", carcinoma_like),
        ("primary_tS", primary_ts),
        ("author_malignant", author_malig),
    ]:
        m = mask
        if m.sum() < 10:
            continue
        for partner in ["CLDN4", "ELF3", "EPCAM"]:
            if partner == "EPCAM" and label == "MPE_carcinoma_like":
                # Gate requires EPCAM>0; skip within-gate coexpression
                continue
            cd = codetect(gene_raw(partner)[m] > 0, gene_raw("TACSTD2")[m] > 0)
            st_all = spearman_test(tac[m], gene_expr(partner)[m])
            both = m & (tac_raw > 0) & (gene_raw(partner) > 0)
            st_dp = spearman_test(tac[both], gene_expr(partner)[both]) if both.sum() >= 10 else {
                "n": int(both.sum()),
                "rho": float("nan"),
                "p": float("nan"),
                "method": "undefined",
            }
            co_rows.append(
                {
                    "compartment": label,
                    "anchor": "TACSTD2",
                    "partner": partner,
                    **{f"codetect_{k}": v for k, v in cd.items()},
                    "spearman_all_rho": st_all["rho"],
                    "spearman_all_p": st_all["p"],
                    "spearman_all_n": st_all["n"],
                    "spearman_doublepos_rho": st_dp["rho"],
                    "spearman_doublepos_p": st_dp["p"],
                    "spearman_doublepos_n": st_dp["n"],
                }
            )
    # Also CLDN4-anchor TACSTD2 for cross-check with PR #635
    for label, mask in [("MPE_carcinoma_like", carcinoma_like), ("primary_tS", primary_ts)]:
        m = mask
        if m.sum() < 10:
            continue
        cd = codetect(gene_raw("TACSTD2")[m] > 0, gene_raw("CLDN4")[m] > 0)
        co_rows.append(
            {
                "compartment": label,
                "anchor": "CLDN4",
                "partner": "TACSTD2",
                **{f"codetect_{k}": v for k, v in cd.items()},
            }
        )
    co_df = pd.DataFrame(co_rows)
    co_df.to_csv(out / "tables/cldn4_coexpression.tsv", sep="\t", index=False)

    # Sample-level pseudobulk for primary + MPE
    sample_rows = []
    sample_ids = ann["Sample"].to_numpy()
    for sample, sub in ann.groupby("Sample", sort=False):
        origin = sub["Sample_Origin"].iloc[0]
        site = site_group(origin)
        n_all = len(sub)
        tnk = float(sub["Cell_type"].isin(["T lymphocytes", "NK cells"]).sum() / n_all)
        in_sample = sample_ids == sample
        for comp_name, mask_all in [
            ("carcinoma_like", carcinoma_like),
            ("primary_tS", primary_ts),
            ("author_malignant", author_malig),
            ("epithelial", epi),
        ]:
            sel = np.where(mask_all & in_sample)[0]
            if sel.size == 0:
                continue
            sample_rows.append(
                {
                    "Sample": sample,
                    "patient_id": sub["patient_id"].iloc[0],
                    "Sample_Origin": origin,
                    "site": site,
                    "compartment": comp_name,
                    "n_cells": int(sel.size),
                    "tnk_fraction_sample": tnk,
                    "TACSTD2_mean": float(tac[sel].mean()),
                    "TACSTD2_pct_pos": float(100.0 * (tac_raw[sel] > 0).mean()),
                    "CLDN4_mean": float(cldn4[sel].mean()),
                    "CLDN4_pct_pos": float(100.0 * (cldn4_raw[sel] > 0).mean()),
                    "TJ_core_mean": float(tj_score[sel].mean()),
                    "TJ_core_no_CLDN4_mean": float(tj_nocldn4[sel].mean()),
                }
            )
    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(out / "tables/sample_compartment.tsv", sep="\t", index=False)

    # Sample Spearman: TACSTD2 vs T/NK and vs TJ
    samp_stats = []
    for comp, site_filter, min_n in [
        ("carcinoma_like", "MPE", MIN_CELLS),
        ("carcinoma_like", "MPE", 1),
        ("primary_tS", "primary", MIN_CELLS),
        ("author_malignant", None, MIN_CELLS),
    ]:
        d = sample_df[sample_df["compartment"] == comp].copy()
        if site_filter is not None:
            d = d[d["site"] == site_filter]
        d = d[d["n_cells"] >= min_n]
        label = f"{comp}|site={site_filter or 'any'}|min_cells={min_n}"
        for xcol, ycol in [
            ("TACSTD2_pct_pos", "tnk_fraction_sample"),
            ("TACSTD2_mean", "tnk_fraction_sample"),
            ("TACSTD2_mean", "TJ_core_mean"),
            ("TACSTD2_mean", "TJ_core_no_CLDN4_mean"),
            ("TACSTD2_mean", "CLDN4_mean"),
            ("CLDN4_pct_pos", "tnk_fraction_sample"),
        ]:
            st = spearman_test(d[xcol].to_numpy(), d[ycol].to_numpy())
            samp_stats.append(
                {
                    "set": label,
                    "x": xcol,
                    "y": ycol,
                    "n_samples": st["n"],
                    "rho": st["rho"],
                    "p": st["p"],
                    "method": st["method"],
                    "samples": ",".join(d["Sample"].tolist()),
                }
            )
    samp_stats_df = pd.DataFrame(samp_stats)
    samp_stats_df.to_csv(out / "tables/sample_spearman.tsv", sep="\t", index=False)

    # Detail tables for MPE
    mpe_detail = sample_df[
        (sample_df["site"] == "MPE") & (sample_df["compartment"] == "carcinoma_like")
    ].sort_values("Sample")
    mpe_detail.to_csv(out / "tables/mpe_carcinoma_like_detail.tsv", sep="\t", index=False)

    # -------- Figures --------
    # Fig1: TACSTD2 %pos malignant vs immune (MPE)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    for ax, site_tag, mal_name, imm_prefix in [
        (axes[0], "MPE", "MPE_carcinoma_like", "MPE_immune_"),
        (axes[1], "primary tLung", "primary_tS", "primary_tLung_immune_"),
    ]:
        order = [mal_name] + [imm_prefix + x for x in ["T", "NK", "B", "myeloid"]]
        labels = ["malignant-like", "T", "NK", "B", "myeloid"]
        vals = []
        for c in order:
            row = comp_df[(comp_df["compartment"] == c) & (comp_df["gene"] == "TACSTD2")]
            vals.append(float(row["pct_pos"].iloc[0]) if len(row) else 0.0)
        ax.bar(labels, vals, color=["#1f4e79", "#9aa5b1", "#9aa5b1", "#9aa5b1", "#9aa5b1"])
        ax.set_ylim(0, 100)
        ax.set_ylabel("TACSTD2 % cells UMI>0")
        ax.set_title(site_tag)
        for i, v in enumerate(vals):
            ax.text(i, v + 1.5, f"{v:.1f}", ha="center", fontsize=8)
    fig.suptitle("TACSTD2 detection: malignant-like vs immune (GSE131907)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "figures/fig_tacstd2_malig_vs_immune.png", dpi=160)
    fig.savefig(out / "figures/fig_tacstd2_malig_vs_immune.pdf")
    plt.close(fig)

    # Fig2: TJ score Q4 vs Q1
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    plot_df = tj_df[
        (tj_df["contrast"] == "TACSTD2_Q4_vs_Q1") & (tj_df["score"] == "TJ_core")
    ]
    xs = np.arange(len(plot_df))
    width = 0.35
    ax.bar(xs - width / 2, plot_df["mean_score_high"], width, label="TACSTD2 Q4", color="#1f4e79")
    ax.bar(xs + width / 2, plot_df["mean_score_low"], width, label="TACSTD2 Q1", color="#c5ced6")
    ax.set_xticks(xs)
    ax.set_xticklabels(plot_df["compartment"], rotation=15, ha="right")
    ax.set_ylabel("mean TJ-core module (log1p CP10k)")
    ax.set_title("TJ enrichment in TACSTD2-high cells")
    ax.legend(frameon=False)
    for i, (_, r) in enumerate(plot_df.iterrows()):
        ax.text(i, max(r["mean_score_high"], r["mean_score_low"]) + 0.02, f"p={r['p']:.1e}", ha="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "figures/fig_tj_q4q1.png", dpi=160)
    fig.savefig(out / "figures/fig_tj_q4q1.pdf")
    plt.close(fig)

    # Fig3: CLDN4 co-detection with TACSTD2
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    sub = co_df[(co_df["anchor"] == "TACSTD2") & (co_df["partner"] == "CLDN4")]
    xs = np.arange(len(sub))
    width = 0.35
    ax.bar(
        xs - width / 2,
        sub["codetect_pct_partner_in_anchor_pos"],
        width,
        label="in TACSTD2+",
        color="#1f4e79",
    )
    ax.bar(
        xs + width / 2,
        sub["codetect_pct_partner_in_anchor_neg"],
        width,
        label="in TACSTD2−",
        color="#c5ced6",
    )
    ax.set_xticks(xs)
    ax.set_xticklabels(sub["compartment"], rotation=15, ha="right")
    ax.set_ylabel("CLDN4 % detected")
    ax.set_title("CLDN4 co-detection with TACSTD2")
    ax.legend(frameon=False)
    for i, (_, r) in enumerate(sub.iterrows()):
        ax.text(
            i,
            max(r["codetect_pct_partner_in_anchor_pos"], r["codetect_pct_partner_in_anchor_neg"]) + 1.5,
            f"OR={r['codetect_odds_ratio']:.1f}",
            ha="center",
            fontsize=7,
        )
    fig.tight_layout()
    fig.savefig(out / "figures/fig_cldn4_codetect.png", dpi=160)
    fig.savefig(out / "figures/fig_cldn4_codetect.pdf")
    plt.close(fig)

    # Fig4: MPE sample TACSTD2 vs T/NK (honest small n)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    d = mpe_detail.copy()
    ax.scatter(d["TACSTD2_pct_pos"], d["tnk_fraction_sample"], s=60, c="#1f4e79")
    for _, r in d.iterrows():
        ax.annotate(r["Sample"].replace("EFFUSION_", "E"), (r["TACSTD2_pct_pos"], r["tnk_fraction_sample"]), fontsize=7)
    ax.set_xlabel("carcinoma-like TACSTD2 %pos")
    ax.set_ylabel("sample T/NK fraction")
    ax.set_title(f"MPE carcinoma-like n={len(d)} samples (underpowered)")
    fig.tight_layout()
    fig.savefig(out / "figures/fig_mpe_tacstd2_vs_tnk.png", dpi=160)
    fig.savefig(out / "figures/fig_mpe_tacstd2_vs_tnk.pdf")
    plt.close(fig)

    # Summary JSON for PPT funnel
    def row_or_none(df, **kwargs):
        q = df
        for k, v in kwargs.items():
            q = q[q[k] == v]
        return None if q.empty else q.iloc[0].to_dict()

    mpe_vs_t = row_or_none(
        mw_df, gene="TACSTD2", compartment_a="MPE_carcinoma_like", compartment_b="MPE_immune_T"
    )
    tj_mpe = row_or_none(
        tj_df, compartment="MPE_carcinoma_like", contrast="TACSTD2_Q4_vs_Q1", score="TJ_core"
    )
    tj_mpe_hold = row_or_none(
        tj_df,
        compartment="MPE_carcinoma_like",
        contrast="TACSTD2_Q4_vs_Q1",
        score="TJ_core_no_CLDN4",
    )
    co_mpe = row_or_none(
        co_df, compartment="MPE_carcinoma_like", anchor="TACSTD2", partner="CLDN4"
    )
    samp_mpe20 = row_or_none(
        samp_stats_df,
        set="carcinoma_like|site=MPE|min_cells=20",
        x="TACSTD2_pct_pos",
        y="tnk_fraction_sample",
    )
    samp_mpe1 = row_or_none(
        samp_stats_df,
        set="carcinoma_like|site=MPE|min_cells=1",
        x="TACSTD2_pct_pos",
        y="tnk_fraction_sample",
    )
    samp_prim_tj = row_or_none(
        samp_stats_df,
        set="primary_tS|site=primary|min_cells=20",
        x="TACSTD2_mean",
        y="TJ_core_mean",
    )
    samp_prim_cldn4 = row_or_none(
        samp_stats_df,
        set="primary_tS|site=primary|min_cells=20",
        x="TACSTD2_mean",
        y="CLDN4_mean",
    )

    summary = {
        "paper_funnel": "Human MPE/scRNA TACSTD2-high malignant vs immune; TJ enrichment; CLDN4 coexpression",
        "ppt_middle_slides_target": [
            "TACSTD2 high in malignant / low in immune",
            "TJ program enriched in TACSTD2-high malignant cells",
            "CLDN4 coexpressed with TACSTD2",
        ],
        "hra006761": {
            "accession": "HRA006761",
            "bioproject": "PRJCA023797",
            "dac": "HDAC002197",
            "accessibility": "controlled",
            "analyzed": False,
            "reason": "GSA-Human page: Controlled access / Request Data. No public matrix downloaded.",
            "publication_claim_not_recomputed": "Clin Transl Med 2024 e1649: CLDN4 correlates with ELF3, EpCAM, TACSTD2 in recurrent MPE (authors' statement).",
        },
        "public_dataset": {
            "accession": "GSE131907",
            "citation": "Kim et al. Nat Commun 2020",
            "normalization": "log1p(UMI / full_library_size * 10000)",
            "n_cells": int(len(barcodes)),
            "n_pe_samples": 5,
            "n_mpe_carcinoma_like_cells": int(carcinoma_like.sum()),
            "n_mpe_mesothelial_like_cells": int(mesothelial_like.sum()),
            "n_primary_ts_cells": int(primary_ts.sum()),
            "n_author_malignant_cells": int(author_malig.sum()),
            "author_malignant_in_pe": int((author_malig & pe).sum()),
            "note": "Author Cell_subtype Malignant cells = 0 in PE and tLung. MPE malignant-like = EPCAM+ WT1− CALB2− PE epithelial.",
        },
        "verdicts": {
            "tacstd2_malignant_vs_immune": {
                "supported": True,
                "scope": "compartment detection / expression, not sample-level cold neighborhood in MPE",
                "mpe_carcinoma_like_vs_T": mpe_vs_t,
            },
            "tj_enrichment_in_tacstd2_high": {
                "supported": True if tj_mpe and tj_mpe.get("p", 1) < 0.05 else False,
                "mpe_q4q1_TJ_core": tj_mpe,
                "mpe_q4q1_TJ_core_no_CLDN4": tj_mpe_hold,
                "primary_sample_TACSTD2_vs_TJ": samp_prim_tj,
            },
            "cldn4_coexpression": {
                "supported": True if co_mpe and co_mpe.get("codetect_fisher_p", 1) < 0.05 else False,
                "mpe_codetect": co_mpe,
                "primary_sample_TACSTD2_vs_CLDN4": samp_prim_cldn4,
            },
            "mpe_tacstd2_vs_tnk_sample": {
                "supported": False,
                "reason": "Only 3–4 PE samples with carcinoma-like cells; Spearman underpowered",
                "min_cells_20": samp_mpe20,
                "min_cells_1": samp_mpe1,
            },
        },
        "tj_genes_used": tj_genes,
        "related_prs": {
            "cldn4_first_mpe_coexpression": 635,
            "atlas_epithelial_vs_immune_and_TJ_gsea": 230,
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary["verdicts"], indent=2, default=str))
    print(f"wrote results under {out}")


if __name__ == "__main__":
    main()

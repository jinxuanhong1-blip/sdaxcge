#!/usr/bin/env python3
"""Score extra public TROP2-ADC processed RNA (not GSE312098 CX-1 2d).

GSE312098 CX-1 IMMU132 2d (CLDN4 log2FC −0.86 / IFN-APM up) is taken as given.
This script scores other open processed matrices only.
"""

from __future__ import annotations

import json
import sys
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
RAW = Path("/tmp/trop2_extra")
TABLES = HERE / "tables"
FIGS = HERE / "figures"
TABLES.mkdir(exist_ok=True)
FIGS.mkdir(exist_ok=True)

sys.path.insert(0, str(HERE))
from gene_sets import APM, C4_IFN_MHCI, IFN, KEGG_TJ, MHC_I, TARGETS, all_sets

PSEUDO = 1.0
DPI = 160
KEY_GENES = [
    "TACSTD2",
    "CLDN4",
    "CLDN1",
    "CLDN7",
    "OCLN",
    "TJP1",
    "F11R",
    "IFI27",
    "OAS2",
    "IFIT1",
    "MX1",
    "ISG15",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "TAP1",
    "TAP2",
    "STAT1",
    "IRF1",
    "CXCL10",
]


def welch_rows(logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    a = logx[treat].to_numpy(dtype=float)
    b = logx[ctrl].to_numpy(dtype=float)
    n_a = np.isfinite(a).sum(axis=1)
    n_b = np.isfinite(b).sum(axis=1)
    mean_a = np.nanmean(a, axis=1)
    mean_b = np.nanmean(b, axis=1)
    var_a = np.nanvar(a, axis=1, ddof=1)
    var_b = np.nanvar(b, axis=1, ddof=1)
    log2fc = mean_a - mean_b
    se2 = var_a / np.maximum(n_a, 1) + var_b / np.maximum(n_b, 1)
    tvals = log2fc / np.sqrt(se2)
    denom = (var_a / np.maximum(n_a, 1)) ** 2 / np.maximum(n_a - 1, 1) + (
        var_b / np.maximum(n_b, 1)
    ) ** 2 / np.maximum(n_b - 1, 1)
    dfree = se2**2 / denom
    pvals = 2 * stats.t.sf(np.abs(tvals), dfree)
    bad = (n_a < 2) | (n_b < 2) | ~np.isfinite(se2) | (se2 == 0)
    tvals[bad] = np.nan
    pvals[bad] = np.nan
    out = pd.DataFrame(
        {
            "gene": logx.index,
            "n_treat": n_a,
            "n_ctrl": n_b,
            "mean_treat_log2": mean_a,
            "mean_ctrl_log2": mean_b,
            "log2FC": log2fc,
            "t": tvals,
            "p": pvals,
        }
    )
    mask = out["p"].notna()
    q = np.full(len(out), np.nan)
    if mask.sum() > 0:
        q[mask.to_numpy()] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    out["q"] = q
    return out


def paired_rows(logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    a = logx[treat].to_numpy(dtype=float)
    b = logx[ctrl].to_numpy(dtype=float)
    d = a - b
    n = np.isfinite(d).sum(axis=1)
    mean_d = np.nanmean(d, axis=1)
    sd = np.nanstd(d, axis=1, ddof=1)
    tvals = mean_d / (sd / np.sqrt(np.maximum(n, 1)))
    pvals = 2 * stats.t.sf(np.abs(tvals), n - 1)
    bad = (n < 3) | ~np.isfinite(sd) | (sd == 0)
    tvals[bad] = np.nan
    pvals[bad] = np.nan
    out = pd.DataFrame(
        {
            "gene": logx.index,
            "n_pairs": n,
            "mean_treat_log2": np.nanmean(a, axis=1),
            "mean_ctrl_log2": np.nanmean(b, axis=1),
            "log2FC": mean_d,
            "t": tvals,
            "p": pvals,
        }
    )
    mask = out["p"].notna()
    q = np.full(len(out), np.nan)
    if mask.sum() > 0:
        q[mask.to_numpy()] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    out["q"] = q
    return out


def geneset_mw(de: pd.DataFrame, genes: list[str]) -> dict:
    present = [g for g in genes if g in set(de["gene"])]
    set_fc = de.loc[de["gene"].isin(present), "log2FC"].dropna()
    bg = de.loc[~de["gene"].isin(present), "log2FC"].dropna()
    rec = {
        "n_in_set": len(genes),
        "n_present": len(present),
        "n_tested": int(len(set_fc)),
        "median_log2FC_set": float(set_fc.median()) if len(set_fc) else None,
        "median_log2FC_bg": float(bg.median()) if len(bg) else None,
        "n_up": int((set_fc > 0).sum()) if len(set_fc) else 0,
        "n_down": int((set_fc < 0).sum()) if len(set_fc) else 0,
        "mw_p": None,
    }
    if len(set_fc) >= 3 and len(bg) >= 20:
        rec["mw_p"] = float(stats.mannwhitneyu(set_fc, bg, alternative="two-sided").pvalue)
    return rec


def load_fpkm_symbol(path: Path, encoding: str) -> pd.DataFrame:
    raw = path.read_bytes()
    import gzip

    text = gzip.decompress(raw).decode(encoding)
    df = pd.read_csv(StringIO(text), sep="\t")
    if "gene_name" not in df.columns:
        raise ValueError(f"no gene_name in {path}")
    value_cols = [c for c in df.columns if c not in {
        "gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
        "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family",
    }]
    # keep highest-mean symbol if duplicated
    mat = df.set_index("gene_name")[value_cols].apply(pd.to_numeric, errors="coerce")
    mat = mat.groupby(mat.index).mean()
    return mat


def score_contrast(name: str, de: pd.DataFrame, expr: pd.DataFrame, treat, ctrl) -> dict:
    key = de[de["gene"].isin(KEY_GENES)].copy()
    sets = {}
    for set_name, genes in all_sets().items():
        rec = geneset_mw(de, genes)
        rec["set"] = set_name
        rec["series"] = name
        sets[set_name] = rec
    out = {
        "series": name,
        "n_genes_tested": int(de["p"].notna().sum()),
        "key_genes": key.to_dict(orient="records"),
        "sets": sets,
    }
    # raw FPKM / count means for CLDN4 / TACSTD2
    for g in TARGETS:
        if g in expr.index:
            out[f"{g}_mean_treat"] = float(np.nanmean(expr.loc[g, treat].to_numpy(dtype=float)))
            out[f"{g}_mean_ctrl"] = float(np.nanmean(expr.loc[g, ctrl].to_numpy(dtype=float)))
    return out


def fig_gene_bars(title: str, de: pd.DataFrame, genes: list[str], out: Path, subtitle: str):
    rows = []
    for g in genes:
        hit = de[de["gene"] == g]
        if hit.empty:
            continue
        r = hit.iloc[0]
        rows.append((g, float(r["log2FC"]), r["p"]))
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    names = [r[0] for r in rows]
    fcs = [r[1] for r in rows]
    colors = ["#b2182b" if v > 0 else "#2166ac" for v in fcs]
    ax.bar(range(len(names)), fcs, color=colors, edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("log2FC (ADC − vehicle)")
    ax.set_title(title)
    ax.text(0.01, 0.98, subtitle, transform=ax.transAxes, va="top", fontsize=8)
    for i, (g, fc, p) in enumerate(rows):
        if pd.notna(p) and p < 0.05:
            ax.text(i, fc + (0.04 if fc >= 0 else -0.04), "*", ha="center", va="bottom" if fc >= 0 else "top")
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=DPI)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_pdx_pairs(expr: pd.DataFrame, models: list[str], genes: list[str], out: Path):
    fig, axes = plt.subplots(1, len(genes), figsize=(3.1 * len(genes), 3.6), sharex=True)
    if len(genes) == 1:
        axes = [axes]
    for ax, g in zip(axes, genes):
        if g not in expr.index:
            ax.set_title(f"{g} missing")
            continue
        xs = [0, 1]
        for m in models:
            c = float(expr.loc[g, f"C_{m}"])
            t = float(expr.loc[g, f"T_{m}"])
            ax.plot(xs, [np.log2(c + PSEUDO), np.log2(t + PSEUDO)], "-o", label=f"PDX-{m}", ms=4)
        ax.set_xticks(xs)
        ax.set_xticklabels(["vehicle", "IMMU132"])
        ax.set_title(g)
        ax.set_ylabel("log2(FPKM+1)")
    axes[-1].legend(fontsize=7, loc="best")
    fig.suptitle("GSE311016 CRC PDX, day 29, paired models (not SuperSeries)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=DPI)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def fig_forest(set_rows: list[dict], out: Path):
    df = pd.DataFrame(set_rows)
    df = df.dropna(subset=["median_log2FC_set"])
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(8.5, 0.42 * len(df) + 1.6))
    y = np.arange(len(df))
    colors = ["#b2182b" if v > 0 else "#2166ac" for v in df["median_log2FC_set"]]
    ax.scatter(df["median_log2FC_set"], y, c=colors, s=40, zorder=3)
    ax.axvline(0, color="black", lw=0.8)
    labels = [f"{r.series} · {r.set} (n={int(r.n_tested)})" for r in df.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("median gene-level log2FC in set vs vehicle")
    ax.set_title("Extra TROP2-ADC series: gene-set median log2FC (not NES; FPKM/pseudobulk)")
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=DPI)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def score_gse311016() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    mat = load_fpkm_symbol(RAW / "GSE311016_gene_fpkm.txt.gz", "utf-16")
    models = ["36", "82", "83", "114", "196"]
    treat = [f"T_{m}" for m in models]
    ctrl = [f"C_{m}" for m in models]
    logx = np.log2(mat + PSEUDO)
    de = paired_rows(logx[treat + ctrl], treat, ctrl)
    rec = score_contrast("GSE311016", de, mat, treat, ctrl)
    rec.update({
        "tissue": "CRC PDX",
        "agent": "IMMU132 (sacituzumab govitecan)",
        "contrast": "IMMU132 vs vehicle, day 29, paired by PDX model",
        "n": "5 pairs (one tumor per model per arm)",
        "statistic": "paired t on log2(FPKM+1)",
        "note": "Regular GEO Series, not a SuperSeries. Companion to given GSE312098 cell-line RNA from the same CRC paper. No organoid RNA deposited.",
    })
    fig_pdx_pairs(mat, models, ["CLDN4", "TACSTD2", "ISG15"], FIGS / "fig_extra1_gse311016_pdx_pairs")
    fig_gene_bars(
        "GSE311016 CRC PDX IMMU132 vs vehicle (paired, day 29)",
        de,
        KEY_GENES,
        FIGS / "fig_extra2_gse311016_keygenes",
        "n=5 pairs. CRC. Not SKB264. Not GSE312098.",
    )
    de.to_csv(TABLES / "de_GSE311016_IMMU132_vs_vehicle.tsv.gz", sep="\t", index=False)
    return rec, de, mat


def score_gse304294() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    mat = load_fpkm_symbol(RAW / "GSE304294_gene_fpkm.txt.gz", "utf-8")
    # Column order matches series-matrix sample order: Control×3, IMMU132×2, IACS×3, Combo×3.
    ctrl = ["OX1_1", "OX1_2", "OX1_3"]
    treat = ["OX2_1", "OX2_2"]
    logx = np.log2(mat + PSEUDO)
    de = welch_rows(logx[treat + ctrl], treat, ctrl)
    rec = score_contrast("GSE304294", de, mat, treat, ctrl)
    rec.update({
        "tissue": "ESCC cell line KYSE30",
        "agent": "IMMU132 (sacituzumab govitecan)",
        "contrast": "IMMU132 vs vehicle, 1 day",
        "n": "2 IMMU132 vs 3 vehicle (underpowered)",
        "statistic": "Welch t on log2(FPKM+1)",
        "note": "Real GEO Series with processed FPKM (4.4 Mb). OX1=Control, OX2=IMMU132 mapped by 3/2/3/3 column counts matching series-matrix sample order. Combination/IACS arms not used.",
    })
    fig_gene_bars(
        "GSE304294 KYSE30 ESCC IMMU132 vs vehicle (1 day)",
        de,
        KEY_GENES,
        FIGS / "fig_extra3_gse304294_keygenes",
        "n=2 vs 3. ESCC, not lung, not SKB264. Underpowered IMMU arm.",
    )
    # strip plot for CLDN4 / TACSTD2
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.6))
    for ax, g in zip(axes, ["CLDN4", "TACSTD2"]):
        cv = np.log2(mat.loc[g, ctrl].to_numpy(dtype=float) + PSEUDO)
        tv = np.log2(mat.loc[g, treat].to_numpy(dtype=float) + PSEUDO)
        ax.scatter(np.zeros(len(cv)), cv, c="#4d4d4d", s=36, label="vehicle")
        ax.scatter(np.ones(len(tv)), tv, c="#b2182b", s=36, label="IMMU132")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["vehicle\nn=3", "IMMU132\nn=2"])
        ax.set_title(g)
        ax.set_ylabel("log2(FPKM+1)")
    axes[1].legend(fontsize=8)
    fig.suptitle("GSE304294 KYSE30, 1 day (ESCC)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra3b_gse304294_points.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra3b_gse304294_points.pdf")
    plt.close(fig)
    de.to_csv(TABLES / "de_GSE304294_IMMU132_vs_vehicle.tsv.gz", sep="\t", index=False)
    return rec, de, mat


def score_emtab16433() -> tuple[dict | None, pd.DataFrame | None]:
    counts_path = Path("/tmp/trop2_extra/emtab16433/treatment_HD4246_SG_counts.txt")
    meta_path = Path("/tmp/trop2_extra/emtab16433/treatment_HD4246_SG_metadata.tsv")
    if not counts_path.exists() or counts_path.stat().st_size < 700_000_000:
        return None, None
    meta = pd.read_csv(meta_path, sep="\t")
    # Stream only needed genes + write mouse pseudobulk for all genes if memory allows.
    # Matrix is genes × cells, first column header is "barcode" (gene symbol).
    usecols = ["barcode"] + meta["barcode"].tolist()
    # pandas will load ~1.5–2 GB; acceptable.
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    # align cells
    cells = [c for c in counts.columns if c in set(meta["barcode"])]
    counts = counts[cells]
    meta = meta.set_index("barcode").loc[cells]
    # mouse-level pseudobulk (sum UMIs)
    pb = {}
    for sid, sub in meta.groupby("sample_id"):
        pb[sid] = counts[sub.index].sum(axis=1)
    pb = pd.DataFrame(pb)
    # CPM-like then log2
    lib = pb.sum(axis=0)
    cpm = pb.div(lib, axis=1) * 1e6
    logx = np.log2(cpm + PSEUDO)
    treat = [c for c in logx.columns if c.startswith("trodelvy_")]
    ctrl = [c for c in logx.columns if c.startswith("vehicle_")]
    de = welch_rows(logx, treat, ctrl)
    rec = score_contrast("E-MTAB-16433", de, cpm, treat, ctrl)
    rec.update({
        "tissue": "CRC PDOX (subcutaneous, HD42466)",
        "agent": "Trodelvy / sacituzumab govitecan",
        "contrast": "SG vs vehicle, 28 days, mouse-level pseudobulk",
        "n": f"{len(treat)} SG mice vs {len(ctrl)} vehicle mice (hashed; 1 library per arm)",
        "statistic": "Welch t on log2(CPM+1) of mouse pseudobulk",
        "n_cells_sg": int((meta["treatment"] == "trodelvy").sum()),
        "n_cells_vehicle": int((meta["treatment"] == "vehicle").sum()),
        "note": (
            "Processed counts 714 Mb (<2 GB). Different paper from GSE311016/GSE312098 "
            "(ORCA-HD TROP2_in_CRC). All vehicle cells are in library S01 and all SG cells "
            "in S02, so library is aliased with treatment. Mouse n=4 vs 4 is hashed within "
            "those two libraries."
        ),
    })
    fig_gene_bars(
        "E-MTAB-16433 CRC PDOX SG vs vehicle (mouse pseudobulk)",
        de,
        KEY_GENES,
        FIGS / "fig_extra4_emtab16433_keygenes",
        f"n={len(treat)} vs {len(ctrl)} mice. CRC PDOX. Library aliased with treatment.",
    )
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 3.6))
    for ax, g in zip(axes, ["CLDN4", "TACSTD2", "ISG15"]):
        if g not in cpm.index:
            continue
        cv = np.log2(cpm.loc[g, ctrl].to_numpy(dtype=float) + PSEUDO)
        tv = np.log2(cpm.loc[g, treat].to_numpy(dtype=float) + PSEUDO)
        ax.scatter(np.zeros(len(cv)), cv, c="#4d4d4d", s=40, label="vehicle")
        ax.scatter(np.ones(len(tv)), tv, c="#b2182b", s=40, label="Trodelvy")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["vehicle\nn=4", "SG\nn=4"])
        ax.set_title(g)
        ax.set_ylabel("log2(CPM+1)")
    axes[2].legend(fontsize=8)
    fig.suptitle("E-MTAB-16433 CRC PDOX mouse pseudobulk (library aliased with treatment)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_extra4b_emtab16433_mouse_points.png", dpi=DPI)
    fig.savefig(FIGS / "fig_extra4b_emtab16433_mouse_points.pdf")
    plt.close(fig)
    de.to_csv(TABLES / "de_EMTAB16433_SG_vs_vehicle_mouse_pseudobulk.tsv.gz", sep="\t", index=False)
    pb.to_csv(TABLES / "EMTAB16433_mouse_pseudobulk_counts.tsv.gz", sep="\t")
    return rec, de


def main():
    scored = []
    set_rows = []
    rec311, de311, _ = score_gse311016()
    scored.append(rec311)
    set_rows.extend(rec311["sets"].values())
    rec304, de304, _ = score_gse304294()
    scored.append(rec304)
    set_rows.extend(rec304["sets"].values())
    rec164, de164 = score_emtab16433()
    if rec164 is not None:
        scored.append(rec164)
        set_rows.extend(rec164["sets"].values())

    fig_forest(set_rows, FIGS / "fig_extra5_geneset_median_log2fc")

    # compact key-gene table
    rows = []
    for rec, de in [
        (rec311, de311),
        (rec304, de304),
        (rec164, de164),
    ]:
        if rec is None or de is None:
            continue
        for g in KEY_GENES:
            hit = de[de["gene"] == g]
            if hit.empty:
                continue
            r = hit.iloc[0]
            rows.append({
                "series": rec["series"],
                "tissue": rec["tissue"],
                "n": rec["n"],
                "gene": g,
                "log2FC": float(r["log2FC"]),
                "p": None if pd.isna(r["p"]) else float(r["p"]),
                "q": None if pd.isna(r["q"]) else float(r["q"]),
            })
    key_tbl = pd.DataFrame(rows)
    key_tbl.to_csv(TABLES / "key_genes_extra_series.tsv", sep="\t", index=False)
    set_tbl = pd.DataFrame(set_rows)
    set_tbl.to_csv(TABLES / "geneset_stats_extra_series.tsv", sep="\t", index=False)

    # JSON without huge key_genes duplication issues
    slim = []
    for rec in scored:
        slim.append({k: v for k, v in rec.items() if k != "key_genes"})
        slim[-1]["key_CLDN4"] = next((x for x in rec["key_genes"] if x["gene"] == "CLDN4"), None)
        slim[-1]["key_TACSTD2"] = next((x for x in rec["key_genes"] if x["gene"] == "TACSTD2"), None)
    (HERE / "key_stats.json").write_text(json.dumps(slim, indent=2, default=str))
    print(json.dumps(slim, indent=2, default=str))


if __name__ == "__main__":
    main()

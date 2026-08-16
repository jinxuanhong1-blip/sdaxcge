#!/usr/bin/env python3
"""Public human lung-tumor scRNA with pre AND post ICI / chemo-IO.

Primary testable cohort: GSE207422 (3 unmatched pre-treatment biopsies +
12 unmatched post-treatment resections). Malignant TACSTD2 and CLDN4 are
compared at the sample unit (Mann-Whitney U).

GSE337519 deposits one unlabeled 10x library despite a paired series title;
pre vs post cannot be scored.

Zhejiang A7 paired IHC is not re-analyzed here.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread

LINEAGE_MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2"],
    "Mast": ["TPSAB1", "CPA3"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
ALVEOLAR = ["SFTPA2", "SFTPA1", "SFTPC", "SFTPB", "AGER"]
CLUB = ["SCGB1A1", "SCGB3A2"]
CILIATED = ["TPPP3", "FOXJ1", "CAPS"]
PANEL = ["TACSTD2", "CLDN4", "PTPRC", "EPCAM"]
MIN_MAL = 10
GENES = ["TACSTD2", "CLDN4"]


def wanted_genes() -> set[str]:
    genes: set[str] = set(PANEL)
    for block in (LINEAGE_MARKERS, {"a": ALVEOLAR, "c": CLUB, "i": CILIATED}):
        for vs in block.values():
            genes.update(vs)
    return genes


def stream_gse207422(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        n_umi = np.zeros(n_cells, dtype=np.float64)
        n_genes = np.zeros(n_cells, dtype=np.int32)
        for i, line in enumerate(handle, start=1):
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            n_umi += arr
            n_genes += arr > 0
            if gene in wanted:
                found[gene] = arr
            if i % 5000 == 0:
                print(f"  streamed {i} genes, stored {len(found)}", flush=True)
    return cell_ids, found, n_umi, n_genes


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best]
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < 0.12] = "Unassigned"
    return labels


def gene_log(log_cp: dict[str, np.ndarray], name: str, n: int) -> np.ndarray:
    return log_cp[name] if name in log_cp else np.zeros(n, dtype=np.float32)


def mw(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    rec = {
        "n_pre": int(len(a)),
        "n_post": int(len(b)),
        "median_pre": float(np.median(a)) if len(a) else None,
        "median_post": float(np.median(b)) if len(b) else None,
        "delta_post_minus_pre": None,
        "U": None,
        "p": None,
        "direction": "NA",
        "test": "Mann-Whitney U (unpaired)",
    }
    if len(a) and len(b):
        rec["delta_post_minus_pre"] = rec["median_post"] - rec["median_pre"]
        rec["direction"] = "up" if rec["delta_post_minus_pre"] > 0 else ("down" if rec["delta_post_minus_pre"] < 0 else "flat")
    if len(a) < 2 or len(b) < 2:
        rec["note"] = "need_n>=2_per_arm"
        return rec
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rec["U"] = float(u)
    rec["p"] = float(p)
    rec["note"] = ""
    return rec


def p_to_z(p: float, direction: str) -> float:
    """Two-sided p → signed Stouffer z (up = positive)."""
    p = min(max(p, 1e-300), 1.0 - 1e-16)
    z = float(stats.norm.isf(p / 2.0))
    if direction == "down":
        z = -z
    elif direction == "flat":
        z = 0.0
    return z


def stouffer(rows: list[dict]) -> dict:
    usable = [r for r in rows if r.get("p") is not None and r.get("direction") in {"up", "down", "flat"}]
    n_up = sum(1 for r in usable if r["direction"] == "up")
    n_down = sum(1 for r in usable if r["direction"] == "down")
    n_flat = sum(1 for r in usable if r["direction"] == "flat")
    if not usable:
        return {
            "n_studies": 0,
            "n_up": 0,
            "n_down": 0,
            "n_flat": 0,
            "stouffer_z": None,
            "stouffer_p": None,
            "note": "no_testable_cohorts",
        }
    zs = np.array([p_to_z(r["p"], r["direction"]) for r in usable], dtype=float)
    z = float(zs.sum() / np.sqrt(len(zs)))
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "n_studies": int(len(usable)),
        "n_up": int(n_up),
        "n_down": int(n_down),
        "n_flat": int(n_flat),
        "stouffer_z": z,
        "stouffer_p": p,
        "note": "one_independent_tumor_cell_cohort" if len(usable) == 1 else "",
    }


def sample_gene_stats(g: pd.DataFrame, gene: str) -> dict:
    if g.empty or g["n_umi"].sum() <= 0:
        return {"n": 0, "pct_pos": np.nan, "mean_log1p_cp10k": np.nan, "pb_cpm": np.nan}
    umi = g[gene].to_numpy(dtype=float)
    scale = 1e4 / g["n_umi"].to_numpy(dtype=float)
    return {
        "n": int(len(g)),
        "pct_pos": float(100.0 * (umi > 0).mean()),
        "mean_log1p_cp10k": float(np.log1p(umi * scale).mean()),
        "pb_cpm": float(umi.sum() / g["n_umi"].sum() * 1e6),
    }


def analyze_gse207422(matrix: Path, metadata: Path, outdir: Path) -> list[dict]:
    meta = pd.read_excel(metadata)
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["response_paper"] = meta["Pathologic Response"].replace({"pCR": "MPR"})
    meta["timepoint"] = np.where(meta["Resource"].astype(str).str.contains("Pre", case=False), "pre", "post")
    meta.to_csv(outdir / "gse207422_sample_metadata.tsv", sep="\t", index=False)

    wanted = wanted_genes()
    print("GSE207422: streaming UMI matrix", flush=True)
    cell_ids, expr, n_umi, n_genes = stream_gse207422(matrix, wanted)
    n = len(cell_ids)
    print(f"  cells={n} stored={len(expr)} median_nUMI={np.median(n_umi):.0f} median_nGene={np.median(n_genes):.0f}", flush=True)

    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    lineage = assign_lineage(scores, gene_log(log_cp, "CD3D", n))
    epi = lineage == "Epithelial"
    alveolar = np.maximum.reduce([gene_log(log_cp, g, n) for g in ALVEOLAR])
    club = np.maximum.reduce([gene_log(log_cp, g, n) for g in CLUB])
    ciliated = np.maximum.reduce([gene_log(log_cp, g, n) for g in CILIATED])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    is_malignant = epi & (~clear_normal)

    sample = np.array([c.rsplit("_", 1)[0] for c in cell_ids])
    cells = pd.DataFrame(
        {
            "barcode": cell_ids,
            "Sample": sample,
            "n_umi": n_umi,
            "n_genes": n_genes,
            "lineage": lineage,
            "is_malignant": is_malignant,
            "is_epithelial": epi,
            "TACSTD2": expr.get("TACSTD2", np.zeros(n)),
            "CLDN4": expr.get("CLDN4", np.zeros(n)),
            "EPCAM": expr.get("EPCAM", np.zeros(n)),
            "PTPRC": expr.get("PTPRC", np.zeros(n)),
        }
    )
    cells = cells.merge(meta, on="Sample", how="left")
    cells.to_csv(outdir / "gse207422_per_cell_labels.tsv.gz", sep="\t", index=False)

    rows = []
    for sample_id, g in cells.groupby("Sample"):
        rec = {
            "accession": "GSE207422",
            "Sample": sample_id,
            "Patient": g["Patient"].iloc[0],
            "timepoint": g["timepoint"].iloc[0],
            "Resource": g["Resource"].iloc[0],
            "Pathology": g["Pathology"].iloc[0],
            "response_paper": g["response_paper"].iloc[0],
            "n_cells": int(len(g)),
            "n_epithelial": int(g["is_epithelial"].sum()),
            "n_malignant": int(g["is_malignant"].sum()),
        }
        for gene in GENES:
            mal = sample_gene_stats(g[g["is_malignant"]], gene)
            epi_s = sample_gene_stats(g[g["is_epithelial"]], gene)
            rec[f"mal_{gene}_n"] = mal["n"]
            rec[f"mal_{gene}_pct_pos"] = mal["pct_pos"]
            rec[f"mal_{gene}_mean_log1p"] = mal["mean_log1p_cp10k"]
            rec[f"mal_{gene}_pb_cpm"] = mal["pb_cpm"]
            rec[f"epi_{gene}_mean_log1p"] = epi_s["mean_log1p_cp10k"]
            rec[f"epi_{gene}_pct_pos"] = epi_s["pct_pos"]
            rec[f"epi_{gene}_pb_cpm"] = epi_s["pb_cpm"]
        rows.append(rec)
    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(outdir / "gse207422_per_sample.tsv", sep="\t", index=False)

    tests = []
    mal_ok = sample_df[sample_df["n_malignant"] >= MIN_MAL].copy()
    epi_ok = sample_df[sample_df["n_epithelial"] >= MIN_MAL].copy()
    contrasts = [
        (
            "epi_all_post_vs_pre",
            epi_ok,
            "epi",
            "unmatched 3 pre vs 12 post; epithelial compartment (CopyKAT not deposited)",
        ),
        (
            "epi_nmpr_post_vs_pre",
            epi_ok[epi_ok["timepoint"].eq("pre") | epi_ok["response_paper"].eq("NMPR")],
            "epi",
            "unmatched pre vs NMPR post; epithelial compartment",
        ),
        (
            "mal_all_post_vs_pre",
            mal_ok,
            "mal",
            "unmatched pre vs post with ≥10 malignant-like cells (epithelial minus alveolar/club/ciliated)",
        ),
        (
            "mal_nmpr_post_vs_pre",
            mal_ok[mal_ok["timepoint"].eq("pre") | mal_ok["response_paper"].eq("NMPR")],
            "mal",
            "unmatched pre vs NMPR post; ≥10 malignant-like cells",
        ),
    ]
    for contrast, sub, prefix, note in contrasts:
        pre = sub[sub["timepoint"] == "pre"]
        post = sub[sub["timepoint"] == "post"]
        for gene in GENES:
            for metric, suffix in [
                ("mean_log1p_cp10k", "mean_log1p"),
                ("pct_pos", "pct_pos"),
                ("pb_cpm", "pb_cpm"),
            ]:
                col = f"{prefix}_{gene}_{suffix}" if suffix != "mean_log1p" else f"{prefix}_{gene}_mean_log1p"
                if suffix == "pct_pos":
                    col = f"{prefix}_{gene}_pct_pos"
                elif suffix == "pb_cpm":
                    col = f"{prefix}_{gene}_pb_cpm"
                if col not in sub.columns:
                    continue
                stat = mw(pre[col], post[col])
                stat.update(
                    {
                        "accession": "GSE207422",
                        "contrast": contrast,
                        "gene": gene,
                        "metric": metric,
                        "compartment": "epithelial" if prefix == "epi" else "malignant-like",
                        "pairing": "unmatched",
                        "n_pre_patients": int(pre["Patient"].nunique()),
                        "n_post_patients": int(post["Patient"].nunique()),
                        "pre_patients": ",".join(pre["Patient"].astype(str)),
                        "post_patients": ",".join(post["Patient"].astype(str)),
                        "design_note": note,
                    }
                )
                tests.append(stat)

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(outdir / "gse207422_pre_vs_post.tsv", sep="\t", index=False)

    # figure: primary epithelial mean log1p (all 15 libraries)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=False)
    plot_df = epi_ok.copy()
    for ax, gene in zip(axes, GENES):
        col = f"epi_{gene}_mean_log1p"
        pre_v = plot_df.loc[plot_df["timepoint"] == "pre", col]
        post_v = plot_df.loc[plot_df["timepoint"] == "post", col]
        ax.boxplot([pre_v, post_v], tick_labels=["pre", "post"], widths=0.45)
        rng = np.random.default_rng(1)
        ax.scatter(1 + rng.uniform(-0.08, 0.08, len(pre_v)), pre_v, c="#1f77b4", s=36, zorder=3)
        colors = plot_df.loc[plot_df["timepoint"] == "post", "response_paper"].map(
            {"MPR": "#2ca02c", "NMPR": "#d62728"}
        )
        ax.scatter(2 + rng.uniform(-0.08, 0.08, len(post_v)), post_v, c=colors, s=36, zorder=3)
        stat = mw(pre_v, post_v)
        ax.set_title(f"{gene}\nMWU p={stat['p']:.3g}" if stat["p"] is not None else gene)
        ax.set_ylabel("epithelial mean log1p(CP10K)")
    fig.suptitle("GSE207422 unmatched pre vs post (sample unit, epithelial, all 15 libraries)", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "fig_gse207422_malignant_pre_post.png", dpi=160)
    plt.close(fig)
    return tests


def analyze_gse337519(raw_dir: Path, outdir: Path) -> dict:
    """One unlabeled 10x library. Cannot assign pre vs post."""
    barcodes = gzip.open(raw_dir / "GSM9856929_barcodes.tsv.gz", "rt").read().splitlines()
    features = [line.split("\t") for line in gzip.open(raw_dir / "GSM9856929_features.tsv.gz", "rt")]
    genes = [f[1] if len(f) > 1 else f[0] for f in features]
    mat = mmread(raw_dir / "GSM9856929_matrix.mtx.gz").tocsr()
    gene_idx = {g: i for i, g in enumerate(genes)}
    n_umi = np.asarray(mat.sum(axis=0)).ravel()
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)

    def vec(name: str) -> np.ndarray:
        if name not in gene_idx:
            return np.zeros(mat.shape[1], dtype=np.float32)
        return np.asarray(mat[gene_idx[name], :].todense()).ravel().astype(np.float32)

    log_cp = {}
    for g in wanted_genes():
        if g in gene_idx:
            log_cp[g] = np.log1p(vec(g) * scale).astype(np.float32)
    n = mat.shape[1]
    scores = {name: module_score(log_cp, gs, n) for name, gs in LINEAGE_MARKERS.items()}
    lineage = assign_lineage(scores, gene_log(log_cp, "CD3D", n))
    epi = lineage == "Epithelial"
    alveolar = np.maximum.reduce([gene_log(log_cp, g, n) for g in ALVEOLAR])
    club = np.maximum.reduce([gene_log(log_cp, g, n) for g in CLUB])
    ciliated = np.maximum.reduce([gene_log(log_cp, g, n) for g in CILIATED])
    clear_normal = (alveolar >= 1.0) | (club >= 1.0) | (ciliated >= 1.0)
    is_mal = epi & (~clear_normal)
    rec = {
        "accession": "GSE337519",
        "n_cells": int(n),
        "n_epithelial": int(epi.sum()),
        "n_malignant": int(is_mal.sum()),
        "pairing_deposited": "series text claims 1 paired patient; GEO deposits 1 unlabeled library (GSM9856929, 'one sample')",
        "pre_vs_post_testable": False,
    }
    for gene in GENES:
        umi = vec(gene)
        mal = is_mal
        rec[f"mal_{gene}_n"] = int(mal.sum())
        rec[f"mal_{gene}_pct_pos"] = float(100.0 * (umi[mal] > 0).mean()) if mal.any() else None
        rec[f"mal_{gene}_mean_log1p"] = float(np.log1p(umi[mal] * scale[mal]).mean()) if mal.any() else None
    pd.DataFrame([rec]).to_csv(outdir / "gse337519_single_library.tsv", sep="\t", index=False)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse207422-matrix", type=Path, required=True)
    ap.add_argument("--gse207422-metadata", type=Path, required=True)
    ap.add_argument("--gse337519-raw", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tests = analyze_gse207422(args.gse207422_matrix, args.gse207422_metadata, args.outdir)
    g337 = analyze_gse337519(args.gse337519_raw, args.outdir)

    # Combined estimate: one independent tumor-cell cohort (GSE207422 primary).
    paper_rows = []
    for gene in GENES:
        primary = next(
            r
            for r in tests
            if r["gene"] == gene and r["contrast"] == "epi_all_post_vs_pre" and r["metric"] == "mean_log1p_cp10k"
        )
        paper_rows.append(
            {
                "accession": "GSE207422",
                "citation": "Hu et al. Genome Med 2023",
                "regimen": "neoadjuvant PD-1 + chemo",
                "pairing": "unmatched (different patients)",
                "n_pre": primary["n_pre"],
                "n_post": primary["n_post"],
                "gene": gene,
                "metric": "epithelial mean log1p(CP10K)",
                "median_pre": primary["median_pre"],
                "median_post": primary["median_post"],
                "direction": primary["direction"],
                "test": primary["test"],
                "U": primary["U"],
                "p": primary["p"],
                "note": "3 pre biopsies vs 12 post resections; sample unit; GEO has no CopyKAT labels",
            }
        )
    paper = pd.DataFrame(paper_rows)
    paper.to_csv(args.outdir / "paper_table_scrna_prepost.tsv", sep="\t", index=False)

    combined = []
    for gene in GENES:
        gene_rows = [r for r in paper_rows if r["gene"] == gene]
        s = stouffer(gene_rows)
        s["gene"] = gene
        s["metric"] = "malignant mean log1p(CP10K)"
        s["studies"] = "GSE207422"
        combined.append(s)
    pd.DataFrame(combined).to_csv(args.outdir / "combined_stouffer.tsv", sep="\t", index=False)

    summary = {
        "testable_tumor_cell_cohorts": ["GSE207422"],
        "claimed_but_not_testable": {
            "GSE337519": g337["pairing_deposited"],
        },
        "min_malignant_cells": MIN_MAL,
        "combined": combined,
        "gse337519": g337,
    }
    (args.outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()

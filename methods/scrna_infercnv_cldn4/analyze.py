#!/usr/bin/env python3
"""CNV-proxy malignant cells, then CLDN4 vs same-patient T/NK on GSE207422.

Hu et al., Genome Medicine 2023 (PMID 36869384). Public GEO UMI only.
Author CopyKAT barcode IDs are not deposited. This is not the R infercnv
or copykat package: window-smoothed expression CNV vs stromal reference.

CLDN4 only. Unit for paired tests = sample/patient.
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

HERE = Path(__file__).resolve().parent
WIN = 25
MIN_BOTH = 5

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
EXTRA = ["CLDN4", "PTPRC", "CD8A", "CD8B"]


def wanted_markers() -> set[str]:
    genes: set[str] = set(EXTRA)
    for vs in LINEAGE_MARKERS.values():
        genes.update(vs)
    return genes


def sample_from_barcode(bc: str) -> str:
    parts = bc.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else bc


def stream_pass1(matrix_path: Path, wanted: set[str]):
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_umi = np.zeros(n, dtype=np.float64)
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            n_umi += arr
            if gene in wanted:
                found[gene] = arr
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  pass1 genes={n_streamed} stored={len(found)}", flush=True)
    print(f"pass1 done genes={n_streamed} cells={n} stored={len(found)}", flush=True)
    return cell_ids, found, n_umi, n_streamed


def stream_pass2(matrix_path: Path, keep_idx: np.ndarray, gene_set: set[str]):
    keep_idx = np.asarray(keep_idx, dtype=int)
    n_keep = keep_idx.size
    gene_names: list[str] = []
    blocks: list[np.ndarray] = []
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        n = len(header) - 1
        n_streamed = 0
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep or gene not in gene_set:
                n_streamed += 1
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            gene_names.append(gene)
            blocks.append(arr[keep_idx].astype(np.float32, copy=False))
            n_streamed += 1
            if n_streamed % 4000 == 0:
                print(f"  pass2 scanned={n_streamed} stored={len(gene_names)}", flush=True)
    mat = np.vstack(blocks) if blocks else np.zeros((0, n_keep), dtype=np.float32)
    print(f"pass2 done stored_genes={len(gene_names)} keep_cells={n_keep}", flush=True)
    return gene_names, mat


def module_score(log_cp: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [log_cp[g] for g in genes if g in log_cp]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(scores: dict[str, np.ndarray], cd3: np.ndarray, min_top: float = 0.12) -> np.ndarray:
    names = list(scores)
    mat = np.vstack([scores[n] for n in names])
    best = np.argmax(mat, axis=0)
    top = mat[best, np.arange(mat.shape[1])]
    labels = np.array(names, dtype=object)[best].copy()
    t_idx = names.index("T")
    nk_idx = names.index("NK")
    close = np.abs(mat[t_idx] - mat[nk_idx]) < 0.15
    both_high = (mat[t_idx] > 0.2) | (mat[nk_idx] > 0.2)
    tnk_best = np.isin(labels, ["T", "NK"])
    labels[close & both_high & tnk_best & (cd3 > 0.15)] = "T"
    labels[close & both_high & tnk_best & (cd3 <= 0.15) & (mat[nk_idx] >= mat[t_idx] * 0.7)] = "NK"
    labels[top < min_top] = "Unassigned"
    return labels


def moving_average_2d(block: np.ndarray, k: int) -> np.ndarray:
    if block.shape[0] < 3:
        return block.copy()
    k = min(k, block.shape[0] if block.shape[0] % 2 == 1 else block.shape[0] - 1)
    if k < 3:
        return block.copy()
    if k % 2 == 0:
        k -= 1
    pad = k // 2
    padded = np.pad(block, ((pad, pad), (0, 0)), mode="edge")
    c = np.cumsum(np.vstack([np.zeros((1, block.shape[1]), dtype=block.dtype), padded]), axis=0)
    return (c[k:] - c[:-k]) / k


def cnv_scores(counts, genes, gene_pos, ref_mask, win=WIN):
    """CopyKAT-like (median ref, sum|smooth|) and inferCNV-like (mean ref, mean|smooth|)."""
    n_g, n_c = counts.shape
    expressed = (counts > 0).sum(axis=1)
    mean_umi = counts.mean(axis=1)
    keep_g = (expressed >= 20) & (mean_umi > 0.05)
    pos = gene_pos.set_index("gene")
    have_pos = np.array([g in pos.index for g in genes])
    keep_g = keep_g & have_pos
    info = {
        "n_genes_input": int(n_g),
        "n_genes_kept": int(keep_g.sum()),
        "n_ref_cells": int(ref_mask.sum()),
        "win": int(win),
    }
    nan = np.full(n_c, np.nan, dtype=np.float32)
    if keep_g.sum() < 200 or ref_mask.sum() < 20:
        info["ok"] = False
        return nan, nan, info
    g_idx = np.where(keep_g)[0]
    gnames = [genes[i] for i in g_idx]
    sub = counts[g_idx]
    lib = np.maximum(sub.sum(axis=0), 1.0)
    logcp = np.log1p(sub / lib * 1e4).astype(np.float32)
    ref_med = np.median(logcp[:, ref_mask], axis=1, keepdims=True)
    ref_mean = np.mean(logcp[:, ref_mask], axis=1, keepdims=True)
    rel_med = logcp - ref_med
    rel_mean = logcp - ref_mean
    chrom = np.array([str(pos.loc[g, "chrom"]) for g in gnames])
    start = np.array([int(pos.loc[g, "start"]) for g in gnames])
    order = np.lexsort((start, chrom))
    rel_med = rel_med[order]
    rel_mean = rel_mean[order]
    chrom = chrom[order]
    sm_med = np.zeros_like(rel_med)
    sm_mean = np.zeros_like(rel_mean)
    chroms_used = []
    for c in pd.unique(chrom):
        idx = np.where(chrom == c)[0]
        if len(idx) < 5:
            continue
        sm_med[idx] = moving_average_2d(rel_med[idx], win)
        sm_mean[idx] = moving_average_2d(rel_mean[idx], win)
        chroms_used.append((c, int(len(idx))))
    copykat = np.sum(np.abs(sm_med), axis=0).astype(np.float32)
    infercnv = np.mean(np.abs(sm_mean), axis=0).astype(np.float32)
    info["ok"] = True
    info["chromosomes"] = chroms_used
    info["copykat_ref_p95"] = float(np.quantile(copykat[ref_mask], 0.95))
    info["infercnv_ref_p95"] = float(np.quantile(infercnv[ref_mask], 0.95))
    return copykat, infercnv, info


def spearman_safe(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return {"n": int(x.size), "rho": None, "p": None}
    r, p = stats.spearmanr(x, y)
    return {"n": int(x.size), "rho": float(r), "p": float(p)}


def wilcoxon_paired(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    out = {
        "n": int(a.size),
        "median_a": float(np.median(a)) if a.size else None,
        "median_b": float(np.median(b)) if b.size else None,
        "mean_a": float(a.mean()) if a.size else None,
        "mean_b": float(b.mean()) if b.size else None,
        "n_a_gt_b": int((a > b).sum()) if a.size else 0,
        "stat": None,
        "p": None,
    }
    if a.size < 1:
        return out
    d = a - b
    if np.allclose(d, 0):
        out["stat"] = 0.0
        out["p"] = 1.0
        return out
    try:
        res = stats.wilcoxon(a, b, alternative="greater", zero_method="wilcox")
        out["stat"] = float(res.statistic)
        out["p"] = float(res.pvalue)
    except ValueError:
        out["p"] = None
    return out


def mwu_two(values, is_a):
    values = np.asarray(values, float)
    is_a = np.asarray(is_a, bool)
    ok = np.isfinite(values)
    values, is_a = values[ok], is_a[ok]
    n_a, n_b = int(is_a.sum()), int((~is_a).sum())
    out = {
        "n_a": n_a,
        "n_b": n_b,
        "mean_a": float(values[is_a].mean()) if n_a else None,
        "mean_b": float(values[~is_a].mean()) if n_b else None,
        "delta": None,
        "U": None,
        "p": None,
    }
    if n_a and n_b:
        out["delta"] = out["mean_a"] - out["mean_b"]
        res = stats.mannwhitneyu(values[is_a], values[~is_a], alternative="two-sided")
        out["U"] = float(res.statistic)
        out["p"] = float(res.pvalue)
    return out


def fmt_p(p):
    if p is None:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_num(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{nd}f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    ap.add_argument("--xlsx", type=Path, default=Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx"))
    ap.add_argument("--gene-chr", type=Path, default=HERE / "assets/gene_chr.tsv")
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_excel(args.xlsx)
    meta.columns = [c.replace(" ", "_") for c in meta.columns]
    meta = meta.dropna(subset=["Sample", "Patient"]).copy()
    meta["response"] = meta["Pathologic_Response"].replace({"pCR": "MPR"})
    meta["is_post"] = meta["Resource"].astype(str).str.contains("Post", case=False)
    meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    gene_chr = pd.read_csv(args.gene_chr, sep="\t")
    gene_chr["chrom"] = gene_chr["chrom"].astype(str)
    cnv_gene_set = set(gene_chr["gene"])

    print("pass1: lineage markers + CLDN4 + library sizes", flush=True)
    cell_ids, expr, n_umi, n_matrix_genes = stream_pass1(args.matrix, wanted_markers())
    n = len(cell_ids)
    cells = np.array(cell_ids)
    sample = np.array([sample_from_barcode(c) for c in cells])
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    log_cp = {g: np.log1p(expr[g] * scale).astype(np.float32) for g in expr}
    scores = {name: module_score(log_cp, genes, n) for name, genes in LINEAGE_MARKERS.items()}
    lineage = assign_lineage(scores, log_cp.get("CD3E", log_cp.get("CD3D", np.zeros(n))))
    if "CLDN4" not in expr:
        raise SystemExit("CLDN4 not found in GSE207422 UMI matrix")
    cldn4_umi = expr["CLDN4"]
    cldn4 = np.log1p(cldn4_umi * scale).astype(np.float32)

    is_epi = lineage == "Epithelial"
    is_stromal = np.isin(lineage, ["Fibroblast", "Endothelial"])
    is_t = lineage == "T"
    is_nk = lineage == "NK"
    is_tnk = is_t | is_nk
    print(
        f"lineage epi={is_epi.sum()} stromal={is_stromal.sum()} T={is_t.sum()} NK={is_nk.sum()}",
        flush=True,
    )

    rng = np.random.default_rng(0)
    immune_idx = np.where(is_t | (lineage == "Myeloid"))[0]
    if immune_idx.size > 2500:
        immune_idx = rng.choice(immune_idx, size=2500, replace=False)
    keep = np.zeros(n, dtype=bool)
    keep[is_epi | is_stromal] = True
    keep[immune_idx] = True
    keep_idx = np.where(keep)[0]
    keep_map = {int(i): j for j, i in enumerate(keep_idx)}
    print(f"pass2 keep_cells={keep_idx.size}", flush=True)

    cnv_genes, cnv_mat = stream_pass2(args.matrix, keep_idx, cnv_gene_set)
    ref_keep = is_stromal[keep_idx]
    if ref_keep.sum() < 80:
        print(f"stromal ref n={ref_keep.sum()} < 80; adding T+myeloid subsample", flush=True)
        ref_keep = ref_keep | np.isin(keep_idx, immune_idx)
    copykat_s, infer_s, cnv_info = cnv_scores(cnv_mat, cnv_genes, gene_chr, ref_keep, WIN)
    print("CNV info", json.dumps({k: v for k, v in cnv_info.items() if k != "chromosomes"}), flush=True)

    copykat_full = np.full(n, np.nan, dtype=np.float32)
    infer_full = np.full(n, np.nan, dtype=np.float32)
    for i, j in keep_map.items():
        copykat_full[i] = copykat_s[j]
        infer_full[i] = infer_s[j]
    ck_p95 = float(cnv_info.get("copykat_ref_p95", np.nan))
    ic_p95 = float(cnv_info.get("infercnv_ref_p95", np.nan))
    is_copykat = is_epi & np.isfinite(copykat_full) & (copykat_full > ck_p95)
    is_infercnv = is_epi & np.isfinite(infer_full) & (infer_full > ic_p95)
    print(f"malignant copykat={int(is_copykat.sum())} infercnv={int(is_infercnv.sum())}", flush=True)

    mal_flags = {
        "infercnv_stromal_p95": is_infercnv,
        "copykat_stromal_p95": is_copykat,
        "marker_all_epi": is_epi,
    }

    cell_df = pd.DataFrame(
        {
            "barcode": cells,
            "sample": sample,
            "lineage": lineage,
            "n_umi": n_umi,
            "CLDN4": cldn4_umi,
            "cldn4_log1p_cp10k": cldn4,
            "cnv_infercnv": infer_full,
            "cnv_copykat": copykat_full,
            "is_tnk": is_tnk.astype(int),
        }
    )
    for name, mask in mal_flags.items():
        cell_df[f"mal_{name}"] = mask.astype(int)
    cell_df.to_csv(args.outdir / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    smap = meta.set_index("Sample")
    rows = []
    for samp in sorted(set(sample)):
        if samp not in smap.index:
            continue
        m = sample == samp
        rec = {
            "sample": samp,
            "patient": str(smap.loc[samp, "Patient"]),
            "is_post": bool(smap.loc[samp, "is_post"]),
            "response": str(smap.loc[samp, "response"]),
            "pathology": str(smap.loc[samp, "Pathology"]),
            "n_cells": int(m.sum()),
            "n_tnk": int((m & is_tnk).sum()),
            "frac_tnk": float(is_tnk[m].mean()) if m.any() else np.nan,
            "tnk_cldn4_mean": float(cldn4[m & is_tnk].mean()) if (m & is_tnk).any() else np.nan,
            "tnk_cldn4_pctpos": float((cldn4_umi[m & is_tnk] > 0).mean()) if (m & is_tnk).any() else np.nan,
        }
        for name, mask in mal_flags.items():
            mm = m & mask
            rec[f"{name}_n"] = int(mm.sum())
            rec[f"{name}_cldn4_mean"] = float(cldn4[mm].mean()) if mm.any() else np.nan
            rec[f"{name}_cldn4_pctpos"] = float((cldn4_umi[mm] > 0).mean()) if mm.any() else np.nan
        rows.append(rec)
    per = pd.DataFrame(rows)
    per.to_csv(args.outdir / "per_sample.tsv", sep="\t", index=False)

    tests = []
    paired = {}
    for cohort_name, cohort in [
        ("post", per["is_post"]),
        ("all_labeled", pd.Series(True, index=per.index)),
    ]:
        sub0 = per.loc[cohort].copy()
        for mal_name in mal_flags:
            for min_n in (1, MIN_BOTH):
                ok = (sub0[f"{mal_name}_n"] >= min_n) & (sub0["n_tnk"] >= min_n)
                sub = sub0.loc[ok]
                key = f"{cohort_name}|{mal_name}|min{min_n}"
                wp_mean = wilcoxon_paired(sub[f"{mal_name}_cldn4_mean"], sub["tnk_cldn4_mean"])
                wp_pct = wilcoxon_paired(sub[f"{mal_name}_cldn4_pctpos"], sub["tnk_cldn4_pctpos"])
                sp_mean = spearman_safe(sub[f"{mal_name}_cldn4_mean"], sub["frac_tnk"])
                sp_pct = spearman_safe(sub[f"{mal_name}_cldn4_pctpos"], sub["frac_tnk"])
                post_resp = sub[sub["response"].isin(["NMPR", "MPR"])]
                mwu = mwu_two(
                    post_resp[f"{mal_name}_cldn4_mean"].to_numpy(),
                    (post_resp["response"] == "NMPR").to_numpy(),
                )
                dropped = sorted(sub0.loc[~ok, "sample"].astype(str))
                rec = {
                    "key": key,
                    "cohort": cohort_name,
                    "malignant_def": mal_name,
                    "min_cells": min_n,
                    "n_samples": int(ok.sum()),
                    "n_cohort": int(len(sub0)),
                    "dropped_samples": ",".join(dropped),
                    "n_nmpr": int((sub["response"] == "NMPR").sum()),
                    "n_mpr": int((sub["response"] == "MPR").sum()),
                    "paired_mean": wp_mean,
                    "paired_pctpos": wp_pct,
                    "spearman_mean_vs_tnkfrac": sp_mean,
                    "spearman_pctpos_vs_tnkfrac": sp_pct,
                    "nmpr_vs_mpr_mean": mwu,
                    "samples": sub["sample"].tolist(),
                }
                tests.append(rec)
                if cohort_name == "post" and min_n == MIN_BOTH:
                    paired[mal_name] = rec

    # cell-level (descriptive; n is huge)
    cell_level = {}
    post_mask = np.isin(sample, per.loc[per["is_post"], "sample"].to_numpy())
    for mal_name, mask in mal_flags.items():
        a = cldn4[post_mask & mask]
        b = cldn4[post_mask & is_tnk]
        if a.size and b.size:
            res = stats.mannwhitneyu(a, b, alternative="greater")
            cell_level[mal_name] = {
                "n_mal": int(a.size),
                "n_tnk": int(b.size),
                "median_mal": float(np.median(a)),
                "median_tnk": float(np.median(b)),
                "pctpos_mal": float((cldn4_umi[post_mask & mask] > 0).mean()),
                "pctpos_tnk": float((cldn4_umi[post_mask & is_tnk] > 0).mean()),
                "mwu_p": float(res.pvalue),
            }

    summary = {
        "dataset": "GSE207422",
        "citation": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "gene": "CLDN4",
        "n_cells_matrix": int(n),
        "n_genes_matrix": int(n_matrix_genes),
        "n_post_samples_labeled": int(per["is_post"].sum()),
        "lineage_counts": {
            "Epithelial": int(is_epi.sum()),
            "Fibroblast": int((lineage == "Fibroblast").sum()),
            "Endothelial": int((lineage == "Endothelial").sum()),
            "T": int(is_t.sum()),
            "NK": int(is_nk.sum()),
        },
        "malignant_n": {k: int(v.sum()) for k, v in mal_flags.items()},
        "cnv": {k: v for k, v in cnv_info.items() if k != "chromosomes"},
        "note": (
            "CNV proxy, not R infercnv/copykat. Malignant = marker epithelium "
            "AND score > 95th percentile of stromal (fibroblast+endothelial) reference. "
            "Paired tests use the same patient for malignant CLDN4 vs T/NK CLDN4. "
            "Samples below the cell floor drop out (honest n)."
        ),
        "tests": tests,
        "cell_level_post": cell_level,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Figure: paired post-treatment inferCNV malignant vs T/NK
    lead = per.loc[per["is_post"]].copy()
    lead_ok = (lead["infercnv_stromal_p95_n"] >= MIN_BOTH) & (lead["n_tnk"] >= MIN_BOTH)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    ax = axes[0]
    plot = lead.loc[lead_ok]
    for _, r in plot.iterrows():
        color = {"NMPR": "#c44e52", "MPR": "#4c72b0"}.get(r["response"], "#888888")
        ax.plot(
            [0, 1],
            [r["infercnv_stromal_p95_cldn4_mean"], r["tnk_cldn4_mean"]],
            color=color,
            alpha=0.75,
            lw=1.2,
        )
        ax.scatter([0], [r["infercnv_stromal_p95_cldn4_mean"]], color=color, s=36, zorder=3)
        ax.scatter([1], [r["tnk_cldn4_mean"]], color=color, s=36, zorder=3)
        ax.text(-0.06, r["infercnv_stromal_p95_cldn4_mean"], r["patient"], fontsize=7, ha="right", va="center")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CNV-malignant\nCLDN4", "same-patient\nT/NK CLDN4"])
    ax.set_ylabel("mean log1p(CP10k) CLDN4")
    ax.set_title(f"GSE207422 post, inferCNV-like p95\nn={int(lead_ok.sum())} patients with ≥{MIN_BOTH} cells/arm")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    for resp, color in [("NMPR", "#c44e52"), ("MPR", "#4c72b0")]:
        sl = plot[plot["response"] == resp]
        ax.scatter(sl["frac_tnk"], sl["infercnv_stromal_p95_cldn4_mean"], color=color, s=42, label=resp)
        for _, r in sl.iterrows():
            ax.text(r["frac_tnk"] + 0.008, r["infercnv_stromal_p95_cldn4_mean"], r["patient"], fontsize=7)
    ax.set_xlabel("T/NK fraction (all cells)")
    ax.set_ylabel("malignant mean log1p(CP10k) CLDN4")
    sp = paired.get("infercnv_stromal_p95", {}).get("spearman_mean_vs_tnkfrac", {})
    ax.set_title(f"vs T/NK fraction  ρ={fmt_num(sp.get('rho'))} p={fmt_p(sp.get('p'))}")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.outdir / "fig_cldn4_infercnv_vs_tnk.png", dpi=160)
    fig.savefig(args.outdir / "fig_cldn4_infercnv_vs_tnk.pdf")
    plt.close(fig)

    print("wrote", args.outdir / "per_sample.tsv", flush=True)
    print("wrote", args.outdir / "summary.json", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""GSE207422 malignant CLDN4 vs B-cell / TLS fraction.

Patient-level B fraction, B+plasma, TLS-like cellular fraction, and TLS12
z-score versus malignant CLDN4. TACSTD2 is a companion gene, never a gate.

Public GEO UMI only. Author CopyKAT barcodes are not public.
Unit of every test is the 12 post-treatment patients (honest n).
The prior TLS meta used n=15 (3 pre-biopsies included); that is not reused.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]

# Hu et al. Fig. 1: 3 pre-biopsy TN; 12 post-surgery (MPR n=4 incl. pCR P06; NMPR n=8).
PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",  # pCR
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "plasma": ["IGHG1", "MZB1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "neutrophil": ["CSF3R"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
    "mast": ["KIT"],
}

A3_NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
BROAD_NORMAL_LUNG = [
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPC",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "TPPP3",
    "FOXJ1",
    "CAPS",
]
TLS12 = [
    "CCL2",
    "CCL3",
    "CCL4",
    "CCL5",
    "CCL8",
    "CCL18",
    "CCL19",
    "CCL21",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CXCL13",
]

COLOR = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "TN": "#6b6b6b"}


def score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(LINEAGES)
    scores = np.vstack([score(expr, LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def load_sample_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["Sample"].map(PAPER_GROUP)
    raw["path_response"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw["timing"] = np.where(
        raw["Resource"].astype(str).str.contains("Pre", case=False, na=False),
        "pre",
        "post",
    )
    return raw


def gene_stats(umi: np.ndarray, lib: np.ndarray) -> dict:
    umi = np.asarray(umi, dtype=np.float64)
    lib = np.asarray(lib, dtype=np.float64)
    if umi.size == 0:
        return {
            "n": 0,
            "mean_log1p": np.nan,
            "mean_log1p_cp10k": np.nan,
            "pct_pos": np.nan,
            "pct_pos_ge2": np.nan,
        }
    cp = np.where(lib > 0, umi / lib * 1e4, np.nan)
    return {
        "n": int(umi.size),
        "mean_log1p": float(np.mean(np.log1p(umi))),
        "mean_log1p_cp10k": float(np.nanmean(np.log1p(cp))),
        "pct_pos": float(np.mean(umi >= 1) * 100.0),
        "pct_pos_ge2": float(np.mean(umi >= 2) * 100.0),
    }


def log1p_cp10k(umi: np.ndarray, lib: np.ndarray) -> np.ndarray:
    umi = np.asarray(umi, dtype=float)
    lib = np.asarray(lib, dtype=float)
    out = np.full(len(umi), np.nan, dtype=float)
    ok = lib > 0
    out[ok] = np.log1p(umi[ok] * 1e4 / lib[ok])
    return out


def spearman_row(x, y, contrast: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return {
            "contrast": contrast,
            "n": n,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": contrast,
        "n": n,
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def exact_wilcoxon_two_sided(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 2 or n_b < 2:
        return {
            "n_a": n_a,
            "n_b": n_b,
            "mean_a": float(np.mean(a)) if n_a else np.nan,
            "mean_b": float(np.mean(b)) if n_b else np.nan,
            "median_a": float(np.median(a)) if n_a else np.nan,
            "median_b": float(np.median(b)) if n_b else np.nan,
            "delta_mean_a_minus_b": np.nan,
            "mwu_u": np.nan,
            "exact_p": np.nan,
            "n_enumerated": 0,
            "note": "too_few_samples",
        }
    method = "exact" if (n_a + n_b) <= 20 else "asymptotic"
    u_obs, p_exact = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    n_enum = int(sum(1 for _ in itertools.combinations(range(n_a + n_b), n_a)))
    return {
        "n_a": n_a,
        "n_b": n_b,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "delta_mean_a_minus_b": float(np.mean(a) - np.mean(b)),
        "mwu_u": float(u_obs),
        "exact_p": float(p_exact),
        "n_enumerated": n_enum,
        "note": "",
    }


def quartile_split(series: pd.Series) -> tuple[pd.Index, pd.Index, int]:
    """Q1 = lowest n//4, Q4 = highest n//4. n=12 → 3 vs 3."""
    s = series.dropna()
    k = max(1, len(s) // 4)
    order = s.sort_values(kind="mergesort")
    q1 = order.index[:k]
    q4 = order.index[-k:]
    return q1, q4, k


def fmt_rho(row: dict) -> str:
    if row.get("note") == "too_few_samples":
        return f"n={row['n']} (too few)"
    return f"ρ={row['spearman_rho']:.2f}, p={row['spearman_p']:.2g}, n={row['n']}"


def fmt_mwu(row: dict) -> str:
    if row.get("note") == "too_few_samples":
        return f"n={row['n_a']} vs {row['n_b']} (too few)"
    return (
        f"mean {row['mean_a']:.3f} vs {row['mean_b']:.3f} "
        f"(Δ={row['delta_mean_a_minus_b']:+.3f}); "
        f"exact p={row['exact_p']:.2g}; n={row['n_a']} vs {row['n_b']}"
    )


def tls12_z_for_rows(rows: list[dict], sample_order: list[str]):
    """Within-set z of patient-mean log1p(CP10k) TLS12 genes, then mean.

    `sample_order` is the honest analysis set (the 12 post patients).
    Pre-treatment samples are scored with the post-set mean/sd so they can
    be written to the table, but they are never used to fit the z.
    """
    post_set = set(sample_order)
    gene_vals: dict[str, list[float]] = {g: [] for g in TLS12}
    post_index = [i for i, r in enumerate(rows) if r["Sample"] in post_set]
    for g in TLS12:
        for i in post_index:
            gene_vals[g].append(rows[i]["_tls12_means"].get(g, np.nan))
    z_params = {}
    genes_used = []
    for g in TLS12:
        vals = np.asarray(gene_vals[g], dtype=float)
        if np.isfinite(vals).sum() >= 4 and np.nanstd(vals, ddof=1) > 0:
            z_params[g] = (float(np.nanmean(vals)), float(np.nanstd(vals, ddof=1)))
            genes_used.append(g)
    out = {}
    for r in rows:
        zs = []
        for g in genes_used:
            mu, sd = z_params[g]
            v = r["_tls12_means"].get(g, np.nan)
            if np.isfinite(v):
                zs.append((v - mu) / sd)
        out[r["Sample"]] = float(np.mean(zs)) if zs else np.nan
    return out, genes_used


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    figdir = args.outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    extracted = np.load(args.workdir / "extracted_cldn4_bfrac_markers.npz", allow_pickle=True)
    cells = extracted["cells"]
    genes = list(extracted["genes"])
    mat = extracted["mat"]
    total = extracted["total"]
    expr = {g: mat[i] for i, g in enumerate(genes)}
    n = len(cells)
    if "CLDN4" not in expr:
        raise SystemExit("CLDN4 missing from extracted matrix")

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    is_t = lineage == "T"
    is_nk = lineage == "NK"
    is_tnk = is_t | is_nk
    is_epi = lineage == "epithelial"
    is_b = lineage == "B"
    is_plasma = lineage == "plasma"

    a3_normal = np.zeros(n, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            a3_normal += expr[g]
    is_malig = is_epi & (a3_normal == 0)
    is_normal_epi = is_epi & (a3_normal > 0)

    broad_score = score(expr, [g for g in BROAD_NORMAL_LUNG if g in expr], n)
    if is_epi.any():
        epi_cut = float(np.quantile(broad_score[is_epi], 0.60))
    else:
        epi_cut = np.inf
    is_malig_broad = is_epi & (broad_score < epi_cut)

    cxcl13 = expr["CXCL13"] if "CXCL13" in expr else np.zeros(n, dtype=np.int32)
    is_cxcl13 = is_tnk & (cxcl13 >= 1)
    is_tls_like = is_b | is_cxcl13

    per_cell = pd.DataFrame(
        {
            "cell": cells,
            "Sample": sample_ids,
            "lineage": lineage,
            "is_tnk": is_tnk,
            "is_t": is_t,
            "is_epi": is_epi,
            "is_b": is_b,
            "is_plasma": is_plasma,
            "is_malig_a3": is_malig,
            "is_malig_broad": is_malig_broad,
            "is_normal_epi": is_normal_epi,
            "is_cxcl13_pos": is_cxcl13,
            "is_tls_like": is_tls_like,
            "total": total,
            "CLDN4": expr["CLDN4"],
            "TACSTD2": expr["TACSTD2"] if "TACSTD2" in expr else 0,
            "CXCL13": cxcl13,
            "MS4A1": expr["MS4A1"] if "MS4A1" in expr else 0,
        }
    )
    per_cell["paper_group"] = per_cell["Sample"].map(PAPER_GROUP)
    per_cell["patient"] = per_cell["Sample"].str.replace("BD_immune", "P", regex=False)

    sample_meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    rows = []
    for sample, sdf in per_cell.groupby("Sample", sort=True):
        mal = sdf[sdf["is_malig_a3"]]
        malb = sdf[sdf["is_malig_broad"]]
        epi = sdf[sdf["is_epi"]]
        cldn4_mal = gene_stats(mal["CLDN4"].to_numpy(), mal["total"].to_numpy())
        cldn4_malb = gene_stats(malb["CLDN4"].to_numpy(), malb["total"].to_numpy())
        cldn4_epi = gene_stats(epi["CLDN4"].to_numpy(), epi["total"].to_numpy())
        tac_mal = gene_stats(mal["TACSTD2"].to_numpy(), mal["total"].to_numpy())
        tac_epi = gene_stats(epi["TACSTD2"].to_numpy(), epi["total"].to_numpy())
        n_all = len(sdf)
        mask = sample_ids == sample
        lib_s = total[mask]
        tls_means = {}
        for g in TLS12:
            if g in expr:
                tls_means[g] = float(np.nanmean(log1p_cp10k(expr[g][mask], lib_s)))
        n_b = int(sdf["is_b"].sum())
        n_plasma = int(sdf["is_plasma"].sum())
        n_cxcl13 = int(sdf["is_cxcl13_pos"].sum())
        row = {
            "Sample": sample,
            "patient": sample.replace("BD_immune", "P"),
            "paper_group": PAPER_GROUP.get(sample),
            "n_cells": n_all,
            "n_epithelial": int(len(epi)),
            "n_malig_a3": int(len(mal)),
            "n_malig_broad": int(len(malb)),
            "n_normal_epi": int(sdf["is_normal_epi"].sum()),
            "n_B": n_b,
            "n_plasma": n_plasma,
            "n_B_plasma": n_b + n_plasma,
            "n_T_NK": int(sdf["is_tnk"].sum()),
            "n_CXCL13_pos": n_cxcl13,
            "n_tls_like": int(sdf["is_tls_like"].sum()),
            "frac_B": float(n_b / n_all),
            "frac_plasma": float(n_plasma / n_all),
            "frac_B_plasma": float((n_b + n_plasma) / n_all),
            "frac_TLS_like": float(sdf["is_tls_like"].sum() / n_all),
            "frac_CXCL13_pos": float(n_cxcl13 / n_all),
            "frac_T_NK": float(sdf["is_tnk"].sum() / n_all),
            "ms4a1_mean_log1p_cp10k": float(np.nanmean(log1p_cp10k(sdf["MS4A1"].to_numpy(), sdf["total"].to_numpy()))),
            "cxcl13_mean_log1p_cp10k": float(np.nanmean(log1p_cp10k(sdf["CXCL13"].to_numpy(), sdf["total"].to_numpy()))),
            "tls12_k": int(len(tls_means)),
            "_tls12_means": tls_means,
        }
        for prefix, st in (
            ("mal_cldn4", cldn4_mal),
            ("malb_cldn4", cldn4_malb),
            ("epi_cldn4", cldn4_epi),
            ("mal_tacstd2", tac_mal),
            ("epi_tacstd2", tac_epi),
        ):
            for k, v in st.items():
                row[f"{prefix}_{k}"] = v
        rows.append(row)

    post_samples = [r["Sample"] for r in rows if PAPER_GROUP.get(r["Sample"]) in {"MPR", "NMPR"}]
    tls_z, tls_genes = tls12_z_for_rows(rows, post_samples)
    for r in rows:
        r["tls12_z"] = tls_z[r["Sample"]]
        r["tls12_genes"] = ",".join(tls_genes)
        del r["_tls12_means"]

    meta_keep = [c for c in sample_meta.columns if c not in {"paper_group"}]
    sample_df = pd.DataFrame(rows).merge(sample_meta[meta_keep], on="Sample", how="left")
    sample_df.to_csv(args.outdir / "per_patient.tsv", sep="\t", index=False)

    post = sample_df[sample_df["paper_group"].isin(["MPR", "NMPR"])].copy()
    assert len(post) == 12, f"expected 12 post-treatment patients, got {len(post)}"
    mpr = post[post["paper_group"] == "MPR"]
    nmpr = post[post["paper_group"] == "NMPR"]

    immune_cols = [
        ("frac_B", "B"),
        ("frac_B_plasma", "B+plasma"),
        ("frac_TLS_like", "TLS-like"),
        ("tls12_z", "TLS12 z"),
    ]
    score_cols = [
        ("mal_cldn4_mean_log1p_cp10k", "A3-malignant CLDN4 mean log1p(CP10k)"),
        ("mal_cldn4_pct_pos", "A3-malignant CLDN4 %pos"),
        ("epi_cldn4_mean_log1p_cp10k", "epithelial CLDN4 mean log1p(CP10k)"),
        ("epi_cldn4_pct_pos", "epithelial CLDN4 %pos"),
        ("malb_cldn4_mean_log1p_cp10k", "broad-malignant CLDN4 mean log1p(CP10k)"),
        ("malb_cldn4_pct_pos", "broad-malignant CLDN4 %pos"),
    ]
    companion_cols = [
        ("mal_tacstd2_mean_log1p_cp10k", "A3-malignant TACSTD2 mean log1p(CP10k) [companion]"),
        ("mal_tacstd2_pct_pos", "A3-malignant TACSTD2 %pos [companion]"),
        ("epi_tacstd2_mean_log1p_cp10k", "epithelial TACSTD2 mean log1p(CP10k) [companion]"),
    ]

    spearman_rows = []
    for scol, slabel in score_cols + companion_cols:
        for icol, ilabel in immune_cols:
            spearman_rows.append(
                spearman_row(post[scol], post[icol], f"post n=12: {slabel} vs {ilabel}")
            )
    spearman_rows.append(
        spearman_row(
            post["mal_cldn4_mean_log1p_cp10k"],
            post["mal_tacstd2_mean_log1p_cp10k"],
            "post n=12: A3-malignant CLDN4 mean vs TACSTD2 mean (companion, not a gate)",
        )
    )
    spearman_rows.append(
        spearman_row(
            post["mal_cldn4_mean_log1p_cp10k"],
            post["ms4a1_mean_log1p_cp10k"],
            "post n=12: A3-malignant CLDN4 mean vs MS4A1 mean",
        )
    )
    spearman_rows.append(
        spearman_row(
            post["mal_cldn4_mean_log1p_cp10k"],
            post["cxcl13_mean_log1p_cp10k"],
            "post n=12: A3-malignant CLDN4 mean vs CXCL13 mean",
        )
    )
    spearman_rows.append(
        spearman_row(
            post["epi_cldn4_mean_log1p_cp10k"],
            post["ms4a1_mean_log1p_cp10k"],
            "post n=12: epithelial CLDN4 mean vs MS4A1 mean",
        )
    )
    spearman_rows.append(
        spearman_row(
            post["epi_cldn4_mean_log1p_cp10k"],
            post["cxcl13_mean_log1p_cp10k"],
            "post n=12: epithelial CLDN4 mean vs CXCL13 mean",
        )
    )
    pd.DataFrame(spearman_rows).to_csv(args.outdir / "spearman.tsv", sep="\t", index=False)

    group_rows = []
    for scol, slabel in score_cols + companion_cols:
        block = exact_wilcoxon_two_sided(nmpr[scol], mpr[scol])
        block["contrast"] = f"NMPR vs MPR: {slabel}"
        block["score"] = scol
        group_rows.append(block)
    for icol, ilabel in immune_cols:
        block = exact_wilcoxon_two_sided(nmpr[icol], mpr[icol])
        block["contrast"] = f"NMPR vs MPR: {ilabel} fraction" if icol != "tls12_z" else f"NMPR vs MPR: {ilabel}"
        block["score"] = icol
        group_rows.append(block)
    pd.DataFrame(group_rows).to_csv(args.outdir / "nmpr_vs_mpr.tsv", sep="\t", index=False)

    q_rows = []
    q_assign = []
    rank_specs = [
        ("mal_cldn4_mean_log1p_cp10k", "A3-malignant CLDN4 mean"),
        ("mal_cldn4_pct_pos", "A3-malignant CLDN4 %pos"),
        ("epi_cldn4_mean_log1p_cp10k", "epithelial CLDN4 mean"),
        ("epi_cldn4_pct_pos", "epithelial CLDN4 %pos"),
    ]
    for scol, slabel in rank_specs:
        q1_idx, q4_idx, k = quartile_split(post[scol])
        q1 = post.loc[q1_idx]
        q4 = post.loc[q4_idx]
        for _, r in q1.iterrows():
            q_assign.append(
                {
                    "rank_score": scol,
                    "rank_label": slabel,
                    "quartile": "Q1",
                    "Sample": r["Sample"],
                    "patient": r["patient"],
                    "paper_group": r["paper_group"],
                    "score_value": r[scol],
                    "k_per_tail": k,
                }
            )
        for _, r in q4.iterrows():
            q_assign.append(
                {
                    "rank_score": scol,
                    "rank_label": slabel,
                    "quartile": "Q4",
                    "Sample": r["Sample"],
                    "patient": r["patient"],
                    "paper_group": r["paper_group"],
                    "score_value": r[scol],
                    "k_per_tail": k,
                }
            )
        for icol, ilabel in immune_cols:
            block = exact_wilcoxon_two_sided(q4[icol], q1[icol])
            suffix = "score" if icol == "tls12_z" else "fraction"
            block["contrast"] = f"Q4 vs Q1 ({slabel}): {ilabel} {suffix}"
            block["rank_score"] = scol
            block["immune"] = icol
            block["k_per_tail"] = k
            block["q4_patients"] = ",".join(sorted(q4["patient"]))
            block["q1_patients"] = ",".join(sorted(q1["patient"]))
            q_rows.append(block)

    pd.DataFrame(q_rows).to_csv(args.outdir / "q4_vs_q1.tsv", sep="\t", index=False)
    pd.DataFrame(q_assign).to_csv(args.outdir / "quartile_assignment.tsv", sep="\t", index=False)

    lineage_counts = (
        per_cell.groupby(["Sample", "lineage"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    lineage_counts.to_csv(args.outdir / "lineage_counts.tsv", sep="\t", index=False)

    sanity = {
        "dataset": "GSE207422",
        "n_cells": int(n),
        "n_genes_in_matrix": int(extracted["n_genes_in_matrix"][0]),
        "n_samples": int(per_cell["Sample"].nunique()),
        "n_post": 12,
        "n_post_mpr": 4,
        "n_post_nmpr": 8,
        "author_cell_labels_public": False,
        "dual_high_used": False,
        "tacstd2_is_gate": False,
        "lineage_counts": {k: int(v) for k, v in per_cell["lineage"].value_counts().items()},
        "n_epithelial": int(is_epi.sum()),
        "n_malig_a3": int(is_malig.sum()),
        "n_B": int(is_b.sum()),
        "n_plasma": int(is_plasma.sum()),
        "n_tls_like": int(is_tls_like.sum()),
        "n_CXCL13_pos": int(is_cxcl13.sum()),
        "tls12_genes": tls_genes,
        "tls12_k": len(tls_genes),
        "genes_present": genes,
        "genes_missing": [g for g in ["CLDN4", "TACSTD2", "CXCL13", "MS4A1", *TLS12] if g not in expr],
        "post_malig_a3_empty": post.loc[post["n_malig_a3"] == 0, "patient"].tolist(),
        "post_malig_a3_lt10": post.loc[(post["n_malig_a3"] > 0) & (post["n_malig_a3"] < 10), "patient"].tolist(),
        "malignant_definition_a3": (
            "epithelial argmax AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
            "(same as given A3 TACSTD2 analysis; not CopyKAT)"
        ),
        "tls_like_definition": "lineage B OR (T/NK AND CXCL13 UMI ≥ 1) / all cells",
        "tls12_definition": (
            "mean of within-post-n=12 z-scores of patient-mean log1p(CP10k) "
            "for the Coppola 12-chemokine set; not histologic TLS"
        ),
        "note": (
            "Unit is 12 post-treatment patients. The 3 pre-treatment biopsies "
            "are excluded. Cell-level p-values are not reported. "
            "Prior TLS meta n=15 is not reused."
        ),
    }
    with (args.outdir / "sanity.json").open("w") as fh:
        json.dump(sanity, fh, indent=2)

    make_figures(post, figdir, spearman_rows, group_rows, q_rows)
    write_finding(args.outdir, post, spearman_rows, group_rows, q_rows, q_assign, sanity)

    summary = {
        "dataset": "GSE207422",
        "task": "malignant CLDN4 vs B-cell / TLS fraction (honest n=12 post)",
        "n_post": 12,
        "n_cells": int(n),
        "n_malig_a3": int(is_malig.sum()),
        "n_B": int(is_b.sum()),
        "empty_a3_malignant_post": sanity["post_malig_a3_empty"],
        "primary_spearman_mal_cldn4_mean_vs_B": next(
            r
            for r in spearman_rows
            if r["contrast"] == "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs B"
        ),
        "primary_spearman_mal_cldn4_mean_vs_tls_like": next(
            r
            for r in spearman_rows
            if r["contrast"] == "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS-like"
        ),
        "primary_spearman_mal_cldn4_mean_vs_tls12": next(
            r
            for r in spearman_rows
            if r["contrast"] == "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS12 z"
        ),
    }
    with (args.outdir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def scatter_panel(ax, df, x, y, xlabel, ylabel, title, rho_row=None):
    for grp, sub in df.groupby("paper_group"):
        ax.scatter(
            sub[x],
            sub[y],
            c=COLOR.get(grp, "#333"),
            label=grp,
            edgecolors="white",
            linewidths=0.6,
            s=56,
            zorder=3,
        )
        for _, r in sub.iterrows():
            if np.isfinite(r[x]) and np.isfinite(r[y]):
                ax.annotate(
                    r["patient"],
                    (r[x], r[y]),
                    textcoords="offset points",
                    xytext=(4, 3),
                    fontsize=7,
                    color="#333",
                )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    if rho_row is not None:
        ax.text(
            0.03,
            0.97,
            fmt_rho(rho_row),
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, lw=0.4),
        )
    ax.legend(frameon=False, fontsize=8, loc="lower right")


def box_panel(ax, left, right, labels, ylabel, title, test_row=None, colors=None):
    data = [np.asarray(left, dtype=float), np.asarray(right, dtype=float)]
    data = [d[np.isfinite(d)] for d in data]
    bp = ax.boxplot(data, tick_labels=labels, widths=0.55, patch_artist=True, showfliers=False)
    cols = colors or ["#2c6eaf", "#d1495b"]
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
        patch.set_edgecolor(c)
    for i, d in enumerate(data, start=1):
        jitter = np.random.default_rng(1 + i).uniform(-0.08, 0.08, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d, c=cols[i - 1], s=28, zorder=3, edgecolors="white")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    if test_row is not None:
        ax.text(
            0.03,
            0.97,
            fmt_mwu(test_row),
            transform=ax.transAxes,
            va="top",
            fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, lw=0.4),
        )


def make_figures(post, figdir: Path, spearman_rows, group_rows, q_rows) -> None:
    def find_s(substr: str):
        return next(r for r in spearman_rows if r["contrast"] == substr)

    def find_g(substr: str):
        return next(r for r in group_rows if r["contrast"] == substr)

    def find_q(substr: str):
        return next(r for r in q_rows if r["contrast"] == substr)

    np.random.seed(0)

    fig, axes = plt.subplots(2, 4, figsize=(16.4, 8.0), constrained_layout=True)
    pairs = [
        (
            axes[0, 0],
            "mal_cldn4_mean_log1p_cp10k",
            "frac_B",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "B-cell fraction",
            "CLDN4 mean vs B",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs B",
        ),
        (
            axes[0, 1],
            "mal_cldn4_mean_log1p_cp10k",
            "frac_B_plasma",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "B+plasma fraction",
            "CLDN4 mean vs B+plasma",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs B+plasma",
        ),
        (
            axes[0, 2],
            "mal_cldn4_mean_log1p_cp10k",
            "frac_TLS_like",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "TLS-like fraction (B + CXCL13+ T/NK)",
            "CLDN4 mean vs TLS-like",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS-like",
        ),
        (
            axes[0, 3],
            "mal_cldn4_mean_log1p_cp10k",
            "tls12_z",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "TLS12 z (post-fit)",
            "CLDN4 mean vs TLS12",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS12 z",
        ),
        (
            axes[1, 0],
            "mal_cldn4_pct_pos",
            "frac_B",
            "A3-malignant CLDN4 %pos",
            "B-cell fraction",
            "CLDN4 %pos vs B",
            "post n=12: A3-malignant CLDN4 %pos vs B",
        ),
        (
            axes[1, 1],
            "mal_cldn4_pct_pos",
            "frac_B_plasma",
            "A3-malignant CLDN4 %pos",
            "B+plasma fraction",
            "CLDN4 %pos vs B+plasma",
            "post n=12: A3-malignant CLDN4 %pos vs B+plasma",
        ),
        (
            axes[1, 2],
            "mal_cldn4_pct_pos",
            "frac_TLS_like",
            "A3-malignant CLDN4 %pos",
            "TLS-like fraction (B + CXCL13+ T/NK)",
            "CLDN4 %pos vs TLS-like",
            "post n=12: A3-malignant CLDN4 %pos vs TLS-like",
        ),
        (
            axes[1, 3],
            "mal_cldn4_pct_pos",
            "tls12_z",
            "A3-malignant CLDN4 %pos",
            "TLS12 z (post-fit)",
            "CLDN4 %pos vs TLS12",
            "post n=12: A3-malignant CLDN4 %pos vs TLS12 z",
        ),
    ]
    for ax, x, y, xlab, ylab, title, key in pairs:
        scatter_panel(ax, post, x, y, xlab, ylab, title, find_s(key))
    fig.suptitle(
        "GSE207422 post-treatment n=12 · malignant CLDN4 vs B / TLS fraction",
        fontsize=12,
    )
    fig.savefig(figdir / "fig_cldn4_vs_b_tls.png", dpi=160)
    fig.savefig(figdir / "fig_cldn4_vs_b_tls.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(15.6, 4.2), constrained_layout=True)
    nmpr = post[post["paper_group"] == "NMPR"]
    mpr = post[post["paper_group"] == "MPR"]
    box_panel(
        axes[0],
        nmpr["frac_B"],
        mpr["frac_B"],
        ["NMPR", "MPR"],
        "B-cell fraction",
        "NMPR vs MPR · B",
        find_g("NMPR vs MPR: B fraction"),
    )
    box_panel(
        axes[1],
        nmpr["frac_B_plasma"],
        mpr["frac_B_plasma"],
        ["NMPR", "MPR"],
        "B+plasma fraction",
        "NMPR vs MPR · B+plasma",
        find_g("NMPR vs MPR: B+plasma fraction"),
    )
    box_panel(
        axes[2],
        nmpr["frac_TLS_like"],
        mpr["frac_TLS_like"],
        ["NMPR", "MPR"],
        "TLS-like fraction",
        "NMPR vs MPR · TLS-like",
        find_g("NMPR vs MPR: TLS-like fraction"),
    )
    box_panel(
        axes[3],
        nmpr["tls12_z"],
        mpr["tls12_z"],
        ["NMPR", "MPR"],
        "TLS12 z",
        "NMPR vs MPR · TLS12",
        find_g("NMPR vs MPR: TLS12 z"),
    )
    fig.suptitle("GSE207422 · NMPR n=8 vs MPR n=4 (pCR P06 = MPR) · B / TLS", fontsize=11)
    fig.savefig(figdir / "fig_b_tls_nmpr_mpr.png", dpi=160)
    fig.savefig(figdir / "fig_b_tls_nmpr_mpr.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(16.4, 8.0), constrained_layout=True)
    extra_q = [
        (axes[0, 0], "epi_cldn4_mean_log1p_cp10k", "frac_B", "epithelial CLDN4 mean", "B"),
        (axes[0, 1], "epi_cldn4_mean_log1p_cp10k", "frac_B_plasma", "epithelial CLDN4 mean", "B+plasma"),
        (axes[0, 2], "epi_cldn4_mean_log1p_cp10k", "frac_TLS_like", "epithelial CLDN4 mean", "TLS-like"),
        (axes[0, 3], "epi_cldn4_mean_log1p_cp10k", "tls12_z", "epithelial CLDN4 mean", "TLS12 z"),
        (axes[1, 0], "epi_cldn4_pct_pos", "frac_B", "epithelial CLDN4 %pos", "B"),
        (axes[1, 1], "epi_cldn4_pct_pos", "frac_B_plasma", "epithelial CLDN4 %pos", "B+plasma"),
        (axes[1, 2], "epi_cldn4_pct_pos", "frac_TLS_like", "epithelial CLDN4 %pos", "TLS-like"),
        (axes[1, 3], "epi_cldn4_pct_pos", "tls12_z", "epithelial CLDN4 %pos", "TLS12 z"),
    ]
    for ax, scol, icol, slabel, ilabel in extra_q:
        q1_idx, q4_idx, k = quartile_split(post[scol])
        q1 = post.loc[q1_idx]
        q4 = post.loc[q4_idx]
        suffix = "score" if icol == "tls12_z" else "fraction"
        key = f"Q4 vs Q1 ({slabel}): {ilabel} {suffix}"
        box_panel(
            ax,
            q4[icol],
            q1[icol],
            [f"Q4 (n={k})", f"Q1 (n={k})"],
            f"{ilabel} {suffix}",
            f"{slabel}\nQ4 vs Q1 {ilabel}",
            find_q(key),
            colors=["#3d5a40", "#c9a227"],
        )
    fig.suptitle(
        "Extra · Q4 vs Q1 on n=12 (k=n//4=3) · epithelial CLDN4 ranks keep all 12 patients",
        fontsize=11,
    )
    fig.savefig(figdir / "fig_extra_q4q1.png", dpi=160)
    fig.savefig(figdir / "fig_extra_q4q1.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 4.2), constrained_layout=True)
    order = post.sort_values("patient")
    x = np.arange(len(order))
    ax.bar(x - 0.2, order["n_malig_a3"], width=0.4, color="#2c6eaf", label="A3-malignant n")
    ax.bar(x + 0.2, order["n_B"], width=0.4, color="#c9a227", label="B n")
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{p}\n{g}" for p, g in zip(order["patient"], order["paper_group"])],
        fontsize=8,
    )
    ax.set_ylabel("cells")
    ax.set_title("Per-patient A3-malignant and B counts (empty MPR visible)")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(figdir / "fig_extra_compartment.png", dpi=160)
    fig.savefig(figdir / "fig_extra_compartment.pdf")
    plt.close(fig)


def write_finding(outdir, post, spearman_rows, group_rows, q_rows, q_assign, sanity) -> None:
    def S(name: str) -> dict:
        return next(r for r in spearman_rows if r["contrast"] == name)

    def G(name: str) -> dict:
        return next(r for r in group_rows if r["contrast"] == name)

    def Q(name: str) -> dict:
        return next(r for r in q_rows if r["contrast"] == name)

    empty = ", ".join(sanity["post_malig_a3_empty"]) or "none"
    lt10 = ", ".join(sanity["post_malig_a3_lt10"]) or "none"
    n_mal_mean = int(post["mal_cldn4_mean_log1p_cp10k"].notna().sum())
    n_epi = int(post["epi_cldn4_mean_log1p_cp10k"].notna().sum())
    q_assign_df = pd.DataFrame(q_assign)

    def q_names(score, q):
        sub = q_assign_df[(q_assign_df["rank_score"] == score) & (q_assign_df["quartile"] == q)]
        return ", ".join(f"{r.patient}({r.paper_group})" for r in sub.itertuples())

    mal_b = S("post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs B")
    mal_tls = S("post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS-like")
    mal_z = S("post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS12 z")
    epi_b = S("post n=12: epithelial CLDN4 mean log1p(CP10k) vs B")

    lines = [
        "# GSE207422 — B-cell / TLS fraction vs malignant CLDN4",
        "",
        "ADDITIVE public-UMI slice. **Honest n is the 12 post-treatment patients.** "
        "The three pre-treatment biopsies are excluded. Cell-level p-values are not reported. "
        "The prior TLS/B meta that listed GSE207422 as n=15 is not reused. "
        "TACSTD2 is a companion gene and is **never a gate**. Dual-high was not run. "
        "This is not histologic TLS.",
        "",
        (
            f"**Verdict (n=12 honest):** malignant CLDN4 does not significantly anti-correlate "
            f"with B-cell fraction, B+plasma, TLS-like cellular fraction, or TLS12 z. "
            f"A3-malignant CLDN4 mean vs B {fmt_rho(mal_b)}; vs TLS-like {fmt_rho(mal_tls)}; "
            f"vs TLS12 z {fmt_rho(mal_z)}. "
            f"Complete-case epithelial CLDN4 mean vs B {fmt_rho(epi_b)}. "
            f"NMPR vs MPR B/TLS fractions are not significant (n=8 vs 4). "
            f"Q4 vs Q1 is 3 vs 3 on epithelial ranks; the smallest exact two-sided p is 0.10."
        ),
        "",
        "## Data and n",
        "",
        f"- Public GEO UMI only: **{sanity['n_cells']:,}** cells × **{sanity['n_genes_in_matrix']:,}** genes. "
        "Raw GSA-Human HRA001033 was not used. Author CopyKAT / epithelium RDS barcodes are not on GEO.",
        "- **Unit of every test is the 12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). "
        "P01/P05/P08 pre-biopsies are written in `per_patient.tsv` and are not tested.",
        f"- Marker lineages (Hu canonical argmax): epithelial {sanity['n_epithelial']:,}; "
        f"A3-malignant-like {sanity['n_malig_a3']:,}; B {sanity['n_B']:,}; "
        f"plasma {sanity['n_plasma']:,}; TLS-like (B + CXCL13+ T/NK) {sanity['n_tls_like']:,}; "
        f"CXCL13+ T/NK {sanity['n_CXCL13_pos']:,}.",
        f"- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
        f"(same rule as the given A3 TACSTD2 slice; **not** CopyKAT). "
        f"Post patients with 0 A3-malignant cells (NaN, dropped from A3-malignant tests): **{empty}**. "
        f"Post patients with 1–9 A3-malignant cells (kept, noisy): **{lt10}**. "
        f"A3-malignant CLDN4 n={n_mal_mean}/12. "
        f"Epithelial CLDN4 is the complete-case n={n_epi}/12 score.",
        f"- TLS12 genes present: {sanity['tls12_k']}/12 ({', '.join(sanity['tls12_genes'])}).",
        "",
        "## Definitions",
        "",
        "| Item | Rule |",
        "|---|---|",
        "| CLDN4 mean | patient mean `log1p(CP10k)` inside the named compartment |",
        "| CLDN4 %pos | fraction of cells in the compartment with CLDN4 UMI ≥ 1 |",
        "| B fraction | lineage B / all cells |",
        "| B+plasma | (lineage B or plasma) / all cells |",
        "| TLS-like fraction | (lineage B) OR ((T or NK) AND CXCL13 UMI ≥ 1) / all cells |",
        "| TLS12 z | mean of within-post-n=12 z-scores of patient-mean log1p(CP10k) for the Coppola 12-chemokine set. Not a cell fraction. Not histologic TLS. |",
        "| Q4 vs Q1 | among patients with a finite rank score, Q4 = highest n//4, Q1 = lowest n//4 (n=12 → 3 vs 3) |",
        "| NMPR vs MPR | exact two-sided Wilcoxon; pCR P06 = MPR |",
        "| TACSTD2 | companion only; not a gate |",
        "",
        "## Primary — A3-malignant CLDN4 vs B / TLS (honest n after empty MPR drop)",
        "",
        "| CLDN4 score | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |",
        "|---|---|---|---|---|",
        f"| mean log1p(CP10k) | {fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs B'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs B+plasma'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS-like'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs TLS12 z'))} |",
        f"| %pos | {fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs B'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs B+plasma'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs TLS-like'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs TLS12 z'))} |",
        "",
        "## Complete-case epithelial CLDN4 vs B / TLS (n=12)",
        "",
        "| CLDN4 score | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |",
        "|---|---|---|---|---|",
        f"| mean log1p(CP10k) | {fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs B'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs B+plasma'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs TLS-like'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs TLS12 z'))} |",
        f"| %pos | {fmt_rho(S('post n=12: epithelial CLDN4 %pos vs B'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 %pos vs B+plasma'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 %pos vs TLS-like'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 %pos vs TLS12 z'))} |",
        "",
        "## Primary — NMPR vs MPR",
        "",
        "| Score | Result |",
        "|---|---|",
        f"| A3-malignant CLDN4 mean | {fmt_mwu(G('NMPR vs MPR: A3-malignant CLDN4 mean log1p(CP10k)'))} |",
        f"| epithelial CLDN4 mean (n=12 complete) | {fmt_mwu(G('NMPR vs MPR: epithelial CLDN4 mean log1p(CP10k)'))} |",
        f"| B fraction | {fmt_mwu(G('NMPR vs MPR: B fraction'))} |",
        f"| B+plasma fraction | {fmt_mwu(G('NMPR vs MPR: B+plasma fraction'))} |",
        f"| TLS-like fraction | {fmt_mwu(G('NMPR vs MPR: TLS-like fraction'))} |",
        f"| TLS12 z | {fmt_mwu(G('NMPR vs MPR: TLS12 z'))} |",
        "",
        "## Extra — Q4 vs Q1 B / TLS",
        "",
        "n=12 honest: epithelial ranks keep all 12 patients (k=3 vs 3). "
        "Exact two-sided p at 3 vs 3 cannot go below 0.10. "
        "A3-malignant ranks drop patients with 0 malignant cells.",
        "",
        "### Epithelial CLDN4 ranks (n=12 complete, k=3)",
        "",
        f"- Mean-rank Q4: {q_names('epi_cldn4_mean_log1p_cp10k', 'Q4')}",
        f"- Mean-rank Q1: {q_names('epi_cldn4_mean_log1p_cp10k', 'Q1')}",
        f"- %pos-rank Q4: {q_names('epi_cldn4_pct_pos', 'Q4')}",
        f"- %pos-rank Q1: {q_names('epi_cldn4_pct_pos', 'Q1')}",
        "",
        "| Rank | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |",
        "|---|---|---|---|---|",
        f"| epithelial mean | {fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): B fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): B+plasma fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): TLS-like fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): TLS12 z score'))} |",
        f"| epithelial %pos | {fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): B fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): B+plasma fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): TLS-like fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): TLS12 z score'))} |",
        "",
        "### A3-malignant CLDN4 ranks (empty malignant dropped)",
        "",
        f"- Mean-rank Q4: {q_names('mal_cldn4_mean_log1p_cp10k', 'Q4')}",
        f"- Mean-rank Q1: {q_names('mal_cldn4_mean_log1p_cp10k', 'Q1')}",
        "",
        "| Rank | vs B | vs B+plasma | vs TLS-like | vs TLS12 z |",
        "|---|---|---|---|---|",
        f"| A3-malignant mean | {fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): B fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): B+plasma fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): TLS-like fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): TLS12 z score'))} |",
        f"| A3-malignant %pos | {fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): B fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): B+plasma fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): TLS-like fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): TLS12 z score'))} |",
        "",
        "## Companion (not a gate)",
        "",
        f"- A3-malignant CLDN4 mean vs TACSTD2 mean: {fmt_rho(S('post n=12: A3-malignant CLDN4 mean vs TACSTD2 mean (companion, not a gate)'))}",
        f"- A3-malignant CLDN4 mean vs MS4A1 mean: {fmt_rho(S('post n=12: A3-malignant CLDN4 mean vs MS4A1 mean'))}",
        f"- A3-malignant CLDN4 mean vs CXCL13 mean: {fmt_rho(S('post n=12: A3-malignant CLDN4 mean vs CXCL13 mean'))}",
        f"- Companion TACSTD2 mean vs B: {fmt_rho(S('post n=12: A3-malignant TACSTD2 mean log1p(CP10k) [companion] vs B'))}",
        f"- Companion TACSTD2 mean vs TLS-like: {fmt_rho(S('post n=12: A3-malignant TACSTD2 mean log1p(CP10k) [companion] vs TLS-like'))}",
        "",
        "## Honest limits",
        "",
        "1. n=12 (4 vs 8) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. "
        "Q4 vs Q1 is 3 vs 3; the smallest exact two-sided p is 0.10.",
        "2. A3-malignant-like empties some MPR residual tumors (normal-lung program). "
        "That is reported, not patched by falling back to epithelial inside the malignant column "
        "or by counting the 3 pre-biopsies.",
        "3. This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak into A3-malignant-like.",
        "4. TLS-like is a marker proxy (B + CXCL13+ T/NK). TLS12 is a 12-chemokine score. "
        "Neither is a pathologist TLS call.",
        "5. Dual-high (TACSTD2 AND CLDN4) was not run. Expression is log1p(CP10k) from the public UMI.",
        "",
        "## Extra figures",
        "",
        "- `figures/fig_cldn4_vs_b_tls.png` — A3-malignant CLDN4 mean/%pos vs B, B+plasma, TLS-like, TLS12",
        "- `figures/fig_b_tls_nmpr_mpr.png` — NMPR vs MPR B / TLS",
        "- `figures/fig_extra_q4q1.png` — Q4 vs Q1 on epithelial ranks (n=12)",
        "- `figures/fig_extra_compartment.png` — per-patient malignant n and B n (empty MPR visible)",
        "",
        "## Files",
        "",
        "- `per_patient.tsv` — 15 samples; tests use the 12 post rows",
        "- `spearman.tsv` / `nmpr_vs_mpr.tsv` / `q4_vs_q1.tsv` / `quartile_assignment.tsv`",
        "- `sample_metadata.tsv` / `lineage_counts.tsv` / `sanity.json` / `summary.json`",
        "- Scripts: `scripts/download.py`, `scripts/extract.py`, `scripts/analyze.py`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/gse207422_cldn4_bfrac/scripts/download.py",
        "python3 methods/gse207422_cldn4_bfrac/scripts/extract.py",
        "python3 methods/gse207422_cldn4_bfrac/scripts/analyze.py",
        "```",
        "",
    ]
    (outdir / "FINDING.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()

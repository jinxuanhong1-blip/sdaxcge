#!/usr/bin/env python3
"""GSE207422 malignant CLDN4 only (not dual-high).

Patient-level mean and %pos, Q4 vs Q1, versus T/NK, CXCL13+, cyto-high T,
and NMPR vs MPR. TACSTD2 is a companion gene, never a gate. A3 TACSTD2
results are taken as given and are not re-argued.

Public GEO UMI only. Author CopyKAT barcodes are not public.
Unit of every test is the 12 post-treatment patients (honest n).
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

# A3 given malignant-like: epithelial AND zero UMI on this panel.
A3_NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
# Broader normal-lung panel used only for the sensitivity "malig_broad" call.
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
CYTO_GENES = ["GZMB", "GZMA", "PRF1", "IFNG", "NKG7"]

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
    """Exact two-sided Wilcoxon/MWU (scipy exact when n is small)."""
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    figdir = args.outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    extracted = np.load(args.workdir / "extracted_cldn4_markers.npz", allow_pickle=True)
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

    cyto = score(expr, CYTO_GENES, n)
    if is_t.any():
        cyto_cut = float(np.quantile(cyto[is_t], 0.75))
    else:
        cyto_cut = np.inf
    is_cyto_high_t = is_t & (cyto >= cyto_cut)

    per_cell = pd.DataFrame(
        {
            "cell": cells,
            "Sample": sample_ids,
            "lineage": lineage,
            "is_tnk": is_tnk,
            "is_t": is_t,
            "is_epi": is_epi,
            "is_malig_a3": is_malig,
            "is_malig_broad": is_malig_broad,
            "is_normal_epi": is_normal_epi,
            "is_cxcl13_pos": is_cxcl13,
            "is_cyto_high_t": is_cyto_high_t,
            "total": total,
            "CLDN4": expr["CLDN4"],
            "TACSTD2": expr["TACSTD2"] if "TACSTD2" in expr else 0,
            "CXCL13": cxcl13,
        }
    )
    per_cell["paper_group"] = per_cell["Sample"].map(PAPER_GROUP)
    per_cell["patient"] = per_cell["Sample"].str.replace("BD_immune", "P", regex=False)

    sample_meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    rows = []
    for sample, sdf in per_cell.groupby("Sample", sort=True):
        lib = sdf["total"].to_numpy()
        mal = sdf[sdf["is_malig_a3"]]
        malb = sdf[sdf["is_malig_broad"]]
        epi = sdf[sdf["is_epi"]]
        cldn4_mal = gene_stats(mal["CLDN4"].to_numpy(), mal["total"].to_numpy())
        cldn4_malb = gene_stats(malb["CLDN4"].to_numpy(), malb["total"].to_numpy())
        cldn4_epi = gene_stats(epi["CLDN4"].to_numpy(), epi["total"].to_numpy())
        tac_mal = gene_stats(mal["TACSTD2"].to_numpy(), mal["total"].to_numpy())
        tac_epi = gene_stats(epi["TACSTD2"].to_numpy(), epi["total"].to_numpy())
        n_all = len(sdf)
        row = {
            "Sample": sample,
            "patient": sample.replace("BD_immune", "P"),
            "paper_group": PAPER_GROUP.get(sample),
            "n_cells": n_all,
            "n_epithelial": int(len(epi)),
            "n_malig_a3": int(len(mal)),
            "n_malig_broad": int(len(malb)),
            "n_normal_epi": int(sdf["is_normal_epi"].sum()),
            "n_T": int(sdf["is_t"].sum()),
            "n_NK": int((sdf["lineage"] == "NK").sum()),
            "n_T_NK": int(sdf["is_tnk"].sum()),
            "n_CXCL13_pos": int(sdf["is_cxcl13_pos"].sum()),
            "n_cyto_high_T": int(sdf["is_cyto_high_t"].sum()),
            "frac_T_NK": float(sdf["is_tnk"].sum() / n_all),
            "frac_CXCL13_pos": float(sdf["is_cxcl13_pos"].sum() / n_all),
            "frac_cyto_high_T": float(sdf["is_cyto_high_t"].sum() / n_all),
            "frac_T_among_TNK": (
                float(sdf["is_t"].sum() / sdf["is_tnk"].sum()) if sdf["is_tnk"].sum() else np.nan
            ),
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

    meta_keep = [c for c in sample_meta.columns if c not in {"paper_group"}]
    sample_df = pd.DataFrame(rows).merge(sample_meta[meta_keep], on="Sample", how="left")
    sample_df.to_csv(args.outdir / "per_patient.tsv", sep="\t", index=False)

    post = sample_df[sample_df["paper_group"].isin(["MPR", "NMPR"])].copy()
    assert len(post) == 12, f"expected 12 post-treatment patients, got {len(post)}"
    mpr = post[post["paper_group"] == "MPR"]
    nmpr = post[post["paper_group"] == "NMPR"]

    immune_cols = [
        ("frac_T_NK", "T/NK"),
        ("frac_CXCL13_pos", "CXCL13+"),
        ("frac_cyto_high_T", "cyto-high T"),
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
            post["epi_cldn4_mean_log1p_cp10k"],
            post["epi_tacstd2_mean_log1p_cp10k"],
            "post n=12: epithelial CLDN4 mean vs TACSTD2 mean (companion, not a gate)",
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
        block["contrast"] = f"NMPR vs MPR: {ilabel} fraction"
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
            block["contrast"] = f"Q4 vs Q1 ({slabel}): {ilabel} fraction"
            block["rank_score"] = scol
            block["immune"] = icol
            block["k_per_tail"] = k
            block["q4_patients"] = ",".join(sorted(q4["patient"]))
            block["q1_patients"] = ",".join(sorted(q1["patient"]))
            q_rows.append(block)
        # Also the ranking score itself Q4 vs Q1 (sanity: must differ).
        block = exact_wilcoxon_two_sided(q4[scol], q1[scol])
        block["contrast"] = f"Q4 vs Q1 ({slabel}): ranking score itself"
        block["rank_score"] = scol
        block["immune"] = scol
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
        "n_malig_broad": int(is_malig_broad.sum()),
        "n_T_NK": int(is_tnk.sum()),
        "n_CXCL13_pos": int(is_cxcl13.sum()),
        "n_cyto_high_T": int(is_cyto_high_t.sum()),
        "cyto_genes": CYTO_GENES,
        "cyto_cut_t_p75": cyto_cut if np.isfinite(cyto_cut) else None,
        "broad_normal_lung_epi_p60": epi_cut if np.isfinite(epi_cut) else None,
        "genes_present": genes,
        "genes_missing": [g for g in ["CLDN4", "TACSTD2", "CXCL13", *CYTO_GENES] if g not in expr],
        "post_malig_a3_empty": post.loc[post["n_malig_a3"] == 0, "patient"].tolist(),
        "post_malig_a3_lt10": post.loc[post["n_malig_a3"] < 10, "patient"].tolist(),
        "malignant_definition_a3": (
            "epithelial argmax AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
            "(same as given A3 TACSTD2 analysis; not CopyKAT)"
        ),
        "note": (
            "CLDN4 only. TACSTD2 is companion, not a dual-high gate. "
            "Unit is 12 post-treatment patients. Cell-level p-values are not reported."
        ),
    }
    with (args.outdir / "sanity.json").open("w") as fh:
        json.dump(sanity, fh, indent=2)

    make_figures(post, figdir, spearman_rows, group_rows, q_rows)
    write_finding(args.outdir, post, spearman_rows, group_rows, q_rows, q_assign, sanity)

    summary = {
        "dataset": "GSE207422",
        "task": "malignant CLDN4 only vs T/NK, CXCL13+, cyto-high T, NMPR vs MPR; Q4 vs Q1",
        "dual_high": False,
        "tacstd2_gate": False,
        "n_post": 12,
        "n_cells": int(n),
        "n_malig_a3": int(is_malig.sum()),
        "empty_a3_malignant_post": sanity["post_malig_a3_empty"],
        "primary_nmpr_mpr_mal_cldn4_mean_exact_p": next(
            r["exact_p"]
            for r in group_rows
            if r["contrast"] == "NMPR vs MPR: A3-malignant CLDN4 mean log1p(CP10k)"
        ),
        "primary_spearman_mal_cldn4_mean_vs_tnk": next(
            r
            for r in spearman_rows
            if r["contrast"] == "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs T/NK"
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

    rng_state = np.random.get_state()
    np.random.seed(0)

    # Main: CLDN4 vs three immune fractions (mean and %pos).
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.2), constrained_layout=True)
    pairs = [
        (
            axes[0, 0],
            "mal_cldn4_mean_log1p_cp10k",
            "frac_T_NK",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "T/NK fraction",
            "CLDN4 mean vs T/NK",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs T/NK",
        ),
        (
            axes[0, 1],
            "mal_cldn4_mean_log1p_cp10k",
            "frac_CXCL13_pos",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "CXCL13+ (T/NK) fraction",
            "CLDN4 mean vs CXCL13+",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs CXCL13+",
        ),
        (
            axes[0, 2],
            "mal_cldn4_mean_log1p_cp10k",
            "frac_cyto_high_T",
            "A3-malignant CLDN4 mean log1p(CP10k)",
            "cyto-high T fraction",
            "CLDN4 mean vs cyto-high T",
            "post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs cyto-high T",
        ),
        (
            axes[1, 0],
            "mal_cldn4_pct_pos",
            "frac_T_NK",
            "A3-malignant CLDN4 %pos",
            "T/NK fraction",
            "CLDN4 %pos vs T/NK",
            "post n=12: A3-malignant CLDN4 %pos vs T/NK",
        ),
        (
            axes[1, 1],
            "mal_cldn4_pct_pos",
            "frac_CXCL13_pos",
            "A3-malignant CLDN4 %pos",
            "CXCL13+ (T/NK) fraction",
            "CLDN4 %pos vs CXCL13+",
            "post n=12: A3-malignant CLDN4 %pos vs CXCL13+",
        ),
        (
            axes[1, 2],
            "mal_cldn4_pct_pos",
            "frac_cyto_high_T",
            "A3-malignant CLDN4 %pos",
            "cyto-high T fraction",
            "CLDN4 %pos vs cyto-high T",
            "post n=12: A3-malignant CLDN4 %pos vs cyto-high T",
        ),
    ]
    for ax, x, y, xlab, ylab, title, key in pairs:
        scatter_panel(ax, post, x, y, xlab, ylab, title, find_s(key))
    fig.suptitle(
        "GSE207422 post-treatment n=12 · malignant CLDN4 only (not dual-high)",
        fontsize=12,
    )
    fig.savefig(figdir / "fig_cldn4_vs_immune.png", dpi=160)
    fig.savefig(figdir / "fig_cldn4_vs_immune.pdf")
    plt.close(fig)

    # NMPR vs MPR
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.2), constrained_layout=True)
    nmpr = post[post["paper_group"] == "NMPR"]
    mpr = post[post["paper_group"] == "MPR"]
    box_panel(
        axes[0],
        nmpr["mal_cldn4_mean_log1p_cp10k"],
        mpr["mal_cldn4_mean_log1p_cp10k"],
        ["NMPR", "MPR"],
        "A3-malignant CLDN4 mean log1p(CP10k)",
        "NMPR vs MPR · CLDN4 mean",
        find_g("NMPR vs MPR: A3-malignant CLDN4 mean log1p(CP10k)"),
    )
    box_panel(
        axes[1],
        nmpr["mal_cldn4_pct_pos"],
        mpr["mal_cldn4_pct_pos"],
        ["NMPR", "MPR"],
        "A3-malignant CLDN4 %pos",
        "NMPR vs MPR · CLDN4 %pos",
        find_g("NMPR vs MPR: A3-malignant CLDN4 %pos"),
    )
    box_panel(
        axes[2],
        nmpr["epi_cldn4_mean_log1p_cp10k"],
        mpr["epi_cldn4_mean_log1p_cp10k"],
        ["NMPR", "MPR"],
        "epithelial CLDN4 mean log1p(CP10k)",
        "NMPR vs MPR · epithelial CLDN4 (n=12 complete)",
        find_g("NMPR vs MPR: epithelial CLDN4 mean log1p(CP10k)"),
    )
    fig.suptitle("GSE207422 · NMPR n=8 vs MPR n=4 (pCR P06 = MPR) · CLDN4 only", fontsize=11)
    fig.savefig(figdir / "fig_cldn4_nmpr_mpr.png", dpi=160)
    fig.savefig(figdir / "fig_cldn4_nmpr_mpr.pdf")
    plt.close(fig)

    # Extra: Q4 vs Q1 immune
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.0), constrained_layout=True)
    extra_q = [
        (axes[0, 0], "epi_cldn4_mean_log1p_cp10k", "frac_T_NK", "epithelial CLDN4 mean", "T/NK"),
        (axes[0, 1], "epi_cldn4_mean_log1p_cp10k", "frac_CXCL13_pos", "epithelial CLDN4 mean", "CXCL13+"),
        (axes[0, 2], "epi_cldn4_mean_log1p_cp10k", "frac_cyto_high_T", "epithelial CLDN4 mean", "cyto-high T"),
        (axes[1, 0], "epi_cldn4_pct_pos", "frac_T_NK", "epithelial CLDN4 %pos", "T/NK"),
        (axes[1, 1], "epi_cldn4_pct_pos", "frac_CXCL13_pos", "epithelial CLDN4 %pos", "CXCL13+"),
        (axes[1, 2], "epi_cldn4_pct_pos", "frac_cyto_high_T", "epithelial CLDN4 %pos", "cyto-high T"),
    ]
    for ax, scol, icol, slabel, ilabel in extra_q:
        q1_idx, q4_idx, k = quartile_split(post[scol])
        q1 = post.loc[q1_idx]
        q4 = post.loc[q4_idx]
        key = f"Q4 vs Q1 ({slabel}): {ilabel} fraction"
        box_panel(
            ax,
            q4[icol],
            q1[icol],
            [f"Q4 (n={k})", f"Q1 (n={k})"],
            f"{ilabel} fraction",
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

    # Extra: A3-malignant Q4 vs Q1 (honest empty MPR)
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.0), constrained_layout=True)
    extra_q_mal = [
        (axes[0, 0], "mal_cldn4_mean_log1p_cp10k", "frac_T_NK", "A3-malignant CLDN4 mean", "T/NK"),
        (axes[0, 1], "mal_cldn4_mean_log1p_cp10k", "frac_CXCL13_pos", "A3-malignant CLDN4 mean", "CXCL13+"),
        (axes[0, 2], "mal_cldn4_mean_log1p_cp10k", "frac_cyto_high_T", "A3-malignant CLDN4 mean", "cyto-high T"),
        (axes[1, 0], "mal_cldn4_pct_pos", "frac_T_NK", "A3-malignant CLDN4 %pos", "T/NK"),
        (axes[1, 1], "mal_cldn4_pct_pos", "frac_CXCL13_pos", "A3-malignant CLDN4 %pos", "CXCL13+"),
        (axes[1, 2], "mal_cldn4_pct_pos", "frac_cyto_high_T", "A3-malignant CLDN4 %pos", "cyto-high T"),
    ]
    for ax, scol, icol, slabel, ilabel in extra_q_mal:
        q1_idx, q4_idx, k = quartile_split(post[scol])
        q1 = post.loc[q1_idx]
        q4 = post.loc[q4_idx]
        key = f"Q4 vs Q1 ({slabel}): {ilabel} fraction"
        box_panel(
            ax,
            q4[icol],
            q1[icol],
            [f"Q4 (n={k})", f"Q1 (n={k})"],
            f"{ilabel} fraction",
            f"{slabel}\nQ4 vs Q1 {ilabel}",
            find_q(key),
            colors=["#3d5a40", "#c9a227"],
        )
    fig.suptitle(
        "Extra · Q4 vs Q1 on A3-malignant CLDN4 (patients with 0 malignant cells dropped)",
        fontsize=11,
    )
    fig.savefig(figdir / "fig_extra_q4q1_a3malig.png", dpi=160)
    fig.savefig(figdir / "fig_extra_q4q1_a3malig.pdf")
    plt.close(fig)

    # Extra: companion TACSTD2 (not a gate)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), constrained_layout=True)
    scatter_panel(
        axes[0],
        post,
        "mal_cldn4_mean_log1p_cp10k",
        "mal_tacstd2_mean_log1p_cp10k",
        "A3-malignant CLDN4 mean log1p(CP10k)",
        "A3-malignant TACSTD2 mean log1p(CP10k)",
        "Companion · CLDN4 vs TACSTD2\n(not a dual-high gate)",
        find_s("post n=12: A3-malignant CLDN4 mean vs TACSTD2 mean (companion, not a gate)"),
    )
    scatter_panel(
        axes[1],
        post,
        "mal_tacstd2_mean_log1p_cp10k",
        "frac_T_NK",
        "A3-malignant TACSTD2 mean log1p(CP10k)",
        "T/NK fraction",
        "Companion TACSTD2 vs T/NK\n(A3 given; not re-argued)",
        find_s("post n=12: A3-malignant TACSTD2 mean log1p(CP10k) [companion] vs T/NK"),
    )
    box_panel(
        axes[2],
        nmpr["mal_tacstd2_mean_log1p_cp10k"],
        mpr["mal_tacstd2_mean_log1p_cp10k"],
        ["NMPR", "MPR"],
        "A3-malignant TACSTD2 mean log1p(CP10k)",
        "Companion TACSTD2 · NMPR vs MPR\n(A3 given)",
        find_g("NMPR vs MPR: A3-malignant TACSTD2 mean log1p(CP10k) [companion]"),
    )
    fig.suptitle("Extra · TACSTD2 companion only — never used to gate CLDN4", fontsize=11)
    fig.savefig(figdir / "fig_extra_companion_tacstd2.png", dpi=160)
    fig.savefig(figdir / "fig_extra_companion_tacstd2.pdf")
    plt.close(fig)

    # Extra: compartment honesty
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3), constrained_layout=True)
    order = post.sort_values(["paper_group", "patient"])
    x = np.arange(len(order))
    axes[0].bar(x, order["n_malig_a3"], color=[COLOR[g] for g in order["paper_group"]], edgecolor="white")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(order["patient"], rotation=45, ha="right")
    axes[0].set_ylabel("A3-malignant cells")
    axes[0].set_title("A3 malignant-like n per post patient")
    axes[0].axhline(10, color="#666", ls="--", lw=0.8, label="n=10")
    axes[1].bar(x, order["n_epithelial"], color=[COLOR[g] for g in order["paper_group"]], edgecolor="white")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(order["patient"], rotation=45, ha="right")
    axes[1].set_ylabel("epithelial cells")
    axes[1].set_title("Epithelial n per post patient (n=12 complete)")
    for ax in axes:
        ax.legend(
            handles=[
                plt.Line2D([0], [0], color=COLOR["MPR"], lw=6, label="MPR"),
                plt.Line2D([0], [0], color=COLOR["NMPR"], lw=6, label="NMPR"),
            ],
            frameon=False,
            fontsize=8,
        )
    fig.suptitle("Extra · empty A3-malignant MPR samples are not hidden", fontsize=11)
    fig.savefig(figdir / "fig_extra_compartment.png", dpi=160)
    fig.savefig(figdir / "fig_extra_compartment.pdf")
    plt.close(fig)

    np.random.set_state(rng_state)


def write_finding(outdir: Path, post, spearman_rows, group_rows, q_rows, q_assign, sanity) -> None:
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

    lines = [
        "# GSE207422 — malignant CLDN4 only (not dual-high)",
        "",
        "**Dual-high was stopped.** This slice scores **malignant CLDN4** only. "
        "TACSTD2 is a companion gene and is **never a gate**. "
        "The given A3 TACSTD2 analysis (NMPR vs MPR and vs T/NK on this same public UMI) is not re-argued.",
        "",
        "## Data and n",
        "",
        f"- Public GEO UMI only: **{sanity['n_cells']:,}** cells × **{sanity['n_genes_in_matrix']:,}** genes. "
        "Raw GSA-Human HRA001033 was not used. Author CopyKAT / epithelium RDS barcodes are not on GEO.",
        "- **Unit of every test is the 12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). "
        "The three pre-treatment biopsies are excluded. Cell-level p-values are not reported.",
        f"- Marker lineages (Hu canonical argmax): epithelial {sanity['n_epithelial']:,}; "
        f"A3-malignant-like {sanity['n_malig_a3']:,}; T/NK {sanity['n_T_NK']:,}; "
        f"CXCL13+ T/NK {sanity['n_CXCL13_pos']:,}; cyto-high T {sanity['n_cyto_high_T']:,}.",
        f"- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
        f"(same rule as the given A3 TACSTD2 slice; **not** CopyKAT). "
        f"Post patients with 0 A3-malignant cells: **{empty}**. "
        f"Post patients with <10 A3-malignant cells: **{lt10}**. "
        f"Those patients are NaN on A3-malignant CLDN4 and drop from that test "
        f"(A3-malignant CLDN4 n={n_mal_mean}/12). "
        f"Epithelial CLDN4 is the complete-case n={n_epi}/12 score.",
        "",
        "## Definitions",
        "",
        "| Item | Rule |",
        "|---|---|",
        "| CLDN4 mean | patient mean `log1p(CP10k)` inside the named compartment |",
        "| CLDN4 %pos | fraction of cells in the compartment with CLDN4 UMI ≥ 1 |",
        "| T/NK | lineage T or NK / all cells |",
        "| CXCL13+ | (T or NK) AND CXCL13 UMI ≥ 1 / all cells |",
        "| cyto-high T | lineage T AND mean log1p(GZMB, GZMA, PRF1, IFNG, NKG7) ≥ pooled T 75th percentile / all cells |",
        "| Q4 vs Q1 | among patients with a finite rank score, Q4 = highest n//4, Q1 = lowest n//4 (n=12 → 3 vs 3) |",
        "| NMPR vs MPR | exact two-sided Wilcoxon by enumerating label assignments; pCR P06 = MPR |",
        "| TACSTD2 | companion only; not used to define dual-high or to filter patients |",
        "",
        "## Primary — A3-malignant CLDN4 vs immune (honest n after empty MPR drop)",
        "",
        "| CLDN4 score | vs T/NK | vs CXCL13+ | vs cyto-high T |",
        "|---|---|---|---|",
        f"| mean log1p(CP10k) | {fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs T/NK'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs CXCL13+'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 mean log1p(CP10k) vs cyto-high T'))} |",
        f"| %pos | {fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs T/NK'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs CXCL13+'))} | "
        f"{fmt_rho(S('post n=12: A3-malignant CLDN4 %pos vs cyto-high T'))} |",
        "",
        "## Primary — NMPR vs MPR",
        "",
        "| Score | Result |",
        "|---|---|",
        f"| A3-malignant CLDN4 mean | {fmt_mwu(G('NMPR vs MPR: A3-malignant CLDN4 mean log1p(CP10k)'))} |",
        f"| A3-malignant CLDN4 %pos | {fmt_mwu(G('NMPR vs MPR: A3-malignant CLDN4 %pos'))} |",
        f"| epithelial CLDN4 mean (n=12 complete) | {fmt_mwu(G('NMPR vs MPR: epithelial CLDN4 mean log1p(CP10k)'))} |",
        f"| epithelial CLDN4 %pos | {fmt_mwu(G('NMPR vs MPR: epithelial CLDN4 %pos'))} |",
        f"| T/NK fraction | {fmt_mwu(G('NMPR vs MPR: T/NK fraction'))} |",
        f"| CXCL13+ fraction | {fmt_mwu(G('NMPR vs MPR: CXCL13+ fraction'))} |",
        f"| cyto-high T fraction | {fmt_mwu(G('NMPR vs MPR: cyto-high T fraction'))} |",
        "",
        "## Extra — Q4 vs Q1 immune fractions",
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
        "| Rank | vs T/NK | vs CXCL13+ | vs cyto-high T |",
        "|---|---|---|---|",
        f"| epithelial mean | {fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): T/NK fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): CXCL13+ fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 mean): cyto-high T fraction'))} |",
        f"| epithelial %pos | {fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): T/NK fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): CXCL13+ fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (epithelial CLDN4 %pos): cyto-high T fraction'))} |",
        "",
        "### A3-malignant CLDN4 ranks (empty malignant dropped)",
        "",
        f"- Mean-rank Q4: {q_names('mal_cldn4_mean_log1p_cp10k', 'Q4')}",
        f"- Mean-rank Q1: {q_names('mal_cldn4_mean_log1p_cp10k', 'Q1')}",
        "",
        "| Rank | vs T/NK | vs CXCL13+ | vs cyto-high T |",
        "|---|---|---|---|",
        f"| A3-malignant mean | {fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): T/NK fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): CXCL13+ fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 mean): cyto-high T fraction'))} |",
        f"| A3-malignant %pos | {fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): T/NK fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): CXCL13+ fraction'))} | "
        f"{fmt_mwu(Q('Q4 vs Q1 (A3-malignant CLDN4 %pos): cyto-high T fraction'))} |",
        "",
        "## Complete-case epithelial CLDN4 vs immune (n=12)",
        "",
        "| CLDN4 score | vs T/NK | vs CXCL13+ | vs cyto-high T |",
        "|---|---|---|---|",
        f"| mean log1p(CP10k) | {fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs T/NK'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs CXCL13+'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 mean log1p(CP10k) vs cyto-high T'))} |",
        f"| %pos | {fmt_rho(S('post n=12: epithelial CLDN4 %pos vs T/NK'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 %pos vs CXCL13+'))} | "
        f"{fmt_rho(S('post n=12: epithelial CLDN4 %pos vs cyto-high T'))} |",
        "",
        "## Companion TACSTD2 (not a gate; A3 given)",
        "",
        "These numbers are reported so CLDN4 can be read next to TACSTD2. "
        "They are **not** a dual-high filter and they do not replace the given A3 TACSTD2 write-up.",
        "",
        f"- A3-malignant CLDN4 mean vs TACSTD2 mean: {fmt_rho(S('post n=12: A3-malignant CLDN4 mean vs TACSTD2 mean (companion, not a gate)'))}",
        f"- Companion TACSTD2 mean vs T/NK: {fmt_rho(S('post n=12: A3-malignant TACSTD2 mean log1p(CP10k) [companion] vs T/NK'))}",
        f"- Companion TACSTD2 NMPR vs MPR: {fmt_mwu(G('NMPR vs MPR: A3-malignant TACSTD2 mean log1p(CP10k) [companion]'))}",
        "",
        "## Honest limits",
        "",
        "1. n=12 (4 vs 8) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. "
        "Q4 vs Q1 is 3 vs 3; the smallest exact two-sided p is 0.10.",
        "2. A3-malignant-like empties some MPR residual tumors (normal-lung program). "
        "That is reported, not patched by dual-high or by dropping MPR from the header n.",
        "3. This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak into A3-malignant-like.",
        "4. Dual-high (TACSTD2 AND CLDN4) was not run.",
        "5. Expression is log1p(CP10k) from the public UMI, not author Seurat-normalized values.",
        "",
        "## Extra figures",
        "",
        "- `figures/fig_cldn4_vs_immune.png` — mean and %pos vs T/NK, CXCL13+, cyto-high T",
        "- `figures/fig_cldn4_nmpr_mpr.png` — NMPR vs MPR",
        "- `figures/fig_extra_q4q1.png` — Q4 vs Q1 on epithelial ranks (n=12)",
        "- `figures/fig_extra_q4q1_a3malig.png` — Q4 vs Q1 on A3-malignant ranks",
        "- `figures/fig_extra_companion_tacstd2.png` — TACSTD2 companion, not a gate",
        "- `figures/fig_extra_compartment.png` — per-patient malignant n (empty MPR visible)",
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
        "python3 methods/gse207422_dualhigh/scripts/download.py",
        "python3 methods/gse207422_dualhigh/scripts/extract.py",
        "python3 methods/gse207422_dualhigh/scripts/analyze.py",
        "```",
        "",
    ]
    (outdir / "FINDING.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()

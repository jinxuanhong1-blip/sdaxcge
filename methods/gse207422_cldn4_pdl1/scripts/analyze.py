#!/usr/bin/env python3
"""GSE207422 malignant CLDN4 vs CD274/HLA in the same cells, plus vs T/NK.

CLDN4 only. TACSTD2 is never a gate. Dual-high is not run.
Author CopyKAT barcodes are not public. Public GEO UMI only.

Honest n:
  - Patient-level tests use the 12 post-treatment patients (empty A3-malignant
    patients are dropped, not patched).
  - Pooled same-cell Spearman is descriptive (cells are not independent).
  - Within-patient Spearman requires ≥20 A3-malignant cells.
  - Patient-pseudobulk Spearman is the honest same-cell inference unit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parents[1]

PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
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
MHC1_GENES = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
HLA_SINGLE = ["HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DRA"]
MIN_CELLS_WITHIN = 20
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


def cp10k(umi: np.ndarray, lib: np.ndarray) -> np.ndarray:
    umi = np.asarray(umi, dtype=np.float64)
    lib = np.asarray(lib, dtype=np.float64)
    return np.where(lib > 0, umi / lib * 1e4, np.nan)


def log1p_cp10k(umi: np.ndarray, lib: np.ndarray) -> np.ndarray:
    return np.log1p(cp10k(umi, lib))


def gene_stats(umi: np.ndarray, lib: np.ndarray) -> dict:
    umi = np.asarray(umi, dtype=np.float64)
    lib = np.asarray(lib, dtype=np.float64)
    if umi.size == 0:
        return {
            "n": 0,
            "mean_log1p_cp10k": np.nan,
            "pct_pos": np.nan,
        }
    return {
        "n": int(umi.size),
        "mean_log1p_cp10k": float(np.nanmean(log1p_cp10k(umi, lib))),
        "pct_pos": float(np.mean(umi >= 1) * 100.0),
    }


def spearman_row(x, y, contrast: str, n_cells: int | None = None, note: str = "") -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 4:
        return {
            "contrast": contrast,
            "n": n,
            "n_cells": n_cells,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": note or "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": contrast,
        "n": n,
        "n_cells": n_cells,
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": note,
    }


def fmt_rho(row: dict, show_cells: bool = False) -> str:
    if row.get("note") == "too_few_samples" or not np.isfinite(row.get("spearman_rho", np.nan)):
        extra = f", n_cells={row['n_cells']}" if show_cells and row.get("n_cells") else ""
        return f"n={row['n']}{extra} (too few)"
    extra = f", n_cells={row['n_cells']}" if show_cells and row.get("n_cells") else ""
    return f"ρ={row['spearman_rho']:.2f}, p={row['spearman_p']:.2g}, n={row['n']}{extra}"


def fmt_rho_nop(row: dict, show_cells: bool = False) -> str:
    """Same-cell pooled: report ρ and n, not p (cells are not independent)."""
    if row.get("note") == "too_few_samples" or not np.isfinite(row.get("spearman_rho", np.nan)):
        extra = f", n_cells={row['n_cells']}" if show_cells and row.get("n_cells") else ""
        return f"n={row['n']}{extra} (too few)"
    extra = f", n_cells={row['n_cells']}" if show_cells and row.get("n_cells") else ""
    return f"ρ={row['spearman_rho']:.2f}, n={row['n']}{extra} (p not used)"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()
    figdir = args.outdir / "figures"
    tabdir = args.outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    extracted = np.load(args.workdir / "extracted_cldn4_pdl1_markers.npz", allow_pickle=True)
    cells = extracted["cells"]
    genes = list(extracted["genes"])
    mat = extracted["mat"]
    total = extracted["total"]
    expr = {g: mat[i] for i, g in enumerate(genes)}
    n = len(cells)
    if "CLDN4" not in expr:
        raise SystemExit("CLDN4 missing from extracted matrix")

    mhc1_used = [g for g in MHC1_GENES if g in expr]
    if len(mhc1_used) < 2:
        raise SystemExit(f"MHC-I cassette too thin: present {mhc1_used}")

    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_lineage(expr, n)
    is_tnk = (lineage == "T") | (lineage == "NK")
    is_epi = lineage == "epithelial"
    a3_normal = np.zeros(n, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        if g in expr:
            a3_normal += expr[g]
    is_malig = is_epi & (a3_normal == 0)

    lib = total.astype(np.float64)
    cldn4_cp = log1p_cp10k(expr["CLDN4"], lib)
    cd274_cp = log1p_cp10k(expr["CD274"], lib) if "CD274" in expr else np.full(n, np.nan)
    hla_cp = {g: log1p_cp10k(expr[g], lib) for g in HLA_SINGLE if g in expr}
    mhc1_cp = np.nanmean(np.vstack([hla_cp[g] for g in mhc1_used]), axis=0)

    per_cell = pd.DataFrame(
        {
            "cell": cells,
            "Sample": sample_ids,
            "lineage": lineage,
            "is_tnk": is_tnk,
            "is_epi": is_epi,
            "is_malig_a3": is_malig,
            "total": total,
            "CLDN4": expr["CLDN4"],
            "CLDN4_log1p_cp10k": cldn4_cp,
            "CD274": expr["CD274"] if "CD274" in expr else 0,
            "CD274_log1p_cp10k": cd274_cp,
            "MHC1_log1p_cp10k": mhc1_cp,
        }
    )
    for g, arr in hla_cp.items():
        per_cell[g] = expr[g]
        per_cell[f"{g}_log1p_cp10k"] = arr
    per_cell["paper_group"] = per_cell["Sample"].map(PAPER_GROUP)
    per_cell["patient"] = per_cell["Sample"].str.replace("BD_immune", "P", regex=False)
    per_cell["timing"] = np.where(per_cell["paper_group"].isin(["MPR", "NMPR"]), "post", "pre")

    sample_meta = load_sample_meta(args.workdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    sample_meta.to_csv(args.outdir / "sample_metadata.tsv", sep="\t", index=False)

    targets = [("CD274", "CD274_log1p_cp10k")]
    for g in HLA_SINGLE:
        if f"{g}_log1p_cp10k" in per_cell.columns:
            targets.append((g, f"{g}_log1p_cp10k"))
    targets.append(("MHC-I mean (HLA-A/B/C + B2M)", "MHC1_log1p_cp10k"))

    rows = []
    for sample, sdf in per_cell.groupby("Sample", sort=True):
        mal = sdf[sdf["is_malig_a3"]]
        epi = sdf[sdf["is_epi"]]
        n_all = len(sdf)
        row = {
            "Sample": sample,
            "patient": sample.replace("BD_immune", "P"),
            "paper_group": PAPER_GROUP.get(sample),
            "n_cells": n_all,
            "n_epithelial": int(len(epi)),
            "n_malig_a3": int(len(mal)),
            "n_T_NK": int(sdf["is_tnk"].sum()),
            "frac_T_NK": float(sdf["is_tnk"].sum() / n_all) if n_all else np.nan,
        }
        for prefix, frame in (("mal", mal), ("epi", epi)):
            st = gene_stats(frame["CLDN4"].to_numpy(), frame["total"].to_numpy())
            row[f"{prefix}_cldn4_n"] = st["n"]
            row[f"{prefix}_cldn4_mean_log1p_cp10k"] = st["mean_log1p_cp10k"]
            row[f"{prefix}_cldn4_pct_pos"] = st["pct_pos"]
            if "CD274" in frame.columns:
                st274 = gene_stats(frame["CD274"].to_numpy(), frame["total"].to_numpy())
                row[f"{prefix}_cd274_mean_log1p_cp10k"] = st274["mean_log1p_cp10k"]
                row[f"{prefix}_cd274_pct_pos"] = st274["pct_pos"]
            if len(frame):
                row[f"{prefix}_mhc1_mean_log1p_cp10k"] = float(
                    np.nanmean(frame["MHC1_log1p_cp10k"].to_numpy())
                )
                for g in HLA_SINGLE:
                    col = f"{g}_log1p_cp10k"
                    if col in frame.columns:
                        row[f"{prefix}_{g}_mean_log1p_cp10k"] = float(np.nanmean(frame[col]))
            else:
                row[f"{prefix}_mhc1_mean_log1p_cp10k"] = np.nan
        rows.append(row)

    sample_df = pd.DataFrame(rows)
    sample_df.to_csv(tabdir / "per_patient.tsv", sep="\t", index=False)
    sample_df.to_csv(args.outdir / "per_patient.tsv", sep="\t", index=False)

    post = sample_df[sample_df["paper_group"].isin(["MPR", "NMPR"])].copy()
    if len(post) != 12:
        raise SystemExit(f"expected 12 post-treatment patients, got {len(post)}")

    post_cells = per_cell[per_cell["timing"] == "post"].copy()
    mal_post = post_cells[post_cells["is_malig_a3"]].copy()
    epi_post = post_cells[post_cells["is_epi"]].copy()

    # --- same-cell pooled (descriptive) ---
    pooled_rows = []
    for label, col in targets:
        pooled_rows.append(
            spearman_row(
                mal_post["CLDN4_log1p_cp10k"],
                mal_post[col],
                f"pooled A3-malignant cells: CLDN4 vs {label}",
                n_cells=int(len(mal_post)),
                note="descriptive_cells_not_independent",
            )
        )
        pooled_rows[-1]["n_patients"] = int(mal_post["Sample"].nunique())
        pooled_rows[-1]["compartment"] = "A3-malignant"
        pooled_rows[-1]["target"] = label
        pooled_rows[-1]["unit"] = "pooled_cells"
    for label, col in targets:
        pooled_rows.append(
            spearman_row(
                epi_post["CLDN4_log1p_cp10k"],
                epi_post[col],
                f"pooled epithelial cells: CLDN4 vs {label}",
                n_cells=int(len(epi_post)),
                note="descriptive_cells_not_independent",
            )
        )
        pooled_rows[-1]["n_patients"] = int(epi_post["Sample"].nunique())
        pooled_rows[-1]["compartment"] = "epithelial"
        pooled_rows[-1]["target"] = label
        pooled_rows[-1]["unit"] = "pooled_cells"
    # overwrite n with n_patients for honesty in the table header; keep n_cells
    for r in pooled_rows:
        r["n_patients_honest"] = r["n_patients"]
    pd.DataFrame(pooled_rows).to_csv(tabdir / "same_cell_pooled.tsv", sep="\t", index=False)

    # --- within-patient Spearman ---
    within_rows = []
    for sample, sdf in mal_post.groupby("Sample", sort=True):
        if len(sdf) < MIN_CELLS_WITHIN:
            continue
        for label, col in targets:
            row = spearman_row(
                sdf["CLDN4_log1p_cp10k"],
                sdf[col],
                f"within {sample}: CLDN4 vs {label}",
                n_cells=int(len(sdf)),
            )
            row["Sample"] = sample
            row["patient"] = sample.replace("BD_immune", "P")
            row["paper_group"] = PAPER_GROUP.get(sample)
            row["target"] = label
            row["compartment"] = "A3-malignant"
            row["unit"] = "within_patient"
            within_rows.append(row)
    within_df = pd.DataFrame(within_rows)
    within_df.to_csv(tabdir / "within_patient.tsv", sep="\t", index=False)

    within_summary = []
    for label, _col in targets:
        sub = within_df[within_df["target"] == label] if len(within_df) else pd.DataFrame()
        rhos = sub["spearman_rho"].to_numpy(dtype=float) if len(sub) else np.array([])
        rhos = rhos[np.isfinite(rhos)]
        n_pt = int(len(rhos))
        if n_pt >= 6:
            w_stat, w_p = stats.wilcoxon(rhos, alternative="two-sided", zero_method="wilcox")
        else:
            w_stat, w_p = np.nan, np.nan
        within_summary.append(
            {
                "contrast": f"within-patient A3-malignant: CLDN4 vs {label}",
                "target": label,
                "n": n_pt,
                "min_cells": MIN_CELLS_WITHIN,
                "median_rho": float(np.median(rhos)) if n_pt else np.nan,
                "q25_rho": float(np.quantile(rhos, 0.25)) if n_pt else np.nan,
                "q75_rho": float(np.quantile(rhos, 0.75)) if n_pt else np.nan,
                "n_negative": int(np.sum(rhos < 0)) if n_pt else 0,
                "n_positive": int(np.sum(rhos > 0)) if n_pt else 0,
                "wilcoxon_stat": float(w_stat) if np.isfinite(w_stat) else np.nan,
                "wilcoxon_p": float(w_p) if np.isfinite(w_p) else np.nan,
                "note": "" if n_pt >= 4 else "too_few_samples",
                "unit": "within_patient_median",
            }
        )
    pd.DataFrame(within_summary).to_csv(tabdir / "within_patient_summary.tsv", sep="\t", index=False)

    # --- patient-pseudobulk (honest same-cell inference) ---
    pb_rows = []
    for prefix, compartment in (("mal", "A3-malignant"), ("epi", "epithelial")):
        for label, colname in (
            [("CD274", f"{prefix}_cd274_mean_log1p_cp10k")]
            + [
                (g, f"{prefix}_{g}_mean_log1p_cp10k")
                for g in HLA_SINGLE
                if f"{prefix}_{g}_mean_log1p_cp10k" in post.columns
            ]
            + [("MHC-I mean (HLA-A/B/C + B2M)", f"{prefix}_mhc1_mean_log1p_cp10k")]
        ):
            row = spearman_row(
                post[f"{prefix}_cldn4_mean_log1p_cp10k"],
                post[colname],
                f"patient-pseudobulk {compartment}: CLDN4 vs {label}",
            )
            row["target"] = label
            row["compartment"] = compartment
            row["unit"] = "patient_pseudobulk"
            pb_rows.append(row)
    pd.DataFrame(pb_rows).to_csv(tabdir / "patient_pseudobulk.tsv", sep="\t", index=False)

    # --- patient-level vs T/NK ---
    tnk_rows = [
        spearman_row(
            post["mal_cldn4_mean_log1p_cp10k"],
            post["frac_T_NK"],
            "patient: A3-malignant CLDN4 mean vs T/NK",
        ),
        spearman_row(
            post["mal_cldn4_pct_pos"],
            post["frac_T_NK"],
            "patient: A3-malignant CLDN4 %pos vs T/NK",
        ),
        spearman_row(
            post["epi_cldn4_mean_log1p_cp10k"],
            post["frac_T_NK"],
            "patient: epithelial CLDN4 mean vs T/NK",
        ),
        spearman_row(
            post["epi_cldn4_pct_pos"],
            post["frac_T_NK"],
            "patient: epithelial CLDN4 %pos vs T/NK",
        ),
    ]
    for r in tnk_rows:
        r["unit"] = "patient"
    pd.DataFrame(tnk_rows).to_csv(tabdir / "vs_tnk.tsv", sep="\t", index=False)

    lineage_counts = (
        per_cell.groupby(["Sample", "lineage"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    lineage_counts.to_csv(args.outdir / "lineage_counts.tsv", sep="\t", index=False)

    empty = post.loc[post["n_malig_a3"] == 0, "patient"].tolist()
    lt10 = post.loc[(post["n_malig_a3"] > 0) & (post["n_malig_a3"] < 10), "patient"].tolist()
    within_pts = sorted(within_df["patient"].unique()) if len(within_df) else []

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
        "n_T_NK": int(is_tnk.sum()),
        "n_malig_a3_post": int(len(mal_post)),
        "n_patients_with_malig_a3_post": int(mal_post["Sample"].nunique()),
        "genes_present": genes,
        "genes_missing": [
            g
            for g in ["CLDN4", "CD274", *MHC1_GENES, "HLA-DRA", *A3_NORMAL_LUNG]
            if g not in expr
        ],
        "mhc1_genes_used": mhc1_used,
        "post_malig_a3_empty": empty,
        "post_malig_a3_lt10": lt10,
        "within_patient_min_cells": MIN_CELLS_WITHIN,
        "within_patient_ids": within_pts,
        "n_within_patient": len(within_pts),
        "cd274_pct_pos_malig_post": float(np.mean(mal_post["CD274"].to_numpy() >= 1) * 100.0)
        if "CD274" in mal_post.columns and len(mal_post)
        else np.nan,
        "cldn4_pct_pos_malig_post": float(np.mean(mal_post["CLDN4"].to_numpy() >= 1) * 100.0)
        if len(mal_post)
        else np.nan,
        "malignant_definition_a3": (
            "epithelial argmax AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
            "(same as given A3 TACSTD2 analysis; not CopyKAT)"
        ),
    }
    with (args.outdir / "sanity.json").open("w") as fh:
        json.dump(sanity, fh, indent=2)

    make_figures(post, mal_post, figdir, pooled_rows, pb_rows, tnk_rows, within_summary)
    write_finding(
        args.outdir,
        post,
        pooled_rows,
        pb_rows,
        tnk_rows,
        within_summary,
        sanity,
    )

    summary = {
        "dataset": "GSE207422",
        "task": "malignant CLDN4 vs CD274/HLA same-cell + patient vs T/NK; CLDN4 only",
        "dual_high": False,
        "tacstd2_gate": False,
        "n_post": 12,
        "n_cells": int(n),
        "n_malig_a3": int(is_malig.sum()),
        "n_malig_a3_post": int(len(mal_post)),
        "empty_a3_malignant_post": empty,
        "n_within_patient": len(within_pts),
        "primary_pb_mal_cldn4_vs_cd274": next(
            r
            for r in pb_rows
            if r["contrast"] == "patient-pseudobulk A3-malignant: CLDN4 vs CD274"
        ),
        "primary_pb_mal_cldn4_vs_mhc1": next(
            r
            for r in pb_rows
            if r["contrast"]
            == "patient-pseudobulk A3-malignant: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)"
        ),
        "primary_tnk_mal_cldn4_mean": next(
            r for r in tnk_rows if r["contrast"] == "patient: A3-malignant CLDN4 mean vs T/NK"
        ),
    }
    with (args.outdir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2, default=str))


def find_row(rows: list[dict], contrast: str) -> dict:
    for r in rows:
        if r["contrast"] == contrast:
            return r
    return {
        "contrast": contrast,
        "n": 0,
        "n_cells": 0,
        "n_patients": 0,
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "note": "missing_contrast",
    }


def scatter_panel(ax, df, x, y, xlabel, ylabel, title, rho_row=None, use_p=True):
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
        txt = fmt_rho(rho_row) if use_p else fmt_rho_nop(rho_row)
        ax.text(
            0.03,
            0.97,
            txt,
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, lw=0.4),
        )
    ax.legend(frameon=False, fontsize=8, loc="lower right")


def hexbin_panel(ax, x, y, xlabel, ylabel, title, rho_row=None):
    hb = ax.hexbin(x, y, gridsize=40, cmap="YlGnBu", mincnt=1, linewidths=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    cb = plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("cells", fontsize=8)
    if rho_row is not None:
        ax.text(
            0.03,
            0.97,
            fmt_rho_nop(rho_row, show_cells=True),
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, lw=0.4),
        )


def make_figures(post, mal_post, figdir, pooled_rows, pb_rows, tnk_rows, within_summary) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4), constrained_layout=True)
    hexbin_panel(
        axes[0],
        mal_post["CLDN4_log1p_cp10k"],
        mal_post["CD274_log1p_cp10k"],
        "A3-malignant CLDN4 log1p(CP10k)",
        "CD274 log1p(CP10k)",
        "Same cells · CLDN4 vs CD274 (descriptive)",
        find_row(pooled_rows, "pooled A3-malignant cells: CLDN4 vs CD274"),
    )
    hexbin_panel(
        axes[1],
        mal_post["CLDN4_log1p_cp10k"],
        mal_post["MHC1_log1p_cp10k"],
        "A3-malignant CLDN4 log1p(CP10k)",
        "MHC-I mean log1p(CP10k)",
        "Same cells · CLDN4 vs MHC-I (descriptive)",
        find_row(
            pooled_rows,
            "pooled A3-malignant cells: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)",
        ),
    )
    scatter_panel(
        axes[2],
        post,
        "mal_cldn4_mean_log1p_cp10k",
        "mal_cd274_mean_log1p_cp10k",
        "A3-malignant CLDN4 mean",
        "A3-malignant CD274 mean",
        "Patient-pseudobulk · CLDN4 vs CD274",
        find_row(pb_rows, "patient-pseudobulk A3-malignant: CLDN4 vs CD274"),
        use_p=True,
    )
    fig.suptitle(
        "GSE207422 · malignant CLDN4 vs CD274/HLA (CLDN4 only; cells not the inference unit)",
        fontsize=11,
    )
    fig.savefig(figdir / "fig_same_cell.png", dpi=160)
    fig.savefig(figdir / "fig_same_cell.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4), constrained_layout=True)
    scatter_panel(
        axes[0],
        post,
        "mal_cldn4_mean_log1p_cp10k",
        "frac_T_NK",
        "A3-malignant CLDN4 mean log1p(CP10k)",
        "T/NK fraction",
        "Patient · CLDN4 mean vs T/NK",
        find_row(tnk_rows, "patient: A3-malignant CLDN4 mean vs T/NK"),
    )
    scatter_panel(
        axes[1],
        post,
        "mal_cldn4_pct_pos",
        "frac_T_NK",
        "A3-malignant CLDN4 %pos",
        "T/NK fraction",
        "Patient · CLDN4 %pos vs T/NK",
        find_row(tnk_rows, "patient: A3-malignant CLDN4 %pos vs T/NK"),
    )
    scatter_panel(
        axes[2],
        post,
        "mal_cldn4_mean_log1p_cp10k",
        "mal_mhc1_mean_log1p_cp10k",
        "A3-malignant CLDN4 mean",
        "A3-malignant MHC-I mean",
        "Patient-pseudobulk · CLDN4 vs MHC-I",
        find_row(
            pb_rows,
            "patient-pseudobulk A3-malignant: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)",
        ),
    )
    fig.suptitle(
        "GSE207422 post-treatment n=12 · CLDN4 only vs T/NK and MHC-I pseudobulk",
        fontsize=11,
    )
    fig.savefig(figdir / "fig_vs_tnk.png", dpi=160)
    fig.savefig(figdir / "fig_vs_tnk.pdf")
    plt.close(fig)

    # within-patient rho strip
    if within_summary:
        fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
        labels = [r["target"] for r in within_summary]
        meds = [r["median_rho"] for r in within_summary]
        q25 = [r["q25_rho"] for r in within_summary]
        q75 = [r["q75_rho"] for r in within_summary]
        xs = np.arange(len(labels))
        yerr = np.vstack(
            [
                np.array(meds) - np.array(q25),
                np.array(q75) - np.array(meds),
            ]
        )
        ax.errorbar(xs, meds, yerr=yerr, fmt="o", color="#2c6eaf", capsize=4)
        ax.axhline(0, color="#888", lw=0.8)
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel("within-patient Spearman ρ")
        n_w = within_summary[0]["n"]
        ax.set_title(
            f"Within-patient A3-malignant CLDN4 vs targets (n={n_w} patients, ≥{MIN_CELLS_WITHIN} cells)"
        )
        fig.savefig(figdir / "fig_within_patient.png", dpi=160)
        fig.savefig(figdir / "fig_within_patient.pdf")
        plt.close(fig)


def fmt_within(row: dict) -> str:
    if row.get("note") == "too_few_samples" or not np.isfinite(row.get("median_rho", np.nan)):
        return f"n={row['n']} (too few)"
    ptxt = f", Wilcoxon p={row['wilcoxon_p']:.2g}" if np.isfinite(row.get("wilcoxon_p", np.nan)) else ""
    return (
        f"median ρ={row['median_rho']:.2f} "
        f"(IQR {row['q25_rho']:.2f}–{row['q75_rho']:.2f}); "
        f"{row['n_negative']}/{row['n']} negative; n={row['n']}{ptxt}"
    )


def write_finding(outdir, post, pooled_rows, pb_rows, tnk_rows, within_summary, sanity) -> None:
    def P(contrast: str) -> dict:
        return find_row(pooled_rows, contrast)

    def B(contrast: str) -> dict:
        return find_row(pb_rows, contrast)

    def T(contrast: str) -> dict:
        return find_row(tnk_rows, contrast)

    def W(target: str) -> dict:
        for r in within_summary:
            if r["target"] == target:
                return r
        return {
            "target": target,
            "n": 0,
            "median_rho": np.nan,
            "q25_rho": np.nan,
            "q75_rho": np.nan,
            "n_negative": 0,
            "wilcoxon_p": np.nan,
            "note": "missing_contrast",
        }

    empty = ", ".join(sanity["post_malig_a3_empty"]) or "none"
    lt10 = ", ".join(sanity["post_malig_a3_lt10"]) or "none"
    within_ids = ", ".join(sanity["within_patient_ids"]) or "none"
    n_mal = int(post["mal_cldn4_mean_log1p_cp10k"].notna().sum())
    n_epi = int(post["epi_cldn4_mean_log1p_cp10k"].notna().sum())
    missing = ", ".join(sanity["genes_missing"]) if sanity["genes_missing"] else "none"
    mhc1 = "/".join(sanity["mhc1_genes_used"])

    pb_cd274 = B("patient-pseudobulk A3-malignant: CLDN4 vs CD274")
    pb_mhc1 = B("patient-pseudobulk A3-malignant: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)")
    tnk_mean = T("patient: A3-malignant CLDN4 mean vs T/NK")
    tnk_epi = T("patient: epithelial CLDN4 mean vs T/NK")

    def sign_word(row, key="spearman_rho"):
        v = row.get(key, np.nan)
        if not np.isfinite(v):
            return "undefined"
        if v > 0.15:
            return "positive"
        if v < -0.15:
            return "negative"
        return "near-null"

    verdict = (
        f"Same-cell CLDN4 vs CD274/HLA is {sign_word(pb_cd274)} for CD274 "
        f"({fmt_rho(pb_cd274)}) and {sign_word(pb_mhc1)} for MHC-I "
        f"({fmt_rho(pb_mhc1)}) on A3-malignant **patient-pseudobulk**. "
        f"Pooled-cell ρ is descriptive only (n_cells={sanity['n_malig_a3_post']:,} "
        f"from {sanity['n_patients_with_malig_a3_post']} patients). "
        f"Patient-level CLDN4 vs T/NK is {sign_word(tnk_mean)} "
        f"({fmt_rho(tnk_mean)}; epithelial complete-case {fmt_rho(tnk_epi)})."
    )

    lines = [
        "# GSE207422 — malignant CLDN4 vs CD274/HLA (same cells) + T/NK",
        "",
        "**CLDN4 only.** TACSTD2 is never a gate. Dual-high was not run. "
        "The given A3 TACSTD2 analysis on this public UMI is not re-argued.",
        "",
        f"**Verdict (honest n):** {verdict}",
        "",
        "## Data and n",
        "",
        f"- Public GEO UMI only: **{sanity['n_cells']:,}** cells × **{sanity['n_genes_in_matrix']:,}** genes. "
        "Raw GSA-Human HRA001033 was not used. Author CopyKAT / epithelium RDS barcodes are not on GEO.",
        "- **Patient is the inference unit** for every test that reports a p-value "
        "(12 post-treatment patients: MPR n=4 including pCR P06; NMPR n=8). "
        "The three pre-treatment biopsies are excluded from tests.",
        f"- Marker lineages (Hu canonical argmax): epithelial {sanity['n_epithelial']:,}; "
        f"A3-malignant-like {sanity['n_malig_a3']:,} "
        f"({sanity['n_malig_a3_post']:,} in the 12 post samples); "
        f"T/NK {sanity['n_T_NK']:,}.",
        f"- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 "
        f"(same rule as the given A3 TACSTD2 slice; **not** CopyKAT). "
        f"Post patients with 0 A3-malignant cells (dropped from malignant tests): **{empty}**. "
        f"Post patients with 1–9 A3-malignant cells (kept, noisy): **{lt10}**. "
        f"A3-malignant n={n_mal}/12. Epithelial complete-case n={n_epi}/12.",
        f"- MHC-I cassette used: **{mhc1}**. Genes missing from the UMI: **{missing}**.",
        f"- CD274 is sparse in A3-malignant post cells: "
        f"{sanity['cd274_pct_pos_malig_post']:.1f}% UMI≥1 "
        f"(CLDN4 {sanity['cldn4_pct_pos_malig_post']:.1f}% UMI≥1).",
        "",
        "## Definitions",
        "",
        "| Item | Rule |",
        "|---|---|",
        "| CLDN4 / CD274 / HLA | `log1p(CP10k)` from public UMI |",
        "| MHC-I | mean of HLA-A, HLA-B, HLA-C, B2M `log1p(CP10k)` (genes present) |",
        "| same-cell pooled | Spearman inside A3-malignant (or epithelial) cells; **p not used** |",
        f"| within-patient | Spearman per post patient with ≥{MIN_CELLS_WITHIN} A3-malignant cells; median ρ |",
        "| patient-pseudobulk | Spearman of per-patient mean `log1p(CP10k)` (honest same-cell n) |",
        "| T/NK | lineage T or NK / all cells, per patient |",
        "| TACSTD2 | not used |",
        "",
        "## Primary table — same-cell CLDN4 vs CD274/HLA (honest n)",
        "",
        "Pooled-cell ρ treats thousands of cells as independent and is **not** the claim. "
        "Patient-pseudobulk and within-patient median are the honest rows.",
        "",
        "| contrast | unit | n | result |",
        "|---|---|---:|---|",
        f"| CLDN4 vs CD274 | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs CD274')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs CD274')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs CD274'), show_cells=True)} |",
        f"| CLDN4 vs MHC-I | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)'), show_cells=True)} |",
        f"| CLDN4 vs HLA-A | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs HLA-A')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs HLA-A')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs HLA-A'), show_cells=True)} |",
        f"| CLDN4 vs HLA-B | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs HLA-B')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs HLA-B')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs HLA-B'), show_cells=True)} |",
        f"| CLDN4 vs HLA-C | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs HLA-C')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs HLA-C')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs HLA-C'), show_cells=True)} |",
        f"| CLDN4 vs B2M | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs B2M')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs B2M')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs B2M'), show_cells=True)} |",
        f"| CLDN4 vs HLA-DRA | pooled A3-malignant cells (descriptive) | {P('pooled A3-malignant cells: CLDN4 vs HLA-DRA')['n_patients']} patients / {P('pooled A3-malignant cells: CLDN4 vs HLA-DRA')['n_cells']} cells | {fmt_rho_nop(P('pooled A3-malignant cells: CLDN4 vs HLA-DRA'), show_cells=True)} |",
        f"| CLDN4 vs CD274 | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('CD274')['n']} | {fmt_within(W('CD274'))} |",
        f"| CLDN4 vs MHC-I | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('MHC-I mean (HLA-A/B/C + B2M)')['n']} | {fmt_within(W('MHC-I mean (HLA-A/B/C + B2M)'))} |",
        f"| CLDN4 vs HLA-A | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('HLA-A')['n']} | {fmt_within(W('HLA-A'))} |",
        f"| CLDN4 vs HLA-B | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('HLA-B')['n']} | {fmt_within(W('HLA-B'))} |",
        f"| CLDN4 vs HLA-C | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('HLA-C')['n']} | {fmt_within(W('HLA-C'))} |",
        f"| CLDN4 vs B2M | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('B2M')['n']} | {fmt_within(W('B2M'))} |",
        f"| CLDN4 vs HLA-DRA | within-patient median (≥{MIN_CELLS_WITHIN} cells) | {W('HLA-DRA')['n']} | {fmt_within(W('HLA-DRA'))} |",
        f"| CLDN4 vs CD274 | patient-pseudobulk A3-malignant | {pb_cd274['n']} | {fmt_rho(pb_cd274)} |",
        f"| CLDN4 vs MHC-I | patient-pseudobulk A3-malignant | {pb_mhc1['n']} | {fmt_rho(pb_mhc1)} |",
        f"| CLDN4 vs HLA-A | patient-pseudobulk A3-malignant | {B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-A')['n']} | {fmt_rho(B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-A'))} |",
        f"| CLDN4 vs HLA-B | patient-pseudobulk A3-malignant | {B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-B')['n']} | {fmt_rho(B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-B'))} |",
        f"| CLDN4 vs HLA-C | patient-pseudobulk A3-malignant | {B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-C')['n']} | {fmt_rho(B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-C'))} |",
        f"| CLDN4 vs B2M | patient-pseudobulk A3-malignant | {B('patient-pseudobulk A3-malignant: CLDN4 vs B2M')['n']} | {fmt_rho(B('patient-pseudobulk A3-malignant: CLDN4 vs B2M'))} |",
        f"| CLDN4 vs HLA-DRA | patient-pseudobulk A3-malignant | {B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-DRA')['n']} | {fmt_rho(B('patient-pseudobulk A3-malignant: CLDN4 vs HLA-DRA'))} |",
        "",
        f"Within-patient patients (≥{MIN_CELLS_WITHIN} A3-malignant cells): **{within_ids}**.",
        "",
        "## Primary table — patient-level CLDN4 vs T/NK",
        "",
        "| CLDN4 score | vs T/NK | n |",
        "|---|---|---:|",
        f"| A3-malignant mean log1p(CP10k) | {fmt_rho(tnk_mean)} | {tnk_mean['n']} |",
        f"| A3-malignant %pos | {fmt_rho(T('patient: A3-malignant CLDN4 %pos vs T/NK'))} | {T('patient: A3-malignant CLDN4 %pos vs T/NK')['n']} |",
        f"| epithelial mean log1p(CP10k) (complete-case) | {fmt_rho(tnk_epi)} | {tnk_epi['n']} |",
        f"| epithelial %pos | {fmt_rho(T('patient: epithelial CLDN4 %pos vs T/NK'))} | {T('patient: epithelial CLDN4 %pos vs T/NK')['n']} |",
        "",
        "## Sensitivity — epithelial patient-pseudobulk (n=12 complete)",
        "",
        "| contrast | n | result |",
        "|---|---:|---|",
        f"| CLDN4 vs CD274 | {B('patient-pseudobulk epithelial: CLDN4 vs CD274')['n']} | {fmt_rho(B('patient-pseudobulk epithelial: CLDN4 vs CD274'))} |",
        f"| CLDN4 vs MHC-I | {B('patient-pseudobulk epithelial: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)')['n']} | {fmt_rho(B('patient-pseudobulk epithelial: CLDN4 vs MHC-I mean (HLA-A/B/C + B2M)'))} |",
        f"| CLDN4 vs HLA-DRA | {B('patient-pseudobulk epithelial: CLDN4 vs HLA-DRA')['n']} | {fmt_rho(B('patient-pseudobulk epithelial: CLDN4 vs HLA-DRA'))} |",
        "",
        "## Honest limits",
        "",
        "1. n=12 (4 vs 8) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. "
        f"A3-malignant tests drop empty residual tumors ({empty}); that n is not patched.",
        "2. Pooled same-cell p-values are not reported as inference. Cells from one patient are not independent.",
        f"3. Within-patient Spearman needs ≥{MIN_CELLS_WITHIN} A3-malignant cells. "
        f"Patients below that floor are excluded from the median-ρ row, not imputed.",
        "4. CD274 (PD-L1) RNA is sparse in scRNA. A near-null same-cell ρ can be dropout, not biology.",
        "5. This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak into A3-malignant-like.",
        "6. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.",
        "7. Expression is log1p(CP10k) from the public UMI, not author Seurat-normalized values, "
        "and not PD-L1 / MHC protein.",
        "",
        "## Figures",
        "",
        "- `figures/fig_same_cell.png` — pooled hexbin (descriptive) + patient-pseudobulk CLDN4 vs CD274",
        "- `figures/fig_vs_tnk.png` — patient CLDN4 vs T/NK and vs MHC-I pseudobulk",
        "- `figures/fig_within_patient.png` — within-patient median ρ strip",
        "",
        "## Files",
        "",
        "- `tables/per_patient.tsv` — 15 samples; tests use the 12 post rows",
        "- `tables/same_cell_pooled.tsv` / `tables/within_patient.tsv` / `tables/within_patient_summary.tsv`",
        "- `tables/patient_pseudobulk.tsv` / `tables/vs_tnk.tsv`",
        "- `sanity.json` / `summary.json` / `file_manifest.tsv`",
        "- Scripts: `scripts/download.py`, `scripts/extract.py`, `scripts/analyze.py`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/gse207422_cldn4_pdl1/scripts/download.py",
        "python3 methods/gse207422_cldn4_pdl1/scripts/extract.py",
        "python3 methods/gse207422_cldn4_pdl1/scripts/analyze.py",
        "```",
        "",
    ]
    (outdir / "FINDING.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()

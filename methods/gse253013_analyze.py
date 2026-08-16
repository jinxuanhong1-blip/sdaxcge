#!/usr/bin/env python3
"""GSE253013 extra scRNA figure: epithelial-restricted TACSTD2/CLDN4 vs T/NK.

Unit of analysis is the patient (n=9 treatment-naive LUAD), not 10x lanes and
not cells. There are no public MPR / R / response labels on this series.
GSE207422 (user A3) is not re-analyzed here.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MIN_EPI = 10
MIN_TNK = 20

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
}
NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
TARGETS = ["TACSTD2", "CLDN4"]


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


def spearman_block(x, y, label: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return {
            "contrast": label,
            "n": int(len(x)),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": label,
        "n": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def parse_patient(text: str) -> str:
    m = re.search(r"(MRC0*\d+)", str(text).replace(" ", ""))
    return m.group(1) if m else "unknown"


def parse_tissue(text: str) -> str:
    t = str(text).upper()
    if t in {"NAT", "ANT", "NORMAL"} or "ANT" in t or "ADJACENT" in t or t == "NAT":
        return "ANT"
    return "Tumor"


def load_inputs(extracted: Path) -> tuple[dict[str, np.ndarray], pd.DataFrame, dict]:
    expr_path = extracted / "gene_panel.npz"
    meta_path = extracted / "cell_metadata.tsv"
    gene_info = json.loads((extracted / "gene_index.json").read_text())
    if not expr_path.exists():
        raise SystemExit(
            f"missing {expr_path}; run methods/gse253013_extract.py and the "
            "CSC assembly step first"
        )
    packed = np.load(expr_path)
    expr = {k: packed[k] for k in packed.files}
    meta = pd.read_csv(meta_path, sep="\t")
    return expr, meta, gene_info


def summarize_compartment(df: pd.DataFrame, gene: str) -> dict:
    if df.empty or gene not in df:
        return {
            f"{gene}_n": 0,
            f"{gene}_mean_log1p": np.nan,
            f"{gene}_mean_log1p_cp10k": np.nan,
            f"{gene}_pct_pos": np.nan,
        }
    umi = df[gene].to_numpy(dtype=float)
    lib = df["total"].to_numpy(dtype=float) if "total" in df.columns else None
    out = {
        f"{gene}_n": int(len(df)),
        f"{gene}_mean_log1p": float(np.mean(np.log1p(umi))),
        f"{gene}_pct_pos": float(np.mean(umi > 0) * 100.0),
    }
    if lib is not None:
        cp = np.where(lib > 0, umi / lib * 1e4, np.nan)
        out[f"{gene}_mean_log1p_cp10k"] = float(np.nanmean(np.log1p(cp)))
    else:
        out[f"{gene}_mean_log1p_cp10k"] = np.nan
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extracted", type=Path, default=Path("data/gse253013/extracted"))
    ap.add_argument("--outdir", type=Path, default=Path("results/gse253013"))
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    audit = {
        "dataset": "GSE253013",
        "pmid": "38335304",
        "n_patients_paper": 9,
        "n_gsm_lanes": 89,
        "treatment": "treatment-naive LUAD tumor + adjacent non-tumor lung",
        "is_blood_only": False,
        "is_ici_response_cohort": False,
        "public_mpr_r_labels": False,
        "gse207422_reanalyzed": False,
        "label_search": {
            "geo_series_matrix_response_fields": [],
            "geo_sample_characteristics": ["group: Tumor | Adjacent Non-Tumor", "tissue: lung"],
            "note": (
                "Sze/Xiang Cancer Research 2024 (PMID 38335304) sequenced 9 "
                "treatment-naive patients. No MPR, pCR, R/NR, DCB, or ICI "
                "outcome is attached to GSE253013."
            ),
        },
    }

    gene_index_path = args.extracted / "gene_index.json"
    if not gene_index_path.exists():
        (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
        (out / "STOP.md").write_text(
            "# GSE253013 extraction incomplete\n\n"
            "Run `python3 methods/gse253013_extract.py` first.\n"
        )
        print("extraction not finished")
        return 2

    gene_info = json.loads(gene_index_path.read_text())
    present = set(gene_info.get("present") or gene_info.get("index", {}).keys())
    audit["genes_present"] = sorted(present)
    audit["TACSTD2_present"] = "TACSTD2" in present
    audit["CLDN4_present"] = "CLDN4" in present

    if "TACSTD2" not in present and "CLDN4" not in present:
        audit["stop_reason"] = "both_target_genes_absent"
        (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
        (out / "STOP.md").write_text(
            "# Stopped: TACSTD2 and CLDN4 both absent\n\n"
            "After a real gene-name check of the GSE253013 processed object, "
            "neither target was present. No figure invented.\n"
        )
        print("STOP: both genes absent")
        return 0

    expr_path = args.extracted / "gene_panel.npz"
    meta_path = args.extracted / "cell_metadata.tsv"
    if not expr_path.exists() or not meta_path.exists():
        audit["stop_reason"] = "matrix_not_assembled"
        (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
        print("gene index exists but CSC assembly is incomplete")
        return 2

    packed = np.load(expr_path)
    expr = {k: packed[k] for k in packed.files}
    meta = pd.read_csv(meta_path, sep="\t")
    n = len(meta)
    audit["n_cells"] = int(n)

    for g in list(LINEAGES["epithelial"]) + list(LINEAGES["T"]) + list(LINEAGES["NK"]) + NORMAL_LUNG + TARGETS:
        if g in expr and g not in meta.columns:
            meta[g] = expr[g]

    if "nCount_RNA" in meta.columns:
        meta["total"] = meta["nCount_RNA"]
        audit["library_size"] = "Seurat nCount_RNA"
    elif "total" not in meta.columns:
        cols = [g for g in expr if g != "total"]
        meta["total"] = np.sum([expr[g] for g in cols], axis=0)
        audit["library_size"] = "sum_of_extracted_panel_only"

    if "patient" not in meta.columns:
        src = meta["orig.ident"] if "orig.ident" in meta.columns else meta.iloc[:, 0]
        meta["patient"] = [parse_patient(x) for x in src]
    if "tissue" not in meta.columns:
        src = meta["orig.ident"] if "orig.ident" in meta.columns else meta.iloc[:, 0]
        meta["tissue"] = [parse_tissue(x) for x in src]

    author_type_col = None
    for cand in (
        "garnett_cluster",
        "garnett_cell_type",
        "cell_type",
        "celltype",
        "CellType",
        "author_cell_type",
        "seurat_clusters",
        "ident",
    ):
        if cand in meta.columns:
            author_type_col = cand
            break
    audit["author_cell_type_column"] = author_type_col

    lineage = assign_lineage(expr, n)
    meta["lineage"] = lineage
    normal_score = score(expr, [g for g in NORMAL_LUNG if g in expr], n)
    meta["malignant_like"] = (meta["lineage"] == "epithelial") & (normal_score <= 0.05)
    meta["tnk"] = meta["lineage"].isin(["T", "NK"])

    if author_type_col:
        lab = meta[author_type_col].astype(str).str.lower()
        epi_author = lab.isin(["epithelial"]) | lab.str.fullmatch("epithelial")
        tnk_author = lab.isin(["t cells", "t cell", "nk", "nk cells"]) | lab.str.contains(
            r"^t cells$|^nk"
        )
        meta["author_epithelial"] = epi_author
        meta["author_tnk"] = tnk_author
    else:
        meta["author_epithelial"] = False
        meta["author_tnk"] = False

    rows = []
    for tissue in ("Tumor", "ANT"):
        sub_t = meta[meta["tissue"] == tissue]
        for patient, pdf in sub_t.groupby("patient"):
            n_cells = len(pdf)
            n_epi = int((pdf["lineage"] == "epithelial").sum())
            n_mal = int(pdf["malignant_like"].sum())
            n_tnk = int(pdf["tnk"].sum())
            rec = {
                "patient": patient,
                "tissue": tissue,
                "n_cells": n_cells,
                "n_epithelial": n_epi,
                "n_malignant_like": n_mal,
                "n_tnk": n_tnk,
                "tnk_fraction": n_tnk / n_cells if n_cells else np.nan,
                "eligible_malig": n_mal >= MIN_EPI and n_tnk >= MIN_TNK,
                "eligible_epi": n_epi >= MIN_EPI and n_tnk >= MIN_TNK,
            }
            for gene in TARGETS:
                if gene not in pdf:
                    continue
                rec.update(summarize_compartment(pdf[pdf["malignant_like"]], gene))
                epi = summarize_compartment(pdf[pdf["lineage"] == "epithelial"], gene)
                rec[f"epi_{gene}_mean_log1p"] = epi[f"{gene}_mean_log1p"]
                rec[f"epi_{gene}_mean_log1p_cp10k"] = epi[f"{gene}_mean_log1p_cp10k"]
                rec[f"epi_{gene}_pct_pos"] = epi[f"{gene}_pct_pos"]
                tnk = summarize_compartment(pdf[pdf["tnk"]], gene)
                rec[f"tnk_{gene}_mean_log1p"] = tnk[f"{gene}_mean_log1p"]
                rec[f"tnk_{gene}_pct_pos"] = tnk[f"{gene}_pct_pos"]
            rows.append(rec)

    per = pd.DataFrame(rows)
    per.to_csv(out / "per_patient_metrics.tsv", sep="\t", index=False)

    tests = []
    tumor = per[per["tissue"] == "Tumor"].copy()
    for gene in TARGETS:
        if gene not in present:
            tests.append(
                {
                    "contrast": f"{gene} absent from processed object",
                    "n": 0,
                    "spearman_rho": np.nan,
                    "spearman_p": np.nan,
                    "note": "gene_absent",
                }
            )
            continue
        for score_col, frac_col, label in (
            (f"{gene}_mean_log1p", "tnk_fraction", f"tumor malig-like {gene} mean_log1p vs T/NK fraction"),
            (f"{gene}_mean_log1p_cp10k", "tnk_fraction", f"tumor malig-like {gene} mean log1p(CP10k) vs T/NK fraction"),
            (f"{gene}_pct_pos", "tnk_fraction", f"tumor malig-like {gene} %pos vs T/NK fraction"),
            (f"epi_{gene}_mean_log1p", "tnk_fraction", f"tumor epithelial {gene} mean_log1p vs T/NK fraction"),
        ):
            elig = tumor[tumor["eligible_malig"]] if "malig" in label else tumor[tumor["eligible_epi"]]
            tests.append(spearman_block(elig[score_col], elig[frac_col], label + " (eligible)"))
            tests.append(spearman_block(tumor[score_col], tumor[frac_col], label + " (all tumor patients)"))

    tests.append(
        {
            "contrast": "response / MPR / R vs TACSTD2 or CLDN4",
            "n": 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "no_public_response_labels",
        }
    )

    # Sensitivity: drop patients with very small malignant-like compartments
    robust = tumor[tumor["n_malignant_like"] >= 50]
    for gene in TARGETS:
        if gene not in present:
            continue
        tests.append(
            spearman_block(
                robust[f"{gene}_mean_log1p"],
                robust["tnk_fraction"],
                f"tumor malig-like {gene} mean_log1p vs T/NK (n_malig>=50)",
            )
        )

    # Author Seurat labels: Epithelial vs T cells (no NK class in this object)
    if "cell_type" in meta.columns:
        arows = []
        for patient, pdf in meta[meta["tissue"] == "Tumor"].groupby("patient"):
            epi = pdf[pdf["cell_type"] == "Epithelial"]
            tnk = pdf[pdf["cell_type"] == "T cells"]
            rec = {
                "patient": patient,
                "n_author_epithelial": int(len(epi)),
                "n_author_tcells": int(len(tnk)),
                "author_t_fraction": len(tnk) / len(pdf) if len(pdf) else np.nan,
            }
            for gene in TARGETS:
                if gene in epi:
                    rec[f"author_epi_{gene}_mean_log1p"] = float(np.mean(np.log1p(epi[gene]))) if len(epi) else np.nan
            arows.append(rec)
        adf = pd.DataFrame(arows)
        adf.to_csv(out / "author_label_per_patient.tsv", sep="\t", index=False)
        for gene in TARGETS:
            if f"author_epi_{gene}_mean_log1p" not in adf:
                continue
            ok = adf[adf["n_author_epithelial"] >= MIN_EPI]
            tests.append(
                spearman_block(
                    ok[f"author_epi_{gene}_mean_log1p"],
                    ok["author_t_fraction"],
                    f"tumor author-Epithelial {gene} vs author T-cell fraction",
                )
            )
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(out / "association_statistics.tsv", sep="\t", index=False)

    # Restriction sanity: epithelial vs T/NK %pos, paired patients
    restriction = []
    for gene in TARGETS:
        if gene not in present:
            continue
        a = tumor[f"{gene}_pct_pos"].to_numpy(dtype=float)
        b = tumor[f"tnk_{gene}_pct_pos"].to_numpy(dtype=float)
        mask = np.isfinite(a) & np.isfinite(b)
        if mask.sum() >= 3:
            try:
                w = stats.wilcoxon(a[mask], b[mask], alternative="greater")
                p = float(w.pvalue)
            except ValueError:
                p = np.nan
            restriction.append(
                {
                    "gene": gene,
                    "n": int(mask.sum()),
                    "median_pct_pos_malig_like": float(np.median(a[mask])),
                    "median_pct_pos_tnk": float(np.median(b[mask])),
                    "wilcoxon_greater_p": p,
                }
            )
    pd.DataFrame(restriction).to_csv(out / "compartment_restriction.tsv", sep="\t", index=False)

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), constrained_layout=True)
    for ax, gene in zip(axes, TARGETS):
        if gene not in present:
            ax.set_axis_off()
            ax.set_title(f"{gene} absent")
            continue
        elig = tumor[tumor["eligible_malig"]]
        x = elig["tnk_fraction"].to_numpy(dtype=float)
        y = elig[f"{gene}_mean_log1p"].to_numpy(dtype=float)
        ax.scatter(x, y, c="#1f4e79", s=42, zorder=3)
        for _, r in elig.iterrows():
            ax.annotate(str(r["patient"]), (r["tnk_fraction"], r[f"{gene}_mean_log1p"]), fontsize=7, xytext=(4, 3), textcoords="offset points")
        blk = spearman_block(x, y, gene)
        ax.set_xlabel("T/NK fraction (tumor cells)")
        ax.set_ylabel(f"Malignant-like {gene}\nmean log1p UMI")
        ax.set_title(f"{gene}  n={blk['n']}  ρ={blk['spearman_rho']:.2f}  p={blk['spearman_p']:.2g}" if blk["n"] >= 4 else f"{gene}  n={blk['n']}  (too few)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("GSE253013 treatment-naïve LUAD (not ICI/MPR) · extra scRNA", fontsize=10)
    fig.savefig(out / "fig_gse253013_tacstd2_cldn4_vs_tnk.png", dpi=180)
    fig.savefig(out / "fig_gse253013_tacstd2_cldn4_vs_tnk.pdf")
    plt.close(fig)

    primary = [t for t in tests if "malig-like TACSTD2 mean_log1p vs T/NK" in t["contrast"] and "eligible" in t["contrast"]]
    primary_c = [t for t in tests if "malig-like CLDN4 mean_log1p vs T/NK" in t["contrast"] and "eligible" in t["contrast"]]
    summary = {
        **audit,
        "n_tumor_patients": int(tumor["patient"].nunique()) if len(tumor) else 0,
        "n_eligible_malig_tumor": int(tumor["eligible_malig"].sum()) if len(tumor) else 0,
        "primary_TACSTD2_vs_tnk": primary[0] if primary else None,
        "primary_CLDN4_vs_tnk": primary_c[0] if primary_c else None,
        "response_test": "not_possible_no_public_labels",
        "honest_unit": "patient",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

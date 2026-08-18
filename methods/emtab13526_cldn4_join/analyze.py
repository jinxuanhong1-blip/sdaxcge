#!/usr/bin/env python3
"""CLDN4-only patient-level join test for E-MTAB-13526.

Question: can this atlas JOIN the concordant pool
(GSE123902+GSE131907+GSE205335+GSE189357)?

Join only if sign matches:
  - T/NK inverse (Spearman ρ < 0 for malignant CLDN4 vs T/NK fraction), and/or
  - IFN/MHC down in CLDN4-high (Q4 vs Q1, or continuous ρ < 0).

No dual-high TACSTD2×CLDN4. Unit = patient. Public CD235a- tumor mtx only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MIN_EPI = 10
MIN_TNK = 20
MIN_Q_TAIL = 3

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


def histolabel(disease: str) -> str:
    d = str(disease).lower()
    if "adenocarcinoma" in d:
        return "LUAD"
    if "squamous" in d:
        return "LUSC"
    return "NSCLC_other"


def q4q1(df: pd.DataFrame, score_col: str, y_col: str, label: str) -> dict:
    work = df[[score_col, y_col]].dropna().copy()
    n = len(work)
    if n < 8:
        return {
            "contrast": label,
            "n": n,
            "n_Q1": 0,
            "n_Q4": 0,
            "mean_Q1": np.nan,
            "mean_Q4": np.nan,
            "delta_Q4_minus_Q1": np.nan,
            "mannwhitney_p": np.nan,
            "note": "n_too_small_for_quartiles",
        }
    try:
        work["q"] = pd.qcut(work[score_col], 4, labels=False, duplicates="drop")
    except ValueError:
        return {
            "contrast": label,
            "n": n,
            "n_Q1": 0,
            "n_Q4": 0,
            "mean_Q1": np.nan,
            "mean_Q4": np.nan,
            "delta_Q4_minus_Q1": np.nan,
            "mannwhitney_p": np.nan,
            "note": "qcut_failed_duplicate_edges",
        }
    qmax = int(work["q"].max())
    low = work[work["q"] == 0]
    high = work[work["q"] == qmax]
    n1, n4 = int(len(low)), int(len(high))
    note = ""
    p = np.nan
    if n1 < MIN_Q_TAIL or n4 < MIN_Q_TAIL:
        note = f"thin_tails_nQ1={n1}_nQ4={n4}"
    elif n1 >= 3 and n4 >= 3:
        try:
            p = float(stats.mannwhitneyu(high[y_col], low[y_col], alternative="two-sided").pvalue)
        except ValueError:
            p = np.nan
            note = "mannwhitney_failed"
    return {
        "contrast": label,
        "n": n,
        "n_Q1": n1,
        "n_Q4": n4,
        "mean_Q1": float(low[y_col].mean()) if n1 else np.nan,
        "mean_Q4": float(high[y_col].mean()) if n4 else np.nan,
        "delta_Q4_minus_Q1": (
            float(high[y_col].mean() - low[y_col].mean()) if n1 and n4 else np.nan
        ),
        "mannwhitney_p": p,
        "note": note,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extracted", type=Path, default=Path("data/emtab13526/extracted"))
    ap.add_argument("--sets", type=Path, default=Path("methods/emtab13526_cldn4_join/data/ifn_mhc_sets.json"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/emtab13526_cldn4_join"))
    args = ap.parse_args()
    out = args.outdir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    gene_index_path = args.extracted / "gene_index.json"
    if not gene_index_path.exists():
        (out / "tables" / "STOP.md").write_text(
            "# Extraction incomplete\n\nRun extract.py first.\n"
        )
        print("extraction not finished")
        return 2

    gene_info = json.loads(gene_index_path.read_text())
    present = set(gene_info.get("present") or [])
    sets = json.loads(args.sets.read_text())
    ifn_genes = [g for g in sets["IFN"] if g in present]
    mhc_genes = [g for g in sets["MHC_I_APM"] if g in present]

    packed = np.load(args.extracted / "gene_panel.npz")
    expr = {k: packed[k] for k in packed.files}
    meta = pd.read_csv(args.extracted / "cell_metadata.tsv", sep="\t")
    n = len(meta)

    lineage = assign_lineage(expr, n)
    meta["lineage"] = lineage
    normal_score = score(expr, [g for g in NORMAL_LUNG if g in expr], n)
    meta["malignant_like"] = (meta["lineage"] == "epithelial") & (normal_score <= 0.05)
    meta["tnk"] = meta["lineage"].isin(["T", "NK"])
    meta["histology"] = meta["disease"].map(histolabel) if "disease" in meta.columns else "NSCLC_other"
    if "CLDN4" in expr:
        meta["CLDN4"] = expr["CLDN4"]
    meta["ifn_log1p"] = score(expr, ifn_genes, n)
    meta["mhc_log1p"] = score(expr, mhc_genes, n)
    meta["ifn_mhc_log1p"] = score(expr, sorted(set(ifn_genes) | set(mhc_genes)), n)

    rows = []
    for patient, pdf in meta.groupby("patient"):
        n_cells = len(pdf)
        n_epi = int((pdf["lineage"] == "epithelial").sum())
        n_mal = int(pdf["malignant_like"].sum())
        n_tnk = int(pdf["tnk"].sum())
        mal = pdf[pdf["malignant_like"]]
        rec = {
            "patient": patient,
            "tissue": "Tumor",
            "histology": pdf["histology"].iloc[0],
            "disease": pdf["disease"].iloc[0] if "disease" in pdf else "",
            "n_lanes": int(pdf["sample"].nunique()),
            "n_cells": n_cells,
            "n_epithelial": n_epi,
            "n_malignant_like": n_mal,
            "n_tnk": n_tnk,
            "n_T": int((pdf["lineage"] == "T").sum()),
            "n_NK": int((pdf["lineage"] == "NK").sum()),
            "tnk_fraction": n_tnk / n_cells if n_cells else np.nan,
            "eligible_malig": bool(n_mal >= MIN_EPI and n_tnk >= MIN_TNK),
            "eligible_epi": bool(n_epi >= MIN_EPI and n_tnk >= MIN_TNK),
        }
        if n_mal:
            rec["CLDN4_mean_log1p"] = float(np.mean(np.log1p(mal["CLDN4"]))) if "CLDN4" in mal else np.nan
            rec["CLDN4_pct_pos"] = float(np.mean(mal["CLDN4"] > 0) * 100.0) if "CLDN4" in mal else np.nan
            rec["IFN_mean_log1p"] = float(mal["ifn_log1p"].mean())
            rec["MHC_mean_log1p"] = float(mal["mhc_log1p"].mean())
            rec["IFN_MHC_mean_log1p"] = float(mal["ifn_mhc_log1p"].mean())
        else:
            rec["CLDN4_mean_log1p"] = np.nan
            rec["CLDN4_pct_pos"] = np.nan
            rec["IFN_mean_log1p"] = np.nan
            rec["MHC_mean_log1p"] = np.nan
            rec["IFN_MHC_mean_log1p"] = np.nan
        rows.append(rec)

    per = pd.DataFrame(rows).sort_values("patient")
    per.to_csv(out / "tables" / "per_patient_metrics.tsv", sep="\t", index=False)

    elig = per[per["eligible_malig"]].copy()
    tests = []
    tests.append(
        spearman_block(
            elig["CLDN4_mean_log1p"],
            elig["tnk_fraction"],
            "tumor malig-like CLDN4 mean_log1p vs T/NK fraction (eligible)",
        )
    )
    tests.append(
        spearman_block(
            elig["CLDN4_pct_pos"],
            elig["tnk_fraction"],
            "tumor malig-like CLDN4 %pos vs T/NK fraction (eligible)",
        )
    )
    tests.append(
        spearman_block(
            per["CLDN4_mean_log1p"],
            per["tnk_fraction"],
            "tumor malig-like CLDN4 mean_log1p vs T/NK (all CD235a- tumor patients)",
        )
    )
    tests.append(
        spearman_block(
            elig["CLDN4_mean_log1p"],
            elig["IFN_mean_log1p"],
            "tumor malig-like CLDN4 vs IFN score (eligible, continuous)",
        )
    )
    tests.append(
        spearman_block(
            elig["CLDN4_mean_log1p"],
            elig["MHC_mean_log1p"],
            "tumor malig-like CLDN4 vs MHC-I/APM score (eligible, continuous)",
        )
    )
    tests.append(
        spearman_block(
            elig["CLDN4_mean_log1p"],
            elig["IFN_MHC_mean_log1p"],
            "tumor malig-like CLDN4 vs IFN/MHC score (eligible, continuous)",
        )
    )

    q_rows = []
    for ycol, ylab in (
        ("IFN_mean_log1p", "IFN"),
        ("MHC_mean_log1p", "MHC-I/APM"),
        ("IFN_MHC_mean_log1p", "IFN/MHC"),
        ("tnk_fraction", "T/NK fraction"),
    ):
        q_rows.append(q4q1(elig, "CLDN4_mean_log1p", ycol, f"CLDN4 mean_log1p Q4 vs Q1 · {ylab}"))
        q_rows.append(q4q1(elig, "CLDN4_pct_pos", ycol, f"CLDN4 %pos Q4 vs Q1 · {ylab}"))

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(out / "tables" / "association_statistics.tsv", sep="\t", index=False)
    q_df = pd.DataFrame(q_rows)
    q_df.to_csv(out / "tables" / "q4q1.tsv", sep="\t", index=False)

    lineage_counts = (
        meta.groupby(["patient", "lineage"]).size().unstack(fill_value=0).reset_index()
    )
    lineage_counts.to_csv(out / "tables" / "lineage_counts.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=True)
    ax = axes[0]
    ax.scatter(elig["tnk_fraction"], elig["CLDN4_mean_log1p"], c="#1f4e79", s=48, zorder=3)
    for _, r in elig.iterrows():
        ax.annotate(str(r["patient"]), (r["tnk_fraction"], r["CLDN4_mean_log1p"]), fontsize=7, xytext=(4, 3), textcoords="offset points")
    tnk = tests[0]
    ax.set_xlabel("T/NK fraction (CD235a− tumor cells)")
    ax.set_ylabel("Malignant-like CLDN4\nmean log1p UMI")
    if tnk["n"] >= 4:
        ax.set_title(f"CLDN4 vs T/NK  n={tnk['n']}  ρ={tnk['spearman_rho']:.2f}  p={tnk['spearman_p']:.2g}")
    else:
        ax.set_title(f"CLDN4 vs T/NK  n={tnk['n']}  (too few)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    ax.scatter(elig["CLDN4_mean_log1p"], elig["IFN_MHC_mean_log1p"], c="#8b1e3f", s=48, zorder=3)
    for _, r in elig.iterrows():
        ax.annotate(str(r["patient"]), (r["CLDN4_mean_log1p"], r["IFN_MHC_mean_log1p"]), fontsize=7, xytext=(4, 3), textcoords="offset points")
    ifnm = [t for t in tests if "IFN/MHC score" in t["contrast"]][0]
    ax.set_xlabel("Malignant-like CLDN4 mean log1p")
    ax.set_ylabel("Malignant-like IFN/MHC\nmean log1p")
    if ifnm["n"] >= 4:
        ax.set_title(f"CLDN4 vs IFN/MHC  n={ifnm['n']}  ρ={ifnm['spearman_rho']:.2f}  p={ifnm['spearman_p']:.2g}")
    else:
        ax.set_title(f"CLDN4 vs IFN/MHC  n={ifnm['n']}  (too few)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle("E-MTAB-13526 De Zuani NSCLC · CLDN4-only join test", fontsize=10)
    fig.savefig(out / "figures" / "fig_cldn4_tnk_ifnmhc.png", dpi=180)
    fig.savefig(out / "figures" / "fig_cldn4_tnk_ifnmhc.pdf")
    plt.close(fig)

    tnk_rho = tnk.get("spearman_rho")
    ifnm_rho = ifnm.get("spearman_rho")
    q_ifnm = next(r for r in q_rows if r["contrast"] == "CLDN4 mean_log1p Q4 vs Q1 · IFN/MHC")
    q_ifn = next(r for r in q_rows if r["contrast"] == "CLDN4 mean_log1p Q4 vs Q1 · IFN")
    q_mhc = next(r for r in q_rows if r["contrast"] == "CLDN4 mean_log1p Q4 vs Q1 · MHC-I/APM")

    tnk_inverse = bool(np.isfinite(tnk_rho) and tnk_rho < 0)
    ifnm_down_cont = bool(np.isfinite(ifnm_rho) and ifnm_rho < 0)
    ifnm_down_q = bool(
        np.isfinite(q_ifnm.get("delta_Q4_minus_Q1")) and q_ifnm["delta_Q4_minus_Q1"] < 0
        and q_ifnm["n_Q1"] >= MIN_Q_TAIL
        and q_ifnm["n_Q4"] >= MIN_Q_TAIL
    )
    ifn_down_q = bool(
        np.isfinite(q_ifn.get("delta_Q4_minus_Q1")) and q_ifn["delta_Q4_minus_Q1"] < 0
        and q_ifn["n_Q1"] >= MIN_Q_TAIL
        and q_ifn["n_Q4"] >= MIN_Q_TAIL
    )
    mhc_down_q = bool(
        np.isfinite(q_mhc.get("delta_Q4_minus_Q1")) and q_mhc["delta_Q4_minus_Q1"] < 0
        and q_mhc["n_Q1"] >= MIN_Q_TAIL
        and q_mhc["n_Q4"] >= MIN_Q_TAIL
    )
    join = bool(tnk_inverse or ifnm_down_cont or ifnm_down_q or ifn_down_q or mhc_down_q)
    decision = "JOIN" if join else "NO-JOIN"

    summary = {
        "dataset": "E-MTAB-13526",
        "citation": "De Zuani, Xue, Park, et al. Nat Commun 2024 (Cvejic NSCLC atlas)",
        "pmid": "38782901",
        "target": "CLDN4-only",
        "dual_high": False,
        "unit": "patient",
        "subset": "CD235a- tumor Cell Ranger mtx (15 lanes)",
        "skipped_h5ads": [
            "10X_Lung_Tumour_Annotated_v2.h5ad (58.74 GB)",
            "10X_Lung_Healthy_Background_Annotated_v2.h5ad (45.51 GB)",
        ],
        "n_cells": int(n),
        "n_tumor_patients": int(per["patient"].nunique()),
        "n_eligible_malig": int(elig["patient"].nunique()),
        "CLDN4_present": "CLDN4" in present,
        "n_IFN_genes": len(ifn_genes),
        "n_MHC_genes": len(mhc_genes),
        "primary_CLDN4_vs_tnk": tnk,
        "primary_CLDN4_vs_ifnmhc": ifnm,
        "q4q1_IFN_MHC": q_ifnm,
        "q4q1_IFN": q_ifn,
        "q4q1_MHC": q_mhc,
        "sign_tnk_inverse": tnk_inverse,
        "sign_ifnmhc_down_continuous": ifnm_down_cont,
        "sign_ifnmhc_down_q4q1": ifnm_down_q,
        "decision": decision,
        "concordant_pool": "GSE123902+GSE131907+GSE205335+GSE189357",
    }
    (out / "tables" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

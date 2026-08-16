#!/usr/bin/env python3
"""Patient-level Spearman: malignant TACSTD2/CLDN4 vs T/NK in E-MTAB-13526.

Extra n only. This atlas is treatment-naive NSCLC (De Zuani / Cvejic,
Nat Commun 2024, PMID 38821935) — not an ICI-response cohort. GSE207422
and other user A3/A* analyses are not re-run here.

Unit of analysis is the patient. CD235a- tumor lanes only (RBC-depleted;
not CD45-sorted), so T/NK fractions are not artificially inflated by
immune FACS. Epithelial restriction is tested as paired %positive
malignant-like vs T/NK.
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


def histolabel(disease: str) -> str:
    d = str(disease).lower()
    if "adenocarcinoma" in d:
        return "LUAD"
    if "squamous" in d:
        return "LUSC"
    return "NSCLC_other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extracted", type=Path, default=Path("data/emtab13526/extracted"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/emtab13526_tacstd2"))
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    audit = {
        "dataset": "E-MTAB-13526",
        "pmid": "38821935",
        "citation": "De Zuani, Xue, Park, et al. Nat Commun 2024 (Cvejic NSCLC atlas)",
        "n_patients_paper": 25,
        "n_patients_arrayexpress": 24,
        "treatment": "treatment-naive NSCLC tumor + matched non-involved lung",
        "is_blood_only": False,
        "is_ici_response_cohort": False,
        "public_mpr_r_labels": False,
        "gse207422_reanalyzed": False,
        "subset_used": "CD235a- tumor Cell Ranger mtx (15 lanes / 12 patients)",
        "skipped_full_atlas": (
            "Author annotated h5ads are 45.5 GB (background/healthy) + 58.7 GB "
            "(tumour). CD45+/MDSC-sorted tumor lanes omitted because they lack "
            "malignant epithelium and inflate T/NK fractions."
        ),
        "label_search": {
            "sdrf_response_fields": [],
            "note": (
                "ArrayExpress SDRF has disease, FACS, sampling site, TNM. "
                "No ICI, MPR, pCR, R/NR, or DCB labels. Paper states "
                "treatment-naive resections."
            ),
        },
    }

    gene_index_path = args.extracted / "gene_index.json"
    if not gene_index_path.exists():
        (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
        (out / "STOP.md").write_text(
            "# E-MTAB-13526 extraction incomplete\n\n"
            "Run `python3 methods/emtab13526_tacstd2/extract.py` first.\n"
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
            "After a real gene-name check of the E-MTAB-13526 features table, "
            "neither target was present. No figure invented.\n"
        )
        print("STOP: both genes absent")
        return 0

    packed = np.load(args.extracted / "gene_panel.npz")
    expr = {k: packed[k] for k in packed.files}
    meta = pd.read_csv(args.extracted / "cell_metadata.tsv", sep="\t")
    n = len(meta)
    audit["n_cells"] = int(n)
    audit["library_size"] = "Cell Ranger total UMI after author-like QC"

    for g in TARGETS + [x for xs in LINEAGES.values() for x in xs] + NORMAL_LUNG:
        if g in expr and g not in meta.columns:
            meta[g] = expr[g]

    lineage = assign_lineage(expr, n)
    meta["lineage"] = lineage
    normal_score = score(expr, [g for g in NORMAL_LUNG if g in expr], n)
    meta["malignant_like"] = (meta["lineage"] == "epithelial") & (normal_score <= 0.05)
    meta["tnk"] = meta["lineage"].isin(["T", "NK"])
    meta["histology"] = meta["disease"].map(histolabel) if "disease" in meta.columns else "NSCLC_other"

    rows = []
    for patient, pdf in meta.groupby("patient"):
        n_cells = len(pdf)
        n_epi = int((pdf["lineage"] == "epithelial").sum())
        n_mal = int(pdf["malignant_like"].sum())
        n_tnk = int(pdf["tnk"].sum())
        rec = {
            "patient": patient,
            "tissue": "Tumor",
            "histology": pdf["histology"].iloc[0] if "histology" in pdf else "",
            "disease": pdf["disease"].iloc[0] if "disease" in pdf else "",
            "n_lanes": int(pdf["sample"].nunique()),
            "n_cells": n_cells,
            "n_epithelial": n_epi,
            "n_malignant_like": n_mal,
            "n_tnk": n_tnk,
            "n_T": int((pdf["lineage"] == "T").sum()),
            "n_NK": int((pdf["lineage"] == "NK").sum()),
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

    per = pd.DataFrame(rows).sort_values("patient")
    per.to_csv(out / "per_patient_metrics.tsv", sep="\t", index=False)

    tests = []
    tumor = per.copy()
    for gene in TARGETS:
        if gene not in present:
            tests.append(
                {
                    "contrast": f"{gene} absent from processed matrices",
                    "n": 0,
                    "spearman_rho": np.nan,
                    "spearman_p": np.nan,
                    "note": "gene_absent",
                }
            )
            continue
        for score_col, label in (
            (f"{gene}_mean_log1p", f"tumor malig-like {gene} mean_log1p vs T/NK fraction"),
            (f"{gene}_mean_log1p_cp10k", f"tumor malig-like {gene} mean log1p(CP10k) vs T/NK fraction"),
            (f"{gene}_pct_pos", f"tumor malig-like {gene} %pos vs T/NK fraction"),
            (f"epi_{gene}_mean_log1p", f"tumor epithelial {gene} mean_log1p vs T/NK fraction"),
        ):
            elig = tumor[tumor["eligible_malig"]] if "malig" in label else tumor[tumor["eligible_epi"]]
            tests.append(spearman_block(elig[score_col], elig["tnk_fraction"], label + " (eligible)"))
            tests.append(spearman_block(tumor[score_col], tumor["tnk_fraction"], label + " (all CD235a- tumor patients)"))

        robust = tumor[tumor["n_malignant_like"] >= 50]
        tests.append(
            spearman_block(
                robust[f"{gene}_mean_log1p"],
                robust["tnk_fraction"],
                f"tumor malig-like {gene} mean_log1p vs T/NK (n_malig>=50)",
            )
        )
        for histo in ("LUAD", "LUSC"):
            sub = tumor[(tumor["histology"] == histo) & tumor["eligible_malig"]]
            tests.append(
                spearman_block(
                    sub[f"{gene}_mean_log1p"],
                    sub["tnk_fraction"],
                    f"tumor malig-like {gene} vs T/NK ({histo} eligible)",
                )
            )

    tests.append(
        {
            "contrast": "response / MPR / R vs TACSTD2 or CLDN4",
            "n": 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "no_public_response_labels_treatment_naive_not_ICI",
        }
    )
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(out / "association_statistics.tsv", sep="\t", index=False)

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

    lineage_counts = (
        meta.groupby(["patient", "lineage"]).size().unstack(fill_value=0).reset_index()
    )
    lineage_counts.to_csv(out / "lineage_counts.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.7), constrained_layout=True)
    for ax, gene in zip(axes, TARGETS):
        if gene not in present:
            ax.set_axis_off()
            ax.set_title(f"{gene} absent")
            continue
        elig = tumor[tumor["eligible_malig"]]
        x = elig["tnk_fraction"].to_numpy(dtype=float)
        y = elig[f"{gene}_mean_log1p"].to_numpy(dtype=float)
        colors = elig["histology"].map({"LUAD": "#1f4e79", "LUSC": "#8b1e3f", "NSCLC_other": "#6b6b6b"})
        ax.scatter(x, y, c=colors.fillna("#6b6b6b"), s=46, zorder=3)
        for _, r in elig.iterrows():
            ax.annotate(
                str(r["patient"]),
                (r["tnk_fraction"], r[f"{gene}_mean_log1p"]),
                fontsize=7,
                xytext=(4, 3),
                textcoords="offset points",
            )
        blk = spearman_block(x, y, gene)
        ax.set_xlabel("T/NK fraction (CD235a− tumor cells)")
        ax.set_ylabel(f"Malignant-like {gene}\nmean log1p UMI")
        if blk["n"] >= 4:
            ax.set_title(f"{gene}  n={blk['n']}  ρ={blk['spearman_rho']:.2f}  p={blk['spearman_p']:.2g}")
        else:
            ax.set_title(f"{gene}  n={blk['n']}  (too few)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("E-MTAB-13526 Cvejic NSCLC atlas (not ICI) · extra scRNA n", fontsize=10)
    fig.savefig(out / "fig_emtab13526_tacstd2_cldn4_vs_tnk.png", dpi=180)
    fig.savefig(out / "fig_emtab13526_tacstd2_cldn4_vs_tnk.pdf")
    plt.close(fig)

    primary = [
        t
        for t in tests
        if "malig-like TACSTD2 mean_log1p vs T/NK" in t["contrast"] and "eligible" in t["contrast"]
    ]
    primary_c = [
        t
        for t in tests
        if "malig-like CLDN4 mean_log1p vs T/NK" in t["contrast"] and "eligible" in t["contrast"]
    ]
    summary = {
        **audit,
        "n_tumor_patients": int(tumor["patient"].nunique()) if len(tumor) else 0,
        "n_eligible_malig_tumor": int(tumor["eligible_malig"].sum()) if len(tumor) else 0,
        "primary_TACSTD2_vs_tnk": primary[0] if primary else None,
        "primary_CLDN4_vs_tnk": primary_c[0] if primary_c else None,
        "response_test": "not_possible_no_public_labels",
        "honest_unit": "patient",
        "restriction": restriction,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

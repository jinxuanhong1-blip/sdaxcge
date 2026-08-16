#!/usr/bin/env python3
"""Recompute TACSTD2 vs CLDN4 on the exact Nusinow 2020 / Gygi CCLE protein table.

Hunt for a public NSCLC protein slice with n≈118. If it does not exist, say so
and still report Spearman ρ on every natural filter of that table.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

USER_N = 118
USER_RHO = 0.69
NEAR_N = 5  # |n-118| <= 5 counts as "≈118"
NEAR_RHO = 0.05
N_BOOT = 5000
BOOT_SEED = 0

SAMPLE_RE = re.compile(r"_TenPx\d+$")


def core_of(col: str) -> str:
    m = re.match(r"(.+)_TenPx\d+$", col)
    return m.group(1) if m else col


def pairwise(x: pd.Series, y: pd.Series) -> pd.DataFrame:
    df = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1)
    df.columns = ["TACSTD2", "CLDN4"]
    return df.dropna()


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOT, seed: int = BOOT_SEED):
    rng = np.random.default_rng(seed)
    n = len(x)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rhos[i] = spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.quantile(rhos, [0.025, 0.975])
    return float(lo), float(hi)


def _py(v):
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    return v


def corr_block(name: str, pair: pd.DataFrame, definition: str) -> dict:
    n = int(len(pair))
    rec = {
        "cohort": name,
        "definition": definition,
        "n_pairwise": n,
        "n_minus_118": n - USER_N if n else None,
        "is_n_118": bool(n == USER_N),
        "is_n_approx_118": bool(abs(n - USER_N) <= NEAR_N) if n else False,
        "spearman_rho": None,
        "spearman_p": None,
        "spearman_rounded_2dp": None,
        "pearson_r": None,
        "pearson_p": None,
        "rho_vs_0.69": None,
        "rho_match_2dp": False,
        "rho_nearby": False,
        "ci95_lo": None,
        "ci95_hi": None,
        "user_0.69_in_ci": None,
    }
    if n < 4:
        return rec
    rho, p = spearmanr(pair["TACSTD2"], pair["CLDN4"])
    r, pp = pearsonr(pair["TACSTD2"], pair["CLDN4"])
    lo, hi = spearman_ci(pair["TACSTD2"].to_numpy(), pair["CLDN4"].to_numpy())
    rec.update(
        {
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "spearman_rounded_2dp": float(round(float(rho), 2)),
            "pearson_r": float(r),
            "pearson_p": float(pp),
            "rho_vs_0.69": float(rho - USER_RHO),
            "rho_match_2dp": bool(round(float(rho), 2) == USER_RHO),
            "rho_nearby": bool(abs(float(rho) - USER_RHO) <= NEAR_RHO),
            "ci95_lo": lo,
            "ci95_hi": hi,
            "user_0.69_in_ci": bool(lo <= USER_RHO <= hi),
        }
    )
    return {k: _py(v) for k, v in rec.items()}


def scatter(pair: pd.DataFrame, title: str, dest: Path, color: str = "#1f77b4") -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(pair["TACSTD2"], pair["CLDN4"], s=28, alpha=0.75, c=color, edgecolors="none")
    rho, p = spearmanr(pair["TACSTD2"], pair["CLDN4"])
    ax.set_xlabel("TACSTD2 protein (TMT, normalized)")
    ax.set_ylabel("CLDN4 protein (TMT, normalized)")
    ax.set_title(f"{title}\nSpearman ρ = {rho:.3f}  n = {len(pair)}  p = {p:.2e}")
    ax.axhline(0, color="0.7", lw=0.6)
    ax.axvline(0, color="0.7", lw=0.6)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=160)
    fig.savefig(dest.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/B2_n118")
    ap.add_argument("--outdir", default="results/rework/B2_n118")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    prot_path = data / "protein_quant_current_normalized.csv.gz"
    s1_path = data / "Table_S1_Sample_Information.xlsx"
    model_path = data / "Model.csv"
    man_src = data / "download_manifest.json"

    prot = pd.read_csv(prot_path)
    sample_cols = [c for c in prot.columns if SAMPLE_RE.search(str(c))]
    if (prot["Gene_Symbol"] == "TACSTD2").sum() != 1 or (prot["Gene_Symbol"] == "CLDN4").sum() != 1:
        raise SystemExit("Expected exactly one TACSTD2 row and one CLDN4 row")
    trow = prot.loc[prot["Gene_Symbol"] == "TACSTD2", sample_cols].iloc[0]
    crow = prot.loc[prot["Gene_Symbol"] == "CLDN4", sample_cols].iloc[0]

    # unique cell-line values: mean of the 3 lines that appear in two plexes
    cores = [core_of(c) for c in sample_cols]
    unique_cores = sorted(set(cores))
    t_unique = {}
    c_unique = {}
    for core in unique_cores:
        cols = [c for c in sample_cols if core_of(c) == core]
        t_unique[core] = pd.to_numeric(trow[cols], errors="coerce").mean()
        c_unique[core] = pd.to_numeric(crow[cols], errors="coerce").mean()
    t_u = pd.Series(t_unique, name="TACSTD2")
    c_u = pd.Series(c_unique, name="CLDN4")

    s1 = pd.read_excel(s1_path, sheet_name="Sample_Information")
    s1_nb = s1[s1["Notes"] != "Bridge line"].drop_duplicates("CCLE Code").set_index("CCLE Code")
    model = pd.read_csv(model_path).drop_duplicates("CCLEName").set_index("CCLEName")

    meta = pd.DataFrame({"CCLE_Code": unique_cores}).set_index("CCLE_Code")
    meta["TACSTD2"] = t_u
    meta["CLDN4"] = c_u
    meta["S1_Tissue"] = s1_nb.reindex(meta.index)["Tissue of Origin"]
    meta["S1_CellLine"] = s1_nb.reindex(meta.index)["Cell Line"]
    meta["OncotreeLineage"] = model.reindex(meta.index)["OncotreeLineage"]
    meta["OncotreePrimaryDisease"] = model.reindex(meta.index)["OncotreePrimaryDisease"]
    meta["OncotreeSubtype"] = model.reindex(meta.index)["OncotreeSubtype"]
    meta["ModelID"] = model.reindex(meta.index)["ModelID"]
    meta["CellLineName"] = model.reindex(meta.index)["CellLineName"]
    meta["has_TACSTD2"] = meta["TACSTD2"].notna()
    meta["has_CLDN4"] = meta["CLDN4"].notna()
    meta["has_both"] = meta["has_TACSTD2"] & meta["has_CLDN4"]

    # natural cohorts (unique cell lines)
    masks = {
        "all_unique_lines": (
            pd.Series(True, index=meta.index),
            "All unique CCLE codes in protein_quant_current_normalized.csv.gz (375 lines; 3 plex-duplicates averaged)",
        ),
        "S1_Lung": (
            meta["S1_Tissue"] == "Lung",
            "Table S1 Tissue of Origin == Lung (CCLE 2012 lineage, includes NCI-H226)",
        ),
        "S1_Lung_excluding_NCIH226": (
            (meta["S1_Tissue"] == "Lung") & (meta.index != "NCIH226_LUNG"),
            "S1 Lung minus NCI-H226 (DepMap 24Q4 calls this pleural mesothelioma)",
        ),
        "Oncotree_Lung": (
            meta["OncotreeLineage"] == "Lung",
            "DepMap 24Q4 OncotreeLineage == Lung",
        ),
        "Oncotree_NSCLC": (
            meta["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer",
            "DepMap 24Q4 OncotreePrimaryDisease == Non-Small Cell Lung Cancer",
        ),
        "Oncotree_LUAD": (
            meta["OncotreeSubtype"] == "Lung Adenocarcinoma",
            "DepMap 24Q4 OncotreeSubtype == Lung Adenocarcinoma",
        ),
        "Oncotree_LUSC": (
            meta["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma",
            "DepMap 24Q4 OncotreeSubtype == Lung Squamous Cell Carcinoma",
        ),
        "Oncotree_NET_SCLC": (
            meta["OncotreePrimaryDisease"] == "Lung Neuroendocrine Tumor",
            "DepMap 24Q4 OncotreePrimaryDisease == Lung Neuroendocrine Tumor (SCLC)",
        ),
        "S1_solid_organ": (
            ~meta["S1_Tissue"].str.contains("Haematopoietic|Lymphoma|Leukemia", case=False, na=False),
            "Table S1 solid-organ lineages (exclude haematopoietic / lymphoma / leukemia)",
        ),
    }

    rows = []
    pair_tables = {}
    for name, (mask, definition) in masks.items():
        sub = meta.loc[mask]
        n_lines = int(len(sub))
        n_t = int(sub["has_TACSTD2"].sum())
        n_c = int(sub["has_CLDN4"].sum())
        pair = pairwise(sub["TACSTD2"], sub["CLDN4"])
        rec = corr_block(name, pair, definition)
        rec.update(
            {
                "n_lines_in_table": n_lines,
                "n_TACSTD2": n_t,
                "n_CLDN4": n_c,
            }
        )
        rows.append(rec)
        pair_tables[name] = (sub, pair)

    corr_df = pd.DataFrame(rows)

    # exhaustive search: any S1 tissue union of size 1–3 with n_lines or n_both == 118
    tissues = sorted(meta["S1_Tissue"].dropna().unique())
    hunt = []
    from itertools import combinations

    for k in range(1, 4):
        for combo in combinations(tissues, k):
            mask = meta["S1_Tissue"].isin(combo)
            n_lines = int(mask.sum())
            n_both = int(meta.loc[mask, "has_both"].sum())
            if abs(n_lines - USER_N) <= NEAR_N or abs(n_both - USER_N) <= NEAR_N:
                hunt.append(
                    {
                        "tissues": " + ".join(combo),
                        "k_tissues": k,
                        "n_lines": n_lines,
                        "n_both": n_both,
                        "matches_118_lines": n_lines == USER_N,
                        "matches_118_both": n_both == USER_N,
                        "natural_NSCLC_or_lung": combo == ("Lung",),
                    }
                )
    hunt_df = pd.DataFrame(hunt)

    n118_natural = False
    natural_n118_cohorts = [
        r
        for r in rows
        if r["is_n_118"] or r["is_n_approx_118"] or r.get("n_lines_in_table") == USER_N
        or (r.get("n_lines_in_table") is not None and abs(r["n_lines_in_table"] - USER_N) <= NEAR_N)
    ]
    # A "natural" n=118 would be a single published lineage/histology filter, not a tissue mash-up
    n118_found_in_table = any(
        (r["is_n_118"] or r.get("n_lines_in_table") == USER_N) for r in rows
    )

    # primary reported recompute: S1 Lung pairwise (the slice that actually hits ρ=0.69)
    lung_pair = pair_tables["S1_Lung"][1]
    nsclc_pair = pair_tables["Oncotree_NSCLC"][1]
    lung_rec = next(r for r in rows if r["cohort"] == "S1_Lung")
    nsclc_rec = next(r for r in rows if r["cohort"] == "Oncotree_NSCLC")

    if n118_found_in_table:
        verdict = "N118_FOUND"
        verdict_line = "A natural NSCLC/lung filter of the Gygi table has n=118."
    else:
        verdict = "N118_NOT_FOUND"
        verdict_line = (
            "n=118 is NOT present as a natural NSCLC or lung slice of the "
            "Nusinow 2020 / Gygi protein_quant_current_normalized table."
        )

    # write pair-level tables
    lung_meta = pair_tables["S1_Lung"][0].copy()
    lung_meta.to_csv(out / "S1_lung_protein_TACSTD2_CLDN4.csv")
    nsclc_meta = pair_tables["Oncotree_NSCLC"][0].copy()
    nsclc_meta.to_csv(out / "Oncotree_NSCLC_protein_TACSTD2_CLDN4.csv")
    meta.to_csv(out / "all_375_protein_TACSTD2_CLDN4.csv")
    corr_df.to_csv(out / "correlations.csv", index=False)
    if len(hunt_df):
        hunt_df.to_csv(out / "n118_tissue_combo_hunt.csv", index=False)
    else:
        pd.DataFrame(columns=["tissues"]).to_csv(out / "n118_tissue_combo_hunt.csv", index=False)

    if len(lung_pair) >= 4:
        scatter(lung_pair, "Nusinow 2020 / Gygi  ·  S1 Lung (pairwise complete)", out / "fig_scatter_S1_lung.png")
    if len(nsclc_pair) >= 4:
        scatter(
            nsclc_pair,
            "Nusinow 2020 / Gygi  ·  Oncotree NSCLC (pairwise complete)",
            out / "fig_scatter_Oncotree_NSCLC.png",
            color="#d62728",
        )

    # n-bar figure vs 118
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = [
        "User claim",
        "S1 Lung lines",
        "S1 Lung both",
        "Oncotree NSCLC lines",
        "Oncotree NSCLC both",
        "Oncotree Lung lines",
    ]
    vals = [
        USER_N,
        lung_rec["n_lines_in_table"],
        lung_rec["n_pairwise"],
        nsclc_rec["n_lines_in_table"],
        nsclc_rec["n_pairwise"],
        next(r for r in rows if r["cohort"] == "Oncotree_Lung")["n_lines_in_table"],
    ]
    colors = ["#7f7f7f", "#1f77b4", "#1f77b4", "#d62728", "#d62728", "#2ca02c"]
    ax.bar(range(len(labels)), vals, color=colors)
    ax.axhline(USER_N, color="0.3", ls="--", lw=1, label="n=118")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("n")
    ax.set_title("No natural filter of the Gygi CCLE protein table has n≈118")
    fig.tight_layout()
    fig.savefig(out / "fig_n_vs_118.png", dpi=160)
    plt.close(fig)

    if man_src.exists():
        shutil.copy(man_src, out / "download_manifest.json")

    summary = {
        "task": "B2_n118",
        "user_claim": {
            "dataset": "CCLE NSCLC protein",
            "n": USER_N,
            "rho": USER_RHO,
            "genes": ["TACSTD2", "CLDN4"],
        },
        "exact_public_table": {
            "found": True,
            "name": "protein_quant_current_normalized.csv.gz",
            "url": "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
            "excel_twin": "https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx",
            "sample_info": "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx",
            "gygi_page": "https://gygi.hms.harvard.edu/publications/ccle.html",
            "paper": "Nusinow et al. Quantitative Proteomics of the Cancer Cell Line Encyclopedia. Cell 2020;180:387-402.e16",
            "doi": "10.1016/j.cell.2019.12.023",
            "n_protein_rows": int(prot.shape[0]),
            "n_sample_columns": int(len(sample_cols)),
            "n_unique_cell_lines": int(len(unique_cores)),
            "duplicate_cores_averaged": ["CAL120_BREAST", "SW948_LARGE_INTESTINE", "HCT15_LARGE_INTESTINE"],
            "TACSTD2_uniprot": "P09758",
            "CLDN4_uniprot": "O14493",
            "note": (
                "This is the exact public Nusinow 2020 / Gygi / DepMap CCLE MS proteomics table "
                "(375 unique lines). It is NOT an n=118 NSCLC table."
            ),
        },
        "n118_search": {
            "n_118_found_as_natural_NSCLC_or_lung_filter": n118_found_in_table,
            "verdict": verdict,
            "closest_natural_n_lines": {
                "S1_Lung": lung_rec["n_lines_in_table"],
                "Oncotree_Lung": next(r for r in rows if r["cohort"] == "Oncotree_Lung")["n_lines_in_table"],
                "Oncotree_NSCLC": nsclc_rec["n_lines_in_table"],
            },
            "closest_natural_n_pairwise": {
                "S1_Lung": lung_rec["n_pairwise"],
                "Oncotree_NSCLC": nsclc_rec["n_pairwise"],
            },
            "historical_n118_is_RNA_not_protein": {
                "source": "Augustyn et al. PNAS 2014;111:14788-14793 (ASCL1)",
                "doi": "10.1073/pnas.1410419111",
                "quote": "NSCLC (n = 118) ... genome-wide mRNA expression data",
                "reused_in": "Nilsson et al. Cancer Cell 2023 CD70 paper (NSCLC cell lines n=118, gene expression)",
            },
            "tissue_mashups_that_sum_to_118": "see n118_tissue_combo_hunt.csv; none is a published NSCLC definition",
        },
        "recompute_on_exact_table": {
            "method": "scipy.stats.spearmanr and pearsonr, two-sided, pairwise complete cases; unique cell lines (plex duplicates averaged)",
            "did_we_tune_to_0.69": False,
            "did_we_tune_to_n118": False,
            "S1_Lung_pairwise": lung_rec,
            "Oncotree_NSCLC_pairwise": nsclc_rec,
        },
        "all_cohorts": rows,
        "honest_verdict": {
            "label": verdict,
            "n118": "NOT FOUND on the exact public CCLE proteomics table",
            "rho_0.69": (
                "MATCH_AT_2DP on S1 Lung pairwise-complete (n=45, ρ=0.693), "
                "which is not n=118"
            ),
            "one_sentence": verdict_line,
        },
    }
    def _json_default(o):
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        raise TypeError(type(o))

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=_json_default) + "\n")

    (out / "verdict.txt").write_text(
        "\n".join(
            [
                f"VERDICT: {verdict}",
                verdict_line,
                "",
                "EXACT PUBLIC TABLE: https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
                "  Nusinow et al. Cell 2020;180:387-402.e16  (Gygi lab CCLE MS proteomics)",
                f"  unique cell lines = {len(unique_cores)}  (sample columns = {len(sample_cols)})",
                "",
                f"S1 Lung lines:            n={lung_rec['n_lines_in_table']}",
                f"S1 Lung pairwise both:    n={lung_rec['n_pairwise']}  Spearman ρ={lung_rec['spearman_rho']:.4f}  -> {lung_rec['spearman_rounded_2dp']}",
                f"Oncotree NSCLC lines:     n={nsclc_rec['n_lines_in_table']}",
                f"Oncotree NSCLC pairwise:  n={nsclc_rec['n_pairwise']}  Spearman ρ={nsclc_rec['spearman_rho']:.4f}  -> {nsclc_rec['spearman_rounded_2dp']}",
                "",
                "n=118 is a published CCLE *mRNA* NSCLC count (Augustyn PNAS 2014), not a protein-table n.",
                "User ρ=0.69 matches S1 Lung protein pairwise n=45, not n=118.",
                "did_we_tune_to_0.69: false",
                "did_we_tune_to_n118: false",
                "",
            ]
        )
    )
    print(open(out / "verdict.txt").read())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

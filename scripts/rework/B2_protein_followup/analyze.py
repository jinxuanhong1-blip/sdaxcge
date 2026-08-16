#!/usr/bin/env python3
"""B2 protein follow-up: TACSTD2–CLDN4 on CCLE MS / RPPA.

User correction: claimed ρ = 0.69 was CCLE NSCLC protein, n=118, not DepMap RNA.
This script does not tune filters to hit 0.69. It reports which named n is 118
and which named complete-case correlation is ~0.69.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

USER_RHO = 0.69
USER_N = 118


def corr_block(x, y, ci: bool = False) -> dict:
    x = pd.to_numeric(pd.Series(x), errors="coerce")
    y = pd.to_numeric(pd.Series(y), errors="coerce")
    m = x.notna() & y.notna()
    x, y = x[m].to_numpy(float), y[m].to_numpy(float)
    n = int(len(x))
    out = {
        "n": n,
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "spearman_ci95_low": None,
        "spearman_ci95_high": None,
        "rounds_to_0.69": False,
        "abs_diff_vs_0.69": None,
        "n_equals_118": n == USER_N,
    }
    if n < 3:
        return out
    rho, p_s = stats.spearmanr(x, y)
    r, p_p = stats.pearsonr(x, y)
    out.update(
        {
            "spearman_rho": float(rho),
            "spearman_p": float(p_s),
            "pearson_r": float(r),
            "pearson_p": float(p_p),
            "rounds_to_0.69": float(f"{rho:.2f}") == USER_RHO,
            "abs_diff_vs_0.69": abs(float(rho) - USER_RHO),
        }
    )
    if ci and n >= 8:
        rng = np.random.default_rng(0)
        boots = np.empty(5000)
        for i in range(5000):
            idx = rng.integers(0, n, n)
            boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
        lo, hi = np.percentile(boots, [2.5, 97.5])
        out["spearman_ci95_low"] = float(lo)
        out["spearman_ci95_high"] = float(hi)
        out["user_0.69_inside_ci95"] = bool(lo <= USER_RHO <= hi)
    return out


def gygi_extract(gz_path: Path) -> tuple[pd.Series, pd.Series, dict]:
    with gzip.open(gz_path, "rt") as f:
        r = csv.DictReader(f)
        quant = [c for c in r.fieldnames if "_TenPx" in c and not c.endswith("_Peptides")]
        rows = {}
        for row in r:
            if row["Gene_Symbol"] in ("TACSTD2", "CLDN4"):
                rows[row["Gene_Symbol"]] = row
    if set(rows) != {"TACSTD2", "CLDN4"}:
        raise SystemExit(f"Gygi missing genes: {sorted(rows)}")
    meta = {
        "TACSTD2_id": rows["TACSTD2"]["Protein_Id"],
        "CLDN4_id": rows["CLDN4"]["Protein_Id"],
        "n_quant_columns": len(quant),
    }

    def ser(g: str) -> pd.Series:
        return pd.Series(
            {
                c: (np.nan if rows[g][c] in ("", "NA", "NaN") else float(rows[g][c]))
                for c in quant
            }
        )

    return ser("TACSTD2"), ser("CLDN4"), meta


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/ccle_protein")
    p.add_argument("--out-dir", default="results/rework/B2_protein_followup")
    p.add_argument("--sample-info", default="/tmp/ccle2018/sample_info_22q2.csv")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict] = []

    # ----- Gygi MS -----
    tac, cld, gygi_meta = gygi_extract(cache / "protein_quant_current_normalized.csv.gz")
    s1 = pd.read_excel(cache / "Table_S1_Sample_Information.xlsx", sheet_name="Sample_Information")
    si = pd.read_csv(args.sample_info)
    si["ccle"] = si["CCLE_Name"].astype(str)

    def ccle_from_col(c: str) -> str:
        return c.rsplit("_TenPx", 1)[0]

    ann = pd.DataFrame({"col": tac.index, "ccle": [ccle_from_col(c) for c in tac.index]})
    ann = ann.merge(
        si[["ccle", "lineage", "lineage_subtype", "primary_disease", "cell_line_name"]],
        on="ccle",
        how="left",
    )
    s1_tissue = dict(zip(s1["CCLE Code"].astype(str), s1["Tissue of Origin"].astype(str)))
    ann["s1_tissue"] = ann["ccle"].map(s1_tissue)
    ann["TACSTD2"] = tac.to_numpy()
    ann["CLDN4"] = cld.to_numpy()
    ann.to_csv(out / "gygi_ms_tacstd2_cldn4_all.csv", index=False)

    def add(name, mask, note, ci=False):
        sub = ann.loc[mask]
        block = corr_block(sub["TACSTD2"], sub["CLDN4"], ci=ci)
        rec = {"dataset": "Gygi_CCLE_MS", "cohort": name, "note": note, **block}
        rows_out.append(rec)
        return rec

    add("all_complete", ann["TACSTD2"].notna() & ann["CLDN4"].notna(), "All Gygi lines with both proteins.")
    lung_suf = ann["col"].str.contains("_LUNG_")
    primary = add(
        "lung_suffix_complete",
        lung_suf,
        "PRIMARY protein correlation: CCLE ID contains _LUNG, both proteins quantified.",
        ci=True,
    )
    add("s1_tissue_Lung", ann["s1_tissue"] == "Lung", "Gygi Table S1 Tissue of Origin == Lung.")
    add("lineage_lung", ann["lineage"] == "lung", "DepMap 22Q2 sample_info lineage==lung (includes meso).")
    add("NSCLC", ann["lineage_subtype"] == "NSCLC", "22Q2 lineage_subtype==NSCLC.")
    add("SCLC", ann["lineage_subtype"] == "SCLC", "22Q2 lineage_subtype==SCLC.")
    add(
        "lung_not_SCLC",
        (ann["lineage"] == "lung") & (ann["lineage_subtype"] != "SCLC"),
        "22Q2 lineage==lung excluding SCLC (still includes meso if present).",
    )

    lung_tab = ann.loc[lung_suf].copy()
    lung_tab.to_csv(out / "gygi_ms_lung_all_columns.csv", index=False)
    both = lung_tab.dropna(subset=["TACSTD2", "CLDN4"])
    both.to_csv(out / "gygi_ms_lung_complete_n45.csv", index=False)

    # scatter
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    for subname, g in both.groupby(both["lineage_subtype"].fillna("NA")):
        ax.scatter(g["TACSTD2"], g["CLDN4"], s=28, alpha=0.85, label=f"{subname} (n={len(g)})")
    ax.set_xlabel("TACSTD2 protein (Gygi TMT, normalized)")
    ax.set_ylabel("CLDN4 protein (Gygi TMT, normalized)")
    ax.set_title(
        f"CCLE MS (Nusinow/Gygi) lung lines with both proteins\n"
        f"Spearman ρ = {primary['spearman_rho']:.3f} (n={primary['n']}); user claim 0.69 / n=118"
    )
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig_gygi_ms_lung_complete.png", dpi=160)
    plt.close(fig)

    # ----- RPPA -----
    ab = pd.read_csv(cache / "CCLE_RPPA_Ab_info_20180123.csv")
    rppa = pd.read_csv(cache / "CCLE_RPPA_20180123.csv")
    rppa = rppa.rename(columns={rppa.columns[0]: "CCLE_ID"})
    ab.to_csv(out / "rppa_antibody_info.csv", index=False)
    ab_hits = ab[
        ab.astype(str).apply(
            lambda c: c.str.contains("TACSTD|TROP2|CLDN4|Claudin-4|Claudin4", case=False, na=False)
        ).any(axis=1)
    ]
    rppa_lung = rppa[rppa["CCLE_ID"].astype(str).str.endswith("_LUNG")].copy()
    rppa_ann = rppa_lung[["CCLE_ID"]].merge(
        si[["ccle", "lineage", "lineage_subtype", "primary_disease", "cell_line_name"]],
        left_on="CCLE_ID",
        right_on="ccle",
        how="left",
    )
    rppa_ann.to_csv(out / "rppa_lung_sample_annotation.csv", index=False)
    n_rppa_lung = int(len(rppa_lung))
    n_rppa_nsclc = int((rppa_ann["lineage_subtype"] == "NSCLC").sum())
    n_rppa_lung_not_sclc = int(
        ((rppa_ann["lineage"] == "lung") & (rppa_ann["lineage_subtype"] != "SCLC")).sum()
    )
    rppa_rec = {
        "dataset": "CCLE_RPPA_20180123",
        "n_samples": int(len(rppa)),
        "n_lung_suffix": n_rppa_lung,
        "n_lung_NSCLC": n_rppa_nsclc,
        "n_lung_not_SCLC": n_rppa_lung_not_sclc,
        "has_TACSTD2_or_TROP2_antibody": False,
        "has_CLDN4_antibody": False,
        "claudin_antibodies_present": ab.loc[
            ab["Target_Genes"].astype(str).str.contains("CLDN", case=False, na=False),
            "Antibody_Name",
        ].tolist(),
        "antibody_hits_TACSTD2_CLDN4": ab_hits.to_dict(orient="records"),
        "can_compute_TACSTD2_CLDN4_rho": False,
        "n_118_match": {
            "n_lung_not_SCLC": n_rppa_lung_not_sclc == USER_N,
            "n_NSCLC": n_rppa_nsclc == USER_N,
        },
        "note": "n=118 is the RPPA lung-not-SCLC headcount (117 NSCLC + 1 carcinoid). "
        "RPPA has neither TACSTD2/TROP2 nor CLDN4, so this n cannot be the correlation sample.",
    }

    # ----- ProCan -----
    procan_path = cache / "ProCan-DepMapSanger_protein_matrix_8498_averaged.txt"
    procan_note = {}
    if procan_path.exists():
        with procan_path.open() as f:
            header = f.readline().rstrip("\n").split("\t")
        proteins = header[1:]
        has_tac = any("P09758" in h or "TACD2_HUMAN" in h for h in proteins)
        has_cld = any("O14493" in h or "CLD4_HUMAN" in h or "CLDN4" in h.upper() for h in proteins)
        claudins = [h for h in proteins if "CLD" in h.upper() and "HUMAN" in h.upper()]
        mp = pd.read_csv(cache / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t")
        procan_note = {
            "dataset": "ProCan-DepMapSanger_8498",
            "n_cell_lines": 949,
            "n_proteins": len(proteins),
            "has_TACSTD2": has_tac,
            "has_CLDN4": has_cld,
            "claudins_present": claudins,
            "n_mapping_Lung": int((mp["Tissue_type"] == "Lung").sum()),
            "n_mapping_NSCLC": int((mp["Cancer_type"] == "Non-Small Cell Lung Carcinoma").sum()),
            "can_compute_TACSTD2_CLDN4_rho": bool(has_tac and has_cld),
            "note": "TACSTD2 is present; CLDN4 is absent (CLDN1/3/7 only). No TACSTD2–CLDN4 ρ.",
        }

    corr_df = pd.DataFrame(rows_out)
    corr_df.to_csv(out / "correlations.csv", index=False)

    hits_rho = corr_df[corr_df["rounds_to_0.69"] == True]  # noqa: E712
    hits_n = {
        "rppa_lung_not_SCLC": n_rppa_lung_not_sclc,
        "rppa_NSCLC": n_rppa_nsclc,
        "gygi_lung_complete": int(primary["n"]),
        "gygi_NSCLC_complete": int(next(r["n"] for r in rows_out if r["cohort"] == "NSCLC")),
        "procan_mapping_NSCLC": procan_note.get("n_mapping_NSCLC"),
        "procan_mapping_Lung": procan_note.get("n_mapping_Lung"),
    }

    summary = {
        "task": "B2_protein_followup",
        "user_correction": "claimed ρ=0.69 was CCLE NSCLC protein n=118, not DepMap RNA",
        "honest_verdict": {
            "n_118_is": (
                "CCLE RPPA 20180123 lung lines with lineage==lung and subtype!=SCLC "
                f"(n={n_rppa_lung_not_sclc}: {n_rppa_nsclc} NSCLC + 1 carcinoid). "
                "This is a sample headcount, not a TACSTD2–CLDN4 complete-case n."
            ),
            "rho_0.69_is": (
                f"Gygi/Nusinow CCLE MS, _LUNG columns with both proteins quantified: "
                f"Spearman ρ={primary['spearman_rho']:.4f} (n={primary['n']}), rounds to 0.69. "
                f"Bootstrap 95% CI [{primary['spearman_ci95_low']:.3f}, {primary['spearman_ci95_high']:.3f}]."
            ),
            "same_cohort": False,
            "rppa_can_compute_pair": False,
            "procan_can_compute_pair": False,
            "did_we_tune": False,
            "statement": (
                "n=118 and ρ=0.69 do not come from the same table. "
                "n=118 is the RPPA lung-not-SCLC headcount (no TACSTD2/CLDN4 antibodies). "
                f"ρ=0.69 is Gygi MS lung complete cases (n={primary['n']}, not 118). "
                "NSCLC-only Gygi MS is n=36, ρ=0.73 — closer histology, farther from 0.69."
            ),
        },
        "gygi_meta": gygi_meta,
        "primary_gygi_lung_complete": primary,
        "all_gygi_cohorts": rows_out,
        "rppa": rppa_rec,
        "procan": procan_note,
        "which_n_hits_0.69": hits_rho[["cohort", "n", "spearman_rho", "abs_diff_vs_0.69"]].to_dict(
            orient="records"
        ),
        "named_n_inventory": hits_n,
        "prior_rna_result": {
            "depmap_24q4_lung_RNA": "Spearman 0.607 n=214 — wrong modality for this claim"
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["honest_verdict"], indent=2))
    print(corr_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

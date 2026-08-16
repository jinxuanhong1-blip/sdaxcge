#!/usr/bin/env python3
"""Recompute epcam_krt_mean with an epithelial-median floor (not all-cell p80)."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze import (
    IMMUNE_DEFS,
    MAL_DEFS,
    SCORE_DEFS,
    draw_grid,
    exact_wilcoxon,
    spearman_safe,
    MAL_LABELS,
    IMM_LABELS,
    SCORE_LABELS,
    RHO_HIGHLIGHT,
    sample_from_barcode,
)

HERE = Path(__file__).resolve().parent.parent
WANTED = ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC"]


def stream_genes(matrix_path: Path, wanted: list[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n = len(cell_ids)
        n_streamed = 0
        need = set(wanted)
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            n_streamed += 1
            if gene not in need:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: {arr.size} != {n}")
            found[gene] = arr
            need.remove(gene)
            print(f"  found {gene} at gene#{n_streamed} left={sorted(need)}", flush=True)
            if not need:
                break
        if need:
            raise SystemExit(f"missing genes: {need}")
    print(f"stream done stored={list(found)}", flush=True)
    return cell_ids, found


def main() -> None:
    matrix = Path("data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
    cells = pd.read_csv(HERE / "cell_calls.tsv.gz", sep="\t")
    print("streaming EPCAM/KRT/PTPRC", flush=True)
    cell_ids, expr = stream_genes(matrix, WANTED)
    if list(cells["barcode"].astype(str)) != list(cell_ids):
        raise SystemExit("barcode order mismatch")
    n_umi = cells["n_umi"].to_numpy(dtype=np.float64)
    scale = np.where(n_umi > 0, 1e4 / n_umi, 0.0)
    mods = []
    for g in ["EPCAM", "KRT8", "KRT18", "KRT19"]:
        mods.append(np.log1p(expr[g] * scale))
    epcam_krt_mod = np.mean(np.vstack(mods), axis=0)
    ptprc = expr["PTPRC"]
    is_epi = cells["lineage_mean"].to_numpy() == "Epithelial"
    ek_cut = float(max(0.40, np.median(epcam_krt_mod[is_epi])))
    is_ek_mean = (epcam_krt_mod >= ek_cut) & (ptprc < 1)
    print(f"ek_cut={ek_cut:.4f} n={int(is_ek_mean.sum())} epi={int(is_epi.sum())}", flush=True)
    cells["mal_epcam_krt_mean"] = is_ek_mean.astype(np.int8)
    cells.to_csv(HERE / "cell_calls.tsv.gz", sep="\t", index=False, compression="gzip")

    per = pd.read_csv(HERE / "per_sample.tsv", sep="\t")
    for i, row in per.iterrows():
        sub = cells.loc[(cells["sample"] == row["sample"]) & (cells["mal_epcam_krt_mean"] == 1)]
        per.at[i, "epcam_krt_mean_n"] = int(len(sub))
        if len(sub) == 0:
            per.at[i, "epcam_krt_mean_mean_log1p"] = np.nan
            per.at[i, "epcam_krt_mean_pct_pos"] = np.nan
            per.at[i, "epcam_krt_mean_ucell_module"] = np.nan
        else:
            per.at[i, "epcam_krt_mean_mean_log1p"] = float(sub["tacstd2_log1p_cp10k"].mean())
            per.at[i, "epcam_krt_mean_pct_pos"] = float((sub["TACSTD2"] > 0).mean())
            per.at[i, "epcam_krt_mean_ucell_module"] = float(sub["ucell_module"].mean())
    per.to_csv(HERE / "per_sample.tsv", sep="\t", index=False)

    post = per[per["is_post"]].copy()
    spear = pd.read_csv(HERE / "grid_spearman.tsv", sep="\t")
    nmpr = pd.read_csv(HERE / "grid_nmpr_mpr.tsv", sep="\t")
    spear = spear[spear["malignant_def"] != "epcam_krt_mean"].copy()
    nmpr = nmpr[nmpr["malignant_def"] != "epcam_krt_mean"].copy()

    mal = "epcam_krt_mean"
    new_nmpr = []
    new_spear = []
    for score in SCORE_DEFS:
        col = f"{mal}_{score}"
        use = post[np.isfinite(post[col])].copy()
        is_nmpr = use["response"].to_numpy() == "NMPR"
        w = exact_wilcoxon(use[col].to_numpy(), is_nmpr)
        new_nmpr.append(
            {
                "family": "NMPR_vs_MPR",
                "malignant_def": mal,
                "malignant_label": MAL_LABELS[mal],
                "score": score,
                "score_label": SCORE_LABELS[score],
                "n": int(len(use)),
                "n_NMPR": int(is_nmpr.sum()) if len(use) else 0,
                "n_MPR": int((~is_nmpr).sum()) if len(use) else 0,
                "mean_NMPR": w["mean_a"],
                "mean_MPR": w["mean_b"],
                "median_NMPR": w["median_a"],
                "median_MPR": w["median_b"],
                "delta_NMPR_minus_MPR": w["delta_a_minus_b"],
                "U": w["U"],
                "p": w["p_exact"],
                "p_onesided_NMPR_gt_MPR": w["p_onesided_A_gt_B"],
                "nmpr_gt_mpr": bool(w["delta_a_minus_b"] is not None and w["delta_a_minus_b"] > 0),
                "highlight": bool(w["delta_a_minus_b"] is not None and w["delta_a_minus_b"] > 0),
                "note": (
                    f"exact C({w['n_a']+w['n_b']},{w['n_a']})={w['n_perm']}; "
                    "pCR P06 counted as MPR; author CopyKAT IDs not public"
                ),
            }
        )
        for imm in IMMUNE_DEFS:
            tcol = f"frac_{imm}"
            u2 = use[np.isfinite(use[tcol])].copy()
            s = spearman_safe(u2[col], u2[tcol])
            rho = s["rho"]
            new_spear.append(
                {
                    "family": "vs_immune",
                    "malignant_def": mal,
                    "malignant_label": MAL_LABELS[mal],
                    "immune_def": imm,
                    "immune_label": IMM_LABELS[imm],
                    "score": score,
                    "score_label": SCORE_LABELS[score],
                    "n": s["n"],
                    "n_NMPR": int((u2["response"] == "NMPR").sum()) if len(u2) else 0,
                    "n_MPR": int((u2["response"] == "MPR").sum()) if len(u2) else 0,
                    "rho": rho,
                    "p": s["p"],
                    "highlight_rho_le_neg035": bool(rho is not None and rho <= RHO_HIGHLIGHT),
                    "note": "Spearman; unit=post-treatment patient; scipy.stats.spearmanr",
                }
            )
    spear = pd.concat([spear, pd.DataFrame(new_spear)], ignore_index=True)
    nmpr = pd.concat([nmpr, pd.DataFrame(new_nmpr)], ignore_index=True)
    # keep original malignant order
    spear["malignant_def"] = pd.Categorical(spear["malignant_def"], MAL_DEFS, ordered=True)
    spear["score"] = pd.Categorical(spear["score"], SCORE_DEFS, ordered=True)
    spear["immune_def"] = pd.Categorical(spear["immune_def"], IMMUNE_DEFS, ordered=True)
    spear = spear.sort_values(["malignant_def", "score", "immune_def"])
    nmpr["malignant_def"] = pd.Categorical(nmpr["malignant_def"], MAL_DEFS, ordered=True)
    nmpr["score"] = pd.Categorical(nmpr["score"], SCORE_DEFS, ordered=True)
    nmpr = nmpr.sort_values(["malignant_def", "score"])
    spear.to_csv(HERE / "grid_spearman.tsv", sep="\t", index=False)
    nmpr.to_csv(HERE / "grid_nmpr_mpr.tsv", sep="\t", index=False)
    hi_rho = spear[spear["highlight_rho_le_neg035"] == True]  # noqa: E712
    hi_nmpr = nmpr[nmpr["highlight"] == True]  # noqa: E712
    hi_rho.to_csv(HERE / "highlight_rho_le_neg035.tsv", sep="\t", index=False)
    hi_nmpr.to_csv(HERE / "highlight_nmpr_gt_mpr.tsv", sep="\t", index=False)

    summary = json.loads((HERE / "summary.json").read_text())
    summary["epcam_krt_mean_cut"] = ek_cut
    summary["malignant_n"]["epcam_krt_mean"] = int(is_ek_mean.sum())
    summary["n_highlight_rho_le_neg035"] = int(len(hi_rho))
    summary["n_highlight_nmpr_gt_mpr"] = int(len(hi_nmpr))
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    draw_grid(spear, nmpr, HERE / "fig_def_grid.png")
    print("updated ek_mean", int(is_ek_mean.sum()), "highlights", len(hi_rho), len(hi_nmpr), flush=True)


if __name__ == "__main__":
    main()

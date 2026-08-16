#!/usr/bin/env python3
"""Combine per-dataset outputs into cross-dataset synthesis tables.

Produces:
  SYNTHESIS_axis_matrix.tsv   dataset x gene_set: direction + competitive p + mean log2FC
  SYNTHESIS_cldn4.tsv         CLDN4 across all contrasts
  SYNTHESIS_perturbation_qc.tsv  TACSTD2/TROP2 knock-down verification per contrast
  SYNTHESIS_overview.json     machine-readable roll-up
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

RES = Path(__file__).resolve().parents[2] / "results" / "fable_tacstd2_kdko"

# contrast label / organism / perturbation for readability
META = {
    "GSE334497": ("Mouse 4T1 breast, Trop2 KO vs WT (tumour)", "mouse", "KO"),
    "GSE289287": ("Human T-47D breast, Trop-2 KO vs WT (xenograft)", "human", "KO"),
    "GSE289287_DSG2tumor": ("Human T-47D, DSG2 KO vs WT (xenograft) [partner]", "human", "KO(DSG2)"),
    "GSE289287_DSG2cell": ("Human T-47D, DSG2 KO vs WT (cells) [partner]", "human", "KO(DSG2)"),
    "GSE245459": ("Human ovarian, shTACSTD2 vs shNC", "human", "KD"),
    "GSE15212": ("Human colorectal RKO, siTACSTD2 vs neg-ctrl", "human", "KD"),
}
ORDER = ["GSE334497", "GSE289287", "GSE245459", "GSE15212",
         "GSE289287_DSG2tumor", "GSE289287_DSG2cell"]


def load_set(ds):
    f = RES / f"{ds}_setlevel.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else None


def load_pergene(ds):
    f = RES / f"{ds}_pergene_panel.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else None


# ---- axis matrix ----
rows = []
for ds in ORDER:
    sd = load_set(ds)
    if sd is None:
        continue
    label, org, pert = META[ds]
    for _, r in sd.iterrows():
        rows.append({
            "dataset": ds, "contrast": label, "organism": org, "perturbation": pert,
            "gene_set": r["gene_set"], "n_detected": r["n_detected"],
            "direction": r["comp_direction"], "comp_rank_biserial": r["comp_rank_biserial"],
            "comp_p": r["comp_p"], "mean_log2FC": r["sc_mean_log2FC"],
            "n_up": r["sc_n_up"], "n_down": r["sc_n_down"], "wilcoxon_p": r["sc_wilcoxon_p"],
        })
axis = pd.DataFrame(rows)
axis.to_csv(RES / "SYNTHESIS_axis_matrix.tsv", sep="\t", index=False, float_format="%.4g")

# compact direction grid for the 3 headline axes
def sig_flag(p):
    if pd.isna(p):
        return ""
    return "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else "ns"))

grid_sets = ["junctions", "tight_junction", "desmosome", "ifn_immune", "ifn_isg",
             "antigen_presentation", "tcell_cytotox"]
grid = {}
for ds in ORDER:
    sd = load_set(ds)
    if sd is None:
        continue
    sd = sd.set_index("gene_set")
    col = {}
    for gs in grid_sets:
        if gs in sd.index:
            d = sd.loc[gs, "comp_direction"]
            p = sd.loc[gs, "comp_p"]
            col[gs] = f"{d} {sig_flag(p)}"
        else:
            col[gs] = ""
    grid[META[ds][0]] = col
grid_df = pd.DataFrame(grid).T[grid_sets]
grid_df.to_csv(RES / "SYNTHESIS_direction_grid.tsv", sep="\t")

# ---- CLDN4 across datasets ----
cl = []
for ds in ORDER:
    pg = load_pergene(ds)
    if pg is None:
        continue
    row = pg[pg["gene"] == "CLDN4"]
    if len(row):
        r = row.iloc[0].to_dict()
        r["dataset"] = ds
        r["contrast"] = META[ds][0]
        cl.append(r)
cldn4 = pd.DataFrame(cl)
cols = [c for c in ["dataset", "contrast", "log2FC", "pvalue", "padj",
                    "t_pvalue", "t_fdr", "mwu_pvalue", "cohens_d"] if c in cldn4.columns]
cldn4 = cldn4[["dataset", "contrast"] + [c for c in cols if c not in ("dataset", "contrast")]]
cldn4.to_csv(RES / "SYNTHESIS_cldn4.tsv", sep="\t", index=False, float_format="%.4g")

# ---- perturbation QC ----
qc = []
for ds in ORDER:
    f = RES / f"{ds}_summary.json"
    if not f.exists():
        continue
    s = json.loads(f.read_text())
    q = s.get("perturbation_qc", {})
    entry = q.get("TACSTD2") or q.get("TROP2") or {}
    qc.append({
        "dataset": ds, "contrast": META[ds][0],
        "gene": "TACSTD2" if "TACSTD2" in q else ("TROP2" if "TROP2" in q else ""),
        "log2FC": entry.get("log2FC"),
        "pvalue": entry.get("pvalue", entry.get("t_pvalue")),
        "padj_or_fdr": entry.get("padj", entry.get("t_fdr")),
    })
qcdf = pd.DataFrame(qc)
qcdf.to_csv(RES / "SYNTHESIS_perturbation_qc.tsv", sep="\t", index=False, float_format="%.4g")

print("=== Perturbation QC (knockdown worked?) ===")
print(qcdf.to_string(index=False))
print("\n=== CLDN4 across datasets ===")
print(cldn4.to_string(index=False))
print("\n=== Direction grid (competitive test; *** p<1e-3, ** p<1e-2, * p<0.05) ===")
print(grid_df.to_string())

overview = {
    "axis_matrix": str(RES / "SYNTHESIS_axis_matrix.tsv"),
    "direction_grid": grid_df.to_dict(orient="index"),
    "cldn4": cldn4.to_dict(orient="records"),
    "perturbation_qc": qcdf.to_dict(orient="records"),
}
(RES / "SYNTHESIS_overview.json").write_text(json.dumps(overview, indent=2, default=str))
print("\n[written] SYNTHESIS_* tables in", RES)

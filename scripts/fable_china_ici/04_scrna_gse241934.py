#!/usr/bin/env python3
"""GSE241934 (NEOTIDE/CTONG2104, Peking University) scRNA-seq:
TACSTD2/CLDN4 in epithelial cells vs pathologic response (MPR incl. pCR vs non-MPR).

Cohorts (both post-neoadjuvant resected tumors, Chinese patients):
  IIT  = trial cohort, EGFR-mutant NSCLC, neoadjuvant sintilimab + chemo (n=11)
  Real = real-world cohort, EGFR-WT LUAD/ASC, neoadjuvant immunochemotherapy (n=34)

Sparse MatrixMarket matrices (genes x cells) are streamed with a filter on the
gene-row index; per-cell normalization uses nCount_RNA from the provided
per-cell metadata (which also contains major.cell.type == 'Epi' annotations).

Outputs -> results/fable_china_ici/{tables,figures}
"""

import gzip
import os
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.environ.get("DATA_DIR", "/tmp/fable_china_ici_data")
ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
TABLES = os.path.join(ROOT, "results", "fable_china_ici", "tables")
FIGS = os.path.join(ROOT, "results", "fable_china_ici", "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

GENES = ["TACSTD2", "CLDN4"]
COHORTS = {
    "IIT": dict(matrix="GSE241934_IIT_Matrix.mtx.gz", features="GSE241934_IIT_features.tsv.gz",
                barcodes="GSE241934_IIT_barcodes.tsv.gz", meta="GSE241934_IIT_Meta.txt.gz",
                label="IIT (EGFR-mut, sintilimab+chemo)"),
    "Real": dict(matrix="GSE241934_Real_Matrix.mtx.gz", features="GSE241934_RWC_features.tsv.gz",
                 barcodes="GSE241934_RWC_barcodes.tsv.gz", meta="GSE241934_Real_Meta.txt.gz",
                 label="Real-world (EGFR-WT, immunochemo)"),
}


def extract_gene_counts(matrix_path, gene_rows):
    """Stream a genes-x-cells .mtx.gz, return {gene: {cell_idx(1-based): count}}."""
    cond = " || ".join(f"$1=={i}" for i in gene_rows.values())
    awk = f"zcat {matrix_path} | tail -n +3 | awk '{cond} {{print $1, $2, $3}}'"
    out = subprocess.run(["bash", "-c", awk], capture_output=True, text=True, check=True).stdout
    row2gene = {v: k for k, v in gene_rows.items()}
    res = {g: {} for g in gene_rows}
    for line in out.splitlines():
        r, c, v = line.split()
        res[row2gene[int(r)]][int(c)] = int(v)
    return res


all_stats = []
patient_tables = []
for name, cfg in COHORTS.items():
    feats = pd.read_csv(f"{DATA}/{cfg['features']}", sep="\t", header=None)
    gene_rows = {}
    for g in GENES:
        idx = feats.index[feats[1] == g]
        assert len(idx) == 1, f"{g} not unique in features"
        gene_rows[g] = int(idx[0]) + 1  # 1-based mtx row
    print(f"[{name}] gene rows: {gene_rows}", flush=True)

    counts = extract_gene_counts(f"{DATA}/{cfg['matrix']}", gene_rows)
    barcodes = pd.read_csv(f"{DATA}/{cfg['barcodes']}", sep="\t", header=None)[0].values
    meta = pd.read_csv(f"{DATA}/{cfg['meta']}", sep="\t")
    meta = meta.set_index("cellID")

    cell_df = pd.DataFrame(index=barcodes)
    for g in GENES:
        arr = np.zeros(len(barcodes))
        for cidx, v in counts[g].items():
            arr[cidx - 1] = v
        cell_df[g] = arr
    cell_df = cell_df.join(meta[["sampleID", "major.cell.type", "nCount_RNA", "Pathological Response",
                                 "PD1", "EGFR", "Histology"]], how="inner")
    print(f"[{name}] cells matched to meta: {len(cell_df)} / {len(barcodes)}", flush=True)

    epi = cell_df[cell_df["major.cell.type"] == "Epi"].copy()
    for g in GENES:
        epi[f"{g}_log"] = np.log1p(epi[g] / epi["nCount_RNA"] * 1e4)
    agg = epi.groupby("sampleID").agg(
        n_epi_cells=("TACSTD2", "size"),
        response=("Pathological Response", "first"),
        pd1_drug=("PD1", "first"),
        egfr=("EGFR", "first"),
        histology=("Histology", "first"),
        TACSTD2_mean_log=("TACSTD2_log", "mean"),
        CLDN4_mean_log=("CLDN4_log", "mean"),
        TACSTD2_pct_pos=("TACSTD2", lambda x: 100 * (x > 0).mean()),
        CLDN4_pct_pos=("CLDN4", lambda x: 100 * (x > 0).mean()),
    ).reset_index()
    agg["cohort"] = name
    agg["mpr"] = np.where(agg["response"].isin(["MPR", "pCR"]), "MPR", "non-MPR")
    patient_tables.append(agg)

    d = agg[agg.n_epi_cells >= 20]
    for g in GENES:
        a = d.loc[d.mpr == "MPR", f"{g}_mean_log"]
        b = d.loc[d.mpr == "non-MPR", f"{g}_mean_log"]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        all_stats.append(dict(dataset=f"GSE241934 {cfg['label']} scRNA epithelial (post-tx)", gene=g,
                              unit="mean log1p CP10K per patient",
                              group_pos="MPR(incl pCR)", n_pos=len(a), median_pos=round(a.median(), 3),
                              group_neg="non-MPR", n_neg=len(b), median_neg=round(b.median(), 3),
                              mannwhitney_p=round(p, 4), auc_pos_high=round(u / (len(a) * len(b)), 3),
                              note="author 'Epi' annotation, patients with >=20 epithelial cells"))

pt = pd.concat(patient_tables, ignore_index=True)
pt.to_csv(f"{TABLES}/gse241934_scrna_epithelial_by_patient.csv", index=False)
sdf = pd.DataFrame(all_stats)
sdf.to_csv(f"{TABLES}/gse241934_scrna_stats.csv", index=False)
print(sdf.to_string(index=False))

fig, axes = plt.subplots(2, 2, figsize=(9, 7.5))
for i, (name, cfg) in enumerate(COHORTS.items()):
    d = pt[(pt.cohort == name) & (pt.n_epi_cells >= 20)]
    for j, g in enumerate(GENES):
        ax = axes[j, i]
        order = ["MPR", "non-MPR"]
        vals = [d.loc[d.mpr == grp, f"{g}_mean_log"].values for grp in order]
        bp = ax.boxplot(vals, tick_labels=[f"{grp}\n(n={len(v)})" for grp, v in zip(order, vals)],
                        widths=0.55, showfliers=False, patch_artist=True)
        for patch, color in zip(bp["boxes"], ["#f6b0a0", "#a8c8e8"]):
            patch.set_facecolor(color)
        for k, v in enumerate(vals):
            ax.scatter(np.random.default_rng(3).normal(k + 1, 0.05, len(v)), v, s=22,
                       color="#333", alpha=0.8, zorder=3)
        row = sdf[(sdf.dataset.str.contains(name if name == "IIT" else "Real")) & (sdf.gene == g)]
        ax.set_title(f"{g} (MW p={row.mannwhitney_p.iloc[0]:.3f})", fontsize=10)
        if j == 0:
            ax.text(0.5, 1.18, cfg["label"], ha="center", va="bottom", transform=ax.transAxes, fontsize=9)
        ax.set_ylabel("mean log1p CP10K (Epi cells)" if i == 0 else "")
fig.suptitle("GSE241934 scRNA: epithelial TACSTD2/CLDN4 in post-treatment tumors, MPR vs non-MPR", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{FIGS}/gse241934_scrna_epithelial.png", dpi=200)
print(f"done -> {FIGS}/gse241934_scrna_epithelial.png")

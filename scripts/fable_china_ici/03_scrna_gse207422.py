#!/usr/bin/env python3
"""GSE207422 scRNA-seq: TACSTD2/CLDN4 in epithelial cells vs pathologic response.

The 175 MB dense UMI matrix (genes x 92,330 cells) is streamed once:
  - per-cell total UMI (for CP10K normalization)
  - rows for TACSTD2, CLDN4, EPCAM, PTPRC, KRT8, KRT18
No per-cell cell-type annotation is provided on GEO, so epithelial cells are
proxied as EPCAM+ (UMI > 0) & PTPRC == 0. Sample-level pathologic response
(MPR incl. pCR vs NMPR) comes from GSE207422_NSCLC_scRNAseq_metadata.xlsx.
Note: 15/20 samples are post-treatment surgical specimens; pre-treatment
biopsies are analyzed separately (small n -> descriptive only).

Outputs -> results/fable_china_ici/{tables,figures}
"""

import gzip
import os

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

WANT = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC", "KRT8", "KRT18"]
MATRIX = f"{DATA}/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"

rows = {}
with gzip.open(MATRIX, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    cells = np.array(header[1:])
    colsum = np.zeros(len(cells))
    for i, line in enumerate(fh):
        tab = line.index("\t")
        vals = np.fromstring(line[tab + 1:], sep="\t")
        colsum += vals
        gene = line[:tab]
        if gene in WANT:
            rows[gene] = vals
        if i % 5000 == 0:
            print(f"  streamed {i} genes...", flush=True)
print(f"cells={len(cells)}, captured={sorted(rows)}")

df = pd.DataFrame({g: rows[g] for g in WANT}, index=cells)
df["total_umi"] = colsum
df["sample"] = ["_".join(c.split("_")[:-1]) for c in cells]
for g in ["TACSTD2", "CLDN4"]:
    df[f"{g}_cp10k_log"] = np.log1p(df[g] / df["total_umi"] * 1e4)
df["is_epi"] = (df["EPCAM"] > 0) & (df["PTPRC"] == 0)

meta = pd.read_excel(f"{DATA}/GSE207422_NSCLC_scRNAseq_metadata.xlsx").dropna(subset=["Patient"])
meta["mpr"] = np.where(meta["Pathologic Response"].astype(str).str.startswith(("MPR", "pCR")), "MPR",
                       np.where(meta["Pathologic Response"].astype(str) == "NE", "NE", "NMPR"))
meta["timepoint"] = np.where(meta["Resource"].str.contains("Pre", case=False), "pre", "post")

epi = df[df["is_epi"]]
agg = epi.groupby("sample").agg(
    n_epi_cells=("TACSTD2", "size"),
    TACSTD2_mean_cp10k_log=("TACSTD2_cp10k_log", "mean"),
    CLDN4_mean_cp10k_log=("CLDN4_cp10k_log", "mean"),
    TACSTD2_pct_pos=("TACSTD2", lambda x: 100 * (x > 0).mean()),
    CLDN4_pct_pos=("CLDN4", lambda x: 100 * (x > 0).mean()),
).reset_index()
agg = agg.merge(meta[["Sample", "Patient", "timepoint", "mpr", "Pathologic Response", "RECIST"]],
                left_on="sample", right_on="Sample").drop(columns="Sample")
agg.to_csv(f"{TABLES}/gse207422_scrna_epithelial_by_sample.csv", index=False)
print(agg.to_string(index=False))

stats_rows = []
for tp in ["post", "pre"]:
    d = agg[(agg.timepoint == tp) & (agg.mpr != "NE") & (agg.n_epi_cells >= 20)]
    for g in ["TACSTD2", "CLDN4"]:
        a = d.loc[d.mpr == "MPR", f"{g}_mean_cp10k_log"]
        b = d.loc[d.mpr == "NMPR", f"{g}_mean_cp10k_log"]
        if len(a) >= 3 and len(b) >= 3:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            auc = u / (len(a) * len(b))
            stats_rows.append(dict(dataset=f"GSE207422 scRNA epithelial ({tp}-treatment)", gene=g,
                                   unit="mean log1p CP10K per sample",
                                   group_pos="MPR", n_pos=len(a), median_pos=round(a.median(), 3),
                                   group_neg="NMPR", n_neg=len(b), median_neg=round(b.median(), 3),
                                   mannwhitney_p=round(p, 4), auc_pos_high=round(auc, 3),
                                   note="EPCAM+/PTPRC- proxy, samples with >=20 epithelial cells"))
        else:
            stats_rows.append(dict(dataset=f"GSE207422 scRNA epithelial ({tp}-treatment)", gene=g,
                                   unit="mean log1p CP10K per sample",
                                   group_pos="MPR", n_pos=len(a), median_pos=np.nan,
                                   group_neg="NMPR", n_neg=len(b), median_neg=np.nan,
                                   mannwhitney_p=np.nan, auc_pos_high=np.nan,
                                   note="too few samples for testing (descriptive only)"))
sdf = pd.DataFrame(stats_rows)
sdf.to_csv(f"{TABLES}/gse207422_scrna_stats.csv", index=False)
print(sdf.to_string(index=False))

d = agg[(agg.timepoint == "post") & (agg.mpr != "NE") & (agg.n_epi_cells >= 20)]
fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))
for ax, g in zip(axes, ["TACSTD2", "CLDN4"]):
    order = ["MPR", "NMPR"]
    vals = [d.loc[d.mpr == grp, f"{g}_mean_cp10k_log"].values for grp in order]
    bp = ax.boxplot(vals, tick_labels=[f"{grp}\n(n={len(v)})" for grp, v in zip(order, vals)],
                    widths=0.55, showfliers=False, patch_artist=True)
    for patch, color in zip(bp["boxes"], ["#f6b0a0", "#a8c8e8"]):
        patch.set_facecolor(color)
    for k, v in enumerate(vals):
        ax.scatter(np.random.default_rng(2).normal(k + 1, 0.05, len(v)), v, s=22, color="#333", alpha=0.8, zorder=3)
    row = sdf[(sdf.dataset.str.contains("post")) & (sdf.gene == g)]
    ax.set_title(f"{g} (MW p={row.mannwhitney_p.iloc[0]:.3f})", fontsize=10)
    ax.set_ylabel("mean log1p CP10K (epithelial)")
fig.suptitle("GSE207422 scRNA: epithelial TACSTD2/CLDN4, post-treatment, MPR vs NMPR", fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIGS}/gse207422_scrna_epithelial.png", dpi=200)
print(f"done -> {FIGS}/gse207422_scrna_epithelial.png")

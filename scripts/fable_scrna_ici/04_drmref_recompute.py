#!/usr/bin/env python3
"""Recompute GSE207422 malignant-only TACSTD2 using DRMref cell labels
(the public annotation that produced the claimed ρ ≈ −0.4 to −0.5).

DRMref (Liu et al. NAR 2024) hosts marker-based Seurat objects
GSE207422_Tor + GSE207422_Sin: 30,877 post-treatment cells, 2,051
labeled 'Malignant cells'. This is NOT Hu et al. CopyKAT.

TACSTD2 UMIs and per-cell totals are taken from our GEO panel extract
(independent of the A3 annotation file's precomputed TACSTD2 columns).
T/NK fraction uses DRMref types: CD8+ T + CD4+ T + NK.

Outputs under results/fable_scrna_ici/.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "/tmp/data_scrna_ici"
RES = "/workspace/results/fable_scrna_ici"
ANN_SRC = "/tmp/a3_annot/drmref_cell_annotation.tsv.gz"
ANN_DST = f"{RES}/drmref_cell_annotation.tsv.gz"
os.makedirs(RES, exist_ok=True)

# keep a local copy so this slice is self-contained
import shutil
if os.path.abspath(ANN_SRC) != os.path.abspath(ANN_DST):
    shutil.copy2(ANN_SRC, ANN_DST)

ann = pd.read_csv(ANN_DST, sep="\t")
panel = pd.read_csv(f"{DATA}/gse207422_panel_cells.tsv.gz", sep="\t")
panel = panel.rename(columns={"barcode": "cell_barcode"})
keep_genes = [c for c in ["TACSTD2", "CLDN4", "GZMB", "GZMA", "GZMK",
                           "PRF1", "IFNG", "NKG7", "PDCD1", "HAVCR2",
                           "LAG3", "TIGIT", "TOX", "CD8A", "CD8B",
                           "CD3D", "CD3E"] if c in panel.columns]
m = ann.merge(panel[["cell_barcode", "total_umi"] + keep_genes],
              on="cell_barcode", how="inner")
print(f"DRMref cells {len(ann)}; joined to GEO panel {len(m)}")
assert len(m) == len(ann), "barcode mismatch vs DRMref"

def log1p_cp10k(umi, total):
    return np.log1p(np.asarray(umi, float) / np.asarray(total, float) * 1e4)

m["tacstd2_ln"] = log1p_cp10k(m["TACSTD2"], m["total_umi"])
m["cldn4_ln"] = log1p_cp10k(m["CLDN4"], m["total_umi"])
cyto = [g for g in ["GZMB", "GZMA", "PRF1", "IFNG", "NKG7"] if g in m.columns]
exh = [g for g in ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"] if g in m.columns]
m["cyto_ln"] = np.mean([log1p_cp10k(m[g], m["total_umi"]) for g in cyto], axis=0)
m["exh_ln"] = np.mean([log1p_cp10k(m[g], m["total_umi"]) for g in exh], axis=0)

TNK_TYPES = {"CD8+ T cells", "CD4+ T cells", "NK cells"}
m["is_mal"] = m["celltype"] == "Malignant cells"
m["is_tnk"] = m["celltype"].isin(TNK_TYPES)
m["response_bin"] = m["response"].map(
    {"MPR": "MPR", "MPR (pCR)": "MPR", "NMPR": "NMPR"})

rows = []
for samp, g in m.groupby("orig.ident"):
    mal = g[g.is_mal]
    tnk = g[g.is_tnk]
    rec = g.iloc[0]
    rows.append(dict(
        sample=samp, patient=rec.patient, cohort=rec.cohort,
        response_raw=rec.response, response=rec.response_bin,
        n_annot=len(g), n_malignant=int(g.is_mal.sum()),
        n_tnk=int(g.is_tnk.sum()),
        frac_tnk=g.is_tnk.mean(),
        malig_TACSTD2_mean=mal["tacstd2_ln"].mean() if len(mal) else np.nan,
        malig_TACSTD2_pct=(mal["TACSTD2"] > 0).mean() if len(mal) else np.nan,
        malig_CLDN4_mean=mal["cldn4_ln"].mean() if len(mal) else np.nan,
        tnk_cyto_mean=tnk["cyto_ln"].mean() if len(tnk) else np.nan,
        tnk_exh_mean=tnk["exh_ln"].mean() if len(tnk) else np.nan,
    ))
samp = pd.DataFrame(rows).sort_values("patient")
samp.to_csv(f"{RES}/gse207422_drmref_sample_metrics.tsv", sep="\t", index=False)
print(samp.round(4).to_string())

def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    return stats.mannwhitneyu(a, b, alternative="two-sided")

def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    r, p = stats.spearmanr(x[mask], y[mask])
    return float(r), float(p), int(mask.sum())

stats_rows = []
nmpr = samp.loc[samp.response == "NMPR", "malig_TACSTD2_mean"]
mpr = samp.loc[samp.response == "MPR", "malig_TACSTD2_mean"]
u, p = mwu(nmpr, mpr)
stats_rows.append(dict(
    dataset="GSE207422_DRMref", analysis="post_NMPR_vs_MPR",
    comparison="NMPR vs MPR", metric="malig_TACSTD2_mean",
    n=f"{len(nmpr)}v{len(mpr)}", stat="MannWhitneyU", value=u, p_value=p,
    note=(f"median NMPR={nmpr.median():.4f} MPR={mpr.median():.4f}; "
          f"mean NMPR={nmpr.mean():.4f} MPR={mpr.mean():.4f}; "
          f"DRMref malignant labels; pCR=MPR; all 12 post-tx samples")))

for yv in ["frac_tnk", "tnk_cyto_mean", "tnk_exh_mean"]:
    r, p, n = spear(samp["malig_TACSTD2_mean"], samp[yv])
    stats_rows.append(dict(
        dataset="GSE207422_DRMref", analysis="spearman_post12",
        comparison=f"malig_TACSTD2_mean~{yv}", metric=yv,
        n=n, stat="spearman_rho", value=r, p_value=p,
        note="DRMref labels; 12 post-tx patients"))

r, p, n = spear(samp["malig_TACSTD2_pct"], samp["frac_tnk"])
stats_rows.append(dict(
    dataset="GSE207422_DRMref", analysis="spearman_post12_secondary",
    comparison="malig_TACSTD2_pct~frac_tnk", metric="frac_tnk",
    n=n, stat="spearman_rho", value=r, p_value=p,
    note="secondary: % TACSTD2+ malignant vs T/NK"))

# drop P06 (<20 malignant) as the A3 sensitivity
sub = samp[samp.n_malignant >= 20]
nmpr2 = sub.loc[sub.response == "NMPR", "malig_TACSTD2_mean"]
mpr2 = sub.loc[sub.response == "MPR", "malig_TACSTD2_mean"]
u, p = mwu(nmpr2, mpr2)
stats_rows.append(dict(
    dataset="GSE207422_DRMref", analysis="post_NMPR_vs_MPR_min20mal",
    comparison="NMPR vs MPR", metric="malig_TACSTD2_mean",
    n=f"{len(nmpr2)}v{len(mpr2)}", stat="MannWhitneyU", value=u, p_value=p,
    note=(f"drops P06 (n_mal={int(samp.loc[samp.patient=='P06','n_malignant'].iloc[0])}); "
          f"mean NMPR={nmpr2.mean():.4f} MPR={mpr2.mean():.4f}; POST-HOC")))
r, p, n = spear(sub["malig_TACSTD2_mean"], sub["frac_tnk"])
stats_rows.append(dict(
    dataset="GSE207422_DRMref", analysis="spearman_min20mal",
    comparison="malig_TACSTD2_mean~frac_tnk", metric="frac_tnk",
    n=n, stat="spearman_rho", value=r, p_value=p,
    note="drops P06; POST-HOC"))

out = pd.DataFrame(stats_rows)
out.to_csv(f"{RES}/gse207422_drmref_stats.tsv", sep="\t", index=False)
print(out.to_string())

# figure
rng = np.random.default_rng(0)
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
# box
ax = axes[0]
data = [mpr.to_numpy(), nmpr.to_numpy()]
ax.boxplot(data, showfliers=False)
ax.set_xticklabels([f"MPR\n(n={len(mpr)})", f"NMPR\n(n={len(nmpr)})"])
for i, d in enumerate(data):
    ax.scatter(rng.normal(i + 1, 0.06, len(d)), d, s=32, zorder=3)
ax.set_ylabel("malignant TACSTD2 mean (log1p CP10K)")
ax.set_title(f"DRMref labels, post-tx\nNMPR mean {nmpr.mean():.2f} vs MPR {mpr.mean():.2f}")
# scatter
ax = axes[1]
ax.scatter(samp["malig_TACSTD2_mean"], samp["frac_tnk"], s=40)
for _, r0 in samp.iterrows():
    ax.annotate(r0.patient, (r0.malig_TACSTD2_mean, r0.frac_tnk),
                fontsize=7, xytext=(3, 3), textcoords="offset points")
r, p, n = spear(samp["malig_TACSTD2_mean"], samp["frac_tnk"])
ax.set_xlabel("malignant TACSTD2 mean (log1p CP10K)")
ax.set_ylabel("T/NK fraction (DRMref CD8+CD4+NK)")
ax.set_title(f"DRMref 12 post-tx patients\nSpearman ρ={r:.3f}  p={p:.3g}  n={n}")
fig.suptitle("GSE207422 DRMref malignant TACSTD2 (not CopyKAT)")
fig.tight_layout()
fig.savefig(f"{RES}/fig3_drmref_tacstd2.png", dpi=150)
print("wrote fig3_drmref_tacstd2.png")

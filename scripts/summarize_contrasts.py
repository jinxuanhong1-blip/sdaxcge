#!/usr/bin/env python3
"""Pooled histology contrasts and summary figures from section_metrics.csv."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

RES = Path("/workspace/results")
MAPS = Path("/workspace/maps")
df = pd.read_csv(RES / "section_metrics.csv")


def stars(p):
    if not np.isfinite(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


rows = []

# paired lepidic vs acinar within GSE273378 sections that have both
g = df[df.dataset == "GSE273378"]
lep = g[g.subset == "lepidic"].set_index("section")
aci = g[g.subset == "acinar"].set_index("section")
both = lep.index.intersection(aci.index)
print("paired lepidic+acinar sections:", list(both))

metrics = [
    ("cldn4_cd8a_spearman_epi", "CLDN4 vs CD8A Spearman (epithelial)"),
    ("q4_minus_q1_nn_cd8high_um", "Q4-Q1 nearest CD8A-high distance (um)"),
    ("q4_minus_q1_krt8resid_neighbor_cd8a", "Q4-Q1 KRT8-residual neighbor CD8A"),
    ("mean_cldn4_subset", "mean CLDN4"),
    ("mean_cd8a_subset", "mean CD8A"),
]
for col, label in metrics:
    a = lep.loc[both, col].astype(float)
    b = aci.loc[both, col].astype(float)
    m = a.notna() & b.notna()
    a, b = a[m], b[m]
    if len(a) >= 3:
        w = stats.wilcoxon(a.values, b.values, alternative="two-sided")
        p = float(w.pvalue)
    else:
        p = np.nan
    rows.append(
        dict(
            contrast="lepidic_vs_acinar_paired_GSE273378",
            metric=col,
            n=int(len(a)),
            lepidic_median=float(a.median()) if len(a) else np.nan,
            acinar_median=float(b.median()) if len(b) else np.nan,
            delta_lep_minus_aci=float((a - b).median()) if len(a) else np.nan,
            p=p,
        )
    )
    print(label, "n", len(a), "lep", a.median() if len(a) else None, "aci", b.median() if len(b) else None, "p", p)

# STAS-containing vs non-STAS sections (epithelial metrics)
e = g[g.subset == "epithelial"].copy()
e["n_stas"] = e.section_histology.str.extract(r"STAS=(\d+)").astype(float)
e["has_stas"] = e["n_stas"] > 0
print("\nSTAS sections", e.loc[e.has_stas, "section"].tolist())
print("nonSTAS sections", e.loc[~e.has_stas, "section"].tolist())
for col, label in metrics:
    a = e.loc[e.has_stas, col].astype(float)
    b = e.loc[~e.has_stas, col].astype(float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) >= 3 and len(b) >= 3:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        p = float(p)
    else:
        p = np.nan
    rows.append(
        dict(
            contrast="STAS_section_vs_nonSTAS_GSE273378",
            metric=col,
            n=int(len(a) + len(b)),
            lepidic_median=float(a.median()) if len(a) else np.nan,  # here = STAS
            acinar_median=float(b.median()) if len(b) else np.nan,  # here = nonSTAS
            delta_lep_minus_aci=float(a.median() - b.median()) if len(a) and len(b) else np.nan,
            p=p,
        )
    )
    print("STAS vs non", label, "STAS med", a.median() if len(a) else None, "non med", b.median() if len(b) else None, "p", p)

# GSE189487 AIS vs IAC
e2 = df[(df.dataset == "GSE189487") & (df.subset == "epithelial")].copy()
for col, label in metrics[:3]:
    a = e2.loc[e2.section_histology == "AIS", col].astype(float)
    b = e2.loc[e2.section_histology == "IAC", col].astype(float)
    print("AIS vs IAC", col, "AIS", list(a), "IAC", list(b))

# sign tests of Spearman < 0
print("\nSign tests Spearman epithelial < 0")
for ds, sub in df[df.subset == "epithelial"].groupby("dataset"):
    x = sub.cldn4_cd8a_spearman_all.dropna()
    n_neg = int((x < 0).sum())
    p = float(stats.binomtest(n_neg, len(x), 0.5, alternative="greater").pvalue) if len(x) else np.nan
    print(ds, f"{n_neg}/{len(x)} negative, one-sided binom p={p:.3g}")
    rows.append(
        dict(
            contrast=f"sign_test_spearman_neg_{ds}",
            metric="cldn4_cd8a_spearman_all",
            n=int(len(x)),
            lepidic_median=float(x.median()) if len(x) else np.nan,
            acinar_median=np.nan,
            delta_lep_minus_aci=n_neg / len(x) if len(x) else np.nan,
            p=p,
        )
    )

out = pd.DataFrame(rows)
out.to_csv(RES / "pooled_contrasts.csv", index=False)

# figures
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), dpi=150)

# 1 paired spearman
ax = axes[0]
a = lep.loc[both, "cldn4_cd8a_spearman_epi"].astype(float)
b = aci.loc[both, "cldn4_cd8a_spearman_epi"].astype(float)
for s in both:
    ax.plot([0, 1], [a.loc[s], b.loc[s]], color="0.7", lw=1)
ax.scatter(np.zeros(len(a)), a, c="#2ca02c", s=40, zorder=3, label="Lepidic")
ax.scatter(np.ones(len(b)), b, c="#d62728", s=40, zorder=3, label="Acinar")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks([0, 1], ["Lepidic", "Acinar"])
ax.set_ylabel("CLDN4–CD8A Spearman")
ax.set_title("GSE273378 paired sections\n(pathologist labels)")
ax.legend(frameon=False, fontsize=8)

# 2 Q4-Q1 nn distance
ax = axes[1]
a = lep.loc[both, "q4_minus_q1_nn_cd8high_um"].astype(float)
b = aci.loc[both, "q4_minus_q1_nn_cd8high_um"].astype(float)
for s in both:
    ax.plot([0, 1], [a.loc[s], b.loc[s]], color="0.7", lw=1)
ax.scatter(np.zeros(len(a)), a, c="#2ca02c", s=40, zorder=3)
ax.scatter(np.ones(len(b)), b, c="#d62728", s=40, zorder=3)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks([0, 1], ["Lepidic", "Acinar"])
ax.set_ylabel("Q4−Q1 nearest CD8A-high (µm)")
ax.set_title("CLDN4-high epithelial spots:\ndistance to CD8A-high")

# 3 STAS vs non STAS spearman
ax = axes[2]
rng = np.random.default_rng(0)
x0 = rng.normal(0, 0.04, size=e.has_stas.sum())
x1 = rng.normal(1, 0.04, size=(~e.has_stas).sum())
ax.scatter(x0, e.loc[e.has_stas, "cldn4_cd8a_spearman_epi"], c="#9467bd", s=40, label="STAS present")
ax.scatter(x1, e.loc[~e.has_stas, "cldn4_cd8a_spearman_epi"], c="#8c564b", s=40, label="No STAS")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks([0, 1], ["STAS+", "STAS−"])
ax.set_ylabel("CLDN4–CD8A Spearman (epithelial)")
ax.set_title("GSE273378 sections\nSTAS vs no STAS")
ax.legend(frameon=False, fontsize=8)

fig.tight_layout()
fig.savefig(MAPS / "summary_lepidic_acinar_stas.png")
plt.close(fig)

# dataset overview forest of epithelial spearman
fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
eall = df[df.subset == "epithelial"].copy()
eall = eall.sort_values(["dataset", "section"])
y = np.arange(len(eall))
colors = {
    "GSE273378": "#1f77b4",
    "GSE300676": "#ff7f0e",
    "GSE189487": "#2ca02c",
    "KERO_AdSpatial2024": "#d62728",
    "VisiumHD_10x_LUAD_STAS": "#9467bd",
}
ax.axvline(0, color="k", lw=0.7)
for i, r in enumerate(eall.itertuples()):
    ax.plot(r.cldn4_cd8a_spearman_all, i, "o", color=colors.get(r.dataset, "k"), ms=5)
ax.set_yticks(y)
ax.set_yticklabels([f"{r.dataset.split('_')[0]} | {r.section[-18:]} | {str(r.section_histology)[:22]}" for r in eall.itertuples()], fontsize=6)
ax.set_xlabel("Section-wide CLDN4 vs CD8A Spearman")
ax.set_title("CLDN4-only vs CD8A (all QC spots)")
fig.tight_layout()
fig.savefig(MAPS / "summary_spearman_all_sections.png")
plt.close(fig)

print("wrote pooled contrasts and summary maps")

#!/usr/bin/env python3
"""CLAIM A1 (LUAD-only): TACSTD2 (TROP2) expression vs immune features in
TCGA-LUAD primary tumors, adjusting for ABSOLUTE tumor purity.

Inputs (downloaded, see provenance.json for URLs/md5):
  data/TCGA.LUAD.HiSeqV2.gz                        Xena TCGA-LUAD RNA-seq, log2(RSEM norm_count + 1)
  data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt PanCanAtlas ABSOLUTE purity/ploidy
  data/TCGA_all_leuk_estimate.masked.20170107.tsv   PanImmune methylation leukocyte fraction
  data/TCGA.Kallisto.fullIDs.cibersort.relative.tsv PanImmune CIBERSORT relative fractions
  data/Scores_160_Signatures.tsv.gz                 PanImmune expression signature scores

Outputs -> results/w200/A1_LUAD/
"""

import gzip
import hashlib
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "w200", "A1_LUAD")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

RNG_NOTE = "deterministic (no sampling)"

# Immune marker genes measured on the same platform as TACSTD2
MARKER_GENES = [
    "CD8A", "CD3E", "CD2", "GZMA", "GZMB", "PRF1", "NKG7", "IFNG",
    "CXCL9", "CXCL10", "CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT",
    "HAVCR2", "IDO1", "FOXP3", "MS4A1", "CD68",
]

WOLF_SIGNATURES = {
    "LIexpression_score": "Wolf_LymphocyteInfiltration",
    "IFNG_score_21050467": "Wolf_IFNgamma",
    "TGFB_score_21050467": "Wolf_TGFbeta",
    "CSF1_response": "Wolf_MacrophageCSF1",
    "CHANG_CORE_SERUM_RESPONSE_UP": "Wolf_WoundHealing_CSR",
}

CIBERSORT_CELLS = [
    "T.cells.CD8", "T.cells.regulatory..Tregs.", "NK.cells.activated",
    "Macrophages.M1", "Macrophages.M2",
]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sample01(barcode, sep="-"):
    """Aliquot barcode -> patient-level primary-tumor sample id, or None."""
    parts = barcode.split(sep)
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3]) + "-01"


def partial_spearman(x, y, z):
    """Spearman correlation of x,y controlling for z (Pearson on ranks of residuals)."""
    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    zr = stats.rankdata(z)
    zc = np.column_stack([np.ones_like(zr), zr])
    bx, *_ = np.linalg.lstsq(zc, xr, rcond=None)
    by, *_ = np.linalg.lstsq(zc, yr, rcond=None)
    rx = xr - zc @ bx
    ry = yr - zc @ by
    r, _ = stats.pearsonr(rx, ry)
    n = len(x)
    # t-test with n - 2 - (number of covariates) df
    df = n - 3
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    p = 2 * stats.t.sf(abs(t), df)
    return r, p


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, p[idx] * n / rank)
        q[idx] = val
        prev = val
    return q


print("Loading expression ...")
expr = pd.read_csv(os.path.join(DATA, "TCGA.LUAD.HiSeqV2.gz"), sep="\t", index_col=0)
# columns already patient-level TCGA-XX-XXXX-NN; keep primary tumors (-01)
expr = expr.loc[:, [c for c in expr.columns if c.endswith("-01")]]
genes_needed = ["TACSTD2"] + MARKER_GENES
missing = [g for g in genes_needed if g not in expr.index]
if missing:
    raise SystemExit(f"missing genes in expression matrix: {missing}")
eg = expr.loc[genes_needed].T  # samples x genes, log2(RSEM+1)
eg["CYT_GZMA_PRF1"] = eg[["GZMA", "PRF1"]].mean(axis=1)

print("Loading ABSOLUTE purity ...")
ab = pd.read_csv(os.path.join(DATA, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), sep="\t")
ab = ab[ab["array"].str.endswith("-01", na=False)][["array", "purity"]].dropna()
purity = ab.groupby("array")["purity"].mean()

print("Loading leukocyte fraction ...")
lf = pd.read_csv(os.path.join(DATA, "TCGA_all_leuk_estimate.masked.20170107.tsv"),
                 sep="\t", header=None, names=["cancer", "aliquot", "leuk"])
lf = lf[lf["cancer"] == "LUAD"].copy()
lf["sample"] = lf["aliquot"].map(sample01)
lf = lf.dropna(subset=["sample"]).groupby("sample")["leuk"].mean()

print("Loading CIBERSORT ...")
cs = pd.read_csv(os.path.join(DATA, "TCGA.Kallisto.fullIDs.cibersort.relative.tsv"), sep="\t")
cs = cs[cs["CancerType"] == "LUAD"].copy()
cs["sample"] = cs["SampleID"].map(lambda b: sample01(b, sep="."))
cs = cs.dropna(subset=["sample"]).groupby("sample")[CIBERSORT_CELLS].mean()
cs.columns = ["CIBERSORT_" + c.replace("T.cells.regulatory..Tregs.", "Tregs")
              .replace("T.cells.CD8", "CD8_T").replace("NK.cells.activated", "NK_activated")
              .replace("Macrophages.M1", "Macrophage_M1").replace("Macrophages.M2", "Macrophage_M2")
              for c in cs.columns]

print("Loading Wolf signature scores ...")
sig = pd.read_csv(os.path.join(DATA, "Scores_160_Signatures.tsv.gz"), sep="\t")
sig = sig[sig["SetName"].isin(WOLF_SIGNATURES)].set_index("SetName")
sig = sig.drop(columns=["Source"]).T
sig.index = [sample01(b) or "" for b in sig.index]
sig = sig[sig.index != ""]
sig = sig.loc[sig.index.isin(eg.index)]  # restrict pan-cancer table to LUAD samples
sig = sig.astype(float).groupby(level=0).mean()
sig.columns = [WOLF_SIGNATURES[c] for c in sig.columns]

print("Merging ...")
df = eg.join(purity.rename("purity"), how="inner")
df = df.join(lf.rename("LeukocyteFraction_methylation"), how="left")
df = df.join(cs, how="left").join(sig, how="left")
df.index.name = "sample"

immune_features = (
    ["LeukocyteFraction_methylation"]
    + list(sig.columns)
    + list(cs.columns)
    + ["CYT_GZMA_PRF1"]
    + MARKER_GENES
)

print(f"n primary tumors with expression: {eg.shape[0]}")
print(f"n with expression + purity (analysis set): {df.shape[0]}")

rows = []
x_all = df["TACSTD2"].values
z_all = df["purity"].values
for feat in immune_features:
    mask = df[feat].notna()
    x, y, z = x_all[mask], df.loc[mask, feat].values, z_all[mask]
    n = int(mask.sum())
    rho, p = stats.spearmanr(x, y)
    prho, pp = partial_spearman(x, y, z)
    fclass = ("gene_marker" if feat in MARKER_GENES + ["CYT_GZMA_PRF1"]
              else "CIBERSORT" if feat.startswith("CIBERSORT")
              else "Wolf_signature" if feat.startswith("Wolf")
              else "leukocyte_fraction")
    rows.append(dict(feature=feat, feature_class=fclass, n=n,
                     spearman_rho=rho, spearman_p=p,
                     partial_rho_purity_adj=prho, partial_p=pp))

res = pd.DataFrame(rows)
res["spearman_fdr"] = bh_fdr(res["spearman_p"])
res["partial_fdr"] = bh_fdr(res["partial_p"])
res["sig_after_purity_adj_fdr05"] = res["partial_fdr"] < 0.05
res = res.sort_values("partial_p").reset_index(drop=True)

# context stats
rho_purity, p_purity = stats.spearmanr(df["TACSTD2"], df["purity"])
lf_mask = df["LeukocyteFraction_methylation"].notna()
rho_lf_pur, _ = stats.spearmanr(df.loc[lf_mask, "purity"],
                                df.loc[lf_mask, "LeukocyteFraction_methylation"])

res.to_csv(os.path.join(OUT, "a1_luad_correlations.tsv"), sep="\t", index=False)
df.reset_index().to_csv(os.path.join(OUT, "a1_luad_analysis_table.tsv"), sep="\t", index=False)

# ---------------- figures ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plot_feats = ["LeukocyteFraction_methylation", "Wolf_LymphocyteInfiltration",
              "CYT_GZMA_PRF1", "CD8A", "Wolf_IFNgamma", "CD274"]
fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
for ax, feat in zip(axes.ravel(), plot_feats):
    mask = df[feat].notna()
    sc = ax.scatter(df.loc[mask, "TACSTD2"], df.loc[mask, feat],
                    c=df.loc[mask, "purity"], cmap="viridis", s=10, alpha=0.7)
    r = res.set_index("feature").loc[feat]
    ax.set_title(f"{feat}\nSpearman rho={r.spearman_rho:.2f}; "
                 f"purity-adj rho={r.partial_rho_purity_adj:.2f} "
                 f"(FDR={r.partial_fdr:.1e}, n={int(r.n)})", fontsize=9)
    ax.set_xlabel("TACSTD2 log2(RSEM+1)")
    ax.set_ylabel(feat, fontsize=8)
fig.colorbar(sc, ax=axes, label="ABSOLUTE purity", shrink=0.6)
fig.suptitle("TCGA-LUAD primary tumors: TACSTD2 vs immune features", fontsize=12)
fig.savefig(os.path.join(FIG, "scatter_tacstd2_vs_immune.png"), dpi=150)
plt.close(fig)

r2 = res.sort_values("partial_rho_purity_adj")
fig, ax = plt.subplots(figsize=(8, 0.32 * len(r2) + 1.5))
colors = ["#c0392b" if s else "#95a5a6" for s in r2["sig_after_purity_adj_fdr05"]]
ax.barh(r2["feature"], r2["partial_rho_purity_adj"], color=colors)
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("Partial Spearman rho (TACSTD2 vs feature | ABSOLUTE purity)")
ax.set_title("TCGA-LUAD: purity-adjusted association of TACSTD2 with immune features\n"
             "(red = FDR < 0.05 after purity adjustment)", fontsize=10)
ax.tick_params(axis="y", labelsize=7)
fig.savefig(os.path.join(FIG, "bar_partial_rho.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------- provenance ----------------
prov = {
    "claim": "A1 (LUAD-only): TACSTD2 expression is associated with immune features in "
             "TCGA-LUAD primary tumors after adjusting for tumor purity",
    "cohort": "TCGA-LUAD primary tumor samples (-01) with RNA-seq and ABSOLUTE purity",
    "n_expression_primary": int(eg.shape[0]),
    "n_analysis": int(df.shape[0]),
    "tacstd2_vs_purity_spearman": {"rho": float(rho_purity), "p": float(p_purity)},
    "purity_vs_leukocyte_fraction_spearman_rho": float(rho_lf_pur),
    "methods": {
        "unadjusted": "Spearman correlation",
        "adjusted": "partial Spearman (Pearson on rank residuals after regressing out "
                    "ranked ABSOLUTE purity), t-test with n-3 df",
        "multiple_testing": "Benjamini-Hochberg FDR across all tested features, computed "
                            "separately for unadjusted and adjusted p-values",
        "duplicates": "replicate aliquots averaged per patient-level -01 sample",
        "randomness": RNG_NOTE,
    },
    "inputs": {
        "expression": {
            "file": "TCGA.LUAD.HiSeqV2.gz",
            "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
            "desc": "UCSC Xena TCGA-LUAD gene expression, log2(RSEM normalized_count + 1)",
        },
        "purity": {
            "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
            "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
            "desc": "PanCanAtlas ABSOLUTE purity/ploidy calls",
        },
        "leukocyte_fraction": {
            "file": "TCGA_all_leuk_estimate.masked.20170107.tsv",
            "url": "https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
            "desc": "PanImmune (Thorsson 2018) DNA-methylation leukocyte fraction",
        },
        "cibersort": {
            "file": "TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
            "url": "https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a",
            "desc": "PanImmune CIBERSORT relative immune-cell fractions",
        },
        "signatures": {
            "file": "Scores_160_Signatures.tsv.gz",
            "url": "https://api.gdc.cancer.gov/data/80a82092-161d-4615-9d96-e858f113618d",
            "desc": "PanImmune 160 expression signature scores (Wolf et al. subset used)",
        },
    },
}
for k, v in prov["inputs"].items():
    v["md5"] = md5(os.path.join(DATA, v["file"]))
with open(os.path.join(OUT, "provenance.json"), "w") as f:
    json.dump(prov, f, indent=2)

print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print(f"\nTACSTD2 vs purity: rho={rho_purity:.3f}, p={p_purity:.2g}")
print(f"purity vs leukocyte fraction: rho={rho_lf_pur:.3f}")
print("done")

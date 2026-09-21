#!/usr/bin/env python3
"""Bar summary of pre-specified family log2FC. Reads the scored TSV."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
fam = pd.read_csv(ROOT / "results" / "tables" / "family_contrasts.tsv", sep="\t")

order = [
    ("GSE237361_UWB_olaparib_96h", "UWB1.289 ola 96h\nBRCA1-mut n=4"),
    ("GSE237361_OVCAR3_olaparib_96h", "OVCAR3 ola 96h\nBRCA1-WT n=4"),
    ("GSE243208_UWB_olaparib_24h", "UWB1.289 ola 24h\nn=3"),
    ("GSE120500_Brca1def_tumor_olaparib", "Brca1-def tumor\nolaparib n=6"),
    ("GSE285827_ovcar3_talazoparib", "OVCAR3 talazoparib\nn=3"),
    ("GSE285827_caov3_talazoparib", "CAOV3 talazoparib\nn=3"),
    ("GSE285827_ovcar3_veliparib", "OVCAR3 veliparib\nn=3"),
    ("GSE233820_SBC5_PARPi_0Gy", "SBC5 SCLC ola\nn=3"),
]
labels = [lab for _, lab in order]
ids = [i for i, _ in order]

def col(family):
    rows = []
    for cid in ids:
        hit = fam[(fam.contrast_id == cid) & (fam.family == family)]
        rows.append(float(hit.mean_gene_log2fc.iloc[0]) if len(hit) else float("nan"))
    return rows

ifn = col("core_IFN")
cl = col("CLDN4")
tj = col("core_TJ")

fig, ax = plt.subplots(figsize=(10.2, 4.6))
import numpy as np
x = np.arange(len(labels))
w = 0.26
ax.bar(x - w, ifn, width=w, color="#b45309", label="core IFN")
ax.bar(x, cl, width=w, color="#1d4e89", label="CLDN4")
ax.bar(x + w, tj, width=w, color="#4a7c59", label="core TJ, CLDN4 out")
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("mean gene log2FC")
ax.set_title("PARP inhibitor RNA-seq: IFN, CLDN4, tight junction")
ax.legend(frameon=False, ncol=3, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
out = ROOT / "results" / "figures"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "parpi_ifn_cldn4_tj.png", dpi=160)
fig.savefig(out / "parpi_ifn_cldn4_tj.pdf")
print(out / "parpi_ifn_cldn4_tj.png")

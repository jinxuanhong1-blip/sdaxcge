"""Generate KM curves + TACSTD2/CD8A scatter for the open cohorts."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

DERIVED = "results/hunt_oak/data/derived"
OUT = "results/hunt_oak/outputs"

# --- IMvigor210 OS KM by TACSTD2 median split ------------------------------
d = pd.read_csv(f"{DERIVED}/imvigor210_tacstd2_persample.csv")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
kmf = KaplanMeierFitter()
for grp, c in (("low", "#2166ac"), ("high", "#b2182b")):
    g = d[(d.TACSTD2_group == grp)].dropna(subset=["os_time", "os_event"])
    kmf.fit(g.os_time, g.os_event, label=f"TACSTD2 {grp} (n={len(g)})")
    kmf.plot_survival_function(ax=ax[0], color=c, ci_show=True)
ax[0].set_title("IMvigor210 (atezolizumab, mUC)\nOverall survival by TACSTD2")
ax[0].set_xlabel("Months"); ax[0].set_ylabel("Survival probability")
ax[0].text(0.02, 0.05, "logrank p=0.32; HR high vs low 0.88 (0.68-1.13)",
           transform=ax[0].transAxes, fontsize=8)

ax[1].scatter(d.log2_TACSTD2, d.log2_CD8A, s=14, alpha=0.6, color="#4d4d4d")
ax[1].set_xlabel("log2 TACSTD2 (norm)"); ax[1].set_ylabel("log2 CD8A (norm)")
ax[1].set_title("IMvigor210: TACSTD2 vs CD8A\nSpearman rho=-0.23, p=2e-5")
fig.tight_layout(); fig.savefig(f"{OUT}/fig_imvigor210.png", dpi=130); plt.close(fig)

# --- GSE135222 PFS KM ------------------------------------------------------
g135 = pd.read_csv(f"{DERIVED}/gse135222_tacstd2_persample.csv")
fig, ax = plt.subplots(figsize=(5.5, 4.4))
for grp, c in (("low", "#2166ac"), ("high", "#b2182b")):
    gg = g135[g135.TACSTD2_group == grp]
    kmf.fit(gg.pfs_time, gg.pfs_event, label=f"TACSTD2 {grp} (n={len(gg)})")
    kmf.plot_survival_function(ax=ax, color=c, ci_show=True)
ax.set_title("GSE135222 (NSCLC, anti-PD-(L)1)\nPFS by TACSTD2 (n=27)")
ax.set_xlabel("Days"); ax.set_ylabel("PFS probability")
ax.text(0.02, 0.05, "logrank p=0.43; median 59 vs 73 d (high vs low)",
        transform=ax.transAxes, fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_gse135222.png", dpi=130); plt.close(fig)
print("figures written to", OUT)

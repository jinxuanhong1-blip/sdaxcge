#!/usr/bin/env python3
"""Figures + expanded per-gene IFN/APM tables + a one-row-per-contrast verdict.

Re-uses the loaders / gene sets / stats from 03 / lib_*.
"""
import importlib.util
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ld = _load("loaders", "03_loaders.py")
S = _load("sets", "lib_sets.py")
ST = _load("stats_lib", "lib_stats.py")
an = _load("analyze", "05_analyze.py")

SETS = {"human": S.gene_sets("human"), "mouse": S.gene_sets("mouse")}
CORE6 = {"human": S.USER_CORE6, "mouse": S.to_mouse(S.USER_CORE6)}

# display order
CONTRAST_ORDER = [
    "GSE207704_T47D", "GSE207704_MCF7", "GSE22493_SKOV3",
    "GSE50927_Cldn4lungWTvsKO", "GSE50927_VILIwtkohi", "GSE50927_VILIwtkolo",
    "GSE274940_EpH4",
    "GSE334497_4T1", "GSE289287_T47Dxeno",
    "GSE245459_SKOV3", "GSE245459_SKOV3_DDP",
]
SHORT = {
    "GSE207704_T47D": "207704 T47D\nCLDN4 KO",
    "GSE207704_MCF7": "207704 MCF7\nCLDN4 KO",
    "GSE22493_SKOV3": "22493 SKOV3\nCLDN4 KD",
    "GSE50927_Cldn4lungWTvsKO": "50927 lung\nCldn4 KO",
    "GSE50927_VILIwtkohi": "50927 VILI-hi\nCldn4 KO",
    "GSE50927_VILIwtkolo": "50927 VILI-lo\nCldn4 KO",
    "GSE274940_EpH4": "274940 EpH4\nCldn-null",
    "GSE334497_4T1": "334497 4T1\nTrop2 KO",
    "GSE289287_T47Dxeno": "289287 T47D xeno\nTrop2 KO",
    "GSE245459_SKOV3": "245459 SKOV3\nshTACSTD2",
    "GSE245459_SKOV3_DDP": "245459 SKOV3\nshTACSTD2+DDP",
}


def build():
    return {c.key: c for c in an.build_contrasts()}


def lfc_of(c):
    return an.contrast_lfc(c)


def verdict_row(c, lfc, p, set_tab):
    """Honest one-line call for this contrast."""
    pv = float(lfc.get(c.perturb_gene, np.nan))
    kd = "YES" if pv == pv and pv <= -0.5 else ("NO / not confirmed" if pv == pv else "gene absent")
    ifn = set_tab[(set_tab.contrast == c.key) &
                  (set_tab.gene_set == "HALLMARK_IFN_ALPHA")].iloc[0]
    core = CORE6[c.species]
    core_lfc = [float(lfc[g]) for g in core if g in lfc.index]
    n_up = sum(x > 0.25 for x in core_lfc)
    n_dn = sum(x < -0.25 for x in core_lfc)
    n_m = len(core_lfc)

    # decision tree
    if c.key.startswith("GSE274940"):
        call = "NOT A CLDN4 TEST (Cldn4 RNA not down; pan-claudin line)"
        ifn_dir = "null"
    elif kd != "YES":
        call = "KD/KO not confirmed in deposited data — do not interpret IFN"
        ifn_dir = "uninterpretable"
    elif c.gse == "GSE50927":
        if ifn["delta_med"] > 0.1 and ifn["p_up"] < 0.05:
            call = "IFN/APM UP, but n=1 whole lung (not a cancer-cell test)"
            ifn_dir = "up (bulk lung, n=1)"
        else:
            call = "no clear IFN up (n=1 whole lung)"
            ifn_dir = "null / mixed"
    elif c.gse == "GSE207704":
        call = "IFN/APM does NOT go up (collapsed FPKM; no replicate test)"
        ifn_dir = "down / null"
    elif c.gse == "GSE22493":
        call = "does NOT open IFN/APM (array QC poor; control is CLDN4-OE)"
        ifn_dir = "not up"
    elif c.gse == "GSE334497":
        call = "CORE6 flat; hallmark IFN weakly up but permutation NS (bulk tumor)"
        ifn_dir = "weak / NS"
    elif c.gse == "GSE289287":
        call = "IFN/ISG UP (cell-intrinsic xenograft in NRG mice)"
        ifn_dir = "up"
    elif c.key == "GSE245459_SKOV3":
        call = "IFN/APM goes DOWN (opposite of the claim)"
        ifn_dir = "down"
    elif c.key == "GSE245459_SKOV3_DDP":
        call = "ISGs up, HLA-A still down; cisplatin-confounded, not a clean test"
        ifn_dir = "mixed"
    else:
        call = "see numbers"
        ifn_dir = "?"

    return dict(
        contrast=c.key, gse=c.gse, arm=c.arm, species=c.species,
        label=c.label, n_test=c.n_test, n_ref=c.n_ref,
        perturb_gene=c.perturb_gene, perturb_log2FC=round(pv, 3) if pv == pv else np.nan,
        knockdown_confirmed=kd,
        core6_n_measured=n_m, core6_n_up=n_up, core6_n_down=n_dn,
        hallmark_ifna_n=int(ifn["n_set"]),
        hallmark_ifna_delta_med=ifn["delta_med"],
        hallmark_ifna_auc=ifn["auc"],
        hallmark_ifna_p_up=ifn["p_up"],
        hallmark_ifna_perm_p_up=ifn["perm_p_up"],
        ifn_direction=ifn_dir,
        honest_call=call,
        caveats=c.caveats,
    )


def main():
    contrasts = build()
    # re-run set tests from the already-written table, then refresh after rebuild
    os.system(f"python3 {os.path.join(HERE, '05_analyze.py')} > {os.path.join(ROOT, 'logs', 'analyze.log')} 2>&1")
    set_tab = pd.read_csv(os.path.join(ROOT, "tables", "geneset_results.tsv"), sep="\t")
    qc = pd.read_csv(os.path.join(ROOT, "tables", "perturbation_qc.tsv"), sep="\t")
    core = pd.read_csv(os.path.join(ROOT, "tables", "core6_per_gene.tsv"), sep="\t")

    # expanded IFN/APM gene table
    rows = []
    for key in CONTRAST_ORDER:
        c = contrasts[key]
        lfc, p = lfc_of(c)
        members = set(SETS[c.species]["ISG_CORE"]) | set(SETS[c.species]["MHC1_APM"]) | set(CORE6[c.species])
        members.add(c.perturb_gene)
        for g in sorted(members):
            if g not in lfc.index:
                rows.append(dict(contrast=key, arm=c.arm, species=c.species,
                                 gene=g, log2FC=np.nan, p=np.nan, present=0))
            else:
                pv = float(p[g]) if p is not None and g in p.index else np.nan
                rows.append(dict(contrast=key, arm=c.arm, species=c.species,
                                 gene=g, log2FC=round(float(lfc[g]), 3),
                                 p=round(pv, 5) if pv == pv else np.nan, present=1))
    wide = pd.DataFrame(rows)
    wide.to_csv(os.path.join(ROOT, "tables", "ifn_apm_per_gene.tsv"),
                sep="\t", index=False)

    # verdict
    vrows = []
    for key in CONTRAST_ORDER:
        c = contrasts[key]
        lfc, p = lfc_of(c)
        vrows.append(verdict_row(c, lfc, p, set_tab))
    verd = pd.DataFrame(vrows)
    verd.to_csv(os.path.join(ROOT, "tables", "verdict.tsv"), sep="\t", index=False)

    # ---------------- figures ----------------
    figdir = os.path.join(ROOT, "figures")
    os.makedirs(figdir, exist_ok=True)

    # 1. perturbation QC
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    x = np.arange(len(CONTRAST_ORDER))
    vals = [float(verd.set_index("contrast").loc[k, "perturb_log2FC"]) for k in CONTRAST_ORDER]
    cols = ["#2a6f97" if v <= -0.5 else "#c1121f" for v in vals]
    ax.bar(x, vals, color=cols, width=0.72, edgecolor="none")
    ax.axhline(-0.5, color="0.4", ls="--", lw=0.8)
    ax.axhline(0, color="0.2", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[k] for k in CONTRAST_ORDER], fontsize=7)
    ax.set_ylabel("perturbation gene log2FC (KD/KO − control)")
    ax.set_title("Did the deposited data actually knock the gene down?")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "01_perturbation_qc.png"), dpi=160)
    plt.close()

    # 2. CORE6 heatmap (human genes; mouse orthologs averaged when 1:many)
    human_core = S.USER_CORE6
    mat = pd.DataFrame(index=human_core, columns=CONTRAST_ORDER, dtype=float)
    for key in CONTRAST_ORDER:
        c = contrasts[key]
        lfc, _ = lfc_of(c)
        if c.species == "human":
            for g in human_core:
                mat.loc[g, key] = float(lfc[g]) if g in lfc.index else np.nan
        else:
            # map each human gene to its mouse list and take median
            for g in human_core:
                mm = S.to_mouse([g])
                vals = [float(lfc[x]) for x in mm if x in lfc.index]
                mat.loc[g, key] = np.median(vals) if vals else np.nan
    fig, ax = plt.subplots(figsize=(11.2, 3.6))
    data = mat.to_numpy(float)
    vmax = 2.0
    im = ax.imshow(data, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(range(len(human_core)))
    ax.set_yticklabels(human_core, fontsize=9)
    ax.set_xticks(range(len(CONTRAST_ORDER)))
    ax.set_xticklabels([SHORT[k] for k in CONTRAST_ORDER], fontsize=7)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, "·", ha="center", va="center", color="0.5", fontsize=8)
            else:
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(v) > 1.1 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cb.set_label("log2FC (loss − control)")
    ax.set_title("User CORE6 (IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A) — mouse = median of orthologs")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "02_core6_heatmap.png"), dpi=160)
    plt.close()

    # 3. Hallmark IFNα competitive effect
    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    sub = set_tab[set_tab.gene_set == "HALLMARK_IFN_ALPHA"].set_index("contrast")
    dmed = [float(sub.loc[k, "delta_med"]) for k in CONTRAST_ORDER]
    pup = [float(sub.loc[k, "p_up"]) for k in CONTRAST_ORDER]
    cols = []
    for k, d, p in zip(CONTRAST_ORDER, dmed, pup):
        if contrasts[k].arm == "CLDN4":
            cols.append("#1d3557" if d > 0 else "#e63939")
        else:
            cols.append("#457b9d" if d > 0 else "#f4a261")
    ax.bar(x, dmed, color=cols, width=0.72, edgecolor="none")
    ax.axhline(0, color="0.2", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[k] for k in CONTRAST_ORDER], fontsize=7)
    ax.set_ylabel("Hallmark IFNα  median log2FC − background")
    ax.set_title("Competitive IFNα shift (positive = set goes UP vs rest of transcriptome)")
    for i, (d, p) in enumerate(zip(dmed, pup)):
        if p == p and p < 0.05:
            ax.text(i, d + (0.03 if d >= 0 else -0.08), f"p={p:.1e}",
                    ha="center", va="bottom" if d >= 0 else "top", fontsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "03_ifna_competitive.png"), dpi=160)
    plt.close()

    # 4. control sets vs IFN (specificity)
    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    width = 0.22
    for i, (gs, lab, col) in enumerate([
        ("HALLMARK_IFN_ALPHA", "IFNα", "#1d3557"),
        ("MHC1_APM", "MHC-I/APM", "#2a9d8f"),
        ("CTRL_HALLMARK_MYC_V1", "MYC (ctrl)", "#adb5bd"),
    ]):
        sub = set_tab[set_tab.gene_set == gs].set_index("contrast")
        vals = [float(sub.loc[k, "delta_med"]) for k in CONTRAST_ORDER]
        ax.bar(x + (i - 1) * width, vals, width=width, color=col, label=lab, edgecolor="none")
    ax.axhline(0, color="0.2", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[k] for k in CONTRAST_ORDER], fontsize=7)
    ax.set_ylabel("median log2FC − background")
    ax.set_title("IFN / APM vs a MYC control set (is everything moving?)")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "04_ifn_vs_myc_control.png"), dpi=160)
    plt.close()

    print("wrote tables/verdict.tsv, ifn_apm_per_gene.tsv and 4 figures")
    print(verd[["contrast", "knockdown_confirmed", "ifn_direction", "honest_call"]].to_string(index=False))


if __name__ == "__main__":
    main()

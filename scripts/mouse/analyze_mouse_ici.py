#!/usr/bin/env python3
"""Analyze Tacstd2 and Cldn4 vs ICI treatment / immune outcomes in mouse lung cohorts.

Loads processed matrices (downloaded under notes/mouse/data), harmonizes each into
a log-scale expression table + sample->group map, then for every treated cohort computes:
  * per-group mean/SD/n of Tacstd2 & Cldn4
  * treated-vs-control log2 fold change, Welch t-test, Mann-Whitney U (BH-FDR across genes)
  * a per-sample cytotoxic/immune score (z-scored immune module) and Spearman correlation
    of Tacstd2 / Cldn4 with that immune score (the "immune outcome" axis)

All statistics are computed from the real downloaded values; nothing is invented.
Outputs land only under results/mouse/.
"""
import glob
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "notes", "mouse", "data")
RES = os.path.join(ROOT, "results", "mouse")
os.makedirs(RES, exist_ok=True)

GMAP = json.load(open(os.path.join(ROOT, "notes", "mouse", "gene_map.json")))
SYM2ENS = GMAP["symbol_to_ensembl"]
TARGETS = ["Tacstd2", "Cldn4"]
# cytotoxic / T-cell inflamed immune module (the immune-outcome axis)
IMMUNE = ["Cd8a", "Cd8b1", "Gzmb", "Gzmk", "Prf1", "Ifng", "Nkg7",
          "Cxcl9", "Cxcl10", "Cd3e", "Pdcd1", "Ptprc"]


# ---------------------------------------------------------------- loaders
def load_emtab():
    df = pd.read_csv(os.path.join(DATA, "E-MTAB-13704_GEMMS_raw_counts.csv"), index_col=0)
    # counts -> CPM -> log2
    cpm = df / df.sum(axis=0) * 1e6
    logexpr = np.log2(cpm + 1)
    # map ensembl -> symbol for target + immune genes
    ens2sym = {v: k for k, v in SYM2ENS.items() if v}
    keep = [e for e in ens2sym if e in logexpr.index]
    sub = logexpr.loc[keep].copy()
    sub.index = [ens2sym[e] for e in keep]
    # sample -> group from sdrf (fastq prefix R# -> stimulus)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(DATA, "E-MTAB-13704.sdrf.txt"))]
    h = rows[0]
    sub_i = h.index("Comment[SUBMITTED_FILE_NAME]")
    stim_i = h.index("Factor Value[stimulus]")
    gmap = {}
    for r in rows[1:]:
        rid = r[sub_i].split("_")[0]  # e.g. R1
        gmap[rid] = r[stim_i].strip()
    group_map = {c: gmap.get(c) for c in logexpr.columns}
    return dict(name="E-MTAB-13704", expr=sub, full=logexpr, ens=True,
                group_map=group_map, control="vehicle",
                platform="bulk RNA-seq (GEMM lung, log2 CPM)")


def load_gse239485():
    df = pd.read_excel(os.path.join(DATA, "GSE239485_Processed_data.xlsx"), sheet_name="DataNorm")
    df = df.dropna(subset=["gene_name"])
    sample_cols = [c for c in df.columns if c[:2] in ("C_", "D_", "T_")]
    expr = df.groupby("gene_name")[sample_cols].mean()
    grp = {"C": "Control_vehicle", "D": "PolyIC_antiPD1", "T": "PolyIC_antiPD1_antiC5aR1"}
    group_map = {c: grp[c[0]] for c in sample_cols}
    return dict(name="GSE239485", expr=expr, full=expr, ens=False,
                group_map=group_map, control="Control_vehicle",
                platform="bulk RNA-seq (LLC tumor, normalized log2)")


def load_gse297630():
    df = pd.read_excel(os.path.join(DATA, "GSE297630_processed_data.xlsx"),
                       sheet_name="Expression", skiprows=4)
    sig_cols = [c for c in df.columns if "Signal" in c]
    ren = {c: c.split("_")[0] for c in sig_cols}  # C-1..P-3
    df = df.dropna(subset=["Gene Symbol"])
    d = df[["Gene Symbol"] + sig_cols].rename(columns=ren)
    d = d.groupby("Gene Symbol").mean()
    group_map = {c: ("Control" if c.startswith("C") else "antiPD1_tolerant") for c in ren.values()}
    return dict(name="GSE297630", expr=d, full=d, ens=False,
                group_map=group_map, control="Control",
                platform="microarray Clariom S (anti-PD-1 tolerant LLC, log2 RMA)")


def load_gse330658():
    frames = {}
    for f in sorted(glob.glob(os.path.join(DATA, "GSE330658", "*.xlsx"))):
        base = os.path.basename(f)
        sname = base.split("_", 1)[1].rsplit("_S1", 1)[0].replace(".xlsx", "")
        sname = sname.replace("_L001", "")
        d = pd.read_excel(f, sheet_name=0)
        s = d.groupby("Name")["TPM"].mean()
        frames[sname] = np.log2(s + 1)
    expr = pd.DataFrame(frames)

    def grp(col):
        c = col.lower()
        if c.startswith("control"):
            return "Control"
        if c.startswith("ptx-anti-vegf"):
            return "PTX_antiVEGF"
        if c.startswith("ptx"):
            return "PTX"
        if c.startswith("anti-vegf"):
            return "antiVEGF"
        return "other"
    group_map = {c: grp(c) for c in expr.columns}
    return dict(name="GSE330658", expr=expr, full=expr, ens=False,
                group_map=group_map, control="Control",
                platform="bulk RNA-seq (Egfr-mut lung, log2 TPM)")


def load_gse197260():
    df = pd.read_csv(os.path.join(DATA, "GSE197260_RNAseqTPM.txt.gz"), sep="\t", index_col=0)
    expr = np.log2(df + 1)
    group_map = {c: c for c in expr.columns}  # n=1 per condition
    return dict(name="GSE197260", expr=expr, full=expr, ens=False,
                group_map=group_map, control="gef_vehicle_d21",
                platform="bulk RNA-seq TPM (EGFR-TKI +/- aPD1/aVEGFR2, n=1/arm)")


# ---------------------------------------------------------------- stats helpers
def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        val = p[idx] * n / (i + 1)
        prev = min(prev, val)
        ranked[idx] = min(prev, 1.0)
    return ranked


def immune_score(full_expr, ens):
    """Per-sample z-scored mean of available immune-module genes."""
    if ens:
        ids = [SYM2ENS[g] for g in IMMUNE if SYM2ENS.get(g)]
        present = [i for i in ids if i in full_expr.index]
        used = [g for g in IMMUNE if SYM2ENS.get(g) in present]
    else:
        idx_upper = {str(i).upper(): i for i in full_expr.index}
        present = [idx_upper[g.upper()] for g in IMMUNE if g.upper() in idx_upper]
        used = [g for g in IMMUNE if g.upper() in idx_upper]
    if len(present) < 3:
        return None, used
    m = full_expr.loc[present]
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), used


def get_gene_row(ds, gene):
    expr = ds["expr"]
    if gene in expr.index:
        r = expr.loc[gene]
        return r if isinstance(r, pd.Series) else r.iloc[0]
    upper = {str(i).upper(): i for i in expr.index}
    if gene.upper() in upper:
        r = expr.loc[upper[gene.upper()]]
        return r if isinstance(r, pd.Series) else r.iloc[0]
    return None


# ---------------------------------------------------------------- main
def analyze(ds):
    out_group = []
    out_cmp = []
    out_corr = []
    gm = pd.Series(ds["group_map"])
    gm = gm[gm.notna()]
    samples = list(gm.index)
    groups = list(pd.unique(gm.values))
    control = ds["control"]

    score, used_imm = immune_score(ds["full"][samples], ds["ens"])

    for gene in TARGETS:
        row = get_gene_row(ds, gene)
        if row is None:
            out_group.append(dict(dataset=ds["name"], gene=gene, group="NA",
                                  n=0, mean=np.nan, sd=np.nan, note="gene not found"))
            continue
        row = row[samples]
        # per-group summary
        for g in groups:
            vals = row[gm[gm == g].index].astype(float)
            out_group.append(dict(dataset=ds["name"], gene=gene, group=g,
                                  n=int(vals.notna().sum()),
                                  mean=round(float(vals.mean()), 4),
                                  sd=round(float(vals.std(ddof=1)), 4) if vals.notna().sum() > 1 else np.nan))
        # treated vs control comparisons
        if control in groups:
            ctrl = row[gm[gm == control].index].astype(float).dropna()
            raw_p = []
            recs = []
            for g in groups:
                if g == control:
                    continue
                t = row[gm[gm == g].index].astype(float).dropna()
                log2fc = float(t.mean() - ctrl.mean())
                if len(t) > 1 and len(ctrl) > 1:
                    tp = stats.ttest_ind(t, ctrl, equal_var=False).pvalue
                    up = stats.mannwhitneyu(t, ctrl, alternative="two-sided").pvalue
                else:
                    tp, up = np.nan, np.nan
                recs.append(dict(dataset=ds["name"], gene=gene, comparison=f"{g}_vs_{control}",
                                 n_treated=len(t), n_control=len(ctrl),
                                 log2FC=round(log2fc, 4),
                                 welch_p=tp, mwu_p=up))
                raw_p.append(tp)
            valid = [p for p in raw_p if p == p]
            fdr = bh_fdr(raw_p) if any(p == p for p in raw_p) else [np.nan] * len(raw_p)
            for r, f in zip(recs, fdr):
                r["welch_fdr"] = float(f) if f == f else np.nan
                out_cmp.append(r)
        # correlation with immune score
        if score is not None:
            common = row.dropna().index.intersection(score.dropna().index)
            if len(common) >= 4:
                rho, pr = stats.spearmanr(row[common].astype(float), score[common].astype(float))
                pear, pp = stats.pearsonr(row[common].astype(float), score[common].astype(float))
                out_corr.append(dict(dataset=ds["name"], gene=gene, n=len(common),
                                     spearman_rho=round(float(rho), 4), spearman_p=float(pr),
                                     pearson_r=round(float(pear), 4), pearson_p=float(pp),
                                     immune_genes_used=";".join(used_imm)))
    return out_group, out_cmp, out_corr


def main():
    loaders = [load_emtab, load_gse239485, load_gse297630, load_gse330658, load_gse197260]
    all_group, all_cmp, all_corr = [], [], []
    meta_rows = []
    for ld in loaders:
        ds = ld()
        g, c, r = analyze(ds)
        all_group += g
        all_cmp += c
        all_corr += r
        gm = pd.Series(ds["group_map"]).dropna()
        meta_rows.append(dict(dataset=ds["name"], platform=ds["platform"],
                              n_samples=len(gm), groups=";".join(f"{k}({v})" for k, v in gm.value_counts().items()),
                              control=ds["control"]))
        print(f"[done] {ds['name']}: {len(gm)} samples, groups={dict(gm.value_counts())}")

    pd.DataFrame(all_group).to_csv(os.path.join(RES, "gene_group_summary.csv"), index=False)
    pd.DataFrame(all_cmp).to_csv(os.path.join(RES, "treated_vs_control_stats.csv"), index=False)
    pd.DataFrame(all_corr).to_csv(os.path.join(RES, "immune_correlation.csv"), index=False)
    pd.DataFrame(meta_rows).to_csv(os.path.join(RES, "cohort_overview.csv"), index=False)
    print("\nWrote CSVs to", RES)
    print("\n== treated vs control (Tacstd2/Cldn4) ==")
    print(pd.DataFrame(all_cmp).to_string(index=False))
    print("\n== immune correlation ==")
    print(pd.DataFrame(all_corr).to_string(index=False))


if __name__ == "__main__":
    main()
